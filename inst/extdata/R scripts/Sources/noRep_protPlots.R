### Coverage maps for proteins of interest
#stopCluster(parClust)
source(parSrc)
nCharLim <- 40
if (prot.list.Cond) {
  setwd(wd)
  evids <- as.integer(unlist(strsplit(PG$`Evidence IDs`[which(PG$`In list` == "+")], ";")))
  if (length(evids)) {
    TMP0 <- ev[which(ev$id %in% evids),]
    if (removeMBR) { TMP0 <- TMP0[grep("-MATCH$", TMP0$Type, invert = TRUE),] }
    TEST0 <- listMelt(strsplit(TMP0$Proteins, ";"), TMP0$id)
    TEST0 <- TEST0[which(TEST0$value %in% unlist(IDs.list)),]
    tmpDB <- db[which(db$`Protein ID` %in% prot.list),
                c("Protein ID", "Common Name", "Sequence")]
    runRat <- FALSE
    l <- length(Exp)
    #if (MakeRatios) {
    if (l > 1) {
      comb <- as.data.frame(gtools::permutations(l, 2, Exp, repeats.allowed = TRUE))
      colnames(comb) <- c("A", "B")
      comb <- comb[which(comb$A != comb$B),]
      myComb <- do.call(paste, c(comb, sep = " / "))
      opt <- vapply(myComb, function(x) { paste(c(x, rep(" ", max(c(1, 250-nchar(x))))), collapse = "") }, "")
      slct <- dlg_list(opt, opt, TRUE, "Ratio plots: select comparisons of interest")$res
      m <- match(slct, opt)
      comb <- comb[m,]
      myComb <- myComb[m]
      runRat <- nrow(comb) > 0
      if (runRat) {
        pepR <- apply(comb, 1, function(x) {
          paste0("R = ", round(cor(pep[[paste0(int.col, " - ", x[[1]])]],
                                   pep[[paste0(int.col, " - ", x[[2]])]]), 3))
        })
        names(pepR) <- apply(comb, 1, function(x) { paste0(x[[1]], " (A) vs ", x[[2]], " (B)") })
      }
    }
    tmpPep <- pep[, grep(topattern(int.col), colnames(pep), value = TRUE)]
    invisible(clusterCall(parClust, function() {
      library(Peptides)
      library(proteoCraft)
      library(magrittr)
      library(plotly)
      library(htmlwidgets)
      return()
    }))
    exports <- list("prot.names", "IDs.list", "TEST0", "tmpDB", "Exp", "TMP0", "nCharLim", "wd", "int.col",
                    "WorkFlow", "dir", "runRat", "Coverage")
    if (runRat) { exports <- append(exports, list("comb", "myComb", "pepR")) }
    clusterExport(parClust, exports, envir = environment())
    #for (prnm in 1:length(prot.names)) { #prnm <- 1
    tstProtPlots <- parLapply(parClust, 1:length(prot.names), function(prnm) { #prnm <- 1
      p <- prot.names[prnm]
      IDs <- IDs.list[[p]]
      if (length(IDs) == 0) {
        p <- toupper(svDialogs::dlgInput(paste0("Sorry, I could not find protein name \"", IDs, "\" in the database, do you want to provide an alternate name?"), "")$res)
        if (p != "") {
          prot.names[prnm] <- p
          IDs <- tmpDB$"Protein ID"[which(tmpDB$"Common Name" == p)]
          if (length(IDs) == 0) { warning("Really sorry, but I really cannot make sense of this protein name, skipping...") }
        } else { warning("Ok, fine, we'll skip then.") }
      }
      if (!length(IDs)) { return() }
      w <- which(TEST0$value %in% IDs)
      if (!length(w)) {
        if (length(IDs) > 1) {
          msg <- paste0("(s) ", paste(IDs[1:(length(IDs)-1)], collapse = ", "), " and ", IDs[length(IDs)])
        } else { msg <- paste0(" ", IDs) }
        warning(paste0("No peptides identified for protein accession", msg, "."))
        return()
      }
      w <- as.numeric(unique(TEST0$L1[w]))
      TMP <- TMP0[match(w, TMP0$id),]
      for (id in IDs) { #id <- IDs[1]
        idMtch <- match(id, tmpDB$"Protein ID")
        nm <- gsub("\\*", "STAR",
                   gsub("/", "-",
                        gsub(" - $", "", paste0(id, "_", tmpDB$"Common Name"[idMtch]))))
        if (nchar(nm) > nCharLim) { nm <- paste0(gsub(" +$", "", substr(nm, 1, nCharLim-3)), "...") }
        nm2 <- gsub(":", "_", gsub("\\.\\.\\.$", "", nm))
        cat("Generating plots for", nm, "\n")
        setwd(wd)
        suppressWarnings({
          drs <- paste0(wd, "/Protein plots/", nm2, c(paste0("/Coverage/", c("Intensity", "PEP")),
                                                      "/Correlation", "/Ratios"))
          for (dr in drs) { if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) } }
        })
        P <- setNames(tmpDB$Sequence[match(id, tmpDB$"Protein ID")], nm)
        P2 <- unlist(strsplit(gsub("I", "L", P), ""))
        s <- data.frame(Seq = unique(TMP$Sequence))
        s$Matches <- lapply(strsplit(gsub("I", "L", s$Seq), ""), function(x) {
          #x <- strsplit(gsub("I", "L", s$Seq), "")[1]
          # Now rewritten to get all possible matches in case a peptide matches multiple times
          # (older versions would've missed overlapping matches)
          rs <- NA
          x <- unlist(x)
          l <- length(x)
          if (l) {
            x <- data.frame(ind = 1:l,
                            seq = x)
            x$Mtch <- lapply(1:l, function(y) {
              which(P2 == x$seq[y]) - x$ind[y] + 1
            })
            x <- unlist(x$Mtch)
            if (length(x)) {
              x <- aggregate(x, list(x), length)
              x <- x[which(x$x == l),]
              if (nrow(x)) { rs <- x$Group.1 }
            }
          }
          return(rs)
        })
        s <- listMelt(s$Matches, s$Seq, c("Matches", "Seq"))
        s <- s[which(!is.na(s$Matches)),]
        if (nrow(s)) {
          s2 <- data.frame("Modified sequence" = unique(TMP$"Modified sequence"), check.names = FALSE)
          s2$Sequence <- TMP$Sequence[match(s2$"Modified sequence", TMP$"Modified sequence")]
          s2$Matches <- s$Matches[match(s2$Sequence, s$Seq)]
          # Important: should use non-imputed intensities
          #(hence why working from ev is good here, as opposed to pep)
          tempev <- setNames(lapply(Exp, function(exp) { #exp <- Exp[1]
            # We are summing over modified sequence for intensities (before log transformation!)
            # and taking the worst (highest) PEP of all PSMs to be conservative.
            wh <- which(TMP$Experiment == exp)
            if (length(wh)) {
              e <- TMP[wh,]
              wh <- which(is.all.good(log10(e$Intensity), 2))
              if (length(wh)) {
                res <- set_colnames(aggregate(e$Intensity[wh], list(e$"Modified sequence"[wh]), sum),
                                    c("Modified sequence", "Intensity"))
                e <- set_colnames(aggregate(e$PEP[wh], list(e$"Modified sequence"[wh]), function(x) {
                  x <- is.all.good(x)
                  if (!length(x)) { x <- NA } else { x <- max(x) }
                  return(x)
                }), c("Modified sequence", "PEP"))
                res$PEP <- e$PEP; rm(e)
                res$"log10(Intensity)" <- log10(res$Intensity)
                res$Intensity <- NULL
                res$Matches <- vapply(s2$Matches[match(res$"Modified sequence", s2$"Modified sequence")],
                                      paste, "", collapse = ";")
              } else { res <- NA }
            } else { res <- NA }
            return(res)
          }), Exp)
          tempev <- tempev[which(vapply(tempev, function(x) { "data.frame" %in% class(x) }, TRUE))]
          if (length(tempev)) {
            mxInt <- ceiling(max(is.all.good(unlist(sapply(names(tempev), function(exp) { #exp <- names(tempev)[1]
              tempev[[exp]]$"log10(Intensity)"
            })))))
            mxPEP <- ceiling(max(is.all.good(unlist(sapply(names(tempev), function(exp) {
              -log10(tempev[[exp]]$PEP)
            })))))
            for (exp in names(tempev)) { #exp <- names(tempev)[1]
              tmp <- tempev[[exp]]
              setwd(paste0(wd, "/Protein plots/", nm2, "/Coverage/Intensity"))
              ttl1a <- paste0("Coverage - ", nm, " - ", exp, ", log10(int.)")
              ttl1b <- paste0("Coverage - ", nm, " - ", exp, ", sum log10(int.)")
              ttl2a <- paste0("Coverage - ", nm, " - ", exp, ", -log10(PEP)")
              Coverage(P, tmp$"Modified sequence", Mode = "Align2", display = FALSE,
                       title = ttl1a, save = c("jpeg", "pdf"), intensities = tmp$`log10(Intensity)`,
                       maxInt = mxInt)
              Coverage(P, tmp$"Modified sequence", Mode = "Heat", display = FALSE,
                       title = ttl1b, save = c("jpeg", "pdf"), intensities = tmp$`log10(Intensity)`)
              setwd(paste0(wd, "/Protein plots/", nm2, "/Coverage/PEP"))
              Coverage(P, tmp$"Modified sequence", Mode = "Align2", display = FALSE,
                       title = ttl2a, save = c("jpeg", "pdf"), intensities = -log10(tmp$PEP),
                       maxInt = mxPEP, colscale = 8)
              # (Doing summed PEP Coverage maps makes no sense)
              # (Note to self: as the number of peptidoforms)
            }
            setwd(wd)
            # Correlation and ratio plots:
            if ((length(tempev) > 1)&&(runRat)) {
              temp1 <- reshape::melt(tempev)
              temp1$Match <- sapply(strsplit(temp1$Matches, ";"), as.integer)
              temp1$Matches <- NULL
              if (class(temp1$Match) == "list") { # Deal with cases with multiple peptide matches in the protein
                temp2 <- listMelt(temp1$Match, 1:nrow(temp1), c("Match", "row"))
                temp2 <- as.data.frame(t(as.data.frame(strsplit(unique(apply(temp2, 1, paste, collapse = "___")), "___"))))
                colnames(temp2) <- c("Match", "row")
                temp2$Match <- as.numeric(temp2$Match)
                temp2$row <- as.numeric(temp2$row)
                temp2[, c("Modified sequence", "variable", "value", "L1")] <- temp1[temp2$row, c("Modified sequence", "variable", "value", "L1")]
                temp1 <- temp2[, c("Modified sequence", "Match", "variable", "value", "L1")]; rm(temp2)
              }
              temp1$variable <- as.character(temp1$variable)
              temp1$Dummy <- apply(temp1[,c("Modified sequence", "L1")], 1, paste, collapse = "---")
              temp11 <- temp1[which(temp1$variable == "log10(Intensity)"),]
              temp12 <- temp1[which(temp1$variable == "PEP"),]
              temp11$PEP <- temp12$value[match(temp11$Dummy, temp12$Dummy)]
              temp1 <- temp11; rm(temp11, temp12)
              comb2 <- as.data.frame(gtools::permutations(length(Exp), 2, Exp, repeats.allowed = TRUE))
              colnames(comb2) <- c("A", "B")
              myComb2 <- do.call(paste, c(comb2, sep = " / "))
              comb2 <- comb2[which(myComb2 %in% myComb),]
              temp2 <- lapply(1:nrow(comb2), function(j) {
                s1 <- temp1[which(temp1$L1 == comb2[j, 1]),]
                s2 <- temp1[which(temp1$L1 == comb2[j, 2]),]
                s2$tmp <- apply(s2[, c("Modified sequence", "Match")], 1, paste, collapse = "___")
                if (!sum(c(nrow(s1), nrow(s2)) > 0) == 2) { return() }
                wtst <- which(s1$"Modified sequence" %in% s2$"Modified sequence")
                if (!length(wtst)) { return() }
                temp3 <- s1[which(s1$"Modified sequence" %in% s2$"Modified sequence"),
                            c("Modified sequence", "Match" , "value")]
                colnames(temp3)[which(colnames(temp3) == "value")] <- "log10(LFQ, A)"
                temp3$"log10(LFQ, B)" <- s2$value[match(apply(temp3[, c("Modified sequence", "Match")],
                                                              1, paste, collapse = "___"),
                                                        s2$tmp)]
                temp3$Comparison <- paste0(comb2[j, 1], " (A) vs ", comb2[j, 2], " (B)")
                return(temp3)
              })
              temp2 <- temp2[which(vapply(temp2, function(x) { "data.frame" %in% class(x) }, TRUE))]
              temp2 <- do.call(rbind, temp2)
              if (nrow(temp2)) {
                test <- apply(temp2[, c("log10(LFQ, A)", "log10(LFQ, B)") ], 1, function(x) { length(is.all.good(x)) }) == 2
                temp2 <- temp2[which(test),]
                temp3 <- as.data.frame(t(sapply(unique(temp2$Comparison), function(x) { #x <- unique(temp2$Comparison)[1]
                  x1 <- temp2[which(temp2$Comparison == unlist(x)), c("log10(LFQ, A)", "log10(LFQ, B)")]
                  x1 <- x1$"log10(LFQ, B)" - x1$"log10(LFQ, A)"
                  return(setNames(c(x, paste0("Median = ", round(median(x1), 3)), paste0("S.D. = ", round(sd(x1), 3))), c("Comparison", "Median", "SD")))
                })))
                temp3$R <- pepR[temp3$Comparison]
                temp2[, c("A", "B")] <- Isapply(strsplit(temp2$Comparison, " vs "), unlist)
                temp3[, c("A", "B")] <- Isapply(strsplit(temp3$Comparison, " vs "), unlist)
                temp2 <- temp2[order(temp2$Comparison),]
                temp3 <- temp3[order(temp3$Comparison),]
                temp2$Sequence <- gsub("\\([^\\)]+\\)|_", "", temp2$"Modified sequence")
                temp2$"C-terminal extent" <- temp2$Match + nchar(temp2$Sequence) - 1
                temp2$"log10(A/B)" <- temp2$`log10(LFQ, A)` - temp2$`log10(LFQ, B)`
                x_min <- min(temp2$`log10(LFQ, A)`)
                y_min <- min(temp2$`log10(LFQ, B)`)
                y_max <- max(temp2$`log10(LFQ, B)`)
                ttl <- paste0("Correlation plot - ", nm)
                plot <- ggplot(temp2) + geom_point(aes(x = `log10(LFQ, A)`, y = `log10(LFQ, B)`, colour = `C-terminal extent`),
                                                   alpha = 1, size = 1, shape = 16) +
                  geom_abline(intercept = 0, slope = 1, colour = "red") + coord_fixed(1) +
                  scale_colour_gradient(low = "green", high = "red") +
                  geom_text(data = temp3, x = x_min, y = y_max - 0.01*(y_max - y_min), aes(label = R), hjust = 0, cex = 2.5) +
                  geom_text(data = temp3, x = x_min, y = y_max - 0.06*(y_max - y_min), aes(label = Median), hjust = 0, cex = 2.5) +
                  geom_text(data = temp3, x = x_min, y = y_max - 0.11*(y_max - y_min), aes(label = SD), hjust = 0, cex = 2.5) +
                  facet_grid(B~A) +
                  theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
                                     strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) +
                  ggtitle(ttl)
                #poplot(plot)
                setwd(paste0(wd, "/Protein plots/", nm2, "/Correlation"))
                suppressMessages({
                  ggsave(paste0(gsub(":|\\*|\\?|<|>|\\||/", "-", ttl), ".jpeg"), plot,
                         dpi = 300, width = 10, height = 10, units = "in")
                  ggsave(paste0(gsub(":|\\*|\\?|<|>|\\||/", "-", ttl), ".pdf"), plot,
                         dpi = 300, width = 10, height = 10, units = "in")
                })
                temp2$A <- paste0("A = ", gsub(" \\(A\\)$", "", temp2$A))
                temp2$B <- paste0("B = ", gsub(" \\(B\\)$", "", temp2$B))
                Aext <- seqL <- nchar(P)
                Bext <- max(temp2$`log10(A/B)`) - min(temp2$`log10(A/B)`)+1
                sk <- Aext/(3*Bext)
                temp2$MW <- vapply(temp2$`C-terminal extent`, function(x) { mw(substr(P, 1, x)) }, 1.5)
                tst <- data.frame(AA = unlist(strsplit(P, "")))
                tst$MW <- vapply(tst$AA, mw, 1.5)
                tst$MW <- cumsum(tst$MW)
                mwScl <- ceiling(max(tst$MW)/1000)*1000
                intSp <- round((mwScl/5)/1000)*1000
                mwScl <- (0:(mwScl/intSp))*intSp
                mwScl <- mwScl[which(mwScl <= max(tst$MW))]
                mwScl <- data.frame(Da = mwScl)
                mwScl$AA <- 0
                mwScl$AA[2:nrow(mwScl)] <- sapply(mwScl$Da[2:nrow(mwScl)], function(x) {
                  #x <- mwScl$Da[2:nrow(mwScl)][1]
                  #x <- mwScl$Da[2:nrow(mwScl)][2]
                  #x <- mwScl$Da[2:nrow(mwScl)][3]
                  #x <- mwScl$Da[2:nrow(mwScl)][4]
                  wh1 <- which(tst$MW <= x)
                  wh2 <- which(tst$MW >= x)
                  if ((length(wh1))&&(length(wh2))) {
                    w1 <- max(wh1)
                    w2 <- min(wh2)
                    x1 <- tst$MW[w1]
                    x2 <- tst$MW[w2]
                    aa <- (x-x1)*(w2-w1)/(x2-x1) + w1
                    return(aa)
                  } else { return(NA) }
                })
                mwScl <- mwScl[which(!is.na(mwScl$AA)),]
                nr <- nrow(mwScl)
                if (nr) {
                  mwScl <- rbind(mwScl,
                                 data.frame("Da" = mwScl$Da[nr]*2-mwScl$Da[nr-1],
                                            "AA" = mwScl$AA[nr]*2-mwScl$AA[nr-1]))
                  mwScl$kDa <- paste0(round(mwScl$Da/1000), " kDa")
                }
                ySum <- summary(is.all.good(temp2$`log10(A/B)`))
                yScl <- ySum["Max."] - ySum["Min."]
                yMin <- ySum["Min."] - 0.1*yScl
                xLim <- c(-10, max(c(seqL, mwScl$AA)))
                yLim <- c(ySum["Min."]-yScl*0.3, ySum["Max."]+yScl*0.05)
                ttl <- paste0("Ratio plot - ", nm)
                tst <- aggregate(temp2$Match, list(temp2$A, temp2$B), function(x) { length(unique(x)) })
                temp2$`Modified sequence` <- gsub("_", "", temp2$`Modified sequence`)
                plot <- ggplot(temp2) +
                  annotate(geom = "rect", xmin = 0, xmax = seqL, ymin = yLim[1], ymax = yLim[2], fill = "lightblue", alpha = 0.1)
                if (min(tst$x) >= 20) {
                  plot <- plot + geom_smooth(aes(x = (`C-terminal extent`+Match)/2, y = `log10(A/B)`), alpha = 0.1, linewidth = 0.5)
                }
                plot <- plot + geom_segment(aes(x = Match, xend = `C-terminal extent`, text = `Modified sequence`,
                                                y = `log10(A/B)`, yend = `log10(A/B)`, colour = `log10(A/B)`), linewidth = 1)
                if (nr) {
                  plot <- plot + geom_point(data = mwScl, aes(x = AA), y = yMin, shape = 17, color = "blue") +
                    geom_text(data = mwScl, aes(x = AA-0.25, label = kDa), y = yMin - 0.025*yScl, hjust = 1, cex = 2, angle = 30, color = "blue")
                }
                plot <- plot + coord_fixed(round(sk*2)) +
                  scale_colour_gradient(low = "green", high = "red") +
                  scale_x_continuous(limits = c(0, seqL), breaks = 50*(1:floor(seqL/50))) +
                  facet_grid(B~A) + ylab("log10(A/B)") +
                  xlim(xLim[1]-5, xLim[2]+5) +
                  ylim(yLim[1], yLim[2]) + xlab("Amino acid") +
                  theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)#,
                                     #strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)
                  ) + ggtitle(ttl)
                setwd(paste0(dir, "/", nm2, "/Ratios"))
                ttl2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", ttl)
                plotLY <- ggplotly(plot)
                tst <- try(saveWidget(partial_bundle(plotLY), paste0(ttl2, ".html")), silent = TRUE)
                if ("try-error" %in% class(tst)) {
                  tst <- try(saveWidget(plotLY, paste0(ttl2, ".html")), silent = TRUE)
                }
                if (WorkFlow == "Band ID") {
                  if (!"try-error" %in% class(tst)) { system(paste0("open \"", ttl2, ".html")) } else {
                    poplot(plot, 12, 20)
                  }
                }
                suppressMessages({
                  ggsave(paste0(ttl2, ".jpeg"), plot, dpi = 600, width = 10, height = 10, units = "in")
                  ggsave(paste0(ttl2, ".pdf"), plot, dpi = 600, width = 10, height = 10, units = "in")
                })
                setwd(wd)
              }
            }
          }
        } else { warning(paste0("No peptides identified for protein accession ", id, ".")) }
      }
      return()
    })
  }
}
# This code seems to damage the cluster: check the source afterwards!
source(parSrc)
