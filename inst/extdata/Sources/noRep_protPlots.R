### Coverage maps for proteins of interest
#stopCluster(parClust)
source(parSrc)
nCharLim <- 40L
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
    if (l > 1L) {
      comb <- as.data.frame(gtools::permutations(l, 2L, Exp, repeats.allowed = TRUE))
      colnames(comb) <- c("A", "B")
      comb <- comb[which(comb$A != comb$B),]
      myComb <- do.call(paste, c(comb, sep = " / "))
      opt <- vapply(myComb, \(x) { paste(c(x, rep(" ", max(c(1L, 250L-nchar(x))))), collapse = "") }, "")
      slct <- dlg_list(opt, opt, TRUE, "Ratio plots: select comparisons of interest")$res
      m <- match(slct, opt)
      comb <- comb[m,]
      myComb <- myComb[m]
      runRat <- nrow(comb) > 0L
      if (runRat) {
        pepR <- apply(comb, 1L, \(x) {
          paste0("R = ", round(cor(pep[[paste0(int.col, " - ", x[[1L]])]],
                                   pep[[paste0(int.col, " - ", x[[2L]])]]), 3L))
        })
        names(pepR) <- apply(comb, 1L, \(x) { paste0(x[[1L]], " (A) vs ", x[[2L]], " (B)") })
      }
    }
    tmpPep <- pep[, grep(topattern(int.col), colnames(pep), value = TRUE)]
    invisible(clusterCall(parClust, \() {
      library(Peptides)
      library(proteoCraft)
      library(magrittr)
      library(plotly)
      library(htmlwidgets)
      return()
    }))
    exports <- list("prot.names", "IDs.list", "TEST0", "tmpDB", "Exp", "TMP0", "nCharLim", "wd", "int.col",
                    "WorkFlow", "runRat", "Coverage", "listMelt", "poplot", "annot_to_tabl")
    if (runRat) { exports <- append(exports, list("comb", "myComb", "pepR")) }
    clusterExport(parClust, exports, envir = environment())
    #for (prnm in 1L:length(prot.names)) { #prnm <- 1L
    tmpPlots <- setNames(parLapply(parClust, 1L:length(prot.names), \(prnm) { #prnm <- 1L
      RES <- list()
      p <- prot.names[prnm]
      IDs <- IDs.list[[p]]
      if (length(IDs) == 0L) {
        p <- toupper(svDialogs::dlgInput(paste0("Sorry, I could not find protein name \"", IDs, "\" in the database, do you want to provide an alternate name?"), "")$res)
        if (p != "") {
          prot.names[prnm] <- p
          IDs <- tmpDB$"Protein ID"[which(tmpDB$"Common Name" == p)]
          if (length(IDs) == 0L) { warning("Really sorry, but I really cannot make sense of this protein name, skipping...") }
        } else { warning("Ok, fine, we'll skip then.") }
      }
      if (!length(IDs)) { return() }
      w <- which(TEST0$value %in% IDs)
      if (!length(w)) {
        msg <- if (length(IDs) > 1L) {
          paste0("(s) ", paste(IDs[1L:(length(IDs)-1L)], collapse = ", "), " and ", IDs[length(IDs)])
        } else { paste0(" ", IDs) }
        warning(paste0("No peptides identified for protein accession", msg, "."))
        return()
      }
      w <- as.numeric(unique(TEST0$L1[w]))
      TMP <- TMP0[match(w, TMP0$id),]
      for (id in IDs) { #id <- IDs[1L]
        idMtch <- match(id, tmpDB$"Protein ID")
        nm <- nm1 <- gsub("\\*", "STAR", gsub("/", "-", gsub(" - $", "", id)))
        nm2 <- gsub("\\*", "STAR", gsub("/", "-", gsub(" - $", "", tmpDB$"Common Name"[idMtch])))
        if (nm1 != nm2) { nm <- paste0(nm, "_", nm2) }
        if (nchar(nm) > nCharLim) { nm <- paste0(gsub(" +$", "", substr(nm, 1L, nCharLim-3L)), "...") }
        nm2 <- gsub(":", "_", gsub("\\.\\.\\.$", "", nm))
        cat("Generating plots for", nm, "\n")
        suppressWarnings({
          drs <- paste0(wd, "/Protein plots/", nm2, c(paste0("/Coverage/", c("Intensity", "PEP")),
                                                      "/Correlation",
                                                      "/Ratios"))
          for (dr in drs) { if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) } }
        })
        P <- setNames(tmpDB$Sequence[match(id, tmpDB$"Protein ID")], nm)
        P2 <- unlist(strsplit(gsub("I", "L", P), ""))
        s <- data.frame(Seq = unique(TMP$Sequence))
        s$Matches <- lapply(strsplit(gsub("I", "L", s$Seq), ""), \(x) {
          #x <- strsplit(gsub("I", "L", s$Seq), "")[1L]
          # Now rewritten to get all possible matches in case a peptide matches multiple times
          # (older versions would've missed overlapping matches)
          rs <- NA
          x <- unlist(x)
          l <- length(x)
          if (l) {
            x <- data.frame(ind = 1L:l,
                            seq = x)
            x$Mtch <- lapply(1L:l, \(y) {
              which(P2 == x$seq[y]) - x$ind[y] + 1L
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
        if (!nrow(s)) {
          warning(paste0("No peptides identified for protein accession ", id, "."))
          return()
        }
        s2 <- data.frame("Modified sequence" = unique(TMP$"Modified sequence"), check.names = FALSE)
        s2$Sequence <- TMP$Sequence[match(s2$"Modified sequence", TMP$"Modified sequence")]
        s2$Matches <- s$Matches[match(s2$Sequence, s$Seq)]
        # Important: should use non-imputed intensities
        #(hence why working from ev is good here, as opposed to pep)
        tempev <- setNames(lapply(Exp, \(exp) { #exp <- Exp[1L]
          # We are summing over modified sequence for intensities (before log transformation!)
          # and taking the worst (highest) PEP of all PSMs to be conservative.
          wh <- which(TMP$Experiment == exp)
          if (!length(wh)) { return(NA) }
          e <- TMP[wh,]
          wh <- which(is.finite(log10(e$Intensity)))
          if (!length(wh)) { return(NA) }
          res <- set_colnames(aggregate(e$Intensity[wh], list(e$"Modified sequence"[wh]), sum),
                              c("Modified sequence", "Intensity"))
          e <- set_colnames(aggregate(e$PEP[wh], list(e$"Modified sequence"[wh]), \(x) {
            x <- x[which(is.finite(x))]
            if (!length(x)) { return(NA) }
            return(max(x))
          }), c("Modified sequence", "PEP"))
          res$PEP <- e$PEP; rm(e)
          res$"log10(Intensity)" <- log10(res$Intensity)
          res$Intensity <- NULL
          res$Matches <- vapply(s2$Matches[match(res$"Modified sequence", s2$"Modified sequence")],
                                paste, "", collapse = ";")
          return(res)
        }), Exp)
        tempev <- tempev[which(vapply(tempev, \(x) { is.data.frame(x) }, TRUE))]
        if (!length(tempev)) { return() }
        mxInt <- unlist(sapply(names(tempev), \(exp) { #exp <- names(tempev)[1L]
          tempev[[exp]]$"log10(Intensity)"
        }))
        mxInt <- ceiling(max(mxInt[which(is.finite(mxInt))]))
        mxPEP <- unlist(sapply(names(tempev), \(exp) {
          -log10(tempev[[exp]]$PEP)
        }))
        mxPEP <- ceiling(max(mxInt[which(is.finite(mxPEP))]))
        RES$Coverage <- list()
        RES$Coverage$logInt <- list()
        RES$Coverage$PEP <- list()
        for (exp in names(tempev)) { #exp <- names(tempev)[1L]
          tmp <- tempev[[exp]]
          ttl1a <- paste0("Coverage - ", nm, " - ", exp, ", log10(int.)")
          ttl1b <- paste0("Coverage - ", nm, " - ", exp, ", sum log10(int.)")
          ttl2a <- paste0("Coverage - ", nm, " - ", exp, ", -log10(PEP)")
          RES$Coverage$logInt[[exp]] <- Coverage(P, tmp$"Modified sequence", Mode = "Align2", display = FALSE,
                                                 title = ttl1a, save = c("jpeg", "pdf", "html"), intensities = tmp$`log10(Intensity)`,
                                                 maxInt = mxInt)
          Coverage(P, tmp$"Modified sequence", Mode = "Heat", display = FALSE,
                   title = ttl1b, save = c("jpeg", "pdf"), intensities = tmp$`log10(Intensity)`)
          dr1 <- paste0(wd, "/Protein plots/", nm2, "/Coverage/Intensity")
          for (ext in c("jpeg", "pdf")) {
            fs::file_move(paste0(ttl1a, ".", ext), dr1)
            fs::file_move(paste0(ttl1b, ".", ext), dr1)
          }
          RES$Coverage$PEP[[exp]] <- Coverage(P, tmp$"Modified sequence", Mode = "Align2", display = FALSE,
                                              title = ttl2a, save = c("jpeg", "pdf", "html"), intensities = -log10(tmp$PEP),
                                              maxInt = mxPEP, colscale = 8L)
          dr2 <- paste0(wd, "/Protein plots/", nm2, "/Coverage/PEP")
          for (ext in c("jpeg", "pdf")) {
            fs::file_move(paste0(ttl2a, ".", ext), dr2)
          }
          # (Doing summed PEP Coverage maps makes no sense)
        }
        # Correlation and ratio plots:
        if ((length(tempev) > 1L) && runRat) {
          tempAB <- reshape::melt(tempev)
          tempAB$Match <- sapply(strsplit(tempAB$Matches, ";"), as.integer)
          tempAB$Matches <- NULL
          if (inherits(tempAB$Match, "list")) { # Deal with cases with multiple peptide matches in the protein
            temp2 <- listMelt(tempAB$Match, 1L:nrow(tempAB), c("Match", "row"))
            temp2 <- as.data.frame(t(as.data.frame(strsplit(unique(apply(temp2, 1L, paste, collapse = "___")), "___"))))
            colnames(temp2) <- c("Match", "row")
            temp2$Match <- as.numeric(temp2$Match)
            temp2$row <- as.numeric(temp2$row)
            temp2[, c("Modified sequence", "variable", "value", "L1")] <- tempAB[temp2$row, c("Modified sequence", "variable", "value", "L1")]
            tempAB <- temp2[, c("Modified sequence", "Match", "variable", "value", "L1")]; rm(temp2)
          }
          tempAB$variable <- as.character(tempAB$variable)
          tempAB$Dummy <- apply(tempAB[,c("Modified sequence", "L1")], 1L, paste, collapse = "---")
          tempAB1 <- tempAB[which(tempAB$variable == "log10(Intensity)"),]
          tempAB2 <- tempAB[which(tempAB$variable == "PEP"),]
          tempAB1$PEP <- tempAB2$value[match(tempAB1$Dummy, tempAB2$Dummy)]
          tempAB <- tempAB1; rm(tempAB1, tempAB2)
          comb2 <- as.data.frame(gtools::permutations(length(Exp), 2, Exp, repeats.allowed = TRUE))
          colnames(comb2) <- c("A", "B")
          myComb2 <- do.call(paste, c(comb2, sep = " / "))
          comb2 <- comb2[which(myComb2 %in% myComb),]
          tempAB$tmp <- as.character(tempAB$Match)
          temp2 <- lapply(1L:nrow(comb2), \(j) { #j <- 1L
            w1 <- which(tempAB$L1 == comb2[j, 1L])
            w2 <- which(tempAB$L1 == comb2[j, 2L])
            if (!length(c(w1, w2))) { return() }
            s1 <- tempAB[w1,]
            s2 <- tempAB[w2,]
            s1$tmp <- apply(s1[, c("Modified sequence", "tmp")], 1L, paste, collapse = "___")
            s2$tmp <- apply(s2[, c("Modified sequence", "tmp")], 1L, paste, collapse = "___")
            s <- data.frame(tmp = union(s1$tmp, s2$tmp),
                            `log10(LFQ, A)` = NA,
                            `log10(LFQ, B)` = NA,
                            check.names = FALSE)
            w1 <- which(s$tmp %in% s1$tmp)
            w2 <- which(s$tmp %in% s2$tmp)
            m1 <- match(s$tmp[w1], s1$tmp)
            m2 <- match(s$tmp[w2], s2$tmp)
            s[w1, c("Modified sequence", "Match", "log10(LFQ, A)")] <- s1[m1, c("Modified sequence", "Match", "value")]
            s[w2, c("Modified sequence", "Match", "log10(LFQ, B)")] <- s2[m2, c("Modified sequence", "Match", "value")]
            s$Comparison <- paste0(comb2[j, 1L], " (A) vs ", comb2[j, 2L], " (B)")
            return(s)
          })
          temp2 <- temp2[which(vapply(temp2, \(x) { is.data.frame(x) }, TRUE))]
          temp2 <- do.call(rbind, temp2)
          if (nrow(temp2)) {
            temp3 <- as.data.frame(t(sapply(unique(temp2$Comparison), \(x) { #x <- unique(temp2$Comparison)[1L]
              x1 <- temp2[which(temp2$Comparison == unlist(x)), c("log10(LFQ, A)", "log10(LFQ, B)")]
              x1 <- x1$"log10(LFQ, B)" - x1$"log10(LFQ, A)"
              return(setNames(c(x,
                                paste0("Median = ", round(median(x1), 3L)),
                                paste0("S.D. = ", round(sd(x1), 3L))), c("Comparison", "Median", "SD")))
            })))
            temp3$R <- pepR[temp3$Comparison]
            temp2[, c("A", "B")] <- Isapply(strsplit(temp2$Comparison, " vs "), unlist)
            temp3[, c("A", "B")] <- Isapply(strsplit(temp3$Comparison, " vs "), unlist)
            temp2 <- temp2[order(temp2$Comparison),]
            temp3 <- temp3[order(temp3$Comparison),]
            temp2$Sequence <- gsub("\\([^\\)]+\\)|_", "", temp2$"Modified sequence")
            temp2$"C-terminal extent" <- temp2$Match + nchar(temp2$Sequence) - 1L
            temp2$"log2(A/B)" <- (temp2$`log10(LFQ, A)` - temp2$`log10(LFQ, B)`)/log10(2L)
            temp2$"avg. log10(intensity)" <- apply(temp2[, c("log10(LFQ, A)", "log10(LFQ, B)")], 1L, \(x){
              mean(x[which(is.finite(x))])
            })
            wInf <- which(!is.finite(temp2$"log2(A/B)"))
            tstInf <- length(wInf)
            tstInfA <- tstInfB <- FALSE
            if (tstInf) {
              logFC_rng <- range(temp2$"log2(A/B)", na.rm = TRUE)
              wA <- wInf[which(is.finite(temp2$`log10(LFQ, A)`[wInf]))]
              wB <- wInf[which(is.finite(temp2$`log10(LFQ, B)`[wInf]))]
              tstInfA <- length(wA)
              tstInfB <- length(wB)
              if (tstInfA) {
                temp2$"log2(A/B)"[wA] <- logFC_rng[2L] + 1
              }
              if (tstInfB) {
                temp2$"log2(A/B)"[wB] <- logFC_rng[1L] - 1
              }
            }
            ttl <- paste0("Correlation plot - ", nm)
            temp2$finiTest <- apply(temp2[, c("log10(LFQ, A)", "log10(LFQ, B)") ], 1L, \(x) { sum(is.finite(x)) }) == 2L
            temp2Corr <- temp2[which(temp2$finiTest),]
            x_min <- min(temp2Corr$`log10(LFQ, A)`)
            y_min <- min(temp2Corr$`log10(LFQ, B)`)
            y_max <- max(temp2Corr$`log10(LFQ, B)`)
            plot <- ggplot(temp2Corr) +
              geom_point(aes(x = `log10(LFQ, A)`, y = `log10(LFQ, B)`, colour = `C-terminal extent`),
                         alpha = 1, size = 1L, shape = 16L) +
              geom_abline(intercept = 0, slope = 1, colour = "red") + coord_fixed(1L) +
              scale_colour_gradient(low = "green", high = "red") +
              geom_text(data = temp3, x = x_min, y = y_max - 0.01*(y_max - y_min), aes(label = R), hjust = 0, cex = 2.5) +
              geom_text(data = temp3, x = x_min, y = y_max - 0.06*(y_max - y_min), aes(label = Median), hjust = 0, cex = 2.5) +
              geom_text(data = temp3, x = x_min, y = y_max - 0.11*(y_max - y_min), aes(label = SD), hjust = 0, cex = 2.5) +
              theme_bw() + theme(strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
                                 strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 0.5)) +
              ggtitle("Correlation plot",
                      subtitle = nm1)
            u <- unique(temp2Corr$Comparison)
            if (length(u) == 1L) {
              v <- unlist(strsplit(sub(" \\(B\\)$", "", u),  " \\(A\\) vs "))
              plot <- plot + xlab(v[1L]) + ylab(v[2L])
            } else {
              plot <- plot + facet_grid(B~A)
            }
            #poplot(plot)
            ttl_ <- gsub(":|\\*|\\?|<|>|\\||/", "-", ttl)
            suppressMessages({
              ggsave(paste0(drs[3L], "/", ttl_, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
              ggsave(paste0(drs[3L], "/", ttl_, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
            })
            temp2$A <- paste0("A = ", gsub(" \\(A\\)$", "", temp2$A))
            temp2$B <- paste0("B = ", gsub(" \\(B\\)$", "", temp2$B))
            Aext <- seqL <- nchar(P)
            Bext <- max(temp2$`log2(A/B)`) - min(temp2$`log2(A/B)`)+1L
            sk <- Aext/(3*Bext)
            temp2$MW <- vapply(temp2$`C-terminal extent`, \(x) { mw(substr(P, 1L, x)) }, 1.5)
            tst <- data.frame(AA = unlist(strsplit(P, "")))
            tst$MW <- vapply(tst$AA, mw, 1.5)
            tst$MW <- cumsum(tst$MW)
            mwScl <- ceiling(max(tst$MW)/1000)*1000
            intSp <- round((mwScl/5)/1000)*1000
            mwScl <- (0L:(mwScl/intSp))*intSp
            mwScl <- mwScl[which(mwScl <= max(tst$MW))]
            mwScl <- data.frame(Da = mwScl)
            mwScl$AA <- 0L
            mwScl$AA[2L:nrow(mwScl)] <- sapply(mwScl$Da[2L:nrow(mwScl)], \(x) {
              #x <- mwScl$Da[2L:nrow(mwScl)][1L]
              #x <- mwScl$Da[2L:nrow(mwScl)][2L]
              #x <- mwScl$Da[2L:nrow(mwScl)][3L]
              #x <- mwScl$Da[2L:nrow(mwScl)][4L]
              wh1 <- which(tst$MW <= x)
              wh2 <- which(tst$MW >= x)
              if ((!length(wh1)) || (!length(wh2))) { return(NA) }
              w1 <- max(wh1)
              w2 <- min(wh2)
              x1 <- tst$MW[w1]
              x2 <- tst$MW[w2]
              aa <- (x-x1)*(w2-w1)/(x2-x1) + w1
              return(aa)
            })
            mwScl <- mwScl[which(!is.na(mwScl$AA)),]
            nr <- nrow(mwScl)
            if (nr) {
              mwScl <- rbind(mwScl,
                             data.frame("Da" = mwScl$Da[nr]*2-mwScl$Da[nr-1L],
                                        "AA" = mwScl$AA[nr]*2-mwScl$AA[nr-1L]))
              mwScl$kDa <- paste0(round(mwScl$Da/1000), " kDa")
            }
            ySum <- summary(temp2$`log2(A/B)`[which(is.finite(temp2$`log2(A/B)`))])
            yScl <- ySum["Max."] - ySum["Min."]
            yMin <- ySum["Min."] - 0.1*yScl
            xLim <- c(-10L, max(c(seqL, mwScl$AA)))
            yLim <- c(ySum["Min."]-yScl*0.3, ySum["Max."]+yScl*0.05)
            ttl <- paste0("Ratio plot - ", nm)
            temp2$`Modified sequence` <- gsub("_", "", temp2$`Modified sequence`)
            myCol <- dfltCol <- "log2(A/B)"
            u <- unique(temp2$Comparison)
            if (length(u) == 1L) {
              v <- unlist(strsplit(sub(" \\(B\\)$", "", u),  " \\(A\\) vs "))
              myCol <- paste0("log2(", v[1L], "/", v[2L], ")")
              if (!myCol %in% colnames(temp2)) { # For if we rerun during testing
                temp2[[myCol]] <- temp2[[dfltCol]]
                temp2[[dfltCol]] <- NULL
              }
              temp2$A <- gsub_Rep(" \\(A\\)$", "", temp2$A)
              temp2$B <- gsub_Rep(" \\(B\\)$", "", temp2$B)
            }
            temp2$"mod. seq." <- apply(temp2[, c("Modified sequence", "Match", "C-terminal extent", myCol, "avg. log10(intensity)")],
                                       1L, \(x) {
                                         x <- gsub("^ +| +$", "", unlist(x)) # because we are mixing characters and numbers of different lengths
                                         paste0(" <i>", x[1L],
                                                "</i>\nstart pos.: <i>", x[2L],
                                                "</i>\nend pos.: <i>", x[3L],
                                                "</i>\n", myCol, "<i>: ", x[4L],
                                                "</i>\navg. log10(int.): <i>", x[5L], "</i>")
                                       })
            plot <- ggplot(temp2) +
              annotate(geom = "rect", xmin = 0L, xmax = seqL, ymin = yLim[1L], ymax = yLim[2L],
                       fill = "lightblue", alpha = 0.1)
            tst <- aggregate(temp2$Match, list(temp2$A, temp2$B), \(x) { length(unique(x)) })
            if (min(tst$x) >= 20L) {
              plot <- plot +
                geom_smooth(data = temp2[which(temp2$finiTest),], aes(x = (`C-terminal extent`+Match)/2, y = .data[[myCol]]),
                            alpha = 0.1, linewidth = 0.5,
                            method = "loess", formula = y ~ x)
            }
            plot <- plot +
              geom_segment(aes(x = Match, xend = `C-terminal extent`, tooltip = `mod. seq.`,
                               y = .data[[myCol]], yend = .data[[myCol]], colour = `avg. log10(intensity)`),
                           linewidth = 1L)
            if (nr) {
              plot <- plot +
                geom_point(data = mwScl, aes(x = AA), y = yMin-0.1, shape = 17L, color = "blue") +
                geom_text(data = mwScl, aes(x = AA-0.25, label = kDa), y = yMin - 0.1*yScl, hjust = 1,
                          cex = 4, angle = 30, color = "blue")
            }
            plot <- plot + coord_fixed(round(sk*2)) +
              scale_colour_gradient(low = "green", high = "red") +
              scale_x_continuous(limits = c(0L, seqL), breaks = 50L*(1L:floor(seqL/50))) +
              xlim(xLim[1L]-5L, xLim[2L]+5L) +
              ylim(yLim[1L], yLim[2L]) +
              xlab("Position") +
              ggtitle("Ratio plot",
                      subtitle = nm1) +
              theme_bw()
            if (length(u) > 1L) {
              plot <- plot + facet_grid(B~A)
            }
            if (tstInf) {
              datInfB <- datInfA <- aggregate(temp2[, c("A", "B")], list(temp2$Comparison), unique)
              if (tstInfA) {
                datInfA$Label <- paste0(sub("^A = ", "", datInfA$A), "-only")
                plot <- plot +
                  geom_hline(yintercept = logFC_rng[2L] + 0.5, color = "darkred", linetype = "dashed") +
                  geom_text(data = datInfA, aes(label = Label),
                            x = 5L, y = logFC_rng[2L] + 0.115*(logFC_rng[2L] - logFC_rng[1L]),
                            color = "darkred", hjust = 0, fontface = "italic")
              }
              if (tstInfB) {
                datInfB$Label <- paste0(sub("^B = ", "", datInfA$B), "-only")
                plot <- plot +
                  geom_hline(yintercept = logFC_rng[1L] - 0.5, color = "darkred", linetype = "dashed") +
                  geom_text(data = datInfB, aes(label = Label),
                            x = 5L, y = logFC_rng[1L] - 0.115*(logFC_rng[2L] - logFC_rng[1L]),
                            color = "darkred", hjust = 0, fontface = "italic")
              }
            }
            plotLY <- ggplotly(plot, tooltip = "tooltip")
            w <- which(vapply(plotLY$x$data, \(x) {
              (x$mode == "text") && (length(x$text) == 1L) && grepl("-only$", x$text)
            }, TRUE))
            for (i in w) {
              plotLY$x$data[[i]]$text <- paste0("<i>", plotLY$x$data[[i]]$text, "</i>")
            }
            if (length(u) > 1L) {
              # Move legend to the right to avoid overlap with y facet label
              w <- which(vapply(1L:length(plotLY$x$data), \(x) { !is.null(plotLY$x$data[[x]]$marker$colorbar) }, TRUE))
              plotLY$x$data[[w]]$marker$colorbar$x <- 1.1
              # Rotate y facet label for ggplot version
              plot <- plot +
                theme(#strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
                  strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0))
            }
            #poplot(plot)
            setwd(drs[4L])
            ttl_ <- gsub(":|\\*|\\?|<|>|\\||/", "-", ttl)
            fl <- paste0(ttl_, ".html")
            tst <- try(saveWidget(partial_bundle(plotLY), fl), silent = TRUE)
            if (inherits(tst, "try-error")) {
              tst <- try(saveWidget(plotLY, fl), silent = TRUE)
            }
            fs::file_move(fl, drs[4L])
            if (WorkFlow == "Band ID") {
              if (!inherits(tst, "try-error")) { system(paste0("open \"", ttl_, ".html")) } else {
                poplot(plot, 12L, 22L)
              }
            }
            suppressMessages({
              ggsave(paste0(drs[4L], "/", ttl_, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
              ggsave(paste0(drs[4L], "/", ttl_, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
            })
            setwd(wd)
            RES$Ratio_plot <- plotLY
          }
        }
      }
      return(RES)
    }), prot.names)
    ratioPlots %<o% setNames(lapply(prot.names, \(nm) { tmpPlots[[nm]]$Ratio_plot }), prot.names)
    covPlots %<o% setNames(lapply(prot.names, \(nm) {
      smpls1 <- names(tmpPlots[[nm]]$Coverage$logInt)
      smpls2 <- names(tmpPlots[[nm]]$Coverage$PEP)
      list(logInt = setNames(lapply(smpls1, \(smpl) { tmpPlots[[nm]]$Coverage$logInt[[smpl]][[1L]][[1L]] }), smpls1),
           PEP = setNames(lapply(smpls2, \(smpl) { tmpPlots[[nm]]$Coverage$PEP[[smpl]][[1L]][[1L]] }), smpls2))
    }), prot.names)
    #ratioPlots[[1L]]
  }
}
# This code seems to damage the cluster: check the source afterwards!
source(parSrc)
