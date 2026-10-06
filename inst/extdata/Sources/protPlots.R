### Coverage maps + ratio plots for proteins of interest
#
# Version re-written by Claude to improve parallelisation!
# The AI also neatly detected a handful of bugs.
# Correlation plots now run for any workflow with >=2 samples/sample groups,
# independently of the ratio-plot selection.
#
#stopCluster(parClust)
source(parSrc)
nCharLim <- 40L
goAhead <- FALSE
if (prot.list.Cond) {
  setwd(wd)
  evids <- as.integer(unlist(strsplit(PG$`Evidence IDs`[PG$`In list` == "+"], ";")))
  goAhead <- length(evids)
}
if (goAhead) {
  TMP0 <- ev[ev$id %in% evids,]
  if (scrptType == "withReps") {
    rfRoot <- paste0("Mean ", pep.ref[length(pep.ref)])
    smplLvl <- VPAL$values
    removeMBR <- FALSE # Could move to parameters as with noReps script
    smplCol <- "Samples group"
    TMP0[[smplCol]] <- Exp.map[match(TMP0$`Parent sample`, Exp.map$`Parent sample`), VPAL$column]
  }
  if (scrptType == "noReps") {
    rfRoot <- paste0(int.col, " - ")
    smplLvl <- Exp
    smplCol <- "Experiment"
  }
  if (removeMBR) {
    TMP0 <- TMP0[grep("-MATCH$", TMP0$Type, invert = TRUE),]
  }
  TEST0 <- listMelt(strsplit(TMP0$Proteins, ";"), TMP0$id)
  TEST0 <- TEST0[TEST0$value %in% unlist(IDs.list),]
  tmpDB <- db[db$`Protein ID` %in% prot.list,
              c("Protein ID", "Common Name", "Sequence")]
  ##############################################################################
  # Sample-pair correlations for the "Correlation" plots: computed for EVERY
  # possible pair of samples/sample groups, whenever there are >=2 of them,
  # independently of scrptType and of the ratio-plot selection below.
  ##############################################################################
  runRat <- FALSE
  myCombAll <- pepR <- NULL
  corrOK <- length(smplLvl) > 1L
  if (corrOK) {
    corrCombAll <- as.data.frame(gtools::permutations(length(smplLvl), 2L, smplLvl, repeats.allowed = TRUE))
    colnames(corrCombAll) <- c("A", "B")
    corrCombAll <- corrCombAll2 <- corrCombAll[corrCombAll$A != corrCombAll$B,]
    if (scrptType == "withReps") {
      corrCombAll2$A <- cleanNms(corrCombAll2$A)
      corrCombAll2$B <- cleanNms(corrCombAll2$B)
    }
    myCombAll <- do.call(paste, c(corrCombAll, sep = " / "))
    pepR <- apply(corrCombAll, 1L, \(x) { # x <- corrCombAll[1L,]
      paste0("R = ", round(cor(pep[[paste0(rfRoot, x[[1L]])]],
                               pep[[paste0(rfRoot, x[[2L]])]],
                               use = "pairwise.complete.obs"), 3L))
    })
    names(pepR) <- apply(corrCombAll2, 1L, \(x) { paste0(x[[1L]], " (A) vs ", x[[2L]], " (B)") })
  }
  ##############################################################################
  # Ratio plots: user-selected subset of comparisons, noReps workflow only
  # (unchanged from before - this is a narrower, separate thing from correlation)
  ##############################################################################
  if (scrptType == "noReps" && corrOK) {
    opt <- vapply(myCombAll, \(x) { paste(c(x, rep(" ", max(c(1L, 250L-nchar(x))))), collapse = "") }, "")
    slct <- dlg_list(opt, opt, TRUE, "Ratio plots: select comparisons of interest")$res
    m <- match(slct, opt)
    comb <- corrCombAll[m,]
    myComb <- myCombAll[m]
    runRat <- nrow(comb) > 0L
  }
  tmpPep <- pep[, grep(topattern(rfRoot), colnames(pep), value = TRUE)]
  mtchTst <- setNames(lapply(1L:length(prot.names), \(iii) { #iii <- 1L
    p <- prot.names[iii]
    IDs <- IDs.list[[p]]
    if (length(IDs) == 0L) {
      p <- toupper(svDialogs::dlgInput(paste0("Sorry, I could not find protein name \"", IDs, "\" in the database, do you want to provide an alternate name?"), "")$res)
      if (p != "") {
        prot.names[iii] <- p
        IDs <- tmpDB$"Protein ID"[tmpDB$"Common Name" == p]
        if (length(IDs) == 0L) { warning("Really sorry, but I really cannot make sense of this protein name, skipping...") }
      } else { warning("Ok, fine, we'll skip then.") }
    }
    w <- which(TEST0$value %in% IDs)
    return(list(IDs = IDs,
                Found = length(w) > 0L))
  }), prot.names)
  IDs_tmp <- setNames(lapply(mtchTst, `[[`, "IDs"), names(mtchTst))
  IDs.list[names(IDs_tmp)] <- IDs_tmp
  mtchTst <- vapply(mtchTst, `[[`, TRUE, "Found")
  whMtch <- which(mtchTst)
  goAhead <- length(whMtch)
}
if (goAhead) {
  invisible(clusterCall(parClust, \() {
    library(Peptides)
    library(proteoCraft)
    library(magrittr)
    library(ggplot2)
    library(plotly)
    library(htmlwidgets)
    return()
  }))
  ##############################################################################
  # 1. MASTER: slim down and serialize the input data, one small file per unit
  ##############################################################################
  dataDir <- tempfile("protPlots_")
  dir.create(dataDir, recursive = TRUE)
  # Only the evidence columns actually needed downstream
  TMPslim <- TMP0[, c("id", "Sequence", "Modified sequence", "Intensity", "PEP", smplCol)]
  colnames(TMPslim)[ncol(TMPslim)] <- "Sample"
  # One row per (protein, accession)
  pidTab <- do.call(rbind, lapply(whMtch, \(iii) {
    p <- prot.names[iii]
    data.frame(Protein = p, ID = IDs.list[[p]], stringsAsFactors = FALSE)
  }))
  idMtch <- match(pidTab$ID, tmpDB$"Protein ID")
  pidTab$Seq <- tmpDB$Sequence[idMtch]
  pidTab$nm1 <- gsub("\\*", "STAR", gsub("/", "-", gsub(" - $", "", pidTab$ID)))
  nm2 <- gsub("\\*", "STAR", gsub("/", "-", gsub(" - $", "", tmpDB$"Common Name"[idMtch])))
  pidTab$nm <- ifelse(pidTab$nm1 != nm2, paste0(pidTab$nm1, "_", nm2), pidTab$nm1)
  wLong <- which(nchar(pidTab$nm) > nCharLim)
  pidTab$nm[wLong] <- paste0(gsub(" +$", "", substr(pidTab$nm[wLong], 1L, nCharLim-3L)), "...")
  pidTab$dirNm <- gsub(":", "_", gsub("\\.\\.\\.$", "", pidTab$nm))
  pidTab$key <- paste0("pid", seq_len(nrow(pidTab)))
  pidTab$dirInt <- paste0(wd, "/Protein plots/", pidTab$dirNm, "/Coverage/Intensity")
  pidTab$dirPEP <- paste0(wd, "/Protein plots/", pidTab$dirNm, "/Coverage/PEP")
  pidTab$dirCorr <- paste0(wd, "/Protein plots/", pidTab$dirNm, "/Correlation")
  pidTab$dirRatio <- paste0(wd, "/Protein plots/", pidTab$dirNm, "/Ratios")
  for (dr in unlist(pidTab[, c("dirInt", "dirPEP", "dirCorr", "dirRatio")])) {
    if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
  }
  # Evidence subset per protein (shared by all its accessions), written once
  protU <- unique(pidTab$Protein)
  pidTab$protFile <- paste0(dataDir, "/prot", match(pidTab$Protein, protU), ".rds")
  for (p in protU) {
    w <- as.numeric(unique(TEST0$L1[which(TEST0$value %in% IDs.list[[p]])]))
    readr::write_rds(TMPslim[match(w, TMPslim$id),],
                     pidTab$protFile[match(p, pidTab$Protein)],
                     compress = "gz")
  }
  rm(TMPslim)
  # Small config object, the only "global" data the workers need
  plotCfg <- list(dataDir = dataDir, scrptType = scrptType, smplLvl = smplLvl, runRat = runRat,
                  corrOK = corrOK,
                  myComb = if (runRat) { myComb } else { NULL },
                  myCombAll = if (corrOK) { myCombAll } else { NULL },
                  pepR = if (corrOK) { pepR } else { NULL })
  
  ##############################################################################
  # 2. WORKER FUNCTIONS
  ##############################################################################
  # Stage 1 (per protein accession): peptide/sequence matching, per-sample tables,
  # comparison tables. Writes one small RDS per sample (coverage) and up to two
  # per accession (correlation / ratio), and returns the list of plot tasks.
  prepPlotData <- function(job) { #job <- pidJobs[[1L]]
    cfg <- plotCfg
    smplLvl <- cfg$smplLvl
    TMP <- readr::read_rds(job$protFile)
    P <- setNames(job$Seq, job$nm)
    P2 <- unlist(strsplit(gsub("I", "L", P), ""))
    s <- data.frame(Seq = unique(TMP$Sequence))
    s$Matches <- lapply(strsplit(gsub("I", "L", s$Seq), ""), \(x) {
      # All possible matches, including overlapping ones
      rs <- NA
      x <- unlist(x)
      l <- length(x)
      if (l) {
        x <- data.frame(ind = 1L:l, seq = x)
        x$Mtch <- lapply(1L:l, \(y) { which(P2 == x$seq[y]) - x$ind[y] + 1L })
        x <- unlist(x$Mtch)
        if (length(x)) {
          x <- aggregate(x, list(x), length)
          x <- x[x$x == l,]
          if (nrow(x)) { rs <- x$Group.1 }
        }
      }
      return(rs)
    })
    s <- listMelt(s$Matches, s$Seq, c("Matches", "Seq"))
    s <- s[!is.na(s$Matches),]
    if (!nrow(s)) {
      return(list(tasks = NULL, msg = paste0("No peptides identified for protein accession ", job$ID, ".")))
    }
    s2 <- data.frame("Modified sequence" = unique(TMP$"Modified sequence"), check.names = FALSE)
    s2$Sequence <- TMP$Sequence[match(s2$"Modified sequence", TMP$"Modified sequence")]
    s2$Matches <- s$Matches[match(s2$Sequence, s$Seq)]
    # Non-imputed intensities: sum over modified sequence (before log), worst PEP
    tempev <- setNames(lapply(smplLvl, \(exp) {
      wh <- which(TMP$Sample == exp)
      if (!length(wh)) { return(NA) }
      e <- TMP[wh,]
      wh <- which(is.finite(log10(e$Intensity)))
      if (!length(wh)) { return(NA) }
      res <- set_colnames(aggregate(e$Intensity[wh], list(e$"Modified sequence"[wh]), sum),
                          c("Modified sequence", "Intensity"))
      e <- set_colnames(aggregate(e$PEP[wh], list(e$"Modified sequence"[wh]), \(x) {
        x <- x[is.finite(x)]
        if (!length(x)) { return(NA) }
        return(max(x))
      }), c("Modified sequence", "PEP"))
      res$PEP <- e$PEP; rm(e)
      res$"log10(Intensity)" <- log10(res$Intensity)
      res$Intensity <- NULL
      res$Matches <- vapply(s2$Matches[match(res$"Modified sequence", s2$"Modified sequence")],
                            paste, "", collapse = ";")
      return(res)
    }), smplLvl)
    tempev <- tempev[vapply(tempev, \(x) { is.data.frame(x) }, TRUE)]
    if (!length(tempev)) {
      return(list(tasks = NULL, msg = paste0("No usable evidence for protein accession ", job$ID, ".")))
    }
    mxInt <- unlist(sapply(names(tempev), \(exp) { tempev[[exp]]$"log10(Intensity)" }))
    mxInt <- ceiling(max(mxInt[is.finite(mxInt)]))
    mxPEP <- unlist(sapply(names(tempev), \(exp) { -log10(tempev[[exp]]$PEP) }))
    mxPEP <- ceiling(max(mxPEP[is.finite(mxPEP)]))
    # One small file per sample for the coverage plots
    covFiles <- vapply(seq_along(tempev), \(k) {
      f <- file.path(cfg$dataDir, paste0(job$key, "_cov", k, ".rds"))
      readr::write_rds(list(P = P, nm = job$nm, exp = names(tempev)[k], tmp = tempev[[k]],
                            mxInt = mxInt, mxPEP = mxPEP, dirInt = job$dirInt, dirPEP = job$dirPEP),
                       f,
                       compress = "gz")
      f
    }, "")
    tasks <- do.call(rbind, lapply(c("CovInt", "CovHeat", "CovPEP"), \(ty) {
      data.frame(Protein = job$Protein, ID = job$ID, Type = ty, Sample = names(tempev),
                 File = covFiles, stringsAsFactors = FALSE)
    }))
    ############################################################################
    # Sample-pair comparison data.
    # Correlation plots use EVERY pair (cfg$myCombAll), gated only on cfg$corrOK.
    # Ratio plots use only the user-selected pairs (cfg$myComb), gated on cfg$runRat.
    ############################################################################
    if ((length(tempev) > 1L) && (cfg$corrOK || cfg$runRat)) {
      tempAB <- reshape::melt(tempev)
      tempAB$Match <- sapply(strsplit(tempAB$Matches, ";"), as.integer)
      tempAB$Matches <- NULL
      if (inherits(tempAB$Match, "list")) { # Multiple peptide matches in the protein
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
      tempAB1 <- tempAB[tempAB$variable == "log10(Intensity)",]
      tempAB2 <- tempAB[tempAB$variable == "PEP",]
      tempAB1$PEP <- tempAB2$value[match(tempAB1$Dummy, tempAB2$Dummy)]
      tempAB <- tempAB1; rm(tempAB1, tempAB2)
      tempAB$tmp <- as.character(tempAB$Match)
      #
      # Builds temp2/temp3 for a given subset of comparisons (e.g. all pairs, or
      # only the user-selected ones); returns NULL if that subset yields nothing.
      buildCompData <- function(myCombSubset) {
        comb2 <- as.data.frame(gtools::permutations(length(smplLvl), 2L, smplLvl, repeats.allowed = TRUE))
        colnames(comb2) <- c("A", "B")
        myComb2 <- do.call(paste, c(comb2, sep = " / "))
        comb2 <- comb2[myComb2 %in% myCombSubset,]
        if (!nrow(comb2)) { return(NULL) }
        temp2 <- lapply(1L:nrow(comb2), \(j) {
          w1 <- which(tempAB$L1 == comb2[j, 1L])
          w2 <- which(tempAB$L1 == comb2[j, 2L])
          if (!length(c(w1, w2))) { return() }
          s1 <- tempAB[w1,]
          s2 <- tempAB[w2,]
          s1$tmp <- apply(s1[, c("Modified sequence", "tmp")], 1L, paste, collapse = "___")
          s2$tmp <- apply(s2[, c("Modified sequence", "tmp")], 1L, paste, collapse = "___")
          s <- data.frame(tmp = union(s1$tmp, s2$tmp),
                          "log10(LFQ, A)" = NA,
                          "log10(LFQ, B)" = NA,
                          check.names = FALSE)
          w1 <- which(s$tmp %in% s1$tmp)
          w2 <- which(s$tmp %in% s2$tmp)
          m1 <- match(s$tmp[w1], s1$tmp)
          m2 <- match(s$tmp[w2], s2$tmp)
          s[w1, c("Modified sequence", "Match", "log10(LFQ, A)")] <- s1[m1, c("Modified sequence", "Match", "value")]
          s[w2, c("Modified sequence", "Match", "log10(LFQ, B)")] <- s2[m2, c("Modified sequence", "Match", "value")]
          if (cfg$scrptType == "withReps") {
            comb2[[1L]] <- cleanNms(comb2[[1L]]) 
            comb2[[2L]] <- cleanNms(comb2[[2L]]) 
          }
          s$Comparison <- paste0(comb2[j, 1L], " (A) vs ", comb2[j, 2L], " (B)")
          return(s)
        })
        temp2 <- temp2[vapply(temp2, \(x) { is.data.frame(x) }, TRUE)]
        temp2 <- do.call(rbind, temp2)
        if (is.null(temp2) || !nrow(temp2)) { return(NULL) }
        temp3 <- as.data.frame(t(sapply(unique(temp2$Comparison), \(x) {
          x1 <- temp2[temp2$Comparison == unlist(x), c("log10(LFQ, A)", "log10(LFQ, B)")]
          x1 <- x1$"log10(LFQ, B)" - x1$"log10(LFQ, A)"
          return(setNames(c(x,
                            paste0("Median = ", round(median(x1, na.rm = TRUE), 3L)),
                            paste0("S.D. = ", round(sd(x1, na.rm = TRUE), 3L))), c("Comparison", "Median", "SD")))
        })))
        temp3$R <- cfg$pepR[temp3$Comparison]
        temp2[, c("A", "B")] <- Isapply(strsplit(temp2$Comparison, " vs "), unlist)
        temp3[, c("A", "B")] <- Isapply(strsplit(temp3$Comparison, " vs "), unlist)
        temp2 <- temp2[order(temp2$Comparison),]
        temp3 <- temp3[order(temp3$Comparison),]
        temp2$Sequence <- gsub("\\([^\\)]+\\)|_", "", temp2$"Modified sequence")
        temp2$"C-terminal extent" <- temp2$Match + nchar(temp2$Sequence) - 1L
        temp2$"log2(A/B)" <- (temp2$"log10(LFQ, A)" - temp2$"log10(LFQ, B)")/log10(2L)
        temp2$"avg. log10(intensity)" <- apply(temp2[, c("log10(LFQ, A)", "log10(LFQ, B)")], 1L, \(x){
          mean(x[is.finite(x)])
        })
        wInf <- which(!is.finite(temp2$"log2(A/B)"))
        tstInf <- length(wInf)
        tstInfA <- tstInfB <- FALSE
        logFC_rng <- NULL
        if (tstInf) {
          logFC_rng <- range(temp2$"log2(A/B)", na.rm = TRUE)
          wA <- wInf[is.finite(temp2$"log10(LFQ, A)"[wInf])]
          wB <- wInf[is.finite(temp2$"log10(LFQ, B)"[wInf])]
          tstInfA <- length(wA)
          tstInfB <- length(wB)
          if (tstInfA) { temp2$"log2(A/B)"[wA] <- logFC_rng[2L] + 1 }
          if (tstInfB) { temp2$"log2(A/B)"[wB] <- logFC_rng[1L] - 1 }
        }
        temp2$finiTest <- apply(temp2[, c("log10(LFQ, A)", "log10(LFQ, B)") ], 1L, \(x) { sum(is.finite(x)) }) == 2L
        return(list(temp2 = temp2, temp3 = temp3, tstInf = tstInf, tstInfA = tstInfA, tstInfB = tstInfB, logFC_rng = logFC_rng))
      }
      
      # Correlation plot: every sample pair, independent of runRat
      if (cfg$corrOK) {
        cd <- buildCompData(cfg$myCombAll)
        if (!is.null(cd)) {
          f <- file.path(cfg$dataDir, paste0(job$key, "_corr.rds"))
          readr::write_rds(list(P = P, nm = job$nm, nm1 = job$nm1, temp2 = cd$temp2, temp3 = cd$temp3,
                                dirCorr = job$dirCorr),
                           f,
                           compress = "gz")
          tasks <- rbind(tasks, data.frame(Protein = job$Protein, ID = job$ID, Type = "Corr",
                                           Sample = NA_character_, File = f, stringsAsFactors = FALSE))
        }
      }
      # Ratio plot: only the user-selected comparisons
      if (cfg$runRat) {
        rd <- buildCompData(cfg$myComb)
        if (!is.null(rd)) {
          f <- file.path(cfg$dataDir, paste0(job$key, "_ratio.rds"))
          readr::write_rds(list(P = P, nm = job$nm, nm1 = job$nm1, temp2 = rd$temp2, temp3 = rd$temp3, tstInf = rd$tstInf,
                                tstInfA = rd$tstInfA, tstInfB = rd$tstInfB, logFC_rng = rd$logFC_rng, dirRatio = job$dirRatio),
                           f,
                           compress = "gz")
          tasks <- rbind(tasks, data.frame(Protein = job$Protein, ID = job$ID, Type = "Ratio",
                                           Sample = NA_character_, File = f, stringsAsFactors = FALSE))
        }
      }
    }
    return(list(tasks = tasks, msg = NULL))
  }
  # Stage 2 helpers: one function per plot type. Each loads only its own file.
  plotCov <- function(d, type) {
    cfg <- plotCfg
    tmp <- d$tmp
    exp_ <- if (cfg$scrptType == "withReps") { cleanNms(d$exp) } else { d$exp }
    ttl <- paste0("Coverage - ", d$nm, " - ", exp_, ", ",
                  switch(type, Int = "log10(int.)", Heat = "sum log10(int.)", PEP = "-log10(PEP)"))
    sv <- paste0(if (type == "PEP") { d$dirPEP } else { d$dirInt }, "/", sub("^Coverage - ", "Cov. ", ttl))
    res <- switch(type,
                  Int = Coverage(d$P, tmp$"Modified sequence", Mode = "Align2", display = FALSE, scale = 100L,
                                 title = ttl, save = c("svg", "html"), save.path = sv,
                                 intensities = tmp$`log10(Intensity)`, maxInt = d$mxInt),
                  Heat = Coverage(d$P, tmp$"Modified sequence", Mode = "Heat", display = FALSE, scale = 100L,
                                  title = ttl, save = "svg", save.path = sv,
                                  intensities = tmp$`log10(Intensity)`),
                  PEP = Coverage(d$P, tmp$"Modified sequence", Mode = "Align2", display = FALSE, scale = 100L,
                                 title = ttl, save = c("svg", "html"), save.path = sv,
                                 intensities = -log10(tmp$PEP), maxInt = d$mxPEP, colscale = 8L))
    if (type == "Heat") { return(NULL) } # Saved to disk only, as before
    return(res[[1L]][[1L]])
  }
  plotCorr <- function(d) {
    cfg <- plotCfg
    temp2 <- d$temp2
    temp3 <- d$temp3
    temp2Corr <- temp2[temp2$finiTest,]
    if (!nrow(temp2Corr)) { return(NULL) }
    x_min <- min(temp2Corr$"log10(LFQ, A)")
    y_min <- min(temp2Corr$"log10(LFQ, B)")
    y_max <- max(temp2Corr$"log10(LFQ, B)")
    if (cfg$scrptType == "withReps") {
      lvls <- cleanNms(union(temp2Corr$A, temp2Corr$B))
      temp2Corr$A <- factor(cleanNms(temp2Corr$A), lvls)
      temp2Corr$B <- factor(cleanNms(temp2Corr$B), lvls)
    }
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
      ggtitle("Correlation plot", subtitle = d$nm1)
    u <- unique(temp2Corr$Comparison)
    if (length(u) == 1L) {
      v <- unlist(strsplit(sub(" \\(B\\)$", "", u),  " \\(A\\) vs "))
      plot <- plot + xlab(v[1L]) + ylab(v[2L])
    } else {
      plot <- plot + facet_grid(B~A)
    }
    ttl_ <- gsub(":|\\*|\\?|<|>|\\||/", "-", paste0("Correlation plot - ", d$nm))
    suppressMessages({
      ggsave(paste0(d$dirCorr, "/", ttl_, ".svg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    })
    return(NULL)
  }
  plotRatio <- function(d) {
    temp2 <- d$temp2
    P <- d$P; nm <- d$nm; nm1 <- d$nm1
    tstInf <- d$tstInf; tstInfA <- d$tstInfA; tstInfB <- d$tstInfB; logFC_rng <- d$logFC_rng
    temp2$A <- paste0("A = ", gsub(" \\(A\\)$", "", temp2$A))
    temp2$B <- paste0("B = ", gsub(" \\(B\\)$", "", temp2$B))
    seqL <- nchar(P)
    Aext <- seqL
    Bext <- max(temp2$`log2(A/B)`) - min(temp2$`log2(A/B)`) + 1L
    sk <- Aext/(3*Bext)
    aaMW <- data.frame(AA = unlist(strsplit(P, "")))
    aaMW$MW <- cumsum(vapply(aaMW$AA, mw, 1.5))
    mwScl <- ceiling(max(aaMW$MW)/1000)*1000
    intSp <- round((mwScl/5)/1000)*1000
    mwScl <- (0L:(mwScl/intSp))*intSp
    mwScl <- mwScl[mwScl <= max(aaMW$MW)]
    mwScl <- data.frame(Da = mwScl)
    mwScl$AA <- 0L
    mwScl$AA[2L:nrow(mwScl)] <- vapply(mwScl$Da[2L:nrow(mwScl)], \(x) {
      wh1 <- which(aaMW$MW <= x)
      wh2 <- which(aaMW$MW >= x)
      if ((!length(wh1)) || (!length(wh2))) { return(NA) }
      w1 <- max(wh1)
      w2 <- min(wh2)
      x1 <- aaMW$MW[w1]
      x2 <- aaMW$MW[w2]
      (x-x1)*(w2-w1)/(x2-x1) + w1
    }, 1)
    mwScl <- mwScl[!is.na(mwScl$AA),]
    nr <- nrow(mwScl)
    if (nr) {
      mwScl <- rbind(mwScl,
                     data.frame("Da" = mwScl$Da[nr]*2-mwScl$Da[nr-1L],
                                "AA" = mwScl$AA[nr]*2-mwScl$AA[nr-1L]))
      mwScl$kDa <- paste0(round(mwScl$Da/1000), " kDa")
    }
    ySum <- summary(temp2$`log2(A/B)`[is.finite(temp2$`log2(A/B)`)])
    yScl <- ySum["Max."] - ySum["Min."]
    yMin <- ySum["Min."] - 0.1*yScl
    xLim <- c(-10L, max(c(seqL, mwScl$AA)))
    yLim <- c(ySum["Min."]-yScl*0.3, ySum["Max."]+yScl*0.05)
    temp2$`Modified sequence` <- gsub("_", "", temp2$`Modified sequence`)
    myCol <- dfltCol <- "log2(A/B)"
    u <- unique(temp2$Comparison)
    if (length(u) == 1L) {
      v <- unlist(strsplit(sub(" \\(B\\)$", "", u),  " \\(A\\) vs "))
      myCol <- paste0("log2(", v[1L], "/", v[2L], ")")
      if (!myCol %in% colnames(temp2)) {
        temp2[[myCol]] <- temp2[[dfltCol]]
        temp2[[dfltCol]] <- NULL
      }
      temp2$A <- gsub_Rep(" \\(A\\)$", "", temp2$A)
      temp2$B <- gsub_Rep(" \\(B\\)$", "", temp2$B)
    }
    temp2$"mod. seq." <- apply(temp2[, c("Modified sequence", "Match", "C-terminal extent", myCol, "avg. log10(intensity)")],
                               1L, \(x) {
                                 x <- gsub("^ +| +$", "", unlist(x))
                                 paste0(" <i>", x[1L],
                                        "</i>\nstart pos.: <i>", x[2L],
                                        "</i>\nend pos.: <i>", x[3L],
                                        "</i>\n", myCol, "<i>: ", x[4L],
                                        "</i>\navg. log10(int.): <i>", x[5L], "</i>")
                               })
    plot <- ggplot(temp2) +
      annotate(geom = "rect", xmin = 0L, xmax = seqL, ymin = yLim[1L], ymax = yLim[2L],
               fill = "lightblue", alpha = 0.1)
    tstN <- aggregate(temp2$Match, list(temp2$A, temp2$B), \(x) { length(unique(x)) })
    if (min(tstN$x) >= 20L) {
      plot <- plot +
        geom_smooth(data = temp2[temp2$finiTest,], aes(x = (`C-terminal extent`+Match)/2, y = .data[[myCol]]),
                    alpha = 0.1, linewidth = 0.5, method = "loess", formula = y ~ x)
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
      ggtitle("Ratio plot", subtitle = nm1) +
      theme_bw()
    if (length(u) > 1L) { plot <- plot + facet_grid(B~A) }
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
        datInfB$Label <- paste0(sub("^B = ", "", datInfB$B), "-only")
        plot <- plot +
          geom_hline(yintercept = logFC_rng[1L] - 0.5, color = "darkred", linetype = "dashed") +
          geom_text(data = datInfB, aes(label = Label),
                    x = 5L, y = logFC_rng[1L] - 0.115*(logFC_rng[2L] - logFC_rng[1L]),
                    color = "darkred", hjust = 0, fontface = "italic")
      }
    }
    plotLY <- ggplotly(plot, tooltip = "tooltip")
    plotLY <- plotly::config(plotLY, modeBarButtonsToRemove = c("select2d", "lasso2d"))
    w <- which(vapply(plotLY$x$data, \(x) {
      isTRUE(x$mode == "text") && (length(x$text) == 1L) && grepl("-only$", x$text)
    }, TRUE))
    for (i in w) {
      plotLY$x$data[[i]]$text <- paste0("<i>", plotLY$x$data[[i]]$text, "</i>")
    }
    if (length(u) > 1L) {
      w <- which(vapply(1L:length(plotLY$x$data), \(x) { !is.null(plotLY$x$data[[x]]$marker$colorbar) }, TRUE))
      plotLY$x$data[[w]]$marker$colorbar$x <- 1.1
      plot <- plot + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0))
    }
    plotLY$x$layout$xaxis$autorange <- TRUE
    plotLY$x$layout$yaxis$autorange <- TRUE
    plotLY <- htmlwidgets::onRender(plotLY, global_autorange)
    plotLY <- plotly::config(plotLY, modeBarButtonsToRemove = c("select2d", "lasso2d"))
    plotLY <- plotly::plotly_build(plotLY)
    ttl_ <- gsub(":|\\*|\\?|<|>|\\||/", "-", paste0("Ratio plot - ", nm))
    fl <- paste0(d$dirRatio, "/", ttl_, ".html")
    tstSave <- try(htmlwidgets::saveWidget(plotly::partial_bundle(plotLY), fl), silent = TRUE)
    okSave <- !inherits(tstSave, "try-error")
    suppressMessages({
      ggsave(paste0(d$dirRatio, "/", ttl_, ".svg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    })
    return(list(plotLY = plotLY, file = fl, okSave = okSave,
                gg = if (okSave) { NULL } else { plot }))
  }
  runPlotTask <- function(task) { #task <- taskList[[1L]]
    tryCatch({
      d <- readr::read_rds(task$File)
      res <- switch(task$Type,
                    CovInt = plotCov(d, "Int"),
                    CovHeat = plotCov(d, "Heat"),
                    CovPEP = plotCov(d, "PEP"),
                    Corr = plotCorr(d),
                    Ratio = plotRatio(d))
      list(ok = TRUE, Result = res, err = NULL)
    }, error = function(e) { list(ok = FALSE, Result = NULL, err = conditionMessage(e)) })
  }
  for (f in c("prepPlotData", "plotCov", "plotCorr", "plotRatio", "runPlotTask")) {
    fn <- get(f); environment(fn) <- .GlobalEnv; assign(f, fn)
  }
  clusterExport(parClust, c("plotCfg", "prepPlotData", "plotCov", "plotCorr", "plotRatio", "runPlotTask",
                            "Coverage", "listMelt", "global_autorange", "Exp", "topattern",
                            "annot_to_tabl", "poplot", "wd"),
                envir = environment())
  if (scrptType == "withReps") {
    clusterExport(parClust, "cleanNms", envir = environment())
  }
  ##############################################################################
  # 3. STAGE 1 (parallel, per accession): prepare data + build the task table
  ##############################################################################
  pidJobs <- lapply(seq_len(nrow(pidTab)), \(i) { as.list(pidTab[i, , drop = FALSE]) })
  prep <- parLapplyLB(parClust, pidJobs, prepPlotData)
  for (m in unlist(lapply(prep, `[[`, "msg"))) { warning(m) }
  tasks <- do.call(rbind, lapply(prep, `[[`, "tasks"))
  rm(prep, pidJobs)
  if ((!is.null(tasks)) && nrow(tasks)) {
    tasks <- tasks[order(match(tasks$Type, c("Ratio", "Corr", "CovInt", "CovPEP", "CovHeat"))),]
    rownames(tasks) <- NULL
    ############################################################################
    # 4. STAGE 2 (parallel, per plot)
    ############################################################################
    taskList <- lapply(seq_len(nrow(tasks)), \(i) { as.list(tasks[i, , drop = FALSE]) })
    results <- parLapplyLB(parClust, taskList, runPlotTask)
    rm(taskList)
    wErr <- which(!vapply(results, \(x) { isTRUE(x$ok) }, TRUE))
    for (i in wErr) {
      warning(paste0("Plot failed: ", tasks$Type[i], " / ", tasks$Protein[i], " / ", tasks$ID[i],
                     if (!is.na(tasks$Sample[i])) { paste0(" / ", tasks$Sample[i]) } else { "" },
                     ": ", results[[i]]$err))
    }
    ############################################################################
    # 5. Re-assemble covPlots/ratioPlots with the original structure
    ############################################################################
    protNames <- prot.names[mtchTst]
    covPlots <- setNames(lapply(protNames, \(p) { list() }), protNames)
    ratioPlots <- setNames(vector("list", length(protNames)), protNames)
    for (p in protNames) {
      tk <- which(tasks$Protein == p)
      if (!length(tk)) { next }
      tk <- tk[tasks$ID[tk] == tail(unique(tasks$ID[tk]), 1L)]
      for (ty in c("CovInt", "CovPEP")) {
        k <- tk[tasks$Type[tk] == ty]
        if (length(k)) {
          smpl <- tasks$Sample[k]
          nms <- if (scrptType == "noReps") { smpl } else { cleanNms(smpl) }
          covPlots[[p]][[if (ty == "CovInt") { "logInt" } else { "PEP" }]] <-
            setNames(lapply(k, \(i) { results[[i]]$Result }), nms)
        }
      }
      k <- tk[tasks$Type[tk] == "Ratio"]
      if (length(k) && results[[k]]$ok) {
        rr <- results[[k]]$Result
        ratioPlots[p] <- list(rr$plotLY)
        if (WorkFlow == "Band ID") {
          if (rr$okSave) { system(paste("open", shQuote(rr$file))) } else { poplot(rr$gg, 12L, 22L) }
        }
      }
    }
    rm(results)
    if (scrptType == "noReps") {
      ratioPlots_fl %<o% paste0(wd, "/Protein plots/ratioPlots.RDS")
      saveFun(ratioPlots, ratioPlots_fl)
    }
    rm(ratioPlots)
    covPlots_fl %<o% paste0(wd, "/Protein plots/covPlots.RDS")
    saveFun(covPlots, covPlots_fl)
    rm(covPlots)
  }
  unlink(dataDir, recursive = TRUE)
}
dirs <- list.dirs(paste0(wd, "/Protein plots"), recursive = TRUE, full.names = TRUE)
tst <- vapply(dirs, \(dr) { length(list.files(dr, recursive = TRUE)) }, 1L)
dirs2rmv <- dirs[which(tst == 0L)]
for (dr in dirs2rmv) { shell(paste0("RMDIR /S /Q \"", dr, "\""), mustWork = FALSE) } # Unlink doesn't work...
# This code seems to damage the cluster: check the source afterwards!
source(parSrc)
