# ROC1 - for excluding PSMs based on GO term annotations
# Only basic in place now, add shiny app
# see https://rpubs.com/Wangzf/pROC
ROCfilt_Pos <- ROCfilt_Neg <- c()
if (length(ROCfilt_GOterms_Pos)) {
  ROCfilt_Pos <- unique(unlist(lapply(ROCfilt_GOterms_Pos, function(x) {
    GO_terms$Offspring[match(x, GO_terms$ID)]
  })))
  if (length(ROCfilt_GOterms_Neg)) {
    ROCfilt_Neg <- unique(unlist(lapply(ROCfilt_GOterms_Neg, function(x) {
      GO_terms$Offspring[match(x, GO_terms$ID)]
    })))
    ROCfilt_Pos <- ROCfilt_Pos[which(!ROCfilt_Pos %in% ROCfilt_Neg)]
    ROCfilt_Neg <- ROCfilt_Neg[which(!ROCfilt_Neg %in% ROCfilt_Pos)]
  }
}
if (length(c(ROCfilt_Pos, ROCfilt_Neg))) {
  pack <- "pROC"
  bioc_req <- unique(c(bioc_req, pack))
  if (!require(pack, character.only = TRUE)) { biocInstall(pack) }
  msg <- "PSMs filtering using ROC analysis"
  ReportCalls <- AddMsg2Report()
  dir <- paste0(wd, "/Workflow control/", evNm, "s/ROC analysis")
  dirlist <- unique(c(dirlist, dir))
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  tmp <- data.table(Proteins = ev$Proteins, Sequence = ev$Sequence)
  tmp <- tmp[, list(Proteins = unique(Proteins)), by = list(Sequence = Sequence)]
  tmp <- as.data.frame(tmp)
  tmp <- listMelt(strsplit(tmp$Proteins, ";"), tmp$Sequence, c("Protein", "Sequence"))
  tmp2 <- data.table(Intensity = ev$Intensity, Sequence = ev$Sequence)
  tmp2 <- tmp2[, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(Sequence = Sequence)]
  tmp2 <- as.data.frame(tmp2)
  tmp2 <- tmp2[which(tmp2$Intensity > 0),]
  tmp2$Predictor <- log10(tmp2$Intensity)
  tmp2 <- tmp2[order(tmp2$Predictor, decreasing = TRUE),]
  #tmp2$Predictor <- nrow(tmp2):1 # Does not have an effect on the result!
  tmp3 <- listMelt(strsplit(db$`GO-ID`, ";"), db$`Protein ID`, c("Term", "Protein"))
  tmp3 <- tmp3[which(tmp3$Term %in% c(ROCfilt_Pos, ROCfilt_Neg)),]
  tst <- rep(FALSE, 2)
  if (length(ROCfilt_Pos)) {
    tmp2$"True Positive" <- FALSE
    p <- tmp3$Protein[which(tmp3$Term %in% ROCfilt_Pos)]
    s <- tmp$Sequence[match(p, tmp$Protein)]
    w <- which(tmp2$Sequence %in% s)
    tmp2$"True Positive"[w] <- TRUE
    if (sum(tmp2$"True Positive") < 10) {
      msg <- "Not enough TRUE positives for ROC analysis (min = 10), skipping!"
      ReportCalls <- AddMsg2Report(Warning = TRUE, Print = FALSE)
    } else {
      rocobj1 <- pROC::roc(tmp2$"True Positive", tmp2$Predictor)
      rocPlot1 <- pROC::ggroc(rocobj1, color = "red") +
        theme_bw()
      tst[1] <- TRUE
    }
  }
  if (length(ROCfilt_Neg)) {
    tmp2$"Not True Negative" <- TRUE
    p <- tmp3$Protein[which(tmp3$Term %in% ROCfilt_Neg)]
    s <- tmp$Sequence[match(p, tmp$Protein)]
    w <- which(tmp2$Sequence %in% s)
    tmp2$"Not True Negative"[w] <- FALSE
    if (sum(!tmp2$"Not True Negative") < 10) {
      msg <- "Not enough TRUE negatives for ROC analysis (min = 10), skipping!"
      ReportCalls <- AddMsg2Report(Warning = TRUE, Print = FALSE)
    } else {
      rocobj2 <- pROC::roc(tmp2$`Not True Negative`, tmp2$Predictor)
      rocPlot2 <- pROC::ggroc(rocobj2, color = "blue") +
        theme_bw()
      tst[2] <- TRUE
    }
  }
  ttl <- paste0("PSMs ROC analysis")
  plotPath <- paste0(dir, "/", ttl, ".png")
  if (sum(tst)) {
    if (sum(tst) == 2) {
      pack <- "gridExtra"
      cran_req <- unique(c(cran_req, pack))
      if (!require(pack, character.only = TRUE)) { pak:pkg_install(pack) }
      rocPlot1 <- as.grob(rocPlot1)
      rocPlot2 <- as.grob(rocPlot2)
      png(plotPath); grid.arrange(rocPlot1, rocPlot2); dev.off()
    } else {
      rocPlot <- get(paste0("rocPlot", which(tst)))
      ggsave(plotPath, rocPlot)
    }
    system(paste0("open \"", plotPath, "\""))
    ROCintPerc %<o% NA
    msg <- "Look at the ROC analysis image then enter the percentage of values to keep..."
    while ((is.na(ROCintPerc))||(ROCintPerc <= 0)||(ROCintPerc > 100)) {
      suppressWarnings({ ROCintPerc <- as.numeric(dlg_input(msg, 100)$res) })
    }
    # Would be much better with, wait for it, a shiny app!
    tmp2Flt <- tmp2[1:round(nrow(tmp2)*ROCintPerc/100),]
    ev <- ev[which(ev$Sequence %in% tmp2Flt$Sequence),]
  }
}
