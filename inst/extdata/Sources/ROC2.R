#### ROC analysis
if (length(ROC_GOterms)) {
  library(ggplot2)
  msg <- "ROC analysis"
  ReportCalls <- AddMsg2Report()
  ROC_GOterms <- unique(unlist(c(ROC_GOterms, GO_terms$Offspring[match(ROC_GOterms, GO_terms$ID)])))
  PG$"ROC - True Positive" <- FALSE
  PG$"ROC - True Positive"[grsep2(ROC_GOterms, PG$`GO-ID`)] <- TRUE
  if (sum(PG$"True Positive") < 10) {
    msg <- "Not enough TRUE positives for ROC analysis (min = 10), skipping!"
    ReportCalls <- AddMsg2Report(Warning = TRUE, Print = FALSE)
  } else {
    dir <- paste0(wd, "/ROC analysis")
    dirlist <- unique(c(dirlist, dir))
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    pack <- "pROC"
    bioc_req <- unique(c(bioc_req, pack))
    biocInstall(pack)
    w1 <- which(paste0(pvalue.col[which(pvalue.use)], VPAL$values) %in% colnames(PG))
    if (!length(w1)) { stop("Where are my P-values?!?!%&$+*!*") } else {
      for (grp in VPAL$values[w1]) { #grp <- VPAL$values[w1[1]]
        grp2 <- cleanNms(grp)
        pkol <- paste0(pvalue.col[which(pvalue.use)], grp)
        w2 <- which(is.all.good(PG[[pkol]], 2))
        PG$"ROC - Predictor" <- 10^(-PG[[pkol]])
        rocobj <- pROC::roc(PG[w2, ], "ROC - True Positive", "ROC - Predictor")
        ttl <- paste0("ROC analysis - ", grp2)
        plot <- pROC::ggroc(rocobj, color = "red") +
          geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dotted") +
          theme_bw() + ggtitle(ttl)
        poplot(plot)
        suppressMessages({
          ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300)
          ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)
        })
        ReportCalls <- AddPlot2Report(Title = gsub(": ?", " - ", ttl))
      }
    }
  }
  PG$"ROC - Predictor" <- NULL
  PG$"ROC - True Positive" <- NULL
}
