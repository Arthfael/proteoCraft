# Impute missing peptide intensities
if (Impute) {
  DatAnalysisTxt <- paste0(DatAnalysisTxt,
                           " Missing values were imputed using two different strategies: i) the KNN (K-Nearest Neighbours) method for Missing-At-Random values within sample groups, and ii) the QRICL (Quantile Regression Imputation of Left-Censored data) method for Missing-Not-At-Random values.")
  msg <- "Imputing missing values..."
  ReportCalls <- AddMsg2Report(Space = FALSE)
  #
  dir <- paste0(wd, "/Workflow control/Peptides/Imputation")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  #
  ref <- pep.ref[length(pep.ref)]
  pat <- topattern(ref)
  if ("Imputation" %in% names(pep.ref)) { ref <- pep.ref[match("Imputation", names(pep.ref))-1] }
  kol <- grep(pat, colnames(pep), value = TRUE)
  groups <- Exp.map[match(gsub(pat, "", kol), Exp.map$Ref.Sample.Aggregate), VPAL$column]
  temp <- Data_Impute2(pep[, kol], groups, is.log = FALSE)
  temp2 <- temp$Imputed_data
  #View(pep[, kol]/temp2)
  # Plot data
  plotTmp <- pep[, kol]
  colnames(plotTmp) <- cleanNms(gsub(pat, "", colnames(plotTmp)))
  plotTmp <- melt(plotTmp, id.vars = NULL)
  colnames(plotTmp) <- c("Sample", "Original")
  plotTmp$Original <- log10(plotTmp$Original)
  plotTmp$Imputed <- log10(melt(temp2, id.vars = NULL)$value)
  plotTmp$Sample <- factor(plotTmp$Sample, levels = cleanNms(RSA$values))
  w <- which(!is.all.good(plotTmp$Original, 2))
  mnX <- floor(min(is.all.good(plotTmp$Original)) - 1)
  plotTmp$Original[w] <- mnX
  ttl <- "Effect of imputation"
  plot <- ggplot(plotTmp) + geom_scattermore(aes(x = Original, y = Imputed, color = Sample), pointsize = 2.5) +
    facet_wrap(~Sample) + theme_bw() + ggtitle(ttl, subtitle = "(non-transformed values)") + coord_fixed() +
    geom_vline(xintercept = mnX + 0.5, color = "darkred", linetype = "dashed") +
    geom_text(x = mnX + 0.525, y = mnX + 0.5, label = "<- imputed", size = 3, color = "darkred", hjust = 0)
  # TO DO:
  # Add a way to look at number of missed values across sample group, so we can distinguish MNAR from MAR on the plots!
  poplot(plot, 12, 22)
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150)
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150)
  #
  colnames(temp2) <- gsub(pat, paste0("imput. ", pep.ref["Original"]), colnames(temp2))
  pep[, colnames(temp2)] <- temp2
  pep.ref["Imputation"] <- paste0("imput. ", pep.ref["Original"])
  rm(list = ls()[which(!ls() %in% .obj)])
  Script <- readLines(ScriptPath)
}