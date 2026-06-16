#### Samples Pearson correlation heatmap
if (scrptTypeFull %in% c("Histones", "withReps_PG_and_PTMs")) {
  # (It's not implemented yet for other types of scripts)
  dir <- paste0(wd, "/Pearson correlation map")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  if (scrptTypeFull == "withReps_PG_and_PTMs") {
    dataType <- "PG"
    prtRfRoot <- if (LocAnalysis) { Prot.Expr.Root2 } else { Prot.Expr.Root }
    prtRfRoot %<o% prtRfRoot
    refRoot <- prtRfRoot
    dirlist <- unique(c(dirlist, dir))
    datNm <- intersect(c("ComBat", "Imputed", #"Filtered",
                         "Original"),
                       names(clustDat[[dataType]]))[1L]
    myData <- clustDat[[dataType]][[datNm]]
    allSamples <- colnames(myData)
    colnames(myData) <- paste0(prtRfRoot, allSamples)
    isLog <- TRUE
    subTtl <- ""
    if (datNm == "ComBat") {
      subTtl <- "batch corrected (ComBat)"
    }
    if (datNm == "Imputed") {
      subTtl <- "Imputed"
    }
  }
  if (scrptTypeFull == "Histones") {
    refRoot <- "Intensity - "
    #allSamples already exists!
    myData <- pep
    isLog <- FALSE
    subTtl<- "Peptide intensities"
  }
  myData <- as.data.frame(myData)
  g <- paste0(refRoot, allSamples)
  corMap <- sapply(seq_along(g), \(x) { #x <- 2L
    vapply(seq_along(g), \(y) { #y <- 1L
      if (x <= y) { return(NA) }
      x <- myData[[g[x]]]
      y <- myData[[g[y]]]
      if (isLog) {
        x <- log10(x)
        y <- log10(y)
      }
      w <- which((is.finite(x))&(is.finite(y)))
      return(cor(x[w], y[w], method = "pearson"))
    }, 1)
  })
  rownames(corMap) <- colnames(corMap) <- allSamples
  corMap <- dfMelt(corMap, ColNames = c("Sample 1", "Pearson corr."))
  corMap$"Sample 2" <- allSamples
  corMap <- corMap[which(!is.na(corMap$"Pearson corr.")),]
  tmp <- myData[, g[2L:length(g)]]
  colnames(tmp) <- allSamples[2L:length(g)]
  tmp <- reshape::melt(tmp, measure.vars = allSamples[2L:length(g)])
  scattrMap <- data.frame("Sample 1" = allSamples[1L],
                          "Sample 2" = tmp$variable,
                          "X" = myData[[g[1L]]],
                          "Y" = tmp$value,
                          check.names = FALSE)
  lAS <- length(allSamples)
  if (lAS > 2L) {
    for (i in 2L:(lAS-1L)) {
      tmp <- tmp[(which(tmp$variable == allSamples[i+1L])[1L]):nrow(tmp),]
      tmp2 <- data.frame("Sample 1" = allSamples[i],
                         "Sample 2" = tmp$variable,
                         "X" = myData[[g[i]]],
                         "Y" = tmp$value,
                         check.names = FALSE)
      scattrMap <- rbind(scattrMap, tmp2)
    }
  }
  scattrMap <- scattrMap[which(!is.na(scattrMap$X)),]
  scattrMap <- scattrMap[which(!is.na(scattrMap$Y)),]
  Mn <- min(scattrMap$X)
  scattrMap$X <- scattrMap$X - Mn
  scattrMap$Y <- scattrMap$Y - Mn
  Mx <- max(scattrMap$X)
  scattrMap$X <- scattrMap$X / Mx
  scattrMap$Y <- scattrMap$Y / Mx
  scattrMap$`Sample 1` <- factor(scattrMap$`Sample 1`, levels = allSamples)
  scattrMap$`Sample 2` <- factor(scattrMap$`Sample 2`, levels = allSamples)
  corMap$`Sample 1` <- factor(corMap$`Sample 1`, levels = allSamples)
  corMap$`Sample 2` <- factor(corMap$`Sample 2`, levels = allSamples)
  ttl <- "Samples Pearson correlation map"
  plot <- ggplot(scattrMap)
  plot <- if (nrow(myData) > 500L) {
    plot + geom_bin_2d(aes(x = X, y = Y), bins = max(c(20L, round(nrow(myData)/100))))
  } else {
    plot + geom_point(aes(x = X, y = Y), size = 0.1)
  }
  plot <- plot +
    scale_fill_viridis() +
    new_scale("fill") +
    geom_tile(data = corMap, aes(fill = `Pearson corr.`, x = 0.5, y = 0.5), width = 1L, height = 1L) +
    scale_fill_viridis(option = "B") + coord_fixed() +
    facet_grid(`Sample 2`~`Sample 1`) + ggtitle(ttl, subtitle = subTtl) +
    theme_minimal() + theme(axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            strip.text.y.right = element_text(angle = 0, size = 5.5),
                            strip.text.x.top = element_text(angle = 45, size = 5.5))
  poplot(plot, 12L, 22L)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 600L)
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 600L)
  })
}
