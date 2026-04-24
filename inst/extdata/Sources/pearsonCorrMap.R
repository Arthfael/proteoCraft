#### Samples Pearson correlation heatmap
if (scrptTypeFull %in% c("Histones", "withReps_PG_and_PTMs")) {
  # (It's not implemented yet for other types of scripts)
  dir <- paste0(wd, "/Pearson correlation map")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  if (scrptTypeFull == "withReps_PG_and_PTMs") {
    if (LocAnalysis) { prtRfRoot %<o% Prot.Expr.Root2 } else { prtRfRoot %<o% Prot.Expr.Root }
    refRoot <- prtRfRoot
    dirlist <- unique(c(dirlist, dir))
    g <- paste0(prtRfRoot, RSA$values)
    w <- which(g %in% colnames(PG))
    g <- g[w]
    allSamples <- cleanNms(RSA$values[w])
    kol <- paste0(prtRfRoot, RSA$values[w])
    myData <- PG[, kol]
    colnames(myData) <- paste0(prtRfRoot, allSamples)
    isLog <- TRUE
  }
  if (scrptTypeFull == "Histones") {
    refRoot <- "Intensity - "
    #allSamples already exists!
    myData <- pep
    isLog <- FALSE
  }
  g <- paste0(refRoot, allSamples)
  corMap <- sapply(seq_along(g), function(x) {
    vapply(seq_along(g), function(y) {
      if (x > y) {
        x <- myData[[g[x]]]
        y <- myData[[g[y]]]
        if (isLog) {
          x <- log10(x)
          y <- log10(y)
        }
        w <- which((is.all.good(x, 2))&(is.all.good(y, 2)))
        return(cor(x[w], y[w], method = "pearson"))
      } else { return(NA) }
    }, 1)
  })
  rownames(corMap) <- colnames(corMap) <- allSamples
  corMap <- reshape2::melt(corMap, measure.vars = allSamples)
  colnames(corMap) <- c("Sample 2", "Sample 1", "Pearson corr.")
  corMap <- corMap[which(!is.na(corMap$"Pearson corr.")),]
  tmp <- myData[, g[2:length(g)]]
  colnames(tmp) <- allSamples[2:length(g)]
  tmp <- reshape::melt(tmp, measure.vars = allSamples[2:length(g)])
  scattrMap <- data.frame("Sample 1" = allSamples[1],
                          "Sample 2" = tmp$variable,
                          "X" = myData[[g[1]]],
                          "Y" = tmp$value,
                          check.names = FALSE)
  if (length(allSamples) > 2) {
    for (i in 2:(length(allSamples)-1)) {
      tmp <- tmp[(which(tmp$variable == allSamples[i+1])[1]):nrow(tmp),]
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
  if (nrow(myData) > 500) {
    plot <- plot + geom_bin_2d(aes(x = X, y = Y), bins = max(c(20, round(nrow(myData)/100))))
  } else {
    plot <- plot + geom_point(aes(x = X, y = Y), size = 0.1)
  }
  plot <- plot +
    scale_fill_viridis() +
    new_scale("fill") +
    geom_tile(data = corMap, aes(fill = `Pearson corr.`, x = 0.5, y = 0.5), width = 1, height = 1) +
    scale_fill_viridis(option = "B") + coord_fixed() +
    facet_grid(`Sample 2`~`Sample 1`) + ggtitle(ttl) +
    theme_minimal() + theme(axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            strip.text.y.right = element_text(angle = 0, size = 5.5),
                            strip.text.x.top = element_text(angle = 45, size = 5.5))
  poplot(plot, 12, 22)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 600)
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 600)
  })
}
