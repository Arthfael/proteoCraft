#### Samples Pearson correlation heatmap
if (scrptTypeFull %in% c("Histones", "withReps_PG_and_PTMs")) {
  # (It's not implemented yet for other types of scripts)
  if (!exists("clustDat")) {
    if (exists("clustDat_fl") && file.exists(clustDat_fl)) {
      loadFun(clustDat_fl)
    } else {
      Src <- paste0(libPath, "/extdata/Sources/cluster_Heatmap_PrepTop.R")
      #rstudioapi::documentOpen(Src)
      source(Src)
    }
  }
  dir <- paste0(wd, "/Pearson correlation map")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  if (scrptTypeFull == "withReps_PG_and_PTMs") {
    dataType <- "PG"
    prtRfRoot %<o% { if (LocAnalysis) { Prot.Expr.Root2 } else { Prot.Expr.Root } }
    refRoot <- prtRfRoot
    dirlist <- unique(c(dirlist, dir))
    datNm <- base::intersect(c("ComBat", "Imputed", #"Filtered",
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
      w <- which(is.finite(x) & is.finite(y))
      return(cor(x[w], y[w], method = "pearson"))
    }, 1)
  })
  rownames(corMap) <- colnames(corMap) <- allSamples
  corMap <- dfMelt(corMap, ColNames = c("Sample 1", "Pearson corr."))
  corMap$"Sample 2" <- allSamples
  corMap <- corMap[!is.na(corMap$"Pearson corr."),]
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
  scattrMap <- scattrMap[!is.na(scattrMap$X),]
  scattrMap <- scattrMap[!is.na(scattrMap$Y),]
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
  plot <- ggplot2::ggplot(scattrMap)
  plot <- if (nrow(myData) > 500L) {
    plot + ggplot2::geom_bin_2d(ggplot2::aes(x = X, y = Y),
                                bins = max(c(20L, round(nrow(myData)/100))))
  } else {
    plot + ggplot2::geom_point(ggplot2::aes(x = X, y = Y),
                               size = 0.1)
  }
  plot <- plot +
    viridis::scale_fill_viridis() +
    ggnewscale::new_scale("fill") +
    ggplot2::geom_tile(data = corMap, width = 1L, height = 1L,
                       ggplot2::aes(fill = `Pearson corr.`, x = 0.5, y = 0.5)) +
    viridis::scale_fill_viridis(option = "B") + ggplot2::coord_fixed() +
    ggplot2::facet_grid(`Sample 2`~`Sample 1`) +
    ggplot2::ggtitle(ttl, subtitle = subTtl) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   strip.text.y.right = ggplot2::element_text(angle = 0, size = 5.5),
                   strip.text.x.top = ggplot2::element_text(angle = 45, size = 5.5))
  poplot(plot, 12L, 22L)
  suppressMessages({
    ggplot2::ggsave(paste0(dir, "/", ttl, ".svg"), plot, dpi = 600L)
  })
}
