#### Dimensionality reduction plots
if (!exists("dimRedPlotLy")) { dimRedPlotLy <- list() }
dimRedPlotLy %<o% dimRedPlotLy
#
Exp.map$RSA_cleaned <- cleanNms(Exp.map$Ref.Sample.Aggregate)
#
# Data to plot
# ------------
# Note on missing values: currently we do not impute but filter by number of valid values.
# This may change.
if (dataType == "PG") {
  # clustPrep <- sum(c(exists("clustDat"),
  #                    exists("clustDatImp"))) == 2L
  Src <- paste0(libPath, "/extdata/R scripts/Sources/cluster_Heatmap_Prep.R") # Run this to generate imputed data which we use for the PCAs
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
  kol <- cleanNms(RSA$values)
  kol <- intersect(kol, colnames(clustDat))
  datMatch <- match(row.names(clustDat), PG$Label)
  dimRedDat <- clustDat
  filt <- which((apply(dimRedDat[, kol], 1L, \(x) { length(is.all.good(x)) }) == ncol(dimRedDat))
                &((is.na(PG$`Potential contaminant`[datMatch]))|(PG$`Potential contaminant`[datMatch] != "+")))
}
if (dataType == "modPeptides") {
  kol <- paste0(ptms.ref[length(ptms.ref)], RSA$values)
  dimRedDat <- ptmpep[, kol]
  colnames(dimRedDat) <- kol <- cleanNms(RSA$values)
  rownames(dimRedDat) <- ptmpep$Name
  datMatch <- match(row.names(dimRedDat), ptmpep$Name) 
  filt <- which((apply(dimRedDat[, kol], 1L, \(x) { length(is.all.good(x)) }) == ncol(dimRedDat))
                &((is.na(ptmpep$`Potential contaminant`[datMatch]))|(ptmpep$`Potential contaminant`[datMatch] != "+")))
}
#
# Plots
# -----
if ((length(filt) > 2L)&&(length(kol) > 2L)) {
  dimRedDat <- dimRedDat[filt, kol]
  # Normalizing properly is crucial:
  datMatch <- match(row.names(dimRedDat), PG$Label)
  # We are now (line below) normalizing by average abundance, just in case, however the effect seems minimal:
  if (dataType == "PG") {
    avgVal <- PG$"Av. log10 abundance"[datMatch]
  }
  if (dataType == "modPeptides") {
    avgVal <- rowMeans(dimRedDat, na.rm = TRUE)
  }
  dimRedDat <- sweep(dimRedDat, 1L, avgVal, "-")
  dimRedDat <- dimRedDat + rnorm(length(unlist(dimRedDat)), 0, 10L^-9L) # To avoid constant/zero columns, add a small random error
  #
  # PCA plots, by sample
  if (dataType == "PG") {
    msg <- "PCA plot, by sample"
    dir <- paste0(wd, "/Dimensionality red. plots/PCA")
    ttl <- "PCA plot - Samples (PG-level)"
  }
  if (dataType == "modPeptides") {
    msg <- paste0(ptm, " PCA plot, by sample")
    dir <- paste0(modDirs[1], "/PCA")
    ttl <- paste0(ptm, " PCA plot - Samples")
  }
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- union(dirlist, dir)
  ReportCalls <- AddMsg2Report(Space = FALSE)
  pc <- prcomp(t(dimRedDat), scale. = TRUE)
  scores <- as.data.frame(pc$x)
  pv <- round(100L*(pc$sdev)^2L / sum(pc$sdev^2L), 0L)
  pv <- pv[which(pv > 5L)]
  pv_ <- paste0("Components: ", paste(vapply(seq_along(pv), \(x) {
    paste0("PC", x, ": ", pv[x], "%")
  }, ""), collapse = ", "))
  scores$Sample <- rownames(scores)
  scores[, RSA$names] <- Isapply(strsplit(scores$Sample, "___"), unlist)
  m <- match(scores$Sample, Exp.map$RSA_cleaned)
  scores$Group <- cleanNms(Exp.map[m, VPAL$column])
  scores$Ratios.Groups <- cleanNms(Exp.map[m, RG$column])
  scores[, RSA$names] <- Exp.map[m, RSA$names]
  if (exists("Tim")) {
    scores$"Time.point" <- as.numeric(scores[[Aggregates[[which(names(Aggregates) == "Tim")]]]])
  }
  form <- if (exists("Tim")) {
    if (length(Exp) > 1L) { "Experiment ~ Time.point" } else { "~Time.point" }
  } else { "~Experiment" }
  if ("PC2" %in% colnames(scores)) {
    xLab <- paste0("PC1 = ", pv[1L], "%")
    yLab <- paste0("PC2 = ", pv[2L], "%")
    plot <- ggplot(scores, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = Group)) +
      ggpubr::stat_conf_ellipse(aes(fill = Group),
                                alpha = 0.2, geom = "polygon", show.legend = FALSE) +
      scale_color_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() +
      xlab(xLab) + ylab(yLab) +
      geom_hline(yintercept = 0, colour = "black", alpha = 0.5) +
      geom_vline(xintercept = 0, colour = "black", alpha = 0.5) +
      ggtitle(ttl#, subtitle = pv_
      ) +
      geom_text_repel(aes(label = Sample), size = 2.5, show.legend = FALSE)
    plot <- if (substr(form, 1L, 1L) == "~") { plot + facet_wrap(form) } else { plot + facet_grid(form) }
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    })
    #
    scores$"Samples group" <- factor(scores$Group)
    Symb <- rep(c("circle", "diamond", "square", "cross", "x"), max(as.numeric(Rep)))[1L:max(as.numeric(Rep))]             
    Symb <- Symb[as.numeric(scores$Replicate)]
    if ("PC3" %in% colnames(scores)) {
      plot_lyPCAProt <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3,
                                color = ~`Samples group`, colors = "viridis",
                                text = ~Sample, type = "scatter3d", mode = "markers",
                                symbol = I(Symb))
      plot_lyPCAProt <- add_trace(plot_lyPCAProt, scores, x = ~PC1, y = ~PC2, z = ~PC3,
                                  type = "scatter3d", mode = "text", showlegend = FALSE)
    } else {
      plot_lyPCAProt <- plot_ly(scores, x = ~PC1, y = ~PC2,
                                color = ~`Samples group`, colors = "viridis",
                                text = ~Sample, type = "scatter", mode = "markers",
                                symbol = I(Symb))
      plot_lyPCAProt <- add_trace(plot_lyPCAProt, scores, x = ~PC1, y = ~PC2,
                                  type = "scatter", mode = "text", showlegend = FALSE)
    }
    plot_lyPCAProt <- layout(plot_lyPCAProt, title = ttl)
    setwd(dir)
    saveWidget(plot_lyPCAProt, paste0(dir, "/", ttl, ".html"))
    setwd(wd)
    dimRedPlotLy[["Samples PCA"]] <- plot_lyPCAProt
    system(paste0("open \"", dir, "/", ttl, ".html"))
    # NB: There is currently no way to create a 3D, faceted plot in plotly for R that I know of) 
    ReportCalls <- AddPlot2Report()
  } else { warning("PCA failed, investigate!") }
  #
  #
  if (dataType == "PG") {
    # PCA, t-SNE and UMAP plots, by protein group
    msg <- "PCA plots, by protein group"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    dir <- paste0(wd, "/Dimensionality red. plots/PCA")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- union(dirlist, dir)
    pc <- prcomp(dimRedDat, scale. = TRUE)
    scores <- as.data.frame(pc$x)
    pv <- round(100L*(pc$sdev)^2L / sum(pc$sdev^2L), 0L)
    pv <- pv[which(pv >= 5L)]
    pv <- paste0("Components: ", paste(vapply(seq_along(pv), \(x) {
      paste0("PC", x, ": ", pv[x], "%")
    }, ""), collapse = ", "))
    datMatch <- match(rownames(scores), PG$Label)
    scores[, c("Protein group", "Av. log10 abundance")] <- PG[datMatch, c(Param$Plot.labels, "Av. log10 abundance")]
    scores$Range <- PG$"Rel. av. log10 abundance"[datMatch] # Useful to check if distribution correlates with expression.
    #
    # We will always plot subcellular localisation markers if we have annotations... and overlay regulation on top of that
    nullVal <- " "
    scores$Classifier <- nullVal
    ClassNm <- nullVal
    if (Annotate) {
      datMatch <- match(row.names(dimRedDat), PG$Label)
      #
      if (!exists("SubCellMark")) {
        # Get or generate subcellular localisation markers - for later use
        Src <- paste0(libPath, "/extdata/R scripts/Sources/SubCellMark.R")
        #rstudioapi::documentOpen(Src)
        source(Src, local = FALSE)
      }
      tst3 <- aggregate(SubCellMark, list(SubCellMark), length)
      tmp <- strsplit(PG$`Leading protein IDs`[datMatch], ";")
      SubCellMark2 <- listMelt(strsplit(PG$`Leading protein IDs`[datMatch], ";"), datMatch, c("ID", "row"))
      SubCellMark2 <- SubCellMark2[which(SubCellMark2$ID %in% names(SubCellMark)),]
      SubCellMark2$Comp <- SubCellMark[match(SubCellMark2$ID, names(SubCellMark))]
      SubCellMark2$Label <- PG$Label[SubCellMark2$row]
      tst <- aggregate(SubCellMark2$row, list(SubCellMark2$Comp), list)
      colnames(tst) <- c("Comp", "row")
      ClassNm <- paste(tst$Comp, collapse = " / ")
      tst <- setNames(tst$row, tst$Comp)
      for (comp in names(tst)) {
        w1 <- which(rownames(scores) %in% PG$Label[unlist(tst[[comp]])])
        scores$Classifier[w1] <- comp
      }
      SubCellMark2 <- aggregate(SubCellMark2$Comp, list(SubCellMark2$Label), unique)
      SubCellMark2 <- SubCellMark2[which(lengths(SubCellMark2$x) == 1L),]
      SubCellMark2$x <- unlist(SubCellMark2$x)
      SubCellMark2 <- setNames(SubCellMark2$x, SubCellMark2$Group.1)
    }
    compVal <- unique(scores$Classifier)
    compVal <- sort(compVal[which(compVal != nullVal)])
    allVal <- c(nullVal, compVal)
    #aggregate(scores$Classifier, list(scores$Classifier), length)
    if (!LocAnalysis) { # -> Class based on regulation, supplemented with compartments
      g <- grep("^Regulated - ", colnames(PG), value = TRUE)
      allVal[1L] <- nullVal <- paste(rep("", length(g)), collapse = " / ")
      w <- which(scores$Classifier == " ")
      scores$Classifier[w] <- nullVal
      datMatch <- match(row.names(dimRedDat), PG$Label)
      if (length(g) <= 6L) {
        ClassNm <- cleanNms(gsub("^Regulated - ", "", g))
        tmp <- do.call(cbind, (lapply(g, \(x) {
          gsub_Rep("^Anti-specific: .*", "down",
                   gsub_Rep("^Specific: .*", "up",
                            gsub_Rep(", FDR = .+|^(non significant)|(too small FC)$", "", PG[datMatch, x])))
        })))
        tmp <- set_colnames(as.data.frame(tmp), ClassNm)
        ClassNm <- paste(ClassNm, collapse = " / ")
        tmp <- do.call(paste, c(tmp, sep = " / "))
        gcl <- grep("^ (/ )+$", tmp, invert = TRUE)
        gcl <- gcl[which(tmp[gcl] != nullVal)]
        scores$Classifier[gcl] <- tmp[gcl]
      } else {
        ClassNm <- "Regulated"
        tmp <- do.call(cbind, (lapply(g, \(x) {
          gsub_Rep("^Anti-specific: .*", "down",
                   gsub_Rep("^Specific: .*", "up",
                            gsub_Rep(", FDR = .+|^(non significant)|(too small FC)$", "", PG[datMatch, x])))
        })))
        tmp <- apply(tmp, 1L, \(x) {
          x <- x[which(x != "")]
          x <- if (length(x)) { "regulated" } else { "" }
        })
        gcl <- which(tmp == "regulated")
        scores$Classifier[gcl] <- tmp[gcl]
      }
      regVal <- unique(tmp[gcl])
      regVal <- regVal[which(regVal != nullVal)]
      allVal <- c(allVal, regVal)
    }
    #aggregate(scores$Classifier, list(scores$Classifier), length)
    scores$Classifier <- factor(scores$Classifier, levels = allVal)
    myColors <- setNames(c("grey", rainbow(length(allVal)-1L)), allVal)
    colScale <- scale_colour_manual(name = "colour", values = myColors)
    scores <- rbind(scores[which(scores$Classifier == nullVal),],
                    scores[which(scores$Classifier %in% compVal),],
                    scores[which(scores$Classifier %in% regVal),])
    wBoring <- which(scores$Classifier == nullVal)
    wComp <- which(scores$Classifier %in% compVal)
    wReg <- which(scores$Classifier %in% regVal)
    scores$Size <- 1L
    scores$Size[wComp] <- 1.2
    scores$Size[wReg] <- 1.6
    ttl <- "PCA plot - Protein groups (PG-level)"
    plot <- ggplot(scores[wBoring,], aes(x = PC1, y = PC2, colour = Classifier, size = Size)) +
      geom_scattermore(shape = 16L)
    if (length(wComp)) {
      plot <- plot + geom_point(data = scores[wComp,], shape = 17L)
    }
    if (length(wReg)) {
      plot <- plot + geom_point(data = scores[wReg,], shape = 15L) +
        geom_text_repel(data = scores[wReg,], aes(label = `Protein group`), show.legend = FALSE)
    }
    plot <- plot + coord_fixed() + colScale +
      geom_hline(yintercept = 0, colour = "black", alpha = 0.5) +
      geom_vline(xintercept = 0, colour = "black", alpha = 0.5) +
      scale_alpha_identity() + scale_size_identity() +
      ggtitle(ttl, subtitle = pv) + theme_bw() +
      guides(alpha = "none", size = "none", colour = guide_legend(title = gsub("/", "/\n", ClassNm)))
    #poplot(plot, 12, 22)
    ReportCalls <- AddPlot2Report()
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    })
    #
    # Plotly
    #   NB: I could never have improved these plotly plots so fast without some invaluable help from chatGPT
    #   Guilty as charged...
    PL_colScale <- myColors
    Categories <- c(" ", "compartment marker", "regulated")
    scores$Category <- " "
    # PL_colScale <- setNames("#808080", " ")
    if (length(wComp)) {
      #   PL_comp <- sort(unique(SubCellMark2))
      #   PL_compScale <- setNames(viridis::plasma(length(PL_comp)), PL_comp)
      #   PL_colScale[names(PL_compScale)] <- PL_compScale
      scores$Category[wComp] <- "compartment marker"
    }
    if (length(wReg)) {
      #   PL_reg <- regVal
      #   PL_regScale <- setNames(rainbow(length(PL_reg)), PL_reg)
      #   PL_colScale[names(PL_regScale)] <- PL_regScale
      scores$Category[wReg] <- "regulated"
    }
    val2Cat <- aggregate(scores$Category, list(scores$Classifier), unique)
    val2Cat <- setNames(val2Cat$x, val2Cat$Group.1)
    catsEyes <- aggregate(scores$Size, list(scores$Category), unique)
    catsEyes <- setNames(catsEyes$x, catsEyes$Group.1)
    # PL_symbolScale <- PL_colScale
    # PL_symbolScale[which(names(PL_symbolScale) == " ")] <- "circle"
    # PL_symbolScale[which(names(PL_symbolScale) %in% compVal)] <- "cross"
    # PL_symbolScale[which(names(PL_symbolScale) %in% regVal)] <- "square"
    PL_symbolScale <- setNames(c("circle", "cross", "square"), Categories)
    # Default camera
    myCam <- list(eye = list(x = 1.6,
                             y = 1.6,
                             z = 1.2))
    #
    tst3D <- ("PC3" %in% colnames(scores))
    base_args <- list(x = ~PC1,
                      y = ~PC2,
                      mode = "markers",
                      type = "scatter",
                      showlegend = TRUE,
                      text = ~paste0("<b>Protein group:</b> ", `Protein group`,
                                     "\n<b>", c("Class", "Compartment")[LocAnalysis+1],
                                     ":</b> ", Classifier))
    if (tst3D) {
      base_args$type <- "scatter3d"
      base_args$z <- ~PC3
    }
    plot_lyPCAProt2 <- plot_ly()
    for (val in allVal) { #val <- allVal[1L]
      myCat <- val2Cat[val]
      sz <- catsEyes[myCat]*1.5
      if (val == nullVal) { subDat <- scores } else {
        subDat <- scores[which(scores$Classifier == val),]
        sz <- sz*3L
      }
      args <- base_args
      args$name <- val
      args$marker <- list(size = sz,
                          color = PL_colScale[val],
                          symbol = PL_symbolScale[myCat])
      plot_lyPCAProt2 <- do.call(add_trace, c(list(plot_lyPCAProt2, data = subDat), args))
    }
    # Dummy layer to avoid random re-zooms
    xr <- range(scores$PC1)
    yr <- range(scores$PC2)
    zr <- if (tst3D) { range(scores$PC3) } else { NULL }
    plot_lyPCAProt2 <- add_trace(plot_lyPCAProt2,
                                 x = xr, y = yr, z = zr,
                                 type = base_args$type,
                                 mode = "markers",
                                 marker = list(size = 0.0001, color = "rgba(0,0,0,0)"),
                                 name = "_bbox_",
                                 showlegend = FALSE,
                                 visible = TRUE)
    #
    plot_lyPCAProt2 <- layout(plot_lyPCAProt2, title = ttl,
                              uirevision = TRUE)
    dimRedPlotLy[["PCA"]] <- plot_lyPCAProt2
    setwd(dir)
    saveWidget(plot_lyPCAProt2, paste0(dir, "/", ttl, ".html"))
    setwd(wd)
    #system(paste0("open \"", dir, "/", ttl, ".html"))
    # NB: There is currently no way to create a 3D, faceted plot in plotly for R that I know of) 
    #
    msg <- "t-SNE plots, by protein group"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    cran_req <- unique(c(cran_req, "Rtsne"))
    if (!require("Rtsne", quietly = TRUE)) { install.packages("Rtsne") }
    require(Rtsne)
    tsne <- try(Rtsne(dimRedDat, dims = 3L, perplexity = 30L, verbose = TRUE, max_iter = 500L), silent = TRUE)
    dir <- paste0(wd, "/Dimensionality red. plots/t-SNE")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- union(dirlist, dir)
    if (!"try-error" %in% class(tsne)) {
      scores2 <- as.data.frame(tsne$Y)[, 1L:3L]
      rownames(scores2) <- rownames(dimRedDat)
      colnames(scores2) <- c("t-SNE Y1", "t-SNE Y2", "t-SNE Y3")
      kol2 <- c("Protein group", "Av. log10 abundance", "Classifier", "Size")
      m21 <- match(rownames(scores2), rownames(scores))
      scores2[, kol2] <- scores[m21, kol2]
      wBoring2 <- which(scores2$Classifier == nullVal)
      wComp2 <- which(scores2$Classifier %in% compVal)
      wReg2 <- which(scores2$Classifier %in% regVal)
      ttl2 <- "t-SNE plot - Protein groups (PG-level)"
      plot <- ggplot(scores2[wBoring2,], aes(x = `t-SNE Y1`, y = `t-SNE Y2`, colour = Classifier, size = Size)) +
        geom_scattermore(shape = 16L)
      if (length(wComp2)) {
        plot <- plot + geom_point(data = scores2[wComp2,], shape = 17L)
      }
      if (length(wReg2)) {
        plot <- plot + geom_point(data = scores2[wReg2,], shape = 15L) +
          geom_text_repel(data = scores2[wReg2,], aes(label = `Protein group`))
      }
      plot <- plot +  coord_fixed() + colScale +
        geom_hline(yintercept = 0, colour = "black", alpha = 0.5) +
        geom_vline(xintercept = 0, colour = "black", alpha = 0.5) +
        scale_alpha_identity() + scale_size_identity() +
        ggtitle(ttl2, subtitle = pv) + theme_bw() +
        guides(alpha = "none", size = "none", colour = guide_legend(title = gsub("/", "/\n", ClassNm)))
      #poplot(plot, 12, 22)
      ReportCalls <- AddPlot2Report()
      suppressMessages({
        ggsave(paste0(dir, "/", ttl2, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
        ggsave(paste0(dir, "/", ttl2, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
      })
      #
      # Plotly
      tst3D <- ("t-SNE Y3" %in% colnames(scores2))
      base_args2 <- base_args
      base_args2$x <- ~`t-SNE Y1`
      base_args2$y <- ~`t-SNE Y2`
      if (tst3D) {
        base_args2$type <- "scatter3d"
        base_args2$z <- ~`t-SNE Y3`
      } else {
        base_args2$z <- NULL
      }
      plot_lytSNE <- plot_ly()
      for (val in allVal) { #val <- allVal[1L]
        myCat <- val2Cat[val]
        sz <- catsEyes[myCat]*1.5
        if (val == nullVal) { subDat <- scores2 } else {
          subDat <- scores2[which(scores2$Classifier == val),]
          sz <- sz*3L
        }
        args <- base_args2
        args$name <- val
        args$marker <- list(size = sz,
                            color = PL_colScale[val],
                            symbol = PL_symbolScale[myCat])
        plot_lytSNE <- do.call(add_trace, c(list(plot_lytSNE, data = subDat), args))
      }
      # Dummy layer to avoid random re-zooms
      xr <- range(scores2$`t-SNE Y1`)
      yr <- range(scores2$`t-SNE Y2`)
      zr <- if (tst3D) { range(scores2$`t-SNE Y3`) } else { NULL }
      plot_lytSNE <- add_trace(plot_lytSNE,
                               x = xr, y = yr, z = zr,
                               type = base_args2$type,
                               mode = "markers",
                               marker = list(size = 0.0001, color = "rgba(0,0,0,0)"),
                               name = "_bbox_",
                               showlegend = FALSE,
                               visible = TRUE)
      #
      plot_lytSNE <- layout(plot_lytSNE, title = ttl2,
                            uirevision = TRUE)
      dimRedPlotLy[["t-SNE"]] <- plot_lytSNE
      setwd(dir)
      saveWidget(plot_lytSNE, paste0(dir, "/", ttl2, ".html"))
      setwd(wd)
      #system(paste0("open \"", dir, "/", ttl2, ".html"))
    } else { warning(tsne) }
    #
    msg <- "UMAP plots, by protein group"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    cran_req <- unique(c(cran_req, "umap"))
    if (!require("umap", quietly = TRUE)) { install.packages("umap") }
    require(umap)
    UMAP <- try(umap(dimRedDat, n_components = 3L), silent = TRUE)
    dir <- paste0(wd, "/Dimensionality red. plots/UMAP")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- union(dirlist, dir)
    if (!"try-error" %in% class(UMAP)) {
      UMAPlayout <- data.frame(UMAP$layout)
      rownames(UMAPlayout) <- rownames(dimRedDat)
      kol3 <- c("Protein group", "Av. log10 abundance", "Classifier", "Size")
      m31 <- match(rownames(UMAPlayout), rownames(scores))
      UMAPlayout[, kol3] <- scores[m31, kol3]
      wBoring3 <- which(UMAPlayout$Classifier == nullVal)
      wComp3 <- which(UMAPlayout$Classifier %in% compVal)
      wReg3 <- which(UMAPlayout$Classifier %in% regVal)
      ttl3 <- "UMAP plot - Protein groups (PG-level)"
      plot <- ggplot(UMAPlayout[wBoring3,], aes(x = X1, y = X2, colour = Classifier, size = Size)) +
        geom_scattermore(shape = 16L)
      if (length(wComp3)) {
        plot <- plot + geom_point(data = UMAPlayout[wComp3,], shape = 17L)
      }
      if (length(wReg3)) {
        plot <- plot + geom_point(data = UMAPlayout[wReg3,], shape = 15L) +
          geom_text_repel(data = UMAPlayout[wReg3,], aes(label = `Protein group`))
      }
      plot <- plot +  coord_fixed() + colScale +
        geom_hline(yintercept = 0, colour = "black", alpha = 0.5) +
        geom_vline(xintercept = 0, colour = "black", alpha = 0.5) +
        scale_alpha_identity() + scale_size_identity() +
        ggtitle(ttl3, subtitle = pv) + theme_bw() +
        guides(alpha = "none", size = "none", colour = guide_legend(title = gsub("/", "/\n", ClassNm)))
      #poplot(plot, 12, 22)
      ReportCalls <- AddPlot2Report()
      suppressMessages({
        ggsave(paste0(dir, "/", ttl3, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
        ggsave(paste0(dir, "/", ttl3, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
      })
      #
      # Plotly
      tst3D <- ("X3" %in% colnames(UMAPlayout))
      base_args3 <- base_args
      base_args3$x <- ~X1
      base_args3$y <- ~X2
      if (tst3D) {
        base_args3$type <- "scatter3d"
        base_args3$z <- ~X3
      } else {
        base_args3$z <- NULL
      }
      plot_lyUMAP <- plot_ly()
      for (val in allVal) { #val <- allVal[1L]
        myCat <- val2Cat[val]
        sz <- catsEyes[myCat]*1.5
        if (val == nullVal) { subDat <- UMAPlayout } else {
          subDat <- UMAPlayout[which(UMAPlayout$Classifier == val),]
          sz <- sz*3L
        }
        args <- base_args3
        args$name <- val
        args$marker <- list(size = sz,
                            color = PL_colScale[val],
                            symbol = PL_symbolScale[myCat])
        plot_lyUMAP <- do.call(add_trace, c(list(plot_lyUMAP, data = subDat), args))
      }
      # Dummy layer to avoid random re-zooms
      xr <- range(UMAPlayout$X1)
      yr <- range(UMAPlayout$X2)
      zr <- if (tst3D) { range(UMAPlayout$X3) } else { NULL }
      plot_lyUMAP <- add_trace(plot_lyUMAP,
                               x = xr, y = yr, z = zr,
                               type = base_args3$type,
                               mode = "markers",
                               marker = list(size = 0.0001, color = "rgba(0,0,0,0)"),
                               name = "_bbox_",
                               showlegend = FALSE,
                               visible = TRUE)
      #
      plot_lyUMAP <- layout(plot_lyUMAP, title = ttl3,
                            uirevision = TRUE)
      dimRedPlotLy[["UMAP"]] <- plot_lyUMAP
      setwd(dir)
      saveWidget(plot_lyUMAP, paste0(dir, "/", ttl3, ".html"))
      setwd(wd)
      #system(paste0("open \"", dir, "/", ttl3, ".html"))
    } else { warning(umap) }
  }
} else { warning("Not enough observations to create Protein Groups-level dimensionality reduction plots!") }
ReportCalls <- AddSpace2Report()
saveFun(dimRedPlotLy, file = paste0(dir, "/DimRedPlots.RData"))
