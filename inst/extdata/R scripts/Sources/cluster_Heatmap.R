#### Heatmaps with clustering at samples and protein groups level, highlighting proteins of interest

if (clustHtMp) {
  if (clustMode == "standard") {
    normTypes <- c("None", "Norm. by row")
    clustDir <- paste0(wd, "/Clustering")
  }
  if (grepl("-tests?$", clustMode)) {
    normTypes <- "Z-scored"
    clustRegFilters <- Reg_filters[[clustMode]]$`By condition`
    clustDir <- paste0(wd, "/Reg. analysis/", clustMode, "/Heatmaps")
  }
  if (!dir.exists(clustDir)) { dir.create(clustDir, recursive = TRUE) }
  if (scrptType == "withReps") {
    dirlist <- unique(c(dirlist, dir))
    mySmpls <- clustMap$Samples <- cleanNms(clustMap$Ref.Sample.Aggregate)
    #
    I <- list(Global = mySmpls)
    # Cluster groups based on Ratios-, GO enrichment- or Normalisation groups
    # Basically any level where we have specified that we have a group smaller than the whole experiment containing sample groups which may be more
    # closely related to each other than the whole experiment.
    # It may make sense to detect automatically factors such as "Tissue", "Developmental.stage" or "Cell.type" here...
    for (myAggr in c("Ratios.Groups", "GO.enrichment.Ref.Aggr", "Norm.Groups")) {
      if (myAggr %in% colnames(Param)) {
        if (gsub(";", "", Param[[myAggr]]) %in% Aggregate.map$Aggregate.Name) {
          ClustGrp <- Param[[myAggr]]
        } else {
          ClustGrp <- "Exp"
        }
        if ((length(Exp) == 1)&&(nchar(ClustGrp) %% 3 > 0)) { ClustGrp <- Param_filter(ClustGrp, "Exp") }
        val <- Aggregate.list[[ClustGrp]]
        nms <- unlist(Aggregate.map$Characteristics[which(Aggregate.map$Aggregate.Name == gsub(";", "", ClustGrp))])
        if (length(nms) == 1) { kol <- nms } else { kol <- ClustGrp }
        ClustGrp <- list(aggregate = ClustGrp,
                         values = val,
                         names = nms,
                         column = kol)
        if (length(ClustGrp$values) > 1) {
          for (i in ClustGrp$values) {
            iNm <- paste0(ClustGrp$column, " = ", cleanNms(i))
            I[[iNm]] <- mySmpls[which(clustMap[[ClustGrp$column]] == i)]
          }
        }
      }
    }
    Fct <- c("Tis", "Cel", "Dev")
    w <- which(Fct %in% Aggregate.map$Aggregate.Name)
    if (length(w)) {
      Fct <- Fct[w]
      for (fct in Fct) {
        col <- Aggregate.map$Characteristics[[match(fct, Aggregate.map$Aggregate.Name)]]
        val <- unique(clustMap[[col]])
        if (length(val) > 1) {
          for (i in val) {
            iNm <- paste0(col, " = ", cleanNms(i))
            I[[iNm]] <- mySmpls[which(clustMap[[col]] == i)]
          }
        }
      }
    }
    I <- I[which(vapply(I, length, 1) > 1)]
  }
  if (scrptType == "noReps") {
    mySmpls <- clustMap$Samples <- cleanNms(clustMap$Experiment)
    I <- list(Global = clustMap$Samples)
    if ((MakeRatios)&&(length(unique(clustMap$`Ratios group`)) > 1)) {
      for (rtGrp in unique(clustMap$`Ratios group`)) {
        tmp <- clustMap$Samples[which(clustMap$`Ratios group` == rtGrp)]
        if (length(tmp) > 1) {
          I[[as.character(rtGrp)]] <- tmp
        }
      }
    }
  }
  Heatmaps %<o% list()
  NHClust %<o% list()
  NVClust %<o% list()
  KlustKols %<o% c()
  MaxHClust %<o% min(c(floor(nrow(PG)/2), 100)) # We want at most 20 clusters
  if (scrptType == "withReps") {
    MaxVClust %<o% length(VPAL$values)
  }
  if (scrptType == "noReps") {
    MaxVClust <- length(Exp)
  }
  VClustScl <- setNames(1:MaxVClust, paste0("Cluster", 1:MaxVClust))
  HClustScl <- setNames(1:MaxHClust, paste0("Cluster", 1:MaxHClust))
  VClustScl <- setNames(rainbow(MaxVClust), VClustScl)
  HClustScl <- setNames(rainbow(MaxHClust), HClustScl)
  VClusters %<o% list()
  # Different options for which proteins to use for vertical clustering (samples level)
  VClustUse <- "All" # Use all
  # VClustUse <- 500 # Max number of proteins to use (starting from ones with highest CV)
  # VClustUse <- 25% # Percentage of proteins to use (starting from ones with highest CV)
  # VClustUse <- "DEP" # Use only Differentially Expressed Proteins (from Reg_filters)
  VClustUse <- toupper(VClustUse)
  KlustRoot %<o% paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ") - ")
  plotLeatMaps %<o% list()
 
  
  
  #vapply(clustXprsKol, function(x) { sum(is.na(PG[clustFilt, x])) }, 1)
  #vapply(clustMap$Samples, function(x) { sum(is.na(clustDat[[x]])) }, 1)
  if (ImputeKlust) {
    if (scrptType == "withReps") {
      Gr <- clustMap[[VPAL$column]] # Here we have replicates and group at samples group level
    }
    if (scrptType == "noReps") {
      # Here we do not have replicates and group at comparison group level
      Gr <- clustMap$`Ratios group`
    }
    Gr <- setNames(match(Gr, unique(Gr)), Gr)
    clustDat <- Data_Impute2(clustDat, Gr)
    clustDatImp <- clustDat$Positions_Imputed
    clustDat <- clustDat$Imputed_data
    rownames(clustDatImp) <- rownames(clustDat)
    colnames(clustDatImp) <- clustMap$Samples
    #vapply(clustMap$Samples, function(x) { sum(is.na(clustDat[[x]])) }, 1) # df
    #vapply(clustMap$Samples, function(x) { sum(clustDatImp[, x]) }, 1) # matrix
  }
  for (i in names(I)) { #i <- names(I)[1] #i <- names(I)[2]
    smpls <- I[[i]]
    smplsMtch <- match(smpls, mySmpls)
    xMp <- clustMap[smplsMtch,]
    xprsKol <- paste0(prtRfRoot, xMp$Samples[smplsMtch])
    lXprs <- length(xprsKol)
    msg <- paste0("Creating ", c("", paste0(i, " "))[(i == "Global")+1], "heatmap",
                  c(paste0(" for ", i), "")[(i == "Global")+1], ".")
    ReportCalls <- AddMsg2Report(Space = FALSE)
    for (normType in normTypes) { #normType <- normTypes[1] #normType <- normTypes[2] #normType <- normTypes[3]
      normTypeInsrt <- paste0(" (", normTypes, ")")
      normTypeInsrt[1] <- ""
      normTypeInsrt <- normTypeInsrt[match(normType, normTypes)]
      temp <- clustDat[, smpls, drop = FALSE]
      #
      temp <- clustDat[, smpls, drop = FALSE]
      if (normType == "Z-scored") {
        # In that case we plot only differentially expressed proteins.
        # Only used for withReps for the time being.
        # We will keep it simple, making one filter only and use any protein which is significant in any test.
        myClustFilter <- PG$Label[sort(unique(unlist(lapply(names(clustRegFilters), function(nm) {
          clustRegFilters[[nm]]$Filter
        }))))]
        temp <- temp[which(rownames(temp) %in% myClustFilter),]
      }
      #tst <- apply(temp, 1, function(x) { length(is.all.good(x)) }) == lXprs
      # Filter to include only rows for which we have at least one valid value
      wAG <- which(apply(clustDat0[, smpls], 1, function(x) { length(is.all.good(x)) }) > 0)
      wAG <- which(rownames(temp) %in% rownames(clustDat0)[wAG])
      temp <- temp[wAG,]
      if (ImputeKlust) {
        whImput <- clustDatImp[wAG,]
      }
      #
      if (normType %in% c("Norm. by row", "Z-scored")) { rwMns <- rowMeans(temp) }
      if (normType == "Norm. by row") {
        temp <- sweep(temp, 1, rwMns, "-")
      }
      if (normType == "Z-scored") {
        SDs <- apply(temp, 1, function(x) { sd(proteoCraft::is.all.good(x)) })
        temp <- sweep(sweep(temp, 1, rwMns, "-"), 1, SDs, "/")
      }
      temp2 <- as.matrix(temp)
      temp3 <- temp2 <- temp2 + runif(length(temp2), min = 0, max = 1e-10) # Small error added to avoid duplicate rows where this breaks
      # Data is now normalized and either imputed or filtered, so ready for clustering
      # 1/ At samples level
      # We perform hierarchical clustering in all cases, because we want to see the dendrogram.
      # But the clusters displayed using colours may be generated using either k-means or hierarchical clustering (default).
      temp2 <- t(temp2)
      tst2 <- apply(temp2, 2, function(x) {
        #x <- proteoCraft::is.all.good(x) # it's been imputed, there are no missing values
        sd(x)/mean(x)
      })
      temp2 <- temp2[, order(tst2, decreasing = TRUE)]
      if (VClustUse != "ALL") {
        if (!is.na(suppressWarnings(as.integer(round(as.numeric(VClustUse)))))) {
          VClustUse <- as.integer(round(as.numeric(VClustUse)))
          if (VClustUse > 0) { temp2 <- temp2[, 1:VClustUse] }
        }
        if (grepl("^[0-9]+\\.?[0-9]*%$", VClustUse)) {
          VClustUse <- as.numeric(gsub("%$", "", VClustUse))/100
          if (VClustUse > 0) { temp2 <- temp2[,1:max(c(1, round(ncol(temp2)*VClustUse)))] }
        }
        if (VClustUse == "DEP") {
          PGflt <- sapply(names(Reg_filters), function(x) { #x <- names(Reg_filters)[1]
            x <- Reg_filters[[x]]$"By condition"
            sapply(names(x), function(y) { x[[y]]$Filter })
          })
          PGflt <- unique(unlist(PGflt))
          temp2 <- temp2[, which(colnames(temp2) %in% PG$Label[PGflt])]
        }
      }
      vcluster <- hclust(dist(temp2))
      vdendro <- as.dendrogram(vcluster)
      # Estimate ideal number of clusters... but ignore it!
      # In fact we know the number of sample groups, so would like to see 1 cluster per group.
      # How well groups and clusters overlap would tell us how well the clustering works, i.e. how different clusters are.
      #
      # Here we use kmeans, but findings apply to any method
      #cat("Estimating optimal number of samples-level clusters...\n")
      vnm <- paste0("Samples-level clusters number analysis - ", i)
      hnm <- paste0("Protein groups-level clusters number analysis - ", i)
      if (scrptType == "withReps") {
        Gr <- xMp[[VPAL$column]] # Here we have replicates and group at samples group level
        Gr <- setNames(match(Gr, unique(Gr)), Gr)
        NVClust[[i]] <- NGr <- length(unique(Gr))
      }
      if (scrptType == "noReps") {
        # Here we do not have replicates and group at comparison group level
        Gr <- xMp$`Ratios group`
        Gr <- setNames(match(Gr, unique(Gr)), Gr)
        NVClust[[i]] <- NGr <- max(c(1, length(mySmpls)))
      }
      if (nrow(temp2) > 2) {
        tst <- cluster::clusGap(temp2, stats::kmeans, max(c(2, min(c(nrow(temp2)-1, NGr-1)))))
        tst2 <- as.data.frame(tst$Tab)  
        yHigh <- max(tst2$gap)
        yLow <- min(tst2$gap)
        yScl <- yHigh-yLow
        tst2 <- sapply(1:NGr, function(x) { tst2$gap[x] >= tst2$gap[x+1] - tst2$SE.sim[x+1] })
        tst2 <- which(tst2)
        # I like to do one more, often these methods feel too conservative
        #if (length(tst2)) { NVClust[[i]] <- tst2[1]+1 } # Not used for now: we use NGr instead
        vplot <- factoextra::fviz_gap_stat(tst)
        tstLy <- capture.output(print(vplot$layers))
        g1 <- grep("geom_vline", tstLy)
        g2 <- grep("\\[\\[[0-9]+\\]\\]", tstLy)
        g2 <- as.numeric(gsub("\\[|\\]", "", tstLy[max(g2[which(g2 < g1)])]))
        vplot$layers[[g2]] <- NULL
        vplot <- vplot +
          geom_vline(xintercept = tst2[1]+1, colour = "red", linetype = "dashed") +
          geom_vline(xintercept = NGr, colour = "orange", linetype = "dashed") +
          geom_text(label = "opt. N. of sample clust.", x = tst2[1]+1-0.2, y = yLow+yScl*0.1,
                    colour = "red", angle = 90, hjust = 1) +
          geom_text(label = "N. of sample groups", x = NGr+0.2, y = yLow+yScl*0.1,
                    colour = "orange", angle = 90, hjust = 1) +
          theme_bw() + ggtitle(vnm)
        #poplot(vplot)
        suppressMessages({
          ggsave(paste0(clustDir, "/", vnm, normTypeInsrt, ".jpeg"), vplot, dpi = 150)
          ggsave(paste0(clustDir, "/", vnm, normTypeInsrt, ".pdf"), vplot, dpi = 150)
        })
        NVClust[[i]] <- max(c(NGr, 2))
      }
      # 2/ At protein groups level
      # As above, we always draw a dendrogram, but colours will be defined by the clustering approach.
      hcluster <- hclust(dist(temp3))
      hdendro <- as.dendrogram(hcluster)
      NHClust[[i]] <- MaxHClust
      Straps <- 10
      # Here we really want to optimize the number of clusters
      # Apply the same method for optimization for any clustering method
      # Number of cluster should not depend on method
      source(parSrc, local = FALSE)
      clusterExport(parClust, list("temp3", "Straps"), envir = environment())
      HClust_rg <- 2:MaxHClust
      tst <- setNames(parLapply(parClust, HClust_rg, function(kl) { #kl <- 2
        try(kmeans(temp3, kl, nstart = Straps)$tot.withinss, silent = TRUE)
      }), HClust_rg)
      tst <- tst[which(vapply(tst, function(x) { !"try-error" %in% class(x) }, TRUE))]
      tst <- setNames(as.numeric(tst)/kmeans(temp3, 1, nstart = 1)$tot.withinss, names(tst))
      yScl2 <- max(tst)
      tst2 <- data.frame("Number of clusters k" = as.integer(names(tst)),
                         "[tot WSS (k)]/[tot WSS (1)]" = tst,
                         check.names = FALSE)
      tst <- parLapply(parClust, HClust_rg, function(kl) {
        try(kmeans(temp3, kl, nstart = Straps)$tot.withinss, silent = TRUE)
      })
      w <- which(vapply(tst, function(x) { !"try-error" %in% class(x) }, TRUE))
      tst <- vapply(tst[w], function(x) { x }, 1)
      tst <- tst/kmeans(temp3, 1, nstart = 1)$tot.withinss
      yScl2 <- max(tst)
      tst2 <- data.frame("Number of clusters k" = HClust_rg[w],
                         "[tot WSS (k)]/[tot WSS (1)]" = tst,
                         check.names = FALSE)
      # For the elbow detection method, we need to normalize to a 1x1 graph which we can rotate by 45 degrees
      tst2$X1 <- tst2$`Number of clusters k`/MaxHClust # divide by max theoretical number of clusters 
      tst2$Y1 <- tst2$"[tot WSS (k)]/[tot WSS (1)]" # Here no need to normalize, this is a ratio
      Angle <- -pi/4
      meanX <- mean(tst2$X1)
      meanY <- mean(tst2$Y1)
      tst2$X2 <- tst2$X1 - meanX
      tst2$Y2 <- tst2$Y1 - meanY
      tst2$X2 <- tst2$X2*cos(Angle)+tst2$Y2*sin(Angle)
      tst2$Y2 <- -tst2$X2*sin(Angle)+tst2$Y2*cos(Angle)
      tst2$X2 <- tst2$X2 - mean(tst2$X2) + meanX
      tst2$Y2 <- tst2$Y2 - mean(tst2$Y2) + meanY
      w <- rev(which(tst2$Y2 == min(tst2$Y2)))[1] # Again, prefer more rather than fewer clusters
      NHClust[[i]] <- tst2$`Number of clusters k`[w]
      tst2$Size <- 1
      tst2$Size[w] <- 2
      xMin <- min(c(tst2$X1, tst2$X2))
      xMax <- max(c(tst2$X1, tst2$X2))
      xScl <- xMax-xMin
      yMin <- min(c(tst2$Y1, tst2$Y2))
      yMax <- max(c(tst2$Y1, tst2$Y2))
      yScl <- yMax-yMin
      hplot <- ggplot(tst2) +
        geom_segment(x = tst2$X2[w], y = tst2$Y2[w],
                     xend = tst2$X1[w], yend = tst2$Y1[w], color = "grey", linetype = "dotted") +
        geom_point(aes(x = X1, y = Y1, size = Size), color = "blue") +
        geom_point(aes(x = X2, y = Y2, size = Size), color = "red") +
        geom_text(x = xScl*0.98+xMin, y = yMin+yScl*0.98, color = "blue", label = "ratio of tot. WSS", hjust = 1) +
        geom_text(x = xScl*0.98+xMin, y = yMin+yScl*0.962, color = "red", label = "ratio of tot. WSS, -pi/4 rotation", hjust = 1) +
        geom_hline(yintercept = tst2$Y2[w], color = "red", linetype = "dashed") +
        geom_vline(xintercept = tst2$X1[w], color = "deepskyblue", linetype = "dashed") +
        geom_text(x = xMin+0.01*xScl, y = tst2$Y2[w]+0.02*yScl, color = "red", label = "Elbow", hjust = 0) +
        geom_text(x = tst2$X1[w]-0.02*xScl, y = yScl*0.9, angle = 90, color = "deepskyblue",
                  label = paste0("Optimal = ", tst2$`Number of clusters k`[w], " clusters"), hjust = 1) +
        ggtitle(hnm) + scale_size_identity() + theme_bw() +
        theme(legend.position = "none") + ylab("Normalised total Within-clusters vs Total Sum of Squares")
      #poplot(hplot)
      suppressMessages({
        ggsave(paste0(clustDir, "/", hnm, normTypeInsrt, ".jpeg"), hplot, dpi = 150)
        ggsave(paste0(clustDir, "/", hnm, normTypeInsrt, ".pdf"), hplot, dpi = 150)
      })
      # Apply cutoffs
      if (KlustMeth == 1) {
        VClusters[[i]] <- kmeans(temp3, NVClust[[i]], 100)$cluster
        tempClust <- kmeans(temp3, NHClust[[i]], 100)$cluster
      }
      if (KlustMeth == 2) {
        VClusters[[i]] <- cutree(vcluster, NVClust[[i]])
        tempClust <- cutree(hcluster, NHClust[[i]])
      }
      KlKol <- paste0(KlustRoot, i)
      KlustKols <- unique(c(KlustKols, KlKol))
      PG[[KlKol]] <- tempClust[match(PG$Label, names(tempClust))]
      #
      Width <- nrow(temp)
      Height <- length(mySmpls)
      #
      if (ImputeKlust) {
        whImps <- which(whImput[match(rownames(temp), rownames(whImput)),
                                colnames(temp)], arr.ind = TRUE)
        temp[whImps] <- NA
      }
      #
      # Get dendrograms
      vddata <- ggdendro::dendro_data(vdendro)
      hddata <- ggdendro::dendro_data(hdendro)
      # vdendro.plot <- ggdendrogram(data = vdendro) +
      #   theme(axis.text.y = element_text(size = 0.1), plot.margin = margin(0, 0, 0, 0, "cm"))
      # poplot(vdendro.plot)
      # Modify dendrograms
      # - Vertical
      vlabs <- ggdendro::label(vddata)
      vlabs$Cluster <- as.factor(VClusters[[i]][match(vlabs$label, smpls)])
      vSeg <- vddata$segments
      # - Horizontal
      hlabs <- ggdendro::label(hddata)
      hlabs$Cluster <- as.factor(tempClust[match(hlabs$label, names(tempClust))])
      hSeg <- hddata$segments
      # Rotate samples-level dendrogram: x -> y, y -> x (for ease of calculation, gets inverted later)
      vlabs$y <- vlabs$x
      vlabs$x <- -0.3
      x <- vSeg$x
      y <- vSeg$y
      vSeg$y <- x
      vSeg$x <- y
      xend <- vSeg$xend
      yend <- vSeg$yend
      vSeg$yend <- xend
      vSeg$xend <- yend
      # Adjust width/height
      # - Vertical dendrogram
      vnch <- Width*0.07*max(nchar(vlabs$label))/12
      xmn <- min(c(vSeg$x, vSeg$xend))
      xmx <- max(c(vSeg$x, vSeg$xend))
      vSeg$x <- -(vSeg$x - xmn)*Width*0.07/xmx - 1.5 - vnch
      vSeg$xend <- -(vSeg$xend - xmn)*Width*0.07/xmx - 1.5 - vnch
      # - Horizontal dendrogram
      hnch <- Height*0.07*max(nchar(hlabs$label))/800
      ymn <- min(c(hSeg$y, hSeg$yend))
      ymx <- max(c(hSeg$y, hSeg$yend))
      hSeg$y <- Height + hnch + 1.5 + (hSeg$y - ymn)*Height*0.07/ymx
      hSeg$yend <- Height + hnch + 1.5 + (hSeg$yend - ymn)*Height*0.07/ymx
      hSeg$x <- hSeg$x - 0.5
      hSeg$xend <- hSeg$xend - 0.5
      # Order labels by order of appearance
      hlabs <- hlabs[order(hlabs$x, decreasing = FALSE),]
      vlabs <- vlabs[order(vlabs$y, decreasing = FALSE),]
      # Re-order our matrix based on extracted dendrogram labels
      temp <- temp[, match(vlabs$label, colnames(temp))]
      temp <- temp[match(hlabs$label, rownames(temp)),]
      if (ImputeKlust) {
        # Just in case: actually we do not use these downstream currently
        whImput <- whImput[, match(vlabs$label, colnames(whImput))]
        whImput <- whImput[match(hlabs$label, rownames(whImput)),]
      }
      # Re-introduce missing values
      MaxChar <- 13
      hlabs$label2 <- hlabs$label
      w <- which(nchar(hlabs$label2) > MaxChar)
      hlabs$label2[w] <- paste0(substr(hlabs$label2[w], 1, MaxChar-3), "...")
      # Create heatmap
      temp$Rowname <- row.names(temp)
      temp2 <- set_colnames(reshape::melt(temp, id.vars = "Rowname"), c("Label", "Sample", "value"))
      temp2$Label <- as.character(temp2$Label)
      temp2$Sample <- as.character(temp2$Sample)
      temp2$"Leading protein IDs" <- PG$"Leading protein IDs"[match(temp2$Label, PG$Label)]
      temp2$Xmax <- match(temp2$Label, hlabs$label) # Explicitly should be the case now!
      temp2$label2 <- hlabs$label2[temp2$Xmax]
      temp2$Xmin <- temp2$Xmax-1
      temp2$Ymax <- vlabs$y[match(temp2$Sample, vlabs$label)]
      temp2$Ymin <- temp2$Ymax-1
      w1 <- which(temp2$Colour == "green")
      w2 <- which((temp2$Ymin == max(temp2$Ymin))&(temp2$Colour == "green"))
      # Color and fill scales
      wV <- round(c(1:NVClust[[i]])*MaxVClust/NVClust[[i]])
      wH <- round(c(1:NHClust[[i]])*MaxHClust/NHClust[[i]])
      vClScl <- setNames(VClustScl[wV], 1:length(wV))
      hClScl <- setNames(HClustScl[wH], 1:length(wH))
      VcolScale <- scale_color_manual(name = "Samples cluster", values = vClScl)
      VfillScale <- scale_fill_manual(name = "Samples cluster", values = vClScl)
      HcolScale <- scale_color_manual(name = "Protein groups cluster", values = hClScl)
      HfillScale <- scale_fill_manual(name = "Protein groups cluster", values = hClScl)
      # Create heatmap plot
      Xlim <- c(NA, Width)
      Ylim <- c(-10, max(c(max(hSeg$y) + Height*0.6), 20))
      #
      prot.list.marks <- FALSE
      if (prot.list.Cond) {
        w0 <- which(temp2$Ymin == 0)
        g <- grsep2(prot.list, temp2$"Leading protein IDs"[w0])
        #tst <- unlist(strsplit(temp2$"Leading protein IDs"[w0], ";"))[1]
        #g <- grsep2(tst, temp2$"Leading protein IDs"[w0])
        if (length(g)) {
          prot.list.marks <- TRUE
          Ylim[1] <- -20
          temp2c <- temp2[w0[g], , drop = FALSE]
        }
      }
      # Main data
      temp2a <- temp2[, c("Xmin", "Ymin", "value", "Label", "Sample")]
      # Colour scale
      temp2b <- data.frame(Xmin = 0:round(Width*0.1),
                           Ymin = Ylim[2]*0.8)
      Mn <- min(temp2a$value, na.rm = TRUE)
      Mx <- max(temp2a$value, na.rm = TRUE)
      temp2b$value <- Mn + temp2b$Xmin*(Mx-Mn)/max(temp2b$Xmin)
      temp2b$Xmin <- temp2b$Xmin-Width*0.15
      temp2b$Label <- temp2b$Sample <- NA
      w2a <- 1:nrow(temp2a)
      w2b <- 1:nrow(temp2b) + max(w2a)
      temp2a <- rbind(temp2a, temp2b)
      #
      # Annotations
      xCntr <- Width*0.6
      yPadUp <- round(rev((1:49 + 0.3)*4/49.3), 2)
      ind <- Height
      if (ind > 50) { ind <- 50 }
      yPadUp <- yPadUp[ind-1]
      nm <- paste0("Clust. heatmap - ", i)
      GrLab <- data.frame(label = c(nm,
                                    normType,
                                    "Protein group",
                                    "Sample",
                                    "Min",
                                    "Max"),
                          x = c(xCntr,
                                xCntr,
                                xCntr,
                                min(vSeg$x)-Width*0.11,
                                min(temp2b$Xmin),
                                max(temp2b$Xmin)),
                          y = c(max(c(hSeg$y, hSeg$yend)) + yPadUp*0.75,
                                max(c(hSeg$y, hSeg$yend)) + yPadUp*0.5,
                                max(c(hSeg$y, hSeg$yend)) + yPadUp*0.25,
                                Height*0.5,
                                temp2b$Ymin[1]-1,
                                temp2b$Ymin[1]-1),
                          angle = c(0, 0, 0, 90, 0, 0),
                          size = c(5, 3.5, 4, 4, 3, 3),
                          fontface = c("bold", "italic", "plain", "plain", "plain", "plain"))
      # Samples-level: how well do clusters fit expectations
      smplsClust <- vlabs[, c("label", "Cluster")]
      if (scrptType == "withReps") {
        smplsClust$Group <- as.factor(xMp[match(smplsClust$label, xMp$Samples), VPAL$column])
      }
      if (scrptType == "noReps") {
        smplsClust$Group <- as.factor(xMp$`Ratios group`[match(smplsClust$label, xMp$Samples)])
      }
      smplsClust$Cluster <- as.factor(paste0("Cluster", as.character(smplsClust$Cluster)))
      if (length(levels(smplsClust$Cluster)) > 1) {
        tstSmplClust <- table(smplsClust$Group, smplsClust$Cluster)
        ClustChiSqTst <- suppressWarnings(chisq.test(tstSmplClust))
        GrLab <- rbind(GrLab,
                       data.frame(label = paste0("Chi-squared contingency test P-value: ", round(ClustChiSqTst$p.value, 5)),
                                  x = min(vSeg$x)-Width*0.095,
                                  y = Height*0.5,
                                  angle = 90,
                                  size = 3,
                                  fontface = "italic"))
      }
      #
      Xlim[1] <- min(c(vSeg$x, vSeg$xend))-Width*0.15
      Ylim[1] <- min(c(Ylim[1], min(temp2a$Ymin)))
      Xlim[2] <- max(c(Xlim[2], max(temp2a$Xmin)+1))
      Ylim[2] <- max(c(Ylim[2], max(temp2a$Ymin)+1))
      xCntr <- mean(Xlim)
      yCntr <- mean(Ylim)
      yScale <- Ylim[2]-Ylim[1]
      #
      # Create heatmap
      heatmap.plot <- ggplot(temp2a[w2a,]) +
        geom_rect(aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value, text = Label)) +
        geom_rect(data = temp2a[w2b,],
                  aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
              panel.background = element_rect(fill = "transparent", color = NA),
              plot.margin = margin(0, 0, 0, 0, "cm")) +
        #scale_fill_gradient2(low = "darkblue", mid = "lightgrey", high = "darkred") +
        scale_fill_viridis(option = "D") +
        xlab(NULL) + ylab(NULL) + theme(legend.position = "none") +
        xlim(Xlim[1], Xlim[2]) + ylim(Ylim[1], Ylim[2])
      #poplot(heatmap.plot, 12, 20)
      # heatmap.plot <- heatmap.plot + 
      #   geom_text(data = temp2a[which(temp2a$Xmin == min(temp2$Xmin)),],
      #             aes(x = Xmin, y = Ymin, label = Sample))
      # Title, axis labels, colour scale annotations
      #
      # Add labels and dendrograms:
      # - Title, scale, chi-squared test
      heatmap.plot <- heatmap.plot +
        new_scale("size") + scale_size_identity() +
        geom_text(data = GrLab, aes(label = label, x = x, y = y, angle = angle, size = size,
                                    fontface = fontface),
                  color = "black", hjust = 0.5)
      #poplot(heatmap.plot, 12, 20)
      # - Cluster boxes
      heatmap.plot <- heatmap.plot + new_scale_color() + HcolScale
      if (KlustMeth == 2) {
        Clutst <- aggregate(hlabs$x, list(hlabs$Cluster), function(x) { min(x)-1 })
        colnames(Clutst)[1] <- "Cluster"
        Clutst$xend <- aggregate(hlabs$x, list(hlabs$Cluster), max)$x
        Clutst$mid <- (Clutst$xend+Clutst$x)/2
        heatmap.plot <- heatmap.plot +
          geom_rect(data = Clutst, aes(xmin = x, xmax = xend, colour = Cluster),
                    ymin = 0, ymax = Height, fill = NA) +
          geom_text(data = Clutst, aes(x = mid, label = Cluster, colour = Cluster),
                    y = Height/2, hjust = 0.5, vjust = 0.5, cex = 5)
        #poplot(heatmap.plot, 12, 20)
      }
      #poplot(heatmap.plot, 12, 20)
      # - Horizontal dendrogram and protein  names
      heatmap.plot <- heatmap.plot +
        geom_segment(data = hSeg, linewidth = 0.25,
                     aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_text(data = hlabs, aes(x = x-0.5, label = label2, colour = Cluster),
                  y = Height+0.05, angle = 90, cex = 0.5, hjust = 0, vjust = 0)
      # - Vertical dendrogram and sample names
      heatmap.plot <- heatmap.plot +
        geom_segment(data = vSeg, linewidth = 0.25,
                     aes(x = x, y = y - 0.5, xend = xend, yend = yend - 0.5)) +
        new_scale_color() + VcolScale +
        geom_text(data = vlabs, aes(y = y - 0.5, label = label, colour = Cluster),
                  x = -0.5, hjust = 1, vjust = 0.5, cex = 2.5)
      # - Proteins of interest
      if (prot.list.marks) {
        heatmap.plot <- heatmap.plot +
          geom_point(data = temp2c, aes(x = Xmin+0.5), y = -0.5, colour = "red", fill = "red", shape = 17) +
          geom_text(data = temp2c, aes(x = Xmin+0.5, label = label2),
                    y = -1, colour = "red", angle = -60, hjust = 0, cex = 2)
      }
      #poplot(heatmap.plot, 12, 20)
      suppressMessages({
        ggsave(paste0(clustDir, "/", nm, normTypeInsrt, ".jpeg"), heatmap.plot)
        ggsave(paste0(clustDir, "/", nm, normTypeInsrt, ".pdf"), heatmap.plot)
      })
      #
      # Plotly version
      tempLy <- temp2a[w2a,]
      tempLy$Sample <- factor(tempLy$Sample, levels = smpls)
      ##tempLy$Ymin <- tempLy$Ymin+0.5
      plotleatmap <- plot_ly(data = tempLy, x = ~Xmin, y = ~Ymin, z = ~value, type = "heatmap", hovertext = tempLy$Label)
      # I cannot find a way to remove tick marks!!!!!
      plLyV <- vSeg
      ##plLyV$x <- -Width*0.2-plLyV$x 
      ##plLyV$xend <- -Width*0.2-plLyV$xend
      plLyV$y <- plLyV$y-1
      plLyV$yend <- plLyV$yend-1 
      plotleatmap <- add_segments(plotleatmap, x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                                  data = plLyV, inherit = FALSE, color = I("black"),
                                  showlegend = FALSE)
      plLyH <- hSeg
      ##plLyH$x <- plLyH$x-0.5
      ##plLyH$xend <- plLyH$xend-0.5
      plLyH$y <- plLyH$y - hnch - 1.5
      plLyH$yend <- plLyH$yend - hnch - 1.5
      plotleatmap <- add_segments(plotleatmap, x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                                  data = plLyH, inherit = FALSE, color = I("black"),
                                  showlegend = FALSE)
      if (KlustMeth == 2) { # Cluster shapes do not appear to work
        plotleatmap <- add_segments(plotleatmap, x = ~ x, xend = ~ xend,
                                    y = I(-0.2*as.numeric(Clutst$Cluster))-1,
                                    yend = I(-0.2*as.numeric(Clutst$Cluster))-1, data = Clutst,
                                    inherit = FALSE, color = ~ Cluster, showlegend = FALSE)
      }
      vlabs2 <- vlabs
      vlabs2$x <- -vnch/2
      vlabs2$y <- vlabs2$y-1
      plotleatmap <- add_trace(plotleatmap, data = vlabs2, y = ~y, x = ~x, text = ~label,
                               color = I("black"), inherit = FALSE, type = "scatter",
                               mode = "text", showlegend = FALSE)
      plotleatmap <- layout(plotleatmap, title = nm,
                            xaxis = list(title = "Protein groups",
                                         tickmode = "array", tickvals = NULL, showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
                            yaxis = list(title = "Samples",
                                         tickmode = "array", tickvals = NULL, showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE))
      if (prot.list.marks) {
        temp2c$X <- (temp2c$Xmin+temp2c$Xmax)/2
        plotleatmap <- add_trace(plotleatmap, data = temp2c, y = I(-0.501), x = ~X,
                                 text = ~Label, color = ~Label, inherit = FALSE,
                                 type = "scatter", mode = "markers", showlegend = FALSE)
      }
      setwd(clustDir)
      saveWidget(plotleatmap, paste0(clustDir, "/", nm, normTypeInsrt, ".html"))
      setwd(wd)
      #system(paste0("open \"", clustDir, "/", nm, normTypeInsrt, ".html\""))
      plotLeatMaps[[i]][[normType]] <- plotleatmap
    }
  }
  saveFun(plotLeatMaps, file = paste0(clustDir, "/HeatMaps.RData"))
  #
  temp <- PG[, c("Leading protein IDs", "Protein names", "Genes", "Mol. weight [kDa]",
                 clustXprsKol, KlustKols)]
  #
  flPath <- paste0(clustDir, "/Protein Groups and Clusters.csv")
  tst <- try(write.csv(temp, file = flPath, row.names = FALSE), silent = TRUE)
  while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
    dlg_message(paste0("File \"", flPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
    tst <- try(write.csv(temp, file = flPath, row.names = FALSE), silent = TRUE)
  }
}
