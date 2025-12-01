#### Heatmaps with clustering at samples and protein groups level, highlighting proteins of interest
Src <- paste0(libPath, "/extdata/R scripts/Sources/cluster_Heatmap_Prep.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
drawPlotly <- TRUE
if (!exists("plotLeatMaps")) { plotLeatMaps <- list() }
plotLeatMaps %<o% plotLeatMaps
heatMaps <- list() # Unlike plotLeatMaps, not persistent
#drawPlotly <- FALSE
if (clustHtMp) {
  if (clustMode == "standard") {
    normTypes <- c("Norm. by row", "None")
    # The order matters: the first is used to generate the PG-level cluster displayed on both heatmaps!
    clustDir <- paste0(wd, "/Clustering")
  }
  if ((grepl("-tests?$", clustMode))||(clustMode %in% c("re-localisation", "SAINTexpress"))) {
    normTypes <- "Z-scored"
    clustRegFilters <- Reg_filters[[clustMode]]$`By condition`
    clustDir <- paste0(wd, "/Reg. analysis/", clustMode, "/Heatmaps")
  }
  verbose <- FALSE
  if (!dir.exists(clustDir)) { dir.create(clustDir, recursive = TRUE) }
  if (scrptType == "withReps") {
    dirlist <- unique(c(dirlist, dir))
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
  VClustScl <- setNames(1:MaxVClust, paste0("Cluster", 1:MaxVClust))
  HClustScl <- setNames(1:MaxHClust, paste0("Cluster", 1:MaxHClust))
  VClustScl <- setNames(rainbow(MaxVClust), VClustScl)
  HClustScl <- setNames(rainbow(MaxHClust), HClustScl)
  #
  # Different options for which proteins to use for vertical clustering (samples level)
  VClustUse <- "All" # Use all
  # VClustUse <- 500 # Max number of proteins to use (starting from ones with highest CV)
  # VClustUse <- 25% # Percentage of proteins to use (starting from ones with highest CV)
  # VClustUse <- "DEP" # Use only Differentially Expressed Proteins (from Reg_filters)
  VClustUse <- toupper(VClustUse)
  if (("Compartment marker" %in% colnames(PG))&&(sum(PG$"Compartment marker" != ""))) {
    uMark <- unique(SubCellMark)
    markColors <- setNames(turbo(length(uMark)), uMark)
  }
  #vapply(clustXprsKol, function(x) { sum(is.na(PG[clustFilt, x])) }, 1)
  #vapply(clustMap$Samples, function(x) { sum(is.na(clustDat[[x]])) }, 1)
  h_clustLst <- v_clustLst <- list()
  for (i in names(I)) { #i <- names(I)[1] #i <- names(I)[2]
    nm <- paste0("Clust. heatmap - ", i)
    smpls <- I[[i]]
    smplsMtch <- match(smpls, mySmpls)
    xMp <- clustMap[smplsMtch,]
    xprsKol <- paste0(prtRfRoot, xMp$Samples[smplsMtch])
    lXprs <- length(xprsKol)
    msg <- paste0("Creating ", c("", paste0(i, " "))[(i == "Global")+1], "heatmap",
                  c(paste0(" for ", i), "")[(i == "Global")+1], ".")
    ReportCalls <- AddMsg2Report(Space = FALSE)
    for (normType in normTypes) { #normType <- normTypes[1] #normType <- normTypes[2] #normType <- normTypes[3]
      clustNm <- paste0(i, " - ", normType)
      if (length(normTypes) > 1) {
        cat(" -", normType, "\n")
      }
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
      if (length(wAG) > 3) {
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
        v_clust <- hclust(dist(temp2))
        v_clustLst[[clustNm]] <- v_clust
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
          NVClust[[clustNm]] <- NGr <- length(unique(Gr))
        }
        if (scrptType == "noReps") {
          # Here we do not have replicates and group at comparison group level
          Gr <- xMp$`Ratios group`
          Gr <- setNames(match(Gr, unique(Gr)), Gr)
          NVClust[[clustNm]] <- NGr <- max(c(1, length(mySmpls)))
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
          #if (length(tst2)) { NVClust[[clustNm]] <- tst2[1]+1 } # Not used for now: we use NGr instead
          if (verbose) {
            vplot <- factoextra::fviz_gap_stat(tst)
            tstLy <- capture.output(print(vplot$layers))
            g1 <- grep("geom_vline", tstLy)
            g2 <- grep("\\[\\[[0-9]+\\]\\]", tstLy)
            if (length(g2)) {
              g2 <- as.numeric(gsub("\\[|\\]", "", tstLy[max(g2[which(g2 < g1)])]))
              vplot$layers[[g2]] <- NULL
            }
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
          }
          NVClust[[clustNm]] <- max(min(c(floor(ncol(temp2)/2), c(NGr, 2))))
        }
        # 2/ At protein groups level
        # As above, we always draw a dendrogram, but colours will be defined by the clustering approach.
        h_clust <- hclust(dist(temp3)) # This is what we will draw, independent of normalisation type
        h_clustLst[[clustNm]] <- h_clust
        if (normType != "None") { #... but we will use cluster colors from "Norm. by row" for "None"
          #NHClust[[clustNm]] <- MaxHClust
          #NHClust[[i]] <- MaxHClust
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
          #NHClust[[clustNm]] <- max(c(2, tst2$`Number of clusters k`[w] + 1), na.rm = TRUE)
          NHClust[[i]] <- max(c(2, tst2$`Number of clusters k`[w] + 1), na.rm = TRUE)
          tst2$Size <- 1
          tst2$Size[w] <- 2
          xMin <- min(c(tst2$X1, tst2$X2))
          xMax <- max(c(tst2$X1, tst2$X2))
          xScl <- xMax-xMin
          yMin <- min(c(tst2$Y1, tst2$Y2))
          yMax <- max(c(tst2$Y1, tst2$Y2))
          yScl <- yMax-yMin
          if (verbose) {
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
          }
        }
        # Apply cutoffs
        if (KlustMeth == 1) {
          #HClusters[[clustNm]] <- kmeans(temp3, NHClust[[clustNm]], 100)$cluster
          if (normType != "None") { # If normType == "None", we inherit this from "Norm. by row"
            HClusters[[i]] <- kmeans(temp3, NHClust[[i]], 100)$cluster
          }
          VClusters[[clustNm]] <- kmeans(t(temp3), NVClust[[clustNm]], 100)$cluster
        }
        if (KlustMeth == 2) {
          #HClusters[[clustNm]] <- cutree(h_clust, NHClust[[clustNm]])
          if (normType != "None") { # If normType == "None", we inherit this from "Norm. by row"
            HClusters[[i]] <- cutree(h_clust, NHClust[[i]])
          }
          VClusters[[clustNm]] <- cutree(v_clust, NVClust[[clustNm]])
        }
        KlKol <- paste0(KlustRoot, i)
        if ((clustMode == "standard")&&(normType == "Norm. by row")) {
          # We use clusters only from normalised by row,
          # because we want to see the effect of relative expression
          KlustKols <- unique(c(KlustKols, KlKol))
          #PG[[KlKol]] <- HClusters[[clustNm]][match(PG$Label, names(HClusters[[clustNm]]))]
          PG[[KlKol]] <- HClusters[[i]][match(PG$Label, names(HClusters[[i]]))]
        }
        #
        Width <- nrow(temp)
        Height <- length(smpls)
        #
        if (ImputeKlust) {
          whImps <- which(whImput[match(rownames(temp), rownames(whImput)),
                                  colnames(temp)], arr.ind = TRUE)
          temp[whImps] <- NA
        }
        #
        # Get dendrograms
        h_ddata <- hclust_to_seg(h_clust)
        v_ddata <- hclust_to_seg(v_clust)
        # Modify dendrograms
        # - Vertical
        v_labs <- v_ddata$labels
        v_labs$Cluster <- as.factor(VClusters[[clustNm]][match(v_labs$label, names(VClusters[[clustNm]]))])
        v_Seg <- v_ddata$segments
        # - Horizontal
        h_labs <- h_ddata$labels
        #h_labs$Cluster <- as.factor(HClusters[[clustNm]][match(h_labs$label, names(HClusters[[clustNm]]))])
        h_labs$Cluster <- as.factor(HClusters[[i]][match(h_labs$label, names(HClusters[[i]]))])
        h_Seg <- h_ddata$segments
        # Rotate samples-level dendrogram: x -> y, y -> x
        v_labs$y <- v_labs$x
        v_labs$x <- -0.3
        x <- v_Seg$x
        y <- v_Seg$y
        v_Seg$y <- x
        v_Seg$x <- -y
        xend <- v_Seg$xend
        yend <- v_Seg$yend
        v_Seg$yend <- xend
        v_Seg$xend <- -yend
        # Re-scale
        # - Vertical dendrogram
        v_labs_wdth <- Width*0.2 # if heatmap is 100 units then reserve 20 units for labels
        v_dendr_wdth <- Width*0.05 # if heatmap is 100 units then reserve 5 units for labels
        xMn <- min(c(v_Seg$x, v_Seg$xend))
        xMx <- max(c(v_Seg$x, v_Seg$xend))  #Should be 0
        xWdth <- xMx - xMn
        v_Seg$x <- v_Seg$x*v_dendr_wdth/xWdth - v_labs_wdth*1.05
        v_Seg$xend <- v_Seg$xend*v_dendr_wdth/xWdth - v_labs_wdth*1.05
        # - Horizontal dendrogram
        hnch <- Height*0.07*max(nchar(h_labs$label))/800
        ymn <- min(c(h_Seg$y, h_Seg$yend))
        ymx <- max(c(h_Seg$y, h_Seg$yend))
        h_Seg$y <- Height + hnch + 1.5 + (h_Seg$y - ymn)*Height*0.07/ymx
        h_Seg$yend <- Height + hnch + 1.5 + (h_Seg$yend - ymn)*Height*0.07/ymx
        h_Seg$x <- h_Seg$x - 0.5
        h_Seg$xend <- h_Seg$xend - 0.5
        # Order labels by order of appearance
        h_labs <- h_labs[order(h_labs$x, decreasing = FALSE),]
        v_labs <- v_labs[order(v_labs$y, decreasing = FALSE),]
        # Re-order our matrix based on extracted dendrogram labels
        temp <- temp[, match(v_labs$label, colnames(temp))]
        temp <- temp[match(h_labs$label, rownames(temp)),]
        if (ImputeKlust) {
          # Just in case: actually we do not use these downstream currently
          whImput <- whImput[, match(v_labs$label, colnames(whImput))]
          whImput <- whImput[match(h_labs$label, rownames(whImput)),]
        }
        # Re-introduce missing values
        MaxChar <- 13
        h_labs$label2 <- gsub("^[^ ]+ - ", "", h_labs$label)
        w <- which(nchar(h_labs$label2) > MaxChar)
        h_labs$label2[w] <- paste0(substr(h_labs$label2[w], 1, MaxChar-3), "...")
        # Create heatmap
        temp$Rowname <- row.names(temp)
        temp2 <- set_colnames(reshape::melt(temp, id.vars = "Rowname"), c("Label", "Sample", "value"))
        temp2$Label <- as.character(temp2$Label)
        temp2$Sample <- as.character(temp2$Sample)
        temp2$"Leading protein IDs" <- PG$"Leading protein IDs"[match(temp2$Label, PG$Label)]
        temp2$Xmax <- match(temp2$Label, h_labs$label) # Explicitly should be the case now!
        temp2$label2 <- h_labs$label2[temp2$Xmax]
        temp2$Xmin <- temp2$Xmax-1
        temp2$Ymax <- v_labs$y[match(temp2$Sample, v_labs$label)] # Also explicit!
        temp2$Ymin <- temp2$Ymax-1
        w1 <- which(temp2$Colour == "green")
        w2 <- which((temp2$Ymin == max(temp2$Ymin))&(temp2$Colour == "green"))
        # Color and fill scales
        wV <- round(c(1:NVClust[[clustNm]])*MaxVClust/NVClust[[clustNm]])
        #wH <- round(c(1:NHClust[[clustNm]])*MaxHClust/NHClust[[clustNm]])
        wH <- round(c(1:NHClust[[i]])*MaxHClust/NHClust[[i]])
        vClScl <- setNames(VClustScl[wV], seq_along(wV))
        hClScl <- setNames(HClustScl[wH], seq_along(wH))
        VcolScale <- scale_color_manual(name = "Samples cluster", values = vClScl)
        VfillScale <- scale_fill_manual(name = "Samples cluster", values = vClScl)
        HcolScale <- scale_color_manual(name = "Protein groups cluster", values = hClScl)
        HfillScale <- scale_fill_manual(name = "Protein groups cluster", values = hClScl)
        # Create heatmap plot
        Ylim <- c(-10, max(c(max(h_Seg$y) + Height*0.6), 20))
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
        sclRng <- 0:20
        sclWdth <- Width*0.1/20
        temp2b <- data.frame(i = sclRng,
                             Xmin = -rev(sclRng)*sclWdth - 1 - Width*0.05,
                             Ymin = Ylim[2]*0.8)
        Mn <- min(temp2a$value, na.rm = TRUE)
        Mx <- max(temp2a$value, na.rm = TRUE)
        temp2b$value <- Mn + temp2b$i*(Mx-Mn)/max(temp2b$i)
        temp2b$Label <- temp2b$Sample <- NA
        w2a <- 1:nrow(temp2a)
        w2b <- 1:nrow(temp2b) + max(w2a)
        temp2b$i <- NULL
        temp2a <- rbind(temp2a, temp2b)
        #
        # Annotations
        xCntr <- Width*0.6
        yPadUp <- round(rev((1:49 + 0.3)*4/49.3), 2)
        #ind <- Height
        ind <- Width # Counter-intuitive... but this padding is to fit labels, whose cex and hence length is Width-dependent
        if (ind > 50) { ind <- 50 }
        yPadUp <- yPadUp[ind-1]
        h_SegMx <- max(c(h_Seg$y, h_Seg$yend))
        lftSpace <- c(v_labs_wdth + v_dendr_wdth) * 1.1
        GrLab <- data.frame(label = c(nm,
                                      normType,
                                      "Protein groups",
                                      "Sample",
                                      "Min",
                                      "Max"),
                            x = c(rep(xCntr, 3),
                                  -lftSpace*1.25,
                                  min(temp2b$Xmin),
                                  max(temp2b$Xmin)+sclWdth),
                            y = c(h_SegMx + yPadUp + c(3.5, 3, 2.5)*Height/10,
                                  Height*0.5,
                                  rep(temp2b$Ymin[1]-Height*0.05, 2)),
                            angle = c(0, 0, 0, 90, 0, 0),
                            size = c(3, 2.5, 2, 2.5, 2, 2),
                            fontface = c("bold", "italic", "plain", "plain", "plain", "plain"))
        # Samples-level: how well do clusters fit expectations
        smplsClust <- v_labs[, c("label", "Cluster")]
        if (scrptType == "withReps") {
          smplsClust$Group <- as.factor(xMp[match(smplsClust$label, xMp$Samples), VPAL$column])
        }
        if (scrptType == "noReps") {
          if (MakeRatios) {
            smplsClust$Group <- as.factor(xMp$`Ratios group`[match(smplsClust$label, xMp$Samples)])
          } else {
            smplsClust$Group <- factor(smplsClust$label, levels = smplsClust$label)
          }
        }
        smplsClust$Cluster <- as.factor(paste0("Cluster", as.character(smplsClust$Cluster)))
        if (length(levels(smplsClust$Cluster)) > 1) {
          tstSmplClust <- table(smplsClust$Group, smplsClust$Cluster)
          ClustChiSqTst <- suppressWarnings(chisq.test(tstSmplClust))
          GrLab <- rbind(GrLab,
                         data.frame(label = paste0("Chi-squared contingency test P-value: ", round(ClustChiSqTst$p.value, 5)),
                                    x = - lftSpace*1.15,
                                    y = Height*0.5,
                                    angle = 90,
                                    size = 3,
                                    fontface = "italic"))
        }
        #
        Xlim <- c(min(GrLab$x)-Width*0.01,
                  Width*1.01)
        Ylim[1] <- min(c(Ylim[1], min(temp2a$Ymin)))
        Ylim[2] <- max(c(Ylim[2], max(c(GrLab$y, temp2b$Ymin+1))))
        xCntr <- mean(Xlim)
        yCntr <- mean(Ylim)
        xScale <- Xlim[2]-Xlim[1]
        yScale <- Ylim[2]-Ylim[1]
        #
        # Create heatmap
        heatmap.plot <- ggplot(temp2a[w2a,]) +
          geom_rect(aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value, text = Label)) +
          geom_rect(data = temp2a[w2b,],
                    aes(xmin = Xmin, xmax = Xmin+sclWdth, ymin = Ymin, ymax = Ymin+Height*0.05, fill = value)) +
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
        # - Clusters
        clustDF <- h_labs[, c("x", "Cluster", "label")]
        #w <- which(!temp2a$Label[w2a] %in% clustDF$label)
        heatmap.plot <- heatmap.plot +
          new_scale_fill() + HfillScale +
          geom_rect(data = clustDF, ymin = -0.75, ymax = -0.25,
                    aes(xmin = x-1, xmax = x, fill = Cluster, text = Cluster))
        #
        #poplot(heatmap.plot, 12, 20)
        # - Horizontal dendrogram and protein names
        heatmap.plot <- heatmap.plot +
          geom_segment(data = h_Seg, linewidth = 0.25,
                       aes(x = x, y = y+yPadUp, xend = xend, yend = yend+yPadUp)) +
          geom_text(data = h_labs, aes(x = x-0.5, label = label2, colour = Cluster),
                    y = Height+0.05, angle = 90, cex = sqrt(yPadUp), hjust = 0, vjust = 0)
        if (nrow(h_labs) < 250) {
          h_labs$label3 <- gsub(" - .*", "", h_labs$label)
          heatmap.plot <- heatmap.plot +
            geom_text(data = h_labs, aes(x = x-0.3, label = label3, colour = Cluster),
                      y = Height+0.05, angle = 90, cex = sqrt(yPadUp)*0.9, hjust = 0, vjust = 0)
        }
        # - Vertical dendrogram and sample names
        heatmap.plot <- heatmap.plot +
          geom_segment(data = v_Seg, linewidth = 0.25,
                       aes(x = x, y = y - 0.5, xend = xend, yend = yend - 0.5)) +
          new_scale_color() + VcolScale +
          geom_text(data = v_labs, aes(y = y - 0.5, label = label, colour = Cluster),
                    x = -xScale/250, hjust = 1, vjust = 0.5, cex = 2.2)
        # - Proteins of interest
        if (prot.list.marks) {
          heatmap.plot <- heatmap.plot +
            geom_point(data = temp2c, aes(x = Xmin+0.5), y = -0.5, colour = "red", fill = "red", shape = 17) +
            geom_text(data = temp2c, aes(x = Xmin+0.5, label = label2),
                      y = -1, colour = "red", angle = -60, hjust = 0, cex = 2)
        }
        # - Subcellular markers
        addSCmarks <- FALSE
        if (("Compartment marker" %in% colnames(PG))&&(sum(PG$"Compartment marker" != ""))) {
          temp2m <- as.data.table(temp2[, c("Leading protein IDs", "Ymin", "Xmin")])
          temp2m <- temp2m[, list(Xmin = min(Xmin, na.rm = TRUE),
                                  Ymin = min(Ymin, na.rm = TRUE)), by = list(`Leading protein IDs` = `Leading protein IDs`)]
          temp2m <- as.data.frame(temp2m)
          m <- match(temp2m$`Leading protein IDs`, PG$`Leading protein IDs`)
          temp2m$Label <- PG$Label[m]
          temp2m$"Compartment marker" <- PG$"Compartment marker"[m]
          temp2m <- temp2m[which(temp2m$`Compartment marker` != ""),]
          addSCmarks <- nrow(temp2m)
          if (addSCmarks) {
            heatmap.plot <- heatmap.plot +
              new_scale_fill() +
              #scale_fill_viridis("turbo", discrete = TRUE) +
              scale_color_manual(values = markColors) +
              geom_rect(data = temp2m, aes(xmin = Xmin, xmax = Xmin+1, fill = `Compartment marker`),
                         ymin = -1.7, ymax = -1.3)
          }
        }
        #poplot(heatmap.plot, 12, 20)
        #
        # Export evaluated version for later parallel saving
        evalPlot <- proteoCraft::plotEval(heatmap.plot)
        heatMaps[[paste0(i, " - ", normType, ".jpeg")]] <- list(Plot = evalPlot,
                                                                Ttl = paste0(clustDir, "/", nm, normTypeInsrt, ".jpeg"),
                                                                Width = 18,
                                                                Height = 9,
                                                                Units = "in")
        heatMaps[[paste0(i, " - ", normType, ".pdf")]] <- list(Plot = evalPlot,
                                                               Ttl = paste0(clustDir, "/", nm, normTypeInsrt, ".pdf"))                                            
        #
        if (drawPlotly) {
          # Plotly version
          tempLy <- temp2a[w2a,]
          tempLy$Sample <- factor(tempLy$Sample, levels = smpls)
          ##tempLy$Ymin <- tempLy$Ymin+0.5
          plotleatmap <- plot_ly(data = tempLy, x = ~Xmin, y = ~Ymin, z = ~value, type = "heatmap",
                                 hovertext = tempLy$Label, hoverinfo = "text")
          # I cannot find a way to remove tick marks!!!!!
          # Vertical dendrogram
          plLyV <- v_Seg
          ##plLyV$x <- -Width*0.2-plLyV$x 
          ##plLyV$xend <- -Width*0.2-plLyV$xend
          plLyV$y <- plLyV$y-1
          plLyV$yend <- plLyV$yend-1
          plotleatmap <- add_segments(plotleatmap, x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                                      data = plLyV, inherit = FALSE, color = I("black"),
                                      showlegend = FALSE, hoverinfo = "none")
          # Horizontal dendrogram
          plLyH <- h_Seg
          ##plLyH$x <- plLyH$x-0.5
          ##plLyH$xend <- plLyH$xend-0.5
          plLyH$y <- plLyH$y - hnch - 1.5
          plLyH$yend <- plLyH$yend - hnch - 1.5
          plotleatmap <- add_segments(plotleatmap, x = ~x-0.5, xend = ~xend-0.5, y = ~y, yend = ~yend,
                                      data = plLyH, inherit = FALSE, color = I("black"),
                                      showlegend = FALSE, hoverinfo = "none")
          # Sample labels
          v_labs2 <- v_labs
          colnames(v_labs2)[which(colnames(v_labs2) == "label")] <- "Sample"
          gap <- abs(max(plLyV$x))
          plotleatmap <- layout(plotleatmap,
                                shapes = list(list(
                                  type = "rect",
                                  fillcolor = "grey",
                                  line = list(color = "darkgrey"),
                                  opacity = 0.3,
                                  x0 = -gap+0.01,
                                  x1 = -0.51,
                                  xref = "x",
                                  y0 = -0.5,
                                  y1 = max(v_labs2$y)-0.5,
                                  yref = "y"
                                )))
          if (scrptType == "withReps") {
            allPal <- c("magma", "inferno", "plasma", "rocket", "turbo", "mako", "cividis")
            nPal <- length(allPal)
            m <- match(v_labs2$Sample, clustMap$Samples)
            v_labs2[, RSA$names] <- clustMap[m, RSA$names]
            tst <- vapply(RSA$names, function(x) { length(unique(v_labs2[[x]])) }, 1) > 1
            w <- which(tst)
            l <- length(w)+1
            sq <- seq_len(l)
            fctNms <- data.frame(Factor = c("Sample", RSA$names[w]),
                                 x = - gap + gap*sq/(l+1))
            for (j in sq) { #j <- 1 #j <- j + 1
              fct <- fctNms$Factor[j]
              vals <- v_labs2[[fct]]
              fvals <- if (is.factor(vals)) { vals } else { factor(vals, levels = unique(vals)) }
              # choose palette function reliably
              pal_nm <- allPal[(j - 1) %% nPal + 1]
              pal_fun  <- get(pal_nm, mode = "function", envir = asNamespace("viridisLite"))
              lvl_nms <- levels(fvals)
              n_lvls <- length(lvl_nms)
              pal_cols <- pal_fun(n_lvls)    # returns character hex vector
              # Map colors to levels (use levels() to ensure stable mapping)
              col_map <- setNames(pal_cols, lvl_nms)
              # Build per-row color column (character vector of hex colors)
              myCol <- col_map[as.character(fvals)]
              # Attach these as columns to a temporary copy of the data (preserve v_labs2)
              vtmp <- v_labs2
              vtmp$.vals_for_hover <- as.character(fvals)  # hover text as character
              vtmp$.col_for_plot <- myCol                 # hex colors per row
              # Sanity checks (print a small table to debug)
              # print(vtmp[, c("Sample", fct, ".vals_for_hover", ".col_for_plot")])
              plotleatmap <- suppressWarnings(
                add_markers(plotleatmap,
                            data = vtmp,
                            x = fctNms$x[j],
                            y = ~y - 1,
                            inherit = FALSE,
                            showlegend = FALSE,
                            text = ~.vals_for_hover,      # hover text bound to data rows
                            hoverinfo = "text",
                            name = fct,
                            marker = list(color = vtmp$.col_for_plot, # explicit vector aligned to vtmp rows
                                          size = 10, showscale = FALSE))
              )
            }
            mxY <- max(v_labs2$y)
            plotleatmap <- layout(plotleatmap,
                                  annotations = lapply(seq_len(nrow(fctNms)), function(j) {
                                    list(x = fctNms$x[j]-0.5,
                                         y = mxY-0.5,
                                         text = fctNms$Factor[j],
                                         textangle = -60,
                                         showarrow = FALSE,
                                         xanchor = "left",
                                         yanchor = "bottom")
                                  }))
          }
          if (scrptType == "noReps") {
            plotleatmap <- add_trace(plotleatmap, data = v_labs2, y = ~y-1, x = -gap/10, text = ~Sample,
                                     color = I("black"), inherit = FALSE, type = "scatter",
                                     mode = "text", showlegend = FALSE,
                                     textposition = "middle left") # ... which means what I would take to be middle right...
          }
          # Cluster segments
          # uCl <- unique(clustDF$Cluster)
          # segmCol <- setNames(rainbow(length(uCl)), uCl)
          # clustDF$Col <- segmCol[clustDF$Cluster]
          # for (j in seq_len(nrow(clustDF))) {
          #   xs <- c(clustDF$x[j]-1.5, clustDF$x[j]-0.5, clustDF$x[j]-0.5, clustDF$x[j]-1.5, clustDF$x[j]-1.5)
          #   ys <- c(-1.5, -1.5, -1, -1, -1.5)
          #   lbl <- rep(clustDF$Cluster[j], length(xs))
          #   plotleatmap <- add_trace(plotleatmap, x = xs, y = ys, type = "scatter", mode = "lines",
          #                            fill = "toself", fillcolor = clustDF$Col[j],
          #                            line = list(width = 0, color = clustDF$Col[j]),
          #                            text = clustDF$Cluster[j],
          #                            hovertemplate = paste0(clustDF$Cluster[j], "<extra></extra>"),
          #                            hoverinfo = "skip", hovertext = clustDF$Cluster[j],
          #                            showlegend = FALSE, inherit = FALSE)
          # }
          plotleatmap <- add_trace(plotleatmap, data = clustDF, x = ~x-1,
                                   y = -1.5, inherit = FALSE, showlegend = FALSE,
                                   hovertext = ~Cluster, hoverinfo = "text",
                                   type = "scatter", mode = "text",
                                   color = clustDF$Cluster,
                                   marker = list(size = 5, showscale = FALSE))
          # Layout
          plotleatmap <- layout(plotleatmap, title = nm,
                                xaxis = list(title = "Protein groups",
                                             tickmode = "array", tickvals = NULL, showticklabels = FALSE,
                                             showgrid = FALSE, zeroline = FALSE),
                                yaxis = list(title = "Samples",
                                             tickmode = "array", tickvals = NULL, showticklabels = FALSE,
                                             showgrid = FALSE, zeroline = FALSE))
          if (prot.list.marks) {
            plotleatmap <- add_trace(plotleatmap, data = temp2c, x = ~Xmin, y = -1,
                                     text = ~Label, color = "red", inherit = FALSE,
                                     type = "scatter", mode = "markers", showlegend = FALSE,
                                     marker = list(size = 5, showscale = FALSE))
          }
          if (addSCmarks) {
            temp2m$.markCol <- markColors[temp2m$`Compartment marker`]
            temp2m$.size <- 5
            plotleatmap <- add_trace(plotleatmap, data = temp2m, x = ~Xmin, y = -2,
                                     text = ~temp2m$`Compartment marker`,
                                     inherit = FALSE,
                                     type = "scatter", mode = "markers", showlegend = FALSE,
                                     marker = list(color = temp2m$.markCol,
                                                   size = temp2m$.size, showscale = FALSE, showlegend = FALSE))
          }
          # setwd(clustDir)
          # saveWidget(plotleatmapa, paste0(clustDir, "/", nm, normTypeInsrt, ".html"))
          # setwd(wd)
          #system(paste0("open \"", clustDir, "/", nm, normTypeInsrt, ".html\""))
          if (clustMode == "standard") {
            plotLeatMaps[[i]][[normType]] <- list(Ttl = paste0(nm, normTypeInsrt),
                                                  Plot = plotleatmap,
                                                  Render = plotly::plotly_build(plotleatmap))
          }
        }
        #poplot(heatmap.plot, 12, 20)
      } else {
        warning(paste0(i, ", ", normType, ": less than 3 rows, skipping..."))
      }
    }
  }
  #
  # Save ggplots
  invisible(parLapply(parClust, heatMaps, function(x) {
    if (grepl("\\.pdf$", x$Ttl)) {
      ggplot2::ggsave(x$Ttl, x$Plot)
    }
    if (grepl("\\.jpeg$", x$Ttl)) {
      ggplot2::ggsave(x$Ttl, x$Plot, width = x$Width, height = x$Height, units = x$Units)
    }
    return()
  }))
  #
  if (clustMode == "standard") {
    saveFun(plotLeatMaps, file = paste0(clustDir, "/HeatMaps.RData"))
    # Save plotly plots
    dr <- clustDir
    myPlotLys <- list()
    for (nm in names(plotLeatMaps)) {
      myPlotLys[paste0(nm, " - ", names(plotLeatMaps[[nm]]))] <- plotLeatMaps[[nm]]
    }
    Src <- paste0(libPath, "/extdata/R scripts/Sources/save_Plotlys.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
  #
  kol <- c("Leading protein IDs", "Protein names", "Genes", "Mol. weight [kDa]")
  kol <- kol[which(kol %in% colnames(PG))]
  temp <- PG[, c(kol, KlustKols)]
  flPath <- paste0(clustDir, "/Protein Groups and Clusters.csv")
  tst <- try(write.csv(temp, file = flPath, row.names = FALSE), silent = TRUE)
  while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
    dlg_message(paste0("File \"", flPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
    tst <- try(write.csv(temp, file = flPath, row.names = FALSE), silent = TRUE)
  }
}
