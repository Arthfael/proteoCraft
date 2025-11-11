# Prepare data for clustering heatmaps script
#
# This script prepares some data which will be re-used when we run the cluster_Heatmap source
if (!exists("clustPrep")) { clustPrep <- FALSE }
clustPrep %<o% clustPrep
if (!clustPrep) {
  if (scrptType == "withReps") {
    clustHtMp %<o% TRUE
    if (LocAnalysis) { prtRfRoot <- Prot.Expr.Root2 } else { prtRfRoot <- Prot.Expr.Root }
  }
  if (scrptType == "noReps") {
    clustHtMp <- (length(Exp) > 1)
    prtRfRoot <- PG.int.col
  }
  clustHtMp %<o% clustHtMp
  prtRfRoot %<o% prtRfRoot
  if (clustHtMp) {
    ImputeKlust %<o% TRUE # Currently always TRUE
    if (scrptType == "withReps") {
      clustXprsKol <- paste0(prtRfRoot, Exp.map$Ref.Sample.Aggregate)
      clustSmpls <- which(clustXprsKol %in% colnames(PG))
      clustXprsKol <- clustXprsKol[clustSmpls]
      clustMap <- Exp.map[clustSmpls,]
      clustMap$Samples <- cleanNms(clustMap$Ref.Sample.Aggregate)
    }
    if (scrptType == "noReps") {
      clustXprsKol <- paste0(prtRfRoot, SamplesMap$Experiment)
      clustSmpls <- which(clustXprsKol %in% colnames(PG))
      clustXprsKol <- clustXprsKol[clustSmpls]
      clustMap <- SamplesMap[clustSmpls,]
      clustMap$Samples <- cleanNms(clustMap$Experiment)
    }
    mySmpls %<o% clustMap$Samples
    clustXprsKol %<o% clustXprsKol
    clustSmpls %<o% clustSmpls
    clustMap %<o% clustMap
    #
    clustFilt %<o% which((apply(PG[, clustXprsKol], 1, function(x) { sum(!is.na(x)) }) > 0)
                         &((is.na(PG$`Potential contaminant`))|(PG$`Potential contaminant` != "+")))
    clustDat %<o% set_rownames(set_colnames(PG[clustFilt, clustXprsKol], clustMap$Samples),
                               PG$Label[clustFilt])
    clustDat0 %<o% clustDat
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
      clustDatImp %<o%  clustDat$Positions_Imputed
      clustDat <- clustDat$Imputed_data
      rownames(clustDatImp) <- rownames(clustDat)
      colnames(clustDatImp) <- clustMap$Samples
      #vapply(clustMap$Samples, function(x) { sum(is.na(clustDat[[x]])) }, 1) # df
      #vapply(clustMap$Samples, function(x) { sum(clustDatImp[, x]) }, 1) # matrix
    }
  }
  #
  KlustRoot %<o% paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ") - ")
  plotLeatMaps %<o% list()
  Heatmaps %<o% list()
  NHClust %<o% list()
  NVClust %<o% list()
  KlustKols %<o% c()
  MaxHClust %<o% min(c(floor(nrow(PG)/2), 100)) # We want at most 20 clusters
  if (scrptType == "withReps") {
    MaxVClust <- length(VPAL$values)
  }
  if (scrptType == "noReps") {
    MaxVClust <- length(Exp)
  }
  MaxVClust %<o% MaxVClust
  VClusters %<o% list()
  #
  # Great function built with the help of chatGPT to extract segments and labels from an hclust efficiently:
  hclust_to_seg %<o% function(hc) {
    stopifnot(inherits(hc, "hclust"))
    #
    n <- length(hc$order)
    merge <- hc$merge
    height <- hc$height
    n_merge <- nrow(merge)
    n_nodes <- n + n_merge
    #
    # Allocate
    x <- numeric(n_nodes)
    y <- numeric(n_nodes)
    #
    # Leaves
    x[hc$order] <- seq_len(n)
    y[1:n] <- 0
    #
    # Precompute absolute indices
    left_idx  <- ifelse(merge[, 1] < 0, -merge[, 1], n + merge[, 1])
    right_idx <- ifelse(merge[, 2] < 0, -merge[, 2], n + merge[, 2])
    #
    # Compute internal node positions (iterative dependency)
    i_seq <- seq_len(n_merge)
    y[n + i_seq] <- hc$height
    for (i in i_seq) {
      x[n + i] <- (x[left_idx[i]] + x[right_idx[i]]) / 2
      y[n + i] <- height[i]
    }
    #
    # ---- build segments ----
    # For each merge, we add 3 segments: left vertical, right vertical, top horizontal
    parent_idx <- n + i_seq
    segs <- data.frame(x = numeric(3 * n_merge),
                       y = numeric(3 * n_merge),
                       xend = numeric(3 * n_merge),
                       yend = numeric(3 * n_merge))
    # vertical left
    segs$x[3*i_seq - 2] <- x[left_idx]
    segs$y[3*i_seq - 2] <- y[left_idx]
    segs$xend[3*i_seq - 2] <- x[left_idx]
    segs$yend[3*i_seq - 2] <- y[parent_idx]
    # vertical right
    segs$x[3*i_seq - 1] <- x[right_idx]
    segs$y[3*i_seq - 1] <- y[right_idx]
    segs$xend[3*i_seq - 1] <- x[right_idx]
    segs$yend[3*i_seq - 1] <- y[parent_idx]
    # horizontal top
    segs$x[3*i_seq]    <- x[left_idx]
    segs$y[3*i_seq]    <- y[parent_idx]
    segs$xend[3*i_seq] <- x[right_idx]
    segs$yend[3*i_seq] <- y[parent_idx]
    # Labels
    labs <- data.frame(x = x[hc$order],
                       y = y[hc$order],
                       label = hc$labels[hc$order])
    #
    list(segments = segs,
         labels = labs)
  }
  #
  clustPrep <- TRUE
}
