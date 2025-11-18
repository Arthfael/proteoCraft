# Prepare data for clustering heatmaps script
#
# This script prepares some data which will be re-used when we run the cluster_Heatmap source
if ((!exists("clustPrep"))||(!is.logical(clustPrep))||(length(clustPrep) != 1)||(is.na(clustPrep))) { clustPrep <- FALSE }
clustPrep %<o% clustPrep
ImputeKlust %<o% TRUE # Currently MUST always be TRUE
#
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
  #
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
  clustPrep <- TRUE
}
