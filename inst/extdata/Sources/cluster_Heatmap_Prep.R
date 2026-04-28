# Prepare data for clustering heatmaps script
#
ImputeKlust %<o% TRUE # Currently MUST always be TRUE
#
if ((!exists("clustDat"))||(!inherits(clustDat, "list"))) { clustDat <- list() }
if ((!exists("clustFilt"))||(!inherits(clustFilt, "list"))) { clustFilt <- list() }
if ((!exists("plotLeatMaps"))||(!inherits(plotLeatMaps, "list"))) { plotLeatMaps <- list() }
if ((!exists("Heatmaps"))||(!inherits(Heatmaps, "list"))) { Heatmaps <- list() }
clustDat %<o% clustDat
clustFilt %<o% clustFilt
plotLeatMaps %<o% plotLeatMaps
Heatmaps %<o% Heatmaps
#
if (scrptType == "withReps") { clustHtMp <- TRUE }
if (scrptType == "noReps") { clustHtMp <- length(Exp) > 1L }
clustHtMp %<o% clustHtMp
if (dataType == "PG") {
  if (scrptType == "withReps") {
    rfRoot <- if (LocAnalysis) { Prot.Expr.Root2 } else { Prot.Expr.Root }
    MaxVClust <- length(VPAL$values)
  }
  if (scrptType == "noReps") {
    rfRoot <- PG.int.col
    MaxVClust <- length(Exp)
  }
  myData <- PG
  rownames(myData) <- myData$Label
}
if (dataType == "peptides") { # Currently only used to prepare data for dim-red plots, not clustering heatmaps!
  # ... but this will likely change
  if (scrptType == "withReps") {
    rfRoot <- pep.ref[length(pep.ref)]
  }
  if (scrptType == "noReps") {
    stop("Not written yet!")
  }
  myData <- pep
  rownames(myData) <- pep$id
}
rfRoot %<o% rfRoot
#
if (scrptType == "withReps") {
  clustXprsKol <- paste0(rfRoot, Exp.map$Ref.Sample.Aggregate)
  clustSmpls <- which(clustXprsKol %in% colnames(myData))
  clustXprsKol <- clustXprsKol[clustSmpls]
  map <- Exp.map[clustSmpls,]
  map$Samples <- cleanNms(map$Ref.Sample.Aggregate)
}
if (scrptType == "noReps") {
  clustXprsKol <- paste0(rfRoot, SamplesMap$Experiment)
  clustSmpls <- which(clustXprsKol %in% colnames(myData))
  clustXprsKol <- clustXprsKol[clustSmpls]
  map <- SamplesMap[clustSmpls,]
  map$Samples <- map$Experiment
}
if (dataType == "PG") {
  MaxVClust %<o% MaxVClust
  clustMap %<o% map
} # (otherwise - for peptides - we re-use those from PG)


# Original data
# -------------
w <- which((apply(myData[, clustXprsKol, drop = FALSE], 1L, \(x) { sum(!is.na(x)) }) > 0L)
           &((is.na(myData$`Potential contaminant`))|(myData$`Potential contaminant` != "+")))
myData <- set_colnames(myData[w, clustXprsKol, drop = FALSE], map$Samples)

filt <- rownames(myData)[which(apply(myData[, map$Samples], 1L, \(x) { length(is.all.good(x)) }) > 0L)]
clustFilt[[dataType]] <- filt
clustDat[[dataType]] <- list()
clustDat[[dataType]]$Original <- myData
#clustDat[[dataType]]$Original -> myData
# Do not update clustFilt within the ImputeKlust loop below!
# The idea is that, yes, we may impute for clustering, but we filter based on pre-imputed data!


# Imputed data
# --------------
# We will always impute here, whether we decide to use imputed data or not: we are preparing data
if (scrptType == "withReps") { # Here we have replicates and group at samples group level
  Grps <- map[[VPAL$column]]
}
if (scrptType == "noReps") { # Here we do not have replicates and group at comparison group level
  Grps <- map$`Ratios group`
}
Grps <- setNames(match(Grps, unique(Grps)), Grps)
tmp <- Data_Impute2(myData, Grps)
clustDat[[dataType]]$Imputed <- myData <- tmp$Imputed_data
myDataImp <- tmp$Positions_Imputed
rownames(myDataImp) <- rownames(myData)
colnames(myDataImp) <- colnames(myData)
clustDat[[dataType]]$Positions_imputed <- myDataImp

# Batch-correct
# -------------
# Sometimes we have identified a batch, but decided no to correct for it in the data we use for statistical tests
# (e.g. since limma can encode it in its design matrix).
# In that case, we may optionally still remove it from the PCA.
if ((scrptType == "withReps")&&(Param$Batch.effect != "")) { # Here we have replicates and group at samples group level
  nms <- Batch.effect$names
  w <- which(vapply(normSequence, \(x) { x$Method }, "") == "ComBat")
  if (length(w)) {
    btchs <- lapply(normSequence[w], \(x) { x$Batch })
    w <- w[which(vapply(btchs, \(x) { sum(nms %in% x) }, 1L) == length(nms))]
  }
  if (length(w)) {
    dcs <- vapply(normSequence[w], \(x) { x$Decision }, TRUE)
    runComBatNow <- sum(!dcs) > 0L
  } else {
    runComBatNow <- TRUE
  }
  if (runComBatNow) {
    m <- match(colnames(myData), map$Samples)
    btchs <- factor(map[m, Batch.effect$column])
    mdlMtr <- designMatr_noBatch[match(gsub("___", "_", as.character(expMap[m, RSA$limmaCol])), row.names(designMatr_noBatch)),]
    #mdlMtr <- model.matrix(~1, data = Exp.map[match(currSamples, Exp.map$Ref.Sample.Aggregate),])
    if (dataType == "PG") {
      tmp <- ComBat(dat = myData,
                    batch = btchs,
                    mod = mdlMtr,
                    par.prior = TRUE)
    }
    if (dataType == "peptides") {
      tmp <- lapply(NormGrps$Group, \(lGrp) { #lGrp <- NormGrps$Group[1L] #lGrp <- NormGrps$Group[2L]
        w <- which(rownames(myData) %in% NormGrps$IDs[[match(lGrp, NormGrps$Group)]])
        if (!length(w)) { return() }
        #
        # For ComBat we only use the longitudinal groups (peptide normalisation group),
        # not the transversal groups (comparison/ratio groups).
        # Indeed, batches will often intersect with the latter
        rs <- ComBat(dat = myData[w,],
                     batch = btchs,
                     mod = mdlMtr,
                     par.prior = TRUE)
        rownames(rs) <- rownames(myData)[w]
        return(rs)
      })
      tmp <- do.call(rbind, tmp)
      tmp <- as.data.frame(tmp[match(rownames(myData), rownames(tmp)),])
    }
    rownames(tmp) <- rownames(myData)
    myDataCorr <- tmp
    clustDat[[dataType]]$ComBat <- myDataCorr
  }
}
#
if (dataType == "PG") {
  KlustRoot %<o% paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ") - ")
  NHClust %<o% list()
  NVClust %<o% list()
  KlustKols %<o% c()
  MaxHClust %<o% min(c(as.integer(floor(nrow(PG)/2)), 100L)) # We want at most 20 clusters
  VClusters %<o% list()
  HClusters %<o% list()
}
