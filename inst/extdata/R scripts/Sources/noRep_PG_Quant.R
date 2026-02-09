# Protein group quantitation
# --------------------------
#
# This script calculates protein group-level quantitative values based on those of individual peptides.
# Exact quantitation parameters (e.g. summarization method) are defined earlier in the workflow.

# You may want to exclude peptides based on how many samples they are found in:
kol <- paste0(int.col, " - ", Exp)
if (!exists("PepFoundInAtLeast")) { PepFoundInAtLeast <- 1 }
if ("PepFoundInAtLeast" %in% names(AnalysisParam)) {
  tmp <- suppressWarnings(as.integer(AnalysisParam$PepFoundInAtLeast))
  if (validIntegPar("tmp")) {
    PepFoundInAtLeast <- tmp
  }
}
PepFoundInAtLeast %<o% PepFoundInAtLeast
tst <- parApply(parClust, pep[, kol, drop = FALSE], 1, function(x) { sum(proteoCraft::is.all.good(x) > 0) })
Pep2Use %<o% which(tst >= PepFoundInAtLeast)

# Map
smplsMap <- SamplesMap
if (MakeRatios) {
  smplsMap$Ratios_group <- paste0("Group", smplsMap$"Ratios group")
  RatGrp <- RefGrp <- list(aggregate = "Rat",
                           values = unique(smplsMap$Ratios_group),
                           names = "Ratios_group",
                           column = "Ratios_group")
  SmplGrp <- list(aggregate = "Exp",
                  values = Exp,
                  names = "Experiment",
                  column = "Experiment")
  Aggregate.map <- list(Aggregate.name = c("Exp", "Rat"), Characteristics = c("Experiment", "Ratios_group"))
  Aggregate.list <- list(Rat = unique(smplsMap$Ratios_group),
                         Exp = unique(smplsMap$Experiment))
  Aggregates <- setNames(c("Ratios_group", "Experiment"), c("Rat", "Exp"))
} else {
  pepRat_col <- Aggregate.map <- Aggregate.list <- Aggregates <- NA
}

# Weights:
# - Higher for peptides with low PEP
# - CV not taken into account since we don't know what a good (regulatory) or a bad CV are, and we may even only have 1 peptide
pep$Weights <- -log10(pep$PEP)
weightsInsrt <- "-log(PEP)"
if ((length(inDirs) == 1)&&(QuantUMS)&&("Quantity Quality" %in% colnames(ev))) {  # DiaNN specific:
  if (!"Quantity Quality" %in% colnames(pep)) {
    tmp <- data.table(Qual = ev$"Quantity Quality", ModSeq = ev$"Modified sequence")
    tmp <- tmp[, .(Qual= mean(Qual)), by = .(ModSeq)]
    tmp <- as.data.frame(tmp)
    pep$"Quantity Quality" <- tmp$Qual[match(pep$"Modified sequence", tmp$ModSeq)]
  }
  weightsInsrt <- "-log10(PEP * (1 - DiaNN_Quantity_Quality))" 
  pep$Weights <- -log10(pep$PEP*(1-pep$"Quantity Quality"))
  # plot <- ggplot(pep) +
  #   geom_density(stat = "density", aes(x = PEP), colour = "red") +
  #   geom_density(stat = "density", aes(x = 1-`Quantity Quality`), colour = "blue") +
  #   theme_bw()
  # poplot(plot)
  # plot <- ggplot(pep) +
  #   geom_density(stat = "density", aes(x = -log10(PEP)), colour = "red") +
  #   geom_density(stat = "density", aes(x = -log10(1-`Quantity Quality`)), colour = "blue") +
  #   theme_bw()
  # poplot(plot)
}
w <- which(!is.all.good(pep$Weights, 2))
pep$Weights[w] <- min(pep$Weights, na.rm = TRUE)
#
m <- max(pep$Weights)
pep$Weights <- pep$Weights/m
pep$Weights[which(pep$Weights < 0.001)] <- 0.001
summary(pep$Weights)
#
Mod.Excl.is.strict %<o% FALSE
if ("Prot.Quant.Mod.Excl.is.strict" %in% names(AnalysisParam)) {
  tmp <- suppressWarnings(as.logical(AnalysisParam$Prot.Quant.Mod.Excl.is.strict))
  if (validLogicPar("tmp")) { Mod.Excl.is.strict <- tmp }
}
Discard.unmod %<o% (Mod.Excl.is.strict+1)
if (Discard.unmod == 1) { Discard.unmod <- as.logical(Discard.unmod) }

#
PG.int.col %<o% "log10(Expr.) - "
PG.int.cols %<o% setNames(PG.int.col, "Original")
PG.rat.col %<o% "log2(Ratio) - " #We need those defaults actually even if !MakeRatios
PG.rat.cols %<o% setNames(PG.rat.col, "Original")
pepInt_col <- paste0(int.cols[max(which(names(int.cols) != "Imputed"))], " - ")
if (MakeRatios) {
  pepRat_col <- paste0(rat.cols[max(which(names(rat.cols) != "Imputed"))], " - ")
} else { pepRat_col <- NA }

# Quantitation
source(parSrc, local = FALSE)
reScalingAlgo <- reScAlgo
if (reScAlgo == "topN") { reScalingAlgo <- paste0("top", topN) }
quantData_list %<o% protQuant(PG,
                              pep[Pep2Use,],
                              Pep4Quant,
                              "Unique peptide IDs",
                              "id",
                              N_unique_Pep,
                              quantAlgo,
                              reScalingAlgo,
                              topN_correct,
                              N_Pep,
                              Inf,
                              Weights = "Weights",
                              useIntWeights = FALSE,
                              Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              skipRatios = !MakeRatios,
                              expMap = smplsMap,
                              expMap_Samples_col = "Experiment",
                              pepInt_Root = pepInt_col,
                              pepRat_root = pepRat_col,
                              pepInt_log = FALSE,
                              pepRat_log = 2,
                              protLFQ_toLog = TRUE,
                              protRat_toLog = TRUE,
                              mods_to_Exclude = Mod2Xclud,
                              mod_Seq = "Modified sequence",
                              discard_unmod = Discard.unmod,
                              prim_Seq = "Sequence (1st accession)",
                              cl = parClust,
                              N.clust = N.clust,
                              N.reserved = 1)
saveFun(quantData_list, file = "quantData.RData")
#loadFun("quantData.RData")
quant.data <- quantData_list$Data
quant.data <- quant.data[, which(!grepl("\\.REF$", colnames(quant.data))), drop = FALSE]
colnames(quant.data) <- gsub("^log10\\(Expr\\.\\) - ", PG.int.col, colnames(quant.data))
if (MakeRatios) {
  colnames(quant.data) <- gsub("^log2\\(Ratio\\) - ", PG.rat.col, colnames(quant.data))
}
w <- which(colnames(quant.data) == "Peptides IDs used for quantitation")
colnames(quant.data)[w] <- paste0("Peptide IDs used for quantitation - ",
                                  names(PG.int.cols)[match(PG.int.col, PG.int.cols)])
PG[, colnames(quant.data)] <- quant.data

if (Impute) {
  # In that case we need to calculate expression values a second time, unfortunately... It doesn't take so long.
  PG.int.cols["Imputed"] <- PG.int.col <- paste0("Imput. ", PG.int.cols["Original"])
  PG.rat.cols["Imputed"] <- PG.rat.col <- paste0("Imput. ", PG.rat.cols["Original"])
  pepImpInt_col <- paste0(int.cols["Imputed"], " - ")
  pepImpRat_col <- paste0(rat.cols["Imputed"], " - ")
  quantData_list_Imput %<o% protQuant(PG,
                                      pep[Pep2Use,],
                                      Pep4Quant,
                                      "Unique peptide IDs",
                                      "id",
                                      N_unique_Pep,
                                      quantAlgo,
                                      reScalingAlgo,
                                      topN_correct,
                                      N_Pep,
                                      Inf,
                                      Weights = "Weights",
                                      useIntWeights = FALSE,
                                      Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                                      skipRatios = !MakeRatios,
                                      expMap = smplsMap,
                                      expMap_Samples_col = "Experiment",
                                      pepInt_Root = pepImpInt_col,
                                      pepRat_root = pepImpRat_col,
                                      pepInt_log = FALSE,
                                      pepRat_log = 2,
                                      protLFQ_toLog = TRUE,
                                      protRat_toLog = TRUE,
                                      mods_to_Exclude = Mod2Xclud,
                                      mod_Seq = "Modified sequence",
                                      discard_unmod = Discard.unmod,
                                      prim_Seq = "Sequence (1st accession)",
                                      cl = parClust,
                                      N.clust = N.clust,
                                      N.reserved = 1)
  saveFun(quantData_list_Imput, file = "quantData_Imput.RData")
  #loadFun("quantData_Imput.RData")
  quant.data2 <- quantData_list_Imput$Data
  quant.data2 <- quant.data2[, which(!grepl("\\.REF$", colnames(quant.data2))), drop = FALSE]
  colnames(quant.data2) <- gsub("^log10\\(Expr\\.\\) - ", PG.int.col, colnames(quant.data2))
  if (MakeRatios) {
    colnames(quant.data2) <- gsub("^log2\\(Ratio\\) - ", PG.rat.col, colnames(quant.data2))
  }
  w2 <- which(colnames(quant.data2) == "Peptides IDs used for quantitation")
  colnames(quant.data2)[w2] <- "Peptide IDs used for quantitation - Imputed"
  PG[, colnames(quant.data2)] <- quant.data2
}
#
rm(Exp.map)
