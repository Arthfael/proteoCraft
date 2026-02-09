# Protein group quantitation
# --------------------------
#
# This script calculates protein group-level quantitative values based on those of individual peptides.
# Exact quantitation parameters (e.g. summarization method) are defined earlier in the workflow.

# You may want to exclude peptides based on how many samples they are found in:
kol <- lapply(VPAL$values, function(x) {
  x <- paste0(pep.ref["Original"], Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)])
  return(x[which(x %in% colnames(pep))])
})
if (!exists("PepFoundInAtLeast")) { PepFoundInAtLeast %<o% 1 }
if ("PepFoundInAtLeast" %in% colnames(Param)) {
  PepFoundInAtLeast %<o% suppressWarnings(as.integer(Param$PepFoundInAtLeast))
  if ((is.na(PepFoundInAtLeast))||(PepFoundInAtLeast < 1)) {
    warning("Invalid \"PepFoundInAtLeast\" parameter, defaulting to 1")
    PepFoundInAtLeast <- 1
  }
}
tst <- parApply(parClust, pep[, unlist(kol)], 1, function(x) { sum(proteoCraft::is.all.good(x) > 0) })
Pep2Use %<o% which(tst >= PepFoundInAtLeast)
maxAllowed <- max(c(2, length(Rep)-1))
if (!exists("PepFoundInAtLeastGrp")) { PepFoundInAtLeastGrp %<o% maxAllowed }
if ("PepFoundInAtLeastGrp" %in% colnames(Param)) {
  PepFoundInAtLeastGrp %<o% suppressWarnings(as.integer(Param$PepFoundInAtLeastGrp))
  if ((is.na(PepFoundInAtLeastGrp))||(PepFoundInAtLeastGrp < 1)||(PepFoundInAtLeastGrp > length(Rep)-1)) {
    warning(paste0("Invalid \"PepFoundInAtLeastGrp\" parameter, defaulting to ", maxAllowed))
    PepFoundInAtLeastGrp <- maxAllowed
  }
}
tst <- pep[, unlist(kol)]
clusterExport(parClust, list("tst", "kol", "PepFoundInAtLeastGrp"), envir = environment())
tst <- parSapply(parClust, 1:nrow(pep), function(x) {
  x <- max(vapply(kol, function(kl) { sum(proteoCraft::is.all.good(unlist(tst[x, kl])) > 0) }, 1) >= PepFoundInAtLeastGrp)
  return(x)
}) > 0
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Pep2Use %<o% which(tst)
#length(Pep2Use)/nrow(pep)

# Weights:
# - Higher for peptides with low intra-sample group CV on average
# - Higher for peptides with low PEP
source(parSrc, local = FALSE)
Kols <- paste0(pep.ref[length(pep.ref)], Exp.map$Ref.Sample.Aggregate)
Kols <- Kols[which(Kols %in% colnames(pep))]
tmp <- pep[, Kols]
clusterExport(parClust, list("Exp.map", "VPAL", "pep.ref", "is.all.good", "tmp"), envir = environment())
CV <- parSapply(parClust, VPAL$values, function(x) { #x <- VPAL$values[1]
  smpls <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)]
  kols <- paste0(pep.ref[length(pep.ref)], smpls)
  kols <- kols[which(kols %in% colnames(tmp))]
  x <- apply(tmp[, kols, drop = FALSE], 1, function(y) {
    y <- is.all.good(log10(unlist(y)))
    if (length(y)) {
      y <- sd(y)/mean(y)
    } else { y <- NA }
    return(y)
  })
  return(x)
})
pep$CV <- rowMeans(CV, na.rm = TRUE)
pep$Weights <- -log10(pep$PEP*pep$CV)
weightsInsrt <- "-log(PEP * CV)"
if ((length(inDirs) == 1)&&(QuantUMS)&&("Quantity Quality" %in% colnames(ev))) {  # DiaNN specific:
  if (!"Quantity Quality" %in% colnames(pep)) {
    tmp <- data.table(Qual = ev$"Quantity Quality", ModSeq = ev$"Modified sequence")
    tmp <- tmp[, .(Qual= mean(Qual)), by = .(ModSeq)]
    tmp <- as.data.frame(tmp)
    pep$"Quantity Quality" <- tmp$Qual[match(pep$"Modified sequence", tmp$ModSeq)]
  }
  weightsInsrt <- "-log10(PEP * CV * (1 - DiaNN_Quantity_Quality))" 
  pep$Weights <- -log10(pep$PEP*pep$CV*(1-pep$"Quantity Quality"))
  # plot <- ggplot(pep) +
  #   geom_density(stat = "density", aes(x = PEP), colour = "red") +
  #   geom_density(stat = "density", aes(x = CV), colour = "green") +
  #   geom_density(stat = "density", aes(x = 1-`Quantity Quality`), colour = "blue") +
  #   theme_bw()
  # poplot(plot)
  # plot <- ggplot(pep) +
  #   geom_density(stat = "density", aes(x = -log10(PEP)), colour = "red") +
  #   geom_density(stat = "density", aes(x = -log10(CV)), colour = "green") +
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
if ("Prot.Quant.Mod.Excl.is.strict" %in% colnames(Param)) {
  tmp <- suppressWarnings(as.logical(Param$Prot.Quant.Mod.Excl.is.strict))
  if (validLogicPar("tmp")) { Mod.Excl.is.strict <- tmp }
}
Discard.unmod %<o% Mod.Excl.is.strict+1
if (Discard.unmod == 1) { Discard.unmod <- as.logical(Discard.unmod) }

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
                              expMap = Exp.map,
                              expMap_Samples_col = "Ref.Sample.Aggregate",
                              param = Param,
                              pepInt_Root = pep.ref[length(pep.ref)],
                              pepRat_root = pep.ratios.ref[length(pep.ratios.ref)],
                              pepInt_log = FALSE,
                              pepRat_log = 2,
                              protLFQ_toLog = TRUE,
                              protRat_toLog = TRUE,
                              mods_to_Exclude = Mod2Xclud,
                              mod_Seq = "Modified sequence",
                              discard_unmod = Discard.unmod,
                              prim_Seq = "Sequence (1st accession)",
                              ref_Mode = RefRat_Mode,
                              cl = parClust,
                              N.clust = N.clust,
                              N.reserved = 1)
insrt <- paste0(", and quantified using ",
                c("limpa's dpcQuant function",
                  paste0("an in-house MaxLFQ-like algorithm which computes a protein group-level profile across samples as the weighted mean of individual, normalized peptidoform profiles (weights = ",
                         weightsInsrt, ")"),
                  "MaxLFQ, as implemented in the iq package",
                  "the aggregateFeatures function form the QFeatures package")[match(quantAlgo, quantAlgoOpt)],
                ".")
if (reScAlgo != quantAlgo) {
  insrt <- paste0(insrt, " Row-normalized relative protein profiles were re-scaled ")
  if (reScAlgo %in% c("limpa", "QFeatures", "MaxLFQ (iq)")) {
    insrt <- paste0(insrt,
                    "to the average row-wise expression scales provided by the output of ",
                    c("limpa's dpcQuant",
                      "QFeatures' aggregateFeatures",
                      "iq's fast_MaxLFQ")[match(reScAlgo, c("limpa", "QFeatures", "MaxLFQ (iq)"))], " function.")
  }
  if (reScAlgo == "topN") {
    if (topN_correct) {
      insrt <- paste0(insrt,
                      "using the top", topN, " method (i.e. to the mean of up to ", topN, " highest intensity peptides).")
    } else {
      insrt <- paste0(insrt,
                      "using a variation on the top", topN, " method, first correcting for each protein peptide intensities ranked by decreasing values for systematic rank-wise intensity bias, then averaging up to ", topN, " peptides.")
    }
  }
  if (reScAlgo == "max") {
    insrt <- paste0(insrt, "to the value of the highest intensity peptides.")
  }
  if (reScAlgo == "weighted.mean") {
    if (quantAlgo == "LM") {
      insrt <- paste0(insrt, "using the weighted mean of each protein group's peptides.")
    } else {
      insrt <- paste0(insrt, "using the mean of each protein group's peptides weighted by ", weightsInsrt, ".")
    }
  }
  if (reScAlgo %in% c("median", "sum")) {
    insrt <- paste0(insrt, "using the ", reScAlgo, " of each protein group's peptides.")
  }
  if (!reScAlgo %in% reScAlgoOpt) {
    insrt <- paste0(insrt, "to the result of summarizing each protein group's peptides with the ", reScAlgo, " function.")
  }
}
DatAnalysisTxt <- gsub("\\.$", insrt, DatAnalysisTxt)
DatAnalysisTxt <- paste0(DatAnalysisTxt, " Estimated expression values were log10-converted...")

saveFun(quantData_list, file = "PG_quant.RData")
quantData <- quantData_list$Data
Prot.Expr.Root %<o% c(Original = "log10(Expr.) - ")
Prot.Rat.Root %<o% c(Original = "log2(Ratio) - ")
stopifnot(length(grep(topattern(Prot.Expr.Root), colnames(quantData))) > 0,
          length(grep(topattern(Prot.Rat.Root), colnames(quantData))) > 0)
