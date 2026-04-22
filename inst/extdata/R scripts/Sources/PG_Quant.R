# Protein group quantitation
# --------------------------
#
# This script calculates protein group-level quantitative values based on those of individual peptides.
# Exact quantitation parameters (e.g. summarization method) are defined earlier in the workflow.

# Default parameter for reps, to be over-written when re-running post-re-normalisation
if (!validLogicPar("post_ReNorm_reRun")) {
  post_ReNorm_reRun <- FALSE
}

# You may want to exclude peptides based on how many samples they are found in:
if (scrptType == "noReps") {
  kol <- paste0(int.col, " - ", Exp)
  tmpPar <- AnalysisParam
}
if (scrptType == "withReps") {
  kol <- lapply(VPAL$values, \(x) {
    x <- paste0(pep.ref["Original"], Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)])
    return(x[which(x %in% colnames(pep))])
  })
  tmpPar <- Param
}
#
# Min number of samples with valid value
if (!exists("PepFoundInAtLeast")) { PepFoundInAtLeast <- 1L }
if ("PepFoundInAtLeast" %in% names(tmpPar)) {
  tmp <- suppressWarnings(as.integer(tmpPar$PepFoundInAtLeast))
  if (validIntegPar("tmp")) {
    PepFoundInAtLeast <- tmp
  }
}
PepFoundInAtLeast %<o% PepFoundInAtLeast
clusterExport(parClust, "is.all.good", envir = environment())
tst <- parApply(parClust, pep[, unlist(kol), drop = FALSE], 1L, \(x) {
  sum(is.all.good(x) > 0)
})
Pep2Use %<o% which(tst >= PepFoundInAtLeast)
#
# Min number of samples with a valid value PER GROUP
if (scrptType == "withReps") {
  maxAllowed <- max(c(2L, length(Rep)-1L))
  if (!exists("PepFoundInAtLeastGrp")) { PepFoundInAtLeastGrp <- 1L }
  if ("PepFoundInAtLeastGrp" %in% names(tmpPar)) {
    tmp <- suppressWarnings(as.integer(tmpPar$PepFoundInAtLeastGrp))
    if (validIntegPar("tmp")) {
      PepFoundInAtLeastGrp <- tmp
    }
  }
  if (PepFoundInAtLeastGrp > maxAllowed) {
    warning(paste0(" Invalid PepFoundInAtLeastGrp parameter, setting to ", maxAllowed))
    PepFoundInAtLeastGrp <- maxAllowed
  }
  tst <- pep[, unlist(kol)]
  clusterExport(parClust, list("tst", "kol", "PepFoundInAtLeastGrp", "is.all.good"), envir = environment())
  tst <- parSapply(parClust, 1L:nrow(pep), \(x) {
    x <- max(vapply(kol, \(kl) {
      sum(is.all.good(unlist(tst[x, kl])) > 0)
    }, 1L) >= PepFoundInAtLeastGrp)
    return(x)
  }) > 0L
  invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
  w <- which(!tst)
  Pep2Use[w] <- FALSE
}
#length(Pep2Use)/nrow(pep)

# Map
bckpNm <-"PG_quant"
if (scrptType == "noReps") {
  smplsMap <- SamplesMap
  smplsKol <- "Experiment"
  PG.int.col %<o% "log10(Expr.) - "
  PG.int.cols %<o% setNames(PG.int.col, "Original")
  PG.rat.col %<o% "log2(Ratio) - " #We need those defaults actually even if !MakeRatios
  PG.rat.cols %<o% setNames(PG.rat.col, "Original")
  pepInt_col <- paste0(int.cols[max(which(names(int.cols) != "Imputed"))], " - ")
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
    pepRat_col <- paste0(rat.cols[max(which(names(rat.cols) != "Imputed"))], " - ")
  } else {
    pepRat_col <- Aggregate.map <- Aggregate.list <- Aggregates <- NA
  }
}
if (scrptType == "withReps") {
  smplsMap <- ExpMap
  smplsKol <- "Ref.Sample.Aggregate"
  if (post_ReNorm_reRun) { bckpNm <- paste0(bckpNm, "_reNorm") }
  refNm <- if (post_ReNorm_reRun) { "Back-norm" } else {
    rev(setdiff(names(pep.ref), "Back-norm"))[1L]
  }
  pepInt_col <- pep.ref[refNm]
  cat(paste0("Input peptide intensities = ", refNm, "\n"))
  #pepRat_col <- pep.ratios.ref[length(pep.ratios.ref)]
  Prot.Expr.Root %<o% c(Original = "log10(Expr.) - ")
  Prot.Rat.Root %<o% c(Original = "log2(Ratio) - ")
}

# Weights:
# - Higher for peptides with low PEP.
# - Reps: Higher for peptides with low intra-sample group CV on average
# - No reps: CV not taken into account since we don't know what a good (regulatory) or bad CV are,
#   and we may even only have 1 peptide.
if (scrptType == "withReps") {
  source(parSrc, local = FALSE)
  Kols <- paste0(pepInt_col, smplsMap$Ref.Sample.Aggregate)
  Kols <- Kols[which(Kols %in% colnames(pep))]
  tmp <- pep[, Kols]
  clusterExport(parClust, list("smplsMap", "VPAL", "pep.ref", "is.all.good", "tmp", "pepInt_col"), envir = environment())
  CV <- parSapply(parClust, VPAL$values, \(x) { #x <- VPAL$values[1L]
    smpls <- smplsMap$Ref.Sample.Aggregate[which(smplsMap[[VPAL$column]] == x)]
    kols <- paste0(pepInt_col, smpls)
    kols <- intersect(kols, colnames(tmp))
    x <- apply(tmp[, kols, drop = FALSE], 1L, \(y) {
      y <- is.all.good(log10(unlist(y)))
      y <- if (length(y)) {
        sd(y)/mean(y)
      } else { NA }
      return(y)
    })
    return(x)
  })
  pep$CV <- rowMeans(CV, na.rm = TRUE)
  pep$Weights <- -log10(pep$PEP*pep$CV)
  weightsInsrt <- "-log(PEP * CV)"
}
if (scrptType == "noReps") {
  pep$Weights <- -log10(pep$PEP)
  weightsInsrt <- "-log(PEP)"
}
if ((length(inDirs) == 1L)&&(QuantUMS)&&("Quantity Quality" %in% colnames(ev))) {  # DiaNN specific:
  if (!"Quantity Quality" %in% colnames(pep)) {
    tmp <- data.table(Qual = ev$"Quantity Quality", ModSeq = ev$"Modified sequence")
    tmp <- tmp[, .(Qual= mean(Qual)), by = .(ModSeq)]
    tmp <- as.data.frame(tmp)
    pep$"Quantity Quality" <- tmp$Qual[match(pep$"Modified sequence", tmp$ModSeq)]
  }
  if (scrptType == "noReps") {
    weightsInsrt <- "-log10(PEP * (1 - DiaNN_Quantity_Quality))" 
    pep$Weights <- -log10(pep$PEP*(1-pep$"Quantity Quality"))
  }
  if (scrptType == "withReps") {
    weightsInsrt <- "-log10(PEP * CV * (1 - DiaNN_Quantity_Quality))"
    pep$Weights <- -log10(pep$PEP*pep$CV*(1-pep$"Quantity Quality"))
  }
  # plot <- ggplot(pep) +
  #   geom_density(stat = "density", aes(x = PEP), colour = "red")
  # if (scrptType == "withReps") {
  #   plot <- plot +
  #     geom_density(stat = "density", aes(x = CV), colour = "green")
  # }
  # plot <- plot +
  #   geom_density(stat = "density", aes(x = 1-`Quantity Quality`), colour = "blue") +
  #   theme_bw()
  # poplot(plot)
  # plot <- ggplot(pep) +
  #   geom_density(stat = "density", aes(x = -log10(PEP)), colour = "red")
  # if (scrptType == "withReps") {
  #   plot <- plot +
  #     geom_density(stat = "density", aes(x = -log10(CV)), colour = "green")
  # }
  # plot <- plot +
  #   geom_density(stat = "density", aes(x = -log10(1-`Quantity Quality`)), colour = "blue") +
  #   theme_bw()
  # poplot(plot, 12L, 22L)
}
w <- which(!is.all.good(pep$Weights, 2L))
pep$Weights[w] <- min(pep$Weights, na.rm = TRUE)
#
m <- max(pep$Weights)
pep$Weights <- pep$Weights/m
pep$Weights[which(pep$Weights < 0.001)] <- 0.001
summary(pep$Weights)
#
Mod.Excl.is.strict %<o% FALSE
if ("Prot.Quant.Mod.Excl.is.strict" %in% names(tmpPar)) {
  tmp <- suppressWarnings(as.logical(tmpPar$Prot.Quant.Mod.Excl.is.strict))
  if (validLogicPar("tmp")) { Mod.Excl.is.strict <- tmp }
}
Discard.unmod %<o% (Mod.Excl.is.strict+1L)
if (Discard.unmod == 1L) { Discard.unmod <- as.logical(Discard.unmod) }

# Quantitation
source(parSrc, local = FALSE)
reScalingAlgo <- reScAlgo
if (reScAlgo == "topN") { reScalingAlgo <- paste0("top", topN) }
quantArgs <- list(Prot = PG,
                  Pep = pep[Pep2Use,],
                  pg_PepIDs = Pep4Quant,
                  pg_PepIDs_unique = "Unique peptide IDs",
                  pep_IDs = "id",
                  N_unique = N_unique_Pep,
                  LFQ_algo = quantAlgo,
                  reScaling = reScalingAlgo,
                  topN_correct = topN_correct,
                  minN = N_Pep,
                  maxN = Inf,
                  Weights = "Weights",
                  useIntWeights = FALSE,
                  Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1L],
                  skipRatios = !MakeRatios,#scrptType != "withReps")&&(!MakeRatios), # For reps, t
                  expMap = smplsMap,
                  expMap_Samples_col = smplsKol,
                  pepInt_Root = pepInt_col,
                  pepRat_root = NULL,
                  pepInt_log = FALSE,
                  pepRat_log = 2L,
                  mods_to_Exclude = Mod2Xclud,
                  mod_Seq = "Modified sequence",
                  discard_unmod = Discard.unmod,
                  prim_Seq = "Sequence (1st accession)",
                  cl = parClust,
                  N.clust = N.clust,
                  N.reserved = 1L,
                  refGroups = NULL,
                  ratGroups = NULL)
if (scrptType == "withReps") { quantArgs$param <- tmpPar }
if ((scrptType == "noReps")&&(MakeRatios)) {
  quantArgs$refGroup <- list(values = unique(smplsMap$Ratios_group),
                                 names = "Ratios_group",
                                 column = "Ratios_group")
  quantArgs$ratGroups <- list(values = unique(smplsMap$Ratios_group),
                              names = "Ratios_group",
                              column = "Ratios_group")
}
# For testing:
#DefArg(protQuant)
#invisible(lapply(names(quantArgs), \(x) { assign(x, quantArgs[[x]], envir = .GlobalEnv); return() }))
#TESTING <- TRUE
quantData_list %<o% do.call(protQuant, quantArgs) # In case this is run after PG reNorm (from back-reNormalized peptides),
quantData_list$reNormalized <- post_ReNorm_reRun
# this replaces the original quantData_list object.
# This is fine, it ensures for limpa we use the latest version.
if (scrptType == "withReps") {
  g <- grep("^log2FC - ", colnames(quantData_list$Data))
  colnames(quantData_list$Data)[g] <- sub("^log2FC - ", Prot.Rat.Root, colnames(quantData_list$Data)[g])
}
myQuantBckpFl <- paste0(wd, "/", bckpNm, ".RData")
saveFun(quantData_list, file = myQuantBckpFl)
#loadFun(myQuantBckpFl)
#
opt <- setNames(c("limpa's dpcQuant function",
                  paste0("an in-house MaxLFQ-like algorithm which computes a protein group-level profile across samples as the weighted mean of individual, normalized peptidoform profiles (weights = ",
                         weightsInsrt, ")"),
                  "MaxLFQ, as implemented in the iq package",
                  "the aggregateFeatures function form the QFeatures package"),
                c("limpa", "LM", "MaxLFQ (iq)", "QFeatures"))
insrt <- paste0(", and quantified using ", opt[quantAlgo], ".")
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
    insrt <- if (topN_correct) {
      paste0(insrt,
             "using the top", topN, " method (i.e. to the mean of up to ", topN, " highest intensity peptides).")
    } else {
      paste0(insrt,
             "using a variation on the top", topN, " method, first correcting for each protein peptide intensities ranked by decreasing values for systematic rank-wise intensity bias, then averaging up to ", topN, " peptides.")
    }
  }
  if (reScAlgo == "max") {
    insrt <- paste0(insrt, "to the value of the highest intensity peptides.")
  }
  if (reScAlgo == "weighted.mean") {
    insrt <- if (quantAlgo == "LM") {
      paste0(insrt, "using the weighted mean of each protein group's peptides.")
    } else {
      paste0(insrt, "using the mean of each protein group's peptides weighted by ", weightsInsrt, ".")
    }
  }
  if (reScAlgo %in% c("median", "sum")) {
    insrt <- paste0(insrt, "using the ", reScAlgo, " of each protein group's peptides.")
  }
  if (!reScAlgo %in% reScAlgoOpt) {
    insrt <- paste0(insrt, "to the result of summarizing each protein group's peptides with the ", reScAlgo, " function.")
  }
}
l <- length(DatAnalysisTxt)
DatAnalysisTxt[l] <- gsub("\\.$", insrt, DatAnalysisTxt[l])
DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l], " Estimated expression values were log10-converted...")

quantData <- quantData_list$Data
if (scrptType == "withReps") {
  stopifnot(length(grep(topattern(Prot.Expr.Root), colnames(quantData))) > 0L)
            #length(grep(topattern(Prot.Rat.Root), colnames(quantData))) > 0L)
  # For script with reps, the data is added later to PG
  g <- grep(topattern(Prot.Rat.Root), colnames(quantData))
  if (!length(g)) {
    # We now calculate ratios separately
    # Add code here if you need to re-calculate them earlier
  }
}
if (scrptType == "noReps") {
  quantData <- quantData[, which(!grepl("\\.REF$", colnames(quantData))), drop = FALSE]
  colnames(quantData) <- gsub("^log10\\(Expr\\.\\) - ", PG.int.col, colnames(quantData))
  if (MakeRatios) {
    colnames(quantData) <- gsub("^log2\\(Ratio\\) - ", PG.rat.col, colnames(quantData))
  }
  m1 <- match("Peptides IDs used for quantitation", colnames(quantData))
  colnames(quantData)[m1] <- paste0("Peptide IDs used for quantitation - ",
                                     names(PG.int.cols)[match(PG.int.col, PG.int.cols)])
  PG[, colnames(quantData)] <- quantData
}
if ((scrptType == "noReps")&&(Impute)) {
  # In that case we need to calculate expression values a second time, unfortunately... It doesn't take so long.
  PG.int.cols["Imputed"] <- PG.int.col <- paste0("Imput. ", PG.int.cols["Original"])
  PG.rat.cols["Imputed"] <- PG.rat.col <- paste0("Imput. ", PG.rat.cols["Original"])
  pepImpInt_col <- paste0(int.cols["Imputed"], " - ")
  pepImpRat_col <- paste0(rat.cols["Imputed"], " - ")
  #
  quantArgs_Imp <- quantArgs
  quantArgs_Imp$pepInt_Root <- pepImpInt_col
  quantArgs_Imp$pepRat_root <- pepImpRat_col
  quantData_list_Imput %<o% do.call(protQuant, quantArgs_Imp)
  myQuantBckpFl2 <- paste0(wd, "/", bckpNm, "2.RData")
  saveFun(quantData_list_Imput, file = myQuantBckpFl2)
  #loadFun(myQuantBckpFl2)
  quantData2 <- quantData_list_Imput$Data
  quantData2 <- quantData2[, which(!grepl("\\.REF$", colnames(quantData2))), drop = FALSE]
  colnames(quantData2) <- gsub("^log10\\(Expr\\.\\) - ", PG.int.col, colnames(quantData2))
  if (MakeRatios) {
    colnames(quantData2) <- gsub("^log2\\(Ratio\\) - ", PG.rat.col, colnames(quantData2))
  }
  m2 <- match("Peptides IDs used for quantitation", colnames(quantData2))
  colnames(quantData2)[m2] <- "Peptide IDs used for quantitation - Imputed"
  PG[, colnames(quantData2)] <- quantData2
}
