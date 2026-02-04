# Currently 9 methods are implemented, but the last 3 are just for testing and do not create all columns required for the script to complete
# (These could be added relatively easily if necessary though)
QuantData <- setNames(paste0("quant.data", seq_along(QuantMethods)), QuantMethods)
QuantMethods_all <- FALSE
.obj <- unique(c("QuantData", "QuantMethods_all", .obj))
exprsCol <- paste0("log10(Expr.) - ", RSA$values)
# Weights:
# - Higher for peptides with low intra-sample group CV on average
# - Higher for peptides with low PEP
if ((length(inDirs) == 1)&&(QuantUMS)) {
  weightsInsrt <- "mean DiaNN \"Quantity Quality\"" 
  pep$Weights <- pep$"Quantity Quality"
} else {
  weightsInsrt <- "-log(PEP)/CV" 
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
  CV <- rowMeans(CV, na.rm = TRUE)
  pep$Weights <- -log10(pep$PEP)/CV
}
summary(pep$Weights)
m <- max(is.all.good(pep$Weights))
pep$Weights <- pep$Weights/m
pep$Weights[which((is.na(pep$Weights))|(pep$Weights < 0.001))] <- 0.001
if ("Prot.Quant.Mod.Excl.is.strict" %in% colnames(Param)) {
  Mod.Excl.is.strict <- Param$Prot.Quant.Mod.Excl.is.strict
  if (!Mod.Excl.is.strict %in% c(1, 0, TRUE, FALSE)) { Mod.Excl.is.strict <- FALSE }
} else { Mod.Excl.is.strict <- FALSE }
Discard.unmod <- Mod.Excl.is.strict+1
if (Discard.unmod == 1) { Discard.unmod <- as.logical(Discard.unmod) }
.obj <- unique(c("Mod.Excl.is.strict", "Discard.unmod", .obj))
if (!grepl("^Prot\\.Quant", Param$QuantMeth)) {
  stop("NB: currently only methods 1 to 6 provide all necessary columns, not just quantitative columns per sample but also reference columns, ratios, etc... Until those are added they cannot be used by this script. Defaulting to method 3!")
  Param$QuantMeth <- "Prot.Quant.Unique"
}

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
