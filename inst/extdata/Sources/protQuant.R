#' protQuant
#' Source version of protQuant because as with many other cases I just can't resolve the slow parallelisation issues.
#'
#' @description
#' A function to calculate estimated protein group Expression (inter-protein group quantitation vectors) from individual peptide values.
#' For legacy reasons, the function can also output Ratios.
#' The input is assumed to already be normalized!\cr
#' 
#' 
#' @param Prot Protein/Protein groups table. A data.frame, which must contain at least peptides IDs and primary sequence of the first accession in each group.
#' @param Pep Peptides table. A data.frame. If it contains unmodified sequences then these are expected to be in column name "Sequence". Argument mod_Seq controls the name of the column used for modified sequence.
#' @param pg_PepIDs Name of the Protein/Protein groups table's peptide IDs column. Default = "Peptide IDs"
#' @param pg_PepIDs_unique Used only if N_unique is larger than 0. Name of the Protein/Protein groups table's unique peptide IDs column. Default = "Unique peptide IDs"
#' @param pep_IDs The name of the Peptides table's IDs column. Default = "id"
#' @param N_unique Logical or numeric. If FALSE or 0, the function just uses the "pg_PepIDs" argument. If non null (default), a second column name of unique peptide IDs should be provided using the "pg_PepIDs_unique" argument. If at least N_unique unique peptides are present, then only those will be used for quantitation, otherwise as many razor/shared peptides as necessary will be added (sorted by decreasing average intensities).
#' @param LFQ_algo Algorithm used to compute average profiles. One of:\cr
#'  - "LM": Levenberg-Marquardt method, backend = [minpack.lm::nls.lm()]; very similar to MaxLFQ. Normalized peptide relative profiles are aligned using Levenberg-Marquardt then summarized.\cr
#'  - "iq": iq's fast implementation of MaxLFQ, backend = [iq::fast_MaxLFQ()]\cr
#'  - "limpa" or "DPC": backend = [limpa::dpcQuant()]; current default\cr
#'  - "QFeatures": backend = [QFeatures::aggregateFeatures()]\cr
#'  - "MSstats": backend = [MSstats::dataProcess()]\cr
#' Note that limpa, QFeatures and MSstats do not "just produce quantitative values": they output specific objects or lists with complex modelling of the dataset, which can be used for downstream analysis (limpa -> limma, DFeatures -> msqrob2, MSstats covers modeling and stats).\cr
#' @param reScaling Optional summary method for re-scaling. May be one of:\cr
#'  - "median",\cr
#'  - "mean",\cr
#'  - "weighted.mean" (requires the "Weights" argument),\cr
#'  - "max",\cr
#'  - "sum" (not recommended),\cr
#'  - "MaxLFQ": like sum, but the value used is the value before any peptides are filtered out (not recommended),\cr
#'  - "topN", where N should be the maximum number of peptides to average (e.g. "top3" - do not use "topN" as there is no default value for N!)\cr
#'  - the name of any valid value of LFQ_algo, which allows re-scaling any LFQ algorithm using the scaling provided by another. Note that this doesn't make sense for all, as some are only meant to model logFC, not estimated abundance scale.\cr
#'  - (you could also use the name of any other available averaging function, this should work in principle assuming similar syntax)
#' @param topN_correct Logical, default = TRUE. In the case where we are using more than one peptide for the re-scaling step, should we correct for systematic peptide intensity biases between peptides of rank 1, 2, 3 and so on and so forth?
#' @param minN Integer, default = 1. How many peptides should at least be present for quantitation? Values lower than 1 are increased to 1. May not be higher than N_unique!
#' @param maxN Integer or Inf, default = 50. Up to how many peptides should we use for the Levenberg-Marquardt procedure (used only for LFQs)? Using too many peptides can be an issue, e.g. with huge proteins like Titin. Default = 50. The most intense peptides will be selected. May not be lower than N_unique!
#' @param LM_fun How should normalized profiles be averaged? One of "median" (default), "mean", or "weighted.mean" (the latter which uses the Weights and useIntWeights arguments).
#' @param Weights Length 1 character, a valid column name of Pep containing user-defined individual peptide weights. Used if LM_fun or reScaling are set to "weighted.mean" (for the former only if LFQ_algo = "LM").
#' @param useIntWeights Logical, default = FALSE. Ignored unless LM_fun or reScaling are set to "mean" or "weighted.mean". If TRUE, will take into account individual peptide intensities when calculating average profile for that step (thus, it will actually be a weighted mean regardless). DOES NOT replace the optional, user-provided Weights, but instead multiplies the former by new intensity-based factors.
#' @param Priority One of "Intensities" (default) or "Ratios" (some flexibility in spelling is allowed). In some rare cases, such as a SILAC dataset processed with MaxQuant, we will want to correct intensities - prior to running the main algorithm - so their ratios reflect the more accurate ratios directly measured by the search engine.
#' @param skipRatios Default = FALSE. If TRUE, ratios will not be calculated.
#' @param experimentMap Map of the experiment map.
#' @param experimentMap_Samples_col Names of the single samples column in the experiment's map. Default = "Sample"
#' @param param This analysis' parameters. If provided, the refGroups argument is not required.
#' @param aggrMap Only needed if ratios (logFCs) are also to be output, lternative way of defining groups when Param is missing. The analysis' aggregate map. Default = Aggregate.map
#' @param aggrList Only needed if ratios (logFCs) are also to be output, lternative way of defining groups when Param is missing. Named list of this analysis' factors aggrNames. Default = Aggregate.list
#' @param aggrNames Only needed if ratios (logFCs) are also to be output, lternative way of defining groups when Param is missing. The experiment's factor aggrNames. Default = Aggregates
#' @param refGroups Only needed if ratios (logFCs) are also to be output. List defining which samples are paired to which references. May alternatively (preferred solution) be provided indirectly through the param argument.
#' @param ratGroups Only needed if ratios (logFCs) are also to be output. List defining groups within which ratios are calculated. May alternatively (preferred solution) be provided indirectly through the param argument.
#' @param pepInt_Root Root of the peptides intensity column(s) names
#' @param pepRat_root Only needed if ratios (logFCs) are also to be output. Root of the peptides ratios column(s) names, used if Priority = "Ratios".
#' @param pepInt_log Set to 0 or FALSE if input peptide intensities are linear scale. If the data is already log scale, set to the relevant scale's base. Default = FALSE
#' @param pepRat_log Only needed if ratios (logFCs) are also to be output. Set to 0 or FALSE if input peptide ratios are linear scale. If the data is already log scale, set to the relevant scale's base. Default = 2
#' @param protLFQ_toLog Should the output protein LFQ values be log-scale or not? Can ve set to the log base desired. If set to TRUE, this will be the same base as the input intensities log scale, or 10 by default.
#' @param protRat_toLog Only needed if ratios (logFCs) are also to be output. Should the output protein ratios be log-scale or not? Can ve set to the log base desired. If set to TRUE, this will be the same base as the input ratios log scale, or its default, 2
#  @param Mods Which modifications (2 lowercase letters PTM code) should be included? If set to FALSE, will not filter any modifications.
#' @param mods_to_Exclude Which modifications should be excluded? (use argument "discard_unmod" to discard unmodified counterpart peptides.) A data.frame with columns "Mark" (2-lettern modification mark) and "Where" (a list, which amino acids are affected, use "Nterm", "Cterm", "protNterm" and "protCterm" for termini). Also see argument "discard_unmod".
#' @param mod_Seq Default = "Modified sequence". The name of the column containing the modified sequence in the peptides table. Can be set to a non-modified sequence if mods_to_Exclude is empty.
#' @param discard_unmod Logical or integer in 0:2. Default = TRUE. Should we discard those unmodified peptides whose primary sequence is the same as that of some modified peptides we will not use? If set to 2, will use the "Where" column in "mods_to_Exclude" to identify (and exclude) peptides which could be modified even if the modified form was not identified. Requires knowledge of protein sequence ("prim_Seq" argument)! Be careful! This will likely massively reduce the number of peptides available for quantitation!
#' @param prim_Seq Default = "Sequence (1st accession)", used if "discard_unmod" is set to 2.
#' @param ref_Mode How are reference ratios calculated?\cr
#'  - If set to "1", only references are considered (i.e. it compares individual references either to each other, or if available to the average reference for the group).\cr
#'  - If set to "2" (default), for each ratios group, reference ratios are based on comparing every possible pair of samples within the group.\cr
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param contrasts Optional contrasts data.frame for which to calculate ratios. It should contain columns of the form "[A-D]_samples" mapping samples to contrasts.
#' @param alsoRun Any valid option for LFQ_algo which was not run as LFQ_algo or reScaling can also be run here, e.g. for comparison between algorithms. The output object(s) will be added to the output list, but will not be returned as main quantitation in the Data slot.
#' 
#' #' @details
#' This function is meant to work from normalised peptide intensities and:\cr
#'  - calculates a quantitative profile for each protein,\cr
#'  - maximizes the usage of peptides-level ratios information.\cr
#' Specifically, it can be broken down into 3 sub-steps:\cr
#'   a) filtering peptides eligible for quantitation\cr
#'   b) estimating/modeling (depending on algorithm) a relative protein profile across samples:\cr
#'    - LM: this uses Levenberg-Marquardt to align peptide profiles, then summarizes them, and is broadly similar to MaxLFQ.\cr
#'    - iq: the iq package's implementation of MaxLFQ, specifically the [iq::fast_MaxLFQ()] function.\cr
#'    - limpa: uses [limpa::dpcQuant()].\cr
#'    - QFeatures: uses [QFeatures::aggregateFeatures()].\cr
#'    - MSstats: uses [MSstats::dataProcess()].\cr
#'   c) Optional re-scaling, i.e. "anchoring" the resulting relative profile to a relative abundance value, such that different proteins are ranked relatively in a manner which - whilst imperfect because of widely different individual peptide detectabilities - would still be correct if the latter could be corrected for, or would otherwise operate under assumptions at least minimizing the expected error. See the reScaling argument.\cr
#'    You may also rescale any quantitation (as specified in LFQ_algo) with the scale provided by another (using the reScaling argument), thus allowing for using the relative quantitation and scaling of any two methods.\cr
#'    To re-scale the same way as MaxLFQ, use "sum" - but see N.B. below!\cr
#'    iBAQ - which would be the equivalent of using "sum" then dividing by the number of observable peptides, but the additional step could easily be added outside this function.\cr
#'\cr
#' NB on re-scaling:\cr
#'   The early - and still occasional - practice of summing peptide intensities to estimate protein abundance (as used in MaxLFQ) is conceptually unsound.\cr
#'   The intensity of a unique peptide with perfect detectability is a function of the peptide's copy number, itself originally (before losses) equal to that of the parent protein.\cr
#'   If detectability could be corrected for, the intensity values of all unique peptides unaffected by PTMs would be a direct unbiased estimator of parent protein copy number.\cr
#'   For this reason, any intensity-based protein abundance estimate much be 'homogenous' with single peptide intensities, either through peptide selection or averaging.\cr
#'   We prefer the topN (1 to 3) methods, which leverage the usually better data S/N ratio of higher intensity peptides, and use the values closer to the original target.\cr
#'   iBAQ's proposed correction of summed intensities by dividing by the number of observable peptides, seems to make little sense: if averaging, one should divide by the number of observed, not observable, values. Hence why this is not implemented.\cr
#'   
#'\cr
#'\cr
#' (Note on refGroups and ratGroups:\cr
#'  These parameters are the method we used for calculating logFCs and together define whether replicates are paired or not.\cr
#'  For an experiment with several Replicates of 2 or more Conditions (incl. one control, aka reference):\cr
#'  - In a paired (= "nested") setup, you would set refGroups to "ExpRep" and ratGroups to "Exp"\cr
#'  - In an unpaired setup, you would set both refGroups and ratGroups to "Exp"\cr
#'  Essentially, refGroups are the groups within which ratios are calculated to all available references in the group,\cr
#'  while ratGroups are comparison groups, i.e. several sample groups to compare to one another.)\cr
#'\cr
#' 
#' @returns
#' A list:
#'  - "Data": output quantitative data.
#'  - "EList_obj": EList object, only present if LFQ_algo or reScaling = "limpa"
#'  - "QFeatures_obj": QFeatures object, only present if LFQ_algo or reScaling  = "QFeatures"
#'  - "MSstats": list created by MSstats::dataProcess() 
#' 
#' 
#' @export


TESTING <- FALSE


DefArg(protQuant)#;TESTING <- TRUE
invisible(lapply(names(quantArgs2), \(x) { assign(x, quantArgs2[[x]], envir = .GlobalEnv); return() }))


if (!exists("ref_Mode")) { ref_Mode <- 2L }
if (!as.numeric(ref_Mode) %in% 1L:2L) { ref_Mode <- 2L }
ref_Mode <- as.character(ref_Mode)
#
# Check arguments
mySmpls <- experimentMap[[experimentMap_Samples_col]]
nPep_0 <- nrow(Pep)
minN <- max(c(1L, suppressWarnings(abs(as.integer(minN)))))
maxN <- max(c(1L, minN, suppressWarnings(abs(as.numeric(maxN)))))
N_unique <- max(c(0L, suppressWarnings(abs(as.numeric(N_unique)))))
if (is.finite(maxN)) { maxN <- suppressWarnings(as.integer(maxN)) }
if (is.finite(N_unique)) { N_unique <- suppressWarnings(as.integer(N_unique)) }
stopifnot(length(N_unique) == 1L,
          length(maxN) == 1L,
          length(minN) == 1L,
          sum(!is.na(c(N_unique, minN, maxN))) == 3L)
cat("Using at least", minN, "and up to", maxN, "peptides...\n")
if (N_unique) {
  stopifnot(pg_PepIDs_unique %in% colnames(Prot))
  if (pg_PepIDs == pg_PepIDs_unique) {
    N_unique <- 0L
  }
}
if (N_unique) {
  cat("If available, up to", N_unique, "unique (= proteotypic) peptides will be used...\n")
}
#
# Algorithms: deal with synonyms and names degeneracy
algoSyn <- data.frame(ALGO = c("LM", "IQ", "LIMPA", "QFEATURES", "QFEATURES", "MSSTATS"),
                      SYNONYM = c(NA_character_, "MAXLFQ", "DPCQUANT", "MSQROB", "MSQROB2", NA_character_))
#
LFQ_algo <- gsub(" -_\\.", "", LFQ_algo)
if (!exists("alsoRun")) { alsoRun <- c() }
alsoRun <- gsub(" -_\\.", "", alsoRun)
if (!exists("reScaling")) { reScaling <- "skip" }
reScaling <- gsub(" ", "", reScaling) # NB: including "-_\\." would be too dangerous, because reScaling also takes function names!
skip_reScaling <- reScaling == "skip"
LFQ_ALGO <- toupper(LFQ_algo)
if (LFQ_ALGO %in% algoSyn$SYNONYM) { LFQ_ALGO <- algoSyn$ALGO[match(LFQ_ALGO, algoSyn$SYNONYM)] }
ALSORUN <- toupper(alsoRun)
w <- which(ALSORUN %in% algoSyn$SYNONYM)
if (length(w)) { ALSORUN[w] <- algoSyn$ALGO[match(ALSORUN[w], algoSyn$SYNONYM)] }
if (!skip_reScaling) {
  RESCALING <- toupper(reScaling)
  if (RESCALING %in% algoSyn$SYNONYM) { RESCALING <- algoSyn$ALGO[match(RESCALING, algoSyn$SYNONYM)] }
}
if (length(mySmpls) == 1L) {
  if (LFQ_ALGO %in% c("LIMPA", "QFEATURES")) {
    warning("Only one sample, LFQ_algo can only be one of LM or iq, defaulting to LM...")
    LFQ_algo <- LFQ_ALGO <- "LM"
  }
  if (RESCALING %in% c("LIMPA", "QFEATURES")) {
    reScaling <- RESCALING <- "LM"
  }
  if (length(ALSORUN)) {
    w <- which(!ALSORUN %in% c("LIMPA", "QFEATURES"))
    alsoRun <- alsoRun[w]
    ALSORUN <- ALSORUN[w]
  }
}
if (scrptTypeFull != "withReps_PG_and_PTMs") {
  if (LFQ_ALGO == "MSSTATS") {
    warning("MSstats necessitates a dataset with multiple sample groups and replicates! Defaulting to LM...")
    LFQ_algo <- LFQ_ALGO <- "LM"
  }
  if (RESCALING == "MSSTATS") {
    reScaling <- RESCALING <- "LM"
  }
  if (length(ALSORUN)) {
    w <- which(ALSORUN != "MSSTATS")
    alsoRun <- alsoRun[w]
    ALSORUN <- ALSORUN[w]
  }
}
stopifnot(nrow(Prot) > 0L,
          nPep_0 > 0L,
          is.character(pg_PepIDs),
          nchar(pg_PepIDs) > 0L,
          length(pg_PepIDs) == 1L,
          pg_PepIDs %in% colnames(Prot),
          inherits(Prot[[pg_PepIDs]], c("character", "integer", "numeric")),
          inherits(N_unique, c("numeric", "integer", "logical")),
          LFQ_ALGO %in% algoSyn$ALGO)
w <- which((ALSORUN %in% algoSyn$ALGO)&(!ALSORUN %in% c(LFQ_ALGO, RESCALING)))
alsoRun <- alsoRun[w]
ALSORUN <- ALSORUN[w]
#
#    I considered including MSstats... but:
#     - MSstats actually does not come with a real protein quantitation algorithm
#     - Its dataProcess() function provides log10 values which... are pre-modelling and DO NOT constitute LFQ
#       (very poor correlation with other methods at both intensity and logFCs level)
#     - MSstats provides average logFCs per contrast, not per sample... (then again, it is the same for limma or MSqRob)
#     - These av. logFCs are based on modelling the data (as part of groupComparison()) and are not consistent with the values provided by dataProcess()
#    -> it is better to export the list of peptides used here for quantitation from this function,
#    and then feed it to a standalone MSstats wrapper further down the road (as part of statistical testing).
#    MSstats would then be best suited for the SAINTexpress treatment, i.e. with its own tab in the report!
#
#    Also note that MSstats is comparatively very, very slow!!!
#    For now the embryo code is still present in the function (see further down below) as reference for when, eventually,
#    a full MSstats workflow is finally added to the stat tests source.
#
# No need to rescale if the algorithm used for LFQ and re-scaling is the same!
if (LFQ_ALGO == RESCALING) {
  skip_reScaling <- TRUE
  reScaling <- "skip"
  RESCALING <- toupper(reScaling)
}
if (!skip_reScaling) {
  reSc_is_topN <- grepl("^TOP[0-9]+$", RESCALING)
}
#
if ("LM" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
  # In-house MaxLFQ-like method (with Levenberg-Marquardt based profiles alignment)
  if ((length(LM_fun) != 1L) ||
      (!is.character(LM_fun)) ||
      is.na(LM_fun) ||
      (!LM_fun %in% c("median", "mean", "weighted.mean"))) {
    warning("Incorrect LM_fun argument, defaulting to \"median\"...")
    LM_fun <- "median"
  }
}
#
# Parameters from optional re-scaling
if (!skip_reScaling) {
  stopifnot(reSc_is_topN ||
              (RESCALING %in% c("LM", "IQ", "LIMPA", "QFEATURES", "MSSTATS")) ||
              (!inherits(try(getFunction(reScaling), silent = TRUE), "try-error")))
  useIntWeights <- FALSE
  reSc_topN <- Inf
  if (reSc_is_topN) {
    reScaling <- "mean"
    reSc_fun <- mean
    reSc_topN <- as.integer(sub("^TOP", "", RESCALING))
    if (!exists("topN_correct") ||
        (length(topN_correct) != 1L) ||
        (!is.logical(topN_correct)) ||
        is.na(topN_correct)) {
      warning("Invalid topN_correct argument, defaulting to TRUE")
      topN_correct <- TRUE
    }
  }
  if (RESCALING == "MAXLFQ") {
    reScaling <- "sum"
  }
  useIntWeights <- reScaling == "weighted.mean"
  if ((!inherits(try(get(reScaling), silent = TRUE), "try-error")) &&
      inherits(get(reScaling), c("standardGeneric", "function"))) {
    reSc_fun <- get(reScaling)
  }
}
#
if (N_unique) {
  # Check that all peptide IDs in Prot[[pg_PepIDs_unique]] are also in Prot[[pg_PepIDs]]
  tmp <- Prot[, c(pg_PepIDs, pg_PepIDs_unique)]
  tmp[[pg_PepIDs]] <- strsplit(tmp[[pg_PepIDs]], ";")
  tmp[[pg_PepIDs_unique]] <- strsplit(tmp[[pg_PepIDs_unique]], ";")
  Prot[[pg_PepIDs]] <- apply(tmp, 1L, \(x) { paste(unique(unlist(x)), collapse = ";") })
}
#
discard_unmod.strict <- FALSE
if (discard_unmod == 2L) {
  discard_unmod.strict <- TRUE
  discard_unmod <- TRUE
}
if (!skipRatios) {
  useContrasts <- (exists("contrasts")) && (!is.function(contrasts))
  if (!useContrasts) {
    if (exists("myContrasts", envir = .GlobalEnv)) {
      contrasts <- myContrasts
      useContrasts <- TRUE
    }
  }
}
if ("Reference" %in% colnames(experimentMap)) {
  experimentMap$Reference <- as.logical(toupper(experimentMap$Reference))
} else {
  if ((!skipRatios) && (!useContrasts)) {
    warning("No reference column was provided, skipping ratios calculation!")
    skipRatios <- TRUE
  }
}
#
noParams <- !exists("param") || (!is.data.frame(param))
# if (noParams && !exists("smplGroups")) {
#   stop("At least one of arguments \"param\" and \"smplGroups\" must be provided!")
# }
if (!skipRatios) {
  if (!noParams) {
    if (!exists("aggrMap")) {
      if (exists("Aggregate.map", envir = .GlobalEnv)) { aggrMap <- Aggregate.map } else {
        warning("Argument \"aggrMap\" must be provided when argument \"param\" is provided!")
        noParams <- TRUE
      }
    }
    if (!exists("aggrList")) {
      if (exists("Aggregate.list", envir = .GlobalEnv)) { aggrList <- Aggregate.list } else {
        warning("Argument \"aggrList\" must be provided when argument \"param\" is provided!")
        noParams <- TRUE
      }
    }
    if (!exists("aggrNames")) {
      if (exists("Aggregates", envir = .GlobalEnv)) { aggrNames <- Aggregates } else {
        warning("Argument \"aggrNames\" must be provided when argument \"param\" is provided!")
        noParams <- TRUE
      }
    }
  }
  if (!noParams) {
    if (!exists("refGroups") && (gsub(";", "", param$Ratios.Ref.Groups) != refGroups$aggregate)) {
      warning("The \"param\" and \"refGroups\" arguments are in disagreement; the former has priority over the latter, so I shall ignore \"refGroups\".")
    }
    if (!exists("ratGroups") && (gsub(";", "", param$Ratios.Groups) != ratGroups$aggregate)) {
      warning("The \"param\" and \"ratGroups\" arguments are in disagreement; the former has priority over the latter, so I shall ignore \"ratGroups\".")
    }
    refGroups <- parse.Param.aggreg(param$Ratios.Ref.Groups,
                                    aggrNames,
                                    aggrMap,
                                    aggrList)
    ratGroups <- parse.Param.aggreg(param$Ratios.Groups,
                                    aggrNames,
                                    aggrMap,
                                    aggrList)
  } else {
    if (!exists("refGroups")) {
      warning("At least one of arguments \"param\" and \"refGroups\" must be provided!")
    }
    if (!exists("ratGroups")) {
      warning("At least one of arguments \"param\" and \"ratGroups\" must be provided!")
    }
    skipRatios <- TRUE
  }
}
Priority <- tolower(substr(Priority, 1L, 3L))
if (!Priority %in% c("rat", "int")) {
  warning("I could not make sense of the value of argument \"Priority\", defaulting to \"Intensities\"")
  Priority <- "int"
}
if (skipRatios && (Priority == "rat")) {
  warning("So, let me get this clear: you want me NOT to calculate ratios BUT to also give ratios priority? Make up your mind! Defaulting to Priority = \"Intensities\"")
  Priority <- "int"
}
#
Pep.Intens.Nms <- paste0(pepInt_Root, mySmpls)
w <- which(Pep.Intens.Nms %in% colnames(Pep))
stopifnot(length(w) > 0L)
mySmpls <- mySmpls[w]
Pep.Intens.Nms <- Pep.Intens.Nms[w]
#
w <- which(apply(Pep[, Pep.Intens.Nms, drop = FALSE], 1L, \(x) {
  sum(is.finite(x))
}) > 0L)
Pep <- Pep[w,]
#
if (!skipRatios) {
  if (useContrasts) {
    Pep.Ratios.Nms <- paste0(pepRat_root, contrasts$Contrast)
  } else {
    Pep.Ratios.Nms <- paste0(pepRat_root, mySmpls)
  }
  Pep.Ratios.Nms <- intersect(Pep.Ratios.Nms, colnames(Pep))
  if (Priority == "rat") { stopifnot(length(Pep.Ratios.Nms) > 0L) }
}
#
Pep$Unmod_Seq <- if ("Sequence" %in% colnames(Pep)) { Pep$Sequence } else { gsub("^_|_$|\\([^\\)]+\\)", "", Pep[[mod_Seq]]) }
#
#
# Make temporary protein-to-peptide list
#  This must be before any filtering of Pep, since we will need to calculate a full peptides intensity sum is re-scaling method = MaxLFQ
#  (Do not re-order)
# - Temporary PG IDs
IDs_vect <- if ("Protein IDs" %in% colnames(Prot)) { Prot$"Protein IDs" } else {
  if ("id" %in% colnames(Prot)) { paste0("ID_", Prot$id) } else {
    paste0("PG#", 1L:nrow(Prot))
  }
}
Prot$temp_IDs <- IDs_vect
# - Full peptides list
quant_pep_IDs <- all_pep_IDss <- setNames(strsplit(Prot[[pg_PepIDs]], ";"),
                                          Prot$temp_IDs)
#
if ((!skip_reScaling) && (RESCALING == "MAXLFQ")) {
  # Equivalent to MaxLFQ
  # Here, we need to know the sum of peptide intensities for each protein group BEFORE any filtering!
  tmp <- listMelt(quant_pep_IDs, ColNames = c("ID", "PG")) 
  tmp$Int <- Pep$Intensity[match(tmp$ID, Pep[[pep_IDs]])]
  tmp <- data.table::data.table(tmp)
  tmp <- tmp[, .(Int = reSc_fun(Int, na.rm = TRUE)), by = .(PG = PG)]
  Prot$"Summed Intensities" <- tmp$Int[match(Prot$temp_IDs,
                                             tmp$PG)]
}
#
# Start filtering peptides
cat(" - Filtering peptides\n")
#
# Use argument mods_to_Exclude to negatively filter peptides:
#  NB:
#   The step controlled by former argument Mods has been commented as redundant.
#   Should this be reverted, the current filtering step, controlled by argument mods_to_Exclude, should remain before the one controlled by argument Mods!
#   Otherwise, the discard_unmod argument may not be used properly.
#   For instance, imagine we first filter peptides to keep only unmodified ones + ones only modified with e.g. "ac" but not "ph"
#   After this, how can we remove unmodified counterparts of observed phospho-peptides, if we do not have the latter's sequences anymore?
#   If changing the order, we must delay the filtering! (as we previously did)
#
modTst <- rep(TRUE, nrow(Pep))
if ((!!exists("mods_to_Exclude")) && nrow(mods_to_Exclude)) {
  if (TESTING) { cat("Identifying peptides with modifications to exclude...\n") }
  mods_to_Exclude$Pattern <- apply(mods_to_Exclude[, c("Mark", "Where")], 1L, \(x) {
    #x <- mods_to_Exclude[1L, c("Mark", "Where")]
    #x <- mods_to_Exclude[2L, c("Mark", "Where")]
    mrk <- x[[1L]]
    wh <- unlist(x[[2L]])
    wh <- grep("^prot[NC]term", wh, invert = TRUE, value = TRUE)
    # For protein N-terminal peptides, we remove them later for those proteins for which the peptide is N-terminal!!!
    #
    wh_A <- grep("^[NC]term", wh, invert = TRUE, value = TRUE) #anywhere
    wh_A <- wh_A[which(wh_A != "")]
    wh_N_trm <- grep("^Nterm$", wh, value = TRUE) # N-terminus, any
    wh_N_trm_sp <- grep("^Nterm_", wh, value = TRUE) # N-terminus, specific
    wh_C_trm <- grep("^Cterm$", wh, value = TRUE) # C-terminus, any
    wh_C_trm_sp <- grep("^Cterm_", wh, value = TRUE) # C-terminus, specific
    res <- c()
    mrkInsrt <- if (discard_unmod.strict) {
      # Optionally, also exclude peptides which could be modified even if they were not found to be...
      # Careful, you will lose A LOT of peptides!
      ""
    } else {
      paste0("\\(([^A-Z]{2},)*", # Changed from "]\\(" to allow catching cases such as "...S(ac,ph)...":
             # Since we are only writing 2 character marks, and excluding capitals from modification marks,
             # this additional optional pattern refinement is sufficient.
             mrk,
             "(,[^A-Z]{2})*\\)") # Same reasoning as above...
    }
    if (length(wh_A)) {
      res <- c(res, paste0("[", paste(wh_A, collapse = ""), "]", mrkInsrt))
    }
    if (length(wh_N_trm)) {
      res <- c(res, paste0("^_", mrkInsrt))
    }
    if (length(wh_N_trm_sp)) {
      wh_N_trm_sp <- sub("^Nterm_", "", wh_N_trm_sp)
      res <- c(res, paste0("^_", mrkInsrt, "[", paste0(wh_N_trm_sp, collapse = ""), "]"))
    }
    if (length(wh_C_trm)) { res <- c(res, paste0(mrkInsrt, "_$")) }
    if (length(wh_C_trm_sp)) {
      wh_C_trm_sp <- sub("^Cterm_", "", wh_C_trm_sp)
      res <- c(res, paste0("[", paste0(wh_C_trm_sp, collapse = ""), "]", mrkInsrt, "_$"))
    }
    if (discard_unmod.strict) {
      # Optionally, also exclude peptides which could be modified even if they were not found to be...
      # Careful, you will lose A LOT of peptides!
      res <- gsub("\\$_", "$", gsub("\\^_", "^", res))
    }
    res <- unique(res)
    return(paste(res, collapse = "|"))
  })
  mods_to_Exclude$"Exclude protein specific" <- vapply(mods_to_Exclude$Where, \(x) {
    sum(grepl("prot[NC]term", unlist(x)))
  }, 1L) > 0L
  fltKol <- mod_Seq
  if (discard_unmod.strict) {
    fltKol <- "Unmod_Seq"
  }
  tst <- nchar(mods_to_Exclude$Pattern)
  mods_to_Exclude$Pattern[which(tst == 0L)] <- NA_character_
  pat <- paste0(mods_to_Exclude$Pattern[which(tst > 0L)], collapse = "|")
  if (nchar(pat)) {
    g <- grep(pat, Pep[[fltKol]])
    modTst[g] <- FALSE
    l <- length(g)
    if (l) {
      if (discard_unmod) {
        # Remove counterpart peptides
        w <- which(Pep$Unmod_Seq %in% unique(Pep$Unmod_Seq[g]))
        l <- length(w)
        if (l) {
          modTst[w] <- FALSE
        }
      }
    }
  }
}
#
Pep <- Pep[which(modTst),]
#
# Filter temp.list based on filtered peptides
quant_pep_IDsA <- listMelt(quant_pep_IDs, ColNames = c("pep", "PG"))
quant_pep_IDsA <- quant_pep_IDsA[which(quant_pep_IDsA$pep %in% Pep[[pep_IDs]]),]
#length(unique(quant_pep_IDsA$pep)) == nrow(Pep)
#
# Now deal with filtering protein-N-terminal mods:
if ((!!exists("mods_to_Exclude")) && nrow(mods_to_Exclude)) {
  Mods2XclTerm <- mods_to_Exclude[which(mods_to_Exclude$"Exclude protein specific"),]
  nrTrm <- nrow(Mods2XclTerm)
  if (nrTrm) {
    Mods2XclTerm$Pattern <- apply(Mods2XclTerm[, c("Mark", "Where")], 1L, \(x) { #x <- Mods2XclTerm[1, c("Mark", "Where")]
      mrk <- x[[1L]]
      wh <- unlist(x[[2L]])
      wh_prt_N_trm <- wh[which(wh == "protNterm")] # protein N-terminus, any
      wh_prt_N_trm_sp <- grep(#"^protNterm_" # (removed the opening "^" in case we have several patterns in one mod)
        "protNterm_", wh, value = TRUE) # protein N-terminus, specific
      wh_prt_C_trm <- wh[which(wh == "protCterm")] # protein C-terminus, any
      wh_prt_C_trm_sp <- grep(#"^protCterm_" # (removed the opening "^" in case we have several patterns in one mod)
        "protCterm_", wh, value = TRUE) # protein C-terminus, specific
      res <- c()
      #
      # !!! Here the effects of discard_unmod.strict can only affect AA-specific cases,
      # at the risk of systematically throwing away every peptide!!!
      mrkInsrt <- paste0("\\(([^A-Z]{2},)*", mrk, "(,[^A-Z]{2})*\\)")
      if (length(wh_prt_N_trm)) {
        res <- c(res, paste0("^_[A-Z]?", mrkInsrt))
        # The "[A-Z]?" here is because some engines may put the mark on the first amino acid or on the "_" N-terminus mark
      }
      if (length(wh_prt_N_trm_sp)) {
        wh_prt_N_trm_sp <- sub("^protNterm_", "", wh_prt_N_trm_sp)
        res <- c(res, paste0("^_",
                             c(mrkInsrt, "")[discard_unmod.strict+1L],
                             "[", paste0(wh_prt_N_trm_sp, collapse = ""), "]"))
      }
      if (length(wh_prt_C_trm)) { res <- c(res, paste0(mrkInsrt, "_$")) }
      if (length(wh_prt_C_trm_sp)) {
        wh_prt_C_trm_sp <- sub("^protCterm_", "", wh_prt_C_trm_sp)
        res <- c(res, paste0("[", paste0(wh_prt_C_trm_sp, collapse = ""), "]",
                             c(mrkInsrt, "")[discard_unmod.strict+1],
                             "_$"))
      }
      res <- unique(res)
      return(paste(res, collapse = "|"))
    })
    Mods2XclTerm <- Mods2XclTerm[which(nchar(Mods2XclTerm$Pattern) > 0L),]
    nrTrm <- nrow(Mods2XclTerm)
  }
  if (nrTrm) {
    modTst2 <- rep(TRUE, nrow(quant_pep_IDsA))
    m1 <- match(quant_pep_IDsA$pep, Pep[[pep_IDs]])
    quant_pep_IDsA$Seq <- Pep$Sequence[m1]
    quant_pep_IDsA$L <- nchar(quant_pep_IDsA$Seq)
    quant_pep_IDsA$PG_seq <- Prot[match(quant_pep_IDsA$PG, Prot$temp_IDs), prim_Seq] # We only filter by first accession!!!
    quant_pep_IDsA$PG_L <- nchar(quant_pep_IDsA$PG_seq)
    # Here I think it makes more sense to deal with each pattern separately
    # It will just be easier, because we need to treat matches to protein N- or C-termini differently
    for (i in 1L:nrTrm) { #i <- 1L
      pat <- paste(Mods2XclTerm$Pattern[i], collapse = "|")
      gy <- grep(pat, Pep[[mod_Seq]])
      if (length(gy)) {
        quant_pep_IDsB <- quant_pep_IDsA[which(quant_pep_IDsA$pep %in% Pep[gy, pep_IDs]),]
        if (grepl("protNterm", Mods2XclTerm$Where[i])) {
          quant_pep_IDsB$Pep_1b <- quant_pep_IDsB$Pep_1a <- substr(quant_pep_IDsB$PG_seq, 1L, quant_pep_IDsB$L)
          w1M <- grep("^M", quant_pep_IDsB$PG_seq) # Allow for N-term methionine loss!
          quant_pep_IDsB$Pep_1b[w1M] <- substr(quant_pep_IDsB$PG_seq[w1M], 2L, quant_pep_IDsB$L[w1M]+1L)
          quant_pep_IDsB$OK <- (quant_pep_IDsB$Seq != quant_pep_IDsB$Pep_1a)&(quant_pep_IDsB$Seq != quant_pep_IDsB$Pep_1b)
        }
        if (grepl("protCterm", Mods2XclTerm$Where[i])) {
          quant_pep_IDsB$Pep_2 <- substr(quant_pep_IDsB$PG_seq, quant_pep_IDsB$PG_L-quant_pep_IDsB$L+1L, quant_pep_IDsB$PG_L)
          quant_pep_IDsB$OK <- quant_pep_IDsB$Seq != quant_pep_IDsB$Pep_2
        }
        modTst2[which(quant_pep_IDsA$pep %in% quant_pep_IDsB$pep[which(!quant_pep_IDsB$OK)])] <- FALSE
      }
    }
    g <- which(!modTst2)
    l2 <- length(g)
    if (l2) {
      w <- which(Pep[[pep_IDs]] %in% quant_pep_IDsA$pep[g])
      if (discard_unmod) {
        # Remove counterpart peptides
        w <- which(Pep$Unmod_Seq %in% quant_pep_IDsA$Seq[g])
        g <- which(quant_pep_IDsA$pep %in% Pep[w, pep_IDs])
        modTst2[g] <- FALSE
      }
      quant_pep_IDsA <- quant_pep_IDsA[which(modTst2),]
      Pep <- Pep[which(Pep[[pep_IDs]] %in% quant_pep_IDsA$pep),]
    }
  }
}
#
# Parse log bases for input and output
# 0 = non-log transformed, 1 = default base (10 for intensities/expression, 2 for ratios))
# Input (peptides-level) intensities
pepInt_log <- as.numeric(pepInt_log)
if (pepInt_log == 1L) { pepInt_log <- 10L }
if (!pepInt_log) {
  # Convert peptide intensities to log
  pepInt_log <- 10L
  if (TESTING) { cat(paste0("Converting input (peptide) intensities to default log", pepInt_log, "...\n")) }
  #origInt <- Pep[, c(pep_IDs, Pep.Intens.Nms)]
  Pep[, Pep.Intens.Nms] <- suppressWarnings(base::log(Pep[, Pep.Intens.Nms],
                                                      pepInt_log))
}
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# Peptide intensities are assumed to be log from now on!!! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
#
# Output (PG-level) expression values
protLFQ_toLog <- as.numeric(protLFQ_toLog)
if (protLFQ_toLog == 1L) { protLFQ_toLog <- 10L }
#
# Input (peptides-level) ratios
if ((!skipRatios) || (Priority == "rat")) {
  pepRat_log <- as.numeric(pepRat_log)
  if (pepRat_log == 1L) { pepRat_log <- 2L }
  if ((Priority == "rat") && (!pepRat_log)) {
    # If we are going to use peptide ratios, and they are not log-transformed, do it now
    pepRat_log <- 2L
    if (TESTING) { cat(paste0("Converting input (peptide) ratios to default log", pepRat_log, "...\n")) }
    Pep[, Pep.Ratios.Nms] <- suppressWarnings(base::log(Pep[, Pep.Ratios.Nms],
                                                        pepRat_log))
  }
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
  # Peptide ratios are assumed to be log from now on!!! #
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
}
#
# Output (PG-level) ratios
if (!skipRatios) {
  protRat_toLog <- as.numeric(protRat_toLog)
  if (protRat_toLog == 1L) { protRat_toLog <- 2L }
}
#
#pepInt_log;protLFQ_toLog;pepRat_log;protRat_toLog
#
Expr.root <- "log10(Expr.)"
Expr.root.full <- paste0(Expr.root, " - ")
#
# Now get peptides average intensity
# When we throw away extra peptides (when we have more than maxN, we want to use in priority the high intensity ones),
# we will use this to decide which to keep.
tmp1 <- data.table::as.data.table(Pep[, c(pep_IDs, Pep.Intens.Nms), drop = FALSE])
tmp1 <- data.table::melt(tmp1, id.vars = pep_IDs)
colnames(tmp1)[1L] <- "pepID"
tmp1 <- tmp1[which(is.finite(tmp1$value)),]
tmp1 <- tmp1[, .(value = mean(value)), by = .(pepID = pepID)]
Pep$avgPepInt <- tmp1$value[match(Pep[[pep_IDs]], tmp1$pepID)]
#
# If N_unique > 0: update list of peptides to use for quantitation
if (N_unique) {
  quant_pep_IDsU <- listMelt(strsplit(Prot[[pg_PepIDs_unique]], ";"), Prot$temp_IDs, c("pep", "PG"))
  quant_pep_IDsU <- quant_pep_IDsU[which(quant_pep_IDsU$pep %in% quant_pep_IDsA$pep),]
  quant_pep_IDsA <- quant_pep_IDsA[which(!quant_pep_IDsA$pep %in% quant_pep_IDsU$pep),]
}
quant_pep_IDsA <- data.table::as.data.table(quant_pep_IDsA[, c("pep", "PG")])
nms <- Prot$temp_IDs
if (N_unique) {
  #
  quant_pep_IDsU <- data.table::as.data.table(quant_pep_IDsU[, c("pep", "PG")])
  #
  # Sort by decreasing intensities and convert to list
  # - Uniques
  quant_pep_IDsU$int <- Pep$avgPepInt[match(quant_pep_IDsU$pep, Pep[[pep_IDs]])]
  quant_pep_IDsU <- quant_pep_IDsU[order(quant_pep_IDsU$int, decreasing = TRUE),]
  quant_pep_IDsU <- quant_pep_IDsU[, list(pep = list(pep)), by = list(PG = PG)]
  quant_pep_IDsU <- setNames(quant_pep_IDsU$pep, quant_pep_IDsU$PG)
  # - Rest
  quant_pep_IDsA$int <- Pep$avgPepInt[match(quant_pep_IDsA$pep, Pep[[pep_IDs]])]
  quant_pep_IDsA <- quant_pep_IDsA[order(quant_pep_IDsA$int, decreasing = TRUE),]
}
quant_pep_IDsA <- quant_pep_IDsA[, list(pep = list(pep)), by = list(PG = PG)]
quant_pep_IDsA <- setNames(quant_pep_IDsA$pep, quant_pep_IDsA$PG)
if (N_unique) {
  # Re-add missing PGs and re-order
  # - Uniques
  w <- which(!nms %in% names(quant_pep_IDsU))
  quant_pep_IDsU[nms[w]] <- lapply(nms[w], \(x) {})
  quant_pep_IDsU <- quant_pep_IDsU[nms]
  # - Rest
  w <- which(!nms %in% names(quant_pep_IDsA))
  quant_pep_IDsA[nms[w]] <- lapply(nms[w], \(x) {})
  quant_pep_IDsA <- quant_pep_IDsA[nms]
  #
  # Test for length
  tstU <- lengths(quant_pep_IDsU)
  tstA <- lengths(quant_pep_IDsA)
  w1 <- which((tstU < N_unique)&(tstA > 0L))
  if (length(w1)) {
    quant_pep_IDsU[w1] <- lapply(w1, \(x) { c(quant_pep_IDsU[[x]], quant_pep_IDsA[[x]]) })
  }
  quant_pep_IDs <- quant_pep_IDsU
} else {
  quant_pep_IDs <- quant_pep_IDsA
}
# Filter by minN and maxN
tst <- lengths(quant_pep_IDs)
w1 <- which(tst < minN)
w2 <- which(tst > maxN)
if (length(w1)) {
  quant_pep_IDs[w1] <- lapply(w1, \(x) { })
}
if (length(w2)) {
  quant_pep_IDs[w2] <- lapply(w2, \(x) {
    quant_pep_IDs[[x]][1L:maxN]
  })
}
#  
#tst1 <- lengths(quant_pep_IDs)
#tst2 <- vapply(quant_pep_IDs, \(x) { length(unique(x)) }, 1L)
#sum(tst2 != tst1)
#
#sum(!Prot$temp_IDs %in% names(quant_pep_IDs))
#
nPep_1 <- nrow(Pep)
nPep_2 <- length(unique(unlist(quant_pep_IDs)))
nRemoved1 <- nPep_0-nPep_1
nRemoved2 <- nPep_1-nPep_2
msg <- ""
if (nRemoved1) {
  msg <- paste0("Excluding ", nRemoved1, " (", round(100*nRemoved1/nPep_0), "%) ")
  msg <- if (discard_unmod.strict) {
    paste0(msg, "whose stoichiometry could be affected by PTMs excluded from quantitation!")
  } else {
    paste0(msg, "bearing PTMs excluded from quantitation",
           c("!", ", as well as the latter's unmodified counterpart peptides!")[discard_unmod+1L])
  }
}
if (nRemoved2) {
  msg <- paste0(msg, "\nAlso discarding ", nRemoved2, " (", round(100*nRemoved2/nPep_0), "%) supernumerary peptides (keeping higher intensity ones).")
}
warning(msg)
#
# Zig-zag order
# Function adapted from https://www.r-bloggers.com/2020/12/going-parallel-understanding-load-balancing-in-r/
nZZ <- max(c(2L, N.clust))
zigzag_ord <- \(x, n = nZZ) {
  #x <- quant_pep_IDs
  ord <- data.frame(Original = seq_along(x),
                    Length = lengths(x))
  ord <- ord[order(ord$Length, decreasing = TRUE),]
  ord$NewOrd <- rep(c(seq(1L, n), seq(n, 1L)), length = length(x))
  ord <- ord[order(ord$NewOrd),]
  ord$NewOrd <- 1L:nrow(ord)
  ord <- ord[order(ord$Original),]
  return(ord)
}
ord <- zigzag_ord(quant_pep_IDs)
ord$ID <- names(quant_pep_IDs)
nuOrd <- order(ord$NewOrd)
quant_pep_IDs <- quant_pep_IDs[nuOrd]
#
Pep <- Pep[which(Pep[[pep_IDs]] %in% unlist(quant_pep_IDs)),] # Update Pep (is this necessary?)
#
# Get summary method:
if (!exists("Weights")) {
  if ("weighted.mean" %in% c(LM_fun, reScaling)) {
    warning("The summary method is \"weighted mean\", but no weights were provided...?\nSetting all weights to 1...")
  }
  Weights <- "Weights"
  Pep[[Weights]] <- 1L
}
#
# Optional intensity weights - so a peptide's weight is a function of its (non log-transformed) intensity
if (!exists("useIntWeights") ||
    (!is.logical(useIntWeights)) ||
    (length(useIntWeights) != 1L) ||
    is.na(useIntWeights)) {
  useIntWeights <- FALSE
}
if (useIntWeights &&
    (sum(c(reScaling, LM_fun) %in% c("mean", "weighted.mean")))) {
  tmp <- data.table::data.table(pepInt_log^Pep[, Pep.Intens.Nms, drop = FALSE])
  tmp$pepID <- Pep[[pep_IDs]]
  tmp <- data.table::melt(tmp, id.vars = "pepID")
  tmp <- tmp[which(is.finite(tmp$value)),]
  tmp <- tmp[which(tmp$value > 0),]
  tmp <- tmp[, .(value = mean(value)), by = .(pepID)]
  Pep$useIntWeights <- tmp$value[match(Pep[[pep_IDs]], tmp$pepID)]
  Pep[[Weights]] <- Pep[[Weights]]*Pep$useIntWeights # Update weights and method
  if (LM_fun == "mean") { LM_fun <- "weighted.mean" }
  if (reScaling == "mean") {
    reScaling <- "weighted.mean"
    reSc_fun <- get(reScaling) # Update
  }
}
#
# If Priority == "Ratios"
#  - we re-calculate pep intensities using ratios, so they exactly - as opposed to approximately - reflect ratios.
#  - we will then use these ratios-based intensities to calculate protein expression downstream.
# This may be used for cases such as SILAC in MaxQuant, where SILAC ratios are more accurate for relative quant
# than the ratio of peak intensities.
if (Priority == "rat") {
  tmpPep <- Pep[, c(pep_IDs, Weights, Pep.Ratios.Nms, Pep.Intens.Nms)]
  intSums <- rowSums(10^tmpPep[, Pep.Intens.Nms], na.rm = TRUE)
  ratGroups$samples <- lapply(ratGroups$values, \(x) {
    experimentMap[which(experimentMap[[ratGroups$column]] == x), experimentMap_Samples_col]
  })
  ratGroups$refSamples <- lapply(ratGroups$values, \(x) {
    experimentMap[which((experimentMap[[ratGroups$column]] == x)&(experimentMap$Reference)), experimentMap_Samples_col]
  })
  ratGroups$newInt <- lapply(1L:length(ratGroups$values), \(x) { #x <- 1L
    allIntCol <- paste0(pepInt_Root, ratGroups$samples[[x]])
    allRatCol <- paste0(pepRat_root, ratGroups$samples[[x]])
    w <- which(allRatCol %in% colnames(tmpPep))
    myRatCol <- allRatCol[w]
    myIntCol <- allIntCol[w]
    refIntCol <- paste0(pepInt_Root, ratGroups$refSamples[[x]])
    refInt <- rowMeans(tmpPep[, refCol, drop = FALSE], na.rm = TRUE)
    allInt <- tmpPep[, allIntCol]
    i1 <- rowMeans(allInt, na.rm = TRUE)
    newInt <- tmpPep[, myRatCol]/base::log(pepInt_log, pepRat_log) + refInt
    colnames(newInt) <- myIntCol
    #View(newInt)
    #View(allInt[, myIntCol])
    allInt[, myIntCol] <- newInt
    i2 <- rowMeans(allInt, na.rm = TRUE)
    allInt <- allInt - i2 + i1
    return(allInt)
  })
  newInt <- do.call(cbind, ratGroups$newInt)
  tmpPep[, colnames(newInt)] <- newInt
  tmpPep <- tmpPep[, which(!colnames(tmpPep) %in% Pep.Ratios.Nms)]
} else {
  tmpPep <- Pep[, c(pep_IDs, Weights, Pep.Intens.Nms)]
}
#
# Re-order by intensity
tmpPep$AvgInt <- rowMeans(tmpPep[, Pep.Intens.Nms, drop = FALSE], na.rm = TRUE)
tmpPep <- tmpPep[order(tmpPep$AvgInt, decreasing = TRUE),]
quant_pep_IDs2 <- listMelt(quant_pep_IDs, ColNames = c("id", "PG"))
quant_pep_IDs2$mtch <- match(quant_pep_IDs2$id, tmpPep[[pep_IDs]])
quant_pep_IDs2 <- quant_pep_IDs2[order(quant_pep_IDs2$mtch, decreasing = FALSE),]
quant_pep_IDs2 <- data.table::as.data.table(quant_pep_IDs2)
quant_pep_IDs2 <- quant_pep_IDs2[, list(IDs = list(id)), by = list(PG = PG)]
quant_pep_IDs2 <- setNames(quant_pep_IDs2$IDs, quant_pep_IDs2$PG)
quant_pep_IDs[names(quant_pep_IDs2)] <- quant_pep_IDs2[names(quant_pep_IDs2)] # There are some empty entries in quant_pep_IDs: proteins with no peptides eligible for quant
rm(quant_pep_IDs2)
#
# Quantitation:
#  - Quantitation step
cat(" - Calculating profile\n")
quntNms <- sub(topattern(pepInt_Root), Expr.root.full, Pep.Intens.Nms)
# Below: object used by most methods and by reScaling
tmpPep2 <- listMelt(quant_pep_IDs, ColNames = c("pepID", "PG"))
tmpPep2[, Pep.Intens.Nms] <- tmpPep[match(tmpPep2$pepID, tmpPep[[pep_IDs]]), Pep.Intens.Nms]
#
allQuants <- list() # Used for reScaling!
cat("   method =", LFQ_algo, "\n")
l <- length(ALSORUN)
if (l) {
  cat(paste0("   (also running alternative method", c("", "s")[(l > 1L)+1L], " ", paste(ALSORUN, collapse = " / "), ")\n"))
}
if ("LM" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
  cat("    ----- running LM method -----\n")
  # Check cluster
  source(parSrc)
  #
  tmpFl1 <- tempfile(fileext = ".rds")
  tmpFl2 <- tempfile(fileext = ".rds")
  readr::write_rds(tmpPep, tmpFl1)
  readr::write_rds(quant_pep_IDs, tmpFl2)
  parallel::clusterExport(parClust,
                          list("diffLog", "LFQ.lm", "tmpFl1", "tmpFl2"),
                          envir = environment()) # (Overwrite the package's versions in case we have local ones)
  invisible(parallel::clusterCall(parClust, \() {
    library(stats)
    library(minpack.lm)
    assign("tmpPep", readr::read_rds(tmpFl1), envir = .GlobalEnv)
    assign("quant_pep_IDs", readr::read_rds(tmpFl2), envir = .GlobalEnv)
    return()
  }))
  cat("            Starting calculations...\n")
  lmDat <- parallel::parLapply(parClust,
                               quant_pep_IDs,
                               LFQ.lm,
                               InputTabl = tmpPep,
                               id = pep_IDs,
                               IntensCol = Pep.Intens.Nms,
                               Summary.method = LM_fun,
                               Summary.weights = Weights,
                               Min.N = minN,
                               Max.N = maxN,
                               reNorm = 1L)
  cat("            Dealing with the results...\n")
  lmDat <- as.data.frame(do.call(rbind, lmDat))
  colnames(lmDat) <- sub(topattern(pepInt_Root), "", colnames(lmDat))
  lmDat <- lmDat[match(names(quant_pep_IDs), row.names(lmDat)),
                 mySmpls]
  rownames(lmDat) <- names(quant_pep_IDs)
  allQuants$LM <- lmDat
  cat("            Done!\n")
}
if (sum(c("IQ", "LIMPA", "QFEATURES", "MSSTATS") %in% c(LFQ_ALGO, RESCALING, ALSORUN))) {
  # These algorithms do not take too long to run,
  # thus they could be run either be run as LFQ, with or without subsequent re-scaling,
  # or used to provide re-scaling for another LFQ method.
  if ("IQ" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
    cat("    ----- running iq MaxLFQ method -----\n")
    tmp4 <- tmpPep2
    w <- which(colnames(tmp4) %in% Pep.Intens.Nms)
    colnames(tmp4)[w] <- sub(topattern(pepInt_Root), "", colnames(tmp4)[w])
    tmp4 <- melt(tmp4, id.vars = c("pepID", "PG"))
    tmp4 <- tmp4[which(is.finite(tmp4$value)),]
    tmp4 <- list(protein_list = tmp4$PG,
                 sample_list = tmp4$variable,
                 id = tmp4$pepID,
                 quant = tmp4$value)
    if (pepInt_log != 2L) {
      tmp4$quant <- tmp4$quant/base::log(2L, pepInt_log)
    }
    l <- unique(lengths(tmp4))
    stopifnot(sum(c("protein_list", "sample_list", "id", "quant") %in% names(tmp4)) == 4L,
              length(l) == 1L,
              l > 0L)
    # tmp4DF <- data.frame(protein_list = tmp4$protein_list,
    #                      sample_list = tmp4$sample_list,
    #                      id = tmp4$id,
    #                      quant = tmp4$quant)
    # View(tmp4DF)
    # sum(is.na(tmp4DF$quant)) == 0L
    # sum(!is.finite(tmp4DF$quant)) == 0L
    iqObj <- iq::fast_MaxLFQ(tmp4)
    iqDat <- iqObj$estimate
    if (pepInt_log != 2L) {
      iqDat <- iqDat/base::log2(pepInt_log)
    }
    # Check columns/row order
    iqDat <- as.data.frame(iqDat)[match(names(quant_pep_IDs), row.names(iqDat)),
                                  mySmpls]
    rownames(iqDat) <- names(quant_pep_IDs)
    allQuants$IQ <- iqDat
  }
  if ("LIMPA" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
    cat("    ----- running limpa method -----\n")
    tmp4 <- as.matrix(tmpPep2[, Pep.Intens.Nms])
    colnames(tmp4) <- sub(topattern(pepInt_Root), "", colnames(tmp4))
    #
    # LIMPA and LIMMA need LOG2!!!...
    if (pepInt_log != 2L) {
      tmp4 <- tmp4/base::log(2L, pepInt_log)
    }
    #
    # Deal with infinites which may break dpcCN!
    w <- which(!is.finite(tmp4))
    tmp4[w] <- NA_real_
    #
    dpcfit <- limpa::dpcCN(tmp4)
    #dpcfit <- limpa::dpc(tmp4)
    #limpa::plotDPC(dpcfit)
    dpcObj <- limpa::dpcQuant(tmp4, tmpPep2$PG, dpc = dpcfit) 
    dpcDat <- dpcObj$E
    #
    # ... but we can work with the log base we want, so convert back!
    if (pepInt_log != 2L) {
      dpcDat <- dpcDat/log2(pepInt_log)
    }
    #
    # Code to replace every value based on 0 observations with NA:
    # w <- which(dpcObj$other$n.observations == 0L)
    # dpcDat[w] <- NA_real_
    # However, Gordon Smyth does not recommend doing it.
    # For a discussion about the proper usage of limpa, why there are no missing values in its output quant matrix,
    # or why it may output repeated values in a row,
    # please see:
    # https://github.com/SmythLab/limpa/issues/5#event-22032865875
    #
    # Check columns/row order
    dpcDat <- as.data.frame(dpcDat)[match(names(quant_pep_IDs), row.names(dpcDat)),
                                    mySmpls]
    #
    rownames(dpcDat) <- names(quant_pep_IDs)
    allQuants$LIMPA <- dpcDat
  }
  if ("QFEATURES" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
    cat("    ----- running QFeatures method -----\n")
    tmp4 <- tmpPep2
    w <- which(colnames(tmp4) %in% Pep.Intens.Nms)
    tmp5 <- as.matrix(tmp4[, w])
    wNotOK <- which(!is.finite(tmp5))
    tmp5[wNotOK] <- NA_real_
    tmp4[, w] <- tmp5
    colnames(tmp4)[w] <- sub(topattern(pepInt_Root), "", colnames(tmp4)[w])
    tmp4$pepID <- 1L:nrow(tmp4) # Otherwise aggregateFeatures throws an error
    QFeatObj <- QFeatures::readQFeatures(assayData = tmp4,
                                         fnames = "pepID",
                                         quantCols = mySmpls,
                                         name = "peptides")
    QFeatObj <- QFeatures::aggregateFeatures(QFeatObj,
                                             i = "peptides",
                                             fcol = "PG",
                                             na.rm = FALSE,
                                             name = "PG")
    qfDat <- SummarizedExperiment::assay(QFeatObj[["PG"]])
    # Check columns/row order
    qfDat <- as.data.frame(qfDat)[match(names(quant_pep_IDs), row.names(qfDat)),
                                  mySmpls]
    rownames(qfDat) <- names(quant_pep_IDs)
    allQuants$QFEATURES <- qfDat
  }
  if ("MSSTATS" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
    tmp4 <- tmpPep2
    tmp4[, Pep.Intens.Nms] <- pepInt_log^tmp4[, Pep.Intens.Nms]
    w <- which(colnames(tmp4) %in% Pep.Intens.Nms)
    colnames(tmp4)[w] <- sub(topattern(pepInt_Root), "", colnames(tmp4)[w])
    tmp4 <- dfMelt(tmp4, c("pepID", "ProteinName", "Run", "Intensity"), c("pepID", "PG"))
    tmp4 <- tmp4[which(is.finite(tmp4$Intensity)),]
    m1 <- match(tmp4$Run, experimentMap$Ref.Sample.Aggregate)
    tmp4$Condition <- experimentMap[m1,
                                    gsub(";", "", param$Volcano.plots.Aggregate.Level)]
    tmp4$BioReplicate <- experimentMap$Replicate[m1]
    m2 <- match(tmp4$pepID, Pep[[pep_IDs]])
    tmp4$PeptideModifiedSequence <- Pep$`Modified sequence`[m2]
    tmp4$PrecursorCharge <- 0 # (we currently do not differentiate between different charge states of the same peptidoforms)
    tmp4$FragmentIon <- "NA"
    tmp4$ProductCharge <- 0
    tmp4$IsotopeLabelType <- "L"
    uID <- unique(tmp4$ProteinName)
    tmp4$tmp <- as.character(match(tmp4$ProteinName, uID))
    tmp4$PeptideSequence  <- do.call(paste, c(tmp4[, c("PeptideModifiedSequence", "tmp")], sep = "pg"))
    # ... not a mistake, cf. how MSstats defines feature!
    #
    tmp4$PeptideModifiedSequence <- NULL
    MSstats_list <- MSstats::dataProcess(tmp4,
                                         normalization = FALSE,
                                         summaryMethod = "TMP"#,
                                         #numberOfCores = N.clust # Unfortunately not available on Windows!
                                         )
    MSstats_quant <- MSstats_list$ProteinLevelData
    MSstats_quant <- reshape::cast(MSstats_quant,
                                   Protein ~ originalRUN,
                                   mean,
                                   value = "LogIntensities",
                                   na.rm = TRUE)
    rownames(MSstats_quant) <- MSstats_quant$Protein
    MSstats_quant$Protein <- NULL
    MSstats_quant <- MSstats_quant[match(names(quant_pep_IDs), rownames(MSstats_quant)),
                                   mySmpls]
    rownames(MSstats_quant) <- names(quant_pep_IDs)
    allQuants$MSSTATS <- MSstats_quant
  }
}
# tst1 <- vapply(names(allQuants), \(nm) { sum(!colnames(allQuants[[nm]]) %in% mySmpls) == 0L }, TRUE)
# tst2 <- vapply(allQuants, nrow, 1L) == nrow(Prot)
# tst3 <- lapply(allQuants, rownames)
# tst3 <- as.data.frame(do.call(cbind, tst3))
# colnames(tst3) <- names(allQuants)
# tst3$Orig <- names(quant_pep_IDs)
# tst3a <- apply(tst3, 1L, \(x) { length(unique(x)) })
# w <- which(tst3a > 1L)
# View(tst3[w,])
res2 <- allQuants[[LFQ_ALGO]]
colnames(res2) <- paste0(Expr.root.full, colnames(res2))
#
#########################################
# Code to compare the different methods #
#########################################
# nr <- nrow(lmDat)
# tst <- lapply(quntNms, \(k) {
#   i <- cleanNms(sub(".* - ", "", k))
#   data.frame(Sample = i,
#              X = c(rep(lmDat[[k]], 2L),
#                    rep(iqDat[[k]], 3L),
#                    dpcDat[[k]]),
#              X_method = c(rep("LM", 3*nr),
#                           rep("iq", 2*nr),
#                           rep("limpa", nr)),
#              Y = c(iqDat[[k]], dpcDat[[k]], qfDat[[k]],
#                    dpcDat[[k]], qfDat[[k]],
#                    qfDat[[k]]),
#              Y_method = c(rep("iq", nr), rep("limpa", nr), rep("QFeatures", nr),
#                           rep("limpa", nr), rep("QFeatures", nr),
#                           rep("QFeatures", nr)))
# })
# tst <- do.call(rbind, tst)
# tst <- tst[which((!is.na(tst$X))&(!is.na(tst$Y))),]
# tst$X_method <- factor(tst$X_method, levels = c("LM", "iq", "limpa", "QFeatures"))
# tst$Y_method <- factor(tst$Y_method, levels = c("LM", "iq", "limpa", "QFeatures"))
# plot <- ggplot(tst) +
#   scattermore::geom_scattermore(aes(x = X, y = Y, color = Sample)) +
#   facet_grid(Sample~X_method*Y_method) + coord_fixed() +
#   theme_bw() + theme(strip.text.y = element_text(angle = 0)) +
#   ggtitle("log10 abundances")
# poplot(plot, 12, 22)
# #
# em <- experimentMap
# em$Replicate <- as.integer(em$Replicate)
# em <- em[order(em$Replicate),]
# f0 <- \(i, y) { #i <- 1L
#   w1 <- which(contr[i,] == 1L)
#   w0 <- which(contr[i,] == -1L)
#   w1 <- which(em[[smplGroups$column]] == colnames(contr)[w1])
#   w0 <- which(em[[smplGroups$column]] == colnames(contr)[w0])
#   if (Nested) {
#     rep <- unique(em$Replicate[c(w1, w0)])
#     res <- lapply(rep, \(r) { #r <- 1L
#       w1_ <- w1[which(em$Replicate[w1] == r)]
#       w0_ <- w0[which(em$Replicate[w0] == r)]
#       if ((length(w1_) != 1L) || (length(w0_) != 1L)) { return() }
#       smpl1_ <- em[w1_, experimentMap_Samples_col]
#       smpl0_ <- em[w0_, experimentMap_Samples_col]
#       return(y[[paste0("log10(Expr.) - ", smpl1_)]] - y[[paste0("log10(Expr.) - ", smpl0_)]])
#     })
#     res <- do.call(cbind, res)
#   } else {
#     smpls1 <- em[w1, experimentMap_Samples_col]
#     smpls0 <- em[w0, experimentMap_Samples_col]
#     rf <- rowMeans(y[, paste0("log10(Expr.) - ", smpls0)], na.rm = TRUE)
#     res <- sweep(y[, paste0("log10(Expr.) - ", smpls1)], 1L, rf, "-")
#   }
#   res <- rowMeans(res, na.rm = TRUE)
#   #res <- res/base::log(pepRat_log, pepInt_log) # res3e is still pepInt_log base
#   return(res)
# }
# res3a <- as.data.frame(do.call(cbind, setNames(lapply(1L:nrow(contr), \(i) { f0(i, y = lmDat) }), colnames(res3e))))
# res3b <- as.data.frame(do.call(cbind, setNames(lapply(1L:nrow(contr), \(i) { f0(i, y = iqDat) }), colnames(res3e))))
# res3c <- as.data.frame(do.call(cbind, setNames(lapply(1L:nrow(contr), \(i) { f0(i, y = dpcDat) }), colnames(res3e))))
# res3d <- as.data.frame(do.call(cbind, setNames(lapply(1L:nrow(contr), \(i) { f0(i, y = qfDat) }), colnames(res3e))))
# tst <- lapply(colnames(res3e), \(k) {
#   data.frame(Group = cleanNms(k),
#              X = c(rep(res3a[[k]], 4L),
#                    rep(res3b[[k]], 3L),
#                    rep(res3c[[k]], 2L),
#                    res3d[[k]]),
#              X_method = c(rep("LM", 4L*nr),
#                           rep("iq", 3L*nr),
#                           rep("limpa", 2L*nr),
#                           rep("QFeatures", nr)),
#              Y = c(res3b[[k]], res3c[[k]], res3d[[k]], res3e[[k]],
#                    res3c[[k]], res3d[[k]], res3e[[k]],
#                    res3d[[k]], res3e[[k]],
#                    res3e[[k]]),
#              Y_method = c(rep("iq", nr), rep("limpa", nr), rep("QFeatures", nr), rep("MSstats", nr),
#                           rep("limpa", nr), rep("QFeatures", nr), rep("MSstats", nr),
#                           rep("QFeatures", nr), rep("MSstats", nr),
#                           rep("MSstats", nr)))
# })
# tst <- do.call(rbind, tst)
# tst <- tst[which((!is.na(tst$X))&(!is.na(tst$Y))),]
# tst$X_method <- factor(tst$X_method, levels = c("LM", "iq", "limpa", "QFeatures", "MSstats"))
# tst$Y_method <- factor(tst$Y_method, levels = c("LM", "iq", "limpa", "QFeatures", "MSstats"))
# plot <- ggplot(tst) +
#   scattermore::geom_scattermore(aes(x = X, y = Y, color = Group)) +
#   facet_grid(Group~X_method*Y_method) + coord_fixed() +
#   theme_bw() + theme(strip.text.y = element_text(angle = 0)) +
#   ggtitle("avg. logFCs")
# poplot(plot, 12, 22)
#######
# End #
#######
#
# Re-scaling step
if (!skip_reScaling) {
  cat(" - Rescaling expression values, method =", RESCALING, "\n")
  if (RESCALING %in% c("MAXLFQ", "IQ", "LIMPA", "QFEATURES", "LM", "MSSTATS")) {
    rescVal <- if (RESCALING == "MAXLFQ") {
      data.frame(PG = Prot$temp_IDs,
                 Value = base::log(Prot$"Summed Intensities", pepInt_log))
    } else {
      data.frame(PG = rownames(allQuants[[RESCALING]]),
                 Value = rowMeans(allQuants[[RESCALING]][, quntNms, drop = FALSE], na.rm = TRUE))
    }
  } else {
    rescVal <- tmpPep2[, c("pepID", "PG")]
    m <- match(rescVal$pepID, Pep[[pep_IDs]])
    rescVal$avgPepInt <- Pep$avgPepInt[m]
    if (reScaling == "weighted.mean") {
      rescVal$Weights <- Pep$Weights[m]
    }
    if (reSc_is_topN || is.finite(reSc_topN)) {
      rescVal$Rank <- as.integer(stats::ave(rescVal$PG, rescVal$PG, FUN = seq_along)) # (thanks chatGPT...)
    }
    if (is.finite(reSc_topN)) {
      rescVal <- rescVal[which(rescVal$Rank <= reSc_topN),]
    }
    rescVal <- data.table::as.data.table(rescVal)
    if (reSc_is_topN && topN_correct) {
      # Not sure about this...
      Md <- median(rescVal$avgPepInt, na.rm = TRUE)
      tst <- as.data.frame(rescVal[, .(Median = median(avgPepInt, na.rm = TRUE)), by = .(Rank = Rank)])
      tst$Rank <- as.integer(tst$Rank)
      for (i in tst$Rank) {
        wi <- which(rescVal$Rank == i)
        rescVal$avgPepInt[wi] <- rescVal$avgPepInt[wi]-tst$Median[match(i, tst$Rank)]
      }
      rescVal$avgPepInt <- rescVal$avgPepInt + Md
    }
    rescVal <- if (reScaling == "weighted.mean") {
      rescVal[, .(Value = weighted.mean(avgPepInt, Weights, na.rm = TRUE)), by = .(PG = PG)]
    } else {
      rescVal[, .(Value = reSc_fun(avgPepInt, na.rm = TRUE)), by = .(PG = PG)] 
    }
    rescVal <- as.data.frame(rescVal)
  }
  rescVal <- rescVal[match(row.names(res2), rescVal$PG),]
  # Re-scale
  wN <- which(!row.names(res2) %in% rescVal$PG)
  wY <- which(row.names(res2) %in% rescVal$PG)
  stopifnot(length(wY) > 0L)
  mY <- match(row.names(res2)[wY], rescVal$PG)
  res2[wN, quntNms] <- NA_real_
  currValY <- rowMeans(res2[wY, quntNms, drop = FALSE], na.rm = TRUE)
  res2[wY, quntNms] <- sweep(res2[wY, quntNms, drop = FALSE], 1L, rescVal$Value[mY]-currValY, "+")
}
# Calculate ratios
if (!skipRatios) {
  if (useContrasts) {
    tmpRt <- make_Rat2(res2,
                       contrasts = contrasts,
                       refGroups = refGroups,
                       pepInt_log,
                       pepRat_log,
                       experiment.map = experimentMap,
                       int.root = Expr.root.full,
                       rat.root = pepRat_root)
    if (is.null(pepRat_root)) { colnames(tmpRt) <- paste0("log", pepRat_log, "FC - ", colnames(tmpRt)) }
  } else {
    tmpRt <- make_Rat(res2,
                      Priority = Priority,
                      pepInt_log,
                      pepRat_log,
                      refGroups = refGroups,
                      experiment.map = experimentMap,
                      int.root = Expr.root.full,
                      rat.root = pepRat_root)
    tmpRt <- tmpRt$Data
  }
  kol <- colnames(tmpRt)
  res2[, kol] <- tmpRt[, kol]
}
# Finally: this should be the very last step to avoid confusions:
# If output should not be log-transformed, de-log!
if (!protLFQ_toLog) { # De-log LFQ values
  w <- grep(topattern(Expr.root.full), colnames(res2))
  g <- colnames(res2)[w]
  res2[, g] <- pepInt_log^(res2[, g])
  colnames(res2)[w] <- sub(topattern(Expr.root.full),
                           "Expr. - ",
                           colnames(res2)[w])
}
if (!skipRatios) {
  w <- grep(topattern(pepRat_root), colnames(res2))
  g <- colnames(res2)[w]
  if (!protRat_toLog) { # If necessary de-log
    res2[, g] <- pepRat_log^(res2[, g])
    colnames(res2)[w] <- sub(topattern(pepRat_root),
                             "Ratio - ",
                             colnames(res2)[w])
  } else {
    # If changing log base from peptides to proteins:
    if (protRat_toLog != pepRat_log) {
      res2[, g] <- res2[, g]/base::log(protRat_toLog, pepRat_log)
      colnames(res2)[w] <- sub(topattern(pepRat_root),
                               paste0("log", protRat_toLog, "FC - "),
                               colnames(res2)[w])
    }
  }
}
for (i in 1L:ncol(res2)) {
  if (!is.numeric(res2[[i]])) {
    stop(paste0("I would expect the class of column ", i, " to be numeric! Investigate!"))
    #res2[[i]] <- as.numeric(res2[[i]])
  }
}
res2$"Peptides IDs used for quantitation" <- vapply(quant_pep_IDs, paste, "", collapse = ";")
ord2 <- ord[nuOrd,]
res2 <- res2[order(ord2$Original),]
#
nC <- dim(res2)[2L]
res3 <- as.data.frame(matrix(rep(NA, nrow(Prot)*nC), ncol = nC)) 
rownames(res3) <- Prot$temp_IDs
colnames(res3) <- colnames(res2)
wY <- which(rownames(res3) %in% rownames(res2))
mY <- match(rownames(res3)[wY], rownames(res2))
res3[wY,] <- res2[mY,]
#
RES <- list(Data = res3)
if ("LM" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
  RES$LM <- lmDat
}
if ("IQ" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
  RES$iq_MaxLFQ <- iqObj
}
if ("LIMPA" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
  RES$EList_obj <- dpcObj
}
if ("QFEATURES" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
  RES$QFeatures_obj <- QFeatObj
}
if ("MSSTATS" %in% c(LFQ_ALGO, RESCALING, ALSORUN)) {
  RES$MSstats_list <- MSstats_list
}
quantData_list %<o% RES
