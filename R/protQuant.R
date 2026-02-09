#' protQuant
#'
#' @description
#' A function to calculate estimated protein group Expression (inter-protein group quantitation vectors) and Ratios (optional, intra-protein group fold changes) from individual peptide values.
#' The input is assumed to already be normalized!\cr
#' 
#' @param Prot Protein/Protein groups table. A data.frame, which must contain at least peptides IDs and primary sequence of the first accession in each group.
#' @param Pep Peptides table. A data.frame. If it contains unmodified sequences then these are expected to be in column name "Sequence". Argument mod_Seq controls the name of the column used for modified sequence.
#' @param pg_PepIDs Name of the Protein/Protein groups table's peptide IDs column. Default = "Peptide IDs"
#' @param pg_PepIDs_unique Used only if mode = "PreferUnique". Name of the Protein/Protein groups table's unique peptide IDs column. Default = "Unique peptide IDs"
#' @param pep_IDs The name of the Peptides table's IDs column. Default = "id"
#' @param N_unique Logical or numeric. If FALSE or 0, the function just uses the "pg_PepIDs" argument. If non null (default), a second column name of unique peptide IDs should be provided using the "pg_PepIDs_unique" argument. If at least N_unique unique peptides are present, then only those will be used for quantitation, otherwise as many razor/shared peptides as necessary will be added (sorted by decreasing average intensities).
#' @param LFQ_algo Algorithm used to compute average profiles. One of:\cr
#'  - "LM": Levenberg-Marquardt method (backend = minpack.lm::nls.lm()); very similar to MaxLFQ. Normalized peptide relative profiles are aligned using Levenberg-Marquardt then summarized.\cr
#'  - "iq": iq's fast implementation of MaxLFQ, backend = iq::fast_MaxLFQ()\cr
#'  - "limpa" or "DPC": backend = limpa::dpcQuant(); current default\cr
#'  - "QFeatures": backend = QFeatures::aggregateFeatures()\cr
#' Note that limpa and QFeatures do not "just produce quantitative values": both output specific objects with additional information which can be used for downstream analysis (respectively with limma or msqrob2).\cr
# @param LFQ_stab Logical. Stabilize ratios? Default = TRUE
#' @param reScaling Optional summary method for re-scaling. May be one of:\cr
#'  - "median",\cr
#'  - "mean",\cr
#'  - "weighted.mean" (requires the "Weights" argument),\cr
#'  - "max",\cr
#'  - "sum" (not recommended),\cr
#'  - "MaxLFQ": like sum, but the value used is the value before any peptides are filtered out (not recommended),\cr
#'  - "topN", where N should be the maximum number of peptides to average (e.g. "top3" - do not use "topN" as there is no default value for N!)
#'  - the name of any valid value of LFQ_algo, which allows re-scaling any LFQ algorithm using the scaling provided by another.
#'  - (you could also use the name of any other available averaging function, this should work in principle assuming similar syntax)
#' @param topN_correct Logical, default = TRUE. In the case where we are using more than one peptide for the re-scaling step, should we correct for systematic peptide intensity biases between peptides of rank 1, 2, 3 and so on and so forth?
#' @param minN Integer, default = 1. How many peptides should at least be present for quantitation? Values lower than 1 are increased to 1. May not be higher than N_unique!
#' @param maxN Integer or Inf, default = 50. Up to how many peptides should we use for the Levenberg-Marquardt procedure (used only for LFQs)? Using too many peptides can be an issue, e.g. with huge proteins like Titin. Default = 50. The most intense peptides will be selected. May not be lower than N_unique!
#' @param LM_fun How should normalized profiles be averaged? One of "median" (default), "mean", or "weighted.mean" (the latter which uses the Weights and useIntWeights arguments).
#' @param Weights Length 1 character, a valid column name of Pep containing user-defined individual peptide weights. Used if LM_fun or reScaling are set to "weighted.mean" (for the former only if LFQ_algo = "LM").
#' @param useIntWeights Logical, default = FALSE. Ignored unless LM_fun or reScaling are set to "mean" or "weighted.mean". If TRUE, will take into account individual peptide intensities when calculating average profile for that step (thus, it will actually be a weighted mean regardless). DOES NOT replace the optional, user-provided Weights, but instead multiplies the former by new intensity-based factors.
#' @param Priority One of "Intensities" (default) or "Ratios" (some flexibility in spelling is allowed). In some rare cases, such as a SILAC dataset processed with MaxQuant, we will want to correct intensities - prior to running the main algorithm - so their ratios reflect the more accurate ratios directly measured by the search engine.
#' @param skipRatios Default = FALSE. If TRUE, ratios will not be calculated. 
#' @param expMap Map of the experiment map.
#' @param expMap_Samples_col Names of the single samples column in the experiment's map. Default = "Sample"
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
#' 
#' #' @details
#' This function is meant to work from normalised peptide intensities and:\cr
#'  - calculates a quantitative profile for each protein,\cr
#'  - maximizes the usage of peptides-level ratios information.\cr
#' Specifically, it can be broken down into 3 sub-steps:\cr
#'   a) filtering peptides eligible for quantitation\cr
#'   b) calculating a relative protein profile across samples:\cr
#'    - LM: this uses Levenberg-Marquardt to align peptide profiles, then summarizes them, and is broadly similar to MaxLFQ.\cr
#'    - iq: the iq package's implementation of MaxLFQ, specifically the iq::fast_MaxLFQ() function.\cr
#'    - limpa: uses the recently introduced limpa package's dpcQuant() function.
#'    - QFeatures: uses the QFeatures::aggregateFeatures() function.
#'   c) Optional re-scaling, i.e. "anchoring" the resulting relative profile to a scale value, such that different proteins are ranked relatively in a manner which - whilst imperfect because of widely different individual peptide detectabilities - would still be correct if the latter could be corrected for, or would otherwise operate under assumptions at least minimizing the expected error. See the reScaling argument.\cr
#'    You may also rescale any quantitation (as specified in LFQ_algo) with the scale provided by another (using the reScaling argument), thus allowing for using the relative quantitation and scaling of any two methods.
#'    To re-scale the same way as MaxLFQ, use "sum" - but see N.B. below!
#'    iBAQ - which would be the equivalent of using "sum" then dividing by the number of observable peptides, but the additional step could easily be added outside this function.
#'\cr
#' NB on re-scaling:
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
#' 
#' @examples
#' temp <- protQuant(Prot = PG, pg_PepIDs = Pep4Quant, Pep = pep, pep_IDs = "New Peptide ID",
#'                   expMap = Exp.map,
#'                   refGroups = Ratios.Ref.Groups,
#'                   pepInt_Root = pep.ref, pepRat_root = pep.ratios.ref,
#'                   pepInt_log = FALSE, pepRat_log = 2,
#'                   protLFQ_toLog = TRUE, protRat_toLog = TRUE,
#'                   mods_to_Exclude = Mod2Xclud, mod_Seq = "Modified sequence",
#'                   minN = 2)
#' 
#' @export

protQuant <- function(Prot,
                      Pep,
                      pg_PepIDs = "Peptide IDs",
                      pg_PepIDs_unique = "Unique peptide IDs",
                      pep_IDs = "id",
                      N_unique,
                      LFQ_algo = "limpa",
                      #LFQ_stab = TRUE,
                      reScaling,
                      topN_correct = TRUE,
                      minN = 1,
                      maxN = 50,
                      LM_fun = "median",
                      Weights,
                      useIntWeights = FALSE,
                      Priority = "Intensities",
                      skipRatios = FALSE,
                      expMap,
                      expMap_Samples_col = "Sample",
                      refGroups,
                      ratGroups,
                      #smplGroups,
                      param,
                      aggrMap,
                      aggrList,
                      aggrNames,
                      pepInt_Root,
                      pepRat_root,
                      pepInt_log = FALSE,
                      pepRat_log = 2,
                      protLFQ_toLog = TRUE,
                      protRat_toLog = TRUE,
                      #Mods,
                      mods_to_Exclude,
                      mod_Seq = "Modified sequence",
                      discard_unmod = TRUE,
                      prim_Seq = "Sequence (1st accession)",
                      ref_Mode = "2",
                      cl,
                      N.clust,
                      N.reserved = 1) {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::protQuant);cl <- parClust;TESTING <- TRUE
  #Prot = PG; pg_PepIDs = "Razor peptide IDs"; pg_PepIDs_unique = "Unique peptide IDs";Pep = pep; pep_IDs = "id";useIntWeights = FALSE;expMap = Exp.map; param = Param;pepInt_Root = pep.ref[length(pep.ref)];pepRat_root = pep.ratios.ref[length(pep.ratios.ref)];pepInt_log = FALSE; pepRat_log = 2;protLFQ_toLog = TRUE; protRat_toLog = TRUE;mods_to_Exclude = Mod2Xclud; discard_unmod = discard_unmod;minN = 1; N.clust = N.clust;Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1]
  #Prot = PG; pg_PepIDs = Pep4Quant; Pep = pep; useIntWeights = FALSE; skipRatios = !MakeRatios; expMap = Exp.map; refGroups = RefGrp; ratGroups = RatGrp; pepInt_Root = paste0(int.col, " - "); pepRat_root = paste0(rat.cols["Original"], " - "); pepInt_log = FALSE; pepRat_log = 2; protLFQ_toLog = TRUE; protRat_toLog = TRUE; mods_to_Exclude = Mod2Xclud; minN = 1; N.clust = N.clust; Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1]
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  if (misFun(ref_Mode)) { ref_Mode <- 2 }
  if (!as.numeric(ref_Mode) %in% 1:2) { ref_Mode <- 2 }
  ref_Mode <- as.character(ref_Mode)
  #
  # Check arguments
  mySmpls <- expMap[[expMap_Samples_col]]
  nPep_0 <- nrow(Pep)
  minN <- max(c(1, suppressWarnings(abs(as.integer(minN)))))
  maxN <- max(c(1, minN, suppressWarnings(abs(as.numeric(maxN)))))
  N_unique <- max(c(0, suppressWarnings(abs(as.numeric(N_unique)))))
  if (is.finite(maxN)) { maxN <- suppressWarnings(as.integer(maxN)) }
  if (is.finite(N_unique)) { N_unique <- suppressWarnings(as.integer(N_unique)) }
  stopifnot(length(N_unique) == 1,
            length(maxN) == 1,
            length(minN) == 1,
            sum(!is.na(c(N_unique, minN, maxN))) == 3)
  cat("Using at least", minN, "and up to", maxN, "peptides...\n")
  if (N_unique) {
    stopifnot(pg_PepIDs_unique %in% colnames(Prot))
    if (pg_PepIDs == pg_PepIDs_unique) {
      N_unique <- 0
    }
  }
  if (N_unique) {
    cat("If available, up to", N_unique, "unique (= proteotypic) peptides will be used...\n")
  }
  LFQ_algo <- gsub(" -_\\.", "", LFQ_algo)
  LFQ_ALGO <- toupper(LFQ_algo)
  stopifnot(nrow(Prot) > 0,
            nPep_0 > 0,
            "character" %in% class(pg_PepIDs),
            nchar(pg_PepIDs) > 0,
            length(pg_PepIDs) == 1,
            pg_PepIDs %in% colnames(Prot),
            sum(c("character", "integer", "numeric") %in% class(Prot[[pg_PepIDs]])) > 0,
            sum(c("numeric", "integer", "logical") %in% class(N_unique)) > 0,
            LFQ_ALGO %in% c("LM",
                            "IQ", #"MAXLFQ",
                            "LIMPA", "DPCQUANT",
                            "QFEATURES", "MSQROB", "MSQROB2"
                            #, "MSSTATS"
                            )
            )
  if (length(mySmpls) == 1) {
    if (LFQ_ALGO %in% c("LIMPA", "DPCQUANT",  "QFEATURES", "MSQROB", "MSQROB2")) {
      warning("Only one sample, LFQ_algo can only be one of LM or iq, defaulting to LM...")
    }
    LFQ_algo <- LFQ_ALGO <- "LM"
  }
  # Resolve synonyms / check names
  if (LFQ_ALGO == "IQ") {
    LFQ_algo <- "iq"
  }
  if (LFQ_ALGO %in% c("LIMPA", "DPCQUANT")) {
    LFQ_algo <- "limpa"
    LFQ_ALGO <- "LIMPA"
  }
  if (LFQ_ALGO %in% c("QFEATURES", "MSQROB", "MSQROB2")) {
    LFQ_algo <- "QFeatures"
    LFQ_ALGO <- "QFEATURES"
  }
  #    I considered including MSstats... but:
  #     - MSstats actually does not come with a real protein quantitation algorithm
  #     - Its dataProcess() function provides log10 values which... are pre-modelling and DO NOT constitute LFQ
  #       (very poor correlation with other methods at both intensity and logFCs level)
  #     - MSstats provides average logFCs per contrast, not per sample.
  #     - The latter are based on modelling the data (as part of groupComparison()) and are not consistent with the values provided by dataProcess()
  #    -> it is better to export the list of peptides used here for quantitation from this function,
  #    and then feed it to a standalone MSstats wrapper further down the road (as part of statistical testing).
  #    MSstats would then best suited for the SAINTexpress treatment, i.e. with its own tab in the report!
  #
  #    Also note that MSstats is comparatively very very slow!!!
  #    For now the embryo code is still present in the function (see further down below)
  #    for when, eventually, a full MSstats workflow is finally added to the stats test available here.
  #
  if (LFQ_ALGO == "LM") {
    # In-house MaxLFQ-like method (with Levenberg-Marquardt based profiles alignment)
    if ((length(LM_fun) != 1)||
        (!"character" %in% class(LM_fun))||
        (is.na(LM_fun))||
        (!LM_fun %in% c("median", "mean", "weighted.mean"))) {
      warning("Incorrect LM_fun argument, defaulting to \"median\"...")
      LM_fun <- "median"
    }
  }
  #
  if (misFun(reScaling)) { reScaling <- "skip" }
  reScaling <- gsub(" ", "", reScaling)
  skip_reScaling <- reScaling == "skip"
  RESCALING <- toupper(reScaling)
  # Resolve synonyms
  if (RESCALING == "IQ") {
    reScaling <- "iq"
  }
  if (RESCALING %in% c("LIMPA", "DPCQUANT")) {
    reScaling <- "limpa"
    RESCALING <- "LIMPA"
  }
  if (RESCALING %in% c("QFEATURES", "MSQROB", "MSQROB2")) {
    reScaling <- "QFeatures"
    RESCALING <- "QFEATURES"
  }
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
  # Parameters from optional re-scaling
  if (!skip_reScaling) {
    stopifnot(reSc_is_topN||
                (RESCALING %in% c("LM", "IQ", "LIMPA", "QFEATURES"))||
                ((!"try-error" %in% class(try(get(reScaling), silent = TRUE)))&&
                   (sum(c("standardGeneric", "function") %in% class(get(reScaling))) > 0)))
    useIntWeights <- FALSE
    reSc_topN <- Inf
    if (reSc_is_topN) {
      reScaling <- "mean"
      reSc_fun <- mean
      reSc_topN <- as.integer(gsub("^TOP", "", RESCALING))
      if ((misFun(topN_correct))||(length(topN_correct) != 1)||(!is.logical(topN_correct))||(is.na(topN_correct))) {
        warning("Invalid topN_correct argument, defaulting to TRUE")
        topN_correct <- TRUE
      }
    }
    if (RESCALING == "MAXLFQ") {
      reScaling <- "sum"
    }
    useIntWeights <- reScaling == "weighted.mean"
    if ((!"try-error" %in% class(try(get(reScaling), silent = TRUE)))&&
        (sum(c("standardGeneric", "function") %in% class(get(reScaling))) > 0)) {
      reSc_fun <- get(reScaling)
    }
  }
  #
  if (N_unique) {
    # Check that all peptide IDs in Prot[[pg_PepIDs_unique]] are also in Prot[[pg_PepIDs]]
    tmp <- Prot[, c(pg_PepIDs, pg_PepIDs_unique)]
    tmp[[pg_PepIDs]] <- strsplit(tmp[[pg_PepIDs]], ";")
    tmp[[pg_PepIDs_unique]] <- strsplit(tmp[[pg_PepIDs_unique]], ";")
    Prot[[pg_PepIDs]] <- apply(tmp, 1 , function(x) { paste(unique(unlist(x)), collapse = ";") })
  }
  #
  discard_unmod.strict <- FALSE
  if (discard_unmod == 2) {
    discard_unmod.strict <- TRUE
    discard_unmod <- TRUE
  }
  if ("Reference" %in% colnames(expMap)) {
    expMap$Reference <- as.logical(toupper(expMap$Reference))
  } else {
    if (!skipRatios) {
      warning("No reference column was provided, skipping ratios calculation!")
      skipRatios <- TRUE
    }
  }
  #
  Priority <- tolower(substr(Priority, 1, 3))
  if (!Priority %in% c("rat", "int")) {
    warning("I could not make sense of the value of argument \"Priority\", defaulting to \"Intensities\"")
    Priority <- "int"
  }
  if ((skipRatios)&&(Priority == "rat")) {
    warning("So, let me get this clear: you want me NOT to calculate ratios BUT to also give ratios priority? Make up your mind! Defaulting to Priority = \"Intensities\"")
    Priority <- "int"
  }
   # if ((misFun(param))&&(misFun(smplGroups))) {
  #   stop("At least one of arguments \"param\" and \"smplGroups\" must be provided!")
  # }
  if ((!misFun(param))&&(!skipRatios)) {
    if ((misFun(param))&&(misFun(refGroups))) {
      stop("At least one of arguments \"param\" and \"refGroups\" must be provided!")
    }
    if ((misFun(param))&&(misFun(ratGroups))) {
      stop("At least one of arguments \"param\" and \"ratGroups\" must be provided!")
    }
    if ((!misFun(refGroups))&&(gsub(";", "", param$Ratios.Ref.Groups) != refGroups$aggregate)) {
      warning("The \"param\" and \"refGroups\" arguments are in disagreement; the former has priority over the latter, so I shall ignore \"refGroups\".")
    }
    # "Soft default" values
    if ((misFun(aggrMap))&&(exists("Aggregate.map", envir = .GlobalEnv))) { aggrMap <- Aggregate.map }
    if ((misFun(aggrList))&&(exists("Aggregate.list", envir = .GlobalEnv))) { aggrList <- Aggregate.list }
    if ((misFun(aggrNames))&&(exists("Aggregates", envir = .GlobalEnv))) { aggrNames <- Aggregates }
    #
    refGroups <- proteoCraft::parse.Param.aggreg(param$Ratios.Ref.Groups,
                                                 aggrNames,
                                                 aggrMap,
                                                 aggrList)
    if ((!misFun(ratGroups))&&(gsub(";", "", param$Ratios.Groups) != ratGroups$aggregate)) {
      warning("The \"param\" and \"ratGroups\" arguments are in disagreement; the former has priority over the latter, so I shall ignore \"ratGroups\".")
    }
    ratGroups <- proteoCraft::parse.Param.aggreg(param$Ratios.Groups,
                                                 aggrNames,
                                                 aggrMap,
                                                 aggrList)
    # smplGroups <- proteoCraft::parse.Param.aggreg(param$Volcano.plots.Aggregate.Level,
    #                                               aggrNames,
    #                                               aggrMap,
    #                                               aggrList)
  }
  #
  Pep.Intens.Nms <- paste0(pepInt_Root, mySmpls)
  w <- which(Pep.Intens.Nms %in% colnames(Pep))
  stopifnot(length(w) > 0)
  mySmpls <- mySmpls[w]
  Pep.Intens.Nms <- Pep.Intens.Nms[w]
  if (!skipRatios) {
    Pep.Ratios.Nms <- paste0(pepRat_root, mySmpls)
    Pep.Ratios.Nms <- Pep.Ratios.Nms[which(Pep.Ratios.Nms %in% colnames(Pep))]
    if (Priority == "rat") { stopifnot(length(Pep.Ratios.Nms) > 0) }
  }
  #
  if ("Sequence" %in% colnames(Pep)) {
    Pep$Unmod_Seq <- Pep$Sequence
  } else {
    Pep$Unmod_Seq <- gsub("^_|_$|\\([^\\)]+\\)", "", Pep[[mod_Seq]])
  }
  #
  #
  # Make temporary protein-to-peptide list
  #  This must be before any filtering of Pep, since we will need to calculate a full peptides intensity sum is re-scaling method = MaxLFQ
  #  (Do not re-order)
  # - Temporary PG IDs
  if ("Protein IDs" %in% colnames(Prot)) { IDs_vect <- Prot$"Protein IDs" } else {
    if ("id" %in% colnames(Prot)) { IDs_vect <- paste0("ID_", Prot$id) } else {
      IDs_vect <- paste0("PG#", 1:nrow(Prot))
    }
  }
  Prot$temp_IDs <- IDs_vect
  # - Full peptides list
  quant_pep_IDss <- all_pep_IDss <- setNames(strsplit(Prot[[pg_PepIDs]], ";"),
                                           Prot$temp_IDs)
  #
  if ((!skip_reScaling)&&(RESCALING == "MAXLFQ")) {
    # Equivalent to MaxLFQ
    # Here, we need to know the sum of peptide intensities for each protein group BEFORE any filtering!
    tmp <- proteoCraft::listMelt(quant_pep_IDss, ColNames = c("ID", "PG")) 
    tmp$Int <- Pep$Intensity[match(tmp$ID, Pep[[pep_IDs]])]
    tmp <- data.table::data.table(tmp)
    tmp <- tmp[, .(Int = reSc_fun(Int, na.rm = TRUE)), by = .(PG = PG)]
    Prot$"Summed Intensities" <- tmp$Int[match(Prot$temp_IDs,
                                               tmp$PG)]
  }
  #
  # Start filtering peptides
  cat(" - filtering peptides\n")
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
  if ((!misFun(mods_to_Exclude))&&(nrow(mods_to_Exclude))) {
    if (TESTING) { cat("Identifying peptides with modifications to exclude...\n") }
    mods_to_Exclude$Pattern <- apply(mods_to_Exclude[, c("Mark", "Where")], 1, function(x) {
      #x <- mods_to_Exclude[1, c("Mark", "Where")]
      #x <- mods_to_Exclude[2, c("Mark", "Where")]
      mrk <- x[[1]]
      wh <- unlist(x[[2]])
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
      if (discard_unmod.strict) {
        # Optionally, also exclude peptides which could be modified even if they were not found to be...
        # Careful, you will lose A LOT of peptides!
        mrkInsrt <- ""
      } else {
        mrkInsrt <- paste0("\\(([^A-Z]{2},)*", # Changed from "]\\(" to allow catching cases such as "...S(ac,ph)...":
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
        wh_N_trm_sp <- gsub("^Nterm_", "", wh_N_trm_sp)
        res <- c(res, paste0("^_", mrkInsrt, "[", paste0(wh_N_trm_sp, collapse = ""), "]"))
      }
      if (length(wh_C_trm)) { res <- c(res, paste0(mrkInsrt, "_$")) }
      if (length(wh_C_trm_sp)) {
        wh_C_trm_sp <- gsub("^Cterm_", "", wh_C_trm_sp)
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
    mods_to_Exclude$"Exclude protein specific" <- vapply(mods_to_Exclude$Where, function(x) {
      sum(grepl("prot[NC]term", unlist(x)))
    }, 1) > 0
    fltKol <- mod_Seq
    if (discard_unmod.strict) {
      fltKol <- "Unmod_Seq"
    }
    tst <- nchar(mods_to_Exclude$Pattern)
    mods_to_Exclude$Pattern[which(tst == 0)] <- NA
    pat <- paste0(mods_to_Exclude$Pattern[which(tst > 0)], collapse = "|")
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
  # # Use argument Mods to positively filter peptides:
  # # Commented: essentially redundant with mods_to_Exclude - the latter which makes more sense
  # #
  # if ((!is.logical(Mods))&&(length(Mods))) {
  #   w <- which(modTst)
  #   if (TESTING) { cat("Filtering out peptidoforms not eligible for quantitation...\n") }
  #   tst <- gsub("\\)\\(", ",", gsub("^\\(|\\)$", "", gsub(paste(c("^_", "_$", AA), collapse = "|"), "", Pep[w, mod_Seq])))
  #   tst <- vapply(strsplit(tst, ","), function(x) {
  #     x <- unlist(x)
  #     return(length(x[which(!x %in% Mods)]))
  #   }, 1) == 0
  #   #sum(!tst)
  #   #Pep[w, mod_Seq][which(!tst)]
  #   l <- sum(!tst)
  #   if (l) {
  #     warning(paste0("Excluding ", l, " (", round(100*l/nPep_0),
  #                    "%) peptides with modifications not specifically included in quantitation parameters!"))
  #     modTst[w[which(!tst)]] <- FALSE
  #   }
  # }
  #
  Pep <- Pep[which(modTst),]
  #
  # Filter temp.list based on filtered peptides
  quant_pep_IDssA <- proteoCraft::listMelt(quant_pep_IDss, ColNames = c("pep", "PG"))
  quant_pep_IDssA <- quant_pep_IDssA[which(quant_pep_IDssA$pep %in% Pep[[pep_IDs]]),]
  #length(unique(quant_pep_IDssA$pep)) == nrow(Pep)
  #
  # Now deal with filtering protein-N-terminal mods:
  if ((!misFun(mods_to_Exclude))&&(nrow(mods_to_Exclude))) {
    Mods2XclTerm <- mods_to_Exclude[which(mods_to_Exclude$"Exclude protein specific"),]
    nrTrm <- nrow(Mods2XclTerm)
    if (nrTrm) {
      Mods2XclTerm$Pattern <- apply(Mods2XclTerm[, c("Mark", "Where")], 1, function(x) { #x <- Mods2XclTerm[1, c("Mark", "Where")]
        mrk <- x[[1]]
        wh <- unlist(x[[2]])
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
          wh_prt_N_trm_sp <- gsub("^protNterm_", "", wh_prt_N_trm_sp)
          res <- c(res, paste0("^_",
                               c(mrkInsrt, "")[discard_unmod.strict+1],
                               "[", paste0(wh_prt_N_trm_sp, collapse = ""), "]"))
        }
        if (length(wh_prt_C_trm)) { res <- c(res, paste0(mrkInsrt, "_$")) }
        if (length(wh_prt_C_trm_sp)) {
          wh_prt_C_trm_sp <- gsub("^protCterm_", "", wh_prt_C_trm_sp)
          res <- c(res, paste0("[", paste0(wh_prt_C_trm_sp, collapse = ""), "]",
                               c(mrkInsrt, "")[discard_unmod.strict+1],
                               "_$"))
        }
        res <- unique(res)
        return(paste(res, collapse = "|"))
      })
      Mods2XclTerm <- Mods2XclTerm[which(nchar(Mods2XclTerm$Pattern) > 0),]
      nrTrm <- nrow(Mods2XclTerm)
    }
    if (nrTrm) {
      modTst2 <- rep(TRUE, nrow(quant_pep_IDssA))
      m1 <- match(quant_pep_IDssA$pep, Pep[[pep_IDs]])
      quant_pep_IDssA$Seq <- Pep$Sequence[m1]
      quant_pep_IDssA$L <- nchar(quant_pep_IDssA$Seq)
      quant_pep_IDssA$PG_seq <- Prot[match(quant_pep_IDssA$PG, Prot$temp_IDs), prim_Seq] # We only filter by first accession!!!
      quant_pep_IDssA$PG_L <- nchar(quant_pep_IDssA$PG_seq)
      # Here I think it makes more sense to deal with each pattern separately
      # It will just be easier, because we need to treat matches to protein N- or C-termini differently
      for (i in 1:nrTrm) { #i <- 1
        pat <- paste(Mods2XclTerm$Pattern[i], collapse = "|")
        gy <- grep(pat, Pep[[mod_Seq]])
        if (length(gy)) {
          quant_pep_IDssB <- quant_pep_IDssA[which(quant_pep_IDssA$pep %in% Pep[gy, pep_IDs]),]
          if (grepl("protNterm", Mods2XclTerm$Where[i])) {
            quant_pep_IDssB$Pep_1b <- quant_pep_IDssB$Pep_1a <- substr(quant_pep_IDssB$PG_seq, 1, quant_pep_IDssB$L)
            w1M <- grep("^M", quant_pep_IDssB$PG_seq) # Allow for N-term methionine loss!
            quant_pep_IDssB$Pep_1b[w1M] <- substr(quant_pep_IDssB$PG_seq[w1M], 2, quant_pep_IDssB$L[w1M]+1)
            quant_pep_IDssB$OK <- (quant_pep_IDssB$Seq != quant_pep_IDssB$Pep_1a)&(quant_pep_IDssB$Seq != quant_pep_IDssB$Pep_1b)
          }
          if (grepl("protCterm", Mods2XclTerm$Where[i])) {
            quant_pep_IDssB$Pep_2 <- substr(quant_pep_IDssB$PG_seq, quant_pep_IDssB$PG_L-quant_pep_IDssB$L+1, quant_pep_IDssB$PG_L)
            quant_pep_IDssB$OK <- quant_pep_IDssB$Seq != quant_pep_IDssB$Pep_2
          }
          modTst2[which(quant_pep_IDssA$pep %in% quant_pep_IDssB$pep[which(!quant_pep_IDssB$OK)])] <- FALSE
        }
      }
      g <- which(!modTst2)
      l2 <- length(g)
      if (l2) {
        w <- which(Pep[[pep_IDs]] %in% quant_pep_IDssA$pep[g])
        if (discard_unmod) {
          # Remove counterpart peptides
          w <- which(Pep$Unmod_Seq %in% quant_pep_IDssA$Seq[g])
          g <- which(quant_pep_IDssA$pep %in% Pep[w, pep_IDs])
          modTst2[g] <- FALSE
        }
        quant_pep_IDssA <- quant_pep_IDssA[which(modTst2),]
        Pep <- Pep[which(Pep[[pep_IDs]] %in% quant_pep_IDssA$pep),]
      }
    }
  }
  #
  # Parse log bases for input and output
  # 0 = non-log transformed, 1 = default base (10 for intensities/expression, 2 for ratios))
  # Input (peptides-level) intensities
  pepInt_log <- as.numeric(pepInt_log)
  if (pepInt_log == 1) { pepInt_log <- 10}
  if (!pepInt_log) {
    # Convert peptide intensities to log
    pepInt_log <- 10
    if (TESTING) { cat(paste0("Converting input (peptide) intensities to default log", pepInt_log, "...\n")) }
    #origInt <- Pep[, c(pep_IDs, Pep.Intens.Nms)]
    Pep[, Pep.Intens.Nms] <- suppressWarnings(log(Pep[, Pep.Intens.Nms],
                                                  pepInt_log))
  }# else {
  #   if (LFQ_ALGO == "MSSTATS") {
  #     origInt <- Pep[, c(pep_IDs, Pep.Intens.Nms)]
  #     origInt[, Pep.Intens.Nms] <- pepInt_log^origInt[, Pep.Intens.Nms]
  #   }
  # }
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
  # Peptide intensities are assumed to be log from now on!!! #
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
  #
  # Output (PG-level) expression values
  protLFQ_toLog <- as.numeric(protLFQ_toLog)
  if (protLFQ_toLog == 1) { protLFQ_toLog <- 10 }
  #
  # Input (peptides-level) ratios
  if (!skipRatios) {
    pepRat_log <- as.numeric(pepRat_log)
    if (pepRat_log == 1) { pepRat_log <- 2 }
    if ((Priority == "rat")&&(!pepRat_log)) {
      # If we are going to use peptide ratios, and they are not log-transformed, do it now
      pepRat_log <- 2
      if (TESTING) { cat(paste0("Converting input (peptide) ratios to default log", pepRat_log, "...\n")) }
      Pep[, Pep.Ratios.Nms] <- suppressWarnings(log(Pep[, Pep.Ratios.Nms],
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
    if (protRat_toLog == 1) { protRat_toLog <- 2 }
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
  colnames(tmp1)[1] <- "pepID"
  tmp1 <- tmp1[which(is.all.good(tmp1$value, 2)),]
  tmp1 <- tmp1[, .(value = mean(value)), by = .(pepID = pepID)]
  Pep$avgPepInt <- tmp1$value[match(Pep[[pep_IDs]], tmp1$pepID)]
  #
  # In PreferUnique mode: update list of peptides to use for quantitation
  if (N_unique) {
    quant_pep_IDssU <- proteoCraft::listMelt(strsplit(Prot[[pg_PepIDs_unique]], ";"), Prot$temp_IDs, c("pep", "PG"))
    quant_pep_IDssU <- quant_pep_IDssU[which(quant_pep_IDssU$pep %in% quant_pep_IDssA$pep),]
    quant_pep_IDssA <- quant_pep_IDssA[which(!quant_pep_IDssA$pep %in% quant_pep_IDssU$pep),]
  }
  quant_pep_IDssA <- data.table::as.data.table(quant_pep_IDssA[, c("pep", "PG")])
  nms <- Prot$temp_IDs
  if (N_unique) {
    #
    quant_pep_IDssU <- data.table::as.data.table(quant_pep_IDssU[, c("pep", "PG")])
    #
    # Sort by decreasing intensities and convert to list
    # - Uniques
    quant_pep_IDssU$int <- Pep$avgPepInt[match(quant_pep_IDssU$pep, Pep[[pep_IDs]])]
    quant_pep_IDssU <- quant_pep_IDssU[order(quant_pep_IDssU$int, decreasing = TRUE),]
    quant_pep_IDssU <- quant_pep_IDssU[, list(pep = list(pep)), by = list(PG = PG)]
    quant_pep_IDssU <- setNames(quant_pep_IDssU$pep, quant_pep_IDssU$PG)
    # - Rest
    quant_pep_IDssA$int <- Pep$avgPepInt[match(quant_pep_IDssA$pep, Pep[[pep_IDs]])]
    quant_pep_IDssA <- quant_pep_IDssA[order(quant_pep_IDssA$int, decreasing = TRUE),]
  }
  quant_pep_IDssA <- quant_pep_IDssA[, list(pep = list(pep)), by = list(PG = PG)]
  quant_pep_IDssA <- setNames(quant_pep_IDssA$pep, quant_pep_IDssA$PG)
  if (N_unique) {
    # Re-add missing PGs and re-order
    # - Uniques
    w <- which(!nms %in% names(quant_pep_IDssU))
    quant_pep_IDssU[nms[w]] <- lapply(nms[w], function(x) {})
    quant_pep_IDssU <- quant_pep_IDssU[nms]
    # - Rest
    w <- which(!nms %in% names(quant_pep_IDssA))
    quant_pep_IDssA[nms[w]] <- lapply(nms[w], function(x) {})
    quant_pep_IDssA <- quant_pep_IDssA[nms]
    #
    # Test for length
    tstU <- vapply(quant_pep_IDssU, length, 1)
    tstA <- vapply(quant_pep_IDssA, length, 1)
    w1 <- which((tstU < N_unique)&(tstA > 0))
    if (length(w1)) {
      quant_pep_IDssU[w1] <- lapply(w1, function(x) { c(quant_pep_IDssU[[x]], quant_pep_IDssA[[x]]) })
    }
    quant_pep_IDss <- quant_pep_IDssU
  } else {
    quant_pep_IDss <- quant_pep_IDssA
  }
  # Filter by minN and maxN
  tst <- vapply(quant_pep_IDss, length, 1)
  w1 <- which(tst < minN)
  w2 <- which(tst > maxN)
  if (length(w1)) {
    quant_pep_IDss[w1] <- lapply(w1, function(x) { })
  }
  if (length(w2)) {
    quant_pep_IDss[w2] <- lapply(w2, function(x) {
      quant_pep_IDss[[x]][1:maxN]
    })
  }
  #  
  #tst1 <- vapply(quant_pep_IDss, length, 1)
  #tst2 <- vapply(quant_pep_IDss, function(x) { length(unique(x)) }, 1)
  #sum(tst2 != tst1)
  #
  #sum(!Prot$temp_IDs %in% names(quant_pep_IDss))
  #
  nPep_1 <- nrow(Pep)
  nPep_2 <- length(unique(unlist(quant_pep_IDss)))
  nRemoved1 <- nPep_0-nPep_1
  nRemoved2 <- nPep_1-nPep_2
  msg <- ""
  if (nRemoved1) {
    msg <- paste0("Excluding ", nRemoved1, " (", round(100*nRemoved1/nPep_0), "%) ")
    if (discard_unmod.strict) {
      msg <- paste0(msg, "whose stoichiometry could be affected by PTMs excluded from quantitation!")
    } else {
      msg <- paste0(msg, "bearing PTMs excluded from quantitation",
                    c("!", ", as well as the latter's unmodified counterpart peptides!")[discard_unmod+1])
    }
  }
  if (nRemoved2) {
    msg <- paste0(msg, "\nAlso discarding ", nRemoved2, " (", round(100*nRemoved2/nPep_0), "%) supernumerary peptides (keeping higher intensity ones).")
  }
  warning(msg)
  #
  # Zig-zag order
  # Function adapted from https://www.r-bloggers.com/2020/12/going-parallel-understanding-load-balancing-in-r/
  zigzag_ord <- function(x, n = length(cl)) {
    #x <- quant_pep_IDss
    ord <- data.frame(Original = seq_along(x),
                      Length = vapply(x, length, 1))
    ord <- ord[order(ord$Length, decreasing = TRUE),]
    ord$NewOrd <- rep(c(seq(1, n), seq(n, 1)), length = length(x))
    ord <- ord[order(ord$NewOrd),]
    ord$NewOrd <- 1:nrow(ord)
    ord <- ord[order(ord$Original),]
    return(ord)
  }
  ord <- zigzag_ord(quant_pep_IDss)
  ord$ID <- names(quant_pep_IDss)
  nuOrd <- order(ord$NewOrd)
  quant_pep_IDss <- quant_pep_IDss[nuOrd]
  #
  Pep <- Pep[which(Pep[[pep_IDs]] %in% unlist(quant_pep_IDss)),] # Update Pep (is this necessary?)
  #
  # Get summary method:
  if (misFun("Weights")) {
    if ("weighted.mean" %in% c(LM_fun, reScaling)) {
      warning("A summary method is \"weighted mean\", but no weights were provided...?\nSetting all weights to 1...")
    }
    Weights <- "Weights"
    Pep[[Weights]] <- 1
  }
  #
  # Optional intensity weights - so a peptide's weight is a function of its (non log-transformed) intensity
  if ((misFun(useIntWeights))||(!is.logical(useIntWeights))||(length(useIntWeights) != 1)||(is.na(useIntWeights))) {
    useIntWeights <- FALSE
  }
  if ((useIntWeights)&&(sum(c(reScaling, LM_fun) %in% c("mean", "weighted.mean")))) {
    tmp <- data.table::data.table(pepInt_log^Pep[, Pep.Intens.Nms, drop = FALSE])
    tmp$pepID <- Pep[[pep_IDs]]
    tmp <- data.table::melt(tmp, id.vars = "pepID")
    tmp <- tmp[which(is.all.good(tmp$value, 2)),]
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
    ratGroups$samples <- lapply(ratGroups$values, function(x) {
      expMap[which(expMap[[ratGroups$column]] == x), expMap_Samples_col]
    })
    ratGroups$refSamples <- lapply(ratGroups$values, function(x) {
      expMap[which((expMap[[ratGroups$column]] == x)&(expMap$Reference)), expMap_Samples_col]
    })
    ratGroups$newInt <- lapply(1:length(ratGroups$values), function(x) { #x <- 1
      allIntCol <- paste0(pepInt_Root, ratGroups$samples[[x]])
      allRatCol <- paste0(pepRat_root, ratGroups$samples[[x]])
      w <- which(allRatCol %in% colnames(tmpPep))
      myRatCol <- allRatCol[w]
      myIntCol <- allIntCol[w]
      refIntCol <- paste0(pepInt_Root, ratGroups$refSamples[[x]])
      refInt <- rowMeans(tmpPep[, refCol, drop = FALSE], na.rm = TRUE)
      allInt <- tmpPep[, allIntCol]
      i1 <- rowMeans(allInt, na.rm = TRUE)
      newInt <- tmpPep[, myRatCol]/log(pepInt_log, pepRat_log) + refInt
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
  quant_pep_IDss2 <- proteoCraft::listMelt(quant_pep_IDss, ColNames = c("id", "PG"))
  quant_pep_IDss2$mtch <- match(quant_pep_IDss2$id, tmpPep[[pep_IDs]])
  quant_pep_IDss2 <- quant_pep_IDss2[order(quant_pep_IDss2$mtch, decreasing = FALSE),]
  quant_pep_IDss2 <- data.table::as.data.table(quant_pep_IDss2)
  quant_pep_IDss2 <- quant_pep_IDss2[, list(IDs = list(id)), by = list(PG = PG)]
  quant_pep_IDss2 <- setNames(quant_pep_IDss2$IDs, quant_pep_IDss2$PG)
  quant_pep_IDss[names(quant_pep_IDss2)] <- quant_pep_IDss2[names(quant_pep_IDss2)] # There are some empty entries in quant_pep_IDss: proteins with no peptides eligible for quant
  rm(quant_pep_IDss2)
  #
  # Quantitation:
  #  - Quantitation step
  cat(" - Calculating profile...\n")
  quntNms <- gsub(proteoCraft::topattern(pepInt_Root), Expr.root.full, Pep.Intens.Nms)
  # Below: object used by most methods and by reScaling
  tmpPep2 <- listMelt(quant_pep_IDss, ColNames = c("pepID", "PG"))
  tmpPep2[, Pep.Intens.Nms] <- tmpPep[match(tmpPep2$pepID, tmpPep[[pep_IDs]]), Pep.Intens.Nms]
  #
  allQuants <- list()
  cat("   method =", LFQ_algo, "\n")
  stopCl <- FALSE
  if (sum(c("INHOUSE_MAXLFQ", "LM") %in% c(LFQ_ALGO, RESCALING))) {
    # Create cluster
    tstCl <- stopCl <- misFun(cl)
    if (!misFun(cl)) {
      tstCl <- suppressWarnings(try({
        a <- 1
        parallel::clusterExport(cl, "a", envir = environment())
      }, silent = TRUE))
      tstCl <- !"try-error" %in% class(tstCl)
    }
    if ((misFun(cl))||(!tstCl)) {
      dc <- parallel::detectCores()
      if (misFun(N.reserved)) { N.reserved <- 1 }
      if (misFun(N.clust)) {
        N.clust <- max(c(dc-N.reserved, 1))
      } else {
        if (N.clust > max(c(dc-N.reserved, 1))) {
          warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
          N.clust <- max(c(dc-N.reserved, 1))
        }
      }
      cl <- parallel::makeCluster(N.clust, type = "SOCK")
    }
    N.clust <- length(cl)
  }
  if ("LM" %in% c(LFQ_ALGO, RESCALING)) {
    tmpFl1 <- tempfile()
    tmpFl2 <- tempfile()
    exports <- list("LFQ_algo", "LFQ_ALGO", "tmpFl1", "tmpFl2", "pep_IDs", "Pep.Intens.Nms", "minN", "maxN", "LM_fun",
                    "Weights" #, "LFQ_stab"
    )
    parallel::clusterExport(cl, exports, envir = environment())
    readr::write_rds(quant_pep_IDss, tmpFl1)
    readr::write_rds(tmpPep, tmpFl2)
    parallel::clusterCall(cl, function() {
      quant_pep_IDss <<- readr::read_rds(tmpFl1)
      tmpPep <<- readr::read_rds(tmpFl2)
      return()
    })
    unlink(tmpFl1)
    unlink(tmpFl2)
    #NB: we should add the option to stabilize large ratios (MaxLFQ-like) to LFQ.lm!!!
    makeProf <- function(ids,
                         fun,
                         pep,
                         intCol,
                         weights,
                         min,
                         max) { #ids <- quant_pep_IDss[[1]]
      proteoCraft::LFQ.lm(ids,
                          InputTabl = pep,
                          id = pep_IDs,
                          IntensCol = intCol,
                          Summary.method = fun,
                          Summary.weights = weights,
                          Min.N = min,
                          Max.N = max,
                          reNorm = 1)
    }
    environment(makeProf) <- .GlobalEnv
    res2a <- try(parallel::parLapply(cl, quant_pep_IDss, function(x) {
      makeProf(ids = x,
               fun = LM_fun,
               pep = tmpPep,
               intCol = Pep.Intens.Nms,
               weights = Weights,
               min = minN,
               max = maxN)
    }), silent = TRUE)
    Done <- !("try-error" %in% class(res2a))
    if (!Done) {
      cat(" - re-running quant algorithm, the slow way...\nsomething clearly went wrong with the cluster\n")
      # This function has occasionally and non-reproducibly failed on weak PCs, this is a back up
      res2a <- lapply(quant_pep_IDss, function(x) {
        makeProf(ids = x,
                 fun = LM_fun,
                 pep = tmpPep,
                 intCol = Pep.Intens.Nms,
                 weights = Weights,
                 min = minN,
                 max = maxN)
      })
    }
    res2a <- as.data.frame(do.call(rbind, res2a))
    allQuants$LM <- res2a
    if (LFQ_ALGO == "LM") {
      res2 <- res2a
    }
  }
  # Unfinished in-house implementation of maxLFQ - still has some issues:
  # in particular does not deal with cases where the ratios tree is disconnected into several independent sub-trees
  # (because of missing values).
  if ("INHOUSE_MAXLFQ" %in% c(LFQ_ALGO, RESCALING)) {
    # Export/serialize to cluster
    # exports <- list("LFQ_algo", "LFQ_ALGO", "wd", "pep_IDs", "Pep.Intens.Nms", "minN", "maxN", "LM_fun",
    #                 "Weights", "mySmpls", "LFQ_stab")
    # parallel::clusterExport(cl, exports, envir = environment())
    # readr::write_rds(quant_pep_IDss, paste0(wd, "/tempIDs.RDS"))
    # readr::write_rds(tmpPep, paste0(wd, "/tempPep.RDS"))
    # invisible(clusterCall(cl, function() {
    #   quant_pep_IDss <<- readr::read_rds(paste0(wd, "/tempIDs.RDS"))
    #   tmpPep <<- readr::read_rds(paste0(wd, "/tempPep.RDS"))
    #   if (LFQ_ALGO == "INHOUSE_MAXLFQ") {
    #     L <<- length(Pep.Intens.Nms)
    #     if (L > 1) {
    #       ratDefs <<- combn(L, 2)
    #       N <<- ncol(ratDefs)
    #       Mtr <- matrix(0, N, L)
    #       for (k in 1:N) {
    #         i <- ratDefs[1, k]
    #         j <- ratDefs[2, k]
    #         Mtr[k, i] <- -1
    #         Mtr[k, j] <-  1
    #       }
    #       Mtr_full <<- Mtr#[, -1] # Fix N1 = 0 by removing that column and solving for the rest // nope, do this later!!!
    #     }
    #   }
    #   return()
    # }))
    # if (LFQ_ALGO == "INHOUSE_MAXLFQ") {
    #   makeProf <- function(ids, stabilize = LFQ_stab) {
    #     #ids <- quant_pep_IDss[[1]]
    #     #ids <- quant_pep_IDss[[2]]
    #     #ids <- quant_pep_IDss[[3]]
    #     nIDs <- length(ids)
    #     if (!nIDs) { return(rep(NA, L)) }
    #     if (L == 1) { return(0) } #!!! We must remain able to output a profile when there is but one sample!!!
    #     m <- match(ids, tmpPep[[pep_IDs]])
    #     qnt <- tmpPep[m, Pep.Intens.Nms]
    #     wghts <- tmpPep[m, Weights]
    #     if (stabilize) {
    #       # See MaxLFQ paper (Eq. 5)
    #       wFeat <- apply(qnt, 2, function(x) { which(!is.na(x)) }) # Per sample
    #       nFeat <- vapply(wFeat, length, 1)
    #       wOverLap <- apply(ratDefs, 2, function(x) { #x <- ratDefs[, 1]
    #         w_ij <- unlist(wFeat[x])
    #         if (!length(w_ij)) { return(0) }
    #         w_ij <- aggregate(w_ij, list(w_ij), length)
    #         return(w_ij$Group.1[which(w_ij$x == 2)])
    #       })
    #       stabTst <- vapply(1:N, function(x) { #x <- 1
    #         nMx <- max(nFeat[ratDefs[, x]]) # Number of features in the sample with the most features
    #         nShrd <- length(wOverLap[[x]]) # Number of shared features
    #         return(nMx/nShrd)
    #       }, 1)
    #     }
    #     # Calculate log10 ratios
    #     pepRats <- lapply(1:nIDs, function(i) { #i <- 1
    #       int <- unlist(qnt[i, ])
    #       ratLst <- lapply(2:L, function(x) {
    #         rs <- int[x]-int[1:(x-1)] 
    #         #if (x < L) { rs <- c(rs, rep(NA, L-x)) }
    #         return(rs)
    #       })
    #       #ratMatr <- do.call(rbind, ratLst)
    #       #colnames(ratMatr) <- mySmpls[1:(L-1)]
    #       #rownames(ratMatr) <- mySmpls[2:L]
    #       #return(ratMatr)
    #       return(setNames(unlist(ratLst), NULL))
    #     })
    #     pepRats <- do.call(rbind, pepRats)
    #     #rowSums(apply(pepRats, 1, is.na))
    #     #tst <- pepRats; View(tst)
    #     prof_fun <- get(LM_fun)
    #     if (LM_fun == "weighted.mean") {
    #       avg_pepRats <- apply(pepRats, 2, prof_fun, wghts, na.rm = TRUE)
    #     } else {
    #       avg_pepRats <- apply(pepRats, 2, prof_fun, na.rm = TRUE)
    #     }
    #     #avg_pepRats
    #     if (stabilize) {
    #       w1 <- which(stabTst > 5)
    #       w2 <- which((stabTst >= 2.5)&(stabTst <= 5))
    #       l1 <- length(w1)
    #       l2 <- length(w2)
    #       if (sum(l1, l2)) {
    #         sum_pepRats <- rep(NA, N)
    #         sum_pepRats[c(w1, w2)] <- vapply(c(w1, w2), function(x) { #x <- c(w1, w2)[1]
    #           j <- ratDefs[, x]
    #           #i <- wOverLap[[x]] # no: the way I read it, it doesn't filter for overlap
    #           #x <- qnt[i, j]
    #           x <- qnt[, j]
    #           x <- 10^x
    #           x <- colSums(x, na.rm = TRUE)
    #           x <- log10(x)
    #           return(x[1]-x[2])
    #         }, 1)
    #         if (l1) {
    #           #View(cbind(avg_pepRats[w1], sum_pepRats[w1]))
    #           wOK1 <- w1[which(!is.na(sum_pepRats[w1]))]
    #           if (length(wOK1)) {
    #             avg_pepRats[wOK1] <- sum_pepRats[wOK1]
    #           }
    #         }
    #         if (l2) {
    #           wOK2 <- w2[which(!is.na(sum_pepRats[w2]))]
    #           if (length(wOK2)) {
    #             w_of_x <- (stabTst[wOK2]-2.5)/2.5
    #             avg_pepRats[wOK2] <- w_of_x*sum_pepRats[wOK2] + (1-w_of_x)*avg_pepRats[wOK2]
    #           }
    #         }
    #       }
    #     }
    #     # Usable ratios
    #     wY <- which(is.finite(avg_pepRats))
    #     #
    #     lY <- length(wY)
    #     A <- matrix(0,
    #                 nrow = lY,
    #                 ncol = L)
    #     b <- numeric(lY)
    #     for (k in 1:lY) { #k <- 1
    #       i <- ratDefs[1, wY[k]]
    #       j <- ratDefs[2, wY[k]]
    #       A[k, i] <- 1
    #       A[k, j] <- -1
    #       b[k] <- avg_pepRats[k]
    #     }
    #     g <- igraph::graph_from_edgelist(t(ratDefs[, wY]), directed = FALSE)
    #
    #     # Here we still need a way to deal with the each independent subtree,
    #     # where missing pairwise ratios result in more than one disconnected subtrees.
    #     # Or we could just skip those - since even outputting a value here is already making an assumption.
    #
    #     igraph::components(g)$no
    #     igraph::is_connected(g)
    #
    #     u <- sort(unique(as.integer(ratDefs[, wY])))
    #     mtr_red <- Mtr_full[wY, u]
    #     mtr_red <- mtr_red[, -1]
    #     #
    #     # 3 methods to solve rank-deficient system
    #     # - QR
    #     qr_m <- qr(mtr_red)
    #     #qr_m$rank
    #     #ncol(mtr_red)
    #     I1 <- qr.coef(qr_m, avg_pepRats[wY])
    #     I1[which(is.na(I1))] <- 0 # Not good: gives us exactly equal values!
    #     # - lm.fit (QR under hood)
    #     fit <- lm.fit(mtr_red, avg_pepRats[wY])
    #     I2 <- fit$coefficients
    #     I2[which(is.na(I2))] <- 0 # Not good: gives us exactly equal values!
    #     # Approximation
    #     library(MASS)
    #     I3 <- MASS::ginv(mtr_red) %*% avg_pepRats[wY]
    #     View(data.frame(QR = I1,
    #                     lm.fit = I2,
    #                     ginv = I3))
    #     #
    #     
    #     I <- qr.solve(mtr_red, avg_pepRats[wY])
    #     I <- c(0, I)  # add N1 = 0 back
    #     I <- I-median(I)
    #     return(I)
    #   }
    # }
    environment(makeProf) <- .GlobalEnv
    res2f <- try(parallel::parLapply(cl, quant_pep_IDss, makeProf), silent = TRUE)
    #
    if ("try-error" %in% class(res2f)) {
      cat(" - re-running quant algorithm, the slow way...\nsomething clearly went wrong with the cluster\n")
      # This function has occasionally and non-reproducibly failed on weak PCs, this is a back up
      res2f <- lapply(quant_pep_IDss, makeProf)
    }
    unlink(paste0(wd, "/tempIDs.RDS"))
    unlink(paste0(wd, "/tempPep.RDS"))
    res2f <- as.data.frame(do.call(rbind, res2f))
    allQuants$INHOUSE_MAXLFQ <- res2f
    if (LFQ_ALGO == "INHOUSE_MAXLFQ") {
      res2 <- res2f
    }
  }
  if (sum(c(LFQ_ALGO, RESCALING) %in% c("IQ", "LIMPA", "QFEATURES" #, "MSSTATS"
                      ))) {
    # These algorithms do not take too long to run,
    # thus they could be run either be run as LFQ, with or without subsequent re-scaling,
    # or used to provide re-scaling for another LFQ method.
    if ("IQ" %in% c(LFQ_ALGO, RESCALING)) {
      tmp4 <- tmpPep2
      colnames(tmp4) <- proteoCraft::cleanNms(gsub(".* - ", "", colnames(tmp4)))
      tmp4 <- melt(tmp4, id.vars = c("pepID", "PG"))
      tmp4 <- tmp4[which(!is.na(tmp4$value)),]
      tmp4 <- list(protein_list = tmp4$PG,
                   sample_list = tmp4$variable,
                   pepID = tmp4$pepID,
                   quant = tmp4$value)
      res2b <- iq::fast_MaxLFQ(tmp4)
      res2b <- as.data.frame(res2b$estimate)
      res2b <- res2b[match(names(quant_pep_IDss), row.names(res2b)),]
      res2 <- res2b
      allQuants$IQ <- res2b
      if (LFQ_ALGO == "IQ") {
        res2 <- res2b
      }
    }
    if ("LIMPA" %in% c(LFQ_ALGO, RESCALING)) {
      tmp4 <- as.matrix(tmpPep2[, Pep.Intens.Nms])
      dpcfit <- limpa::dpcCN(tmp4)
      #dpcfit <- limpa::dpc(tmp4)
      #limpa::plotDPC(dpcfit)
      dpcRes <- limpa::dpcQuant(tmp4, tmpPep2$PG, dpc = dpcfit)
      res2c <- dpcRes$E
      #
      # Code to replace every value based on 0 observations with NA:
      # w <- which(dpcRes$other$n.observations == 0)
      # res2c[w] <- NA
      # However, Gordon Smyth does not recommend doing it.
      # For a discussion about the proper usage of limpa, why there are no missing values in its output quant matrix,
      # or why it may output repeated values in a row,
      # please see:
      # https://github.com/SmythLab/limpa/issues/5#event-22032865875
      #
      res2c <- as.data.frame(res2c)
      colnames(res2c) <- mySmpls
      res2c <- res2c[match(names(quant_pep_IDss), row.names(res2c)),]
      allQuants$LIMPA <- res2c
      if (LFQ_ALGO == "LIMPA") {
        res2 <- res2c
      }
    }
    if ("QFEATURES" %in% c(LFQ_ALGO, RESCALING)) {
      tmp4 <- tmpPep2
      tmp4$pepID <- 1:nrow(tmp4) # Otherwise aggregateFeatures throws an error
      QFeat_obj <- QFeatures::readQFeatures(assayData = tmp4,
                                            fnames = "pepID",
                                            quantCols = Pep.Intens.Nms,
                                            name = "peptides")
      QFeat_obj <- QFeatures::aggregateFeatures(QFeat_obj,
                                                i = "peptides",
                                                fcol = "PG",
                                                na.rm = TRUE,
                                                name = "PG")
      res2d <- as.data.frame(SummarizedExperiment::assay(QFeat_obj[["PG"]]))
      res2d <- res2d[match(names(quant_pep_IDss), row.names(res2d)),]
      allQuants$QFEATURES <- res2d
      if (LFQ_ALGO == "QFEATURES") {
        res2 <- res2d
      }
    }
    # if ("MSSTATS" %in% c(LFQ_ALGO, RESCALING)) {
    #   tmp4 <- tmpPep2[, c("pepID", "PG")]
    #   tmp4[, gsub(".* - ", "", Pep.Intens.Nms)] <- origInt[match(tmp4$pepID, origInt[[pep_IDs]]), Pep.Intens.Nms]
    #   #
    #   tmp4 <- melt(tmp4, id.vars = c("pepID", "PG"))
    #   colnames(tmp4)[2:4] <- c("ProteinName", "Run", "Intensity")
    #   tmp4[, c("PeptideModifiedSequence", "PeptideSequence")] <- Pep[match(tmp4$pepID, Pep[[pep_IDs]]),
    #                                                                  c("Modified sequence", "Sequence")]
    #   tmp4$PeptideModifiedSequence <- gsub("^_|_$", "", tmp4$PeptideModifiedSequence)
    #   tmp4$IsotopeLabelType <- "light"
    #   tmp4$ProductCharge <- tmp4$FragmentIon <- tmp4$PrecursorCharge <- NA
    #   tmp4$Fraction <- 1
    #   u <- as.character(unique(tmp4$Run))
    #   mu <- match(tmp4$Run, u)
    #   grp <- expMap[match(u, expMap[[expMap_Samples_col]]), smplGroups$column]
    #   rep <- expMap$Replicate[match(u, expMap[[expMap_Samples_col]])]
    #   tmp4$Condition <- grp[mu]
    #   tmp4$BioReplicate <- rep[mu]
    #   #
    #   # Make MSstats contrast matrix
    #   grps1 <- unique(expMap[which(!expMap$Reference), smplGroups$column])
    #   grps0 <- unique(expMap[which(expMap$Reference), smplGroups$column])
    #   grps <- c(grps1, grps0)
    #   contr <- do.call(rbind, lapply(ratGroups$values, function(x) { #x <- ratGroups$values[1]
    #     em <- expMap[which(expMap[[ratGroups$column]] == x),]
    #     if ((!nrow(em))||(length(unique(em$Reference)) < 2)) { return() }
    #     grp1 <- unique(em[which(!em$Reference), smplGroups$column])
    #     grp0 <- unique(em[which(em$Reference), smplGroups$column])
    #     m1 <- match(grp1, grps)
    #     m0 <- match(grp0, grps)
    #     return(do.call(rbind, lapply(grp0, function(y) { #y <- grp0[1]
    #       mtr <- matrix(rep(0, length(grps)*length(grp1)), ncol = length(grps))
    #       mtr[, match(y, grps)] <- -1
    #       for (i in 1:length(grp1)) {
    #         mtr[i, match(grp1[i], grps)] <- 1
    #       }
    #       colnames(mtr) <- grps
    #       rownames(mtr) <- paste0(proteoCraft::cleanNms(grp1, rep = "_"), "_vs_", proteoCraft::cleanNms(y, rep = "_"))
    #       return(mtr)
    #     })))
    #   }))
    #   # Slowww...
    #   MSstats_summarized <- MSstats::dataProcess(tmp4,
    #                                              logTrans = pepInt_log,
    #                                              normalization = FALSE, # The data is already pre-normalized
    #                                              numberOfCores = N.clust,
    #                                              featureSubset = "all", # IMPORTANT: we already filtered peptides
    #                                              summaryMethod = "TMP",
    #                                              censoredInt = "NA")
    #   # setequal(colnames(contr),
    #   #          unique(as.character(MSstats_summarized$ProteinLevelData$GROUP)))
    #   #
    #   #usedPep <- unique(gsub("_NA_NA_NA$", "", MSstats_summarized$FeatureLevelData$FEATURE))
    #   #origPep <- unique(tmp4$PeptideModifiedSequence)
    #   #length(setdiff(origPep, usedPep)) # Number of peptides not used by MSstats
    #   # Important caveat here: this doesn't already contain the logFCs which MSstats would use.
    #   # So we cannot use MSstats_summarized$ProteinLevelData, we must proceed to the next step:
    #   # Slowww too...
    #   MSstats_model <- MSstats::groupComparison(contr,
    #                                             MSstats_summarized,
    #                                             log_base = pepInt_log,
    #                                             numberOfCores = N.clust)
    #   
    #   table(summarized$FeatureLevelData$GROUP,
    #         summarized$FeatureLevelData$RUN)
    #   
    #   
    #   #View(MSstats_model$ComparisonResult)
    #   res3e <- MSstats_model$ComparisonResult[, c("Protein", "Label", paste0("log", pepInt_log, "FC"))]
    #   colnames(res3e)[3] <- "logFC"
    #   res3e <- cast(res3e, Protein~Label, mean, value = "logFC")
    #   res3e <- as.data.frame(res3e)
    #   rg <- 2:ncol(res3e)
    #   colnames(res3e)[rg] <- gsub("_vs_.*", "", colnames(res3e)[rg])
    #   colnames(res3e)[rg] <- expMap[match(colnames(res3e)[rg],
    #                                                proteoCraft::cleanNms(expMap[[smplGroups$column]], rep = "_")),
    #                                          smplGroups$column]
    #   res3e <- res3e[match(names(quant_pep_IDss), res3e$Protein), rg, drop = FALSE]
    # }
  }
  colnames(res2) <- quntNms
  rownames(res2) <- names(quant_pep_IDss)
  #
  #########################################
  # Code to compare the different methods #
  #########################################
  # colnames(res2a) <- colnames(res2b) <- colnames(res2c) <- colnames(res2d) <- quntNms
  # nr <- nrow(res2a)
  # tst <- lapply(quntNms, function(k) {
  #   i <- proteoCraft::cleanNms(gsub(".* - ", "", k))
  #   data.frame(Sample = i,
  #              X = c(rep(res2a[[k]], 2),
  #                    rep(res2b[[k]], 3),
  #                    res2c[[k]]),
  #              X_method = c(rep("LM", 3*nr),
  #                           rep("iq", 2*nr),
  #                           rep("limpa", nr)),
  #              Y = c(res2b[[k]], res2c[[k]], res2d[[k]],
  #                    res2c[[k]], res2d[[k]],
  #                    res2d[[k]]),
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
  # proteoCraft::poplot(plot, 12, 22)
  # #
  # em <- expMap
  # em$Replicate <- as.integer(em$Replicate)
  # em <- em[order(em$Replicate),]
  # f0 <- function(i, y) { #i <- 1
  #   w1 <- which(contr[i,] == 1)
  #   w0 <- which(contr[i,] == -1)
  #   w1 <- which(em[[smplGroups$column]] == colnames(contr)[w1])
  #   w0 <- which(em[[smplGroups$column]] == colnames(contr)[w0])
  #   if (Nested) {
  #     rep <- unique(em$Replicate[c(w1, w0)])
  #     res <- lapply(rep, function(r) { #r <- 1
  #       w1_ <- w1[which(em$Replicate[w1] == r)]
  #       w0_ <- w0[which(em$Replicate[w0] == r)]
  #       if ((length(w1_) != 1)||(length(w0_) != 1)) { return() }
  #       smpl1_ <- em[w1_, expMap_Samples_col]
  #       smpl0_ <- em[w0_, expMap_Samples_col]
  #       return(y[[paste0("log10(Expr.) - ", smpl1_)]] - y[[paste0("log10(Expr.) - ", smpl0_)]])
  #     })
  #     res <- do.call(cbind, res)
  #   } else {
  #     smpls1 <- em[w1, expMap_Samples_col]
  #     smpls0 <- em[w0, expMap_Samples_col]
  #     rf <- rowMeans(y[, paste0("log10(Expr.) - ", smpls0)], na.rm = TRUE)
  #     res <- sweep(y[, paste0("log10(Expr.) - ", smpls1)], 1, rf, "-")
  #   }
  #   res <- rowMeans(res, na.rm = TRUE)
  #   #res <- res/log(pepRat_log, pepInt_log) # res3e is still pepInt_log base
  #   return(res)
  # }
  # res3a <- as.data.frame(do.call(cbind, setNames(lapply(1:nrow(contr), function(i) { f0(i, y = res2a) }), colnames(res3e))))
  # res3b <- as.data.frame(do.call(cbind, setNames(lapply(1:nrow(contr), function(i) { f0(i, y = res2b) }), colnames(res3e))))
  # res3c <- as.data.frame(do.call(cbind, setNames(lapply(1:nrow(contr), function(i) { f0(i, y = res2c) }), colnames(res3e))))
  # res3d <- as.data.frame(do.call(cbind, setNames(lapply(1:nrow(contr), function(i) { f0(i, y = res2d) }), colnames(res3e))))
  # tst <- lapply(colnames(res3e), function(k) {
  #   data.frame(Group = proteoCraft::cleanNms(k),
  #              X = c(rep(res3a[[k]], 4),
  #                    rep(res3b[[k]], 3),
  #                    rep(res3c[[k]], 2),
  #                    res3d[[k]]),
  #              X_method = c(rep("LM", 4*nr),
  #                           rep("iq", 3*nr),
  #                           rep("limpa", 2*nr),
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
  # proteoCraft::poplot(plot, 12, 22)
  #######
  # End #
  #######
  #
  # Re-scaling step
  if (!skip_reScaling) {
    cat("Rescaling expression values...\n")
    if (RESCALING %in% c("MAXLFQ", "IQ", "LIMPA", "QFEATURES", "LM")) {
      if (RESCALING == "MAXLFQ") {
        rescVal <- data.frame(PG = Prot$temp_IDs,
                              Value = log(Prot$"Summed Intensities", pepInt_log))
      } else {
        rescVal <- data.frame(PG = rownames(allQuants[[RESCALING]]),
                              Value = rowMeans(allQuants[[RESCALING]][, quntNms, drop = FALSE], na.rm = TRUE))
      }
    } else {
      rescVal <- tmpPep2[, c("pepID", "PG")]
      m <- match(rescVal$pepID, Pep[[pep_IDs]])
      rescVal$avgPepInt <- Pep$avgPepInt[m]
      if (reScaling == "weighted.mean") {
        rescVal$Weights <- Pep$Weights[m]
      }
      if ((reSc_is_topN)||(is.finite(reSc_topN))) {
        rescVal$Rank <- stats::ave(rescVal$PG, rescVal$PG, FUN = seq_along) # (thanks chatGPT...)
      }
      if (is.finite(reSc_topN)) {
        rescVal <- rescVal[which(rescVal$Rank <= reSc_topN),]
      }
      rescVal <- data.table::as.data.table(rescVal)
      if (reSc_is_topN&&topN_correct) {
        # Not sure about this...
        Md <- median(rescVal$avgPepInt, na.rm = TRUE)
        tst <- as.data.frame(rescVal[, .(Median = median(avgPepInt, na.rm = TRUE)), by = .(Rank = Rank)])
        for (i in tst$Rank) {
          wi <- which(rescVal$Rank == i)
          rescVal$avgPepInt[wi] <- rescVal$avgPepInt[wi]-tst$Median[match(i, tst$Rank)]
        }
        rescVal$avgPepInt <- rescVal$avgPepInt + Md
      }
      if (reScaling == "weighted.mean") {
        rescVal <- rescVal[, .(Value = weighted.mean(avgPepInt, Weights, na.rm = TRUE)), by = .(PG = PG)]
      } else {
        rescVal <- rescVal[, .(Value = reSc_fun(avgPepInt, na.rm = TRUE)), by = .(PG = PG)] 
      }
      rescVal <- as.data.frame(rescVal)
    }
    rescVal <- rescVal[match(row.names(res2), rescVal$PG),]
    # Re-scale
    currVal <- rowMeans(res2[, quntNms, drop = FALSE], na.rm = TRUE)
    res2[, quntNms] <- sweep(res2[, quntNms, drop = FALSE], 1, rescVal$Value-currVal, "+")
  }
  # Homogenize ratios and expression values, and create reference ratios
  if (!skipRatios) {
    tmpRt <- proteoCraft::make_Rat(res2,
                                   Priority = Priority,
                                   pepInt_log,
                                   pepRat_log,
                                   refGroups = refGroups,
                                   experiment.map = expMap,
                                   int.root = Expr.root.full,
                                   rat.root = pepRat_root)
    kol <- colnames(tmpRt$Data)
    res2[, kol] <- tmpRt$Data[, kol]
    # Calculate Reference-to-Reference ratio vectors:
    # Mirrored from the code after "# Create peptide-level Ref-to-Ref ratios (useful for PTMs analysis):" in the main replicates-workflow
    # The main difference is that the input data is log-transformed (not for peptides).
    if (!misFun(param)) {
      rat_cont_grps <- param$Ratios.Contaminant.Groups
    } else { rat_cont_grps <- "Ratio groups" } # Default
    res3 <- try(proteoCraft::make_RefRat(data = res2,
                                         experiment.map = expMap,
                                         int.root = Expr.root.full,
                                         rat.root = pepRat_root,
                                         rat.con.grps = rat_cont_grps,
                                         mode = ref_Mode,
                                         parameters = param,
                                         logInt = pepInt_log,
                                         logRat = pepRat_log), silent = TRUE)
    if (!"try-error" %in% class(res3)) {
      res2[, colnames(res3)] <- res3
    } else { warning("No Ref to Ref columns were generated!") }
  }
  # Finally: this should be the very last step to avoid confusions:
  # If output should not be log-transformed, de-log!
  if (!protLFQ_toLog) { # De-log LFQ values
    w <- grep(proteoCraft::topattern(Expr.root.full), colnames(res2))
    g <- colnames(res2)[w]
    res2[, g] <- pepInt_log^(res2[, g])
    colnames(res2)[w] <- gsub(proteoCraft::topattern(Expr.root.full),
                              "Expr. - ",
                              colnames(res2)[w])
  }
  if (!skipRatios) {
    w <- grep(proteoCraft::topattern(pepRat_root), colnames(res2))
    g <- colnames(res2)[w]
    if (!protRat_toLog) { # If necessary de-log
      res2[, g] <- pepRat_log^(res2[, g])
      colnames(res2)[w] <- gsub(proteoCraft::topattern(pepRat_root),
                                "Ratio - ",
                                colnames(res2)[w])
    } else {
      # If changing log base from peptides to proteins:
      if (protRat_toLog != pepRat_log) {
        res2[, g] <- res2[, g]/log(protRat_toLog, pepRat_log)
        colnames(res2)[w] <- gsub(proteoCraft::topattern(pepRat_root),
                                  paste0("log", protRat_toLog, "(Ratio) - "), colnames(res2)[w])
      }
    }
  }
  for (i in 1:ncol(res2)) {
    if (!"numeric" %in% class(res2[[i]])) {
      stop(paste0("I would expect the class of column ", i, " to be numeric! Investigate!"))
      #res2[[i]] <- as.numeric(res2[[i]])
    }
  }
  res2$"Peptides IDs used for quantitation" <- vapply(quant_pep_IDss, paste, "", collapse = ";")
  ord2 <- ord[nuOrd,]
  res2 <- res2[order(ord2$Original),]
  rownames(res2) <- rownames(Prot)
  #
  if (stopCl) { parallel::stopCluster(cl) }
  RES <- list(Data = res2)
  if ("LIMPA" %in% c(LFQ_ALGO, RESCALING)) {
    RES$EList_obj <- dpcRes
  }
  if ("QFEATURES" %in% c(LFQ_ALGO, RESCALING)) {
    RES$QFeatures_obj <- QFeat_obj
  }
  return(RES)
}
