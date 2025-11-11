#' Prot.Quant
#'
#' @description 
#' A function to calculate estimated protein group Expression (Absolute Quant) and Ratios (optional) from individual peptide values.
#' The input is assumed to be normalized.
#' 
#' A note on ref.groups and ratio.groups:
#' These parameters are the method I have used to define if replicates are paired or not.
#' For an experiment with one master Experiment, several Conditions (incl. one reference) and Replicates:
#' - in a paired (= "nested") setup, you would set ref.groups to "ExpRep" and ratio.groups to "Exp"
#' - in an unpaired setup, you would set ref.groups to "Exp" to and ratio.groups to "Exp"
#' Essentially, ref.groups are the groups within which ratios are calculated to all available references,
#' while ratio.groups are groups within which reference-to-reference values are calculated:
#' - between individual references in the paired setup
#' - from individual references to the average reference in the unpaired setup
#' 
#' @param Prot Protein/Protein groups table. A data.frame.
#' @param Peptide.IDs Name of the Protein/Protein groups table's peptide IDs column. Default = "Peptide IDs"
#' @param Pep Peptides table. A data.frame.
#' @param Mode In "Classic" mode (current default), the function just uses the "Peptide.IDs" argument. In "PreferUnique" mode, the "Peptide.IDs" column should contain razor or shared peptide IDs, and a second column name of unique peptides should be provided using the "Unique.peptide.IDs" argument. If at least "Min.Unique" unique peptides are present, then only those will be used for quantitation, otherwise as many razor/shared peptides as necessary will be added (sorted by decreasing average intensities).
#' @param id The name of the Peptides table's IDs column. Default = "id"
#' @param Summary.method The summary method used for ratios (the Levenberg-Marquardt algorithm is used for intensities). One of "mean", "median", or "weighted.mean".
#' @param Summary.weights If a "weighted.mean" summary method is chosen, then a vector of weights must be provided. This should be the name of a column of the peptides table
#' @param Intensity.weights Logical, default = FALSE; if TRUE, will take into account individual peptide intensities when calculating average profile if "Summary.method" is either mean "mean" or "weighted.mean" (thus, it will actually be a weighted mean regardless). 
#' @param Priority One of "Ratios" or "Intensities" (default). Some flexibility in spelling is allowed. You want to prioritize ratios for SILAC because in this case MaxQuant measures peptides ratios more precisely than intensities. Otherwise, you want to prioritize Intensities and re-calculate ratios from them.
#' @param Skip.ratios Logical, default = FALSE. If TRUE, ratios will not be calculated. 
#' @param experiments.map Map of the experiment, default = Exp.map
#' @param param The experiment's parameters object. If provided, the ref.groups argument is not required.
#' @param aggregate.map The aggregate map. Default = Aggregate.map
#' @param aggregate.list The named list of aggregates. Default = Aggregate.list
#' @param aggregates The aggregates themselves. Default = Aggregates
#' @param ref.groups Defines which samples are paired to which references. May alternatively (preferred solution) be provided indirectly through the param argument.
#' @param ratio.groups Defines groups within which ratios are calculated. May alternatively (preferred solution) be provided indirectly through the param argument.
#' @param sample.groups Defines groups of samples which are replicates of the same category. Used to calculate within group reference ratios.
#' @param Pep.Intens.root Root of the peptides intensity column(s) names
#' @param Pep.Ratios.root Root of the peptides ratios column(s) names
#' @param log.Pep.Intens Set to 0 or FALSE if input peptide intensities are linear scale. If the data is already log scale, set to the relevant scale's base. Default = FALSE
#' @param log.Pep.Ratios Set to 0 or FALSE if input peptide ratios are linear scale. If the data is already log scale, set to the relevant scale's base. Default = 2
#' @param Prot.LFQ.to.log Should the output protein LFQ values be log-scale or not? Can ve set to the log base desired. If set to TRUE, this will be the same base as the input intensities log scale, or 10 by default.
#' @param Prot.Ratios.to.log Should the output protein ratios be log-scale or not? Can ve set to the log base desired. If set to TRUE, this will be the same base as the input ratios log scale, or its default, 2
#' @param Mods Which modifications (2 lowercase letters PTM code) should be included? If set to FALSE, will not filter any modifications.
#' @param Mods.to.Exclude Alternative/complementary to Mods. Which modifications should be excluded? (use argument "Discard.unmod" to discard unmodified counterpart peptides.) A data.frame with columns "Mark" (code) and "Where", see argument "Discard.unmod".
#' @param Mod.Nms Default = "Modified sequence". The name of the column containing the modified sequence in the peptides table. Can be set to a non-modified sequence if "Mods" = FALSE and "Mods.to.Exclude" is empty.
#' @param Discard.unmod Default = TRUE. Should we discard those unmodified peptides whose primary sequence is the same as that of some modified peptides we will not use? If set to 2, will use the "Where" column in "Mods.to.Exclude" to identify (and exclude) peptides which could be modified even if the modified form was not identified. Requires knowledge of protein sequence ("Prot.Seq" argument)! Be careful! This will likely massively reduce the number of peptides available for quantitation!
#' @param Prot.Seq Default = "Sequence (1st accession)", used if "Discard.unmod" is set to 2
#' @param Min.N How many peptides should at least be present? Should be at the very least 1 (default = 2).
#' @param Max.N How many peptides can we use at most for the Levenberg-Marquardt procedure (used only for LFQs)? Using too many peptides can be an issue, e.g. with huge proteins like Titin. Default = 50. The most intense peptides will be selected.
#' @param Ratios.SD Deprecated (Should the output ratios also include peptide level standard deviations? Default = FALSE)
#' @param Ratios.Pvalue Deprecated (Should the Output also include peptide level P-values? Default = FALSE)
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param Unique.peptide.IDs Used only if mode = "PreferUnique". Name of the Protein/Protein groups table's unique peptide IDs column. Default = "Unique peptide IDs"
#' @param Min.Unique Used only if mode = "PreferUnique". Minimum number of unique peptides to consider before adding shared peptides. Default = 3
#' @param Refs_Mode How are reference ratios calculated?\cr
#'  - If set to "1", only references are considered (i.e. it compares individual references either to each other, or if available to the average reference for the group).\cr
#'  - If set to "2" (default), for each ratios group, reference ratios are based on comparing every possible pair of samples within the group.\cr
#' @param Refs_AllGroups Deprecated, use Refs_Mode instead.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' 
#' @examples
#' temp <- Prot.Quant(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep, id = "New Peptide ID",
#'                    experiments.map = Exp.map,
#'                    ref.groups = Ratios.Ref.Groups,
#'                    Pep.Intens.root = pep.ref, Pep.Ratios.root = pep.ratios.ref,
#'                    log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
#'                    Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
#'                    Mods = Mod4Quant, Mod.Nms = "Modified sequence",
#'                    Min.N = 2, Ratios.SD = FALSE, Ratios.Pvalue = FALSE)
#' 
#' @export

Prot.Quant <- function(Prot,
                       Peptide.IDs = "Peptide IDs",
                       Mode = "Classic",
                       Pep,
                       id = "id",
                       Summary.method = "median",
                       Summary.weights,
                       Intensity.weights = FALSE,
                       Priority = "Intensities",
                       Skip.ratios = FALSE,
                       experiments.map,
                       ref.groups,
                       ratio.groups,
                       sample.groups,
                       param,
                       aggregate.map = Aggregate.map,
                       aggregate.list = Aggregate.list,
                       aggregates = Aggregates,
                       Pep.Intens.root,
                       Pep.Ratios.root,
                       log.Pep.Intens = FALSE,
                       log.Pep.Ratios = 2,
                       Prot.LFQ.to.log = TRUE,
                       Prot.Ratios.to.log = TRUE,
                       Mods,
                       Mods.to.Exclude,
                       Mod.Nms = "Modified sequence",
                       Discard.unmod = TRUE,
                       Prot.Seq = "Sequence (1st accession)",
                       Min.N = 2,
                       Max.N = 50,
                       Ratios.SD = FALSE,
                       Ratios.Pvalue = FALSE,
                       N.clust,
                       N.reserved = 1,
                       Unique.peptide.IDs = "Unique peptide IDs",
                       Min.Unique = 3,
                       Refs_Mode = "2",
                       Refs_AllGroups = NULL,
                       cl) {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::Prot.Quant);TESTING <- TRUE; cl <- parClust
  #Prot = PG; Mode = "PreferUnique";Peptide.IDs = "Razor peptide IDs"; Unique.peptide.IDs = "Unique peptide IDs";Pep = pep; id = "id";Summary.method = "mean";Intensity.weights = FALSE;experiments.map = Exp.map; param = Param;Pep.Intens.root = pep.ref[length(pep.ref)];Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)];log.Pep.Intens = FALSE; log.Pep.Ratios = 2;Prot.LFQ.to.log = TRUE; Prot.Ratios.to.log = TRUE;Mods = Mod4Quant; Mods.to.Exclude = Mod2Xclud; Discard.unmod = Discard.unmod;Min.N = 1; N.clust = N.clust;Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1]
  #Prot = PG; Peptide.IDs = Pep4Quant; Mode = "PreferUnique"; Pep = pep; Summary.method = "mean"; Intensity.weights = FALSE; Skip.ratios = !MakeRatios; experiments.map = Exp.map; ref.groups = RefGrp; ratio.groups = RatGrp; sample.groups = SmplGrp; Pep.Intens.root = paste0(int.col, " - "); Pep.Ratios.root = paste0(rat.cols["Original"], " - "); log.Pep.Intens = FALSE; log.Pep.Ratios = 2; Prot.LFQ.to.log = TRUE; Prot.Ratios.to.log = TRUE; Mods = Mod4Quant; Mods.to.Exclude = Mod2Xclud; Min.N = 1; N.clust = N.clust; Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1]
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  if ((exists("Refs_AllGroups"))&&(!is.null(Refs_AllGroups))) {
    warning("Argument \"Refs_AllGroups\" is deprecated, use argument Refs_Mode instead.")
  }
  if (misFun(Refs_Mode)) { Refs_Mode <- 2 }
  if (!as.numeric(Refs_Mode) %in% 1:2) { Refs_Mode <- 2 }
  Refs_Mode <- as.character(Refs_Mode)
  #
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
  #
  # Adapted from https://www.r-bloggers.com/2020/12/going-parallel-understanding-load-balancing-in-r/
  zigzag_ord <- function(x, n = length(cl)) {
    #x <- temp.ids
    ord <- data.frame(Original = seq_along(x),
                      Length = sapply(x, length))
    ord <- ord[order(ord$Length, decreasing = TRUE),]
    ord$NewOrd <- rep(c(seq(1, n), seq(n, 1)), length = length(x))
    ord <- ord[order(ord$NewOrd),]
    ord$NewOrd <- 1:nrow(ord)
    ord <- ord[order(ord$Original),]
    return(ord)
  }
  #
  Mode <- toupper(Mode)
  if (!Mode %in% c("CLASSIC", "PREFERUNIQUE")) {
    warning("I could not understand the \"Mode\" argument, defaulting to classic behavior.")
    Mode <- "CLASSIC"
  }
  if (!is.integer(suppressWarnings(as.integer(Min.N)))) { Min.N = 2 }
  if (!is.integer(suppressWarnings(as.integer(Max.N)))) { Max.N = 50 }
  Min.N <- as.integer(Min.N)
  Max.N <- as.integer(Max.N)
  stopifnot(!is.null(Peptide.IDs),
            Peptide.IDs %in% colnames(Prot),
            is.integer(Min.N),
            is.integer(Max.N))
  if (Mode == "PREFERUNIQUE") {
    stopifnot(!is.null(Unique.peptide.IDs),
              Unique.peptide.IDs %in% colnames(Prot))
    if (Peptide.IDs == Unique.peptide.IDs) {
      Mode <- "CLASSIC"
    }
  }
  if (Mode == "PREFERUNIQUE") {
    # Check that all peptide IDs in Prot[[Unique.peptide.IDs]] are also in Prot[[Peptide.IDs]]
    tmp <- Prot[, c(Peptide.IDs, Unique.peptide.IDs)]
    tmp[[Peptide.IDs]] <- strsplit(tmp[[Peptide.IDs]], ";")
    tmp[[Unique.peptide.IDs]] <- strsplit(tmp[[Unique.peptide.IDs]], ";")
    Prot[[Peptide.IDs]] <- apply(tmp, 1 , function(x) { paste(unique(unlist(x)), collapse = ";") })
    if ((Min.Unique < 1)||(!is.integer(as.integer(Min.Unique)))) {
      warning("Bad value of \"Min.Unique\" provided, defaulting to 1!")
      Min.Unique <- 1
    }
  }
  Discard.unmod.strict <- FALSE
  if (Discard.unmod == 2) {
    Discard.unmod.strict <- TRUE
    Discard.unmod <- TRUE
  }
  if ("Reference" %in% colnames(experiments.map)) {
    experiments.map$Reference <- as.logical(toupper(experiments.map$Reference))
  } else {
    if (!Skip.ratios) {
      warning("No reference column was provided, skipping ratios calculation!")
      Skip.ratios <- TRUE
    }
  }
  if ((Min.N < 1)||(!is.integer(as.integer(Min.N)))) {
    warning("Bad value of \"Min.N\" provided, defaulting to 1!")
    Min.N <- 1
  }
  if (Max.N < Min.N) {
    warning("For obvious reasons, \"Max.N\" cannot be smaller than \"Min.N\", defaulting to \"Min.N\"'s value.")
    Max.N <- Min.N
  }
  if ((!Skip.ratios)&&(Ratios.SD)) {
    warning("Argument \"Ratios.SD\" is deprecated!")
    Ratios.SD <- FALSE
  }
  if ((!Skip.ratios)&&(Ratios.Pvalue)) {
    warning("Argument \"Ratios.Pvalue\" is deprecated!")
    Ratios.Pvalue <- FALSE
  }
  Priority <- tolower(substr(Priority, 1, 3))
  if ((Skip.ratios)&&(Priority == "rat")) {
    warning("So, let me get this clear: you want me not to calculate ratios BUT to also give ratios priority? Make up your mind! Defaulting to Priority = \"Intensities\"")
    Priority <- "int"
  }
  if (!Priority %in% c("rat", "int")) {
    warning("I could not make sense of the value of argument \"Priority\", defaulting to \"Intensities\"")
    Priority <- "int"
  }
  if ((misFun(param))&&(misFun(ref.groups))) {
    stop("At least one of arguments \"param\" and \"ref.groups\" must be provided!")
  }
  if ((misFun(param))&&(misFun(ratio.groups))) {
    stop("At least one of arguments \"param\" and \"ratio.groups\" must be provided!")
  }
  if ((misFun(param))&&(misFun(sample.groups))) {
    stop("At least one of arguments \"param\" and \"sample.groups\" must be provided!")
  }
  if (!misFun(param)) {
    if ((!misFun(ref.groups))&&(gsub(";", "", param$Ratios.Ref.Groups) != ref.groups$aggregate)) {
      warning("The \"param\" and \"ref.groups\" arguments are in disagreement; the former has priority over the latter, so I shall ignore \"ref.groups\".")
    }
    ref.groups <- proteoCraft::parse.Param.aggreg(param$Ratios.Ref.Groups,
                                            aggregates,
                                            aggregate.map,
                                            aggregate.list)
    if ((!misFun(ratio.groups))&&(gsub(";", "", param$Ratios.Groups) != ratio.groups$aggregate)) {
      warning("The \"param\" and \"ratio.groups\" arguments are in disagreement; the former has priority over the latter, so I shall ignore \"ratio.groups\".")
    }
    ratio.groups <- proteoCraft::parse.Param.aggreg(param$Ratios.Groups,
                                              aggregates,
                                              aggregate.map,
                                              aggregate.list)
    sample.groups <- proteoCraft::parse.Param.aggreg(param$Volcano.plots.Aggregate.Level,
                                               aggregates,
                                               aggregate.map,
                                               aggregate.list)
  }
  #
  Pep.Intens.Nms <- paste0(Pep.Intens.root, experiments.map$Ref.Sample.Aggregate)
  Pep.Intens.Nms <- Pep.Intens.Nms[which(Pep.Intens.Nms %in% colnames(Pep))]
  stopifnot(length(Pep.Intens.Nms) > 0)
  if (!Skip.ratios) {
    Pep.Ratios.Nms <- paste0(Pep.Ratios.root, experiments.map$Ref.Sample.Aggregate)
    Pep.Ratios.Nms <- Pep.Ratios.Nms[which(Pep.Ratios.Nms %in% colnames(Pep))]
    stopifnot(length(Pep.Ratios.Nms) > 0)
  }
  #
  if ("Sequence" %in% colnames(Pep)) {
    Pep$UnmodSeq <- Pep$Sequence
  } else {
    Pep$UnmodSeq <- gsub("^_|_$|\\([^\\)]+\\)", "", Pep[[Mod.Nms]])
  }
  test <- rep(TRUE, nrow(Pep))
  # Identify peptides with modifications not included in "Mods"
  if (sum(Mods != FALSE)) {
    if (TESTING) { cat("Filtering out peptidoforms not eligible for quantitation...\n") }
    #a <- paste(c(paste0(c(AA, "_"), collapse = "|"), sapply(Mods, function(x) {paste("\\(", x, "\\)")})), collapse = "|")
    #test <- sapply(gsub(a, "", Pep[[Mod.Nms]]), nchar) == 0
    test <- gsub("\\)\\(", ",", gsub("^\\(|\\)$", "", gsub(paste(c("^_", "_$", AA), collapse = "|"), "", Pep[[Mod.Nms]])))
    test <- sapply(strsplit(test, ","), function(x) {
      x <- unlist(x)
      return(length(x[which(!x %in% Mods)]))
    }) == 0
    if (sum(!test)) {
      warning(paste0("Excluding ", sum(!test), " (", round(100*sum(!test)/nrow(Pep)),
                     "%) peptides with modifications not specifically included in quantitation parameters!"))      
    }
  }
  # Identify peptides with modifications included in "Mods.to.Exclude"
  if ((!misFun(Mods.to.Exclude))&&(nrow(Mods.to.Exclude))) {
    if (TESTING) { cat("Identifying peptides with modifications to exclude...\n") }
    Mods.to.Exclude$Pattern <- apply(Mods.to.Exclude[ , c("Mark", "Where")], 1, function(x) {
      x1 <- x[[1]]
      x2 <- x[[2]]
      x2 <- grep("^prot[NC]term", x2, invert = TRUE, value = TRUE)
      # For protein N-terminal peptides, we remove them later for those proteins for which the peptide is N-terminal!!!
      x2a <- grep("^[NC]term", x2, invert = TRUE, value = TRUE) #anywhere
      x2na <- grep("^Nterm$", x2, value = TRUE) # N-terminus, any
      x2ns <- grep("^Nterm_", x2, value = TRUE) # N-terminus, specific
      x2ca <- grep("^Cterm$", x2, value = TRUE) # C-terminus, any
      x2cs <- grep("^Cterm_", x2, value = TRUE) # C-terminus, specific
      res <- c()
      if (length(x2a)) { res <- c(res, paste0("[", paste0(x2a, collapse = ""), "]\\(", x1, "\\)")) }
      if (length(x2na)) { res <- c(res, paste0("^_\\(", x1, "\\)")) }
      if (length(x2ns)) {
        x2ns <- gsub("^Nterm_", "", x2ns)
        res <- c(res, paste0("^_\\(", x1, "\\)[", paste0(x2ns, collapse = ""), "]"))
      }
      if (length(x2ca)) { res <- c(res, paste0("\\(", x1, "\\)_$")) }
      if (length(x2cs)) {
        x2cs <- gsub("^Cterm_", "", x2cs)
        res <- c(res, paste0("[", paste0(x2cs, collapse = ""), "]\\(", x1, "\\)_$"))
      }
      return(paste(res, collapse = "|"))
    })
    Mods.to.Exclude$"Exclude protein specific" <- sapply(Mods.to.Exclude$Where, function(x) {
      sum(grepl("^prot[NC]term", unlist(x))) > 0
    })
    kol <- Mod.Nms
    if (Discard.unmod.strict) {
      # Optionally, also exclude peptides which could be modified even if they were not found to be...
      # Careful, you will lose A LOT of peptides!
      Mods.to.Exclude$Pattern <- apply(Mods.to.Exclude[, c("Mark", "Pattern")], 1, function(x) {
        gsub("\\$_", "$", gsub("\\^_", "^", gsub(paste0("\\\\\\(", x[[1]], "\\\\\\)"), "", x[[2]])))
      })
      kol <- "UnmodSeq"
    }
    tst <- nchar(Mods.to.Exclude$Pattern)
    Mods.to.Exclude$Pattern[which(tst == 0)] <- NA
    pat <- paste0(Mods.to.Exclude$Pattern[which(tst > 0)], collapse = "|")
    g <- grep(pat, Pep[[kol]])
    l <- sum(test[g])
    if (l) {
      if (Discard.unmod.strict) {
        warning(paste0("Also excluding ", l, " (", round(100*l/nrow(Pep)),
                       "%) peptides which could at least be modified with modifications to exclude for quantitation!"))
      } else {
        warning(paste0("Excluding ", l, " (", round(100*l/nrow(Pep)),
                       "%) peptides with modifications to exclude for quantitation!"))
      }
    }
    test[g] <- FALSE
  }
  # Remove counterpart peptides
  if (Discard.unmod) {
    w <- which(!test)
    if (length(w)) {
      if (TESTING) { cat("Identifying counterpart peptides...\n") }
      w <- which(Pep$UnmodSeq %in% unique(Pep$UnmodSeq[w]))
      l <- sum(test[w])
      if (l) {
        warning(paste0("Also excluding ", l, " (", round(100*l/nrow(Pep)),
                     "%) unmodified counterpart peptides!"))
      }
      test[w] <- FALSE
    }
  }
  Pep <- Pep[which(test),]
  # Do not re-order
  if ("Protein IDs" %in% colnames(Prot)) { idV <- Prot$"Protein IDs" } else {
    if ("id" %in% colnames(Prot)) { idV <- paste0("ID_", Prot$id) } else {
      idV <- paste0("Row_", 1:nrow(Prot))
    }
  }
  temp.ids <- setNames(strsplit(Prot[[Peptide.IDs]], ";"), idV)
  ord <- zigzag_ord(temp.ids)
  ord$ID <- idV
  nuOrd <- order(ord$NewOrd)
  temp.ids <- temp.ids[nuOrd]
  #
  cat(" - filtering peptides\n")
  tmp1 <- Pep[[id]]
  exports <- list("tmp1", "temp.ids")
  parallel::clusterExport(cl, exports, envir = environment())
  f0Flt <- function(x) { x[which(x %in% tmp1)] }
  environment(f0Flt) <- .GlobalEnv
  temp.ids <- setNames(parallel::parLapply(cl, temp.ids, f0Flt), names(temp.ids))
  if ((!misFun(Mods.to.Exclude))&&(nrow(Mods.to.Exclude))) {
    Mods2XclTerm <- Mods.to.Exclude[which(Mods.to.Exclude$"Exclude protein specific"),]
    Mods2XclTerm$Pattern <- apply(Mods2XclTerm[, c("Mark", "Where")], 1, function(x) {
      x1 <- x[[1]]
      x2 <- x[[2]]
      x2na <- grep("^protNterm$", x2, value = TRUE) # protein N-terminus, any
      x2ns <- grep("^protNterm_", x2, value = TRUE) # protein N-terminus, specific
      x2ca <- grep("^protCterm$", x2, value = TRUE) # protein C-terminus, any
      x2cs <- grep("^protCterm_", x2, value = TRUE) # protein C-terminus, specific
      res <- c()
      if (length(x2na)) { res <- c(res, paste0("^_\\(", x1, "\\)")) }
      if (length(x2ns)) {
        x2ns <- gsub("^protNterm_", "", x2ns)
        res <- c(res, paste0("^_\\(", x1, "\\)[", paste0(x2ns, collapse = ""), "]"))
      }
      if (length(x2ca)) { res <- c(res, paste0("\\(", x1, "\\)_$")) }
      if (length(x2cs)) {
        x2cs <- gsub("^protCterm_", "", x2cs)
        res <- c(res, paste0("[", paste0(x2cs, collapse = ""), "]\\(", x1, "\\)_$"))
      }
      return(paste(res, collapse = "|"))
    })
    pat <- paste(Mods2XclTerm$Pattern, collapse = "|")
    gy <- grep(pat, Pep[[Mod.Nms]])
    if (length(gy)) {
      tmp1 <- setNames(Prot[[Prot.Seq]][nuOrd], idV[nuOrd]) # Already re-ordered
      tmp2 <- Pep[, c(id, Mod.Nms, "UnmodSeq")]
      exports <- list("temp.ids", "tmp1", "tmp2", "id", "Mod.Nms",
                      "gy", "Discard.unmod")
      parallel::clusterExport(cl, exports, envir = environment())
      f0 <- function(x) { #x <- 1
        PrSeq <- tmp1[[x]]
        PrSeq <- c(PrSeq, gsub("^M", "", PrSeq))
        ti <- temp.ids[[x]]
        tiy <- ti[which(ti %in% tmp2[gy, id])]
        if (length(tiy)) {
          PpMdSeq <- tmp2[match(tiy, tmp2[gy, id]), Mod.Nms]
          PpMdSeqPat <- gsub("\\([^\\)]{2}\\)|^_|_$", "", gsub("\\([^\\)]{2}\\)_$", "$", gsub("^_\\([^\\)]{2}\\)", "^", PpMdSeq)))
          tst <- sapply(PpMdSeqPat, function(x) { sum(grepl(x, PrSeq)) }) > 0
          tiyR <- tiy[which(tst)]
          if (Discard.unmod) {
            PpSeq <- tmp2$UnmodSeq[match(ti, tmp2[, id])]
            PpSeqR <- tmp2$UnmodSeq[match(tiyR, tmp2[, id])]
            tiyR <- ti[which(PpSeq %in% PpSeqR)]
          }
          ti[which(!ti %in% tiyR)]
        }
        return(ti)
      }
      environment(f0) <- .GlobalEnv
      temp.ids2 <- setNames(parallel::parLapply(cl, seq_along(tmp1), f0), names(tmp1)) # Already re-ordered
      l1 <- length(unlist(temp.ids))
      l2 <- length(unlist(temp.ids2))
      if (l2 < l1) {
        warning(paste0("Also removing ", l1-l2, " cases of use of peptides with protein-terminal modifications to exclude (protein sequence-specific process, they may still be used for other protein groups if internal...)"))
      }
      temp.ids <- temp.ids2[match(names(temp.ids), names(temp.ids2))]; rm(temp.ids2)
    }
    if (Discard.unmod.strict) {
      # Exclude mods which are just on any N/C terminus regardless of associated amino-acid
      # (otherwise all peptides would go!)
      g <- grep("^\\^_\\\\\\([^\\\\)]{2}\\\\\\)$|^\\\\\\([^\\\\)]{2}\\\\\\)_\\$$", Mods2XclTerm$Pattern, invert = TRUE)
      if (length(g)) {
        Mods2XclTerm2 <- Mods2XclTerm[g,]
        Mods2XclTerm2$Pattern <- gsub("^\\^_\\\\\\([^\\\\)]{2}\\\\\\)", "^",
                                      gsub("\\\\\\([^\\\\)]{2}\\\\\\)\\$$", "$", Mods2XclTerm2$Pattern))
        pat <- paste(Mods2XclTerm2$Pattern, collapse = "|")
        gy <- grep(pat, Pep$UnmodSeq)
        if (length(gy)) {
          tmp1 <- setNames(Prot[[Prot.Seq]][nuOrd], idV[nuOrd]) # Already re-ordered
          tmp2 <- Pep[gy, id]
          exports <- list("temp.ids", "tmp1", "tmp2")
          parallel::clusterExport(cl, exports, envir = environment())
          f0 <- function(x) { #x <- 1
            PrSeq <- tmp1[[x]]
            PrSeq <- c(PrSeq, gsub("^M", "", PrSeq))
            ti <- temp.ids[[x]]
            return(ti[which(!ti %in% tmp2)])
          }
          environment(f0) <- .GlobalEnv
          temp.ids3 <- setNames(parallel::parLapply(cl, seq_along(tmp1), f0), names(tmp1)) # Already re-ordered
          l1 <- length(unlist(temp.ids))
          l3 <- length(unlist(temp.ids3))
          if (l3 < l1) {
            warning(paste0("Also removing ", l1-l3, " cases of use of peptides which could be affected by modifications to exclude (protein sequence-specific process, they may still be used for other protein groups if internal...)"))
          }
          temp.ids <- temp.ids3[match(names(temp.ids), names(temp.ids3))]; rm(temp.ids3)
        }
      }
    }
  }
  Pep <- Pep[which(Pep$id %in% unlist(temp.ids)),] # Update Pep
  #
  # Parse log bases for input and output
  # 0 = non-log transformed, 1 = default base (10 for intensities/expression, 2 for ratios))
  # Input (peptides-level) intensities
  log.Pep.Intens <- as.numeric(log.Pep.Intens)
  if (log.Pep.Intens == 1) { log.Pep.Intens <- 10}
  if (!log.Pep.Intens) {
    # Convert peptide intensities to log - they are assumed to be log from now on
    log.Pep.Intens <- 10
    if (TESTING) { cat(paste0("Converting input (peptide) intensities to default log", log.Pep.Intens, "...\n")) }
    Pep[, Pep.Intens.Nms] <- suppressWarnings(log(Pep[, Pep.Intens.Nms],
                                                  log.Pep.Intens))
  }
  # Output (PG-level) expression values
  Prot.LFQ.to.log <- as.numeric(Prot.LFQ.to.log)
  if (Prot.LFQ.to.log == 1) { Prot.LFQ.to.log <- 10 }
  #
  # Input (peptides-level) ratios
  log.Pep.Ratios <- as.numeric(log.Pep.Ratios)
  if (log.Pep.Ratios == 1) { log.Pep.Ratios <- 2 }
  if (!log.Pep.Ratios) {
    # Convert peptide ratios to log - they are assumed to be log from now on
    log.Pep.Ratios <- 2
    if (TESTING) { cat(paste0("Converting input (peptide) ratios to default log", log.Pep.Ratios, "...\n")) }
    Pep[, Pep.Ratios.Nms] <- suppressWarnings(log(Pep[, Pep.Ratios.Nms],
                                                  log.Pep.Ratios))
  }
  # Output (PG-level) ratios
  Prot.Ratios.to.log <- as.numeric(Prot.Ratios.to.log)
  if (Prot.Ratios.to.log == 1) { Prot.Ratios.to.log <- 2 }
  #
  #log.Pep.Intens;Prot.LFQ.to.log;log.Pep.Ratios;Prot.Ratios.to.log
  #
  Expr.root <- "log10(Expr.)"
  Expr.root.full <- paste0(Expr.root, " - ")
  #
  # In PreferUnique mode: update list of peptides to use for quantitation
  if (Mode == "PREFERUNIQUE") {
    if (TESTING) { cat("Updating peptides list to preferentially use unique peptides...\n") }
    tmp1 <- Pep[[id]] # Update
    temp.ids2 <- setNames(strsplit(Prot[[Unique.peptide.IDs]][nuOrd], ";"), idV[nuOrd]) # Already re-ordered
    exports <- list("tmp1", "temp.ids2")
    parallel::clusterExport(cl, exports, envir = environment())
    temp.ids2 <- setNames(parallel::parLapply(cl, temp.ids2, f0Flt), names(temp.ids2)) # List of unique peptides, remove all peptides which may have been already excluded
    tmp2 <- unlist(temp.ids2)
    exports <- list("temp.ids", "tmp2")
    parallel::clusterExport(cl, exports, envir = environment())
    f0Flt2 <- function(x) { x[which(!x %in% tmp2)] }
    environment(f0Flt2) <- .GlobalEnv
    temp.ids <- setNames(parallel::parLapply(cl, temp.ids, f0Flt2), names(temp.ids)) # List of razor/shared peptides, not including unique ones!!!
    # Calculate average intensities
    tmp1 <- Pep[, Pep.Intens.Nms, drop = FALSE]
    parallel::clusterExport(cl, "tmp1", envir = environment())
    f0 <- function(x) { mean(proteoCraft::is.all.good(x)) }
    environment(f0) <- .GlobalEnv
    AvPepInt <- parallel::parApply(cl, tmp1, 1, f0)
    # NB: we have already made sure that all unique peptides are in the razor/shared list as well.
    # Fill the gaps
    l2 <- sapply(temp.ids2, length)
    l1 <- sapply(temp.ids, length)
    wh <- which((l2 < Min.Unique)&(l1 > 0)) # No point doing this if there are no shared ones to use
    if (length(wh)) {
      tmp1 <- Pep$id
      exports <- list("temp.ids", "temp.ids2", "tmp1", "wh", "Min.Unique", "AvPepInt", "l2")
      parallel::clusterExport(cl, exports, envir = environment())
      f0 <- function(x) { #x <- wh[2] 
        l <- Min.Unique-l2[x]
        x2 <- unlist(temp.ids[x])
        avint <- AvPepInt[match(x2, tmp1)]
        w <- which(!is.na(avint))
        if (length(w)) {
          x2 <- x2[w]; avint <- avint[w]
          x2 <- x2[order(avint, decreasing = TRUE)]
          x2 <- x2[1:l]
        }
        return(c(temp.ids2[[x]], x2))
      }
      environment(f0) <- .GlobalEnv
      temp.ids2[wh] <- parallel::parLapply(cl, wh,  f0)
    } 
    temp.ids <- temp.ids2[match(names(temp.ids), names(temp.ids2))]; rm(temp.ids2)
  }
  Pep <- Pep[which(Pep$id %in% unlist(temp.ids)),] # Update Pep
  # Get summary method:
  if (misFun("Summary.weights")) {
    if (Summary.method == "weighted.mean") {
      warning("Summary method is weighted mean, but no weights were provided...?")
    }
    Summary.weights <- "Summary.weights"
    Pep[[Summary.weights]] <- 1
  }
  #
  # Optional intensity weights (so a peptide's weight is a function of its non-log intensity)
  if (Intensity.weights) {
    tmp <- log.Pep.Intens^Pep[, Pep.Intens.Nms, drop = FALSE]
    parallel::clusterExport(cl, "tmp", envir = environment())
    f0 <- function(x) {
      x <- proteoCraft::is.all.good(x)
      x <- x[which(x > 0)]
      return(mean(x))
    }
    environment(f0) <- .GlobalEnv
    Pep$Intensity.weights <- parallel::parApply(cl, tmp, 1, f0)
    Pep[[Summary.weights]] <- Pep[[Summary.weights]]*Pep$Intensity.weights # Update weights and method
    Summary.method <- "weighted.mean"
  }
  sum.func <- get(Summary.method)
  #
  # Calculate ratios (method used if Priority = "Ratios"):
  if (Priority == "rat") {
    if (TESTING) { cat("Calculating ratios...\n") }
    tmpPep <- Pep[, c(id, Summary.weights, Pep.Ratios.Nms, Pep.Intens.Nms)]
    ratios.fun <- function(x) { #x <- temp.ids[1]
      res <- rep(NA, length(Pep.Ratios.Nms)*(1+Ratios.SD+Ratios.Pvalue)) # This is the default vector of NAs to replace if the data is complete enough
      mtch <- match(as.numeric(unlist(x)), tmpPep[[id]])
      if (length(mtch)) {
        temp <- tmpPep[mtch, c(Pep.Ratios.Nms, Pep.Intens.Nms), drop = FALSE]
        if (Summary.method == "weighted.mean") { temp[[Summary.weights]] <- tmpPep[mtch, Summary.weights] }
        if (length(Pep.Ratios.Nms) >= Min.N) {
          tst1 <- apply(temp[, Pep.Ratios.Nms, drop = FALSE], 1, function(y) { length(proteoCraft::is.all.good(y)) }) > 0
          if (sum(tst1)) {
            temp <- temp[which(tst1), , drop = FALSE]
            if (sum(tst1) > Max.N) { # Remove excess peptides, for cases where we have much more than we need and including all would slow down the calculations (Titin!!!)
              temp$Rank <- nrow(temp) + 1 - rank(apply(temp[, Pep.Intens.Nms, drop = FALSE], 1, function(y) { sum(proteoCraft::is.all.good(y)) }))
              temp <- temp[which(temp$Rank <= Max.N), , drop = FALSE]
              temp$Rank <- NULL
            }
            res <- sapply(Pep.Ratios.Nms, function(y) { #y <- Pep.Ratios.Nms[[1]]
              y1 <- temp[[y]]
              length(y1)
              if (Summary.method == "weighted.mean") { rs <- sum.func(y1, temp[[Summary.weights]]) } else {
                rs <- sum.func(y1)
              }
              if (Ratios.SD) { # (legacy code, currently ignored)
                rs <- c(rs, sd(y1))
              }
              if (Ratios.Pvalue) { # (legacy code, currently ignored)
                if (length(unique(y1)) > 1) {
                  tt <- try(-log10(t.test(x = y1, y = NULL, alternative = "two.sided")$p.value), silent = TRUE)
                  rs <- c(rs, c(tt, NA)[(class(tt) == "try-error")+1])
                } else { rs <- c(rs, NA) }
              }
              return(rs)
            })
          }
        }
      }
      return(res)
    }
    exports <- list("temp.ids", "tmpPep", "Pep.Ratios.Nms", "Min.N", "id", "Ratios.SD", "Ratios.Pvalue", "ratios.fun", 
                    "Summary.method", "sum.func", "Max.N", "Pep.Intens.Nms", "Summary.weights", "t.test")
    parallel::clusterExport(cl, exports, envir = environment())
    f0 <- function(x) { ratios.fun(unlist(x)) }
    environment(f0) <- .GlobalEnv
    res <- parallel::parSapply(cl, temp.ids, f0, USE.NAMES = TRUE)
    if ((length(Pep.Ratios.Nms) > 1)||(Ratios.SD+Ratios.Pvalue)) { res <- as.data.frame(t(res)) } else { res <- as.data.frame(res) }
    rownames(res) <- names(temp.ids)
    # Column names:
    cola <- ""
    if (Ratios.SD) { cola <- c(cola,  ": SD") } # (legacy code, currently ignored)
    if (Ratios.Pvalue) { cola <- c(cola, ": -log10(peptides Pvalue)") } # (legacy code, currently ignored)
    colnames(res) <- sapply(Pep.Ratios.Nms, function(x) {
      sapply(cola, function(y) { gsub(proteoCraft::topattern(gsub(" - $", "", Pep.Ratios.root)), paste0(gsub(" - $", "", Pep.Ratios.root), y), x) })
    })
    # Reorder:
    res <- res[, sapply(cola, function(x) {
      gsub(proteoCraft::topattern(gsub(" - $", "", Pep.Ratios.root)), paste0(gsub(" - $", "", Pep.Ratios.root), x), Pep.Ratios.Nms)
    })]
    #tst <- setNames(apply(res[, grep(": (SD|-log10\\(peptides Pvalue\\))", colnames(res), invert = TRUE)], 2, median, na.rm = TRUE), Pep.Ratios.Nms)
    #tst
  }
  # Calculate absolute abundance values for all conditions:
  ## This uses the Levenberg-Marquardt procedure to align peptide profiles and is based on the "best flyer" hypothesis.
  #if (TESTING) {
    cat(" - exporting to cluster...\n")
  #}
  Viz <- FALSE
  tmpPep <- Pep[, c(id, Summary.weights, Pep.Intens.Nms)]
  exports <- list("temp.ids", "tmpPep", "id", "Pep.Intens.Nms", "Min.N", "Max.N",
                  "Summary.method", "Summary.weights", "Viz")
  parallel::clusterExport(cl, exports, envir = environment())
  #db$`Common Name`[match(prot.list, db$`Protein ID`)]
  #match(prot.list, db$`Protein ID`)
  #proteoCraft::grsep(prot.list[2], x = Prot$`Leading protein IDs`)
  #View(ord)
  f0 <- function(ids) {
    proteoCraft::LFQ.lm(ids,
                        InputTabl = tmpPep,
                        id = id,
                        IntensCol = Pep.Intens.Nms,
                        Summary.method = Summary.method,
                        Summary.weights = Summary.weights,
                        Min.N = Min.N,
                        Max.N = Max.N,
                        Viz = Viz)
  }
  environment(f0) <- .GlobalEnv
  #if (TESTING) {
    cat(" - running quant algorithm...\n")
  #}
  res2 <- try(parallel::parSapply(cl, temp.ids, f0, USE.NAMES = TRUE), silent = TRUE)
  if ("try-error" %in% class(res2)) {
    cat(" - re-running quant algorithm, the slow way...\nsomething clearly went wrong with the cluster\n")
    # This function has occasionally and non-reproducibly failed on weak PCs, this is a back up
    res2 <- sapply(temp.ids, f0, USE.NAMES = TRUE)
  }
  if (length(Pep.Intens.Nms) > 1) { res2 <- as.data.frame(t(res2)) } else { res2 <- as.data.frame(res2) } #as.data.frame is weird...
  colnames(res2) <- gsub(proteoCraft::topattern(Pep.Intens.root), Expr.root.full, Pep.Intens.Nms)
  rownames(res2) <- names(temp.ids)
  #
  if (Priority == "rat") { res2[, colnames(res)] <- res[match(rownames(res2), rownames(res)),] }
  if (!Skip.ratios) {
    # Calculate ratios references
    kount <- 0
    for (i in ref.groups$values) { #i <- ref.groups$values[1]
      j <- setNames(unlist(strsplit(i, "___")), ref.groups$names)
      temp <- sapply(ref.groups$names, function(x) { list(which(experiments.map[[x]] == j[[x]])) })
      temp2 <- sort(unique(unlist(temp)))
      test <- sapply(temp2, function(x) { sum(sapply(temp, function(y) {x %in% unlist(y)})) })
      temp2 <- temp2[which(test == length(temp))]
      temp3 <- experiments.map[temp2,]
      temp3 <- temp3[which(temp3$Reference),]
      b <- temp3$Ref.Sample.Aggregate
      if (length(b)) {
        kol <- paste0(Expr.root.full, b)
        w <- which(kol %in% colnames(res2))
        if (length(w)) {
          b1 <- as.numeric(apply(res2[, kol[w], drop = FALSE], 1, proteoCraft::log_ratio_av)) # (Expression data is already log-transformed for Prot.Quant)
          b2 <- paste0(Expr.root.full, i, ".REF")
          #print(b2)
          res2[[b2]] <- b1 # log-transformed
          kount <- kount + 1
          if (kount == 1) {
            rat.2.ref <- data.frame(Name = b2, Source = paste(b, collapse = ";"))
          } else { rat.2.ref <- rbind(rat.2.ref, c(b2, paste(b, collapse = ";"))) }
        } else { warning(paste0("Empty group: ", i, ", skipping!")) }
      } else { warning(paste0("There is no reference for level ", i)) }
    }
    #View(res2[,grep("\\.REF$", colnames(res2))])
    rat.2.ref$Source <- strsplit(rat.2.ref$Source, ";")
    rat.2.ref$Used_by <- list(NA)
    # Re-calculate individual ratios/intensities to relevant reference, depending on priority
    # If the reference is an average, we will also calculate individual ref ratios to it;
    # this will be useful further down the line.
    for (i in ref.groups$values) { #i <- ref.groups$values[1]
      j <- setNames(unlist(strsplit(i, "___")), ref.groups$names)
      temp <- sapply(ref.groups$names, function(x) { list(which(experiments.map[[x]] == j[[x]])) })
      temp2 <- sort(unique(unlist(temp)))
      test <- sapply(temp2, function(x) { sum(sapply(temp, function(y) {x %in% unlist(y)})) })
      temp2 <- temp2[which(test == length(temp))]
      temp3 <- experiments.map[temp2,]
      # Get reference
      k <- j[which(names(j) %in% ref.groups$names)]
      b <- paste0(Expr.root.full, paste(k, collapse = "___"), ".REF")
      b1 <- res2[[b]] # Log-transformed!!!
      A <- a <- temp3$Ref.Sample.Aggregate
      if (length(which(temp3$Reference)) == 1) {
        # If there is only one ref in the group, remove it as there is no point calculating a ratio to itself,
        # or re-calculating its expression to itself.
        a <- temp3$Ref.Sample.Aggregate[which(!temp3$Reference)]
      }
      a <- a[which(paste0(Expr.root.full, a) %in% colnames(res2))]
      if (!length(a)) {
        warning(paste0(paste0("There are no ", c("ratios", "expression values")[match(Priority, c("int", "rat"))], " to calculate for level ", i)))
      } else {
        if (Priority == "int") { # Here we calculate log ratios so they precisely reflect log expression
          res2[, paste0(Pep.Ratios.root, a)] <- sweep(res2[, paste0(Expr.root.full, a), drop = FALSE],
                                                      1, b1, "-")/log(log.Pep.Ratios, log.Pep.Intens)
        }
        if (Priority == "rat") { # Here we re-calculate log expression to reflect log ratios
          res2[, paste0(Expr.root.full, a)] <- sweep(res2[, paste0(Pep.Ratios.root, a), drop = FALSE]*log(log.Pep.Ratios, log.Pep.Intens),
                                                     1, b1, "+")
        }
        PrintSummary <- TRUE
        if (PrintSummary) {
          for (a1 in A) { #a1 <- A[1]
            tmp <- proteoCraft::is.all.good(res2[[paste0(Expr.root.full, a1)]])
            if (!log.Pep.Intens) { tmp <- tmp[which(tmp > 0)] }
            emd <- signif(median(tmp), 3)
            esd <- signif(sd(tmp), 3)
            msg <- paste0("-> Sample ", cleanNms(a1), ":\n", paste0("   - log", log.Pep.Intens, "(expression): median = ", emd, ", SD = ", esd, "\n"))
            if (paste0(Pep.Ratios.root, a1) %in% colnames(res2)) {
              tmp <- proteoCraft::is.all.good(res2[[paste0(Pep.Ratios.root, a1)]])
              rmd <- signif(median(tmp), 3)
              rsd <- signif(sd(tmp), 3)
              msg <- paste0(msg, paste0("   - log", log.Pep.Ratios, "(ratio): median = ", rmd, ", SD = ", rsd, "\n"))
            } else { msg <- paste0(msg, "     (no ratios calculated)\n") }
            cat(msg)
          }
        }
      }
      rat.2.ref[["Used_by"]][which(rat.2.ref$Name == b)] <- list(a)
    }
  }
  # Calculate Reference-to-Reference ratio vectors:
  # Mirrored from the code after "# Create peptide-level Ref-to-Ref ratios (useful for PTMs analysis):" in the main replicates-workflow
  # The main difference is that the input data is log-transformed (not for peptides).
  if (!Skip.ratios) {
    if (!misFun(param)) {
      rat_cont_grps <- param$Ratios.Contaminant.Groups
    } else { rat_cont_grps <- "Ratio groups" } # Default
    res3 <- try(proteoCraft::make_RefRat(data = res2,
                                   experiment.map = experiment.map,
                                   int.root = Expr.root.full,
                                   rat.root = Pep.Ratios.root,
                                   rat.con.grps = rat_cont_grps,
                                   mode = Refs_Mode,
                                   parameters = param,
                                   logInt = log.Pep.Intens,
                                   logRat = log.Pep.Ratios), silent = TRUE)
    if (!"try-error" %in% class(res3)) {
      res2[, colnames(res3)] <- res3
    } else { warning("No Ref to Ref columns were generated!") }
  }
  # Finally: this should be the very last step to avoid confusions:
  # If output should be log-transformed, transform log and add log base to column names
  if (!Prot.LFQ.to.log) { # De-log LFQ values
    g <- grep(proteoCraft::topattern(Expr.root.full), colnames(res2), value = TRUE)
    res2[, g] <- log.Pep.Intens^(res2[, g])
    colnames(res2) <- gsub(proteoCraft::topattern(Expr.root.full),
                           "Expr.",
                           colnames(res2))
  }
  if (!Skip.ratios) {
    # For ratio values, either... 
    I <- colnames(res2)[which(!grepl(": SD$|: -log10\\(peptides Pvalue\\)$", colnames(res2)))]
    if (Prot.Ratios.to.log) {
      if (Prot.Ratios.to.log != log.Pep.Ratios) {# ... if necessary change to the desired base:
        res2[, I] <- res2[, I]/log(Prot.Ratios.to.log, log.Pep.Ratios) }
      colnames(res2) <- gsub(proteoCraft::topattern(Pep.Ratios.root),
                             paste0("log", Prot.Ratios.to.log, "(Ratio) - "), colnames(res2))
    } else { #... or remove log from name and de-log:
      res2[, I] <- log.Pep.Ratios^res2[, I]
      colnames(res2) <- gsub(proteoCraft::topattern(Pep.Ratios.root), "Ratio - ", colnames(res2))
    }
  }
  for (i in 1:ncol(res2)) {
    if (!"numeric" %in% class(res2[[i]])) {
      stop(paste0("I would expect the class of column ", i, " to be numeric! Investigate!"))
      #res2[[i]] <- as.numeric(res2[[i]])
    }
  }
  res2$"Peptides IDs used for quantitation" <- sapply(temp.ids, paste, collapse = ";")
  ord2 <- ord[nuOrd,]
  res2 <- res2[order(ord2$Original),]
  rownames(res2) <- rownames(Prot)
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(res2)
}
