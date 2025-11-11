#' Prot.Quant2
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
#' @param id The name of the Peptides table's IDs column. Default = "id"
#' @param Summary.method The summary method used for ratios (the Levenberg-Marquardt algorithm is used for intensities). One of "mean", "median", or "weighted.mean".
#' @param Summary.weights If a "weighted.mean" summary method is chosen, then a vector of weights must be provided. This should be the name of a column of the peptides table
#' @param Priority One of "Ratios" or "Intensities" (default). Some flexibility in spelling is allowed. You want to prioritize ratios for SILAC because in this case MaxQuant measures peptides ratios more precisely than intensities. Otherwise, you want to prioritize Intensities and re-calculate ratios from them.
#' @param Skip.ratios Logical, default = FALSE. If TRUE, ratios will not be calculated. 
#' @param experiments.map Map of the experiment, default = Exp.map
#' @param param The experiment's parameters object If provided, the ref.groups argument is not required.
#' @param aggregate.map The aggregate map. Default = Aggregate.map
#' @param aggregate.list The named list of aggregates. Default = Aggregate.list
#' @param aggregates The aggregates themselves. Default = Aggregates
#' @param ref.groups Defines which samples are paired to which references. May alternatively (preferred solution) be provided indirectly through the param argument.
#' @param ratio.groups Defines groups within which ratios are calculated. May alternatively (preferred solution) be provided indirectly through the param argument.
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
#' 
#' @examples
#' temp <- Prot.Quant2(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep, id = "New Peptide ID",
#'                    experiments.map = Exp.map,
#'                    ref.groups = Ratios.Ref.Groups,
#'                    Pep.Intens.root = pep.ref, Pep.Ratios.root = pep.ratios.ref,
#'                    log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
#'                    Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
#'                    Mods = Mod4Quant, Mod.Nms = "Modified sequence",
#'                    Min.N = 2)
#' 
#' @export

Prot.Quant2 <- function(Prot,
                        Peptide.IDs = "Peptide IDs",
                        Pep,
                        id = "id",
                        Summary.method = "median",
                        Summary.weights,
                        Priority = "Intensities",
                        Skip.ratios = FALSE,
                        experiments.map,
                        ref.groups,
                        ratio.groups,
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
                        Min.N = 2) {
  #proteoCraft::DefArg(proteoCraft::Prot.Quant2); ref.groups <- proteoCraft::parse.Param.aggreg(Param$Ratios.Ref.Groups, aggregates, aggregate.map, aggregate.list); ratio.groups <- proteoCraft::parse.Param.aggreg(Param$Ratios.Groups, aggregates, aggregate.map, aggregate.list)
  #Prot = PG; Peptide.IDs = Pep4Quant; Pep = pep; Summary.method = "weighted.mean"; Summary.weights = "Weights"; experiments.map = Exp.map; param = Param; Pep.Intens.root = pep.ref[length(pep.ref)]; Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)]; Mods = Mod4Quant; Mods.to.Exclude = Mod2Xclud; Min.N = 1; Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1]
  #Prot = PG; Peptide.IDs = Pep4Quant; Pep = pep; ref.groups = RefGrp; ratio.groups = RatGrp; Skip.ratios = !MakeRatios; experiments.map = Exp.map; Pep.Intens.root = "Intensity - "; Mods = Mod4Quant; Mods.to.Exclude = Mod2Xclud; log.Pep.Intens = FALSE; log.Pep.Ratios = 2; Prot.LFQ.to.log = TRUE; Min.N = 1; if (!Skip.ratios) { Pep.Ratios.root = rat.col }
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
  Priority <- tolower(substr(Priority, 1, 3))
  if ((Skip.ratios)&&(Priority == "rat")) {
    warning("So, let me get this clear: you want me not to calculate ratios BUT to also give ratios priority? Make up your mind! Defaulting to Priority = \"Intensities\"")
    Priority <- "int"
  }
  if (!Priority %in% c("rat", "int")) {
    warning("I could not make sense of the value of argument \"Priority\", defaulting to \"Intensities\"")
    Priority <- "int"
  }
  if ((missing("param"))&&(missing("ref.groups"))) {
    stop("At least one of arguments \"param\" and \"ref.groups\" must be provided!")
  }
  if ((missing("param"))&&(missing("ratio.groups"))) {
    stop("At least one of arguments \"param\" and \"ratio.groups\" must be provided!")
  }
  if (!missing("param")) {
    if ((!missing("ref.groups"))&&(param$Ratios.Ref.Groups != ref.groups)) {
      warning("The \"param\" and \"ref.groups\" arguments are in disagreement; the former has priority over the latter, so I shall ignore \"ref.groups\".")
    }
    ref.groups <- proteoCraft::parse.Param.aggreg(param$Ratios.Ref.Groups,
                                            aggregates,
                                            aggregate.map,
                                            aggregate.list)
    if ((!missing("ratio.groups"))&&(param$Ratios.Groups != ratio.groups)) {
      warning("The \"param\" and \"ratio.groups\" arguments are in disagreement; the former has priority over the latter, so I shall ignore \"ratio.groups\".")
    }
    ratio.groups <- proteoCraft::parse.Param.aggreg(param$Ratios.Groups,
                                              aggregates,
                                              aggregate.map,
                                              aggregate.list)
  }
  # Get summary method:
  if (missing("Summary.weights")) {
    if (Summary.method == "weighted.mean") {
      warning("Summary method is weighted mean, but no weights were provided...?")
    }
    Summary.weights <- "Summary.weights"
    Pep$Summary.weights <- 1
  }
  sum.func <- get(Summary.method)
  #
  Pep.Intens.Nms <- paste0(Pep.Intens.root, experiments.map$Ref.Sample.Aggregate)
  Pep.Intens.Nms <- Pep.Intens.Nms[which(Pep.Intens.Nms %in% colnames(Pep))]
  stopifnot(length(Pep.Intens.Nms) > 0)
  if (!Skip.ratios) {
    Pep.Ratios.Nms <- paste0(Pep.Ratios.root, experiments.map$Ref.Sample.Aggregate)
    Pep.Ratios.Nms <- Pep.Ratios.Nms[which(Pep.Ratios.Nms %in% colnames(Pep))]
    stopifnot(length(Pep.Ratios.Nms) > 0)
  }
  Pep$UnmodSeq <- gsub("^_|_$|\\([^\\)]+\\)", "", Pep[[Mod.Nms]])
  test <- rep(TRUE, nrow(Pep))
  # Identify peptides with modifications not included in "Mods"
  if (sum(Mods != FALSE)) {
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
  if ((!missing(Mods.to.Exclude))&&(nrow(Mods.to.Exclude))) {
    Mods.to.Exclude$Pattern <- apply(Mods.to.Exclude[, c("Mark", "Where")], 1, function(x) {
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
  temp.ids <- strsplit(Prot[[Peptide.IDs]], ";")
  temp.ids <- lapply(temp.ids, function(x) { x[which(x %in% Pep[[id]])] })
  if ((!missing(Mods.to.Exclude))&&(nrow(Mods.to.Exclude))) {
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
      temp.ids2 <- lapply(1:nrow(Prot), function(x) { #x <- 1
        PrSeq <- Prot[x, Prot.Seq]
        PrSeq <- c(PrSeq, gsub("^M", "", PrSeq))
        ti <- temp.ids[[x]]
        tiy <- ti[which(ti %in% Pep[gy, id])]
        if (length(tiy)) {
          PpMdSeq <- Pep[match(tiy, Pep[gy, id]), Mod.Nms]
          PpMdSeqPat <- gsub("\\([^\\)]{2}\\)|^_|_$", "", gsub("\\([^\\)]{2}\\)_$", "$", gsub("^_\\([^\\)]{2}\\)", "^", PpMdSeq)))
          tst <- sapply(PpMdSeqPat, function(x) { sum(grepl(x, PrSeq)) }) > 0
          tiyR <- tiy[which(tst)]
          if (Discard.unmod) {
            PpSeq <- Pep$UnmodSeq[match(ti, Pep[, id])]
            PpSeqR <- Pep$UnmodSeq[match(tiyR, Pep[, id])]
            tiyR <- ti[which(PpSeq %in% PpSeqR)]
          }
          ti[which(!ti %in% tiyR)]
        }
        return(ti)
      })
      l1 <- length(unlist(temp.ids))
      l2 <- length(unlist(temp.ids2))
      if (l2 < l1) {
        warning(paste0("Also removing ", l1-l2, " cases of use of peptides with protein-terminal modifications to exclude (protein sequence-specific process, they may still be used for other protein groups if internal...)"))
      }
      temp.ids <- temp.ids2; rm(temp.ids2)
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
          temp.ids3 <- lapply(1:nrow(Prot), function(x) { #x <- 1
            PrSeq <- Prot[x, Prot.Seq]
            PrSeq <- c(PrSeq, gsub("^M", "", PrSeq))
            ti <- temp.ids[[x]]
            return(ti[which(!ti %in% Pep[gy, id])])
          })
          l1 <- length(unlist(temp.ids))
          l3 <- length(unlist(temp.ids3))
          if (l3 < l1) {
            warning(paste0("Also removing ", l1-l3, " cases of use of peptides which could be affected by modifications to exclude (protein sequence-specific process, they may still be used for other protein groups if internal...)"))
          }
          temp.ids <- temp.ids3; rm(temp.ids3)
        }
      }
    }
  }
  # Parse log bases for input and output
  # 0 = non-log transformed, 1 = default base (10 for intensities/expression, 2 for ratios))
  # Input (peptides-level) intensities
  log.Pep.Intens <- as.numeric(log.Pep.Intens)
  if (log.Pep.Intens == 1) { log.Pep.Intens <- 10 }
  if (log.Pep.Intens) {
    # De-log peptide intensities - they are assumed to be non-log transformed from now on
    message("Reverting the log-transformation on input (peptide) intensities...")
    Pep[, Pep.Intens.Nms] <- suppressWarnings(log.Pep.Intens^Pep[, Pep.Intens.Nms])
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
    message(paste0("Converting input (peptide) ratios to default log", log.Pep.Ratios, "..."))
    Pep[, Pep.Ratios.Nms] <- suppressWarnings(log(Pep[, Pep.Ratios.Nms], log.Pep.Ratios))
  }
  # Output (PG-level) ratios
  Prot.Ratios.to.log <- as.numeric(Prot.Ratios.to.log)
  if (Prot.Ratios.to.log == 1) { Prot.Ratios.to.log <- 2 }
  #
  #log.Pep.Intens;Prot.LFQ.to.log;log.Pep.Ratios;Prot.Ratios.to.log
  #
  Expr.root <- "Expr."
  Expr.root.full <- paste0(Expr.root, " - ")
  # Calculate ratios (method used if Priority = "Ratios"):
  if (Priority == "rat") {
    ratios.fun <- function(x) { #x <- temp.ids[1]
      res <- rep(NA, length(Pep.Ratios.Nms)) # This is the default vector of NAs to replace if the data is complete enough
      mtch <- match(as.numeric(unlist(x)), Pep[[id]])
      if (length(mtch)) {
        temp <- Pep[mtch, c(Pep.Ratios.Nms, Pep.Intens.Nms, Summary.weights), drop = FALSE]
        if (length(Pep.Ratios.Nms) >= Min.N) {
          tst1 <- apply(temp[,Pep.Ratios.Nms, drop = FALSE], 1, function(y) {
            length(proteoCraft::is.all.good(y))
          }) > 0
          if (sum(tst1)) {
            temp <- temp[which(tst1), , drop = FALSE]
            res <- sapply(Pep.Ratios.Nms, function(y) { #y <- Pep.Ratios.Nms[[1]]
              y1 <- temp[[y]]
              length(y1)
              if (Summary.method == "weighted.mean") { rs <- sum.func(y1, temp[[Summary.weights]]) } else {
                rs <- sum.func(y1)
              }
              return(rs)
            })
          }
        }
      }
      return(res)
    }
    res <- sapply(temp.ids, function(x) { ratios.fun(unlist(x)) })
    if (length(Pep.Ratios.Nms) > 1) { res <- as.data.frame(t(res)) } else { res <- as.data.frame(res) }
    # Column names:
    cola <- ""
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
  LFQ2 <- function(i) {
    #for (i in 1:length(temp.ids)) {# i <- 1
    #print(i)
    rs <- rep(NA, length(Pep.Intens.Nms)) # This is the default vector of NAs to replace if the data is complete enough
    x <- temp.ids[[unlist(i)]]
    # Check that we have at least enough peptidoforms:
    # NB: Currently done at peptidoforms level, i.e. accepts two peptidoforms of the same primary sequence as different.
    # If I ever decide to change it, this would be the place.
    mtch <- match(unlist(x), Pep[[id]])
    if (length(mtch) >= Min.N) {
      temp2 <- Pep[mtch, Pep.Intens.Nms, drop = FALSE]
      extWeights <- Pep[mtch, Summary.weights]
      w1 <- which(apply(temp2[, Pep.Intens.Nms, drop = FALSE], 1, function(y) {
        length(proteoCraft::is.all.good(log(y)))
      }) > 0) # (log to also exclude 0s)
      if (length(w1)) {
        # Remove peptides with only non-valid or missing values
        temp2 <- temp2[w1, , drop = FALSE]; extWeights <- extWeights[w1]
        tst2 <- sapply(Pep.Intens.Nms, function(y) {
          length(proteoCraft::is.all.good(log(temp2[[y]])))
        }) # (log to also exclude 0s)
        # We should average unscaled peptide values, using weighted mean,
        # which presents a double advantage:
        # - This approach will give more weight to higher intensity (less noisy) peptides.
        # - We may not have the same peptides in all samples; w weighted mean using individual peptide intensity range
        #   allows us to correct for this.
        #
        # Estimate individual peptide average intensity range
        intWeights <- apply(temp2, 1, function(x) {
          x <- proteoCraft::is.all.good(x); x <- x[which(x > 0)] # Filter bad values
          mean(x)
        })
        w2 <- which(intWeights == max(intWeights)) # Identify highest intensity value(s)
        #w1 <- which(tst2 == max(tst2)) # (identify longest vector)
        intScale <- intWeights[w2[1]] # Scale for absolute quant (best flyer hypothesis)
        rs <- apply(temp2, 2, function(x) {
          x <- as.numeric(x)
          w3 <- which(proteoCraft::is.all.good(log(x), 2)) # Exclude nulls and non-valid values
          if (length(w3)) { x <- weighted.mean(x[w3], intWeights[w3]*extWeights[w3]) } else {
            x <- 0
          }
          return(x)
        })
        rs <- rs*intScale/mean(rs[which(rs > 0)])
      }
    }
    return(rs)
  }
  res2 <- sapply(1:length(temp.ids), LFQ2)
  if (length(Pep.Intens.Nms) > 1) { res2 <- as.data.frame(t(res2)) } else { res2 <- as.data.frame(res2) } #as.data.frame is weird...
  colnames(res2) <- gsub(proteoCraft::topattern(Pep.Intens.root), Expr.root.full, Pep.Intens.Nms)
  #
  if (Priority == "rat") { res2[, colnames(res)] <- res }
  if (!Skip.ratios) {
    # Calculate ratios references
    kount <- 0
    for (i in ref.groups$values) { #i <- ref.groups$values[1]
      j <- setNames(unlist(strsplit(i, "___")), ref.groups$names)
      temp <- sapply(ref.groups$names, function(x) {
        list(which(experiments.map[[x]] == j[[x]]))
      })
      temp2 <- sort(unique(unlist(temp)))
      test <- sapply(temp2, function(x) {
        sum(sapply(temp, function(y) {
          x %in% unlist(y)
        }))
      })
      temp2 <- temp2[which(test == length(temp))]
      temp3 <- experiments.map[temp2,]
      temp3 <- temp3[which(temp3$Reference),]
      b <- temp3$Ref.Sample.Aggregate
      if (length(b)) {
        kol <- paste0(Expr.root.full, b)
        w <- which(kol %in% colnames(res2))
        if (length(w)) {
          b1 <- 10^as.numeric(apply(res2[, kol[w], drop = FALSE], 1, function(x) {
            proteoCraft::log_ratio_av(log10(x)) # (Expression data is non log-transformed for Prot.Quant2)
          }))
          b2 <- paste0(Expr.root.full, i, ".REF")
          #print(b2)
          res2[[b2]] <- b1 # non log-transformed
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
      b1 <- res2[[b]] # Not logged
      A <- a <- temp3$Ref.Sample.Aggregate
      if (length(which(temp3$Reference)) == 1) {
        # If there is only one ref in the group, remove it as there is no point calculating a ratio to itself,
        # or re-calculating its expression to itself.
        a <- temp3$Ref.Sample.Aggregate[which(!temp3$Reference)]
      }
      if ((!length(a))||(!sum(paste0(Expr.root.full, a) %in% colnames(res2)))) {
        warning(paste0(paste0("There are no ", c("ratios", "expression values")[match(Priority, c("int", "rat"))], " to calculate for level ", i)))
      } else {
        if (Priority == "int") { # Here we calculate log ratios so they precisely reflect expression
          res2[, paste0(Pep.Ratios.root, a)] <- log(sweep(res2[, paste0(Expr.root.full, a), drop = FALSE], 1, b1, "/"),
                                                    log.Pep.Ratios)
        }
        if (Priority == "rat") { # Here we re-calculate expression to reflect log ratios
          res2[, paste0(Expr.root.full, a)] <- sweep(2^res2[, paste0(Pep.Ratios.root, a), drop = FALSE],
                                                     1, b1, "*")
        }
        PrintSummary <- TRUE
        if (PrintSummary) {
          for (a1 in A) { #a1 <- A[1]
            tmp <- proteoCraft::is.all.good(res2[[paste0(Expr.root.full, a1)]])
            if (!log.Pep.Intens) { tmp <- tmp[which(tmp > 0)] }
            emd <- signif(median(tmp), 3)
            esd <- signif(sd(tmp), 3)
            msg <- paste0("-> Sample ", a1, ":\n", paste0("   - expression: median = ", emd, ", SD = ", esd, "\n"))
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
  if (!Skip.ratios) {
    # Calculate Reference-to-Reference ratio vectors:
    res3 <- lapply(ratio.groups$values, function(x) { #x <- ratio.groups$values[1] # Within each ratios group
      m <- experiments.map[which((experiments.map[[ratio.groups$column]] == x)&(experiments.map$Reference)),]
      rs <- NA
      if (nrow(m)) {
        if (nrow(m) > 1) {
          x0 <- unique(m[[ref.groups$column]])
          xr0 <- paste0(Expr.root.full, x0, ".REF")
          w <- which(xr0 %in% colnames(res2))
          if (length(w)) {
            x0 <- x0[w]
            xr0 <- xr0[w]
            if (ratio.groups$aggregate != ref.groups$aggregate) {
              # This is the case where there is no average reference and we need to do permutations
              if (length(unique(m$Ref.Sample.Aggregate)) > 1) {
                perm <- gtools::permutations(length(xr0), 2, xr0)
                rs <- apply(perm, 1, function(y) {
                  log(res2[[y[[1]]]]/res2[[y[[2]]]], log.Pep.Ratios)
                })
              } else { stop("Unexpected situation: there should always be replicates of reference conditions per ratio group!") }
            } else {
              # This is the case where there is an average reference
              rs <- setNames(lapply(1:length(x0), function(y) { #y <- 1
                y2 <- unique(m$Ref.Sample.Aggregate[which(m[[ref.groups$column]] == x0[y])])
                xry <- paste0(Expr.root.full, y2)
                wy <- which(xry %in% colnames(res2))
                xry <- xry[wy]
                y <- log(sweep(res2[, xry, drop = FALSE], 1, res2[[xr0[y]]], "/"), log.Pep.Ratios)
              }), x0)
            }
          }
        } else { warning(paste0("Only one reference column for group \"", x,"\", no \"ref-to-ref\" ratios will be calculated!")) }
      } else { warning(paste0("Not a single reference column for group \"", x,"\", yet I am supposed to be calculating ratios!? Investigate!!!")) }
      if ((length(rs) > 1)||(!is.na(rs))) {
        rs <- as.data.frame(rs)
        rs <- suppressWarnings(log(rs, log.Pep.Ratios))
        colnames(rs) <- paste0(Pep.Ratios.root, x, "_REF.to.REF_", 1:ncol(rs))
      }
      return(rs)
    })
    res3 <- res3[which(sapply(res3, function(x) { sum(class(x) %in% c("data.frame", "matrix", "array")) }) > 0)]
    L <- length(res3)
    if (L) {
      res4 <- res3[[1]] # Not the most elegant solution but avoids some weird name corruption effects!
      if (L > 1) { for (l in 2:L) { res4 <- cbind(res4, res3[[l]]) } }
      for (k in colnames(res4)) { res2[[k]] <- as.numeric(res4[,k]) }
    }
  }
  # Finally: this should be the very last step to avoid confusions:
  # If output should be log-transformed, transform log and add log base to column names
  if (Prot.LFQ.to.log) {
    g <- grep(proteoCraft::topattern(Expr.root.full), colnames(res2), value = TRUE)
    res2[, g] <- log(res2[, g], Prot.LFQ.to.log)
    colnames(res2) <- gsub(proteoCraft::topattern(Expr.root.full),
                           paste0("log", Prot.LFQ.to.log, "(", Expr.root, ") - "),
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
  return(res2)
}
