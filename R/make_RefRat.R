#' make_RefRat
#'
#' @description
#' A function to create reference ratios according to a variety of methods, for the purpose of a threshold for extreme fold changes.
#' 
#' @param data A data.frame containing expression values. Default = pep
#' @param experiments.map The experiments map. Default = Exp.map
#' @param int.root Name root for expression columns. Default = pep.ref[length(pep.ref)]
#' @param rat.root Name root for ratio columns. Default = pep.ratios.ref[1]
#' @param rat.con.grps Which sample groups are we considering. Default = RatConGrps
#' @param mode One of 1 or 2
#' @param parameters The experiment's parameters file. Default = Param
#' @param logInt Set to 0 or FALSE if input peptide intensities are linear scale. If the data is already log scale, set to the relevant scale's base. Default = FALSE
#' @param log.Ratios Set to 0 or FALSE if input peptide ratios are linear scale. If the data is already log scale, set to the relevant scale's base. Default = 2
#'
#' @examples
#' pep.Ref.Ratios %<o% make_RefRat()
#' pep[, colnames(pep.Ref.Ratios)] <- pep.Ref.Ratios              
#' 
#' @export

make_RefRat <- function(data = pep,
                        experiment.map = Exp.map,
                        int.root = pep.ref[length(pep.ref)],
                        rat.root = pep.ratios.ref[1],
                        rat.con.grps = RatConGrps,
                        mode = RefRat_Mode,
                        parameters = Param,
                        logInt = FALSE,
                        logRat = 2) {
  #proteoCraft::DefArg(proteoCraft::make_RefRat)
  mode <- as.character(mode)
  stopifnot(mode %in% c("1", "2"))
  isLog <- TRUE
  if ((!is.numeric(logInt))||(logInt <= 0)) {
    isLog <- FALSE
  }
  if ((!is.numeric(logRat))||(logRat <= 0)) {
    warning("Invalid \"logRat\" argument, must be valid log base!")
    logRat <- 2
  }
  rGrp <- proteoCraft::parse.Param.aggreg(parameters$Ratios.Groups)
  rrGrp <- proteoCraft::parse.Param.aggreg(parameters$Ratios.Ref.Groups)
  smplsGrp <- proteoCraft::parse.Param.aggreg(parameters$Volcano.plots.Aggregate.Level)
  #
  if (mode == "1") {
    # Here we are measuring intrinsic variability by comparing intra-sample group variation for control/reference sample group(s) only
    RES <- setNames(lapply(rGrp$values, function(x) { #x <- rGrp$values[1]
      if (rat.con.grps == "Ratio groups") {
        x1 <- unique(experiment.map[which(experiment.map[[rGrp$column]] == x), rrGrp$column])
      }
      if (rat.con.grps == "Experiments") {
        x1 <- unique(experiment.map$Experiment[which(experiment.map[[rGrp$column]] == x)])
        x1 <- unique(experiment.map[which(experiment.map$Experiment == x1), rrGrp$column])
      }
      if (rat.con.grps == "Whole dataset") {
        x1 <- unique(experiment.map[[rrGrp$column]])
      }
      xr1 <- paste0(int.root, x1, ".REF")
      xr1 <- xr1[which(xr1 %in% colnames(data))]
      if (length(xr1) > 1) {
        perm <- gtools::permutations(length(xr1), 2, xr1)
        if (isLog) {
          tmp <- apply(perm, 1, function(y) { (data[[y[[1]]]]-data[[y[[2]]]])/log(logRat, logInt) })
        } else {
          tmp <- apply(perm, 1, function(y) { log(data[[y[[1]]]]/data[[y[[2]]]], logRat) })
        }
      } else {
        x2 <- unique(experiment.map$Ref.Sample.Aggregate[which((experiment.map[[rGrp$column]] == x)&(experiment.map$Reference))])
        xr2 <- paste0(int.root, x2)
        if (isLog) {
          tmp <- sapply(xr2, function(y) { (data[[y]]-data[[xr1]])/log(logRat, logInt) })
        } else {
          tmp <- sapply(xr2, function(y) { log(data[[y]]/data[[xr1]], logRat) })
        }
      }
      tmp <- as.data.frame(tmp)
      colnames(tmp) <- paste0(rat.root, x, "_REF.to.REF_", 1:ncol(tmp))
      return(tmp)
    }), rGrp$values)
  }
  if (mode == "2") {
    # Here we are measuring intrinsic variability by comparing intra-sample group variation (replicate to replicate of the same stuff) for all sample groups,
    # not just control/reference sample group(s)!
    RES <- setNames(lapply(rGrp$values, function(rtGrp) { #rtGrp <- rGrp$values[1]
      if (rat.con.grps == "Ratio groups") {
        em <- experiment.map[which(experiment.map[[rGrp$column]] == rtGrp),]
      }
      if (rat.con.grps == "Experiments") {
        xp <- unique(experiment.map$Experiment[which(experiment.map[[rGrp$column]] == rtGrp)])
        em <- experiment.map[which(experiment.map$Experiment == xp),]
      }
      if (rat.con.grps == "Whole dataset") {
        em <- experiment.map
      }
      grps <- unique(em[[smplsGrp$column]])
      x <- setNames(lapply(grps, function(grp) {  #grp <- grps[1]
        y <- unique(em$Ref.Sample.Aggregate[which(em[[smplsGrp$column]] == grp)])
        y <- paste0(int.root, y)
        y <- y[which(y %in% colnames(data))]
        if (length(y)) {
          y <- gtools::permutations(length(y), 2, y) # Important to use permutations and not combinations since we want symmetry!
          y <- apply(y, 1, function(z) {
            if (isLog) {
              z <- (data[[z[[1]]]]-data[[z[[2]]]])/log(logRat, logInt)
            } else {
              z <- suppressWarnings(log(data[[z[[1]]]]/data[[z[[2]]]], logRat))
            }
            return(z)
          })
          y <- list(Outcome = TRUE, Res = y)
        } else { y <- list(Outcome = FALSE) }
        return(y)
      }), grps)
      x <- x[which(sapply(x, function(y) { y$Outcome }))]
      if (length(x)) {
        x <- lapply(x, function(y) { y$Res })
        x <- as.data.frame(x)
        colnames(x) <- paste0(rat.root, rtGrp, "_REF.to.REF_", 1:ncol(x))
        x <- list(Outcome = TRUE, Res = x)
      } else { x <- list(Outcome = FALSE) }
      return(x)
    }), rGrp$values)
    w <- which(sapply(RES, function(x) { x$Outcome }))
    if (length(w)) {
      RES <- lapply(RES[w], function(x) { x$Res })
    }
  }
  L <- length(RES)
  if (L) {
    names(RES) <- NULL
    RES <- data.frame(RES, check.names = FALSE)
  } else { stop("No Ref to Ref columns were generated, investigate!!!") }
  return(RES)
}
