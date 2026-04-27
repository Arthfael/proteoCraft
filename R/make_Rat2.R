#' make_Rat2
#'
#' @description
#' Internal function to output average log ratios.
#' This partially replaces the slowly deprecated make_Rat() function.
#' 
#' @param myData Input data.
#' @param contrasts Contrasts data.frame for which to calculate ratios. It should contain columns of the form "[A-D]_samples" mapping samples to contrasts.
#' @param refGroups Factors aggregate defining blocks.
#' @param int.log Intensities log base, default = 10. Set to FALSE if input intensities are not log-transformed.
#' @param rat.log Ratios log base, default = 2
#' @param experiment.map Entire experiment map (unfiltered).
#' @param experiment.map_col Name of the sample column in the experiment map.
#' @param int.root Root of input intensity column names.
#' @param rat.root Root of output ratios column names.
#' 
#' @returns
#' Data supplemented with new ratio columns. Unlike make_Rat, no intensities are output.
#' 
#' @export

make_Rat2 <- function(myData = pep,
                      contrasts = myContrasts,
                      refGroups = RRG,
                      int.log = 10L,
                      rat.log = 2L,
                      experiment.map,
                      experiment.map_col = "Ref.Sample.Aggregate",
                      int.root,
                      rat.root) {
  #DefArg(make_Rat2)
  #DefArg(make_Rat2, silent = FALSE)
  #myData = pep; contrasts = myContrasts; refGroups = RRG; int.log = FALSE; rat.log = 2L; experiment.map = Exp.map; int.root = pep.ref[length(pep.ref)]; rat.root = pep.ratios.ref
  #
  # Process int.log argument:
  #  - Input intensity data may be log-transformed or not
  #  - We will output data in the same scale.
  if ((!is.logical(int.log))&&(!is.numeric(int.log))) { stop("Invalid \"int.log\" argument!") }
  if (is.logical(int.log)) {
    if (!is.na(int.log)&&int.log) {
      logTrans <- TRUE
      warning("Assuming default log10 intensities.")
      int.log <- 10L
    } else {
      int.log <- FALSE
      logTrans <- FALSE
    }
  }
  if (is.numeric(int.log)) {
    if (int.log <= 0) { stop("Invalid \"int.log\" argument!") }
    logTrans <- TRUE
  }
  #
  # Process rat.log argument
  #  - Here we just want the desired log base.
  #  - We will always output log-transformed ratios.
  if ((!is.numeric(rat.log))||(rat.log <= 0)) { stop("Invalid \"rat.log\" argument!") }
  #
  if (logTrans) {
    if (rat.log == int.log) {
      cat("make_Ratios2: input data is log-transformed.\n")
      log_Int2Rat <- 1L
    } else {
      cat(paste0("make_Ratios2: input data is log-transformed and its base (", int.log, 
                 ") is different from that of the desired ratios (", rat.log, ")\n"))
      log_Int2Rat <- base::log(rat.log, int.log)
    }
  } else {
    cat("make_Ratios2: input data is not log-transformed.\n")
  }
  #
  ABkol <- c("A_samples", "B_samples")
  CDkol <- intersect(c("C_samples", "D_samples"), colnames(contrasts))
  CDkol_tst <- length(CDkol) == 2L
  myData <- myData[, paste0(int.root, unique(unlist(contrasts[, c(ABkol, CDkol)])))]
  #
  # If input was not log, we log-transform to ratios's log base, otherwise we adjust base.
  myData <- if (!logTrans) { base::log(myData, rat.log) } else { myData/log_Int2Rat }
  # -> myData's base is that of rat.log
  myRats <- lapply(1L:nrow(contrasts), \(i) { #i <- 1L #i <- 3L #i <- 5L
    rat <- setNames(lapply(refGroups$values, \(grp) { #grp <- refGroups$values[1L]
      em <- experiment.map[which(experiment.map[[refGroups$column]] == grp),]
      rs <- setNames(lapply(ABkol, \(x) { #x <- ABkol[1L]
        x <- contrasts[[x]][[i]]
        x <- intersect(x, em[[experiment.map_col]])
        if (!length(x)) { return() }
        kol <- paste0(int.root, x)
        dat <- myData[, kol, drop = FALSE]
        rowMeans(dat, na.rm = TRUE)
      }), c("A", "B"))
      rs <- as.data.frame(do.call(cbind, rs))
      return(rs$A - rs$B)
    }), refGroups$values)
    rat <- as.data.frame(do.call(cbind, rat))
    rat <- rowMeans(rat, na.rm = TRUE)
    if (CDkol_tst) {
      second <- setNames(lapply(refGroups$values, \(grp) { #grp <- refGroups$values[1L]
        em <- experiment.map[which(experiment.map[[refGroups$column]] == grp),]
        rs <- setNames(lapply(CDkol, \(x) { #x <- CDkol[1L]
          x <- contrasts[[x]][[i]]
          x <- intersect(x, em[[experiment.map_col]])
          if (!length(x)) { return() }
          kol <- paste0(int.root, x)
          dat <- myData[, kol, drop = FALSE]
          rowMeans(dat, na.rm = TRUE)
        }), c("C", "D"))
        rs <- as.data.frame(do.call(cbind, rs))
        return(rs$C - rs$D)
      }), refGroups$values)
      second <- as.data.frame(do.call(cbind, second))
      if (nrow(second)) {
        second <- rowMeans(second, na.rm = TRUE)
        rat <- rat-second
      }
    }
    return(rat)
  })
  myRats <- do.call(cbind, myRats)
  colnames(myRats) <- paste0(rat.root, contrasts$Contrast)
  return(myRats)
}
