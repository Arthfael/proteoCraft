#' make_Rat2
#'
#' @description
#' Internal function to output consistent ratios and intensity values.
#' This partially replaces the now deprecated make_Rat() function.
#' 
#' @param myData Input data.
#' @param contrasts Contrasts data.frame for which to calculate ratios. It should contain columns of the form "[A-D]_samples" mapping samples to contrasts.
#' @param refGroups Factors aggregate defining blocks.
#' @param int.log Intensities log base, default = 10. Set to FALSE if intensities should not be log-transformed.
#' @param rat.log Ratios log base, default = 2
#' @param experiment.map Entire experiment map (unfiltered).
#' @param experiment.map_col Name of the sample column in the experiment map.
#' @param int.root Root of input intensity column names.
#' @param rat.root Root of output ratios column names.
#' 
#' @returns
#' Data supplemented with new ratios and intensities columns
#' 
#' @export

make_Rat2 <- function(myData = pep,
                      contrasts = myContrasts,
                      refGroups = RRG,
                      int.log = 10,
                      rat.log = 2,
                      experiment.map,
                      experiment.map_col = "Ref.Sample.Aggregate",
                      int.root,
                      rat.root) {
  #DefArg(make_Rat2)
  #DefArg(make_Rat2, silent = FALSE)
  #myData = pep; contrasts = myContrasts; refGroups = RRG; int.log = FALSE; rat.log = 2; experiment.map = Exp.map; int.root = pep.ref[length(pep.ref)]; rat.root = pep.ratios.ref
  #
  if (is.logical(int.log)) {
    if (!is.na(int.log)&&int.log) {
      logTrans <- TRUE
      warning("Assuming default log10 intensities.")
      int.log <- 10
    } else {
      int.log <- FALSE
      logTrans <- FALSE
    }
  }
  if (is.numeric(int.log)) {
    stopifnot(int.log > 0)
    logTrans <- TRUE
  }
  stopifnot(is.numeric(rat.log),
            rat.log > 0)
  if (logTrans) {
    log_Int2Rat <- log(rat.log, int.log)
  }
  #
  if (length(unique(experiment.map$Experiment)) == 1) {
    smplGroups$column <- sub("Exp", "", smplGroups$column)
  }
  ABkol <- paste0(c("A", "B"), "_samples")
  CDkol <- paste0(c("C", "D"), "_samples")
  CDkol <- intersect(CDkol, colnames(contrasts))
  CDkol_tst <- length(CDkol) == 2
  myData <- myData[, paste0(int.root, unique(unlist(contrasts[, c(ABkol, CDkol)])))]
  if (!logTrans) { myData <- log(myData, rat.log) } else {
    myData <- myData/log_Int2Rat
  }
  myRats <- lapply(1:nrow(contrasts), function(i) { #i <- 1 #i <- 3
    rat <- setNames(lapply(refGroups$values, function(grp) { #grp <- refGroups$values[1]
      em <- experiment.map[which(experiment.map[[refGroups$column]] == grp),]
      rs <- setNames(lapply(ABkol, function(x) { #x <- ABkol[1]
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
      second <- setNames(lapply(refGroups$values, function(grp) { #grp <- refGroups$values[1]
        em <- experiment.map[which(experiment.map[[refGroups$column]] == grp),]
        rs <- setNames(lapply(CDkol, function(x) { #x <- CDkol[1]
          x <- contrasts[[x]][[i]]
          x <- intersect(x, em[[experiment.map_col]])
          if (!length(x)) { return() }
          kol <- paste0(int.root, x)
          dat <- myData[, kol, drop = FALSE]
          if (logTrans) { dat <- log(dat, rat.log) } else {
            dat <- dat/log_Int2Rat
          }
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
