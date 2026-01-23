#' Data_Impute2
#'
#' @description
#' A function to impute missing values. Uses the imputeLCMD package.
#' 
#' @param quant_data A matrix or data.frame with more than 1 column, containing missing or invalid (e.g. infinite).
#' @param groups A vector of groups of length equal to the number of columns of the input data frame.
#' @param is.log Default = TRUE. Are the input values log transformed? If set to FALSE, the data will be log-transformed prior to processing, then de-logged.
#' @param seed Use a fixed seed for reproducibility purposes.
#' @param min_misses Minimum number of misses for assuming MAR - this is a way to bypass imputeLCMD's model.Selector(), whose predictions as to which rows are MAR (as opposed to left censored MNAR) seem sometimes rather optimistic. The current default is 1 (i.e. no effect of this parameter), but where the data will be subjected to statistical analysis it may be a good idea to increase this to 2. 
#' 
#' @details
#' Works in two steps:
#'  - Step 1: within user-defined groups of user-defined related samples (replicates from a treatment),
#'            on which it performs MAR/MCAR imputation using KNN.
#'  - Step 2: globally on all remaining missing values, which it imputes as MNAR using QRILC.
#'
#' @returns
#' A list of two arrays:
#'  - "Imputed_data" = new quantitative values
#'  - "Positions_Imputed" = positions of imputed values
#' 
#' @examples
#' New_Data <- Data_Impute2(data)
#' 
#' @export

Data_Impute2 <- function(quant_data,
                         groups,
                         is.log = TRUE,
                         seed,
                         min_misses = 1) {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::Data_Impute2);TESTING <- TRUE
  #quant_data <- temp; groups = grps
  #quant_data <- tmpDat; groups = Gr
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  if (misFun(seed)) { if (exists("mySeed")) { seed <- mySeed } else { seed <- 1234567 } }
  set.seed(seed)
  norm::rngseed(seed)
  #
  if ((!is.numeric(min_misses))||(is.na(min_misses))||(min_misses < 1)) { min_misses <- 1 }
  min_misses <- round(min_misses)
  #
  kol <- colnames(quant_data)
  # Get class
  Klass <- class(quant_data)
  if (!sum(Klass %in% c("data.frame", "matrix"))) { stop("Input \"quant_data\" must be a matrix or data.frame!") }
  if (!"matrix" %in% Klass) { quant_data <- as.matrix(quant_data) }
  if (misFun(groups)) { groups <- rep("Group 1", length(kol)) }
  imputed_data <- quant_data
  # If necessary log-transform
  if (!is.log) { imputed_data <- log10(imputed_data) }
  for (kk in kol) { imputed_data[which(!proteoCraft::is.all.good(imputed_data[,kk], 2)), kk] <- NA }
  wNA <- which(is.na(imputed_data), arr.ind = TRUE)
  if (!is.log) {
    wOK <- which(!is.na(imputed_data), arr.ind = TRUE)
  }
  misses <- matrix(rep(FALSE, nrow(imputed_data)*ncol(imputed_data)),
                   ncol = ncol(imputed_data))
  misses[wNA] <- TRUE
  if (nrow(wNA)) {
    message("Missing values found, correcting (method = knn for MAR, QRILC for MNAR)!")
    # Step one: work within groups to impute MAR/MAR
    for (grp in unique(groups)) { #grp <- unique(groups)[1] #grp <- unique(groups)[2]
      # Groups, where we have at least one non-missing value: MAR
      m <- which(groups == grp)
      NAs_in_Grp <- rowSums(is.na(imputed_data[, m, drop = FALSE]))
      #aggregate(NAs_in_Grp, list(NAs_in_Grp), length)
      #w <- which(NAs_in_Grp < length(m))
      #w1 <- which((NAs_in_Grp < length(m))&(NAs_in_Grp >= 1))
      #w2 <- which((NAs_in_Grp < length(m))&(NAs_in_Grp >= 2))
      w <- which((NAs_in_Grp < length(m))&(NAs_in_Grp >= min_misses))
      if (length(w)) {
        ms <- imputeLCMD::model.Selector(imputed_data[w, m, drop = FALSE])
        a <- capture.output({ # Using capture.output to suppress the messages of impute.MAR which are really painful!
          tst <- try(imputeLCMD::impute.MAR(imputed_data[w, m, drop = FALSE], ms, "KNN"), silent = TRUE)
        })
        if (!"try-error" %in% class(tst)) {
          imputed_data[w, m] <- tst 
        } else {
          warning(tst[1])
        }
      }
    }
    # Step 2: apply MNAR correction to remaining misses within imputed data
    w <- which(is.na(imputed_data), arr.ind = TRUE)
    if (nrow(w)) {
      tst <- capture.output({
        imputed_data <- imputeLCMD::impute.QRILC(imputed_data)[[1]]
      })
    }
  }
  # De-log if necessary
  if (!is.log) {
    imputed_data <- 10^(imputed_data)
    # In case the log/delog cycle would have introduced small changes:
    imputed_data[wOK] <- quant_data[wOK]
  }
  #
  if (!sum(Klass %in% class(imputed_data))) { imputed_data <- get(paste0("as.", Klass))(imputed_data) }
  res <- list("Imputed_data" = imputed_data,
              "Positions_Imputed" = misses)
  return(res)
}
