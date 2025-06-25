#' Data_Impute2
#'
#' @description
#' Better function to impute missing values than my former version, "Data_Impute". Uses the imputeLCMD package.
#' Works in two steps:
#'  - Step 1: within user-defined groups of user-defined related samples (replicates from a treatment),
#'            on which it performs MAR imputation using KNN.
#'  - Step 2: globally on all remaining missing values, which it imputed as MNAR using QRILC.
#' The output is a list of the two arrays: new quantitative values and positions which were imputed.
#' 
#' @param quant_data A matrix or data.frame with more than 1 column, containing missing or invalid (e.g. infinite).
#' @param groups A vector of groups of length equal to the number of columns of the input data frame.
#' @param is.log Default = TRUE. Are the input values log transformed? If set to FALSE, the data will be log-transformed prior to processing, then de-logged.
#' @param seed Use a fixed seed for reproducibility purposes.
#' 
#' @examples
#' New_Data <- Data_Impute2(data)
#' 
#' @export


Data_Impute2 <- function(quant_data,
                         groups,
                         is.log = TRUE,
                         seed) {
  #proteoCraft::DefArg(proteoCraft::Data_Impute2)
  #quant_data <- temp; groups = grps
  if (missing(seed)) { if (exists("mySeed")) { seed <- mySeed } else { seed <- 1234567 } }
  set.seed(seed)
  norm::rngseed(seed)
  #
  kol <- colnames(quant_data)
  # Get class
  Klass <- class(quant_data)
  if (!sum(Klass %in% c("data.frame", "matrix"))) { stop("Input \"quant_data\" must be a matrix or data.frame!") }
  if (!"matrix" %in% Klass) { quant_data <- as.matrix(quant_data) }
  if (missing(groups)) { groups <- rep("Group 1", length(kol)) }
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
    # Step one: work within groups to impute MAR/MCAR
    for (grp in unique(groups)) { #grp <- unique(groups)[1]
      # Groups, where we have at least one non-missing value: MAR
      m <- which(groups == grp)
      rwSms <- rowSums(is.na(imputed_data[, m, drop = FALSE]))
      #aggregate(rwSms, list(rwSms), length)
      w <- which() < length(m))
      if (length(w)) {
        ms <- imputeLCMD::model.Selector(imputed_data[w, m, drop = FALSE])
        imputed_data[w, m] <- imputeLCMD::impute.MAR(imputed_data[w, m, drop = FALSE], ms, "KNN")
      }
    }
    # Step 2: apply MNAR correction to remaining misses within imputed data
    w <- which(is.na(imputed_data), arr.ind = TRUE)
    if (nrow(w)) { imputed_data <- imputeLCMD::impute.QRILC(imputed_data)[[1]] }
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
