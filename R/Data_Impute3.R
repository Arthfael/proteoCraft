#' Data_Impute3
#'
#' @description
#' Work in Progress for a better wrapper function to impute missing values than my former versions, "Data_Impute" and "Data_Impute2" (the latter is the one used as of 29/06/2021).
#' Distinguishes between Missing-At-Random (MAR) and Missing-Not-At-Random, using knn for the former, QRILC for the latter.
#' The output is a list of the two arrays: new quantitative values and positions which were imputed.
#' 
#' @param quant_data A matrix or data.frame with more than 1 column, containing missing or invalid (e.g. infinite).
#' @param groups A vector of groups of length equal to the number of columns of the input data frame.
#' @param is.log Default = TRUE. Are the input values log transformed? If set to FALSE, the data will be log-transformed prior to processing, then de-logged.
#' @param seed Use a fixed seed for reproducibility purposes.
#' 
#' @examples
#' New_Data <- Data_Impute3(data)
#' 
#' @export


Data_Impute3 <- function(quant_data,
                         groups,
                         is.log = TRUE,
                         seed) {
  #proteoCraft::DefArg(proteoCraft::Data_Impute3); quant_data <- PG[, paste0(Prot.Expr.Root, Ref.Sample.Aggregate$values)]; groups = Exp.map[match(Ref.Sample.Aggregate$values, Exp.map$Ref.Sample.Aggregate), Volcano.plots.Aggregate.Level$column]
  if (missing(seed)) { if (exists("mySeed")) { seed <- mySeed } else { seed <- 1234567 } }
  set.seed(seed)
  norm::rngseed(seed)
  #
  Klass <- class(quant_data)
  if (!Klass %in% c("data.frame", "matrix")) { stop("Input \"quant_data\" must be a matrix or data.frame!") }
  if (!Klass != "matrix") { quant_data <- as.matrix(quant_data) }
  kol <- colnames(quant_data)
  imputed_data <- quant_data
  for (kk in kol) {
    if (!is.log) { imputed_data[[kk]] <- log10(imputed_data[[kk]]) }
    imputed_data[which(!proteoCraft::is.all.good(imputed_data[[kk]], 2)), kk] <- NA
  }
  W <- which(is.na(imputed_data), arr.ind = TRUE)
  misses <- as.data.frame(matrix(rep(FALSE, nrow(imputed_data)*ncol(imputed_data)), ncol = ncol(imputed_data)))
  misses[W] <- TRUE
  res <- list("Positions_Imputed" = misses)
  if (nrow(W)) {
    message("Missing values found, correcting (method = knn for MAR, QRILC for MNAR)!")
    for (grp in unique(groups)) { #grp <- unique(groups)[1]
      # Groups, where we have at least one missing value: MAR
      m <- which(groups == grp)
      w <- which(apply(imputed_data[, m], 1, function(x) { sum(is.na(x)) }) < length(m))
      ms <- imputeLCMD::model.Selector(imputed_data[w, m])
      tst <- imputeLCMD::impute.MAR(imputed_data[w, m], ms, method = "KNN")
      imputed_data[w, m] <- tst[[1]]
    }
    # View(quant_data)
    # View(imputed_data)
    W <- which(is.na(imputed_data), arr.ind = TRUE)
    if (nrow(W)) {
      # The rest we assume is MNAR
      imputed_data2 <- imputeLCMD::impute.QRILC(imputed_data)
    }
    # View(quant_data)
    # View(imputed_data2)
  } else { warning("No missing data!") }
  if (!is.log) { for (kk in kol) { imputed_data[[kk]] <- 10^(imputed_data[[kk]]) } }
  if (class(imputed_data) != Klass) { imputed_data <- get(paste0("as.", Klass))(imputed_data) }
  res[["Imputed_data"]] <- imputed_data
  return(res)
}
