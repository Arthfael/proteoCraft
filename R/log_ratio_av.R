#' log_ratio_av
#'
#' @description 
#' Averages log-transformed ratio values.\cr
#'  - The mean is used when some valid values are present\cr
#'  - If only infinite values are present, if these are all the same (-Inf or +Inf, not both), then this consensus value is returned, otherwise NA is returned\cr
#'  - All NA or NaN values return NA\cr
#'
#' @param x Log-transformed ratio values to be averaged.
#' 
#' @export

log_ratio_av <- function(x) {
  x1 <- proteoCraft::is.all.good(as.numeric(unlist(x)))
  if (length(x1)) { x1 <- mean(x1) } else {
    x1 <- unique(x[which(!is.na(x))])
    if (length(x1) != 1) { x1 <- NA }
  }
  return(x1)
}
