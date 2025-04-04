#' log_ratio_av
#'
#' @description 
#' Averages log-transformed ratio values.
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
