#' geomean
#'
#' @description 
#' Calculates the geometric mean of linear-scale data.
#' NB: I don't think that this function is being used anymore, is it?
#' 
#' @param x The numeric data which will be averaged.
#' @param filter If set to TRUE (default), this will remove all non numeric, NA, NaN and +/-Inf values.
#' 
#' @examples
#' > geomean(c(1, 3))
#'   1.732051
#'   
#' @export

geomean <- function(x, filter = TRUE) {
  exp(mean(proteoCraft::is.all.good(log(unlist(x)))))
}
