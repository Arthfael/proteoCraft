#' diff.log
#'
#' @description 
#' A function whose sole purpose is to be turned into another whose sum will be minimized by nls.lm within the LFQ function.
#' 
#' @param p Starting estimates.
#' @param dat Data frame of log-transformed, row-normalized quantitative profiles.
#' 
#' @examples
#' temp3 <- ... # Data frame of log-transformed, row-normalized quantitative values
#' diff.log.v <- function(...) {
#'   p <- list(...)
#'   res <- diff.log(p, dat = temp3)
#'   return(res)
#' }
#' 
#' @export diff.log

diff.log <- function(p, dat) {
  p <- c(0, unlist(p))
  dat <- sweep(dat, 1, p, "-")
  dat <- as.matrix(dat)
  nr <- nrow(dat)
  res <- as.data.frame(t(apply(dat, 2, function(z) {
    z <- as.numeric(z)
    m <- sapply(1:(nr-1), function(n) {
      proteoCraft::is.all.good(z[n]-z[(n+1):nr])
    }) # creates all pairwise differences - the order does not matter since the sum square of the vector will then be calculated (and minimized)
    return(proteoCraft::is.all.good(unlist(m)))
  })))
  res <- unlist(res)
  return(res)
}
