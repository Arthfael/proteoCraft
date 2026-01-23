#' diffLog
#'
#' @description 
#' Internal function whose sole purpose is to be fed to another whose sum will be minimized by minpack.lm::nls.lm() within the LFQ.lm() function.
#' 
#' @param p Starting estimates, and N-1 length numeric vector.
#' @param dat Data.frame of log-transformed, row-normalized quantitative profiles with N columns.
#' 
#' @examples
#' temp3 <- ... # Data frame of log-transformed, row-normalized quantitative values
#' diffLog_v <- function(...) {
#'   p <- list(...)
#'   res <- diffLog(p, dat = temp3)
#'   return(res)
#' }
#' 
#' @export

diffLog <- function(p, 
                     dat) {
  p <- c(0, unlist(p))
  dat <- sweep(dat, 1, p, "-")
  dat <- as.matrix(dat)
  nr <- nrow(dat)
  res <- apply(dat, 2, function(z) {
    z <- as.numeric(z)
    return(
      lapply(1:(nr-1), function(n) { # creates all pairwise differences - the order does not matter since the sum square of the vector will then be calculated (and minimized)
        z[n]-z[(n+1):nr]
      })
    )
  })
  res <- proteoCraft::is.all.good(unlist(res))
  return(res)
}
