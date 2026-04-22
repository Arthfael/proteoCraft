#' Av_SE_fun
#' 
#' @description 
#' A simple function to calculate the arithmetic mean and standard error of a numeric vector.
#' 
#' @param vect Input numeric vector.
#' 
#' @returns
#' Returns the arithmetic mean and standard error of a numeric vector.
#' 
#' @export

Av_SE_fun <- function(vect) {
  res <- is.all.good(as.numeric(vect))
  l <- length(res)
  if (l) { res <- c(mean(res), sd(res)/sqrt(l)) } else {
    res <- unique(vect[which(!is.na(vect))])
    res <- if (l == 1L) { c(res, NA) } else { c(NA, NA) }
  }
  return(res)
}
