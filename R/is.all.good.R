#' is.all.good
#'
#' @description 
#' A function to filter numeric vectors to just include real numbers.
#' 
#' @param x Input numeric or integer vector, or all-numeric/integer matrix.
#' @param mode The function can work in two ways; if set to "values", "1" or 1 (default), removes all +Inf, -Inf, NA and NaN values. If set to "logical", "2" or 2, returns TRUE for real finite values and FALSE for others. If input is non-numeric, will throw a warning and attempt to convert to numeric.
#' 
#' @examples
#' is.all.good(c(1, NA, NaN, 2, -Inf, +Inf))
#' # [1] 1 2
#' is.all.good(c(1, NA, NaN, 2, -Inf, +Inf), mode = "logical")
#' # [1]  TRUE FALSE FALSE  TRUE FALSE FALSE
#' @export

is.all.good <- function(x, mode = "values") {
  mode <- c(0, 0, 1, 1)[match(mode, c("1", "values", "2", "logical"))]
  stopifnot(!is.na(mode))
  klass <- class(x)
  if (!sum(c("integer", "integer64", "numeric", "matrix") %in% klass)) {
    x <- as.numeric(unlist(x))
    warning("I had to convert data to numerics, you may want to check your input!")
  }
  res <- (!is.na(x))&(!is.nan(x))&(is.finite(x))
  if (!mode) {
    if ("matrix" %in% klass) {
      warning("Converting matrix data to numerics to filter it...")
      x <- as.numeric(x)
    }
    res <- x[which(res)]
  }
  return(res)
}
