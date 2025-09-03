#' logical.op
#'
#' @description 
#' A function to call logical operators by name.
#' 
#' @param name Name of the operator.
#' @param x Left-hand side term.
#' @param y Right-hand side term.
#' @param ties Logical: are ties accepted?
#' 
#' @examples
#' logical.op("equals", x = 1, y = 2)
#' 
#' @export

logical.op <- function(x, test, y, ties = TRUE) {
  lX <- length(x)
  lY <- length(y)
  if (!lX) { stop("x must be a vector of non-null length!") }
  if (!lY) { stop("y must be a vector of non-null length!") }
  if (lX != lY) { stop("x and y must be the same length!") }
  test <- toupper(gsub(" |_|-|\\.", "", test))
  EQUAL <- c("EQUAL", "EQUALS", "EQUALTO", "ISEQUALTO", "=", "==")
  DIFFERENT <- c("DIFFERENT", "DIFFERENTFROM", "ISDIFFERENT", "ISDIFFERENTFROM", "!=")
  SUPERIOR <- c("SUPERIOR", "SUPERIOUR", "ISSUPERIOR", "ISSUPERIOUR", "ISSUPERIORTO", "ISSUPERIOURTO", ">", ">=")
  INFERIOR <- c("INFERIOR", "INFERIOUR", "ISINFERIOR", "ISINFERIOUR", "ISINFERIORTO", "ISINFERIOURTO", "<", "<=")
  stopifnot(test %in% c(EQUAL, DIFFERENT, SUPERIOR, INFERIOR))
  if (test %in% EQUAL) { return(x == y) }
  if (test %in% DIFFERENT) { return(x != y) }
  if (test %in% SUPERIOR) {
    if ((ties)||(test == ">=")) { return(x >= y) } else { return(x > y) }
  }
  if (test %in% INFERIOR) {
    if ((ties)||(test == "<=")) { return(x <= y) } else { return(x < y) }
  }
}
