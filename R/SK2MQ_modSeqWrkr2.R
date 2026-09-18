#' .SK2MQ_modSeqWrkr2
#' 
#' Worker function used by Skyline_to_MQ().
#'
#' @param x Temporary modified sequence
#' 
#' @export

.SK2MQ_modSeqWrkr2 <- function(x) {
  if (!length(x)) { return("Unmodified") }
  x <- stats::aggregate(x, list(x), length)
  x <- x[order(x$Group.1, decreasing = FALSE), 2L:1L]
  x$x <- gsub("1 ", "", paste0(as.character(x$x), " "))
  x <- paste0(apply(x, 1L, paste, collapse = ""), collapse = ",")
  return(x)
}
