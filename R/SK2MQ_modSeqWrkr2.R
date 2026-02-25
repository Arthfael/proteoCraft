#' .SK2MQ_modSeqWrkr2
#' 
#' Worker function used by Skyline_to_MQ().
#'
#' @param x Temporary modified sequence

.SK2MQ_modSeqWrkr2 <- function(x) {
  if (!length(x)) { return("Unmodified") }
  x <- aggregate(x, list(x), length)
  x <- x[order(x$Group.1, decreasing = FALSE), 2:1]
  x$x <- gsub("1 ", "", paste0(as.character(x$x), " "))
  x <- paste0(apply(x, 1, paste, collapse = ""), collapse = ",")
  return(x)
}
