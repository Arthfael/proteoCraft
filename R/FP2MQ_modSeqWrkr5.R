#' .FP2MQ_modSeqWrkr5
#' 
#' Worker function used by FP_to_MQ().
#'
#' @param x Temp data, see FP_to_MQ code...
#' @param mods Modifications table.
#' 
#' @export

.FP2MQ_modSeqWrkr5 <- function(x,
                               mods) {
  #x <- tmp[1L]
  x <- unlist(x)
  x <- stats::aggregate(x, list(x), length)
  x <- x[order(x$Group.1, decreasing = FALSE),]
  x$Group.1 <- mods$"Full name"[match(x$Group.1, mods$Mark)]
  return(paste(apply(x, 1L, function(y) {
    gsub("^1 ", "", paste(rev(y), collapse = " "))
  }), collapse = ","))
}
