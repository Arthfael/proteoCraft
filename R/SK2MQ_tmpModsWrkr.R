#' .SK2MQ_tmpModsWrkr
#' 
#' Worker function used by Skyline_to_MQ().
#'
#' @param x 
#' 
#' @export

.SK2MQ_tmpModsWrkr <- function(x) {
  #x <- tmp$ModSeq[1]
  x <- annot_to_tabl(x)[[1]]
  w <- which(x$Annotations != "")
  x$Sequence[2] <- paste0("_", x$Sequence[2])
  do.call(paste, c(x[w, , drop = FALSE], sep = "mod_"))
}
