#' .SK2MQ_tmpModsWrkr
#' 
#' @description
#' Worker function used by Skyline_to_MQ().
#'
#' @param x Modified sequence.

.SK2MQ_tmpModsWrkr <- function(x) {
  #x <- tmp$ModSeq[1L]
  x <- annot_to_tabl(x)[[1L]]
  w <- which(x$Annotations != "")
  x$Sequence[2L] <- paste0("_", x$Sequence[2L])
  do.call(paste, c(x[w, , drop = FALSE], sep = "mod_"))
}
