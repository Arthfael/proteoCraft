#' .FP2MQ_modSeqWrkr4
#' 
#' Worker function used by FP_to_MQ().
#'
#' @param x Temp data, see FP_to_MQ code...
#' @param mods Modifications table.

.FP2MQ_modSeqWrkr4 <- function(x,
                               mods) {
  #x <- tmp[1,]
  x1 <- unlist(x[[1L]])
  x2 <- unlist(x[[2L]])
  x1 <- x1[which(x1 != "!m")]
  if (length(x1)) {
    mds <- stats::aggregate(x1, list(x1), length)
    mds$Group.1 <- mods$"Full name"[match(mds$Group.1, mods$Mark)]
  }
  if (length(x2)) {
    x2 <- stats::aggregate(x2, list(x2), length)
    mds <- { if (length(x1)) { rbind(mds, x2) } else { x2 } }
  }
  mds <- mds[order(mds$Group.1, decreasing = FALSE),]
  return(paste(apply(mds, 1L, function(y) {
    gsub("^1 ", "", paste(rev(y), collapse = " "))
  }), collapse = ","))
}
