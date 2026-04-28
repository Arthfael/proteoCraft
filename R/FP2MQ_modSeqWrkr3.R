#' .FP2MQ_modSeqWrkr3
#' 
#' Worker function used by FP_to_MQ().
#'
#' @param x Temp data, see FP_to_MQ code...
#' @param mods Modifications table.
#' 
#' @export

.FP2MQ_modSeqWrkr3 <- function(x,
                               mods) {
  #x <- strsplit(EV$"Modified sequence_verbose"[wMdSq2[1]], "\\(|\\)")
  x <- unlist(x)
  l <- length(x)
  wmds <- (1L:((l - 1L)/2L)) * 2L
  mds <- strsplit(x[wmds], ",")
  w1 <- grep("^[0-9]+ ", mds, invert = TRUE)
  w2 <- grep("^[0-9]+ ", mds)
  mds2 <- list()
  if (length(w1)) { mds2[w1] <- mds[w1] }
  if (length(w2)) {
    mds2[w2] <- lapply(mds[w2], function(y) {
      y <- unlist(strsplit(y, " "))
      rep(y[2L], as.integer(y[1L]))
    })
  }
  mds <- paste0("(", vapply(mds2, function(y) {
    y <- Modifs$Mark[match(y, Modifs$"Full name")]
    w <- which(y == "!m")
    if (length(w) > 1L) { y <- y[c(which(y != "!m"), w[1L])] }
    return(paste(y, collapse = ","))
  }, ""), ")")
  x[wmds] <- mds
  return(paste(x, collapse = ""))
}
