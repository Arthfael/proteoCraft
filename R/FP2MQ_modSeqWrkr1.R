#' .FP2MQ_modSeqWrkr1
#' 
#' Worker function used by FP_to_MQ().
#'
#' @param x List, comma-split assigned modifications column from FragPipe.
#' @param mods Modifications table.
#' @param pat Regex for amino acids.

.FP2MQ_modSeqWrkr1 <- function(x,
                               mods,
                               pat) {
  #x <- strsplit(PSMs$`Assigned Modifications`[wMdSq[1]], ", ?")
  x <- gsub("\\)$", "", unlist(x))
  ptms <- as.data.frame(t(sapply(x, function(y) { unlist(strsplit(y, split = "\\(")) })))
  colnames(ptms) <- c("Site", "Mass shift")
  rownames(ptms) <- NULL
  ptms[, c("Site", "Position")] <- t(sapply(ptms$Site, function(aa) {
    if (aa == "N-term") { ps <- 0L }
    if (aa == "C-term") { ps <- l+1L }
    if (!aa %in% c("N-term", "C-term")) {
      ps <- gsub(paste0(pat, "$"), "", aa)
      aa <- gsub("^[0-9]+", "", aa)
    }
    return(c(aa, ps))
  }))
  ptms$"Mass shift" <- as.numeric(ptms$"Mass shift")
  ptms$Position <- as.numeric(ptms$Position)
  ptms$"Full name" <- NA
  w <- which(is.na(ptms$"Full name"))
  k <- 5L
  while ((length(w))&&(k > 2L)) {
    ptms$"Full name"[w] <- mods$"Full name"[match(round(ptms$`Mass shift`[w], k),
                                                  round(mods$`Mass delta`, k))]
    w <- which(is.na(ptms$"Full name"))
    k <- k-1L
  }
  if (length(w)) {
    stop("I am an error and I mean that some assigned PTM was not recognized... but actually you should probably convert me to a warning!")
  }
  return(ptms)
}
