#' .summaryWrkr2
#' 
#' Worker function used by MQ_summary().
#'
#' @param x Raw file
#' @param tempDat Temp data.
#' @param size Type of chromatogram to export.

.summaryWrkr2 <- function(x,
                          tempDat,
                          size) { #x <- rawFls[1L]
  w <- data.frame(Wh = which(tempDat$"Raw file path" == x))
  w$tst <- vapply(1L:length(w$Wh), \(x) {
    m <- max(tempDat$"Retention time"[w$Wh])
    rng <- c(max(c(0, tempDat$"Retention time"[w$Wh[x]]-size/2)),
             min(c(tempDat$"Retention time"[w$Wh[x]]+size/2, m)))
    w2 <- which((tempDat$"Retention time"[w$Wh] >= rng[1L])&(tempDat$"Retention time"[w$Wh] <= rng[2L]))
    return(tempDat$Intensity[w$Wh[x]] == max(tempDat$Intensity[w$Wh][w2]))
  }, TRUE)
  return(w)
}
