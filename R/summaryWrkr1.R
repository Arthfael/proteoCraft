#' .summaryWrkr1
#' 
#' Worker function used by MQ_summary().
#'
#' @param x File.
#' @param raw Table of MS raw files.
#' @param type Type of chromatogram to export.
#' 
#' @export

.summaryWrkr1 <- function(x,
                          raw,
                          type) {
  x2 <- rawrr::readChromatogram(raw$Path[x],
                                type = type)
  return(data.frame("Raw file" = raw$`Raw file`[x],
                    "Raw file path" = raw$Path[x],
                    "Retention time" = as.numeric(x2$times),
                    "Intensity" = as.numeric(x2$intensities),
                    check.names = FALSE))
}
