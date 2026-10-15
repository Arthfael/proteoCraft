#' .labelEdit_worker
#' 
#' Worker function within Volcano.plot() for editing labels to make them multi-row.
#'
#' @param x Character label split at spaces
#' @param colchar How many characters wide a column of text may be. Default = 25
#' 
#' @export

.labelEdit_worker <- function(x,
                              colchar = 25L) {
  #x <- strsplit(Prot$Labels[weech], "  ?")[1L]
  x <- unlist(x)
  if (length(x) > 8L) { x <- c(x[1L:8L], "...") } # This is because:
  # a) sentences of too many words can cause problems with subsequent combinatorial steps, slowing down the script or even causing it to fail
  # and
  # b) do you really expect a protein label more than 8 words long to be helpful?
  # Usually this will not be an issue with Uniprot but can be for other databases, e.g. TAIR, where protein names contain additional information for unknown proteins.
  nc <- min(c(ceiling((sum(nchar(x)) + length(x) - 1L)/colchar),
              length(x)))
  if (nc > 1L) {
    tstbrk <- cbind(0L, gtools::combinations(length(x)-1L, nc-1L), length(x))
    tstbrk <- apply(tstbrk, 1L, \(y) {
      vapply(1L:nc, \(z) { paste(x[(y[z]+1L):y[z+1L]], collapse = " ") }, "")
    })
    tstsd <- apply(tstbrk, 2L, \(y) { sum(nchar(y)^2L) })
    label <- paste(tstbrk[, which(tstsd == min(tstsd))[1L]], collapse = "\n")
    rm(tstbrk, tstsd)
  } else { label <- paste(x, collapse = " ") }
  return(label)
}
