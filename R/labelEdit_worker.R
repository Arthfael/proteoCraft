#' .labelEdit_worker
#' 
#' Worker function within Volcano.plot() for editing labels to make them multi-row.
#'
#' @param x Character label split at spaces
#' @param colchar How many characters wide a column of text may be. Default = 25

.labelEdit_worker <- function(x,
                              colchar = 25) {
  #x <- strsplit(Prot$Labels[weech], "  ?")[1]
  x <- unlist(x)
  if (length(x) > 8) { x <- c(x[1:8], "...") } # This is because:
  # a) sentences of too many words can cause problems with subsequent combinatorial steps, slowing down the script or even causing it to fail
  # and
  # b) do you really expect a protein label more than 8 words long to be helpful?
  # Usually this will not be an issue with Uniprot but can be for other databases, e.g. TAIR, where protein names contain additional information for unknown proteins.
  nc <- min(c(ceiling((sum(nchar(x)) + length(x) - 1)/colchar),
              length(x)))
  if (nc > 1) {
    tstbrk <- cbind(0, gtools::combinations(length(x)-1, nc-1), length(x))
    tstbrk <- apply(tstbrk, 1, function(y) {
      vapply(1:nc, function(z) { paste(x[(y[z]+1):y[z+1]], collapse = " ") }, "")
    })
    tstsd <- apply(tstbrk, 2, function(y) { sum(nchar(y)^2) })
    label <- paste(tstbrk[, which(tstsd == min(tstsd))[1]], collapse = "\n")
    rm(tstbrk, tstsd)
  } else { label <- paste(x, collapse = " ") }
  return(label)
}
