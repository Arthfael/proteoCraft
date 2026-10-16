#' .SK2MQ_modSeqWrkr1
#' 
#' @description
#' Worker function used by Skyline_to_MQ().
#'
#' @param x Temporary modified sequence
#' @param tempMods Temporary modifications table

.SK2MQ_modSeqWrkr1 <- function(x,
                               tempMods) {
  #x <- temp1[[1L]]
  x1 <- x2 <- unlist(x)
  l <- length(x1)
  rg <- c(1L:((l-1L)/2L))*2L
  nc <- nchar(x1[rg-1L])
  aa <- substr(x1[rg-1L], nc, nc)
  if (nc[1L] == 1L) { aa[1L] <- substr(x1[3L], 1L, 1L) }
  m <- match(x1[rg], tempMods$Old)
  x1[rg] <- paste0("(", tempMods$Mark[m], ")")
  x2[rg] <- paste0("(", tempMods$Name[m], ")")
  return(c(paste(x1, collapse = ""),
           paste(x2, collapse = "")))
}
