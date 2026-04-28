#' .SK2MQ_modSeqWrkr1
#' 
#' Worker function used by Skyline_to_MQ().
#'
#' @param x Temporary modified sequence
#' @param tempMods Temporary modifications table
#' 
#' @export

.SK2MQ_modSeqWrkr1 <- function(x,
                               tempMods) {
  #x <- temp1[[1]]
  x1 <- x2 <- unlist(x)
  l <- length(x1)
  rg <- c(1:((l-1)/2))*2
  nc <- nchar(x1[rg-1])
  aa <- substr(x1[rg-1], nc, nc)
  if (nc[1] == 1) { aa[1] <- substr(x1[3], 1, 1) }
  m <- match(x1[rg], tempMods$Old)
  x1[rg] <- paste0("(", tempMods$Mark[m], ")")
  x2[rg] <- paste0("(", tempMods$Name[m], ")")
  return(c(paste(x1, collapse = ""),
           paste(x2, collapse = "")))
}
