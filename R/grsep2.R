#' grsep2
#'
#' @description 
#' A function to look for IDs in a separated character list.
#' Unlike "grsep", it is not a wrapper for "grep/grepl". Unlike it, it should not break for large number of IDs.
#' Also unlike grsep it does not have "gsub-like" capabilities (yet?)
#' This function is based on reshape2; it is significantly faster than grsep on large vectors, but can be slower on smaller ones.
#' 
#' @param IDs The IDs to search for.
#' @param x A character vector listing IDs in which to search.
#' @param sep The characters list's separator. Default = ";"
#' @param ignore.case The base parameter (default = FALSE)
#' @param value Return values from matches (TRUE) or just match ids (FALSE, default)
#' @param invert The base parameter (default = FALSE)
#' 
#' @examples
#' tst <- grsep(c("MYFAVOURITEPROTEINSACCESSION", "ALSOTHISONE", "ANDHERESANOTHERONEFORGOODMEASURE"), PG$"Protein IDs")
#' 
#' @export

grsep2 <- function(IDs, x, sep = ";", ignore.case = FALSE, value = FALSE, invert = FALSE) {
  #proteoCraft::DefArg(proteoCraft::grsep2)
  stopifnot(length(x) > 0, is.null(dim(as.character(x))))
  Wh <- 1:length(x)
  x <- as.character(x)
  g1 <- grep(sep, x)
  g2 <- grep(sep, x, invert = TRUE)
  if (length(g1)) {
    w <- which(sapply(x[g1], length) > 0)
    temp <- proteoCraft::listMelt(strsplit(x[g1[w]], sep), g1[w])
  }
  if (length(g2)) {
    temp2 <- data.frame(value = x[g2], L1 = g2)
    if (length(g1)) { temp <- rbind(temp, temp2) } else { temp <- temp2 }
  }
  if (ignore.case) {
    IDs <- tolower(IDs)
    temp$value <- tolower(temp$value)
  }
  res <- sort(as.numeric(unique(temp$L1[which(temp$value %in% IDs)])))
  if (invert) { res <- Wh[which(!Wh %in% res)] }
  if (value) { res <- x[res] }
  return(res)
}
