#' insertElems
#' 
#' @description 
#' Nifty little function to insert elements in a vector after the indicated positions.
#' 
#' @param vect The target vector of length N.
#' @param pos The k positions where k elements should be inserted; pos can range from 0 (before the first item) to N (after the last).
#' @param elems The k elements to insert.
#' @param warn Should we display a warning when we have no choice but to convert a vector from numeric to character? Default = FALSE.
#' 
#' @examples
#' insertElems(c(1:4), 3, "A")
#' # [1] "1" "2" "A" "3" "4"
#' 
#' @export

insertElems <- function(vect, pos, elems, warn = FALSE) {
  tt <- setNames(c("logical", "integer", "integer", "numeric", "character", "list"),
                 c("logical", "integer", "integer64", "numeric", "character", "list"))
  klasses <- unique(c(class(vect), class(elems)))
  if (sum(which(!klasses %in% names(tt)) > 0L)) {
    stop("This function can only deal with vectors of class \"logical\", \"integer\", \"numeric\", \"character\" and \"list\"!")
  }
  konvert <- FALSE
  if (length(klasses) > 1L) {
    konvert <- TRUE
    ok <- FALSE
    for (i in rev(names(tt))) {
      if ((!ok)&&(i %in% klasses)) {
        if (warn) { warning(paste0("The output will be converted to the type \"", i, "\".")) }
        f <- get(paste0("as.", tt[i]))
        ok <- TRUE
      }
    }
  }
  lPos <- length(pos)
  if (lPos != length(unique(pos))) { stop("Duplicate positions provided! Only one insertion is allowed for each position!") }
  if (lPos != length(elems)) { stop("Different number of elements to insert and positions provided!") }
  l1 <- length(vect)
  if (sum(!pos %in% 0L:l1)) { stop("The supplied positions do not make sense") }
  elems <- elems[order(pos)]
  pos <- sort(pos)
  pos0 <- FALSE
  if (0L %in% pos) {
    pos0 <- TRUE
    elem0 <- elems[1L]
    if (lPos > 1L) {
      pos <- pos[2L:lPos]
      elems <- elems[2L:lPos]
      lPos <- length(pos)
    }
  }
  res <- rep(NA, l1+lPos)
  if (lPos) {
    pos <- pos+c(1L:lPos)
    res[pos] <- elems
    res[which(is.na(res))] <- vect
  }
  if (pos0) { res <- c(elem0, res) }
  if (konvert) { res <- f(res) }
  if (("list" %in% klasses)&&(!inherits(res, "list"))) { res <- as.list(res) }
  return(res)
}
