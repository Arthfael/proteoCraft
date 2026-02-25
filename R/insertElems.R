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

insertElems = function(vect, pos, elems, warn = FALSE) {
  tt <- setNames(c("logical", "integer", "integer", "numeric", "character", "list"),
                 c("logical", "integer", "integer64", "numeric", "character", "list"))
  klasses <- unique(c(class(vect), class(elems)))
  if (sum(which(!klasses %in% names(tt)) > 0)) {
    stop("This function can only deal with vectors of class \"logical\", \"integer\", \"numeric\", \"character\" and \"list\"!")
  }
  konvert <- FALSE
  if (length(klasses) > 1) {
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
  if (length(pos) != length(unique(pos))) { stop("Duplicate positions provided! Only one insertion is allowed for each position!") }
  if (length(pos) != length(elems)) { stop("Different number of elements to insert and positions provided!") }
  l1 <- length(vect)
  if (length(which(!pos %in% c(0:(l1))))) { stop("The supplied positions do not make sense") }
  elems <- elems[order(pos)]
  pos <- sort(pos)
  l2 <- length(pos)
  pos0 <- FALSE
  if (0 %in% pos) {
    pos0 <- TRUE
    elem0 <- elems[1]
    if (l2 > 1) {
      pos <- pos[2:l2]
      elems <- elems[2:l2]
      l2 <- length(pos)
    }
  }
  res <- rep(NA, l1+l2)
  if (l2 > 0) {
    pos <- pos+c(1:l2)
    res[pos] <- elems
    res[which(is.na(res))] <- vect
  }
  if (pos0) { res <- c(elem0, res) }
  if (konvert) { res <- f(res) }
  if (("list" %in% klasses)&&(class(res) != "list")) { res <- as.list(res) }
  return(res)
}
