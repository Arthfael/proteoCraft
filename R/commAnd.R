#' .commAnd
#'
#' @description 
#' \preformatted{
#'               comma           and
#'                comma         and
#'                 comma       and
#'                  comma     and
#'                   comma   and
#'                    comma and
#'                     commAnd
#' }
#' ...
#' and now I hope you see why this function is named like this...?
#'
#' Anyhow...
#' 
#' This weirdly named function is meant to simplify a step I have had to write waaaaay too many times:\cr
#' It collapse a vector (e.g. c("A", "B", "C")) onto a single character string of the shape "A, B and C" for the purpose of writing human-readable text.
#' 
#' @returns
#' A single character vector.
#' 
#' @param x A vector of any length. Unlisted and converted to character if necessary.
#' @param sort Logical, default = FALSE. If TRUE, values are sorted before collapsing.
#' @param unique Logical, default = FALSE. If TRUE, only unique values are used before collapsing.
#' @param sep1 1st separator, default = ", "
#' @param sep2 2nd separator, default = " and "
#' @param end Final character, default = "". Use "." for final dot.
#' @param maxN Round number (numeric or integer). How many items should be collapsed? ; if more than the max number of objects are present,  
#' 
#' @examples
#' .commAnd(c("Charles", "Alice", "Bob", "Bob"))
#' .commAnd(c("Charles", "Alice", "Bob", "Bob"), TRUE)
#' .commAnd(c("Charles", "Alice", "Bob", "Bob"), TRUE, TRUE)
#' .commAnd(list(A = c("David", "Emilia"), B = c("Charles", "Alice", "Bob", "Bob")), FALSE, TRUE)
#' .commAnd(c("Charles", "Alice", "Bob"), TRUE, sep1 = " ", sep2 = " or ", end = ".")
#' .commAnd(list(A = c("David", "Emilia", "Fiona"),
#'               B = c("Charles", "Alice", "Bob", "Bob")), maxN = 5L)
#' .commAnd(list(A = c("David", "Emilia", "Fiona"),
#'               B = c("Charles", "Alice", "Bob", "Bob")), TRUE, maxN = "5")
#' .commAnd(list(A = c("David", "Emilia", "Fiona"),
#'               B = c("Charles", "Alice", "Bob", "Bob")), TRUE, TRUE, maxN = "5")
#' 
#' @export

.commAnd <- function(x,
                     sort = FALSE,
                     unique = FALSE,
                     sep1 = ", ",
                     sep2 = " and ",
                     end = "",
                     maxN) {
  #DefArg(.commAnd); x <- c("Charles", "Alice", "Bob", "Bob")
  x <- as.character(unlist(x)) # As.character is important to avoid issues with numerics, otherwise additional " " may be added to shorter numbers!
  if (unique) { x <- unique(x) }
  l <- length(x)
  if (l == 0L) { x <- "" }
  if (l <= 1L) { return(paste0(x, end)) }
  if (sort) { x <- sort(x) }
  if (!missing(maxN)) {
    maxN <- suppressWarnings(as.integer(maxN))
    if ((!is.finite(maxN)) || (maxN < 1L)) { warning("Invalid maxN parameter, ignoring...") } else {
      if (maxN < l) {
        return(paste0(paste(x[1L:maxN], collapse = sep1), ", ..."))
      }
    }
  }
  paste0(paste(x[1L:(l-1L)], collapse = sep1), sep2, x[l], end)
}
