#' topattern
#'
#' @description 
#' A function to convert character strings to regex patterns
#' 
#' @param x The character vector to be turned into a regex pattern.
#' @param start Default = TRUE, which means the pattern will require a match to the start of the string.
#' @param end Default = FALSE. Set to TRUE for the pattern require a match to the end of the string.
#' @param collapse For vectors of length > 1, character to collapse the vector to in order to create a single pattern. Default = "|". Set to FALSE to turn it off.
#'
#' @details
#' Internal function to convert a chain of characters into a regex pattern usable by gsub() or grep().
#' Note: this may not work for absolutely every string, but is a very useful shortcut for the controlled usage we make of it within this package.
#' 
#' @returns
#' A regex pattern.
#' 
#' @examples
#' topattern(pep.ref[length(pep.ref)])
#' 
#' @export

topattern <- function(x, start = TRUE, end = FALSE, collapse = "|") {
  #
  # Check arguments
  stopifnot(length(start) == 1,
            length(end) == 1,
            length(collapse) == 1,
            is.logical(start),
            is.logical(end))
  if (is.na(start)) { start <- TRUE }
  if (is.na(end)) { start <- FALSE }
  if (is.logical(collapse)) {
    if (is.na(collapse)||collapse) { collapse <- "|" }
  }
  #
  x <- gsub("\\\\", "\\\\\\\\", as.character(x)) # Should be first so as not to escape escape signs!
  x <- gsub("\\.", "\\\\.", x)
  x <- gsub("\\*", "\\\\*", x)
  x <- gsub("\\$", "\\\\$", x)
  x <- gsub("\\^", "\\\\^", x)
  x <- gsub("\\+", "\\\\+", x)
  x <- gsub("\\?", "\\\\?", x)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  x <- gsub("\\[", "\\\\[", x)
  x <- gsub("\\]", "\\\\]", x)
  x <- gsub("\\(", "\\\\(", x)
  x <- gsub("\\)", "\\\\)", x)
  x <- gsub("\\|", "\\\\|", x)
  if (start) { x <- paste0("^", x) }
  if (end) { x <- paste0(x, "$") }
  if ((length(x) > 1)&&(is.character(collapse))) { x <- paste(x, collapse = collapse) }
  return(x)
}
