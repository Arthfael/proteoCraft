#' topattern
#'
#' @description 
#' A function to convert character strings to regex patterns, not necessarily universal but based on my common usage for my script.
#' 
#' @param x The character vector to be turned into a regex pattern.
#' @param start Default = TRUE, which means the pattern will require a match to the start of the string.
#' @param end Default = FALSE. Set to TRUE for the pattern require a match to the end of the string.
#' @param collapse For vectors of length > 1, character to collapse the vector to in order to create a single pattern. Default = "|". Set to FALSE to turn it off.
#'
#' @examples
#' topattern(pep.ref[length(pep.ref)])
#' 
#' @export

topattern <- function(x, start = TRUE, end = FALSE, collapse = "|") {
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
  if ((length(x) > 1)&&(collapse != FALSE)) { x <- paste(x, collapse = collapse) }
  return(x)
}
