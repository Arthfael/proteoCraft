#' gsub_Rep
#'
#' @description 
#' A function which in some specific situations will yield a performance increase over classic gsub.
#' This is only useful for large vectors of few repeated values.
#' All arguments are passed to base::gsub
#' The difference with normal gsub is that first a vector of unique values is created, on which gsub is called.
#' The input values are then matched to the unique ones and replaced with the corresponding edited result.
#' 
#' @param pattern	A character string containing a regular expression (or character string for fixed = TRUE) to be matched in the given character vector. Coerced by as.character to a character string if possible. If a character vector of length 2 or more is supplied, the first element is used with a warning. Missing values are allowed.
#' @param replacement A replacement for matched pattern in sub and gsub. Coerced to character if possible. For fixed = FALSE this can include backreferences "\1" to "\9" to parenthesized subexpressions of pattern. For perl = TRUE only, it can also contain "\U" or "\L" to convert the rest of the replacement to upper or lower case and "\E" to end case conversion. If a character vector of length 2 or more is supplied, the first element is used with a warning. If NA, all elements in the result corresponding to matches will be set to NA.
#' @param x a character vector where matches are sought, or an object which can be coerced by as.character to a character vector. Long vectors are supported.
#' @param ignore.case	If FALSE, the pattern matching is case sensitive and if TRUE, case is ignored during matching.
#' @param perl Logical. Should Perl-compatible regexps be used?
#' @param fixed	Logical. If TRUE, pattern is a string to be matched as is. Overrides all conflicting arguments.
#' @param useBytes Logical. If TRUE the matching is done byte-by-byte rather than character-by-character. See ‘Details’.
#' 
#' @examples
#' tA1 <- sys.time()
#' PSMs$Sample <- gsub(".*/|\\.d$", "", PSMs$"Raw file")
#' tA2 <- sys.time()
#' tB1 <- sys.time()
#' PSMs$Sample <- gsub_Long(".*/|\\.d$", "", PSMs$"Raw file")
#' tB2 <- sys.time()
#' tA2 - tA1
#' tB2 - tB1
#' 
#' @export

gsub_Rep <- function(pattern,
                      replacement,
                      x,
                      ignore.case = FALSE,
                      perl = FALSE,
                      fixed = FALSE,
                      useBytes = FALSE) {
  u <- unique(x)
  nu <- gsub(pattern, replacement, u, ignore.case, perl, fixed, useBytes)
  nu[match(x, u)]
}
