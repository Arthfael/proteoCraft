#' grsep
#'
#' @description 
#' A wrapper function for grep/grepl/gsub to search IDs in a separated characters vector listing a set of IDs.
#' N.B.: Other modes could easily be added to the code without major rewrites of this or other scripts using it.
#' 
#' @param pattern The IDs to search for.
#' @param replacement The replacement, only used if mode = gsub.
#' @param x A character vector listing IDs in which to search.
#' @param sep The characters list's separator. Default = ";"
#' @param mode Should this work in grep, gsub or grepl mode?
#' @param ignore.case The base parameter (default = FALSE)
#' @param perl The base parameter (default = FALSE)
#' @param value The base parameter (default = FALSE)
#' @param fixed The base parameter (default = FALSE)
#' @param useBytes The base parameter (default = FALSE)
#' @param invert The base parameter (default = FALSE)
#' 
#' @examples
#' tst <- grsep(c("MYFAVOURITEPROTEINSACCESSION", "ALSOTHISONE", "ANDHERESANOTHERONEFORGOODMEASURE"), x = PG$"Protein IDs")
#' 
#' @export

grsep <- function (pattern,
                   replacement = "",
                   x,
                   sep = ";", 
                   mode = "grep",
                   ignore.case = FALSE,
                   perl = FALSE,
                   value = FALSE, 
                   fixed = FALSE,
                   useBytes = FALSE,
                   invert = FALSE) {
  stopifnot(mode %in% c("grep", "grepl", "gsub"))
  x <- paste0(sep, x, sep)
  pattern <- paste(paste0(sep, pattern, sep), collapse = "|")
  if (mode == "grep") {
    res <- grep(pattern, x, ignore.case = ignore.case, perl = perl, value = value, fixed = fixed, useBytes = useBytes,
                invert = invert)
    if ((nchar(sep) > 0) && (value == TRUE)) {
      res <- substr(res, start = nchar(sep) + 1, stop = nchar(res) - nchar(sep))
    }
  }
  if (mode == "grepl") {
    res <- grepl(pattern, x, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
  }
  if (mode == "gsub") {
    res <- gsub(pattern, replacement, x, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
  }
  return(res)
}
