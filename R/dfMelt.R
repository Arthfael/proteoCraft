#' dfMelt
#'
#' @description 
#' Wrapper function for utils::stack(), allowing much faster melting of simple - and only simple! - data.frames!
#' 
#' @param df A data.frame (matrices get silently converted to data.frame).
#' @param ColNames Column names of the output data frame. Default = c("variable", "value")
#' @param id.vars Id variables - column name only (indices are not supported!) If missing, all columns are assumed to be measure variables.
#' 
#' @details
#' This should only be used for what is essentially a "matrix in a data.frame"!
#' It will assume that everything in the table is a value and that column names are variable names.
#' In this case, it should be much faster than reshape::melt.data.frame()
#' It cannot be thought of as a full substitute for the latter though. Never use it for complex data.frames, e.g. where columns are lists. 
#' 
#' @returns
#' A melted (molten?) data.frame.
#' 
#' @examples
#' df <- as.data.frame(matrix(rnorm(1000), ncol = 10))
#' df2 <- dfMelt(df)
#' 
#' @export

dfMelt <- function(df,
                   ColNames = c("variable", "value"),
                   id.vars) {
  #df <- as.data.frame(matrix(rnorm(1000), ncol = 10))
  if ("matrix" %in% class(df)) { df <- as.data.frame(df) } # Stack needs a data.frame
  stopifnot("data.frame" %in% class(df),
            length(ColNames) == 2)
  if (!missing(id.vars)) {
    stopifnot(sum(!id.vars %in% colnames(df)) == 0)
    ids <- df[, id.vars]
    df <- df[, which(!colnames(df) %in% id.vars)]
  }
  df <- utils::stack(df)
  df <- df[, 2:1]
  colnames(df) <- ColNames
  if (!missing(id.vars)) {
    df[, id.vars] <- ids
    df <- df[, c(id.vars, ColNames)]
  }
  return(df)
}
