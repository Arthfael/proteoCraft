#' dfMelt
#'
#' @description 
#' Wrapper function for utils::stack(), allowing much faster melting of simple - and only simple! - data.frames!
#' This should only be used for what is essentially a "matrix in a data.frame"!
#' It will assume that everything in the table is a value and that column names are variable names.
#' In this case, it should be much faster than reshape::melt.data.frame()
#' 
#' @param df A data.frame (matrices get silently converted to data.frame).
#' @param ColNames Column names of the output data frame. Default = c("variable", "value")
#' 
#' @examples
#' df <- as.data.frame(matrix(rnorm(1000), ncol = 10))
#' df2 <- dfMelt(df)
#' 
#' @export

dfMelt <- function(df,
                   ColNames = c("variable", "value")) {
  #df <- as.data.frame(matrix(rnorm(1000), ncol = 10))
  if ("matrix" %in% class(df)) { df <- as.data.frame(df) } # Stack needs a data.frame
  stopifnot("data.frame" %in% class(df),
            length(ColNames) == 2)
  df <- utils::stack(df)
  df <- df[, 2:1]
  colnames(df) <- ColNames
  return(df)
}
