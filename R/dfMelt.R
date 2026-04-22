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
#' df <- as.data.frame(matrix(rnorm(1000), ncol = 10L))
#' df2 <- dfMelt(df)
#' 
#' @export

dfMelt <- function(df,
                   ColNames = c("variable", "value"),
                   id.vars) {
  #df <- as.data.frame(matrix(rnorm(1000), ncol = 10L))
  #df <- temp; ColNames <- c("Label", "Sample", "value"); id.vars <- "Rowname"
  if (is.matrix(df)) { df <- as.data.frame(df) } # Stack needs a data.frame
  stopifnot(is.data.frame(df),
            nrow(df) > 0L)
  id.vars_tst <- (!missing(id.vars))&&(!is.null(id.vars))
  if (id.vars_tst) {
    l_IDvars <- length(id.vars)
    stopifnot(sum(!id.vars %in% colnames(df)) == 0L)
    ids <- df[, id.vars, drop = FALSE]
    df <- df[, setdiff(colnames(df), id.vars), drop = FALSE]
    if (missing(ColNames)) {
      ColNames <- c(id.vars, ColNames)
    } else { stopifnot(length(ColNames) == 2L + l_IDvars) }
  } else {
    stopifnot(length(ColNames) == 2L)
    id.vars <- c()
    l_IDvars <- 0L
  }
  df <- utils::stack(df)
  df <- df[, 2L:1L]
  colnames(df) <- ColNames[(1L + l_IDvars):length(ColNames)]
  if (l_IDvars) {
    colnames(ids) <- ColNames[1L:l_IDvars]
    rownames(ids) <- NULL
    df <- cbind(ids, df) # recycling ids which now has fewer rows
  }
  return(df)
}
