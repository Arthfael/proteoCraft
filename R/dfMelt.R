#' dfMelt
#'
#' @description 
#' Wrapper function for utils::stack(), allowing much faster melting of data.frames or lists of data.frames!
#' 
#' @param df A data.frame (matrices get silently converted to data.frame). May also be a list of data.frames, which is then traversed recursively until a data.frame is hit.
#' @param ColNames Column names of the output data frame. Default = c("variable", "value")
#' @param id.vars Id variables - column name only, column indices are not supported! If missing, all columns are assumed to be measure variables.
#' 
#' @details
#' Meant for "matrix-in-a-data.frame", does not handle columns of type list!
#' Should be much faster than reshape::melt.data.frame(), and frees us from this dependency.
#' But this is an in-house solution so may not fully substitute.
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
  if (inherits(df, "list")) {
    # If we have a list
    warning("dfMelt: list input detected!")
    if (is.null(names(df))) { names(df) <- as.character(1L:length(df)) }
    res <- lapply(names(df), \(nm) { #nm <- names(df)[1L]
      x <- dfMelt(df[[nm]], id.vars = id.vars)
      g <- grep("^L[0-9]+$", colnames(x), value = TRUE)
      kol <- if (length(g)) { paste0("L", max(as.integer(sub("^L", "", g)))+1L) } else { "L1" }
      x[[kol]] <- nm
      return(x)
    })
    return(do.call(plyr::rbind.fill, res))
  }
  if (is.matrix(df)) { df <- as.data.frame(df) } # Stack needs a data.frame
  stopifnot(is.data.frame(df),
            nrow(df) > 0L)
  id.vars_tst <- (!missing(id.vars))&&(!is.null(id.vars))
  if (id.vars_tst) {
    l_IDvars <- length(id.vars)
    stopifnot(sum(!id.vars %in% colnames(df)) == 0L)
    ids <- df[, id.vars, drop = FALSE]
    df <- df[, base::setdiff(colnames(df), id.vars), drop = FALSE]
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
