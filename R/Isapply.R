#' Isapply
#' 
#' @description 
#' I love the sapply function, but its output can be difficult to deal with:
#' It can require different treatment depending on whether results per x are expected to be length 1 or more.
#' This is a wrapper so I don't have to rewrite the same stuff all the time, hopefully it will help.
#' It does not always work currently... 
#' 
#' @inheritParams base::sapply
#' @param col.names Optional column names to provide to the output, if it will be 2D. Ignored (with a warning) if output is 1D.
#' @param convert.to By default ("df" or 1), converts 2D output to a data.frame. If set as "m" or "2", converts it to a matrix.
#' 
#' @export

Isapply <- function(x, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, col.names, convert.to = "df") {
  res <- sapply(x, FUN = FUN, ..., simplify = simplify, USE.NAMES = USE.NAMES)
  if (sum(c("data.frame", "matrix") %in% class(res)) > 0) {
    convert.to <- tolower(as.character(convert.to))
    if (!convert.to %in%  c("df", "dataframe", "data.frame", "1", "matrix", "m", "2")) {
      warning("Argument \"convert.to\" unrecognized, defaulting to converting 2D output to a data.frame!")
      convert.to <- "df"
    }
    nr <- nrow(res)
    nc <- ncol(res)
    if (nc == length(x)) {
      if (convert.to %in% c("df", "dataframe", "data.frame", "1")) {
        res <- as.data.frame(t(res))
      }
      if (convert.to %in% c("matrix", "m", "2")) {
        res <- as.matrix(t(res))
      }
    } else {
      warning("I did not expect this to happen! You may want to check the \"Isapply\" function's code.")
    }
  }
  if (!missing("col.names")) {
    if  (sum(c("data.frame", "matrix") %in% class(res)) > 0) {
      colnames(res) <- col.names
    } else {
      warning("Argument \"col.names\" will be ignored as the result is a vector, not tabular data!")
    }
  }
  return(res)
}
