#' listMelt
#'
#' @description 
#' Wrapper function for utils::stack(), allowing much faster and efficient melting of simple lists, i.e. lists where each element is a 1-dimension vector.
#' This won't work on more complex lists!!!
#' 
#' @param List A list.
#' @param Names Names to apply to each element, if missing. Can be convenient in case you want to preserve Name's original non-character data type (e.g. numeric).
#' @param ColNames Column names of the output data frame. Default = c("value", "L1"); note the reverse order compared with classic melt!
#' 
#' @details
#' This function allows much faster melting than reshape::melt.list(). However, it may only be used on "simple" lists, that is, lists where each element is a vector.
#' 
#' @returns
#' A data.frame with two columns. The default column names are "value" and "L1" as in classic melt(). The column order is inverted, i.e. "value" then "L1" - a stupid oversight, unfortunately changing this would be difficult now.
#' 
#' @examples
#' temp <- listMelt(strsplit(PG$"Protein IDs", ";"), PG$id, c("Protein accession", "PG id"))
#' 
#' @export

listMelt <- function(List,
                     Names,
                     ColNames = c("value", "L1")) {
  #List <- prjcts
  toInt <- FALSE
  stopifnot(is.atomic(List)|is.list(List),
            length(ColNames) == 2L)
  lL <- length(List)
  w <- which(lengths(List) > 0L) # Filter necessary not only to avoid warning for empty elements ("In stack.default(List) : non-vector elements will be ignored")
  # But also a bug which would cause improper names assignment if the list contains empty elements!!!
  if (!missing(Names)) {
    stopifnot(length(Names) == lL)
    # Important to catch improper use!!!
  } else {
    if (!is.null(names(List))) {
      Names <- names(List)
    } else {
      Names <- as.character(1L:lL)
      toInt <- TRUE
    }
  }
  # if (length(Names) != lL) {
  #   Names <- as.character(1L:lL)
  #   toInt <- TRUE
  # }
  List <- setNames(List[w], as.character(w))
  List <- setNames(lapply(List, unlist), names(List)) # Very unsafe!
  List <- utils::stack(List)
  Names <- Names[w] # Do NOT forget that filter!
  colnames(List) <- c(ColNames[1L], "ind")
  List$ind <- as.integer(List$ind)
  List[[ColNames[2]]] <- Names[List$ind]
  if (toInt) { List[[ColNames[2]]] <- as.integer(List[[ColNames[2]]]) }
  List$ind <- NULL
  return(List)
}
