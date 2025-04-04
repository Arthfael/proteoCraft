#' listMelt
#'
#' @description 
#' Wrapper function for utils::stack(), allowing much faster and efficient melting of simple lists, i.e. lists where each element is a 1-dimension vector.
#' This won't work on more complex lists!!!
#' 
#' Allows the following:
#'  - Much faster melting than reshape::melt.list() for large lists
#'  - Avoids using equally slow reshape::melt(), which is apparently deprecated and causing issues with data.tables
#'  - Simplifies the syntax of the rewrites
#' 
#' @param List A list.
#' @param Names Names to apply to each element, if missing. Can be convenient in case you want to preserve Name's original non-character data type (e.g. numeric).
#' @param ColNames Column names of the output data frame. Default = c("value", "L1"); note the reverse order compared with classic melt!
#' 
#' @examples
#' temp <- listMelt(strsplit(PG$"Protein IDs", ";"), PG$id, c("Protein accession", "PG id"))
#' 
#' @export

listMelt <- function(List,
                     Names = NULL,
                     ColNames = c("value", "L1")) {
  #List <- prjcts
  stopifnot(is.atomic(List)|is.list(List),
            length(ColNames) == 2)
  lL <- length(List)
  w <- which(vapply(List, length, integer(1)) > 0) # Filter necessary not only to avoid warning for empty elements ("In stack.default(List) : non-vector elements will be ignored")
  # But also a bug which would cause improper names assignment if the list contains empty elements!!!
  if (!is.null(Names)) {
    stopifnot(length(Names) == lL)
    # Important to catch improper use!!!
  } else {
    if (length(names(List)) == lL) {
      Names <- names(List)
    } else {
      Names <- 1:lL
    }
  }
  List <- setNames(List[w], w)
  List <- utils::stack(List)
  #if (!is.null(Names)) { # Now with this rewrite names are never NULL at this stage
    Names <- Names[w] # Do NOT forget that filter!
    colnames(List) <- c(ColNames[1], "ind")
    List$ind <- as.integer(List$ind)
    List[[ColNames[2]]] <- Names[List$ind]
    List$ind <- NULL
  # } else {
  #   colnames(List) <- ColNames
  #   List[[ColNames[2]]] <- as.integer(List[[ColNames[2]]])
  # }
  return(List)
}
