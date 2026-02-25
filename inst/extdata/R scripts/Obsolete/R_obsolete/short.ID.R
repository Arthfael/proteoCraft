#' short.ID
#'
#' @description 
#' A convenient function to create a meaningful tag for proteins in graphs, using the gene and protein columns.
#' This is not being used at all by any of the functions -> retired on 03/09/2025!
#' 
#' @param Prot The protein/protein groups file.
#' @param Proteins.colname The name of the protein IDs column.
#' @param Genes.colname The name of the gene IDs column.
#' 
#' @examples
#' PG$Short.ID <- function(PG, Proteins.colname = "Protein.IDs", Genes.colname = "Genes")
#' 
#' @export

short.ID <- function(Prot, Proteins.colname, Genes.colname) {
  a <- strsplit(Prot[[Proteins.colname]], ";")
  b <- strsplit(Prot[[Genes.colname]], ";")
  Res <- apply(cbind(a, b), 1, function(x) {
    x1 <- unlist(x[1])
    x2 <- unlist(x[2])
    if (length(x2)) {
      if (length(x2) > 1) { res <- paste0("gn_", x2[1], "...") } else { res <- paste0("gn_", x2) }
    } else {
      if (!length(x1)) { res <- "pr_unnamed protein" } else {
        if (length(x1) > 1) { res <- paste0("pr_", x1[1], "...") } else { res <- paste0("pr_", x1) }
      }
    }
    return(res)
  })
  return(Res)
}
