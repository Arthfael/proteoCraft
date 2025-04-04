#' GO_annotate_PG
#'
#' @description 
#' A function to process the results of the GO_annotate_DB function into new columns added to a protein groups file.
#'
#' @param Prot The protein groups file.
#' @param annotations The annotations file; must be the output of the GO_annotate_DB function.
#' @param ids The name of the column with the protein IDs; usually "Protein IDs" (default).
#' 
#' @examples
#' annot <- GO_annotate_DB(db$Protein.ID, "Homo sapiens")
#' GO_annotate_PG(PG, annotations = annot)
#'
#' @export

GO_annotate_PG <- function(Prot,
                           annotations,
                           ids = "Protein IDs") {
  GO <- setNames(c("BP", "CC", "MF"), nm = c("biological_process", "cellular_component", "molecular_function"))
  GO.col <- as.character(sapply(GO, function(x) {paste("go_", x, c("", "_name"), sep = "")}))
  a <- strsplit(Prot[[ids]], split = ";")
  temp <- annotations$Annotated_Proteins
  temp2 <- as.data.frame(t(sapply(a, function(x) {
    x <- which(temp$"Protein ID" %in% unlist(x))
    if (length(x) == 0) { x <- rep("", length(GO.col)) } else {
      if (length(x) > 1) {
        x <- temp[x, GO.col]
        x <- apply(x, 2, function(y) {
          paste(unique(unlist(strsplit(unlist(y), split = ";"))), collapse = ";")
        })
      } else { x <- temp[x, GO.col] }
    }
    return(x)
  })))
  for (i in colnames(temp2)) { temp2[[i]] <- unlist(temp2[[i]]) }
  Prot[, colnames(temp2)] <- temp2
  return(Prot)
}
