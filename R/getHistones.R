#' getHistones
#'
#' @description
#' This function takes a data.frame resulting from the parsing of a fasta file by Format.DB() and text mines its header to identify histones.
#' 
#' @param DB A data.frame with at least "Header" and "Protein ID" columns.
#' 
#' @returns
#' A list of 2 vectors:
#'  - All = accessions of all histones
#'  - Core = accessions of core histones
#' 
#' @examples
#' histIDs <- getHistones(db)
#' 
#' @export

getHistones <- function(DB) {
  Hist <- grep("histone", DB$Header, ignore.case = TRUE, value = TRUE)
  Hist <- grep("ase(,|\\.| )", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
  Hist <- grep("-(binding|like)", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
  Hist <- grep("chaperone", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
  Hist <- grep("ethyl", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
  Hist <- grep("factor", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
  Hist <- grep("non-histone", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
  Hist <- grep("\\(fragment\\)", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
  coreHist <- grep("H1|H5|linker", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
  HistIDs <- DB$`Protein ID`[match(Hist, DB$Header)]
  coreHistIDs <- DB$`Protein ID`[match(coreHist, DB$Header)]
  #histDB <- DB[match(HistIDs, DB$`Protein ID`),]
  return(list(All = HistIDs,
              Core = coreHistIDs))
}
