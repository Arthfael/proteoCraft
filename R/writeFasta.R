#' writeFasta
#'
#' @description 
#' A function to write a data table of protein sequences, parsed from a fasta by Format.DB, as a fasta.
#' The fasta is:
#'  - written down if a destination file is provided using destFl,
#'  - or, if newFl is TRUE, a modal dialog is opened to interactively choose it.
#'  - The fasta can also optionally be returned as a character string.
#' 
#' @param DB A data base of protein sequences. Must contain "Header" and "Sequence" column.
#' @param destFl Destination file path.
#' @param filter A character of protein accessions (IDs) to use as (positive) filter for the output. Requires that DB contain a "Protein ID" column.
#' @param return Should the resulting character string be returned (TRUE), or only written to destFl (default)?
#' @param newFl New file? If so the user will be prompted for the path of the destination file (and destFl will be ignored). Default = FALSE
#' 
#' @examples
#' DB <- Format.DB(".../My fasta database.fasta")
#' path <- ".../My fasta database's backup.fasta"
#' writeFasta(DB, path)
#' fasta <- writeFasta(DB, return = TRUE) # Assign as object in environment without writing to file.
#' fasta <- writeFasta(DB, path, return = TRUE) # Write to file & assign as object in environment.
#' # Filter
#' filt <- c("Protein 1", "Protein 2")
#' path <- ".../My filtered fasta database.fasta"
#' writeFasta(DB, path, filt)
#' 
#' @export

writeFasta <- function(DB, destFl, filter, return = FALSE, newFl = FALSE) {
  nr <- nrow(DB)
  stopifnot(nrow(DB) > 0,
            length(unique(DB$"Protein ID")) == nrow(DB))
  if ((!missing("filter"))&&("Protein. ID" %in% colnames(DB))) {
    DB <- DB[which(DB$"Protein ID" %in% filter),]
    nr <- nrow(DB)
    stopifnot(nrow(DB) > 0)
  }
  tmp <- rep("", nr*3)
  tmp[(1:nr)*3-2] <- DB$Header
  tmp[(1:nr)*3-1] <- DB$Sequence
  if (newFl) { # Ignores existing value of destFl
    filt <- matrix(data = c("fasta", "*.fasta;*.fas;*.fa;*.fasta.fas"), ncol = 2,
                   dimnames = list("Fasta"))
    destFl <- svDialogs::dlg_save("proteins of interest.fasta", "Select destination file", filters = filt)$res
  }
  if (!missing("destFl")) { writeLines(tmp, destFl) }
  if (return) { return(tmp) }
}
