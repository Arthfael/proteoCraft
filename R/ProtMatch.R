#' ProtMatch
#'
#' @description
#' 
#' This function is very slow and should not be used! It is only provided for legacy purposes. Instead use ProtMatch2, which is ~50-to-1000 times faster!
#' 
#' Some search software do not report all proteins from the database which match a specific sequence.
#' Others, such as MaxQuant, do report weird matches for a minority of PSMs, or miss a few perfectly valid matches.
#' While this may be correct from their point of view, this function allows re-checking matches for consistency between software.
#' This will return a data frame of matches between observed peptides and protein IDs from a formatted protein sequences database.
#' 
#' @param Seq Character vector of peptide sequences.
#' @param DB The formated protein database, with proteins IDs and sequences.
#' @param Cut For digestion. The character vector of amino acid(s) after which to cut. For each, provide one letter followed or preceded by underscore indicating whether the cut is C- or N-terminal, respectively. Default (for Trypsin) = c("K_", "R_")
#' @param strict.avoid For digestion. The pattern(s) that will strictly block cleavage when overlapping with cleavage sites. Write cleavage site as underscore. Default = ""
#' @param loose.avoid For digestion. The pattern(s) that will partially block cleavage when overlapping with cleavage sites. Write cleavage site as underscore. Default =  c("K_P", "R_P")
#' @param min For digestion. The minimum length of peptides to report. Default = 7
#' @param max For digestion. The maximum length of peptides to report. Set to 100 (default) to allow even very large peptides. If one of the peptides in Seq in longer, than that length will be used instead.
#' @param missed Let's say we have observations with up to 6 missed cleavages. This will probably break the function. In this case, it is better to specify a smaller number of missed cleavages, e.g. 2.

ProtMatch <- function(Seq,
                      DB,
					  Cut = c("K_", "R_"),
					  strict.avoid = "",
					  loose.avoid = c("K_P", "R_P"),
                      min = 7,
					  max = 100,
					  missed = 2) {
  #proteoCraft::DefArg(proteoCraft::ProtMatch); w <- 1:nrow(ev); Seq = unique(ev$Sequence[w]); DB = db; min = MinPepSz
  stopifnot( min(nchar(Cut)) == 2, max(nchar(Cut)) == 2)
  if (!is.logical(max)) { max <- max(c(max, nchar(Seq))) }
  Seq <- data.frame(Sequence = Seq)
  #n0 <- nchar(Seq$Sequence)
  #for (C in gsub("_", "", Cut)) { Seq[[paste0(C, ".Count")]] <- n0 - nchar(gsub(C, "", Seq$Sequence)) }
  if (missing("missed")) {
    n1 <- nchar(gsub(paste(paste0(gsub("_", "", Cut), "$"), collapse = "|"), "", Seq$Sequence))
    n2 <- nchar(gsub(paste(gsub("_", "", Cut), collapse = "|"), "", Seq$Sequence))
    Seq$"Missed cleavages" <- n1 - n2
    missed <- max(Seq$"Missed cleavages")
  }
  # Identify matching proteins:
  # Digest database and match to predicted peptides:
  dig <- proteoCraft::Digest(setNames(DB$Sequence, DB$`Protein ID`),
                             Cut,
                             strict.avoid,
                             loose.avoid,
                             missed,
                             min,
                             max,
                             FALSE)
  lTst <- vapply(dig, length, 1)
  w <- which(lTst == 0)
  lW <- length(w)
  if (lW) {
    warning(paste0(lW, " protein digest", c("", "s")[(lW > 1)+1], " did not result in any peptide sequences, investigate!!!"))
    dig <- dig[which(lTst > 0)]
  }
  temp2 <- magrittr::set_colnames(reshape2::melt(dig), c("Pep", "Prot"))
  #temp2$Sequence <- paste0("_", DB$Sequence[match(temp2$Prot, DB$`Protein ID`)], "_")
  #tst <- vapply(1:nrow(temp2), function(x) { length(unlist(strsplit(temp2$Sequence[x], temp2$Pep[x]))) }, 1)
  #min(tst) > 2
  #temp2a <- temp2[which(temp2$Pep %in% Seq$Sequence),]
  #temp2b <- temp2[which(gsub("I", "L", temp2$Pep) %in% gsub("I", "L", Seq$Sequence)),] # There is a peptides assignment loss, though less than 1%
  temp2 <- temp2[which(gsub("I", "L", temp2$Pep) %in% gsub("I", "L", Seq$Sequence)),]
  # NB: The "unique" in the call below should be required for cases where we have sequence repetitions.
  temp3 <- magrittr::set_colnames(aggregate(temp2$Prot, list(temp2$Pep), function(x) { paste(unique(x), collapse = ";") }), c("Pep", "Prot"))
  w1 <- which(Seq$Sequence %in% temp3$Pep)
  w2 <- which(!Seq$Sequence %in% temp3$Pep)
  if (length(w2) > 0) { warning(paste0(length(w2), " observed peptides are not found in the canonical digest. Investigations of previous cases showed that rare reported peptides are non-canonical (e.g. internal breaks?), probably when spectra allow unambiguous identification(?)")) }
  Seq$Proteins <- ""
  Seq$Proteins[w1] <- temp3$Prot[match(Seq$Sequence[w1], temp3$Pep)]
  return(Seq)
}
