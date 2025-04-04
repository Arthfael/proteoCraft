#' ProtFormula
#' 
#' @description
#' Calculate the chemical formula of a protein.
#' 
#' @param Seq Protein single letter code sequence.
#' @param collapse If TRUE, collapse as single character chain, otherwise provided as table of individual element counts.
#' @param reduced If TRUE, then disulfide bonds are assumed to be reduced. Note: We are assuming that all cysteines are parts of disulfide bonds... but we do not know that! One solution would be to provide a new way to map known disulfide bonds to the sequence. For now the function cannot handle this case other than by assuming that all cysteines are oxidized.
#' @param alkylated If one of "IAM", "NEM", "CAM", "MMTS" or "iST_NHS" (alkylating agent used in the iST-NHS kit), then the appropriate modification will be applied on cysteines.
#' 
#' @examples
#' Form <- ProtFormula(Seq, collapse = FALSE, reduced = FALSE)
#' Form <- ProtFormula(Seq, collapse = FALSE, alkylated = FALSE)
#' Form <- ProtFormula(Seq, collapse = FALSE)
#' 
#' @export

ProtFormula <- function(Seq, collapse = TRUE, reduced = TRUE, alkylated = "IAM") {
  Seq <- unlist(strsplit(toupper(Seq), ""))
  L <- length(Seq)
  Tbl <- proteoCraft::AA_table
  alk <- (!is.logical(alkylated))&&(!is.na(alkylated))&&(nchar(as.character(alkylated)) > 0)
  if ((!reduced)&&(alk)) {
    warning("Did you really alkylate but did not reduce? Ignoring alkylations.")
    alk <- FALSE
  }
  if ((alk)&&(alkylated %in% proteoCraft::alkTbl$AA)) {
    warning(paste0("\"", alkylated, "\" is not recognized as a valid alkylation and will be ignored!"))
  }
  wC <- which(Tbl$AA == "C")
  if (!reduced) {
    warning("We are assuming that all cysteines are parts of disulfide bonds... but we do not know that! One solution would be to provide a list of known disulfide bonds, if available. For now the function cannot handle this case other than by assuming that all cysteines are oxidized.")
    # Create mass shift from disulfide bonds
    Tbl$H[wC] <- Tbl$H[wC]-1
  }
  if (alk) {
    w <- which(proteoCraft::alkTbl$Alk == alkylated)
    for (at in grep("^[A-Z][a-z]?$", colnames(Tbl), value = TRUE)) {
      if (at %in% colnames(proteoCraft::alkTbl)) { Tbl[wC, at] <- Tbl[wC, at] + proteoCraft::alkTbl[w, at] }
    }
  }
  w <- which(!Seq %in% Tbl$AA)
  if (length(w)) { warning(paste0("Bypassing ", w, " un-recognised amino acid",
                                  c("", "s")[(length(w) > 1)+1], "!")) }
  Seq <- Seq[which(Seq %in% Tbl$AA)]
  atoms <- c("C", "H", "O", "N", "S", "Se") 
  Seq <- unlist(lapply(Seq, function(x) { #x <- Seq[1]
    m <- match(x, Tbl$AA)
    unlist(sapply(atoms, function(y) { rep(y, Tbl[m, y]) }))
  }))
  Seq <- aggregate(Seq, list(Seq), length)
  colnames(Seq) <- c("Element", "Count")
  w <- which(atoms %in% Seq$Element)
  Seq <- Seq[match(atoms[w], Seq$Element),]
  # Effect of peptide bonds:
  Seq$Count[which(Seq$Element == "O")] <- Seq$Count[which(Seq$Element == "O")]-(L-1)
  Seq$Count[which(Seq$Element == "H")] <- Seq$Count[which(Seq$Element == "H")]-(L-1)*2
  #
  if (collapse) {
    Seq$Count <- as.character(Seq$Count)
    Seq <- paste(do.call(paste, c(Seq[, c("Element", "Count")], sep = "")), collapse = "")
  }
  return(Seq)
}
