#' phos_QC
#'
#' @description 
#' A function to create a column saying whether a phospho-peptide is high quality or not.
#' Useful for creating a table of, e.g. sequence with amino-acid modifications, site probabilities, etc...
#' 
#' @param PhPep The phospho-peptides table.
#' @param min_P The minimum site probability acceptable. Default = 0.75
#' @param min_delta_sc The minimum delta score acceptable. Default = 5
#' @param max_PEP The maximum Posterior Error Probability acceptable. Default = 0.01
#' @param ModSeq_col The name of the modified sequence column. Default = "Modified sequence"
#' @param P_col The name of the site probabilities column. Default = "Phospho (STY) Probabilities"
#' @param ScD_col The name of the score differences column. Default = "Phospho (STY) Score Diffs"
#' @param PEP_col The name of the Posterior Error Probability column. Default = "PEP"
#' 
#' @examples
#' phospep$High_Quality <- phos_QC(phospep)
#' 
#' @export

phos_QC <- function(PhPep,
                    min_P = 0.75,
                    min_delta_sc = 5,
                    max_PEP = 0.01,
                    ModSeq_col = "Modified sequence",
                    P_col = "Phospho (STY) Probabilities",
                    ScD_col = "Phospho (STY) Score Diffs",
                    PEP_col = "PEP") {
  if (!ModSeq_col %in% colnames(PhPep)) { stop("The modified sequence column was not found!") }
  if (!P_col %in% colnames(PhPep)) { stop("The site probabilities column was not found!") }
  if (!ScD_col %in% colnames(PhPep)) { stop("The score differences column was not found!") }
  l <- nrow(PhPep)
  if (l == 0) { stop("There are no rows to analyse!") }
  M <- proteoCraft::annot_to_tabl(PhPep[[ModSeq_col]])
  P <- proteoCraft::annot_to_tabl(PhPep[[P_col]], Nterm = FALSE, Cterm = FALSE, numeric_data = TRUE)
  ScD <- proteoCraft::annot_to_tabl(PhPep[[ScD_col]], Nterm = FALSE, Cterm = FALSE, numeric_data = TRUE)
  res <- sapply(c(1:l), function(x) {
    m <- M[[x]]
    p <- P[[x]]
    scd <- ScD[[x]]
    if (class(m) == "data.frame") {
      stopifnot(class(p) == "data.frame", class(scd) == "data.frame")
      m$P <- as.numeric(p$Annotations)
      m$ScD <- as.numeric(scd$Annotations)
      m <- m[2:(nrow(m)-1),]
      w <- which(m$Annotations == "ph")
      test <- sapply(w, function(y) {
        (m$P[[y]] >= min_P)&(m$ScD[[y]] >= min_delta_sc)
      })
      r <- (sum(test) == length(w))&(PhPep[x,PEP_col] <= max_PEP)
    } else { r <- FALSE }
    return(r)
  })
  return(res)
}
