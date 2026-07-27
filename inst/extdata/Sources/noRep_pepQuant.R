# No replicates script:
# Calculate per-sample peptide intensities from PSM intensities
for (e in Exp) { #e <- Exp[1L]
  temp <- ev[which(ev$Experiment == e),]
  temp <-  magrittr::set_colnames(aggregate(temp[[int.col]], list(temp$"Modified sequence"), \(x) { sum(x[which(is.finite(x))]) }),
                                  c("Modified sequence", int.col))
  pep[, paste0(int.col, " - ", e)] <- 0
  w <- which(pep$"Modified sequence" %in% temp$"Modified sequence")
  pep[w, paste0(int.col, " - ", e)] <- temp[match(pep$"Modified sequence"[w], temp$"Modified sequence"), int.col]
}
ev_to_pep <- match(pep$`Modified sequence`, ev$`Modified sequence`)
if (!exists("PTMriched")) { PTMriched <- FALSE }
if (PTMriched) {
  pep[, EnrichedPTMs] <- ev[ev_to_pep, EnrichedPTMs]
}
