# No replicates workflow only:
# Calculate per-sample peptide intensities from PSM intensities
for (e in Exp) { #e <- Exp[1L]
  temp <- ev[ev$Experiment == e,]
  temp <- aggregate(temp[[int.col]], list(temp$"Modified sequence"), \(x) { sum(x[is.finite(x)]) })
  colnames(temp) <- c("Modified sequence", int.col)
  k <- paste0(int.col, " - ", e)
  pep[, k] <- 0
  w <- which(pep$"Modified sequence" %in% temp$"Modified sequence")
  pep[w, k] <- temp[match(pep$"Modified sequence"[w], temp$"Modified sequence"), int.col]
}
ev_to_pep <- match(pep$`Modified sequence`, ev$`Modified sequence`)
if (!exists("PTMriched")) { PTMriched <- FALSE }
if (PTMriched) {
  pep[, EnrichedPTMs] <- ev[ev_to_pep, EnrichedPTMs]
}
