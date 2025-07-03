# Create experiment Factors shortcuts
Aggregates %<o% Factors
a <- substr(Aggregates, 1, 3)
if (length(unique(a)) != length(a)) {
  stop("Factors must start with a unique 3 letter tag! Edit the experiment map columns and restart.")
} else { names(Aggregates) <- a }
for (i in Aggregates) {
  temp <- sort(unique(Exp.map[[i]]))
  tst <- sum(as.character(suppressWarnings(as.numeric(temp))) != temp) # Are they numerics?
  if ((!is.na(tst))&&(!tst)) { temp <- as.numeric(temp) }
  tst <- sum(as.character(suppressWarnings(as.integer(temp))) != temp) # Are they integers?
  if ((!is.na(tst))&&(!tst)) { temp <- as.integer(temp) }
  substr(i, 1, 3) %<c% temp
}
if (exists("Rep")) { names(Rep) <- gsub("rep", "", Rep) }
if (LabelType == "Isobaric") {
  if (!exists("Iso")) {
    stop("See older versions of this script - but in the newer version I would always expect Iso to exist here.")
  } else { test.iso.set %<o% TRUE }
}
MQ.Frac %<o% sort(suppressWarnings(unique(as.integer(unlist(strsplit(as.character(Exp.map$Fractions), ";"))))), na.last = TRUE)

# Test
if ((SearchSoft == "MAXQUANT")&&(LabelType == "Isobaric")) {
  tst <- sum(!paste0("Reporter intensity ", Exp.map$"Isobaric label") %in% colnames(ev))
  stopifnot(tst == 0)
  # If this is not correct, then we should consider re-introducing older "Code chunk - Isobarically-labelled samples only! Check whether 1 must be subtracted from MaxQuant's isobaric labels index"
}
