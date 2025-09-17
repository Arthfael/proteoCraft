# Annotate the protein groups table
p <- strsplit(PG$"Leading protein IDs", ";") #Here taking just the minimum set of protein IDs to explain the observed dataset.
db$Observed <- db$"Protein ID" %in% unique(unlist(p))
if (globalGO) {
  temp <- listMelt(strsplit(PG$"Leading protein IDs", ";"), PG$id, c("Accession", "id"))
  kol <- annot.col[which(annot.col %in% colnames(db))]
  if ("Taxonomy" %in% kol) { # Taxonomy can be dealt with differently
    PG$Taxonomy <- db$Taxonomy[match(gsub(";.*", "", PG$`Leading protein IDs`), db$`Protein ID`)]
  }
  kol2 <- annot.col[which(!annot.col %in% "Taxonomy")]
  kol2 <- kol2[which(kol2 %in% colnames(db))]
  temp[, kol2] <- db[match(temp$Accession, db$"Protein ID"), kol2]
  tst1 <- unlist(strsplit(temp$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(temp$GO, ";"))
  tst2a <- gsub(".*\\[|\\]$", "", tst2)
  # Test that...
  tst3 <- data.table(ID = tst1, Name = tst2)
  # - ... each ID has a single name, and that names are unique
  tst3a <- tst3[, list(Name = list(unique(Name))), by = list(ID = ID)]
  tst3a <- as.data.frame(tst3a)
  tst3a$N_names <- vapply(tst3a$Name, length, 1)
  # - ... each name is unique to one ID
  tst3b <- tst3[, list(ID = list(unique(ID))), by = list(Name = Name)]
  tst3b <- as.data.frame(tst3b)
  tst3b$N_IDs <- vapply(tst3b$ID, length, 1)
  stopifnot(max(tst3b$N_IDs) == 1)
  #
  for (i in kol2) { temp[[i]] <- strsplit(as.character(temp[[i]]), ";") }
  f0 <- function(x) { list(unique(unlist(x))) }
  temp <- aggregate(temp[, kol2], list(temp$id), f0)
  #
  # Below commented data.table aggregation code... which is slower so not used.
  #
  # temp2 <- as.data.table(temp[, c("id", kol2)])
  # temp2 <- temp2[, lapply(.SD, f0), by = list(Group.1 = id), .SDcols = kol2]
  # temp2 <- as.data.frame(temp2)
  #
  for (i in kol2) {
    temp[[i]] <- parSapply(parClust, temp[[i]], function(x) { paste(unique(unlist(x)), collapse = ";") }) # Do not use sort here or it will break the matching between "GO" and "GO-ID"
  }
  tst1 <- unlist(strsplit(temp$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(temp$GO, ";"))
  tst3 <- data.table(ID = tst1, Name = tst2)
  tst3 <- tst3[, list(Name = list(unique(Name))), by = list(ID = ID)]
  tst3 <- as.data.frame(tst3)
  tst3b$N_IDs <- vapply(tst3b$ID, length, 1)
  stopifnot(max(tst3b$N_IDs) == 1)
  #
  PG[, kol2] <- ""
  w <- which(PG$id %in% temp$Group.1)
  PG[w, kol2] <- temp[match(PG$id[w], temp$Group.1), kol2]
  #
  # Also peptides (minor approximation: use first protein group)
  pep[, kol] <- PG[match(as.integer(gsub(";.*", "", pep$`Protein group ID`)), PG$id), kol]
  #
  PG$Ontology <- NULL # Temporary fix for now, this column is broken
  #
  stopCluster(parClust)
  source(parSrc, local = FALSE)
  Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_prepare.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}
