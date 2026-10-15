# Annotate our data
setwd(wd)
if (!Annotate) { globalGO <- FALSE }
p <- unique(unlist(strsplit(pep$Proteins, ";"))) #Here taking just the minimum set of protein IDs to explain the observed dataset.
db$Observed <- db$"Protein ID" %in% p
if (globalGO) {
  if (scrptTypeFull == "withReps_PTMs_only") {
    myData <- PG
    myCol <- "Leading protein IDs"
    # NB: I used to get functional annotations for all proteins in the protein group.
    #     However we are now - and I think with reason - only using annotations from the leading protein(s)!
  } else {
    myData <- pep
    myCol <- "Proteins"
  }
  #
  Src <- paste0(libPath, "/extdata/Sources/GO_prepare.R")
  #rstudioapi::documentOpen(Src)
  source(Src)
  #
  temp <- listMelt(strsplit(myData[[myCol]], ";"), myData$id, c("Accession", "id"))
  kol <- intersect(annot.col, colnames(db))
  if ("Taxonomy" %in% kol) { # Taxonomy can be dealt with differently
    myData$Taxonomy <- db$Taxonomy[match(sub(";.*", "", myData[[myCol]]), db$`Protein ID`)]
  }
  kol2 <- setdiff(annot.col, "Taxonomy")
  kol2 <- intersect(kol2, colnames(db))
  temp[, kol2] <- db[match(temp$Accession, db$"Protein ID"), kol2]
  tst1 <- unlist(strsplit(temp$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(temp$GO, ";"))
  tst2a <- gsub(".*\\[|\\]$", "", tst2)
  # Test that...
  tst3 <- data.table(ID = tst1, Name = tst2)
  # - ... each ID has a single name, and that names are unique
  tst3a <- tst3[, .(Name = list(unique(Name))), by = .(ID = ID)]
  tst3a <- as.data.frame(tst3a)
  tst3a$N_names <- lengths(tst3a$Name)
  # - ... each name is unique to one ID
  tst3b <- tst3[, .(ID = list(unique(ID))), by = .(Name = Name)]
  tst3b <- as.data.frame(tst3b)
  tst3b$N_IDs <- lengths(tst3b$ID)
  stopifnot(max(tst3b$N_IDs) == 1L)
  #
  for (i in kol2) { temp[[i]] <- strsplit(as.character(temp[[i]]), ";") }
  f0 <- \(x) { list(unique(unlist(x))) }
  temp <- aggregate(temp[, kol2], list(temp$id), f0)
  #
  # Below commented data.table aggregation code... which is slower so not used.
  #
  # temp2 <- as.data.table(temp[, c("id", kol2)])
  # temp2 <- temp2[, lapply(.SD, f0), by = .(Group.1 = id), .SDcols = kol2]
  # temp2 <- as.data.frame(temp2)
  #
  for (i in kol2) {
    temp[[i]] <- parSapply(parClust, temp[[i]], \(x) { paste(unique(unlist(x)), collapse = ";") }) # Do not use sort here or it will break the matching between "GO" and "GO-ID"
  }
  tst1 <- unlist(strsplit(temp$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(temp$GO, ";"))
  tst3 <- data.table(ID = tst1, Name = tst2)
  tst3 <- tst3[, .(Name = list(unique(Name))), by = .(ID = ID)]
  tst3 <- as.data.frame(tst3)
  tst3b$N_IDs <- lengths(tst3b$ID)
  stopifnot(max(tst3b$N_IDs) == 1L)
  #
  myData[, kol2] <- ""
  w <- which(myData$id %in% temp$Group.1)
  myData[w, kol2] <- temp[match(myData$id[w], temp$Group.1), kol2]
  #
  myData$Ontology <- NULL # Temporary fix for now, this column is broken
  if (scrptTypeFull != "withReps_PTMs_only") {
    pep <- myData
  } else {
    PG <- myData
    # Also peptides (minor approximation: use first protein group)
    pep[, kol] <- PG[match(as.integer(sub(";.*", "", pep$`Protein group ID`)), PG$id), kol]
  }
  #
  stopCluster(parClust)
  source(parSrc)
}
