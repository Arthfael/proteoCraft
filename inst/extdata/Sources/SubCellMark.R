if (Annotate) {
  require(pRoloc)
  # Define subcellular "Markers" useable by various steps of the workflow
  CompGOTerms %<o% setNames(c("GO:0005634", "GO:0005654", "GO:0000785", "GO:0005730", "GO:0005635", "GO:0005737", "GO:0005829",
                              "GO:0005783", "GO:0005794", "GO:0031988", "GO:0005739", "GO:0009536", "GO:0005886"),
                            c("Nucleus", "Nucleoplasm", "Chromatin", "Nucleolus", "Nuclear envelope", "Cytoplasm", "Cytosol",
                              "ER", "Golgi", "Vesicle", "Mitochondrion", "Plastid", "Plasma membrane"))
  tstSp <- if ((exists("Org"))
      &&(is.data.frame(Org))
      &&(nrow(Org) >= 1L)
      &&("Organism" %in% colnames(Org))
      &&(!is.na(Org$Organism[1L]))) {
    Org$Organism[1L]
  } else {
    unlist(strsplit(Param$Search.DB.species, ";"))[1L]
  }
  tstSp <- tolower(unlist(strsplit(tstSp, " ")))
  tstSp <- paste0(substr(tstSp[1L], 1L, 1L), substr(tstSp[2L], 1L, 3L))
  tst2 <- capture.output(pRolocmarkers())
  tst2 <- tst2[2L:length(tst2)]
  l2 <- length(tst2)
  tst2 <- data.frame(Species = tst2[c(1L:(l2/2L))*2L-1L],
                     IDsType = tst2[c(1L:(l2/2L))*2L])
  tst2$SpTag <- gsub("^[^\\[]+\\[|(\\]|_).*:$", "", tst2$Species)
  useProloc <- (tstSp %in% tst2$SpTag)
  if (useProloc) { # Use pRoloc markers if our species is available in pRoloc...
    SubCellMark <- pRolocmarkers(tstSp)
    tst <- data.frame(Comp = SubCellMark,
                      Mark = names(SubCellMark) %in% db$`Protein ID`)
    tst <- aggregate(tst$Mark, list(tst$Comp), sum)
    tst <- tst[which(tst$x > 3L),]
    useProloc <- nrow(tst) >= 5L
  }
  if (!useProloc) { # ... or generate them automatically
    tst <- setNames(lapply(CompGOTerms, \(x) {
      unique(unlist(strsplit(PG$`Leading protein IDs`[grsep(x, x = PG$"GO-ID")], ";")))
    }), names(CompGOTerms))
    tst <- listMelt(tst, ColNames = c("Accession", "Compartment"))
    tst2 <- aggregate(tst$Compartment, list(tst$Accession), \(x) {
      length(unique(x))
    })
    tst <- tst[which(tst$Accession %in% tst2$Group.1[which(tst2$x == 1L)]),]
    SubCellMark <- setNames(tst$Compartment, tst$Accession)
  }
  SubCellMark %<o% SubCellMark
  #
  # Annotate PG
  if ((exists("PG"))&&(length(SubCellMark))) {
    tst <- listMelt(strsplit(PG$`Leading protein IDs`, ";"), ColNames = c("Accession", "Row"))
    tst <- tst[which(tst$Accession %in% names(SubCellMark)),]
    tst$Comp <- SubCellMark[match(tst$Accession, names(SubCellMark))]
    tst <- as.data.table(tst)
    tst <- tst[, .(Comp = list(unique(Comp))), by = .(Row = as.integer(Row))]
    tst <- as.data.frame(tst)
    tst$L <- lengths(tst$Comp)
    tst <- tst[which(tst$L == 1L),]
    tst$Comp <- vapply(tst$Comp, unlist, "")
    PG$"Compartment marker" <- ""
    PG$"Compartment marker"[tst$Row] <- tst$Comp
  }
}
