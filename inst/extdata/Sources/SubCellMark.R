if (Annotate) {
  
  # If scrptTypeFull == "custom", you must have created the following variables:
  # - myData = the data to be annotation, which must have a "GO-ID" column!
  # - myCol = the column with protein IDs in myData
  # - useProloc = logical, use pRoloc markers automatically or an alternative solution?
  # before running this source...
  #
  if (scrptTypeFull == "withReps_PTMs_only") {
    myData <- PG
    myCol <- "Leading protein IDs"
    # NB: I used to get functional annotations for all proteins in the protein group.
    #     However we are now - and I think with reason - only using annotations from the leading protein(s)!
  } else if (scrptTypeFull != "custom") {
    myData <- pep
    myCol <- "Proteins"
  }
  #
  # Define subcellular "Markers" useable by various steps of the workflow
  CompGOTerms %<o% setNames(c("GO:0005634", "GO:0005654", "GO:0000785", "GO:0005730", "GO:0005635", "GO:0005737", "GO:0005829",
                              "GO:0005783", "GO:0005794", "GO:0031988", "GO:0005739", "GO:0009536", "GO:0005886"),
                            c("Nucleus", "Nucleoplasm", "Chromatin", "Nucleolus", "Nuclear envelope", "Cytoplasm", "Cytosol",
                              "ER", "Golgi", "Vesicle", "Mitochondrion", "Plastid", "Plasma membrane"))
  if (scrptTypeFull != "custom") {
    tstSp <- if (exists("Org")
                 && is.data.frame(Org)
                 && (nrow(Org) >= 1L)
                 && ("Organism" %in% colnames(Org))
                 && (!is.na(Org$Organism[1L]))) {
      Org$Organism[1L]
    } else {
      if (exists("Param") && ("Search.DB.species" %in% colnames(Param))) {
        unlist(strsplit(Param$Search.DB.species, ";"))[1L]
      } else {
        svDialogs::dlg_input("Enter the dataset's main species (Linnean name, e.g. 'Homo sapiens')")$res
      }
    }
    tstSp <- tolower(unlist(strsplit(tstSp, " ")))
    tstSp <- paste0(substr(tstSp[1L], 1L, 1L), substr(tstSp[2L], 1L, 3L))
    tst2 <- capture.output(pRoloc::pRolocmarkers())
    tst2 <- tst2[2L:length(tst2)]
    l2 <- length(tst2)
    tst2 <- data.frame(Species = tst2[c(1L:(l2/2L))*2L-1L],
                       IDsType = tst2[c(1L:(l2/2L))*2L])
    tst2$SpTag <- gsub("^[^[]+\\[|(\\]|_).*:$", "", tst2$Species)
    useProloc <- (tstSp %in% tst2$SpTag)
    if (useProloc) { # Use pRoloc markers if our species is available in pRoloc...
      SubCellMark <- pRoloc::pRolocmarkers(tstSp)
      tst <- data.frame(Comp = SubCellMark,
                        Mark = names(SubCellMark) %in% db$`Protein ID`)
      tst <- aggregate(tst$Mark, list(tst$Comp), sum)
      tst <- tst[tst$x > 3L,]
      useProloc <- nrow(tst) >= 5L
    }
  }
  if (!useProloc) { # ... or generate them automatically
    tst <- setNames(lapply(CompGOTerms, \(x) {
      unique(unlist(strsplit(myData[grsep(x, x = myData$"GO-ID"), myCol], ";")))
    }), names(CompGOTerms))
    tst <- listMelt(tst, ColNames = c("Accession", "Compartment"))
    tst2 <- aggregate(tst$Compartment, list(tst$Accession), \(x) {
      length(unique(x))
    })
    tst <- tst[tst$Accession %in% tst2$Group.1[tst2$x == 1L],]
    SubCellMark <- setNames(tst$Compartment, tst$Accession)
  }
  SubCellMark %<o% SubCellMark
  #
  # Annotate myData
  if (length(SubCellMark)) {
    tst <- listMelt(strsplit(myData[[myCol]], ";"), ColNames = c("Accession", "Row"))
    tst <- tst[tst$Accession %in% names(SubCellMark),]
    tst$Comp <- SubCellMark[match(tst$Accession, names(SubCellMark))]
    requireNamespace("data.table")
    tst <- data.table::as.data.table(tst)
    tst <- tst[, .(Comp = list(unique(Comp))), by = .(Row = as.integer(Row))]
    tst <- as.data.frame(tst)
    tst$L <- lengths(tst$Comp)
    tst <- tst[tst$L == 1L,]
    tst$Comp <- vapply(tst$Comp, unlist, "")
    myData$"Compartment marker" <- ""
    myData$"Compartment marker"[tst$Row] <- tst$Comp
  }
  #
  if (scrptTypeFull == "withReps_PTMs_only") {
    PG <- myData
    # Also peptides (minor approximation: use first protein group)
    pep[, kol] <- PG[match(as.integer(sub(";.*", "", pep$`Protein group ID`)), PG$id), kol]
  } else if (scrptTypeFull != "custom") {
    pep <- myData
  }
}
