### Load and process annotations
dfltLocsFl <- paste0(homePath, "/Default_locations.xlsx")
dfltLocs <- openxlsx2::read_xlsx(dfltLocsFl)
fastaLoc <- dfltLocs$Path[match("Fasta files", dfltLocs$Folder)]
#
GO.col %<o% c("GO", "GO-ID")
if (Annotate) {
  source(parSrc, local = FALSE)
  Parsed_annotations_lst <- lapply(1:nrow(AnnotFlsTbl), \(i) { #i <- 1L
    fl <- AnnotFlsTbl$Path[i]
    tp <- AnnotFlsTbl$Type[i]
    # If the annotation file is not present locally, make a local copy
    if (!file.exists(paste0(wd, "/", basename(fl)))) { fs::file_copy(fl, wd) }
    # Now parse it
    if (tp == "UniProtKB .txt") {
      x <- Format.DB_txt(fl, usePar = TRUE, cl = parClust)
    }
    if (tp == "NCBI .gtf") {
      x <- annot_from_GTF(fl, mode = "GTF")
    }
    if (tp == "NCBI .gff") {
      x <- annot_from_GTF(fl, mode = "GFF")
    }
    return(x)
  })
  Parsed_annotations <- dplyr::bind_rows(Parsed_annotations_lst)
}
Annotate <- ((exists("Parsed_annotations"))
             &&(is.data.frame(Parsed_annotations))
             &&(nrow(Parsed_annotations) > 0L))
if (Annotate) {
  if ("Sequence" %in% colnames(Parsed_annotations)) {
    mType <- "Sequence"
    mAnnot <- match(db$Sequence, Parsed_annotations$Sequence)
    Annotate <- sum(!is.na(mAnnot)) > 0L
  } else {
    mType <- "Accession"
    mAnnot <- match(db$`Protein ID`, Parsed_annotations$Accession)
    Annotate <- sum(!is.na(mAnnot)) > 0L
  }
}
if (Annotate) {
  Parsed_annotations %<o% Parsed_annotations
  #
  # Check GO for name degeneracies
  # (the same term has had different names in different files, and we allow for multiple files)
  # Match Parsed_annotations rows to GO name(s)
  wAnnot <- which(nchar(Parsed_annotations$GO) > 0L)
  tst1 <- listMelt(strsplit(Parsed_annotations$GO[wAnnot], ";"),
                   wAnnot,
                   c("name", "row"))[, c("row", "name")]
  # Match GO name/IDs to Parsed_annotations row(s)
  tst2 <- set_colnames(aggregate(tst1$row, list(tst1$name), list), c("name", "rows"))
  tst2$ID <- gsub(".*\\[|\\]$", "", tst2$name)
  tst1$ID <- tst2$ID[match(tst1$name, tst2$name)]
  # Match GO IDs to their name(s)
  tst3 <- set_colnames(aggregate(tst2$name, list(tst2$ID), unique), c("ID", "names"))
  tst3$L <- lengths(tst3$names)
  wDegen <- which(tst3$L > 1L)
  if (length(wDegen)) { # Indicates degeneracy and the need to fix names
    degen <- tst3[which(tst3$L > 1L),] # Degenerate ID-to-names mappings
    degen$name <- vapply(degen$names, \(x) { x[[1L]] }, "") # Take the first name
    rws1 <- unique(unlist(tst2$rows[which(tst2$ID %in% degen$ID)])) # Rows where the degenerate IDs are present
    tst1 <- tst1[which(tst1$row %in% rws1),] # Keep ALL annotations for those rows
    tst2 <- tst2[which(tst2$ID %in% degen$ID),]
    tst2$name_from_3 <- degen$name[match(tst2$ID, degen$ID)] # Get replacement names
    #View(tst2[which(tst2$name != tst2$name_from_3),]) # Check before updating
    # If those are ok, replace:
    tst2$name <- tst2$name_from_3
    tst2$name_from_3 <- NULL
    #
    #tst3b <- set_colnames(aggregate(tst2$name, list(tst2$ID), unique), c("ID", "names"))
    #tst3b$L <- lengths(tst3b$names)
    #max(tst3b$L) # Should be 1 now
    #
    # Apply new names to tst1
    w2 <- which(tst1$ID %in% tst2$ID)
    tst1$name[w2] <- tst2$name[match(tst1$ID[w2], tst2$ID)]
    #tst3c <- set_colnames(aggregate(tst1$name, list(tst1$ID), unique), c("ID", "names"))
    #tst3c$L <- lengths(tst3c$names)
    #max(tst3c$L) # Should be 1 now
    #
    tst1 <- as.data.table(tst1)
    tst1 <- tst1[, list(GO = paste(name, collapse = ";"),
                        ID = paste(ID, collapse = ";")),
                 by = list(row = tst1$row)]
    tst1 <- as.data.frame(tst1)
    # Check results
    a1 <- unlist(strsplit(tst1$ID, ";"))
    b1 <- unlist(strsplit(tst1$GO, ";"))
    c1 <- gsub(".*\\[|\\]$", "", b1)
    stopifnot(length(a1) == length(b1),
              sum(a1 != c1) == 0L)
    #
    w <- which(1L:nrow(Parsed_annotations) %in% tst1$row)
    Parsed_annotations[w, c("GO", "GO-ID")] <- tst1[match(w, tst1$row), c("GO", "ID")]
    #View(tst1[match(w, tst1$tst1), c("GO", "ID")])
  }
  tst1 <- unlist(strsplit(Parsed_annotations$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(Parsed_annotations$GO, ";"))
  tst3 <- data.table(ID = tst1, Name = tst2)
  tst3 <- tst3[, list(Name = list(unique(Name))), by = list(ID = ID)]
  tst3 <- as.data.frame(tst3)
  tst3$N_names <- lengths(tst3$Name)
  #max(tst3$N_names)
  stopifnot(max(tst3$N_names) == 1L)
  #View(tst3[which(tst3$N_names > 1L),])
  #View(tst3[which(tst3$N_names == 0L),])
  #View(Parsed_annotations[, c("GO", "GO-ID")])
  #db <- db[, which(!colnames(db) %in% annot.col)]
  kol <- colnames(Parsed_annotations)
  annot.col %<o% kol[which(!kol %in% c("Accession", "id", "ID", "Names", "Sequence", "MW (Da)", "Common Name"))]
  annot.col2 %<o% annot.col[which(!annot.col %in% colnames(db))]
  annot.col3 %<o% annot.col[which(annot.col %in% colnames(db))]
  if (length(annot.col2)) { db[, annot.col2] <- NA }
  w <- which(db$`Protein ID` %in% Parsed_annotations$Accession)
  if (mType == "Sequence") {
    mtch <- match(db$Sequence[w], Parsed_annotations$Sequence)
  }
  if (mType == "Accession") {
    mtch <- match(db$`Protein ID`[w], Parsed_annotations$Accession)
  }
  db[, annot.col] <- ""
  db[w, annot.col] <- Parsed_annotations[mtch, annot.col]
  if ("Common Name" %in% kol) {
    w2 <- which(!nchar(db$"Common Name"))
    if (length(w2)) {
      if (mType == "Sequence") {
        mtch2 <- match(db$Sequence[w2], Parsed_annotations$Sequence)
      }
      if (mType == "Accession") {
        mtch2 <- match(db$`Protein ID`[w2], Parsed_annotations$Accession)
      }
      db$"Common Name"[w] <- Parsed_annotations$"Common Name"[mtch2]
    }
  }
  #View(db)
  tst1 <- unlist(strsplit(db$`GO-ID`, ";"))
  tst4 <- unique(tst1)
  tst4 <- tst4[which(nchar(tst4) > 0L)]
  Annotate <- length(tst4) > 0L
}
if (Annotate) {
  tst2 <- unlist(strsplit(db$GO, ";"))
  tst3 <- data.table(ID = tst1, Name = tst2)
  tst3 <- tst3[, list(Name = list(unique(Name))), by = list(ID = ID)]
  tst3 <- as.data.frame(tst3)
  tst3$N_names <- lengths(tst3$Name)
  #max(tst3$N_names)
  stopifnot(max(tst3$N_names) == 1L)
  if (length(annot.col3)) {
    for (kol in annot.col3) { #kol <- annot.col3[1L]
      w <- which(!is.na(mtch)) # check that there is a valid match...
      w <- w[which((is.na(db[w, kol]))|(db[w, kol] %in% c("", "NA", "NaN")))] #... and that it is useful!
      db[w, kol] <- Parsed_annotations[mtch[w], kol]
      w <- which((is.na(db[[kol]]))|(db[[kol]] %in% c("", "NA", "NaN"))) #... and that it is useful!
      db[w, kol] <- ""
    }
  }
  tst1 <- unlist(strsplit(db$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(db$GO, ";"))
  tst3 <- data.table(A1 = tst1, A2 = tst2)
  tst3 <- tst3[, list(x = unique(A2)), by = list(Group.1 = A1)]
  tst3 <- as.data.frame(tst3)
  #View(tst3[which(lengths(tst3$x) > 1L),])
  #View(tst3[which(lengths(tst3$x) == 0L),])
  stopifnot(length(tst3$x) == length(unique(tst3$x)),
            is.character(tst3$x))
  db$Ontology <- NULL # Temporary fix for now, this column is broken
  #
  w <- which(vapply(colnames(db), \(x) { inherits(db[[x]], "list") }, TRUE))
  if (length(w)) {
    for (i in w) {
      db[[i]] <- parSapply(parClust, db[[i]], paste, collapse = ";")
    }
  }
  saveFun(Parsed_annotations, file = "Parsed_annotations.RData")
  pth <- paste0(wd, "/Parsed, annotated search db.csv")
  f0 <- \() {
    data.table::fwrite(db, pth,
                       quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE, na = "NA")
  }
  tst <- try(f0(), silent = TRUE)
  while ((inherits(tst, "try-error"))&&(grepl("cannot open the connection", tst[1L]))) {
    dlg_message(paste0("File \"", pth, "\" appears to be locked for editing, close the file then click ok..."), "ok")
    tst <- try(f0(), silent = TRUE)
  }
  #db <- data.table::fread(paste0(wd, "/Parsed, annotated search db.csv"), integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  #write.csv(db, "Parsed, annotated search db.csv", row.names = FALSE)
  #db <- read.csv("Parsed, annotated search db.csv", check.names = FALSE, colClasses = "character")
}
