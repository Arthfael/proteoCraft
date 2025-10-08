### Load and process annotations
# This includes a QC step in case the database differs slightly from the one used by MQ, or if somehow some IDs have not been properly parsed.
dfltLocsFl <- paste0(homePath, "/Default_locations.xlsx")
dfltLocs <- openxlsx2::read_xlsx(dfltLocsFl)
fastaLoc <- dfltLocs$Path[match("Fasta files", dfltLocs$Folder)]
#
GO.col %<o% c("GO", "GO-ID")
ObjNm <- "Annotate"
if ((scrptType == "withReps")&&(ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) {
  ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]]
} else {
  tst <- "Parsed functional annotations" %in% reloadedBckps$Role
  if (!tst) {
    msg <- "Can you provide functional annotations? (required for GO analysis)"
    tmp <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
    ObjNm %<c% tmp
    if (scrptType == "withReps") {
      AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
      tmp <- AllAnsw[1,]
      tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
      tmp$Value <- list(get(ObjNm))
      m <- match(ObjNm, AllAnsw$Parameter)
      if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
    }
  } else { Annotate %<o% TRUE }
}
if (Annotate) {
  if (!exists("Parsed_annotations")) {
    tmpFls <- gsub("\\.fa((s(ta(\\.fas)?)?)|a?)?$", ".txt", fastasTbl$Full)
    AnnotFls <- vapply(tmpFls, function(x) { #x <- tmpFls[1]
      x2 <- gsub(".+/", paste0(fastaLoc, "/"), x)
      if (!file.exists(x)) {
        if (file.exists(x2)) {
          fs::file_copy(x2, wd)
          x <- x2
        } else { x <- NA }
      }
      return(as.character(x))
    }, "")
    AnnotFls %<o% AnnotFls[which(!is.na(AnnotFls))]
    l <- length(AnnotFls)
    if (!l) {
      moar <- TRUE
      while (moar) {
        if (l == 0) {
          msg <- "No functional annotation file detected. Select one?"
        } else { msg <- "Select more?" }
        moar <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
        if (moar) {
          msg <- "Select annotation file"
          #filt <- matrix(c("Annotations txt file", "*.txt"), ncol = 2)
          #AnnotFls <- c(AnnotFls, normalizePath(choose.files(paste0(fastaLoc, "/*.txt"), msg, TRUE, filt), winslah = "/"))
          AnnotFls <- unique(c(AnnotFls,
                               rstudioapi::selectFile(msg,
                                                      path = paste0(fastaLoc, "/*.txt"),
                                                      filter = "UniProtKB txt annotations file (*.txt)")))
          AnnotFls <- AnnotFls[which(!is.na(AnnotFls))]
          l <- length(AnnotFls)
        }
      }
    }
    if (!length(AnnotFls)) {
      warning("No annotations file(s) provided, skipping annotations!")
      Annotate <- FALSE
    } else {
      source(parSrc, local = FALSE)
      Parsed_annotations_lst <- lapply(AnnotFls, function(x) { #x <- AnnotFls[1]
        # If the annotations is not present locally, make a local copy
        if (!file.exists(basename(x))) { fs::file_copy(x, wd) }
        # Parse it
        #tst <- Format.DB_txt(x, usePar = TRUE, cl = parClust)
        return(Format.DB_txt(x, usePar = TRUE, cl = parClust))
      })
      Parsed_annotations <- dplyr::bind_rows(Parsed_annotations_lst)
    }
  }
}
Annotate <- ((exists("Parsed_annotations"))
             &&("data.frame" %in% class(Parsed_annotations))
             &&(nrow(Parsed_annotations) > 0))
if (Annotate) {
  Parsed_annotations %<o% Parsed_annotations
  #
  # Check GO for name degeneracies
  # (the same term has had different names in different files, and we allow for multiple files)
  # Match Parsed_annotations rows to GO name(s)
  wAnnot <- which(nchar(Parsed_annotations$GO) > 0)
  tst1 <- listMelt(strsplit(Parsed_annotations$GO[wAnnot], ";"),
                   wAnnot,
                   c("name", "row"))[, c("row", "name")]
  # Match GO name/IDs to Parsed_annotations row(s)
  tst2 <- set_colnames(aggregate(tst1$row, list(tst1$name), list), c("name", "rows"))
  tst2$ID <- gsub(".*\\[|\\]$", "", tst2$name)
  tst1$ID <- tst2$ID[match(tst1$name, tst2$name)]
  # Match GO IDs to their name(s)
  tst3 <- set_colnames(aggregate(tst2$name, list(tst2$ID), unique), c("ID", "names"))
  tst3$L <- vapply(tst3$names, length, 1)
  wDegen <- which(tst3$L > 1)
  if (length(wDegen)) { # Indicates degeneracy and the need to fix names
    degen <- tst3[which(tst3$L > 1),] # Degenerate ID-to-names mappings
    degen$name <- vapply(degen$names, function(x) { x[[1]] }, "") # Take the first name
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
    #tst3b$L <- vapply(tst3b$names, length, 1)
    #max(tst3b$L) # Should be 1 now
    #
    # Apply new names to tst1
    w2 <- which(tst1$ID %in% tst2$ID)
    tst1$name[w2] <- tst2$name[match(tst1$ID[w2], tst2$ID)]
    #tst3c <- set_colnames(aggregate(tst1$name, list(tst1$ID), unique), c("ID", "names"))
    #tst3c$L <- vapply(tst3c$names, length, 1)
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
              sum(a1 != c1) == 0)
    #
    w <- which(1:nrow(Parsed_annotations) %in% tst1$row)
    Parsed_annotations[w, c("GO", "GO-ID")] <- tst1[match(w, tst1$row), c("GO", "ID")]
    #View(tst1[match(w, tst1$tst1), c("GO", "ID")])
  }
  tst1 <- unlist(strsplit(Parsed_annotations$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(Parsed_annotations$GO, ";"))
  tst3 <- data.table(ID = tst1, Name = tst2)
  tst3 <- tst3[, list(Name = list(unique(Name))), by = list(ID = ID)]
  tst3 <- as.data.frame(tst3)
  tst3$N_names <- vapply(tst3$Name, length, 1)
  #max(tst3$N_names)
  stopifnot(max(tst3$N_names) == 1)
  #View(tst3[which(tst3$N_names > 1),])
  #View(tst3[which(tst3$N_names == 0),])
  #View(Parsed_annotations[, c("GO", "GO-ID")])
  #db <- db[, which(!colnames(db) %in% annot.col)]
  kol <- colnames(Parsed_annotations)
  annot.col %<o% kol[which(!kol %in% c("Accession", "id", "ID", "Names", "Sequence", "MW (Da)"))]
  annot.col2 %<o% annot.col[which(!annot.col %in% colnames(db))]
  annot.col3 %<o% annot.col[which(annot.col %in% colnames(db))]
  if (length(annot.col2)) { db[, annot.col2] <- NA }
  w <- which(db$`Protein ID` %in% Parsed_annotations$Accession)
  mtch <- match(db$`Protein ID`[w], Parsed_annotations$Accession)
  db[, annot.col] <- ""
  db[w, annot.col] <- Parsed_annotations[mtch, annot.col]
  tst1 <- unlist(strsplit(db$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(db$GO, ";"))
  tst3 <- data.table(ID = tst1, Name = tst2)
  tst3 <- tst3[, list(Name = list(unique(Name))), by = list(ID = ID)]
  tst3 <- as.data.frame(tst3)
  tst3$N_names <- vapply(tst3$Name, length, 1)
  #max(tst3$N_names)
  stopifnot(max(tst3$N_names) == 1)
  if (length(annot.col3)) {
    for (kol in annot.col3) { #kol <- annot.col3[1]
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
  #View(tst3[which(vapply(tst3$x, length, 1) > 1),])
  #View(tst3[which(vapply(tst3$x, length, 1) == 0),])
  stopifnot(length(tst3$x) == length(unique(tst3$x)), "character" %in% class(tst3$x))
  db$Ontology <- NULL # Temporary fix for now, this column is broken
  #
  w <- which(vapply(colnames(db), function(x) { "list" %in% class(db[[x]]) }, TRUE))
  if (length(w)) {
    for (i in w) {
      db[[i]] <- parSapply(parClust, db[[i]], paste, collapse = ";")
    }
  }
  saveFun(Parsed_annotations, file = "Parsed_annotations.RData")
  pth <- paste0(wd, "/Parsed, annotated search db.csv")
  f0 <- function() {
    data.table::fwrite(db, pth,
                       quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE, na = "NA")
  }
  tst <- try(f0(), silent = TRUE)
  while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
    dlg_message(paste0("File \"", pth, "\" appears to be locked for editing, close the file then click ok..."), "ok")
    tst <- try(f0(), silent = TRUE)
  }
  #db <- data.table::fread(paste0(wd, "/Parsed, annotated search db.csv"), integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  #write.csv(db, "Parsed, annotated search db.csv", row.names = FALSE)
  #db <- read.csv("Parsed, annotated search db.csv", check.names = FALSE, colClasses = "character")
}
