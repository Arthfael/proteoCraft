### Load and process annotations
# This includes a QC step in case the database differs slightly from the one used by MQ, or if somehow some IDs have not been properly parsed.
GO.col %<o% c("GO", "GO-ID")
ObjNm <- "Annotate"
if ((scrptType == "withReps")&&(ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) {
  ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]]
} else {
  tmp <- "Parsed functional annotations" %in% reloadedBckps$Role
  if (!tmp) {
    msg <- "Can you provide functional annotations? (required for GO analysis)"
    tmp <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  }
  ObjNm %<c% tmp
  if (scrptType == "withReps") {
    AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
    tmp <- AllAnsw[1,]
    tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
    tmp$Value <- list(get(ObjNm))
    m <- match(ObjNm, AllAnsw$Parameter)
    if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
  }
}
if (Annotate) {
  if (!exists("Parsed_annotations")) {
    tmpFls <- gsub("\\.fa((s(ta(\\.fas)?)?)|a?)?$", ".txt", fastasTbl$Full)
    AnnotFls <- vapply(tmpFls, function(x) { #x <- tmpFls[1]
      x2 <- gsub(".+/", "D:/Fasta_databases/", x)
      if (!file.exists(x)) {
        if (file.exists(x2)) {
          fs::file_copy(x2, wd)
          x <- x2
        } else { x <- NA }
      }
      return(as.character(x))
    }, "")
    AnnotFls %<o% AnnotFls[which(!is.na(AnnotFls))]
    if (!length(AnnotFls)) {
      moar <- TRUE
      kount <- 0
      while (moar) {
        if (kount == 0) {
          msg <- "No functional annotation file detected. Select one (or more)?"
        } else { msg <- "Select more?" }
        moar <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
        if (moar) {
          msg <- "Choose annotation file(s):"
          filt <- matrix(c("Annotations txt file", "*.txt"), ncol = 2)
          AnnotFls <- c(AnnotFls, choose.files("D:/Fasta_databases/*.txt", msg, TRUE, filt))
          AnnotFls <- normalizePath(AnnotFls, winslash = "/")
          AnnotFls <- AnnotFls[which(!is.na(AnnotFls))]
          kount <- kount + 1
        }
      }
    }
    if (!length(AnnotFls)) {
      warning("No annotations file(s) provided, skipping annotations!")
      Annotate <- FALSE
    } else {
      source(parSrc, local = FALSE)
      Parsed_annotations <- lapply(AnnotFls, function(x) { #x <- AnnotFls[1]
        # If the annotations is not present locally, make a local copy
        if (!file.exists(basename(x))) { fs::file_copy(x, wd) }
        # Parse it
        return(Format.DB_txt(x, usePar = TRUE, cl = parClust))
      })
      Parsed_annotations <- dplyr::bind_rows(Parsed_annotations)
      tst1 <- unlist(strsplit(Parsed_annotations$`GO-ID`, ";"))
      tst2 <- unlist(strsplit(Parsed_annotations$GO, ";"))
      tst3 <- data.table(A1 = tst1, A2 = tst2)
      tst3 <- tst3[, list(x = unique(A2)), by = list(Group.1 = A1)]
      tst3 <- as.data.frame(tst3)
      stopifnot(length(tst3$x) == length(unique(tst3$x)), "character" %in% class(tst3$x))
      #View(tst3[which(vapply(tst3$x, length, 1) > 1),])
      #View(tst3[which(vapply(tst3$x, length, 1) == 0),])
      #View(Parsed_annotations[, c("GO", "GO-ID")])
    }
  }
  # Check GO for name degeneracies
  # (the same term has had different names in different files, and we allow for multiple files)
  w <- which(nchar(Parsed_annotations$GO) > 0)
  # rows-to-names
  tst1 <- listMelt(strsplit(Parsed_annotations$GO[w], ";"), w, c("name", "row"))
  # names-to-rows-to-IDs
  tst2 <- set_colnames(aggregate(tst1$row, list(tst1$name), list), c("name", "rows"))
  tst2$ID <- gsub(".*\\[|\\]$", "", tst2$name)
  tst1$ID <- tst2$ID[match(tst1$name, tst2$name)]
  tst3 <- set_colnames(aggregate(tst2$name, list(tst2$ID), unique), c("ID", "names"))
  tst3$L <- vapply(tst3$names, length, 1)
  if (max(tst3$L) > 1) { # this would indicate degeneracy and the need to fix names
    tst3 <- tst3[which(tst3$L > 1),]
    tst3$name <- vapply(tst3$names, function(x) { x[[1]] }, "")
    rws <- unique(unlist(tst2$rows[which(tst2$ID %in% tst3$ID)]))
    tst1 <- tst1[which(tst1$row %in% rws),]
    tst2 <- tst2[which(tst2$ID %in% tst3$ID),]
    tst2$name <- tst3$name[match(tst2$ID, tst3$ID)]
    w2 <- which(tst1$ID %in% tst2$ID)
    tst1$name[w2] <- tst2$name[match(tst1$ID[w2], tst2$ID)]
    tst1 <- as.data.table(tst1)
    tst1 <- tst1[, list(GO = paste(name, collapse = ";"),
                        ID = paste(ID, collapse = ";")),
                 by = list(tst1$row)]
    tst1 <- as.data.frame(tst1)
    Parsed_annotations[w, c("GO", "GO-ID")] <- tst1[match(w, tst1$tst1), c("GO", "ID")]
  }
}
Annotate <- exists("Parsed_annotations") # Update condition
if (Annotate) {
  Parsed_annotations %<o% Parsed_annotations
  saveFun(Parsed_annotations, file = "Parsed_annotations.RData")
  #db <- db[, which(!colnames(db) %in% annot.col)]
  kol <- colnames(Parsed_annotations)
  annot.col %<o% kol[which(!kol %in% c("Accession", "id", "ID", "Names", "Sequence", "MW (Da)"))]
  annot.col2 %<o% annot.col[which(!annot.col %in% colnames(db))]
  annot.col3 %<o% annot.col[which(annot.col %in% colnames(db))]
  if (length(annot.col2)) { db[, annot.col2] <- NA }
  mtch <- match(db$`Protein ID`, Parsed_annotations$Accession)
  db[, annot.col2] <- Parsed_annotations[mtch, annot.col2]
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
  data.table::fwrite(db, paste0(wd, "/Parsed, annotated search db.csv"),
                     quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE, na = "NA")
  #db <- data.table::fread(paste0(wd, "/Parsed, annotated search db.csv"), integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  #write.csv(db, "Parsed, annotated search db.csv", row.names = FALSE)
  #db <- read.csv("Parsed, annotated search db.csv", check.names = FALSE, colClasses = "character")
}
