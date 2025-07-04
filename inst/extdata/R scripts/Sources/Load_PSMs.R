# Load PSMs
QuantUMS %<o% FALSE
LabelType %<o% "LFQ" # Default
if (SearchSoft == "MAXQUANT") {
  ObjNm <- "mqparFile"
  if ((scrptType == "withReps")&&(ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) {
    ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]]
  } else {
    filt <- matrix(c("MaxQuant parameters xml file", "*.xml"), ncol = 2)
    msg <- "Select MaxQuant mqpar.xml file"
    tmp <- choose.files(paste0(gsub("(/combined)?/txt.*", "", indir), "/*.xml"), msg, FALSE, filt, 1)
    tmp <- gsub("\\\\", "/", tmp)
    ObjNm %<c% tmp
    indir2 <- dirname(get(ObjNm))
    if (scrptType == "withReps") {
      AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
      tmp <- AllAnsw[1,]
      tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
      tmp$Value <- list(get(ObjNm))
      m <- match(ObjNm, AllAnsw$Parameter)
      if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
    }
    if (scrptType == "noReps") {
      AnalysisParam$"MaxQuant mqpar.xml file" <- get(ObjNm)
    }
  }
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> mqpar.xml file: \", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(paste0(\"     \", mqparFile), prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))")
  mqpar %<o% readLines(mqparFile)
  fastas %<o% gsub(" *</?fastaFilePath> *", "", grep("<fastaFilePath>", mqpar, value = TRUE))
  fastas <- unlist(strsplit(gsub("\\\\", "/", gsub("^Fasta file\t", "", fastas)), ";"))
  MQvers %<o% gsub(" *</?maxQuantVersion> *", "", grep("<maxQuantVersion>", mqpar, value = TRUE))
  MQFold %<o% paste0("C:/MaxQuant/MaxQuant_", c("", "v"), MQvers)
  MQFold <- MQFold[which(dir.exists(MQFold))]
  dr <- gsub("\\\\", "/", gsub(" *</?fixed((Combined)|(Search))Folder> *", "", grep("<fixed((Combined)|(Search))Folder>", mqpar, value = TRUE)))
  if ((!length(dr))||(dr == "")) { dr <- indir }
  PSMsFl %<o% paste0(dr, c("", "/txt", "/combined/txt"), "/evidence.txt")
  PSMsFl <- PSMsFl[which(file.exists(PSMsFl))]
  if (length(PSMsFl) > 1) {
    msg <- "I could not identify MaxQuant's evidence.txt file, please select it manually:"
    PSMsFl <- choose.files(paste0(indir, "/*.txt"), msg)
  }
  tmpMQ <- data.table::fread(PSMsFl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  Cor_MQ <- cor_mod_seq(tmpMQ)
  ev %<o% Cor_MQ$Peptides
  Modifs %<o% Cor_MQ$PTMs
  #
  # Raw files
  g1 <- grep("</?filePaths>", mqpar)
  rawFiles <- gsub("\\\\", "/", gsub(" *</?string> *", "", mqpar[(g1[1]+1):(g1[2]-1)]))
  rawFilesExt <- gsub(".*\\.", "", rawFiles)
  wNtFnd <- which(!file.exists(rawFiles))
  updtFls <- FALSE
  if (length(wNtFnd)) {
    msg <- paste0(c("Some", "All")[(length(wNtFnd) == length(rawFiles))+1], " MS files are missing at the expected location.\n")
    tbl <- data.frame(path = rawFiles[wNtFnd], file = gsub(".*/", "", rawFiles[wNtFnd]), ext = rawFilesExt[wNtFnd])
    locs <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
    dirs <- unique(unlist(lapply(unique(c(indir, dirname(mqparFile))), function(dir) {
      dir <- c(dir, dirname(dir))
      return(c(dir,
               gsub(".*/", paste0(locs$Path[match("Archive folder", locs$Folder)], "/"), dir[1]),
               gsub(".*/", paste0(locs$Path[match("Archive folder", locs$Folder)], "/"), dir[2])))
    })))
    dirs <- dirs[which(dir.exists(dirs))]
    dirs <- dirs[which(dirs != ".")]
    dirFls <- setNames(lapply(dirs, function(dir) { list.files(dir, recursive = TRUE, full.names = TRUE) }), dirs)
    dirDFls <- setNames(lapply(dirs, function(dir) { grep("\\.d$", list.dirs(dir, recursive = TRUE, full.names = TRUE), value = TRUE) }), dirs)
    dirFls2 <- setNames(lapply(names(dirFls), function(nm) { basename(dirFls[[nm]]) }), dirs)
    dirDFls2 <- setNames(lapply(names(dirDFls), function(nm) { gsub(".*/", "", dirDFls[[nm]]) }), dirs)
    for (dir in dirs) {
      tbl[[dir]] <- apply(tbl[, c("file", "ext")], 1, function(x) { #x <- tbl[1, c("file", "ext")]
        rs <- NA
        if (x[[2]] == "d") {
          m <- match(x[[1]], dirDFls2[[dir]])
          if (!is.na(m)) { rs <- dirDFls[[dir]][m] }
        } else {
          m <- match(x[[1]], dirFls2[[dir]])
          if (!is.na(m)) { rs <- dirFls[[dir]][m] }
        }
        return(rs)
      })
    }
    tbl$nuLoc <- apply(tbl[, dirs, drop = FALSE], 1, function(x) {
      x <- x[which(!is.na(x))]
      if (!length(x)) { x <- NA }
      return(x)
    })
    wY <- which(!is.na(tbl$nuLoc))
    wN <- which(is.na(tbl$nuLoc))
    lY <- length(wY)
    lN <- length(wN)
    tstY1 <- (lY > 0)+1 # Did we manage to locate some files automatically?
    tstY2 <- (lY > 1)+1 # ... more than 1?
    tstN1 <- (lN > 0)+1 # Did we fail to for some?
    tstN2 <- (lN > 1)+1 # ... more than 1?
    if (tstY1 == 2) {
      updtFls <- TRUE
      msg <- paste0(msg, "However, the script automatically detected ",
                    c(c("the", "all")[tstY2], "the following")[tstN1], " file", c("", "s")[tstY2],
                    " at the following location", c("", "s")[tstY2], ":\n",
                    paste(paste0(" - ", tbl$nuLoc[wY], "\n"), collapse = ""))
      nuRawFiles <- rawFiles
      nuRawFiles[wNtFnd[wY]] <- tbl$nuLoc[wY]
    }
    cat(msg)
    if (tstN1 == 2) {
      msg2 <- paste0("Select the location of the missing file(s) (or cancel if files are unavailable):")
      newDir <- selectDirectory(msg2, path = wd)
      newFls <- list.files(newDir, full.names = TRUE, include.dirs = TRUE, recursive = TRUE)
      newFlsTst <- gsub(".*/", "", newFls)
      wNA <- which(is.na(tbl$nuLoc))
      tst <- vapply(tbl$file[wNA], function(x) {
        x <- newFls[which(newFlsTst == x)]
        if (length(x) > 1) {
          nc <- nchar(x)
          x <- x[which(nc == min(nc))]
        }
        return(x)
      }, "")
      tst <- tst[which(vapply(tst, length, 1) == 1)]
      tst <- setNames(unlist(tst), names(tst))
      wY2 <- which(is.na(tbl$nuLoc)&(tbl$file %in% names(tst)))
      wN2 <- which(is.na(tbl$nuLoc)&(!tbl$file %in% names(tst)))
      tstY2 <- (length(wY2) > 1)+1
      tstN2 <- (length(wN2) > 1)+1
      if (length(wY2)) {
        m <- match(tbl$file[wY2], names(tst))
        tbl$nuLoc[wY2] <- tst[m]
        msg <- "The user was able to locate the "
        if (length(wN2) == 0) {
          msg <- paste0(msg, "missing file", c("", "s")[tstY2],
                        " in directory ", newDir, "\n\n")
        } else {
          msg <- paste0(msg, "following missing file", c("", "s")[tstY2],
                        " in directory ", newDir, ":\n", paste0(" - ", tbl$file[wY2], collapse = "\n"), "\n\n")
        }
        cat(msg)
      } else {
        msg <- paste0("The following file", c("", "s")[tstN2], " could not be located:\n",
                      paste(paste0(" - ", tbl$path[wN], "\n\n"), collapse = ""))
      }
      nuRawFiles <- rawFiles
      w <- which(!is.na(tbl$nuLoc))
      updtFls <- length(w) > 0
      nuRawFiles[w] <- tbl$nuLoc[w]
    }
  }
  #rawFiles2 %<o% gsub(".*[\\\\/]|\\.((raw)|(mzX?ML)|(d))$", "", rawFiles, ignore.case = TRUE)
  rawFiles2 %<o% unique(ev$"Raw file") # For MaxQuant only!
  tmp <- gsub(".*/|\\.((raw)|(mzX?ML)|(d))$", "", rawFiles, ignore.case = TRUE)
  stopifnot(sum(!rawFiles2 %in% tmp) == 0)
  rawFiles2 <- tmp
  # We also need a full path column!
  if (updtFls) {
    m <- lapply(rawFiles2, function(fl) { which(ev$`Raw file` == fl) })
    m <- listMelt(m, nuRawFiles)
    stopifnot(class(m$value) == "integer",
              sum(!1:nrow(ev) %in% m$value) == 0)
    m <- m[order(m$value),]
    ev$"Raw file path" <- m$L1[match(1:nrow(ev), m$value)]
    rawFiles <- nuRawFiles
  } else {
    ev$"Raw file path" <- rawFiles[match(ev$"Raw file", rawFiles2)]
  }
  #
  g2 <- grep("</?experiments>", mqpar)
  g2 <- gsub(" *</?string> *", "", mqpar[(g2[1]+1):(g2[2]-1)])
  if (length(g2[which(g2 != "")]) == 0) { g2[which(g2 == "")] <- "Exp1" }
  g3 <- grep("</?fractions>", mqpar)
  g3 <- as.integer(gsub("^32767$", "1", gsub(" *</?short> *", "", mqpar[(g3[1]+1):(g3[2]-1)])))
  g4 <- grep("</?ptms>", mqpar)
  g4 <- as.logical(toupper(gsub(" *</?boolean> *", "", mqpar[(g4[1]+1):(g4[2]-1)])))
  g5 <- grep("</?paramGroupIndices>", mqpar)
  g5 <- as.integer(gsub(" *</?int> *", "", mqpar[(g5[1]+1):(g5[2]-1)]))
  FracMap %<o% data.frame("Raw file" = rawFiles,
                          "Raw files name" = rawFiles2,
                          "MQ.Exp" = g2,
                          "Fraction" = g3,
                          "PTMs" = g4,
                          "Parameter group" = g5,
                          "Use" = TRUE,
                          check.names = FALSE)
  # Labeling
  tmp <- gsub(" |_|-", "", toupper(unique(gsub(" *</?lcmsRunType> *", "", grep("<lcmsRunType>", mqpar, value = TRUE)))))
  if (length(tmp) > 1) {
    stop("There appear to be different parameter groups with different labeling methods in this MaxQuant search.\tWe have never tested this case, but it would certainly break the code as currently written.")
  }
  isDIA %<o% grepl("MaxDIA", tmp)
  if (tmp == "STANDARD") {
    tmp2 <- as.integer(gsub(" *</?multiplicity> *", "", grep("<multiplicity>", mqpar, value = TRUE)))
    if (tmp2 != 1) {
      stop("SILAC-labeling (or analogous multiple channel labeling MS1-based workflows) is not currently supported by this workflow!")
    }
  }
  if (tmp == "NEUCODE") { stop("NeuCode SILAC-labeling is not supported by this workflow!") }
  m <- match(tmp,
             c("STANDARD", "REPORTERMS2", "REPORTERMS3", "REPORTERIONMS2", "REPORTERIONMS3", "BOXCAR", "TIMSDDA", "TIMSMAXDIA", "BOXCARMAXDIA"))
  LabelType <- c("LFQ", "Isobaric", "Isobaric", "Isobaric", "Isobaric", "LFQ", "LFQ", "LFQ", "LFQ")[m]
  if ((scrptType == "noReps")&&(LabelType != "LFQ")) { stop("Currently, this workflow only supports LFQ experiments!") }
  if (LabelType == "Isobaric") { # Currently used only with Reps
    FracMap$Isobaric.set <- c(1, "?")[(nrow(FracMap) > 1)+1]
    lab <- gsub(" *</?((inter)|(termi))nalLabel> *", "", grep("<((inter)|(termi))nalLabel>", mqpar, value = TRUE))
    lab <- unique(gsub("-((Nter)|(Lys))", "-", lab))
    IsobarLabDet %<o% gsub(".*-", "", lab)
    IsobarLabPrec %<o% unique(gsub("-.*", "", lab))
    IsobarLab %<o% gsub("[0-9]*plex$", "", IsobarLabPrec)
    kol <- grep("^Reporter intensity [0-9]+$", colnames(ev), value = TRUE)
    stopifnot(length(kol) > 0)
    val <- as.integer(gsub("^Reporter intensity ", "", kol))
    IsobarLab %<c% val
  } else {
    FracMap$"Parent sample" <- FracMap$MQ.Exp
    k <- colnames(FracMap)[which(colnames(FracMap) != "Parent sample")]
    FracMap <- FracMap[, c("Parent sample", k)]
  }
  #
  MinPepSz %<o% as.integer(gsub(" *</?min((PepLen)|(PeptideLength))> *", "", grep(" *</?min((PepLen)|(PeptideLength))> *", mqpar, value = TRUE)))
  Missed %<o% as.integer(gsub(" *</?maxMissedCleavages> *", "", grep(" *</?maxMissedCleavages> *", mqpar, value = TRUE)))
}
if (SearchSoft == "DIANN") {
  isDIA %<o% TRUE
  ObjNm <- "DIANN_log"
  if ((scrptType == "withReps")&&(ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) {
    ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]]
  } else {
    filt <- matrix(c("DiaNN log file", "*.log.txt"), ncol = 2)
    msg <- "Select DiaNN log file"
    ObjNm %<c% choose.files(paste0(indir, "/*.log.txt"), msg, multi = FALSE, filt, 1)
    if (scrptType == "withReps") {
      AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
      tmp <- AllAnsw[1,]
      tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
      tmp$Value <- list(get(ObjNm))
      m <- match(ObjNm, AllAnsw$Parameter)
      if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
    }
    if (scrptType == "noReps") {
      AnalysisParam$"DiaNN .log.txt file" <- get(ObjNm)
    }
  }
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> DiaNN log file: \", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(paste0(\"     \", DIANN_log), prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))")
  DIANNlog %<o% readLines(DIANN_log)
  DIANNCall %<o% grep("^diann.exe ", DIANNlog, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(DIANNCall)) {
    # Sometimes, the diaNN call does not include .exe (when run on Linux), so we can identify the proper row indirectly, based on the fact that afaik it is always immediately after the logical cores
    DIANNCall <- DIANNlog[grep("^Logical CPU cores:", DIANNlog)+1]
    if (!grepl(" --(cfg|f) ", DIANNCall)) {
      # In case the other failed, we try to identify the call by the present of an input files or files list arguments
      DIANNCall <- grep(" --(cfg|f) ", DIANNlog, ignore.case = TRUE, value = TRUE)[1]
      if (is.na(DIANNCall)) { stop("Unhandled exception: could not detect diaNN call in this log!") }
    }
  }
  Args <- unlist(strsplit(DIANNCall, " +--"))
  #
  PSMsFl %<o% gsub("^\"|\"$", "", gsub("\\\\", "/", gsub("^out +", "", grep("^out ", Args, value = TRUE))))
  if (grepl("\\.parquet$", PSMsFl)) {
    PSMsFls <- c(gsub("\\.parquet$", ".tsv", PSMsFl),
                 PSMsFl)
  }
  if (grepl("\\.tsv$", PSMsFl)) {
    PSMsFls <- c(PSMsFl,
                 gsub("\\.tsv$", ".parquet", PSMsFl))
  }
  w <- which(file.exists(PSMsFls))
  if (!length(w)) {
    PSMsFls <- paste0(dirname(DIANN_log), gsub(".*/", "/", PSMsFls))
    w <- which(file.exists(PSMsFls))
    if (length(w)) {
      message("The DiaNN output folder has been renamed or moved since the search was run, but could be located automatically.")
    } else {
      warning("The DiaNN output folder has been renamed or moved since DiaNN was run.\nThe psms report could not be located automatically.\nPrompting user...")
      PSMsFls <- rstudioapi::selectFile("Select DiaNN report file (.tsv or .parquet)", path = indir, filter = "DiaNN report file s (*.tsv|*.parquet)")
      w <- which(file.exists(PSMsFls))
    }
  }
  PSMsFls <- PSMsFls[w]
  l <- length(PSMsFls)
  stopifnot(l > 0)
  #
  fastas %<o% gsub("^\"|\"$", "", gsub("\\\\", "/", gsub("^fasta ", "", grep("^fasta ", Args, value = TRUE))))
  if (!length(fastas)) { # In the case of DiaNN, sometimes there is no fasta, only a library!
    msg <- "Select fasta file(s) used to create the library, or a compatible fasta file"
    dflt <- paste0(wd, "/*.fasta")
    filt <- matrix(data = c("fasta", "*.fasta;*.fas;*.fa;*.fasta.fas"), ncol = 2,
                   dimnames = list("Fasta"))
    fastas <- normalizePath(choose.files(dflt, filters = filt), winslash = "/")
  }
  w <- which(!file.exists(fastas))
  if (length(w)) {
    for (x in w) {
      if (gsub("/[^/]+$", "", fastas[x]) == dirname(DIANN_log)) {
        tmp <- paste0(dirname(DIANN_log), gsub(".*/", "/", fastas[x]))
        if (file.exists(tmp)) {
          warning(paste0("Fasta file ", fastas[x], " has been relocated since DiaNN was run, but could be located automatically."))
          fastas[x] <- tmp
        } else {
          warning(paste0("Fasta file ", fastas[x], " has been relocated since DiaNN was run, and could not be located automatically. Prompting user..."))
          fastas[x] <- normalizePath(choose.files(gsub("/[^/]+$", "/*.tsv", fastas[x]), "Select Fasta report file", FALSE), winslash = "/")
        }
      }
    }
  }
  tmp <- unlist(strsplit(DIANNCall, " +--"))
  rawFiles %<o% unique(gsub("^\"|\"$", "", gsub("\\\\", "/", gsub("^f +", "", grep("^f +", tmp, value = TRUE))))) # Raw files existence is not critical for the script!
  if (!length(rawFiles)) {
    # This is the case where a list of files was provided as argument - usually by FP
    tmp <- gsub("\\\\", "/", gsub("^cfg +", "",  grep("^cfg +", Args, value = TRUE)))
    tmp <- gsub("-- $", "", tmp) # Don't ask me why this is required...
    if (!file.exists(tmp)) { tmp <- paste0(dirname(DIANN_log), gsub(".*/", "/", tmp)) }
    if (!file.exists(tmp)) { tmp <- paste0(gsub("/[^/]+$", "", dirname(DIANN_log)), gsub(".*/", "/", tmp)) } # The DiaNN-run-through-FragPipe case
    if (file.exists(tmp)) {
      rawFiles <- gsub("^ *--f +", "", grep("^ *--f +", readLines(tmp), value = TRUE))
    }
  }
  rawFiles <- grep("\\.exe$", rawFiles, value = TRUE, invert = TRUE) # Yes, that happened!
  rawFiles <- gsub("^ *[\"'] *| *[\"'] *$", "", rawFiles)
  rawFilesExt <- gsub(".*\\.", "", rawFiles)
  wNtFnd <- which(!file.exists(rawFiles))
  updtFls <- FALSE
  if (length(wNtFnd)) {
    msg <- paste0(c("Some", "All")[(length(wNtFnd) == length(rawFiles))+1], " MS files are missing at the expected location.\n")
    tbl <- data.frame(path = rawFiles[wNtFnd], file = gsub(".*/", "", rawFiles[wNtFnd]), ext = rawFilesExt[wNtFnd])
    locs <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
    dirs <- unique(unlist(lapply(unique(c(indir, dirname(DIANN_log))), function(dir) {
      dir <- c(dir, dirname(dir))
      return(c(dir,
               gsub(".*/", paste0(locs$Path[match("Archive folder", locs$Folder)], "/"), dir[1]),
               gsub(".*/", paste0(locs$Path[match("Archive folder", locs$Folder)], "/"), dir[2])))
    })))
    dirs <- dirs[which(dir.exists(dirs))]
    dirs <- dirs[which(dirs != ".")]
    dirFls <- setNames(lapply(dirs, function(dir) { list.files(dir, recursive = TRUE, full.names = TRUE) }), dirs)
    dirDFls <- setNames(lapply(dirs, function(dir) { grep("\\.d$", list.dirs(dir, recursive = TRUE, full.names = TRUE), value = TRUE) }), dirs)
    dirFls2 <- setNames(lapply(names(dirFls), function(nm) { basename(dirFls[[nm]]) }), dirs)
    dirDFls2 <- setNames(lapply(names(dirDFls), function(nm) { gsub(".*/", "", dirDFls[[nm]]) }), dirs)
    for (dir in dirs) {
      tbl[[dir]] <- apply(tbl[, c("file", "ext")], 1, function(x) { #x <- tbl[1, c("file", "ext")]
        rs <- NA
        if (x[[2]] == "d") {
          m <- match(x[[1]], dirDFls2[[dir]])
          if (!is.na(m)) { rs <- dirDFls[[dir]][m] }
        } else {
          m <- match(x[[1]], dirFls2[[dir]])
          if (!is.na(m)) { rs <- dirFls[[dir]][m] }
        }
        return(rs)
      })
    }
    tbl$nuLoc <- apply(tbl[, dirs, drop = FALSE], 1, function(x) {
      x <- x[which(!is.na(x))]
      if (!length(x)) { x <- NA }
      return(x)
    })
    wY <- which(!is.na(tbl$nuLoc))
    wN <- which(is.na(tbl$nuLoc))
    lY <- length(wY)
    lN <- length(wN)
    tstY1 <- (lY > 0)+1 # Did we manage to locate some files automatically?
    tstY2 <- (lY > 1)+1 # ... more than 1?
    tstN1 <- (lN > 0)+1 # Did we fail to locate some?
    tstN2 <- (lN > 1)+1 # ... more than 1?
    if (tstY1 == 2) {
      updtFls <- TRUE
      msg <- paste0(msg, "However, the script automatically detected ",
                    c(c("the", "all")[tstY2], "the following")[tstN1], " file", c("", "s")[tstY2],
                    " at the following location", c("", "s")[tstY2], ":\n",
                    paste(paste0(" - ", tbl$nuLoc[wY], "\n"), collapse = ""))
      nuRawFiles <- rawFiles
      nuRawFiles[wNtFnd[wY]] <- tbl$nuLoc[wY]
    }
    cat(msg)
    if (tstN1 == 2) {
      msg2 <- paste0("Select the location of the missing file(s) (or cancel if files are unavailable):")
      newDir <- selectDirectory(msg2, path = wd)
      newFls <- list.files(newDir, full.names = TRUE, include.dirs = TRUE, recursive = TRUE)
      newFlsTst <- gsub(".*/", "", newFls)
      wNA <- which(is.na(tbl$nuLoc))
      tst <- vapply(tbl$file[wNA], function(x) {
        x <- newFls[which((newFlsTst == x)|(gsub(" ", "", newFlsTst) == x))]
        if (length(x) > 1) {
          nc <- nchar(x)
          x <- x[which(nc == min(nc))]
        }
        return(x)
      }, "")
      tst <- tst[which(vapply(tst, length, 1) == 1)]
      if (length(tst)) {
        tst <- setNames(unlist(tst), names(tst))
        wY2 <- which(is.na(tbl$nuLoc)&(tbl$file %in% names(tst)))
        wN2 <- which(is.na(tbl$nuLoc)&(!tbl$file %in% names(tst)))
        tstY2 <- (length(wY2) > 1)+1
        tstN2 <- (length(wN2) > 1)+1
        if (length(wY2)) {
          m <- match(tbl$file[wY2], names(tst))
          tbl$nuLoc[wY2] <- tst[m]
          msg <- "The user was able to locate the "
          if (length(wN2) == 0) {
            msg <- paste0(msg, "missing file", c("", "s")[tstY2],
                          " in directory ", newDir, "\n\n")
          } else {
            msg <- paste0(msg, "following missing file", c("", "s")[tstY2],
                          " in directory ", newDir, ":\n", paste0(" - ", tbl$file[wY2], collapse = "\n"), "\n\n")
          }
          cat(msg)
        } else {
          msg <- paste0("The following file", c("", "s")[tstN2], " could not be located:\n",
                        paste(paste0(" - ", tbl$path[wN], "\n\n"), collapse = ""))
        }
        nuRawFiles <- rawFiles
        w <- which(!is.na(tbl$nuLoc))
        updtFls <- length(w) > 0
        nuRawFiles[w] <- tbl$nuLoc[w]
      }
    }
  }
  if ((ReLoadPSMsBckp)&&(file.exists(PSMsBckp))) {
    loadFun(PSMsBckp)
  }
  if (!exists("ev_DIANN2MQ")) {
    source(parSrc, local = FALSE)
    ev_DIANN2MQ <- DIANN_to_MQ(PSMsFls,
                               cl = parClust)
    saveFun(ev_DIANN2MQ, file = PSMsBckp)
  }
  ev %<o% ev_DIANN2MQ$Evidence
  if ("QuantUMS" %in% names(ev_DIANN2MQ)) {
    QuantUMS <- ev_DIANN2MQ$QuantUMS
  }
  ev$`Raw file path` <- gsub_Rep("\\\\", "/", ev$`Raw file path`)
  if (updtFls) {
    w <- which(rawFiles != nuRawFiles)
    m <- lapply(rawFiles[w], function(fl) { which(ev$`Raw file path` == fl) })
    lTst <- vapply(m, length, 1)
    wL <- which(lTst > 0)
    if (length(wL)) {
      m <- listMelt(m[wL], nuRawFiles[w[wL]])
      stopifnot(length(which(is.na(m$L1))) == 0)
      w2 <- which(ev$`Raw file path` %in% rawFiles[w])
      m2 <- match(w2, m$value)
      stopifnot(length(which(is.na(m2))) == 0)
      ev$`Raw file path`[w2] <- m$L1[match(w2, m$value)]
    }
  }
  Modifs %<o% ev_DIANN2MQ$PTMs
  ev$"Raw file path" <- gsub_Rep("\\\\", "/", ev$"Raw file path")
  rawFiles %<o% unique(ev$"Raw file path")
  rawFiles2 %<o% gsub(".*/|\\.((raw)|(mzX?ML)|(d))$", "", rawFiles, ignore.case = TRUE)
  FracMap %<o% data.frame("Raw file" = rawFiles,
                          "Raw files name" = rawFiles2,
                          "Parent sample" = rawFiles2,
                          "Fraction" = 1,
                          "Use" = TRUE,
                          check.names = FALSE)
  temp <- unlist(strsplit(DIANNCall, " --min-pep-len "))
  if (length(temp) == 2) { MinPepSz %<o% as.integer(gsub(" --.+", "", unlist(strsplit(DIANNCall, " --min-pep-len "))[2])) } else {
    MinPepSz %<o% min(nchar(ev$Sequence))
  }
  temp <- unlist(strsplit(DIANNCall, " --missed-cleavages "))
  if (length(temp) == 2) { Missed %<o% as.integer(gsub(" --.+", "", unlist(strsplit(DIANNCall, " --missed-cleavages "))[2])) } else {
    Missed %<o% (max(vapply(strsplit(ev$Sequence, "R|K"), length, 1))-1)
  }
}
if (SearchSoft == "FRAGPIPE") {
  ObjNm <- "FP_WorkflowFl"
  if ((scrptType == "withReps")&&(ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) {
    ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]]
  } else {
    filt <- matrix(c("FragPipe workflow", "*.workflow"), ncol = 2)
    msg <- "Select FragPipe .workflow file"
    ObjNm %<c% normalizePath(choose.files(paste0(indir, "/*.workflow"), msg, multi = FALSE, filt, 1), winslash = "/")
    if (scrptType == "withReps") {
      AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
      tmp <- AllAnsw[1,]
      tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
      tmp$Value <- list(get(ObjNm))
      m <- match(ObjNm, AllAnsw$Parameter)
      if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
    }
    if (scrptType == "noReps") {
      AnalysisParam$"FragPipe .workflow file" <- get(ObjNm)
    }
  }
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> FragPipe workflow file: \", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(paste0(\"     \", FP_WorkflowFl), prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))")
  ObjNm <- "FP_ManifestFl"
  if ((scrptType == "withReps")&&(ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) {
    ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]]
  } else {
    filt <- matrix(c("FragPipe samples manifest", "*.fp-manifest"), ncol = 2)
    msg <- "Select FragPipe .fp-manifest file"
    ObjNm %<c% normalizePath(choose.files(paste0(indir, "/*.fp-manifest"), msg, multi = FALSE, filt, 1), winslash = "/")
    if (scrptType == "withReps") {
      AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
      tmp <- AllAnsw[1,]
      tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
      tmp$Value <- list(get(ObjNm))
      m <- match(ObjNm, AllAnsw$Parameter)
      if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
    }
    if (scrptType == "noReps") {
      AnalysisParam$"FragPipe .fp-manifest file" <- get(ObjNm)
    }
  }
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> FragPipe manifest file: \", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(paste0(\"     \", FP_ManifestFl), prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))")
  if ((ReLoadPSMsBckp)&&(file.exists(PSMsBckp))) {
    loadFun(PSMsBckp)
  }
  if (!exists("ev_FP2MQ")) {
    source(parSrc, local = FALSE)
    ev_FP2MQ <- FP_to_MQ(FP_WorkflowFl,
                         FP_ManifestFl,
                         FailIfNoQuant = TRUE,
                         cl = parClust)
    saveFun(ev_FP2MQ, file = PSMsBckp)
  }
  ev %<o% ev_FP2MQ$Evidence
  Modifs %<o% ev_FP2MQ$PTMs
  FP_Manifest %<o% ev_FP2MQ$FracMap
  FP_Workflow %<o% ev_FP2MQ$WorkFlow
  if ("TMT_annotations" %in% names(ev_FP2MQ)) {
    FP_TMT_annot %<o% ev_FP2MQ$TMT_annotations
    LabelType <- "Isobaric"
    IsobarLab %<o% "TMT"
    IsobarLabPrec %<o% paste0(IsobarLab, unique(FP_TMT_annot$plex), "plex")
    IsobarLab %<c% unique(FP_TMT_annot$channel_code)
    IsobarLabDet %<o% unique(FP_TMT_annot$channel)
  }
  rawFiles %<o% FP_Manifest$Path
  rawFilesExt <- gsub(".*\\.", "", rawFiles)
  wNtFnd <- which(!file.exists(rawFiles))
  updtFls <- FALSE
  if (length(wNtFnd)) {
    msg <- paste0(c("Some", "All")[(length(wNtFnd) == length(rawFiles))+1], " MS files are missing at the expected location.\n")
    tbl <- data.frame(path = rawFiles[wNtFnd], file = gsub(".*/", "", rawFiles[wNtFnd]), ext = rawFilesExt[wNtFnd])
    locs <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
    dirs <- unique(unlist(lapply(unique(c(indir, FP_Manifest$Path)), function(dir) {
      dir <- c(dir, dirname(dir))
      return(c(dir,
               gsub(".*/", paste0(locs$Path[match("Archive folder", locs$Folder)], "/"), dir[1]),
               gsub(".*/", paste0(locs$Path[match("Archive folder", locs$Folder)], "/"), dir[2])))
    })))
    dirs <- dirs[which(dir.exists(dirs))]
    dirs <- dirs[which(dirs != ".")]
    dirFls <- setNames(lapply(dirs, function(dir) { list.files(dir, recursive = TRUE, full.names = TRUE) }), dirs)
    dirDFls <- setNames(lapply(dirs, function(dir) { grep("\\.d$", list.dirs(dir, recursive = TRUE, full.names = TRUE), value = TRUE) }), dirs)
    dirFls2 <- setNames(lapply(names(dirFls), function(nm) { basename(dirFls[[nm]]) }), dirs)
    dirDFls2 <- setNames(lapply(names(dirDFls), function(nm) { gsub(".*/", "", dirDFls[[nm]]) }), dirs)
    for (dir in dirs) {
      tbl[[dir]] <- apply(tbl[, c("file", "ext")], 1, function(x) { #x <- tbl[1, c("file", "ext")]
        rs <- NA
        if (x[[2]] == "d") {
          m <- match(x[[1]], dirDFls2[[dir]])
          if (!is.na(m)) { rs <- dirDFls[[dir]][m] }
        } else {
          m <- match(x[[1]], dirFls2[[dir]])
          if (!is.na(m)) { rs <- dirFls[[dir]][m] }
        }
        return(rs)
      })
    }
    tbl$nuLoc <- apply(tbl[, dirs, drop = FALSE], 1, function(x) {
      x <- x[which(!is.na(x))]
      if (!length(x)) { x <- NA }
      return(x)
    })
    wY <- which(!is.na(tbl$nuLoc))
    wN <- which(is.na(tbl$nuLoc))
    lY <- length(wY)
    lN <- length(wN)
    tstY1 <- (lY > 0)+1 # Did we manage to locate some files automatically?
    tstY2 <- (lY > 1)+1 # ... more than 1?
    tstN1 <- (lN > 0)+1 # Did we fail to for some?
    tstN2 <- (lN > 1)+1 # ... more than 1?
    if (tstY1 == 2) {
      updtFls <- TRUE
      msg <- paste0(msg, "However, the script automatically detected ",
                    c(c("the", "all")[tstY2], "the following")[tstN1], " file", c("", "s")[tstY2],
                    " at the following location", c("", "s")[tstY2], ":\n",
                    paste(paste0(" - ", tbl$nuLoc[wY], "\n"), collapse = ""))
      nuRawFiles <- rawFiles
      nuRawFiles[wNtFnd[wY]] <- tbl$nuLoc[wY]
    }
    cat(msg)
    if (tstN1 == 2) {
      msg2 <- paste0("Select the location of the missing file(s) (or cancel if files are unavailable):")
      newDir <- rstudioapi::selectDirectory(msg2, path = wd)
      if (!is.null(newDir)) {
        newFls <- list.files(newDir, full.names = TRUE, include.dirs = TRUE, recursive = TRUE)
        newFlsTst <- gsub(".*/", "", newFls)
        wNA <- which(is.na(tbl$nuLoc))
        tst <- vapply(tbl$file[wNA], function(x) {
          x <- newFls[which(newFlsTst == x)]
          if (length(x) > 1) {
            nc <- nchar(x)
            x <- x[which(nc == min(nc))]
          }
          return(x)
        }, "")
        tst <- tst[which(vapply(tst, length, 1) == 1)]
        tst <- setNames(unlist(tst), names(tst))
        wY2 <- which(is.na(tbl$nuLoc)&(tbl$file %in% names(tst)))
        wN2 <- which(is.na(tbl$nuLoc)&(!tbl$file %in% names(tst)))
        tstY2 <- (length(wY2) > 1)+1
        tstN2 <- (length(wN2) > 1)+1
        if (length(wY2)) {
          m <- match(tbl$file[wY2], names(tst))
          tbl$nuLoc[wY2] <- tst[m]
          msg <- "The user was able to locate the "
          if (length(wN2) == 0) {
            msg <- paste0(msg, "missing file", c("", "s")[tstY2],
                          " in directory ", newDir, "\n\n")
          } else {
            msg <- paste0(msg, "following missing file", c("", "s")[tstY2],
                          " in directory ", newDir, ":\n", paste0(" - ", tbl$file[wY2], collapse = "\n"), "\n\n")
          }
          cat(msg)
        } else {
          msg <- paste0("The following file", c("", "s")[tstN2], " could not be located:\n",
                        paste(paste0(" - ", tbl$path[wN], "\n\n"), collapse = ""))
        }
        nuRawFiles <- rawFiles
        w <- which(!is.na(tbl$nuLoc))
        updtFls <- length(w) > 0
        nuRawFiles[w] <- tbl$nuLoc[w]
      }
    }
  }
  if (updtFls) {
    w <- which(rawFiles != nuRawFiles)
    m <- lapply(rawFiles[w], function(fl) { which(ev$`Raw file path` == fl) })
    m <- listMelt(m, nuRawFiles[w])
    stopifnot(length(which(is.na(m$L1))) == 0)
    w2 <- which(ev$`Raw file path` %in% rawFiles[w])
    m2 <- match(w2, m$value)
    stopifnot(length(which(is.na(m2))) == 0)
    ev$"Raw file path"[w2] <- m$L1[match(w2, m$value)]
  }
  rawFiles <- unique(ev$`Raw file path`)
  rawFiles2 %<o% gsub(".*[\\\\/]|\\.((raw)|(mzX?ML)|(d))$", "", basename(rawFiles), ignore.case = TRUE)
  fastas %<o% gsub("^database\\.db-path=", "", grep("^database\\.db-path=", FP_Workflow, value = TRUE))
  fastas <- gsub("\\\\", "", gsub("\\\\\\\\", "/", fastas))
  PSMsFl %<o% gsub("/+", "/", paste0(indir, "/", FP_Manifest$Samples, "/psm.tsv"))
  Samples %<o% FP_Manifest$Samples
  FracMap %<o% data.frame("Raw file" = rawFiles,
                          "Raw files name" = rawFiles2,
                          "MQ.Exp" = Samples,
                          "Fraction" = 1,
                          "Use" = TRUE,
                          check.names = FALSE)
  if (LabelType == "Isobaric") {
    FracMap$Isobaric.set <- c(1, "?")[(nrow(FracMap) > 1)+1]
  } else {
    FracMap$"Parent sample" <- FracMap$MQ.Exp
  }
  pat <- topattern("msfragger.digest_min_length=")
  MinPepSz %<o% as.integer(gsub(pat, "", grep(pat, FP_Workflow, value = TRUE)))
  Ump %<o% as.logical(toupper(gsub(topattern("diaumpire.run-diaumpire="), "", grep(topattern("diaumpire.run-diaumpire="), FP_Workflow, value = TRUE))))
  Diane %<o% as.logical(toupper(gsub(topattern("diann.run-dia-nn="), "", grep(topattern("diann.run-dia-nn="), FP_Workflow, value = TRUE))))
  isDIA %<o% (Diane|Ump) # I think theoretically you could use either without the other... although probably for any DIA dataset you will use DiaNN
  pat <- topattern("msfragger.allowed_missed_cleavage_1=")
  Missed %<o% as.integer(gsub(pat, "", grep(pat, FP_Workflow, value = TRUE)))
}
if (SearchSoft == "PROTEOMEDISCOVERER") { stop("This part has not yet been re-written for PD!") }
stopifnot(nrow(ev) > 0)
#
if (exists("FracMap_reloaded")) {
  m <- match(FracMap$`Raw file`, FracMap_reloaded$`Raw file`)
  tst <- sum(is.na(m))
  gs <- FALSE
  if (tst) {
    m <- match(gsub(" ", "", FracMap$`Raw file`), gsub(" ", "", FracMap_reloaded$`Raw file`))
    tst <- sum(is.na(m))
    gs <- TRUE
  }
  if (tst) {
    m <- match(FracMap$`Raw files name`, FracMap_reloaded$`Raw files name`)
    tst <- sum(is.na(m))
    gs <- FALSE
  }
  if (tst) {
    m <- match(gsub(" ", "", FracMap$`Raw files name`), gsub(" ", "", FracMap_reloaded$`Raw files name`))
    tst <- sum(is.na(m))
    gs <- TRUE
  }
  if (!tst) {
    k <- c("Raw file", "Raw files name")
    FracMap[, k] <- FracMap_reloaded[m, k]
    if (gs) {
      mEv <- match(gsub(" ", "", ev$`Raw file path`), gsub(" ", "", FracMap$`Raw file`))
    } else {
      mEv <- match(ev$`Raw file path`, FracMap$`Raw file`)
    }
    ev$`Raw file` <- FracMap$`Raw files name`[mEv]
  } else {
    rm(FracMap_reloaded)
  }
}
if (!file.exists(FracMapPath)) {
  write.csv(FracMap, file = FracMapPath, row.names = FALSE)
}
if (exists("fastas_reloaded")) {
  fastas <- unique(c(fastas, fastas_reloaded))
}
