# Load PSMs
source(parSrc)
#
require(unimod)
UniMod <- unimod::modifications
#
searchOutputs %<o% list()
l_inDirs <- length(inDirs)
for (dir_i in 1:l_inDirs) { #dir_i <- 1 #dir_i <- 2
  cat(paste0("Processing input folder",
             c("", paste0(" #", dir_i))[(l_inDirs > 1)+1],
             ":\n\t", inDirs[dir_i], "\n"))
  # Better to run this with a for loop:
  # This will at times ask the user some questions, or will assign some values to the global environment in various ways, etc...
  #
  if (SearchSoft[dir_i] == "MAXQUANT") {
    mqparFl_i <- c()
    if (grepl("/combined/txt$", inDirs[dir_i])) {
      mqparFl_i <- list.files(gsub("/txt$", "", inDirs[dir_i]), "\\.xml$", full.names = TRUE)
    }
    if (grepl("/combined$", inDirs[dir_i])) {
      mqparFl_i <- list.files(inDirs[dir_i], "\\.xml$", full.names = TRUE)
    }
    if (length(mqparFl_i) == 1) {
      cat(" - MaxQuant mqpar.xml file detected automatically\n")
    } else {
      mqparFl_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select MaxQuant mqpar.xml file"),
                                          path = paste0(gsub("/txt[^/]*$", "", inDirs[dir_i]), "/*.xml"),
                                          filter = "XML file (*.xml)")
    }
    #
    ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> mqpar.xml file: \", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
    ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_fpar(Report, fpar(ftext(paste0(\"     \", ", mqparFl_i,
                                                          "), prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))"))
    mqpar_i <- readLines(mqparFl_i)
    fastas_i <- gsub(" *</?fastaFilePath> *", "", grep("<fastaFilePath>", mqpar_i, value = TRUE))
    fastas_i <- unlist(strsplit(gsub("\\\\", "/", gsub("^Fasta file\t", "", fastas_i)), ";"))
    mqVers_i <- gsub(" *</?maxQuantVersion> *", "", grep("<maxQuantVersion>", mqpar_i, value = TRUE))
    mqFold_i <- paste0("C:/MaxQuant/MaxQuant_", c("", "v"), mqVers_i)
    mqFold_i <- mqFold_i[which(dir.exists(mqFold_i))]
    dr <- gsub("\\\\", "/", gsub(" *</?fixed((Combined)|(Search))Folder> *", "", grep("<fixed((Combined)|(Search))Folder>", mqpar_i, value = TRUE)))
    if ((!length(dr))||(dr == "")) { dr <- inDirs[dir_i] }
    psmFls_i <- paste0(dr, c("", "/txt", "/combined/txt"), "/evidence.txt")
    psmFls_i <- psmFls_i[which(file.exists(psmFls_i))]
    if (length(psmFls_i) > 1) {
      psmFls_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": could not find MaxQuant's evidence.txt file, please select it manually"),
                                         path = paste0(inDirs[dir_i], "/*.txt"),
                                         filter = "txt files (*.txt)")
    }
    tmpMQ <- data.table::fread(psmFls_i, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
    Cor_MQ <- cor_mod_seq(tmpMQ)
    ev_i <- Cor_MQ$Peptides
    ev_i$Search_ID <- inDirs[dir_i]
    mods_i <- Cor_MQ$PTMs
    #
    # Raw files
    g1 <- grep("</?filePaths>", mqpar_i)
    rawFiles_i <- gsub("\\\\", "/", gsub(" *</?string> *", "", mqpar_i[(g1[1]+1):(g1[2]-1)]))
    rawFiles_i_ext <- gsub(".*\\.", "", rawFiles_i)
    wNtFnd <- which(!file.exists(rawFiles_i))
    updtFls <- FALSE
    if (length(wNtFnd)) {
      msg <- paste0(" - ", c("Some", "All")[(length(wNtFnd) == length(rawFiles_i))+1], " MS files are missing at the expected location.\n")
      tbl <- data.frame(path = rawFiles_i[wNtFnd], file = gsub(".*/", "", rawFiles_i[wNtFnd]), ext = rawFiles_i_ext[wNtFnd])
      locs <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
      dirs <- unique(unlist(lapply(unique(c(inDirs[dir_i], dirname(mqparFl_i))), function(dir) {
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
        msg <- paste0(msg, "   However, the script automatically detected ",
                      c(c("the", "all")[tstY2], "the following")[tstN1], " file", c("", "s")[tstY2],
                      " at the following location", c("", "s")[tstY2], ":\n",
                      paste(paste0(" - ", tbl$nuLoc[wY], "\n"), collapse = ""))
        rawFiles_i_nu <- rawFiles_i
        rawFiles_i_nu[wNtFnd[wY]] <- tbl$nuLoc[wY]
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
          msg <- "   The user was able to locate the "
          if (length(wN2) == 0) {
            msg <- paste0(msg, "missing file", c("", "s")[tstY2],
                          " in directory ", newDir, "\n")
          } else {
            msg <- paste0(msg, "following missing file", c("", "s")[tstY2],
                          " in directory ", newDir, ":\n", paste0(" - ", tbl$file[wY2], collapse = "\n"), "\n")
          }
          cat(msg)
        } else {
          msg <- paste0("   The following file", c("", "s")[tstN2], " could not be located:\n",
                        paste(paste0(" - ", tbl$path[wN], "\n\n"), collapse = ""))
        }
        rawFiles_i_nu <- rawFiles_i
        w <- which(!is.na(tbl$nuLoc))
        updtFls <- length(w) > 0
        rawFiles_i_nu[w] <- tbl$nuLoc[w]
      }
    }
    #rawFiles_i_2 <- gsub(".*[\\\\/]|\\.((raw)|(mzX?ML)|(d))$", "", rawFiles_i, ignore.case = TRUE)
    rawFiles_i_2 <- unique(ev_i$"Raw file") # For MaxQuant only!
    tmp <- gsub(".*/|\\.((raw)|(mzX?ML)|(d))$", "", rawFiles_i, ignore.case = TRUE)
    stopifnot(sum(!rawFiles_i_2 %in% tmp) == 0)
    rawFiles_i_2 <- tmp
    # We also need a full path column!
    if (updtFls) {
      m <- lapply(rawFiles_i_2, function(fl) { which(ev_i$`Raw file` == fl) })
      m <- listMelt(m, rawFiles_i_nu)
      stopifnot(class(m$value) == "integer",
                sum(!1:nrow(ev_i) %in% m$value) == 0)
      m <- m[order(m$value),]
      ev_i$"Raw file path" <- m$L1[match(1:nrow(ev_i), m$value)]
      rawFiles_i <- rawFiles_i_nu
    } else {
      ev_i$"Raw file path" <- rawFiles_i[match(ev_i$"Raw file", rawFiles_i_2)]
    }
    #
    g2 <- grep("</?experiments>", mqpar_i)
    g2 <- gsub(" *</?string> *", "", mqpar_i[(g2[1]+1):(g2[2]-1)])
    if (length(g2[which(g2 != "")]) == 0) { g2[which(g2 == "")] <- "Exp1" }
    g3 <- grep("</?fractions>", mqpar_i)
    g3 <- as.integer(gsub("^32767$", "1", gsub(" *</?short> *", "", mqpar_i[(g3[1]+1):(g3[2]-1)])))
    g4 <- grep("</?ptms>", mqpar_i)
    g4 <- as.logical(toupper(gsub(" *</?boolean> *", "", mqpar_i[(g4[1]+1):(g4[2]-1)])))
    g5 <- grep("</?paramGroupIndices>", mqpar_i)
    g5 <- as.integer(gsub(" *</?int> *", "", mqpar_i[(g5[1]+1):(g5[2]-1)]))
    fracMap_i <- data.frame("Raw file" = rawFiles_i,
                            "Raw files name" = rawFiles_i_2,
                            "MQ.Exp" = g2,
                            "Fraction" = g3,
                            "PTMs" = g4,
                            "Parameter group" = g5,
                            "Use" = TRUE,
                            check.names = FALSE)
    # Labeling
    tmp <- gsub(" |_|-", "", toupper(unique(gsub(" *</?lcmsRunType> *", "", grep("<lcmsRunType>", mqpar_i, value = TRUE)))))
    if (length(tmp) > 1) {
      stop("There appear to be different parameter groups with different labeling methods in this MaxQuant search.\tWe have never tested this case, but it would certainly break the code as currently written.")
    }
    isDIA_i <- grepl("MaxDIA", tmp)
    if (tmp == "STANDARD") {
      tmp2 <- as.integer(gsub(" *</?multiplicity> *", "", grep("<multiplicity>", mqpar_i, value = TRUE)))
      if (tmp2 != 1) {
        stop("SILAC-labeling (or analogous multiple channel labeling MS1-based workflows) is not currently supported by this workflow!")
      }
    }
    if (tmp == "NEUCODE") { stop("NeuCode SILAC-labeling is not supported by this workflow!") }
    m <- match(tmp,
               c("STANDARD", "REPORTERMS2", "REPORTERMS3", "REPORTERIONMS2", "REPORTERIONMS3", "BOXCAR", "TIMSDDA", "TIMSMAXDIA", "BOXCARMAXDIA"))
    labelType_i <- c("LFQ", "Isobaric", "Isobaric", "Isobaric", "Isobaric", "LFQ", "LFQ", "LFQ", "LFQ")[m]
    if ((scrptType == "noReps")&&(labelType_i != "LFQ")) { stop("Currently, this workflow only supports LFQ experiments!") }
    if (labelType_i == "Isobaric") { # Currently used only with Reps
      fracMap_i$Isobaric.set <- c(1, "?")[(nrow(fracMap_i) > 1)+1]
      lab <- gsub(" *</?((inter)|(termi))nalLabel> *", "", grep("<((inter)|(termi))nalLabel>", mqpar_i, value = TRUE))
      lab <- unique(gsub("-((Nter)|(Lys))", "-", lab))
      isobarLabDet_i <- gsub(".*-", "", lab)
      isobarLabPrec_i <- unique(gsub("-.*", "", lab))
      isobarLab_i <- gsub("[0-9]*plex$", "", isobarLabPrec_i)
      kol <- grep("^Reporter intensity [0-9]+$", colnames(ev_i), value = TRUE)
      stopifnot(length(kol) > 0)
      val <- as.integer(gsub("^Reporter intensity ", "", kol))
      assign(isobarLab_i, val)
    } else {
      fracMap_i$"Parent sample" <- fracMap_i$MQ.Exp
      k <- colnames(fracMap_i)[which(colnames(fracMap_i) != "Parent sample")]
      fracMap_i <- fracMap_i[, c("Parent sample", k)]
    }
    #
    minPepSz_i <- as.integer(gsub(" *</?min((PepLen)|(PeptideLength))> *", "", grep(" *</?min((PepLen)|(PeptideLength))> *", mqpar_i, value = TRUE)))
    missed_i <- as.integer(gsub(" *</?maxMissedCleavages> *", "", grep(" *</?maxMissedCleavages> *", mqpar_i, value = TRUE)))
    #
    # MatMet text template
    moult <- (length(rawFiles_i_2) > 1)+1
    gFx <- (grep("<fixedModifications>", mqpar_i)+1):(grep("</fixedModifications>", mqpar_i)-1)
    gVar <- (grep("<variableModifications>", mqpar_i)+1):(grep("</variableModifications>", mqpar_i)-1)
    FxMd <- gsub(" *</?string> *", "", mqpar_i[gFx])
    FxMdC <- grep("\\(C\\)", FxMd, value = TRUE)
    FxMd <- FxMd[which(!FxMd %in% FxMdC)]
    VarMd <- gsub(" *</?string> *", "", mqpar_i[gVar])
    searchTxt_i <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in MaxQuant")
    if ((length(mqVers_i))&&(nchar(mqVers_i))) {
      searchTxt_i <- paste0(searchTxt_i, " version ", mqVers_i)
    }
    searchTxt_i <- paste0(searchTxt_i, " against TEMPLATELIBTEXT")
    MBR <- gsub(" *</?matchBetweenRuns> *", "", grep("<matchBetweenRuns>", mqpar_i, value = TRUE))
    SecPep <- gsub(" *</?secondPeptide> *", "", grep("<secondPeptide>", mqpar_i, value = TRUE))
    DPep <- gsub(" *</?dependentPeptides> *", "", grep("<dependentPeptides>", mqpar_i, value = TRUE))
    tstFDRs <- grep("Fdr", mqpar_i, value = TRUE)
    FDRs <- gsub(" *</?[A-Z,a-z]*Fdr[A-Z,a-z]*> *", "", tstFDRs)
    w <- which(!FDRs %in% c("False", "True"))
    FDRs <- setNames(as.numeric(FDRs[w]), gsub(" *</?|>.+", "", tstFDRs[w]))
    FDRs <- FDRs[which(names(FDRs) != "psmFdrCrosslink")]
    if (DPep == "False") { FDRs <- FDRs[which(names(FDRs) != "dependentPeptideFdr")] }
    if (length(FxMdC)) {
      if (length(FxMdC) == 1) {
        txt <- paste0("Fixed cysteine modification was set to ", FxMdC, ".")
      } else {
        txt <- paste0(paste(FxMdC[1:(length(FxMdC)-1)], collapse = ", "), " and ", rev(FxMdC)[1], " were included as fixed cysteine modifications.")
      }
      searchTxt_i <- paste0(searchTxt_i, " ", txt)
    }
    if (length(FxMd)) {
      if (length(FxMd) == 1) {
        txt <- paste0("An additional fixed modifications, ", FxMd, ", was also included.")
      } else {
        txt <- paste0("In addition, ", paste(FxMd[1:(length(FxMd)-1)], collapse = ", "), " and ", rev(FxMd)[1], " were included as fixed modifications.")
      }
      searchTxt_i <- paste0(searchTxt_i, " ", txt)
    }
    if (length(VarMd)) {
      if (length(VarMd) == 1) {
        txt <- paste0(VarMd, " was set as variable modification.")
      } else {
        txt <- paste0("Variable modifications were set to ", paste(VarMd[1:(length(VarMd)-1)], collapse = ", "), " and ", rev(VarMd)[1], ".")
      }
      searchTxt_i <- paste0(searchTxt_i, " ", txt)
    }
    tsts <- setNames(as.logical(toupper(c(MBR, SecPep, DPep))), c("Match-between-Runs", "Second Peptide Search", "Dependent Peptides"))
    y <- which(tsts); ly <- length(y)
    n <- which(!tsts); ln <- length(n)
    if (ly) {
      if (ly == 1) { txt <- paste0(names(tsts)[y], " was set to active") } else {
        txt <- paste0(paste(names(tsts)[y[1:(ly-1)]], collapse = ", "), " and ", names(tsts)[y[ly]], " were set to active")
      }
      if (ln) {
        if (ln == 1) { txt <- paste0(txt, ", while ", names(tsts)[n], " was turned off.") } else {
          txt <- paste0(txt, ", while ", paste(names(tsts)[n[1:(ln-1)]], collapse = ", "), " and ", names(tsts)[n[ln]], " were turned off.")
        }
      } else { txt <- paste0(txt, ".") }
    } else {
      txt <- paste0(paste(names(tsts)[n[1:(ln-1)]], collapse = ", "), " and ", names(tsts)[n[ln]], " were all turned off.")
    }
    searchTxt_i <- paste0(searchTxt_i, " ", txt)
    tst <- unique(FDRs)
    lf <- length(FDRs)
    nms <- gsub("Fdr$", "s", names(FDRs))
    if (length(tst) == 1) { txt <- paste0("All FDRs were set to ", tst*100, "%.") } else {
      txt <- paste0("FDRs were set to: ", paste(paste0(FDRs[1:(lf-1)]*100, "% (", nms[1:(lf-1)], ")"), collapse = ", "),
                    " and ", paste0(FDRs[lf]*100, "% (", nms[lf], ")"), ".")
    }
    searchTxt_i <- paste0(searchTxt_i, " ", txt)
    #
    searchOutputs[[dir_i]] <- list(ev = ev_i,
                                   Software_param = list(Software = SearchSoft[dir_i],
                                                         Param_file = mqparFl_i,
                                                         Param = mqpar_i,
                                                         MinPepSz = minPepSz_i,
                                                         Missed = missed_i,
                                                         Vers = mqVers_i,
                                                         Dir = mqFold_i,
                                                         LabelType = labelType_i),
                                   PSMsFls = psmFls_i,
                                   Mods = mods_i,
                                   Fastas = fastas_i,
                                   isDIA = isDIA_i,
                                   Raw_files = list(Full = rawFiles_i,
                                                    Names = rawFiles_i_2),
                                   FracMap = fracMap_i,
                                   MatMet_txt = searchTxt_i)
    if (labelType_i == "Isobaric") {
      tmp <- list(IsobarLab = isobarLab_i,
                  IsobarLabDet = isobarLabDet_i,
                  IsobarLabPrec = isobarLabPrec_i)
      tmp[[isobarLab_i]] <- get(isobarLab_i)
      searchOutputs[[dir_i]]$Software_param[names(tmp)] <- tmp
    }
  }
  if (SearchSoft[dir_i] == "DIANN") {
    labelType_i <- "LFQ" # Could be SILAC in future, will need to be detected!
    isDIA_i <- TRUE
    diaNN_logFl_i <- list.files(inDirs[dir_i], "\\.log\\.txt", full.names = TRUE)
    if (length(diaNN_logFl_i) == 1) {
      cat(" - DiaNN .log.txt file detected automatically\n")
    } else {
      diaNN_logFl_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select DiaNN .log.txt file"),
                                              path = paste0(inDirs[dir_i], "/*.log.txt"),
                                              filter = "DiaNN .log.txt file (*.log.txt)")
    }
    ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> DiaNN log file: \", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
    ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_fpar(Report, fpar(ftext(paste0(\"     \", ", diaNN_logFl_i,
                                                          "), prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))"))
    diaNN_logFlDir_i <- dirname(diaNN_logFl_i)
    diannLog_i <- readLines(diaNN_logFl_i)
    diannCall_i <- grep("^diann.exe ", diannLog_i, ignore.case = TRUE, value = TRUE)[1]
    if (is.na(diannCall_i)) {
      # Sometimes, the diaNN call does not include .exe (when run on Linux), so we can identify the proper row indirectly, based on the fact that afaik it is always immediately after the logical cores
      diannCall_i <- diannLog_i[grep("^Logical CPU cores:", diannLog_i)+1]
      if (!grepl(" --(cfg|f) ", diannCall_i)) {
        # In case the other failed, we try to identify the call by the present of an input files or files list arguments
        diannCall_i <- grep(" --(cfg|f) ", diannLog_i, ignore.case = TRUE, value = TRUE)[1]
        if (is.na(diannCall_i)) { stop("Unhandled exception: could not detect diaNN call in this log!") }
      }
    }
    Args <- unlist(strsplit(diannCall_i, " +--"))
    #
    psmFls_i <- gsub("^\"|\"$", "", gsub("\\\\", "/", gsub("^out +", "", grep("^out ", Args, value = TRUE))))
    tstParquet <- grepl("\\.parquet$", psmFls_i)
    tstTSV <- grepl("\\.tsv$", psmFls_i)
    if (tstParquet) {
      psmFls_i <- c(gsub("\\.parquet$", ".tsv", psmFls_i),
                    psmFls_i)
    }
    if (tstTSV) {
      psmFls_i <- c(psmFls_i,
                    gsub("\\.tsv$", ".parquet", psmFls_i))
    }
    w <- which(file.exists(psmFls_i))
    if (!length(w)) {
      psmFls_i <- paste0(diaNN_logFlDir_i, gsub(".*/", "/", psmFls_i))
      w <- which(file.exists(psmFls_i))
      if (length(w)) {
        cat(" - The DiaNN output folder has been renamed or moved since the search was run, but could be located automatically.\n")
      } else {
        warning(" - The DiaNN output folder has been renamed or moved since DiaNN was run.\nThe psms report could not be located automatically.\nPrompting user...\n")
        psmFls_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select DiaNN report file (.tsv or .parquet)"),
                                           path = inDirs[dir_i],
                                           filter = "DiaNN report file s (*.tsv|*.parquet)")
        w <- which(file.exists(psmFls_i))
      }
    }
    psmFls_i <- psmFls_i[w]
    l <- length(psmFls_i)
    stopifnot(l > 0)
    #
    fastas_i <- gsub("^\"|\"$", "", gsub("\\\\", "/", gsub("^fasta ", "", grep("^fasta ", Args, value = TRUE))))
    if (!length(fastas_i)) { # In the case of DiaNN, sometimes there is no fasta, only a library!
      msg <- "Select fasta file(s) used to create the library, or a compatible fasta file"
      dflt <- paste0(wd, "/*.fasta")
      filt <- matrix(data = c("fasta", "*.fasta;*.fas;*.fa;*.fasta.fas"), ncol = 2,
                     dimnames = list("Fasta"))
      fastas_i <- normalizePath(choose.files(dflt, filters = filt), winslash = "/")
    }
    w <- which(!file.exists(fastas_i))
    if (length(w)) {
      for (x in w) {
        if (gsub("/[^/]+$", "", fastas_i[x]) == diaNN_logFlDir_i) {
          tmp <- paste0(diaNN_logFlDir_i, gsub(".*/", "/", fastas_i[x]))
          if (file.exists(tmp)) {
            warning(paste0(" - ", "Fasta file ", fastas_i[x], " has been relocated since DiaNN was run, but could be located automatically."))
            fastas_i[x] <- tmp
          }
        }
      }
    }
    tmp <- unlist(strsplit(diannCall_i, " +--"))
    rawFiles_i <- unique(gsub("^\"|\"$", "", gsub("\\\\", "/", gsub("^f +", "", grep("^f +", tmp, value = TRUE))))) # Raw files existence is not critical for the script!
    if (!length(rawFiles_i)) {
      # This is the case where a list of files was provided as argument - usually by FP
      tmp <- gsub("\\\\", "/", gsub("^cfg +", "",  grep("^cfg +", Args, value = TRUE)))
      tmp <- gsub("-- $", "", tmp) # Don't ask me why this is required...
      if (!file.exists(tmp)) { tmp <- paste0(diaNN_logFlDir_i, gsub(".*/", "/", tmp)) }
      if (!file.exists(tmp)) { tmp <- paste0(gsub("/[^/]+$", "", diaNN_logFlDir_i), gsub(".*/", "/", tmp)) } # The DiaNN-run-through-FragPipe case
      if (file.exists(tmp)) {
        rawFiles_i <- gsub("^ *--f +", "", grep("^ *--f +", readLines(tmp), value = TRUE))
      }
    }
    rawFiles_i <- grep("\\.exe$", rawFiles_i, value = TRUE, invert = TRUE) # Yes, that happened!
    rawFiles_i <- gsub("^ *[\"'] *| *[\"'] *$", "", rawFiles_i)
    rawFiles_i_ext <- gsub(".*\\.", "", rawFiles_i)
    wNtFnd <- which(!file.exists(rawFiles_i))
    updtFls <- FALSE
    if (length(wNtFnd)) {
      msg <- paste0(" - ", c("Some", "All")[(length(wNtFnd) == length(rawFiles_i))+1], " MS files are missing at the expected location.\n")
      tbl <- data.frame(path = rawFiles_i[wNtFnd], file = gsub(".*/", "", rawFiles_i[wNtFnd]), ext = rawFiles_i_ext[wNtFnd])
      locs <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
      dirs <- unique(unlist(lapply(unique(c(inDirs[dir_i], diaNN_logFlDir_i)), function(dir) {
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
        msg <- paste0(msg, "   However, the script automatically detected ",
                      c(c("the", "all")[tstY2], "the following")[tstN1], " file", c("", "s")[tstY2],
                      " at the following location", c("", "s")[tstY2], ":\n",
                      paste(paste0(" - ", tbl$nuLoc[wY], "\n"), collapse = ""))
        rawFiles_i_nu <- rawFiles_i
        rawFiles_i_nu[wNtFnd[wY]] <- tbl$nuLoc[wY]
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
            msg <- "   The user was able to locate the "
            if (length(wN2) == 0) {
              msg <- paste0(msg, "missing file", c("", "s")[tstY2],
                            " in directory ", newDir, "\n")
            } else {
              msg <- paste0(msg, "following missing file", c("", "s")[tstY2],
                            " in directory ", newDir, ":\n", paste0(" - ", tbl$file[wY2], collapse = "\n"), "\n")
            }
            cat(msg)
          } else {
            msg <- paste0("   The following file", c("", "s")[tstN2], " could not be located:\n",
                          paste(paste0(" - ", tbl$path[wN], "\n\n"), collapse = ""))
          }
          rawFiles_i_nu <- rawFiles_i
          w <- which(!is.na(tbl$nuLoc))
          updtFls <- length(w) > 0
          rawFiles_i_nu[w] <- tbl$nuLoc[w]
        }
      }
    }
    if (exists("ev_DIANN2MQ")) { rm(ev_DIANN2MQ) }
    psmsBckpFl_i <- paste0("diaNN PSMs converted to MQ-like format_", dir_i, ".RData")
    if ((psmsBckpFl_i %in% reloadedBckps)&&(file.exists(psmsBckpFl_i))) {
      loadFun(psmsBckpFl_i)
    }
    if (!exists("ev_DIANN2MQ")) {
      cat(" - Processing PSMs...\n")
      source(parSrc, local = FALSE)
      ev_DIANN2MQ <- DIANN_to_MQ(psmFls_i,
                                 cl = parClust)
      cat(" - Saving...\n")
      saveFun(ev_DIANN2MQ, file = psmsBckpFl_i)
      cat(" -> Done!\n")
    }
    ev_i <- ev_DIANN2MQ$Evidence
    ev_i$Search_ID <- inDirs[dir_i]
    if ("QuantUMS" %in% names(ev_DIANN2MQ)) {
      quantUMS_i <- ev_DIANN2MQ$QuantUMS
    }
    ev_i$`Raw file path` <- gsub_Rep("\\\\", "/", ev_i$`Raw file path`)
    if (updtFls) {
      w <- which(rawFiles_i != rawFiles_i_nu)
      m <- lapply(rawFiles_i[w], function(fl) { which(ev_i$`Raw file path` == fl) })
      lTst <- vapply(m, length, 1)
      wL <- which(lTst > 0)
      if (length(wL)) {
        m <- listMelt(m[wL], rawFiles_i_nu[w[wL]])
        stopifnot(length(which(is.na(m$L1))) == 0)
        w2 <- which(ev_i$`Raw file path` %in% rawFiles_i[w])
        m2 <- match(w2, m$value)
        stopifnot(sum(is.na(m2)) == 0)
        ev_i$`Raw file path`[w2] <- m$L1[match(w2, m$value)]
      }
    }
    mods_i <- ev_DIANN2MQ$PTMs
    ev_i$"Raw file path" <- gsub_Rep("\\\\", "/", ev_i$"Raw file path")
    rawFiles_i <- unique(ev_i$"Raw file path")
    rawFiles_i_2 <- gsub(".*/|\\.((raw)|(mzX?ML)|(d))$", "", rawFiles_i, ignore.case = TRUE)
    fracMap_i <- data.frame("Raw file" = rawFiles_i,
                            "Raw files name" = rawFiles_i_2,
                            "Parent sample" = rawFiles_i_2,
                            "Fraction" = 1,
                            "Use" = TRUE,
                            check.names = FALSE)
    temp <- unlist(strsplit(diannCall_i, " --min-pep-len "))
    if (length(temp) == 2) { minPepSz_i <- as.integer(gsub(" --.+", "", unlist(strsplit(diannCall_i, " --min-pep-len "))[2])) } else {
      minPepSz_i <- min(nchar(ev_i$Sequence))
    }
    temp <- unlist(strsplit(diannCall_i, " --missed-cleavages "))
    if (length(temp) == 2) { missed_i <- as.integer(gsub(" --.+", "", unlist(strsplit(diannCall_i, " --missed-cleavages "))[2])) } else {
      missed_i <- (max(vapply(strsplit(ev_i$Sequence, "R|K"), length, 1))-1)
    }
    diannVers_i <- gsub("^DIA-NN | .+$", "", grep("^DIA-NN [0-9]+\\.?[0-9]*", diannLog_i, value = TRUE))
    #
    # MatMet text template
    moult <- (length(rawFiles_i_2) > 1)+1
    diannCall_i2 <- unlist(strsplit(diannCall_i, " +--"))
    FxMd <- gsub(" enabled as a fixed modification$", "", grep(" enabled as a fixed modification$", diannLog_i, value = TRUE))
    FxMdC <- grep("Cysteine", FxMd, value = TRUE)
    FxMd <- FxMd[which(!FxMd %in% FxMdC)]
    FxMdC <- gsub("cysteine[ ,-]?", "", FxMdC, ignore.case = TRUE)
    FxMdC <- paste0(toupper(substr(FxMdC, 1, 1)), substr(FxMdC, 2, nchar(FxMdC)))
    VarMd <- gsub(" will be considered as variable$", "", grep(" will be considered as variable$", diannLog_i, value = TRUE))
    if (length(VarMd)) {
      VarMd <- set_colnames(Isapply(strsplit(gsub("^Modification ", "", VarMd), " with mass delta | at "), unlist), c("UniMod", "Delta mass", "AA"))
      VarMd <- set_colnames(aggregate(VarMd$AA, list(VarMd$UniMod, VarMd$`Delta mass`), list), c("UniMod", "Delta mass", "AA"))
      VarMd$UniMod <-  gsub("^UniMod:", "", VarMd$UniMod)
      VarMd$Type <- "Variable"
      VarMd$"Full name" <- apply(VarMd[, c("UniMod", "Delta mass")], 1, function(x) { #x <- VarMd[1, c("UniMod", "Delta mass")]
        unique(UniMod$Name[which((UniMod$UnimodId == x[[1]])&(round((UniMod$MonoMass == x[[2]])*2)/2 == 0))])
      })
      w <- which(lapply(VarMd$`Full name`, length) == 0)
      if (length(w)) {
        VarMd$`Full name`[w] <- apply(VarMd[w, c("Delta mass", "AA")], 1, function(x) {
          pos <- sort(unlist(x[[2]]))
          pos[which(pos == "*n")] <- "_" # Check, this may not be correct...
          # in fact, check whether this part of DIANN_to_MQ() should not be improved.
          Modifs$`Full name`[which((Modifs$`Mass shift` == x[[1]])&(Modifs$AA %in% pos))]
        })
      }
      VarMd$Text <- apply(VarMd[, c("Full name", "AA")], 1, function(x) {
        x[[2]][which(x[[2]] == "n")] <- "N-term"
        x[[2]][which(x[[2]] == "*n")] <- "protein N-term"
        x[[2]][which(x[[2]] == "c")] <- "C-term"
        x[[2]][which(x[[2]] == "*c")] <- "protein C-term"
        paste0(x[[1]], " (", paste(x[[2]], collapse = ""), ")")
      })
    }
    # Note that it is possible that a modification was searched but not found (so is present in VarMd but not Modifs)
    MBR <- grepl("--reanaly[s,z]e", diannCall_i)+1
    Libraries <- gsub("^lib +", "", grep("^lib", diannCall_i2, value = TRUE))
    nL <- length(Libraries)
    LibFree <- (nL == 0)+1
    searchTxt_i <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in DiaNN")
    if ((length(diannVers_i))&&(nchar(diannVers_i))) { searchTxt_i <- paste0(searchTxt_i, " version ", diannVers_i) }
    searchTxt_i <- paste0(searchTxt_i,
                        c(paste0(" against the following spectral librar", c("y", "ies")[(nL > 1)+1], ": ",
                                 paste(Libraries, collapse = "/"), " (see supplementary material), using "),
                          " in library-free mode against ")[LibFree],
                        "TEMPLATELIBTEXT")
    searchTxt_i <- paste0(searchTxt_i, " Match-Between-Runs was turned ", c("off", "on")[MBR], ".")
    if (length(FxMdC)) {
      if (length(FxMdC) == 1) {
        txt <- paste0("Fixed cysteine modification was set to ", FxMdC, ".")
      } else {
        txt <- paste0(paste(FxMdC[1:(length(FxMdC)-1)], collapse = ", "), " and ", rev(FxMdC)[1], " were included as fixed cysteine modifications.")
      }
      searchTxt_i <- paste0(searchTxt_i, " ", txt)
    }
    if (length(FxMd)) {
      if (length(FxMd) == 1) {
        txt <- paste0("An additional fixed modifications, ", FxMd, ", was also included.")
      } else {
        txt <- paste0("In addition, ", paste(FxMd[1:(length(FxMd)-1)], collapse = ", "), " and ", rev(FxMd)[1], " were included as fixed modifications.")
      }
      searchTxt_i <- paste0(searchTxt_i, " ", txt)
    }
    if (("data.frame" %in% class(VarMd))&&(nrow(VarMd))) {
      if (nrow(VarMd) == 1) {
        txt <- paste0(VarMd$Text, " was set as variable modification.")
      } else {
        txt <- paste0("Variable modifications were set to ", paste(VarMd$Text[1:(nrow(VarMd)-1)], collapse = ", "), " and ", rev(VarMd$Text)[1], ".")
      }
      searchTxt_i <- paste0(searchTxt_i, " ", txt)
    }
    tstFDR <- grep("^Output will be filtered at .+ FDR$", diannLog_i, value = TRUE)
    FDR <- gsub("^Output will be filtered at | FDR$", "", tstFDR)
    if (!grepl(" ?%$", FDR)) { FDR <- as.numeric(FDR)*100 }
    if ((FDR >= 0)&&(FDR <= 1)) { searchTxt_i <- paste0(searchTxt_i, " Data was filtered at ", FDR, "% FDR.") }
    #
    searchOutputs[[dir_i]] <- list(ev = ev_i,
                                   Software_param = list(Software = SearchSoft[dir_i],
                                                         Param_file = diaNN_logFl_i,
                                                         Param = diannLog_i,
                                                         Call = diannCall_i,
                                                         MinPepSz = minPepSz_i,
                                                         Missed = missed_i,
                                                         Vers = diannVers_i,
                                                         LabelType = labelType_i,
                                                         QuantUMS = quantUMS_i),
                                   PSMsFls = psmFls_i,
                                   Mods = mods_i,
                                   Fastas = fastas_i,
                                   isDIA = isDIA_i,
                                   Raw_files = list(Full = rawFiles_i,
                                                    Names = rawFiles_i_2),
                                   FracMap = fracMap_i,
                                   MatMet_txt = searchTxt_i)
  }
  if (SearchSoft[dir_i] == "FRAGPIPE") {
    fpWorkflowFl_i <- list.files(inDirs[dir_i], "\\.workflow$", full.names = TRUE)
    if (length(fpWorkflowFl_i) == 1) {
      cat(" - FragPipe workflow file detected automatically\n")
    } else {
      fpWorkflowFl_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select FragPipe workflow file"),
                                               path = paste0(inDirs[dir_i], "/*.workflow"),
                                               filter = "FragPipe workflow file (*.workflow)")
    }
    fpManifestFl_i <- list.files(inDirs[dir_i], "\\.fp-manifest$", full.names = TRUE)
    if (length(fpManifestFl_i) == 1) {
      cat(" - FragPipe manifest file detected automatically\n")
    } else {
      fpManifestFl_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select FragPipe manifest file"),
                                               path = paste0(inDirs[dir_i], "/*.fp-manifest"),
                                               filter = "FragPipe manifest file (*.fp-manifest)")
    }
    ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> FragPipe workflow file: \", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
    ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_fpar(Report, fpar(ftext(paste0(\"     \", ", fpWorkflowFl_i,
                                                          "), prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))"))
    ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> FragPipe manifest file: \", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
    ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_fpar(Report, fpar(ftext(paste0(\"     \", ", fpManifestFl_i,
                                                          "), prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))"))
    if (exists("ev_FP2MQ")) { rm(ev_FP2MQ) }
    psmsBckpFl_i <- paste0("FragPipe PSMs converted to MQ-like format_", dir_i, ".RData")
    if ((psmsBckpFl_i %in% reloadedBckps)&&(file.exists(psmsBckpFl_i))) {
      loadFun(psmsBckpFl_i)
    }
    if (!exists("ev_FP2MQ")) {
      cat(" - Processing PSMs...\n")
      source(parSrc, local = FALSE)
      ev_FP2MQ <- FP_to_MQ(fpWorkflowFl_i,
                           fpManifestFl_i,
                           FailIfNoQuant = TRUE,
                           cl = parClust)
      cat(" - Saving...\n")
      saveFun(ev_FP2MQ, file = psmsBckpFl_i)
      cat(" -> Done!\n")
    }
    ev_i <- ev_FP2MQ$Evidence
    ev_i$Search_ID <- inDirs[dir_i]
    mods_i <- ev_FP2MQ$PTMs
    fpManifest_i <- ev_FP2MQ$FracMap
    fpWorkflow_i <- ev_FP2MQ$WorkFlow
    psmFls_i <- ev_FP2MQ$PTMs
    labelType_i <- "LFQ"
    if ("TMT_annotations" %in% names(ev_FP2MQ)) {
      fpTMT_annot_i <- ev_FP2MQ$TMT_annotations
      labelType_i <- "Isobaric"
      isobarLab_i <- "TMT"
      isobarLabPrec_i <- paste0(isobarLab_i, unique(fpTMT_annot_i$plex), "plex")
      assign(isobarLab_i, unique(fpTMT_annot_i$channel_code))
      isobarLabDet_i <- unique(fpTMT_annot_i$channel)
    }
    rawFiles_i <- fpManifest_i$Path
    rawFiles_i_ext <- gsub(".*\\.", "", rawFiles_i)
    wNtFnd <- which(!file.exists(rawFiles_i))
    updtFls <- FALSE
    if (length(wNtFnd)) {
      msg <- paste0(" - ", c("Some", "All")[(length(wNtFnd) == length(rawFiles_i))+1], " MS files are missing at the expected location.\n")
      tbl <- data.frame(path = rawFiles_i[wNtFnd], file = gsub(".*/", "", rawFiles_i[wNtFnd]), ext = rawFiles_i_ext[wNtFnd])
      locs <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
      dirs <- unique(unlist(lapply(unique(c(inDirs[dir_i], fpManifest_i$Path)), function(dir) {
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
        msg <- paste0(msg, "   However, the script automatically detected ",
                      c(c("the", "all")[tstY2], "the following")[tstN1], " file", c("", "s")[tstY2],
                      " at the following location", c("", "s")[tstY2], ":\n",
                      paste(paste0(" - ", tbl$nuLoc[wY], "\n"), collapse = ""))
        rawFiles_i_nu <- rawFiles_i
        rawFiles_i_nu[wNtFnd[wY]] <- tbl$nuLoc[wY]
      }
      cat(msg)
      if (tstN1 == 2) {
        msg2 <- paste0("Select the location of the missing file(s) (or cancel if files are unavailable):")
        newDir <- rstudioapi::selectDirectory(msg2, path = wd)
        if (!is.null(newDir)) {
          newFls <- list.files(newDir, full.names = TRUE, include.dirs = TRUE, recursive = TRUE)
          newFlsTst <- gsub(".*/", "", newFls)
          wNA <- which(is.na(tbl$nuLoc))
          tst <- setNames(lapply(tbl$file[wNA], function(x) { #x <- tbl$file[wNA[1]]
            x <- newFls[which(newFlsTst == x)]
            if (length(x) > 1) {
              nc <- nchar(x)
              x <- x[which(nc == min(nc))]
            }
            return(x)
          }), tbl$file[wNA])
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
              msg <- "   The user was able to locate the "
              if (length(wN2) == 0) {
                msg <- paste0(msg, "missing file", c("", "s")[tstY2],
                              " in directory ", newDir, "\n")
              } else {
                msg <- paste0(msg, "following missing file", c("", "s")[tstY2],
                              " in directory ", newDir, ":\n", paste0(" - ", tbl$file[wY2], collapse = "\n"), "\n")
              }
              cat(msg)
            } else {
              msg <- paste0("   The following file", c("", "s")[tstN2], " could not be located:\n",
                            paste(paste0(" - ", tbl$path[wN], "\n\n"), collapse = ""))
            }
          }
          rawFiles_i_nu <- rawFiles_i
          w <- which(!is.na(tbl$nuLoc))
          updtFls <- length(w) > 0
          rawFiles_i_nu[w] <- tbl$nuLoc[w]
        }
      }
    }
    if (updtFls) {
      w <- which(rawFiles_i != rawFiles_i_nu)
      m <- lapply(rawFiles_i[w], function(fl) { which(ev_i$`Raw file path` == fl) })
      m <- listMelt(m, rawFiles_i_nu[w])
      stopifnot(length(which(is.na(m$L1))) == 0)
      w2 <- which(ev_i$`Raw file path` %in% rawFiles_i[w])
      m2 <- match(w2, m$value)
      stopifnot(length(which(is.na(m2))) == 0)
      ev_i$"Raw file path"[w2] <- m$L1[match(w2, m$value)]
    }
    rawFiles_i <- unique(ev_i$`Raw file path`)
    rawFiles_i_2 <- gsub(".*[\\\\/]|\\.((raw)|(mzX?ML)|(d))$", "", basename(rawFiles_i), ignore.case = TRUE)
    fastas_i <- gsub("^database\\.db-path=", "", grep("^database\\.db-path=", fpWorkflow_i, value = TRUE))
    fastas_i <- gsub("\\\\", "", gsub("\\\\\\\\", "/", fastas_i))
    psmFls_i <- gsub("/+", "/", paste0(inDirs[dir_i], "/", fpManifest_i$Samples, "/psm.tsv"))
    samples_i <- fpManifest_i$Samples
    fracMap_i <- data.frame("Raw file" = rawFiles_i,
                            "Raw files name" = rawFiles_i_2,
                            "MQ.Exp" = samples_i,
                            "Fraction" = 1,
                            "Use" = TRUE,
                            check.names = FALSE)
    if (labelType_i == "Isobaric") {
      fracMap_i$Isobaric.set <- c(1, "?")[(nrow(fracMap_i) > 1)+1]
    } else {
      fracMap_i$"Parent sample" <- fracMap_i$MQ.Exp
    }
    pat <- topattern("msfragger.digest_min_length=")
    minPepSz_i <- as.integer(gsub(pat, "", grep(pat, fpWorkflow_i, value = TRUE)))
    pat <- topattern("msfragger.allowed_missed_cleavage_1=")
    missed_i <- as.integer(gsub(pat, "", grep(pat, fpWorkflow_i, value = TRUE)))
    fpVers_i <- gsub("^# FragPipe \\(|\\).+$", "", grep("^# FragPipe \\(", fpWorkflow_i, value = TRUE))
    #
    # MatMet text template
    moult <- (length(rawFiles_i_2) > 1)+1
    FxMdTbl <- gsub(topattern("msfragger.table.fix-mods="), "", grep(topattern("msfragger.table.fix-mods="), fpWorkflow_i, value = TRUE))
    FxMdTbl <- unlist(strsplit(FxMdTbl, ";"))
    FxMdTbl <- Isapply(strsplit(FxMdTbl, ","), unlist)
    colnames(FxMdTbl) <- c("Shift", "Site", "Enabled", "MaxOccurences")
    FxMdTbl <- FxMdTbl[which(FxMdTbl$Enabled == "true"),]
    FxMdTbl$AAName <- gsub(".+\\(|\\)", "", FxMdTbl$Site)
    FxMdTbl$AAName <- sapply(FxMdTbl$AAName, function(x) {
      paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    })
    FxMdTbl$AA <- gsub(" \\(.+\\)", "", FxMdTbl$Site)
    FxMdTbl$Shift <- as.numeric(FxMdTbl$Shift)
    FxMdTbl <- FxMdTbl[which(FxMdTbl$Shift != 0),]
    FxMdTbl$Shift <- as.character(FxMdTbl$Shift)
    g <- grep("^-", FxMdTbl$Shift, invert = TRUE)
    FxMdTbl$Shift[g] <- paste0("+", FxMdTbl$Shift[g])
    FxMdTbl$Name <- apply(FxMdTbl[, c("Shift", "AAName")], 1, function(x) { paste0(x[[1]], " (", x[[2]], ")") })
    FxMdC <- FxMdTbl$Name[which(FxMdTbl$AA == "C")]
    FxMd <- FxMdTbl$Name[which(FxMdTbl$AA != "C")]
    VarMd <- gsub(topattern("msfragger.table.var-mods="), "", grep(topattern("msfragger.table.var-mods="), fpWorkflow_i, value = TRUE))
    VarMd <- unlist(strsplit(VarMd, ";"))
    VarMd <- Isapply(strsplit(VarMd, ","), unlist)
    colnames(VarMd) <- c("Shift", "Site", "Enabled", "MaxOccurences")
    VarMd <- VarMd[which(VarMd$Enabled == "true"),]
    for (i in 1:nrow(AA_table)) { VarMd$Site <- gsub(AA_table$AA[i], paste0(";", AA_table$AA[i], ";"), VarMd$Site) }
    VarMd$Site <- gsub("n", ";peptide N-term;", VarMd$Site)
    VarMd$Site <- gsub("c", ";peptide C-term;", VarMd$Site)
    VarMd$Site <- gsub("\\[\\^", ";protein N-term;", VarMd$Site)
    VarMd$Site <- gsub("\\^\\]", ";protein C-term;", VarMd$Site)
    VarMd$Site <- strsplit(gsub("^;|;$", "", VarMd$Site), ";;")
    VarMd$AAName <- sapply(VarMd$Site, function(x) { #
      x <- unlist(x)
      w <- which(nchar(x) == 1)
      x[w] <- as.character(AA_table$Name[match(x[w], as.character(AA_table$AA))])
      return(x)
    })
    VarMd$Shift <- gsub(" +", "", VarMd$Shift)
    g <- grep("^-", VarMd$Shift, invert = TRUE)
    VarMd$Shift[g] <- paste0("+", VarMd$Shift[g])
    VarMd$Text <- apply(VarMd[, c("Shift", "AAName", "Site")], 1, function(x) {
      x2 <- x[[2]]
      if (length(x2) > 1) { x2 <- paste(x[[3]], collapse = "") }
      return(paste0(x[[1]], " (", x2, ")"))
    })
    # Note that it is possible that a modification was searched but not found (so is present in VarMd but not Modifs)
    #
    # FragPipe is very versatile... which will not make this easy...
    # Unfortunately, the workflow file does not appear to store the type of Workflow which was preselected in the 2nd tab.
    Ump <- as.logical(toupper(gsub(topattern("diaumpire.run-diaumpire="), "", grep(topattern("diaumpire.run-diaumpire="), fpWorkflow_i, value = TRUE))))
    Diane <- as.logical(toupper(gsub(topattern("diann.run-dia-nn="), "", grep(topattern("diann.run-dia-nn="), fpWorkflow_i, value = TRUE))))
    isDIA_i <- Diane|Ump # I think theoretically you could use either without the other... although probably for any DIA dataset you will use DiaNN
    #stopifnot(((Param$Label == "DIA")&isDIA)|((Param$Label != "DIA")&!isDIA)) # Sanity check
    # OpenSearch is indirectly identified as cases where unusually high values for tolerances are used
    PrTolUp <- gsub(topattern("msfragger.precursor_mass_upper="), "", grep(topattern("msfragger.precursor_mass_upper="), fpWorkflow_i, value = TRUE))
    PrTolDwn <- gsub(topattern("msfragger.precursor_mass_lower="), "", grep(topattern("msfragger.precursor_mass_lower="), fpWorkflow_i, value = TRUE))
    # Validations
    Chris <- as.logical(toupper(gsub(topattern("crystalc.run-crystalc="), "", grep(topattern("crystalc.run-crystalc="), fpWorkflow_i, value = TRUE))))
    Moses <- as.logical(toupper(gsub(topattern("peptide-prophet.run-peptide-prophet="), "", grep(topattern("peptide-prophet.run-peptide-prophet="), fpWorkflow_i, value = TRUE))))
    Illy <- as.logical(toupper(gsub(topattern("percolator.run-percolator="), "", grep(topattern("percolator.run-percolator="), fpWorkflow_i, value = TRUE))))
    Maud <- as.logical(toupper(gsub(topattern("ptmprophet.run-ptmprophet="), "", grep(topattern("ptmprophet.run-ptmprophet="), fpWorkflow_i, value = TRUE))))
    Camus <- as.logical(toupper(gsub(topattern("phi-report.run-report="), "", grep(topattern("phi-report.run-report="), fpWorkflow_i, value = TRUE))))
    David <- as.logical(toupper(gsub(topattern("ptmshepherd.run-shepherd="), "", grep(topattern("ptmshepherd.run-shepherd="), fpWorkflow_i, value = TRUE))))
    Kant <- as.logical(toupper(gsub(topattern("quantitation.run-label-free-quant="), "", grep(topattern("quantitation.run-label-free-quant="), fpWorkflow_i, value = TRUE))))
    Jon <- as.logical(toupper(gsub(topattern("ionquant.run-ionquant="), "", grep(topattern("ionquant.run-ionquant="), fpWorkflow_i, value = TRUE))))
    Matt <- as.logical(as.numeric(gsub(topattern("ionquant.mbr="), "", grep(topattern("ionquant.mbr="), fpWorkflow_i, value = TRUE))))
    Fritz <- as.logical(toupper(gsub(topattern("freequant.run-freequant="), "", grep(topattern("freequant.run-freequant="), fpWorkflow_i, value = TRUE))))
    Timothy <- as.logical(toupper(gsub(topattern("tmtintegrator.run-tmtintegrator="), "", grep(topattern("tmtintegrator.run-tmtintegrator="), fpWorkflow_i, value = TRUE))))
    # (I know, this is really bad, but I was tired, ok!)
    #
    searchTxt_i <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in FragPipe")
    if ((length(fpVers_i))&&(nchar(fpVers_i))) {
      searchTxt_i <- paste0(searchTxt_i, " version ", fpVers_i)
    }
    searchTxt_i <- paste0(searchTxt_i, " against TEMPLATELIBTEXT")
    if (Ump == 1) {
      searchTxt_i <- gsub("\\.$",
                        " As a preliminary step, the DIAUmpire module was run to extract pseudo-MS/MS spectra from DIA spectra files.",
                        searchTxt_i)
    }
    if (length(FxMdC)) {
      if (length(FxMdC) == 1) {
        txt <- paste0("Fixed cysteine modification was set to ", FxMdC, ".")
      } else {
        txt <- paste0(paste(FxMdC[1:(length(FxMdC)-1)], collapse = ", "), " and ", rev(FxMdC)[1], " were included as fixed cysteine modifications.")
      }
      searchTxt_i <- paste0(searchTxt_i, " ", txt)
    }
    if (length(FxMd)) {
      if (length(FxMd) == 1) {
        txt <- paste0("An additional fixed modifications, ", FxMd, ", was also included.")
      } else {
        txt <- paste0("In addition, ", paste(FxMd[1:(length(FxMd)-1)], collapse = ", "), " and ", rev(FxMd)[1], " were included as fixed modifications.")
      }
      searchTxt_i <- paste0(searchTxt_i, " ", txt)
    }
    if (nrow(VarMd)) {
      if (nrow(VarMd) == 1) {
        txt <- paste0(VarMd$Text, " was set as variable modification.")
      } else {
        txt <- paste0("Variable modifications were set to ", paste(VarMd$Text[1:(nrow(VarMd)-1)], collapse = ", "), " and ", rev(VarMd$Text)[1], ".")
      }
      searchTxt_i <- paste0(searchTxt_i, " ", txt)
    }
    w <- which(c(PrTolUp != "20", PrTolDwn != "-20"))
    l <- length(w)
    if (l) {
      PrTol <- c("lower", "upper")[w]
      PrTol[1] <- paste0(toupper(substr(PrTol[1], 1, 1)), substr(PrTol[1], 2 , nchar(PrTol[w])))
      PrTol <- paste(PrTol, collapse = " and ")
      searchTxt_i <- paste0(searchTxt_i, " ", PrTol, " precursor mass tolerance", c("", "s")[l], " w", c("as", "ere")[l], " set to ",
                          paste(c(PrTolDwn, PrTolUp)[w], collapse = " and "), " ppm", c("", ", respectively")[l], ".")
    }
    # Validations
    w <- which(c(Chris, Moses, Illy)) # Normally should only ever be length = 1
    l <- length(w)
    if (l) {
      txt <- c("Crystal-C", "PeptideProphet", "Percolator")[w]
      if (l > 1) { txt <- paste0(paste(txt[1:(l-1)], collapse = ", "), " and ", txt[l]) }
      searchTxt_i <- paste0(searchTxt_i, " Peptide identifications were validated using ", txt, ".")
      if (Maud) { searchTxt_i <- gsub("\\.$", ", and", searchTxt_i) }
    }
    if (Maud) { searchTxt_i <- paste0(searchTxt_i, " PTM sites", c(" were", "")[(l>0)+1], " localized using PTMProphet.") }
    # FDR levels
    if (Camus) {
      CamusPar <- toupper(gsub(topattern("phi-report.filter="), "", grep(topattern("phi-report.filter="), fpWorkflow_i, value = TRUE)))
      CamusPar <- unlist(strsplit(CamusPar, "--"))
      CamusPar <- CamusPar[2:length(CamusPar)]
      CamusPar <- data.frame(Par = CamusPar)
      CamusPar$Val <- gsub("[^ ]+ ", "", CamusPar$Par)
      CamusPar$Par <- gsub(" .+", "", CamusPar$Par)
      lvls <- setNames(c("ion", "psm", "peptide", "protein"), c("ION", "PSM", "PEP", "PROT"))
      w <- which(CamusPar$Par %in% names(lvls))
      l <- length(w)
      if (length(w)) {
        txt1 <- lvls[CamusPar$Par[w]]
        txt2 <- 100*as.numeric(CamusPar$Val[w])
        if (l > 1) {
          txt1 <- paste0(paste(txt1[1:(l-1)], collapse = ", "), " and ", txt1[l])
          txt2 <- paste0(paste(txt2[1:(l-1)], collapse = ", "), " and ", txt2[l])
        }
        searchTxt_i <- paste0(searchTxt_i, " Results were filtered in Philosopher at ", txt1, " level", c("", "s")[(l>1)+1],
                            " at FDR", c("", "s")[(l>1)+1], " ", txt2, "%", c("", ", respectively")[(l>1)+1], ".")
      }
    }
    # PTM Shepherd
    if (David) {
      searchTxt_i <- gsub("\\.$", ", then PTM_Shepperd was run to localize Open search identified mass shifts." , searchTxt_i)
      # Basic for now, of course eventually this should be expended.
    }
    # MS1 Quant
    if (Kant) {
      txt <- c()
      if (Jon) { txt <- c(txt, paste0("IonQuant", c("", " with match-between-runs turned on")[Matt+1])) }
      if (Fritz) { txt <- c(txt, "FreeQuant") }
      if (length(txt) > 1) { txt <- paste0("both ", paste(txt, collapse = paste0(" and ", c("", "with ")[Matt+1]))) }
      searchTxt_i <- paste0(searchTxt_i, " MS1-level peptide quantitation was performed using ", txt, ".")
    }
    # MS2 Quant
    if (Timothy) { searchTxt_i <- paste0(searchTxt_i, " MS2 reporter intensities were measured using TMT-integrator.") }
    # DiaNN
    if (Diane) {
      Libby <- gsub(topattern("diann.library="), "", grep(topattern("diann.library="), fpWorkflow_i, value = TRUE))
      Libby <- gsub("\\\\", "", gsub("\\\\\\\\", "/", Libby))
      searchTxt_i <- paste0(searchTxt_i, " DIA quantitation was performed using DiaNN ",
                          c("in library-free mode",
                            paste0("using the following library: \"", Libby, "\" (see supplementary information)"))[(nchar(Libby) > 0)+1], ".")
    }
    #
    searchOutputs[[dir_i]] <- list(ev = ev_i,
                                   Samples = samples_i,
                                   Software_param = list(Software = SearchSoft[dir_i],
                                                         Param_file = fpWorkflowFl_i,
                                                         Param = fpWorkflow_i,
                                                         Manifest_file = fpManifestFl_i,
                                                         MinPepSz = minPepSz_i,
                                                         Missed = missed_i,
                                                         Vers = fpVers_i,
                                                         LabelType = labelType_i,
                                                         Umpire = Ump,
                                                         runDiaNN = Diane),
                                   PSMsFls = psmFls_i,
                                   Mods = mods_i,
                                   Fastas = fastas_i,
                                   isDIA = isDIA_i,
                                   Raw_files = list(Full = rawFiles_i,
                                                    Names = rawFiles_i_2),
                                   FracMap = fracMap_i,
                                   MatMet_txt = searchTxt_i)
    if (labelType_i == "Isobaric") {
      tmp <- list(IsobarLab = isobarLab_i,
                  IsobarLabDet = isobarLabDet_i,
                  IsobarLabPrec = isobarLabPrec_i,
                  TMT_annot = fpTMT_annot_i)
      tmp[[isobarLab_i]] <- get(isobarLab_i)
      searchOutputs[[dir_i]]$Software_param[names(tmp)] <- tmp
    }
  }
  if (SearchSoft[dir_i] == "PROTEOMEDISCOVERER") { stop("This part has not yet been re-written for Proteome Discoverer!") }
  if (SearchSoft[dir_i] == "ALPHADIA") { stop("This part has not yet been re-written for alphaDIA!") }
  cat("\n\n")
}
#
# Now combine the outputs of the different search directories
#  - Modifications table
Modifs <- lapply(searchOutputs, function(x) { x$Mods })
modsTst1 <- modsTst2 <- do.call(rbind, lapply(Modifs, function(x) { x[, c("Mark", "Mass shift")] }))
modsTst1 <- aggregate(round(modsTst1$`Mass shift`, 4), list(modsTst1$Mark), unique)
colnames(modsTst1) <- c("Mark", "Mass shifts")
modsTst1$N <- vapply(modsTst1$"Mass shifts", length, 1)
w <- which(modsTst1$N > 1)
if (length(w)) {
  stop("We need to add code to deal with ambiguous cases where the same mark is used in different searches for different PTMs!!!")
}
modsTst2 <- aggregate(modsTst2$Mark, list(round(modsTst2$`Mass shift`, 4)), unique)
colnames(modsTst2) <- c("Mass shifts", "Mark")
modsTst2$N <- vapply(modsTst2$Mark, length, 1)
w <- which(modsTst2$N > 1)
if (length(w)) {
  stop("We need to add code deal with these ambiguous cases where different marks are assigned to the same PTM by different searches!!!")
}
modsTst <- do.call(plyr::rbind.fill, Modifs)
modsTst_a <- aggregate(modsTst$"Full name", list(modsTst$Mark, round(modsTst$"Mass shift", 4)), function(x) { x[[1]] })
colnames(modsTst_a) <- c("Mark", "Mass shifts", "Full name")
k <- c("Site", "Site_long", "AA", "UniMod", "Position")
w <- which(k %in% colnames(modsTst))
modsTst_b <- aggregate(modsTst[, k[w]], list(modsTst$Mark, round(modsTst$"Mass shift", 4)), function(x) { unique(unlist(x)) })
colnames(modsTst_b) <- c("Mark", "Mass shifts", k[w])
modsTst_a[, k[w]] <- modsTst_b[, k[w]]
modsTst_a$Type <- aggregate(modsTst$Type, list(modsTst$Mark, round(modsTst$"Mass shift", 4)), function(x) {
  x <- unique(x)
  x <- x[which(!is.na(x))]
  x <- x[which(x != "NA")]
  if ((length(x) > 1)&&("Variable" %in% x)) {
    x <- "Variable"
  }
  return(x)
})$x
if ("Mass delta" %in% colnames(modsTst)) {
  modsTst_a$"Mass delta" <- aggregate(1:nrow(modsTst), list(modsTst$Mark, round(modsTst$"Mass shift", 4)), function(x) {
    rs <- unique(unlist(modsTst$`Mass delta`[x]))
    rs <- rs[which(!is.na(rs))]
    if (!length(rs)) {
      rs <- unique(unlist(modsTst$`Mass shift`[x]))
      rs <- rs[which(!is.na(rs))]
    }
    if (length(rs)) { rs <- mean(rs) } else { rs <- NA}
    return(rs)
  })$x
}
if ("Max occurences" %in% colnames(modsTst)) {
  modsTst_a$"Max occurences" <- aggregate(modsTst$"Max occurences", list(modsTst$Mark, round(modsTst$"Mass shift", 4)), function(x) {
    max(x, na.rm = TRUE)
  })$x
}
Modifs %<o% modsTst_a
#  - MatMet text templates
MatMet_Search %<o% setNames(lapply(1:l_inDirs, function(x) {
  paste0(c("", paste0("Search #", x, " = '", gsub(".*/", "", inDirs[x]), "':\n"))[(l_inDirs > 1)+1],
         searchOutputs[[x]]$MatMet_txt)
}), inDirs)
# These will need to be updated later
#cat(paste(unlist(MatMet_Search), collapse = "\n"), "\n")
#  - Fastas
fastas_map %<o% list(Original = setNames(lapply(searchOutputs, function(x) { x$Fastas }),
                                         inDirs))
orig_fastas <- fastas <- unique(unlist(fastas_map))
w <- which(!file.exists(fastas))
if (length(w)) {
  dflt <- openxlsx2::read_xlsx(paste0(homePath, "/Default_locations.xlsx"))
  dflt <- dflt$Path[match("Fasta files", dflt$Folder)]
  for (i in w) {
    cat(paste0(" - Fasta file ", fastas[i], " has been relocated since the search",
               c(" was", "es were")[(l_inDirs > 1)+1],
               " run and could not be located automatically.\n   Prompting user...\n"))
    tmp <- c()
    while (!length(tmp)) {
      tmp <- rstudioapi::selectFile(paste0("Select missing fasta file ", fastas[i]),
                                    path = dflt,
                                    filter = "fasta file (*.fasta|*.fas|*.fa|*.faa|*.fasta.fas|*.txt)")
    }
    fastas[i] <- tmp
  }
}
fastas_map$Actual <- setNames(lapply(names(fastas_map$Original), function(nm) {
  fastas[match(fastas_map$Original[[nm]], orig_fastas)]
}), names(fastas_map$Original))
fastas %<o% unique(fastas)
#  - FracMap
FracMap %<o% lapply(searchOutputs, function(x) { x$FracMap })
FracMap <- do.call(plyr::rbind.fill, c(FracMap))
tst1 <- "Parent sample" %in% colnames(FracMap)
tst2 <- "MQ.Exp" %in% colnames(FracMap)
if (tst2) {
  if (!tst1) { FracMap$"Parent sample" <- FracMap$`Raw files name` }
  FracMap$MQ.Exp <- FracMap$"Parent sample"
}
if (tst1) {
  if (!tst2) { FracMap$MQ.Exp <- FracMap$`Raw files name` }
  FracMap$"Parent sample" <- FracMap$MQ.Exp
}
w <- which(is.na(FracMap$"Parent sample"))
if (length(w)) { FracMap$"Parent sample"[w] <- FracMap$MQ.Exp[w] }
w <- which(is.na(FracMap$"Parent sample"))
if (length(w)) { FracMap$"Parent sample"[w] <- FracMap$"Raw files name"[w] }
w <- which(is.na(FracMap$MQ.Exp))
if (length(w)) { FracMap$MQ.Exp[w] <- FracMap$"Parent sample"[w] }
#  - LabelType
#    For now we must enforce single values for LabelType!
LabelType %<o% unique(vapply(searchOutputs, function(x) { x$Software_param$LabelType }, ""))
if (length(LabelType) > 1) {
  stop("It is ok to combine several searches but do it in a responsible manner!\nFor now, it is not possible to combine different searches with different labelling methods.\nWe're not sure that it even makes sense to do it so this is unlikely to change in the future!")
}
#  - MinPepSz
MinPepSz %<o% min(vapply(searchOutputs, function(x) { x$Software_param$MinPepSz }, 1))
#  - Missed
Missed %<o% max(vapply(searchOutputs, function(x) { x$Software_param$Missed }, 1))
#  - Raw files
rawFiles %<o% unique(unlist(lapply(searchOutputs, function(x) { x$Raw_files$Full })))
rawFiles2 %<o% unique(unlist(lapply(searchOutputs, function(x) { x$Raw_files$Names })))
#  - PSM files
PSMsFls_2_Engine %<o% setNames(lapply(1:l_inDirs, function(x) {
  data.frame(File = searchOutputs[[x]]$PSMsFls,
             Search_ID = inDirs[x],
             Search_Engine = paste0(setNames(searchOutputs[[x]]$Software_param$Software, NULL),
                                    "_v", searchOutputs[[x]]$Software_param$Vers))
}), inDirs)
PSMsFls %<o% unlist(lapply(PSMsFls_2_Engine, function(x) { x$File }))
searchEngines %<o% setNames(unlist(lapply(PSMsFls_2_Engine, function(x) { x$Search_Engine })), PSMsFls)
# - isDIA:
isDIA %<o% setNames(vapply(searchOutputs, function(x) { x$isDIA }, TRUE),
                    inDirs)
# - QuantUMS:
QuantUMS %<o% setNames(vapply(searchOutputs, function(x) {
  if ("QuantUMS" %in% names(x$Software_param)) { x <- x$Software_param$QuantUMS } else { x <- FALSE }
  return(x)
}, TRUE), inDirs)
#  - PSMs: 
ev %<o% do.call(plyr::rbind.fill, lapply(searchOutputs, function(x) { x$ev }));ev$id <- 1:nrow(ev) ### THE 2ND PART IS SO IMPORTANT I AM PUTTING THEM ON 1 ROW

#tstEvs <- do.call(plyr::rbind.fill, lapply(searchOutputs, function(x) { x$ev[1:10,] }));View(tstEvs)

if (exists("FracMap_reloaded")) {
  m <- match(FracMap$`Raw file`, FracMap_reloaded$`Raw file`)
  tst <- sum(is.na(m))
  gs <- FALSE # Should we remove spaces?
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
    k <- c("Raw file", "Raw files name", "Parent sample", "Fraction", "Use", "PTM-enriched")
    k <- k[which(k %in% colnames(FracMap_reloaded))]
    FracMap[, k] <- FracMap_reloaded[m, k]
    if (LabelType == "LFQ") {
      if ("MQ.Exp" %in% colnames(FracMap_reloaded)) {
        FracMap$"Parent sample" <- FracMap_reloaded$MQ.Exp
      }
    }
    if (gs) {
      mEv <- match(gsub(" ", "", ev$`Raw file path`), gsub(" ", "", FracMap$`Raw file`))
    } else {
      mEv <- match(ev$`Raw file path`, FracMap$`Raw file`)
    }
    tst <- sum(is.na(mEv))
    if (!tst) {
      ev$`Raw file` <- FracMap$`Raw files name`[mEv]
    } else {
      rm(FracMap_reloaded)
    }
  } else {
    rm(FracMap_reloaded)
  }
}
#if (!file.exists(FracMapPath)) {
  write.csv(FracMap, file = FracMapPath, row.names = FALSE)
#}
if (exists("fastas_reloaded")) {
  fastas <- unique(c(fastas, fastas_reloaded))
}

