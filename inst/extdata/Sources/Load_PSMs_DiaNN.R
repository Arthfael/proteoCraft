# Load and process DiaNN PSMs (from .parquet or .tsv report file)
#
labelType_i <- "LFQ" # Could be SILAC in future, will need to be detected!
isDIA_i <- TRUE
diaNN_logFl_i <- list.files(inDirs[dir_i], "\\.log\\.txt", full.names = TRUE)
if (length(diaNN_logFl_i) == 1L) {
  cat(" - DiaNN .log.txt file detected automatically\n")
} else {
  diaNN_logFl_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select DiaNN .log.txt file"),
                                          path = paste0(inDirs[dir_i], "/*.log.txt"),
                                          filter = "DiaNN .log.txt file (*.log.txt)")
}
cat(paste0(" -> DiaNN log file: ", diaNN_logFl_i, "\n"))
diaNN_logFlDir_i <- dirname(diaNN_logFl_i)
diannLog_i <- readr::read_lines(diaNN_logFl_i)
diannCall_i <- grep("^diann.exe ", diannLog_i, ignore.case = TRUE, value = TRUE)[1L]
if (is.na(diannCall_i)) {
  # Sometimes, the diaNN call does not include .exe (when run on Linux), so we can identify the proper row indirectly, based on the fact that afaik it is always immediately after the logical cores
  diannCall_i <- diannLog_i[grep("^Logical CPU cores:", diannLog_i)+1L]
  if (!grepl(" --(cfg|f) ", diannCall_i)) {
    # In case the other failed, we try to identify the call by the present of an input files or files list arguments
    diannCall_i <- grep(" --(cfg|f) ", diannLog_i, ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(diannCall_i)) { stop("Unhandled exception: could not detect diaNN call in this log!") }
  }
}
Args <- unlist(strsplit(diannCall_i, " +--"))
#
psmFls_i <- gsub("^\"|\"$", "", gsub("\\\\", "/", gsub("^out +", "", grep("^out ", Args, value = TRUE))))
tstParquet <- grepl("\\.parquet$", psmFls_i)
tstTSV <- grepl("\\.tsv$", psmFls_i)
if (tstParquet) {
  psmFls_i <- c(sub("\\.parquet$", ".tsv", psmFls_i),
                psmFls_i)
}
if (tstTSV) {
  psmFls_i <- c(psmFls_i,
                sub("\\.tsv$", ".parquet", psmFls_i))
}
w <- which(file.exists(psmFls_i))
if (!length(w)) {
  psmFls_i <- paste0(diaNN_logFlDir_i, sub(".*/", "/", psmFls_i))
  w <- which(file.exists(psmFls_i))
  if (length(w)) {
    cat(" - The DiaNN output folder has been renamed or moved since the search was run, but the PSMs report could be located automatically.\n")
  } else {
    warning(" - The DiaNN output folder has been renamed or moved since DiaNN was run.\nThe PSMs report could not be located automatically.\nPrompting user...\n")
    psmFls_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select DiaNN report file (.tsv or .parquet)"),
                                       path = inDirs[dir_i],
                                       filter = "DiaNN report file s (*.tsv|*.parquet)")
    w <- which(file.exists(psmFls_i))
  }
}
psmFls_i <- psmFls_i[w]
l <- length(psmFls_i)
stopifnot(l > 0L)
#
fastas_i <- gsub("^\"|\"$", "", gsub("\\\\", "/", sub("^fasta ", "", grep("^fasta ", Args, value = TRUE))))
if (!length(fastas_i)) { # In the case of DiaNN, sometimes there is no fasta, only a library!
  msg <- "Select fasta file(s) used to create the library, or a compatible fasta file"
  dflt <- paste0(wd, "/*.fasta")
  filt <- matrix(data = c("fasta", "*.fasta;*.fas;*.fa;*.fasta.fas"), ncol = 2L,
                 dimnames = list("Fasta"))
  fastas_i <- normalizePath(choose.files(dflt, filters = filt), winslash = "/")
}
w <- which(!file.exists(fastas_i))
if (length(w)) {
  for (x in w) {
    if (sub("/[^/]+$", "", fastas_i[x]) == diaNN_logFlDir_i) {
      tmp <- paste0(diaNN_logFlDir_i, sub(".*/", "/", fastas_i[x]))
      if (file.exists(tmp)) {
        warning(paste0(" - ", "Fasta file ", fastas_i[x], " has been relocated since DiaNN was run, but could be located automatically."))
        fastas_i[x] <- tmp
      }
    }
  }
}
tmp <- unlist(strsplit(diannCall_i, " +--"))
rawFiles_i <- unique(gsub("^\"|\"$", "", gsub("\\\\", "/", sub("^f +", "", grep("^f +", tmp, value = TRUE))))) # Raw files existence is not critical for the script!
if (!length(rawFiles_i)) {
  # This is the case where a list of files was provided as argument - usually by FP
  tmp <- gsub("\\\\", "/", sub("^cfg +", "",  grep("^cfg +", Args, value = TRUE)))
  tmp <- sub("-- $", "", tmp) # Don't ask me why this is required...
  if (!file.exists(tmp)) { tmp <- paste0(diaNN_logFlDir_i, sub(".*/", "/", tmp)) }
  if (!file.exists(tmp)) { tmp <- paste0(sub("/[^/]+$", "", diaNN_logFlDir_i), sub(".*/", "/", tmp)) } # The DiaNN-run-through-FragPipe case
  if (file.exists(tmp)) {
    rawFiles_i <- sub("^ *--f +", "", grep("^ *--f +", readr::read_lines(tmp), value = TRUE))
  }
}
rawFiles_i <- grep("\\.exe$", rawFiles_i, value = TRUE, invert = TRUE) # Yes, that happened!
rawFiles_i <- gsub("^ *[\"'] *| *[\"'] *$", "", rawFiles_i)
rawFiles_i_ext <- sub(".*\\.", "", rawFiles_i)
wNtFnd <- which(!file.exists(rawFiles_i))
updtFls <- FALSE
if (ProcessedByUs && length(wNtFnd)) {
  msg <- paste0(" - ", c("Some", "All")[(length(wNtFnd) == length(rawFiles_i))+1L], " MS raw files are missing at the expected location.\n")
  tbl <- data.frame(path = rawFiles_i[wNtFnd], file = sub(".*/", "", rawFiles_i[wNtFnd]), ext = rawFiles_i_ext[wNtFnd])
  dirs <- unique(unlist(lapply(unique(c(inDirs[dir_i], diaNN_logFlDir_i)), \(dir) {
    dir <- c(dir, dirname(dir))
    return(c(dir,
             sub(".*/", paste0(archDir, "/"), dir[1L:2L]),
             archDir))
  })))
  dirs <- dirs[dir.exists(dirs)]
  dirs <- base::setdiff(dirs, ".")
  if (length(dirs) > 1L) {
    dirTst <- lapply(dirs, \(x) { grep(topattern(x), base::setdiff(dirs, x), value = TRUE) })
    dirs <- base::setdiff(dirs, unlist(dirTst))
  }
  if (length(dirs)) {
    dirFls <- setNames(lapply(dirs, .safe_listFls), dirs)
    dirDFls <- setNames(lapply(dirs, \(dr) { grep("\\.d$", .safe_listFls(dr, "dir"), value = TRUE) }), dirs)
    dirFls2 <- setNames(lapply(names(dirFls), \(nm) { basename(dirFls[[nm]]) }), dirs)
    dirDFls2 <- setNames(lapply(names(dirDFls), \(nm) { sub(".*/", "", dirDFls[[nm]]) }), dirs)
    for (dir in dirs) {
      tbl[[dir]] <- apply(tbl[, c("file", "ext")], 1L, \(x) { #x <- tbl[1, c("file", "ext")]
        rs <- NA_character_
        if (x[[2L]] == "d") {
          m <- match(x[[1L]], dirDFls2[[dir]])
          if (!is.na(m)) { rs <- dirDFls[[dir]][m] }
        } else {
          m <- match(x[[1L]], dirFls2[[dir]])
          if (!is.na(m)) { rs <- dirFls[[dir]][m] }
        }
        return(rs)
      })
    }
    tbl$nuLoc <- apply(tbl[, dirs, drop = FALSE], 1L, \(x) {
      x <- x[!is.na(x)]
      if (!length(x)) { x <- NA_character_ }
      return(x)
    })
  } else { tbl$nuLoc <- NA_character_ }
  wY <- which(!is.na(tbl$nuLoc))
  wN <- which(is.na(tbl$nuLoc))
  lY <- length(wY)
  lN <- length(wN)
  tstY1 <- (lY > 0L)+1L # Did we manage to locate some files automatically?
  tstY2 <- (lY > 1L)+1L # ... more than 1?
  tstN1 <- (lN > 0L)+1L # Did we fail to locate some?
  tstN2 <- (lN > 1L)+1L # ... more than 1?
  if (tstY1 == 2L) {
    updtFls <- TRUE
    msg <- paste0(msg, "   However, the script automatically detected ",
                  c(c("this", "all")[tstY2], "the following")[tstN1], " file", c("", "s")[tstY2],
                  " at the following location", c("", "s")[tstY2], ":\n",
                  paste(paste0(" - ", tbl$nuLoc[wY], "\n"), collapse = ""))
    rawFiles_i_nu <- rawFiles_i
    rawFiles_i_nu[wNtFnd[wY]] <- tbl$nuLoc[wY]
  }
  cat(msg)
  if (tstN1 == 2L) {
    msg2 <- paste0("Select the location of the missing MS file(s) (or cancel if files are unavailable):")
    newDir <- rstudioapi::selectDirectory(msg2, path = archDir)
    if (!is.null(newDir)) {
      newFls <- list.files(newDir, full.names = TRUE, include.dirs = TRUE, recursive = TRUE)
      newFlsTst <- sub(".*/", "", newFls)
      wNA <- which(is.na(tbl$nuLoc))
      tst <- setNames(lapply(tbl$file[wNA], \(fl) { #fl <- tbl$file[wNA[1L]]
        x <- newFls[newFlsTst == fl]
        if (length(x) > 1L) {
          nc <- nchar(x)
          x <- x[nc == min(nc)]
        } else {
          # Backup solution in case we replaced spaces with something else - e.g. when working on Linux
          x <- newFls[(gsub(" ", "-", newFlsTst) == gsub(" ", "-", fl))
                      | (gsub(" ", "_", newFlsTst) == gsub(" ", "_", fl))
                      | (gsub(" ", "", newFlsTst) == gsub(" ", "", fl))]
          if (length(x) > 1L) {
            nc <- nchar(x)
            x <- x[nc == min(nc)]
          }
        }
        return(x)
      }), tbl$file[wNA])
      tst <- tst[lengths(tst) == 1L]
      if (length(tst)) {
        tst <- setNames(unlist(tst), names(tst))
        wY2 <- which(is.na(tbl$nuLoc) & (tbl$file %in% names(tst)))
        wN2 <- which(is.na(tbl$nuLoc) & (!tbl$file %in% names(tst)))
        tstY2 <- (length(wY2) > 1L)+1L
        tstN2 <- (length(wN2) > 1L)+1L
        if (length(wY2)) {
          m <- match(tbl$file[wY2], names(tst))
          tbl$nuLoc[wY2] <- tst[m]
          msg <- "   The user was able to locate "
          msg <- if (length(wN2) == 0L) {
            paste0(msg, c("", "all ")[tstY2], "the missing file", c("", "s")[tstY2],
                   " in directory ", newDir, "\n")
          } else {
            paste0(msg, "the following missing file", c("", "s")[tstY2],
                   " in directory ", newDir, ":\n", paste0(" - ", tbl$file[wY2], collapse = "\n"), "\n")
          }
          cat(msg)
        } else {
          msg <- paste0("   The following file", c("", "s")[tstN2], " could not be located:\n",
                        paste(paste0(" - ", tbl$path[wN], "\n\n"), collapse = ""))
        }
        rawFiles_i_nu <- rawFiles_i
        w <- which(!is.na(tbl$nuLoc))
        updtFls <- length(w) > 0L
        rawFiles_i_nu[w] <- tbl$nuLoc[w]
      }
    }
  }
}
if (exists("ev_DIANN2MQ")) { rm(ev_DIANN2MQ) }
psmsBckpFl_i <- paste0("diaNN PSMs converted to MQ-like format_", dir_i, ".RDS")
if (exists("reloadedBckps") && (psmsBckpFl_i %in% reloadedBckps$File) && file.exists(psmsBckpFl_i)) {
  loadFun(psmsBckpFl_i)
}
if (!exists("ev_DIANN2MQ")) {
  cat(" - Processing PSMs...\n")
  source(parSrc)
  ev_DIANN2MQ <- DIANN_to_MQ(psmFls_i,
                             cl = parClust,
                             log_Fl = diaNN_logFl_i)
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
  m <- lapply(rawFiles_i[w], \(fl) { which(ev_i$`Raw file path` == fl) })
  lTst <- lengths(m)
  wL <- which(lTst > 0L)
  if (length(wL)) {
    m <- listMelt(m[wL], rawFiles_i_nu[w[wL]])
    stopifnot(sum(is.na(m$L1)) == 0L)
    w2 <- which(ev_i$`Raw file path` %in% rawFiles_i[w])
    m2 <- match(w2, m$value)
    stopifnot(sum(is.na(m2)) == 0L)
    ev_i$`Raw file path`[w2] <- m$L1[match(w2, m$value)]
  }
}
mods_i <- ev_DIANN2MQ$PTMs
ev_i$"Raw file path" <- gsub_Rep("\\\\", "/", ev_i$"Raw file path")
rawFiles_i <- unique(ev_i$"Raw file path")
rawFiles_i_2 <- gsub(".*/|\\.(raw|mzX?ML|d)$", "", rawFiles_i, ignore.case = TRUE)
fracMap_i <- data.frame("Raw file" = rawFiles_i,
                        "Raw files name" = rawFiles_i_2,
                        "Parent sample" = rawFiles_i_2,
                        "Fraction" = 1L,
                        "Use" = TRUE,
                        check.names = FALSE)
temp <- unlist(strsplit(diannCall_i, " --min-pep-len "))
minPepSz_i <- if (length(temp) == 2L) { as.integer(gsub(" --.+", "", unlist(strsplit(diannCall_i, " --min-pep-len "))[2L])) } else {
  min(nchar(ev_i$Sequence))
}
temp <- unlist(strsplit(diannCall_i, " --missed-cleavages "))
missed_i <- if (length(temp) == 2L) { as.integer(gsub(" --.+", "", unlist(strsplit(diannCall_i, " --missed-cleavages "))[2L])) } else {
  max(lengths(strsplit(ev_i$Sequence, "R|K")))-1L
}
diannVers_i <- gsub("^DIA-NN | .+$", "", grep("^DIA-NN [0-9]+\\.?[0-9]*", diannLog_i, value = TRUE))
#
# MatMet text template
moult <- (length(rawFiles_i_2) > 1L)+1L
diannCall_i2 <- unlist(strsplit(diannCall_i, " +--"))
FxMd <- gsub(" enabled as a fixed modification$", "", grep(" enabled as a fixed modification$", diannLog_i, value = TRUE))
FxMdC <- grep("Cysteine", FxMd, value = TRUE)
FxMd <- base::setdiff(FxMd, FxMdC)
FxMdC <- gsub("cysteine[ ,-]?", "", FxMdC, ignore.case = TRUE)
FxMdC <- paste0(toupper(substr(FxMdC, 1L, 1L)), substr(FxMdC, 2L, nchar(FxMdC)))
VarMd <- sub("^Modification ", "", sub(" will be considered as variable$", "", grep(" will be considered as variable$", diannLog_i, value = TRUE)))
if (length(VarMd)) {
  VarMd <- Isapply(strsplit(VarMd, " with mass delta | at "), unlist)
  colnames(VarMd) <- c("UniMod", "Delta mass", "AA")
  VarMd$UniMod <- sub("^Uni[Mm]od:?", "", VarMd$UniMod)
  VarMd$AA <- lapply(strsplit(VarMd$AA, ""), \(aa) {
    w <- which(aa == "*")
    lW <- length(w)
    if (lW) {
      stopifnot(lW == 1L, aa[w+1L] == "n") # check assumptions
      aa <- aa[-(w+1L)]
      aa[w] <- "*n"
    }
    return(aa)
  })
  tmp <- listMelt(VarMd$AA, ColNames = c("AA", "row"))
  tmp[, c("UniMod", "Delta mass")] <- VarMd[tmp$row, c("UniMod", "Delta mass")]
  VarMd <- aggregate(tmp$AA, list(tmp$UniMod, tmp$`Delta mass`), list)
  colnames(VarMd) <- c("UniMod", "Delta mass", "AA")
  VarMd$Type <- "Variable"
  VarMd$"Full name" <- apply(VarMd[, c("UniMod", "Delta mass")], 1L, \(x) { #x <- VarMd[1, c("UniMod", "Delta mass")]
    unique(UniMod$Name[(UniMod$UnimodId == x[[1L]]) & (round((UniMod$MonoMass == x[[2L]])*2L)/2 == 0L)])
  })
  w <- which(lapply(VarMd$`Full name`, length) == 0L)
  if (length(w)) {
    VarMd$`Full name`[w] <- apply(VarMd[w, c("Delta mass", "AA")], 1L, \(x) { #x <- VarMd[w[1L], c("Delta mass", "AA")]
      pos <- sort(unlist(x[[2L]]))
      pos[pos == "*n"] <- "_"
      delta <- as.numeric(x[[1L]])
      mods_i <- mods_i[order(abs(mods_i$`Mass shift` - delta)),]
      # in fact, check whether this part of DIANN_to_MQ() should not be improved.
      res <- mods_i$`Full name`[(abs(mods_i$`Mass shift`) - delta < 0.0001) & (vapply(mods_i$AA, \(y) { sum(y %in% pos) > 0L }, TRUE))]
      if (!length(res)) { res <- "unknown_PTM" }
      return(res)
    })
    w <- which(VarMd$`Full name` ==  "unknown_PTM")
    lW <- length(w)
    if (lW) {
      VarMd$`Full name`[w] <- paste0(VarMd$`Full name`, "_", as.character(1L:lW))
    }
  }
  VarMd$Text <- apply(VarMd[, c("Full name", "AA")], 1L, \(x) {
    x[[2L]][x[[2L]] == "n"] <- "N-term"
    x[[2L]][x[[2L]] == "*n"] <- "protein N-term"
    x[[2L]][x[[2L]] == "c"] <- "C-term"
    x[[2L]][x[[2L]] == "*c"] <- "protein C-term"
    paste0(x[[1L]], " (", paste(x[[2L]], collapse = ", "), ")")
  })
}
# Note that it is possible that a modification was searched but not found (so is present in VarMd but not mods_i)
MBR <- grepl("--reanaly[s,z]e", diannCall_i)+1L
Libraries <- sub("^lib +", "", grep("^lib", diannCall_i2, value = TRUE))
nL <- length(Libraries)
LibFree <- (nL == 0L)+1L
searchTxt_i <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in DiaNN")
if (length(diannVers_i) && nchar(diannVers_i)) { searchTxt_i <- paste0(searchTxt_i, " version ", diannVers_i) }
searchTxt_i <- paste0(searchTxt_i,
                      c(paste0(" against the following spectral librar", c("y", "ies")[(nL > 1L)+1L], ": ",
                               paste(Libraries, collapse = "/"), " (see supplementary material), using "),
                        " in library-free mode against ")[LibFree],
                      "TEMPLATELIBTEXT")
searchTxt_i <- paste0(searchTxt_i, " Match-Between-Runs was turned ", c("off", "on")[MBR], ".")
if (length(FxMdC)) {
  txt <- if (length(FxMdC) == 1L) {
    paste0("Fixed cysteine modification was set to ", FxMdC, ".")
  } else {
    paste0(.commAnd(FxMdC), " were included as fixed cysteine modifications.")
  }
  searchTxt_i <- paste0(searchTxt_i, " ", txt)
}
if (length(FxMd)) {
  txt <- if (length(FxMd) == 1L) {
    paste0("An additional fixed modification, ", FxMd, ", was also included.")
  } else {
    paste0("In addition, ", .commAnd(FxMd), " were also included as fixed modifications.")
  }
  searchTxt_i <- paste0(searchTxt_i, " ", txt)
}
if (is.data.frame(VarMd) && nrow(VarMd)) {
  txt <- if (nrow(VarMd) == 1L) {
    paste0(VarMd$Text, " was set as variable modification.")
  } else {
    paste0("Variable modifications were set to ", .commAnd(VarMd$Text), ".")
  }
  searchTxt_i <- paste0(searchTxt_i, " ", txt)
}
tstFDR <- grep("^Output will be filtered at .+ FDR$", diannLog_i, value = TRUE)
FDR <- gsub("^Output will be filtered at | FDR$", "", tstFDR)
if (!grepl(" ?%$", FDR)) { FDR <- as.numeric(FDR)*100 }
if ((FDR >= 0) && (FDR <= 1)) { searchTxt_i <- paste0(searchTxt_i, " Data was filtered at ", FDR, "% FDR.") }
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
