# Load and process MaxQuant PSMs (evidence.txt file)
#
mqparFl_i <- c()
if (grepl("/combined/txt$", inDirs[dir_i])) {
  mqparFl_i <- list.files(sub("/txt$", "", inDirs[dir_i]), "\\.xml$", full.names = TRUE)
}
if (grepl("/combined$", inDirs[dir_i])) {
  mqparFl_i <- list.files(inDirs[dir_i], "\\.xml$", full.names = TRUE)
}
if (length(mqparFl_i) == 1L) {
  cat(" - MaxQuant mqpar.xml file detected automatically\n")
} else {
  mqparFl_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select MaxQuant mqpar.xml file"),
                                      path = paste0(sub("/txt[^/]*$", "", inDirs[dir_i]), "/*.xml"),
                                      filter = "XML file (*.xml)")
}
#
cat(paste0(" -> mqpar.xml file: ", mqparFl_i, "\n"))
mqpar_i <- readr::read_lines(mqparFl_i)
fastas_i <- gsub(" *</?fastaFilePath> *", "", grep("<fastaFilePath>", mqpar_i, value = TRUE))
fastas_i <- unlist(strsplit(gsub("\\\\", "/", sub("^Fasta file\t", "", fastas_i)), ";"))
mqVers_i <- gsub(" *</?maxQuantVersion> *", "", grep("<maxQuantVersion>", mqpar_i, value = TRUE))
mqFold_i <- paste0("C:/MaxQuant/MaxQuant_", c("", "v"), mqVers_i)
mqFold_i <- mqFold_i[dir.exists(mqFold_i)]
dr <- gsub("\\\\", "/", gsub(" *</?fixed(Combined|Search)Folder> *", "", grep("<fixed(Combined|Search)Folder>", mqpar_i, value = TRUE)))
if ((length(dr) != 1L) || (dr == "")) { dr <- inDirs[dir_i] }
psmFls_i <- paste0(dr, c("", "/txt", "/combined/txt"), "/evidence.txt")
psmFls_i <- psmFls_i[file.exists(psmFls_i)]
if (length(psmFls_i) > 1L) {
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
rawFiles_i <- gsub("\\\\", "/", gsub(" *</?string> *", "", mqpar_i[(g1[1L]+1L):(g1[2L]-1L)]))
rawFiles_i_ext <- sub(".*\\.", "", rawFiles_i)
wNtFnd <- which(!file.exists(rawFiles_i))
updtFls <- FALSE
if (ProcessedByUs && length(wNtFnd)) {
  msg <- paste0(" - ", c("Some", "All")[(length(wNtFnd) == length(rawFiles_i))+1L], " MS raw files are missing at the expected location.\n")
  tbl <- data.frame(path = rawFiles_i[wNtFnd], file = sub(".*/", "", rawFiles_i[wNtFnd]), ext = rawFiles_i_ext[wNtFnd])
  dirs <- unique(unlist(lapply(unique(c(inDirs[dir_i], dirname(mqparFl_i))), \(dir) {
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
      tbl[[dir]] <- apply(tbl[, c("file", "ext")], 1L, \(x) { #x <- tbl[1L, c("file", "ext")]
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
#tmp <- gsub(".*[\\\\/]|\\.(raw|mzX?ML|d)$", "", rawFiles_i, ignore.case = TRUE)
tmp <- unique(ev_i$"Raw file")
rawFiles_i_2 <- gsub(".*/|\\.(raw|mzX?ML|d)$", "", rawFiles_i, ignore.case = TRUE)
stopifnot(sum(!tmp %in% rawFiles_i_2) == 0L)
# We also need a full path column:
if (updtFls) {
  m <- lapply(rawFiles_i_2, \(fl) { which(ev_i$`Raw file` == fl) })
  m <- listMelt(m, rawFiles_i_nu)
  stopifnot(is.integer(m$value),
            sum(!1L:nrow(ev_i) %in% m$value) == 0L)
  m <- m[order(m$value),]
  ev_i$"Raw file path" <- m$L1[match(1L:nrow(ev_i), m$value)]
  rawFiles_i <- rawFiles_i_nu
} else {
  ev_i$"Raw file path" <- rawFiles_i[match(ev_i$"Raw file", rawFiles_i_2)]
}
#
g2 <- grep("</?experiments>", mqpar_i)
g2 <- gsub(" *</?string> *", "", mqpar_i[(g2[1L]+1L):(g2[2L]-1L)])
if (!sum(g2 != "")) { g2[g2 == ""] <- "Exp1" }
g3 <- grep("</?fractions>", mqpar_i)
g3 <- as.integer(sub("^32767$", "1", gsub(" *</?short> *", "", mqpar_i[(g3[1L]+1L):(g3[2L]-1L)])))
g4 <- grep("</?ptms>", mqpar_i)
g4 <- as.logical(toupper(gsub(" *</?boolean> *", "", mqpar_i[(g4[1L]+1L):(g4[2L]-1L)])))
g5 <- grep("</?paramGroupIndices>", mqpar_i)
g5 <- as.integer(gsub(" *</?int> *", "", mqpar_i[(g5[1L]+1L):(g5[2L]-1L)]))
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
if (length(tmp) > 1L) {
  stop("There appear to be different parameter groups with different labeling methods in this MaxQuant search.\tWe have never tested this case, but it would certainly break the code as currently written.")
}
isDIA_i <- grepl("MaxDIA", tmp)
if (tmp == "STANDARD") {
  tmp2 <- as.integer(gsub(" *</?multiplicity> *", "", grep("<multiplicity>", mqpar_i, value = TRUE)))
  if (tmp2 != 1L) {
    stop("SILAC-labeling (or analogous multiple channel labeling MS1-based workflows) is not currently supported by this workflow!")
  }
}
if (tmp == "NEUCODE") { stop("NeuCode SILAC-labeling is not supported by this workflow!") }
m <- match(tmp,
           c("STANDARD", "REPORTERMS2", "REPORTERMS3", "REPORTERIONMS2", "REPORTERIONMS3", "BOXCAR", "TIMSDDA", "TIMSMAXDIA", "BOXCARMAXDIA"))
labelType_i <- c("LFQ", "Isobaric", "Isobaric", "Isobaric", "Isobaric", "LFQ", "LFQ", "LFQ", "LFQ")[m]
if ((scrptType == "noReps") && (labelType_i != "LFQ")) { stop("Currently, this workflow only supports LFQ experiments!") }
if (labelType_i == "Isobaric") { # Currently used only with Reps
  fracMap_i$Isobaric.set <- c(as.character(1), "?")[(nrow(fracMap_i) > 1L)+1L]
  lab <- gsub(" *</?(inter|termi)nalLabel> *", "", grep("<(inter|termi)nalLabel>", mqpar_i, value = TRUE))
  lab <- unique(gsub("-(Nter|Lys)", "-", lab))
  isobarLabDet_i <- gsub(".*-", "", lab)
  isobarLabPrec_i <- unique(gsub("-.*", "", lab))
  isobarLab_i <- sub("[0-9]*plex$", "", isobarLabPrec_i)
  kol <- grep("^Reporter intensity [0-9]+$", colnames(ev_i), value = TRUE)
  stopifnot(length(kol) > 0L)
  val <- as.integer(sub("^Reporter intensity ", "", kol))
  assign(isobarLab_i, val)
} else {
  fracMap_i$"Parent sample" <- fracMap_i$MQ.Exp
  k <- base::setdiff(colnames(fracMap_i), "Parent sample")
  fracMap_i <- fracMap_i[, c("Parent sample", k)]
}
#
minPepSz_i <- as.integer(gsub(" *</?min(PepLen|PeptideLength)> *", "", grep(" *</?min(PepLen|PeptideLength)> *", mqpar_i, value = TRUE)))
missed_i <- as.integer(gsub(" *</?maxMissedCleavages> *", "", grep(" *</?maxMissedCleavages> *", mqpar_i, value = TRUE)))
#
# MatMet text template
moult <- (length(rawFiles_i_2) > 1L)+1L
gFx <- (grep("<fixedModifications>", mqpar_i)+1L):(grep("</fixedModifications>", mqpar_i)-1L)
gVar <- (grep("<variableModifications>", mqpar_i)+1L):(grep("</variableModifications>", mqpar_i)-1L)
FxMd <- gsub(" *</?string> *", "", mqpar_i[gFx])
FxMdC <- grep("\\(C\\)", FxMd, value = TRUE)
FxMd <- base::setdiff(FxMd, FxMdC)
VarMd <- gsub(" *</?string> *", "", mqpar_i[gVar])
searchTxt_i <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in MaxQuant")
if (length(mqVers_i) && nchar(mqVers_i)) {
  searchTxt_i <- paste0(searchTxt_i, " version ", mqVers_i)
}
searchTxt_i <- paste0(searchTxt_i, " against TEMPLATELIBTEXT")
MBR <- gsub(" *</?matchBetweenRuns> *", "", grep("<matchBetweenRuns>", mqpar_i, value = TRUE))
SecPep <- gsub(" *</?secondPeptide> *", "", grep("<secondPeptide>", mqpar_i, value = TRUE))
DPep <- gsub(" *</?dependentPeptides> *", "", grep("<dependentPeptides>", mqpar_i, value = TRUE))
tstFDRs <- grep("Fdr", mqpar_i, value = TRUE)
FDRs <- gsub(" *</?[A-Za-z]*Fdr[A-Za-z]*> *", "", tstFDRs)
w <- which(!FDRs %in% c("False", "True"))
FDRs <- setNames(as.numeric(FDRs[w]), gsub(" *</?|>.+", "", tstFDRs[w]))
FDRs <- FDRs[names(FDRs) != "psmFdrCrosslink"]
if (DPep == "False") { FDRs <- FDRs[names(FDRs) != "dependentPeptideFdr"] }
lFxMdC <- length(FxMdC)
if (lFxMdC) {
  txt <- if (length(FxMdC) == 1L) {
    paste0("Fixed cysteine modification was set to ", FxMdC, ".")
  } else {
    paste0(.commAnd(FxMdC), " were included as fixed cysteine modifications.")
  }
  searchTxt_i <- paste0(searchTxt_i, " ", txt)
}
lFxMd <- length(FxMd)
if (lFxMd) {
  txt <- if (lFxMd == 1L) {
    paste0("An additional fixed modification, ", FxMd, ", was also included.")
  } else {
    paste0("In addition, ", .commAnd(FxMd), " were also included as fixed modifications.")
  }
  searchTxt_i <- paste0(searchTxt_i, " ", txt)
}
lVarMd <- length(VarMd)
if (lVarMd) {
  txt <- if (lVarMd == 1L) {
    paste0(VarMd, " was set as variable modification.")
  } else {
    paste0("Variable modifications were set to ", .commAnd(VarMd), ".")
  }
  searchTxt_i <- paste0(searchTxt_i, " ", txt)
}
tsts <- setNames(as.logical(toupper(c(MBR, SecPep, DPep))), c("Match-between-Runs", "Second Peptide Search", "Dependent Peptides"))
y <- which(tsts); ly <- length(y)
n <- which(!tsts); ln <- length(n)
if (ly) {
  txt <- if (ly == 1L) { paste0(names(tsts)[y], " was set to active") } else {
    paste0(.commAnd(names(tsts)[y]), " were set to active")
  }
  if (ln) {
    txt <- if (ln == 1L) { paste0(txt, ", while ", names(tsts)[n], " was turned off.") } else {
      paste0(txt, ", while ", .commAnd(names(tsts)[n]), " were turned off.")
    }
  } else { txt <- paste0(txt, ".") }
} else {
  txt <- paste0(.commAnd(names(tsts)[n]), " were all turned off.")
}
searchTxt_i <- paste0(searchTxt_i, " ", txt)
tst <- unique(FDRs)
lf <- length(FDRs)
nms <- sub("Fdr$", "s", names(FDRs))
txt <- if (length(tst) == 1L) { paste0("All FDRs were set to ", tst*100, "%.") } else {
  paste0("FDRs were set to: ", .commAnd(paste0(FDRs*100, "% (", nms)), ".")
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
