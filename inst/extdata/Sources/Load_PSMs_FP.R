# Load and process FragPipe PSMs
#
fpWorkflowFl_i <- list.files(inDirs[dir_i], "\\.workflow$", full.names = TRUE)
if (length(fpWorkflowFl_i) == 1L) {
  cat(" - FragPipe workflow file detected automatically\n")
} else {
  fpWorkflowFl_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select FragPipe workflow file"),
                                           path = paste0(inDirs[dir_i], "/*.workflow"),
                                           filter = "FragPipe workflow file (*.workflow)")
}
fpManifestFl_i <- list.files(inDirs[dir_i], "\\.fp-manifest$", full.names = TRUE)
if (length(fpManifestFl_i) == 1L) {
  cat(" - FragPipe manifest file detected automatically\n")
} else {
  fpManifestFl_i <- rstudioapi::selectFile(paste0(inDirs[dir_i], ": select FragPipe manifest file"),
                                           path = paste0(inDirs[dir_i], "/*.fp-manifest"),
                                           filter = "FragPipe manifest file (*.fp-manifest)")
}
cat(paste0(" -> FragPipe workflow file: ", fpWorkflowFl_i, "\n"))
cat(paste0(" -> FragPipe manifest file: ", fpManifestFl_i, "\n"))
if (exists("ev_FP2MQ")) { rm(ev_FP2MQ) }
psmsBckpFl_i <- paste0("FragPipe PSMs converted to MQ-like format_", dir_i, ".RDS")
if (exists("reloadedBckps") && (psmsBckpFl_i %in% reloadedBckps$File) && file.exists(psmsBckpFl_i)) {
  loadFun(psmsBckpFl_i)
}
if (!exists("ev_FP2MQ")) {
  cat(" - Processing PSMs...\n")
  source(parSrc)
  ev_FP2MQ <- FP_to_MQ(fpWorkflowFl_i,
                       fpManifestFl_i,
                       FailIfNoQuant = TRUE,
                       cl = parClust)
  cat(" - Saving...\n")
  saveFun(ev_FP2MQ, file = psmsBckpFl_i)
  cat(" -> Done!\n")
}
ev_i <- ev_FP2MQ$Evidence
if ((!"Intensity" %in% colnames(ev_i)) || (!sum(!is.na(ev_i$Intensity))) || (!sum(ev_i$Intensity, na.rm = TRUE))) {
  stop(paste0("You do not appear to have ticked any compatible quantitation method in the Fragpipe search in folder \"",
              inDirs[dir_i], "\"!\nIf you know that the search did include quantitation, then function FP_to_MQ() needs an update.\nOtherwise, rerun the search with quantitation turned on and try again."))
}
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
rawFiles_i_ext <- sub(".*\\.", "", rawFiles_i)
wNtFnd <- which(!file.exists(rawFiles_i))
updtFls <- FALSE
if (ProcessedByUs && length(wNtFnd)) {
  msg <- paste0(" - ", c("Some", "All")[(length(wNtFnd) == length(rawFiles_i))+1L], " MS raw files are missing at the expected location.\n")
  tbl <- data.frame(path = rawFiles_i[wNtFnd], file = sub(".*/", "", rawFiles_i[wNtFnd]), ext = rawFiles_i_ext[wNtFnd])
  dirs <- unique(unlist(lapply(unique(c(inDirs[dir_i], fpManifest_i$Path)), \(dir) {
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
          msg <- if (!length(wN2)) {
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
      }
      rawFiles_i_nu <- rawFiles_i
      w <- which(!is.na(tbl$nuLoc))
      updtFls <- length(w) > 0L
      rawFiles_i_nu[w] <- tbl$nuLoc[w]
    }
  }
}
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
rawFiles_i <- unique(ev_i$`Raw file path`)
rawFiles_i_2 <- gsub(".*[\\\\/]|\\.(raw|mzX?ML|d)$", "", basename(rawFiles_i), ignore.case = TRUE)
fastas_i <- sub("^database\\.db-path=", "", grep("^database\\.db-path=", fpWorkflow_i, value = TRUE))
fastas_i <- gsub("\\\\", "", gsub("\\\\\\\\", "/", fastas_i))
psmFls_i <- gsub("/+", "/", paste0(inDirs[dir_i], "/", fpManifest_i$Samples, "/psm.tsv"))
samples_i <- fpManifest_i$Samples
fracMap_i <- data.frame("Raw file" = rawFiles_i,
                        "Raw files name" = rawFiles_i_2,
                        "MQ.Exp" = samples_i,
                        "Fraction" = 1L,
                        "Use" = TRUE,
                        check.names = FALSE)
if (labelType_i == "Isobaric") {
  fracMap_i$Isobaric.set <- c(as.character(1L), "?")[(nrow(fracMap_i) > 1L)+1L]
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
moult <- (length(rawFiles_i_2) > 1L)+1L
FxMdTbl <- sub(topattern("msfragger.table.fix-mods="), "", grep(topattern("msfragger.table.fix-mods="), fpWorkflow_i, value = TRUE))
FxMdTbl <- unlist(strsplit(FxMdTbl, ";"))
FxMdTbl <- Isapply(strsplit(FxMdTbl, ","), unlist)
colnames(FxMdTbl) <- c("Shift", "Site", "Enabled", "MaxOccurences")
FxMdTbl <- FxMdTbl[FxMdTbl$Enabled == "true",]
FxMdTbl$AAName <- gsub(".+\\(|\\)", "", FxMdTbl$Site)
FxMdTbl$AAName <- sapply(FxMdTbl$AAName, \(x) {
  paste0(toupper(substr(x, 1L, 1L)), substr(x, 2L, nchar(x)))
})
FxMdTbl$AA <- gsub(" \\(.+\\)", "", FxMdTbl$Site)
FxMdTbl$Shift <- as.numeric(FxMdTbl$Shift)
FxMdTbl <- FxMdTbl[FxMdTbl$Shift != 0,]
FxMdTbl$Shift <- as.character(FxMdTbl$Shift)
g <- grep("^-", FxMdTbl$Shift, invert = TRUE)
FxMdTbl$Shift[g] <- paste0("+", FxMdTbl$Shift[g])
FxMdTbl$Name <- apply(FxMdTbl[, c("Shift", "AAName")], 1L, \(x) { paste0(x[[1L]], " (", x[[2L]], ")") })
FxMdC <- FxMdTbl$Name[FxMdTbl$AA == "C"]
FxMd <- FxMdTbl$Name[FxMdTbl$AA != "C"]
VarMd <- sub(topattern("msfragger.table.var-mods="), "", grep(topattern("msfragger.table.var-mods="), fpWorkflow_i, value = TRUE))
VarMd <- unlist(strsplit(VarMd, ";"))
VarMd <- Isapply(strsplit(VarMd, ","), unlist)
colnames(VarMd) <- c("Shift", "Site", "Enabled", "MaxOccurences")
VarMd <- VarMd[VarMd$Enabled == "true",]
for (i in 1L:nrow(AA_table)) { VarMd$Site <- gsub(AA_table$AA[i], paste0(";", AA_table$AA[i], ";"), VarMd$Site) }
VarMd$Site <- gsub("n", ";peptide N-term;", VarMd$Site)
VarMd$Site <- gsub("c", ";peptide C-term;", VarMd$Site)
VarMd$Site <- gsub("\\[\\^", ";protein N-term;", VarMd$Site)
VarMd$Site <- gsub("\\^\\]", ";protein C-term;", VarMd$Site)
VarMd$Site <- strsplit(gsub("^;|;$", "", VarMd$Site), ";;")
VarMd$AAName <- sapply(VarMd$Site, \(x) { #
  x <- unlist(x)
  w <- which(nchar(x) == 1L)
  x[w] <- as.character(AA_table$Name[match(x[w], as.character(AA_table$AA))])
  return(x)
})
VarMd$Shift <- gsub(" +", "", VarMd$Shift)
g <- grep("^-", VarMd$Shift, invert = TRUE)
VarMd$Shift[g] <- paste0("+", VarMd$Shift[g])
VarMd$Text <- apply(VarMd[, c("Shift", "AAName", "Site")], 1L, \(x) {
  x2 <- x[[2L]]
  if (length(x2) > 1L) { x2 <- paste(x[[3L]], collapse = "") }
  return(paste0(x[[1L]], " (", x2, ")"))
})
# Note that it is possible that a modification was searched but not found
#
# FragPipe is very versatile... which will not make this easy...
# Unfortunately, the workflow file does not appear to store the type of Workflow which was preselected in the 2nd tab.
myPat <- c(Ump = "diaumpire.run-diaumpire=",
           Diane = "diann.run-dia-nn=",
           PrTolUp = "msfragger.precursor_mass_upper=",
           PrTolDwn = "msfragger.precursor_mass_lower=",
           Chris = "crystalc.run-crystalc=",
           Moses = "peptide-prophet.run-peptide-prophet=",
           Illy = "percolator.run-percolator=",
           Maud = "ptmprophet.run-ptmprophet=",
           Camus = "phi-report.run-report=",
           CamusPar = "phi-report.filter=",
           David = "ptmshepherd.run-shepherd=",
           Kant = "quantitation.run-label-free-quant=",
           Jon = "ionquant.run-ionquant=",
           Matt = "ionquant.mbr=",
           Fritz = "freequant.run-freequant=",
           Timothy = "tmtintegrator.run-tmtintegrator=",
           Libby = "diann.library=")
myPat <- vapply(myPat, topattern, "")
Ump <- as.logical(toupper(sub(myPat["Ump"], "", grep(myPat["Ump"], fpWorkflow_i, value = TRUE))))
Diane <- as.logical(toupper(sub(myPat["Diane"], "", grep(myPat["Diane"], fpWorkflow_i, value = TRUE))))
isDIA_i <- Diane|Ump # I think theoretically you could use either without the other... although probably for any DIA dataset you will use DiaNN
#stopifnot(((Param$Label == "DIA") & isDIA) | ((Param$Label != "DIA") & (!isDIA))) # Sanity check
# OpenSearch is indirectly identified as cases where unusually high values for tolerances are used
PrTolUp <- sub(myPat["PrTolUp"], "", grep(myPat["PrTolUp"], fpWorkflow_i, value = TRUE))
PrTolDwn <- sub(myPat["PrTolDwn"], "", grep(myPat["PrTolDwn"], fpWorkflow_i, value = TRUE))
# Validations
Chris <- as.logical(toupper(sub(myPat["Chris"], "", grep(myPat["Chris"], fpWorkflow_i, value = TRUE))))
Moses <- as.logical(toupper(sub(myPat["Moses"], "", grep(myPat["Moses"], fpWorkflow_i, value = TRUE))))
Illy <- as.logical(toupper(sub(myPat["Illy"], "", grep(myPat["Illy"], fpWorkflow_i, value = TRUE))))
Maud <- as.logical(toupper(sub(myPat["Maud"], "", grep(myPat["Maud"], fpWorkflow_i, value = TRUE))))
Camus <- as.logical(toupper(sub(myPat["Camus"], "", grep(myPat["Camus"], fpWorkflow_i, value = TRUE))))
David <- as.logical(toupper(sub(myPat["David"], "", grep(myPat["David"], fpWorkflow_i, value = TRUE))))
Kant <- as.logical(toupper(sub(myPat["Kant"], "", grep(myPat["Kant"], fpWorkflow_i, value = TRUE))))
Jon <- as.logical(toupper(sub(myPat["Jon"], "", grep(myPat["Jon"], fpWorkflow_i, value = TRUE))))
Matt <- as.logical(as.numeric(sub(myPat["Matt"], "", grep(myPat["Matt"], fpWorkflow_i, value = TRUE))))
Fritz <- as.logical(toupper(sub(myPat["Fritz"], "", grep(myPat["Fritz"], fpWorkflow_i, value = TRUE))))
Timothy <- as.logical(toupper(sub(myPat["Timothy"], "", grep(myPat["Timothy"], fpWorkflow_i, value = TRUE))))
# (I know, this is really bad, but I was tired, ok!)
#
searchTxt_i <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in FragPipe")
if (length(fpVers_i) && nchar(fpVers_i)) {
  searchTxt_i <- paste0(searchTxt_i, " version ", fpVers_i)
}
searchTxt_i <- paste0(searchTxt_i, " against TEMPLATELIBTEXT")
if (Ump == 1L) {
  searchTxt_i <- sub("\\.$",
                     " As a preliminary step, the DIAUmpire module was run to extract pseudo-MS/MS spectra from DIA spectra files.",
                     searchTxt_i)
}
lFxMdC <- length(FxMdC)
if (lFxMdC) {
  txt <- if (lFxMdC == 1L) {
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
nVarMd <- nrow(VarMd)
if (nVarMd) {
  txt <- if (nVarMd == 1L) {
    paste0(VarMd$Text, " was set as variable modification.")
  } else {
    paste0("Variable modifications were set to ", .commAnd(VarMd$Text), ".")
  }
  searchTxt_i <- paste0(searchTxt_i, " ", txt)
}
w <- which(c(PrTolUp != "20", PrTolDwn != "-20"))
l <- length(w)
if (l) {
  PrTol <- c("lower", "upper")[w]
  PrTol[1L] <- paste0(toupper(substr(PrTol[1L], 1L, 1L)), substr(PrTol[1L], 2L, nchar(PrTol[w])))
  PrTol <- paste(PrTol, collapse = " and ")
  searchTxt_i <- paste0(searchTxt_i, " ", PrTol, " precursor mass tolerance", c("", "s")[l], " w", c("as", "ere")[l], " set to ",
                        paste(c(PrTolDwn, PrTolUp)[w], collapse = " and "), " ppm", c("", ", respectively")[l], ".")
}
# Validations
w <- which(c(Chris, Moses, Illy)) # Normally should only ever be length = 1
l <- length(w)
if (l) {
  txt <- c("Crystal-C", "PeptideProphet", "Percolator")[w]
  if (l > 1L) { txt <- .commAnd(txt) }
  searchTxt_i <- paste0(searchTxt_i, " Peptide identifications were validated using ", txt, ".")
  if (Maud) { searchTxt_i <- sub("\\.$", ", and", searchTxt_i) }
}
if (Maud) { searchTxt_i <- paste0(searchTxt_i, " PTM sites", c(" were", "")[(l > 0L) + 1L], " localized using PTMProphet.") }
# FDR levels
if (Camus) {
  CamusPar <- toupper(sub(myPat["CamusPar"], "", grep(myPat["CamusPar"], fpWorkflow_i, value = TRUE)))
  CamusPar <- unlist(strsplit(CamusPar, "--"))
  CamusPar <- CamusPar[2L:length(CamusPar)]
  CamusPar <- data.frame(Par = CamusPar)
  CamusPar$Val <- gsub("[^ ]+ ", "", CamusPar$Par)
  CamusPar$Par <- gsub(" .+", "", CamusPar$Par)
  lvls <- setNames(c("ion", "psm", "peptide", "protein"), c("ION", "PSM", "PEP", "PROT"))
  w <- which(CamusPar$Par %in% names(lvls))
  l <- length(w)
  if (length(w)) {
    txt1 <- lvls[CamusPar$Par[w]]
    txt2 <- 100*as.numeric(CamusPar$Val[w])
    ind <- (l>1L)+1L
    if (l > 1L) {
      rg <- 1L:(l-1L)
      txt1 <- .commAnd(txt1)
      txt2 <- .commAnd(txt2)
    }
    searchTxt_i <- paste0(searchTxt_i, " Results were filtered in Philosopher at ", txt1, " level", c("", "s")[ind],
                          " at FDR", c("", "s")[ind], " ", txt2, "%", c("", ", respectively")[ind], ".")
  }
}
# PTM Shepherd
if (David) {
  searchTxt_i <- sub("\\.$", ", then PTM_Shepperd was run to localize Open search identified mass shifts." , searchTxt_i)
  # Basic for now, of course eventually this should be expended.
}
# MS1 Quant
if (Kant) {
  txt <- c()
  if (Jon) { txt <- c(txt, paste0("IonQuant", c("", " with match-between-runs turned on")[Matt+1L])) }
  if (Fritz) { txt <- c(txt, "FreeQuant") }
  if (length(txt) > 1L) { txt <- paste0("both ", paste(txt, collapse = paste0(" and ", c("", "with ")[Matt+1L]))) }
  searchTxt_i <- paste0(searchTxt_i, " MS1-level peptide quantitation was performed using ", txt, ".")
}
# MS2 Quant
if (Timothy) { searchTxt_i <- paste0(searchTxt_i, " MS2 reporter intensities were measured using TMT-integrator.") }
# DiaNN
if (Diane) {
  Libby <- sub(myPat["Libby"], "", grep(myPat["Libby"], fpWorkflow_i, value = TRUE))
  Libby <- gsub("\\\\", "", gsub("\\\\\\\\", "/", Libby))
  searchTxt_i <- paste0(searchTxt_i, " DIA quantitation was performed using DiaNN ",
                        c("in library-free mode",
                          paste0("using the following library: \"", Libby, "\" (see supplementary information)"))[(nchar(Libby) > 0L)+1L], ".")
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
