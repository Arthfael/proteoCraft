############################################
#                                          #
# Histones tail modified peptides analysis #
#                                          #
############################################
#
# Can deal with DiaNN, FragPipe, Skyline or alphaDIA input.
ScriptPath <- normalizePath(gtools::script_file(), winslash = "/")
scrptType <- scrptTypeFull <- "Histones"

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
if (!require(pak)) {
  install.packages("pak")
  library(pak)
}
pkgs <- c("rstudioapi", "qs2", "plyr", "data.table", "openxlsx2", "svDialogs", "parallel", "ggplot2", "viridis",
          "limma", "plotly", "sva", "plotly", "shiny", "htmlwidgets", "shinyWidgets", "shinyjs",
          "shinycssloaders", "DT", "ggrepel", "proteoCraft")
tst <- vapply(pkgs, function(pkg) { require(pkg, character.only = TRUE) }, TRUE)
if (sum(!tst)) {
  w <- pkgs[which(!tst)]
  pkgs_to_inst <- pkgs[w]
  wP <- which(pkgs_to_inst == "proteoCraft")
  if (length(wP)) {
    pkgs_to_inst[wP] <- "Arthfael/proteoCraft"
  }
  pak::pkg_install(pkgs_to_inst)
}
for (pkg in pkgs) {
  library(pkg, character.only = TRUE)
}
homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")

#proteoCraft::load_Bckp()

# For now I always need that one,
# but this could be easily modified
dbFl <- "D:/Fasta_databases/Mus_musculus/Mus_musculus_(C57BL-6J)_UP_20230201_Iso_noDupl_cont.fasta"

#updateMe <- c(TRUE, FALSE)[match(dlg_message("Update the proteoCraft package?", "yesno")$res, c("yes", "no"))]
updateMe <- FALSE
if (updateMe) {
  try({
    unloadNamespace("proteoCraft")
    remove.packages("proteoCraft")
  }, silent = TRUE)
  pak::pkg_install("Arthfael/proteoCraft")
}

# saveFun2 <- function(x, file) {
#   tmp <- paste0("qs::qsavem(", deparse(substitute(x)),
#                 ", file = '", file, "', nthreads = max(c(parallel::detectCores()-1, 1)))")
#   #cat(tmp)
#   eval(parse(text = tmp))
# }
# saveImgFun <- function(file) { # More elegant rewriting - I think
#   args2 <- list(file = file)
#   obj <- objects(.GlobalEnv)
#   obj <- setNames(lapply(obj, get, envir = .GlobalEnv), obj)
#   args2$x <- obj
#   do.call(qs::qsave, args2)
# }
# loadFun <- function(file) {
#   qs::qload(file, env = globalenv(), nthreads = max(c(parallel::detectCores()-1, 1)))
# }
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")

# Set Shiny options, load functions for creating a Word report, create Excel styles
Src <- paste0(libPath, "/extdata/R scripts/Sources/ShinyOpt_Styles_and_Report.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Create parallel processing cluster
parSrc <- paste0(libPath, "/extdata/R scripts/Sources/make_check_Cluster.R")
#rstudioapi::documentOpen(Src)
source(parSrc, local = FALSE)

# Mode
inputTypes <- c("DiaNN", "FragPipe", "Skyline", "alphaDIA")
opt <- sapply(inputTypes, function(x) { paste(c(x, rep(" ", 250 - nchar(x))), collapse = "") })
inputType <- gsub(" +$", "", dlg_list(opt, opt[1], title = "Select type of input (NB: if wanting to combine the output of multiple folders, first run the corresponding combination script)")$res)

if ((exists("inDir"))&&(!is.null(inDir))&&(dir.exists(inDir))) { dfltDir <- inDir } else {
  if ((exists("parDir"))&&(!is.null(parDir))&&(dir.exists(parDir))) { dfltDir <- parDir } else {
    fl <- paste0(homePath, "/Default_locations.xlsx")
    if (file.exists(fl)) {
      tmp <- openxlsx2::read_xlsx(fl)
      dfltDir <- tmp$Path[match("Search folder", tmp$Folder)]
    } else {
      dfltDir <- "D:/groups_temp"
    }
  }
}
if (inputType == "FragPipe") {
  msg <- "Select folder in which FragPipe search result folder(s) are located"
  dflt <- inDir <- rstudioapi::selectDirectory(msg, path = dfltDir)
  parDir <- dirname(inDir)
}
if (inputType == "DiaNN") {
  msg <- "Select DiaNN log file"
  #filt <- matrix(c("DiaNN log file", "*.log.txt"), ncol = 2)
  #diaNN_log_fl <- normalizePath(choose.files(paste0(dfltDir, "/*.log.txt"), msg, multi = FALSE, filt, 1), winslash = "/")
  diaNN_log_fl <- rstudioapi::selectFile(msg,
                                         path = paste0(dfltDir, "/*.log.txt"),
                                         filter = "DiaNN log file (*.log.txt)")
  dflt <- inDir <- dirname(diaNN_log_fl)
  parDir <- dirname(inDir)
}
if (inputType == "Skyline") {
  msg <- "Select Skyline export .tsv file"
  #filt <- matrix(c("Skyline .tsv file", "Skyline .csv file", "*.tsv", "*.csv"), ncol = 2)
  #skyline_fl <- normalizePath(choose.files(paste0(dfltDir, "/*.tsv"), msg, multi = FALSE, filt, 1), winslash = "/")
  skyline_fl <- rstudioapi::selectFile(msg,
                                       path = paste0(dfltDir, "/*.tsv"),
                                       filter = "tsv file exported from Skyline (*.tsv)")
  dflt <- inDir <- dirname(skyline_fl)
  parDir <- dirname(inDir)
}
if (inputType == "alphaDIA") {
  msg <- "Select alphaDIA precursors .tsv file"
  #filt <- matrix(c("alphaDIA .tsv precursors file", "alphaDIA .parquet precursors file", "*.tsv", "*.parquet"), ncol = 2)
  #alphaDIA_fl <- normalizePath(choose.files(paste0(dfltDir, "/*.tsv"), msg, multi = FALSE, filt, 1), winslash = "/")
  alphaDIA_fl <- rstudioapi::selectFile(msg,
                                        path = paste0(dfltDir, "/*.tsv"),
                                        filter = "alphaDIA output tsv file (*.tsv)")
  dflt <- inDir <- dirname(alphaDIA_fl)
  parDir <- dirname(inDir)
}
kount <- 0
while ((kount == 0)||
       (!exists("dstDir"))||
       (!is.character(dstDir))||
       (length(dstDir) != 1)||
       (is.na(dstDir))||
       (!dir.exists(dstDir))) {
  dstDir <- rstudioapi::selectDirectory("Select output directory", path = dflt)
  kount <- kount+1
}
setwd(dstDir)
dtstNm <- gsub(topattern(paste0(dirname(parDir), "/")), "", parDir)

backupFl <- paste0(dstDir, "/Backup.RData")
write(inDir, paste0(dstDir, "/Input search directory.txt")) # In case I reprocess and do not have the backup file
#saveImgFun(backupFl)
#loadFun(backupFl)
if (inputType == "FragPipe") {
  fls <- list.files(inDir, full.names = TRUE, recursive = FALSE)
  FP_WorkflowFl <- grep("/fragpipe\\.workflow$", fls, value = TRUE)
  FP_ManifestFl <- grep("/fragpipe-files\\.fp-manifest$", fls, value = TRUE)
  stopifnot(length(FP_WorkflowFl) > 0, length(FP_ManifestFl) > 0)
  if (length(FP_WorkflowFl) > 1) {
    msg <- "Select FragPipe workflow files in the folder"
    opt <- setNames(sapply(FP_WorkflowFl, function(x) { paste(c(x, rep(" ", 250 - nchar(x))), collapse = "") }),
                    FP_WorkflowFl)
    FP_WorkflowFl <- gsub(" +$", "", dlg_list(opt, opt[1], title = msg)$res)
  }
  if (length(FP_ManifestFl) > 1) {
    msg <- "Select FragPipe manifest files in the folder"
    opt <- setNames(sapply(FP_ManifestFl, function(x) { paste(c(x, rep(" ", 250 - nchar(x))), collapse = "") }),
                    FP_ManifestFl)
    FP_ManifestFl <- gsub(" +$", "", dlg_list(opt, opt[1], title = msg)$res)
  }
  FP2MQ <- FP_to_MQ(FP_WorkflowFl, FP_ManifestFl, FailIfNoQuant = TRUE, cl = parClust)
  ev <- FP2MQ$Evidence
  ev$`Raw file name` <- ev$`Raw file`
  #
  rawFiles <- unique(ev$"Raw file path")
  rawFiles2 <- unique(ev$`Raw file`)
  dbFl <- gsub("\\\\", "", gsub("\\\\\\\\", "/", gsub("^database.db-path=", "", grep("^database.db-path=", FP2MQ$WorkFlow, value = TRUE))))
  #
  ev$Bruker_runID <- as.integer(gsub(".*_|\\.d$", "", ev$`Raw file name`))
  ev$Scan_ID <- do.call(paste, c(ev[, c("Bruker_runID", "MS/MS scan number")], sep = "___"))
  ev$tmp <- round(ev$"m/z", 3)
  ev$Scan_ID_MZ <- do.call(paste, c(ev[, c("Scan_ID", "tmp")], sep = "___"))
  ev$tmp <- NULL
  Exp <- unique(ev$Experiment)
  Modifs <- FP2MQ$PTMs
  FracMap <- FP2MQ$FracMap
}
if (inputType == "DiaNN") {
  diaNN_log <- readLines(diaNN_log_fl)
  diaNN_Call <- grep("diann.exe ", diaNN_log, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(diaNN_Call)) {
    diaNN_Call <- diaNN_log[grep("^Logical CPU cores:", diaNN_log)+1]
  }
  tmpCall <- unlist(strsplit(diaNN_Call, " +--"))
  PSMsFl <- gsub("\\\\", "/", gsub("^out +", "", grep("^out ", tmpCall, value = TRUE)))
  if (!file.exists(PSMsFl)) {
    PSMsFl2 <- paste0(dirname(diaNN_log_fl), gsub(".*/", "/", PSMsFl))
    if (file.exists(PSMsFl2)) {
      warning("The input folder has been renamed since DiaNN was run, but the report file could be located automatically.")
      PSMsFl <- PSMsFl2
    } else {
      warning("The input folder has been renamed since DiaNN was run, and the report file could not be located automatically; prompting user...")
      msg <- "Select DiaNN report file"
      xt <- gsub(".*\\.", "", PSMsFl)
      #PSMsFl <- choose.files(paste0(inDir, "/*.", xt), msg, FALSE)
      PSMsFl <- rstudioapi::selectFile(msg,
                                       path = paste0(inDir, "/*.", xt),
                                       filter = paste0("DiaNN ", xt, " report file (*.", xt, ")"))
    }
  }
  #
  PSMsFl <- normalizePath(PSMsFl, winslash = "/")
  #
  DiaNN2MQ <- try(DIANN_to_MQ(PSMsFl, cl = parClust), silent = TRUE)
  if ("try-error" %in% class(DiaNN2MQ)) { stop("Something went wrong!") }
  #
  Exp <- unique(DiaNN2MQ$Evidence$Experiment)
  if (!length(Exp)) {
    Exp <- unique(DiaNN2MQ$Evidence$"Raw file")
    DiaNN2MQ$Evidence$Experiment <- DiaNN2MQ$Evidence$"Raw file"
  }
  dbFl <- gsub("\\\\", "/", gsub("^fasta ", "", grep("^fasta ", tmpCall, value = TRUE)))
  if (!length(dbFl)) { # In the case of DiaNN, sometimes there is no fasta, only a library!
    dflt <- paste0(inDir, "/*.fasta")
    #msg <- "Select fasta file(s) used to create the library, or a compatible fasta file"
    #filt <- matrix(data = c("fasta", "*.fasta;*.fas;*.fa;*.fasta.fas"), ncol = 2, dimnames = list("Fasta"))
    #dbFl <- normalizePath(choose.files(dflt, filters = filt), winslash = "/")
    msg <- "Select fasta file used to create the library, or a compatible fasta file"
    dbFl <- rstudioapi::selectFile(msg,
                                   path = dflt,
                                   filter = "fasta file (*.fasta|*.fas|*.fa|*.faa|*.fasta.fas|*.txt)")
  }
  DiaNN2MQ$Exp <- Exp
  DiaNN2MQ$dbFl <- dbFl
  ev <- DiaNN2MQ$Evidence
  Modifs <- DiaNN2MQ$PTMs
  ev$"Raw file path" <- gsub("\\\\", "/", ev$"Raw file path")
  rawFiles <- unique(ev$"Raw file path")
  # nuFiles <- rawFiles <- unique(ev$"Raw file path")
  # wN <- which(!file.exists(nuFiles))
  # tmpDr <- inDir
  # while ((length(wN))&&(grepl("/", tmpDr))) {
  #   tmpFls <- paste0(tmpDr, gsub(".*/", "/", nuFiles[wN]))
  #   wY <- which(file.exists(tmpFls))
  #   nuFiles[wN][wY] <- tmpFls[wY]
  #   wN <- which(!file.exists(nuFiles))
  #   tmpDr <- gsub("/[^/]+$", "", tmpDr)  
  # }
  rawFilesExt <- gsub(".*\\.", "", rawFiles)
  wNtFnd <- which(!file.exists(rawFiles))
  updtFls <- FALSE
  if (length(wNtFnd)) {
    msg <- paste0(c("Some", "All")[(length(wNtFnd) == length(rawFiles))+1], " MS files are missing at the expected location.\n")
    tbl <- data.frame(path = rawFiles[wNtFnd], file = gsub(".*/", "", rawFiles[wNtFnd]), ext = rawFilesExt[wNtFnd])
    locs <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
    dirs <- unique(unlist(lapply(unique(c(inDir, dirname(diaNN_log_fl))), function(dir) {
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
      cat(msg)
    }
    cat(msg)
    if (tstN1 == 2) {
      msg2 <- paste0("Select the location of the missing file(s) (or cancel if files are unavailable):")
      newDir <- selectDirectory(msg2, path = wd)
      if (!is.null(newDir)) {
        newFls <- list.files(newDir, full.names = TRUE, include.dirs = TRUE, recursive = TRUE)
        newFlsTst <- gsub(".*/", "", newFls)
        wNA <- which(is.na(tbl$nuLoc))
        tst <- sapply(tbl$file[wNA], function(x) {
          x <- newFls[which(newFlsTst == x)]
          if (length(x) > 1) {
            nc <- nchar(x)
            x <- x[which(nc == min(nc))]
          }
          x
        })
        tst <- tst[which(sapply(tst, length) == 1)]
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
  rawFiles2 <- unique(ev$`Raw file`)
  FracMap <- data.frame("Raw file" = rawFiles,
                        "Raw files name" = rawFiles2,
                        "Parent sample" = rawFiles2,
                        "Fraction" = 1,
                        "Use" = TRUE,
                        check.names = FALSE)
}
if (inputType %in% c("alphaDIA", "Skyline")) {
  if (inputType == "alphaDIA") {
    #evTmp <- data.table::fread(alphaDIA_fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
    alpha2MQ <- proteoCraft::alphaDIA_to_MQ(alphaDIA_fl, cl = parClust)
    #alpha2MQ <- alphaDIA_to_MQ(alphaDIA_fl, cl = parClust)
    ev <- alpha2MQ$Evidence
    Modifs <- alpha2MQ$PTMs
  }
  if (inputType == "Skyline") {
    #evTmp <- data.table::fread(skyline_fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
    Sky2MQ <- proteoCraft::Skyline_to_MQ(skyline_fl, cl = parClust)
    #Sky2MQ <- Skyline_to_MQ(skyline_fl, cl = parClust)
    ev <- Sky2MQ$Evidence
    Modifs <- Sky2MQ$PTMs
  }
  if ("Raw file path" %in% colnames(ev)) {
    rawFiles <- unique(ev$"Raw file path")
    rawFiles2 <- gsub("\\.[^\\.]+$", "", gsub(".*/", "", rawFiles))
  } else {
    rawFiles2 <- rawFiles <- unique(ev$"Raw file")
    ev$"Raw file path" <- ev$"Raw file"
  }
  if (!"Experiment" %in% colnames(ev)) {
    ev$Experiment <- ev$"Raw file path"
  }
  FracMap <- data.frame("Raw file" = rawFiles,
                        "Raw files name" = rawFiles2,
                        "Parent sample" = rawFiles2,
                        "Fraction" = 1,
                        "Use" = TRUE,
                        check.names = FALSE)
}

# Remove reverse database hits
if ("Reverse" %in% colnames(ev)) {
  ev <- ev[which((is.na(ev$Reverse))|(ev$Reverse != "+")),]
}

# Fasta database and annotations 
w <- which(!file.exists(dbFl))
if (length(w)) {
  # Locate files!!!
  for (i in w) { #i <- w[1]
    tmp <- selectFile(paste0("Select missing fasta ", dbFl[i]), path = inDir)
    dbFl[i] <- gsub("^~", normalizePath(Sys.getenv("HOME"), winslash = "/"), tmp)
  }
}
annotFl_tst <- function() {
  !((!exists("annot_Fl", envir = .GlobalEnv))||(length(annot_Fl) != 1)||(!file.exists(annot_Fl))||(gsub(".*\\.", "", annot_Fl) != "txt"))
}
if (!annotFl_tst()) {
  annot_Fl <- gsub("\\.fas.*$", ".txt", dbFl)
}
while (!annotFl_tst()) {
  annot_Fl <- selectFile("Select UniProtKB annotation txt file...", path = dirname(dbFl))
}
prsDB_Fl <- paste0(dstDir, "/Parsed DB.csv")
if (!file.exists(prsDB_Fl)) {
  dbs <- lapply(dbFl, function(x) { Format.DB(x, cl = parClust) })
  db <- plyr::rbind.fill(dbs)
  data.table::fwrite(db, prsDB_Fl, sep = ",", row.names = FALSE, na = "NA")
} else {
  db <- data.table::fread(prsDB_Fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE) 
}
db <- db[grep("^>rev_", db$Header, invert = TRUE),] # Remove reverse database entries
parsedAnnot_Fl <- paste0(dstDir, "/Parsed Annot.rds")
reUseAnnotBckp <- file.exists(parsedAnnot_Fl)
if (reUseAnnotBckp) {
  reUseAnnotBckp <- c(TRUE, FALSE)[match(dlg_message("Parsed annotation backup found in folder. Re-use?", "yesno")$res,
                                         c("yes", "no"))]
}
if (reUseAnnotBckp) {
  annot <- readr::read_rds(parsedAnnot_Fl)
} else {
  annot <- Format.DB_txt(annot_Fl, Features = TRUE, cl = parClust)
  readr::write_rds(annot, parsedAnnot_Fl)
}

# Check peptide-to-protein assignment, and isolate histones peptides
### NB: This doesn't take into consideration N-terminal methionines!!!
msg <- paste0("Update ", inputType, "'s original protein-to-peptides assignments?")
Update_Prot_matches <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
if (Update_Prot_matches) {
  I_eq_L %<o% TRUE
  fl <- paste0(dstDir, "/evmatch.RData")
  if (file.exists(fl)) { loadFun(fl) } else {
    Seq <- unique(ev$Sequence)
    DB <- db
    Src <- paste0(libPath, "/extdata/R scripts/Sources/ProtMatch.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    saveFun(evmatch, fl)
  }
  ev$Search_engine_proteins <- ev$Proteins
  ev$Proteins <- evmatch$Proteins[match(ev$Sequence, evmatch$Sequence)]
  #ev$Proteins <- ev$Search_engine_proteins
  #View(ev[, c("Search_engine_proteins", "Proteins", "Histone(s)", "Modified sequence_verbose")])
  #sum(ev$Search_engine_proteins != ev$Proteins)
  #sum(ev$Search_engine_proteins == ev$Proteins)
}
w <- which((ev$Proteins == "")|(is.na(ev$Proteins)))
lw <- length(w)
if (lw) {
  warning(paste0("Removing ", lw, " PSMs with no match to the database (", round(100*lw/nrow(ev), 1), "%)."))
  ev <- ev[which((ev$Proteins != "")&(!is.na(ev$Proteins))),]
}

#
# The part below is commented as it doesn't work!!!
# The idea was to ensure consistency when combining different searches...
# But when analyzing a single search, it was still removing PSMs!
# The reason: one MS2 can of course yield several PSMs, with close or even identical m/z...
# and all can be valid!!!
# One would need to also check which peaks were matched:
# - get all sequences and all MS2 spectra involved
# - calculate new scores
# - take the best match
# This would be very work intensive to write so has been dropped for now.
#
# When combining different searches, we want to make sure the same spectrum is not assigned
# different peptidoforms.
# Note: this can only work in DDA mode!
# tmp <- as.data.table(ev[, c("Scan_ID_MZ", "Modified sequence_verbose")])
# tmp <- tmp[, list(ModSeq = list(`Modified sequence_verbose`)), by = list(Scan_ID_MZ = Scan_ID_MZ)]
# tmp <- as.data.frame(tmp)
# tmp$L <- sapply(tmp$ModSeq, length)
# mx <- max(tmp$L)
# if (mx > 1) {
#   tstAmbig <- aggregate(tmp$L, list(tmp$L), length)
#   colnames(tstAmbig) <- c("Matches reported per m/z", "Scans")
#   View(tstAmbig)
#   Ambig <- tmp[which(tmp$L > 1),]
#   wAmbig <- which(ev$Scan_ID_MZ %in% Ambig$Scan_ID_MZ)
#   wClear <- which(!ev$Scan_ID_MZ %in% Ambig$Scan_ID_MZ)
#   kol <- c("Scan_ID", "Modified sequence_verbose", "Theoretical m/z", "m/z", "Score", "PEP")
#   stopifnot(sum(!kol %in% colnames(ev)) == 0)
#   tmp2 <- as.data.table(ev[wAmbig, kol])
#   tmp2 <- tmp2[, list(Score = list(Score),
#                       PEP = list(PEP),
#                       `m/z` = list(`m/z`),
#                       `theor m/z` = list(`Theoretical m/z`),
#                       ModSeq = list(`Modified sequence_verbose`)), by = list(Scan_ID = Scan_ID)]
#   tmp2 <- as.data.frame(tmp2)
#   tmp2$BestModSeq <- sapply(1:nrow(tmp2), function(x) {
#     Sc <- unlist(tmp2$Score[x])
#     mdSq <- unlist(tmp2$ModSeq[x])
#     PEP <- unlist(tmp2$PEP[x])
#     delta <- unlist(tmp2$`m/z`[x])-unlist(tmp2$`theor m/z`[x])
#     w <- which(Sc == max(Sc))
#     if (length(w) > 1) {
#       w <- w[which(delta[w] == min(delta[w]))]
#       if (length(w) > 1) {
#         w <- w[which(PEP[w] == min(PEP[w]))]
#       }
#     }
#     mdSq[w]
#   })
#   #class(tmp2$BestModSeq)
#   wU <- which(sapply(tmp2$BestModSeq, length) == 1)
#   wNU <- which(sapply(tmp2$BestModSeq, length) > 1)
#   if (length(wU)) {
#     tmpA <- do.call(paste, c(ev[wAmbig, c("Scan_ID", "Modified sequence_verbose")], sep = "___"))
#     tmpB <- do.call(paste, c(tmp2[wU, c("Scan_ID", "BestModSeq")], sep = "___"))
#     w <- which(tmpA %in% tmpB)
#     rmv <- nrow(ev)-length(c(wClear, w))
#     cat(paste0("Removing ", rmv, " PSMs in favour of a higher quality PSM to the same spectrum but from another search.\n"))
#     ev <- ev[c(wClear, w),]
#   }
#   if (length(wNU)) {
#     warning(paste0(length(wNU), " ambiguous spectra could not be resolved and were removed."))
#     ev <- ev[which(!ev$Scan_ID %in% tmp2$Scan_ID[wNU]),]
#   }
# }

mods <- Modifs
if (!"Delta" %in% colnames(mods)) {
  # Let's allow for my inconsistencies in naming the delta mass columns:
  if ("Mass delta" %in% colnames(mods)) { mods$Delta <- mods$"Mass delta" } else {
    if ("Mass shift" %in% colnames(mods)) { mods$Delta <- mods$"Mass shift" } else {
      if ("Delta mass" %in% colnames(mods)) { mods$Delta <- mods$"Delta mass" } else {
        if ("MW" %in% colnames(mods)) { mods$Delta <- mods$MW }
      }
    }
  }
}
mods <- mods[c(which(is.na(mods$`Old mark`)), which(!is.na(mods$`Old mark`))),]
mods$Delta <- as.numeric(mods$Delta)
mods <- mods[order(mods$Delta, decreasing = FALSE),]
kol <- colnames(mods)
kol <- kol[which(kol != "Old mark")]
modsTst <- do.call(paste, c(mods[, kol], sep = ""))
modsTst <- aggregate(1:nrow(mods), list(modsTst), min)
mods <- mods[modsTst$x,]
mods <- mods[c(which(is.na(mods$`Old mark`)), which(!is.na(mods$`Old mark`))),]
mods <- mods[order(mods$Delta, decreasing = FALSE),]
saveImgFun(backupFl)
#loadFun(backupFl)

# We now need to resolve cases where we got the PTM's name wrong
# - not that we made a mistake, but automation can only go so far.
# Here the issue is that FP takes as input just a delta mass, not a PTM,
# and to convert to MQ-like format we got a name from UniMod which may be incorrect.
# We also need to resolve ambiguous cases, where a single "mark" is assigned to 2 or more distinct PTMs.
# Let's do all of this at once.
msg <- paste0("Do you want to load an external table of PTMs names to force use of specific names?")
opt <- c("Yes                                                                                                                                                           ",
         "No                                                                                                                                                            ")
useExtTbl <- c(TRUE, FALSE)[match(dlg_list(opt, opt[2], title = msg)$res, opt)]
if (useExtTbl) {
  # This file is hard-coded in here: this is the table of mass-shift-to-PTM-name matches.
  # You can edit it but PLEASE stick to the current layout/logic or the script will BREAK!
  modTblFl <- rstudioapi::selectFile("histone_modification_masses.xlsx")
  stopifnot(file.exists(modTblFl))
  #
  useExtTbl <- dlg_message() 
  modTbl <- openxlsx2::read_xlsx(modTblFl)
  #openxlsx2::xl_open(modTblFl)
  colnames(modTbl)[1] <- "Priority"
  modTbl <- modTbl[which(!is.na(modTbl$"molecular weight")),]
  modTbl$`molecular weight` <- as.numeric(modTbl$`molecular weight`)
  modTbl$Priority[which(is.na(modTbl$Priority))] <- ""
  modTbl <- modTbl[which(modTbl$Priority != "Not analysed"),]
  tmp <- listMelt(strsplit(modTbl$`modified residues`, ","), 1:nrow(modTbl), c("Pos", "row"))
  tmp[, c("MW", "Name")] <- modTbl[tmp$row, c("molecular weight", "Modification")]
  tst <- aggregate(tmp$Name, list(tmp$Pos, tmp$MW), function(x) { length(unique(x)) })
  stopifnot(max(tst$x) == 1) # We assume that there is only one mod name for each combination of position and mass shift
  if (inputType %in% c("DiaNN", "Skyline")) {
    if ("Pos" %in% colnames(mods)) {
      mods$Site <- mods$Pos
    } else {
      if ("AA" %in% colnames(mods)) {
        mods$Site <- mods$A
      } else { stop() }
    }
    mods$"Mass delta" <- mods$Delta
  }
  mods2 <- listMelt(mods$Site, 1:nrow(mods), c("Pos", "row"))
  mods2[, c("MW", "Name", "Mark")] <- mods[mods2$row, c("Mass delta", "Full name", "Mark")]
  kol <- c("Pos", "MW")
  tmp$aggr <- do.call(paste, c(tmp[, kol], sep = "___"))
  mods2$aggr <- do.call(paste, c(mods2[, kol], sep = "___"))
  mods2$newName <- ""
  mods2$newMark <- mods2$Mark
  w <- which(mods2$aggr %in% tmp$aggr)
  mods2$newName[w] <- tmp$Name[match(mods2$aggr[w], tmp$aggr)]
  g1 <- grep("-", mods2$newName[w])
  g2 <- grep("-", mods2$newName[w], invert = TRUE)
  if (length(g1)) {
    mods2$newMark[w[g1]] <- sapply(strsplit(tolower(mods2$newName[w[g1]]), "-"), function(x) {
      paste0(substr(x[[1]], 1, 1), substr(x[[2]], 1, 1))
    })
  }
  if (length(g2)) {
    mods2$newMark[w[g2]] <- substr(tolower(mods2$newName[w[g2]]), 1, 2)
  }
  tst <- aggregate(mods2$MW, list(mods2$newMark), unique)
  wMult <- which(sapply(tst$x, length) > 1)
  if (length(wMult)) {
    mods2$newMark_tmp <- mods2$newMark
    for (i in wMult) {
      #i <- wMult[1]
      w <- which(mods2$newMark_tmp == tst$Group.1[i])
      m <- mods2[w,]
      m$Pos[which(sapply(m$Pos, length) == 0)] <- "X"
      r <- 1
      s <- 1:nrow(m); s <- s[which(s != r)]
      tst <- lapply(s, function(x) {
        paste0(tolower(m$Pos[[x]]), substr(m$newMark[[x]], 1, 1))
      })
      tst <- lapply(seq_along(tst), function(x) {
        rs <- tst[[x]]
        rs[which(!rs %in% mods2$newMark)]
      })
      l <- length(tst)
      if (l > 1) {
        for (i in 2:l) {
          tst[[i]] <- tst[[i]][which(!tst[[i]] %in% unlist(tst[1:(i-1)]))]
        }
      }
      tst <- sapply(tst, function(x) {
        x <- unlist(x)
        if (length(x)) { x <- x[1] } else { x <- "That didnae work, did it?" }
        return(x)
      })
      tst2 <- ((tst %in% mods2$newMark)|(tst == "That didnae work, did it?"))
      w2 <- which(!tst2)
      m$newMark[s][w2] <- tst[w2]
      w1 <- which(tst2)
      if (length(w1)) {
        # not tested
        s1 <- s[w1]
        rs <- c()
        kount <- 1
        char <- c(0:9, letters)
        taken <- unique(c(mods2$newMark, m$newMark))
        for (j in s1) {
          tst <- paste0(tolower(m$Pos[s1]), char[kount])
          while (((tst) %in% taken)&&(kount < length(char))) {
            kount <- kount+1
            tst <- paste0(tolower(m$Pos[s1]), char[kount])
          }
          if (kount == length(char)) {
            stop("I am really out of options here, never thought this would go this far! Check the code just in case, this should never happen.")
          } else {
            rs <- c(rs, tst)
          }
        }
        m$newMark[[s1]] <- rs
      }
      mods2[w,] <- m
    }
  }
  tst <- aggregate(mods2$MW, list(mods2$newMark), unique)
  stopifnot(is.numeric(tst$x))
  mods <- mods2; rm(mods2)
  w <- which(mods$newName == "")
  mods$newName[w] <- mods$Name[w]
  mods$newName <- gsub("\\)", ">", gsub("\\(", "<", mods$newName))
  nmKol <- "Name"
  if (!"Modified sequence_verbose" %in% colnames(ev)) {
    ev$"Modified sequence_verbose" <- ev$`Modified sequence`
    nmKol <- "Mark"
  }
  unq <- grep("\\(", unique(ev$`Modified sequence_verbose`), value = TRUE)
  mods2 <- mods
  mods2$newMark <- paste0("(", mods2$newMark, ")")
  mods2$newName <- paste0("(", mods2$newName, ")")
  clusterExport(parClust, list("mods2", "nmKol"), envir = environment())
  unq2 <- as.data.frame(t(parSapply(parClust, unq, function(x) {
    #x <- unq[1]
    x <- y <- proteoCraft::annot_to_tabl(x)[[1]]
    w <- which(x$Annotations != "")
    m <- match(x$Annotations[w], mods2[[nmKol]])
    x$Annotations[w] <- mods2$newMark[m]
    y$Annotations[w] <- mods2$newName[m]
    x <- do.call(paste, c(x, sep = "", collapse = ""))
    y <- do.call(paste, c(y, sep = "", collapse = ""))
    return(c(x, y))
  })))
  tst1 <- gsub("\\)[A-Z]*\\(", "___", gsub("_[A-Z]+_|_[A-Z]*\\(|\\)[A-Z]*_", "", unq2[[1]]))
  tst2 <- gsub("\\)[A-Z]*\\(", "___", gsub("_[A-Z]+_|_[A-Z]*\\(|\\)[A-Z]*_", "", unq2[[2]]))
  tst1 <- unique(unlist(strsplit(tst1, "___")))
  tst2 <- unique(unlist(strsplit(tst2, "___")))
  stopifnot(sum(!tst1 %in% mods$newMark) == 0,
            sum(!tst2 %in% mods$newName) == 0)
  w <- which(ev$`Modified sequence_verbose` %in% unq)
  m <- match(ev$`Modified sequence_verbose`[w], unq)
  ev$`Modified sequence`[w] <- unq2[m, 1]
  ev$`Modified sequence_verbose`[w] <- unq2[m, 2]
  # Done!!!
  # Pfeewwwww...
}

# Edit samples map
Src <- paste0(libPath, "/extdata/R scripts/Sources/hist_Fractions_Map_editor.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

if (!"Old_Experiment" %in% colnames(ev)) { ev$Old_Experiment <- ev$Experiment }
w <- which(!ev$Old_Experiment %in% samplesMap$Sample)
if (length(w)) { 
  ev <- ev[which(ev$Old_Experiment %in% samplesMap$Sample),]
}
m <- match(ev$Old_Experiment, samplesMap$Sample)
ev$Experiment <- samplesMap$"Sample name"[m]
# u <- unique(ev$Experiment)
# if ((length(u) == 1)&&(is.na(u))) {
#   # For when rerunning in bits!
#   ev$Experiment <- ev$Old_Experiment
#   ev$Old_Experiment <- samplesMap$Sample[match(ev$Experiment, samplesMap$"Sample name")]
# }
Exp <- unique(ev$Experiment)
Groups <- unique(samplesMap$Group)

# Edit groups map
Src <- paste0(libPath, "/extdata/R scripts/Sources/hist_Groups_Map_editor.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

statTsts_tst <- aggregate(groupsMap$Reference, list(groupsMap$Comparison_group), function(x) {
  sum(c(FALSE %in% x,
        TRUE %in% x,
        length(unique(x)) == 2))
})
statTsts <- (3 %in% statTsts_tst$x)
if (!statTsts) {
  warning("No statistical tests will be performed!")
}
samplesMap$Comparison_group <- groupsMap$Comparison_group[match(samplesMap$Group, groupsMap$Group)]
samplesMap$Reference <- groupsMap$Reference[match(samplesMap$Group, groupsMap$Group)]
compGrps <- unique(groupsMap$Comparison_group)

isQuant <- ("Intensity" %in% colnames(ev))
if (!isQuant) {
  stop("No Intensity column detected, this is necessary for going further with this workflow...")
}

Nested <- FALSE
if (statTsts) {
  # Nested (paired) replicates?
  opt <- c("Yes                                                                                                                                                               ",
           "No                                                                                                                                                                ")
  Nested <- c(TRUE, FALSE)[match(dlg_list(opt, opt[2], title = "Are the replicates paired?")$res, opt)]
}

Remove0Int <- FALSE
#Remove0Int <- c(FALSE, TRUE)[match(dlg_message("Should we keep peptides with 0 intensity values?", "yesno")$res, c("yes", "no"))]
if (Remove0Int) { ev <- ev[which(ev$Intensity > 0),] }

histDir <- paste0(dstDir, "/Histones_coverage_maps")
if (!dir.exists(histDir)) { dir.create(histDir) } else {
  ls <- list.files(histDir, full.names = TRUE)
  if (length(ls)) {
    msg <- "The coverage maps folder is not empty, clear it?"
    clearDir <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
    if (clearDir) {
      for (fl in ls) { unlink(fl) }
    }
  }
}

# Histone IDs
Hist <- grep("histone", db$Header, ignore.case = TRUE, value = TRUE)
Hist <- grep("ase(,|\\.| )", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("-(binding|like)", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("chaperone", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("ethyl", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("factor", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("non-histone", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("\\(fragment\\)", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
coreHist <- grep("H1|H5|linker", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
HistIDs <- db$`Protein ID`[match(Hist, db$Header)]
coreHistIDs <- db$`Protein ID`[match(coreHist, db$Header)]
histDB <- db[match(HistIDs, db$`Protein ID`),]
#unique(db$`Common Name`[match(Hist, db$Header)])
# Default organism may not be correct:
writeFasta(histDB, paste0(dstDir, "/M_musculus_all_histones.fasta"))
histDB$ID_Nm <- apply(histDB[, c("Protein ID", "Common Name")], 1, paste, collapse = " - ")
histDB$ID_Nm <- gsub("/", "_", histDB$ID_Nm)
clusterExport(parClust, "histDB", envir = environment())

###################################################################################
# Identify the boundary between each histone's ordered domain and disordered tail #
###################################################################################
#
# We will want to normalize to the relatively unmodified, histone-fold domains.
# That way we normalize to what we are interested in (histones) and not affected by tail-borne variable PTMs (what we want to test).
# Note: yes I know there are PTMs also on the structured part. We will either ignore them because they are less many/less variable,
# or we will remove them and their counterpart peptide when identified. TO BE DECIDED!!!
# The annotation tables from UniProtKB identify part of histone sequences as "Disordered": this is the tail, and will be how we identify them.
#
coreHistAnnot <- setNames(lapply(coreHistIDs, function(x) { annot[grsep(x, x = annot$Accession),] }), coreHistIDs)
coreHistAnnot <- coreHistAnnot[which(sapply(coreHistAnnot, nrow) > 0)]
#
#View(coreHistAnnot$Features[[1]])
Features <- lapply(names(coreHistAnnot), function(hist) { #hist <- names(coreHistAnnot)[1] #hist <- "G3UX40"
  x <- coreHistAnnot[[hist]]$Features[[1]]
  #View(x)
  rs <- list(Histone = hist,
             Sequence = histDB$Sequence[match(hist, histDB$`Protein ID`)])
  rs$Ordered <- list(rs$Sequence)
  if ("Note" %in% colnames(x)) {
    w <- which(sapply(x$Note, function(y) { "Disordered" %in% y }))
    lW <- length(w)
    if (lW) {
      sq <- unlist(strsplit(rs$Sequence, ""))
      rmv <- x[w, c("Start", "End"), drop = FALSE]
      rmv <- apply(rmv, 1, function(x) { x[[1]]:x[[2]] })
      rmv <- as.numeric(unlist(rmv))
      sq[rmv] <- "_"
      sq <- unlist(strsplit(gsub("^_+|_+$", "", paste(sq, collapse = "")), "_+"))
      rs$Ordered <- sq
    }
  }
  return(rs)
})
Ordered <- data.frame(Accession = names(coreHistAnnot),
                      Name = histDB$`Common Name`[match(names(coreHistAnnot), histDB$`Protein ID`)],
                      Sequence = histDB$Sequence[match(names(coreHistAnnot), histDB$`Protein ID`)])
Ordered$Ordered <- lapply(Features, function(x) { x$Ordered })
Ordered2 <- Ordered
Ordered2$Ordered <- sapply(Ordered2$Ordered, paste, collapse = " / ")
write.csv(Ordered2, file = paste0(dstDir, "/Ordered_histone_domains.csv"), row.names = FALSE)

if (useExtTbl) {
  tmp <- data.frame(ModSeq =  ev$"Modified sequence_verbose")
  tmp$tst <- gsub("^_[A-Z]+_$|^_[A-Z]*\\(|\\)[A-Z]*_$", "",
                  gsub("\\)[A-Z]+\\(", "___", tmp$ModSeq))
  tst <- unlist(strsplit(tmp$tst, "___"))
  tst <- aggregate(tst, list(tst), length)
  nrow(tst) # Number of identified PTMs
  length(unique(mods$newName)) # Number of searched PTMs
  stopifnot(sum(!tst$Group.1 %in% mods$newName) == 0) # Checking that the names match
}

tmp <- strsplit(ev$Proteins, ";")
tmp <- proteoCraft::listMelt(tmp, 1:nrow(ev), c("ID", "row"))
tmp <- tmp[which(tmp$ID %in% histDB$`Protein ID`),]
tmp <- data.table::as.data.table(tmp)
tmp <- tmp[, list(IDs = list(ID)), by = list(row = row)]
tmp <- as.data.frame(tmp)
tmp$IDs <- parSapply(parClust, tmp$IDs, paste, collapse = ";")
ev$"Histone(s)" <- ""
ev$"Histone(s)"[tmp$row] <- tmp$IDs
#View(ev[tmp$row, c("Proteins", "Histone(s)")])
sum(ev$Proteins != ev$`Histone(s)`)
sum(ev$Proteins == ev$`Histone(s)`)

tmp <- ev[, c("Histone(s)", "Proteins")]
tmp$"Proteins" <- strsplit(tmp$"Proteins", ";")
tmp$"Histone(s)" <- strsplit(tmp$"Histone(s)", ";")
tst <- apply(tmp[, c("Histone(s)", "Proteins")], 1, function(x) {
  sum(!x[[1]] %in% x[[2]])
})
w <- which(tst > 0)
#View(ev[w, c("Modified sequence_verbose", "Histone(s)", "Proteins")])
#View(ev[grep("phospho", ev$`Modified sequence_verbose`, ignore.case = TRUE), c("Modified sequence_verbose", "Histone(s)", "Proteins")])

invisible(parLapply(parClust, 1:N.clust, function(x) { rm(list = ls());gc() }))
saveImgFun(backupFl)
#loadFun(backupFl)

# Summarize over charge
if (!"Raw file path" %in% colnames(ev)) {
  stopifnot("Raw file" %in% colnames(ev))
  ev$"Raw file path" <- ev$`Raw file`
}
if (!"Potential contaminant" %in% colnames(ev)) { ev$"Potential contaminant" <- "" }
if (!"Reverse" %in% colnames(ev)) { ev$"Reverse" <- "" }
kol1 <- c("Modified sequence", "Raw file path", "Experiment")
kol2 <- c("Intensity", "Charge", "m/z", "Score")
#c(kol1, kol2)[which(!c(kol1, kol2) %in% colnames(ev))]
if (sum(!c(kol1, kol2) %in% colnames(ev))) {
  stop("Some columns necessary for going further with this workflow are missing...")
}
# We will allow summing if several raw files are assigned the same Experiment
evSum <- as.data.table(ev[, c(kol1, kol2)])
evSum <- evSum[,
               list(Intensity = sum(Intensity, na.rm = TRUE),
                    Charge = paste(Charge, collapse = ";"),
                    `m/z` = paste(`m/z`, collapse = ";"),
                    Score = max(Score, na.rm = TRUE)),
               by = list(`Modified sequence` = `Modified sequence`,
                         #`Raw file path` = `Raw file path`,
                         Experiment = Experiment)]
evSum <- as.data.frame(evSum)
kol3 <- c("Modified sequence_verbose", "Modifications", "Proteins", "Histone(s)", "Mass",
          "Missed cleavages", "Potential contaminant", "Reverse")
kol3 <- kol3[which(kol3 %in% colnames(ev))]
evSum[, kol3] <- ev[match(evSum[["Modified sequence"]], ev[["Modified sequence"]]), kol3]
kol4 <- c("Modified sequence", "Modified sequence_verbose", "Modifications", "Intensity",
          "Raw file path", "Experiment", "Proteins", "Histone(s)", "Mass", "Charge", "m/z", "Score",
          "Missed cleavages", "Potential contaminant", "Reverse")
kol4 <- kol4[which(kol4 %in% colnames(evSum))]
evSum <- evSum[, kol4]

Exp <- unique(evSum$Experiment)

PSMs_per_Histone <- aggregate(evSum$"Histone(s)", list(evSum$"Histone(s)"), length)
colnames(PSMs_per_Histone) <- c("Histone matches", "Peptides")
PSMs_per_Histone <- PSMs_per_Histone[which(PSMs_per_Histone$`Histone matches` != ""),]
PSMs_per_Histone <- PSMs_per_Histone[order(PSMs_per_Histone$Peptides, decreasing = TRUE),]
View(PSMs_per_Histone)
write.csv(PSMs_per_Histone, file = paste0(dstDir, "/PSMs per Histone.csv"), row.names = FALSE)

w <- which(parSapply(parClust, strsplit(evSum$`Histone(s)`, ";"), function(x) {
  x <- unlist(x)
  x <- x[which(x != "")]
  (length(x) > 1)||(!is.na(x))
}))
histEv <- evSum[w,]
tst <- gsub("^_[A-Z]+_$|^_[A-Z]*\\(|\\)[A-Z]*_$", "", gsub("\\)[A-Z]+\\(", "___",
                                                           histEv$`Modified sequence_verbose`))
tst <- unlist(strsplit(tst, "___"))
tst <- aggregate(tst, list(tst), length)
nrow(tst) # Number of identified PTMs
stopifnot(sum(!unique(unlist(strsplit(histEv$`Histone(s)`, ";"))) %in% histDB$`Protein ID`) == 0)
histEv$Sequence <- ev$Sequence[match(histEv$`Modified sequence`, ev$`Modified sequence`)]

# Check PTMs
wMod <- grep("\\(", histEv$`Modified sequence_verbose`)
tmpMds <- gsub("^_?[A-Z]+_?$|^_?[A-Z]*\\(|\\)[A-Z]*_?$", "",
               gsub("\\)[A-Z]*\\(", "___", unique(histEv$`Modified sequence_verbose`[wMod])))
ptmTst <- unlist(strsplit(tmpMds, "___"))
ptmTst <- aggregate(ptmTst, list(ptmTst), length)
colnames(ptmTst) <- c("PTM", "Nb. of Histone peptidoforms")
ptmTst <- ptmTst[order(ptmTst$`Nb. of Histone peptidoforms`, decreasing = TRUE),]
View(ptmTst)
data.table::fwrite(ptmTst, paste0(histDir, "/Histone PTMs summary.csv"), sep = ",", row.names = FALSE, na = "NA")

# Histones PTM sites count
kol <- c("Sequence", "Modified sequence_verbose", "Proteins")
tmp <- histEv[wMod, kol]
tst <- do.call(paste, c(tmp, sep = "___"))
tst <- unique(tst)
tmp <- Isapply(strsplit(tst, split = "___"), unlist)
colnames(tmp) <- kol
tmp2 <- unlist(strsplit(tmp$Proteins, ";"))
tmp2 <- aggregate(tmp2, list(tmp2), length)
tmp2 <- tmp2[order(tmp2$x, decreasing = TRUE),]
clusterExport(parClust, "tmp2", envir = environment())
tmp$Proteins_by_pep <- parSapply(parClust, strsplit(tmp$Proteins, ";"), function(x) { #x <- strsplit(tmp$Proteins[5], ";")
  x <- unlist(x)
  if (length(x) > 1) {
    o <- tmp2$x[match(x, tmp2$Group.1)]
    x <- x[order(o, decreasing = TRUE)]
    x <- paste(sort(unlist(x)), collapse = ";")
  }
  return(x)
})
tmp$First_Protein <- gsub(";.*", "", tmp$Proteins_by_pep)
tmpSq <- paste0("_", db$Sequence[match(tmp$First_Protein, db$`Protein ID`)])
tmp$Match <- sapply(1:nrow(tmp), function(x) {
  nchar(unlist(strsplit(tmpSq[x], tmp$Sequence[x]))[1])
})
tmp$Tbl <- lapply(tmp$`Modified sequence_verbose`, proteoCraft::annot_to_tabl)
tmp$Sites <- lapply(1:nrow(tmp), function(x) { #x <- 1
  mtch <- tmp$Match[x]
  pr <- tmp$First_Protein[x]
  tbl <- tmp$Tbl[[x]]
  tbl <- tbl[[1]]
  nr <- nrow(tbl)-2
  tbl$Dist <- c(1, 1:nr, nr)
  tbl$Pos <- tbl$Dist + mtch - 1
  tbl <- tbl[which(tbl$Annotations != ""), , drop = FALSE]
  tbl$Pos <- as.character(tbl$Pos)
  tbl$Site <- apply(tbl[, c("Sequence", "Pos", "Annotations")], 1, function(y) {
    paste0(pr, " ", y[[1]], y[[2]], "(", y[[3]], ")")
  })
  return(tbl$Site)
})
hist_Sites <- unique(unlist(tmp$Sites))
hist_Sites <- data.frame(Site = hist_Sites)
hist_Sites$Type <- gsub(".*\\(|\\)", "", hist_Sites$Site)
nSites <- aggregate(hist_Sites$Type, list(hist_Sites$Type), length)
nSites <- nSites[order(nSites$x, decreasing = TRUE),]
colnames(nSites) <- c("PTM", "Nb. of sites")
View(nSites)
data.table::fwrite(nSites, paste0(histDir, "/Histone PTM sites summary.csv"), sep = ",", row.names = FALSE, na = "NA")
cat(paste0("Number of PTM sites = ",
           sum(nSites$`Nb. of sites`[which(!nSites$PTM %in% c("Carbamidomethyl", "Oxidation"))]), "\n"))

saveImgFun(backupFl)
#loadFun(backupFl)

# Coverage maps
kol <- c("Experiment", "Histone(s)", "Intensity", "Modified sequence", "Modified sequence_verbose")
kol %in% colnames(histEv)
tmpEv <- histEv[, kol]
clusterExport(parClust, list("tmpEv", "histDir", "Hist", "Coverage", "histDB"), envir = environment())
#clusterExport(parClust, "Coverage", envir = environment())
for (e in Exp) { #e <- Exp[1]
  w <- which(tmpEv$Experiment == e)
  #w <- which(tmpEv$Experiment %in% Exp)
  if (length(w)) {
    pp <- tmpEv[w,]
    pp <- aggregate(pp$Intensity,
                    list(pp$`Modified sequence`, pp$`Modified sequence_verbose`),
                    sum, na.rm = TRUE)
    colnames(pp) <- c("Modified sequence", "Modified sequence_verbose", "Intensity")
    pp$"Histone(s)" <- tmpEv$"Histone(s)"[match(pp$`Modified sequence`, tmpEv$`Modified sequence`)]
    tmpPP <- listMelt(strsplit(pp$"Histone(s)", ";"), 1:nrow(pp), c("Histone", "row"))
    tmpPP <- aggregate(tmpPP$row, list(tmpPP$"Histone"), c)
    tmpPP <- setNames(lapply(tmpPP$x, function(x) {
      pp[x, c("Modified sequence_verbose", "Intensity")]
    }), tmpPP$Group.1)
    e <- gsub("[\\\\/:\\*\\?\"<>\\|]", "_", e)
    clusterExport(parClust, list("e", "tmpPP"), envir = environment())
    # tst <- sapply(names(tmpPP), function(x) {
    #   g <- grep("nyl", tmpPP[[x]]$`Modified sequence_verbose`, value = TRUE)
    #   if (length(g)) { stop(cat(x, "\n", paste(g, collapse = "\n"), "\n")) }
    # })
    # View(tmpPP[["A0A2R8Y619"]])
    w <- which(histDB$`Protein ID`[match(Hist, histDB$Header)] %in% names(tmpPP))
    if (length(w)) {
      a <- parSapply(parClust, Hist[w], function(hist) {
        #hist <- Hist[w][1]
        #hist <- grep("G3UX40", Hist[w], value = TRUE)
        #hist <- grep("G3UWL7", Hist[w], value = TRUE)
        #hist <- grep("A0A0N4SV66", Hist[w], value = TRUE)
        #hist <- grep("A0A8I4SYN6", Hist[w], value = TRUE)
        #for (hist in Hist[w]) {
        #grep("Q07133", Hist[w])
        m <- match(hist, histDB$Header)
        if (length(m)) {
          nm <- histDB$ID_Nm[m]
          sq <- setNames(histDB$Sequence[m], histDB$ID_Nm[m])
          tmpSq <- tmpPP[[histDB$`Protein ID`[m]]]
          Coverage(sq, tmpSq$`Modified sequence_verbose`, "Heat", FALSE, FALSE,
                   title = paste0(nm, " coverage - ", e),
                   intensities = tmpSq$Intensity,
                   save.path = paste0(histDir, "/", nm, " - ", e),
                   save = c("pdf", "jpeg"))
          fl <- paste0(histDir, "/", nm, " - ", e)
          if (file.exists(fl)) { unlink(fl) }
        } else { stop() }
      }) 
    }
  }
}

allSamples <- unique(samplesMap$"Sample name")

#openwd(histDir)
data.table::fwrite(histEv, paste0(histDir, "/Hist_peptidoforms.tsv"),
                   sep = "\t", row.names = FALSE, na = "NA")
#openxlsx2::xl_open(paste0(histDir, "/Hist_peptidoforms.tsv"))

# Wide peptides table
kol <- colnames(evSum)
kol2 <- c("id", "Raw file", "Raw file2", "Raw file path", "Retention time (start)", "Retention time (end)", "Retention length", "Score", "q-value", "Charge", "Library index", "Library rt",
          "MS2 intensities", "MS/MS scan number", "IM", "IM (predicted)", "iIM", "iIM (predicted)", "Experiment", "Old_Experiment", "Search_engine_proteins", "Retention time",
          "Retention time (predicted)", "iRT", "iRT (predicted)", "Intensity", "PEP", "Theoretical m/z", "m/z", "Mass", "Delta score")
kol1 <- kol[which(!kol %in% kol2)]
tmp <- do.call(paste, c(evSum[, kol1], sep = "_<>_"))
tmp <- aggregate(1:nrow(evSum), list(tmp), min)
pep <- evSum[tmp$x, kol1]
#
# Quantitative values
intType <- "Original"
intRoot <- setNames("Intensity - ", intType)
logIntRoot <- setNames("log10(int.) - ", intType)
ratRoot <- "log2(Ratio) - "
#
tmp <- evSum[, c("Experiment", "Modified sequence", "Intensity")]
exports <- list("allSamples", "samplesMap", "tmp", "intRoot", "is.all.good")
clusterExport(parClust, exports, envir = environment())
clusterCall(parClust, function() library(data.table))
tmp4 <- setNames(parLapply(parClust, allSamples, function(smpl) { #smpl <- allSamples[1] #smpl <- allSamples[4]
  w <- which(tmp$Experiment == smpl)
  tmp2 <- data.frame(mod = NA, Intensity = NA)
  if (length(w)) {
    tmp2 <- data.table(mod = tmp[w, "Modified sequence"],
                       Intensity = tmp[w, "Intensity"])
    tmp2$Intensity[which(!is.all.good(tmp2$Intensity, 2))] <- NA
    tmp2$Intensity[which(tmp2$Intensity <= 0)] <- NA
    tmp2 <- tmp2[, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(mod)]
    tmp2 <- as.data.frame(tmp2)
  }
  return(tmp2)
}), allSamples)
for (smpl in allSamples) { #smpl <- allSamples[1] #smpl <- allSamples[4]
  tmp <- tmp4[[smpl]]
  pep[[paste0(intRoot[intType], smpl)]] <- NA
  w <- which(pep$"Modified sequence" %in% tmp$mod)
  m <- match(pep$"Modified sequence"[w], tmp$mod)
  pep[w, paste0(intRoot[intType], smpl)] <- tmp$Intensity[m]
}
pep$id <- 1:nrow(pep)
#
labCol <- c("ModSeq", "Proteins", "Name", "Gene names")
w <- which((labCol %in% colnames(ev))&(!labCol %in% colnames(evSum)))
if (length(w)) { evSum[, labCol[w]] <- ev[match(evSum$`Modified sequence`, ev$`Modified sequence`), labCol[w]] }
w <- which((labCol %in% colnames(evSum))&(!labCol %in% colnames(pep)))
if (length(w)) { pep[, labCol[w]] <- evSum[match(pep$`Modified sequence`, evSum$`Modified sequence`), labCol[w]] }
pep$ModSeq <- gsub("^_|_$", "", pep$`Modified sequence`)
pep$Name <- db$`Common Name`[match(gsub(";.*", "", pep$Proteins), db$`Protein ID`)]
w <- which(is.na(pep$Name))
pep$Name[w] <- "No match!"
labCol <- labCol[which(labCol %in% colnames(pep))]
clusterExport(parClust, "labCol", envir = environment())
pep$Label <- parApply(parClust, pep[, labCol, drop = FALSE], 1, function(x) {
  paste0(labCol, ": ", x, collapse = "\n")
})
pep$Sequence <- ev$Sequence[match(pep$`Modified sequence_verbose`, ev$`Modified sequence_verbose`)]
intType <- "Original"
normSummary <- data.frame(Sample = allSamples, Original = sapply(allSamples, function(x) { median(pep[[paste0(intRoot[intType], x)]], na.rm = TRUE) }))

#
saveImgFun(backupFl)
#loadFun(backupFl)

# Batch correction:
kol <- colnames(samplesMap)
kol <- kol[which(!kol %in% c("Sample", "Sample name", "Group", "Comparison_group",  "Reference"))]
kol <- kol[which(vapply(kol, function(k) {
  length(unique(samplesMap[[k]])) > 1
}, TRUE))]
if (length(kol)) {
  dflt <- kol
  if (!Nested) {
    dflt <- dflt[which(dflt != "Replicate")]
  }
  if (exists("myBatches")) {
    dflt <- myBatches
  }
  kol <- setNames(vapply(kol, function(k) { paste(c(k, rep(" ", max(c(200, nchar(k))))), collapse = "") }, ""),
                  kol)
  dflt <- setNames(kol[match(dflt, names(kol))], dflt)
  myBatches <- dlg_list(kol, dflt,
                        title = "Select any known batch variables to explore sequentially correcting against",
                        multiple = TRUE)$res
  myBatches <- names(kol)[match(myBatches, kol)]
  combatNorm <- length(myBatches)
  if (combatNorm) {
    #
    # Check for synonymous batches
    tmp <- samplesMap[, myBatches, drop = FALSE]
    for (btch in myBatches) {
      tmp[[btch]] <- as.numeric(as.factor(tmp[[btch]]))
    }
    lBat <- length(myBatches)
    if (lBat > 1) {
      comb <- gtools::combinations(length(myBatches), 2, myBatches)
      comb <- as.data.frame(comb)
      tst <- apply(comb, 1, function(x) {
        identical(tmp[[x[1]]], tmp[[x[2]]])
      })
      w <- which(tst)
      if (length(w)) {
        btchs2Remove <- unique(comb[w, 2])
        warning(paste0("The following batches are synonymous: ",
                       paste(do.call(paste, c(comb[w,], sep = " and ")), collapse = ", "),
                       ".\nThe following batches will be removed: ", btchs2Remove))
        myBatches <- myBatches[which(!myBatches %in% btchs2Remove)]
      }
    }
    #
    btchDir <- paste0(dstDir, "/Batch correction")
    if (!dir.exists(btchDir)) { dir.create(btchDir, recursive = TRUE) }
    #
    scoresLst <- PCAlyLst <- PCsLst <- list()
    orig <- "original"
    prevBatch <- sapply(myBatches, function(btch) {
      m <- match(btch, myBatches)
      if (m == 1) {
        rs <- orig
      } else {
        rs <- myBatches[m - 1]
      }
      return(rs)
    })
    #
    # Prepare data
    intCols <- paste0(intRoot[intType], allSamples)
    logIntCols <- gsub(proteoCraft::topattern(intRoot[intType]), logIntRoot[intType], intCols)
    edata <- pep[, intCols]
    for (k in intCols) {
      edata[[k]] <- log10(edata[[k]])
    }
    colnames(edata) <- logIntCols
    grps <- samplesMap$Group[match(allSamples, samplesMap$"Sample name")]
    #
    # Impute missing values
    tempImp <- proteoCraft::Data_Impute2(edata, grps)
    impEdata <- tempImp$Imputed_data
    impVal <- tempImp$Positions_Imputed
    #
    # First let's look at the pre-batch correction data with PCA plots:
    tst <- try(pcaBatchPlots, silent = TRUE)
    if ("try-error" %in% class(tst)) {
      source("H:/aRmel_package/proteoCraft/R/pcaBatchPlots.R")
    }
    tmp <- pcaBatchPlots(impEdata,
                         orig,
                         intRoot = logIntRoot[intType],
                         map = samplesMap,
                         SamplesCol = "Sample name")
    PCAlyLst[[orig]] <- tmp$PlotLy
    scoresLst[[orig]] <- tmp$Scores
    PCsLst[[orig]] <- tmp$PCs
    #
    mod0 <- model.matrix(~1, data = samplesMap)
    mod <- model.matrix(~as.factor(samplesMap$Group), data = samplesMap)
    n.sv <- sva::num.sv(impEdata, mod, method = "leek")
    KeepComBatResDflt <- (n.sv > 0)
    if (n.sv == 0) {
      msg <- "We do not estimate that there are any surrogate variables hidden in the data -> no batch correction required. Try correcting nonetheless?"
      opt <- c("No                                                                                                                                                                                                                                                                                                                ",
               "Yes                                                                                                                                                                                                                                                                                                               ")
      combatNorm <- c(FALSE, TRUE)[match(dlg_list(opt, opt[1], title = msg)$res, opt)]
    }
  }
  if (combatNorm) {
    ComBat_Data <- list()
    KeepComBatRes <- c()
    for (btch in myBatches) { #btch <- myBatches[1] #btch <- myBatches[2]
      keepCmBtRs <- KeepComBatResDflt
      combat_edata <- list()
      for (cmpGrp in compGrps) { #cmpGrp <- compGrps[1]  #cmpGrp <- compGrps[2]
        last <- btch
        prev <- prevBatch[btch]
        w <- which(samplesMap$Comparison_group == cmpGrp)
        k <- paste0(logIntRoot[intType], samplesMap$"Sample name"[w])
        map <- samplesMap[w,]
        mod0a <- model.matrix(~1, data = map)
        btchVect <- samplesMap[[btch]][w]
        tst <- aggregate(btchVect, list(btchVect), length)
        if (min(tst$x == 1)) {
          warning(paste0("Could not correct for group ", cmpGrp, ": not enough values"))
          combat_edata[[cmpGrp]] <- impEdata[, k]
        } else {
          combat_edata[[cmpGrp]] <- ComBat(dat = impEdata[, k],
                                           batch = btchVect,
                                           mod = mod0a,
                                           par.prior = TRUE) 
        }
      }
      combat_edata <- do.call(cbind, combat_edata)
      combat_edata <- as.data.frame(combat_edata)
      colnames(combat_edata) <- paste0(logIntRoot[intType], gsub(".* - ", "", colnames(combat_edata)))
      # Plot
      tmp <- try(pcaBatchPlots(combat_edata,
                               btch,
                               intRoot = logIntRoot[intType],
                               map = samplesMap,
                               SamplesCol = "Sample name"), silent = TRUE)
      if (!"try-error" %in% class(tmp)) {
        PCAlyLst[[btch]] <- tmp$PlotLy
        scoresLst[[btch]] <- tmp$Scores
        PCsLst[[btch]] <- tmp$PCs
        #
        appNm <- paste0("Batch corr.: ", prev, " -> ", btch)
        msg <- "Keep results from ComBat batch correction? (untick to cancel correction)"
        PCs <- data.frame("Component" = paste0("PC", as.character(seq_along(PCsLst[[prev]]$sdev))),
                          "Before (%)" = round(100*(PCsLst[[prev]]$sdev)^2 / sum(PCsLst[[prev]]$sdev^2), 0),
                          "After (%)" = round(100*(PCsLst[[btch]]$sdev)^2 / sum(PCsLst[[btch]]$sdev^2), 0))
        if (exists("IHAVERUN")) { rm(IHAVERUN) }
        ui <- fluidPage(
          useShinyjs(),
          setBackgroundColor( # Doesn't work
            color = c(#"#F8F8FF",
              "#EEFAE6"),
            gradient = "linear",
            direction = "bottom"
          ),
          extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
          tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
          titlePanel(tag("u", appNm),
                     appNm),
          br(),
          fluidRow(column(5,
                          checkboxInput("KeepResults", msg, keepCmBtRs),
                          actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                          h4("Recommended criteria:"),
                          h5(HTML("&nbsp;Does the original grouping follow known batches?")),
                          h5(HTML("&nbsp;&nbsp;-> If no: only accept the correction if it improves the apparent grouping of expected sample groups.")),
                          h5(HTML("&nbsp;&nbsp;-> If yes: accept the correction if...")),
                          h5(HTML("&nbsp;&nbsp;&nbsp;- ... it removes the original grouping by batches...")),
                          h5(HTML("&nbsp;&nbsp;&nbsp;- ... or it improves the apparent grouping of expected sample groups.")),
                          withSpinner(DTOutput("PCs")),
                          br(),
                          br(),
                          br()),
                   column(7,
                          withSpinner(plotlyOutput("Before", height = "550px")),
                          withSpinner(plotlyOutput("After", height = "550px")))),
          br(),
          br()
        )
        server <- function(input, output, session) {
          output$Before <- renderPlotly(PCAlyLst[[prev]][[btch]])
          output$After <- renderPlotly(PCAlyLst[[btch]][[btch]])
          output$PCs <- renderDT({ PCs },
                                 FALSE,
                                 escape = FALSE,
                                 class = "compact",
                                 selection = "none",
                                 editable = FALSE,
                                 rownames = FALSE,
                                 options = list(
                                   dom = 't',
                                   paging = FALSE,
                                   ordering = FALSE
                                 ),
                                 callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
          observeEvent(input[["KeepResults"]], {
            assign("keepCmBtRs", as.logical(input[["KeepResults"]]), envir = .GlobalEnv)
          })
          observeEvent(input$saveBtn, {
            assign("keepCmBtRs", as.logical(input[["KeepResults"]]), envir = .GlobalEnv)
            assign("IHAVERUN", TRUE, .GlobalEnv)
            stopApp()
          })
          #observeEvent(input$cancel, { stopApp() })
          session$onSessionEnded(function() { stopApp() })
        }
        runKount <- 0
        while ((!runKount)||(!exists("IHAVERUN"))) {
          eval(parse(text = runApp), envir = .GlobalEnv)
          runKount <- runKount+1
        }
        #
        msg <- paste0(" -> correction of ", btch, "-batch associated effect from ", prev, " ", c("rejec", "accep")[keepCmBtRs+1], "ted.\n")
        if (!keepCmBtRs) {
          last <- prev
          m <- match(btch, myBatches)
          if (m < length(myBatches)) {
            prevBatch[m+1] <- prev
          }
        }
        cat(msg)
        ComBat_Data[[btch]] <- combat_edata
        KeepComBatRes[btch] <- keepCmBtRs
      } else {
        msg <- paste0(" -> correction of ", btch, "-batch associated effect from ", prev, " failed!\n")
        cat(msg)
        ComBat_Data[[btch]] <- NA
        KeepComBatRes[btch] <- keepCmBtRs
      }
    }
    w <- which(KeepComBatRes)
    if (length(w)) {
      # If accepted batch correction, we must now put the data back into our pep object
      #
      btch <- names(KeepComBatRes)[max(which(KeepComBatRes))]
      btchCorrData <- ComBat_Data[[btch]]
      # New expression column root names
      intType <- "ComBat"
      intRoot["ComBat"] <- "ComBat corr. int. - "
      logIntRoot["ComBat"] <- "ComBat corr. log10(int.) - "
      # De-log
      for (k in logIntCols) {
        btchCorrData[[k]] <- 10^(btchCorrData[[k]])
      }
      colnames(btchCorrData) <- gsub(proteoCraft::topattern(logIntRoot["Original"]),
                                     intRoot["ComBat"],
                                     colnames(btchCorrData))
      # Re-introduce missing values
      w <- which(impVal, arr.ind = TRUE)
      btchCorrData[w] <- pep[, intCols][w]
      #View(pep[, intCols]);View(btchCorrData)
      #
      comBatCols <- colnames(btchCorrData)
      pep[, comBatCols] <- btchCorrData[, comBatCols]
      normSummary[[intType]] <- sapply(allSamples, function(x) { median(pep[[paste0(intRoot[intType], x)]], na.rm = TRUE) })
    }
  }
}

# Define Core peptides as peptides matching the non-tail ordered region of histones 
pep$isCore <- FALSE
w <- which(nchar(pep$`Histone(s)`) > 0)
tmp <- listMelt(Ordered$Ordered, Ordered$Accession, c("Ordered", "Accession"))
tmp$Ordered <- paste0("_", unlist(tmp$Ordered), "_")
# We will define a peptide as Core if, for each core histone it matches, it's match falls entirely within the ordered region
pep$tmp <- strsplit(pep$`Histone(s)`, ";")
pep$isCore[w] <- apply(pep[w, c("Sequence", "tmp")], 1, function(x) { #x <- pep[w[1], c("Sequence", "tmp")]
  sq <- tmp$Ordered[match(unlist(x[[2]]), tmp$Accession)]
  sq <- strsplit(sq, unlist(x[[1]]))
  tst <- sapply(sq, length)
  return(min(tst) > 1)
})
pep$tmp <- NULL

# Normalize
quntCol <- paste0(intRoot[intType], allSamples)
if (length(Exp) > 1) {
  normOpt <- c("Histone-fold region of all core histones",
               "Histone-fold region of all core histones, per histone (-> non-histone peptides excluded from analysis!)",
               "All proteins",
               "Do not normalize")
  opt2 <- setNames(sapply(normOpt, function(x) { paste(c(x, rep(" ", 150 - nchar(x))), collapse = "") }), normOpt)
  normMeth <- dlg_list(opt2, opt2[1], title = "Select which normalization method you wish to use")$res
  normMeth <- names(opt2[match(normMeth, opt2)])
  if (normMeth %in% normOpt[1:3]) {
    intRoot["Normalized"] <- "norm. int. - "
    logIntRoot["Normalized"] <- "norm. log10(int.) - "
    quntCol2 <- paste0(intRoot["Normalized"], allSamples)
    if (normMeth %in% normOpt[1:2]) {
      w <- which(pep$isCore)
      if (!length(w)) { stop("Not enough peptides to normalize") }
      pep$isModified <- gsub("^_|_$", "", pep$`Modified sequence`) != pep$Sequence
      tst <- sum(!pep$isModified[w])
      if (tst) {
        msg <- paste0("Should we exclude modified peptides from the set used for normalization?\n(out of a total of ", length(w), " histone ordered region peptidoforms, ", tst, " are unmodified)")
        removeMods <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
        if (removeMods) {
          w <- w[which(!pep$isModified[w])]
        }
      }
      # Global normalisation to histones
      #View(pep[w, c("id", quntCol)])
      tmp <- proteoCraft::AdvNorm.IL(pep[w, c("id", quntCol)], "id", quntCol)
      #View(tmp)
      #tst <- sapply(paste0("AdvNorm.", quntCol), function(x) { summary(tmp[[x]]) });View(tst)
      x <- unlist(pep[w, quntCol])
      x <- x[which(x > 0)]
      M <- median(x, na.rm = TRUE)
      m <- sapply(paste0("AdvNorm.", quntCol), function(x) {
        x <- unlist(tmp[[x]])
        x <- x[which(x > 0)]
        median(is.all.good(as.numeric(x)))
      })
      pep[, quntCol2] <- sweep(pep[, quntCol], 2, m, "/")*M
      if (normMeth == normOpt[2]) {
        # Here we want to normalize each peptide to the parent protein(s)
        # I guess we can calculate normalization factors for each unique protein,
        # then average those if a peptide matches several proteins?
        w <- which(pep$`Histone(s)` != "")
        prt2pep <- proteoCraft::listMelt(strsplit(pep$Proteins[w], ";"), w, c("Protein", "pepRow"))
        prt2pep <- aggregate(prt2pep$pepRow, list(prt2pep$Protein), list)
        prt2pep$id <- sapply(prt2pep$x, function(x) { pep$id[unlist(x)] })
        prt2pep$Core <- sapply(prt2pep$id, function(x) {
          x[which(pep$isCore[match(x, pep$id)])]
        })
        wL <- which(sapply(prt2pep$Core, length) > 0)
        prt2pep <- prt2pep[wL,]
        # Calculate Core peptides-only re-normalization factor at protein level:
        #prt2pep$Norm <- lapply(1:nrow(prt2pep), function(x) { c() })
        tmp <- pep[, c("id", quntCol2)]
        clusterExport(parClust, list("tmp", "quntCol2", "allSamples"), envir = environment())
        prt2pep[, allSamples] <- as.data.frame(t(parSapply(parClust, prt2pep$Core, function(x) { #x <- prt2pep$Core[2]
          tmp1 <- tmp[match(unlist(x), tmp$id),]
          M <- unlist(tmp1[, quntCol2])
          M <- M[which(M > 0)]
          M <- median(M, na.rm = TRUE)
          tmp1 <- proteoCraft::AdvNorm.IL(tmp1, "id", quntCol2)
          #View(tmp)
          #tst <- sapply(paste0("AdvNorm.", quntCol2), function(x) { summary(tmp[[x]]) });View(tst)
          m <- sapply(paste0("AdvNorm.", quntCol2), function(x) {
            x <- unlist(tmp1[[x]])
            x <- x[which(x > 0)]
            median(proteoCraft::is.all.good(as.numeric(x)))
          })
          return(setNames(M/m, allSamples)) 
        })))
        #
        # Turn it into peptide-level normalization factors, averaging as needed
        pep2prt <- listMelt(prt2pep$x, prt2pep$Group.1, c("pep", "prot"))
        pep2prt[, allSamples] <- prt2pep[match(pep2prt$prot, prt2pep$Group.1[wL]), allSamples]
        Norm <- aggregate(pep2prt[, allSamples], list(pep2prt$pep), function(x) { exp(mean(log(x), na.rm = TRUE)) })
        # If we cannot calculate a normalization factor, we need to remove the peptide
        tst <- parApply(parClust, Norm[, allSamples], 1, function(x) { length(proteoCraft::is.all.good(x)) })
        Norm <- Norm[which(tst > 0),]
        #
        # Now we need to remove peptides which cannot be normalized with this method
        wN <- which(!pep$id %in% Norm$Group.1)
        wY <- which(pep$id %in% Norm$Group.1)
        warning(paste0("Removing ", length(w), " non-normalizable peptides (keeping ", length(wY), ")"))
        pep <- pep[wY,]
        Norm <- Norm[match(pep$id, Norm$Group.1),]
        #
        # Finally we apply the normalization
        for (i in seq_along(allSamples)) {
          pep[[quntCol2[i]]] <- pep[[quntCol2[i]]]*Norm[[allSamples[i]]]
        }
      }
    }
    if (normMeth == normOpt[3]) {
      tmp <- proteoCraft::AdvNorm.IL(pep[, c("id", quntCol)], "id", quntCol)
      pep[, quntCol2] <- tmp[, paste0("AdvNorm.", quntCol)]
    }
    intType <- "Normalized"
    normSummary[[intType]] <- sapply(allSamples, function(x) { median(pep[[paste0(intRoot[intType], x)]], na.rm = TRUE) })
    quntCol <- paste0(intRoot[intType], allSamples) # Update
  }
}

# Look at shift induced by normalization for all samples
normTst <- data.frame(Sample = allSamples, NormFact = round(normSummary[[intType]]/normSummary$Original, 3))
max(normTst$NormFact)/min(normTst$NormFact)
write.csv(normTst, paste0(dstDir, "/Norm. summary.csv"), row.names = FALSE)
#tmp <- paste0(do.call(paste, c(normTst, sep = " ->\t\t")), "\n")
#writeClipboard(capture.output(cat(tmp)))

### PCA dimensionality reduction plots
dimRedPlotLy %<o% list()
pep$"Av. log10 abundance" <- rowMeans(pep[, quntCol], na.rm = TRUE)
Rep <- unique(samplesMap$Replicate)
compGrps <- unique(groupsMap$Comparison_group)
pca_types <- c("Global", compGrps)
pca_types2 <- "Histones"
tstHist <- pep$`Histone(s)` > "" # Sometimes we only have histone peptide options, typically when working with Skyline  /n
if (sum(!tstHist)) {
  pca_types2 <- c("Global", pca_types2)
}
for (pca_type in pca_types) {
  for (pca_type2 in pca_types2) {
    kol <- quntCol
    if (pca_type != "Global") {
      smpls <- samplesMap$"Sample name"[which(samplesMap$Comparison_group == pca_type)] 
      tst <- gsub(".* - ", "", kol)
      kol <- kol[which(tst %in% smpls)]
    }
    if (length(kol) > 1) {
      temp <- pep[, kol, drop = FALSE]
      colnames(temp) <- gsub(topattern(intRoot[intType]), "", colnames(temp))
      wOK <- which((apply(temp, 1, function(x) {length(proteoCraft::is.all.good(x))}) == ncol(temp))
                   &((is.na(pep$"Potential contaminant"))|(pep$"Potential contaminant" != "+")))
      temp <- temp[wOK,]
      if (pca_type2 == "Histones") {
        temp <- temp[which(tstHist[wOK]),]
      }
      temp <- sweep(temp, 1, pep$"Av. log10 abundance"[w], "-") # Minimal effect:
      # (the idea is that the final result is not affected by peptide intensity but by peptide logFCs only)
      temp <- temp + rnorm(length(unlist(temp)), 0, 10^-9) # To avoid constant/zero columns, add a small random error
      
      pc <- prcomp(t(temp), scale. = TRUE)
      scores <- as.data.frame(pc$x)
      pv <- round(100*(pc$sdev)^2 / sum(pc$sdev^2), 0)
      pv <- pv[which(pv > 0)]
      pv <- paste0("Components: ", paste(sapply(1:length(pv), function(x) {
        paste0("PC", x, ": ", pv[x], "%")
      }), collapse = ", "))
      scores$Sample <- rownames(scores)
      m <- match(scores$Sample, samplesMap$"Sample name")
      scores$Group <- samplesMap$Group[m]
      scores$Replicate <- samplesMap$Replicate[m]
      #
      scores$"comparison group" <- groupsMap$Comparison_group[match(scores$Group, groupsMap$Group)]
      form <- "~`comparison group`"
      ttl <- "Samples PCA plot"
      if (pca_type2 == "Histones") {
        ttl <- paste0(ttl, " (Hist. only)")
      }
      if (pca_type != "Global") {
        ttl <- paste0(ttl, " - ", pca_type)
      }
      if ("PC2" %in% colnames(scores)) {
        plot <- ggplot(scores, aes(x = PC1, y = PC2)) +
          geom_point(aes(color = Group)) +
          ggpubr::stat_conf_ellipse(aes(fill = Group, color = Group), alpha = 0.1, geom = "polygon") +
          scale_color_viridis_d(begin = 0.25) +
          coord_fixed() +
          geom_hline(yintercept = 0, colour = "black", alpha = 0.5) +
          geom_vline(xintercept = 0, colour = "black", alpha = 0.5) +
          ggtitle(ttl, subtitle = pv) +
          geom_text_repel(aes(label = Sample), size = 2.5, show.legend = FALSE) + 
          theme_bw()
        if (substr(form, 1, 1) == "~") { plot <- plot + facet_wrap(form) } else { plot <- plot + facet_grid(form) }
        #
        scores$Group <- factor(scores$Group)
        if ("PC3" %in% colnames(scores)) {
          Symb <- rep(c("circle", "diamond", "square", "cross", "x"), max(as.numeric(Rep)))[1:max(as.numeric(Rep))]             
          Symb <- Symb[as.numeric(scores$Replicate)]
          plot_lyPCAProt <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3,
                                    color = ~Group, colors = "viridis",
                                    text = ~Sample, type = "scatter3d", mode = "markers",
                                    symbol = I(Symb))
          plot_lyPCAProt <- add_trace(plot_lyPCAProt, scores, x = ~PC1, y = ~PC2, z = ~PC3,
                                      type = "scatter3d", mode = "text", showlegend = FALSE)
          plot_lyPCAProt <- layout(plot_lyPCAProt, title = ttl)
          setwd(dstDir)
          saveWidget(plot_lyPCAProt, paste0(dstDir, "/", ttl, ".html"))
          dimRedPlotLy[["Samples PCA"]] <- plot_lyPCAProt
          #system(paste0("open \"", dstDir, "/", ttl, ".html"))
          # NB: There is currently no way to create a 3D, faceted plot in plotly for R that I know of) 
        } else { poplot(plot, width = 18) }
        suppressWarnings({
          ggsave(paste0(dstDir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
          ggsave(paste0(dstDir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
        })
      } else { warning(paste0(pca_type, " PCA failed, investigate!")) }
    }
  }
}
#openwd(dstDir)

# Table of PTM sites
kol <- c("Sequence", "Modified sequence_verbose", "Proteins")
w <- which((kol %in% colnames(ev))&(!kol %in% colnames(evSum)))
if (length(w)) { evSum[, kol[w]] <- ev[match(evSum$`Modified sequence`, ev$`Modified sequence`), kol[w]] }
w <- which((kol %in% colnames(evSum))&(!kol %in% colnames(pep)))
if (length(w)) { pep[, kol[w]] <- evSum[match(pep$`Modified sequence`, evSum$`Modified sequence`), kol[w]] }
tmp <- pep[, kol]
tmp$`Modified sequence_verbose` <- gsub("^_|_$", "", tmp$`Modified sequence_verbose`)
w <- which(tmp$Proteins == "")
tmp$Proteins[w] <- "No match!"
tst <- do.call(paste, c(tmp, sep = "///"))
tst <- unique(tst)
tmp <- Isapply(strsplit(tst, split = "///"), unlist)
colnames(tmp) <- kol
tmp$pepID <- pep$id
tmp2 <- unlist(strsplit(tmp$Proteins, ";"))
tmp2 <- aggregate(tmp2, list(tmp2), length)
tmp2 <- tmp2[order(tmp2$x, decreasing = TRUE),]
clusterExport(parClust, "tmp2", envir = environment())
tmp$Proteins_by_pep <- parSapply(parClust, strsplit(tmp$Proteins, ";"), function(x) { #x <- strsplit(tmp$Proteins[5], ";")
  x <- unlist(x)
  if (length(x) > 1) {
    o <- tmp2$x[match(x, tmp2$Group.1)]
    x <- x[order(o, decreasing = TRUE)]
    x <- paste(sort(unlist(x)), collapse = ";")
  }
  return(x)
})
tmp$First_Protein <- gsub(";.*", "", tmp$Proteins_by_pep)
m <- match(tmp$First_Protein, db$`Protein ID`)
w <- which(!is.na(m))
tmp <- tmp[w,]; m <- m[w]
tmpSq <- paste0("_", db$Sequence[m])
tmp$Match <- sapply(1:nrow(tmp), function(x) {
  nchar(unlist(strsplit(tmpSq[x], tmp$Sequence[x]))[1])
})
tmp$Tbl <- parLapply(parClust, paste0("_", tmp$`Modified sequence_verbose`, "_"), proteoCraft::annot_to_tabl)
tmp$Sites <- lapply(1:nrow(tmp), function(x) { #x <- 1
  # Do not parallelize me: it's slower!
  mtch <- tmp$Match[x]
  pr <- tmp$First_Protein[x]
  tbl <- tmp$Tbl[[x]]
  tbl <- tbl[[1]]
  nr <- nrow(tbl)-2
  tbl$Dist <- c(1, 1:nr, nr)
  tbl$Pos <- tbl$Dist + mtch - 1
  tbl <- tbl[which(tbl$Annotations != ""), , drop = FALSE]
  tbl$Pos <- as.character(tbl$Pos)
  tbl$Site <- apply(tbl[, c("Sequence", "Pos", "Annotations")], 1, function(y) {
    paste0(pr, " ", y[[1]], y[[2]], "(", y[[3]], ")")
  })
  return(tbl$Site)
})
Sites <- listMelt(tmp$Sites, tmp$pepID, c("Site", "pepID"))
quntCol <- paste0(intRoot[intType], allSamples) # Update
Sites[, quntCol] <- pep[match(Sites$pepID, pep$id), quntCol]
Sites2 <- Sites
Sites2$pepID <- NULL
Sites2 <- as.data.table(Sites2)
Sites2 <- Sites2[, lapply(.SD, sum, na.rm = TRUE), keyby = Site]
Sites2 <- as.data.table(Sites2)
Sites2 <- Sites2[, lapply(.SD, sum, na.rm = TRUE), keyby = Site] # NB: converts all NAs to 0... equivalent but we had NAs, I will re-introduce them
Sites2 <- as.data.frame(Sites2)
w <- which(as.matrix(Sites2[, quntCol]) == 0, arr.ind = TRUE)
Sites2[, quntCol][w] <- NA
Sites2$Type <- gsub(".*\\(|\\)", "", Sites2$Site)
Sites2$Protein <- gsub(" .*", "", Sites2$Site)
Sites$`Histone(s)` <- pep$`Histone(s)`[match(Sites$pepID, pep$id)]
Sites2$`Histone(s)` <- Sites$`Histone(s)`[match(Sites2$Site, Sites$Site)]
Sites2$Protein_name <- db$`Common Name`[match(Sites2$Protein, db$`Protein ID`)]
Sites <- Sites2; rm(Sites2)
labCol <- c("Site", "Protein_name")
Sites$Label <- apply(Sites[, labCol], 1, function(x) {
  paste0(labCol, ": ", x, collapse = "\n")
})

if (statTsts) {
  # Statistics
  #
  source(parSrc)
  clusterExport(parClust, list("Nested"), envir = environment())
  #reNorm <- FALSE
  #
  # From samplesMap to design matrix
  groupsMap2 <- groupsMap[which(groupsMap$Comparison_group %in% statTsts_tst$Group.1[which(statTsts_tst$x == 3)]),]
  samplesMap2 <- samplesMap[which(samplesMap$Group %in% groupsMap2$Group),]
  samplesMap2$Sample <- samplesMap2$"Sample name"
  compGrps <- unique(groupsMap2$Comparison_group)
  Factors2 <- c("Group")
  Factors <- c("Sample", "Replicate", Factors2)
  for (fct in Factors) {
    assign(paste0(fct, "_"), as.factor(samplesMap2[[fct]]))
  }
  if (Nested) {
    designCall <- paste0("designMatr <- model.matrix(~0 + Replicate_ + ", paste0(Factors2, "_", collapse = " + "), ")")
  } else {
    designCall <- paste0("designMatr <- model.matrix(~0 + ", paste0(Factors2, "_", collapse = " + "), ")")
  }
  #cat(designCall, "\n")
  eval(parse(text = designCall))
  #
  # Define contrasts
  expContrasts <- list()
  for (ratGrp in compGrps) {
    em <- samplesMap2[which(samplesMap2$Comparison_group == ratGrp),]
    grp1 <- unique(em$Group[which(!em$Reference)])
    grp0 <- unique(em$Group[which(em$Reference)])
    expContrasts[[ratGrp]] <- plyr::rbind.fill(lapply(grp0, function(g0) { data.frame(x1 = paste0("Group_", grp1), x0 = paste0("Group_", g0)) }))
  }
  expContrasts <- plyr::rbind.fill(expContrasts)
  tmp <- expContrasts
  w <- which(!tmp$x1 %in% colnames(designMatr))
  if (length(w)) { tmp$x1[w] <- "" }
  w <- which(!tmp$x2 %in% colnames(designMatr))
  if (length(w)) { tmp$x2[w] <- "" }
  tmp <- apply(tmp, 1, function(x) {
    #x <- x[which(x != "")]
    if (length(x) == 2) { x <- paste(x, collapse = " - ") }
    return(x)
  })
  contrCall <- paste0("contrMatr <- makeContrasts(",
                      paste(tmp, collapse = ", "),
                      ", levels = designMatr)")
  # (NB: these could be renamed to something shorter, e.g. makeContrasts(Comp1 = A - B, Comp2 = A - C)
  #cat(contrCall, "\n")
  eval(parse(text = contrCall))
  expContrasts$Name <- colnames(contrMatr)
  #
  colMatch <- data.frame(x = c("up", "down", "n.s.", "n.t."),
                         y = c(1, -1, 0, NA))
  myColors <- setNames(c("limegreen","red3", "black", "grey"),
                       c("up", "down", "n.s.", "n.t."))
  colScale <- scale_colour_manual(name = "colour", values = myColors)
  #
  FDR <- 0.3
  FC <- 1.5
  l10FC <- log10(FC) # For decideTests, to which we feed log10 data and so which needs a log10 FC!
  l2FC <- log2(FC) # For plotting
  statsXprs <- expression({
    intCols <- grep(intRoot[intType], colnames(myData), value = TRUE)
    w <- which(myData[, intCols] == 0, arr.ind = TRUE)
    if (nrow(w)) {
      myData[, intCols][w] <- NA
    }
    logIntCols <- gsub(topattern(intRoot[intType]),
                       logIntRoot[intType],
                       intCols)
    myData[, logIntCols] <- log10(myData[, intCols])
    #
    # First calculate mean per group
    for (grp in Groups) { #grp <- Groups[1] #grp <- Groups[2]
      #cat(" - Mean intensities ", grp, "\n")
      w <- which(samplesMap2$Group == grp)
      e1 <- samplesMap2[w,]
      stopifnot(length(unique(e1$"Sample name")) == length(e1$"Sample name"))
      allSamples1 <- e1$"Sample name"
      col1 <- paste0(logIntRoot[intType], allSamples1)
      w1 <- which(col1 %in% colnames(myData))
      if (length(w1)) { # Calculate average expression
        e1 <- e1[w1,]; col1 <- col1[w1]
        ke1 <- paste0("Mean ", logIntRoot[intType], grp)
        tempVal <- myData[, col1, drop = FALSE]
        row.names(tempVal) <- myData[[namesCol]]
        if (length(w1) == 1) {
          cat(paste0("Warning: Only 1 replicate for ", grp, "!\n"))
          myData[[ke1]] <- tempVal[[1]]
        } else {
          myData[[ke1]] <- parApply(parClust, tempVal, 1, mean, na.rm = TRUE)
        }
      }
    }
    rm(grp)
    #
    # Average lfc
    if (Nested) {
      for (contr in colnames(contrMatr)) { #contr <- colnames(contrMatr)[1] #contr <- colnames(contrMatr)[2]
        tmp <- gsub("Group_", "", expContrasts[match(contr, expContrasts$Name), c("x1", "x0")])
        e <- samplesMap2[which(samplesMap2$Group %in% tmp),]
        uRps <- unique(e$Replicate)
        rps <- lapply(uRps, function(x) {
          x0 <- e$"Sample name"[which((e$Replicate == x)&(e$Reference))]
          x1 <- e$"Sample name"[which((e$Replicate == x)&(!e$Reference))]
          if ((length(x0) == 1)&&(length(x1) == 1)) {
            rs <- data.frame(Rep = x,
                             Ref = x0,
                             Sample = x1)
          } else {
            rs <- data.frame(Rep = x,
                             Ref = NA,
                             Sample = NA)
          }
          return(rs)
        })
        rps <- plyr::rbind.fill(rps)
        rps <- rps[which(!is.na(rps$Ref)),]
        rps <- rps[which(!is.na(rps$Sample)),]
        myData[, paste0(ratRoot, nm, " ", uRps)] <- cbind(lapply(uRps, function(x) { #x <- uRps[1]
          m <- match(x, rps$Rep)
          (myData[[paste0(logIntRoot[intType], rps$Sample[m])]] - myData[[paste0(logIntRoot[intType], rps$Ref[m])]])/log10(2)
        }))
        myData[[paste0(ratRoot, nm)]] <- apply(myData[, paste0(ratRoot, nm, " ", uRps)], 1, mean, na.rm = TRUE)
      }
    } else {
      myData[, paste0(ratRoot, gsub("Group_", "", colnames(contrMatr)))] <- cbind(lapply(1:nrow(expContrasts), function(x) { #x <- 1
        x <- gsub("^Group_", "", expContrasts[x,])
        (myData[[paste0("Mean ", logIntRoot[intType], x[1])]] - myData[[paste0("Mean ", logIntRoot[intType], x[2])]])/log10(2)
      }))
    }
    #View(myData[, c(paste0(logIntRoot[intType], samplesMap2$Sample), paste0(ratRoot, gsub("Group_", "", colnames(contrMatr))))])
    #
    tempVal <- myData[, paste0(logIntRoot[intType], samplesMap2$Sample)]
    fit <- lmFit(tempVal, designMatr) # It's a nested design
    fit <- contrasts.fit(fit, contrMatr)
    fit <- eBayes(fit) # Assumes 1% of genes change!!!
    #View(fit$p.value)
    pVal <- fit$p.value
    kol2 <- sapply(colnames(pVal), function(contr) { #contr <- colnames(pVal)[2]
      tmp <- gsub("Group_", "", expContrasts[match(contr, expContrasts$Name), c("x1", "x0")])
      nm <- paste(tmp, collapse = " - ")
      return(paste0("Mod. P-value - ", nm))
    })
    myData[, kol2] <- pVal
    FDR_thresh <- list()
    for (i in 1:ncol(pVal)) {
      contr <- colnames(contrMatr)[i]
      tmp <- gsub("Group_", "", expContrasts[match(contr, expContrasts$Name), c("x1", "x0")])
      nm <- paste(tmp, collapse = " - ")
      sig <- proteoCraft::FDR(as.data.frame(pVal),
                              pvalue_col = colnames(pVal)[i],
                              fdr = 30,
                              returns = c(TRUE, TRUE))
      #View(sig$`Significance vector`)
      myData[[paste0("Signif. - ", nm)]] <- sig$`Significance vector`
      FDR_thresh[[nm]] <- sig$Thresholds
    }
    decisions <- decideTests(fit, adjust.method = "BH", p.value = FDR, lfc = l10FC) # Remember! The log base of the lfc threshold must be the same as that of the input data!!!
    decisions <- as.data.frame(decisions@.Data)
    nms <- sapply(colnames(decisions), function(contr) { #contr <- colnames(decisions)[2]
      tmp <- gsub("Group_", "", expContrasts[match(contr, expContrasts$Name), c("x1", "x0")])
      nm <- paste(tmp, collapse = " - ")
      return(nm)
    })
    kol2 <- paste0("Decision - ", nms)
    colnames(decisions) <- nms
    myData[, kol2] <- decisions
    venn <- vennCounts(decisions)
    print(venn)
    #
    for (contr in colnames(contrMatr)) { #contr <- colnames(contrMatr)[2]
      tmp <- gsub("Group_", "", expContrasts[match(contr, expContrasts$Name), c("x1", "x0")])
      nm <- paste(tmp, collapse = " - ")
      dc <- paste0("Decision - ", nm)
      myData[[dc]] <- colMatch$x[match(myData[[dc]], colMatch$y)]
    }
    myData$"is Histone" <- c("-", "+")[(myData$`Histone(s)` != "") + 1]
    myData$"is Histone"[which(is.na(myData$"is Histone"))] <- "-"
    myData$Alpha <- c(0.2, 1)[match(myData$"is Histone", c("-", "+"))]
    myData$Size <- c(0.1, 1.5)[match(myData$"is Histone", c("-", "+"))]
    for (contr in colnames(contrMatr)) { #contr <- colnames(contrMatr)[2]
      tmp <- gsub("Group_", "", expContrasts[match(contr, expContrasts$Name), c("x1", "x0")])
      nm <- paste(tmp, collapse = " - ")
      thresh <- FDR_thresh[[nm]]
      ttl <- paste0(namesRoot, " ", nm)
      fc <- paste0(ratRoot, nm)
      pv <- paste0("Mod. P-value - ", nm)
      dc <- paste0("Decision - ", nm)
      kl <- c("is Histone", "Alpha", "Size", "Label")
      myDat <- myData[, c(fc, pv, dc, kl)]
      myDat[[pv]] <- -log10(myDat[[pv]])
      colnames(myDat) <- c("log2 FC", "-log10 mod. P-value", "Decision", kl)
      plot <- ggplot(myDat) +
        geom_point(aes(x = `log2 FC`, y = `-log10 mod. P-value`, colour = `Decision`,
                       shape = `is Histone`, alpha = Alpha, size = Size, text = Label)) +
        theme_bw() + ggtitle(ttl) + colScale +
        geom_vline(xintercept = -l2FC, color = "red3", linetype = "dashed") +
        geom_vline(xintercept = l2FC, color = "limegreen", linetype = "dashed") +
        geom_hline(yintercept = -log10(thresh), color = "purple", linetype = "dashed") +
        scale_alpha_identity(guide = "none") + scale_size_identity(guide = "none") +
        ylab("-log10(moderated P-value") +
        scale_y_continuous(expand = c(0, 0))
      #poplot(plot, 12, 22)
      suppressWarnings({
        ggsave(paste0(dstDir, "/", ttl, ".jpeg"), plot, dpi = 300)
      })
      plotLy <- ggplotly(plot, tooltip = c("x", "y", "text"))
      htmlwidgets::saveWidget(plotLy, paste0(dstDir, "/", ttl, ".html"), selfcontained = TRUE)
      system(paste0("open \"", dstDir, "/", ttl, ".html\""))
    }
  })
  #
  for (iii in 1:2) {
    #iii <- 1
    #iii <- 2
    if (iii == 1) {
      myData <- pep
      namesCol <- "Label"
      namesRoot <- "Peptides"
    }
    if (iii == 2) {
      myData <- Sites
      namesCol <- "Site"
      namesRoot <- "Sites"
    }
    eval(statsXprs)
  }
}

# Specific peptides/sites
for (cmpGrp in compGrps) { #cmpGrp <- compGrps[1]
  grps <- unique(groupsMap2$Group[which(groupsMap2$Comparison_group == cmpGrp)])
  grps2smpls <- setNames(lapply(grps, function(grp) { unique(samplesMap2$"Sample name"[which(samplesMap2$Group == grp)]) }), grps)
  grps2kol <- setNames(lapply(grps2smpls, function(x) { paste0(intRoot[intType], x) }), grps)
  stopifnot(sum(!unlist(grps2kol) %in% colnames(pep)) == 0,
            sum(!unlist(grps2kol) %in% colnames(Sites)) == 0)
  for (grp1 in grps) { #grp1 <- grps[1]
    k1 <- grps2kol[[grp1]]
    tst <- groupsMap2$Reference[match(grp1, groupsMap2$Group)]
    # This will mean comparing to reference
    grp0 <- groupsMap2$Group[which((groupsMap2$Comparison_group == cmpGrp)&(groupsMap2$Reference != tst))]
    k0 <- unlist(grps2kol[grp0])
    pep[[paste0("Specific - ", grp1)]] <- vapply(1:nrow(pep), function(i) {
      x1 <- is.all.good(as.numeric(pep[i, k1]))
      x0 <- is.all.good(as.numeric(pep[i, k0]))
      res <- c("", "+")[(((!length(x0))&(length(x1) >= 3))|((length(x0) == 1)&(length(x1) >= 4))) + 1]
      return(res)
    }, "a")
    Sites[[paste0("Specific - ", grp1)]] <- vapply(1:nrow(Sites), function(i) {
      x1 <- is.all.good(as.numeric(Sites[i, k1]))
      x0 <- is.all.good(as.numeric(Sites[i, k0]))
      res <- c("", "+")[(((!length(x0))&(length(x1) >= 3))|((length(x0) == 1)&(length(x1) >= 4))) + 1]
      return(res)
    }, "a")
  }
}

# Write tables
# - Peptides
samplesMap3 <- samplesMap[order(samplesMap$Comparison_group,
                            samplesMap$Group,
                            samplesMap$Replicate),]
allSamples3 <- unique(samplesMap3$"Sample name")
Groups3 <- unique(samplesMap3$Group)
kol <- c(paste0(intRoot[intType], allSamples3),
         paste0("Specific - ", Groups3))
kol <- kol[which(kol %in% colnames(pep))]
tmp <- pep[which(pep$`Histone(s)` != ""),]
tmp$Modifications <- gsub(",", ";", tmp$`Modified sequence`)
tmp <- tmp[, c(colnames(tmp)[which(!colnames(tmp) %in% kol)], kol)]
tmp$Label <- gsub("\n", " ", tmp$Label)
fl <- paste0(dstDir, "/Histone peptides.csv")
data.table::fwrite(tmp, fl, quote = FALSE, na = "NA", sep = ",")
#openxlsx2::xl_open(fl)
# - Sites
tmp <- Sites[which(Sites$`Histone(s)` != ""),]
tmp <- tmp[, c(colnames(tmp)[which(!colnames(tmp) %in% kol)], kol)]
tmp$Label <- gsub("\n", " ", tmp$Label)
fl <- paste0(dstDir, "/Histone sites.csv")
data.table::fwrite(tmp, fl, quote = FALSE, sep = ",", na = "NA")
#openxlsx2::xl_open(fl)
#
saveImgFun(backupFl)

#### Heatmaps with clustering
if (length(Exp) > 1) {
  pkgs <- unique(c(pkgs, "ggdendro", "gridExtra", "ggpubr", "colorspace", "ggnewscale", "factoextra", "NbClust"))
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      pak::pkg_install(pkg, ask = FALSE)
    }
  }
  for (pkg in pkgs) {
    library(pkg, character.only = TRUE)
  }
  dir <- paste0(dstDir, "/Clustering")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  xMap <- samplesMap
  #
  ImputeKlust <- TRUE
  MaxHClust <- min(c(floor(nrow(pep)/2), 100)) # We want at most 20 clusters
  MaxVClust <- nrow(groupsMap2)
  VClustScl <- setNames(1:MaxVClust, paste0("Cluster", 1:MaxVClust))
  HClustScl <- setNames(1:MaxHClust, paste0("Cluster", 1:MaxHClust))
  VClustScl <- setNames(rainbow(MaxVClust), VClustScl)
  HClustScl <- setNames(rainbow(MaxHClust), HClustScl)
  # Different options for which proteins to use for vertical clustering (samples level)
  KlustMeth <- 2
  KlustRoot <- paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ") - ")
  normTypes <- c("Norm. by row", "Z-scored")
  kol <- paste0(intRoot[intType], allSamples)
  plotLeatMaps <- list()
  prot.list.Cond <- TRUE
  prot.list <- histDB$`Protein ID`
  KlustKols <- c()
  #
  Filter <- list("Global" = 1:nrow(pep),
                 "Histones only" = which(pep$`Histone(s)` != ""))
  for (normType in normTypes) { #normType <- normTypes[1]
    normTypeInsrt <- c("", paste0(" (", normType, ")"))[match(normType, normTypes)]
    plotLeatMaps[[normType]] <- list()
    for (flt in names(Filter)) {
      fltInsrt <- c("", " Histones only")[match(flt, names(Filter))]
      temp <- magrittr::set_rownames(magrittr::set_colnames(pep[Filter[[flt]], kol], allSamples), pep$Label[Filter[[flt]]])
      temp <- log10(temp)
      Gr <- xMap$Comparison_group
      Gr <- setNames(match(Gr, unique(Gr)), Gr)
      tst <- apply(temp, 1, function(x) { length(is.all.good(x)) })
      filt <- rep(FALSE, nrow(temp))
      temp <- temp[which(tst > 0),]
      if (ImputeKlust) {
        temp2 <- Data_Impute2(temp, Gr)
        imputed <- temp2$Positions_Imputed
        temp <- temp2$Imputed_data
        rownames(imputed) <- rownames(temp)
        colnames(imputed) <- colnames(temp)
      }
      w <- which(apply(temp, 1, function(x) { length(is.all.good(x)) }) == length(kol))
      filt[which(tst > 0)[w]] <- TRUE
      temp <- temp[w,]
      imputed <- imputed[w,]
      rwMns <- rowMeans(temp)
      if (normType == "Norm. by row") {
        temp <- sweep(temp, 1, rwMns, "-")
      }
      if (normType == "Z-scored") {
        SDs <- apply(temp, 1, function(x) { sd(is.all.good(x)) })
        temp <- sweep(sweep(temp, 1, rwMns, "-"), 1, SDs, "/")
      }
      temp2 <- temp3 <- as.matrix(temp) + runif(length(temp), min = 0, max = 1e-10) # Small error added to avoid duplicate rows where this breaks
      # Data is now normalized and either imputed or filtered, so ready for clustering
      # 1/ At samples level
      # We perform hierarchical clustering in all cases, because we want to see the dendrogram.
      # But the clusters displayed using colours may be generated using either k-means or hierarchical clustering (default).
      temp2 <- t(temp2)
      tst2 <- apply(temp2, 2, function(x) {
        #x <- is.all.good(x) # it's been imputed, there are no missing values
        sd(x)/mean(x)
      })
      temp2 <- temp2[, order(tst2, decreasing = TRUE)]
      vcluster <- hclust(dist(temp2))
      vdendro <- as.dendrogram(vcluster)
      # Estimate ideal number of clusters... but ignore it!
      # In fact we know the number of sample groups, so would like to see 1 cluster per group.
      # How well groups and clusters overlap would tell us how well the clustering works, i.e. how different clusters are.
      #
      # Here we use kmeans, but findings apply to any method
      #cat("Estimating optimal number of samples-level clusters...\n")
      vnm <- paste0("Samples-level clusters number analysis", fltInsrt)
      NVClust <- NGr <- length(unique(Gr))
      tst <- cluster::clusGap(temp2, kmeans, min(c(nrow(temp2), NGr+1)))
      tst2 <- as.data.frame(tst$Tab)
      yScl <- max(tst2$gap)
      tst2 <- sapply(1:NGr, function(x) { tst2$gap[x] >= tst2$gap[x+1] - tst2$SE.sim[x+1] })
      tst2 <- which(tst2)
      # I like to do one more, often these methods feel too conservative
      #if (length(tst2)) { NVClust <- tst2[1]+1 } # Not used for now: we use NGr instead
      vplot <- fviz_gap_stat(tst)
      tstLy <- capture.output(print(vplot$layers))
      g1 <- grep("geom_vline", tstLy)
      g2 <- grep("\\[\\[[0-9]+\\]\\]", tstLy)
      g2 <- as.numeric(gsub("\\[|\\]", "", tstLy[max(g2[which(g2 < g1)])]))
      vplot$layers[[g2]] <- NULL
      vplot <- vplot +
        geom_vline(xintercept = tst2[1]+1, colour = "red", linetype = "dashed") +
        geom_vline(xintercept = NGr, colour = "orange", linetype = "dashed") +
        geom_text(label = "Optimal number of sample clusters", x = tst2[1]+1-0.2, y = yScl*0.9,
                  colour = "red", angle = 90, hjust = 1) +
        geom_text(label = "Number of sample groups", x = NGr+0.2, y = yScl*0.9,
                  colour = "orange", angle = 90, hjust = 1) +
        theme_bw() + ggtitle(vnm)
      #poplot(vplot)
      suppressWarnings({
        ggsave(paste0(dir, "/", vnm, normTypeInsrt, ".jpeg"), vplot, dpi = 150)
        ggsave(paste0(dir, "/", vnm, normTypeInsrt, ".pdf"), vplot, dpi = 150)
      })
      NVClust <- max(c(NGr, 2))
      # 2/ At protein groups level
      # As above, we always draw a dendrogram, but colours will be defined by the clustering approach.
      hcluster <- hclust(dist(temp3))
      hdendro <- as.dendrogram(hcluster)
      hnm <- paste0("Protein groups-level clusters number analysis", fltInsrt)
      NHClust <- MaxHClust
      Straps <- 10
      # Here we really want to optimize the number of clusters
      # Apply the same method for optimization for any clustering method
      # Number of cluster should not depend on method
      clusterExport(parClust, list("temp3", "Straps"), envir = environment())
      tst <- parSapply(parClust, 2:MaxHClust, function(kl) {
        kmeans(temp3, kl, nstart = Straps)$tot.withinss
      })/kmeans(temp3, 1, nstart = 1)$tot.withinss
      yScl2 <- max(tst)
      tst2 <- data.frame("Number of clusters k" = 2:MaxHClust,
                         "[tot WSS (k)]/[tot WSS (1)]" = tst,
                         check.names = FALSE)
      # For the elbow detection method, we need to normalize to a 1x1 graph which we can rotate by 45 degrees
      tst2$X1 <- tst2$`Number of clusters k`/MaxHClust # divide by max theoretical number of clusters 
      tst2$Y1 <- tst2$"[tot WSS (k)]/[tot WSS (1)]" # Here no need to normalize, this is a ratio
      Angle <- -pi/4
      meanX <- mean(tst2$X1)
      meanY <- mean(tst2$Y1)
      tst2$X2 <- tst2$X1 - meanX
      tst2$Y2 <- tst2$Y1 - meanY
      tst2$X2 <- tst2$X2*cos(Angle)+tst2$Y2*sin(Angle)
      tst2$Y2 <- -tst2$X2*sin(Angle)+tst2$Y2*cos(Angle)
      tst2$X2 <- tst2$X2 - mean(tst2$X2) + meanX
      tst2$Y2 <- tst2$Y2 - mean(tst2$Y2) + meanY
      w <- rev(which(tst2$Y2 == min(tst2$Y2)))[1] # Again, prefer more rather than fewer clusters
      NHClust <- tst2$`Number of clusters k`[w]
      tst2$Size <- 1
      tst2$Size[w] <- 2
      xMin <- min(c(tst2$X1, tst2$X2))
      xMax <- max(c(tst2$X1, tst2$X2))
      xScl <- xMax-xMin
      yMin <- min(c(tst2$Y1, tst2$Y2))
      yMax <- max(c(tst2$Y1, tst2$Y2))
      yScl <- yMax-yMin
      hplot <- ggplot(tst2) +
        geom_segment(x = tst2$X2[w], y = tst2$Y2[w],
                     xend = tst2$X1[w], yend = tst2$Y1[w], color = "grey", linetype = "dotted") +
        geom_point(aes(x = X1, y = Y1, size = Size), color = "blue") +
        geom_point(aes(x = X2, y = Y2, size = Size), color = "red") +
        geom_text(x = xScl*0.98+xMin, y = yMin+yScl*0.98, color = "blue", label = "ratio of tot. WSS", hjust = 1) +
        geom_text(x = xScl*0.98+xMin, y = yMin+yScl*0.962, color = "red", label = "ratio of tot. WSS, -pi/4 rotation", hjust = 1) +
        geom_hline(yintercept = tst2$Y2[w], color = "red", linetype = "dashed") +
        geom_vline(xintercept = tst2$X1[w], color = "deepskyblue", linetype = "dashed") +
        geom_text(x = xMin+0.01*xScl, y = tst2$Y2[w]+0.02*yScl, color = "red", label = "Elbow", hjust = 0) +
        geom_text(x = tst2$X1[w]-0.02*xScl, y = yScl*0.9, angle = 90, color = "deepskyblue",
                  label = paste0("Optimal = ", tst2$`Number of clusters k`[w], " clusters"), hjust = 1) +
        ggtitle(hnm) + scale_size_identity() + theme_bw() +
        theme(legend.position = "none") + ylab("Normalised total Within-clusters vs Total Sum of Squares")
      #poplot(hplot)
      suppressWarnings({
        ggsave(paste0(dir, "/", hnm, normTypeInsrt, ".jpeg"), hplot, dpi = 150)
        ggsave(paste0(dir, "/", hnm, normTypeInsrt, ".pdf"), hplot, dpi = 150)
      })
      # Apply cutoffs
      if (KlustMeth == 1) {
        VClusters <- kmeans(t(temp3), NVClust, 100)$cluster
        tempClust <- kmeans(temp3, NHClust, 100)$cluster
      }
      if (KlustMeth == 2) {
        VClusters <- cutree(vcluster, NVClust)
        tempClust <- cutree(hcluster, NHClust)
      }
      KlKol <- KlustRoot
      KlustKols <- unique(c(KlustKols, KlKol))
      pep[[KlKol]] <- tempClust[match(pep$Label, names(tempClust))]
      #
      Width <- nrow(temp)
      Height <- ncol(temp)
      #
      whImps <- which(imputed, arr.ind = TRUE)
      temp[whImps] <- NA
      #
      # Get dendrograms
      vddata <- dendro_data(vdendro)
      hddata <- dendro_data(hdendro)
      # vdendro.plot <- ggdendrogram(data = vdendro) +
      #   theme(axis.text.y = element_text(size = 0.1), plot.margin = margin(0, 0, 0, 0, "cm"))
      # poplot(vdendro.plot)
      # Modify dendrograms
      # - Vertical
      vlabs <- ggdendro::label(vddata)
      vlabs$Cluster <- as.factor(VClusters[match(vlabs$label, allSamples)])
      vSeg <- vddata$segments
      # - Horizontal
      hlabs <- ggdendro::label(hddata)
      hlabs$Cluster <- as.factor(tempClust[match(hlabs$label, names(tempClust))])
      hSeg <- hddata$segments
      # Rotate samples-level dendrogram: x -> y, y -> x (for ease of calculation, gets inverted later)
      vlabs$y <- vlabs$x
      vlabs$x <- -0.3
      x <- vSeg$x
      y <- vSeg$y
      vSeg$y <- x
      vSeg$x <- y
      xend <- vSeg$xend
      yend <- vSeg$yend
      vSeg$yend <- xend
      vSeg$xend <- yend
      # Adjust width/height
      # - Vertical dendrogram
      vnch <- Width*0.07*max(nchar(vlabs$label))/12
      xmn <- min(c(vSeg$x, vSeg$xend))
      xmx <- max(c(vSeg$x, vSeg$xend))
      vSeg$x <- -(vSeg$x - xmn)*Width*0.07/xmx - 1.5 - vnch
      vSeg$xend <- -(vSeg$xend - xmn)*Width*0.07/xmx - 1.5 - vnch
      # - Horizontal dendrogram
      hnch <- Height*0.07*max(nchar(hlabs$label))/800
      ymn <- min(c(hSeg$y, hSeg$yend))
      ymx <- max(c(hSeg$y, hSeg$yend))
      hSeg$y <- Height + hnch + 1.5 + (hSeg$y - ymn)*Height*0.07/ymx
      hSeg$yend <- Height + hnch + 1.5 + (hSeg$yend - ymn)*Height*0.07/ymx
      # Order labels by order of appearance
      hlabs <- hlabs[order(hlabs$x, decreasing = FALSE),]
      vlabs <- vlabs[order(vlabs$y, decreasing = FALSE),]
      # Re-order our matrix based on extracted dendrogram labels
      temp <- temp[, match(vlabs$label, colnames(temp))]
      temp <- temp[match(hlabs$label, rownames(temp)),]
      imputed <- imputed[, match(vlabs$label, colnames(imputed))]
      imputed <- imputed[match(hlabs$label, rownames(imputed)),]
      # Re-introduce missing values
      MaxChar <- 13
      hlabs$label2 <- hlabs$label
      w <- which(nchar(hlabs$label2) > MaxChar)
      hlabs$label2[w] <- paste0(substr(hlabs$label2[w], 1, MaxChar-3), "...")
      # Create heatmap
      temp$Rowname <- row.names(temp)
      temp2 <- magrittr::set_colnames(reshape2::melt(temp, id.vars = "Rowname"), c("Label", "Sample", "value"))
      temp2$Label <- as.character(temp2$Label)
      temp2$Sample <- as.character(temp2$Sample)
      temp2$Xmax <- match(temp2$Label, hlabs$label) # Explicitly should be the case now!
      temp2$label2 <- hlabs$label2[temp2$Xmax]
      temp2$Xmin <- temp2$Xmax-1
      temp2$Ymax <- vlabs$y[match(temp2$Sample, vlabs$label)]
      temp2$Ymin <- temp2$Ymax-1
      w1 <- which(temp2$Colour == "green")
      w2 <- which((temp2$Ymin == max(temp2$Ymin))&(temp2$Colour == "green"))
      # Color and fill scales
      wV <- unique(ceiling(c(1:NVClust)*MaxVClust/NVClust))
      wH <- unique(ceiling(c(1:NHClust)*MaxHClust/NHClust))
      vClScl <- setNames(VClustScl[wV], seq_along(wV))
      hClScl <- setNames(HClustScl[wH], seq_along(wH))
      VcolScale <- scale_color_manual(name = "Samples cluster", values = vClScl)
      VfillScale <- scale_fill_manual(name = "Samples cluster", values = vClScl)
      HcolScale <- scale_color_manual(name = "Protein groups cluster", values = hClScl)
      HfillScale <- scale_fill_manual(name = "Protein groups cluster", values = hClScl)
      # Create heatmap plot
      Xlim <- c(NA, Width)
      Ylim <- c(-10, max(c(max(hSeg$y) + Height*0.6), 20))
      #
      prot.list.marks <- FALSE
      if (prot.list.Cond) {
        w0 <- which(temp2$Ymin == 0)
        temp2$`Histone(s)` <- pep$`Histone(s)`[match(temp2$Label, pep$Label)]
        g <- which(temp2$`Histone(s)` != "")
        #tst <- unlist(strsplit(temp2$"Leading protein IDs"[w0], ";"))[1]
        #g <- grsep(tst, x = temp2$"Leading protein IDs"[w0])
        if (length(g)) {
          prot.list.marks <- TRUE
          Ylim[1] <- -20
          temp2c <- temp2[w0[g], , drop = FALSE]
        }
      }
      # Main data
      temp2a <- temp2[, c("Xmin", "Ymin", "value", "Label", "Sample")]
      # Colour scale
      temp2b <- data.frame(Xmin = 0:round(Width*0.1),
                           Ymin = Ylim[2]*0.8)
      Mn <- min(temp2a$value, na.rm = TRUE)
      Mx <- max(temp2a$value, na.rm = TRUE)
      temp2b$value <- Mn + temp2b$Xmin*(Mx-Mn)/max(temp2b$Xmin)
      temp2b$Xmin <- temp2b$Xmin-Width*0.15
      temp2b$Label <- temp2b$Sample <- NA
      w2a <- 1:nrow(temp2a)
      w2b <- 1:nrow(temp2b) + max(w2a)
      temp2a <- rbind(temp2a, temp2b)
      # Samples-level: how well do clusters fit expectations
      SamplesClust <- vlabs[, c("label", "Cluster")]
      SamplesClust$Group <- as.factor(samplesMap$Group[match(SamplesClust$label, samplesMap$"Sample name")])
      SamplesClust$Cluster <- as.factor(paste0("Cluster", as.character(SamplesClust$Cluster)))
      tstSmplClust <- table(SamplesClust$Group, SamplesClust$Cluster)
      ClustChiSqTst <- suppressWarnings(chisq.test(tstSmplClust))
      #
      xCntr <- Width*0.6
      yScale <- Height*2
      nm <- paste0("Clust. heatmap", normTypeInsrt, fltInsrt)
      GrLab <- data.frame(label = c(nm,
                                    "Sample",
                                    paste0("Chi-squared contingency test P-value: ", round(ClustChiSqTst$p.value, 5)),
                                    "Protein group",
                                    "Min",
                                    "Max"),
                          x = c(xCntr,
                                min(vSeg$x)-Width*0.11,
                                min(vSeg$x)-Width*0.095,
                                xCntr,
                                min(temp2b$Xmin),
                                max(temp2b$Xmin)),
                          y = c(max(c(hSeg$y, hSeg$yend)) + yScale*0.15,
                                Height*0.5,
                                Height*0.5,
                                max(c(hSeg$y, hSeg$yend)) + yScale*0.1,
                                temp2b$Ymin[1]-1,
                                temp2b$Ymin[1]-1),
                          angle = c(0, 90, 90, 0, 0, 0),
                          size = c(5, 4, 3, 4, 3, 3),
                          fontface = c("bold", "plain", "italic", "plain", "plain", "plain"))
      Xlim[1] <- min(c(vSeg$x, vSeg$xend))-Width*0.15
      Ylim[1] <- min(c(Ylim[1], min(temp2a$Ymin)))
      Xlim[2] <- max(c(Xlim[2], max(temp2a$Xmin)+1))
      Ylim[2] <- max(c(Ylim[2], max(temp2a$Ymin)+1))
      xCntr <- mean(Xlim)
      yCntr <- mean(Ylim)
      yScale <- Ylim[2]-Ylim[1]
      #
      # Create graph
      heatmap.plot <- ggplot(temp2a[w2a,]) +
        geom_rect(aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value, text = Label)) +
        geom_rect(data = temp2a[w2b,], aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
              panel.background = element_rect(fill = "transparent", color = NA),
              plot.margin = margin(0, 0, 0, 0, "cm")) +
        #scale_fill_gradient2(low = "darkblue", mid = "lightgrey", high = "darkred") +
        scale_fill_viridis(option = "D") +
        xlab(NULL) + ylab(NULL) + theme(legend.position = "none") +
        xlim(Xlim[1], Xlim[2]) + ylim(Ylim[1], Ylim[2])
      #poplot(heatmap.plot, 12, 20)
      # heatmap.plot <- heatmap.plot + 
      #   geom_text(data = temp2a[which(temp2a$Xmin == min(temp2$Xmin)),],
      #             aes(x = Xmin, y = Ymin, label = Sample))
      # Title, axis labels, colour scale annotations
      #
      # Add labels and dendrograms:
      # - Title, scale, chi-squared test
      heatmap.plot <- heatmap.plot +
        new_scale("size") + scale_size_identity() +
        geom_text(data = GrLab, aes(label = label, x = x, y = y, angle = angle, size = size, fontface = fontface),
                  color = "black", hjust = 0.5)
      #poplot(heatmap.plot, 12, 20)
      # - Cluster boxes
      heatmap.plot <- heatmap.plot + new_scale_color() + HcolScale
      if (KlustMeth == 2) {
        Clutst <- aggregate(hlabs$x, list(hlabs$Cluster), function(x) { min(x)-1 })
        colnames(Clutst)[1] <- "Cluster"
        Clutst$xend <- aggregate(hlabs$x, list(hlabs$Cluster), max)$x
        Clutst$mid <- (Clutst$xend+Clutst$x)/2
        heatmap.plot <- heatmap.plot +
          geom_rect(data = Clutst, aes(xmin = x, xmax = xend, colour = Cluster),
                    ymin = 0, ymax = Height, fill = NA) +
          geom_text(data = Clutst, aes(x = mid, label = Cluster, colour = Cluster),
                    y = Height/2, hjust = 0.5, vjust = 0.5, cex = 5)
        #poplot(heatmap.plot2, 12, 20)
      }
      #poplot(heatmap.plot, 12, 20)
      # - Horizontal dendrogram and protein  names
      heatmap.plot <- heatmap.plot +
        geom_segment(data = hSeg, linewidth = 0.5,
                     aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_text(data = hlabs, aes(x = x-0.5, label = label2, colour = Cluster),
                  y = Height+0.05, angle = 90, cex = 0.5, hjust = 0, vjust = 0)
      # - Vertical dendrogram and sample names
      heatmap.plot <- heatmap.plot +
        geom_segment(data = vSeg, linewidth = 0.5,
                     aes(x = x, y = y - 0.5, xend = xend, yend = yend - 0.5)) +
        new_scale_color() + VcolScale +
        geom_text(data = vlabs, aes(y = y - 0.5, label = label, colour = Cluster),
                  x = -0.5, hjust = 1, vjust = 0.5, cex = 2.5)
      # - Proteins of interest
      if (prot.list.marks) {
        heatmap.plot <- heatmap.plot +
          geom_point(data = temp2c, aes(x = Xmin+0.5), y = -0.5, colour = "red", fill = "red", shape = 17) +
          geom_text(data = temp2c, aes(x = Xmin+0.5, label = label2),
                    y = -1, colour = "red", angle = -60, hjust = 0, cex = 2)
      }
      #poplot(heatmap.plot, 12, 20)
      suppressWarnings({
        ggsave(paste0(dir, "/", nm, ".jpeg"), heatmap.plot)
        ggsave(paste0(dir, "/", nm, ".pdf"), heatmap.plot)
      })
      #
      # Plotly version
      tempLy <- temp2a[w2a,]
      tempLy$Sample <- factor(temp2$Sample, levels = allSamples)
      ##tempLy$Ymin <- tempLy$Ymin+0.5
      plotleatmap <- plot_ly(data = tempLy, x = ~Xmin, y = ~Ymin, z = ~value, type = "heatmap", hovertext = tempLy$Label)
      # I cannot find a way to remove tick marks!!!!!
      plLyV <- vSeg
      ##plLyV$x <- -Width*0.2-plLyV$x 
      ##plLyV$xend <- -Width*0.2-plLyV$xend
      plLyV$y <- plLyV$y-1
      plLyV$yend <- plLyV$yend-1 
      plotleatmap <- add_segments(plotleatmap, x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                                  data = plLyV, inherit = FALSE, color = I("black"),
                                  showlegend = FALSE)
      plLyH <- hSeg
      ##plLyH$x <- plLyH$x-0.5
      ##plLyH$xend <- plLyH$xend-0.5
      plLyH$y <- plLyH$y - hnch - 1.5
      plLyH$yend <- plLyH$yend - hnch - 1.5
      plotleatmap <- add_segments(plotleatmap, x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                                  data = plLyH, inherit = FALSE, color = I("black"),
                                  showlegend = FALSE)
      if (KlustMeth == 2) { # Cluster shapes do not appear to work
        plotleatmap <- add_segments(plotleatmap, x = ~ x, xend = ~ xend,
                                    y = I(-0.2*as.numeric(Clutst$Cluster))-1,
                                    yend = I(-0.2*as.numeric(Clutst$Cluster))-1, data = Clutst,
                                    inherit = FALSE, color = ~ Cluster, showlegend = FALSE)
      }
      vlabs2 <- vlabs
      vlabs2$x <- -vnch/2
      vlabs2$y <- vlabs2$y-1
      plotleatmap <- add_trace(plotleatmap, data = vlabs2, y = ~y, x = ~x, text = ~label,
                               color = I("black"), inherit = FALSE, type = "scatter",
                               mode = "text", showlegend = FALSE)
      plotleatmap <- layout(plotleatmap, title = nm,
                            xaxis = list(title = "Protein groups",
                                         tickmode = "array", tickvals = NULL, showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
                            yaxis = list(title = "Samples",
                                         tickmode = "array", tickvals = NULL, showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE))
      if ((prot.list.marks)&&(flt == "Global")) {
        temp2c$X <- (temp2c$Xmin+temp2c$Xmax)/2
        plotleatmap <- add_trace(plotleatmap, data = temp2c, y = I(-0.501), x = ~X,
                                 text = ~Label, color = ~Label, inherit = FALSE,
                                 type = "scatter", mode = "markers", showlegend = FALSE)
      }
      htmlwidgets::saveWidget(plotleatmap, paste0(dir, "/", nm, ".html"))
      system(paste0("open \"", dir, "/", nm, ".html\""))
      plotLeatMaps[[normType]][[flt]] <- plotleatmap
    }
  }
}
invisible(parLapply(parClust, 1:N.clust, function(x) { rm(list = ls());gc() }))
# Save modifications table
tmp <- Modifs
tmp$Site_long <- NULL
tmp$Site <- NULL
tmp$origMark <- NULL
for (k in colnames(tmp)) {
  if (typeof(tmp[[k]]) == "list") { tmp[[k]] <- vapply(tmp[[k]], paste, "", collapse = ";") }
}
write.csv(tmp, paste0(dstDir, "/Mods table.csv"), row.names = FALSE)

if (dirname(ScriptPath) != dstDir) { file.copy(ScriptPath, dstDir, overwrite = TRUE) }
saveImgFun(backupFl)
#openwd(dstDir)

###################################################################################
#                                 Et voilà, done!                                 #
###################################################################################
