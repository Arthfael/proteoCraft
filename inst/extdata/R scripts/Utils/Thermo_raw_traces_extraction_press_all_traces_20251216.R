#
options(stringsAsFactors = FALSE)
getInt <- FALSE
if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
if(!require("data.table", quietly = TRUE)) { install.packages("data.table") }
if(!require("XML", quietly = TRUE)) { install.packages("XML") }
if(!require("plyr", quietly = TRUE)) { install.packages("plyr") }
if(!require("svDialogs", quietly = TRUE)) { install.packages("svDialogs") }
if(!require("parallel", quietly = TRUE)) { install.packages("parallel") }
if(!require("snow", quietly = TRUE)) { install.packages("snow") }
if(!require("ggplot2", quietly = TRUE)) { install.packages("ggplot2") }
require(svDialogs)
require(parallel)
require(ggplot2)

rm(list = ls()[!ls() %in% c("fls")])

# Create cluster
if (exists("parClust")) { try(stopCluster(parClust), silent = TRUE) } # A fresh start if re-running
N.clust <- detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(stopCluster(parClust), silent = TRUE)
  parClust <- makeCluster(N.clust, type = "SOCK")
}

# Some useful functions
cleanRawNm <- function(rawFileNm) { gsub("\\.raw$", "", rawFileNm, ignore.case = TRUE) }
environment(cleanRawNm) <- .GlobalEnv
clusterExport(parClust, "cleanRawNm", envir = environment())
listMelt <- function (List, Names = NULL, ColNames = c("value", "L1")) {
  List <- setNames(List, 1:length(List))
  List <- stack(List)
  if (!is.null(Names)) {
    colnames(List) <- c(ColNames[1], "ind")
    List$ind <- as.integer(List$ind)
    List[[ColNames[2]]] <- Names[List$ind]
    List$ind <- NULL
  }
  else {
    colnames(List) <- ColNames
    List[[ColNames[2]]] <- as.integer(List[[ColNames[2]]])
  }
  return(List)
}

# Install/load packages
rawrrVers <- c("github",
               "bioc",
               "bioc_1.11.14")
tst <- try({
  if (!require(rawrr, quietly = TRUE)) {
    rawrrVers <- dlg_list(rawrrVers, rawrrVers[3], title = "Select which rawrr version should be installed")$res
    if (rawrrVers == "github") {
      if(!require(devtools)) { install.packages("devtools") }
      devtools::install_github("cpanse/rawrr")
      #install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawrr_0.2.1.tar.gz', repo = NULL) # Old address, now on github
    }
    if (rawrrVers == "bioc") {
      if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
      BiocManager::install("rawrr", version = "1.11")
      BiocManager::install("rawrr")
    }
    if (rawrrVers == "bioc_1.11.14") {
      url <- "https://bioconductor.org/packages/3.19/bioc/src/contrib/Archive/rawrr/rawrr_1.11.14.tar.gz"
      require(curl)
      dstfl <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/rawrr_1.11.14.tar.gz")
      curl_download(url, dstfl)
      install.packages(dstfl)
    }
    # 
    require(rawrr)
    yesRawFileReaderLicenseIsAccepted <- function () {
      licenseFile <- file.path(system.file(package = "rawrr"), 
                               "rawrrassembly", "RawFileReaderLicense.txt")
      stopifnot(file.exists(licenseFile))
      eulaFile <- file.path(rawrrAssemblyPath(), "eula.txt")
      msg <- c("# By changing the setting below to TRUE you are accepting ", 
               "the Thermo License agreement.")
      if (!file.exists(eulaFile)) {
        file.show(licenseFile)
        response <- "y"
        if (tolower(response) == "y") {
          if (isFALSE(dir.exists(dirname(eulaFile)))) {
            dir.create(dirname(eulaFile), recursive = TRUE)
          }
          fileConn <- file(eulaFile)
          writeLines(paste(msg, paste0("# ", date()), "eula=true", 
                           sep = "\n"), fileConn)
          close(fileConn)
          return(TRUE %in% grepl("eula=true", tolower(readLines(eulaFile))))
        }
      } else { return(TRUE %in% grepl("eula=true", tolower(readLines(eulaFile)))) }
      msg <- "Yes, we accept Thermo's License agreement, get on with it!"
      cat(msg, "\n")
    }
    installRawFileReaderDLLsNoAcpt <- function (sourceUrl = .thermofisherlsmsUrl(), ...) {
      rawfileReaderDLLsPath <- rawrrAssemblyPath()
      if (isTRUE(dir.exists(rawfileReaderDLLsPath))) {
        msg <- sprintf("removing files in directory '%s'", rawfileReaderDLLsPath)
        message(msg)
        file.remove(file.path(rawrrAssemblyPath(), list.files(rawrrAssemblyPath())))
      }
      if (isFALSE(dir.exists(rawfileReaderDLLsPath))) {
        dir.create(rawfileReaderDLLsPath, recursive = TRUE)
      }
      #
      stopifnot(yesRawFileReaderLicenseIsAccepted())
      .rawfileReaderDLLs <- getAnywhere(.rawfileReaderDLLs)
      .rawfileReaderDLLs <- .rawfileReaderDLLs$objs[[1]]
      rv <- vapply(.rawfileReaderDLLs(), function(dll) {
        destfile <- file.path(rawfileReaderDLLsPath, dll)
        download.file(file.path(sourceUrl, dll), destfile = destfile, 
                      mode = "wb", ...)
      }, 0)
      rv
    }
    installRawFileReaderDLLsNoAcpt(sourceUrl = .thermofisherlsmsUrl())
    #rawrr::installRawFileReaderDLLs(sourceUrl = .thermofisherlsmsUrl())
    #require(rawrr); getAnywhere(".isRawFileReaderLicenseAccepted")
    rawrr::installRawrrExe()
  }
  require(rawrr)
}, silent = TRUE)
getInt <- !("try-error" %in% class(tst))
#
if (getInt) {
  ticFun <- function(fl) { #fl <- fls[1]
    x <- try(readChromatogram(fl, type = "tic"), silent = TRUE)
    if (!"try-error" %in% class(x)) {
      res <- list(Outcome = TRUE,
                  Output = data.frame("Raw file" = fl,
                                      "Raw file name" = cleanRawNm(basename(fl)),
                                      "Retention time" = as.numeric(x$times),
                                      "Intensity" = as.numeric(x$intensities),
                                      check.names = FALSE))
    } else { res <- list(Outcome = FALSE,
                         Error = x) }
    return(res)
  }
  bpcFun <- function(fl) { #fl <- fls[1]
    x <- try(readChromatogram(fl, type = "bpc"), silent = TRUE)
    if (!"try-error" %in% class(x)) {
      res <- list(Outcome = TRUE,
                  Output = data.frame("Raw file" = fl,
                                      "Raw file name" = cleanRawNm(basename(fl)),
                                      "Retention time" = as.numeric(x$times),
                                      "Intensity" = as.numeric(x$intensities),
                                      check.names = FALSE))
    } else { res <- list(Outcome = FALSE,
                         Error = x) }
    return(res)
  }
  xicFun <- function(fl, mass, tol, filter) { #fl <- fls[1]
    x <- try(readChromatogram(fl, type = "xic", mass, tol, filter), silent = TRUE)
    if (!"try-error" %in% class(x)) {
      res <- list(Outcome = TRUE,
                  Output = data.frame("Raw file" = fl,
                                      "Raw file name" = cleanRawNm(basename(fl)),
                                      "M/Z" = mass,
                                      "Tolerance" = tol,
                                      "Retention time" = as.numeric(x[[1]]$times),
                                      "Intensity" = as.numeric(x[[1]]$intensities),
                                      check.names = FALSE))
    } else { res <- list(Outcome = FALSE,
                         Error = x) }
    return(res)
  }
  environment(ticFun) <- .GlobalEnv
  environment(bpcFun) <- .GlobalEnv
  environment(xicFun) <- .GlobalEnv
}
#require(RaMS)
dflt <- "D:/groups_temp"
if (!dir.exists(dflt)) { dflt <- "D:/Data" }
if (!dir.exists(dflt)) { dflt <- "B:/group/lsfgrp/Mass_Spec/Acquired data_v2" }
if (!dir.exists(dflt)) { dflt <- "B:/archive/lsfgrp/MS/Acquired data" }

# Raw files
if (exists("fls")) {
  fls <- fls[which(file.exists(fls))]
} else { fls <- c() }
l <- length(fls)
if (l) {
  msg <- paste0(l, " input Raw files already present in environment from a previous run:")
  opt <- c("Remove them                                                                                                             ",
           "Reprocess them                                                                                                          ")
  startFresh <- c(TRUE, FALSE)[match(dlg_list(opt, opt[1], title = msg)$res, opt)]
  if (startFresh) { fls <- c() }
}
l <- length(fls)
if (!l) {
  filt <- matrix(data = c("Thermo raw file", "*.raw;*.RAW"), ncol = 2,
                 dimnames = list("raw file"))
  #fls <- "B:/archive/lsfgrp/MS/Acquired data/frimlgrp/LCMS_JiFrATeplova1_m"
  fls <- normalizePath(choose.files(paste0(dflt, "/*.raw"), filters = filt), winslash = "/")
  fls <- fls[order(file.info(fls)$mtime, decreasing = FALSE)]
  if (getInt) {
    tst <- parSapply(parClust, fls, function(fl) {
      rawrr::readFileHeader(fl)$`Creation date`
    }) # Get created time directly from raw file!
    tst <- as.POSIXct.default(tst, tz = Sys.timezone(), format = "%d/%m/%Y %H:%M:%S")
    fls <- fls[order(tst, decreasing = FALSE)]
  }
}
#a <- rawrr::readFileHeader(rawfile = fl)

# Work directory
wd <- unique(dirname(fls))[1]
dtst <- gsub(".*/", "", wd)
tst <- try(suppressWarnings(write("Test", paste0(wd, "/test.txt"))), silent = TRUE)
while ("try-error" %in% class(tst)) {
  wd <- rstudioapi::selectDirectory("Choose a work directory where we have write permission!", path = "D:/")
  tst <- try(suppressWarnings(write("Test", paste0(wd, "/test.txt"))), silent = TRUE)
}
unlink(paste0(wd, "/test.txt"))
setwd(wd)
clusterExport(parClust, "wd", envir = environment())

# Convert to mzML
deer <- list()
ParsDirs <- grep("/ThermoRawFileParser", c(list.dirs("C:/Program Files", full.names = TRUE, recursive = FALSE),
                                           list.dirs(paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads"), full.names = TRUE, recursive = FALSE)),
                 value = TRUE)
if (!length(ParsDirs)) {
  url <- "https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.4/ThermoRawFileParser1.4.4.zip"
  require(curl)
  dstfl <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/ThermoRawFileParser1.4.4.zip")
  curl_download(url, dstfl)
  ParsDirs <- "C:/ThermoRawFileParser/V_1.4.4"
  if (!dir.exists(ParsDirs)) { dir.create(ParsDirs, recursive = TRUE) }
  unzip(dstfl, exdir = ParsDirs)
}
if (length(ParsDirs) > 1) {
  w <- which(vapply(ParsDirs, function(x) { "ThermoRawFileParser.exe" %in% list.files(x) }, TRUE))
  ParsDirs <- ParsDirs[w]
}
if (length(ParsDirs) > 1) {
  tst <- sapply(ParsDirs, function(x) { file.info(x)$ctime })
  ParsDirs <- ParsDirs[which(tst == max(tst))[1]]
}
deer$ParsDir <- ParsDirs
MSConvertInst <- ("C:/Program Files/ProteoWizard"%in% list.dirs("C:/Program Files", recursive = FALSE))
if (!MSConvertInst) {
  tst <- grep("ProteoWizard ", list.dirs("C:/Users/Thermo/AppData/Local/Apps", full.names = FALSE), value = TRUE)
  if (length(tst)) {
    MSConvertInst <- TRUE
    MSConvertDir <- paste0("C:/Users/Thermo/AppData/Local/Apps/", tst[1])
    deer$MSConvertDir <- MSConvertDir
  }
} else {
  MSConvertDirs <- grep("/ProteoWizard/ProteoWizard [^/]+",
                        list.dirs("C:/Program Files/ProteoWizard", recursive = FALSE), value = TRUE)
  if (!length(MSConvertDirs)) { MSConvertInst <- FALSE } else {
    #MSConvertDirs <- c("C:/Program Files/ProteoWizard/ProteoWizard 3.0.19172.57d620127",
    #                   "C:/Program Files/ProteoWizard/ProteoWizard 3.0.22099.89b871a",
    #                  "C:/Program Files/ProteoWizard/ProteoWizard 3.0.22317.1e024d4")
    MSConvertVers <- as.data.frame(t(sapply(strsplit(gsub(".*/ProteoWizard ", "", MSConvertDirs), "\\."), unlist)))
    MSConvertVers <- MSConvertVers[order(MSConvertVers$V1, MSConvertVers$V2, MSConvertVers$V3, MSConvertVers$V4, decreasing = TRUE),]
    MSConvertVers <- MSConvertVers[1,]
    MSConvertDir <- paste0("C:/Program Files/ProteoWizard/ProteoWizard ", paste(MSConvertVers, collapse = "."))
    deer$MSConvertDir <- MSConvertDir
  }
}
Convert_mode <- "thermorawfileparser"
zlib <- TRUE
PeakPicking <- TRUE
mzMLs <- paste0(wd, "/", gsub("\\.raw$", ".mzML", basename(fls), ignore.case = TRUE))
pressFls <- gsub("\\.mzML$", ".csv", mzMLs)
clusterExport(parClust, list("fls", "mzMLs", "pressFls"), envir = environment())
w <- which((!file.exists(mzMLs))|(file.size(mzMLs) <= 2000))
if (length(w)) {
  if (tolower(Convert_mode) == "thermorawfileparser") { # Mode 1: using ThermoRawFileParser
    clusterExport(parClust, list("deer", "zlib", "PeakPicking"), envir = environment())
    tst <- parSapply(parClust, w, function(i) { #i <- w[1]
      cmd <- paste0("\"", deer$ParsDir, "/ThermoRawFileParser.exe\" -i=\"",
                    fls[i], "\" -b=\"", gsub(".*/", paste0(wd, "/"), mzMLs[i]), "\" -f=2 -a",
                    c(" -z", "")[zlib+1], c(" -p", "")[PeakPicking+1])
      write(cmd, paste0(wd, "/tmp", i, ".bat"))
      cmd2 <- paste0("\"", wd, "/tmp", i, ".bat\"") # I have to go through an intermediate batch file to run cmd,
      # because it is one of those which works in Windows command line but not when passed to system() or shell() in R.
      #cat(cmd)
      #writeClipboard(cmd)
      system(cmd2)
      unlink(paste0(wd, "/tmp", i, ".bat"))
    })
  }
  if (tolower(Convert_mode) == "msconvert") { # Mode 2: using msconvert
    write(fls[w], file = paste0(wd, "/tmp_MS_files.txt"))
    precRecal <- FALSE
    cmd <- paste0("\"", deer$MSConvertDir, "/msconvert.exe\" -f \"", wd, "/tmp_MS_files.txt\" -o \"",
                  wd, "\" --mzML --64 -z --filter \"peakPicking true 1-\"")
    #cat(cmd)
    system(cmd)
    unlink(paste0(wd, "/tmp_MS_files.txt"))
  }
}

# Files range
nFls <- length(fls)
sz <- 24
N <- ceiling(nFls/sz)
Ranges <- lapply(1:N, function(x) {
  x <- (1:sz)+(x-1)*sz
  x[which(x <= nFls)]
})

# Getting the pressure profile is unfortunately not doable using rawrr... but with xml2, using the mzML, it is feasible!
# The R functions here are based on the Python script at
# https://gist.githubusercontent.com/caetera/0921b33f0c6201a538436906cc965fff/raw/d1af134fe228ce6a23be3e5bc3c49d20a8447ab2/mzML_pressure_to_csv.py
# translated to R with sypport from chagPT:
mzML_scan <- function(file # mzML file
                      ,
                      ont = "MS:1003019" # Which ontology to scan for, the default defines pressure
                      ) {
  library(xml2)
  # parse mzML
  doc <- xml2::read_xml(file)
  #
  # read namespace & normalize prefix
  ns <- xml2::xml_ns(doc)
  if ("d1" %in% names(ns)) {
    ns <- c(mzml = ns[["d1"]])
  } else if ("" %in% names(ns)) {
    ns <- c(mzml = ns[[""]])
  }
  #
  # find pressure chromatograms
  xpath_chrom <- paste0(".//mzml:chromatogram[mzml:cvParam[@accession='", ont, "']]")
  chroms <- xml2::xml_find_all(doc, xpath_chrom, ns)
  # XPath: chromatograms containing a cvParam with accession MS:1003019
  xpath <- paste0(".//", names(ns)[1], ":chromatogram[", names(ns)[1], ":cvParam[@accession='", ont, "']]")
  ontChroms <- xml2::xml_find_all(doc, xpath, ns = ns)
  ids <- xml2::xml_attr(ontChroms, "id")
  return(ids)
}
mzML_getChrom <- function(file # mzML file
                          ,
                          chroms # Chromatograms to extract. If missing, all chromatograms matching the ontology are extracted.
                          ,
                          ont = "MS:1003019" # Which ontology to scan for, the default defines pressure
                          ) {
  library(xml2)
  library(base64enc)
  # parse mzML
  doc <- xml2::read_xml(file)
  #
  # read namespace & normalize prefix
  ns <- xml2::xml_ns(doc)
  if ("d1" %in% names(ns)) {
    ns <- c(mzml = ns[["d1"]])
  } else if ("" %in% names(ns)) {
    ns <- c(mzml = ns[[""]])
  }
  #
  # find pressure chromatograms
  xpath_chrom <- paste0(".//mzml:chromatogram[mzml:cvParam[@accession='", ont, "']]")
  ontChroms <- xml2::xml_find_all(doc, xpath_chrom, ns)
  ids <- xml2::xml_attr(ontChroms, "id")
  if (!missing(chroms)) {
    ontChroms <- ontChroms[match(chroms, ids)]
  }
  Res <- lapply(ontChroms, function(Crom) { #Crom <- chroms[1]
    Crom_id <- xml2::xml_attr(Crom, "id")
    # find binary arrays for this chromatogram
    arrays <- xml2::xml_find_all(Crom, "mzml:binaryDataArrayList/mzml:binaryDataArray", ns)
    decoded <- list()
    for (array in arrays) {
      # extract base64 text
      bin_node <- xml2::xml_find_first(array, "mzml:binary", ns)
      bin_raw  <- base64enc::base64decode(xml2::xml_text(bin_node))
      # parse cvParams
      params <- xml2::xml_find_all(array, "mzml:cvParam", ns)
      kind <- NULL
      compression <- NULL
      size <- NULL
      for (p in params) {
        name_parts <- strsplit(xml2::xml_attr(p, "name"), " ")[[1]]
        last <- tail(name_parts, 1)
        if (last == "array") {
          kind <- name_parts[1]
        } else if (last == "compression") {
          compression <- name_parts[1]
        } else if (last == "float") {
          size <- name_parts[1]
        }
      }
      if (is.null(kind)) { stop("Array kind is not set") }
      # float size
      dtype <- switch(size,
                      "32-bit" = list(type = "float", bytes = 4),
                      "64-bit" = list(type = "double", bytes = 8),
                      stop(paste("Unknown float size:", size))
      )
      # decompress if needed
      if (compression == "zlib") {
        bin_raw <- memDecompress(bin_raw, type = "gzip")
      } else {
        if (!identical(compression, "no")) { stop(paste("Unknown compression method:", compression)) }
      }
      # decode array
      n_vals <- length(bin_raw)/dtype$bytes
      vec <- readBin(bin_raw, what = dtype$type, n = n_vals, size = dtype$bytes, endian = "little")
      decoded[[kind]] <- vec
    }
    # Output CSV like in Python
    rs <- as.data.frame(decoded)
    colnames(rs) <- c("Retention time", "Pressure")
    return(rs)
  })
  names(Res) <- chroms
  #View(Res[[1]])
  return(Res)
}
clusterExport(parClust, list("mzML_scan", "mzML_getChrom", "mzMLs"), envir = environment())
allPressChrom <- parLapply(parClust, mzMLs, mzML_scan)
whichPressChrom <- allPressChrom <- unique(unlist(allPressChrom))
whichPressChromDum <- vapply(whichPressChrom, function(x) {
  paste(c(x, rep(" ", 250-nchar(x))), collapse = "")
}, "")
if (length(whichPressChrom) > 1) {
  whichPressChromDum <- dlg_list(whichPressChromDum, whichPressChromDum[1], title = "Select pressure chromatogram(s) to extract", multiple = TRUE)$res
  whichPressChrom <- gsub(" +$", "", whichPressChromDum)
}
if (length(whichPressChrom)) {
  myPressChrom <- parLapply(parClust, mzMLs, mzML_getChrom, chroms = whichPressChrom)
  names(myPressChrom) <- mzMLs
}
pressChrom <- setNames(lapply(whichPressChrom, function(Crom) { #Crom <- whichPressChrom[1]
  a <- setNames(lapply(1:length(fls), function(i) { #i <- 1
    mzML <- mzMLs[i]
    fl <- fls[i]
    x <- myPressChrom[[mzML]][[Crom]]
    x$"Raw file" <- factor(fl, levels = fls)
    x$"Raw file name" <- factor(cleanRawNm(basename(fl)), levels = cleanRawNm(basename(fls)))
    return(x)
  }), fls)
}), whichPressChrom)
allChroms <- list(Pressure = pressChrom)
#allChroms$Pressure <- press
#
#chromtypes <- setNames(c("tic", "bpc", "xic"), c("TIC", "Base peak", "XIC"))
if (getInt) {
  tol <- as.numeric(dlg_input("Enter mass tolerance (ppm)", 20)$res)
  clusterCall(parClust, function() library(rawrr))
  MS2s <- parLapply(parClust, fls, function(fl) { #fl <- fls[1]
    x <- rawrr::readIndex(fl)
    x <- x[which(x$MSOrder == "Ms2"),]
    return(unique(gsub("@.*", "", unique(x$scanType))))
  })
  ms2s <- unique(unlist(MS2s))
  if (length(ms2s) == 1) {
    cat(paste0("Single MS2 filter detected:\"", ms2s, "\"\n"))
  } else {
    tmp <- sapply(ms2s, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") })
    tmp <- dlg_list(tmp, tmp[1], TRUE, "Select filter(s) for which you want to get an XIC")$res
    ms2s <- gsub(" *$", "", tmp)
  }
  masses <- setNames(lapply(ms2s, function(ms2) { #ms2 <- ms2s[1]
    msg <- paste0(ms2, " XIC: enter M/Z to extract (if multiple M/Zs, separate them with \" \")")
    mass <- as.numeric(unlist(strsplit(dlg_input(msg, 134.0472)$res, " ")))
    return(mass[which(!is.na(mass))])
  }), ms2s)
  # TIC
  tic <- setNames(parLapply(parClust, fls, ticFun), fls)
  w <- sapply(tic, function(x) { x$Outcome })
  tic <- lapply(tic[w], function(x) { x$Output })
  if (length(tic)) { allChroms$TIC <- tic }
  # BPC
  bpc <- setNames(parLapply(parClust, fls, bpcFun), fls)
  w <- sapply(bpc, function(x) { x$Outcome })
  bpc <- lapply(bpc[w], function(x) { x$Output })
  if (length(bpc)) { allChroms$BPC <- bpc }
  # XICs
  clusterExport(parClust, list("ms2s", "tol", "masses", "xicFun"), envir = environment())
  fltms2 <- listMelt(masses, names(masses))
  clusterExport(parClust, list("tol", "xicFun"), envir = environment())
  xicNms <- do.call(paste, c(fltms2, sep = " for "))
  xic <- setNames(apply(fltms2, 1, function(x) {
    mass <- x[[1]]
    ms2 <- x[[2]]
    clusterExport(parClust, list("ms2", "mass"), envir = environment())
    setNames(parLapply(parClust, fls, function(fl) { xicFun(fl, mass, tol, ms2) }), fls)
  }), xicNms)
  w <- setNames(lapply(xicNms, function(x) { #x <- names(xic)[1]
    which(sapply(xic[[x]], function(y) { y$Outcome }))
  }), xicNms)
  xic <- setNames(lapply(xicNms, function(x) { #x <- names(xic)[1]
    setNames(lapply(fls[w[[x]]], function(fl) { xic[[x]][[fl]]$Output }), fls[w[[x]]])
  }), xicNms)
  w <- which(setNames(lapply(xicNms, length), xicNms) > 0)
  xic <- xic[w]
  for (nm in names(xic)) { allChroms[[nm]] <- xic[[nm]] }
}
#
for (nm in names(allChroms)) { #nm <- names(allChroms)[1]
  if (nm == "Pressure") { #nm <- "Pressure"
    nms2 <- names(allChroms[[nm]])
    allChroms[[nm]] <- setNames(lapply(nms2, function(nm2) { #nm2 <- nms2[1]
      temp <- lapply(fls, function(fl) { #fl <- fls[1]
        tmp <- as.data.table(allChroms[[nm]][[nm2]][[fl]])
        tmp$Bin <- round(tmp$`Retention time`, 1)
        tmp <- tmp[, list(Pressure = mean(Pressure)),
                   by = list(`Raw file` = `Raw file`,
                             `Raw file name` = `Raw file name`,
                             `Retention time` = Bin)]
        return(tmp)
      })
      temp <- do.call(rbind, temp)
      #temp <- data.table::dcast(temp, `Retention time`~`Raw file name`)
      temp <- as.data.frame(temp)
      return(temp)
    }), nms2)
  } else {
    allChroms[[nm]] <- plyr::rbind.fill(allChroms[[nm]])
  }
}
#
fullRTRange <- unlist(lapply(names(allChroms), function(nm) {
  if (nm == "Pressure") {
    return(unlist(lapply(allChroms[[nm]], function(x) {
      x$`Retention time`
    })))
  }
  return(allChroms[[nm]]$`Retention time`)
}))
fullRTRange <- c(min(fullRTRange), max(fullRTRange))
getRTRange <- TRUE; kount <- 1
while (getRTRange) {
  for (nm in names(allChroms)) { #nm <- names(allChroms)[1] #nm <- names(allChroms)[2]
    if ("data.frame" %in% class(allChroms[[nm]])) {
      tmpDat <- list()
      tmpDat[[nm]] <- allChroms[[nm]]
    }
    if ("list" %in% class(allChroms[[nm]])) {
      tmpDat <- allChroms[[nm]]
    }
    for (nm2 in names(tmpDat)) { #nm2 <- names(tmpDat)[1]
      x <- tmpDat[[nm2]]
      if (kount > 1) {
        x <- x[which((x$`Retention time` >= RTRange[1])&(x$`Retention time` <= RTRange[2])),]
      }
      xLim <- c(min(x$`Retention time`, na.rm = TRUE), max(x$`Retention time`, na.rm = TRUE))
      rgX <- xLim[2]-xLim[1]
      Ykol <- c("Intensity", "Pressure")[(nm == "Pressure")+1]
      yLim <- c(min(x[[Ykol]], na.rm = TRUE), max(c(x[[Ykol]], 100), na.rm = TRUE)) # Hard-coded minimum 100 max Y to avoid issues when values are very low!
      rgY <- yLim[2]-yLim[1]
      if ("data.frame" %in% class(x)) {
        if (!nm %in% c("TIC", "BPC", "Pressure")) {
          tmp <- unlist(strsplit(nm2, " for "))
          mass <- as.numeric(tmp[1])
          ms2 <- tmp[2]
          ms2Trg <- as.numeric(gsub(".* ", "", ms2))
        }
        x$`Raw file` <- factor(x$`Raw file`, levels = fls)
        x$`Raw file name` <- factor(x$`Raw file name`, levels = cleanRawNm(basename(fls)))
        #offset <- max(x$Int)/3
        #nc <- nchar(floor(offset))
        #offset <- ceiling(offset/10^(nc-1))*10^(nc-1)
        #if (!is.finite(offset)) { offset <- 0 }
        #x$Int <- x[[Ykol]] + offset*(match(x$`Raw file`, fls)-1)
        #unique(x$`Raw file`)
        #
        for (i in 1:length(Ranges)) { #i <- 1
          rg <- Ranges[[i]]
          w <- which(x$`Raw file` %in% fls[rg])
          if (length(w)) {
            if (nm %in% c("TIC", "BPC", "Pressure")) { ttl <- paste0(dtst, " - ", nm2) } else {
              ttl <- paste0(dtst, "- XIC, transition ", ms2Trg, " to ", mass)
            }
            if (kount > 1) { ttl <- paste0(ttl, ", RT ", RTRange[1], "-", RTRange[2], " min") }
            if (length(Ranges) > 1) { ttl <- paste0(ttl, ", files ", min(rg), "-", max(rg)) }
            plot <- ggplot(x[w,], aes(x = `Retention time`, y = .data[[Ykol]], colour = `Raw file name`)) + ggtitle(ttl) +
              coord_fixed(rgX/(rgY*3)) +
              facet_wrap(~`Raw file name`) +
              #geom_point() +
              geom_line(aes(group = `Raw file name`)) +
              theme_bw()  +
              theme(strip.text.y = element_text(angle = 0)) +
              xlim(xLim[1], xLim[2]) + ylim(yLim[1], yLim[2])
            # Only pop up pressure and XIC plots
            if (!nm %in% c("TIC", "BPC")) { windows(22, 12); print(plot) }
            ggsave(paste0(wd, "/", ttl, ".jpeg"), plot, dpi = 100)
          }
        }
      } else {
        if (nm %in% c("TIC", "BPC")) {
          warning(paste0("Failure to extract ", nm2))
        } else {
          warning(paste0("Failure to extract XICs @ mass = ", mass, " Da, tol. = ", tol, ", filter = '", ms2, "'"))
        }
      }
    }
  }
  kount <- kount + 1
  msg <- paste0("Extract a", c("", "nother")[(kount > 2)+1], " RT sub-range?")
  getRTRange <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  if (getRTRange) {
    if (exists("RTRange")) { dflt <- RTRange } else { dflt <- fullRTRange }
    RTRange <- as.numeric(unlist(strsplit(dlg_input("Enter min. and max. RT divided by \"-\"",
                                                    paste(dflt, collapse = "-"))$res, " *- *")))
    RTRange <- RTRange[which(!is.na(RTRange))]
    if (length(RTRange) != 2) {
      warning("Invalid RT range! Cancelling...")
      getRTRange <- FALSE
    } else {
      RTRange <- sort(RTRange)
      if (RTRange[1] < fullRTRange[1]) {
        warning("Too low RT range!")
        RTRange[1] <- fullRTRange[1]
      }
      if (RTRange[2] > fullRTRange[2]) {
        warning("Too high RT range!")
        RTRange[2] <- fullRTRange[2]
      }
    }
  }
}
#nms <- names(allChroms)[which(names(allChroms) != "Pressure")] # No need to write pressure, we have it already
nms <- names(allChroms) # Actually let's write it, it's spread over different files
sapply(nms, function(nm) { #nm <- "BPC"
  nm2 <- paste0(nm, ".csv")
  if (!nm %in% c("TIC", "BPC", "Pressure")) { nm2 <- paste0("XIC ", nm2) }
  x <- allChroms[[nm]]
  x$`Raw file name` <- NULL
  x <- x[, c("Raw file", colnames(x)[which(colnames(x) != "Raw file")])]
  data.table::fwrite(x, paste0(wd, "/", nm2), row.names = FALSE, na = "NA")
})

# To do:
# - Parallelize plots printing bit
# - Shiny-based app to display single graph for file of interest
# - Propagate lessons learned in this script to "Convert Raw files..." script
