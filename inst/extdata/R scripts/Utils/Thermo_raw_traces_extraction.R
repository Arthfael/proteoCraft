#
# NB: this script assumes a python installation with numpy and pandas is available.
# If not, pressure profile extraction will fail, since it relies on a python script!
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
parXprs <- expression({
  if (exists("parClust")) { try(parallel::stopCluster(parClust), silent = TRUE) } # A fresh start if re-running
  N.clust <- parallel::detectCores()-1
  a <- 1
  tst <- try(parallel::clusterExport(parClust, "a", envir = environment()), silent = TRUE)
  if ("try-error" %in% class(tst)) {
    try(parallel::stopCluster(parClust), silent = TRUE)
    parClust <- parallel::makeCluster(N.clust, type = "SOCK")
  }
})
eval(parXprs, .GlobalEnv)

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
    cleanRawNm <- function(rawFileNm) { gsub("\\.raw$", "", rawFileNm, ignore.case = TRUE) }
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
    cleanRawNm <- function(rawFileNm) { gsub("\\.raw$", "", rawFileNm, ignore.case = TRUE) }
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
    cleanRawNm <- function(rawFileNm) { gsub("\\.raw$", "", rawFileNm, ignore.case = TRUE) }
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
dflt <- "...Search_Folder"
if (!dir.exists(dflt)) { dflt <- "D:/Data" }
if (!dir.exists(dflt)) { dflt <- "...Projects_Doc_Folder/Mass_Spec/Acquired data_v2" }
if (!dir.exists(dflt)) { dflt <- "...MS_File_Archive" }

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
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
  wd <- rstudioapi::selectDirectory("Choose a work directory where we have write permission!", path = "D:/")
  tst <- try(suppressWarnings(write("Test", paste0(wd, "/test.txt"))), silent = TRUE)
}
unlink(paste0(wd, "/test.txt"))
setwd(wd)
clusterExport(parClust, "wd", envir = environment())

# Convert to mzML
deer <- list()
ParsDirs <- grep("/ThermoRawFileParser1\\.4", # It seems I should restrict to specific sub-versions as not all appear to work
                 c(list.dirs("C:/ThermoRawFileParser", full.names = TRUE, recursive = FALSE),
                   list.dirs("C:/Program Files", full.names = TRUE, recursive = FALSE),
                   list.dirs(paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads"), full.names = TRUE, recursive = FALSE)),
                 value = TRUE)
tst <- vapply(ParsDirs, function(x) { "ThermoRawFileParser.exe" %in% list.files(x) }, TRUE)
ParsDirs <- ParsDirs[which(tst)]
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
  tst <- vapply(ParsDirs, function(x) { file.info(x)$ctime }, 1)
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
    clusterExport(parClust, list("wd", "deer", "zlib", "PeakPicking"), envir = environment())
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

# Getting the pressure profile is unfortunately not doable using rawrr
# This must be done from mzMLs
require(XML)
#tst <- grabMzmlData(mzMLs[1], "everything")
#View(tst$metadata)
pressScript <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/mzML_pressure_to_csv.py")
if (!file.exists(pressScript)) {
  download.file("https://gist.githubusercontent.com/caetera/0921b33f0c6201a538436906cc965fff/raw/d1af134fe228ce6a23be3e5bc3c49d20a8447ab2/mzML_pressure_to_csv.py",
                pressScript)
}
stopifnot(file.exists(pressScript))
# Create temporary python pressure script edited to create a more predictably named output
tmp <- suppressWarnings(readLines(pressScript))
w <- which(tmp == "        csv_path = file.replace('.mzML', f'_{convert_to_safe_filename(chrom_name)}.csv')")
tmp[w] <- "        csv_path = file.replace('.mzML', f'.csv')"
pressScript2 <- paste0(wd, "/pressScript.py")
write(tmp, pressScript2)
w <- which(file.exists(mzMLs))
clusterExport(parClust, list("pressScript2", "mzMLs"), envir = environment())
tst <- parSapply(parClust, w, function(i) { #i <- w[1]
  cmd <- paste0("python \"", pressScript2, "\" ", paste0("\"", mzMLs[i], "\"", collapse = " "))
  #cat(cmd)
  system(cmd)
})
#
cmd <- paste0("open \"", wd, "\"")
#cat(cmd)
#system(cmd)
tstPress <- file.exists(pressFls)
#stopifnot(length(w) == 0)
if (sum(!tstPress)) {
  warning(paste0("Couldn't get pressure profile from the following file(s):",
                 paste0("\n - ", fls[which(!tstPress)], collapse = ""), "\n"))
}
unlink(pressScript2) # Remove temporary python script 
#
clusterExport(parClust, "cleanRawNm", envir = environment())
press <- setNames(parLapply(parClust, which(tstPress), function(i) {
  x <- data.table::fread(pressFls[i], integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  fl <- fls[i]
  colnames(x) <- c("Retention time", "Pressure")
  x$"Raw file" <- factor(fl, , levels = fls)
  x$"Raw file name" <- factor(cleanRawNm(basename(fl)), levels = cleanRawNm(basename(fls)))
  return(x)
}), fls[which(tstPress)])
allChroms <- list(Pressure = press)
#allChroms$Pressure <- press
#
#chromtypes <- setNames(c("tic", "bpc", "xic"), c("TIC", "Base peak", "XIC"))
if (getInt) {
  tol <- as.numeric(dlg_input("Enter mass tolerance (ppm)", 20)$res)
  #stopCluster(parClust)
  eval(parXprs, .GlobalEnv)
  clusterCall(parClust, function() library(rawrr))
  invisible(clusterCall(parClust, function() {
    library(rawrr)
    return()
  }))
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
  #stopCluster(parClust)
  eval(parXprs, .GlobalEnv)
  invisible(clusterCall(parClust, function() {
    library(rawrr)
    return()
  }))
  tic <- setNames(parLapply(parClust, fls, ticFun), fls)
  w <- sapply(tic, function(x) { x$Outcome })
  tic <- lapply(tic[w], function(x) { x$Output })
  if (length(tic)) { allChroms$TIC <- tic }
  # BPC
  #stopCluster(parClust)
  eval(parXprs, .GlobalEnv)
  invisible(clusterCall(parClust, function() {
    library(rawrr)
    return()
  }))
  bpc <- setNames(parLapply(parClust, fls, bpcFun), fls)
  w <- sapply(bpc, function(x) { x$Outcome })
  bpc <- lapply(bpc[w], function(x) { x$Output })
  if (length(bpc)) { allChroms$BPC <- bpc }
  # XICs
  #stopCluster(parClust)
  eval(parXprs, .GlobalEnv)
  invisible(clusterCall(parClust, function() {
    library(rawrr)
    return()
  }))
  fltms2 <- listMelt(masses, names(masses))
  colnames(fltms2) <- c("Mass", "Filter")
  fltms2 <- lapply(fls, function(fl) {
    x <- fltms2
    x$File <- fl
    return(x)
  })
  fltms2 <- do.call(rbind, fltms2)
  fltms2$Nms <- do.call(paste, c(fltms2[, c("Mass", "Filter")], sep = " for "))
  clusterExport(parClust, list("tol", "xicFun"), envir = environment())
  fltms2$xic <- parApply(parClust, fltms2[, c("Mass", "Filter", "File")], 1, function(x) { #x <- fltms2[1,]
    tst <- try({ xicFun(x[[3]], x[[1]], tol, x[[2]]) }, silent = TRUE)
    if ("try-error" %in% class(tst)) { return() }
    return(tst)
  })
  uNms <- unique(fltms2$Nms)
  xic <- setNames(lapply(uNms, function(nm) { #nm <- uNms[1]
    w <- which((fltms2$Nms %in% nm)&(vapply(fltms2$xic, function(x) { #x <- fltms2$xic[[1]]
      ("list" %in% class(x))&&
        (sum(c("Outcome", "Output") %in% names(x)) == 2)&&
        (x$Outcome)&&
        ("data.frame" %in% class(x$Output))&&
        (nrow(x$Output))
    }, TRUE)))
    setNames(lapply(w, function(x) { fltms2$xic[[x]]$Output }), fltms2$File[w])
  }), uNms)
  for (nm in names(xic)) { allChroms[[nm]] <- xic[[nm]] }
}
#
for (nm in names(allChroms)) { #nm <- names(allChroms)[1]
  allChroms[[nm]] <- plyr::rbind.fill(allChroms[[nm]])
  if (nm == "Pressure") {
    tmp <- as.data.table(allChroms[[nm]])
    tmp$Bin <- round(as.numeric(tmp$`Retention time`, 1))
    tmp <- tmp[, list(Pressure = mean(Pressure)),
               by = list(`Raw file` = `Raw file`,
                         `Raw file name` = `Raw file name`,
                         `Retention time` = Bin)]
    allChroms[[nm]] <- as.data.frame(tmp)
  }
}
#
fullRTRange <- unlist(lapply(allChroms, function(x) { x$`Retention time` }))
fullRTRange <- c(min(fullRTRange), max(fullRTRange))
getRTRange <- TRUE; kount <- 1
while (getRTRange) {
  for (nm in names(allChroms)) { #nm <- names(allChroms)[1]
    x <- allChroms[[nm]]
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
        tmp <- unlist(strsplit(nm, " for "))
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
          if (nm %in% c("TIC", "BPC", "Pressure")) { ttl <- paste0(dtst, " - ", nm) } else {
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
        warning(paste0("Failure to extract ", nm))
      } else {
        warning(paste0("Failure to extract XICs @ mass = ", mass, " Da, tol. = ", tol, ", filter = '", ms2, "'"))
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
