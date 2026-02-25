# Convert Thermo .raw files or Bruker .d folders to mzML (or other) as a preliminary to a FragPipe or OpenMS search:
# NB: all files are expected to be in a single location!
options(stringsAsFactors = FALSE)
library(parallel)
#library(proteoCraft)
library(svDialogs)

# Pattern replacement function
if (!require(proteoCraft)) {
  topattern <- function (x, start = TRUE, end = FALSE, collapse = "|") {
    x <- gsub("\\\\", "\\\\\\\\", as.character(x))
    x <- gsub("\\.", "\\\\.", x)
    x <- gsub("\\*", "\\\\*", x)
    x <- gsub("\\$", "\\\\$", x)
    x <- gsub("\\^", "\\\\^", x)
    x <- gsub("\\+", "\\\\+", x)
    x <- gsub("\\?", "\\\\?", x)
    x <- gsub("\\{", "\\\\{", x)
    x <- gsub("\\}", "\\\\}", x)
    x <- gsub("\\[", "\\\\[", x)
    x <- gsub("\\]", "\\\\]", x)
    x <- gsub("\\(", "\\\\(", x)
    x <- gsub("\\)", "\\\\)", x)
    x <- gsub("\\|", "\\\\|", x)
    if (start) { x <- paste0("^", x) }
    if (end) { x <- paste0(x, "$") }
    if ((length(x) > 1) && (collapse != FALSE)) { x <- paste(x, collapse = collapse) }
    return(x)
  }
}
# Convenient function to "fix" directory formats
dirfix <- function(dir, win = FALSE, quotes = FALSE) {
  dir <- gsub("/$", "", gsub("/+", "/", dir))
  if (win) {
    dir <- gsub("Program Files(x86)", "PROGRA~2",
                gsub("Program files", "PROGRA~1", dir, ignore.case = TRUE),
                ignore.case = TRUE)
  }
  if (quotes) { dir <- paste0("\"", dir, "\"") }
  return(dir)
}
# Default parameters:
deer <- list(RawConvertDir = "C:/Program Files/RawConverter_x64")
ParsDirs <- grep("/ThermoRawFileParser", c(list.dirs("C:/Program Files", full.names = TRUE, recursive = FALSE),
                                           list.dirs(paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads"), full.names = TRUE, recursive = FALSE)),
                 value = TRUE)
if (!length(ParsDirs)) {
  url <- "https://github.com/compomics/ThermoRawFileParser/archive/refs/heads/master.zip"
  require(curl)
  dstfl <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/ThermoRawFileParser-master.zip")
  curl_download(url, dstfl)
  ParsDirs <- gsub("\\.zip$", "", dstfl)
  unzip(dstfl, exdir = ParsDirs)
}
if (length(ParsDirs) > 1) {
  tst <- sapply(ParsDirs, function(x) {
    file.info(x)$ctime
  })
  ParsDirs <- ParsDirs[which(tst == max(tst))[1]]
}
deer$ParsDir <- ParsDirs
MSConvertInst <- ("C:/Program Files/ProteoWizard"%in% list.dirs("C:/Program Files", recursive = FALSE))
if (MSConvertInst) {
  MSConvertDirs <- grep("/ProteoWizard/ProteoWizard [^/]+", list.dirs("C:/Program Files/ProteoWizard", recursive = FALSE), value = TRUE)
  if (!length(MSConvertDirs)) { MSConvertInst <- FALSE } else {
    
    
    
    MSConvertVers <- as.data.frame(t(sapply(strsplit(gsub(".*/ProteoWizard ", "", MSConvertDirs), "\\."), unlist)))
    MSConvertVers <- MSConvertVers[order(MSConvertVers$V1, MSConvertVers$V2, MSConvertVers$V3, MSConvertVers$V4, decreasing = TRUE),]
    MSConvertVers <- MSConvertVers[1,]
    MSConvertDir <- paste0("C:/Program Files/ProteoWizard/ProteoWizard ", paste(MSConvertVers, collapse = "."))
    deer$MSConvertDir <- MSConvertDir
  }
}




formatz <- c("mgf", "mzml", "indexedmzml", "parquet", "mzxml", paste0("ms", 1:2), "mz5", "text",
             paste0("cms", 1:2), "log")
extensionz2 <- c("mgf", "mzML", "indexed-mzML", "parquet", "mzXML", paste0("ms", 1:2), "mz5", "txt",
                 paste0("cms", 1:2), "log")
extensionz <- gsub("^indexed-", "", extensionz2)
#alks <- c("Iodoacetamide (IAA)", "N-ethylmaleimide (NEM)", "Methyl Methanethiosulfonate (MMTS)")
#alks2 <- c("iaa", "nem", "mmts")
converters <- c("ThermoRawFileParser", "msconvert", "RawConverter")
convertersexe <- c("ThermoRawFileParser.exe", "msconvert.exe", "RawConverter.exe")
convertersdir <- c("ParsDir", "MSConvertDir", "RawConvertDir")
Software <- c("FragPipe", "OpenMS", "MSGF+", "Other")
Modes <- setNames(c("Bruker .d files (folders)                                                 ",
                    "Thermo .raw file                                                          "),
                    c("Bruker", "Thermo"))
Mode <- dlg_list(Modes, Modes[1], FALSE, "Which file types do you want to convert?")$res
Mode <- names(Modes)[match(Mode, Modes)]
WD_backup <- getwd()
if (Mode == names(Modes)[1]) {
  FILES <- c()
  Moar <- TRUE
  dflt <- getwd()
  while (Moar) {
    msg <- c("Select input Bruker .d folder.", "Select another Bruker .d folder, or escape to stop adding more.")[length(FILES > 0)+1]
    fl <- rstudioapi::selectDirectory(msg, path = dflt)
    Moar <- length(fl) > 0
    if (Moar) {
      dflt <- gsub("/[^/]+$", "", fl)
      FILES <- unique(c(FILES, fl))
    }
  }
}
if (Mode == names(Modes)[2]) {
  FILES <- gsub("\\\\", "/", choose.files("D:\\groups_temp\\."))
}
check.ext <- sum(!tolower(sapply(strsplit(FILES, "\\."), function(x) { rev(x)[1] })) %in% c("raw", "d")) == 0
if (!check.ext) { stop("Only raw files are accepted!") } else {
  WD <- unique(gsub("/[^/]+$", "", FILES))
  if (length(WD) > 1) { stop("Only one local folder is allowed for now!") } else {
    # For simplicity, this will only accept a single folder location.
    # There can be issues with some of the converters otherwise.
    # Since their syntax is unique to each, it's simpler to just request all files to be in one location.
    setwd(WD)
    deer$Dir = WD
    FILES <- gsub(paste0("^", topattern(WD), "/?"), "", FILES)
    msg <- paste0("Enter the number corresponding to the output file type you want to convert to:",
                paste(sapply(1:length(extensionz2), function(x) { paste0(x, " = ", extensionz2[x]) }), collapse = "; "))
    format <- formatz[as.numeric(dlgInput(msg, 3)$res)]
    wf <- which(formatz == format)
    ext <- extensionz[wf]
    #msg <- paste0("Enter the number of alkylating agent used:",
    #            paste(sapply(1:length(alks), function(x) { paste0(x, " = ", alks[x]) }), collapse = "; "))
    #alk <- alks2[as.numeric(dlgInput(msg, 2)$res)]
    #alk.mass <- c(57.0214637236, 125.0476784741, 45.9877207542)[which(c("iaa", "nem", "mmts") == alk)]
    msg <- "Which software will the files be searched with?\n(Choose \"Other\" to ignore default settings)"
    SearchSoft <- dlg_list(Software, "FragPipe", title = msg)$res
    stopifnot(length(SearchSoft %in% Software) == 1, SearchSoft %in% Software)
    if (SearchSoft == "FragPipe") {
      if (MSConvertInst) { Convert_mode <- "msconvert" } else {
        stop("Could not find any valid installation of msconvert!\n ")
      }
    }
    if (Mode == names(Modes)[1]) { Convert_mode <- "msconvert" }
    if (SearchSoft == "OpenMS") { Convert_mode <- "ThermoRawFileParser" }
    if ((Mode == names(Modes)[2])&&(SearchSoft == "Other")) {
      if (!MSConvertInst) {
        wY <- which(converters != "msconvert")
        converters <- converters[wY]
        convertersexe <- convertersexe[wY]
      }
      msg <- "Select a converter:"
      Convert_mode <- dlg_list(converters, 1, title = msg)$res
    }
    w <- which(converters == Convert_mode)
    if (!convertersexe[w] %in% list.files(deer[[convertersdir[w]]])) {
      stop(paste0("Could not find \"", convertersexe[w], "\" at the specified location!"))
    } else {
      for (n in names(deer)) { deer[[n]] <- dirfix(deer[[n]], TRUE, FALSE) }
      PeakPicking <- TRUE
      zlib <- FALSE
      if (SearchSoft == "FragPipe") { zlib <- TRUE }
      if ((SearchSoft == "Other")&&(Convert_mode %in% c("msconvert", "ThermoRawFileParser"))) {
        msg <- "Should we perform zlib compression (y/n)?"
        zlib <- c(FALSE, TRUE)[match(dlgInput(msg, "y")$res, c("n", "y"))]
        msg <- "Should we perform peak-picking (y/n)?"
        PeakPicking <- c(FALSE, TRUE)[match(dlgInput(msg, "y")$res, c("n", "y"))]
      }
      files <- paste0(gsub("\\.raw$", "", FILES, ignore.case = TRUE), ".", ext)
      tst <- files %in% list.files(deer$Dir)
      if (sum(tst)) {
        msg <- "Some output files already exist in the directory, should we overwrite them (y/n)?"
        Convert_overwrite <- c(FALSE, TRUE)[which(c("n", "y") == dlgInput(msg, "y")$res)]
        if (!Convert_overwrite) { FILES <- FILES[which(!tst)] }
      }
      if (length(FILES)) {
        if (tolower(Convert_mode) == "thermorawfileparser") { # Mode 1: using ThermoRawFileParser
          if (!ext %in% c("mgf", "mzML", "parquet")) {
            stop(paste0("ThermoRawFileParser cannot convert Raw files to ", ext, "!"))
          } else {
            wf2 <- which(c("mgf", "mzml", "indexedmzml", "parquet") == format) - 1
            for (i in 1:length(FILES)) { #i <- 1
              cmd <- paste0("\"", deer$ParsDir, "/ThermoRawFileParser.exe\" -i=\"",
                            deer$Dir, "/", FILES[i], "\" -b=\"", deer$Dir, "/", files[i], "\" -f=", wf2,
                            c(" -z", "")[zlib+1], c(" -p", "")[PeakPicking+1])
              #cat(cmd)
              system(cmd)
            }
          }
        }
        if (tolower(Convert_mode) == "msconvert") { # Mode 2: using msconvert
          if (!ext %in% c("mgf", "mzML", "mzXML", "ms1", "ms2", "mz5", "txt", "cms1", "cms2")) {
            stop(paste0("msconvert cannot convert Raw files to ", ext, "!"))
          } else {
            # Using msconvert:
            #write(FILES, file = "tmp_MS_files.txt")
            write(paste0(WD, "/", FILES), file = "tmp_MS_files.txt")
            if (ext == "txt") { x <- "text" } else { x <- ext }
            # For some reason I have to set the R working directory to the one where the files are
            # for this to work, I cannot get it to work by appending the folder locations to the files!
            if ((SearchSoft == "FragPipe")||(Mode == names(Modes)[1])) {
              precRecal <- FALSE
            } else {
              msg <- "Should we perform Precursor Recalculation (y/n)? - note: only on Orbitrap or FT data, and can cause issues, so use with caution."
              precRecal <- c(FALSE, TRUE)[which(c("n", "y") == dlgInput(msg, "n")$res)]
            }
            cmd <- paste0("\"", deer$MSConvertDir, "/msconvert.exe\" -f tmp_MS_files.txt -o \"",
                          deer$Dir, "\" --", ext, " --64",
                          c("", " --noindex")[(format == "mzml")+1], c("", " -z")[zlib+1],
                          c("", " --filter \"peakPicking true 1-\"")[PeakPicking+1],
                          c("", " --filter \"precursorRecalculation\"")[precRecal+1])
            #cat(cmd)
            system(cmd)
            unlink("tmp_MS_files.txt")
          }
        }
        
        if (tolower(Convert_mode) == "rawconverter") { # Mode 3: using RawConverter
          
          # Needs fixing!!!
          
          if (!ext %in% c("mgf", "mzML", "mzXML", paste0("ms", 1:3), "log")) {
            stop(paste0("RawConverter cannot convert Raw files to ", ext, "!"))
          } else {
            for (i in 1:length(FILES)) { #i <- 1
              setwd(deer$Dir)
              system(paste0("cd ", dirfix(deer$Dir, TRUE, TRUE)))
              cmd <- paste0(dirfix(paste0(deer$RawConvertDir, "/RawConverter.exe"), TRUE, TRUE), " ", ext,
                            dirfix(paste0(deer$Dir, "/", FILES[i]), TRUE, TRUE), 
                            " --", ext, " --select_mono_prec")
              #cat(cmd)
              system(cmd)
            }
          }
        }
      } else {
        warning("No files were converted!")
      }
    }
  }
}
#proteoCraft::openwd(WD)
setwd(WD_backup)
#rm(list = ls())
