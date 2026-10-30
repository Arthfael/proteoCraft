# Convert Thermo .raw files or Bruker .d folders to mzML (or other) as a preliminary to a FragPipe or OpenMS search:
# NB: all files are expected to be in a single location!
options(stringsAsFactors = FALSE)
#library(proteoCraft)

# Pattern replacement function
if (!require(proteoCraft)) {
  topattern <- \(x, start = TRUE, end = FALSE, collapse = "|") {
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
    if ((length(x) > 1L) && (collapse != FALSE)) { x <- paste(x, collapse = collapse) }
    return(x)
  }
} else {
  topattern <- proteoCraft::topattern
}
# Convenient function to "fix" directory formats
dirfix <- \(dir, win = FALSE, quotes = FALSE) {
  dir <- sub("/$", "", gsub("/+", "/", dir))
  if (win) {
    dir <- sub("/Program Files(x86)/", "/PROGRA~2/",
               sub("/Program files/", "/PROGRA~1/", dir, ignore.case = TRUE),
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
  dstfl <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/ThermoRawFileParser-master.zip")
  curl::curl_download(url, dstfl)
  ParsDirs <- sub("\\.zip$", "", dstfl)
  utils::unzip(dstfl, exdir = ParsDirs)
}
if (length(ParsDirs) > 1L) {
  tst <- sapply(ParsDirs, \(x) {
    file.info(x)$ctime
  })
  ParsDirs <- ParsDirs[which(tst == max(tst))[1L]]
}
deer$ParsDir <- ParsDirs
MSConvertInst <- ("C:/Program Files/ProteoWizard"%in% list.dirs("C:/Program Files", recursive = FALSE))
if (MSConvertInst) {
  MSConvertDirs <- grep("/ProteoWizard/ProteoWizard [^/]+", list.dirs("C:/Program Files/ProteoWizard", recursive = FALSE), value = TRUE)
  if (!length(MSConvertDirs)) { MSConvertInst <- FALSE } else {
    
    
    
    MSConvertVers <- as.data.frame(t(sapply(strsplit(sub(".*/ProteoWizard ", "", MSConvertDirs), "\\."), unlist)))
    MSConvertVers <- MSConvertVers[order(MSConvertVers$V1, MSConvertVers$V2, MSConvertVers$V3, MSConvertVers$V4, decreasing = TRUE),]
    MSConvertVers <- MSConvertVers[1L,]
    MSConvertDir <- paste0("C:/Program Files/ProteoWizard/ProteoWizard ", paste(MSConvertVers, collapse = "."))
    deer$MSConvertDir <- MSConvertDir
  }
}




formatz <- c("mgf", "mzml", "indexedmzml", "parquet", "mzxml", paste0("ms", 1L:2L), "mz5", "text",
             paste0("cms", 1L:2L), "log")
extensionz2 <- c("mgf", "mzML", "indexed-mzML", "parquet", "mzXML", paste0("ms", 1L:2L), "mz5", "txt",
                 paste0("cms", 1L:2L), "log")
extensionz <- sub("^indexed-", "", extensionz2)
#alks <- c("Iodoacetamide (IAA)", "N-ethylmaleimide (NEM)", "Methyl Methanethiosulfonate (MMTS)")
#alks2 <- c("iaa", "nem", "mmts")
converters <- c("ThermoRawFileParser", "msconvert", "RawConverter")
convertersexe <- c("ThermoRawFileParser.exe", "msconvert.exe", "RawConverter.exe")
convertersdir <- c("ParsDir", "MSConvertDir", "RawConvertDir")
Software <- c("FragPipe", "OpenMS", "MSGF+", "Other")
Modes <- setNames(c("Bruker .d files (folders)                                                 ",
                    "Thermo .raw file                                                          "),
                  c("Bruker", "Thermo"))
Mode <- svDialogs::dlg_list(Modes, Modes[1L], FALSE, "Which file types do you want to convert?")$res
Mode <- names(Modes)[match(Mode, Modes)]
WD_backup <- getwd()
if (Mode == names(Modes)[1L]) {
  FILES <- c()
  Moar <- TRUE
  dflt <- getwd()
  while (Moar) {
    msg <- c("Select input Bruker .d folder.", "Select another Bruker .d folder, or escape to stop adding more.")[length(FILES > 0)+1]
    fl <- rstudioapi::selectDirectory(msg, path = dflt)
    Moar <- length(fl) > 0L
    if (Moar) {
      dflt <- sub("/[^/]+$", "", fl)
      FILES <- unique(c(FILES, fl))
    }
  }
}
if (Mode == names(Modes)[2L]) {
  FILES <- gsub("\\\\", "/", choose.files("D:\\groups_temp\\."))
}
check.ext <- sum(!tolower(sapply(strsplit(FILES, "\\."), \(x) { rev(x)[1L] })) %in% c("raw", "d")) == 0L
if (!check.ext) { stop("Only raw files are accepted!") } else {
  WD <- unique(sub("/[^/]+$", "", FILES))
  if (length(WD) > 1L) { stop("Only one local folder is allowed for now!") } else {
    # For simplicity, this will only accept a single folder location.
    # There can be issues with some of the converters otherwise.
    # Since their syntax is unique to each, it's simpler to just request all files to be in one location.
    setwd(WD)
    deer$Dir <- WD
    FILES <- sub(paste0(topattern(WD), "/?"), "", FILES)
    msg <- paste0("Enter the number corresponding to the output file type you want to convert to:",
                  paste(vapply(1L:length(extensionz2), \(x) { paste0(x, " = ", extensionz2[x]) }, ""), collapse = "; "))
    format <- formatz[as.numeric(svDialogs::dlg_input(msg, 3L)$res)]
    wf <- which(formatz == format)
    ext <- extensionz[wf]
    #msg <- paste0("Enter the number of alkylating agent used:",
    #            paste(sapply(1L:length(alks), \(x) { paste0(x, " = ", alks[x]) }), collapse = "; "))
    #alk <- alks2[as.numeric(svDialogs::dlg_input(msg, 2)$res)]
    #alk.mass <- c(57.0214637236, 125.0476784741, 45.9877207542)[which(c("iaa", "nem", "mmts") == alk)]
    msg <- "Which software will the files be searched with?\n(Choose \"Other\" to ignore default settings)"
    SearchSoft <- svDialogs::dlg_list(Software, "FragPipe", title = msg)$res
    stopifnot(length(SearchSoft %in% Software) == 1L, SearchSoft %in% Software)
    if (SearchSoft == "FragPipe") {
      if (MSConvertInst) { Convert_mode <- "msconvert" } else {
        stop("Could not find any valid installation of msconvert!\n ")
      }
    }
    if (Mode == names(Modes)[1L]) { Convert_mode <- "msconvert" }
    if (SearchSoft == "OpenMS") { Convert_mode <- "ThermoRawFileParser" }
    if ((Mode == names(Modes)[2L]) && (SearchSoft == "Other")) {
      if (!MSConvertInst) {
        wY <- which(converters != "msconvert")
        converters <- converters[wY]
        convertersexe <- convertersexe[wY]
      }
      msg <- "Select a converter:"
      Convert_mode <- svDialogs::dlg_list(converters, 1L, title = msg)$res
    }
    w <- which(converters == Convert_mode)
    if (!convertersexe[w] %in% list.files(deer[[convertersdir[w]]])) {
      stop(paste0("Could not find \"", convertersexe[w], "\" at the specified location!"))
    } else {
      for (n in names(deer)) { deer[[n]] <- dirfix(deer[[n]], TRUE, FALSE) }
      PeakPicking <- TRUE
      zlib <- FALSE
      if (SearchSoft == "FragPipe") { zlib <- TRUE }
      if ((SearchSoft == "Other") && (Convert_mode %in% c("msconvert", "ThermoRawFileParser"))) {
        msg <- "Should we perform zlib compression (y/n)?"
        zlib <- c(FALSE, TRUE)[match(svDialogs::dlg_input(msg, "y")$res, c("n", "y"))]
        msg <- "Should we perform peak-picking (y/n)?"
        PeakPicking <- c(FALSE, TRUE)[match(svDialogs::dlg_input(msg, "y")$res, c("n", "y"))]
      }
      files <- paste0(sub("\\.raw$", "", FILES, ignore.case = TRUE), ".", ext)
      tst <- files %in% list.files(deer$Dir)
      if (sum(tst)) {
        msg <- "Some output files already exist in the directory, should we overwrite them (y/n)?"
        Convert_overwrite <- c(FALSE, TRUE)[which(c("n", "y") == svDialogs::dlg_input(msg, "y")$res)]
        if (!Convert_overwrite) { FILES <- FILES[which(!tst)] }
      }
      if (length(FILES)) {
        if (tolower(Convert_mode) == "thermorawfileparser") { # Mode 1: using ThermoRawFileParser
          if (!ext %in% c("mgf", "mzML", "parquet")) {
            stop(paste0("ThermoRawFileParser cannot convert Raw files to ", ext, "!"))
          } else {
            wf2 <- match(format, c("mgf", "mzml", "indexedmzml", "parquet")) - 1L
            for (i in 1L:length(FILES)) { #i <- 1L
              cmd <- paste0("\"", deer$ParsDir, "/ThermoRawFileParser.exe\" -i=\"",
                            deer$Dir, "/", FILES[i], "\" -b=\"", deer$Dir, "/", files[i], "\" -f=", wf2,
                            c(" -z", "")[zlib+1L], c(" -p", "")[PeakPicking+1L])
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
            x <- if (ext == "txt") { "text" } else { ext }
            # For some reason I have to set the R working directory to the one where the files are
            # for this to work, I cannot get it to work by appending the folder locations to the files!
            if ((SearchSoft == "FragPipe") || (Mode == names(Modes)[1L])) {
              precRecal <- FALSE
            } else {
              msg <- "Should we perform Precursor Recalculation (y/n)? - note: only on Orbitrap or FT data, and can cause issues, so use with caution."
              precRecal <- c(FALSE, TRUE)[which(c("n", "y") == svDialogs::dlg_input(msg, "n")$res)]
            }
            cmd <- paste0("\"", deer$MSConvertDir, "/msconvert.exe\" -f tmp_MS_files.txt -o \"",
                          deer$Dir, "\" --", ext, " --64",
                          c("", " --noindex")[(format == "mzml")+1L], c("", " -z")[zlib+1L],
                          c("", " --filter \"peakPicking true 1-\"")[PeakPicking+1L],
                          c("", " --filter \"precursorRecalculation\"")[precRecal+1L])
            #cat(cmd)
            system(cmd)
            unlink("tmp_MS_files.txt")
          }
        }
        
        if (tolower(Convert_mode) == "rawconverter") { # Mode 3: using RawConverter
          
          # Needs fixing!!!
          
          if (!ext %in% c("mgf", "mzML", "mzXML", paste0("ms", 1L:3L), "log")) {
            stop(paste0("RawConverter cannot convert Raw files to ", ext, "!"))
          } else {
            for (i in 1L:length(FILES)) { #i <- 1L
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
