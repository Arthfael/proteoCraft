# Clone folder securely to a new location
options(stringsAsFactors = FALSE)

if (!require(tools)) { install.packages("tools") }
if (!require(rstudioapi)) { install.packages("rstudioapi") }
if (!require(svDialogs)) { install.packages("svDialogs") }
if (!require(fs)) { install.packages("fs") }
if (!require(snow)) { install.packages("snow") }
if (!require(parallel)) { install.packages("parallel") }
require(rstudioapi)
require(svDialogs)
require(fs)
require(snow)
require(parallel)
N.clust <- detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  cat("Creating parallel cluster for faster operations...\n")
  parClust <- makeCluster(N.clust, type = "SOCK")
}
MaxAttempts <- 5

# Start with clean slate
#rm(list = ls())

if (!require(proteoCraft)) {
  openwd <- function(Dir) {
    cmd <- paste0("explorer \"", normalizePath(Dir), "\"")
    #cat(cmd)
    suppressMessages(suppressWarnings(shell(cmd, intern = TRUE)))
  }
  topattern <- function(x, start = TRUE, end = FALSE, collapse = "|") {
    x <- gsub("\\\\", "\\\\\\\\", as.character(x)) # Should be first so as not to escape escape signs!
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
    if ((length(x) > 1)&&(collapse != FALSE)) { x <- paste(x, collapse = collapse) }
    return(x)
  }
} else { openwd <- proteoCraft::openwd }

if ((exists("TargDir"))&&(!is.null(TargDir))&&(dir.exists(TargDir))) {
  dflt <- TargDir
} else {
  dflt <- c("D:/Data/Projects", "...Search_Folder", "C:/", getwd())
  dflt <- dflt[which(dir.exists(dflt))[1]]
}
#
nrmPath4PS <- function(Path) {
  Path <- gsub("\\)", "`)", gsub("\\(", "`(", Path))
  Path <- gsub("\\]", "``]", gsub("\\[", "``[", Path))
  Path <- gsub("\\.", "`.", Path)
  Path <- gsub(",", "`,", Path)
  Path <- gsub("-", "`-", Path)
  Path <- gsub(" ", "` ", Path)
  Path <- gsub("\'", "`'", Path)
  return(Path)
}
# Currently time adjustment is off by default
adjstTimes <- FALSE

#
# Capture operator's name
Whodunnit <- ""
while(!nchar(Whodunnit)) {
  Whodunnit <- dlg_input("Enter your first and last names:\n(For now this script will operate under the assumption that users are honest.\nDo not make me change this policy!)")$res
}
Whodunnit <- paste(sapply(strsplit(Whodunnit, " +"), function(x) {
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2, nchar(x))))
}), collapse = " ")
#
TargDir <- selectDirectory("Select folder to clone", path = dflt)
if (!is.null(TargDir)) {
  msg <- "Please confirm that this folder is not currently being modified?
(e.g. do not run this if raw files are stil being recorded here!)"
  OngoingAcq <- !c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  if (!OngoingAcq) {
    Fls <- list.files(TargDir, recursive = TRUE, full.names = TRUE)
    if (length(Fls)) {
      #
      parDir <- gsub(".*/", "", TargDir)
      rootDir <- gsub("/[^/]+$", "", TargDir)
      #
      # Create hash function
      HashFun <- function(x) { as.character(tools::md5sum(x)) } #tools::md5sum takes a character path, not a connection
      clusterExport(parClust, "HashFun", envir = environment())
      #
      tmp <- as.data.frame(t(parSapply(parClust, Fls, file.info)))
      tmp$isdir <- NULL
      tmp$exe <- NULL
      tmp$atime <- NULL
      tmp$mode <- NULL
      FilesDF <- data.frame(Files = Fls)
      FilesDF[, colnames(tmp)] <- tmp
      FilesDF$size <- unlist(FilesDF$size)
      FilesDF$ctime <- unlist(FilesDF$ctime)
      FilesDF$mtime <- unlist(FilesDF$mtime)
      FilesDF$flsCrTm <- as.POSIXct(FilesDF$ctime, origin = "1970-01-01 00:00.00 UTC")
      FilesDF$flsMdTm <- as.POSIXct(FilesDF$mtime, origin = "1970-01-01 00:00.00 UTC")
      FilesDF$TimeDiff <- difftime(FilesDF$flsMdTm, FilesDF$flsCrTm, units = "hours")
      #
      DestDir <- NULL
      while (is.null(DestDir)) { DestDir <- selectDirectory("Choose location where files will be cloned and/or the log(s) saved", path = "C:/") }
      #
      # Check available space
      # Adapted from https://stackoverflow.com/questions/32200879/how-to-get-disk-space-of-windows-machine-with-r
      disks <- system("wmic logicaldisk get size,freespace,caption", inter = TRUE)
      disks <- strsplit(disks[1:(length(disks)-1)], " +")
      disks <- disks[which(sapply(disks, length) == 4)]
      k <- unlist(disks[1])
      disks <- as.data.frame(t(sapply(disks[2:length(disks)], function(x) { x[1:3] })))
      colnames(disks) <- k[1:3]
      rownames(disks) <- NULL
      disks$FreeSpace <- as.numeric(as.character(disks$FreeSpace))
      disks$Size <- as.numeric(disks$Size)
      disks$FreeSpace[which(disks$Caption == "B:")] <- Inf # Hack because archive space is not being properly read and we will always have enough space there
      FrSp <- as.numeric(disks$FreeSpace[match(gsub("/.*", "", DestDir), disks$Caption)])
      #
      pat <- gsub("/+", "/", paste0(rootDir, "/", parDir, "/"))
      if (grepl("^[A-Z]:/", pat)) { pat <- paste0("^", pat) }
      FilesDF$New <- gsub(pat, paste0(DestDir, "/", parDir, "/"), FilesDF$Files)
      #
      # Test for path length
      # Windows has a default max path length of 260
      # Long paths remove the limitation
      mx <- max(nchar(FilesDF$New))
      abort <- FALSE
      if (mx > 260) { 
        abort <- dlg_message("Long paths (> 260 characters) predicted for destination files. Please check in the registry that HKEY_LOCAL_MACHINE/SYSTEM/CurrentControlSet/Control/FileSystem/LongPathsEnabled=1; If yes, we can continue...",
                             "yesno")$res
        abort <- abort != "yes"
      }
      if (abort) {
        warning("Illegal path length (> 260), set HKEY_LOCAL_MACHINE/SYSTEM/CurrentControlSet/Control/FileSystem/LongPathsEnabled to 1 and restart or choose a less deep target location...")
      } else {
        # Now, py_diAID is creating files with impossibly long names.
        # I can check their existence but cannot move or copy them programmatically with either:
        # - file.copy()
        # - fs::file_copy()
        # - copy in Windows command line
        # Instead, we will make an exception for those files rename the copy to a more sensible file name
        FilesDF$py_diAID <- grepl("/Optimization_plot_A1_ *[0-9]+\\.[0-9]+_A2_ *[0-9]+\\.[0-9]+_B1_ *[0-9]+\\.[0-9]+_B2_ *[0-9]+\\.[0-9]+_result_ *[0-9]+\\.[0-9]+\\.png",
                                  FilesDF$Files)
        if (sum(FilesDF$py_diAID)) {
          FilesDF$New[which(FilesDF$py_diAID)] <- parSapply(parClust, FilesDF$Files[which(FilesDF$py_diAID)], function(x) {
            #x <- FilesDF$Files[py_diAID[1]]
            x <- unlist(strsplit(x, "/"))
            dr <- paste(x[1:(length(x)-1)], collapse = "/")
            x <- rev(x)[1]
            x <- unlist(strsplit(gsub("^Optimization_plot_|\\.png", "", x), "_"))
            nms <- grep("^[A-B][1-2]$|^result$", x)
            vals <- round(as.numeric(x[nms[1:4]+1]), 3)
            res <- round(as.numeric(x[nms[5]+1]))
            x <- paste0(dr, "/Optimization_plot_A1_", vals[1], "_A2_ ", vals[2], "_B1_ ", vals[3], "_B2_ ", vals[4], "_result_", res, ".png")
            return(x)
          })
        }
        #
        FilesDF$Hash <- "Bleh!"
        FilesDF$NewHash <- "Blah!" # The two hashes' dummy values should be different!!!
        dstDirs <- unique(dirname(FilesDF$New))
        for (dstDir in dstDirs) { if (!dir.exists(dstDir)) { dir.create(dstDir, recursive = TRUE) } }
        FilesDF <- FilesDF[which(!is.na(FilesDF$size)),] # Exclude some system files if at root of drive
        FilesDF$Invoice <- grepl("\\.raw$|\\.d", FilesDF$Files, ignore.case = TRUE) # Files to invoice
        FilesDF$IsINI <- gsub(".*\\.", "", FilesDF$Files) == "ini"
        #
        #
        # Decide whether to invoice and if so, which files to use to calculate MS time
        # (e.g. exclude some blanks, reruns, failed injections etc...)
        #msgRoot <- paste0("Cloning from \n       ", gsub(".*/", ".../", TargDir), "\nto \n       ", gsub(".*/", ".../", DestDir), "\n")
        msgRoot <- paste0("Cloning from \n       ", TargDir, "\nto\n       ",
                          DestDir, "\n")
        logTm <- Sys.time()
        Invoicing <- FALSE
        FilesDF$ThermoRaw <- FilesDF$BrukerD <- FALSE
        wI <- which(FilesDF$Invoice)
        FilesDF$ThermoRaw[wI] <- grepl("\\.raw$", FilesDF$Files[wI], ignore.case = TRUE) # Thermo .raw files
        FilesDF$BrukerD[wI] <- grepl("[^/]\\.d", dirname(FilesDF$Files[wI]), ignore.case = TRUE) # Bruker .d folders
        isThermo <- sum(FilesDF$ThermoRaw) > 0
        isBruker <- sum(FilesDF$BrukerD) > 0
        if (isThermo) {
          # We will try to read instrument time from the raw file
          # Typically, MS spectra do not cover the whole chromatogram, but the pressure profile does.
          # So we will focus on that:
          # 1 - Download pressure extracted script
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
          pressScript2 <- paste0(DestDir, "/pressScript.py")
          write(tmp, pressScript2)
          # 2 - Convert files to mzML
          rawFiles <- FilesDF$Files[which(FilesDF$ThermoRaw)]
          mzMLs <- gsub("\\.raw$", ".mzML", rawFiles, ignore.case = TRUE)
          wDel <- which(!file.exists(mzMLs))
          pressFls <- gsub("\\.mzML$", ".csv", mzMLs)
          w <- which(!file.exists(mzMLs))
          if (length(w)) {
            # Convert to mzML
            deer <- list()
            ParsDirs <- grep("/ThermoRawFileParser", c(list.dirs("C:/Program Files", full.names = TRUE, recursive = FALSE),
                                                       list.dirs(paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads"), full.names = TRUE, recursive = FALSE)),
                             value = TRUE)
            tst <- sapply(ParsDirs, function(x) { "ThermoRawFileParser.exe" %in% list.files(x) })
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
            clusterExport(parClust, list("rawFiles", "mzMLs", "pressFls", "w"), envir = environment())
            cat("Temporarily converting", length(w), "Thermo raw files to mzML to extract run time...\n")
            # (Unfortunately, it does not seem to work if I include only the pressure chromatograms...)
            if (tolower(Convert_mode) == "thermorawfileparser") { # Mode 1: using ThermoRawFileParser
              clusterExport(parClust, list("DestDir", "deer", "zlib", "PeakPicking"), envir = environment())
              tst <- parSapply(parClust, w, function(i) { #i <- w[1]
                cmd <- paste0("\"", deer$ParsDir, "/ThermoRawFileParser.exe\" -i=\"",
                              rawFiles[i], "\" -b=\"", gsub(".*/", paste0(DestDir, "/"), mzMLs[i]), "\" -f=2 -a",
                              c(" -z", "")[zlib+1], c(" -p", "")[PeakPicking+1]#, " -L=1"
                              )
                write(cmd, paste0(DestDir, "/tmp", i, ".bat"))
                cmd2 <- paste0("\"", DestDir, "/tmp", i, ".bat\"") # I have to go through an intermediate batch file to run cmd,
                # because it is one of those which works in Windows command line but not when passed to system() or shell() in R.
                #cat(cmd)
                #writeClipboard(cmd)
                system(cmd2)
                unlink(paste0(DestDir, "/tmp", i, ".bat"))
              })
            }
            if (tolower(Convert_mode) == "msconvert") { # Mode 2: using msconvert
              write(rawFiles[w], file = paste0(DestDir, "/tmp_MS_files.txt"))
              precRecal <- FALSE
              cmd <- paste0("\"", deer$MSConvertDir, "/msconvert.exe\" -f \"", DestDir, "/tmp_MS_files.txt\" -o \"",
                            DestDir, "\" --mzML --64 -z --filter \"peakPicking true 1-\""#, " --filter \"msLevel 1-1\""
                            )
              #cat(cmd)
              system(cmd)
              unlink(paste0(DestDir, "/tmp_MS_files.txt"))
            }
          }
          # 3 - get pressure traces
          clusterExport(parClust, list("pressScript2", "mzMLs"), envir = environment())
          tst <- parSapply(parClust, w, function(i) { #i <- w[1]
            cmd <- paste0("python \"", pressScript2, "\" ", paste0("\"", mzMLs[i], "\"", collapse = " "))
            #cat(cmd)
            system(cmd)
          })
          
          # 4 - read run times from pressure trace

          # 5 - Cleanup
          unlink(mzMLs[wDel])
          
        }
        if (isBruker) {
          dDrs <- unique(grep("[^/]+\\.d$", dirname(FilesDF$Files[which(FilesDF$BrukerD)]), value = TRUE))
          tst <- parSapply(parClust, dDrs, function(dr) {
            length(grep("/[0-9]+\\.m$", list.dirs(dr))) == 1
            # If not specific enough, add more tests:
            # Our d folders typically always contain the same files, but which specifically?
            # Can we trust that these will always be there, or are some only present with some options?
          })
          dDrs <- dDrs[which(tst)]
        }
        isMS <- isThermo+isBruker
        if (isMS == 2) {
          warning("This script detected both Thermo .raw files and Bruker .d folders in this folder. Check script and only proceed with the next steps with caution if you know what you are doing!")
        }
        justSimul <- FALSE
        if (isMS) {
          modes <- c("Copy data to the destination folder                                                                    ",
                     "Run invoicing part only                                                                                ")
          justSimul <- c(FALSE, TRUE)[match(dlg_list(modes, modes[1], title = "MS acquisition folder detected. What do you want to do?")$res, modes)]
          ## Thermo case
          if (isThermo) {
            msg <- paste0(msgRoot, "\n-> Write a log of which Thermo .raw files will be invoiced?\n\n")
            Invoicing <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res,
                                              c("yes", "no"))]
            if (Invoicing) {
              wRw <- which(FilesDF$ThermoRaw)
              rwFls <- FilesDF[wRw,]
              rwFls <- rwFls[order(rwFls$flsCrTm, decreasing = FALSE),]
              rwFlsSz <- as.character(rwFls$size)
              tmp <- nchar(rwFls$Files)+nchar(rwFlsSz)
              clusterExport(parClust, "tmp", env = environment())
              tmp <- parSapply(parClust, tmp, function(x) {
                paste(rep(" ", (max(tmp) - x + 3)*2), collapse = "")
              })
              rwFlstmp <- paste0(rwFls$Files, tmp, rwFlsSz, " bytes")
              msg <- paste0(msgRoot,
                            "\n-> Select Thermo .raw files to include in calculations of useful instrument time (to charge):")
              Fls2Chrg <- dlg_list(rwFlstmp, rwFlstmp, TRUE, msg)$res
              Fls2Chrg <- rwFls[match(Fls2Chrg, rwFlstmp),]
              m <- match(Fls2Chrg$Files, FilesDF$Files)
              TimeDiff <- difftime(FilesDF$flsMdTm[m], FilesDF$flsCrTm[m], units = "hours")
              # In that case, you do not want to use it, because "file created date" will now be lated than file modified date!
              Log1 <- c("THERMO INVOICING LOG",
                        "-----------------------",
                        "",
                        paste0("Folder: ", TargDir),
                        paste0(logTm, ": ", Whodunnit, " ran this script. The following files were selected for invoicing:"),
                        "",
                        paste0(" - ", Fls2Chrg$Files),
                        "",
                        paste0("Selected instrument run time to charge: ", ceiling(sum(TimeDiff)*10)/10, " h"),
                        "")
              Log1Fl <- paste0(TargDir, "/Thermo invoicing log_", gsub("[^0-9]", "", logTm), ".txt")
              write(Log1, Log1Fl)
              file.copy(Log1Fl, paste0(gsub("/$", "", DestDir), "/", parDir))
              system(paste0("open \"", Log1Fl, "\""))
            }
          }
          ## Bruker case
          if (isBruker) {
            msg <- paste0(msgRoot, "\n-> Write a log of which Bruker .d folders will be invoiced?\n\n")
            Invoicing <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res,
                                              c("yes", "no"))]
            if (Invoicing) {
              back2basics <- FALSE
              dCrTm <- parSapply(parClust, dDrs, function(dFl) { #dFl <- dDrs[1]
                fls <- list.files(dFl, full.names = TRUE)
                min(file.info(fls)$ctime, na.rm = TRUE)
              })
              dDrsTbl <- data.frame(d = dDrs, Created = as.POSIXct(dCrTm, origin = "1970-01-01 00:00.00 UTC"))
              dDrsTbl <- dDrsTbl[order(dDrsTbl$Created, decreasing = FALSE),]
              dDrs <- dDrsTbl$d
              pat <- topattern(gsub("/+", "/", paste0(rootDir, "/", parDir, "/")))
              tst <- lapply(dDrs, function(dFl) { #dFl <- dDrs[1]
                sbdr <- grep("\\.m$", list.dirs(dFl), invert = TRUE, value = TRUE)
                sbdr <- sbdr[which(sbdr != dFl)]
                lgFl <- grep("execution-log", list.files(sbdr, full.names = TRUE), value = TRUE)
                if (length(lgFl) != 1) {
                  if (length(lgFl) > 1) { res <- list(outcome = FALSE) } else {
                    nuFl <- gsub(pat, paste0(DestDir, "/", parDir, "/"), dFl)
                    sbdr <- grep("\\.m$", list.dirs(nuFl), invert = TRUE, value = TRUE)
                    sbdr <- sbdr[which(sbdr != nuFl)]
                    lgFl <- grep("execution-log", list.files(sbdr, full.names = TRUE), value = TRUE)
                    if (length(lgFl) != 1) { res <- list(outcome = FALSE) } else{
                      res <- list(outcome = TRUE, file = lgFl)
                    }
                  }
                } else {
                  res <- list(outcome = TRUE, file = lgFl)
                }
                return(res)
              })
              w <- which(sapply(tst, function(x) { x$outcome }))
              goOn <- (length(w) > 0)
              if (goOn) {
                dDrs <- dDrs[w]; dDrsTbl <- dDrsTbl[w,]; tst <- tst[w]
                dDrsTbl$Log_file <- sapply(tst, function(x) { x$file })
                clusterExport(parClust, list("pat", "DestDir", "parDir"), envir = environment())
                dDrsTbl$"Length (h)" <- parApply(parClust, dDrsTbl[, c("d", "Log_file")], 1, function(x) { #x <- dDrsTbl[1, c("d", "Log_file")]
                  lgFl <- x[[2]]
                  lg <- readLines(lgFl)
                  #system(paste0("open \"", lgFl, "\""))
                  stp <- grep("COMPLETED", lg, value = TRUE)
                  if (!length(stp)) { stp <- grep("ABORTED", lg, value = TRUE) }
                  if (!length(stp)) { res <- list(outcome = FALSE) } else {
                    stp <- unlist(strsplit(gsub("\t.*", "", stp), " - "))[2]
                    if (!grep("^000\\.", stp)) { stop("This doesn't (yet) support methods this long!") }
                    tm <- unlist(strsplit(gsub("^000\\.", "", stp), ":"))
                    res <- list(outcome = TRUE, value = sum(as.numeric(tm) * c(1, 1/60, 1/3600)))
                  }
                  return(res)
                })
                w <- which(sapply(dDrsTbl$"Length (h)", function(x) { x$outcome }))
                dDrs <- dDrs[w]; dDrsTbl <- dDrsTbl[w,]
                dDrsTbl$"Length (h)" <- sapply(dDrsTbl$"Length (h)", function(x) { x$value })
              } else {
                msg <- "Looks like this was not run with a Bruker LC, is this correct?"
                back2basics <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
              }
              if (back2basics) {
                goOn <- TRUE
                dDrsTbl$Files <- parLapply(parClust, dDrsTbl$d, list.files, full.names = TRUE)
                dDrsTbl$Created <- parSapply(parClust, dDrsTbl$Files, function(x) { min(file.info(x)$ctime) })
                dDrsTbl$Closed <- parSapply(parClust, dDrsTbl$Files, function(x) { max(file.info(x)$mtime) })
                dDrsTbl$"Length (h)" <- (dDrsTbl$Closed - dDrsTbl$Created)/3600
              }
              if (goOn) {
                dDrsTbl <- dDrsTbl[order(dDrsTbl$Created),]
                tmp <- nchar(dDrs) + nchar(dDrsTbl$Created)
                clusterExport(parClust, "tmp", env = environment())
                tmp <- parSapply(parClust, tmp, function(x) {
                  paste(rep(" ", (max(tmp) - x + 3)*2), collapse = "")
                })
                dDrstmp <- paste0(dDrs, tmp, dDrsTbl$Created)
                msg <- paste0(msgRoot,
                              "\n-> Select Bruker .d folders to include in calculations of useful instrument time (to charge):")
                Fls2Chrg <- dlg_list(dDrstmp, dDrstmp, TRUE, msg)$res
                Fls2Chrg <- dDrs[match(Fls2Chrg, dDrstmp)]
                #View(cbind(dDrsCrTm, dDrsMdTm, dDrsMdTm-dDrsCrTm))
                #View(cbind(dDrsCrTmPos, dDrsMdTmPos, dDrsMdTmPos-dDrsCrTmPos))
                TimeDiff <- sum(dDrsTbl$`Length (h)`[match(Fls2Chrg, dDrsTbl$d)])
                stopifnot(min(TimeDiff) >= 0) # This would indicate that some of these files are duplicates made without using the script
                # In that case, you do not want to use it, because "file created date" will now be later than file modified date!
                Log1 <- c("BRUKER INVOICING LOG",
                          "-----------------------",
                          "",
                          paste0("Folder: ", TargDir),
                          paste0(logTm, ": ", Whodunnit, " ran this script. The following files were selected for invoicing:"),
                          "",
                          paste0(" - ", Fls2Chrg),
                          "",
                          paste0("Selected instrument run time to charge: ", ceiling(sum(TimeDiff)*10)/10, " h"),
                          "")
                Log1Fl <- paste0(TargDir, "/Bruker invoicing log_", gsub("[^0-9]", "", logTm), ".txt")
                write(Log1, Log1Fl)
                file.copy(Log1Fl, paste0(gsub("/$", "", DestDir), "/", parDir))
                system(paste0("open \"", Log1Fl, "\""))
              } else { warning("Nothing to invoice: there are no valid/complete Bruker .d folders in the target folder.") }
            }
          }
        }
        if (!justSimul) {
          FilesDF$Copy <- TRUE
          FilesDF$NewSize <- NA
          FilesDF$NewExists <- file.exists(FilesDF$New)
          if (sum(FilesDF$NewExists)) {
            wExst <- which(FilesDF$NewExists)
            FilesDF$NewSize[wExst] <- file.size(FilesDF$New[wExst])
            # Calculate hash for overlapping files
            cat("Some matching files found in the target directory,\nchecking hashes for overlapping files:
(hang in there, this will take a while...)\n")
            cat(" - Original files...\n")
            FilesDF$Hash[wExst] <- parSapply(parClust, FilesDF$Files[wExst], HashFun)
            w <- which(is.na(FilesDF$Hash[wExst]))
            if (length(w)) {
              FilesDF$Hash[wExst[w]] <- parSapply(parClust, FilesDF$Files[wExst[w]], HashFun)
            }
            w <- which(is.na(FilesDF$Hash[wExst]))
            if (length(w)) {
              a <- (length(w) > 1)+1
              msg <- paste0("Md5sum calculations failed for the following file", c("", "s")[a], " after 2 attempts:\n",
                            paste(paste0(" - ", FilesDF$Files[wExst[w]]), collapse = "\n"), "\nCheck th", c("is", "ose")[a], " file", c("", "s")[a], ", ",
                            c("it", "they")[a], " may be corrupt.")
              #cat(msg)
              warning(msg)
            }
            cat(" - Destination files...\n")
            FilesDF$NewHash[wExst] <- parSapply(parClust, FilesDF$New[wExst], HashFun)
            BdHsh <- (is.na(FilesDF$NewHash[wExst]))|(FilesDF$Hash[wExst] != FilesDF$NewHash[wExst])
            BdHsh[which(is.na(BdHsh))] <- TRUE
            tstBdHsh <- (sum(BdHsh)>0)+1
            cat(paste0(" -> Result: ", c("no ", "")[tstBdHsh], "discrepancies found, ", c("ignor", "overwrit")[tstBdHsh], "ing overlapping files.\n"))
            SzOk <- FilesDF$NewSize[wExst] == FilesDF$size[wExst]
            FilesDF$Copy[wExst] <- (!SzOk)|(BdHsh)
            cat("\n")
          }
          kount <- 0 # Attempts
          CopyFls <- sum(FilesDF$Copy) > 0
          FilesDF$Copied <- FALSE
          #cat(msgRoot)
          wC <- which(FilesDF$Copy) # Initiate
          while ((CopyFls)&&(kount < MaxAttempts)) {
            kount <- kount+1
            cat(paste0("Copying attempt #", kount, "\n"))
            wC <- wC[which((!FilesDF$Copied[wC])
                           |(!FilesDF$NewExists[wC])
                           |((FilesDF$NewHash[wC] == "Blah!")|(is.na(FilesDF$NewHash[wC])))
                           |((!is.na(FilesDF$NewSize[wC]))&(FilesDF$size[wC] != FilesDF$NewSize[wC]))
                           |(is.na(FilesDF$NewHash[wC]))
                           |(is.na(FilesDF$Hash[wC]))
                           |((!is.na(FilesDF$NewHash[wC]))&(!is.na(FilesDF$Hash[wC]))&(FilesDF$Hash[wC] != FilesDF$NewHash[wC])))]
            # Check for space
            if (sum(FilesDF$size[wC], na.rm = TRUE) > FrSp) { 
              msg <- "Insufficient free space in target volume! Copy anyway? (only do this for network resources which may not report free space correctly)"
              CopyFls <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
            }
            if (CopyFls) {
              # Calculate original hash
              wBlh <- wC[which(FilesDF$Hash[wC] == "Bleh!")]
              if (length(wBlh)) {
                cat(" - Calculating original file hashes...\n")
                FilesDF$Hash[wBlh] <- parSapply(parClust, FilesDF$Files[wBlh], HashFun)
                w <- which(is.na(FilesDF$Hash[wBlh]))
                if (length(w)) {
                  FilesDF$Hash[wBlh[w]] <- parSapply(parClust, FilesDF$Files[wBlh[w]], HashFun)
                }
                w <- which(is.na(FilesDF$Hash[wBlh]))
                if (length(w)) {
                  a <- (length(w) > 1)+1
                  msg <- paste0("Md5sum calculations failed for the following file", c("", "s")[a], " after 2 attempts:\n",
                                paste(paste0(" - ", FilesDF$Files[wExst[w]]), collapse = "\n"), "\nCheck th", c("is", "ose")[a], " file", c("", "s")[a], ", ",
                                c("it", "they")[a], " may be corrupt.")
                  #cat(msg)
                  warning(msg)
                }
              }
              cat(" - Copying files...\n")
              # Only export what we need
              Fls_wC <- FilesDF$Files[wC]
              nuFls_wC <- FilesDF$New[wC]
              clusterExport(parClust, list("Fls_wC", "nuFls_wC"), env = environment())
              tst <- try(parSapply(parClust, 1:length(wC), function(x) {
                fs::file_copy(Fls_wC[x], nuFls_wC[x], overwrite = TRUE)
              }), silent = TRUE)
              FilesDF$NewExists[wC] <- file.exists(FilesDF$New[wC])
              wY <- which(FilesDF$NewExists[wC]) # Indices of wC, exceptionally and for practical reasons
              wN <- which(!FilesDF$NewExists[wC]) # Indices of wC, exceptionally and for practical reasons
              FilesDF$Copied[wC[wY]] <- TRUE
              # Potential issue and its fix:
              # We are copying in one call a lot of files. In some cases, it has occurred that this would fail on  a large (spectra binaries) file,
              # or for files with absurdly long paths.
              # Solution: try copying this file individually.
              if (length(wN)) {
                TRYFIX <- TRUE
                while ((length(wN))&&(TRYFIX)) {
                  tst <- try(file.copy(Fls_wC[wN[1]], nuFls_wC[wN[1]], overwrite = TRUE, recursive = TRUE), silent = TRUE)
                  if (!"try-error" %in% class(tst)) {
                    FilesDF$NewExists[wC[wN[1]]] <- file.exists(FilesDF$New[wC[wN[1]]])
                    wN <- wN[which(wN != wN[1])]
                    if (length(wN)) {
                      tst <- try(parSapply(parClust, wN, function(x) {
                        fs::file_copy(Fls_wC[wC[x]], nuFls_wC[wC[x]], overwrite = TRUE)
                      }), silent = TRUE)
                      FilesDF$NewExists[wC[wN]] <- file.exists(FilesDF$New[wC[wN]])
                      wY <- which(FilesDF$NewExists[wC]) # Indices of wC
                      wN <- which(!FilesDF$NewExists[wC]) # Indices of wC
                      FilesDF$Copied[wC[wY]] <- TRUE
                    }
                  } else { TRYFIX <- FALSE }
                }
              }
              # Check success based on whether:
              # - the file exists
              # - its size is correct
              # - its hash is correct
              wExst <- wC[which(FilesDF$NewExists[wC])] # Normally wExst should be identical to wC
              FilesDF$NewSize[wExst] <- file.size(FilesDF$New[wExst])
              cat(" - Calculating cloned file hashes...\n")
              FilesDF$NewHash[wExst] <- parSapply(parClust, FilesDF$New[wExst], HashFun)
              BdHsh <- (is.na(FilesDF$NewHash[wExst]))|(FilesDF$Hash[wExst] != FilesDF$NewHash[wExst])
              BdHsh[which(is.na(BdHsh))] <- TRUE
              tstBdHsh <- (sum(BdHsh)>0)+1
              if (tstBdHsh == 2) {
                cat("   checking discrepancies...\n")
                wBdHsh <- wExst[which(BdHsh)]
                FilesDF$size[wBdHsh] <- file.size(FilesDF$Files[wBdHsh])
                FilesDF$Hash[wBdHsh] <- parSapply(parClust, FilesDF$Files[wBdHsh], HashFun)
                BdHsh <- (is.na(FilesDF$NewHash[wExst]))|(FilesDF$Hash[wExst] != FilesDF$NewHash[wExst])
                BdHsh[which(is.na(BdHsh))] <- TRUE
                tstBdHsh <- (sum(BdHsh)>0)+1
                wC <- unique(c(wC[which(!wC %in% wExst)],
                               wExst[which((FilesDF$size[wExst] != FilesDF$NewSize[wExst])
                                           |BdHsh)]))
                if (tstBdHsh == 2) {
                  cat(paste0(" -> Result: ", sum(BdHsh), " confirmed issue, cloning attempt #", kount,
                             " was unsuccessful.\n   Is it possible that the input files are actually still being modified?\n"))
                } else {
                  CopyFls <- FALSE
                  cat(paste0(" -> Result: no discrepancies found, cloning was successful", c("", paste0(" after ", kount, " attempts"))[(kount > 1)+1],
                             ".\nHowever, original file hashes needed re-calculating, indicating that they were modified since the script was starting. Are you exploring the files (e.g. in Bruker's Data Analysis software)? If not, check again whether data acquisition is not still ongoing in the folder!\n"))
                }
              } else {
                wC <- unique(c(wC[which(!wC %in% wExst)],
                               wExst[which((FilesDF$size[wExst] != FilesDF$NewSize[wExst])
                                           |BdHsh)]))
                CopyFls <- FALSE
                cat(paste0(" -> Result: no discrepancies found, cloning was successful", c("", paste0(" after ", kount, " attempts"))[(kount > 1)+1], ".\n"))
              }
            }
          }
          tsts <- c(length(wC) == 0, # Do we still have files to copy?
                    sum(!file.exists(FilesDF$New)) == 0, # Do all of the "new" (= destination) files exist?
                    sum(FilesDF$size != FilesDF$NewSize) == 0, # Do they have the expected size?
                    sum(is.na(FilesDF$NewHash)) == 0, # Do we have a hash for all of those new files?
                    sum(FilesDF$Hash != FilesDF$NewHash) == 0) # And are those valid hashes?
          tsts[which(is.na(tsts))] <- FALSE # Fix NAs: they are failures!
          if (!sum(!tsts)) { # Are all tests successes?
            FilesDF$nuFlsCrTm <- parSapply(parClust, FilesDF$New, function(fl) { file.info(fl)$ctime })
            FilesDF$nuFlsMdTm <- parSapply(parClust, FilesDF$New, function(fl) { file.info(fl)$mtime })
            # Adjust created time (modified time is normally good)
            # This is currently turned off by default, for 2 reasons:
            # - Data is currently mostly being transferred to the archive... where we do not want to modify it once written. Does this count as a modification?
            # - Invoicing (for which we are using instrument time) - is now being done on original files in most cases.
            if (adjstTimes) {
              wJ <- which((!FilesDF$IsINI)&(FilesDF$Copy))
              if (length(wJ)) {
                # Get time created and modified
                wCr <- which(FilesDF$nuFlsCrTm[wJ] != FilesDF$flsCrTm[wJ]) # Indices of wJ
                wMd <- which(FilesDF$nuFlsMdTm[wJ] != FilesDF$flsMdTm[wJ]) # Indices of wJ
                if (length(c(wCr, wMd))) {
                  nuFls_wJ <- FilesDF$New[wJ]
                  flsCrTm_wJ <- FilesDF$flsCrTm[wJ]
                  flsMdTm_wJ <- FilesDF$flsMdTm[wJ]
                  clusterExport(parClust, list("nuFls_wJ", "nrmPath4PS", "flsCrTm_wJ", "flsMdTm_wJ"), env = environment())
                  if (length(wCr)) {
                    cat("Adjusting file created times...\n")
                    parSapply(parClust, wCr, function(i) {
                      nufl <- nuFls_wJ[i]
                      flnm4pwrshll <- nrmPath4PS(nufl)
                      tm <- flsCrTm_wJ[i]
                      tmY <- as.integer(format(tm, "%Y"))
                      tmM <- as.integer(format(tm, "%m"))
                      tmD <- as.integer(format(tm, "%d"))
                      tmHr <- as.integer(format(tm, "%H"))
                      tmMn <- as.integer(format(tm, "%M"))
                      tmSc <- as.integer(format(tm, "%S"))
                      cmd <- paste0('powershell -command "$(Get-Item \"', flnm4pwrshll, '\").CreationTime=$(Get-Date -Year ', tmY,
                                    ' -Month ', tmM, ' -Day ', tmD, ' -Hour ', tmHr, ' -Minute ', tmMn, ' -Second ', tmSc, ')"')
                      #cat(cmd)
                      resCr <- system(cmd, intern = TRUE)
                    })
                  }
                  if (length(wMd)) {
                    cat("Adjusting file modified times...\n")
                    parSapply(parClust, wCr, function(i) {
                      nufl <- nuFls_wJ[i]
                      flnm4pwrshll <- nrmPath4PS(nufl)
                      tm <- flsMdTm_wJ[i]
                      tmY <- as.integer(format(tm, "%Y"))
                      tmM <- as.integer(format(tm, "%m"))
                      tmD <- as.integer(format(tm, "%d"))
                      tmHr <- as.integer(format(tm, "%H"))
                      tmMn <- as.integer(format(tm, "%M"))
                      tmSc <- as.integer(format(tm, "%S"))
                      cmd <- paste0('powershell -command "$(Get-Item \"', flnm4pwrshll, '\").LastWriteTime=$(Get-Date -Year ', tmY,
                                    ' -Month ', tmM, ' -Day ', tmD, ' -Hour ', tmHr, ' -Minute ', tmMn, ' -Second ', tmSc, ')"')
                      #cat(cmd)
                      resMd <- system(cmd, intern = TRUE)
                    })
                  } 
                }
              } 
            }
            # Logs
            if (isMS) {
              # - Optionally remove local MS files (it will leave other ones behind!)
              Remove <- c(TRUE, FALSE)[match(dlg_message(paste0("Do you want to delete local MS files?\n(", msgRoot, ")"),
                                                         "yesno")$res, c("yes", "no"))]
              if (Remove) {
                openwd(TargDir)
                openwd(paste0(gsub("/$", "", DestDir), "/", parDir))
                # Filter to only keep files to remove
                rawFls2Rmv <- dDrs2Rmv <- c()
                if (isThermo) {
                  rawFls2Rmv <- FilesDF$Files[which((FilesDF$NewExists)
                                                    &(FilesDF$size == FilesDF$NewSize)
                                                    &(FilesDF$NewHash == FilesDF$Hash)
                                                    &(FilesDF$ThermoRaw))]
                }
                if (isBruker) {
                  dDrs2Rmv <- FilesDF$Files[which((FilesDF$NewExists)
                                                  &(FilesDF$size == FilesDF$NewSize)
                                                  &(FilesDF$NewHash == FilesDF$Hash)
                                                  &(FilesDF$BrukerD))]
                  dDrs2Rmv <- unique(grep("[^/]+\\.d$", dirname(dDrs2Rmv), value = TRUE))
                }
                Fls2Rmv <- c(rawFls2Rmv, dDrs2Rmv)
                if (length(Fls2Rmv)) {
                  RmvTsts <- setNames(rep(FALSE, length(Fls2Rmv)), Fls2Rmv)
                  if (isThermo) {
                    for (fl in rawFls2Rmv) { #fl <- rawFls2Rmv[1]
                      tst <- try(unlink(fl, recursive = TRUE), silent = TRUE)
                      RmvTsts[fl] <- (!"try-error" %in% class(tst))
                    }
                  }
                  if (isBruker) {
                    for (dr in dDrs2Rmv) { #dr <- dDrs2Rmv[1]
                      # I had to use Powershell to remove the folders!!!
                      cmd <- paste0('powershell -command "$(Remove-Item -Path \"', nrmPath4PS(dr), '\" -Recurse -Force )"')
                      #cat(cmd, "\n")
                      #writeClipboard(cmd)
                      system(cmd)
                      RmvTsts[dr] <- !dir.exists(dr)
                    }
                  }
                  Log2 <- c("MS FILES DELETION LOG",
                            "-----------------------",
                            "",
                            paste0("Folder: ", TargDir),
                            paste0(logTm, ": ", Whodunnit, " authorized the deletion, following archiving, of the following local MS files:"),
                            "",
                            paste0(" - ", Fls2Rmv),
                            "")
                  if (sum(!RmvTsts)) {
                    Log2 <- c(Log2, "However, deletion failed for the following files:", paste0(" - ", Fls2Rmv[which(!RmvTsts)]), "")
                  }
                  if ((Invoicing)&&(length(Fls2Rmv))) {
                    Log2 <- c(Log2, paste0("Selected instrument run time to charge: ", ceiling(sum(TimeDiff)*10)/10, " h"), "")
                  }
                  Log2Fl <- paste0(TargDir, "/MS files deletion log_", gsub("[^0-9]", "", logTm), ".txt")
                  write(Log2, Log2Fl)
                  file.copy(Log2Fl, paste0(gsub("/$", "", DestDir), "/", parDir))
                  system(paste0("open \"", Log2Fl, "\""))
                } else { message("No files were deleted.") }
              }
            } 
          } else {
            if (!CopyFls) {
              msg <- paste0(msgRoot,
                            "\nCloning the files failed, investigate (not enough space? stable connection? bug? ...)\n")
              warning(msg)
            } else {
              msg <- paste0(msgRoot,
                            "\nCloning the files failed after ", MaxAttempts,
                            " attempts, investigate!\n")
              warning(msg)
            }
          }
        }
      }
    }
  }
} else { cat("No target directory selected...\nGOODBYE!\n") }

# NB:
# Sometimes when running the script several times, the cluster gets corrupted. It can be a good idea to stop the cluster in between by running:
# stopCluster(parClust) # (This may become the default behaviour: re-creating the cluster from scratch takes maybe 10s, not a huge deal...)
