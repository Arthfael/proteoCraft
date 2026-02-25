# Move one QC file
options(stringsAsFactors = FALSE) # Needed because we will stick to an older version

# Default behavior:
# - We only want to load files which are younger than this time
TimeWindow <- 30 # In days
# - We only want to load files from the main folder, not the test folders (method development)
Rec <- FALSE
# We only want to move one file every 6h
Every <- 0.01 # in hours

topattern <- function (x, start = TRUE, end = FALSE, collapse = "|") {
  x <- gsub("\\.", "\\\\.", as.character(x))
  x <- gsub("\\*", "\\\\*", as.character(x))
  x <- gsub("\\$", "\\\\$", as.character(x))
  x <- gsub("\\^", "\\\\^", as.character(x))
  x <- gsub("\\+", "\\\\+", as.character(x))
  x <- gsub("\\?", "\\\\?", as.character(x))
  x <- gsub("\\{", "\\\\{", as.character(x))
  x <- gsub("\\}", "\\\\}", as.character(x))
  x <- gsub("\\[", "\\\\[", as.character(x))
  x <- gsub("\\]", "\\\\]", as.character(x))
  x <- gsub("\\(", "\\\\(", as.character(x))
  x <- gsub("\\)", "\\\\)", as.character(x))
  x <- gsub("\\|", "\\\\|", as.character(x))
  if (start) { x <- paste0("^", x) }
  if (end) { x <- paste0(x, "$") }
  if ((length(x) > 1) && (collapse != FALSE)) {
    x <- paste(x, collapse = collapse)
  }
  return(x)
}

LocQCDir <- "D:/Data/Standard runs"
TargQCDir <- "...QC_File_Archive/QCloud"
dir.exists(TargQCDir)
LogFl <- paste0(TargQCDir, "/R_Script_file_transfer_log.txt")
if (dir.exists(TargQCDir)) {
  # Detect raw files in QC directory
  fls <- paste0(LocQCDir, "/", list.files(LocQCDir, recursive = Rec))
  raws <- grep("\\.raw$", fls, value = TRUE, ignore.case = TRUE)
  if (length(raws)) {
    QC1 <- grep("qc0?1", raws, ignore.case = TRUE, value = TRUE)
    QC2 <- grep("qc0?2", raws, ignore.case = TRUE, value = TRUE)
    QC3 <- grep("qc0?3", raws, ignore.case = TRUE, value = TRUE)
    if (length(c(QC1, QC2))) {
      QCs <- data.frame(Path = c(QC1, QC2),
                        check.names = FALSE)
      QCs$`Raw file` <- basename(QCs$Path)
      QCs$Dir <- dirname(QCs$Path)
      # Create new name with time stamp (code re-used for archiving script)
      QCs$"Modified time" <- sapply(QCs$Path, function(x) { as.character(file.info(x)$mtime) })
      QCs$"Creation time" <- sapply(QCs$Path, function(x) { as.character(file.info(x)$ctime) })
      QCs[, c("Year", "Month", "Day", "Date")] <- as.data.frame(t(sapply(strsplit(QCs$"Creation time", "-| +"),
                                                                         function(x) {
                                                                           c(x[[1]], x[[2]], x[[3]], paste0(paste(c(x[[1]], x[[2]], x[[3]]), collapse = "-"), " ", x[[4]], " CEST"))
                                                                         })))
      for (i in c("Year", "Month", "Day")) { QCs[[i]] <- as.numeric(QCs[[i]]) }
      QCs$Date <- as.POSIXct(QCs$Date, origin = "1970-01-01 00:00.00 UTC")
      QCs$`New name` <- QCs$`Raw file`
      QCs$Extension <- sapply(strsplit(QCs$`Raw file`, "\\."), function(x) { rev(unlist(x))[1] }) #To accept all combinations of upper/lower case
      QCs$Size <- file.size(QCs$Path)
      #
      # Deal with time-stamps
      # - 1) Files with a time stamp
      g1 <- grep("_[0-9]{14}\\.raw$", QCs$`Raw file`, ignore.case = TRUE)
      #View(Raw[g1,])
      if (length(g1)) {
        QCs$`New name`[g1] <- paste0(gsub("_[0-9]+\\.raw$", "", QCs$`Raw file`[g1], ignore.case = TRUE),
                                     "_", gsub("[^0-9]", "", QCs$Date[g1]), ".", QCs$Extension[g1])
        w <- which(QCs$`New name`[g1] != QCs$`Raw file`[g1])
        if (length(w)) {
          #View(QCs[g1[w], c("Raw file", "New name")])
          warning("Some predicted time stamps do not match the ones created by XCalibur, our assumptions as to how it created them may not hold! Investigate!")
          # For 1 project, there was a small difference of 1s, it is ok to stick to the old file name then
          QCs$`New name`[g1][w] <- QCs$`Raw file`[g1][w]
        }
      }
      # - 2) Files without a time stamp
      g2 <- c(1:nrow(QCs))[which(!c(1:nrow(QCs)) %in% g1)]
      if (length(g2)) {
        QCs$`New name`[g2] <- paste0(gsub("\\.raw$", "", QCs$`Raw file`[g2], ignore.case = TRUE),
                                     "_", gsub("[^0-9]", "", QCs$Date[g2]), ".", QCs$Extension[g2])
      }
      #
      QCs$`New path` <- apply(QCs[, c("Dir", "New name")], 1, paste, collapse = "/")
      # Are the files old and large enough?
      QCs$Date_created <- sapply(QCs$Path, function(x) { as.numeric(file.info(x)$ctime) })
      QCs$Copy <- (((as.numeric(Sys.time()) - QCs$Date_created) <= TimeWindow*24*60*60)&
                     (QCs$Size > 1024*50))
      
      #Niklas: we cannot us sum(copy), right? as QCloud only processes one file, it uploads all, 
      # but processes just one. we had that issue in the past. do you know more? So just copy, if there is no file.
      # and still the question if the qcloud workflow is working, if one file is uploaded but still processing, 
      #and you copy and upload the next file
      
      if (sum(QCs$Copy)) { # Is there at least one to copy?
        # Are the files still being acquired?
        Sys.sleep(5)
        whC <- which(QCs$Copy)
        QCs$Copy[whC] <- QCs$Size[whC] == file.size(QCs$Path[whC])
        whC <- which(QCs$Copy)
        if (length(whC)) {
          # Are the files already present in one of the processed folders?
          QCs2 <- grep("\\.raw$", list.files(paste0(TargQCDir, c("", "/processed")), recursive = TRUE, full.names = TRUE),
                       value = TRUE, ignore.case = TRUE)
          if (length(QCs2)) {
            QCs2 <- data.frame(Path = QCs2, check.names = FALSE)
            QCs2$Size <- file.size(QCs2$Path)
            QCs2$"Raw file" <- basename(QCs2$Path)
            QCs2$Dir <- dirname(QCs2$Path)
            wRdt <- whC[which(QCs$"New name"[whC] %in% QCs2$`Raw file`)] # Redundant files
            # Either they are not identical and we need to investigate (since there should be a time stamp, all conflicts  should be identical),
            # or they are identical and we just remove the local one.
            if (length(wRdt)) {
              tst <- apply(QCs[wRdt, c("New name", "Size")], 1, function(x) {
                paste(c(x[[1]], as.character(x[[2]])), collapse = "___")
              })
              tst2 <- apply(QCs2[match(QCs$"New name"[wRdt], QCs2$`Raw file`),
                                 c("Raw file", "Size")], 1, function(x) {
                paste(c(x[[1]], as.character(x[[2]])), collapse = "___")
              })
              wY <- which(tst %in% tst2)
              wN <- which(!tst %in% tst2)
              if (length(wY)) {
                # If a file is already there, we need to delete it locally to make space
                # Check hash
                HashFun <- function(x) { as.character(tools::md5sum(x)) } #tools::md5sum takes a character path, not a connection
                HashTst1 <- sapply(QCs$Path[wRdt[wY]], HashFun) 
                HashTst2 <- sapply(QCs2$Path[match(tst[wY], tst2)], HashFun) 
                wD <- which(HashTst1 == HashTst2)
                if (length(wD)) {
                  cat("Deleting redundant files\n")
                  for (w in wD) {
                    cat(paste0("   ", QCs$Path[wRdt[wY]][w], "\n"))
                    unlink(QCs$Path[wRdt[wY]][w])
                  }
                }
                QCs$Copy[wRdt] <- !(tst %in% tst2)
              }
              if (length(wN)) {
                shell(paste0("open \"", LocQCDir, "\""))
                shell(paste0("open \"", TargQCDir, "\""))
                print(QCs$`New name`[wRdt[wN]])
                stop("Unexpected conflict, investigate!")
              }
            }
          }
          whC <- which(QCs$Copy) # Update whC
          if (length(whC)) {
            # We will allow for up to one file to be copied at a time
            Witch <- whC[1]
            if (file.exists(LogFl)) {
              LastTime <- readLines(LogFl)
              LastTime <- grep("^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}\\:[0-9]{2}\\:[0-9]{2} - ", LastTime, value = TRUE)
              if (length(LastTime)) {
                LastTime <- unlist(strsplit(rev(LastTime)[1], " - "))
                LastTime <- LastTime[[1]]
                LastTime <- as.numeric(as.POSIXct(LastTime))
              } else { LastTime <- NA }
            } else {
              LastTime <- Sys.time()
              Log <- paste0(LastTime, " - Created log file")
              write(Log, LogFl)
              LastTime <- as.numeric(LastTime)
            }
            if (is.na(LastTime)) { LastTime <- as.numeric(Sys.time()) }
            #LastTime <- LastTime - Every*3600
            CurrTime <- as.numeric(Sys.time())
            if ((CurrTime - LastTime)/3600 >= Every) {
              # Rename file to apply time stamp
              tstRnm <- FALSE
              if (QCs$"New name"[Witch] != QCs$"Raw file"[Witch]) {
                if (file.exists(QCs$"New path"[Witch])) { stop() } else {
                  tstRnm <- file.rename(QCs$Path[Witch], QCs$"New path"[Witch])
                }
              }
              if (tstRnm) {
                tst <- try(fs::file_copy(QCs$"New path"[Witch],
                                         TargQCDir),
                           silent = TRUE)
              } else {
                # If we cannot rename the local file (e.g. locked for viewing), we copy+rename in one go
                tst <- try(fs::file_copy(QCs$Path[Witch],
                                         paste0(TargQCDir, "/", QCs$"New name"[Witch])),
                           silent = TRUE)
              }
              if (
                (!"try-error" %in% class(tst))
                &&
                (tst == gsub(topattern(LocQCDir), TargQCDir, QCs$"New path"[Witch]))
                &&
                (QCs$"New name"[Witch] %in% gsub(".*/", "", list.files(TargQCDir, recursive = TRUE)))
              ) {
                if (tstRnm) { tst2 <- try(unlink(QCs$"New path"[Witch]), silent = TRUE) } else {
                  tst2 <- try(unlink(QCs$"Raw file"[Witch]), silent = TRUE)
                }
                CurrTime0 <- Sys.time()
                if (file.exists(LogFl)) { Log <- readLines(LogFl) } else { Log <- c() }
                Log <- c(Log, paste0(as.character(CurrTime0), " - copied file \"", QCs$"New name"[Witch], "\""))
                write(Log, LogFl)
              } else {
                warning("Unexpected results, investigate!")
                write("Unexpected results, investigate!", paste0("error_log", gsub(":", "-", Sys.time()), ".txt"))
              }
            }
          }
        }
      }
    }
    # Also move QC03 to archive
    if (length(QC3)) {
      QC3s <- data.frame(Path = QC3, check.names = FALSE)
      # Deal with time-stamps
      # - 1) Files with a time stamp
      g1 <- grep("_[0-9]{14}\\.raw$", QC3s$`Raw file`, ignore.case = TRUE)
      #View(Raw[g1,])
      if (length(g1)) {
        QC3s$`New name`[g1] <- paste0(gsub("_[0-9]+\\.raw$", "", QC3s$`Raw file`[g1], ignore.case = TRUE),
                                      "_", gsub("[^0-9]", "", QC3s$Date[g1]), ".", QC3s$Extension[g1])
        w <- which(QC3s$`New name`[g1] != QC3s$`Raw file`[g1])
        if (length(w)) {
          #View(QC3s[g1[w], c("Raw file", "New name")])
          warning("Some predicted time stamps do not match the ones created by XCalibur, our assumptions as to how it created them may not hold! Investigate!")
          # For 1 project, there was a small difference of 1s, it is ok to stick to the old file name then
          QC3s$`New name`[g1][w] <- QC3s$`Raw file`[g1][w]
        }
      }
      # - 2) Files without a time stamp
      g2 <- c(1:nrow(QC3s))[which(!c(1:nrow(QC3s)) %in% g1)]
      if (length(g2)) {
        QC3s$`New name`[g2] <- paste0(gsub("\\.raw$", "", QC3s$`Raw file`[g2], ignore.case = TRUE),
                                      "_", gsub("[^0-9]", "", QC3s$Date[g2]), ".", QC3s$Extension[g2])
      }
      # Copy without further ado
    }
    # Delete blanks older than 2 weeks
    Blnks <- grep("blank|zig|zag", raws, value = TRUE, ignore.case = TRUE)
    if (length(Blnks)) {
      tst <- Sys.time() - file.info(Blnks)$mtime > 14*24*3600
      Blnks <- Blnks[which(tst)]
      for (blnk in Blnks) { unlink(Blnks) }
    }
  }
}
