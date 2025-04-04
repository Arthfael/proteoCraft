# Archive acquired files
options(stringsAsFactors = FALSE) # Needed because we will stick to an older version

if (!require(svDialogs)) { install.packages("svDialogs") }
require(svDialogs)
if (!require(tools)) { install.packages("tools") }
require(tools)
if (!require(Rcpp)) { install.packages("fs") }
require(fs)
if (!require(fs)) { install.packages("fs") }
require(fs)
#if (!require(openssl)) { install.packages("openssl") }
#require(openssl)
if (!require(digest)) { install.packages("digest") }
require(digest)

# Start with clean slate
rm(list = ls())

if (!require(proteoCraft)) {
  openwd <- function(Dir) {
    cmd <- paste0("explorer \"", normalizePath(Dir), "\"")
    #cat(cmd)
    suppressMessages(suppressWarnings(shell(cmd, intern = TRUE)))
  }
} else { openwd <- proteoCraft::openwd() }

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

# Select target directory
# The script will identify all raw files in the directory and transfer/back them up to mirrored folders on msmonster01 and the LSF archive
TargDir <- choose.dir("D:/Data/", "Select a completed project folder")
if (!is.na(TargDir)) {
  TargDir <- normalizePath(choose.dir("D:/Data/", "Select a completed project folder"), winslash = "/")
  if (dir.exists(TargDir)) {
    if (grepl("D:/Data/", TargDir)) {
      msg <- "Please confirm: has data acquisition in this folder completed?
(This script will fail if raw files are stil being recorded!)"
      OngoingAcq <- !c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
      if (!OngoingAcq) {
        #openwd(TargDir)
        # 
        msg <- "Should we copy files to the analysis PC = msmonster01? (y/n)"
        dflt <- c("y", "n")[grepl("Solgate", TargDir) + 1]
        SkipMSMonster <- tolower(dlg_input(msg, dflt)$res)
        while (!SkipMSMonster %in% c("y", "n")) { SkipMSMonster <- tolower(dlg_input(msg, dflt)$res) }
        SkipMSMonster <- !c(TRUE, FALSE)[match(SkipMSMonster, c("y", "n"))]
        #
        if (!SkipMSMonster) {
          msg <- "Do you want to transfer blanks to msmonster01 too?"
          BlanksToo <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
        } else { BlanksToo <- FALSE }
        #
        TargSubDir <- gsub("^D:/Data/", "", TargDir)
        DestDir1 <- "T:/" # Top level destination directory on msmonster01
        DestDir2 <- "//archive3.ist.local/archive-lsfgrp/MS/Acquired data/" # Top level destination directory on LSF archive
        Fls <- list.files(TargDir, recursive = TRUE)
        Raw <- grep("\\.raw$", Fls, ignore.case = TRUE, value = TRUE)
        Seq <- grep("\\.sld", Fls, ignore.case = TRUE, value = TRUE)
        #
        if (!length(Raw)) {
          warning("No raw files found in the folder!")
        } else {
          # Raw files
          Raw <- data.frame("Raw file" = Raw,
                            "Path" = paste0(TargDir, "/", Raw),
                            "Copy" = TRUE,
                            "In there" = FALSE,
                            check.names = FALSE)
          # Also sequence files
          Seq <- grep(".sld", list.files(TargDir), value = TRUE)
          # Deal with files in sub-directory
          Raw$Dir <- gsub("/[^/]+$", "", Raw$Path)
          Raw$`Raw file` <- gsub(".+/", "", Raw$`Raw file`)
          Raw$"Creation time (POSIXct)" <- as.POSIXct(sapply(Raw$Path, function(x) { file.info(x)$ctime }), origin = "1970-01-01 00:00.00 UTC")
          Raw$"Modified time (POSIXct)" <- as.POSIXct(sapply(Raw$Path, function(x) { file.info(x)$mtime }), origin = "1970-01-01 00:00.00 UTC")
          Raw$"Creation time" <- sapply(Raw$Path, function(x) { as.character(file.info(x)$ctime) })
          Raw$"Modified time" <- sapply(Raw$Path, function(x) { as.character(file.info(x)$mtime) })
          Raw$Year <- as.integer(format(Raw$`Creation time (POSIXct)`, "%Y"))
          Raw$Month <- as.numeric(format(Raw$`Creation time (POSIXct)`, "%m"))
          Raw$Day <- as.integer(format(Raw$`Creation time (POSIXct)`, "%d"))
          Raw$`New name` <- Raw$`Raw file`
          Raw$Extension <- sapply(strsplit(Raw$`Raw file`, "\\."), function(x) { rev(unlist(x))[1] }) #To accept all combinations of upper/lower case
          Raw$Size <- file.size(Raw$Path)
          Raw <- Raw[order(Raw$`Creation time (POSIXct)`, decreasing = FALSE),]
          #
          # Deal with time-stamps
          # - 1) Files with a time stamp
          g1 <- grep("_[0-9]{14}\\.raw$", Raw$`Raw file`, ignore.case = TRUE)
          #View(Raw[g1,])
          if (length(g1)) {
            Raw$`New name`[g1] <- paste0(gsub("_[0-9]+\\.raw$", "", Raw$`Raw file`[g1], ignore.case = TRUE),
                                         "_", gsub("[^0-9]", "", Raw$`Creation time (POSIXct)`[g1]), ".", Raw$Extension[g1])
            w <- which(Raw$`New name`[g1] != Raw$`Raw file`[g1])
            if (length(w)) {
              # TO DO: edit this to accept the occasional 1s difference - this type of error is acceptable
              View(Raw[g1[w], c("Raw file", "New name")])
              warning("Some predicted time stamps do not match the ones created by XCalibur, our assumptions as to how it created them may not hold! Investigate!")
              # For 1 project, there was a small difference of 1s, it is ok to stick to the old file name then.
              # Sometimes, the issue is also that in the past I have sometimes added time stamps manually but did not have seconds-level information.
              Raw$`New name`[g1][w] <- Raw$`Raw file`[g1][w]
            }
          }
          # - 2) Files without a time stamp
          g2 <- c(1:nrow(Raw))[which(!c(1:nrow(Raw)) %in% g1)]
          if (length(g2)) {
            Raw$`New name`[g2] <- paste0(gsub("\\.raw$", "", Raw$`Raw file`[g2], ignore.case = TRUE),
                                         "_", gsub("[^0-9]", "", Raw$`Creation time (POSIXct)`[g2]), ".", Raw$Extension[g2])
          }
          # Check for name uniqueness
          if (length(unique(Raw$`New name`)) < nrow(Raw)) { stop("Check the files renaming code!!!") }
          Raw$"New (no ext.)" <- gsub("\\.raw$", "", Raw$`New name`, ignore.case = TRUE)
          # Hash
          HashFun <- function(x) { as.character(tools::md5sum(x)) } #tools::md5sum takes a character path, not a connection
          #HashFun <- function(x) { as.character(openssl::md5(file(x))) }
          #HashFun <- function(x) { digest(file = file(x), algo = "md5") }
          # md5 seems faster than sha1; both are "usually sufficient for collision-resistant identifiers" but "no longer considered secure for cryptographic purposes."
          # The digest implementation of md5 seems ~20% faster than the openssl one
          # but digest syntax is trickier, until I get it right I will stick to openssl
          # I ended up switching to tools::md5sum, which is about as fast as digest but easier to write.
          cat("Calculating raw file hashes...\n")
          Raw$Hash <- sapply(Raw$Path, HashFun) 
          # NB: It is intentional that I check/modify names sequentially, not independently per directory.
          # (but do check that this works!)
          Raw2 <- list()
          Dirs <- setNames(c(DestDir1, DestDir2), c("msmonster01", "Archive"))
          Copied <- Skip <- setNames(rep(FALSE, length(Dirs)), names(Dirs))
          if (sum(Raw$Copy)) {
            for (i in 1:2) { #i <- 1
              Dir <- Dirs[i]
              DirNm <- names(Dirs)[i]
              Raw2[[DirNm]] <- Raw
              #Dir <- DestDir1
              #Dir <- DestDir2
              if (Dir == DestDir1) { Skip[DirNm] <- SkipMSMonster }
              if (!Skip[DirNm]) {
                if (!dir.exists(Dir)) {
                  if (Dir == DestDir1) {
                    warning("\"groups_temp\" (aka \"msmonster01\" or \"T:\") is unavailable, please run the connection script on the desktop of the acquisition PC!")
                  }
                  if (Dir == DestDir2) {
                    warning("\"//archive3.ist.local/archive-lsfgrp/MS/Acquired data\" is unavailable, check that the external drive is correctly mapped!")
                  }
                } else {
                  if (Dir == DestDir1) {
                    # Check available space on msmonster01
                    # From https://stackoverflow.com/questions/32200879/how-to-get-disk-space-of-windows-machine-with-r
                    disks <- system("wmic logicaldisk get size,freespace,caption", inter=TRUE)
                    disks <- read.fwf(textConnection(disks[1:(length(disks)-1)]), 
                                      widths=c(9, 13, 13), strip.white=TRUE, stringsAsFactors=FALSE)
                    colnames(disks) <- disks[1,]
                    disks <- disks[-1,]
                    rownames(disks) <- NULL
                    FrSp <- as.numeric(disks$FreeSpace[match("T:", disks$Caption)])
                    # Identify blank samples - ask if they should be removed
                    Blks <- grep("(([Bb]lank)|([Zz]ig-?[Zz]ag))", Raw2[[DirNm]]$"Raw file")
                    if (length(Blks)) { Raw2[[DirNm]]$Copy[Blks] <- BlanksToo }
                  }
                  #
                  dir <- paste0(Dir, TargSubDir) # Mirror of our target directory (may have sub-directories)
                  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
                  if (!"DestDir" %in% colnames(Raw2[[DirNm]])) { Raw2[[DirNm]]$DestDir <- NA }
                  Raw2[[DirNm]]$DestDir <- paste0(dir, substr(Raw2[[DirNm]]$Dir, nchar(TargDir)+1,
                                                              nchar(Raw2[[DirNm]]$Dir))) # Proper destination directory, may include sub-directories
                  for (dr in unique(Raw2[[DirNm]]$DestDir)) {
                    if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
                  }
                  NewFls <- paste0(Raw2[[DirNm]]$DestDir, "/", Raw2[[DirNm]]$`New name`)
                  DestFls <- list.files(dir, recursive = TRUE) # Files in destination directory
                  OLRaw <- grep("\\.raw$", DestFls, ignore.case = TRUE, value = TRUE)
                  SzTst <- sapply(paste0(dir, "/", OLRaw), file.size)
                  OLRaw <- OLRaw[which(SzTst > 1024)] # Remove dummies, abortive files, etc...
                  OLRaw2 <- gsub("\\.raw$", "", OLRaw, ignore.case = TRUE)
                  w <- which(Raw2[[DirNm]]$"New (no ext.)" %in% OLRaw2)
                  if (length(w)) {
                    Raw2[[DirNm]]$Copy[w] <- FALSE
                    m <- match(Raw2[[DirNm]]$"New (no ext.)"[w], OLRaw2)
                    # Calculate file checksums - this is slow but very safe
                    w2 <- which(is.na(Raw2[[DirNm]]$Hash[w]))
                    if (length(w2)) {
                      Raw2[[DirNm]]$Hash[w[w2]] <- sapply(Raw2[[DirNm]]$Path[w[w2]], HashFun)
                      Raw$Hash[w[w2]] <- Raw2[[DirNm]]$Hash[w[w2]] # So we don' do it twice for the same files!
                    }
                    tst <- sapply(paste0(dir, "/", OLRaw[m]), HashFun)
                    w <- w[which(Raw2[[DirNm]]$Hash[w] != tst)] # Filter by Hash: we rename if Hash is different, otherwise we skip
                    if (length(w)) {
                      Raw2[[DirNm]]$Copy[w] <- TRUE
                      kount <- rep(1, length(w))
                      g <- grep("_rerun", Raw2[[DirNm]]$"New (no ext.)"[w])
                      if (length(g)) {
                        kount[g] <- as.numeric(gsub("^.+_rerun", "", Raw2[[DirNm]]$"New (no ext.)"[w][g]))+1
                      }
                      Raw2[[DirNm]]$"New (no ext.)"[w] <- paste0(Raw2[[DirNm]]$"New (no ext.)"[w], "_rerun", kount)
                      Raw2[[DirNm]]$"New name"[w] <- apply(Raw2[[DirNm]][w, c("New (no ext.)", "Extension")], 1, paste, collapse = ".")
                    }
                  }
                  CopyFls <- FALSE
                  w <- which(Raw2[[DirNm]]$Copy)
                  #View(Raw2[[DirNm]])
                  if (length(w)) {
                    CopyFls <- TRUE
                    # If on the LSF archive, space is not an issue.
                    # On msmonster01, however, it is necessary to check
                    if (Dir == DestDir1) {
                      sztst <- sum(Raw2[[DirNm]]$Size[w])
                      if (sztst > FrSp) { 
                        warning("No files will be copied as there is insufficient free space on \"msmonster01\"!")
                        CopyFls <- FALSE
                      } else {
                        if (sztst >= FrSp*0.8) { 
                          warning("Free space on \"msmonster01\" is running out, time for a clean-up!")
                        }
                      }
                    }
                    if (CopyFls) {
                      tst <- substr(Raw2[[DirNm]]$Dir[w], 1, nchar(TargDir))
                      w2 <- which(tst != TargDir)
                      if (length(w2)) { stop("Something went wrong!") }
                      for (diri in Raw2[[DirNm]]$DestDir[w]) { if (!dir.exists(diri)) { dir.create(diri, recursive = TRUE) } }
                      w3 <- which(!apply(Raw2[[DirNm]][w, c("DestDir", "New name")], 1, function(x) {
                        x[[2]] %in% list.files(x[[1]])
                      }))
                      if (length(w3)) {
                        cat(paste0("Copying these files to \"", dir, "\":\n",
                                   paste0(" - ", Raw2[[DirNm]]$"New name"[w[w3]], collapse = "\n"), "\n\n"))
                        tst <- try(fs::file_copy(Raw2[[DirNm]]$Path[w[w3]], NewFls[w[w3]]), silent = TRUE)
                        #Sys.setFileTime(NewFls, Raw2[[DirNm]]$`Modified time`[w[w3]]) # We actually want to set both created and modified time
                        if (!"try-error" %in% class(tst)) {
                          # Adjust created and modified time
                          for (w4 in w3) {
                            fl <- NewFls[w[w4]]
                            flnm4pwrshll <- nrmPath4PS(nufl)
                            # Creation time
                            tm <- Raw2[[DirNm]]$`Creation time (POSIXct)`[w[w4]]
                            tmY <- as.integer(format(tm, "%Y"))
                            tmM <- as.integer(format(tm, "%m"))
                            tmD <- as.integer(format(tm, "%d"))
                            tmHr <- as.integer(format(tm, "%H"))
                            tmMn <- as.integer(format(tm, "%M"))
                            tmSc <- as.integer(format(tm, "%S"))
                            cmd <- paste0('powershell -command "$(Get-Item \"', flnm4pwrshll, '\").CreationTime=$(Get-Date -Year ', tmY,
                                          ' -Month ', tmM, ' -Day ', tmD, ' -Hour ', tmHr, ' -Minute ', tmMn, ' -Second ', tmSc, ')"')
                            #cat(cmd)
                            res <- system(cmd, intern=TRUE)
                            # Modified time (just a precaution: should not be different)
                            tm <- Raw2[[DirNm]]$`Modified time (POSIXct)`[w[w4]]
                            tmY <- as.integer(format(tm, "%Y"))
                            tmM <- as.integer(format(tm, "%m"))
                            tmD <- as.integer(format(tm, "%d"))
                            tmHr <- as.integer(format(tm, "%H"))
                            tmMn <- as.integer(format(tm, "%M"))
                            tmSc <- as.integer(format(tm, "%S"))
                            cmd <- paste0('powershell -command "$(Get-Item \"', flnm4pwrshll, '\").LastWriteTime=$(Get-Date -Year ', tmY,
                                          ' -Month ', tmM, ' -Day ', tmD, ' -Hour ', tmHr, ' -Minute ', tmMn, ' -Second ', tmSc, ')"')
                            #cat(cmd)
                            res <- system(cmd, intern=TRUE)
                          }
                        }
                      }
                      # Also copy sequences
                      if (length(Seq)) {
                        SeqFull <- paste0(TargDir, "/", Seq)
                        SeqDir <- dirname(SeqFull)
                        SeqDir2 <- paste0(dir, substr(SeqDir, nchar(TargDir)+1, nchar(SeqDir)))
                        file.copy(SeqFull, paste0(SeqDir2, "/", Seq))
                      }
                    }
                  }
                  #Copied[DirNm] <- (length(tst) > 1)||(class(tst) != "try-error")
                  # Re-check copied material for inconsistencies:
                  #if (Copied[DirNm]) { ## Check that all copied files are present
                  wh <- which(!apply(Raw2[[DirNm]][, c("New name", "DestDir")], 1, function(x) {
                    x[[1]] %in% list.files(x[[2]])
                  }))
                  Copied[DirNm] <- length(wh) == 0
                  if (Copied[DirNm]) { ## Check that copied files have the right size
                    wh <- which(Raw2[[DirNm]]$Size - file.size(NewFls) != 0)
                    Copied[DirNm] <- length(wh) == 0
                    if (Copied[DirNm]) { # Finally, calculate and check hashes
                      cat("Checking copied raw file hashes...\n")
                      FinalHashes <- sapply(NewFls, HashFun)
                      wh <- which(Raw2[[DirNm]]$Hash != FinalHashes)
                      Copied[DirNm] <- length(wh) == 0
                      if ((Copied[DirNm])&&(CopyFls)) {
                        cat(paste0("Files successfully copied to \"", dir, "\"!"))
                        Raw$"In there" <- TRUE
                        message("\n")
                      } else {
                        WarnMsg <- paste0("The md5 hash of some original raw files and their copy in\" ", dir,
                                          "\" differ, investigate!")
                      }
                    } else {
                      WarnMsg <- paste0("The size of some of the raw files copied to \"", dir,
                                        "\" differs from that of the originals, investigate!")
                    }
                  } else {
                    WarnMsg <- paste0("Some raw files missing from \"", dir, "\", investigate!")
                  }
                  #} else {
                  #  wh <- 1:nrow(Raw2[[DirNm]])
                  #  WarnMsg <- paste0("Error during copying process to ", Dir, ", investigate!")
                  #}
                  if (Copied[DirNm]) {
                    write.csv(Raw2[[DirNm]], paste0(dir, "/File names map.csv"), row.names = FALSE)
                  } else {
                    FaultyFls <- paste0(Raw2[[DirNm]]$`New name`[wh])
                    warning(paste0(WarnMsg, "\nFiles affected:\n", paste0(" - ", FaultyFls, "\n", collapse = "")))
                    # Note: the warning about missing raw files on msmonster01 can be ignored if affecting only blanks and the user has chosen not to transfer blanks!
                  }
                }
              }
            }
            # Write copying log and map from old names to new ones
            if (sum(Copied)) {
              Tm <- Sys.time()
              temp <- c("ARCHIVING SCRIPT LOG",
                        "--------------------",
                        "",
                        paste0("Archiving script run on ", Tm),
                        "")
              for (i in which(Copied)) { #i <- 1
                Dir <- paste0(Dirs[i], TargSubDir)
                DirNm <- names(Dirs)[i]
                raw2 <- Raw2[[DirNm]]
                wY <- which(raw2$Copy)
                wN <- which(!raw2$Copy)
                a <- paste0("Files successfully copied to ", DirNm, ":")
                temp <- c(temp,
                          "",
                          "##########",
                          "",
                          a,
                          paste(rep("-", nchar(a)), collapse = ""),
                          "   Old name\t\t\t\t\tNew name\t\t\t\t\tDirectory",
                          paste0(" - ", apply(raw2[wY, c("Raw file", "New name", "DestDir")], 1, paste0, collapse = "\t\t\t\t\t")),
                          "")
                if (length(wN)) {
                  temp <- c(temp,
                            "Files not copied:",
                            paste0(" - ", raw2$"Raw file"[wN]),
                            "")
                }
              }
              write(temp, paste0(TargDir, "/Archiving script log_", gsub("[^0-9]", "", Tm), ".txt"))
            } else { cat("Nothing was copied!\\n") }
          }
          #
          # Manually inspect the results; if you are happy, delete (forever, they are too large) the local files!
          if (Copied["Archive"]) {
            msg <- "Do you want to delete local files now?"
            SearchNDestroy <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
            if (SearchNDestroy) {
              Whodunnit <- dlg_input("Enter your first and last names:\n(For now this script will operate under the assumption that users are honest.\nDo not make me change this policy!)")$res
              Whodunnit <- paste(sapply(strsplit(Whodunnit, " +"), function(x) {
                paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2, nchar(x))))
              }), collapse = " ")
              msg <- "I will now open both the original acquisition directory and its mirror in the Archive.\nYou will need to check the files manually in both, then you will be prompted again to authorize deletion of the local files.\nOk?"
              TwoMinutesToMidnight <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
              if (TwoMinutesToMidnight) {
                openwd(TargDir)
                #cat(normalizePath(TargDir))
                openwd(paste0(DestDir2, TargSubDir))
                #cat(normalizePath(paste0(DestDir2, TargSubDir)))
                msg <- "Have you checked the files?\nDue to their size, they will be erased without sending them to the recycle bin.\nCarry on?\nBY DOING THIS, YOU ARE ENGAGING YOUR RESPONSIBILITY IN CASE SOMETHING GOES WRONG!!!"
                NukeMFromOrbit <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
                if (NukeMFromOrbit) {
                  Tm2 <- Sys.time()
                  wh <- which(file.exists(paste0(DestDir2, TargSubDir, "/", Raw$`New name`)))
                  fls <- dlg_list(Raw$`Raw file`[wh], Raw$`Raw file`[wh], TRUE, "Select files to include in calculations of useful instrument time (to charge):")$res
                  m <- match(fls, Raw$`Raw file`[wh])
                  TimeDiff <- difftime(Raw$`Modified time (POSIXct)`[wh[m]], Raw$`Creation time (POSIXct)`[wh[m]], units = "hours")
                  temp2 <- c("RAW FILES DELETION LOG",
                             "-----------------------",
                             "",
                             paste0("Folder: ", normalizePath(TargDir)),
                             paste0(Tm2, ": ", Whodunnit, " authorized the deletion, following archiving, of the following local raw files:"),
                             "",
                             paste0(" - ", Raw$"Raw file"[wh]),
                             "",
                             paste0("Selected instrument run time to charge: ", round(sum(TimeDiff), 2), " h"),
                             "")
                  write(temp2, paste0(TargDir, "/Raw files deletion log_", gsub("[^0-9]", "", Tm2), ".txt"))
                  system(paste0("open \"", TargDir, "/Raw files deletion log_", gsub("[^0-9]", "", Tm2), ".txt\""))
                  unlink(Raw$Path[wh])
                }
              }
            }
          }
        }
      }
    } else { warning("Raw files subfolder must be in D:/Data/! Skipping...") }
  }  
}

