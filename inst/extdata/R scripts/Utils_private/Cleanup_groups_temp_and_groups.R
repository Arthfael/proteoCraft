###########
# Cleanup #
###########
#
# Removes all Thermo Raw and Bruker D folders older than 30 days as long as they are properly backed up on the Archive
#
require(svDialogs)
wd <- "...Search_Folder"
arch <- "...MS_File_Archive"
now <- Sys.time()

raws <- list.files(wd, "\\.raw$", recursive = TRUE, full.names = TRUE)
dirs <- list.dirs(wd, recursive = TRUE, full.names = TRUE)
dDirs <- grep("\\.d$", dirs, value = TRUE)

# Raw files
raws <- data.frame(Local = raws)
raws$Archive <- gsub(proteoCraft::topattern(wd), arch, raws$Local)
raws$BackedUp <- file.exists(raws$Archive)
raws$Created <- file.info(raws$Local)$ctime
raws$"Size (local)" <- file.size(raws$Local)
raws$"Size (archive)" <- NA
raws$"Size (archive)"[which(raws$BackedUp)] <- file.size(raws$Archive[raws$BackedUp])
sum(raws$BackedUp)
rawsRmv <- raws[which((raws$BackedUp) # Are they backed up ?
                      &(!is.na(raws$Created)) # Do they have a valid created date...
                      &(now - raws$Created > 30) # ... and is it at least 30 days old?
                      &(!is.na(raws$"Size (archive)")) # Can we measure they size in the archive...
                      &(raws$"Size (local)" == raws$"Size (archive)")) # ... and is it the same ?
                ,]
nr <- nrow(rawsRmv)
if (nr) {
  msg <- paste0("Found ", nr, " archived Thermo .raw files. Erase them?")
  ok2del <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  if (ok2del) {
    tst <- try(unlink(rawsRmv$Local), silent = TRUE)
    if ("try-error" %in% class(tst)) { warning("This didn't work, investigate!") } else { cat("Success!\n") }
  }
}
#unlink(rawsRmv$Local[w])

# Bruker D folders
dDirs <- data.frame(Local = dDirs)
dDirs$Archive <- gsub(proteoCraft::topattern(wd), arch, dDirs$Local)
dDirs$BackedUp <- dir.exists(dDirs$Archive)
dDirs$Created <- file.info(dDirs$Local)$ctime
dDirs$Files <- lapply(dDirs$Local, function(dr) { list.files(dr, recursive = TRUE, full.names = TRUE) })
dDirs$"Size (local)" <- sapply(dDirs$Files, function(fls) {
  sum(file.size(fls), na.rm = TRUE)
})
dDirs$nFiles <- sapply(dDirs$Files, length)
wBckp <- which(dDirs$BackedUp)
dDirs$"Files (archive)" <- dDirs$"Size (archive)" <- NA
dDirs$"Files (archive)"[wBckp] <- lapply(dDirs$Archive[wBckp], function(dr) {
  list.files(dr, recursive = TRUE, full.names = TRUE)
})
dDirs$"nFiles (archive)" <- sapply(dDirs$"Files (archive)", length)
dDirs$"Size (archive)"[wBckp] <- sapply(dDirs$"Files (archive)"[wBckp], function(fls) {
  sum(file.size(fls), na.rm = TRUE)
})
sum(dDirs$BackedUp)
dDirsRmv <- dDirs[which((dDirs$BackedUp) # Are they backed up ?
                        &(!is.na(dDirs$Created)) # Do they have a valid created date...
                        &(now - dDirs$Created > 30) # ... and is it at least 30 days old?
                        &(dDirs$nFiles == dDirs$"nFiles (archive)") # Is there the same number of files in both?
                        &(!is.na(dDirs$"Size (archive)")) # Can we measure they size in the archive...
                        &(dDirs$"Size (local)" == dDirs$"Size (archive)")) # ... and is it the same ?
                  ,]
nr <- nrow(dDirsRmv)
if (nr) {
  msg <- paste0("Found ", nr, " archived Bruker .d files. Erase them?")
  ok2del <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  if (ok2del) {
    for (x in 1:nr) { #x <- 1
      Ddr <- dDirsRmv$Local[[x]]
      print(Ddr)
      fls <- dDirsRmv$Files[[x]]
      tst <- try(unlink(fls, recursive = TRUE), silent = TRUE)
      drs <- c(Ddr, list.dirs(Ddr, recursive = TRUE))
      drs <- drs[order(nchar(drs), decreasing = TRUE)]
      tst <- list()
      for (dr in drs) {
        tst[[dr]] <- try(fs::dir_delete(dr), silent = TRUE)
        if ("try-error" %in% class(tst[[dr]])) {
          cmd <- paste0("rmdir \"", dr, "\"")
          system(cmd)
        }
      }
    }
  }
}
#proteoCraft::openwd(dDirs$Local[1])
