require(rstudioapi)
require(openxlsx2)
require(parallel)
require(proteoCraft)
require(svDialogs)
require(data.table)

# Create parallel processing cluster
N.clust <- detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(stopCluster(parClust), silent = TRUE)
  parClust %<o% makeCluster(N.clust, type = "SOCK")
}
#stopCluster(parClust)

homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
fl <- paste0(homePath, "/Default_locations.xlsx")
flTst <- file.exists(fl)
if (!flTst) { file.copy(paste0(RPath, "/proteoCraft/extdata/Default_locations.xlsx"), fl) }
inRoot <- read_xlsx(fl)

# Target directory
targDir %<o% inRoot$Path[match("Search folder", inRoot$Folder)]
while(!dir.exists(targDir)) {
  targDir <- selectDirectory("Select groups_temp archive folder", path = "C:/")
}

# Where we will archive it
archDir %<o% inRoot$Path[match("Archive folder", inRoot$Folder)]
archDir <- selectDirectory("Select groups_temp archive folder", path = archDir)

# MS raw file archive(s)
MS_archDirs %<o% c(inRoot$Path[match("Archive folder", inRoot$Folder)],
                   gsub("/[^/]+$", "/Acquired data_v2", archDir))
MS_archDirs <- MS_archDirs[which(dir.exists(MS_archDirs))]

grpDirs <- list.dirs(targDir, recursive = FALSE, full.names = TRUE)
dirNms <- gsub(".*/", "", grpDirs)
names(grpDirs) <- dirNms
prjcts <- setNames(lapply(grpDirs, list.dirs, recursive = FALSE, full.names = TRUE),
                   dirNms)
prjcts <- prjcts[which(sapply(prjcts, length) > 0)]
prjcts <- listMelt(prjcts, names(prjcts), c("Project", "Group"))
prjcts$Files <- parLapply(parClust, prjcts$Project, list.files, full.names = TRUE, recursive = TRUE)
fls <- listMelt(prjcts$Files, prjcts$Project, c("File", "Project"))
tmp <- parSapply(parClust, fls$File, file.info)
tmp <- as.data.frame(t(tmp))
fls[, colnames(tmp)] <- tmp
#
# Now a first step would be to check .raw files and .d folders: are they in the archive?
rwFls <- fls[grep("\\.raw$", fls$File, ignore.case = TRUE),]
tmp <- gsub(topattern(paste0(targDir, "/")), "", rwFls$Project)
tmp <- strsplit(tmp, "/")
tmp <- Isapply(tmp, function(x) { unlist(x)[1:2] })
rwFls[, c("ResGrp", "Project_name")] <- tmp
rwFls$Name <- basename(rwFls$File)
rwFls$Archived <- apply(rwFls[, c("Name", "ResGrp", "Project_name", "size", "mtime")], 1, function(x) {
  #x <- rwFls[1, c("Name", "ResGrp", "Project_name", "size", "mtime")]
  #
  x <- unlist(x)
  # Try to match directories
  dirs <- paste0(MS_archDirs, "/", x[[2]], "/", x[[3]])
  dirs <- dirs[which(dir.exists(dirs))]
  if ((!length(dirs))&&(grepl("_[a-z]$", x[[3]]))) {
    dirs <- paste0(MS_archDirs, "/", x[[2]], "/", gsub("_[a-z]$", "", x[[3]]))
    dirs <- dirs[which(dir.exists(dirs))]
  }
  if ((!length(dirs))&&(!grepl("_[a-z]$", x[[3]]))) {
    w <- which(sapply(paste0(MS_archDirs, "/", x[[2]]), function(y) {
      #y <- paste0(MS_archDirs[1], "/", x[[2]])
      x[[3]] %in% gsub("_[a-z]$", "", list.dirs(y, recursive = FALSE, full.names = FALSE))
    }))
    if (length(w)) {
      dirs <- unlist(sapply(paste0(MS_archDirs[w], "/", x[[2]]), function(y) {
        #y <- paste0(MS_archDirs[1], "/", x[[2]])
        drs <- list.dirs(y, recursive = FALSE, full.names = FALSE)
        drs <- drs[which(gsub("_[a-z]$", "", drs) == x[[3]])]
        if (length(drs)) { drs <- paste0(y, "/", drs) }
        return(drs)
      }))
    }
  }
  # Search for the file
  if (length(dirs)) {
    fl <- paste0(dirs, "/", x[[1]])
    fl <- fl[which(file.exists(fl))]
    if (!length(fl)) {
      fl <- unlist(lapply(dirs, function(dr) { #dr <- dirs[1]
        list.files(dirs, "\\.raw$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
      }))
      fl <- fl[which(basename(fl) == x[[1]])]
    }
    if (length(fl)) {
      inf <- file.info(fl)
      res <- (inf$size == x[[4]])&(as.numeric(inf$mtime) == x[[5]])
    } else { res <- FALSE }
  } else { res <- FALSE }
  return(res)
})
#
tmp <- as.data.table(fls)
tmp <- tmp[, list(mtime = max(mtime)), by = list(Project = Project)]
tmp <- as.data.frame(tmp)
tmp$mtime <- as.POSIXct(tmp$mtime, origin = "1970-01-01 00:00.00 UTC")
wN <- which(!prjcts$Project %in% tmp$Project)
wY <- which(prjcts$Project %in% tmp$Project)
# Remove empty directories
empties <- prjcts$Project[wN]
if (length(empties)) {
  msg <- paste0("The following empty folders have been identified, which can we delete?")
  m <- max(c(250, nchar(empties)))
  empties2 <- sapply(empties, function(x) { paste(c(x, rep(" ", m-nchar(x))), collapse = "") })
  m <- match(dlg_list(empties2, empties2, TRUE, msg)$res, empties2)
  if (length(empties)) {
    empties <- empties[m]; rm(empties2)
    try(unlink(empties, recursive = TRUE), silent = TRUE)
  }
}
# Archive old directories
prjcts$mtime <- tmp$mtime[match(prjcts$Project, tmp$Project)]
prjcts$mtime <- as.POSIXct(prjcts$mtime, origin = "1970-01-01 00:00.00 UTC")
#View(prjcts[, c("Project", "mtime")])
now <- Sys.time()
ThreeMonthsAgo <- as.POSIXct(Sys.time() - 3*30*24*60*60, origin = "1970-01-01 00:00.00 UTC")
wO <- which(prjcts$mtime < ThreeMonthsAgo)
if (length(wO)) {
  oldies <- prjcts$Project[wO]
  msg <- paste0("The following folders haven't been modified since 3 months, can we archive them?")
  m <- max(c(250, nchar(oldies)))
  oldies2 <- sapply(oldies, function(x) { paste(c(x, rep(" ", m-nchar(x))), collapse = "") })
  m <- match(dlg_list(oldies2, oldies2, TRUE, msg)$res, oldies2)
  oldies <- oldies[m]; rm(oldies2)
  if (length(oldies)) {
    oldies <- data.frame(Project = oldies,
                         Outcome = NA)
    oldies$DestDir <- paste0(archDir, "/", gsub(".*/", "", dirname(oldies$Project)))
    for (i in 1:nrow(oldies)) { #i <- 1
      dir <- oldies$DestDir[i]
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      cmd <- paste0("mv \"", oldies$Project[i], "\" \"", dir, "\"")
      tst <- try(system(cmd), silent = TRUE)
      oldies$Outcome[i] <- (!"try-error" %in% class(tst))
      if (!oldies$Outcome[i]) {
        warning(paste0("Automated move to archive failed for project \"", oldies$Project[i], "\""))
      } else {
        cat(paste0("Project \"", oldies$Project[i], "\" successfully moved to archive\n"))
      }
    }
  }
}
