# Transfer data from old archives to new one

require(proteoCraft)
require(tools)
require(fs)
wd <- "C:/Users/anicolas/Documents"

# Better file.exists, dir.exists and file.size functions for dealing with the LONGPATH issue...
lp <- function(p) {
  p <- normalizePath(p, winslash = "\\", mustWork = FALSE)
  paste0("\\\\?\\", p)
}
safe_file_exists <- function(files) {
  res <- file.exists(files)
  #res <- as.logical(file.access(files) + 1) # Alternative
  w <- which(!res)
  if (length(w)) {
    # That's the GODDAM Windows long path issue!!!
    # Luckily, chatGPT has a great solution!
    p <- normalizePath(files[w], winslash = "\\")
    p2 <- paste0("\\\\?\\", p)
    res[w] <- file.exists(p2)
  }
  return(res)
}
safe_dir_exists <- function(dirs) { # Not tested!
  res <- dir.exists(dirs)
  w <- which(!res)
  if (length(w)) {
    # That's the GODDAM Windows long path issue!!!
    # Luckily, chatGPT has a great solution!
    res[w] <- dir.exists(lp(dirs[w]))
  }
  return(res)
}
safe_file_size <- function(files) {
  res <- file.size(files)
  w <- which(is.na(res))
  if (length(w)) {
    res[w] <- file.size(lp(files[w]))
  }
  return(res)
}
is_system_or_hidden <- function(file) {
  out <- suppressWarnings(system2("attrib", shQuote(file), stdout = TRUE, stderr = TRUE))
  if (!length(out)) { return(FALSE) }
  out <- gsub(proteoCraft::topattern(file, start = FALSE), "", gsub("\\\\", "/", out))
  return(grepl("[SH]", out))
}
safe_file_copy <- function(from, to, copy.date = TRUE, overwrite = TRUE) {
  ok <- file.copy(lp(from), lp(to), copy.date = copy.date, overwrite = overwrite)
  if (!ok) { ok <- is_system_or_hidden(from) }
  return(ok)
}
safe_dir_create <- function(paths) { # Recursive by default, which is good
  paths <- lp(paths)
  paths <- shQuote(paths)
  for (p in paths) {
    system2("cmd", c("/c", "mkdir", p))
  }
}

# Create hash function
HashFun <- function(x) { as.character(tools::md5sum(x)) } #tools::md5sum takes a character path, not a connection

RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
parSrc <- paste0(libPath, "/extdata/R scripts/Sources/make_check_Cluster.R")
source(parSrc)
clusterExport(parClust, "HashFun", envir = environment())

oldArch <- c("B:/archive/mspecgrp/MS/Acquired_data",
             "B:/group/lsfgrp/Mass_Spec/Acquired_data_v2")
newArch <- "B:/group/mspecgrp/_Archive/Raw"
stopifnot(sum(safe_dir_exists(oldArch)) == length(oldArch),
          safe_dir_exists(newArch),
          nchar(newArch) <= max(nchar(oldArch)))

# We will take into account the groups/project sub-structure.
# Also, because of the size of the data to transfer, we will take our sweet time and check at each step.

file_Tables <- list()
for (oldDir in oldArch) { #oldDir <- oldArch[1]
  # Step 1:
  # Identify subfolder structure
  grpDirs <- list.dirs(oldDir, recursive = FALSE, full.names = TRUE)
  projDirs <- setNames(lapply(grpDirs, function(grpDir) {
    list.dirs(grpDir, recursive = FALSE, full.names = TRUE)
  }), grpDirs)
  projDirs <- listMelt(projDirs, ColNames = c("Path", "Group"))[, c("Group", "Path")]
  projDirs$Group <- gsub(".*/", "", projDirs$Group)
  projDirs$Project <- gsub(".*/", "", projDirs$Path)
  projDirs$New_path <- paste0(newArch, "/", do.call(paste, c(projDirs[, c("Group", "Project")], sep = "/")))
  #
  file_Tables[[oldDir]] <- list()
  for (i in 1:nrow(projDirs)) { #i <- 1
    cat(paste0("Project ", projDirs$Path[i], "...\n"))
    fls <- list.files(projDirs$Path[i], recursive = TRUE)
    if (length(fls)) {
      flsTbl <- data.frame(File = fls,
                           Old_path = paste0(projDirs$Path[i], "/", fls),
                           New_path = paste0(projDirs$New_path[i], "/", fls))
      #flsTbl$Exists <- as.logical(file.access(flsTbl$Old_path) + 1)
      flsTbl$Exists <- safe_file_exists(flsTbl$Old_path)
      stopifnot(sum(!flsTbl$Exists) == 0)
      flsTbl$Old_dir <- dirname(flsTbl$Old_path)
      flsTbl$New_exists <- safe_file_exists(flsTbl$New_path)
      flsTbl$New_dir <- gsub("/[^/+]$", "", flsTbl$New_path)
      flsTbl$Copy_me <- !flsTbl$New_exists
      #sum(flsTbl$Copy_me)
      flsTbl$Old_size <- safe_file_size(flsTbl$Old_path)
      stopifnot(sum(is.na(flsTbl$Old_size)) == 0)
      flsTbl$New_size <- NA
      wNwE <- which(flsTbl$New_exists)
      if (length(wNwE)) {
        flsTbl$New_size[wNwE] <- safe_file_size(flsTbl$New_path[wNwE])
        flsTbl$Copy_me[wNwE] <- flsTbl$Old_size[wNwE] != flsTbl$New_size[wNwE]
        #sum(flsTbl$Copy_me)
      }
      wCopy <- which(flsTbl$Copy_me)
      #wCopy <- 1:nrow(flsTbl)
      if (length(wCopy)) {
        cat("   Copying...\n")
        #
        # First create new directories
        old_drs <- unique(gsub("/[^/]+$", "", flsTbl$Old_path[wCopy]))
        old_drs <- old_drs[order(nchar(old_drs), decreasing = TRUE)]
        drs <- substr(old_drs, nchar(projDirs$Path[i])+1, nchar(old_drs))
        new_drs <- paste0(projDirs$New_path[i], drs)
        w <- which(!safe_dir_exists(new_drs))
        if (length(w)) {
          safe_dir_create(new_drs[w])
        }
        w <- which(!safe_dir_exists(new_drs))
        stopifnot(length(w) == 0)
        #new_drs[w]
        #nchar(new_drs[w])
        #
        wCopy2 <- wCopy[order(flsTbl$Old_size[wCopy], decreasing = TRUE)]
        clusterExport(parClust, list("flsTbl", "wCopy2", "wd", "lp", "is_system_or_hidden", "safe_file_copy"), envir = environment())
        # tst <- try({
        #   Map(copyFun,
        #       file = flsTbl$Old_path[wCopy2],
        #       from = flsTbl$Old_dir[wCopy2],
        #       to = flsTbl$New_dir[wCopy2])
        # })
        tst <- parLapply(parClust, wCopy2, function(x) { #x <- wCopy2[1]
          currWD <- getwd()
          on.exit(setwd(currWD))
          setwd(flsTbl$Old_dir[x])
          # I am going for file.copy because it seems to have some advantages over fs::file_copy:
          # - Seems to fail less frequently
          # - Can in theory preserve file dates metadata... at least some
          safe_file_copy(flsTbl$Old_path[x], flsTbl$New_path[x], copy.date = TRUE, overwrite = TRUE)
        })
        flsTbl$New_exists[wCopy] <- safe_file_exists(flsTbl$New_path[wCopy])
        flsTbl$New_size[wCopy] <- safe_file_size(flsTbl$New_path[wCopy])
      }
      flsTbl$Success <- (flsTbl$New_exists)&(flsTbl$New_size == flsTbl$Old_size)
      w <- which(!flsTbl$Success)
      if (length(w)) {
        # We skip locked/hidden system files
        flsTbl$Success[w] <- vapply(flsTbl$Old_path[w], is_system_or_hidden, TRUE)
      }
      Success <- sum(!flsTbl$Success) == 0
      stopifnot(Success)
    } else {
      cat("   Empty project, skipping.../n")
    }
    cat("\n")
    file_Tables[[projDirs$Project[i]]] <- flsTbl
  }
  #
  # There could also be files at the top level of the project folder or group folders...
  projFls <- setNames(lapply(c(oldDir, grpDirs), function(dr) {
    list.files(dr, recursive = FALSE, full.names = TRUE, include.dirs = FALSE)
  }), c(oldDir, grpDirs))
  
  for (grpDir in grpDirs) { #grpDir <- grpDirs[1]
    
  }
  
  
  
}

tst <- lapply(oldArch, list.files, recursive = TRUE, full.names = TRUE)

