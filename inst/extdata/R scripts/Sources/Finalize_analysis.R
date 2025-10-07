# Write Materials and Method to text file
setwd(wd)
tmp <- paste0("MatMet <- ", unlist(MatMetCalls$Calls))
tmpSrc <- paste0(wd, "/tmp.R")
write(tmp, tmpSrc)
MatMetFl <- paste0(wd, "/Materials and methods_WIP.docx")
tst <- try({
  source(tmpSrc, local = FALSE)
  #rstudioapi::documentOpen(tmpSrc)
  MatMet %<o% MatMet
  print(MatMet, target = MatMetFl)
}, silent = TRUE)
if ("try-error" %in% class(tst)) {
  warning("Couldn't write materials and methods template, investigate...")
} else {
  system(paste0("open \"", MatMetFl, "\""))
  dlg_message("Check materials and methods text, make any necessary edits, then save and click ok", "ok")
}
unlink(tmpSrc)

# Write pdf report
#if (scrptType == "withReps") {
  tst <- try({
    tmp <- paste0("Report <- ", unlist(ReportCalls$Calls))
    tmpSrc <- paste0(wd, "/tmp.R")
    write(tmp, tmpSrc)
    source(tmpSrc, local = FALSE)
    Report %<o% Report
    # Write report to Word file
    print(Report, target = paste0(wd, "/Workflow control/Analysis report.docx"))
    #system(paste0("open \"", wd, "/Workflow control/Analysis report.docx"))
    unlink(tmpSrc)
  }, silent = TRUE)
  if ("try-error" %in% class(tst)) { warning("Couldn't write pdf report, investigate...") }
#}

# Save session info
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
if (scrptType == "withReps") { dirlist <- unique(c(dirlist, dir)) }
writeLines(capture.output(sessionInfo()), paste0(dir, "/sessionInfo.txt"))

# Copy script itself and guide to the results there:
#file.copy(system.file("extdata/R scripts", "Regulation analysis - master script.R", package = "proteoCraft"), wd, overwrite = TRUE)
if (dirname(ScriptPath) != wd) { file.copy(ScriptPath, wd, overwrite = TRUE) }
file.copy(system.file("extdata", paste0("Guide to the Results (", c("", "no ")[match(scrptType,
                                                                                     c("withReps", "noReps"))],
                                        "replicates).docx"), package = "proteoCraft"), wd, overwrite = TRUE)

## Copy final results to the delivery folder
RPath %<o% as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath %<o% paste0(RPath, "/proteoCraft")
homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
dirLocsFl <- paste0(homePath, "/Default_locations.xlsx")
dirLocs <- openxlsx2::read_xlsx(dirLocsFl)
if ((exists("outDir"))&&(length(outDir) == 1)&&(!is.na(outDir))&&(dir.exists(outDir))) {
  dflt <- outDir
} else {
  dflt <- dirLocs$Path[match("Results delivery folder", dirLocs$Folder)]
}
if (!dir.exists(dflt)) { dflt <- wd }
outDir %<o% rstudioapi::selectDirectory("Select data delivery folder",
                                        path = dflt)
ok2Deliver %<o% ((exists("outDir"))&&(length(outDir) == 1)&&(!is.na(outDir))&&(dir.exists(outDir)))
dataDeliveryOk %<o% FALSE
if (ok2Deliver) {
  Tsts <- list()
  # - 1/ MS files
  dir <- paste0(outDir, "/1_MS_files")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  if (writeRaws) {
    Tsts$"MS files" <- list()
    xst <- sapply(rawFiles, file.exists)
    if (sum(!xst)) {
      if (sum(!xst) == length(rawFiles)) {
        msg <- "None of your MS files could be located. We will not be able to copy them to the delivery folder!"
        Tsts$"MS files" <- 1
      } else { msg <- "Some of your MS files could be located, they will not be copied to the delivery folder!" }
      warning(msg)
    }
    w <- which(xst)
    ext <- gsub(".*\\.", "", rawFiles[w])
    for (xt in unique(ext)) {
      w2 <- which(ext == xt)
      eval(parse(text = paste0("copFun <- fs::", c("file", "dir")[(xt == "d")+1], "_copy")), envir = .GlobalEnv)
      cat(" - Copying MS raw files to\n\t\t", dir, "\n\n")
      for (fl in rawFiles[w[w2]]) {
        Tsts$"MS files"[[fl]] <- try(copFun(fl, dir, overwrite = FALSE), silent = TRUE)
      }
    }
  } else {
    if (!length(list.files(dir))) {
      file.copy(system.file("extdata", "MS_files.txt", package = "proteoCraft"), dir, overwrite = TRUE)
    }
  }
  # - 2 Input directories
  topSrchDir <- paste0(outDir, "/2_Search(es)")
  if (!dir.exists(topSrchDir)) { dir.create(topSrchDir, recursive = TRUE) }
  Tsts$"Search folder" <- list()
  if (writeSearch) {
    cat(" - Copying search results to raw files to\n\t\t", topSrchDir, "\n\n")
    lapply(1:length(inDirs), function(dir_i) { #dir_i <- 1 #dir_i <- 2
      rs <- FALSE
      indir <- inDirs[dir_i]
      indirNm <- gsub(".*/" , "", indir)
      fls1 <- list.files(indir, recursive = TRUE, full.names = TRUE) # Compare with what we want to copy
      doYouCopy <- length(fls1) > 0
      if (!doYouCopy) {
        warning(paste0("     NB: input directory \"", indir, "\" is currently empty!\n         did you clear it since starting this analysis?"))
      } else {
        flsTbl1 <- data.frame(File = fls1,
                              Name = basename(fls1),
                              Size = file.size(fls1))
        fls1 <- flsTbl1$Name
        kount <- 0
        destDirOK <- FALSE
        while (!destDirOK) {
          destDir <- paste0(topSrchDir, "/", indirNm, c("", paste0("v", as.character(kount)))[(kount > 0)+1])
          if (dir.exists(destDir)) {
            fls2 <- list.files(destDir, recursive = TRUE, full.names = TRUE) # Does it already contain stuff?
            if (length(fls2)) { # Yes? Then...
              flsTbl2 <- data.frame(File = fls2,
                                    Name = basename(fls1),
                                    Size = file.size(fls2))
              fls2 <- flsTbl2$Name
              destDirOK <- sum(!fls1 %in% fls2) == 0 # All the files should be here...
              if (destDirOK) {
                m12 <- match(fls1, fls2)
                destDirOK <- sum(file.size(flsTbl1$File) != file.size(flsTbl2$File[m12]),
                                 na.rm = TRUE) == 0 # ... and they should have the same size
              }
              if (!destDirOK) { kount <- kount+1 }
            } else {
              destDirOK <- doYouCopy <- TRUE
            }
          } else { destDirOK <- doYouCopy <- TRUE }
        }
        if (doYouCopy) {
          flsTbl1$New <- paste0(destDir, "/", flsTbl1$Name)
          flsTbl1$NewExists <- file.exists(flsTbl1$New)
          w <- which(!flsTbl1$NewExists)
          if (length(w)) {
            cat(" - Copying input data from\n\t\t", indir, "\n   to\n\t\t", destDir, "\n\n")
            tmpDir <- paste0(outDir, "/tmp_", dir_i)
            if (!dir.exists(tmpDir)) { dir.create(tmpDir, recursive = TRUE) }
            tst <- try(fs::dir_copy(indir, tmpDir), silent = TRUE)
            rs <- !"try-error" %in% class(Tsts$"Search folder")
            if (rs) { file.rename(paste0(tmpDir, "/", indirNm), destDir) }
          }
        }
      }
      return(rs)
    })
    Tsts$"Search folder - outcome" <- sum(vapply(Tsts$"Search folder", function(x) { "try-error" %in% class(x) }, TRUE)) == 0
    if (Tsts$"Search folder - outcome") {
      tmp <- c(paste0("These are the search engine", c("'s", "s'")[(length(unique(SearchSoft)) > 1)+1], " output files."),
               "Please save them to one of your own groups' shares then delete them from here to free up space.",
               "Once this is done please send us a confirmation email.",
               "",
               "Keep these files preciously as you may need to upload them as supporting data for publications.",
               "")
      write(tmp, paste0(topSrchDir, "/Search_results.txt"))
      lapply(1:length(inDirs), function(dir_i) { #dir_i <- 1
        tmpDir <- paste0(outDir, "/tmp_", dir_i)
        unlink(tmpDir, TRUE, TRUE)
      })
    }
  } else {
    if (!"Search_results.txt" %in% list.files(topSrchDir)) {
      tmp <- c("Output files from the search(es), as well as details about search parameters, are available from the facility upon request.",
               "")
      write(tmp, paste0(topSrchDir, "/Search_results.txt"))
    }
  }
  #
  # - 3 Post-processing directory
  # - 3a) remove empty directories
  tstDrs <- list.dirs(wd, full.names = TRUE, recursive = TRUE)
  tstDrs2 <- unique(dirname(list.files(wd, full.names = TRUE, recursive = TRUE)))
  tstDrs2 <- unique(unlist(lapply(strsplit(tstDrs2, "/"), function(x) {
    x <- unlist(x) 
    vapply(1:length(x), function(y) { paste(x[1:y], collapse = "/") }, "")
  })))
  tstDrs <- tstDrs[which(!tstDrs %in% tstDrs2)]
  if (length(tstDrs)) { for (dr in tstDrs) { unlink(dr) } } # Remove empty directories
  # - 3b) create final output directory
  procDir %<o% paste0(outDir, "/3_Post_processing_", Sys.Date())
  g <- grep("^3(\\.[0-9]+)?_", list.dirs(outDir, FALSE, FALSE), value = TRUE)
  kount <- length(g)
  if (kount) {
    procDir <- paste0(outDir, "/3.", kount, "_Post_processing_", Sys.Date())
  }
  #unloadNamespace("devtools")
  #unloadNamespace("usethis")
  #unloadNamespace("fs")
  # pak::pak("r-lib/fs")
  # Above ^ = Great function which seems to simplify package installation
  # compared with install.packages() and devtools::install_github()
  # - 3c) temporarily move .RData files out of the way
  tmpRDat <- grep("\\.RData$", list.files(wd, full.names = TRUE), value = TRUE) # Move all RData files out of the way for now
  tmpRDat <- grep("/Backup\\.RData$", tmpRDat, value = TRUE, invert = TRUE) # We do want to export the final Backup.RData file: it is large, but useful to have
  tmpRDat2 <- gsub(topattern(wd), projDir, tmpRDat)
  if (length(tmpRDat)) {
    fs::file_move(tmpRDat, tmpRDat2)
  }
  # - 3d) copy analysis results
  cat(" - Copying processing results from\n\t\t", wd, "\n   to\n\t\t", procDir, "\n\n")
  Tsts$"Data analysis" <- try(fs::dir_copy(wd, procDir, overwrite = FALSE), silent = TRUE)
  if ("try-error" %in% class(Tsts$"Data analysis")) {
    warning("Analysis results not copied to destination - usually this is a path length issue.\nYou will have to copy them manually.")
    Tsts$"Data analysis" <- FALSE
  } else {
    tmpDr <- paste0(outDir, gsub(".*/" , "/", wd))
    if (("character" %in% class(Tsts$"Data analysis"))&&(Tsts$"Data analysis" == tmpDr)) {
      Tsts$"Data analysis" <- file.rename(tmpDr, procDir)
      Tsts$"Data analysis" <- Tsts$"Data analysis" == procDir
      tmp <- grep("\\.RData$", list.files(procDir, all.files = TRUE, full.names = TRUE), value = TRUE)
      tmp <- grep("/Backup\\.RData$", tmp, value = TRUE, invert = TRUE) # We want to export the final Backup.RData file: it is large, but useful to have
      unlink(tmp) # Unlink the other RData files
      unlink(paste0(procDir, "/.RHistory"))
    }
  }
  # - 3e) move .RData files back
  if (length(tmpRDat)) { # Now put them back in their original place
    fs::file_move(tmpRDat2, tmpRDat)
  }
  dataDeliveryOk <- (is.logical(Tsts$"Data analysis"))&&(!is.na(Tsts$"Data analysis"))&&(Tsts$"Data analysis" == TRUE)
}
# - 4 Cleanup!
if (ok2Deliver&&dataDeliveryOk) {
  cleanUp <- c(TRUE, FALSE)[match(dlg_message("Should we cleanup the temporary folder? (parameter files will remain)", "yesno")$res, c("yes", "no"))]
  if (cleanUp) {
    drs <- list.dirs(wd, recursive = FALSE, full.names = TRUE)
    unlink(drs, TRUE, TRUE)
    #openwd(wd)
  }
}

# End logging:
if (scrptType == "withReps") { sink(NULL, type = "message") }
#close(logcon)
rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
invisible(parLapply(parClust, 1:N.clust, function(x) { rm(list = ls());gc() }))
rm(ReportCalls) # Temporary fix until I figure out how to fix the grphtype bug - I thought I had
setwd(wd); saveImgFun(BckUpFl) # Leave an ultimate backup in the temporary folder
#loadFun(BckUpFl)

# Also save a citations report
if (dataDeliveryOk) { setwd(procDir) } else { setwd(wd) }
if (!require(grateful)) {
  pak::pkg_install("grateful")
}
invisible(suppressMessages({
  a <- captureOutput(grateful::cite_packages(out.dir = ".", pkgs = "Session", quiet = TRUE))
})) # This thing won't shut up!

# Cleanup
locsFl <- paste0(homePath, "/Default_locations.xlsx")
locs <- openxlsx2::read_xlsx(locsFl)
archDirDflt <- locs$Path[match("Archive folder (searches)", locs$Folder)]
archDirDflt <- archDirDflt[which(dir.exists(archDirDflt))]
if ((dataDeliveryOk)&&(length(archDirDflt) == 1)) {
  inDirs2 <- grep(topattern(archDirDflt), inDirs, value = TRUE, invert = TRUE)
  L <- length(inDirs2)
  inDirs2Arch <- c()
  if (L == 1) {
    msg <- "Archive input (= search) folder?"
    archiveIndirs <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
    inDirs2Arch <- inDirs2[which(archiveIndirs)]
  }
  if (L > 1) { 
    msg <- paste0("Which input (= search) folder(s) would you like to archive?")
    opt <- gsub(".*/", ".../", inDirs2)
    nc <- max(c(nchar(opt), 200))
    opt <- vapply(opt, function(x) { paste(c(x, rep(" ", nc-nchar(x))), collapse = "") }, "")
    archiveIndirs <- dlg_list(opt, opt, TRUE, title = msg)$res
    inDirs2Arch <- inDirs2[match(archiveIndirs, opt)]
  }
  if (length(inDirs2Arch)) {
    lapply(inDirs2Arch, function(indir) { #indir <- inDirs2Arch[1]
      indirArch <- paste0(archDir, "/", gsub(".*/", "", indir))
      if (!dir.exists(indirArch)) { dir.create(indirArch, recursive = TRUE) }
      fls <- list.files(indir, recursive = TRUE#, all.files = TRUE
      )
      if (length(fls)) { # Checking because we may already have archived...
        clusterExport(parClust, list("fls", "indir", "indirArch"), envir = environment())
        parLapply(parClust, 1:length(fls), function(x) { #x <- 1
          fl <- fls[x]
          oldFl <- paste0(indir, "/", fl)
          nuFl <- paste0(indirArch, "/", fl)
          nuDr <- gsub("/[^/]+$", "", nuFl)
          if (!dir.exists(nuDr)) { dir.create(nuDr, recursive = TRUE) }
          fs::file_move(oldFl, nuFl)
        })
        write(c(paste0("The data in this folder was archived on ", Sys.Date(), "; its new location is:"),
                paste0("\t", indirArch),
                ""),
              paste0(indir, "/Archiving_log.txt"))
        drs <- list.dirs(indir, recursive = TRUE, full.names = TRUE)
        drs <- drs[which(drs != indir)]
        fls <- list.files(indir, recursive = TRUE)
        if ((length(drs))&&(length(fls) == 1)&&(fls == "Archiving_log.txt")) {
          for (dr in drs) { unlink(dr, TRUE, TRUE) }
        }
      }
      return(TRUE)
    })
  }
}

# Save final state of the environment
# This is done within the destination folder (outDir) because it will restart the session so has to be done last
# (this will interrupt the script flow so all commands queued after that are gone)
if (dataDeliveryOk) { dr <- procDir } else { dr <- wd }
setwd(dr)
pkgs <- gtools::loadedPackages()
dscrptFl <- paste0(dr, "/DESCRIPTION")
tmp <- paste0(do.call(paste, c(pkgs[, c("Name", "Version")], sep = " (")), ")")
tmp <- paste0("Depends: ", paste(tmp, collapse = ", "))
write(tmp, dscrptFl)
renv::snapshot(force = TRUE, prompt = FALSE, type = "explicit")
if ((exists("renv"))&&(renv)) { try(renv::deactivate(), silent = TRUE) }
setwd(wd)

# Done!
