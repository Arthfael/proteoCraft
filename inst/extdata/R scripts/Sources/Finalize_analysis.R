# Write Materials and Method to text file
setwd(wd)
tmp <- paste0("MatMet <- ", unlist(MatMetCalls$Calls))
Src <- paste0(wd, "/tmp.R")
write(tmp, Src)
MatMetFl <- paste0(wd, "/Materials and methods_WIP.docx")
tst <- try({
  source(Src)
  MatMet %<o% MatMet
  print(MatMet, target = MatMetFl)
}, silent = TRUE)
if ("try-error" %in% class(tst)) {
  warning("Couldn't write materials and methods template, investigate...")
} else {
  system(paste0("open \"", MatMetFl, "\""))
  dlg_message("Check materials and methods text, make any necessary edits, then save and click ok", "ok")
}

# Write pdf report
if (scrptType == "withReps") {
  # Write report to Word file
  tmp <- paste0("Report <- ", unlist(ReportCalls$Calls))
  Src <- paste0(wd, "/tmp.R")
  write(tmp, Src)
  tst <- try({
    source(Src)
    Report %<o% Report
    print(Report, target = paste0(wd, "/Workflow control/Analysis report.docx"))
  }, silent = TRUE)
  if ("try-error" %in% class(tst)) { warning("Couldn't write pdf report, investigate...") }
  unlink(Src)
  #system(paste0("open \"", wd, "/Workflow control/Analysis report.docx"))
}

# Save session info
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
if (scrptType == "withReps") { dirlist <- unique(c(dirlist, dir)) }
writeLines(capture.output(sessionInfo()), paste0(dir, "/sessionInfo.txt"))

# Copy script itself and guide to the results there:
#file.copy(system.file("extdata/R scripts", "Regulation analysis - master script.R", package = "proteoCraft"), wd, overwrite = TRUE)
file.copy(ScriptPath, wd, overwrite = TRUE)
file.copy(system.file("extdata", paste0("Guide to the Results (", c("", "no ")[match(scrptType,
                                                                                     c("withReps", "noReps"))],
                                        "replicates).docx"), package = "proteoCraft"), wd, overwrite = TRUE)

## Copy final results to the delivery folder
Tsts <- list()
# - 1/ MS files
dir <- paste0(outdir, "/1_MS_files")
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
    eval(parse(text = paste0("copFun <- fs::", c("file", "dir")[(xt == "d")+1], "_copy")))
    for (fl in rawFiles[w[w2]]) {
      Tsts$"MS files"[[fl]] <- try(copFun(fl, dir, overwrite = FALSE), silent = TRUE)
    }
  }
} else {
  file.copy(system.file("extdata", "MS_files.txt", package = "proteoCraft"), dir, overwrite = TRUE)
}
# - 2 Input directory
dir <- paste0(outdir, "/2_", names(SearchSoft), "_search")
nuDir <- dir.exists(dir)
doICopy <- writeSearch
if (writeSearch) {
  kount <- 0
  while (nuDir) {
    tst1 <- list.files(dir, recursive = TRUE)
    nuDir <- length(tst1) > 0
    if (nuDir) {
      tst2 <- list.files(indir, recursive = TRUE)
      if ((sum(!tst2 %in% tst1) == 0)&&(length(tst2) == length(tst1))) {
        m <- match(tst2, tst1)
        nuDir <- (sum(file.size(paste0(indir, "/", tst2)) != file.size(paste0(dir, "/", tst1[m]))) != 0)
        doICopy <- nuDir
      }
    }
    if (nuDir) {
      kount <- kount + 1
      dir <- paste0(outdir, "/2.", kount, "_", names(SearchSoft), "_search")
      nuDir <- dir.exists(dir)
    }
  }
  if (doICopy) {
    Tsts$"Search folder" <- try(fs::dir_copy(indir, outdir, overwrite = FALSE), silent = TRUE)
    if (!"try-error" %in% class(Tsts$"Search folder")) { file.rename(paste0(outdir, gsub(".*/" , "/", indir)), dir) }
  }
} else {
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  tst1 <- list.files(dir, recursive = TRUE)
  if (!length(tst1)) {
    tmp <- paste0("Output files from the ", names(SearchSoft),
                  " search, as well as details about search parameters, are available from the facility upon request.\n")
    write(tmp, paste0(dir, "/Search_results.txt"))
  } else {
    tmp <- paste0("These are ", names(SearchSoft),
                  "' O's output files. Please save them to one of your own groups' shares then delete them to free up space.\n")
    write(tmp, paste0(dir, "/Search_results.txt"))
  }
}
# - 3 Post-processing directory
# - 3a) remove empty directories
tstDrs <- list.dirs(wd, full.names = TRUE, recursive = TRUE)
tstDrs2 <- unique(dirname(list.files(wd, full.names = TRUE, recursive = TRUE)))
tstDrs2 <- unique(unlist(lapply(strsplit(tstDrs2, "/"), function(x) {
  x <- unlist(x) 
  sapply(1:length(x), function(y) { paste(x[1:y], collapse = "/") })
})))
tstDrs <- tstDrs[which(!tstDrs %in% tstDrs2)]
if (length(tstDrs)) { for (dr in tstDrs) { unlink(dr) } } # Remove empty directories
# - 3b) create final output directory
procdir %<o% paste0(outdir, "/3_Post_processing_", Sys.Date())
g <- grep("^3(\\.[0-9]+)?_", list.dirs(outdir, FALSE, FALSE), value = TRUE)
kount <- length(g)
if (kount) {
  procdir <- paste0(outdir, "/3.", kount, "_Post_processing_", Sys.Date())
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
Tsts$"Data analysis" <- try(fs::dir_copy(wd, procdir, overwrite = FALSE), silent = TRUE)
if ("try-error" %in% class(Tsts$"Data analysis")) {
  warning("Analysis results not copied to destination - usually this is a path length issue.\nYou will have to copy them manually.")
  Tsts$"Data analysis" <- FALSE
} else {
  tmpDr <- paste0(outdir, gsub(".*/" , "/", wd))
  if (("character" %in% class(Tsts$"Data analysis"))&&(Tsts$"Data analysis" == tmpDr)) {
    Tsts$"Data analysis" <- file.rename(tmpDr, procdir)
    tmp <- grep("\\.RData$", list.files(procdir, full.names = TRUE), value = TRUE)
    tmp <- grep("/Backup\\.RData$", tmp, value = TRUE, invert = TRUE) # We want to export the final Backup.RData file: it is large, but useful to have
    unlink(tmp) # Unlink the other RData files
    unlink(paste0(procdir, "/.RHistory"))
  }
}
# - 3e) move .RData files back
if (length(tmpRDat)) { # Now put them back in their original place
  fs::file_move(tmpRDat2, tmpRDat)
}
# - 4 Cleanup!
if ((is.logical(Tsts$"Data analysis"))&&(!is.na(Tsts$"Data analysis"))&&(Tsts$"Data analysis" == TRUE)) {
  cleanUp <- c(TRUE, FALSE)[match(dlg_message("Should we cleanup the temporary folder? (parameter files will remain)", "yesno")$res, c("yes", "no"))]
  if (cleanUp) {
    drs <- list.dirs(wd, recursive = FALSE, full.names = TRUE)
    unlink(drs, TRUE, TRUE)
    #openwd(wd)
  }
}
