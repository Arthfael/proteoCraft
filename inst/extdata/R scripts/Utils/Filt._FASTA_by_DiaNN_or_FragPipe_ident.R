# Filter a Fasta with the results from a 1st unfiltered DIA-NN search, then predict a library from it.
# NB: Un-necessary for FragPipe: that fasta is automatically saved in the during the search!

require(reshape2)
require(svDialogs)
require(proteoCraft)
wd <- "...Search_Folder"

Inputs <- c("DIA-NN", "FragPipe")
Input <- dlg_list(Inputs, Inputs[1], title = "Which software search was used?")$res
Modes <- c("Write output fasta", "Just get a report on number of PSMs/proteins (no)")
Mode <- dlg_list(Modes, Modes[1], title = "What do you want to do?")$res
# DIA-NN
if (Input == Inputs[1]) {
  # Detect DIA-NN exe
  DIANNdir <- grep("DIA-NN", list.dirs("C:/", recursive = FALSE), value = TRUE)
  msg <- "Select DiaNN.exe executable"
  if (length(DIANNdir) == 1) {
    DIANNdir2 <- grep("DIA-NN", list.dirs(DIANNdir, recursive = FALSE), value = TRUE)
    stopifnot(length(DIANNdir2) > 0)
    DIANNexe <- paste0(normalizePath(DIANNdir2), "\\DiaNN.exe")
    stopifnot(sum(file.exists(DIANNexe)) > 0)
    if (length(DIANNexe) > 1) {
      #DIANNexe <- DIANNexe[order(file.mtime(DIANNexe), decreasing = TRUE)[1]]
      #DIANNexe <- normalizePath(choose.files(paste0(normalizePath(DIANNdir), "\\*.exe"), msg, multi = FALSE))
      DIANNexe <- rstudioapi::selectFile(msg,
                                         path = paste0(DIANNdir, "/*.exe"),
                                         filter = "DiaNN exe (*.exe)")
    }
  } else {
    #DIANNexe <- normalizePath(choose.files("C:/*.exe", msg, multi = FALSE))
    DIANNexe <- rstudioapi::selectFile(msg,
                                       path = "C:/*.exe",
                                       filter = "DiaNN exe (*.exe)")
  }
  if ((length(DIANNexe) == 1)&&(file.exists(DIANNexe))) {
    DIANNexe <- gsub("\\\\DIA-NN\\.exe$", "\\\\DiaNN.exe", DIANNexe) # In case you selected the wrong executable in the folder
    DIANNexe <- paste0("\"", DIANNexe, "\"")
    #
    msg <- "Select a DIA-NN log file"
    #DIANNlogFl <- choose.files(paste0(wd, "/*.log.txt"), msg, multi = FALSE)
    DIANNlogFl <- rstudioapi::selectFile(msg,
                                         path = paste0(wd, "/*.log.txt"),
                                         filter = "DiaNN log file (*.log.txt)")
    wd <- dirname(DIANNlogFl)
    setwd(wd)
    #
    cat("Getting list of identified proteins...\n")
    DIANNlog <- readLines(DIANNlogFl)
    DIANNcall <- grep("^diann.exe ", DIANNlog, value = TRUE)
    DIANNargs <- grep("^--", unlist(strsplit(gsub(" --", " Arg--", DIANNcall), " Arg")), value = TRUE)
    rprtFl <- gsub("\\\\", "/", gsub("^--out ", "", grep("^--out ", DIANNargs, value = TRUE)))
    rprt <- try(read.delim(rprtFl), silent = TRUE)
    if ("try-error" %in% class(rprt)) {
      cat("DiaNN file not found at the expected location, trying relative file path\n")
      rprtFl <- paste0(dirname(DIANNlogFl), "/", basename(rprtFl))
      if (file.exists(rprtFl)) {
        cat("Success!\n")
        rprt <- read.delim(rprtFl)
      } else { stop("Nope... I have no idea where that report went!") }
    }
    #libFl <- gsub("\\\\", "/", gsub("^--out-lib ", "", grep("^--out-lib ", DIANNargs, value = TRUE)))
    FastaFl <- gsub("\\\\", "/", gsub("^--fasta ", "", grep("^--fasta ", DIANNargs, value = TRUE)))
    #FastaFl <- gsub(" ", "_", FastaFl) # Can be useful in old cases - I have changed some paths to remove spaces
    prt <- unique(unlist(strsplit(rprt$Protein.Ids, ";")))
    if (!length(FastaFl)) { stop("The selected DIA-NN run was not set to use a FASTA!") } else {
      cat("Parsing fasta file...\n")
      dbs <- lapply(FastaFl, Format.DB, Unique = FALSE) # We do not want to filter duplicates here!
      # Otherwise we risk missing some accessions!
      db <- plyr::rbind.fill(dbs)
      dbFlt <- db[which(db$`Protein ID` %in% prt),]
      # Instead, we will filter for duplicates here
      seq <- unique(dbFlt$Sequence)
      if (length(seq) < nrow(dbFlt)) {
        wY <- match(seq, dbFlt$Sequence)
        wN <- which(!c(1:nrow(dbFlt)) %in% wY)
        cat(paste0("", length(wN), " duplicate sequences were removed.\n"))
        dbFlt <- dbFlt[wY, , drop = FALSE]
      }
      cat(paste0("DIA-NN identified ", nrow(rprt), " PSMs covering ", nrow(dbFlt), " proteins\n"))
      if (Mode == Modes[1]) {
        fasta2 <- rep("", nrow(dbFlt)*3)
        fasta2[(1:nrow(dbFlt))*3-2] <- dbFlt$Header
        fasta2[(1:nrow(dbFlt))*3-1] <- dbFlt$Sequence
        nm <- grep("cRAPome", basename(FastaFl), invert = TRUE, value = TRUE)
        nm <- dlg_input("Enter a file name for the new fasta database:", gsub("\\.fas(ta(\\.fas)?)?$", "_filt", nm))$res
        write(fasta2, paste0(wd, "/", nm, ".fasta"))
        cat(paste0("Filtered fasta saved to \"", wd, "/", nm, ".fasta\"\nIf using in DIA-NN, do not forget to add contaminants again!"))
        #
        DIANNargs2 <- grep("^--((f)|(lib)) ", DIANNargs, invert = TRUE, value = TRUE)
        gO <- grep("^--out ", DIANNargs2)
        tmpO <- paste0("--out \"", gsub("/", "\\\\", wd), "\\Report.tsv\"")
        if (length(gO)) { DIANNargs2[gO] <- tmpO } else { DIANNargs2 <- c(DIANNargs2, tmpO) }
        gOL <- grep("^--out-lib ", DIANNargs2)
        tmpOL <- paste0("--out-lib \"", gsub("/", "\\\\", wd), "\\Report_filt_lib.tsv\"")
        if (length(gOL)) { DIANNargs2[gOL] <- tmpOL } else { DIANNargs2 <- c(DIANNargs2, tmpOL) }
        gF <- grep("^--fasta ", DIANNargs2)
        DIANNargs2[gF] <- paste0("--fasta \"", gsub("/", "\\\\", wd), "\\", nm, ".fasta", "\"")
        if (!"--fasta-search" %in% DIANNargs2) { DIANNargs2 <- c(DIANNargs2, "--fasta-search")}
        if (!"--predictor" %in% DIANNargs2) { DIANNargs2 <- c(DIANNargs2, "--predictor")}
        cmd <- c(DIANNexe, DIANNargs2)
        cmd <- paste(cmd, collapse = " ")
        #cat(cmd)
        Run <- c(TRUE, FALSE)[match(dlg_message("Do you want to run the DIA-NN command now?\nThis will create a library from the filtered fasta.",
                                                "yesno")$res, c("yes", "no"))]
        if (Run) { system(cmd) }
      }
    }
  }
}
# FragPipe
if (Input == Inputs[2]) {
  msg <- "Choose FragPipe workflow file"
  #FP_Workflow <- normalizePath(choose.files(paste0(wd, "/*.workflow"), msg, multi = FALSE), winslash = "/")
  FP_Workflow <- rstudioapi::selectFile(msg,
                                        path = paste0(wd, "/*.workflow"),
                                        filter = "FragPipe workflow file (*.workflow)")
  stopifnot(length(FP_Workflow) > 0)
  wd <- dirname(FP_Workflow)
  msg <- "Choose FragPipe manifest"
  #FP_Manifest <- normalizePath(choose.files(paste0(wd, "/*.fp-manifest"), , multi = FALSE), winslash = "/")
  FP_Manifest <- rstudioapi::selectFile(msg,
                                        path = paste0(wd, "/*.fp-manifest"),
                                        filter = "FragPipe MANIFEST file (*.fp-manifest)")
  stopifnot(length(FP_Manifest) > 0)
  FP_Workflow <- readLines(FP_Workflow)
  FP_Dir <- gsub("\\\\", "", gsub("\\\\\\\\", "/", gsub("^workdir=", "", grep("^workdir=", FP_Workflow, value = TRUE))))
  if (!dir.exists(FP_Dir)) {
    stop(paste0("FragPipe's output directory was not found: \"", FP_Dir, "\""))
  }
  FastaFl <- grep("^database\\.db-path=", FP_Workflow, value = TRUE)
  FastaFl <- gsub("^database\\.db-path=", "", FastaFl) 
  FastaFl <- gsub("\\\\", "", gsub("\\\\\\\\", "/", FastaFl))
  #
  # Parse FragPipe's samples manifest
  FP_Manifest <- suppressWarnings(read.delim(FP_Manifest, header = FALSE))
  colnames(FP_Manifest) <- c("Path", "Experiment", "Replicate", "Data type")
  stopifnot(sum(file.exists(FP_Manifest$Path)) == nrow(FP_Manifest))
  FP_Manifest$Path <- gsub("\\\\", "/", FP_Manifest$Path)
  FP_Manifest$"File name" <- gsub("\\.(raw|mz(X?ML|BIN)|mgf|d)$", "", gsub(".+/", "", FP_Manifest$Path))
  Exp <- Samples <- FP_Manifest$Experiment
  Rep <- FP_Manifest$Replicate
  if (length(Rep)) {
    Samples <- paste0(Exp, "_", Rep, "_", Rep)
    # Rep seems to be duplicated somehow... I am testing if this is always the case...
  }
  FP_PSMs <- paste0(FP_Dir, "/", Samples, "/psm.tsv")
  FP_Fas <- paste0(FP_Dir, "/", Samples, "/protein.fas")
  tst <- file.exists(FP_PSMs)
  w <- which(!tst); l <- length(w)
  if (l) {
    SamplesNo <- unique(SamplesNo[w]); l <- length(SamplesNo)
    if (l > 1) { SamplesNo <- paste0(paste(SamplesNo[1:(l-1)], collapse = ", "), " and ", SamplesNo[l]) }
    plur <- c("", "s")[(l > 1)+1]
    warning(paste0("The \"psm.tsv\" file", plur, " for sample", plur, " ", SamplesNo,
                   " could not be found, skipping missing sample", plur, "!"))
    FP_PSMs <- FP_PSMs[which(tst)]
    FP_Fas <- FP_Fas[which(tst)]
    Samples <- Samples[which(tst)]
    FP_Manifest <- FP_Manifest[which(tst),]
    Exp <- Exp[which(Exp %in% FP_Manifest$Experiment)]
    Rep <- Rep[which(Rep %in% FP_Manifest$Replicate)]
  }
  tst <- file.exists(FP_Fas)
  w <- which(!tst); l <- length(w)
  if (l) {
    SamplesNo <- unique(SamplesNo[w]); l <- length(SamplesNo)
    if (l > 1) { SamplesNo <- paste0(paste(SamplesNo[1:(l-1)], collapse = ", "), " and ", SamplesNo[l]) }
    plur <- c("", "s")[(l > 1)+1]
    warning(paste0("The \"protein.fas\" file", plur, " for sample", plur, " ", SamplesNo,
                   " could not be found, skipping missing sample", plur, "!"))
    FP_PSMs <- FP_PSMs[which(tst)]
    FP_Fas <- FP_Fas[which(tst)]
    Samples <- Samples[which(tst)]
    FP_Manifest <- FP_Manifest[which(tst),]
    Exp <- Exp[which(Exp %in% FP_Manifest$Experiment)]
    Rep <- Rep[which(Rep %in% FP_Manifest$Replicate)]
  }
  Exp <- unique(Exp)
  Samples <- unique(Samples)
  Rep <- unique(Rep)
  if ((Mode == Modes[1])&&(length(FP_Fas) == 1)) {
    cat(paste0("A single \"protein.fas\" file was written by FragPipe, no processing required:\n -> ",
               FP_Fas, "\n"))
    Mode <- Modes[2]
  }
  cat("Parsing PSM files...\n")
  PSMs <- lapply(FP_PSMs, function(FP_PSM) {
    x <- read.delim(FP_PSM, check.names = FALSE)
    x$"PSM file" <- FP_PSM
    return(x)
  })
  PSMs <- plyr::rbind.fill(PSMs)
  cat("Parsing filtered fasta files...\n")
  dbFlts <- lapply(FP_Fas, Format.DB, Unique = FALSE) # We do not want to filter duplicates here!
  # Otherwise we risk missing some accessions!
  # Filter for unicity
  dbFlt <- plyr::rbind.fill(dbFlts)
  tst <- aggregate(1:nrow(dbFlt), list(dbFlt$Sequence), function(x) { min(x) })
  dbFlt <- dbFlt[tst$x,]
  #
  cat(paste0("FragPipe identified ", nrow(PSMs), " PSMs covering ", nrow(dbFlt), " proteins\n"))
  if (Mode == Modes[1]) {
    nm <- grep("cRAPome", basename(FastaFl), invert = TRUE, value = TRUE)
    nm <- dlg_input("Enter a file name for the new fasta database:", gsub("\\.fas(ta(\\.fas)?)?$", "_filt", nm))$res
    write(fasta2, paste0(wd, "/", nm, ".fasta"))
    cat(paste0("Filtered fasta saved to \"", wd, "/", nm, ".fasta\"\nIf using in FragPipe, do not forget to add contaminants and decoys again!"))
  }
}  

