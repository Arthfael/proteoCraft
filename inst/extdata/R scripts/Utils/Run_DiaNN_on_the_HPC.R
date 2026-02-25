# Header
# This script will
# - interactively create a SLURM script for running DiaNN jobs on the HPC
# - start the job (this may send a "Houston, we have a problem" error,
#   in which case you can run:
#      ssh_exec_wait(sshsess, batchcmd)
#   to start the job)
#

packs <- c("svDialogs", "ssh")
for (pack in packs) {
  if (!require(pack, character.only = TRUE)) { install.packages(pack) }
  require(pack, character.only = TRUE)
}
if (!require(unimod)) { suppressMessages(devtools::install_github("rformassspectrometry/unimod")) }
require(unimod)

# A version of topattern is given here (so installing the proteoCraft package is not necessary):
if (!require(proteoCraft)) {
  topattern <- function(x, start = TRUE, end = FALSE, collapse = "|") {
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
    if ((length(x) > 1)&&(collapse != FALSE)) { x <- paste(x, collapse = collapse) }
    return(x)
  }
}

# Initiate ssh session and check existing jobs:
if ((!exists("sshsess"))||(class(sshsess) != "ssh_session")||(!exists("usrNm"))) {
  hstNm <- dlg_input("Enter host name")$res
  usrNm <- dlg_input("Enter user name")$res
  sshsess <- try(ssh_connect(paste0(usrNm, "@", hstNm)), silent = TRUE)
  stopifnot(class(sshsess) == "ssh_session")
  user <- gsub("^connected: |@.+", "", as.character(capture.output(sshsess))[2])
  stopifnot(user == usrNm)
  print(sshsess)
  ssh_exec_wait(sshsess, "source /etc/profile.d/modules.sh")
}
#shell <- "source /etc/profile.d/10_lsb-release-env-vars.sh ; source /etc/profile.d/99_load_slurm.sh ; source /etc/profile.d/modules.sh ; "
shell <- "source /etc/profile.d/10_lsb-release-env-vars.sh ; source /etc/profile.d/modules.sh ; "
print("Checking ongoing jobs:")
queuecmd <- paste0(shell, "squeue --user=", user, " -o \"%.7i %.8u %.8j %.9a %.9P %.9N %.7T %.10M %.9L %.4C %.7m %.6D %R\"")
#queuecmd <- paste0(shell, "squeue --user=", user, " -l")
#cat(queuecmd)
ssh_exec_wait(sshsess, queuecmd)
# Function to disconnect the session
sshDisc <- function(session = sshsess) {
  try(ssh_disconnect(sshsess), silent = TRUE)
  if (exists("sshsess", .GlobalEnv)) { rm(sshsess, .GlobalEnv) }
}
rm(list = ls()[which(!ls() %in% c("sshsess", "shell", "queuecmd", "sshDisc", "BATCH", "topattern", "user"))])

EMAIL <- dlg_input("Enter email to which SLURM will be sending messages", paste0(user, "@"))$res
# Function from https://www.r-bloggers.com/2012/07/validating-email-adresses-in-r/
isValidEmail <- function(x) { grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x), ignore.case = TRUE) }
while(!isValidEmail(EMAIL)) { EMAIL <- dlg_input("Invalid email address, try again!", paste0(user, "@"))$res }

# Choose the name of the batch folder
if (exists("BATCH")) { def <- BATCH } else { def <- "My_batch" }
BATCH <- gsub(" ", "_", dlg_input("Enter the name of your batch folder", def)$res)
# Job name (currently the same as that of the batch folder)
JOBNAME <- BATCH #gsub(" ", "_", dlg_input("Enter the name of your job", BATCH)$res)
# Get the location where the script will be saved:
LINLOC <- capture.output(ssh_exec_wait(sshsess, "pwd"))[1]
LINLOC <- gsub("/+$", "", gsub("/+", "/", paste0("/", LINLOC)))
BATCHDIR <- paste0(LINLOC, "/", BATCH)

#ssh_exec_wait(sshsess, paste0("env > ", BATCHDIR, "/e0.txt"))
ssh_exec_wait(sshsess, paste0("cd \"", LINLOC, "\""))
out <- capture.output(ssh_exec_wait(sshsess, "ls -l | grep '^d'"))
if (length(grep(paste0(" ", BATCH, "$"), out)) == 0) { ssh_exec_wait(sshsess, paste0("mkdir -p ", BATCH)) }
#
Modes <- setNames(c("Bruker d files (folders)                                   ",
                    "Any other type of MS file                                  ",
                    "None (library building)                                    "),
                  c("Bruker", "MS files", "Library"))
Mode <- dlg_list(Modes, Modes[1], FALSE, "Which file types do you want to analyse?")$res
Mode <- names(Modes)[match(Mode, Modes)]
if (Mode == names(Modes)[1]) {
  FILES <- c()
  Moar <- TRUE
  dflt <- getwd()
  while (Moar) {
    msg <- c("Select input Bruker .d folder.", "Select another Bruker .d folder, or escape to stop adding more.")[length(FILES > 0)+1]
    fl <- rstudioapi::selectDirectory(msg, path = dflt)
    Moar <- length(fl) > 0
    if (Moar) {
      dflt <- gsub("/[^/]+$", "", fl)
      FILES <- unique(c(FILES, fl))
    }
  }
}
if (Mode == names(Modes)[2]) {
  FILES <- choose.files("", "Select input MS files (or escape for library-building only mode - no file search)")
}
NTHREADS <- NA
while (is.na(suppressWarnings(as.integer(NTHREADS)))) { NTHREADS <- dlg_input("How many threads would you like to request?", max(c(12, length(FILES))))$res }
#
filt <- matrix(c("fasta file", "fas file", "fa file", "faa file", "fasta.fas file", "txt (fasta) file", ".fasta", ".fas", ".fa", ".faa", "fasta.fas", ".txt"),
               ncol = 2)
Moar <- TRUE
FASTAS <- c()
dflt2 <- gsub("/[^/]+$", "/*.fasta", dflt)
msg <- "Select fasta files (or Esc. for lib-build)"
while (Moar) {
  fl <- normalizePath(choose.files(dflt2, msg, filters = filt), winslash = "/")
  Moar <- length(fl)
  if (Moar) {
    FASTAS <- unique(c(FASTAS, fl))
    msg <- "Select more fasta files? (Esc. to cancel)"
    dflt2 <- gsub("/[^/]+$", "/*.fasta", fl)
  }
}
reAnnot <- length(FASTAS) > 0
if (reAnnot) { FASTAS <- normalizePath(FASTAS, winslash = "/") }
#
filt <- matrix(c("speclib library file", "tsv library file", "*.speclib", "*.tsv"), ncol = 2)
LIBRARIES <- choose.files(gsub("/[^/]+$", "/*.tsv", dflt), "Select .speclib or .tsv library files (or escape for library-free mode)", filters = filt)
LibFree <- length(LIBRARIES) == 0
if (!LibFree) { LIBRARIES <- normalizePath(LIBRARIES, winslash = "/") }
FDR <- NA
while (is.na(suppressWarnings(as.numeric(FDR)))) { FDR <- dlg_input("FDR threshold?", 0.01)$res }

# Outputs - will be saved in the same folder
hpcFls <- capture.output(ssh_exec_wait(sshsess, paste0("ls ", BATCHDIR)))
i <- 0
while ((i == 0)||(RPRTNM %in% hpcFls)) {
  i <- i +1
  RPRTNM <- paste0(dlg_input("Enter DiaNN output report name", paste0("report", i))$res)
}
i <- 0; RPRTLIBNM <- NA
while ((i == 0)||((length(RPRTLIBNM) > 0)&&(RPRTLIBNM %in% hpcFls))) {
  i <- i +1
  RPRTLIBNM <- paste0(dlg_input("Enter DiaNN output library name", paste0(RPRTNM, "-lib"))$res)
}
# Ranges and mass accuracies
MZ1 <- MZ2 <- PEPL <- MS1Acc <- MS2Acc <- Z <- MISSES <- NA
while (sum(is.na(suppressWarnings(as.numeric(MZ1))))) { MZ1 <- as.character(as.numeric(unlist(strsplit(dlg_input("Enter precursor m/z range:", "200 - 1600")$res, " *- *")))) }
while (sum(is.na(suppressWarnings(as.numeric(MZ2))))) { MZ2 <- as.character(as.numeric(unlist(strsplit(dlg_input("Enter fragment m/z range:", "200 - 1800")$res, " *- *"))))}
while (sum(is.na(suppressWarnings(as.integer(Z))))) { Z <- as.character(as.numeric(unlist(strsplit(dlg_input("Enter precursor charge range:", "1 - 4")$res, " *- *"))))}
while (sum(is.na(suppressWarnings(as.integer(PEPL))))) { PEPL <- as.character(as.numeric(unlist(strsplit(dlg_input("Enter peptide length range:", "7 - 30")$res, " *- *")))) }
while (sum(is.na(suppressWarnings(as.numeric(MS1Acc))))) { MS1Acc <- as.character(as.numeric(dlg_input("Enter MS1 mass accuracy:", "20")$res)) }
while (sum(is.na(suppressWarnings(as.numeric(MS2Acc))))) { MS2Acc <- as.character(as.numeric(dlg_input("Enter MS2 mass accuracy:", MS1Acc)$res)) }
# Missed cleavages
while (sum(is.na(suppressWarnings(as.integer(MISSES))))) { MISSES <- dlg_input("Enter max. missed cleavages:", 1)$res }
# PTMs
UniMod <- unimod::modifications
CommonMods <- data.frame(Name = c("Acetylation (protein N-term)", "M oxidation", "deamidation NQ", "Gln->pyroGlu", "Phospho (STY)", "GlyGly"),
                         Site = c("*n", "M", "NQ", "Q", "STY", "K"),
                         UniMod = paste0("UniMod:", c(1, 35, 7, 28, 21, 121)))
opts <- c(apply(CommonMods, 1, paste, collapse = " / "), "... add other(s) not in this list")
VARMODS <- dlg_list(opts, opts[1:4], TRUE, "Select variable modifications to include:")$res
Moar <- ("... add other(s) not in this list" %in% VARMODS)
VARMODS <- VARMODS[which(VARMODS != "... add other(s) not in this list")]
if (length(VARMODS)) {
  VARMODS <- CommonMods[match(VARMODS, opts),]
}
if (Moar) {
  VARMODS2 <- dlg_input("Enter UniMod IDs and sites (Amino acid single-letter, n and *n for peptide and protein N-term), separate different mods with \",\"", "1 K, 34 K")$res
  VARMODS2 <- unlist(strsplit(VARMODS2, " *, *"))
  VARMODS2 <- data.frame(UniMod = paste0("UniMod:", gsub(" .*", "", VARMODS2)),
                         Site = gsub("[^ ]+ +", "", VARMODS2))
  VARMODS2$Name <- UniMod$Name[match(gsub("^UniMod:", "", VARMODS2$UniMod), UniMod$UnimodId)]
  if (class(VARMODS) == "data.frame") { VARMODS <- rbind(VARMODS, VARMODS2) } else { VARMODS <- VARMODS2 }
}
inclVarMods <- class(VARMODS) == "data.frame"
if (inclVarMods) {
  VARMODS$MassShift <- as.character(apply(VARMODS[, c("UniMod", "Site")], 1, function(x) { #x <- VARMODS[1, c("UniMod", "Site")]
    x[[1]] <- gsub("^UniMod:", "", x[[1]])
    if (x[[2]] %in% c("n", "*n")) {
      if (x[[2]] == "n") { x[[2]] <- "N-term"; pos <- "Any N-term" }
      if (x[[2]] == "*n") {x[[2]] <- "N-term";  pos <- "Protein N-term" }
      m <- which((UniMod$UnimodId == x[[1]])&(UniMod$Site == x[[2]])&(UniMod$Position == pos))
    } else {
      # If there are several amino acids and they are grouped, there is only one mass shift
      x[[2]] <- substr(x[[2]], 1, 1)
      m <- which((UniMod$UnimodId == x[[1]])&(UniMod$Site == x[[2]]))
    }
    return(UniMod$MonoMass[m[1]])
  }))
  VARMODS <- apply(VARMODS[, c("UniMod", "MassShift", "Site")], 1, paste, collapse = ",")
  NVARMODS <- dlg_input("Max. variable modifications per peptide:", 3)$res
  while (is.na(suppressWarnings(as.integer(NVARMODS)))) { NVARMODS <- dlg_input("Max. variable modifications per peptide (enter a valid number!):", 3)$res }
}

# Create DiaNN call
DIANNcall <- paste0("srun --cpu_bind=verbose /usr/bin/time -v diann")
if (Mode != Modes[3]) { DIANNcall <- paste0(DIANNcall, paste0(" --f \"${running_dir}/", basename(FILES), "\"", collapse = "")) }
#if (!LibFree) {
  # Somehow if doing library-free it appears you should still provide ' --lib "" ' as argument...
  DIANNcall <- paste0(DIANNcall, paste0(" --lib ", c("\"\"", paste0("\"${running_dir}/", basename(LIBRARIES), "\""))[(length(LIBRARIES)>0)+1], collapse = ""))
#}
DIANNcall <- paste0(DIANNcall, " --threads $SLURM_CPUS_PER_TASK --verbose 1 --out \"${running_dir}/", RPRTNM, ".tsv\" --qvalue ", FDR, " --matrices")
if (length(RPRTLIBNM)) { DIANNcall <- paste0(DIANNcall, "--out-lib \"${running_dir}/", RPRTLIBNM, ".tsv\"") }
DIANNcall <- paste0(DIANNcall, " --gen-spec-lib --predictor")
if (reAnnot) { DIANNcall <- paste0(DIANNcall, paste0(" --fasta \"${running_dir}/", basename(FASTAS), "\"", collapse = ""), " --fasta-search") }
DIANNcall <- paste0(DIANNcall, " --min-fr-mz ", MZ2[1], " --max-fr-mz ", MZ2[2], " --met-excision --cut K*,R* --missed-cleavages ", MISSES,
                    " --min-pep-len ", PEPL[1], " --max-pep-len ", PEPL[2], " --min-pr-mz ", MZ1[1], " --max-pr-mz ", MZ1[2],
                    " --min-pr-charge ", Z[1], " --max-pr-charge ", Z[2])
tstV1 <- length(VARMODS)
tstV2 <- grepl("UniMod:7,", VARMODS)+1
if (tstV1) {
  DIANNcall <- paste0(DIANNcall, " --var-mods ", NVARMODS,
                      paste0(" --var-mod ", VARMODS, c(paste0(" --monitor-mod ", gsub(",.*", "", VARMODS)), "")[tstV2], collapse = ""))
}
DIANNcall <- paste0(DIANNcall, " --reanalyse --relaxed-prot-inf --smart-profiling --peak-center --no-ifs-removal --min-fr-mz ", MZ2[1], " --threads ", NTHREADS)
if (length(MS1Acc)) { DIANNcall <- paste0(DIANNcall, " --mass-acc-ms1 ", MS1Acc) }
if (length(MS2Acc)) { DIANNcall <- paste0(DIANNcall, " --mass-acc ", MS2Acc) }

# Upload MS files + libraries
if (Mode %in% names(Modes)[1:2]) {
  flsSizes <- list()
  if (Mode == names(Modes)[1]) {
    require(RCurl)
    fls <- FILES[which(dir.exists(FILES))]
    tst <- sapply(fls, function(fl) { #fl <- FILES[1]
      res <- gsub(".*/", "", fl) %in% hpcFls
      if (res) {
        hpcFls2 <- capture.output(ssh_exec_wait(sshsess, paste0("ls ", gsub(topattern(dirname(fl)), BATCHDIR, fl), " -R")))
        w <- which(hpcFls2 %in% c("", "[1] 0"))
        l <- length(w); stopifnot(l > 0)
        hpcFls2 <- sapply(1:length(w), function(x) {
          paste0(gsub(":$", "/", hpcFls2[c(1, w[1:(l-1)]+1)][x]), hpcFls2[c(2, w[1:(l-1)]+2)[x]:(w[x]-1)])
        })
        hpcFls2 <- unlist(hpcFls2)
        fls1 <- data.frame(File = list.files(fl, full.names = TRUE, recursive = TRUE))
        fls1$Size <- file.info(fls1$File)$size
        flsSizes[[fl]] <<- sum(fls1$Size, na.rm = TRUE)
        fls1$HPC <- gsub(topattern(dirname(fl)), BATCHDIR, fls1$File)
        fls1$HPCxst <- fls1$HPC %in% hpcFls2
        res <- !sum(!fls1$HPCxst)
        if (res) {
          fls1$HPCsz <- sapply(fls1$HPC, function(fl) { as.numeric(capture.output(ssh_exec_wait(sshsess, paste0("ls -lrt ", fl, " | nawk '{print $5}'")))[1]) })
          res <- !sum(fls1$Size != fls1$HPCsz)
        }
      }
      return(res)
    })
    fls <- fls[which(!tst)]
  }
  if (Mode == names(Modes)[2]) {
    fls <- FILES[which(file.exists(FILES))]
    hpcFls2 <- capture.output(ssh_exec_wait(sshsess, paste0("ls ", BATCHDIR)))
    tst <- sapply(fls, function(fl) { #fl <- fls[1]
      res <- basename(fl) %in% hpcFls2
      if (res) {
        HPCsz <- as.numeric(capture.output(ssh_exec_wait(sshsess, paste0("ls -lrt ", gsub(topattern(dirname(fl)), BATCHDIR, fl), " | nawk '{print $5}'")))[1])
        flsSizes[[fl]] <<- file.size(fl)
        res <- flsSizes[[fl]] == HPCsz
      }
      return(res)
    })
    fls <- fls[which(!tst)]
  }
  flsSizes <- flsSizes[fls]
  fls <- c(fls, LIBRARIES) # Those we can overwrite for now
  if (length(fls)) {
    for (fl in fls) { #fl <- fls[1]
      if (!fl %in% names(flsSizes)) { flsSizes[[fl]] <- file.size(fl) }
      if (Mode == names(Modes)[1]) {
        allFls <- list.files(fl, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
        lrgFls <- grep("\\.tdf_bin$|/linespectra$", allFls, value = TRUE)
        allFls <- allFls[which(!allFls %in% lrgFls)]
        for (aFl in allFls) { #aFl <- allFls[1]
          btchDr <- gsub(topattern(dirname(fl)), BATCHDIR, dirname(aFl))
          if (!dir.exists(btchDr)) { dir.create(btchDr, recursive = TRUE) }
          tst <- try(scp_upload(sshsess, aFl, btchDr, verbose = TRUE), silent = TRUE)
        }
        if (length(lrgFls)) {
          for (lrgfl in lrgFls) {
            btchDr <- gsub(topattern(dirname(fl)), BATCHDIR, dirname(lrgfl))
            cmd <- paste0("scp \"", lrgfl, "\"  ", usrNm, "@", hstNm, ":\"", btchDr, "/", basename(lrgfl), "\"")
            #cat(cmd)
            tst <- try(system(cmd), silent = TRUE)
          }
        }
      }
    }
    if (Mode == names(Modes)[2]) {
      # Upload in one batch
      scp_upload(sshsess, fls, BATCHDIR, verbose = TRUE)
    }
    #scp_upload(sshsess, fls, BATCHDIR, verbose = TRUE) # Could be done like this, but for too large data amounts the function seems to sometimes fail
  }
}

# Upload human readable files
# For these we should make sure they have Unix line returns
fls2 <- FASTAS
if (length(fls2)) {
  for (fl in fls2) { #fl <- fls2[1]
    if (!fl %in% names(flsSizes)) { flsSizes[[fl]] <- file.size(fl) }
    flnm <- basename(fl)
    FILE <- file(flnm, "wb")
    dat <- readLines(fl)
    write(dat, FILE)
    close.connection(FILE)
    scp_upload(sshsess, flnm, BATCHDIR, verbose = TRUE)
    unlink(flnm)
  }
}

# Create and upload SLURM script:
# NB: In the SLURM script:
# - lines without a starting # are system commands
# - lines with one starting # are SLURM commands
# - lines starting with 2 or more starting # are ignored by SLURM (commented), but kept for future reference.
hpcFls2 <- capture.output(ssh_exec_wait(sshsess, paste0("ls ", BATCHDIR)))
ITER <- 0
SLURMnm <- NA
while ((is.na(SLURMnm))||(SLURMnm %in% hpcFls2)) {
  ITER <- ITER+1
  SLURMnm <- paste0("SLURM_script", as.character(ITER), ".sh")
}
flsSizes <- sum(unlist(flsSizes), na.rm = TRUE)
SLURM <- c("#!/bin/bash",
           paste0("#SBATCH --job-name=\"", JOBNAME, "\""),
           "#SBATCH --account=lsfgrp",
           "#SBATCH --time=240:00:00",
           #paste0("#SBATCH --output=\"", BATCHDIR, "/log_", ITER,"_out.txt\""), # No need to write what is essentially a duplicate of DiaNN's log at log level = 5
           "###Define the amount of RAM used by your job in MegaBytes",
           #"#SBATCH --mem=0", #To allow access to all available memory on each node
           paste0("#SBATCH --mem=", ceiling(flsSizes*4/1000000), "G"), #... but we actually need this, otherwise it may run on nodes with too little memory!
           #"#SBATCH --exclusive", # I do not think that we need this
           "##Send emails when a job starts, it is finished or it exits",
           paste0("#SBATCH --mail-user=", EMAIL),
           "#SBATCH --mail-type=ALL",
           paste0("#SBATCH -c ", NTHREADS),
           #"source /etc/profile.d/*.sh", 
           "",
           "unset SLURMEXPORT_ENV",
           "",
           "module purge",
           "module load diann",
           #"export LC_ALL=C", # Only if locale is not English
           #"unset LANGUAGE", # Only if locale is not English
           "ulimit -S -n 131072",
           "",
           "ahost=$(hostname)",
           "",
           "echo \"Running on hostname ... ${ahost}\"",
           "echo \" \"",
           "",
           "echo \" \"",
           "echo \"PATH=$PATH\"",
           "echo \" \"",
           "",
           paste0("running_dir=\"", BATCHDIR, "\""),
           "",
           "export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK",
           #"wait $!",
           DIANNcall,
           "cat /proc/sys/vm/max_map_count",
           "",
           "echo \"End of the processing ...\"",
           "date")
FILE <- file(SLURMnm, "wb")
write(SLURM, FILE)
close.connection(FILE)
scp_upload(sshsess, SLURMnm, BATCHDIR, verbose = TRUE)
unlink(SLURMnm)

# Run
batchcmd <- paste0(shell, "sbatch ", LINLOC, "/", BATCH, "/", SLURMnm)
ssh_exec_wait(sshsess, batchcmd)
ssh_exec_wait(sshsess, queuecmd)

#rm(list = ls()[which(!ls() %in% c("sshsess", "shell", "queuecmd", "sshDisc", "BATCH", "topattern"))])
#sshDisc()
