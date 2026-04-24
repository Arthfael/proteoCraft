# This script will
# - interactively create a SLURM script for running MaxQuant jobs on the HPC
# - edit the mqpar.xml selected by the user
# - transfer all data to the HPC and install the correct MaxQuant version there
# - start the job (this may send a "Houston, we have a problem" error,
#   in which case you can run:
#      ssh::ssh_exec_wait(sshsess, batchcmd)
#   to start the job)

packs <- c("svDialogs", "ssh")
for (pack in packs) {
  if (!require(pack, character.only = TRUE)) { install.packages(pack) }
  require(pack, character.only = TRUE)
}

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
require(tools)
SSH_on <- FALSE
if (!SSH_on) {
  msg <- "SSH session required\nEnter host name!"
  sshost <- svDialogs::dlg_input(msg)$res
  sshsess <- c()
  kount <- 1
  while ((kount < 5)&&(!SSH_on)) {
    sshsess <- try(ssh::ssh_connect(sshost), silent = TRUE)
    SSH_on <- "ssh_session" %in% class(sshsess)
    kount <- kount + 1
  }
  #print(sshsess)
}
if (!SSH_on) { stop() }
user <- gsub("^connected: |@.+", "", as.character(capture.output(sshsess))[2])
print(sshsess)
ssh::ssh_exec_wait(sshsess, "source /etc/profile.d/modules.sh")

shell <- "source /etc/profile.d/10_lsb-release-env-vars.sh ; source /etc/profile.d/99_load_slurm.sh ; source /etc/profile.d/modules.sh ; "
print("Checking ongoing jobs:")
queuecmd <- paste0(shell, "squeue --user=", user, " -l -o \"%.7i %.8u %.8j %.9a %.9P %.9N %.7T %.10M %.9L %.4C %.7m %.6D %R\"")
#queuecmd <- paste0(shell, "squeue --user=", user, " -l")
#cat(queuecmd)
ssh::ssh_exec_wait(sshsess, queuecmd)
# Function to disconnect the session
sshDisc <- function(session = sshsess) {
  try(ssh::ssh_disconnect(sshsess), silent = TRUE)
  if (exists("sshsess", .GlobalEnv)) { rm(sshsess, .GlobalEnv) }
}

MYEMAIL <- "e.g. myemail@domain.com"
while (!grepl("^[^ ]+@[^ ]+\\.[^ ]+$", MYEMAIL)) {
  MYEMAIL <- svDialogs::dlg_input("Enter your email address", MYEMAIL)$res
}

rm(list = ls()[which(!ls() %in% c("sshsess", "shell", "queuecmd", "sshDisc", "BATCH", "topattern"))])
# Create basic SLURM script, which we will edit as we go:
# NB: In the SLURM script:
# - lines without a starting # are system commands
# - lines with one starting # are SLURM commands
# - lines starting with 2 or more starting # are ignored by SLURM (commented), but kept for future reference.
SLURM <- c("#!/bin/bash",
           "#SBATCH --job-name=\"JOBNAME\"",
           "#SBATCH --nodes=1",
           "#SBATCH --account=lsfgrp",
           "#SBATCH --time=240:00:00",
           "#SBATCH --output=\"LINLOC/BATCH/log_ITER.out\"",
           "###Define the amount of RAM used by your job in MegaBytes",
           "#SBATCH --mem=0", #To access all available memory
           "#SBATCH --partition=MYPARTITION",
           "#SBATCH --exclusive",
           "##Send emails when a job starts, it is finished or it exits",
           paste0("#SBATCH --mail-user=", MYEMAIL),
           "#SBATCH --mail-type=ALL",
           "source /etc/profile.d/*.sh", 
           "",
           "unset SLURMEXPORT_ENV",
           "env > BATCHDIR/e1.txt",
           "",
           "module purge",
           "module load mono",
           "module load --ignore_cache dotnet/2.1.817",
           "",
           "ahost=$(hostname)",
           "",
           "echo \"Running on hostname ... ${ahost}\"",
           "echo \" \"",
           "",
           "echo \" \"",
           "echo \"PATH= $PATH\"",
           "echo \" \"",
           "",
           "running_dir=\"MQ_LOC\"",
           "input_file_dir=\"BATCHDIR\"",
           "MQPAR=\"${input_file_dir}/mqpar_tmp.xml\"",
           "",
           "NumThreads=$(($SLURM_CPUS_ON_NODE - 1))",
           "sed 's#<numThreads>.*</numThreads>#<numThreads>'$NumThreads'</numThreads>#' $MQPAR > ${input_file_dir}/mqpar_ITER.xml",
           #"wait $!",
           "srun --cpu_bind=verbose /usr/bin/time -v $(which mono) AOT_YESNO ${running_dir}/MaxQuantCmd.exe ${input_file_dir}/mqpar_ITER.xml",
           "",
           "cat /proc/sys/vm/max_map_count",
           "",
           "echo \"End of the processing ...\"",
           "date")
#
# Choose the name of the batch folder
if (exists("BATCH")) { def <- BATCH } else { def <- "My_batch" }
BATCH <- gsub(" ", "_", svDialogs::dlg_input("Enter the name of your batch folder", def)$res)
#
# Job name (currently the same as that of the batch folder)
JobName <- BATCH #gsub(" ", "_", svDialogs::dlg_input("Enter the name of your job", BATCH)$res)
SLURM[which(SLURM == "#SBATCH --job-name=\"JOBNAME\"")] <- paste0("#SBATCH --job-name=\"", JobName, "\"")
#
# Get the location where the script will be saved:
LINLOC <- capture.output(ssh::ssh_exec_wait(sshsess, "pwd"))[1]
LINLOC <- gsub("/+$", "", gsub("/+", "/", paste0("/", LINLOC)))
BATCHDIR <- paste0(LINLOC, "/", BATCH)
#ssh::ssh_exec_wait(sshsess, paste0("env > ", BATCHDIR, "/e0.txt"))
SLURM <- gsub("BATCHDIR", BATCHDIR, SLURM)
ssh::ssh_exec_wait(sshsess, paste0("cd \"", LINLOC, "\""))
out <- capture.output(ssh::ssh_exec_wait(sshsess, "ls -l | grep '^d'"))
if (length(grep(paste0(" ", BATCH, "$"), out)) == 0) { ssh::ssh_exec_wait(sshsess, paste0("mkdir -p ", BATCH)) }
#
# Select local, pre-created MaxQuant parameters file:
# This should have been created in MaxQuant and contain all the required information, e.g. files, search database, mods, etc...
# Make sure it is from the same MaxQuant version you will be using.
# Currently this does NOT check whether the raw/mzXML files actually exist!!!
msg <- "Select the MaxQuant parameters file to edit"
#filt <- matrix(c("MaxQuant parameters xml file", "*.xml"), ncol = 2)
#xmlnm <- normalizePath(choose.files("D:\\groups_temp\\mqpar.xml", msg, multi = FALSE, filt, 1), winslash = "/")
xmlnm <- rstudioapi::selectFile(msg,
                                path = "D:\\groups_temp\\mqpar.xml",
                                filter = "XML file (*.xml)")
xml <- readLines(xmlnm)
xmlnm <- basename(xmlnm)
xmlmqvers <- unlist(strsplit(grep("maxQuantVersion", xml, value = TRUE), "> ?| ?<"))[3]
#

MQvers <- grep("^MaxQuant_.+\\.zip$", list.files("...MaxQuant_vers_Archive"), value = TRUE)
MQvers <- data.frame(File = MQvers, Version = gsub("^MaxQuant_|\\.zip$", "", MQvers))
MQvers[,paste0("V", 1:4)] <- as.data.frame(t(sapply(strsplit(MQvers$Version, "\\."), function(x) {
  x <- as.numeric(unlist(x))
  l <- length(x)
  if (l < 4) { x <- c(x, rep(0, 4-l)) }
  return(x)
})))
MQvers <- MQvers[order(MQvers$V1, MQvers$V2, MQvers$V3, MQvers$V4, decreasing = TRUE), , drop = FALSE]
m <- match(xmlmqvers, MQvers$Version)
if (is.na(m)) {
  warning("The MaxQuant version used to create this parameter file does not appear to be available in our MaxQuant versions archive folder.")
  m <- 1
}
wmq <- as.numeric(svDialogs::dlg_input(paste(c("Which version of MaxQuant should we use?",
                                               paste0(" - ", 1:nrow(MQvers), ": ", MQvers$Version)), collapse = "\n"), m)$res)
mqvers <- MQvers$Version[wmq]
#
# Copy the selected MaxQuant zip to the HPC and extract
MQloc <- paste0(LINLOC, "/", BATCH, "/MaxQuant_", mqvers)
MQ_LOC <- paste0(MQloc, "/MaxQuant/bin")
out <- ssh::ssh_exec_wait(sshsess, paste0("[ -f \"", MQ_LOC, "/MaxQuantCmd.exe\" ]"))
if (out != 0) {
  # Aaaahhh... consistency between releases...
  MQ_LOC <- paste0(MQloc, "/MaxQuant_", mqvers, "/bin")
  out <- ssh::ssh_exec_wait(sshsess, paste0("[ -f \"", MQ_LOC, "/MaxQuantCmd.exe\" ]"))
}
if (out != 0) {
  MQ_LOC <- paste0(MQloc, "/MaxQuant/bin") # It wasn't that after all!
  ssh::scp_upload(sshsess, paste0("...MaxQuant_vers_Archive/MaxQuant_", mqvers, ".zip"), paste0(LINLOC, "/", BATCH), verbose = TRUE)
  ssh::ssh_exec_wait(sshsess, paste0("unzip \"", MQloc, ".zip\" -d ", MQloc))
  ssh::ssh_exec_wait(sshsess, paste0("unlink \"", MQloc, ".zip\""))
  #
}
out <- ssh::ssh_exec_wait(sshsess, paste0("[ -f \"", MQ_LOC, "/MaxQuantCmd.exe\" ]"))
if (out != 0) {
  # Same as above (both are needed)
  MQ_LOC <- paste0(MQloc, "/MaxQuant_", mqvers, "/bin")
  out <- ssh::ssh_exec_wait(sshsess, paste0("[ -f \"", MQ_LOC, "/MaxQuantCmd.exe\" ]"))
}
runnow <- FALSE
dryrun <- FALSE
if (out != 0) { stop("MaxQuant extraction failed!") } else {
  # Also ask for local modifications file
  locmod <- c(FALSE, TRUE)[match(svDialogs::dlg_message("Do you want to transfer a custom modifications file?",
                                                        "yesno")$res, c("no", "yes"))]
  if (locmod) {
    msg <- "Select the local modifications file to upload"
    #filt <- matrix(c("MaxQuant modifications xml file", "*.xml"), ncol = 2)
    #locmod <- normalizePath(choose.files(paste0("C:/MaxQuant/MaxQuant_", MQvers$Version[1], "/bin/conf/modifications.local.xml"), msg, multi = FALSE, filt, 1), winslash = "/")
    locmod <- rstudioapi::selectFile(msg,
                                     path = paste0("C:/MaxQuant/MaxQuant_", MQvers$Version[1], "/bin/conf/modifications.local.xml"),
                                     filter = "XML file (*.xml)")
    ssh::scp_upload(sshsess, locmod, paste0(MQ_LOC, "/conf/"), verbose = TRUE)
  }
  SLURM <- gsub("MQ_LOC", MQ_LOC, SLURM)
  continue <- TRUE
  if (mqvers != xmlmqvers) {
    continue <- svDialogs::dlg_message(paste0("This mqpar.xml file was created with a different version of MaxQuant (",
                                              xmlmqvers,
                                              ") than the one you chose to use (", mqvers, ")!\nThis can cause issues or crashed. Do you want to continue?"),
                                       "yesno")$res
    continue <- c(TRUE, FALSE)[match(continue, c("yes", "no"))]
  } 
  if (continue) {
    #
    # Edit raw (or mzXML) file paths and transfer them to the cluster
    # !!!Note: if files are in different subfolders, on the cluster they will all be in a single folder!!!
    r1 <- grep("<filePaths>", xml, ignore.case = TRUE)+1
    r2 <- grep("</filePaths>", xml, ignore.case = TRUE)-1
    if ((length(r1) > 1)||((length(r2) > 1))) { stop("Your MaxQuant parameters file appears to be corrupted, there can only be one \"filepaths\" field!")}
    if ((length(r1) == 0)||((length(r2) == 0))) { stop("Your MaxQuant parameters file  appears to be corrupted, check your \"filepaths\" field!")}
    nraw <- r2-r1+1
    raw <- xml[r1:r2]
    Raw <- as.data.frame(t(sapply(strsplit(raw, "<string>|</string>"), unlist)))
    sp <- unique(Raw[, 1])[1]
    Raw <- data.frame(Full = Raw[,2], Path = dirname(Raw[,2]), Name = basename(Raw[,2]))
    print(Raw$Name)
    ext <- unique(gsub(".+\\.", "", Raw$Name))
    if (length(ext) > 1) { stop("Only one file type is allowed!") }
    if (!tolower(ext) %in% c("raw", "mzxml")) { warning("Only Thermo Raw files (.raw) or mzXML files have been tested so far, proceed at your own risk!") }
    Raw$Size <- file.info(Raw$Full)$size
    out <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", BATCHDIR, " -l")))
    out <- grep(paste0("\\.", ext, "$"), out, value = TRUE, ignore.case = TRUE)
    if (length(out)) {
      out <- as.data.frame(t(sapply(strsplit(out, " "), function(x) {
        x <- unlist(x)
        tt <- which(nchar(x) > 0)
        x <- x[c(tt[1:7], tt[8]:length(x))]
        return(c(x[1:8], paste(x[9:length(x)], collapse = " ")))
      })))
      colnames(out) <- c("Permissions", "Dir.&links", "User", "Group", "bytes", "Mod_mm", "Mod_dd", "Mod_tt", "Name")
      out <- out[which(out$Name %in% Raw$Name),]
      out$bytes <- as.numeric(out$bytes)
      Raw$Copied_bytes <- out$bytes[match(Raw$Name, out$Name)]
      w <- which((!Raw$Name %in% out$Name)|(Raw$Copied_bytes != Raw$Size))
    } else { w <- 1:nrow(Raw) }
    while (length(w) > 0) {
      ssh::scp_upload(sshsess, Raw$Full[w], BATCHDIR, verbose = TRUE)
      out <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", BATCHDIR, " -l")))
      out <- grep(paste0("\\.", ext, "$"), out, value = TRUE, ignore.case = TRUE)
      out <- as.data.frame(t(sapply(strsplit(out, " "), function(x) {
        x <- unlist(x)
        tt <- which(nchar(x) > 0)
        x <- x[c(tt[1:7], tt[8]:length(x))]
        return(c(x[1:8], paste(x[9:length(x)], collapse = " ")))
      })))
      colnames(out) <- c("Permissions", "Dir.&links", "User", "Group", "bytes", "Mod_mm", "Mod_dd", "Mod_tt", "Name")
      out$bytes <- as.numeric(out$bytes)
      out <- out[which(out$Name %in% Raw$Name),]
      Raw$Copied_bytes <- out$bytes[match(Raw$Name, out$Name)]
      w <- which((!Raw$Name %in% out$Name)|(Raw$Copied_bytes != Raw$Size))
    }
    xml[r1:r2] <- paste0(sp, "<string>", LINLOC, "/", BATCH, "/", Raw$Name, "</string>")
    rm(raw, r1, r2, sp)
    #
    # Edit database file paths:
    f <- grep("<fastaFilePath>", xml, ignore.case = TRUE)
    if (length(f) == 0) { stop("Your MaxQuant parameters file looks corrupted, at least one \"fastaFilePath\" field should be present!")}
    fasta <- xml[f]
    Fasta <- as.data.frame(t(sapply(strsplit(fasta, "<fastaFilePath>|</fastaFilePath>"), unlist)))
    sp <- unique(Fasta[,1])[1]
    Fasta <- data.frame(Full = Fasta[, 2], Path = dirname(Fasta[,2]), Name = basename(Fasta[,2]))
    print(Fasta$Name)
    Fasta$Size <- file.info(Fasta$Full)$size
    out <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", BATCHDIR, " -l")))
    out <- grep("\\.fasta$", out, value = TRUE, ignore.case = TRUE)
    if (length(out)) {
      out <- as.data.frame(t(sapply(strsplit(out, " "), function(x) {
        x <- unlist(x)
        tt <- which(nchar(x) > 0)
        x <- x[c(tt[1:7], tt[8]:length(x))]
        return(c(x[1:8], paste(x[9:length(x)], collapse = " ")))
      })))
      colnames(out) <- c("Permissions", "Dir.&links", "User", "Group", "bytes", "Mod_mm", "Mod_dd", "Mod_tt", "Name")
      out$bytes <- as.numeric(out$bytes)
      out <- out[which(out$Name %in% Fasta$Name),]
      Fasta$Copied_bytes <- out$bytes[match(Fasta$Name, out$Name)]
      w <- which((!Fasta$Name %in% out$Name)|(Fasta$Copied_bytes != Fasta$Size))
    } else { w <- 1:nrow(Fasta) }
    while (length(w)) {
      ssh::scp_upload(sshsess, Fasta$Full[w], BATCHDIR, verbose = TRUE)
      out <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", BATCHDIR, " -l")))
      out <- grep("\\.fasta$", out, value = TRUE, ignore.case = TRUE)
      out <- as.data.frame(t(sapply(strsplit(out, " "), function(x) {
        x <- unlist(x)
        tt <- which(nchar(x) > 0)
        x <- x[c(tt[1:7], tt[8]:length(x))]
        return(c(x[1:8], paste(x[9:length(x)], collapse = " ")))
      })))
      colnames(out) <- c("Permissions", "Dir.&links", "User", "Group", "bytes", "Mod_mm", "Mod_dd", "Mod_tt", "Name")
      out$bytes <- as.numeric(out$bytes)
      out <- out[which(out$Name %in% Fasta$Name),]
      Fasta$Copied_bytes <- out$bytes[match(Fasta$Name, out$Name)]
      w <- which((!Fasta$Name %in% out$Name)|(Fasta$Copied_bytes != Fasta$Size))
    }
    xml[f] <- paste0(sp, "<fastaFilePath>", LINLOC, "/", BATCH, "/", Fasta$Name, "</fastaFilePath>")
    rm(fasta, f, sp)
    #
    # For MaxQuant >= 2.0.0.0, DIA search is available
    # If so, we need to transfer library files
    type <- gsub(" *</?lcmsRunType> *", "", grep("^ *<lcmsRunType>", xml, value = TRUE))
    if (grepl("MaxDIA", type)) {
      libType <- as.numeric(gsub(" *</?diaLibraryType> *", "", grep("^ *<diaLibraryType>", xml, value = TRUE)))
      if (libType != 0) {
        stop("This script currently only supports MaxQuant txt libraries... for now. You know what to do, do it with style!")
      } else {
        grPep <- grep(" *</?diaPeptidePaths> *", xml)
        grEv <- grep(" *</?diaEvidencePaths> *", xml)
        grMSMS <- grep(" *</?diaMsmsPaths> *", xml)
        PepLib <- normalizePath(gsub(" *</?string> *", "", xml[(grPep[1]+1):(grPep[2]-1)]), winslash = "/")
        EvLib <- normalizePath(gsub(" *</?string> *", "", xml[(grEv[1]+1):(grEv[2]-1)]), winslash = "/")
        MSMSLib <- normalizePath(gsub(" *</?string> *", "", xml[(grMSMS[1]+1):(grMSMS[2]-1)]), winslash = "/")
        Lib <- c(PepLib, EvLib, MSMSLib)
        tst <- vapply(Lib, function(x) { file.info(x)$size }, 1)
        # Create destination folder if it doesn't exist
        capture.output(ssh::ssh_exec_wait(sshsess, paste0("mkdir -p ", BATCHDIR, "/DIA_library")))
        # Upload the files
        LibNms <- gsub(".+/", "", Lib)
        w <- 1:length(Lib)
        while (length(w)) {
          out <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", BATCHDIR, "/DIA_library", " -l")))
          out <- grep("\\.txt$", out, value = TRUE, ignore.case = TRUE)
          if (length(out)) {
            out <- as.data.frame(t(sapply(strsplit(out, " "), function(x) {
              x <- unlist(x)
              tt <- which(nchar(x) > 0)
              x <- x[c(tt[1:7], tt[8]:length(x))]
              return(c(x[1:8], paste(x[9:length(x)], collapse = " ")))
            })))
            colnames(out) <- c("Permissions", "Dir.&links", "User", "Group", "bytes", "Mod_mm", "Mod_dd", "Mod_tt", "Name")
            out$bytes <- as.numeric(out$bytes)
            w <- which(!LibNms %in% out$Name)
            w2 <- which(LibNms %in% out$Name)
            if (length(w2)) {
              m <- match(LibNms[w2], out$Name)
              w <- c(w, w2[which(tst[w2] != out$bytes[m])])
            }
          }
          if (length(w)) {
            print("Uploading DIA library files (MaxQuant txt format):")
            ssh::scp_upload(sshsess, Lib[w], paste0(BATCHDIR, "/DIA_library"))
          }
        }
        xml[(grPep[1]+1):(grPep[2]-1)] <- paste0("            <string>", BATCHDIR, "/DIA_library/", LibNms[1], "</string>")
        xml[(grEv[1]+1):(grEv[2]-1)] <- paste0("            <string>", BATCHDIR, "/DIA_library/", LibNms[2], "</string>")
        xml[(grMSMS[1]+1):(grMSMS[2]-1)] <- paste0("            <string>", BATCHDIR, "/DIA_library/", LibNms[3], "</string>")
      }
    }
    #
    # Only valid for MaxQuant up to 1.6.12.0 (not 1.6.14, I don't know about 1.6.13) 
    # Edit combined folder location:
    # If we have enough space, we would want to have it running in the local node's home folder for faster I/O.
    v <- as.numeric(unlist(strsplit(mqvers, "\\.")))
    if ((v[1] < 1)||((v[1] == 1)&&(v[2] < 6))||((v[1] == 1)&&(v[2] == 6)&&(v[3] <= 13))) {
      localcomb <- svDialogs::dlg_message("Do you want to use a local temporary folder?", "yesno")$res
      localcomb <- c("", paste0("/localhome/", user))[which(c("no", "yes") == localcomb)]
      cf1 <- grep("<fixedCombinedFolder>", xml, ignore.case = TRUE)
      cf2 <- grep("</fixedCombinedFolder>", xml, ignore.case = TRUE)
      if ((length(cf1) != 1)||(length(cf2) != 1)) {
        stop("Your MaxQuant parameters file looks corrupted, no more than one \"fixedCombinedFolder\" field should be present!")
      } else { xml[cf1:cf2] <- paste0("   <fixedCombinedFolder>", localcomb, "</fixedCombinedFolder>") }
      rm(cf1, cf2)
    }
    #
    # Use .Net core?
    g <- grep("useDotNetCore", xml)
    if (length(g)) {
      useDotNet <- c(FALSE, TRUE)[match(svDialogs::dlg_message("Do you want to use .NET Core?", "yesno")$res,
                                        c("no", "yes"))]
      xml[g] <- paste0("   <useDotNetCore>", c("False", "True")[useDotNet+1], "</useDotNetCore>")
      w1 <- which(SLURM == "module load --ignore_cache dotnet/2.1.817")
      w2 <- c(1:length(SLURM))[which(!c(1:length(SLURM)) %in% w1)]
      if (useDotNet) {
        v <- as.numeric(unlist(strsplit(mqvers, "\\.")))
        if  (((v[1] == 2)&&(v[2] == 0)&&(v[3] >= 3))||((v[1] == 2)&&(v[2] > 0))||(v[1] > 2)) {
          SLURM[w1] <- "module load --ignore_cache dotnet/3.1.414"
          # NB: do not use 3.1.406 (or earlier?) or MQ will fail!!!
        }
      } else { SLURM <- SLURM[w2] }
    }
    #
    # Choose partition and node:
    # The Cluster contains different partitions, each containing one or more nodes.
    # Nodes from a same partition may have different numbers of vCPUs.
    # Because we need to enter a fixed number of threads in the MaxQuant parameters file, we need to choose a specific
    # node with known number of vCPUs.
    # We apparently still need to specify the partition too!
    #partcmd <- paste0(shell, "sinfo -l --user=", user)
    partcmd <- paste0(shell, "sinfo -l")
    partitions <- capture.output(ssh::ssh_exec_wait(sshsess, partcmd))
    partitions <- partitions[2:(length(partitions)-1)]
    tmp <- partitions[1]
    partitions <- partitions[2:length(partitions)]
    partitions <- as.data.frame(t(sapply(strsplit(partitions, " +"), unlist)))
    colnames(partitions) <- unlist(strsplit(tmp, " +"))
    partitions <- partitions[c(which(partitions$STATE %in% c("idle", "idle*")),
                               which(partitions$STATE %in% c("mix", "idle*")),
                               which(partitions$STATE %in% c("alloc", "alloc*")),
                               which(partitions$STATE %in% c("drain", "drain*")),
                               which(partitions$STATE %in% c("down", "down*")),
                               which(!partitions$STATE %in% c("idle", "idle*", "mix", "idle*", "alloc", "alloc*", "drain", "drain*", "down", "down*"))),]
    partitions$PARTITION <- gsub("\\*", "", partitions$PARTITION)
    partitionsmssg <- unique(sort(partitions$PARTITION))
    PARTITION <- as.numeric(svDialogs::dlg_input(paste(c("Which of the cluster's partitions do you want to use?",
                                                         paste0(" - ", partitionsmssg, ": ", 1:length(partitionsmssg))), collapse = "\n"),
                                                 match("defaultp", partitionsmssg))$res)
    PARTITION <- partitionsmssg[PARTITION]
    SLURM[which(SLURM == "#SBATCH --partition=MYPARTITION")] <- paste0("#SBATCH --partition=", PARTITION)
    nodecmd <- paste0(shell, "sinfo -N --partition=", PARTITION,
                      " --format=\"%.5a %.10l %.6D %.6t %N %C %m %e\"")
    nodes <- capture.output(ssh::ssh_exec_wait(sshsess, nodecmd))
    nodes <- nodes[1:(length(nodes)-1)]
    nodes <- gsub("^ +", "", nodes)
    tmp <- nodes[1]
    nodes <- nodes[2:length(nodes)]
    nodes <- magrittr::set_colnames(as.data.frame(t(sapply(strsplit(nodes, " +"), unlist))),
                                    unlist(strsplit(tmp, " +")))
    if (PARTITION == "defaultp") {
      nodes <- nodes[grep("^(bea|bigterra)[0-9]+$", nodes$NODELIST, invert = TRUE),]
    }
    nodes$Total_CPUs <- vapply(strsplit(nodes$`CPUS(A/I/O/T)`, "/"), function(x) { rev(unlist(x))[1] }, 1)
    nodes <- nodes[order(as.numeric(nodes$Total_CPUs), decreasing = TRUE),]
    nodes <- nodes[c(which(nodes$STATE %in% c("idle", "idle*")),
                     which(nodes$STATE %in% c("mix", "idle*")),
                     which(nodes$STATE %in% c("alloc", "alloc*")),
                     which(nodes$STATE %in% c("drain", "drain*")),
                     which(nodes$STATE %in% c("down", "down*")),
                     which(!nodes$STATE %in% c("idle", "idle*", "mix", "idle*", "alloc", "alloc*", "drain", "drain*", "down", "down*"))),]
    nodesmess <- apply(nodes[,c("NODELIST", "STATE", "Total_CPUs")], 1, function(x) {
      paste0(x[[1]], " (", x[[2]], ", ", x[[3]], " CPUs)")
    })
    nodes$MEMORY <- as.numeric(nodes$MEMORY)
    nodes$FREE_MEM <- suppressWarnings(as.numeric(nodes$FREE_MEM))
    #nodes <- nodes[which(!is.na(nodes$FREE_MEM)),]
    # Also specifically ask for bjoern's nodes for the purpose of the dry run (same code as above)
    bjoernnodecmd <- gsub(PARTITION, "defaultp", nodecmd)
    bjoern <- capture.output(ssh::ssh_exec_wait(sshsess, bjoernnodecmd))
    bjoern <- gsub("^ +", "", bjoern)
    tmp <- bjoern[1]
    bjoern <- bjoern[2:max(c(2, (length(bjoern)-1)))]
    bjoern <- magrittr::set_colnames(as.data.frame(t(sapply(strsplit(bjoern, " +"), unlist))), unlist(strsplit(tmp, " +")))
    bjoern$Total_CPUs <- vapply(strsplit(bjoern$`CPUS(A/I/O/T)`, "/"), function(x) { rev(unlist(x))[1] }, 1)
    bjoern <- bjoern[order(as.numeric(bjoern$Total_CPUs), decreasing = TRUE),]
    bjoern <- bjoern[c(which(bjoern$STATE %in% c("idle", "idle*")),
                       which(bjoern$STATE %in% c("mix", "idle*")),
                       which(bjoern$STATE %in% c("alloc", "alloc*")),
                       which(bjoern$STATE %in% c("drain", "drain*")),
                       which(bjoern$STATE %in% c("down", "down*")),
                       which(!bjoern$STATE %in% c("idle", "idle*", "mix", "idle*", "alloc", "alloc*", "drain", "drain*", "down", "down*"))),]
    bjoern <- bjoern[grep("^bjoern", bjoern$NODELIST),]
    bjoern$MEMORY <- as.numeric(bjoern$MEMORY)
    bjoern$FREE_MEM <- suppressWarnings(as.numeric(bjoern$FREE_MEM))
    bjoern <- bjoern[which(!is.na(bjoern$FREE_MEM)),]
    if (length(unique(gsub("/.+", "", nodes$`CPUS(A/I/O/T)`))) > 1) {
      # Choose node or constraint:
      msg <- "Not all nodes on this partition are equivalent.\nDo you want to:\n - 1: use any node,\n - 2: restrict the job to a constraint,\n - 3: or to a specific node?"
      constraints <- as.numeric(svDialogs::dlg_input(msg, 1)$res)
      while (!constraints %in% 1:3) { constraints <- as.numeric(svDialogs::dlg_input(msg, 1)$res) }
      if (constraints == 2) {
        tmp <- gsub("[0-9]+", "", nodes$NODELIST)
        msg <- aggregate(tmp, list(tmp), length)
        msg$Idles <- aggregate(nodes$STATE, list(tmp), function(x) { length(x[which(x == "idle")])})$x
        msg <- msg[order(msg$x, decreasing = TRUE), ]
        msg <- msg[order(msg$Idles, decreasing = TRUE), ]
        msg <- msg[which(msg$x > 1),]
        tmp <- msg$Group.1
        msg <- apply(msg[, c("Group.1", "Idles")], 1, function(x) { paste0(" - ", x[[1]], " (", x[[2]], " idle nodes)") })
        msg <- paste(c("Which constraint do you want to use? Suggested names (NB: detected, some or all may be invalid)", msg), collapse = "\n")
        CONSTRAINT <- svDialogs::dlg_input(msg, tmp[1])$res
        w <- grep("^#SBATCH --constraint=", SLURM)
        if (length(w)) {
          SLURM[w] <- paste0("#SBATCH --constraint=", CONSTRAINT)
        } else {
          w <- grep("^#SBATCH --partition=", SLURM)
          SLURM <- c(SLURM[1:w], paste0("#SBATCH --constraint=", CONSTRAINT), SLURM[(w+1):length(SLURM)])
        }
        #RAMavail <- min(nodes$MEMORY[grep(topattern(CONSTRAINT), nodes$NODELIST)]) - 1000
        #if (is.na(RAMavail)) { RAMavail <- min(nodes$MEMORY) - 1000 }
        #w <- grep("^#SBATCH --constraint=", SLURM)
        #SLURM <- c(SLURM[1:w], paste0("#SBATCH --mem=", RAMavail), SLURM[(w+1):length(SLURM)])
      }
      if (constraints == 3) {
        NODENAME <- as.numeric(svDialogs::dlg_input(paste(c("Which of the nodes do you want to use?",
                                                            paste0(" - ", nodesmess, ": ", 1:length(nodesmess))),
                                                          collapse = "\n"),
                                                    1)$res)
        NODENAME <- nodes$NODELIST[NODENAME]
        w <- grep("^#SBATCH --nodelist=", SLURM)
        if (length(w)) {
          SLURM[w] <- paste0("#SBATCH --nodelist=", NODENAME)
        } else {
          w <- grep("^#SBATCH --partition=", SLURM)
          SLURM <- c(SLURM[1:w], paste0("#SBATCH --nodelist=", NODENAME), SLURM[(w+1):length(SLURM)])
        }
        #RAMavail <- nodes$MEMORY[match(NODENAME, nodes$NODELIST)] - 1000
        #w <- grep("^#SBATCH --nodelist=", SLURM)
        #SLURM <- c(SLURM[1:w], paste0("#SBATCH --mem=", RAMavail), SLURM[(w+1):length(SLURM)])
      }
    }
    #
    #SLURM[grep("^#SBATCH --mem=", SLURM)] <- paste0("#SBATCH --mem=", RAMavail, "M")
    #
    # Ask about aot compilation:
    aot <- tolower(svDialogs::dlg_input("Should we use ahead of time (aot) compilation? (y/n)", "n")$res)
    while (!aot %in% c("y", "n")) { aot <- tolower(svDialogs::dlg_input("I repeat, dumbass: Should. We. Use. Ahead. Of. Time. Compilation? (y/n)", "n")$res) }
    aot <- c(TRUE, FALSE)[match(aot, c("y", "n"))]
    # Detect local SLURM scripts to identify iteration level:
    ITER <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", LINLOC, "/", BATCH, "/*")))
    ITER <- grep(topattern(paste0(LINLOC, "/", BATCH, "/")), ITER, value = TRUE)
    ITER <- gsub(topattern(paste0(LINLOC, "/", BATCH, "/")), "", ITER)
    ITER <- grep("^SLURM_script_", ITER, value = TRUE)
    ITER <- ITER[which(!grepl("_dryrun\\.sh$|_partial\\.sh", ITER))]
    ITER <- length(ITER) + 1
    SLURM <- gsub("ITER", ITER, SLURM)
    SLURM[grep("^#SBATCH --output=", SLURM)] <- paste0("#SBATCH --output=\"", BATCHDIR, "/log_", ITER, ".out\"")
    # Create mono call:
    xmlnm2 <- "mqpar_tmp.xml"
    g <- grep(" AOT_YESNO ", SLURM)
    SLURM[g] <- gsub(" AOT_YESNO ", c(" ", " --aot=full ")[aot+1], SLURM[g])
    SLURMnm <- paste0("SLURM_script_", ITER, ".sh")
    # Also create dry-run script and partial processing template script
    w <- which(bjoern$STATE == "idle")
    dryrun <- length(w) > 0
    if (dryrun) {
      w <- w[length(w)]
      SLURMdr <- SLURM
      SLURMdr <- SLURMdr[which(!grepl("#SBATCH --mail-user=", SLURMdr))] # For the dryrun, we can dispense with emails...
      SLURMdr[grep("^#SBATCH --partition=", SLURMdr)]  <- "#SBATCH --partition=defaultp"
      SLURMdr[grep("^#SBATCH --constraint=", SLURMdr)]  <- "#SBATCH --constraint=bjoern"
      #SLURMdr[grep("^#SBATCH --nodelist=", SLURMdr)]  <- paste0("#SBATCH --nodelist=", bjoern$NODELIST[w])
      SLURMdr[grep("^#SBATCH --mem=", SLURMdr)] <- paste0("#SBATCH --mem=",  bjoern$MEMORY[w], "M")
      SLURMdr[grep("sed", SLURMdr)] <- gsub("\\.xml", "_dryrun.xml", SLURMdr[grep("sed", SLURMdr)])
      SLURMdr[grep("srun", SLURMdr)] <- paste0(gsub("\\.xml", "_dryrun.xml", SLURMdr[grep("srun", SLURMdr)]), " -n")
      SLURMdr[grep("^#SBATCH --output=", SLURMdr)] <- paste0("#SBATCH --output=\"", BATCHDIR, "/log_", ITER, "_dryrun.out\"")
      xmldr <- xml
      cpudr <- as.numeric(unlist(strsplit(bjoern$`CPUS(A/I/O/T)`[w], "/"))[1])
      xmldr[grep("Threads", xmldr)] <- paste0("   <numThreads>",  bjoern$Total_CPUs[w]-1, "</numThreads>")
      SLURMdrnm <- gsub("\\.sh$", "_dryrun.sh", SLURMnm)
    }
    SLURMpp <- SLURM
    SLURMpp[grep("srun", SLURMpp)]  <- paste0(SLURMpp[grep("srun", SLURMpp)], " --partial-processing=PARTIALPROCESSINGSTEP")
    # 
    #
    # This is it, editing done!
    #
    # Write the modified files to a folder in Q: then move them to scistore13 and run the SLURM script:
    setwd("...Archive/")
    if (!dir.exists(BATCH)) { dir.create(BATCH, recursive = TRUE) }
    write(xml, paste0(BATCH, "/", gsub(" ", "_", xmlnm2)))
    FILE <- file(paste0("./", BATCH, "/", SLURMnm), "wb")
    write(SLURM, FILE)
    close.connection(FILE)
    ssh::scp_upload(sshsess, paste0("...Archive/", BATCH, "/", c(gsub(" ", "_", xmlnm2), SLURMnm)), BATCHDIR, verbose = TRUE)
    if (dryrun) {
      FILEdr <- file(paste0("./", BATCH, "/", SLURMdrnm), "wb")
      write(SLURMdr, FILEdr)
      close.connection(FILEdr)
      write(xmldr, paste0(BATCH, "/", gsub("\\.xml", "_dryrun.xml", gsub(" ", "_", xmlnm2))))
      ssh::scp_upload(sshsess, paste0("...Archive/", BATCH, "/", c(gsub("\\.xml", "_dryrun.xml", gsub(" ", "_", xmlnm2)), SLURMdrnm)),
                      BATCHDIR, verbose = TRUE)
    }
    Sys.sleep(1)
    setwd("...Archive/") # Just in case
    unlink(paste0("...Archive/", BATCH), recursive = TRUE, force = TRUE)
    runnow <- c(FALSE, TRUE)[match(svDialogs::dlg_message("Do you want to run the analysis now?", "yesno")$res,
                                   c("no", "yes"))]
  }
}
if (runnow) {
  partialrun <- FALSE
  out <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", LINLOC, "/", BATCH, "/*")))
  out <- grep(topattern(paste0(LINLOC, "/", BATCH, "/")), out, value = TRUE)
  if (paste0(LINLOC, "/", BATCH, "/log_", ITER, ".out") %in% out) {
    ssh::ssh_exec_wait(sshsess, paste0("unlink \"", LINLOC, "/", BATCH, "/log_", ITER, ".out\""))
  }
  if (dryrun) {
    print("Dry run:")
    drbatchcmd <- paste0(shell, "sbatch ", LINLOC, "/", BATCH, "/", SLURMdrnm)
    ssh::ssh_exec_wait(sshsess, drbatchcmd)
    Sys.sleep(1)
    out <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", LINLOC, "/", BATCH, "/*")))
    out <- grep(topattern(paste0(LINLOC, "/", BATCH, "/")), out, value = TRUE)
    while (!paste0(LINLOC, "/", BATCH, "/log_", ITER, "_dryrun.out") %in% out) {
      Sys.sleep(1)
      out <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", LINLOC, "/", BATCH, "/*")))
      out <- grep(topattern(paste0(LINLOC, "/", BATCH, "/")), out, value = TRUE)
    }
    suppressWarnings(rm(tasks))
    sleeptm <- 1
    tasks <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("cat ", LINLOC, "/", BATCH, "/log_", ITER, "_dryrun.out")))
    Sys.sleep(sleeptm)
    while ((!length(tasks))&&(sleeptm < 60)) {
      sleeptm <- sleeptm + 1
      tasks <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("cat ", LINLOC, "/", BATCH, "/log_", ITER, "_dryrun.out")))
      Sys.sleep(sleeptm)
    }
    Sys.sleep(1)
    tasks <- grep("[0-9]+\t", tasks, value = TRUE)
    if (!length(tasks)) {
      tasks <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("cat ", LINLOC, "/", BATCH, "/log_", ITER, "_dryrun.out")))
      tasks <- grep("[0-9]+\t", tasks, value = TRUE)
    }
  }
}
# Breaking the script like this feels absurd, but somehow otherwise the tasks are just not captured properly!
if (runnow) {
  if (dryrun) {
    tasks <- grep("[0-9]+\t", tasks, value = TRUE)
    if (!length(tasks)) {
      tasks <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("cat ", LINLOC, "/", BATCH, "/log_", ITER, "_dryrun.out")))
      tasks <- grep("[0-9]+\t", tasks, value = TRUE)
    }
    if (length(tasks)) {
      tasks <- as.data.frame(t(sapply(strsplit(tasks, "\t"), unlist)))
      colnames(tasks) <- c("ID", "Task")
      tasks$ID <- as.numeric(tasks$ID)
      out <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls ", LINLOC, "/", BATCH, "/combined/proc/")))
      out <- grep(topattern(paste0(LINLOC, "/", BATCH, "/combined/proc/")), out, value = TRUE)
      if (paste0(LINLOC, "/", BATCH, "/combined/proc/#runningTimes.txt") %in% out) {
        proc <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("cat ", LINLOC, "/", BATCH, "/combined/proc/#runningTimes.txt")))
        proc <- as.data.frame(t(as.data.frame(sapply(strsplit(proc, "\t"), unlist))))
        rownames(proc) <- NULL
        colnames(proc) <- proc[1,]
        proc <- proc[2:(nrow(proc)-1),]
        proc$TaskID <- tasks$ID[match(proc$Job, tasks$Task)]
        partialrun <- c(FALSE, TRUE)[match(svDialogs::dlg_message("Do you want to perform partial processing?",
                                                                  "yesno")$res, c("no", "yes"))]
      } else { partialrun <- FALSE }
    } else { stop("Houston we have a problem!!!") }
  }
  print("Submitting actual job...")
  if (partialrun) {
    print(paste0(apply(tasks[which(tasks$ID %in% proc$TaskID),], 1, paste, collapse = ": "), collapse = "\n"))
    partialrun_step <- svDialogs::dlg_input(paste0("From which task (1 to ", max(proc$TaskID),
                                                   ") do you want to perform partial processing?"),
                                            max(proc$TaskID) - 1)$res
    SLURMpp <- gsub("PARTIALPROCESSINGSTEP", partialrun_step, SLURMpp)
    SLURMppnm <- gsub("\\.sh$", "_partial.sh", SLURMnm)
    FILEpp <- file(paste0("./", BATCH, "/", SLURMppnm), "wb")
    write(SLURMpp, FILEpp)
    close.connection(FILEpp)
    ssh::scp_upload(sshsess, paste0("...Archive/", BATCH, "/", SLURMppnm), BATCHDIR, verbose = TRUE)
    ssh::ssh_exec_wait(sshsess, paste0("sbatch ", LINLOC, "/", BATCH, "/", SLURMppnm))
  } else {
    dirs <- capture.output(ssh::ssh_exec_wait(sshsess, paste0("ls -d ", LINLOC, "/", BATCH, "/*/")))
    if (paste0(LINLOC, "/", BATCH, "/combined/") %in% dirs) { ssh::ssh_exec_wait(sshsess, paste0("rm -r \"", LINLOC, "/", BATCH, "/combined/\"")) }
    for (i in Raw$Name) {
      if (paste0(LINLOC, "/", BATCH, "/", i, "/") %in% dirs) {
        ssh::ssh_exec_wait(sshsess, paste0("rm -r \"", LINLOC, "/", BATCH, "/", i, "/\""))
        ssh::ssh_exec_wait(sshsess, paste0("unlink \"", LINLOC, "/", BATCH, "/", i, ".index\""))
      }
    }
    batchcmd <- paste0(shell, "sbatch ", LINLOC, "/", BATCH, "/", SLURMnm)
    ssh::ssh_exec_wait(sshsess, batchcmd)
  }
}
# If "Error: Houston we have a problem!!!" is displayed then run the two commented lines below to start the job:
batchcmd <- paste0(shell, "sbatch ", LINLOC, "/", BATCH, "/", SLURMnm)
#ssh::ssh_exec_wait(sshsess, batchcmd)
ssh::ssh_exec_wait(sshsess, queuecmd)

#rm(list = ls()[which(!ls() %in% c("sshsess", "shell", "queuecmd", "sshDisc", "BATCH", "topattern"))])
#sshDisc()
