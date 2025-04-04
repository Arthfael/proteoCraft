# A script to combine n FragPipe output folders into one,
# for the purpose of combining the results into a single data analysis using the scripts from this package.
# This can then be processed by FP_to_MQ().
#
# Useful if you searched different files from the same experiment with different parameters but want to combine the outputs.
#
# Obviously, this is not meant to be used for just any FragPipe searches. This will handle specific limited cases only where it makes sense.
# It will break in some unauthorized, un-recommended or unsupported cases (e.g. different fixed modifications).

require(rstudioapi)
require(data.table)
try(setDTthreads(threads = parallel::detectCores()-1), silent = TRUE)

dflt <- "C:"
if ((exists("dirs"))&&(length(dirs))) { dflt <- gsub("/[^/]+$", "", rev(dirs)[1]) }
goOn <- TRUE
dirs <- c()
kount <- 0
while (goOn) {
  dir <- rstudioapi::selectDirectory(paste0("Select a", c("", "nother")[(kount > 0)+1], " FragPipe results directory"), path = dflt)
  if (is.null(dir)) { goOn <- FALSE } else {
    kount <- kount + 1
    dirs <- unique(c(dirs, dir))
    dflt <- gsub("/[^/]+$", "", dir)
  }
}

# Select destination directory
dstDir <- NA
while ((is.na(dstDir))||(is.null(dstDir))||(dstDir%in% dirs)) {
  dstDir <- selectDirectory("Select destination directory", path = gsub("/[^/]+$", "/.*", dirs[1]))
}
setwd(dstDir)

FPs <- setNames(lapply(dirs, function(dir) { #dir <- dirs[1]
  print(dir)
  fls <- list.files(dir, full.names = TRUE, recursive = FALSE)
  FP_WorkflowFl <- grep("/fragpipe\\.workflow$", fls, value = TRUE)
  FP_ManifestFl <- grep("/fragpipe-files\\.fp-manifest$", fls, value = TRUE)
  tst <- ((!is.na(FP_WorkflowFl))&&(length(FP_WorkflowFl) == 1)&&(!is.na(FP_ManifestFl))&&(length(FP_ManifestFl) == 1))
  if (tst) {
    FP_Mnfst <- suppressWarnings(read.delim(FP_ManifestFl, header = FALSE))
    colnames(FP_Mnfst) <- c("Path", "Experiment", "Replicate", "Data type")
    FP_Mnfst$Path <- gsub("\\\\", "/", FP_Mnfst$Path)
    return(list(Outcome = TRUE,
                Workflow = readLines(FP_WorkflowFl),
                Manifest = FP_Mnfst))
  } else {
    return(list(Outcome = FALSE))
  }
}), dirs)
w <- which(sapply(FPs, function(x) { x$Outcome }))
FPs <- FPs[w]; dirs <- dirs[w]
#
# Process manifests
Manifests <- setNames(lapply(FPs, function(x) { x$Manifest }), dirs)
tst <- sapply(dirs, function(x) {
  fls <- list.files(x)
  sum(dir.exists(paste0(x, "/", Manifests[[x]]$Experiment))) == nrow(Manifests[[x]])
})
if (sum(!tst)) {
  warning("We have encountered a non-supported case (no value entered in the Experiment column of the manifest)!\nWe will assume that all samples are from one experiment.")
  nuManifest <- Manifests[[1]]
} else {
  nuManifest <- plyr::rbind.fill(Manifests)
  tst <- do.call(paste, c(nuManifest, sep = "_|_"))
  tst <- aggregate(1:length(tst), list(tst), min)
  nuManifest <- nuManifest[tst$x,]
}
nuManifest_fl <- paste0(dstDir, "/fragpipe-files.fp-manifest")
write.table(nuManifest, nuManifest_fl, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, na = "")
# Obviously this does not deal with the case of the same file being saved in different locations.
# This can be dealt with during data post-processing.
#
#
# Process Workflows
Workflows <- setNames(lapply(FPs, function(x) { x$Workflow }), dirs)
#l <- lapply(Workflows, length) # Nope, they are not all the same length
# Version
fpVers <- unique(sapply(Workflows, function(x) { grep("^# FragPipe \\([0-9]+(\\.[0-9]+)?\\) runtime properties", x, value = TRUE) }))
if (length(fpVers) > 1) { fpVers <- "# FragPipe (different versions) runtime properties" }
# Fasta(s)
fastas <- #gsub("\\\\", "", gsub("\\\\\\\\", "/", 
  unique(sapply(Workflows, function(x) { grep("^database.db-path=", x, value = TRUE) }))#))
tst <- unique(gsub(".*/[0-9]{4}-[0-9]{2}-[0-9]{2}-", "", gsub("\\\\", "", gsub("\\\\\\\\", "/", fastas))))
if (length(tst) > 1) {
  stop("What should I do? FragPipe cannot deal with multiple fastas, yet it looks like you used different ones.")
  # Suggested ways to deal with this if this occurs:
  #  - Locate those fastas, if necessary prompting users to locate them.
  #  - Then combine them, save them to destDir and enter new database location in FragPipe pseudo-workflow
}
tmp <- grep("database.db-path=/nfs/", fastas, invert = TRUE, value = TRUE) # Prefer local to remote located databases
if (length(tmp)) { fastas <- tmp }
fastas <- fastas[1]

# Fixed mods
fixMods <- unique(sapply(Workflows, function(x) { grep("^msfragger\\.table\\.fix-mods=", x, value = TRUE) }))
if (length(fixMods) > 1) { stop("Multiple searches with different fixed modifications are currently unsupported!") } # This could change, since we are only parsing those to extract PTM information
# Var mods
varMods <- unique(sapply(Workflows, function(x) { grep("^msfragger\\.table\\.var-mods=", x, value = TRUE) }))
if (length(varMods) > 1) {
  # We are combining them, because we want to have a table with all those mods
  tmp <- unique(unlist(strsplit(gsub("^msfragger\\.table\\.var-mods=", "", varMods), "; ")))
  tmp <- grep("^0\\.0+,", tmp, invert = TRUE, value = TRUE)
  varMods <- paste0("msfragger.table.var-mods=", paste(tmp, collapse = "; "))
}
# Was DiaNN run?
pat <- proteoCraft::topattern("diann.run-dia-nn=")
isActuallyDIANN <- unique(sapply(Workflows, function(x) { as.logical(toupper(gsub(pat, "", grep(pat, x, value = TRUE)))) }))
if (length(isActuallyDIANN) > 1) {
  isActuallyDIANN <- "diann.run-dia-nn=true"
} else {
  isActuallyDIANN <- paste0("diann.run-dia-nn=", tolower(as.character(isActuallyDIANN)))
}
# Is TMT
pat <- proteoCraft::topattern("tmtintegrator.run-tmtintegrator=")
isTMT <- unique(sapply(Workflows, function(x) { as.logical(toupper(gsub(pat, "", grep(pat, x, value = TRUE)))) }))
if (length(isTMT) > 1) { stop("This is not supported!") } else {
  TMT <- paste0("tmtintegrator.run-tmtintegrator=", tolower(as.character(isTMT)))
  if (isTMT) {
    pat <- proteoCraft::topattern("tmtintegrator.channel_num=")
    TMTvers <- unique(sapply(Workflows, function(x) { grep(pat, x, value = TRUE) }))
    if (length(TMTvers) > 1) { stop("This is not supported!") }
  }
}
# WorkDir
FP_Dir <- paste0("workdir=", dstDir)
#
# Create and save our fake workflow
nuWorkflow <- c(fpVers, "", "", FP_Dir, fastas, fixMods, varMods, isActuallyDIANN, TMT)
if (isTMT) { nuWorkflow <- c(nuWorkflow, TMTvers) }
nuWorkflow <- c(nuWorkflow, "")
nuWorkflow_fl <- paste0(dstDir, "/fragpipe.workflow")
write(nuWorkflow, nuWorkflow_fl)
#
#
# PSMs
Exp <- unique(nuManifest$Experiment)
Exp <- Exp[which(!is.na(Exp))]
if (length(Exp)) {
  for (exp in Exp) { #exp <- Exp[1]
    dir <- paste0(dstDir, "/", exp)
    if (!dir.exists(dir)) { dir.create(dir) }
    psms <- paste0(dirs, "/", exp, "/psm.tsv")
    w <- which(file.exists(psms))
    if (length(w)) {
      psms <- lapply(psms[w], function(x) {
        data.table::fread(x, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
      })
      psms <- plyr::rbind.fill(psms)
      data.table::fwrite(psms, paste0(dir, "/psm.tsv"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE, na = "")
    }
  }
} else {
  psms <- paste0(dirs, "/psm.tsv")
  w <- which(file.exists(psms))
  if (length(w)) {
    psms <- lapply(psms[w], function(x) {
      data.table::fread(x, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
    })
    psms <- plyr::rbind.fill(psms)
    data.table::fwrite(psms, paste0(dstDir, "/psm.tsv"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE, na = "")
  }
}

#
# Test
N.clust <- parallel::detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(parallel::stopCluster(parClust), silent = TRUE)
  parClust <- parallel::makeCluster(N.clust, type = "SOCK")
}
FP_Workflow <- nuWorkflow_fl
FP_Manifest <- nuManifest_fl
tst <- try(proteoCraft::FP_to_MQ(nuWorkflow_fl, nuManifest_fl, cl = parClust), silent = TRUE)
if (!"try-error" %in% class(tst)) {
  msg <- paste0("Done!\nNow you can treat folder:\n\t", normalizePath(dstDir), "\nas a (minimal) FragPipe pseudo-output folder for the purpose of running this package's data analysis scripts,\nand ONLY FOR THAT PURPOSE!\n")
  cat(msg)
} else {
  warning("This didn't work...")
}
