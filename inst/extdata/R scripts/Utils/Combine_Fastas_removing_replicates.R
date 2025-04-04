#####################################################################################################
#                                                                                                   #
# This script combines all fastas in the destination folder below, removing any duplicate sequences #
#                                                                                                   #
#####################################################################################################

# Update the path below as needed:
wd <- "...My_Dataset/FRAGPIPE_F_theoretical_amount_inj/combine_fasta"

setwd(wd)
destFl <- paste0(wd, "/Combined_fasta.fasta") # The destination file

# Create parallel cluster
require(parallel)
N.clust <- parallel::detectCores()
a <- 1
tst <- try(parallel::clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(parallel::stopCluster(parClust), silent = TRUE)
  parClust <- parallel::makeCluster(N.clust, type = "SOCK")
}

# Scan folder for fasta files
fastaFls <- list.files(wd, "\\.fasta$", full.names = TRUE)
fastaFls <- fastaFls[which(fastaFls != destFl)]

DBs <- parSapply(parClust, fastaFls, readLines) # Read files,...
DBs <- unlist(DBs) # ... combine them as one character vector,...
db <- proteoCraft::Format.DB(DBs, in.env = TRUE, cl = parClust) # ... then process them as a data.frame, filtering out repeat sequences in the process.
# (not the only way to do this, but the fastest to write)

nrow(db) == length(unique(db$Sequence)) # Check that the result contains each sequence only once
proteoCraft::writeFasta(db, destFl) # Write results as fasta
