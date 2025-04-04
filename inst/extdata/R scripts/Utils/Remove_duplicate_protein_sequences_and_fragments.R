# An R script designed to remove all duplicates, as well as X-containing 
require(parallel)
# Create parallel processing cluster
N.clust <- detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(stopCluster(parClust), silent = TRUE)
  parClust <- makeCluster(N.clust, type = "SOCK")
}
#
AA <- proteoCraft::AA
AA <- AA[which(!AA %in% c("O"))]
dbFl <- choose.files("D:/Fasta_databases/*.fasta", multi = FALSE) # Select fasta file
ext <- rev(unlist(strsplit(dbFl, "\\.")))[1] # In case the extension is different, e.g. .txt or .fa
db <- proteoCraft::Format.DB(dbFl, cl = parClust) # loads database as table and removes duplicate entries by default
tst <- nchar(gsub(paste(AA, collapse = "|"), "", db$Sequence))
wY <- which(tst == 0)
wN <- which(tst != 0)
l <- length(wN)
if (l) { warning(paste0("Removing ", l, " entries with non-standard amino acids.")) }
db <- db[wY,]
g <- grepl(" \\(Fragment\\) OS=", db$Header)
if (sum(g)) { warning(paste0("Removing ", sum(g), " fragmentary entries.")) }
db <- db[which(!g),]
db2 <- proteoCraft::Format.DB("D:/Fasta_databases/CCP_cRAPome.fasta") # loads contaminants
#db <- plyr::rbind.fill(db, db2)
db <- rbind(db, db2)
proteoCraft::writeFasta(db, gsub(paste0("\\.", ext, "$"), paste0("_noDupl_cont.", ext), dbFl)) # Write fixed database with original extension
proteoCraft::openwd(dirname(dbFl))
