# An R script designed to process a fasta to:
# - remove all sequence duplicates
# - remove entries with invalid amino acids
# - remove entries marqued as fragments
# - add common contaminants (CCP)
#
require(parallel)
require(proteoCraft)

# Create parallel processing cluster
# Cluster
RPath %<o% as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath %<o% paste0(RPath, "/proteoCraft")
parSrc %<o% paste0(libPath, "/extdata/R scripts/Sources/make_check_Cluster.R")
source(parSrc)

homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
locFl <- paste0(homePath, "/Default_locations.xlsx")
loc <- openxlsx2::read_xlsx(locFl)
fastaDir <- loc$Path[which(loc$Folder == "Fasta files")]

#
AA <- proteoCraft::AA
AA <- AA[which(!AA %in% c("O"))]
dbFl <- choose.files(paste0(fastaDir, "/*.fasta"), multi = FALSE) # Select fasta file
ext <- rev(unlist(strsplit(dbFl, "\\.")))[1] # In case the extension is different, e.g. .txt or .fa
db <- proteoCraft::Format.DB(dbFl, cl = parClust) # loads database as table and removes duplicate entries by default
tst <- nchar(gsub(paste(AA, collapse = "|"), "", db$Sequence))
wY <- which(tst == 0)
wN <- which(tst != 0)
l <- length(wN)
if (l) { warning(paste0("Removing ", l, " entries with non-standard amino acids.")) }
db <- db[wY,]
g <- grepl(" \\([Ff]ragment\\) OS=", db$Header)
if (sum(g)) { warning(paste0("Removing ", sum(g), " fragmentary entries.")) }
db <- db[which(!g),]
db2 <- proteoCraft::Format.DB(paste0(fastaDir, "/CCP_cRAPome.fasta")) # loads contaminants
#db <- plyr::rbind.fill(db, db2)
db <- rbind(db, db2)
proteoCraft::writeFasta(db, gsub(paste0("\\.", ext, "$"), paste0("_noDupl_cont.", ext), dbFl)) # Write fixed database with original extension
proteoCraft::openwd(dirname(dbFl))
