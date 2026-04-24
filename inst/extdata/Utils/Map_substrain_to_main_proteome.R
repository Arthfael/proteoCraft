# This script 
require(pwalign)
require(proteoCraft)
wd <- "D:/Fasta_databases/Marchantia_polymorpha"
setwd(wd)
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
# Create parallel processing cluster
parSrc <- paste0(libPath, "/extdata/Sources/make_check_Cluster.R")
#rstudioapi::documentOpen(parSrc)
source(parSrc, local = FALSE)

dbFl1 <- paste0(wd, "/MpTak1v5.1_r1.protein.fasta")
dbFl2 <- paste0(wd, "/Marchantia_polymorpha_UP_20260121_Iso.fasta")
# We updated these functions when writing this script, so we can load the latest versions
source("H:/aRmel_package/proteoCraft/R/Format.DB.R")
source("H:/aRmel_package/proteoCraft/R/Format.DB_worker.R")
#
db1 <- Format.DB(dbFl1, cl = parClust)
db2 <- Format.DB(dbFl2, cl = parClust)
stopifnot(length(unique(db1$Sequence)) == nrow(db1),
          length(unique(db2$Sequence)) == nrow(db2))
# (If this fails, this means we used an older version of Format.DB and have to add code to drop duplicates)

seq <- unique(db1$Sequence)
if (length(seq) < nrow(db1)) {
  wY <- match(seq, db1$Sequence)
  wN <- setdiff(1L:nrow(db1), wY)
  lN <- length(wN)
  cat(paste0(lN, " duplicate sequence", c(" was", "s were")[(lN > 1L)+1L], " removed.\n"))
  temp1 <- db1[wY, , drop = FALSE]
}

tst <- gsub(paste(AA, collapse = "|"), "", db1$Sequence)
#max(nchar(tst))
#unique(tst)
db1 <- db1[which(nchar(tst) == 0L),]
stopifnot(length(unique(db1$`Protein ID`)) == nrow(db1))

wY <- which(db1$Sequence %in% db2$Sequence)
wN <- which(!db1$Sequence %in% db2$Sequence)
tmp <- db2[match(db1$Sequence[wY], db2$Sequence),]
length(unique(db1$Sequence[wY]))
length(wY)

#View(db1[wN,])
stopifnot(length(unique(tmp$`Protein ID`)) == nrow(tmp))

db1[wY, ] <- tmp
#View(db1[wN,])
stopifnot(length(unique(db1$`Protein ID`)) == nrow(db1))

library(data.table)

# We want to find likely hits, but this will be difficult knowing that there are a lot of sequences to search
# Faster euristic:
# - Digest sequences
# - Only compare proteins sharing at least one peptide
dig1 <- Digest(setNames(db1$Sequence[wN], db1$`Protein ID`[wN]), cl = parClust)
dig2 <- Digest(setNames(db2$Sequence, db2$`Protein ID`), cl = parClust)
dig1DF <- listMelt(dig1, ColNames = c("Seq", "ID"))
dig2DF <- listMelt(dig2, ColNames = c("Seq", "ID"))
dig2DF <- dig2DF[which(dig2DF$Seq %in% dig1DF$Seq),]
dig2DF <- as.data.table(dig2DF)
dig2DF <- dig2DF[, list(ID = list(unique(ID))), by = list(Seq = Seq)]
dig2DF <- as.data.frame(dig2DF)
dig1DF$Candidates <- dig2DF$ID[match(dig1DF$Seq, dig2DF$Seq)]
dig1DF <- listMelt(setNames(dig1DF$Candidates, dig1DF$ID), ColNames =  c("ID2", "ID1"))
# dig1DF <- as.data.table(dig1DF)
# dig1DF <- dig1DF[, list(ID2 = list(unique(ID2))), by = list(ID1 = ID1)]
# dig1DF <- as.data.frame(dig1DF)
dig1DF$ID1_Seq <- db1$Sequence[match(dig1DF$ID1, db1$`Protein ID`)]
#dig1DF$ID2_Seq <- lapply(dig1DF$ID2, \(x) { db2$Sequence[match(x, db2$`Protein ID`)] })
dig1DF$ID2_Seq <- db2$Sequence[match(dig1DF$ID2, db2$`Protein ID`)]
# dig1DFLst <- apply(dig1DF[, c("ID1", "ID2", "ID1_Seq", "ID2_Seq")], 1, \(x) {
#   nms <- unlist(x[1:2])
#   seq <- unlist(x[3:4])
#   setNames(seq, nms)
# })
selfAlign <- rbind(as.data.frame(aggregate(dig1DF$ID1_Seq, list(dig1DF$ID1), unique)),
                   as.data.frame(aggregate(dig1DF$ID2_Seq, list(dig1DF$ID2), unique)))
colnames(selfAlign) <- c("ID", "Sequence")
dig1DF$ID1_Seq <- Biostrings::AAStringSetList(as.list(dig1DF$ID1_Seq))
dig1DF$ID2_Seq <- Biostrings::AAStringSetList(as.list(dig1DF$ID2_Seq))
selfAlign$Seq <- Biostrings::AAStringSetList(as.list(selfAlign$Sequence))
clusterEvalQ(parClust, { # Prevent nested parallelism (otherwise protr crashes or the computer becomes very, very sad - and slow!)
  options(mc.cores = 1L)
  Sys.setenv(
    OMP_NUM_THREADS = 1L,
    OPENBLAS_NUM_THREADS = 1L,
    MKL_NUM_THREADS = 1L,
    VECLIB_MAXIMUM_THREADS = 1L,
    NUMEXPR_NUM_THREADS = 1L
  )
})
tmpFl <- tempfile(".rds")
clusterExport(parClust, "tmpFl", envir = environment())
readr::write_rds(dig1DF, tmpFl)
clusterCall(parClust, \(x) {
  dig1DF <<- readr::read_rds(tmpFl)
  return(NULL)
})
dig1DF$"Similarity score" <- parSapply(parClust, 1L:nrow(dig1DF), \(x) { #x <- 1L
  pwalign::pairwiseAlignment(dig1DF$ID1_Seq[[x]],
                             dig1DF$ID2_Seq[[x]],
                             substitutionMatrix = "BLOSUM62",
                             gapOpening = 10L,
                             gapExtension = 0.5,
                             scoreOnly = TRUE)
})
readr::write_rds(selfAlign, tmpFl)
clusterCall(parClust, \(x) {
  selfAlign <<- readr::read_rds(tmpFl)
  return(NULL)
})
selfAlign$"Self similarity score" <- parSapply(parClust, 1L:nrow(selfAlign), \(x) { #x <- 1L
  pwalign::pairwiseAlignment(selfAlign$Seq[[x]],
                             selfAlign$Seq[[x]],
                             substitutionMatrix = "BLOSUM62",
                             gapOpening = 10L,
                             gapExtension = 0.5,
                             scoreOnly = TRUE)
})
stopCluster(parClust)
# Below: if parallelizing this you will need to re-export the objects, because the score columns do not exist in the versions currently stored on the cluster!
dig1DF$"Normalised similarity score" <- vapply(1L:nrow(dig1DF), \(x) { #x <- 1L
  sc <- dig1DF$`Similarity score`[x]
  id1 <- dig1DF$ID1[x]
  id2 <- dig1DF$ID2[x]
  slfSc1 <- selfAlign$`Self similarity score`[match(id1, selfAlign$ID)]
  slfSc2 <- selfAlign$`Self similarity score`[match(id2, selfAlign$ID)]
  return(sc/sqrt(slfSc1*slfSc2))
}, 1)
bestMatch <- aggregate(1L:nrow(dig1DF), list(dig1DF$ID1), \(x) {
  x[which.max(dig1DF$"Normalised similarity score"[x])]
})
colnames(bestMatch) <- c("ID", "best match row")
bestMatch$"best match" <- dig1DF$ID2[bestMatch$`best match row`]
bestMatch$"best match score" <- dig1DF$`Normalised similarity score`[bestMatch$`best match row`]
# Interpretation
# - 0.8       | same protein / isoform |
# - 0.6–0.8   | very likely same gene |
# - 0.45–0.6  | homolog — ambiguous mapping |
# - 0.3–0.45  | weak similarity (peptide collision possible) |
# - < 0.3     | probably unrelated
bestMatch$Decision <- "unrelated"
bestMatch$Decision[which(bestMatch$`best match score` >= 0.3)] <- "similar"
bestMatch$Decision[which(bestMatch$`best match score` >= 0.45)] <- "homolog"
bestMatch$Decision[which(bestMatch$`best match score` >= 0.6)] <- "same gene"
bestMatch$Decision[which(bestMatch$`best match score` >= 0.8)] <- "isoform"
bestMatch$Decision[which(bestMatch$`best match score` == 1)] <- "same protein" # (shouldn't be here, we have excluded identical sequences)
w0 <- which(bestMatch$Decision == "unrelated")
w1 <- which(bestMatch$Decision == "isoform")
w2 <- which(bestMatch$Decision == "same gene")
w3 <- which(bestMatch$Decision == "homolog")
w4 <- which(bestMatch$Decision == "similar")
bestMatch$`best match db2 row` <- match(bestMatch$`best match`, db2$`Protein ID`)
bestMatch$upID <- toupper(bestMatch$ID)
bestMatch$Header <- paste0(">mp|", bestMatch$upID, "|", bestMatch$upID, "_MARPO ")
bestMatch$Header[w0] <- paste0(bestMatch$Header[w0], "unknown protein")
bestMatch$Header[c(w1, w2)] <- paste0(bestMatch$Header[c(w1, w2)], "variant of ", db2$`Common Name`[bestMatch$`best match db2 row`[c(w1, w2)]])
bestMatch$Header[w3] <- paste0(bestMatch$Header[w3], "putative homolog of ", db2$`Common Name`[bestMatch$`best match db2 row`[w3]])
bestMatch$Header[w4] <- paste0(bestMatch$Header[w4], "similar to ", db2$`Common Name`[bestMatch$`best match db2 row`[w4]])
bestMatch$Header <- paste0(bestMatch$Header, " OS=Marchantia polymorpha subsp. ruderalis OX=1480154 ")
bestMatch$Header[c(w1, w2)] <- paste0(bestMatch$Header[c(w1, w2)], "GN=", db2$Gene[bestMatch$`best match db2 row`[c(w1, w2)]])
bestMatch$Header[c(w0, w3, w4)] <- paste0(bestMatch$Header[c(w0, w3, 24L == w4)], "GN=", bestMatch$ID[c(w0, w3, w4)])
bestMatch$Header <- paste0(bestMatch$Header, " PE=4 SV=1")
#bestMatch$Header[1:10]

wY <- which(db1$`Protein ID` %in% bestMatch$ID)
wN <- which((!db1$`Protein ID` %in% bestMatch$ID)&(!db1$Sequence %in% db2$Sequence))
db1$Header[wY] <- bestMatch$Header[match(db1$`Protein ID`[wY], bestMatch$ID)]
db1$Header[wN] <- paste0(">mp|", toupper(db1$`Protein ID`[wN]), "|", toupper(db1$`Protein ID`[wN]), "_MARPO unknown protein OS=Marchantia polymorpha subsp. ruderalis OX=1480154 GN=",
                         db1$`Protein ID`[wN], " PE=4 SV=1")

dbFl1_UP <- paste0(dirname(dbFl1), "/UPlike_", basename(dbFl1))
writeFasta(db1, dbFl1_UP)
