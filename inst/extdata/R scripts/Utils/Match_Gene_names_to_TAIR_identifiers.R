options(stringsAsFactors = FALSE)
require(proteoCraft)
setwd("D:/Fasta_databases")
list.files()
UniProt <- "Arabidopsis_thaliana_-_uniprot-proteome%3AUP000006548.fasta"
TAIR <- "TAIR10_pep_20101214.fasta"
UniProt <- Format.DB(UniProt)
TAIR <- Format.DB(TAIR, mode = "TAIR")
temp <- aggregate(UniProt$Gene, list(UniProt$Sequence), function(x) {paste(sort(unique(x)), collapse = ";") })
TAIR$"Gene (from UniProt)" <- ""
w <- which(TAIR$Sequence %in% temp$Group.1)
TAIR$"Gene (from UniProt)"[w] <- toupper(temp$x[match(TAIR$Sequence[w], temp$Group.1)])
temp <- TAIR[,c("Gene", "Gene (from UniProt)")]
for (i in colnames(temp)) { temp[[i]] <- strsplit(temp[[i]], ";") }
TAIR$"Gene (all)" <- toupper(apply(temp, 1, function(x) {
  paste(sort(unique(unlist(x))), collapse = ";")
}))

write.csv(TAIR, "TAIR as table with gene names.csv", row.names = FALSE)
openwd()
