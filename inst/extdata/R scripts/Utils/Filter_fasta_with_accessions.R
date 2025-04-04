# Create fasta from protein accession list
fastaFl <- selectFile("Select parent fasta database", path = "D:/Fasta_databases")
fastaFl <- gsub("^~", normalizePath(Sys.getenv("HOME"), winslash = "/"), fastaFl)
db <- proteoCraft::Format.DB(fastaFl)
svDialogs::dlg_message("Now copy to the clipboard the accessions you want to write a fasta for, then click ok", "ok")
intPrt <- readClipboard()
intPrt <- unlist(strsplit(intPrt, " *[,;] *"))
intPrt <- intPrt[which(intPrt != "")]
L1 <- length(intPrt)
intPrt <- intPrt[which(intPrt %in% db$`Protein ID`)]
L2 <- length(intPrt)
cat("Read ", L1, " accessions from the clipboard, ", c(L2, "all")[(L1 == L2)+1], " of which were found in the parent database.\n")
if (L2) {
  db <- db[match(intPrt, db$`Protein ID`),]
  writeFasta(db, newFl = TRUE)
} else { warning("No fasta file created!") }

