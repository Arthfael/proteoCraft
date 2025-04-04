library(Biostrings)
devtools::install_github("CambridgeCentreForProteomics/camprotR", dependencies = TRUE)
library(camprotR)
library(httr)
dir <- "D:/Fasta_databases"
if (!dir.exists(dir)) { dir.create(dir) }
ccp_tmp <- tempfile(fileext = ".fasta")
download_ccp_crap(ccp_tmp, is_crap = TRUE, verbose = TRUE)
file.copy(ccp_tmp, paste0(dir, "/CCP_cRAPome.fasta"), overwrite = TRUE)

# Despite CCP claims to the contrary, the headers are actually slightly different than would be expected from a standard fasta. Let's fix that!
library(proteoCraft)
db <- Format.DB(paste0(dir, "/CCP_cRAPome.fasta"))
db$`Full ID`[1:10]
db$`Protein ID` <- paste0("CON__", gsub("\\|", "", gsub("^[a-z]{2}\\||\\|[^\\|]+$", "", db$`Full ID`)))
db$`Protein ID`[1:10]
stopifnot(length(unique(db$`Protein ID`)) == nrow(db))
db$Header <- apply(db[, c("Header", "Protein ID")], 1, function(x) {
  y <- unlist(strsplit(x[[1]], "\\|[^\\|]+\\|[^\\|]+\\|"))
  return(paste0(y[1], "|", x[[2]], "|", y[2]))
})
writeFasta(db, paste0(dir, "/CCP_cRAPome.fasta"))

