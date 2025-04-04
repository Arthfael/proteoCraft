#
libFl <- rstudioapi::selectFile("Select DiaNN library file",
                                path = getwd(),
                                filter = "DiaNN library (*.tsv)")
wd <- dirname(libFl)
setwd(wd)

# Check for SILAC mods:
psmsFl <- paste0(wd, "/report.tsv")
psms <- data.table::fread(psmsFl, integer64 = "numeric", check.names = FALSE,
                          data.table = FALSE)
mods <- grep("^UniMod", unique(unlist(strsplit(psms$Modified.Sequence, "\\(|\\)"))), value = TRUE)
mods

# Load library
#libFl <- paste0(wd, "/library.tsv")
libFl <- paste0(wd, "/report-lib.tsv")
lib <- data.table::fread(libFl, integer64 = "numeric", check.names = FALSE,
                         data.table = FALSE)
libKols <- colnames(lib)

# Load DIA-NN parameters and identify fasta database used
#require(proteoCraft)
#DIANNparam <- readLines(paste0(wd, "/report.log.txt"))
#DIANNcall <- grep("^diann.exe ", DIANNparam, value = TRUE)
#fasta <- gsub("\\.fasta.+", ".fasta", unlist(strsplit(DIANNcall, "--fasta "))[2])
#fasta <- normalizePath(fasta, winslash = "/")
#db <- Format.DB(fasta) # Load and parse search database

# R and K count (parent peptide)
for (i in c("R", "K")) {
  lib[[paste0(i, " count (peptide)")]] <- nchar(lib$PeptideSequence) - nchar(gsub(i, "", lib$PeptideSequence))
  cat(paste0("Max number of ", i, ": ", max(lib[[paste0(i, " count (peptide)")]]), "\n"))
}
lib$"RK count" <- rowSums(lib[, paste0(c("R", "K"), " count (peptide)")])
cat(paste0("Max number of R+K per peptide: ", max(lib$"RK count"), "\n"))
# R and K count (fragment)
lib$"Fragment sequence" <- ""
wb <- which(lib$FragmentType == "b")
wy <- which(lib$FragmentType == "y")
nc <- nchar(lib$PeptideSequence)
lib$"Fragment sequence"[wb] <- substr(lib$PeptideSequence[wb], 1, lib$FragmentSeriesNumber[wb])
lib$"Fragment sequence"[wy] <- substr(lib$PeptideSequence[wy], nc[wy]-lib$FragmentSeriesNumber[wy]+1, nc[wy])
for (i in c("R", "K")) {
  lib[[paste0(i, " count")]] <- nchar(lib$"Fragment sequence") - nchar(gsub(i, "", lib$"Fragment sequence"))
}  

# Define SILAC parameters:
require(unimod)
UniMods <- unimod::modifications
SILAC <- list("Medium" = setNames(c(6.020129, 4.025107), c("R", "K")),
              "Heavy" = setNames(c(10.008269, 8.014199), c("R", "K")))
SILAC2 <- lapply(SILAC, function(silac) { #silac <- SILAC[[1]]
  silac <- sapply(names(silac), function(i) { #i <- "R"
    UM <- UniMods[which(UniMods$Site == i),]
    w <- which(UM$MonoMass == silac[i])
    stopifnot(length(w) == 1)
    return(paste0("UniMod:", UM$UnimodId[w]))
  })
})
Libs <- setNames(lapply(names(SILAC), function(silac) { #silac <- "Medium"
  silac2 <- SILAC2[[silac]]
  silac <- SILAC[[silac]]
  tmp <- lib[which(lib$"RK count" > 0),]
  tmp$PrecursorMz <- tmp$PrecursorMz + silac["R"] * tmp$"R count (peptide)" + silac["K"] * tmp$"K count (peptide)"
  tmp$ProductMz <- tmp$ProductMz + silac["R"] * tmp$"R count" + silac["K"] * tmp$"K count"
  for (i in c("R", "K")) {
    tmp$ModifiedPeptide <- gsub("UniMod", "unimod", tmp$ModifiedPeptide)
    tmp$ModifiedPeptide <- gsub(i, paste0(i, "(", silac2[i], ")"), tmp$ModifiedPeptide)
    tmp$ModifiedPeptide <- gsub("unimod", "UniMod", tmp$ModifiedPeptide)
  }
  tmp$FullUniModPeptideName <- tmp$ModifiedPeptide # I checked, these are duplicate columns
  return(tmp)
}), names(SILAC))
Libs[["Light"]] <- lib
Libs <- Libs[c("Light", "Medium", "Heavy")]
Libs <- setNames(lapply(Libs, function(lb) {
  lb[, libKols]
}), c("Light", "Medium", "Heavy"))
SILACLib <- plyr::rbind.fill(Libs)
data.table::fwrite(SILACLib, paste0(wd, "/SILAC_library.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
                   na = "NA")
