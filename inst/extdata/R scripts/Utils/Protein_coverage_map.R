# A simplified script for mapping
packs <- c("ggplot2", "svDialogs")
for (pack in packs) {
  if (!require(pack, character.only = TRUE)) { install.packages(pack) }
  require(pack, character.only = TRUE)
}
require(proteoCraft)

dir <- choose.dir("...Search_Folder/", "Select destination directory")
dir <- normalizePath(dir, winslash = "/")

Seq <- dlg_input("Enter protein sequence (single-letter code, capital letters only)")$res
#Seq <- "TESTPEPTIDERANDSOMESTUFF"
Nm <- dlg_input("Enter protein name")$res
#Nm <- "My test protein"
names(Seq) <- Nm

filt <- as.matrix(data.frame(Type = c("csv or tsv file"), Extension = c("*.csv;*.tsv")))
PepFl <- dlg_open(paste0(dir, "/*"), title = "Select csv file containing peptide modified sequences (\"Modified sequence\" column, \"_PEPT(ph)IDER_\" format) and optionally intensity values (\"Intensity\" column)",
                  multiple = FALSE, filters = filt)$res
ext <- gsub(".+\\.", "", PepFl)
if (ext == "tsv") { Pep <- read.csv(PepFl, check.names = FALSE, sep = "\t") }
if (ext == "csv") { Pep <- read.csv(PepFl, check.names = FALSE) }
if (ext == "xlsx") {
  if (!require("openxlsx", character.only = TRUE)) { install.packages("openxlsx") }
  require("openxlsx", character.only = TRUE)
  Pep <- read.xlsx(PepFl, check.names = FALSE)
}
Kol <- colnames(Pep)
dflt <- "Modified sequence"
ModSeqKol <- c()
while (!length(ModSeqKol)) { ModSeqKol <- dlg_list(Kol, dflt, title = "Select the modified sequence column:")$res }
Pep2 <- data.frame("Modified sequence" = Pep[[ModSeqKol]], check.names = FALSE)
dflt <- "Intensity"
IntKol <- dlg_list(Kol, dflt, title = "Select the Intensity column:")$res
if (length(IntKol)) { Pep2$Intensity <- Pep[[IntKol]] }
ttl <- paste0("Peptide coverage - ", Nm)
path <- paste0(dir, "/", gsub(paste(c("\\\\", "/", ":", "\\*", "?", "<", ">", "\\|"), collapse = "|"), "", ttl), ".jpg")
if ("Intensity" %in% colnames(Pep2)) {
  Coverage(Seq, Pep2$`Modified sequence`, Mode = "Align2", intensities = Pep2$Intensity, title = ttl, save = "jpg", save.path = path)
  Coverage(Seq, Pep2$`Modified sequence`, Mode = "Align2", intensities = Pep2$Intensity, title = ttl, save = "pdf", save.path = path)
} else {
  Coverage(Seq, Pep2$`Modified sequence`, Mode = "Align2", title = ttl, save = "jpg", save.path = path)
  Coverage(Seq, Pep2$`Modified sequence`, Mode = "Align2", title = ttl, save = "pdf", save.path = path)
}
