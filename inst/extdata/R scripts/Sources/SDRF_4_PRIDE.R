# This is an embryo script to write SDRF files for PRIDE submissions.
# For now it is probably incomplete, and only works for label-free datasets
# There will probably have to be additions for phospho-enriched samples too
# This could get included in the two main workflows, though I think it's better to keep it separate for now...
# Also the path to the validation tool is hard coded for now
# Assumes Trypsin
# Assumes single instrument type for all
# Assumes...
# Just... check the output!!!

require(svDialogs)
require(proteoCraft)

#load_Bckp()

stopifnot(LabelType %in% c("LFQ", "DIA")) # That part would be easy to write but isn't done yet
w <- which(c("Organism_Full", "Organism") %in% colnames(db))
tstorg %<o% (length(w) > 0)
if (tstorg) {
  OrgKol %<o% c("Organism_Full", "Organism")[w[1]]
  tst <- gsub(" *(\\(|\\[).*", "", db[[OrgKol]])
  tst <- aggregate(tst, list(tst), length)
  tst <- tst[order(tst$x, decreasing = TRUE),]
  mainOrg %<o% tst$Group.1[1]
} else { dlg_input("What is the parent organism of the samples?", "not available")$res }

setwd(wd)
if (scrptType == "withReps") {
  L <- nrow(Frac.map)
  Frac.map2 <- Frac.map
  kol <- colnames(Exp.map)[which(!colnames(Exp.map) %in% names(Aggregate.list))]
  kol <- kol[which(!kol %in% c("Fractions", "Reference", "MQ.Exp", "Sample.name", "Use", "Ref.Sample.Aggregate"))]
  if (length(unique(Exp.map$Experiment)) == 1) { kol <- kol[which(kol != "Experiment")] }
  tmp <- listMelt(Exp.map$MQ.Exp, 1:nrow(Exp.map))
  tmp[, kol] <- Exp.map[tmp$L1, kol]
  Frac.map2[, kol] <- tmp[match(Frac.map2$MQ.Exp, tmp$value), kol]
  SDRF <- data.frame("source name" = paste0("sample ", 1:L),
                     "characteristics[organism]" = tolower(mainOrg),
                     check.names = FALSE)
}
if (scrptType == "noReps") {
  L <- nrow(FracMap)
  Frac.map2 <- FracMap
  kol <- colnames(SamplesMap)
  kol <- kol[which(!kol %in% c("Fractions", "Reference", "MQ.Exp", "Sample.name", "Use"))]
  kol <- grep("__$", kol, value = TRUE, invert = TRUE)
  tmp <- listMelt(SamplesMap$Experiment, 1:nrow(SamplesMap))
  tmp[, kol] <- SamplesMap[tmp$L1, kol]
  Frac.map2[, kol] <- tmp[match(Frac.map2$`Parent sample`, tmp$value), kol]
  SDRF <- data.frame("source name" = paste0("sample ", 1:L),
                     "characteristics[organism]" = tolower(mainOrg),
                     check.names = FALSE)
}
#
if ("Tissue" %in% kol) { SDRF$"characteristics[organism part]" <- Frac.map2$Tissue } else {
  SDRF$"characteristics[organism part]" <- dlg_input("Which organism part was analysed in the experiment?", "not available")$res
}
kol <- kol[which(!kol %in% c("Tissue"))]
if ("Replicate" %in% kol) { SDRF$"characteristics[biological replicate]" <- Frac.map2$Replicate } else {
  SDRF$"characteristics[biological replicate]" <- 1
}
kol <- kol[which(!kol %in% c("Replicate"))]
if ("Disease" %in% kol) { SDRF$"characteristics[disease]" <- Frac.map2$Disease } else {
  SDRF$"characteristics[disease]" <- dlg_input("Which disease is of interest in the experiment?", "not available")$res
}
kol <- kol[which(!kol %in% c("Disease"))]
if ("Cell type" %in% kol) { SDRF$"characteristics[cell type]" <- Frac.map2$"Cell type" } else {
  SDRF$"characteristics[cell type]" <- dlg_input("From which cell type are the samples?", "not available")$res
}
kol <- kol[which(!kol %in% c("Cell type"))]
if ("Developmental stage" %in% kol) { SDRF$"characteristics[developmental stage]" <- Frac.map2$"Developmental stage" } else {
  SDRF$"characteristics[developmental stage]" <- dlg_input("At which development stage are the samples?", "not available")$res
}
kol <- kol[which(!kol %in% c("Developmental stage"))]
SDRF$"assay name" <- paste0("Run ", 1:L)
SDRF$"comment[technical replicate]" <- 1
SDRF$"comment[fraction identifier]" <- Frac.map2$Fraction
SDRF$"comment[label]" <- "label free"
SDRF$"comment[data file]" <- paste0(Frac.map2$`Raw files name`, ".d.zip")
if (!exists("homePath")) {
  homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
}
msModelDefFl <- paste0(homePath, "/MS_models.csv")
if (file.exists(msModelDefFl)) {
  msModelDef <- read.csv(msModelDefFl, check.names = FALSE)
  msModelDef$Vendor <- gsub(" +$", "", msModelDef$Vendor)
  vOpts <- unique(msModelDef$Vendor)
  vOpts <- c(grep("bruker", vOpts, value = TRUE, ignore.case = TRUE),
             grep("thermo", vOpts, value = TRUE, ignore.case = TRUE),
             grep("sciex", vOpts, value = TRUE, ignore.case = TRUE),
             grep("waters", vOpts, value = TRUE, ignore.case = TRUE),
             grep("agilent", vOpts, value = TRUE, ignore.case = TRUE),
             grep("shimadzu", vOpts, value = TRUE, ignore.case = TRUE),
             grep("((bruker)|(thermo)|(sciex)|(waters)|(agilent)|(shimadzu))", vOpts, value = TRUE, ignore.case = TRUE, invert = TRUE))
  vOpts <- vapply(vOpts, function(x) { paste(c(x, rep(" ", max(c(1, 250-nchar(x))))), collapse = "") }, "")
  vendor <- dlg_list(vOpts, vOpts[1], "Select MS instrument vendor from list")$res
  vendor <- gsub(" +$", "", vendor)
  iOpts <- msModelDef$Instrument[which(msModelDef$Vendor == vendor)]
  iOpts <- vapply(iOpts, function(x) { paste(c(x, rep(" ", max(c(1, 250-nchar(x))))), collapse = "") }, "")
  SDRF$"comment[instrument]" <- dlg_list(iOpts, iOpts[1], title = "Select MS instrument model")$res
  # See https://www.ebi.ac.uk/ols4/ontologies/ms/classes?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMS_1000031&viewMode=All&siblings=false for values
}
SDRF$"comment[cleavage agent details]" <- "NT=Trypsin; AC=MS:1001251; CS=(?⇐[KR])(?!P)"
if (length(kol)) { for (k in kol) { SDRF[[paste0("factor value[", tolower(k), "]")]] <- Frac.map2[[k]] } }
#View(SDRF)
sdrfPth <- paste0(outdir, "/SDRF.tsv")
#writeClipboard(sdrfPth)
write.table(SDRF, sdrfPth, sep = "\t", quote = FALSE, row.names = FALSE)
#openxlsx::openXL(sdrfPth)

fl <- paste0("", gsub("/Documents$", "", gsub("\\\\", "/", Sys.getenv("HOME"))),
             "/AppData/Roaming/Python/Python310/Scripts/parse_sdrf.exe")
if (file.exists(fl)) {
  cat("Validating SDRF file...\n")
  cmd <- paste0("\"", fl, "\" validate-sdrf --sdrf_file \"", sdrfPth, "\"")
  #cat(cmd)
  system(cmd)
  #openwd(dirname(sdrfPth))
} else {
  cat("To allow this workflow to automatically validate this file, make sure that a) python is installed and b) follow installation instructions for https://github.com/bigbio/sdrf-pipelines - the current script will then be able to run the validation tool directly...")
}
cat("\n... but remember: it's always a good idea to also check it manually!\n")
