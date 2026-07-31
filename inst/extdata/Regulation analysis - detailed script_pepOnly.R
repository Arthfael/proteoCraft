#### Code chunk - Initialization
if (!interactive()) { stop("This script should only be run within an interactive R session!") }
options(stringsAsFactors = FALSE)
options(install.packages.compile.from.source = "never")
options(svDialogs.rstudio = TRUE)
#rm(list = ls()[which(!ls() %in% c("dtstNm", "wd", "inDirs", "outDir"))])
closeAllConnections()

## Load proteoCraft
if (exists(".obj")) { rm(".obj") }
library(proteoCraft)
dirlist %<o% c() # This should go!!!
# Get local work directory:
ScriptPath %<o% normalizePath(gtools::script_file(), winslash = "/")
RunByMaster %<o% grepl(" - master script\\.R$", ScriptPath)
if (RunByMaster) { ScriptPath <- BehindTheScenes$ScriptFile }
Script %<o% readr::read_lines(ScriptPath)

RPath %<o% as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath %<o% paste0(RPath, "/proteoCraft")
homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")

fls <- paste0(homePath, "/", c(#"Regulation analysis - master script.R",
  "Regulation analysis - detailed script.R",
  "Regulation analysis - detailed script_pepOnly.R",
  "No replicates analysis - detailed script.R",
  "Reload_renv_from_lock_file.R",
  "Default_locations.xlsx",
  "LC_columns.xlsx"))
tst <- sum(!file.exists(fls))
if (tst) { proteoCraft::Configure() }
xplorSrc %<o% paste0(libPath, "/extdata/Sources/xplorData.R")
locDirs_fl %<o% paste0(homePath, "/Default_locations.xlsx")
locDirs %<o% openxlsx2::read_xlsx(locDirs_fl)

# Load backup?
load_a_Bckp %<o% c(TRUE, FALSE)[match(svDialogs::dlg_message("Re-load a backup?", "yesno")$res, c("yes", "no"))]
if (load_a_Bckp) {
  tst <- try({
    locDirs %<o% openxlsx2::read_xlsx(locDirs_fl)
    load_Bckp(startDir = locDirs$Path[match("Temporary folder", locDirs$Folder)])
  }, silent = TRUE)
  # Update values!
  xplorSrc %<o% paste0(libPath, "/extdata/Sources/xplorData.R")
  locDirs_fl %<o% paste0(homePath, "/Default_locations.xlsx")
  locDirs %<o% openxlsx2::read_xlsx(locDirs_fl)
}

if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1L, 1L)) }
parSrc %<o% paste0(libPath, "/extdata/Sources/make_check_Cluster.R")
bckpSrc %<o% paste0(libPath, "/extdata/Sources/updateBackup.R")


# Boolean functions to check parameter values
Src <- paste0(libPath, "/extdata/Sources/parBooleans.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
ReUseAnsw %<o% FALSE
scrptType %<o% "withReps"
scrptTypeFull %<o% "withReps_PTMs_only"
ExcelMax %<o% 32767L
MakeRatios %<o% TRUE

# Parameters used by the start analysis script:
###-|-### Workflows: setNames(c("Differential Peptide Expression analysis", "Pull-Down (e.g. co-IP)", "Biotin-based Pull-Down (BioID, TurboID, APEX...)", "Time Course", "SubCellular Localisation analysis"), c("REGULATION", "PULLDOWN", "BIOID", "TIMECOURSE", "LOCALISATION"))
###-|-### Replicates? TRUE
###-|-### External dependencies: Excel (loose); ScanHeadsman (loose); Cytoscape (loose); saintExpress (auto)

### Packages
## For convenience all (or most) of the packages used are loaded or installed here:
## CRAN packages:
if(!exists("cran_req")) { cran_req <- "pak" }
cran_req %<o% cran_req
if(!exists("bioc_req")) { bioc_req <- c() } 
bioc_req %<o% bioc_req
cran_req <- unique(c(cran_req, "pak", "fs", "shiny", "renv", "R.utils", "data.table", "devtools", "qs2", "shinyWidgets", "DT", "shinyBS", "stringr",
                     "gplots", "ggplot2", "ggpubr", "gtools", "reshape", "reshape2", "compiler", "stats", "rgl", "ggrepel", "rstudioapi", "modeest",
                     "minpack.lm", "snow", "viridis", "pcaMethods", "impute", "imputeLCMD", "parallel", "coin", "openxlsx", "openxlsx2", "plotly",
                     "Peptides", "xml2", "pdftools", "statmod", "ggpolypath", "venn", "gridExtra", "svDialogs", "htmlwidgets", "magrittr", "tibble",
                     "officer", "hexbin", "igraph", "matlib", "umap", "plyr", "ggnewscale", "shinyjs", "shinyFiles", "TeachingDemos", "shinycssloaders",
                     "tidyr", "ggplotify", "jpeg", "scattermore", "rpanel", "stringi", "lmtest", "ssh", "taxize", "arrow", "PTMods",
                     "ggdendro", "colorspace", "factoextra", "NbClust", "BH", "DEP", "iq", "Rtsne"))
bioc_req <- unique(c(bioc_req, "biomaRt", "GO.db", "UniProt.ws", "limma", "sva", "qvalue", "MSnbase",
                     "Rgraphviz", "RCy3", "siggenes", "DEqMS", "limpa", "QFeatures"))
inst <- as.data.frame(installed.packages())
for (pack in cran_req) {
  if (!pack %in% inst$Package) {
    if (pack %in% c("pak", "uchardet", "taxize")) {
      # Exceptions where for now we want a specific version to be installed,
      # or have to help the installer so it finds the right location
      if (pack == "pak") {
        install.packages("pak", dependencies = TRUE)
      }
      if (pack == "uchardet") {
        url <- "https://cran.r-project.org/src/contrib/Archive/uchardet/uchardet_1.1.1.tar.gz"
        destfile <- "uchardet_1.1.1.tar.gz"
        tst <- try(download.file(url, destfile, "curl"), silent = TRUE)
        if (inherits(tst, "try-error")) { try(download.file(url, destfile, "wget"), silent = TRUE) }
        install.packages(destfile, dependencies = TRUE)
        unlink(destfile)
      }
      if (pack == "taxize") {
        pak::pak("ropensci/bold", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
        pak::pak("ropensci/taxize", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
      }
    } else {
      tst <- try(pak::pak(pack, ask = FALSE, upgrade = TRUE, dependencies = TRUE), silent = TRUE)
      if (inherits(tst, "try-error")) {
        tst <- try(install.packages(pack, dependencies = TRUE), silent = TRUE)
      }
      if (inherits(tst, "try-error")) {
        warning(paste0("Package ", pack, " wasn't installed properly, skipping..."))
        cran_req <- cran_req[which(cran_req != pack)]
        bioc_req <- bioc_req[which(bioc_req != pack)]
      }
    }
    inst <- as.data.frame(installed.packages())
  }
}
## Bioconductor packages:
biocInstall %<o% function(pack, load = TRUE) {
  inst <- as.data.frame(installed.packages())
  if (!pack %in% inst$Package) {
    tst <- try(pak::pak(pack, ask = FALSE, upgrade = TRUE, dependencies = TRUE), silent = TRUE)
    if ("try-error"%in% class(tst)) {
      tst <- try(pak::pak(pack, ask = FALSE, upgrade = TRUE, dependencies = FALSE), silent = TRUE)
    }
    if ("try-error"%in% class(tst)) {
      tst <- try(pak::pak(pack, ask = FALSE, upgrade = FALSE, dependencies = FALSE), silent = TRUE)
    }
    if ("try-error"%in% class(tst)) {
      stop(tst)
    }
  }
  if (load) { library(pack, character.only = TRUE) }
}
for (pack in bioc_req) { biocInstall(pack, load = FALSE) }

# Run local scripts at startup - keep this after loading the backup!
locScrptSrc %<o% paste0(libPath, "/extdata/Sources/runLocScrpts.R")
#rstudioapi::documentOpen(locScrptSrc)
source(locScrptSrc)

# Set Shiny options, load functions for creating a Word report, create Excel styles
Src <- paste0(libPath, "/extdata/Sources/ShinyOpt_Styles_and_Report.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Select input/output folders and define experimental structure
# Reuse answers?
# The script sometimes pauses to ask the user a question in a popup. These answers are stored.
# When re-running the script, it can be useful to just answer whether one would like to re-use the answers to all of these?
if (file.exists("All_decisions.RData")) {
  msg <- "This script sometimes pauses to ask the user a question in a popup.
These answers are saved locally.
Answers from a previous run could be found in the analysis folder, do you want to re-use them?
"
  ReUseAnsw <- c(TRUE, FALSE)[match(svDialogs::dlg_message(msg, "yesno")$res, c("yes", "no"))]
  if (is.na(ReUseAnsw)) { ReUseAnsw <- FALSE }
}
if (ReUseAnsw) {
  load("All_decisions.RData")
  AllAnsw <- AllAnsw
}
if (!exists("AllAnsw")) {
  AllAnsw <- data.frame(Parameter = "Which question is it?", Message = "Message of the question")
  AllAnsw$Value <- list("Value of the answer")
}
AllAnsw %<o% AllAnsw

# Update the proteoCraft package?
# msg <- "Should we update the proteoCraft package?"
# updt_proteoCraft %<o% c(TRUE, FALSE)[match(svDialogs::dlg_message(msg, "yesno")$res, c("yes", "no"))]
updt_proteoCraft %<o% FALSE

# Define input, output, project folder etc...
Src <- paste0(libPath, "/extdata/Sources/Start_analysis.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

setwd(wd)
#loadFun(BckUpFl)
#openwd()

# Log the current analysis:
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
start_date %<o% gsub(":", "-", gsub(" ", "_", Sys.time()))
if (!"Current analysis start data.RData" %in% list.files(paste0(wd, "/Workflow control"))) {
  save(start_date, file = paste0(wd, "/Workflow control/Current analysis start data.RData"))
}
#load(paste0(wd, "/"Workflow control/Current analysis start data.RData"))
if (! paste0(wd, "/Workflow control/Data_analysis_log_", start_date, ".txt") %in% list.files()) {
  write(c(paste0("Data_analysis_log_", start_date), "__________________________", ""),
        file = paste0(wd, "/Workflow control/Data_analysis_log_", start_date, ".txt"))
}
#logcon %<o% file(paste0("Workflow control/Data_analysis_log_", start_date, ".txt"), open = "a")
#sink(logcon, type = "message", split = TRUE)
sink(paste0("Workflow control/Data_analysis_log_", start_date, ".txt"), type = "output", append = TRUE, split = TRUE)

# Create parallel processing cluster
source(parSrc, local = FALSE)
setDTthreads(threads = N.clust)

# Load PSMs
Src <- paste0(libPath, "/extdata/Sources/Load_PSMs.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Install and/or load rawrr only if we have .raw files
g <- grep("\\.raw$", rawFiles, ignore.case = TRUE)
if (length(g)) {
  rawrrSrc <- paste0(libPath, "/extdata/Sources/install_rawrr.R")
  #rstudioapi::documentOpen(rawrrSrc)
  source(rawrrSrc)
}

# MS raw files map
Src <- paste0(libPath, "/extdata/Sources/Fractions_Map_editor.R")
#rstudioapi::documentOpen(Src)
tstFrMp <- FALSE
while (!tstFrMp) {
  source(Src, local = FALSE)
}

#### Code chunk - Edit Experimental Factors
Src <- paste0(libPath, "/extdata/Sources/Experimental_Factors_editor.R")
#rstudioapi::documentOpen(Src)
tstXpFct <- FALSE
while (!tstXpFct) {
  source(Src, local = FALSE)
}
#

#### Code chunk - Edit Experiment map
Src <- paste0(libPath, "/extdata/Sources/Experiment_Map_editor.R")
#rstudioapi::documentOpen(Src)
tstXpMp <- FALSE
while (!tstXpMp) {
  source(Src, local = FALSE)
}
#

#### Code chunk - Load and process search database(s)
Src <- paste0(libPath, "/extdata/Sources/Process_Fasta_DBs.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#evNm %<o% c("PSM", "Evidence")[(SearchSoft == "MAXQUANT")+1L]
evNm %<o% "PSM"

#### Code chunk - Load and process annotations
Src <- paste0(libPath, "/extdata/Sources/Load_Annotations.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
source(parSrc, local = FALSE)
Src <- paste0(libPath, "/extdata/Sources/GO_prepare.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Create experiment Factors shortcuts
Src <- paste0(libPath, "/extdata/Sources/XpFact_shortcuts.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

# Check and process Fractions map
Src <- paste0(libPath, "/extdata/Sources/Fractions_Map_check.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - MA plots //ask
Src <- paste0(libPath, "/extdata/Sources/MA_plots.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Define analysis parameters
# Define analysis parameters
paramSrc <- paste0(libPath, "/extdata/Sources/rep_Parameters_editor_Main.R")
#rstudioapi::documentOpen(paramSrc)
source(paramSrc, local = FALSE)
#
Src <- paste0(libPath, "/extdata/Sources/rep_Parameters_editor_Contr.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#
Src <- paste0(libPath, "/extdata/Sources/PTMs_check.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Start writing Materials and Methods
Src <- paste0(libPath, "/extdata/Sources/autoMatMet.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Start processing the PSMs table
ReportCalls <- AddSpace2Report()
ReportCalls <- AddTxt2Report("Processing PSMs...")
# Remove reverse database hits
ev <- ev[which(ev$Reverse == ""),]

# Optionally remove charge 1 PSMs - off for now, but may become either user decision or parameter controlled
RemovZ1 <- FALSE
w1 <- which(ev$Charge == 1L)
wHt1 <- which(ev$Charge > 1L)
if ((RemovZ1)&&(length(w1))) {
  AmIBogus <- paste(unique(ev$"Modified sequence"[w1]), collapse = "\n")
  #cat(AmIBogus)
  cat("Removing the following presumably bogus identifications with Z=1:\n", AmIBogus, "\n")
  ev <- ev[wHt1,]
}

w <- grep("CONTAMINANT", colnames(ev), ignore.case = TRUE)
if (length(w) > 1L) { warning("Hmmm..., you might wanna check what is happening here...") } else {
  colnames(ev)[w] <- "Potential contaminant"
}
for (i in c("Potential contaminant", "Reverse")) {
  w <- which(is.na(ev[[i]]))
  ev[w, i] <- ""
}

#### Code chunk - Update peptide-to-protein mappings
Src <- paste0(libPath, "/extdata/Sources/checkPep2Prot.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

# Filter to keep only PSMs with valid quantitative values:
Src <- paste0(libPath, "/extdata/Sources/filtPSMs.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

rm(list = ls()[which(!ls() %in% .obj)])

#### Code chunk - Summary table and QC plots
Src <- paste0(libPath, "/extdata/Sources/rep_Summary.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Test for amino acid biases:
Src <- paste0(libPath, "/extdata/Sources/AA_biases_test.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Isobaric label purity correction
Src <- paste0(libPath, "/extdata/Sources/isobarCorr.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

ev$"Unique State" <- do.call(paste, c(ev[, c("Modified sequence", "Charge")], sep = ""))

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

# DIA-only: MS2-based correction of MS1-based quantitative values
Src <- paste0(libPath, "/extdata/Sources/MS2corr2MS1.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Optional - Normalize PSM MS1 intensities, then, if applicable, MS2 reporter (Isobaric labelling) or fragment (DIA) intensities
Src <- paste0(libPath, "/extdata/Sources/evNorm.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - ROC analysis for optimizing a threshold to include/exclude peptides mapped to specific GO terms
Src <- paste0(libPath, "/extdata/Sources/ROC1.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Create modified peptides table
Src <- paste0(libPath, "/extdata/Sources/pepMake.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)


#### Code chunk - Peptidoforms-level, calculate quantitation and test for outliers
Src <- paste0(libPath, "/extdata/Sources/pepQuant.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

####################################################
# Optional: choose whether to remove any outliers  #
####################################################
Src <- paste0(libPath, "/extdata/Sources/remove_Outliers.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

DatAnalysisTxt %<o% DatAnalysisTxt # Just in case...
l <- length(DatAnalysisTxt)
if (length(inDirs) > 1L) {
  DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                              " PSM tables ")
  DatAnalysisTxt[l] <- if (length(unique(SearchSoft)) > 1L) {
    paste0(DatAnalysisTxt[l],
           "from the different search engines were converted to a similar format and ")
  } else {
    paste0(DatAnalysisTxt[l],
           "were ")
  }
  DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                              "combined into a single table, then this ")
} else {
  DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                              " The long format PSMs table ")
}
DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                            "was transformed into a wide format peptidoforms table, summing up quantitative values where necessary.")

g <- grep(topattern(pep.ref["Original"]), colnames(pep), value = TRUE)
# View(pep[, g])
test <- rowSums(pep[, g])
l <- length(which(test == 0))
if (l) {
  msg <- paste0("Removing ", l, " peptide", c("", "s")[(l > 1L)+1L], " with invalid expression values - this is unexpected, investigate!")
  ReportCalls <- AddMsg2Report(Space = FALSE, Warning = TRUE)
  pep <- pep[which(test > 0),]
  w <- which(ev$id %in% unique(as.integer(unlist(strsplit(pep$"Evidence IDs", ";")))))
  ev <- ev[w,]
}

LocAnalysis2 %<o% FALSE

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Optionally impute missing peptide intensities
Src <- paste0(libPath, "/extdata/Sources/pep_Impute.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Re-normalize peptide intensities
rfnm <- c("Original", "Imputation")[Impute+1L]
Src <- paste0(libPath, "/extdata/Sources/pepNorm_VarPlot.R") # This may go, it currently isn't used much, and is kinda redundant with the MA plots...
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# - Run normalisations
if (Param$Norma.Pep.Intens) {
  Src <- paste0(libPath, "/extdata/Sources/pepNorm_Main.R")
  #rstudioapi::documentOpen(Src)
  #rstudioapi::documentOpen(nrmSrc)
  source(Src, local = FALSE)
}
#
rfnm <- names(pep.ref)[length(pep.ref)]
Src <- paste0(libPath, "/extdata/Sources/pepNorm_VarPlot.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#View(pep[, grep(topattern(pep.ref[length(pep.ref)]), colnames(pep), value = TRUE)]) # Check final data visually

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readr::read_lines(ScriptPath)

# Calculate peptide ratios - currently off
makePepRat %<o% FALSE
if (makePepRat) {
  Src <- paste0(libPath, "/extdata/Sources/rep_pepRat_Calc.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

# Backup data/update cluster
stopClust <- TRUE
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Annotations
# (still required for dataset enrichment analysis, and also for output tables; not for the GO terms enrichment analysis)
# I used to get functional annotations for all proteins in the protein group.
# However we are now - and I think with reason - only using annotations from the leading protein(s)!
p <- strsplit(PG$"Leading protein IDs", ";") #Here taking just the minimum set of protein IDs to explain the observed dataset.
db$Observed <- db$"Protein ID" %in% unique(unlist(p))
if (globalGO) {
  temp <- listMelt(strsplit(PG$"Leading protein IDs", ";"), PG$id)
  kol <- annot.col[which(annot.col %in% colnames(db))]
  if ("Taxonomy" %in% kol) {
    PG$Taxonomy <- db$Taxonomy[match(gsub(";.*", "", PG$`Leading protein IDs`), db$`Protein ID`)]
  }
  kol2 <- annot.col[which(!annot.col %in% "Taxonomy")]
  kol2 <- annot.col[which(annot.col %in% colnames(db))]
  temp[, kol2] <- db[match(temp$value, db$"Protein ID"), kol2]
  tst1 <- unlist(strsplit(temp$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(temp$GO, ";"))
  tst2a <- gsub(".*\\[|\\]$", "", tst2)
  tst3 <- data.table(A1 = tst1, A2 = tst2)
  tst3 <- tst3[, list(x = unique(A2)), by = list(Group.1 = A1)]
  tst3 <- as.data.frame(tst3)
  stopifnot(length(tst3$x) == length(unique(tst3$x)),
            is.character(tst3$x),
            length(which(temp$value != temp$Accession)) == 0L)
  for (i in kol2) { temp[[i]] <- strsplit(as.character(temp[[i]]), ";") }
  f0 <- function(x) { list(unique(unlist(x))) }
  temp <- aggregate(temp[, kol2], list(temp$L1), f0)
  #
  # Below commented data.table aggregation code... which is slower so not used.
  #
  # temp2 <- as.data.table(temp[, c("L1", kol2)])
  # temp2 <- temp2[, lapply(.SD, f0), by = list(Group.1 = L1), .SDcols = kol2]
  # temp2 <- as.data.frame(temp2)
  #
  for (i in kol2) {
    temp[[i]] <- parSapply(parClust, temp[[i]], \(x) { paste(unique(unlist(x)), collapse = ";") }) # Do not use sort here or it will break the correspondance between "GO" and "GO-ID"
  }
  tst1 <- unlist(strsplit(temp$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(temp$GO, ";"))
  tst3 <- data.table(A1 = tst1, A2 = tst2)
  tst3 <- tst3[, list(x = unique(A2)), by = list(Group.1 = A1)]
  tst3 <- as.data.frame(tst3)
  stopifnot(length(tst3$x) == length(unique(tst3$x)),
            is.character(tst3$x))
  #
  PG[, kol2] <- temp[match(PG$id, temp$Group.1), kol2]
  #
  #View(tst3[which(lengths(tst3$x) > 1L),])
  #View(tst3[which(lengths(tst3$x) == 0L),])
  #
  # Also peptides (minor approximation: use first protein group)
  pep[, kol] <- PG[match(as.integer(gsub(";.*", "", pep$`Protein group ID`)), PG$id), kol]
  #
  PG$Ontology <- NULL # Temporary fix for now, this column is broken
  #
  # It makes sense to close/re-create parallel clusters regularly to reduce memory usage + avoid corruption
  stopCluster(parClust)
  source(parSrc, local = FALSE)
  #
  Src <- paste0(libPath, "/extdata/Sources/GO_prepare.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

# Contrast log2FC columns
pep.rat.ref %<o% "log2(FC) - "
tmp <- make_Rat2(experiment.map = Exp.map,
                 int.root = pep.ref[length(pep.ref)],
                 rat.root = pep.rat.ref)
pep[, colnames(tmp)] <- tmp

#### Code chunk - Perform statistical tests
dataType <- "peptides"
Src <- paste0(libPath, "/extdata/Sources/Stat_tests.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

### Visualize and check P-values
dataType <- "peptides"
Src <- paste0(libPath, "/extdata/Sources/pVal_check.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

useSAM %<o% ((names(pvalue.col)[which(pvalue.use)] == "Student")&&(useSAM_thresh))

# Create list of control ratio values for the purpose of identifying vertical thresholds for plots:
Src <- paste0(libPath, "/extdata/Sources/ratThresh.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

create_plotly %<o% TRUE
create_plotly_local %<o% TRUE # No need for a licence when I can save local htmls! Still, old legacy code kept below.
# Arbitrary thresholds
arbitrary.thr %<o% data.frame(yintercept = -log10(c(0.05, 0.01)),
                              slope = c(0, 0),
                              xintercept = c(NA, NA),
                              colour = c("orange", "red"),
                              label = c("5% P-value", "1% P-value"))
# Default volcano plot arguments
Src <- paste0(libPath, "/extdata/Sources/dfltVolcPlotArgs.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

volcano.plots %<o% list()
filter_types %<o% tolower(unlist(strsplit(Param$Filters.type, ";")))
filter_types[grep("^dat.+2$", filter_types, invert = TRUE)] <- substr(filter_types[which(!grepl("^dat.+2$", filter_types))], 1L, 3L)
filter_types[grep("^dat.+2$", filter_types)] <- "dat2"
filter_types <- unique(c("con", filter_types))
if ("ref" %in% filter_types) {
  if ((RRG$aggregate != RG$aggregate)||(Nested)) {
    warning("Grouping filter by reference is not feasible if replicates are paired!")
    filter_types <- filter_types[which(filter_types != "ref")]
  } else {
    if (sum(vapply(RG$names, \(x) {! x %in% RSA$names }, TRUE)) > 0L) {
      warning("Grouping filter by reference is not feasible if the factors used for \"RG\" are not included in those used for \"RRG\"!")
      filter_types <- filter_types[which(filter_types != "ref")]
    }
  }
}
F_Root %<o% "mod. F-test -log10(Pvalue)"
GO.enrich.MultiRefs %<o% (("GO.enrichment.Ref.Aggr" %in% colnames(Param))&&(!Param$GO.enrichment.Ref.Aggr %in% c("", "NA", NA)))
F.test %<o% Param$F.test

if (F.test) {
  EM <- Exp.map[match(rownames(designMatr), Exp.map$Ref.Sample.Aggregate),]
  Group_ <- do.call(paste, c(EM[, Coefficients, drop = FALSE], sep = "_"))
  Group_ <- as.factor(Group_)
  EM$Group_ <- Group_
  expContrasts_F %<o% expContrasts
  expContrasts_F$Type <- "Simple"
  expContrasts_F$Contrasts <- tmp <- apply(expContrasts_F[, c("x1", "x0")], 1L, \(x) {
    x <- x[which(x %in% colnames(designMatr))]
    paste(x, collapse = " - ")
  })
  expContrasts_F$Map <- NULL
  # Double contrasts - for now we create all by default
  # Eventually we will have an app with all contrasts:
  # - Select sample group reference(s?) per ratio group
  # - Generate automatically all contrasts from references (click button to update contrasts)
  # - Choose which to run for (all) t-tests and F-tests
  l <- length(tmp)
  if (l > 1L) {
    tmp2 <- unlist(vapply(1L:(l-1L), \(x) {
      paste0("(", tmp[x], ") - (", tmp[(x+1L):l], ")")
    }, ""))
    tmp2Nm <- unlist(vapply(1L:(l-1L), \(x) {
      paste0("(", expContrasts_F$name[x], ") - (", expContrasts_F$name[(x+1L):l], ")")
    }, ""))
    tmp2Tbl <- data.frame(Contrasts = tmp2,
                          x1 = gsub("^\\(|\\) - \\(.*", "", tmp2),
                          x0 = gsub(".*\\) - \\(|\\)$", "", tmp2),
                          name = tmp2Nm,
                          Type = "Double")
    tmp2Tbl$All <- lapply(1L:nrow(tmp2Tbl), \(x) { gsub("^Group_", "", unlist(strsplit(unlist(tmp2Tbl[x, c("x1", "x0")]), " - "))) })
    tmp <- c(tmp, tmp2)
    expContrasts_F <- rbind(expContrasts_F, tmp2Tbl)
  }
}

#### Code chunk - Modified peptides analysis
modPepSrc <- paste0(libPath, "/extdata/Sources/modPeptides.R")
#rstudioapi::documentOpen(modPepSrc)
source(modPepSrc, local = FALSE)

# Backup data/update cluster
stopClust <- TRUE
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Create output tables
## PSMs
dir <- paste0(wd, "/Tables")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
w <- which(vapply(colnames(ev), \(x) { inherits(ev[[x]], "list") }, TRUE))
if (length(w)) { for (i in w) { ev[[i]] <- parSapply(parClust, ev[[i]], paste, collapse = ";") } }
data.table::fwrite(ev, paste0(dir, "/evidence.tsv"), sep = "\t", row.names = FALSE, na = "NA")
#
## Main peptidoforms- and protein groups-level, multi-tabs report
xlSrc <- paste0(libPath, "/extdata/Sources/rep_Write_Excel_pepOnly.R")
#rstudioapi::documentOpen(xlSrc)
source(xlSrc, local = FALSE)
#xl_open(repFl)

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
Script <- readr::read_lines(ScriptPath)

#### Code chunk - Venn diagrams
vennTst <- try({
  Src <- paste0(libPath, "/extdata/Sources/rep_Venn.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}, silent = TRUE)
if (inherits(vennTst, "try-error")) {
  warning("Venn diagrams failed ---> investigate!")
}

#### Code chunk - Finalize analysis and export results
# Remove empty directories:
#dirlist <- list.dirs()
dirlist <- dirlist[order(nchar(dirlist), decreasing = TRUE)]
for (dir in dirlist) { #d <- dirlist[1L]
  if (!length(list.files(dir))) {
    unlink(dir, recursive = TRUE)
    dirlist <- dirlist[which(dirlist != dir)]
  }
}
# Save decisions
save(AllAnsw, file = "All_decisions.RData")

# Finalize reports
#
MatMetCalls$Texts$DatAnalysis <- c(MatMetCalls$Texts$DatAnalysis, DatAnalysisTxt)
for (i in seq_along(MatMetCalls$Texts$DatAnalysis)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls,
                              paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$DatAnalysis[", i,"], prop = WrdFrmt$",
                                     c("Body", "Template_text")[(MatMetCalls$Texts$DatAnalysis[i] == "TEMPLATE")+1L],
                                     "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#

# Write SDRF file in case you want to submit to PRIDE
Src <- paste0(libPath, "/extdata/Sources/SDRF_4_PRIDE.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Finalize analysis
Src <- paste0(libPath, "/extdata/Sources/Finalize_analysis.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

### That's it, done!
#openwd(outDir)
#rm(list = ls())
