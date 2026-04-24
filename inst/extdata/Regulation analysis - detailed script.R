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
ReUseAnsw %<o% FALSE
scrptType %<o% "withReps"
scrptTypeFull %<o% "withReps_PG_and_PTMs"
ExcelMax %<o% 32767L
MakeRatios %<o% TRUE

RPath %<o% as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath %<o% paste0(RPath, "/proteoCraft")
homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
parSrc %<o% paste0(libPath, "/extdata/Sources/make_check_Cluster.R")
bckpSrc %<o% paste0(libPath, "/extdata/Sources/updateBackup.R")
# Boolean functions to check parameter values
Src <- paste0(libPath, "/extdata/Sources/parBooleans.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

fls <- paste0(homePath, "/", c(#"Regulation analysis - master script.R",
  "Regulation analysis - detailed script.R",
  "Regulation analysis - detailed script_pepOnly.R",
  "No replicates analysis - detailed script.R",
  "Reload_renv_from_lock_file.R",
  "Default_locations.xlsx",
  "LC_columns.xlsx"))
tst <- sum(!file.exists(fls))
if (tst) { Configure() }
xplorSrc %<o% paste0(libPath, "/extdata/Sources/xplorData.R")
locDirs_fl %<o% paste0(homePath, "/Default_locations.xlsx")
locDirs %<o% openxlsx2::read_xlsx(locDirs_fl)

# Parameters used by the start analysis script:
###-|-### Workflows: setNames(c("Differential Protein Expression analysis", "Pull-Down (e.g. co-IP)", "Biotin-based Pull-Down (BioID, TurboID, APEX...)", "Time Course", "SubCellular Localisation analysis"), c("REGULATION", "PULLDOWN", "BIOID", "TIMECOURSE", "LOCALISATION"))
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
                     "ggdendro", "colorspace", "factoextra", "NbClust", "BH", "plogr", "iq", "Rtsne"))
bioc_req <- unique(c(bioc_req, "biomaRt", "GO.db", "UniProt.ws", "limma", "sva", "qvalue", "MSnbase", "DEP",
                     "Rgraphviz", "RCy3", "siggenes", "DEqMS", "pRoloc", "pRolocGUI", "rbioapi", "png", "Rhdf5lib", "limpa", "QFeatures"))
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
    if (inherits(tst, "try-error")) {
      tst <- try(pak::pak(pack, ask = FALSE, upgrade = TRUE, dependencies = FALSE), silent = TRUE)
    }
    if (inherits(tst, "try-error")) {
      tst <- try(pak::pak(pack, ask = FALSE, upgrade = FALSE, dependencies = FALSE), silent = TRUE)
    }
    if (inherits(tst, "try-error")) {
      stop(tst)
    }
  }
  if (load) { library(pack, character.only = TRUE) }
}
for (pack in bioc_req) { biocInstall(pack, load = FALSE) }

# Load backup?
load_a_Bckp %<o% c(TRUE, FALSE)[match(svDialogs::dlg_message("Do you want to load a backup?", "yesno")$res, c("yes", "no"))]
if (load_a_Bckp) {
  tst <- try({
    locDirs %<o% openxlsx2::read_xlsx(locDirs_fl)
    load_Bckp(startDir = locDirs$Path[match("Temporary folder", locDirs$Folder)])
  }, silent = TRUE)
}

# Run local scripts at startup - keep this after loading the backup!
locScrptSrc %<o% paste0(libPath, "/extdata/Sources/runLocScrpts.R")
source(locScrptSrc)

# Set Shiny options, load functions for creating a Word report, create Excel styles
Src <- paste0(libPath, "/extdata/Sources/ShinyOpt_Styles_and_Report.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Select input/output folders and define experimental structure
# Get local work directory:
ScriptPath %<o% normalizePath(gtools::script_file(), winslash = "/")
RunByMaster %<o% grepl(" - master script\\.R$", ScriptPath)
if (RunByMaster) { ScriptPath <- BehindTheScenes$ScriptFile }
Script %<o% readLines(ScriptPath)

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
Src <- paste0(libPath, "/extdata/Sources/GO_prepare.R") # Doing this earlier but also keep latter instance for now
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
if (LabelType == "LFQ") {
  source(parSrc, local = FALSE)
  clusterExport(parClust, "is.all.good", envir = environment()) # Use this (rather than loading library/calling by package::function syntax), in case we use a modified version of the function
  if ((Param$Label == "DIA")&&("MS2 intensities" %in% colnames(ev))) {
    ev$MS2_intensities <- strsplit(ev$"MS2 intensities", ";")
    ev$MS2_intensities <- parLapply(parClust, ev$MS2_intensities, as.numeric) # (Let's keep this as a numeric list)
    temp <- ev[, c(ev.col["Original"], "MS2_intensities")]
    temp$SumS2 <- parSapply(parClust, temp$MS2_intensities, \(x) { sum(is.all.good(x)) })
    # While we're at it, let's estimate missing MS1 intensities if we only have MS2:
    # (sum of MS2 intensities * median ratio of precursor intensity to sum of MS2 intensities)
    temp2 <- temp$Intensity/temp$SumS2
    m <- median(is.all.good(temp2))
    #sd(is.all.good(temp2))
    #plot <- ggplot(temp) + geom_point(aes(x = log10(Intensity), y = log10(SumS2))) + theme_bw() + geom_abline(intercept = log10(1/m), slope = 1, colour = "red")
    #poplot(plot)
    w <- which(((!is.all.good(ev[[ev.col["Original"]]], 2L))|(ev[[ev.col["Original"]]] <= 0))&(temp$SumS2 > 0))
    if (length(w)) { ev[w, ev.col["Original"]] <- temp$SumS2[w]*m }
    test <- parApply(parClust, temp[, c(ev.col["Original"], "SumS2")], 1L, sum)
  } else {
    temp <- ev[, ev.col["Original"], drop = FALSE]
    test <- parApply(parClust, temp, 1L, \(x) { sum(is.all.good(x)) })
  }
  l <- length(which(test == 0))
  if (l) {
    msg <- paste0("Removing ", l, " (", signif(100L*l/nrow(ev), 2L), "%) PSMs with invalid expression values!")
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
    w <- which(test > 0)
    ev <- ev[w,]
  }
}
if (LabelType == "Isobaric") { # If isobaric
  kol <- paste0(ev.ref["Original"], get(IsobarLab))
  w <- which(kol %in% colnames(ev))
  # Remove unused channels (if applicable)
  tmpIso <- get(IsobarLab)[w]
  u <- sort(as.numeric(unique(Exp.map$"Isobaric label")))
  w1 <- which(tmpIso %in% u)
  w2 <- which(!tmpIso %in% u)
  if (length(w2)) {
    kol <- lapply(tmpIso[w2], \(x) {
      paste0(c("Reporter intensity corrected ", "Reporter intensity ", "Reporter intensity count "), x)
    })
    w <- which(!colnames(ev) %in% unlist(kol))
    ev <- ev[, w]
  }
  assign(IsobarLab, tmpIso[w1])
  #
  kol <- paste0(ev.ref["Original"], get(IsobarLab))
  tst <- temp <- ev[, kol, drop = FALSE]
  tst$MS1 <- ev[[ev.col["Original"]]]
  tst$Reporter <- rowSums(temp, na.rm = TRUE)
  # Check dependency: there should be one in log space
  temp2 <- tst$MS1/tst$Reporter
  m <- median(is.all.good(temp2))
  #sd(is.all.good(temp2))
  #plot <- ggplot(tst) + geom_point(aes(x = log10(MS1), y = log10(Reporter))) + theme_bw() + geom_abline(intercept = log10(1/m), slope = 1, colour = "red")
  #poplot(plot)
  #
  #View()
  # If precursor intensity is missing, replace by estimate (sum of reporter intensities * median ratio of precursor intensity to sum of reporter intensities)
  w <- which((!is.all.good(ev[[ev.col["Original"]]], 2L))|(ev[[ev.col["Original"]]] <= 0))
  if (length(w)) { ev[w, ev.col["Original"]] <- tst$Reporter[w]*m }
  # Now the reverse scenario: no reporters, but we have precursor intensities; these are throw-away stuff 
  w <- which(!is.all.good(tst$Reporter, 2L)|(tst$Reporter <= 0))
  l <- length(w)
  if (l) {
    RemEv %<o% ev[w,]
    #View(RemEv[, kol])
    msg <- paste0("Removing ", l, " (", signif(100L*l/nrow(ev), 2L), "%) PSMs with invalid expression values!")
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
    w <- which(is.all.good(tst$Reporter, 2L)&(tst$Reporter > 0))
    ev <- ev[w,]
  }
}
rm(list = ls()[which(!ls() %in% .obj)])

# Filter by intensity
if ((!is.na(minInt))&&(is.numeric(minInt))&&(is.finite(minInt))&&(minInt >= 0)) {
  wY <- which(ev$Intensity >= minInt)
  wN <- which(ev$Intensity < minInt)
  lN <- length(wN)
  if (lN) {
    warning(paste0("Removing ", lN, " PSMs with intensity lower than ", minInt, " minimum threshold...\n"))
    ev <- ev[which(ev$Intensity >= minInt),]
  }
}

#### Code chunk - MA plots //ask
Src <- paste0(libPath, "/extdata/Sources/MA_plots.R")
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
DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                            "was transformed into a wide format peptidoforms table, summing up quantitative values where necessary.")

LocAnalysis2 %<o% FALSE

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Impute missing peptide intensities
Src <- paste0(libPath, "/extdata/Sources/pep_Impute.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Re-normalize peptide intensities
rfnm <- c("Original", "Imputation")[Impute+1L]
Src <- paste0(libPath, "/extdata/Sources/pepNorm_VarPlot.R")
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
# - Check variance/intensity dependency before or after normalisation
Src <- paste0(libPath, "/extdata/Sources/pepNorm_VarPlot.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#View(pep[, grep(topattern(pep.ref[length(pep.ref)]), colnames(pep), value = TRUE)]) # Check final data visually

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)

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

#### Code chunk - Assemble protein groups
ReportCalls <- AddSpace2Report()
ReportCalls <- AddTxt2Report("Starting protein groups assembly...")
if ("N. of peptidoforms for quantitation" %in% colnames(Param)) {
  N_Pep <- as.integer(Param$"N. of peptidoforms for quantitation")
  if ((is.na(N_Pep))||(N_Pep <= 0L)) {
    warning("Invalid `\"N. of peptidoforms for quantitation\" parameter, defaulting to 1!")
    N_Pep <- 1L
  }
} else { N_Pep <- 1L }
.obj <- unique(c("N_Pep", .obj))
#
tm1 <- Sys.time()
source(parSrc, local = FALSE)
Src <- paste0(libPath, "/extdata/Sources/PG_Assemble.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#loadFun("PG_assembly.RData")
#
PG %<o% PG_assembly$Protein.groups
pep <- PG_assembly$Peptides
db <- PG_assembly$Database
if ("Evidences" %in% names(PG_assembly)) { ev <- PG_assembly$Evidences }

# Here would be a good place to check protein taxonomy and [if necessary/as per parameters] split hybrid groups!
# Should be controlled by a parameter only showing up if taxonomy is present in db and has more than one value!
warning("(TO DO: add 'split-by-taxonomy?' here!)")
# Default = TRUE
# Don't forget to update PG IDs in PG, ev and pep afterwards! Also check potential contaminant column!


# Check those rare proteins IDs which are not in the search DB (should be contaminants)
tst <- unlist(strsplit(pep$Proteins, ";"))
if (length(tst)) {
  msg <- paste0("These protein accessions in peptides are not in the database: ", paste(tst[which(!tst %in% db$`Protein ID`)], collapse = " - "))
  ReportCalls <- AddMsg2Report(Space = FALSE, Print = FALSE)
}

# Basic fix, because I do not like the way I was doing Quality filters up to now
g <- grep("^Quality filter: ", colnames(PG), value = TRUE)
if (length(g)) {
  for (h in g) { #h <- g[1L]
    PG[[h]] <- c("no -> dubious!", "")[match(PG[[h]], c("", "Keep"))]
  }
}
#
if (tstOrg) {
  test <- vapply(strsplit(PG$`Protein IDs`, ";"), \(x) {
    paste(sort(unique(c(db[match(x, db$`Protein ID`), dbOrgKol]))), collapse = ";")
  }, "")
  pgOrgKol %<o% c("Organism", "Organism(s)")[(sum(grepl(";", test)) > 0L)+1L]
  PG[[pgOrgKol]] <- test
}
l <- length(DatAnalysisTxt)
DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                            " Protein groups were inferred from observed peptides.")

# Some stats on protein groups
tmp <- aggregate(PG$id, list(PG$`Peptides count`), length)
colnames(tmp) <- c("Peptides count", "Protein groups")
tmp$"log10(Protein groups count)" <- log10(tmp$"Protein groups")
pal <- colorRampPalette(c("brown", "yellow"))(max(tmp$"Peptides count")-1L)
tmp$Colour <- c("blue", pal)[tmp$`Peptides count`]
tmp2 <- summary(PG$`Peptides count`)
tmp2 <- data.frame(Variable = c(names(tmp2), "", "Protein groups", "Protein groups with 2+ peptidoforms"),
                   Value = c(as.character(signif(as.numeric(tmp2), 3L)),
                             "",
                             as.character(c(nrow(PG), sum(PG$"Peptides count" >= 2)))))
tmp2$Txt <- apply(tmp2[, c("Variable", "Value")], 1L, \(x) {
  x <- x[which(x != "")]
  x <- if (length(x)) { paste(x, collapse = ": ") } else { "" }
  return(x)
})
tmp2$X <- max(as.numeric(tmp2$Value[match("Max.", tmp2$Variable)]))*0.98
tmp2$Y <- max(tmp$"log10(Protein groups count)")*(0.98-(0L:(nrow(tmp2) - 1L))*0.02)
ttl <- "Peptidoforms per PG"
plot <- ggplot(tmp) + geom_col(aes(x = `Peptides count`, y = `log10(Protein groups count)`, fill = Colour),
                               colour = NA) +
  geom_text(data = tmp2, aes(x = X, y = Y, label = Txt), hjust = 1L, size = 3L) +
  scale_fill_identity() + theme_bw() + ggtitle(ttl)
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
dir <- paste0(wd, "/Summary plots")
dirlist<- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300L)
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L)
})

source(parSrc, local = FALSE)
tmp1 <- strsplit(pep$Proteins, ";")
tmp2 <- db[which(db$"Protein ID" %in% unlist(tmp1)), c("Protein ID", "Gene", "Name", "Common Name")]
exports <- list("tmp1", "tmp2")
clusterExport(parClust, "tmp2", envir = environment())
tmp <- parSapply(parClust, tmp1, \(x) {
  m <- match(unlist(x), tmp2$"Protein ID")
  x1 <- tmp2$Gene[m]
  x2 <- tmp2$Name[m]
  x3 <- paste(x1, collapse = ";")
  x4 <- paste(x2, collapse = ";")
  x1 <- paste(unique(x1), collapse = ";")
  x2 <- paste(unique(x2), collapse = ";")
  x5 <- paste(tmp2$`Common Name`[m], collapse = ";")
  rm(m)
  return(c(x1, x2, x3, x4, x5))
})
pep[, c("Gene names", "Protein names", "Gene names (all)", "Protein names (all)", "Common protein names")] <- t(tmp)
kol <- c("Leading proteins", "Leading razor proteins", "Gene names", "Protein names", "Gene names (all)",
         "Protein names (all)", "Common protein names", "Protein group IDs")
ev[, kol] <- pep[match(ev$"Modified sequence", pep$"Modified sequence"), kol]

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)

# Some more columns
tmp <- strsplit(PG$"Leading protein IDs", ";")
tmp2 <- db$"Protein ID"
tstFllID <- "Full ID" %in% colnames(db)
exports <- list("tmp", "tmp2", "tstFllID")
if (tstFllID) {
  tmp3 <- db$"Full ID"
  exports <- append(exports, "tmp3")
}
source(parSrc, local = FALSE)
clusterExport(parClust, exports, envir = environment())
for (i in c("No Isoforms", "Names", "Genes")) { #i <- "No Isoforms"
  j <- if (i == "No Isoforms") { i } else { sub("s$", "", i) }
  if (!j %in% colnames(db)) {
    j <- paste0(sub("s$", "", i), c(" ID", " IDs", ""))
    w <- which(j %in% colnames(db))
    if (length(w)) { j <- j[w[1L]] } else {
      warning(paste0("No near matching column name found for \"", i, "\" in the protein data base table."))
    }
  }
  if (length(j) == 1L) {
    tmp4 <- db[[j]]
    exports <- list("i", "tmp4")
    clusterExport(parClust, exports, envir = environment())
    PG[[i]] <- parSapply(parClust, tmp, \(x) {
      x <- unlist(x)
      m1 <- match(x, tmp2)
      if (tstFllID) {
        m1 <- data.frame(m1 = m1, m2 = match(x, tmp3))
        m1 <- apply(m1, 1L, \(y) {
          y <- unique(y[which(!is.na(y))])
          if (!length(y)) { y <- "" }
          return(y)
        })
        m1 <- unlist(m1)
      }
      x <- tmp4[m1]
      x[which(x %in% c("", " ", "NA", NA))] <- ""
      if (!length(x)) { x <- "" }
      if (i == "Genes") { x <- unique(x) }
      x <- x[which(x != "")]
      x <- paste(x, collapse = ";")
      return(x)
    })
  }
}
# Simplify Gene columns
genkol <- c("Genes", "Gene names")
w <- which(genkol %in% colnames(PG))
if (length(w) == 2L) {
  temp <- PG[, genkol]
  for (i in genkol) { temp[[i]] <- strsplit(temp[[i]], ";") }
  PG$Genes <- apply(temp, 1L, \(x) { paste(sort(unique(unlist(x))), collapse = ";") })
  PG$"Gene names" <- NULL
} else { if (length(w) == 1L) { colnames(PG)[which(colnames(PG) %in% genkol)] <- "Genes" } }
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
#
# If Arabidopsis:
if (("ARATH" %in% db$Organism)||(3702L %in% db$TaxID)) {
  fl <- system.file("extdata", "Uniprot2AGI.txt", package = "proteoCraft")
  tmp <- read.delim(fl, header = FALSE)
  colnames(tmp) <- c("UniProt", "TAIR")
  klnm <- "TAIR"
  if (klnm %in% colnames(db)) { klnm <- "TAIR_v2" }
  db[[klnm]] <- ""
  w <- which(db$`Protein ID` %in% tmp$UniProt)
  db[w, klnm] <- gsub("\\.[0-9]+$", "",
                      gsub("\\.[0-9]+;", ";", toupper(tmp$TAIR[match(db$`Protein ID`[w], tmp$UniProt)])))
  kol <- c("TAIR", "TAIR_v2")
  if (sum(kol %in% colnames(db)) == 2L) {
    #View(db[, kol])
    w <- which((nchar(db$TAIR) > 0L)&(nchar(db$TAIR_v2) > 0L))
    tmp4 <- db[w, kol]
    tmp4$TAIR <- strsplit(tmp4$TAIR, ";")
    tmp4$TAIR_v2 <- strsplit(tmp4$TAIR_v2, ";")
    tmp4 <- apply(tmp4, 1L, unique)
    tmp4 <- parSapply(parClust, tmp4, \(x) { paste(unlist(x), collapse = ";") })
    db$TAIR[w] <- tmp4
    w <- which(db$TAIR == "")
    db$TAIR[w] <- db$TAIR_v2[w]
    db$TAIR_v2 <- NULL
    db$TAIR[which(is.na(db$TAIR))] <- ""
  }
  tmp2 <- listMelt(strsplit(PG$`Leading protein IDs`, ";"), PG$id)
  tmp2$TAIR <- db$TAIR[match(tmp2$value, db$`Protein ID`)]
  tmp3 <- listMelt(strsplit(tmp2$TAIR, ";"), tmp2$L1)
  tmp3 <- aggregate(tmp3$value, list(tmp3$L1), \(x) { paste(unique(x), collapse = ";") })
  PG$TAIR <- ""
  w <- which(PG$id %in% tmp3$Group.1)
  PG$TAIR[w] <- tmp3$x[match(PG$id[w], tmp3$Group.1)]
}

# Here, if this is a BioID type experiment, we also want to mark protein groups which have Biotin peptides:
IsBioID2 %<o% FALSE
if (IsBioID) {
  wbiot %<o% grep("biot", Modifs$"Full name", ignore.case = TRUE)
  l <- length(wbiot)
  if (l) {
    tmp <- if (l == 1L) { Modifs$"Full name"[wbiot] } else {
      paste0(paste(Modifs$"Full name"[wbiot[1L:(l-1L)]], collapse = "\", \""), "\" and \"", Modifs$"Full name"[wbiot[l]])
    }
    warning(paste0("Modifications \"", tmp, "\" were detected as biotinylations, check that this is correct!"))
    g <- grep(topattern(Modifs$Mark[wbiot], start = FALSE), pep$"Modified sequence")
    if (length(g)) {
      wpg <- unique(as.integer(unlist(strsplit(pep$"Protein group IDs"[g], ";"))))
      wpg <- which(PG$id %in% wpg)
      PG[["Biot. peptide IDs"]] <- ""
      temp <- setNames(lapply(strsplit(PG$"Peptide IDs", ";"), as.integer), PG$id)
      temp <- listMelt(temp)
      temp <- temp[which(temp$value %in% pep$id[g]),]
      temp <- aggregate(temp$value, list(temp$L1), \(x) { paste(sort(x), collapse = ";") })
      PG[wpg, "Biot. peptide IDs"] <- temp$x[match(PG$id[wpg], temp$Group.1)]
      PG[["Biot. peptides count"]] <- lengths(strsplit(PG[["Biot. peptide IDs"]], ";"))
      PG[["Biot. peptides [%]"]] <- round(100*PG[["Biot. peptides count"]]/PG$"Peptides count", 1L)
      IsBioID2 <- TRUE
    } else { warning("I could not find any biotinylated peptides!") }
  } else { warning("I could not identify any biotinylated PTMs in the modifications table, did you include them in the search?") }
}
# First sequence
m <- match(vapply(strsplit(PG$"Leading protein IDs", ";"), \(x) { x[[1L]] }, ""),
           db$"Protein ID")
PG$"Protein ID (1st accession)" <- db$`Protein ID`[m]
PG$"Sequence (1st accession)" <- db$Sequence[m]

# Number of spectra, PSMs and peptides per sample:
source(parSrc, local = FALSE)
invisible(clusterCall(parClust, \() {
  library(proteoCraft)
  library(reshape)
  library(data.table)
  return()
}))
temp_PG <- data.frame(id = PG$id,
                      Accession1 = PG$`Protein ID (1st accession)`,
                      Pep = PG$"Peptide IDs")
tmp <- listMelt(strsplit(temp_PG$Pep, ";"), PG$id, c("id", "PG"))
tmp$id <- as.integer(tmp$id)
tmp$Seq <- pep$Sequence[match(tmp$id, pep$id)]
tmp <- as.data.table(tmp)
tmp <- tmp[, .(Seq = list(Seq)), by = .(PG = PG)]
temp_PG$Pep <- tmp$Seq[match(PG$id, tmp$PG)]
temp_PG$Seq <- db$Sequence[match(temp_PG$Accession1, db$"Protein ID")]
exports <-
  if (!"Sequence coverage [%]" %in% colnames(PG)) {
    exports <- list("Coverage")
    clusterExport(parClust, exports, envir = environment())
    PG$"Sequence coverage [%]" <- round(100*parApply(parClust, temp_PG[, c("Seq", "Pep")], 1L, \(x) {
      Coverage(x[[1L]], x[[2L]])
    }), 1L)
  }
CreateMSMSKol %<o% (("MS/MS IDs" %in% colnames(ev))&&(class(ev$"MS/MS IDs") %in% c("integer", "character")))
if (CreateMSMSKol) {
  # There appear to be no MSMS IDs for DIA in MaxQuant.
  ev$temp <- parLapply(parClust, strsplit(as.character(ev$"MS/MS IDs"), ";"), as.integer)
  #PG[, paste0("Spectr", c("al count", "um IDs"))]
  temp <- listMelt(lapply(strsplit(PG$`Evidence IDs`, ";"), as.integer), PG$id, c("Ev_id", "PG_id"))
  temp$MSMSIDs <- ev$temp[match(temp$Ev_id, ev$id)]
  temp <- temp[which(lengths(temp$MSMSIDs) > 0L),] # Remove Match-Between-Runs evidences (no MS/MS)
  temp <- listMelt(temp$MSMSIDs, temp$PG_id, c("MSMSIDs", "PG_id"))
  temp <- do.call(data.frame, aggregate(temp$MSMSIDs, list(temp$PG_id), \(x) {
    x <- unique(x)
    return(c(Count = length(x), List = list(x)))
  }))
  temp$x.Count <- unlist(temp$x.Count)
  temp$Pasted <- vapply(temp$x.List, \(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") }, "")
  PG[, paste0("Spectr", c("al count", "um IDs"))] <- temp[match(PG$id, temp$Group.1), c("x.Count", "Pasted")]
  ev$temp <- NULL
}
temp_ev <- ev[, c("id", "MQ.Exp", "Protein group IDs", "Peptide ID")]
if (CreateMSMSKol) { temp_ev$"MS/MS IDs" <- ev$"MS/MS IDs" }
temp_pep <- pep[, c("id", "Sequence")]
clusterExport(parClust, exports, envir = environment())
Samplez <- list(Sample = c(), Group = c())
MQ_Exp <- list()
for (gr in VPAL$values) { #gr <- VPAL$values[1L]
  wh <- which(Exp.map[[VPAL$column]] == gr)
  smplz <- Exp.map$Ref.Sample.Aggregate[wh]
  mqexp <- setNames(lapply(smplz, \(x) { Exp.map$MQ.Exp[match(x, Exp.map$Ref.Sample.Aggregate)] }), smplz)
  mqexp[[gr]] <- Exp.map$MQ.Exp[wh]
  Samplez$Sample <- c(Samplez$Sample, smplz, gr)
  Samplez$Group <- c(Samplez$Group, rep(gr, length(smplz)), gr)
  MQ_Exp[[gr]] <- mqexp
}
Samplez <- data.frame(Sample = Samplez$Sample,
                      Group = Samplez$Group)
tmpFl1 <- tempfile(fileext = ".rds")
tmpFl2 <- tempfile(fileext = ".rds")
tmpFl3 <- tempfile(fileext = ".rds")
readr::write_rds(temp_ev, tmpFl1)
readr::write_rds(temp_pep, tmpFl2)
readr::write_rds(temp_PG, tmpFl3)
exports <- list("Exp.map", "tmpFl1", "tmpFl2", "tmpFl3", "IsBioID2", "MQ_Exp", "Samplez", "Modifs", "CreateMSMSKol")
if (IsBioID2) { exports <- append(exports, "wbiot") }
clusterExport(parClust, exports, envir = environment())
invisible(clusterCall(parClust, \() {
  temp_ev <- readr::read_rds(tmpFl1)
  temp_pep <- readr::read_rds(tmpFl2)
  temp_PG <- readr::read_rds(tmpFl3)
  assign("temp_ev", temp_ev, envir = .GlobalEnv)
  assign("temp_pep", temp_pep, envir = .GlobalEnv)
  assign("temp_PG", temp_PG, envir = .GlobalEnv)
}))
unlink(tmpFl1)
unlink(tmpFl2)
unlink(tmpFl3)
temp <- parApply(parClust, Samplez, 1L, \(Smpl) { #Smpl <- unlist(Samplez[1,])
  smpl <- Smpl[[1L]]
  gr <- Smpl[[2L]]
  res <- temp_PG[, "id", drop = FALSE]
  kol <- c()
  if (CreateMSMSKol) {
    kols <- paste0("Spectr", c("al count", "um IDs"), " - ", smpl)
    kol <- c(kol, kols)
  }
  kole <- paste0("Evidence", c("s count", " IDs"), " - ", smpl)
  kolp <- paste0("Peptide", c("s count", " IDs"), " - ", smpl)
  kol <- c(kol, kole, kolp)
  kolk <- grep(" count - ", kol, value = TRUE)
  koli <- grep(" IDs - ", kol, value = TRUE)
  res[, kolk] <- 0L
  res[, koli] <- ""
  res[[paste0("Sequence coverage [%] - ", smpl)]] <- 0
  mqexp <- MQ_Exp[[gr]]
  w <- which(temp_ev$MQ.Exp %in% unique(unlist(mqexp[[smpl]])))
  if (length(w)) {
    e <- temp_ev[w, , drop = FALSE]
    temp1 <- lapply(strsplit(e$"Protein group IDs", ";"), as.integer)
    temp1 <- listMelt(temp1, e$"Peptide ID")
    temp1 <- do.call(data.frame, aggregate(temp1$L1, list(temp1$value), \(x) {
      x <- unique(x)
      return(c(Count = length(x), List = list(x)))
    }))
    temp1$x.Count <- unlist(temp1$x.Count)
    temp1$Pasted <- vapply(temp1$x.List, \(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") }, "")
    tmp <- temp1$x.List
    tmp <- listMelt(tmp)
    tmp$Seq <- temp_pep$Sequence[match(as.numeric(tmp$value), temp_pep$id)]
    tmp <- data.table(Seq = tmp$Seq, Row = tmp$L1)
    tmp <- tmp[, list(Seq = list(unique(Seq))), by = Row]
    tmp <- as.data.frame(tmp)
    temp1$Pepseq <- tmp$Seq[match(tmp$Row, 1L:nrow(temp1))]
    temp2 <- lapply(strsplit(e$"Protein group IDs", ";"), as.integer)
    temp2 <- listMelt(temp2, e$id)
    temp2 <- do.call(data.frame, aggregate(temp2$L1, list(temp2$value), \(x) {
      c(Count = length(x), IDs = paste(sort(as.numeric(x)), collapse = ";"))
    }))
    temp2$x.Count <- as.integer(temp2$x.Count)
    if (CreateMSMSKol) {
      w3 <- which(e$"MS/MS IDs" != "")
      temp3 <- lapply(strsplit(e$"Protein group IDs"[w3], ";"), as.integer)
      temp3 <- listMelt(temp3, e$"MS/MS IDs"[w3])
      temp3$L1 <- lapply(strsplit(temp3$L1, ";"), as.integer)
      temp3 <- listMelt(temp3$L1, temp3$value)
      temp3 <- do.call(data.frame, aggregate(temp3$value, list(temp3$L1), \(x) {
        x <- unique(x)
        return(c(Count = length(x), List = list(x)))
      }))
      temp3$x.Count <- unlist(temp3$x.Count)
      temp3$Pasted <- vapply(temp3$x.List, \(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") }, "")
      w <- which(res$id %in% temp3$Group.1)
      m <- match(res$id[w], temp3$Group.1)
      res[w, kols] <- temp3[m, c("x.Count", "Pasted")]
    }
    w <- which(res$id %in% temp1$Group.1)
    m <- match(res$id[w], temp1$Group.1)
    res[w, kolp] <- temp1[m, c("x.Count", "Pasted")]
    temp_PG$Pep <- NA
    temp_PG$Pep[w] <- temp1$Pepseq[m]
    res[w, paste0("Sequence coverage [%] - ", smpl)] <- round(100*apply(temp_PG[w, c("Seq", "Pep")], 1L, \(x) {
      Coverage(x[[1L]], x[[2L]])
    }), 1L)
    w <- which(res$id %in% temp2$Group.1)
    m <- match(res$id[w], temp2$Group.1)
    res[w, kole] <- temp2[m, c("x.Count", "x.IDs")]
    if (IsBioID2) {
      kolB <- paste0("Biot. ", paste0(tolower(substr(kol, 1L, 1L)), substr(kol, 2L, nchar(kol))))
      kolBe <- grep("^Biot\\. evidence", kolB, value = TRUE)
      kolBp <- grep("^Biot\\. peptide", kolB, value = TRUE)
      if (CreateMSMSKol) {
        kolBs <- grep("^Biot\\. spectr", kolB, value = TRUE)
        kolB <- kolB[which(!kolB %in% kolBs)]; rm(kolBs) #I don't think we need those columns now... too many is too many
      }
      kolBk <- grep(" count - ", kolB, value = TRUE)
      kolBi <- grep(" IDs - ", kolB, value = TRUE)
      res[, kolBk] <- 0L
      res[, kolBi] <- ""
      g <- grep(topattern(Modifs$Mark[wbiot], start = FALSE), e$"Modified sequence")
      if (length(g)) {
        eB <- e[g, , drop = FALSE]
        temp1 <- listMelt(lapply(strsplit(eB$"Protein group IDs", ";"), as.integer), eB$"Peptide ID")
        temp1 <- do.call(data.frame, aggregate(temp1$L1, list(temp1$value), \(x) {
          x <- unique(x)
          return(c(Count = length(x), IDs = paste(sort(as.numeric(x)), collapse = ";")))
        }))
        temp1$x.Count <- as.integer(temp1$x.Count)
        temp1 <- do.call(data.frame, temp1)
        temp2 <- strsplit(eB$"Protein group IDs", ";")
        temp2 <- listMelt(temp2, eB$id)
        temp2 <- do.call(data.frame, aggregate(temp2$L1, list(temp2$value), \(x) {
          c(Count = length(x), IDs = paste(sort(x), collapse = ";"))
        }))
        temp2$x.Count <- as.integer(temp2$x.Count)
        w <- which(res$id %in% temp1$Group.1)
        m <- match(res$id[w], temp1$Group.1)
        res[w, kolBp] <- temp1[m, c("x.Count", "x.IDs")]
      }
    }
  }
  res$id <- NULL
  return(res)
})
temp <- do.call(cbind, temp)
PG[, colnames(temp)] <- temp
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
#View(PG[, grep("^Spectr|^Peptide|^Evidence", colnames(PG))])

# CRAPome
if (CRAPome) {
  tst <- strsplit(PG$`Protein IDs`, ";")
  tst <- listMelt(tst)
  tst <- tst[which(tst$value %in% CRAPomeProteins),]
  PG$`Potential contaminant`[which(PG$id %in% tst$L1)] <- "+"
}

# Proteins in list
if (prot.list.Cond) {
  PG$"In list" <- ""
  g1 <- grsep2(prot.list, PG$"Protein IDs")
  PG$`In list`[g1] <- "+"
  PG$"Potential contaminant"[g1] <- ""
  pep$"In list" <- ""
  g2 <- grsep2(prot.list, pep$"Proteins")
  pep$`In list`[g2] <- "+"
  pep$"Potential contaminant"[g2] <- ""
  ev$"Potential contaminant"[grsep2(prot.list, ev$Proteins)] <- ""
  pep$"Potential contaminant"[grsep2(prot.list, pep$Proteins)] <- ""
  db$"Potential contaminant"[which(db$`Protein ID` %in% prot.list)] <- ""
}

# Backup data/update cluster
stopClust <- TRUE
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Calculate protein group-level quantitative values
post_ReNorm_reRun <- FALSE
quntSrc %<o% paste0(libPath, "/extdata/Sources/PG_Quant.R")
#rstudioapi::documentOpen(quntSrc)
source(quntSrc, local = FALSE)

#### Code chunk - Re-normalize protein group expression values
Src <- paste0(libPath, "/extdata/Sources/PG_ReNorm.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Optional "mock" back correction
# -------------------------------
# In cases where we do have a batch, but did not correct for it,
# we would still do a batch correction here so we would be able to draw PCA plots without the batch effect.
# The idea is that we can see a good approximation of the data's structure if the batch is removed.


# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

# Test expression values:
quantData <- quantData_list$Data
g <- grep(topattern(Prot.Expr.Root), colnames(quantData), value = TRUE)
g <- grep(": SD$", g, value = TRUE, invert = TRUE)
g <- grep("\\.REF$", g, value = TRUE, invert = TRUE)
test <- quantData[, g]
colnames(test) <- gsub(topattern(Prot.Expr.Root), "", colnames(test))
test <- test[which(apply(test, 1L, \(x) { length(is.all.good(x)) }) > 0L),]
test <- suppressMessages(dfMelt(test))
test$variable <- as.character(test$variable)
test[, RSA$names] <- ""
w <- rep(FALSE, nrow(test))
test[which(!w), RSA$names] <- Isapply(strsplit(test$variable[which(!w)], "___"), unlist)
a <- RSA$names
w <- which(vapply(a, \(x) { length(unique(test[[x]])) }, 1L) > 1L)
if (length(w)) { a <- a[w] }
test[[a[1L]]] <- factor(test[[a[1L]]], levels = sort(unique(test[[a[1L]]])))
test <- test[which(is.all.good(test$value, 2L)),]
test2 <- set_colnames(aggregate(test$value, list(test$variable), median), c("variable", "value"))
test2[, a] <- test[match(test2$variable, test$variable), a]
MinMax <- c(min(test$value), max(test$value))
nbinz <- ceiling((MinMax[2L]-MinMax[1L])/0.1)
binz <- c(0L:nbinz)/nbinz
binz <- binz*(MinMax[2L]-MinMax[1L])+MinMax[1L]
binz[1L] <- binz[1L]-0.000001
testI <- data.frame(Intensity = (binz[2L:(nbinz+1L)]+binz[1L:nbinz])/2)
for (v in unique(test$variable)) {
  wv <- which(test$variable == v)
  testI[[v]] <- vapply(1L:nbinz, \(x) {
    sum((test$value[wv] > binz[x])&(test$value[wv] <= binz[x+1L]))
  }, 1L)
}
testI <- reshape2::melt(testI, id.vars = "Intensity")
testI[, a] <- test[match(testI$variable, test$variable), a]
testI$variable <- cleanNms(testI$variable)
dir <- paste0(wd, "/Workflow control/Protein groups/Expression")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
ttl <- "Protein groups - distribution of Expression values"
plot <- ggplot(testI) +
  geom_area(aes(x = Intensity, y = value, fill = variable, group = variable,
                colour = variable), alpha = 0.25) +
  scale_color_viridis_d(begin = 0.25) +
  scale_fill_viridis_d(begin = 0.25) +
  geom_vline(data = test2, aes(xintercept = value), linetype = "dashed", color = "grey") +
  ggtitle(ttl) + theme_bw() + theme(legend.position = "none", strip.text.y = element_text(angle = 0)) +
  scale_y_continuous(limits = c(0, max(testI$value)*1.1), expand = c(0L, 0L))
plot <- if (length(a) == 1L) { plot + facet_wrap(as.formula(paste0("~", a))) } else {
  plot + facet_grid(as.formula(paste0(a[1L], "~", paste(a[2L:length(a)], collapse = "+"))))
}
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150L, width = 10L, height = 10L, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150L, width = 10L, height = 10L, units = "in")
})
ReportCalls <- AddPlot2Report()
#
# Test ratio values:
g <- grep(topattern(Prot.Rat.Root), colnames(quantData), value = TRUE)
g <- grep(": SD$", g, value = TRUE, invert = TRUE)
g <- grep("REF\\.to\\.REF", g, value = TRUE, invert = TRUE)
test <- quantData[, g]
colnames(test) <- gsub(topattern(Prot.Rat.Root), "", colnames(test))
test <- test[which(apply(test, 1L, \(x) { length(is.all.good(x)) }) > 0L),]
test <- suppressMessages(dfMelt(test))
test$Contrast <- gsub_Rep(" - ", " -\n", as.character(test$variable))
allContr <- unique(test$Contrast)
test$Contrast <- factor(test$Contrast, levels = allContr)
test <- test[which(is.all.good(test$value, 2L)),]
test2 <- set_colnames(aggregate(test$value, list(test$Contrast), median), c("Contrast", "value"))
MinMax <- c(min(test$value), max(test$value))
nbinz <- ceiling((MinMax[2L]-MinMax[1L])/0.1)
binz <- c(0L:nbinz)/nbinz
binz <- binz*(MinMax[2L]-MinMax[1L])+MinMax[1L]
binz[1L] <- binz[1L]-0.000001
testR <- data.frame(Intensity = (binz[2L:(nbinz+1L)]+binz[1L:nbinz])/2)
for (ctr in allContr) {
  wv <- which(test$Contrast == ctr)
  testR[[ctr]] <- vapply(1L:nbinz, \(x) {
    sum((test$value[wv] > binz[x])&(test$value[wv] <= binz[x+1L]))
  }, 1L)
}
testR <- dfMelt(testR, id.vars = "Intensity")
colnames(testR)[which(colnames(testR) == "variable")] <- "Contrast"
dir <- paste0(wd, "/Workflow control/Protein groups/Expression")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
ttl <- "Protein groups - distribution of Ratios"
plot <- ggplot(testR) +
  geom_area(aes(x = Intensity, y = value, fill = Contrast, group = Contrast,
                colour = Contrast), alpha = 0.25) +
  scale_color_viridis(begin = 0.25, discrete = TRUE, option = "D") +
  scale_fill_viridis(begin = 0.25, discrete = TRUE, option = "D") +
  geom_vline(data = test2, aes(xintercept = value), linetype = "dashed", color = "grey") +
  ggtitle(ttl) + theme_bw() + theme(legend.position = "none", strip.text.y = element_text(angle = 0)) +
  scale_y_continuous(limits = c(0, max(testR$value)*1.1), expand = c(0L, 0L)) + facet_wrap(~Contrast)
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150L, width = 10L, height = 10L, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150L, width = 10L, height = 10L, units = "in")
})
ReportCalls <- AddPlot2Report()

#### SubCellular Localisation analysis: re-scale to known proportions
Src <- paste0(libPath, "/extdata/Sources/SubCellResc.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Code chunk - Add quant data to PG:
PG <- PG[, which(!colnames(PG) %in% colnames(quantData))]
PG[, colnames(quantData)] <- quantData
if (Param$Prot.Only.with.Quant) {
  colnames(quantData)
  test1 <- apply(quantData[,grep(topattern(Prot.Expr.Root), colnames(quantData), value = TRUE)],
                 1L, \(x) { length(is.all.good(x)) })
  a <- grep(topattern(Prot.Rat.Root), colnames(quantData), value = TRUE)
  a <- a[which(!grepl(": SD$|: -log10\\(peptides Pvalue\\)$", a))]
  test2 <- apply(quantData[,a],
                 1L, \(x) { length(is.all.good(x)) })
  PG <- PG[which((test1 > 0L)|(test2 > 0L)),]
}
if (!"Peptides count" %in% colnames(PG)) {
  PG$"Peptides count" <- lengths(strsplit(PG$"Peptide IDs", ";"))
}

if (!Param$Plot.labels %in% colnames(PG)) {
  tmp <- gsub("\\.", " ", Param$Plot.labels)
  if (tmp %in% colnames(PG)) {
    warning(paste0("Protein groups table column \"", Param$Plot.labels, "\" not found, the column is called \"", tmp, "\" (check parameter \"Plot.labels\")"))
    Param$Plot.labels <- tmp
  } else {
    tmp <- c("Common Name (short)", "Common Names", "Names", "Protein IDs", "Common.Names.short", "Common.Names", "Protein.IDs")
    w <- which(tmp %in% colnames(PG))
    warning(paste0("Protein groups table column \"", Param$Plot.labels, "\" not found (check parameter \"Plot.labels\"), defaulting to \"", tmp[w[1L]], "\""))
    Param$Plot.labels <- tmp[w[1L]]
  }
}

#### Code chunk - Optional - Apply True-Discovery or negative filter
Src <- paste0(libPath, "/extdata/Sources/PG_Filters.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - samples Pearson correlation heatmap
Src <- paste0(libPath, "/extdata/Sources/pearsonCorrMap.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

# Average expression columns
clusterExport(parClust, "is.all.good", envir = environment())
for (grp in VPAL$values) { #grp <- VPAL$values[1L] #grp <- VPAL$values[3L]
  em <- Exp.map[which(Exp.map[[VPAL$column]] == grp),]
  # PG
  kol <- paste0(Prot.Expr.Root, em$Ref.Sample.Aggregate)
  kol <- intersect(kol, colnames(PG))
  if (length(kol)) {
    tmp <- as.data.frame(t(parApply(parClust, PG[, kol, drop = FALSE], 1L, Av_SE_fun)))
    PG[[paste0("Mean ", Prot.Expr.Root, grp)]]  <- tmp[, 1L]
  }
  # Pep
  kol <- paste0(pep.ref[length(pep.ref)], em$Ref.Sample.Aggregate)
  kol <- intersect(kol, colnames(pep))
  if (length(kol)) {
    tmp <- as.data.frame(t(parApply(parClust, log10(pep[, kol, drop = FALSE]), 1L, Av_SE_fun)))
    pep[[paste0("Mean ", pep.ref[length(pep.ref)], grp)]]  <- 10^tmp[, 1L]
  }
}

#### Code chunk - Perform statistical tests
dataType <- "PG"
Src <- paste0(libPath, "/extdata/Sources/Stat_tests.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

### Visualize and check P-values
Src <- paste0(libPath, "/extdata/Sources/pVal_check.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

useSAM %<o% ((names(pvalue.col)[which(pvalue.use)] == "Student")&&(useSAM_thresh))
if (useSAM) {
  # In this case, we bypass the original decision and base it off SAM even though we plot Student's P-values
  for (i in names(SAM_thresh)) { #i <- names(SAM_thresh)[1L]
    dec <- SAM_thresh[[i]]$decision
    mKol <- rev(colnames(dec))[1L]
    FCkol <- paste0("Mean ", Prot.Rat.Root, i)
    stopifnot(FCkol %in% names(PG))
    regKol <- paste0("Regulated - ", i)
    PG[[regKol]] <- "non significant"
    fdrs <- as.numeric(gsub("FDR$", "", colnames(dec)[which(colnames(dec) != mKol)]))
    fdrs <- sort(fdrs, decreasing = TRUE)
    for (f in fdrs) { #f <- fdrs[1L]
      w <- which(PG[[mKol]] %in% dec[which(dec[[paste0(f, "FDR")]] == "+"), mKol])
      if (length(w)) {
        PG[which(PG[w, FCkol] > 0), regKol] <- paste0("up, FDR = ", f*100, "%")
        PG[which(PG[w, FCkol] < 0), regKol] <- paste0("down, FDR = ", f*100, "%")
      }
    }
  }
}
#
rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)

#### Code chunk - Prepare Annotations and (if applicable) GO terms
# NB: I used to get functional annotations for all proteins in the protein group.
#     However we are now - and I think with reason - only using annotations from the leading protein(s)!
setwd(wd)
Src <- paste0(libPath, "/extdata/Sources/Annotate_me.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Get or generate subcellular localisation markers - for later use
Src <- paste0(libPath, "/extdata/Sources/SubCellMark.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)

#### Code chunk - ROC analysis
Src <- paste0(libPath, "/extdata/Sources/ROC2.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Estimate P-value significance for a set of accepted FDRs
## NB: For graphical reasons (volcano plots), there is only support for 4 different FDR values. This should suffice anyway.
a <- sapply(strsplit(Param$Plot.metrics, ";"), \(x) { strsplit(x, ":") })
a[[2L]][2L] <- gsub("\\.$", "", pvalue.col[which(pvalue.use)])
Param$Plot.metrics <- paste(vapply(a, paste, "", collapse = ":"), collapse = ";")
FDR.thresholds %<o% c()

A <- myContrasts$Contrast
test <- vapply(A, \(x) { #x <- A[6]
  x <- paste0(pvalue.col[which(pvalue.use)], x)
  r <- x %in% colnames(PG)
  if (r) { r <- length(is.all.good(as.numeric(PG[[x]]))) > 0L }
  return(r)
}, TRUE)
A <- A[which(test)]
PG <- PG[, grep("^Significant-FDR=", colnames(PG), invert = TRUE)]
for (a in A) { #a <- A[1L]
  temp <- FDR(data = PG,
              aggregate = a,
              pvalue_root = pvalue.col[which(pvalue.use)],
              fdr = BH.FDR, returns = rep(TRUE, 3L), method = "BH")
  PG[, colnames(temp$`Significance vector`)] <- temp$`Significance vector`
  FDR.thresholds <- c(FDR.thresholds, temp$Thresholds)
}
#View(PG[, grep("^Significant-FDR=", colnames(PG))])

#### Code chunk - Optional: adjust P-values
Adj_Pval %<o% (("Adjust.P.values" %in% colnames(Param))&&(Param$Adjust.P.values))
if (Adj_Pval) {
  pval.adjust %<o% TRUE
  w <- grep("^adj. ", names(pvalue.col), invert = TRUE)
  pvalue.col <- pvalue.col[w]
  pvalue.use <- pvalue.use[w]
  pkol <- grep(topattern(pvalue.col), colnames(PG), value = TRUE)
  if (!length(pkol)) { stop("There should be P-value columns in the protein groups table at this stage!") }
  msg <- "Adjusting P-values..."
  ReportCalls <- AddMsg2Report(Space = FALSE, Print = FALSE)
  for (pk in pkol) { #pk <- pkol[1L]
    pk2 <- gsub("-log10\\(Pvalue\\) ", "-log10(adj. Pvalue) ", pk)
    if (pk2 == pk) { stop("Bug!!!") } else {
      PG[[pk2]] <- -log10(p.adjust(10^(-PG[[pk]]), method = "BH"))
    }
  }
  pvalue.col2 <- gsub("-log10\\(Pvalue\\)\\.", "-log10(adj. Pvalue).", pvalue.col)
  names(pvalue.col2) <- paste0("adj. ", names(pvalue.col2))
  pvalue.col <- c(pvalue.col, pvalue.col2)
  pvalue.use <- c(pvalue.use, rep(FALSE, length(pvalue.use))) # For now we are providing adjusted P-value columns but still plotting raw P-values on volcano plots.
  #adj.pval.thresh %<o% data.frame(yintercept = -log10(c(0.01, 0.05)),
  #                              slope = rep(0, 2L),
  #                              xintercept = rep(NA, 2L),
  #                              colour = colorRampPalette(c("orange", "red"))(2L),
  #                              label = paste0(c(0.01, 0.05)*100, "% P-value"))
  l <- length(DatAnalysisTxt)
  DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                              " P-values were adjusted using the Benjamini-Hochberg (FDR) method.")
}

# Create list of control ratio values for the purpose of identifying vertical thresholds for plots:
Src <- paste0(libPath, "/extdata/Sources/ratThresh.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Heatmaps with clustering at samples and protein groups level, highlighting proteins of interest
clustMode <- "standard"
Src <- paste0(libPath, "/extdata/Sources/cluster_Heatmap_Main.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - t-test volcano plot(s)
#
PG$"1-PEP" <- 1 - PG$PEP
PG$"log10(1-PEP)" <- log10(PG$"1-PEP")
PG$"log10(Peptides count)" <- log10(PG$"Peptides count")
a <- grep(topattern(Prot.Expr.Root), colnames(PG), value = TRUE)
a <- grep("\\.REF$", a, value = TRUE, invert = TRUE)
PG$"Av. log10 abundance" <- apply(PG[, a], 1L, \(x) { mean(is.all.good(unlist(x))) })
PG$"Rel. av. log10 abundance" <- PG$"Av. log10 abundance"/max(is.all.good(PG$"Av. log10 abundance"))
PG$"Rel. log10(Peptides count)" <- PG$"log10(Peptides count)"/max(is.all.good(PG$"log10(Peptides count)"))
# Plotly
create_plotly %<o% TRUE
create_plotly_local %<o% TRUE # No need for a licence when I can save local htmls! Still, old legacy code kept below.
# create_plotly <- !((as.character(Param$Plotly_user_name) %in% c("", "NA", " "))&(as.character(Param$Plotly_API_key) %in% c("", "NA", " ")))
# if ((create_plotly)&&(!create_plotly_local)) {
#   plotly_subfolder %<o% gsub(":|\\*|\\?|<|>|\\|", "-", Param$Project)
#   plotly_subfolder <- paste0(gsub("/+$", "", plotly_subfolder), "/")
#   Sys.setenv("plotly_username" = Param$Plotly_user_name)
#   Sys.setenv("plotly_api_key" = Param$Plotly_API_key)
#   plot_ly %<o% list()
#   plot_ly$"t-tests" <- list()
# }
# Arbitrary thresholds
arbitrary.thr %<o% data.frame(yintercept = -log10(c(0.05, 0.01)),
                              slope = c(0, 0),
                              xintercept = c(NA, NA),
                              colour = c("orange", "red"),
                              label = c("5% P-value", "1% P-value"))
volcano.plots %<o% list()
# For now, we are plotting non adjusted P-values because in bad cases adjusting causes all to collapse on 1 (0 as -log10)
# Instead, we are sticking to plotting raw P-values with FDR thresholds
#
PrLabKol %<o% setNames(c("Common Name (short)", "Protein IDs", "Genes", "PEP"),
                       c("Protein name", "Protein ID(s)", "Gene(s)", "PEP"))
subDr <- "Reg. analysis/t-tests"
setwd(wd)
# Default volcano plot arguments
# (I've given this poor mess of an ever-evolving function so many arguments over the years!)
volcPlot_args %<o% list(mode = "custom",
                        experiments.map = Exp.map,
                        contrasts = myContrasts,
                        X.root = Prot.Rat.Root,
                        Y.root = pvalue.col[which(pvalue.use)],
                        aggregate.map = Aggregate.map,
                        aggregate.name = VPAL$aggregate,
                        aggregate.list = Aggregate.list,
                        parameters = Param,
                        save = c("jpeg", "pdf"),
                        labels = c("FDR", "both")[useSAM+1L],
                        Ref.Ratio.method = paste0("obs", RefRat_Mode),
                        ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                        FDR.thresh = FDR.thresholds,
                        arbitrary.lines = arbitrary.thr,
                        proteins = prot.list,
                        proteins_split = protsplit,
                        return = TRUE,
                        return.plot = TRUE,
                        title = "Volcano plot ",
                        subfolder = subDr,
                        subfolderpertype = FALSE,
                        Alpha = "Rel. log10(Peptides count)",
                        Size = "Rel. av. log10 abundance",
                        Size.max = 2L,
                        plotly = create_plotly,
                        plotly_local = create_plotly_local,
                        plotly_labels = PrLabKol,
                        SAM = useSAM,
                        curved_Thresh = SAM_thresh,
                        saveData = TRUE)
volcPlot_args2 <- volcPlot_args
volcPlot_args2$Prot <- PG
volcPlot_args2$cl <- parClust
tempVP <- try(do.call(Volcano.plot, volcPlot_args2), silent = TRUE)
if ((inherits(tempVP, "try-error"))||(is.character(tempVP))) {
  stop("MAJOR ERROR: No volcano plots were created, investigate!")
}
#
# Save plotly plots
dr <- paste0(wd, "/", subDr)
myPlotLys <- tempVP$"Plotly plots"
Src <- paste0(libPath, "/extdata/Sources/save_Plotlys.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
VP_list <- tempVP
insrt <- ""
Src <- paste0(libPath, "/extdata/Sources/thresholds_Excel.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
thresh <- lapply(names(tempVP$Thresholds$Absolute), \(x) {
  y <- tempVP$Thresholds$Absolute[[x]]
  x <- data.frame(Test = rep(cleanNms(x), nrow(y)))
  return(cbind(x, y))
})
thresh <- plyr::rbind.fill(thresh)
thresh$Name <- NULL
thresh$Root <- gsub(" - $", "", thresh$Root)
thresh$Value <- thresh$Text.value
thresh$Text.value <- NULL
fdrThresh <- tempVP$Thresholds$FDR
fdrThresh$Test <- cleanNms(fdrThresh$Sample)
fdrThresh$Sample <- NULL
colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.up")] <- "Colour (up)"
colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.down")] <- "Colour (down)"
colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.line")] <- "Colour (line)"
fdrThresh <- fdrThresh[, c("Test", colnames(fdrThresh)[which(colnames(fdrThresh) != "Test")])]
fl <- paste0(wd, "/", subDr, "/Thresholds.xlsx")
wb <- wb_workbook()
wb <- wb_set_creators(wb, "Me")
wb <- wb_add_worksheet(wb, "Thresholds")
dms <- wb_dims(2L, 1L)
wb <- wb_add_data(wb, "Thresholds", "Absolute thresholds", dms)
wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                  italic = "true", underline = "single")
dms <- wb_dims(3L, 2L)
wb <- wb_add_data_table(wb, "Thresholds", thresh, dms,
                        col_names = TRUE, table_style = "TableStyleMedium2",
                        banded_rows = TRUE, banded_cols = FALSE)
dms <- wb_dims(3L+nrow(thresh)+3L, 1L)
wb <- wb_add_data(wb, "Thresholds", "FDR thresholds", dms)
wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                  italic = "true", underline = "single")
dms <- wb_dims(3L+nrow(thresh)+4L, 2L)
wb <- wb_add_data_table(wb, "Thresholds", fdrThresh, dms,
                        col_names = TRUE, table_style = "TableStyleMedium2",
                        banded_rows = TRUE, banded_cols = FALSE)
wb <- wb_set_col_widths(wb, "Thresholds", 1L, 3L)
tmp1 <- rbind(colnames(thresh), thresh)
colnames(tmp1) <- paste0("V", 1L:ncol(tmp1))
tmp2 <- rbind(colnames(fdrThresh), fdrThresh)
colnames(tmp2) <- paste0("V", 1L:ncol(tmp2))
tst <- plyr::rbind.fill(tmp1, tmp2)
tst <- setNames(apply(tst, 2L, \(x) { max(nchar(x), na.rm = TRUE) }), NULL)
wb <- wb_set_col_widths(wb, "Thresholds", 1L:(length(tst)+1L), c(3L, tst))
wb_save(wb, fl)
#xl_open(fl)
#
temp3 <- tempVP$Protein_groups_file
temp4 <- tempVP$Plots
g <- grep("Regulated - ", colnames(temp3), value = TRUE)
#View(temp3[,g])
PG[,g] <- temp3[,g]
n2 <- names(temp4$Labelled)
volcano.plots$Unlabelled <- temp4$Unlabelled
volcano.plots$Labelled <- temp4$Labelled
dir <- paste0(wd, "/Reg. analysis/t-tests")
for (ttl in n2) {
  plot <- volcano.plots$Labelled[[ttl]]
  ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
}

# Also calculate Q-values - for now, the plot is created but not saved!
if (("Q.values" %in% colnames(Param))&&(is.logical(Param$Q.values))&&(Param$Q.values)) {
  require(qvalue)
  pkol <- grep(topattern(pvalue.col[which(pvalue.use)]), colnames(PG), value = TRUE)
  if (length(pkol)) {
    msg <- "Calculating Q-values..."
    ReportCalls <- AddMsg2Report(Space = FALSE, Print = FALSE)
    for (pk in pkol) { #pk <- pkol[2L]
      temp <- 10^(-PG[[pk]])
      wag <- which(is.all.good(temp, 2L))
      pi0 <- qvalue::pi0est(temp)$pi0[1L]
      temp <- try(qvalue::qvalue(temp[wag]), silent = TRUE) # For now we do not explicitly set pi0
      if (inherits(temp, "try-error")) {
        temp <- try(qvalue::qvalue(temp[wag], pi0 = pi0), silent = TRUE)
        while ((inherits(temp, "try-error"))&&(pi0 <= 1)) {
          pi0 <- pi0 + 0.05
          temp <- try(qvalue::qvalue(temp[wag], pi0 = pi0), silent = TRUE)
        }
      }
      if (!inherits(temp, "try-error")) {
        PG[[gsub(topattern(pvalue.col[which(pvalue.use)]), "-log10(Qvalue) - ", pk)]] <- NA
        PG[[gsub(topattern(pvalue.col[which(pvalue.use)]), "local FDR ", pk)]] <- NA
        PG[wag, gsub(topattern(pvalue.col[which(pvalue.use)]), "-log10(Qvalue) - ", pk)] <- -log10(temp$qvalues)
        PG[wag, gsub(topattern(pvalue.col[which(pvalue.use)]), "local FDR ", pk)] <- temp$lfdr
      }
    }
    qval.thresh %<o% data.frame(yintercept = -log10(BH.FDR),
                                slope = rep(0, length(BH.FDR)),
                                xintercept = rep(NA, length(BH.FDR)),
                                colour = colorRampPalette(c("orange", "red"))(length(BH.FDR)),
                                label = paste0(BH.FDR*100, "% FDR"))
    # Probably deprecated... check arguments before running
    subDr <- "Reg. analysis/t-tests"
    volcPlot_args2 <- volcPlot_args
    volcPlot_args2$Prot <- PG
    volcPlot_args2$Y.root <- "-log10(Qvalue) - "
    volcPlot_args2$arbitrary.thresh <- qval.thresh
    volcPlot_args2$title <- "Q-values volcano plot "
    volcPlot_args2$cl <- parClust
    tempVP2 <- do.call(Volcano.plot, volcPlot_args2)
    #
    # Save plotly plots
    dr <- paste0(wd, "/", subDr)
    myPlotLys <- tempVP2$"Plotly plots"
    Src <- paste0(libPath, "/extdata/Sources/save_Plotlys.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    VP_list <- tempVP2
    insrt <- "_Qvalues"
    Src <- paste0(libPath, "/extdata/Sources/thresholds_Excel.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    l <- length(DatAnalysisTxt)
    DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                                " Q-values were computed using package qvalue.")
  }
}

# Specificity mark for untested proteins
# We can assume that proteins with PSMs only in the specific pull-down samples are actually specifically enriched!
# -> Label them as such!
# NB: This used to be for pull-downs only but I am extending it to the whole workflow for now.
#     This could also be parameter-controlled.
# NB: For the F-test, this is done within the function
#if (IsPullDown) {
MinPep4Spec %<o% 2L
if ("Min.pep.per.sample.for.spec" %in% colnames(Param)) {
  MinPep4Spec <- Param$Min.pep.per.sample.for.spec
  if (!is.numeric(MinPep4Spec)) {
    warning("Invalid value for parameter \"Min.pep.per.sample.for.spec\", defaulting to 2!")
    MinPep4Spec <- 2L
  }
}
kolTR <- paste0("Regulated - ", myContrasts$Contrast)
w <- which(kolTR %in% colnames(PG))
if (length(w)) {
  kolTR <- kolTR[w]
  for (i in w) { #i <- w[1L]
    koleA <- paste0("Evidences count - ", myContrasts$A_samples[[i]])
    kolpA <- paste0("Peptides count - ", myContrasts$A_samples[[i]])
    koleB <- paste0("Evidences count - ", myContrasts$B_samples[[i]])
    kolpB <- paste0("Peptides count - ", myContrasts$B_samples[[i]])
    tstA <- rowSums(PG[, kolpA] >= MinPep4Spec, na.rm = TRUE) == length(kolpA) # Are there at least MinPep4Spec peptidoforms in A?...
    tstB <- rowSums(PG[, kolpB] == 0L, na.rm = TRUE) == length(kolpB) #... and concurrently none in B?
    wA <- which(tstA & tstB)
    if (length(wA)) {
      pepmin <- apply(PG[wA, kolpA, drop = FALSE], 1L, min)
      evcount <- apply(PG[wA, koleA, drop = FALSE], 1L, sum)
      txtup <- paste0("Specific: at least ", pepmin, " pep./sample (", evcount, " PSMs tot.)")
      PG[wA, kolTR[i]] <- txtup
    }
  }
} else { stop("I would expect \"Regulated ...\" columns in the PG table by this stage!") }

# Create t-test filters:
## These can then be used for further steps down the line, such as volcano plots, etc...
Reg_filters %<o% list()
filter_types %<o% tolower(unlist(strsplit(Param$Filters.type, ";")))
filter_types[grep("^dat.+2$", filter_types, invert = TRUE)] <- substr(filter_types[which(!grepl("^dat.+2$", filter_types))], 1L, 3L)
filter_types[grep("^dat.+2$", filter_types)] <- "dat2"
filter_types <- unique(c("con", filter_types))
if ("ref" %in% filter_types) {
  if (Nested) {
    warning("Grouping filter by reference is not feasible if replicates are paired!")
    filter_types <- filter_types[which(filter_types != "ref")]
  } else {
    if (sum(vapply(RG$names, \(x) { !x %in% RSA$names }, TRUE))) {
      warning("Grouping filter by reference is not feasible if the factors used for \"RG\" are not included in those used for \"RRG\"!")
      filter_types <- filter_types[which(filter_types != "ref")]
    }
  }
}
g <- grep("^Regulated - ", colnames(PG), value = TRUE)
g1 <- gsub("^Regulated - ", "", g)
up <- grep("^up|^Specific", unique(unlist(PG[, g])), value = TRUE)
down <- grep("^down|^Anti-specific", unique(unlist(PG[, g])), value = TRUE) # The "anti-specific" part will only become relevant if in future I disconnect symmetry and/or adding specific tags from pull-down experiments
Reg_filters$"t-tests" <- list()
if ("con" %in% filter_types) {
  Reg_filters$"t-tests"$"By condition" <- list()
  for (i in seq_along(g)) {
    Reg_filters$"t-tests"$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                          Filter_up = sort(which(PG[[g[i]]] %in% up)),
                                                          Filter_down = sort(which(PG[[g[i]]] %in% down)),
                                                          Filter = sort(which(PG[[g[i]]] %in% c(up, down))))
  }
}
if ("ref" %in% filter_types) {
  Reg_filters$"t-tests"$"By reference" <- list()
  g2 <- as.data.frame(t(as.data.frame(strsplit(g1, "___"))))
  colnames(g2) <- VPAL$names
  tst <- do.call(paste, c(Exp.map[, RRG$names, drop = FALSE], sep = "___"))
  tmp <- do.call(paste, c(g2[, RRG$names, drop = FALSE], sep = "___"))
  g2$Ref <- vapply(tmp, \(x) { unique(tst[which((Exp.map$Reference)&(tst == x))]) }, "")
  for (i in unique(g2$Ref)) {
    w <- which(g2$Ref == i)
    u <- grep("^up|^Specific", unique(as.character(PG[, g[w]])), value = TRUE)
    d <- grep("^down", unique(as.character(PG[, g[w]])), value = TRUE)
    Reg_filters$"t-tests"$"By reference"[[i]] <- list(Columns = g[w],
                                                      Filter_up = sort(which(apply(PG[, g[w], drop = FALSE], 1L, \(x) {
                                                        length(which(x %in% up))
                                                      }) > 0L)),
                                                      Filter_down = sort(which(apply(PG[, g[w], drop = FALSE], 1L, \(x) {
                                                        length(which(x %in% down))
                                                      }) > 0L)),
                                                      Filter = sort(which(apply(PG[, g[w], drop = FALSE], 1L, \(x) {
                                                        length(which(x %in% c(up, down)))
                                                      }) > 0L)))
  }
}
if (sum(c("dat", "dat2") %in% filter_types)) {
  Reg_filters$"t-tests"$"Whole dataset" <- list(Columns = g,
                                                Filter_up = sort(which(apply(PG[, g, drop = FALSE], 1L, \(x) {
                                                  length(which(x %in% up))
                                                }) > 0L)),
                                                Filter_down = sort(which(apply(PG[, g, drop = FALSE], 1L, \(x) {
                                                  length(which(x %in% down))
                                                }) > 0L)),
                                                Filter = sort(which(apply(PG[, g, drop = FALSE], 1L, \(x) {
                                                  length(which(x %in% c(up, down)))
                                                }) > 0L)))
}
#
# Z-scored clustering heatmaps of regulated proteins
clustersTest <- try({
  clustMode <- "t-tests"
  clstSrc <- paste0(libPath, "/extdata/Sources/cluster_Heatmap_Main.R")
  #rstudioapi::documentOpen(clstSrc)
  source(clstSrc, local = FALSE)
}, silent = TRUE) # Allowed to fail, but with a warning!
if (inherits(clustersTest, "try-error")) {
  warning("Could not draw heatmap for t-test results!")
}
#

# Backup data/update cluster
stopClust <- TRUE
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - F-test
#Param <- Param.load()
F.test %<o% FALSE
if (("F.test" %in% colnames(Param))&&(is.logical(Param$F.test))&&(length(Param$F.test) == 1L)&&(!is.na(Param$F.test))&&(Param$F.test)) {
  F.test <- !((length(VPAL$values) == 2L)&&(pvalue.col[pvalue.use] == "Moderated t-test -log10(Pvalue) - "))
}
if (F.test) {
  dir <- paste0(wd, "/", c("Reg. analysis/F-tests"#,
                           #"Reg. analysis/F-tests/pdf",
                           #"Reg. analysis/F-tests/jpeg",
                           #"Reg. analysis/F-tests/html"
  ))
  for (d in dir) { if (!dir.exists(d)) { dir.create(d, recursive = TRUE) }}
  dirlist <- unique(c(dirlist, dir))
  msg <- "Performing F-test"
  ReportCalls <- AddMsg2Report(Space = FALSE, Print = FALSE)
  #
  if (("F.test_within" %in% colnames(Param))&&(Param$F.test_within != "")) {
    warning("Parameter \"F.test_within\" is deprecated!")
  }
  if (("F.test_factors" %in% colnames(Param))&&(Param$F.test_factors != "")) {
    warning("Parameter \"F.test_factors\" is deprecated!")
  }
  if (("F.test_factors_ref" %in% colnames(Param))&&(Param$F.test_factors_ref != "")) {
    warning("Parameter \"F.test_factors_ref\" is deprecated!")
  }
  #
  F_Root %<o% "mod. F-test -log10(Pvalue)"
  dataType <- "PG"
  #
  FSrc %<o% paste0(libPath, "/extdata/Sources/run_F_test.R")
  #rstudioapi::documentOpen(FSrc)
  tstFtst <- try(source(FSrc, local = FALSE), silent = TRUE)
  #
  if (!inherits(tstFtst, "try-error")) {
    F.test <- TRUE
    #F_test_ref_ratios %<o% F_volc$`Reference ratios` # Not needed
    volcano.plots$"F-tests_Unlabelled" <- F_volc$Plots$"Unlabelled"
    volcano.plots$"F-tests_Labelled" <- F_volc$Plots$"Labelled"
    n2 <- names(volcano.plots$"F-tests_Labelled")
    dir <- paste0(dir, "/Reg. analysis/F-tests")
    for (ttl in n2) {
      plot <- volcano.plots$"F-tests_Labelled"[[ttl]]
      ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
    }
    #
    # Create F-test filters:
    Freg_Root <- gsub(" -log10\\(Pvalue\\)", " Regulated", F_Root)
    pat <- topattern(paste0(Freg_Root, " - "))
    g <- grep(pat, colnames(F_test_data), value = TRUE)
    g1 <- gsub(pat, "", g)
    up <- grep("^up|^Specific", unique(unlist(F_test_data[, g])), value = TRUE)
    down <- grep("^down|^Anti-specific", unique(unlist(F_test_data[, g])), value = TRUE) # The "anti-specific" part will only become relevant if in future I disconnect symmetry and/or adding specific tags from pull-down experiments
    Reg_filters$"F-tests" <- list()
    if ("con" %in% filter_types) {
      Reg_filters$"F-tests"$"By condition" <- list()
      for (i in seq_along(g)) { #i <- 1
        Reg_filters$"F-tests"$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                              Filter_up = sort(which(F_test_data[[g[i]]] %in% up)),
                                                              Filter_down = sort(which(F_test_data[[g[i]]] %in% down)),
                                                              Filter = sort(which(F_test_data[[g[i]]] %in% c(up, down))))
      }
    }
    if ("dat2" %in% filter_types) {
      Reg_filters$"F-tests"$"Whole dataset" <- list(Columns = g,
                                                    Filter_up = sort(which(apply(F_test_data[, g, drop = FALSE], 1L, \(x) {
                                                      length(which(x %in% up))
                                                    }) > 0L)),
                                                    Filter_down = sort(which(apply(F_test_data[, g, drop = FALSE], 1L, \(x) {
                                                      length(which(x %in% down))
                                                    }) > 0L)),
                                                    Filter = sort(which(apply(F_test_data[, g, drop = FALSE], 1L, \(x) {
                                                      length(which(x %in% c(up, down)))
                                                    }) > 0L)))
    }
    if (("Q.values" %in% colnames(Param))&&(Param$Q.values)) {
      require(qvalue)
      if (F_Root %in% colnames(F_test_data)) {
        temp <- 10^(-F_test_data[[F_Root]])
        wag <- which(is.all.good(temp, 2L))
        pi0 <- qvalue::pi0est(temp)$pi0[1L]
        temp <- try(qvalue::qvalue(temp[wag]), silent = TRUE) # For now we do not explicitly set pi0
        if (inherits(temp, "try-error")) {
          temp <- try(qvalue::qvalue(temp[wag], pi0 = pi0), silent = TRUE)
          while ((inherits(temp, "try-error"))&&(pi0 <= 1)) {
            pi0 <- pi0 + 0.05
            temp <- try(qvalue(temp[wag], pi0 = pi0), silent = TRUE)
          }
        }
        if (!inherits(temp, "try-error")) {
          F_test_data[["-log10(Qvalue)"]] <- NA
          F_test_data[["local FDR"]] <- NA
          F_test_data[wag, "-log10(Qvalue)"] <- -log10(temp$qvalues)
          F_test_data[wag, "local FDR"] <- temp$lfdr
        } else { warning(paste0("F-test: Q-values calculation failed for ", F_Root, "! No q-values column will be created.")) }
      }
    }
    #
    # Z-scored clustering heatmaps of regulated proteins
    clustersTest <- try({
      clustMode <- "F-tests"
      clstSrc <- paste0(libPath, "/extdata/Sources/cluster_Heatmap_Main.R")
      #rstudioapi::documentOpen(clstSrc)
      source(clstSrc, local = FALSE)
    }, silent = TRUE) # Allowed to fail, but with a warning!
    if (inherits(clustersTest, "try-error")) {
      warning("Could not draw heatmap for F-test results!")
    }
    #
    # Backup data/update cluster
    #rstudioapi::documentOpen(bckpSrc)
    source(bckpSrc, local = FALSE)
    #loadFun(BckUpFl)
    #
  } else { warning("F-test analysis failed, check your parameters!")}
}
# Mat-meth text
tmp <- BH.FDR*100
l <- length(tmp)
if (l > 1L) { tmp <- paste0(paste(tmp[1L:(l-1L)], collapse = "%, "), " and ", tmp[l], "%") }
tmp2 <- Param$Ratios.Contamination.Rates
tmpPVal <- gsub(" -log10\\(pvalue\\) - ", "",
                gsub("welch", "Welch",
                     gsub("student", "Student", tolower(pvalue.col[which(pvalue.use)]))))
if (grepl("^moderated", tmpPVal, ignore.case = TRUE)) { tmpPVal <- paste0(tmpPVal, " (limma)") }
if (grepl("^DEqMS mod\\.", tmpPVal, ignore.case = TRUE)) { tmpPVal <- paste0(gsub("deqms mod\\. ", "DEqMS moderated ", tmpPVal), " (limma + edge)") }
if (grepl("^((EBA)|(S))AM ", tmpPVal, ignore.case = TRUE)) { tmpPVal <- paste0(tmpPVal, " (siggenes)") }
if (grepl("^permutations ", tmpPVal, ignore.case = TRUE)) { tmpPVal <- paste0(tmpPVal, " (coin)") }
if (grepl("^(ODP)|(LRT) ", tmpPVal, ignore.case = TRUE)) { tmpPVal <- paste0(tmpPVal, " (edge)") }
l <- length(DatAnalysisTxt)
DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                            " Average log10 expression values were tested for significance using a ",
                            #c("two", "one")[match(AltHyp, c("two.sided", "greater", "lower"))], "-sided ", # Nope!
                            tmpPVal, " per samples group",
                            c("", " and a moderated ANOVA (limma, F-test run using voomaLmFit for heteroskedasticity correction with individual moderated t-tests as post-hoc tests)")[F.test+1L],
                            ". Significance thresholds were calculated using the Benjamini-Hochberg procedure for False Discovery Rate (FDR) values of ", tmp,
                            ". For all tests, differentially expressed protein groups were defined as those with a significant P-value and a",
                            c("n absolute", "")[IsPullDown+1L],
                            " log2 average ratio greater than ",
                            c(paste0(Param$Ratios.Contamination.Rates*100, "% of ",
                                     c(paste0("control-to", c("-average", "")[Nested+1L],
                                              "-control"),
                                       "intra-sample groups")[match(RefRat_Mode, c("1", "2"))],
                                     " ratios"),
                              Param$Ratios.Contamination.Rates)[match(Param$Ratios.Thresholds,
                                                                      threshOpt)], ".")

#### Code chunk - SAINTexpress
Src <- paste0(libPath, "/extdata/Sources/SAINTexpress.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Now let's create a table of regulated protein groups per test made:
g <- grep("^Regulated - ", colnames(PG), value = TRUE)
regPG_TTest <- data.frame(Test = gsub("^Regulated - ", "", g))
dir <- c("up", "down")
kolstms <- c("count", "PG IDs", "Leading Protein IDs", "Genes")
DF <- data.frame(dir = rep(dir, length(BH.FDR)),
                 FDR = unlist(lapply(BH.FDR, function(f) { rep(f, length(dir)) })))
regPG_TTest <- lapply(1L:nrow(DF), \(ii) { #ii <- 1L
  d <- DF$dir[[ii]]
  i <- match(DF$FDR[[ii]], BH.FDR)
  tmp <- paste0(d, ", FDR = ", BH.FDR[1L:i]*100, "%")
  kolnms <- paste0(tmp[i], " - ", kolstms)
  tmp2 <- set_colnames(Isapply(g, \(x) {
    w <- which(PG[[x]] %in% tmp)
    x1 <- length(w)
    x2 <- paste0(PG$id[w], collapse = ", ")
    x3 <- paste0(PG$"Leading protein IDs"[w], collapse = ", ")
    x4 <- paste0(PG$Genes[w], collapse = ", ")
    return(c(x1, x2, x3, x4))
  }), kolnms)
  tmp2[[kolnms[1L]]] <- as.numeric(tmp2[[kolnms[1L]]])
  tst <- unique(tmp2[[kolnms[1L]]])
  tst <- tst[which(tst > 0L)]
  if (!length(tst)) { regPG_TTest[, kolnms] <- tmp2 }
  return(tmp2)
})
regPG_TTest <- cbind(gsub("^Regulated - ", "", g), do.call(cbind, regPG_TTest))
colnames(regPG_TTest)[1L] <- "Contrast"
if (IsPullDown) {
  tmp <- grep("^Specific: ", unique(unlist(PG[, g])), value = TRUE)
  if (length(tmp)) {
    kolnms <- paste0("Specific - ", kolstms)
    tmp2 <- set_colnames(Isapply(g, \(x) {
      w <- which(PG[[x]] %in% tmp)
      x1 <- length(w)
      x2 <- paste0(PG$id[w], collapse = ", ")
      x3 <- paste0(PG$"Leading protein IDs"[w], collapse = ", ")
      x4 <- paste0(PG$Genes[w], collapse = ", ")
      return(c(x1, x2, x3, x4))
    }), kolnms)
    tmp2[[kolnms[1L]]] <- as.numeric(tmp2[[kolnms[1L]]])
    tst <- unique(tmp2[[kolnms[1L]]])
    tst <- tst[which(tst > 0)]
    if (length(tst)) { regPG_TTest[, kolnms] <- tmp2 }
  }
}
if (ncol(regPG_TTest) > 1L) {
  dir <- paste0(wd, "/Tables")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  regPG_TTest <- regPG_TTest[, c("Contrast", as.character(sapply(kolstms, \(x) {
    grep(paste0(" - ", x, "$"), colnames(regPG_TTest), value = TRUE)
  })))]
  write.csv(regPG_TTest, file = paste0(dir, "/Reg. PGs - t-test.csv"), row.names = FALSE)
} else { warning("The t-test(s) did not identify any regulated protein groups!") }
if (F.test) {
  g <- grep("^mod\\. F-test Regulated - ", colnames(F_test_data), value = TRUE)
  # For the F-test we always want to include both directions because there can be up and down for a pull-down if doing a secondary comparison
  regPG_FTest <- lapply(1L:nrow(DF), \(ii) { #ii <- 1L
    d <- DF$dir[[ii]]
    i <- match(DF$FDR[[ii]], BH.FDR)
    tmp <- paste0(d, ", FDR = ", BH.FDR[1L:i]*100, "%")
    kolnms <- paste0(tmp[i], " - ", kolstms)
    tmp2 <- set_colnames(Isapply(g, \(x) {
      w <- which(F_test_data[[x]] %in% tmp)
      x1 <- length(w)
      x2 <- paste0(PG$id[w], collapse = ", ")
      x3 <- paste0(PG$"Leading protein IDs"[w], collapse = ", ")
      x4 <- paste0(PG$Genes[w], collapse = ", ")
      return(c(x1, x2, x3, x4))
    }), kolnms)
    tmp2[[kolnms[1L]]] <- as.numeric(tmp2[[kolnms[1L]]])
    tst <- unique(tmp2[[kolnms[1L]]])
    tst <- tst[which(tst > 0L)]
    if (!length(tst)) { regPG_TTest[, kolnms] <- tmp2 }
    return(tmp2)
  })
  regPG_FTest <- cbind(gsub("^mod\\. F-test Regulated - ", "", g), do.call(cbind, regPG_FTest))
  colnames(regPG_FTest)[1L] <- "Contrast"
  if (IsPullDown) {
    tmp <- grep("^Specific: ", unique(unlist(F_test_data[, g])), value = TRUE)
    if (length(tmp)) {
      kolnms <- paste0("Specific - ", kolstms)
      tmp2 <- set_colnames(Isapply(g, \(x) {
        w <- which(F_test_data[[x]] %in% tmp)
        x1 <- length(w)
        x2 <- paste0(PG$id[w], collapse = ", ")
        x3 <- paste0(PG$"Leading protein IDs"[w], collapse = ", ")
        x4 <- paste0(PG$Genes[w], collapse = ", ")
        return(c(x1, x2, x3, x4))
      }), kolnms)
      tmp2[[kolnms[1L]]] <- as.numeric(tmp2[[kolnms[1L]]])
      tst <- unique(tmp2[[kolnms[1L]]])
      tst <- tst[which(tst > 0L)]
      if (length(tst)) { regPG_FTest[, kolnms] <- tmp2 }
    }
  }
  if (ncol(regPG_FTest) > 1L) {
    dir <- paste0(wd, "/Tables")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    regPG_FTest <- regPG_FTest[, c("Contrast", as.character(sapply(kolstms, \(x) {
      grep(paste0(" - ", x, "$"), colnames(regPG_FTest), value = TRUE)
    })))]
    write.csv(regPG_FTest, file = paste0(dir, "/Reg. PGs - F-test.csv"), row.names = FALSE)
  } else { warning("The F-test did not identify any regulated protein groups!") }
}
# To do: SAINTexpress!

# Summary table and heatmap of number of regulated protein groups
Tsts <- c("t-tests", "F-tests")
WhTsts <- which(Tsts %in% names(Reg_filters))
for (tt in WhTsts) { #tt <- WhTsts[1L]
  tstrt <- Tsts[tt]
  stopifnot(!is.na(tstrt))
  dir <- paste0(wd, "/Reg. analysis/", tstrt)
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  filt <- Reg_filters[[tstrt]]
  By <- c("By condition", "By reference", "By analysis")
  By <- By[which(By %in% names(filt))]
  if (length(By)) {
    for (bee in By) { #bee <- By[1L]
      flt <- filt[[bee]]
      if (length(flt) >= 2L) {
        #flt <- flt[order(names(flt))]
        flt <- rev(flt)
        N <- length(flt)
        temp <- as.data.frame(matrix(rep("", (N+1L)^2L), ncol = N+1L))
        temp[2L:(N+1L), 1L] <- temp[1L, 2L:(N+1L)] <- names(flt)
        for (i in 2L:(N+1L)) {
          x <- flt[[temp[i, 1L]]]$Filter
          temp[i, 2L:(N+1)] <- vapply(temp[1L, 2L:(N+1L)], \(y) { sum(x %in% flt[[y]]$Filter) }, 1L)
        }
        nms <- names(flt)
        #nms <- cleanNms(nms)
        # tst <- lengths(strsplit(nms, " - "))
        # tst <- (min(tst) > 1L)&(length(unique(tst)) == 1L)
        # if (tst) {
        #   tst <- as.data.frame(t(sapply(strsplit(nms, " - "), unlist)))
        #   l <- apply(tst, 2L, \(x) { length(unique(x)) })
        #   tst <- tst[, which(l > 1L), drop = FALSE]
        #   nms <- do.call(paste, c(tst, sep = " - "))
        # }
        temp[2L:(N+1L), 1L] <- temp[1L, 2L:(N+1L)] <- nms
        nm <- paste0("N. of co-regulated PGs\n", tstrt, "\n(", tolower(bee), ")")
        write.csv(temp, file = paste0(dir, "/", gsub("\n", " - ", gsub("\n\\(", " (", nm)), ".csv"), row.names = FALSE)
        temp2 <- temp[2L:(N+1L), 2L:(N+1L)]
        colnames(temp2) <- temp[1L, 2L:(N+1L)]
        rownames(temp2) <-  temp[2L:(N+1L), 1L]
        for (i in 1L:nrow(temp2)) { temp2[[i]] <- as.numeric(temp2[[i]]) }
        if (max(is.all.good(unlist(temp2)))) {
          temp2 <- as.matrix(temp2)
          basic.heatmap(temp2,
                        "N. of co-regulated PGs",
                        paste0(tstrt, "\n(", tolower(bee), ")"),
                        save = c("pdf", "jpeg"),
                        folder = dir)
        } else { warning(paste0("Not a single regulated protein group in any of the ", tstrt, " performed, skipping.")) }
      } else {
        bb <- gsub("By ", "", bee)
        msg <- if (bb == "condition") { paste0(tstrt, " analysis: only one ", bb, " tested, skipping.") } else {
          paste0(tstrt, ": only one ", bb, " ", c("performed", "used")[(bb == "reference")+1L], " -> skipping.")
        }
        warning(msg)
      }
    }
  } else { warning(paste0(tstrt, " filters only defined for dataset, skipping.")) }
}

# Gene-Set Enrichment Analysis (GSEA)
if (!exists("runGSEA")) { runGSEA %<o% FALSE }
if (runGSEA) {
  dataType <- "PG"
  GSEAmode <- "standard"
  Src <- paste0(libPath, "/extdata/Sources/GSEA.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

# Weighted Gene Correlation Networks Analysis (WGCNA)
if (!exists("runWGCNA")) { runWGCNA %<o% FALSE }
if (runWGCNA) {
  Src <- paste0(libPath, "/extdata/Sources/WGCNA.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

# Backup data/update cluster
stopClust <- TRUE
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Dimensionality reduction plots
dataType <- "PG"
dmrdSrc <- paste0(libPath, "/extdata/Sources/dimRed_plots.R")
#rstudioapi::documentOpen(dmrdSrc)
source(dmrdSrc, local = FALSE)

#### Code chunk - Protein group profile plots and ranked abundance plots
PrfRASrc %<o% paste0(libPath, "/extdata/Sources/profile_and_rankedAbund_plots.R")
#rstudioapi::documentOpen(PrfRASrc)
source(PrfRASrc, local = FALSE)

# Visualize results
if (!exists("xplorSrc")) {
  xplorSrc <- paste0(libPath, "/extdata/Sources/xplorData.R") # Backwards compatibility
}
xplorSrc %<o% xplorSrc
#rstudioapi::documentOpen(xplorSrc)
source(xplorSrc, local = FALSE)

#### Code chunk - Optional: if there are time points, plot the curve of the ratios of one or all protein(s) over time
if (exists("Tim")) {
  msg <- "Time profile plots"
  ReportCalls <- AddMsg2Report(Space = FALSE)
  dir <- paste0(wd, "/Time profile plots")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  r <- paste0("Mean ", Prot.Rat.Root)
  p <- pvalue.col[which(pvalue.use)]
  a <- gsub("Tim", "", VPAL$aggregate)
  o <- parse.Param.aggreg(Param_filter(Param$Volcano.plots.Aggregate.Level, "Tim"))$names
  A <- get(a)
  a1 <- VPAL$aggregate
  o1 <- VPAL$names
  A1 <- get(a1)
  ylim <- paste0(r, A1)
  ylim <- ylim[which(ylim %in% colnames(PG))]
  ylim <- max(is.all.good(unlist(PG[,ylim])))*1.05
  temp <- list()
  for (i in A) { #i <- A[1L]
    i1 <- unlist(strsplit(i, "___"))
    e <- lapply(seq_along(o), \(x) { which(Exp.map[[o[x]]] == i1[x]) })
    l <- unique(unlist(e))
    t <- vapply(l, \(x) { length(which(unlist(e) == x)) == length(o) }, TRUE)
    l <- l[which(t)]
    e <- Exp.map[l,]
    t1 <- paste0(r, e[[a1]])
    t2 <- paste0(p, e[[a1]])
    w <- which((t1 %in% colnames(PG))&(t2 %in% colnames(PG)))
    if (length(w)) {
      e <- e[w,]
      t1 <- unique(t1[w])
      t2 <- unique(t2[w])
      tp <- unique(e[[Aggregates[which(names(Aggregates) == "Tim")]]])
      tp <- sort(as.numeric(tp))
      if (length(t1) <= 1L) {
        if (!length(t1)) { cat("   There is no valid data for aggregate", i, "\n")
        } else { cat("   There is only a single time point for aggregate", i, "\n") }
      } else {
        test <- apply(PG[,c(t1, t2)], 1L, \(x) { length(is.all.good(x)) == length(tp)*2L })
        col <- c("Protein IDs", "Names", "ID")
        col <- col[which(col %in% colnames(PG))]
        temp1 <- PG[which(test), c(col, Param$Plot.labels, t1)]
        temp1$IDs <- as.character(1L:nrow(temp1))
        temp2 <- PG[which(test), c(col, Param$Plot.labels, t2)]
        temp1 <- reshape2::melt(temp1, id.vars = c(col, "IDs", Param$Plot.labels))
        colnames(temp1)[which(colnames(temp1) == "value")] <- "log2(Ratio)"
        temp1$variable <- gsub_Rep(topattern(r, start = FALSE), "", as.character(temp1$variable))
        temp2 <- reshape2::melt(temp2, id.vars = c(col, Param$Plot.labels))
        temp1$"-log10(Pvalue)" <- temp2$value[match(temp1$"Protein IDs", temp2$"Protein IDs")]
        temp1[,o1] <- Isapply(strsplit(temp1$variable, "___"), unlist)
        temp[[i]] <- temp1
      }
    }
  }
  if (length(temp)) {
    kount <- 0L
    for (tp in A) {
      if (tp %in% names(temp)) {
        kount <- kount + 1L
        if (kount == 1L) {
          tmp <- temp[[tp]]
          tmp$Aggregate <- tp
        } else {
          tmp2 <- temp[[tp]]
          tmp2$Aggregate <- tp
          tmp <- rbind(tmp, tmp2)
        }
      }
    }
    tmp$Aggregate <- cleanNms(tmp$Aggregate)
    tp2 <- cleanNms(tp)
    tmp$Label <- tmp[[Param$Plot.labels]]
    tmp$"Time point" <- as.numeric(tmp[[Aggregates[which(names(Aggregates) == "Tim")]]])
    tmp$"Time point" <- factor(tmp$"Time point", levels = as.character(sort(as.numeric(unique(tmp$"Time point")))))
    levels(tmp$`Time point`) <- sort(as.numeric(levels(tmp$`Time point`)))
    ttl <- paste0("Global time profile - ", tp2)
    plot <- ggplot(tmp) +
      geom_line(aes(x = `Time point`, y = `log2(Ratio)`, group = IDs, color = IDs)) +
      scale_color_viridis_d(begin = 0.25) +
      ggtitle(ttl) + guides(color = "none", group = "none") +
      facet_wrap(~Aggregate) +
      ylim(c(-ylim, ylim)) + theme_bw()
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 600L, width = 10L, height = 10L, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 600L, width = 10L, height = 10L, units = "in")
    })
    ReportCalls <- AddPlot2Report()
    if (create_plotly) {
      #test <- aggregate(tmp$`log2(Ratio)`, list(tmp$IDs), \(x) {length(is.all.good(x)) == length(Tim)-1})
      tmp2 <- aggregate(tmp$`log2(Ratio)`, list(tmp$IDs), \(x) { max(abs(is.all.good(x)))})
      tmp2 <- tmp2[order(tmp2$x, decreasing = TRUE),]
      tmp2 <- tmp2$Group.1[1L:min(c(1000L, nrow(tmp2)))]
      tmp2 <- tmp[which(tmp$IDs %in% tmp2),]
      plot <- ggplot(tmp2) +
        geom_line(aes(x = `Time point`, y = `log2(Ratio)`, group = IDs, color = IDs, text = Label)) +
        scale_color_viridis_d(begin = 0.25) +
        ggtitle(ttl) + guides(color = "none", group = "none") +
        facet_wrap(~Aggregate) + ylim(c(-ylim, ylim)) + theme_bw()
      plot.ly <- ggplotly(plot, tooltip = "text")
      if (create_plotly_local) {
        saveWidget(plot.ly, paste0(dir, "/", ttl, ".html"))
        system(paste0("open \"", dir, "/", ttl, ".html"))
      } else {
        poplot(plot)
        # Legacy code for web-hosted plotly plots:
        plot_ly_object <- api_create(x = plot.ly, filename = paste0(plotly_subfolder, ttl),
                                     fileopt = "overwrite", sharing = "secret")
        plot_ly$"Time profile" <- list()
        plot_ly$"Time profile"$"Global time profile" <- plot_ly_object  
      }
    } else { poplot(plot) }
  }
}

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)

#### Code chunk - Sub-Cellular localisation analysis
# Includes:
# - pRoloc-based prediction of localisation
# - Re-localisation analysis based on testing the Sums of Squared Differences of Profiles
Src <- paste0(libPath, "/extdata/Sources/SubCellLoc.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

### Check that Cytoscape is installed and can run, then launch it.
Src <- paste0(libPath, "/extdata/Sources/Cytoscape_init.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
runClueGO <- runClueGO&CytoScape

#### Code chunk - Gene Ontology terms enrichment analysis
goSrc <- paste0(libPath, "/extdata/Sources/rep_GO.R")
#rstudioapi::documentOpen(goSrc)
source(goSrc, local = FALSE)

# Backup data/update cluster
stopClust <- TRUE
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Modified peptides analysis
Src <- paste0(libPath, "/extdata/Sources/Cytoscape_init.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
modPepSrc <- paste0(libPath, "/extdata/Sources/modPeptides.R")
#rstudioapi::documentOpen(modPepSrc)
source(modPepSrc, local = FALSE)

# Backup data/update cluster
stopClust <- TRUE
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Proteomic ruler
if (protrul) {
  ProtRulRoot %<o% "log10(est. copies/cell) - "
  tempPG <- PG
  exprsRt <- paste0("Mean ", prtRfRoot)
  if (LocAnalysis) {
    tempPG <- tempPG[, grep(topattern(exprsRt), colnames(tempPG), invert = TRUE)]
    for (grp2 in SubCellFracAggr2$values) { #grp2 <- SubCellFracAggr2$values[1L]
      em2 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == grp2),]
      for (grp in unique(em2[[SubCellFracAggr$column]])) { #grp <- unique(em2[[SubCellFracAggr$column]])[1L]
        em <- em2[which(em2[[SubCellFracAggr$column]] == grp),]
        kol <- paste0(prtRfRoot, em$Ref.Sample.Aggregate)
        tempPG[[paste0(prtRfRoot, grp)]] <- apply(10^tempPG[, kol], 1L, \(x) {
          log10(sum(is.all.good(x)))
        })
      }
      tempPG[[paste0(exprsRt, grp2)]] <- apply(tempPG[, paste0(prtRfRoot, unique(em2[[SubCellFracAggr$column]]))], 1L, \(x) {
        mean(is.all.good(x))
      })
    }
  }
  temp <- try(Prot.Ruler(tempPG, db, exprsRt, NuclL = ProtRulNuclL), silent = TRUE)
  if (is.list(temp)) {
    db <- temp$Database
    temp <- temp$Protein.groups
    kol <- c(grep(topattern(exprsRt), colnames(temp), value = TRUE), grep(topattern(ProtRulRoot), colnames(temp), value = TRUE))
    PG[, kol] <- temp[, kol]
    protrul <- TRUE
    l <- length(DatAnalysisTxt)
    DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                                " Protein group copy numbers per cell were estimated using a variant of the proteome ruler logic, normalizing to scaled values of all identified histones.")
  } else {
    warning("Failed to run Prot.Ruler function; is the remote NCBI server available?")
    protrul <- FALSE
  }
  rm(temp)
}

#### Code chunk - Summary table and QC plots
if ((create_plotly)&&(!create_plotly_local)) { # This code is so old it's auld!
  plot_ly_addresses %<o% data.frame(`Plot type` = NA, `Name` = NA, `Address` = NA)
  kount <- 1L
  for (n in seq_along(plot_ly)) { #n <- 1L
    i <- plot_ly[[n]]
    if (length(i)) {
      plot_ly_addresses[kount,] <- c(names(plot_ly)[n], "", "")
      for (j in seq_along(i)) { #j <- 1L
        kount <- kount + 1L
        k <- i[[j]]
        plot_ly_addresses[kount,] <- c("",
                                       names(i)[j],
                                       paste0(gsub("\\.embed$", "", k$embed_url), "/?share_key=", k$share_key))
      }
      kount <- kount + 1L
    }
  }
  class(plot_ly_addresses$Address) <- "hyperlink"
  wb <- createWorkbook()
  sheet  <- addWorksheet(wb, "Plotly plot addresses")
  writeData(wb, sheet = "Plotly plot addresses", x = plot_ly_addresses)
  style1 <- createStyle(textDecoration = "bold")
  addStyle(wb, sheet = "Plotly plot addresses", style = style1, rows = 2L:nrow(plot_ly_addresses)+1L, cols = 1L, gridExpand = FALSE, stack = FALSE)
  style2 <- createStyle()
  addStyle(wb, sheet = "Plotly plot addresses", style = style2, rows = 2L:nrow(plot_ly_addresses)+1L, cols = 2L, gridExpand = FALSE, stack = FALSE)
  style3 <- createStyle(textDecoration = "italic")
  addStyle(wb, sheet = "Plotly plot addresses", style = style3, rows = 2L:nrow(plot_ly_addresses)+1L, cols = 3L, gridExpand = FALSE, stack = FALSE)
  setColWidths(wb, sheet = "Plotly plot addresses", cols = 1L, widths = 25L)
  setColWidths(wb, sheet = "Plotly plot addresses", cols = 2L, widths = 50L)
  setColWidths(wb, sheet = "Plotly plot addresses", cols = 3L, widths = 60L)
  freezePane(wb, 1L, firstRow = TRUE)
  saveWorkbook(wb, file = "Tables/Plotly plot addresses.xlsx", overwrite = TRUE)
}
mods <- setNames(Modifs$Mark[which(Modifs$Type == "Variable")],
                 nm = Modifs$"Full name"[which(Modifs$Type == "Variable")])
tmp <- aggregate(Frac.map$"Raw file", list(Frac.map$MQ.Exp), length)
tmp <- round(mean(tmp$x)) # Size of a full fraction set, rounding for cases where we removed some fractions
defSc <- 60L # (Non-strict) default maximum number of files to look at per plot
if (tmp > defSc) {
  # If one set of fractions is larger than defSc
  sc <- tmp
} else {
  # What is closest to default: n or n+1 full sets of fractions?
  tst <- (defSc %% tmp) >= defSc/2L
  # Identify N = fixed number of files from full fraction sets we can fit in one plot
  N <- c(floor(defSc/tmp), ceiling(defSc/tmp))[tst+1L]*tmp
  # If we divide the total number of files by that number, how many plots do we have?
  Nplts <- ceiling(nrow(Frac.map)/N)
  # So that makes how many files per plot:
  sc <- ceiling(nrow(Frac.map)/Nplts)
}
sc <- max(c(sc, 1L))
source(parSrc, local = FALSE)
Exp_summary %<o% MQ.summary(ev = ev, pg = PG, wd = wd, mods = mods, save = "pdf",
                            raw.files = rawFiles, sc = sc, cl = parClust,
                            MQtxt = inDirs[which(SearchSoft == "MAXQUANT")])
Exp_summary$"Biological sample" <- ""
if ("Parent sample" %in% colnames(Frac.map)) {
  Exp_summary$"Biological sample" <- Frac.map$"Parent sample"[match(Exp_summary$Sample, Frac.map$"Raw file")]
}
if ("MQ.Exp" %in% colnames(Frac.map)) {
  Exp_summary$"Biological sample" <- Frac.map$MQ.Exp[match(Exp_summary$Sample, Frac.map$"Raw file")]
}
Exp_summary$"Biological sample"[1L] <- "All samples"
Exp_summary <- Exp_summary[, c("Sample", "Biological sample",
                               colnames(Exp_summary)[which(!colnames(Exp_summary) %in% c("Sample", "Biological sample"))])]
write.csv(Exp_summary, paste0(wd, "/Workflow control/Summary.csv"), row.names = FALSE)
#Exp_summary <- read.csv(paste0(wd, "/Workflow control/Summary.csv"), check.names = FALSE)

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)

#### Code chunk - XML coverage columns
Src <- paste0(libPath, "/extdata/Sources/xml_Coverage_columns.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
# To do: also PTMs in a different color (one for all)

#### Code chunk - GO term columns
GO_PG_col %<o% unique(unlist(strsplit(Param$GO.tabs, ";")))
GO_filt %<o% length(GO_PG_col) > 0L
if (GO_filt) {
  if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) { loadFun("GO_terms.RData") }
  GO_PG_col <- GO_PG_col[which(GO_PG_col %in% GO_terms$ID)]
  GO_filt <- length(GO_PG_col) > 0L
}
if (GO_filt) {
  tmp <- listMelt(strsplit(PG$`GO-ID`, ";"), 1L:nrow(PG), c("Term", "Row"))
  Offspring <- setNames(lapply(GO_PG_col, \(x) { #x <- "GO:0009725"
    ont <- Ontology(x)
    x <- c(x, get(paste0("GO", ont, "OFFSPRING"))[[x]])
    x <- x[which(!is.na(x))]
    return(x)
  }), GO_PG_col)
  tmp <- tmp[which(tmp$Term %in% unlist(Offspring)),]
  tmp <- aggregate(tmp$Row, list(tmp$Term), c)
  colnames(tmp) <- c("Term", "Rows")
  w <- which(vapply(GO_PG_col, \(x) { sum(Offspring[[x]] %in% tmp$Term) }, 1L) == 0L)
  if (length(w)) {
    msg <- stringi::stri_join("No proteins found for the following GO terms:",
                              stringi::stri_join("\n - ", GO_terms$Term[match(GO_PG_col[w], GO_terms$ID)],
                                                 collapse = ""),
                              "\n\n")
    cat(msg)
    GO_PG_col <- GO_PG_col[which(GO_PG_col %in% tmp$Term)]
  }
  if (length(GO_PG_col)) {
    GO_PG_col2 %<o% setNames(GO_terms$Term[match(GO_PG_col, GO_terms$ID)],
                             GO_PG_col)
    PG[, GO_PG_col2] <- ""
    for (go in GO_PG_col) { #go <- GO_PG_col[1L]
      w <- which(tmp$Term %in% Offspring[[go]])
      w2 <- unique(unlist(tmp$Rows[w]))
      PG[w2, GO_PG_col2[go]] <- "+"
    }
    #View(PG[, GO_PG_col2])
    tst <- setNames(vapply(GO_PG_col2, \(x) {
      sum(PG[[x]] == "+")
    }, 1L), GO_PG_col2)
    tst <- paste0(GO_PG_col2, " -> ", tst, collapse = "\n - ")
    tst <- paste0("Number of protein groups per GO term of interest\n - ", tst, "\n\n")
    cat(tst)
    write(tst, paste0(wd, "/Reg. analysis/GO enrich/Dataset/GO_terms_of_interest.txt"))
  }
}

#### Code chunk - Create output tables
## PSMs
dir <- paste0(wd, "/Tables")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
w <- which(vapply(colnames(ev), \(x) { is.list(ev[[x]]) }, TRUE))
if (length(w)) { for (i in w) { ev[[i]] <- parSapply(parClust, ev[[i]], paste, collapse = ";") } }
data.table::fwrite(ev, paste0(dir, "/evidence.tsv"), sep = "\t", row.names = FALSE, na = "NA")
#
## Main peptidoforms- and protein groups-level, multi-tabs report
xlSrc <- paste0(libPath, "/extdata/Sources/rep_Write_Excel.R")
#rstudioapi::documentOpen(xlSrc)
source(xlSrc, local = FALSE)
#xl_open(repFl)

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
source(parSrc, local = FALSE)#### Code chunk - Amica input tables
# Write tables for Amica input:
## PG table
if (Param$Amica) {
  cat("Also writing Amica table...\n")
  w <- which(as.logical(Exp.map$Use))
  AmicaDesign <- data.frame(groups = cleanNms(Exp.map[w, VPAL$column], rep = "."),
                            samples = cleanNms(Exp.map[w, RSA$column], rep = "."))
  AmicTbl <- data.frame(Majority.protein.IDs = PG$"Leading protein IDs",
                        Gene.names = PG$Genes,
                        razorUniqueCount = PG$"Razor + unique peptides",
                        Potential.contaminant = PG$"Potential contaminant")
  tmp <- data.frame(IDs = PG$"Peptide IDs", Razor = PG$"Peptide is razor")
  tmp$Razor <- lapply(strsplit(tmp$Razor, ";"), \(x) {
    as.logical(toupper(x))
  })
  tmp$IDs <- lapply(strsplit(tmp$IDs, ";"), as.numeric)
  tmp$RazorIDs <- apply(tmp[, c("IDs", "Razor")], 1L, \(x) {
    x[[1L]][which(x[[2L]])]
  })
  for (i in Exp.map$Ref.Sample.Aggregate[which(as.logical(Exp.map$Use))]) { #i <- Exp.map$Ref.Sample.Aggregate[which(as.logical(Exp.map$Use))][1L]
    i2 <- cleanNms(i, rep = ".")
    kol <- paste0("LFQIntensity_", i2)
    AmicTbl[[kol]] <- PG[[paste0(prtRfRoot, i)]]/log10(2L)
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2L)), kol] <- NaN
    kol <- paste0("razorUniqueCount_", i2)
    tmp$Tmp <- strsplit(PG[[paste0("Peptide IDs - ", i)]], ";")
    AmicTbl[[kol]] <- apply(tmp[, c("IDs", "Tmp")], 1L, \(x) {
      sum(x[[2L]] %in% x[[1L]])
    })
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2L)), kol] <- 0L
  }
  tmp <- Exp.map$Ref.Sample.Aggregate[which(as.logical(Exp.map$Use))]
  kol <- paste0("LFQIntensity_", cleanNms(tmp, rep = "."))
  temp <- Data_Impute2(AmicTbl[, kol],
                       Exp.map[match(tmp, Exp.map$Ref.Sample.Aggregate), VPAL$column])
  temp <- temp$Imputed_data
  colnames(temp) <- gsub("^LFQIntensity_", "ImputedIntensity_", colnames(temp))
  AmicTbl[, colnames(temp)] <- temp
  w <- which(!Exp.map$Reference)
  grps <- unique(Exp.map[w, VPAL$column])
  for (g in grps) { #g <- grps[1L]
    gEd <- cleanNms(g, rep = ".")
    m <- Exp.map[which(Exp.map[[VPAL$column]] == g),]
    ratgrps <- unique(m[[RG$column]])
    g0 <- unique(Exp.map[which((Exp.map[[RG$column]] == ratgrps)&(Exp.map$Reference)),
                         VPAL$column])
    gEd0 <- cleanNms(g0, rep = ".")
    kol <- paste0("P.Value_", gEd, "__vs__", paste(gEd0, collapse = "&"))
    AmicTbl[[kol]] <- 10^(-PG[[paste0(pvalue.col[which(pvalue.use)], g)]])
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2L)), kol] <- NaN
    kol <- paste0("adj.P.Val_", gEd, "__vs__", paste(gEd0, collapse = "&"))
    PVkol <- paste0(pvalue.col[which(pvalue.use)], g)
    AmicTbl[[kol]] <- p.adjust(10^(-PG[[PVkol]]), method = "BH")
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2L)), kol] <- NaN
    kol <- paste0("logFC_", gEd, "__vs__", paste(gEd0, collapse = "&"))
    AmicTbl[[kol]] <- PG[[paste0("Mean ", Prot.Rat.Root, g)]]
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2L)), kol] <- NaN
    kol <- paste0("AveExpr_", gEd, "__vs__", paste(gEd0, collapse = "&"))
    AmicTbl[[kol]] <- PG[[paste0("Mean ", prtRfRoot, g)]]/log10(2)
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2L)), kol] <- NaN
  }
  tst <- apply(AmicTbl[, grep("^AveExpr_", colnames(AmicTbl), value = TRUE), drop = FALSE], 1L, \(x) {
    length(is.all.good(x))
  }) > 0L
  AmicTbl$quantified <- c("", "+")[((AmicTbl$razorUniqueCount >= 2L)&(tst))+1L]
  AmicTbl <- AmicTbl[which(AmicTbl$quantified == "+"),]
  dir <- paste0(wd, "/Amica")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  data.table::fwrite(AmicTbl, paste0(wd, "/Amica/Amica_file.csv"), row.names = FALSE, na = "NaN", sep = "\t", quote = FALSE)
  data.table::fwrite(AmicaDesign, paste0(wd, "/Amica/Experimental_design.csv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")
  #data.table::fwrite(AmicTbl[1L:500L,], paste0(wd, "/Amica/Amica_file_short.csv"), row.names = FALSE, na = "NaN", quote = FALSE)
}

# Backup data/update cluster
stopClust <- TRUE
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Venn diagrams
Src <- paste0(libPath, "/extdata/Sources/Venn_diagrams.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Coverage maps, XICs and heatmaps for proteins of interest
protlspep <- prot.list_pep
if (length(protlspep)) {
  test <- vapply(protlspep, \(i) { length(grsep2(i, PG$"Leading protein IDs")) }, 1L)
  if (0 %in% test) {
    w <- which(test == 0L)
    for (w1 in w) {
      m <- match(protlspep[w1], db$"Protein ID")
      nm <- paste0(db$"Protein ID"[m], " - ", db$"Common Name"[m])
      warning(paste0("Protein of interest ",nm, " was not found in the dataset!"))
    }
    protlspep <- protlspep[which(test > 0L)]
  }
}
if (length(protlspep)) { # Coverage
  setwd(wd) # To make sure we are in the working directory
  dir <- paste0(wd, "/Coverage")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  xKol <- paste0(pep.ref[length(pep.ref)], Exp.map$Ref.Sample.Aggregate) 
  tmpDB <- db[match(protlspep, db$"Protein ID"), c("Common Name", "Protein ID", "Sequence")]
  tst <- lapply(protlspep, \(x) { grsep2(x, pep$Proteins) })
  w <- which(lengths(tst) > 0L)
  prots <- protlspep[w]
  tst <- unique(unlist(tst))
  tmpPep <- pep[tst, c("Proteins", "Modified sequence", xKol)]
  source(parSrc, local = FALSE)
  clusterExport(parClust, list("tmpDB", "tmpPep", "pep.ref", "xKol", "wd", "VPAL", "Exp.map", "Exp"), envir = environment())
  lst <- parLapply(parClust, prots, \(i) { #i <- prots[1L]
    nm <- tmpDB$"Common Name"[match(i, tmpDB$"Protein ID")]
    nm <- gsub("[<>:\"/\\\\\\|\\?]", "-", nm)
    if (nchar(nm) > 20L) { nm <- paste0(gsub(" $", "", substr(nm, 1L, 17L)), "...") }
    seq <- setNames(tmpDB$Sequence[which(tmpDB$"Protein ID" == i)],
                    paste(tmpDB[which(tmpDB$"Protein ID" == i), c("Protein ID", "Common Name")], collapse = " - "))
    grs <- grsep2(i, tmpPep$Proteins)
    drLst <- dir <- paste0(wd, "/Coverage/", i)
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    p <- tmpPep[grs,]
    m <- apply(p[, xKol], 1L, \(x) {
      10^mean(is.all.good(log10(x)))
    })
    w <- which(m == 0)
    if (length(w)) { stop("I didn't expect this, investigate!")}
    p1 <- data.frame(Sequence = p$"Modified sequence",
                     Intensity = m)
    print(Coverage(seq, p1$Sequence))
    ttl <- gsub(":|/", "-", names(seq))
    setwd(dir) # To control precisely where it is saved
    Coverage(seq, p1$Sequence, Mode = "Align2", title = paste0("Coverage map - ", nm), save = c("jpeg", "pdf"),
             intensities = p1$Intensity, display = FALSE)
    setwd(wd) # To make sure I return to the working directory
    for (j in VPAL$values) { #j <- VPAL$values[1L]
      sm <- Exp.map[which(Exp.map[[VPAL$column]] == j),]
      m <- apply(p[, paste0(pep.ref[length(pep.ref)], sm$Ref.Sample.Aggregate)], 1L, \(x) {
        10L^mean(is.all.good(log10(x)))
      })
      w <- which(m > 0)
      if (length(w)) {
        p1 <- data.frame(Sequence = p$"Modified sequence"[w],
                         Intensity = m[w])
        dir2 <- paste0(dir, "/", gsub(":|\\*|\\?|<|>|\\|", "-", cleanNms(j, rep = "_")))
        if (!dir.exists(dir2)) { dir.create(dir2, recursive = TRUE) }
        drLst <- unique(c(drLst, dir2))
        ttl <- paste0(gsub(":|/", "-", names(seq)), " - ", cleanNms(j, rep = "_"))
        setwd(dir2) # To control precisely where it is saved
        Coverage(seq, p1$Sequence, Mode = "Align2", save = c("jpeg", "pdf"), title = paste0("Coverage map - ", nm),
                 intensities = p1$Intensity, display = FALSE)
        setwd(wd) # To make sure I return to the working directory
        for (k in sm$Ref.Sample.Aggregate) { #k <- sm$Ref.Sample.Aggregate[1L]
          m <- p[[paste0(pep.ref[length(pep.ref)], k)]]
          w <- which(m > 0)
          if (length(w)) {
            p1 <- data.frame(Sequence = p$"Modified sequence"[w],
                             Intensity = m[w])
            ttl <- paste0(gsub(":|/", "-", names(seq)), " - ", cleanNms(k, rep = "_"))
            setwd(dir2) # To control precisely where it is saved
            Coverage(seq, p1$Sequence, Mode = "Align2", save = c("jpeg", "pdf"),
                     title = paste0("Coverage map - ", nm),
                     intensities = p1$Intensity, display = FALSE)
            setwd(wd) # To make sure I return to the working directory
          }
        }
        setwd(wd) # To make sure I return to the working directory
      }
    }
    return(drLst)
  })
  dirlist <- unique(c(dirlist, unlist(lst)))
  setwd(wd) # To make sure I return to the working directory
}
if (length(protlspep)) { # XICs
  for (indir in inDirs[which(SearchSoft == "DIANN")]) {
    xicDir <- paste0(indir, "/report_xic")
    if (dir.exists(xicDir)) {
      XIC_fls <- list.files(xicDir, "\\.parquet$", full.names = TRUE)
      if (length(XIC_fls)) {
        require(arrow)
        source(parSrc, local = FALSE)
        g <- grsep2(protlspep, ev$Proteins)
        u <- ev$"Mod. seq. (DiaNN format)"[g]
        tmp <- Frac.map$`Raw files name`
        clusterExport(parClust, list("tmp", "g", "u"), envir = environment())
        XICs <- parLapply(parClust, XIC_fls, \(x) {
          res <- arrow::read_parquet(x)
          res$"Mod. seq." <- gsub_Rep("[0-9]+$", "", res$pr)
          res <- res[which(res$"Mod. seq." %in% u),]
          nm <- gsub(".*/|\\.xic\\.parquet$", "", x)
          res$File <- nm
          res$"Seq_Run" <- do.call(paste, c(res[, c("pr", "File")], sep = ">>>"))
          res$File <- factor(res$File, levels = tmp)
          res <- res[which(res$feature != "index"),]
          return(res)
        })
        #View(XICs[[1L]])
        XICs <- plyr::rbind.fill(XICs)
        #View(XICs[1L:100L,])
        #
        m <- match(XICs$"Mod. seq.", ev$"Mod. seq. (DiaNN format)")
        myKol <- c("Proteins", "Sequence", "Modified sequence", "PEP", "Quantity Quality")
        XICs[, myKol] <- ev[m, myKol]
        ev$tmp <- ">>>"
        Boundaries <- do.call(paste0, c(ev[, c("Mod. seq. (DiaNN format)", "Charge", "tmp", "Raw file")]))
        ev$tmp <- NULL
        w <- which(Boundaries %in% XICs$Seq_Run)
        Boundaries <- data.frame(Seq_Run = Boundaries[w],
                                 RT = ev$`Retention time`[w],
                                 `RT (start)` = ev$`Retention time (start)`[w],
                                 `RT (end)` = ev$`Retention time (end)`[w],
                                 check.names = FALSE)
        Boundaries$File <- gsub(".*>>>", "", Boundaries$Seq_Run)
        for (pr in protlspep) { #pr <- protlspep[1L]
          xicDir2 <- paste0(wd, "/XIC/", pr)
          if (!dir.exists(xicDir2)) { dir.create(xicDir2, recursive = TRUE) }
          dirlist <- unique(c(dirlist, xicDir2))
          g <- grsep2(pr, XICs$Proteins)
          if (length(g)) {
            XIC <- XICs[g,]
            pkBnds <- Boundaries[which(Boundaries$Seq_Run %in% XIC$Seq_Run),]
            u <- unique(XIC$"Modified sequence")
            clusterExport(parClust, list("XIC", "pkBnds", "xicDir2", "pr"), envir = environment())
            invisible(parLapply(parClust, u, \(sq) { #sq <- u[1L] #sq <- u[2L]
              sq2 <- gsub("^_|_$", "", sq)
              ppXIC <- XIC[which(XIC$"Modified sequence" == sq),]
              yMax <- aggregate(ppXIC$value, list(ppXIC$File), max)
              xMin <- min(ppXIC$rt)
              bnds <- pkBnds[which(pkBnds$Seq_Run %in% ppXIC$Seq_Run),]
              bnds$yMax <- yMax$x[match(bnds$File, yMax$Group.1)]
              wMS1 <- which(ppXIC$feature == "ms1")
              wMS2 <- which(ppXIC$feature != "ms1")
              aNNOt <- aggregate(ppXIC[, c("PEP", "Quantity Quality")], list(ppXIC$File), \(x) { signif(mean(x, na.rm = TRUE), 3L) })
              colnames(aNNOt)[1L] <- "File"
              aNNOt$PEP <- paste0("PEP = ", aNNOt$PEP)
              aNNOt$"Quantity Quality" <- paste0("Quantity Quality = ", aNNOt$"Quantity Quality")
              aNNOt$Text <- do.call(paste, c(aNNOt[, c("PEP", "Quantity Quality")], sep = "\n"))
              aNNOt$y <- yMax$x[match(aNNOt$File, yMax$Group.1)]
              plot <- ggplot2::ggplot() + ggplot2::scale_y_continuous(expand = c(0L, 10L))
              if (nrow(bnds)) {
                plot <- plot +
                  ggplot2::geom_rect(data = bnds, ggplot2::aes(xmin = `RT (start)`, ymin = 0, xmax = `RT (end)`, ymax = yMax),
                                     fill = "lightblue", alpha = 0.2) +
                  ggplot2::geom_vline(data = bnds, ggplot2::aes(xintercept = RT),
                                      color = "darkblue", linewidth = 0.5) +
                  ggplot2::geom_vline(data = bnds, ggplot2::aes(xintercept = `RT (start)`),
                                      color = "darkblue", linewidth = 0.5, linetype = "dashed") +
                  ggplot2::geom_vline(data = bnds, ggplot2::aes(xintercept = `RT (end)`),
                                      color = "darkblue", linewidth = 0.5, linetype = "dashed")
              }
              if (length(wMS1)) {
                plot <- plot +
                  ggplot2::geom_line(data = ppXIC[wMS1,], ggplot2::aes(x = rt, y = value), color = "red",
                                     linewidth = 0.5)
              }
              if (length(wMS2)) {
                plot <- plot +
                  ggplot2::geom_line(data = ppXIC[wMS2,], ggplot2::aes(x = rt, y = value, color = feature),
                                     linewidth = 0.3)
              }
              plot <- plot +
                ggplot2::geom_text(data = aNNOt, ggplot2::aes(label = Text, y = y), x = xMin, hjust = 0, vjust = 1, size = 3) +
                ggplot2::scale_colour_viridis_d() +
                ggplot2::facet_wrap(~File, scales = "free_y") + ggplot2::theme_bw() +
                ggplot2::ggtitle(paste0(pr, " ", sq2)) +
                ggplot2::xlab("Retention time") + ggplot2::ylab("Intensity")
              #poplot(plot, 12L, 22L)
              #
              suppressMessages({
                ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".jpeg"), plot, dpi = 450L, height = 10L, width = 10L)
                ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".pdf"), plot, height = 10L, width = 10L)
              })
              return()
            }))
          }
        }
      }
    }
  }
}
if (length(protlspep)) {
  source(parSrc, local = FALSE)
  dir <- paste0(wd, "/Heatmaps")
  pepHtmp(protlspep,
          pep,
          pep.ref[length(pep.ref)],
          dir,
          cl = parClust)
}
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))

#### Code chunk - peptide tables for visualizing the coverage of proteins of interest in 3D using SCV
# Edit me to also just use the newer cov3D function!!!
if (length(protlspep)) {
  # From https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r
  valid_url <- function(url_in, t = 2){
    con <- url(url_in)
    check <- suppressWarnings(try(open.connection(con, open = "rt", timeout = t), silent = TRUE)[1L])
    suppressWarnings(try(close.connection(con), silent = TRUE))
    ifelse(is.null(check), TRUE, FALSE)
  }
  # Mods and their mass shifts
  SCV_PTMs <- TRUE
  if (!"Mass shift" %in% colnames(Modifs)) {
    if ("UniMod" %in% colnames(Modifs)) {
      if (!require("PTMods", quietly = TRUE)) { pak::pak("rformassspectrometry/PTMods") }
      require(PTMods)
      data(modifications, package = "PTMods")
      UniMod <- modifications
      Modifs$"Mass shift" <- UniMod$MonoMass[match(Modifs$UniMod, UniMod$UnimodId)]
    } else {
      if ("MAXQUANT" %in% SearchSoft) {
        if ((length(MQFold) == 1L)&&(dir.exists(MQFold))) {
          modFls <- paste0(MQFold, "/bin/conf/modifications", c("", ".local"), ".xml")
          modFls <- modFls[which(file.exists(modFls))]
        } else {
          dflt <- if ((exists("mqFld"))&&(dir.exists(mqFld))) { mqFld } else { "C:" }
          dflt <- paste0(dflt, "/*.xml")
          #modFls <- choose.files(dflt, "Select MaxQuant modifications file(s) as source of PTMs mass shifts")
          #if ((length(modFls) > 1)||(!is.na(modFls))) { mqFld <- unique(dirname(modFls)) }
          modFls <- rstudioapi::selectFile("Select a MaxQuant modifications file as source of PTMs mass shifts",
                                           path = dflt,
                                           filter = "XML file (*.xml)")
          modFls <- unique(c(modFls, rstudioapi::selectFile("Optionally also select a local MaxQuant modifications file",
                                                            path = dflt,
                                                            filter = "XML file (*.xml)")))
        }
        modFls <- modFls[which(!is.na(modFls))]
        if (length(modFls)) {
          cran_req <- unique(c(cran_req, "xml2"))
          if (!require("xml2", quietly = TRUE)) { install.packages(xml2) }
          require(xml2)
          modFls <- lapply(modFls, \(modFl) { #modFl <- modFls[1L]
            xml_lst <- as_list(read_xml(modFl))
            xml_lst <- xml_lst[[1L]]
            xml_lst <- as.data.frame(t(sapply(xml_lst, \(x) {
              #x <- xml_lst[[1L]]
              return(c(attr(x, "title"), attr(x, "composition")))
            })))
          })
          modFls <- plyr::rbind.fill(modFls)
          colnames(modFls) <- c("Name", "Composition")
          Modifs$Composition <- modFls$Composition[match(Modifs$"Full name", modFls$Name)]
          w <- which(is.na(Modifs$Composition))
          if (length(w)) {
            for (i in w) {
              msg <- paste0("Enter the chemical formula of modification ", Modifs$"Full name"[i], " (Example: enter \"H(-1) N(-1) O\" for \"deamidation (NQ)\")")
              Modifs$Composition[i] <- dlg_input(msg, "")$res
            }
            w <- which(vapply(colnames(Modifs), \(x) { is.list(Modifs[[x]]) }, 1L))
            for (i in w) { temp[[i]] <- vapply(temp[[i]], paste, 1L, collapse = ", ") }
            write.csv(temp, "Workflow control/Modifications.csv", row.names = FALSE)
          }
          Modifs$"Mass shift" <- vapply(strsplit(Modifs$Composition, " "), \(x) {
            #x <- strsplit(Modifs$Composition, " ")[1L]
            x <- unlist(x)
            x <- as.data.frame(t(sapply(strsplit(gsub("\\)$", "", x), "\\("), \(y) {
              if (length(y) == 1L) { y <- c(y, 1L) }
              return(y)
            })))
            m <- match(x[[1L]], IsotopeProbs$Atom)
            stopifnot(sum(is.na(m)) == 0L)
            # For now the code above throws an error if an elements is missing from the table
            # If it ever does, I should expand the table to add isotopic probabilities for more elements!!!
            x <- sum(as.numeric(gsub("_.+", "", IsotopeProbs$Monoisotopic[m]))*as.integer(x[[2L]]))
            return(x)
          }, 1L)
        } else {
          warning("Could not map PTMs to mass shifts, these will be ignored from the SCV visualisations.")
          SCV_PTMs <- FALSE
        }
      } else {
        warning("Could not map PTMs to mass shifts, these will be ignored from the SCV visualisations.")
        SCV_PTMs <- FALSE
      }
    }
  }
  g1 <- grsep2(protlspep, pep$Proteins)
  prVect <- pep$Proteins[g1]
  modSq <- pep$"Modified sequence"[g1]
  PDB_in_DB <- ("PDB" %in% colnames(db))
  g2 <- grsep2(protlspep, db$`Protein ID`)
  dbPDB <- db$PDB[g2]
  dbPID <- db$"Protein ID"[g2]
  source(parSrc, local = FALSE)
  clusterExport(parClust,
                list("prVect", "wd", "PDB_in_DB", "SCV_PTMs", "dbPDB", "dbPID", "Modifs", "valid_url", "dirlist", "modSq"),
                envir = environment())
  #for (plp in protlspep) { #plp <- protlspep[1L]
  Tst <- parSapply(parClust, protlspep, \(plp) { #plp <- protlspep[1L]
    grs <- grsep2(plp, prVect)
    OutCome <- FALSE
    if (length(grs)) {
      dir <- paste0(wd, "/Coverage/", plp)
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      # For each protein we want to download:
      # - all models for all fragments.
      # - latest version only!
      # PDB: we get PDB IDs from parsing the txt file
      kPBD <- 0L
      if (PDB_in_DB) {
        tmp <- unlist(unlist(strsplit(dbPDB[match(plp, dbPID)], ";")))
        tmp <- tmp[which(tmp != "")]
        kPBD <- length(tmp)
        if (kPBD) {
          for (i in tmp) {
            url <- paste0("https://files.rcsb.org/download/", i, ".pdb")
            kPBD <- valid_url(url)
            if (kPBD) { download.file(url, paste0(dir, "/PDB-", i, ".pdb")) }
          }
        }
      }
      # Then AlphaFold
      kAlpha <- 0L
      tstF <- TRUE # Continue looking for the next fragment?
      while (tstF) {
        kAlpha <- kAlpha+1L
        kV <- 0L
        while ((!kV)||(tstV)) {
          kV <- kV + 1L
          mdlNm <- paste0("AF-", plp,"-F", kAlpha, "-model_v", kV, ".pdb")
          url <- paste0("https://alphafold.ebi.ac.uk/files/", mdlNm)
          tstV <- valid_url(url) # We want to find out which is the latest v version of a model for that protein
        }
        kV <- kV - 1L # The last is always a failure
        if (kV) { # Did we find a valid url?
          mdlNm <- paste0("AF-", plp,"-F", kAlpha, "-model_v", kV, ".pdb")
          url <- paste0("https://alphafold.ebi.ac.uk/files/", mdlNm)
          download.file(url, paste0(dir, "/", mdlNm))
        }
        tstF <- kV > 0L
      }
      kAlpha <- kAlpha-1L # The last is always a failure
      if (kPBD||kAlpha) {
        # We have found at least one model model which can be used to visualize coverage
        # Let's write peptidoforms
        tmp <- modSq[grs]
        if (SCV_PTMs) {
          tmp <- gsub("_", "", tmp)
          tmp <- gsub("\\)", "]_",gsub("\\(", "_[", tmp))
          tmp <- strsplit(tmp, "_")
          tmp <- vapply(tmp, \(x) { #x <- tmp[1L]
            x <- unlist(x)
            w <- grep("\\[.+\\]", x)
            x[w] <- paste0("[", round(Modifs$"Mass shift"[match(x[w], paste0("[", Modifs$Mark, "]"))], 0L), "]")
            return(paste(x, collapse = ""))
          }, "")
        } else { tmp <- gsub("[^A-Z]", "", tmp) }
        write(tmp, paste0(dir, "/SCV - observed peptides.txt"))
        OutCome <- TRUE
      }
    }
    return(OutCome)
  })
  # If this worked, we write a guide in the Coverage folder
  if (sum(Tst)) {
    Guide <- c("Visualing protein coverage in 3D using SCV",
               "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
               "",
               "If, for protein accession of interest \"PROTEIN\", a 3D model is available, then the following files are created in its subfolder:",
               " - \"SCV - expected peptides.txt\": sequences of identified peptides",
               " - For each fragment number i for which a structure is available (ranging from 1 to n), \"AF-PROTEIN-F1-model_v#.pdb\", where \"#\" is the latest valid version.",
               "",
               "To visualise the peptides onto the folded protein structure, navigate to https://scv.lab.gy and:",
               " - paste the peptides into the \"PSM/peptide list\" field",
               " - load the pdb file",
               " - Have fun!",
               "")
    write(Guide, paste0(wd, "/Coverage/SCV - how to visualise protein coverage in 3D.txt"))
  }
}

#### Code chunk - Create STRINGdb graph(s) and Cytoscape networks for regulated proteins
strngSrc <- paste0(libPath, "/extdata/Sources/STRINGdb.R")
#rstudioapi::documentOpen(strngSrc)
source(strngSrc, local = FALSE)

#### Code chunk - For pull-downs: create table summarizing types of evidence for all proteins of interest
# if (IsPullDown) {
#   g <- paste0("Regulated - ", unique(Exp.map[which(!Exp.map$Reference), VPAL$column]))
#   test <- apply(PG[, g, drop = FALSE], 1L, \(x) {
#     length(which(!x %in% c("", NA, "NA", "non significant", "too small FC")))
#   })
#   prot <- unique(c(prot.list, unlist(strsplit(PG$"Leading protein IDs"[which(test > 0)], ";"))))
#   if ((!is.null(prot.list_pep))&&(length(prot.list_pep))) { prot <- unique(c(prot, prot.list_pep)) }
#   if (length(prot)) {
#     dir <- paste0(wd, "/Evidences type tables")
#     if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
#     dirlist <- unique(c(dirlist, dir))
#     for (i in prot) {
#       temp <- ev[grsep2(i, ev$"Leading proteins"),]
#       temp <- sapply(unique(ev$Type), \(x) {
#         vapply(unique(ev$"Raw file"), \(y) { length(which((temp$Type == x)&(temp$"Raw file" == y))) }, 1L)
#       })
#       write.csv(temp, paste0(dir, "/Ev table - ", i, ".csv"))
#     }
#   }
# }

#### Code chunk - Finalize analysis and export results
# Remove empty directories:
#dirlist <- list.dirs()
dirlist <- as.character(unlist(dirlist))
dirlist <- dirlist[order(nchar(dirlist), decreasing = TRUE)]
for (dir in dirlist) { #d <- dirlist[1L]
  if (!length(list.files(dir))) {
    unlink(dir, recursive = TRUE)
    dirlist <- dirlist[which(dirlist != dir)]
  }
}
# Save decisions - should this go?
save(AllAnsw, file = "All_decisions.RData")

# Finalize reports
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

# Visualize results - already run earlier, could be expanded here to add more visualizations at further stages
#rstudioapi::documentOpen(xplorSrc)
#source(xplorSrc, local = FALSE)
