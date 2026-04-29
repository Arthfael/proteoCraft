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
load_a_Bckp %<o% c(TRUE, FALSE)[match(svDialogs::dlg_message("Do you want to load a backup?", "yesno")$res, c("yes", "no"))]
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
###-|-### Workflows: setNames(c("Differential Protein Expression analysis", "Pull-Down (e.g. co-IP)", "Biotin-based Pull-Down (BioID, TurboID, APEX...)", "Time Course","SubCellular Localisation analysis"), c("REGULATION", "PULLDOWN", "BIOID", "TIMECOURSE", "LOCALISATION"))
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
ReportCalls$Calls <- AddTxt2Report("Processing PSMs...")
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
    msg <- paste0("Removing ", l, " (", signif(100*l/nrow(ev), 2L), "%) PSMs with invalid expression values!")
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
    msg <- paste0("Removing ", l, " (", signif(100*l/nrow(ev), 2L), "%) PSMs with invalid expression values!")
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

useSAM %<o% ((names(pvalue.col)[which(pvalue.use)] == "Student")&&(useSAM_thresh))

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

#### Code chunk - Perform statistical tests
dataType <- "PG"
Src <- paste0(libPath, "/extdata/Sources/Stat_tests.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

### Visualize and check P-values
# Needs rewriting to accept peptide-level input
# Src <- paste0(libPath, "/extdata/Sources/pVal_check.R")
# #rstudioapi::documentOpen(Src)
# source(Src, local = FALSE)

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
w <- which(vapply(colnames(ev), \(x) { is.list(ev[[x]]) }, TRUE))
if (length(w)) { for (i in w) { ev[[i]] <- parSapply(parClust, ev[[i]], paste, collapse = ";") } }
data.table::fwrite(ev, paste0(dir, "/evidence.tsv"), sep = "\t", row.names = FALSE, na = "NA")
#
## Main peptidoforms-level, multi-tabs report
# Create openxlsx2 styles
#   It may make sense from the way Excel works, but I HATE how openxlsx2 deals with styles!
#   Anyway... 
#   So. We will. CHEAT!
#   I have saved a dummy tab with my old openxlsx styles,
#   which I will load in openxlsx2 to get and copy the styles from.
intNms <- function(nms, topLvl = FALSE, type = "PG") {
  m <- match(type, c("pep", "PG"))
  root <- c("Intensity", "Expression")[m]
  mode <- topLvl+1L
  sapply(nms, \(nm) {
    if (nm %in% c("Original", "Intensity", "Expression", "Original int.", "Intensity int.", "Expression int.")) {
      nm <- c(c("int.", "expr.")[m], root)[mode]
    } else {
      if (nm %in% c("ReNorm.", "ReNorm. int.")) { nm <- "re-norm" } else { nm <- substr(nm, 1L, min(c(3L, nchar(nm)))) }
      nm <- paste0(nm, ". ", c(c("int.", "expr.")[m], root)[mode])
    }
    paste0("log10(", nm, ")")
  })
}
ratNms <- function(nms, topLvl = FALSE) {
  mode <- topLvl+1L
  sapply(nms, \(nm) {
    if (nm %in% c("Original", "Ratios", "Original rat.", "Ratios rat.")) {
      nm <- c("rat.", "Ratio")[mode]
    } else {
      if (nm %in% c("ReNorm.", "ReNorm. rat.")) { nm <- "Re-norm" } else { nm <- substr(nm, 1L, min(c(3L, nchar(nm)))) }
      nm <- paste0(nm, ". ", c("rat.", "ratios")[mode])
    }
    paste0("log2(", nm, ")")
  })
}
for (nm in names(pep.ref)) { #nm <- names(pep.ref[1L])
  rpl <- intNms(nm, type = "pep")
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
fl <- system.file("extdata", "Report - column names - with replicates.xlsx", package = "proteoCraft")
styleNms <- openxlsx2::read_xlsx(fl, "tmp", colNames = FALSE)[, 1L]
WorkBook %<o% wb_load(fl)
repFl <- paste0(wd, "/Tables/Report_", dtstNm, ".xlsx")
WorkBook <- wb_add_data(WorkBook, "Description", dtstNm, wb_dims(2L, 5L))
WorkBook <- wb_add_data(WorkBook, "Description", format(Sys.Date(), "%d/%m/%Y"), wb_dims(3L, 5L))
WorkBook <- wb_add_data(WorkBook, "Description", WhoAmI, wb_dims(4L, 5L))
tmp <- loadedPackages(TRUE)
WorkBook <- wb_add_data(WorkBook, "Description", tmp$Version[grep("proteoCraft", tmp$Name)], wb_dims(5L, 5L))
WorkBook <- wb_set_base_font(WorkBook, 11L, font_name = "Calibri")
cat(" - Writing Excel report...\n")
#
# Function for editing our header
KolEdit <- function(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl) {
  #KolNames <- ColumnsTbl$Col; intTbl = intColsTbl; ratTbl = ratColsTbl
  klnms <- KolNames
  KolNames <- gsub("Peptides?", "Pep.", KolNames)
  KolNames <- gsub("Evidences?", c("PSMs", "Ev.")[(SearchSoft == "MAXQUANT")+1L], KolNames)
  KolNames <- gsub("Spectr((al)|(um))", "Spec.", KolNames)
  KolNames <- gsub("Razor", "Raz.", KolNames)
  KolNames <- gsub("Unique", "Uniq.", KolNames)
  KolNames <- gsub("MS/MS", "MS2", KolNames)
  # This would be the place to edit PER sample evidence and MS/MS columns,
  # currently not needed because those do not exist yet
  for (nm in names(intTbl)) { #nm <- names(intTbl)[1L] #nm <- names(intTbl)[2]
    m <- match(intTbl[[nm]]$Log, KolNames)
    w <- which(!is.na(m))
    if (length(w)) {
      rpl <- intNms(nm, type = "pep")
      KolNames[m[w]] <- paste0(rpl, " ", intTbl[[nm]]$Sample[w])
    }
  }
  if (!missing("ratTbl")) { # Actually, should never be missing in this workflow!
    for (nm in names(ratTbl)) {
      m <- match(ratTbl[[nm]]$Log, KolNames)
      w <- which(!is.na(m))
      if (length(w)) {
        rpl <- ratNms(nm)
        KolNames[m[w]] <- paste0(rpl, " ", ratTbl[[nm]]$Sample[w])
      }
    }
  }
  wNF <- grep("^mod\\. F-test ", KolNames, invert = TRUE)
  if (F.test) {
    wF <- grep("^mod\\. F-test ", KolNames)
    KolNames[wF] <- gsub("^mod\\. F-test +", "F-test ", KolNames[wF]) # Shorter F-test tag
    KolNames[wF] <- gsub(" +-log10\\(Pvalue\\)( - )?", " -log10 pval. ", KolNames[wF])
    KolNames[wF] <- gsub(" +Regulated - ", " reg. ", KolNames[wF])
    KolNames[wF] <- gsub(" +Significant-", " signif. ", KolNames[wF])
  }
  KolNames[wNF] <- gsub(".*-log10\\(Pvalue\\)( - )?", "-log10 pval. ", KolNames[wNF])
  KolNames[wNF] <- gsub(".*Significant-", "signif. ", KolNames[wNF])
  KolNames[wNF] <- gsub(".*Regulated - ", "reg. ", KolNames[wNF])
  KolNames[wNF] <- gsub(".*Significant-", "signif. ", KolNames[wNF])
  KntKol <- paste0(AA, " Count")
  KolNames[which(KolNames %in% KntKol)] <- gsub(" Count$", "", KolNames[which(KolNames %in% KntKol)])
  #
  KolNames <- gsub("( - )|(___)", " ", KolNames)
  #
  # F-test
  g <- grep("F-test: ", klnms)
  if (length(g)) {
    KolNames[g] <- paste0("F-test ", KolNames[g])
    g <- grep("F-test .*F(-| |\\.|_)?test", KolNames)
    stopifnot(length(g) == 0L)
  }
  #
  # Those names must be unique if the data is to be written as a table!
  # Which is annoying, because this limits how much fat we can cut
  tst <- aggregate(KolNames, list(KolNames), c)
  tst$L <- lengths(tst$x)
  tst <- tst[which(tst$Group.1 != ""),]
  stopifnot(max(tst$L) == 1L)
  #tst$x[which(tst$L > 1L)]
  #
  KolNames <- as.data.frame(t(KolNames))
  colnames(KolNames) <- klnms
  return(KolNames)
}
#KolEdit(xlTabs[[sheetnm]], intColsTbl, ratColsTbl)
if ((prot.list.Cond)&&(!"In list" %in% colnames(ev))) {
  g <- grsep2(prot.list, ev$Proteins)
  w <- rep(FALSE, nrow(ev))
  w[g] <- TRUE
  ev$"In list" <- c("", "+")[w+1L]
}
QualFilt %<o% c("In list", "Potential contaminant", "Only identified by site")
if ((DiscFilt)&&(DiscFiltMode == DiscFiltModes[3])) { QualFilt <- c(QualFilt, DiscFiltCols) }
II <- setNames(1, "All peptidoforms")
if ((exists("PTMs_pep"))&&(length(PTMs_pep))) {
  stop("Really? There is no PTM-modified class of peptides to analyze? Why did you run this workflow then?")
} else {
  Mod2Write <- names(PTMs_pep)
  II[paste0(Mod2Write, "-mod. pept.")] <- 1L+(seq_along(length(Mod2Write)))
}
for (ii in II) { #ii <- II[1L] #ii <- II[2]
  tblMode <- tblMode2 <- "pep"
  TbNm <- names(II)[ii]
  tempData <- get(tblMode)
  intRf <- pep.ref
  ratRf <- pep.ratios.ref
  if (ii > 1L) {
    Ptm <- Mod2Write[[ii-1L]]
    ptm <- Modifs$Mark[match(Ptm, Modifs$`Full name`)]
    tempData <- PTMs_pep[[Ptm]]
    tempData$Name <- gsub("\n", " ", tempData$Name)
    intRf <- PTMs_int.ref[[Ptm]]
    ratRf <- PTMs_rat.ref[[Ptm]]
    tblMode2 <- paste0(Ptm, "-modified pep")
  }
  names(intRf) <- paste0(names(intRf), " int.")
  names(ratRf) <- paste0(names(ratRf), " rat.")
  if (nrow(tempData)) {
    CoreCol %<o% c("id", "Modified sequence")
    if ("Modified sequence_verbose" %in% colnames(tempData)) { CoreCol <- c(CoreCol, "Modified sequence_verbose") }
    CoreCol <- c(CoreCol, "Sequence", "Proteins")
    if (ii > 1L) { CoreCol <- c(CoreCol, "Name", paste0(Ptm, "-site(s)")) }
    CoreCol2 %<o% c("Leading razor proteins",
                    "Protein names", "Gene names", "Protein group IDs", "Razor protein group ID",
                    grep(" (Probabilities|Score Diffs)$", colnames(tempData), value = TRUE),
                    "Normalisation group")
    evcol <- "Evidence IDs"
    spcol <- "MS/MS count"
    gel <- setNames(lapply(intRf, \(rf) {
      x <- c(paste0(rf, RSA$values),
             paste0("Mean ", rf, VPAL$values))
      return(x[which(x %in% colnames(tempData))])
    }), intRf)
    if (ii == 1) {
      # Log transform for normal tables - not necessary for PTM-modified tables as we already transformed
      gel2 <- setNames(lapply(intRf, \(rf) {
        x <- c(paste0(rf, RSA$values),
               paste0("Mean ", rf, VPAL$values))
        w <- which(x %in% colnames(tempData))
        rf2 <- sub("Int\\. - ", "log10(Int.) - ", rf)
        x <- c(paste0(rf2, RSA$values),
               paste0("Mean ", rf, VPAL$values))
        return(x[w])
      }), intRf)
      for (rf in intRf) { #rf <- intRf[1L]
        kol1 <- gel[[rf]]
        kol2 <- gel2[[rf]]
        if (length(kol1)) {
          tempData[, kol2] <- suppressWarnings(log10(tempData[, kol1]))
          for (kl in kol2) {
            w <- which(is.infinite(tempData[[kl]]))
            tempData[w, kl] <- NA
          }
        }
      }
      intRf <- sapply(intRf, \(rf) {
        sub("Int\\. - ", "log10(Int.) - ", rf)
      })
      gel <- unlist(gel2)
    }
    #
    smpls <- c(VPAL$values, RSA$values)
    if (length(Exp) == 1L) { smpls <- gsub(topattern(paste0(Exp, "___")), "", smpls) }
    smpls <- gsub("___", " ", smpls)
    intColsTbl <- setNames(lapply(names(intRf), \(nm) { #nm <- names(intRf)[1L]
      res <- data.frame(Log = c(paste0("Mean ", intRf[nm], VPAL$values),
                                paste0(intRf[nm], RSA$values)),
                        Type = c(rep("Average", length(VPAL$values)),
                                 rep("Individual", length(RSA$values))),
                        Sample = smpls)
      w <- which(res$Log %in% colnames(tempData))
      return(res[w,])
    }), names(intRf))
    w <- which(vapply(intColsTbl, nrow, 1L) > 0L)
    intColsTbl <- intColsTbl[w]; intRf <- intRf[w]
    quantCols <- intCols <- lapply(intColsTbl, \(x) { x$Log })
    ratColsTbl <- setNames(lapply(names(ratRf), \(nm) {
      res <- data.frame(Log = c(paste0("Mean ", ratRf[nm], VPAL$values),
                                paste0(ratRf[nm], RSA$values)),
                        Type = c(rep("Average", length(VPAL$values)),
                                 rep("Individual", length(RSA$values))),
                        Sample = smpls)
      w <- which(res$Log %in% colnames(tempData))
      return(res[w,])
    }), names(ratRf))
    w <- which(vapply(ratColsTbl, nrow, 1L) > 0L)
    ratColsTbl <- ratColsTbl[w]; ratRf <- ratRf[w]
    ratCols <- lapply(ratColsTbl, \(x) { x$Log })
    grl <- unlist(ratCols)
    for (gr in grl) {
      w <- which(is.infinite(tempData[[gr]]))
      tempData[w, gr] <- NA
    }
    quantCols[names(ratRf)] <- ratCols
    regcol <- grep("^((Enriched)|(Regulated)) - ", colnames(tempData), value = TRUE)
    signcol <- grep("^Significant-FDR=[1-9][0-9]*\\.*[0-9]*% - ", colnames(tempData), value = TRUE)
    signcol <- grep(" - Analysis_[0-9]+", signcol, invert = TRUE, value = TRUE)
    quantcol <- unlist(quantCols)
    PepColList %<o% c("gel", "grl", "quantcol", "signcol", "regcol") # These are any column for which we want to gsub "___" to " "
    .obj <- unique(c(PepColList, .obj)) # Here easier than using a custom operator
    if (ii > 1L) {
      gpl <- grep(topattern(pvalue.col[which(pvalue.use)]), colnames(tempData), value = TRUE)
      quantcol <- c(quantcol, gpl)
      PepColList <- c(PepColList, "gpl")
    }
    aacol <- paste0(AA, " Count")
    qualFlt <- QualFilt[which(QualFilt %in% colnames(ev))]
    w <- which(!qualFlt %in% colnames(tempData))
    if (length(w)) {
      tempData[, qualFlt[w]] <- ev[match(tempData$"Modified sequence", ev$"Modified sequence"), qualFlt[w]]
    }
    kol <- c(CoreCol, CoreCol2, evcol, spcol, "PEP", quantcol, signcol, regcol, qualFlt, aacol)
    if (ii > 1L) { kol <- c(kol, "Code") }
    if (Annotate) {
      PepAnnotCol %<o% annot.col[which(annot.col %in% colnames(tempData))]
      kol <- c(kol, PepAnnotCol)
    }
    #tst <- data.frame(Names = names(kol), Column = setNames(kol, NULL), Found = kol %in% colnames(tempData));View(tst)
    kol <- kol[which(kol %in% colnames(tempData))]
    #kol[which(!kol %in% colnames(tempData))]
    #colnames(tempData)[which(!colnames(tempData) %in% kol)]
    tempData <- tempData[, kol]
    # If there is only one experiment, remove it from the names here...
    colnames(tempData) <- cleanNms(colnames(tempData), start = FALSE)
    for (i in PepColList) { assign(i, cleanNms(get(i), start = FALSE)) }
    intColsTbl <- lapply(intColsTbl, \(x) {
      x$Log <- cleanNms(x$Log, start = FALSE)
      x
    })
    ratColsTbl <- lapply(ratColsTbl, \(x) {
      x$Log <- cleanNms(x$Log, start = FALSE)
      x
    })
    intCols <- lapply(intCols, cleanNms, start = FALSE)
    ratCols <- lapply(ratCols, cleanNms, start = FALSE)
    quantCols <- lapply(quantCols, cleanNms, start = FALSE)
    colnames(tempData) <- gsub("_names$", " names", colnames(tempData))
    if (Annotate) { PepAnnotCol <- gsub("_names$", " names", PepAnnotCol) }
    for (k in regcol) {
      tempData[which(tempData[[k]] == "non significant"), k] <- "n.s."
      tempData[which(tempData[[k]] == ""), k] <- "n.t."
    }
    if ((ii > 1L)&&(F.test)) {
      tempPepF <- PTMs_F_test_data[[Ptm]]
      tempPepF <- tempPepF[, which(!colnames(tempPepF) %in% c(Param$Plot.labels, "Rel. log10(Peptides count)", "Av. log10 abundance"))]
      colnames(tempPepF) <- cleanNms(colnames(tempPepF), start = FALSE)
      mnratcolF %<o% grep("Mean log2\\(Ratio\\) - ", colnames(tempPepF), value = TRUE)
      m <- match(mnratcolF, colnames(tempPepF))
      colnames(tempPepF)[m] <- paste0("mod. F-test ", mnratcolF)
      mnratcolF <- paste0("mod. F-test ", mnratcolF)
      pvalcolF %<o% F_Root
      signcolF %<o% grep("^mod\\. F-test Significant", colnames(tempPepF), value = TRUE)
      regcolF %<o% grep("^mod\\. F-test Regulated", colnames(tempPepF), value = TRUE)
      Fkol %<o% c(regcolF, mnratcolF, pvalcolF, signcolF)
      for (k in regcolF) {
        tempPepF[which(tempPepF[[k]] == "non significant"), k] <- "n.s."
        tempPepF[which(tempPepF[[k]] == ""), k] <- "n.t."
      }
      tempData[, Fkol] <- tempPepF[match(tempData$"Modified sequence", tempPepF$"Modified sequence"), Fkol]
    }
    dir <- paste0(wd, "/Tables")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    data.table::fwrite(tempData, paste0(dir, "/", TbNm, ".tsv"), sep = "\t", row.names = FALSE, na = "NA")
    w <- grsep2(prot.list, tempData$Proteins)
    if (length(w)) {
      data.table::fwrite(tempData[w,], paste0(wd, "/Tables/", TbNm, " - Proteins in list.tsv"),
                         sep = "\t", row.names = FALSE, na = "NA")
    }
    # Which columns are affected by each style
    # - IDs
    ColumnsTbl %<o% list(IDs = c(CoreCol, CoreCol2))
    # - Counts
    ColumnsTbl[["AA counts"]] <- aacol
    # - Evidence counts and IDs
    ColumnsTbl[["Global Ev. IDs"]] <- "Evidence IDs"
    ColumnsTbl[["Global Spec. counts"]] <- "MS/MS count"
    # - Individual Expr
    ColumnsTbl[["Individual Expr"]] <- grep("^Mean ", unlist(intCols), value = TRUE, invert = TRUE)
    # - Summary Expr
    ColumnsTbl[["Summary Expr"]] <- grep("^Mean ", unlist(intCols), value = TRUE)
    # - Individual Ratios
    ColumnsTbl[["Individual Ratios"]] <- grep("^Mean ", unlist(ratCols), value = TRUE, invert = TRUE)
    # - Summary Ratios
    ColumnsTbl[["Summary Ratios"]] <- grep("^Mean ", unlist(ratCols), value = TRUE)
    if (ii > 1L) {
      # - Summary Ratios: P-values and significance
      ColumnsTbl[["P-values"]] <- gpl
      # - Significant
      ColumnsTbl[["Significant"]] <- signcol
      # - Regulated
      ColumnsTbl[["Regulated"]] <- regcol
      # F-test
      if ((ii > 1L)&&(F.test)) {
        ColumnsTbl[["F-test summary Ratios"]] <- mnratcolF
        ColumnsTbl[["F-test P-values"]] <- pvalcolF
        ColumnsTbl[["F-test significant"]] <- signcolF
        ColumnsTbl[["F-test regulated"]] <- regcolF
      }
    }
    # - Annotations
    if (Annotate) {
      AnnotTbl$Columns <- list(c("GO", "GO-ID"), c("Taxonomy", "TaxID"), NA, NA, NA, NA, "EMBL", NA)
      for (i in annot) { AnnotTbl$Columns[match(i, AnnotTbl$Name)] <- list(c(i, paste0(i, " names"))) }
      AnnotTbl$Columns[match("Other", AnnotTbl$Name)] <- list(annot.col2[which(!annot.col2 %in% unlist(AnnotTbl$Columns))])
      for (i in 1L:nrow(AnnotTbl)) { ColumnsTbl[[paste0(AnnotTbl$Name[i], " annotations")]] <- AnnotTbl$Columns[[i]] }
    }
    # - PEP
    ColumnsTbl[["PEP"]] <- "PEP"
    # - Filters
    ColumnsTbl[["Filters"]] <- qualFlt
    # Melt
    ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, \(x) { length(x[which(!is.na(x))]) }, 1L) > 0L)]
    ColumnsTbl <- set_colnames(reshape::melt.list(ColumnsTbl), c("Col", "Grp"))
    #tst <- aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), length); View(tst)
    #tst <- aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), unique); View(tst)
    #tst <- aggregate(1L:nrow(ColumnsTbl), list(ColumnsTbl$Col), unique); w <- which(lengths(tst$x) > 1L); setNames(tst$x[w], tst$Group.1[w])
    stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
    ColumnsTbl$Class <- ""
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Peptides information"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evcol))] <- "Evidence IDs"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(spcol))] <- "Spectral count"
    for (nm in names(intRf)) { #nm <- names(intRf)[1L]
      rpl <- intNms(nm, TRUE, type = "pep")
      ColumnsTbl$Class[grep(topattern(intRf[nm]), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[grep(topattern(paste0("Mean ", intRf[nm])), ColumnsTbl$Col)] <- rpl
    }
    for (nm in names(ratRf)) { #nm <- names(ratRf)[1L]
      rpl <- ratNms(nm, TRUE)
      ColumnsTbl$Class[grep(topattern(ratRf[nm]), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[grep(topattern(paste0("Mean ", ratRf[nm])), ColumnsTbl$Col)] <- rpl
    }
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "P-values")] <- gsub(" - $", "", pvalue.col[which(pvalue.use)])
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% signcol)] <- "Significant"
    if ((ii > 1L)&&(F.test)) {
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test summary Ratios")] <- "F-test"
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test P-values")] <- "F-test"
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test significant")] <- "F-test"
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test regulated")] <- "F-test"
    }
    ColumnsTbl$Class[grep("[Aa]nnotations", ColumnsTbl$Grp)] <- "Annotations"
    ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters"))] <- "QC filters"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% aacol)] <- "Amino Acid counts"
    ColumnsTbl$Hide <- ColumnsTbl$Class %in% c("Spectral count", "Spectrum IDs", "Amino Acid counts", "Annotations", "Cluster (hierarch.)")
    #
    if (MakeRatios) { a <- KolEdit(ColumnsTbl$Col, intColsTbl, ratColsTbl) } else { a <- KolEdit(ColumnsTbl$Col, intColsTbl) }
    ColumnsTbl$edit_Col <- unlist(a)
    #
    Src <- paste0(libPath, "/extdata/Sources/fstWrite_Excel_core_script.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
}
#
Src <- paste0(libPath, "/extdata/Sources/Write_Excel_end_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#WorkBook$get_active_sheet()
#xl_open(repFl)

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)

#### Code chunk - Venn diagrams
Src <- paste0(libPath, "/extdata/Sources/Venn_diagrams.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

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
