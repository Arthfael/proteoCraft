#### Code chunk - Initialization
if (!interactive()) { stop("This script should only be run within an interactive R session!") }
options(stringsAsFactors = FALSE)
options(install.packages.compile.from.source = "never")
options(svDialogs.rstudio = TRUE)
#rm(list = ls()[which(!ls() %in% c("dtstNm", "wd", "inDirs", "outDir"))])

## The proteoCraft package can be re-installed at any time in the workflow (there is a specific script for this in the package's library folder),
## or just load it here:
if (exists(".obj")) { rm(".obj") }
#myPackNm %<o% "proteoCraft" # Bad idea
library(proteoCraft)
dirlist %<o% c() # This should go!!!
ReUseAnsw %<o% FALSE
scrptType %<o% "withReps"
scrptTypeFull %<o% "withReps_PTMs_only"

RPath %<o% as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath %<o% paste0(RPath, "/proteoCraft")
homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1, 1)) }
parSrc %<o% paste0(libPath, "/extdata/R scripts/Sources/make_check_Cluster.R")
fls <- paste0(homePath, "/", c(#"Regulation analysis - master script.R",
                               "Regulation analysis - detailed script.R",
                               "Regulation analysis - detailed script_pepOnly.R",
                               "No replicates analysis - detailed script.R",
                               "Reload_renv_from_lock_file.R",
                               "Default_locations.xlsx",
                               "LC_columns.xlsx"))
tst <- sum(!file.exists(fls))
if (tst) { proteoCraft::Configure() }
#xplorSrc %<o% paste0(libPath, "/extdata/R scripts/Sources/xplorData.R")

# Parameters used by the master script:
###-|-### Workflows: setNames(c("Differential Protein Expression analysis", "Pull-Down (e.g. co-IP)", "Biotin-based Pull-Down (BioID, TurboID, APEX...)", "Time Course","SubCellular Localisation analysis"), c("REGULATION", "PULLDOWN", "BIOID", "TIMECOURSE", "LOCALISATION"))
###-|-### Replicates? TRUE
###-|-### External dependencies: Excel (loose); ScanHeadsman (loose); Cytoscape (loose); saintExpress (auto)

### Packages
## For convenience all (or most) of the packages used are loaded or installed here:
## CRAN packages:
if(!exists("cran_req")) { cran_req %<o% "pak" } else { cran_req %<o% cran_req }
if(!exists("bioc_req")) { bioc_req %<o% c() } else { bioc_req %<o% bioc_req }
cran_req <- unique(c(cran_req, "pak", "fs", "shiny", "renv", "R.utils", "data.table", "devtools", "qs2", "shinyWidgets", "DT", "shinyBS", "stringr",
                     "gplots", "ggplot2", "ggpubr", "gtools", "reshape", "reshape2", "compiler", "stats", "rgl", "ggrepel", "rstudioapi", "modeest",
                     "minpack.lm", "snow", "viridis", "pcaMethods", "impute", "imputeLCMD", "parallel", "coin", "openxlsx", "openxlsx2", "plotly",
                     "Peptides", "xml2", "pdftools", "statmod", "ggpolypath", "venn", "gridExtra", "svDialogs", "htmlwidgets", "magrittr", "tibble",
                     "officer", "hexbin", "igraph", "matlib", "umap", "plyr", "ggnewscale", "shinyjs", "shinyFiles", "TeachingDemos", "shinycssloaders",
                     "tidyr", "ggplotify", "jpeg", "scattermore", "rpanel", "stringi", "lmtest", "ssh", "taxize", "arrow", "unimod"))
bioc_req <- unique(c(bioc_req, "biomaRt", "GO.db", "UniProt.ws", "limma", "sva", "qvalue", "MSnbase", "DEP",
                     "Rgraphviz", "RCy3", "siggenes", "DEqMS"))
inst <- as.data.frame(installed.packages())
for (pack in cran_req) {
  if (!pack %in% inst$Package) {
    if (pack %in% c("pak", #"shiny",
                    "uchardet", #"openxlsx2",
                    "taxize", "unimod")) {
      # Exceptions where for now we want a specific version to be installed,
      # or have to help the installer so it finds the right location
      if (pack == "pak") {
        install.packages("pak", dependencies = TRUE)
      }
      # if (pack == "shiny") { # Should be fixed now
      #   install.packages("https://cran.r-project.org/src/contrib/Archive/shiny/shiny_1.7.5.tar.gz", dependencies = TRUE)
      # }
      if (pack == "uchardet") {
        url <- "https://cran.r-project.org/src/contrib/Archive/uchardet/uchardet_1.1.1.tar.gz"
        destfile <- "uchardet_1.1.1.tar.gz"
        tst <- try(download.file(url, destfile, "curl"), silent = TRUE)
        if ("try-error" %in% class(tst)) { try(download.file(url, destfile, "wget"), silent = TRUE) }
        install.packages(destfile, dependencies = TRUE)
        unlink(destfile)
      }
      # if (pack == "openxlsx2") {
      #   pak::pkg_install("JanMarvin/openxlsx2@v1.10", ask = FALSE, upgrade = TRUE, dependencies = TRUE) # ... until I can figure out what is happening...
      # }
      # if (pack == "myTAI") {
      #   pak::pkg_install("drostlab/myTAI@v0.9.3", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
      # }
      if (pack == "taxize") {
        pak::pkg_install("ropensci/bold", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
        pak::pkg_install("ropensci/taxize", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
      }
      if (pack == "unimod") {
        pak::pkg_install("rformassspectrometry/unimod", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
      }
    } else {
      tst <- try(pak::pkg_install(pack, ask = FALSE, upgrade = TRUE, dependencies = TRUE), silent = TRUE)
      if ("try-error" %in% class(tst)) {
        tst <- try(install.packages(pack, dependencies = TRUE), silent = TRUE)
      }
      if ("try-error" %in% class(tst)) {
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
    tst <- try(pak::pkg_install(pack, ask = FALSE, upgrade = TRUE, dependencies = TRUE), silent = TRUE)
    if ("try-error"%in% class(tst)) {
      tst <- try(pak::pkg_install(pack, ask = FALSE, upgrade = TRUE, dependencies = FALSE), silent = TRUE)
    }
    if ("try-error"%in% class(tst)) {
      tst <- try(pak::pkg_install(pack, ask = FALSE, upgrade = FALSE, dependencies = FALSE), silent = TRUE)
    }
    if ("try-error"%in% class(tst)) {
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
    tmp <- openxlsx2::read_xlsx(paste0(homePath, "/Default_locations.xlsx"))
    load_Bckp(startDir = tmp$Path[which(tmp$Folder == "Temporary folder")])
  }, silent = TRUE)
}

# Set Shiny options, load functions for creating a Word report, create Excel styles
Src <- paste0(libPath, "/extdata/R scripts/Sources/ShinyOpt_Styles_and_Report.R")
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
Src <- paste0(libPath, "/extdata/R scripts/Sources/Start_analysis.R")
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
Src <- paste0(libPath, "/extdata/R scripts/Sources/Load_PSMs.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Install and/or load rawrr only if we have .raw files
g <- grep("\\.raw$", rawFiles, ignore.case = TRUE)
if (length(g)) {
  rawrrSrc <- paste0(libPath, "/extdata/R scripts/Sources/install_rawrr.R")
  #rstudioapi::documentOpen(rawrrSrc)
  source(rawrrSrc)
}

# MS raw files map
Src <- paste0(libPath, "/extdata/R scripts/Sources/Fractions_Map_editor.R")
#rstudioapi::documentOpen(Src)
tstFrMp <- FALSE
while (!tstFrMp) {
  source(Src, local = FALSE)
}

#### Code chunk - Edit Experimental Factors
Src <- paste0(libPath, "/extdata/R scripts/Sources/Experimental_Factors_editor.R")
#rstudioapi::documentOpen(Src)
tstXpFct <- FALSE
while (!tstXpFct) {
  source(Src, local = FALSE)
}
#

#### Code chunk - Edit Experiment map
Src <- paste0(libPath, "/extdata/R scripts/Sources/Experiment_Map_editor.R")
#rstudioapi::documentOpen(Src)
tstXpMp <- FALSE
while (!tstXpMp) {
  source(Src, local = FALSE)
}
#

#### Code chunk - Load and process search database(s)
Src <- paste0(libPath, "/extdata/R scripts/Sources/Process_Fasta_DBs.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#evNm %<o% c("PSM", "Evidence")[(SearchSoft == "MAXQUANT")+1]
evNm %<o% "PSM"

#### Code chunk - Load and process annotations
Src <- paste0(libPath, "/extdata/R scripts/Sources/Load_Annotations.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
source(parSrc, local = FALSE)
Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_prepare.R") # Doing this earlier but also keep latter instance for now
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Create experiment Factors shortcuts
Src <- paste0(libPath, "/extdata/R scripts/Sources/XpFact_shortcuts.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

# Check and process Fractions map
Src <- paste0(libPath, "/extdata/R scripts/Sources/Fractions_Map_check.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Define analysis parameters
# Define analysis parameters
paramSrc <- paste0(libPath, "/extdata/R scripts/Sources/rep_Parameters_editor_Main.R")
#rstudioapi::documentOpen(paramSrc)
source(paramSrc, local = FALSE)
# Temporary solution to the cross-app contamination issue: unload-reload packages
unloadNamespace("pRolocGUI")
unloadNamespace("colourpicker")
unloadNamespace("devtools")
unloadNamespace("miniUI")
unloadNamespace("shinyWidgets")
unloadNamespace("shinyFiles")
unloadNamespace("DEP")
unloadNamespace("shiny")
unloadNamespace("shinyjs")
unloadNamespace("shinyWidgets")
unloadNamespace("shinyhelper")
unloadNamespace("shinyFiles")
unloadNamespace("shinydashboardPlus")
unloadNamespace("shinydashboard")
unloadNamespace("shinycssloaders")
unloadNamespace("shinyBS")
unloadNamespace("DT")
library("pRolocGUI", character.only = TRUE)
library("colourpicker", character.only = TRUE)
library("devtools", character.only = TRUE)
library("miniUI", character.only = TRUE)
library("shinyWidgets", character.only = TRUE)
library("shinyFiles", character.only = TRUE)
library("DEP", character.only = TRUE)
library("shiny", character.only = TRUE)
library("shinyjs", character.only = TRUE)
library("shinyWidgets", character.only = TRUE)
library("shinyhelper", character.only = TRUE)
library("shinyFiles", character.only = TRUE)
library("shinydashboard", character.only = TRUE)
library("shinydashboardPlus", character.only = TRUE)
library("shinycssloaders", character.only = TRUE)
library("shinyBS", character.only = TRUE)
library("DT", character.only = TRUE)
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/rep_Parameters_editor_Stats.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Start writing Materials and Methods
Src <- paste0(libPath, "/extdata/R scripts/Sources/autoMatMet.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Start processing the PSMs table
ReportCalls <- AddSpace2Report()
ReportCalls$Calls <- AddTxt2Report("Processing PSMs...")
# Remove reverse database hits
ev <- ev[which(ev$Reverse == ""),]

# Optionally remove charge 1 PSMs - off for now, but may become either user decision or parameter controlled
RemovZ1 <- FALSE
w1 <- which(ev$Charge == 1)
wHt1 <- which(ev$Charge > 1)
if ((RemovZ1)&&(length(w1))) {
  AmIBogus <- paste(unique(ev$"Modified sequence"[w1]), collapse = "\n")
  #cat(AmIBogus)
  cat("Removing the following presumably bogus identifications with Z=1:\n", AmIBogus, "\n")
  ev <- ev[wHt1,]
}

w <- grep("CONTAMINANT", colnames(ev), ignore.case = TRUE)
if (length(w) > 1) { warning("Hmmm..., you might wanna check what is happening here...") } else {
  colnames(ev)[w] <- "Potential contaminant"
}
for (i in c("Potential contaminant", "Reverse")) {
  w <- which(is.na(ev[[i]]))
  ev[w, i] <- ""
}

#### Code chunk - Update peptide-to-protein mappings
Src <- paste0(libPath, "/extdata/R scripts/Sources/checkPep2Prot.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

# Filter to keep only PSMs with valid quantitative values:
if (LabelType == "LFQ") {
  source(parSrc, local = FALSE)
  if ((Param$Label == "DIA")&&("MS2 intensities" %in% colnames(ev))) {
    ev$MS2_intensities <- strsplit(ev$"MS2 intensities", ";")
    ev$MS2_intensities <- parLapply(parClust, ev$MS2_intensities, as.numeric) # (Let's keep this as a numeric list)
    temp <- ev[, c(ev.col["Original"], "MS2_intensities")]
    temp$SumS2 <- parSapply(parClust, temp$MS2_intensities, function(x) { sum(proteoCraft::is.all.good(x)) })
    # While we're at it, let's estimate missing MS1 intensities if we only have MS2:
    # (sum of MS2 intensities * median ratio of precursor intensity to sum of MS2 intensities)
    temp2 <- temp$Intensity/temp$SumS2
    m <- median(is.all.good(temp2))
    #sd(is.all.good(temp2))
    #plot <- ggplot(temp) + geom_point(aes(x = log10(Intensity), y = log10(SumS2))) + theme_bw() + geom_abline(intercept = log10(1/m), slope = 1, colour = "red")
    #poplot(plot)
    w <- which(((!is.all.good(ev[[ev.col["Original"]]], 2))|(ev[[ev.col["Original"]]] <= 0))&(temp$SumS2 > 0))
    if (length(w)) { ev[w, ev.col["Original"]] <- temp$SumS2[w]*m }
    test <- parApply(parClust, temp[, c(ev.col["Original"], "SumS2")], 1, sum)
  } else {
    temp <- ev[, ev.col["Original"], drop = FALSE]
    test <- parApply(parClust, temp, 1, function(x) { sum(proteoCraft::is.all.good(x)) })
  }
  l <- length(which(test == 0))
  if (l) {
    msg <- paste0("Removing ", l, " (", signif(100*l/nrow(ev), 2), "%) PSMs with invalid expression values!")
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
    kol <- lapply(tmpIso[w2], function(x) {
      paste0(c("Reporter intensity corrected ", "Reporter intensity ", "Reporter intensity count "), x)
    })
    w <- which(!colnames(ev) %in% unlist(kol))
    ev <- ev[, w]
  }
  assign(IsobarLab, tmpIso[w1])
  #
  kol <- paste0(ev.ref["Original"], get(IsobarLab))
  tst <- temp <- ev[, kol, drop = FALSE]
  tst$MS1 = ev[[ev.col["Original"]]]
  tst$Reporter = rowSums(temp, na.rm = TRUE)
  # Check dependency: there should be one in log space
  temp2 <- tst$MS1/tst$Reporter
  m <- median(is.all.good(temp2))
  #sd(is.all.good(temp2))
  #plot <- ggplot(tst) + geom_point(aes(x = log10(MS1), y = log10(Reporter))) + theme_bw() + geom_abline(intercept = log10(1/m), slope = 1, colour = "red")
  #poplot(plot)
  #
  #View()
  # If precursor intensity is missing, replace by estimate (sum of reporter intensities * median ratio of precursor intensity to sum of reporter intensities)
  w <- which((!is.all.good(ev[[ev.col["Original"]]], 2))|(ev[[ev.col["Original"]]] <= 0))
  if (length(w)) { ev[w, ev.col["Original"]] <- tst$Reporter[w]*m }
  # Now the reverse scenario: no reporters, but we have precursor intensities; these are throw-away stuff 
  w <- which(!is.all.good(tst$Reporter, 2)|(tst$Reporter <= 0))
  l <- length(w)
  if (l) {
    RemEv %<o% ev[w,]
    #View(RemEv[, kol])
    msg <- paste0("Removing ", l, " (", signif(100*l/nrow(ev), 2), "%) PSMs with invalid expression values!")
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
    w <- which(is.all.good(tst$Reporter, 2)&(tst$Reporter > 0))
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
Src <- paste0(libPath, "/extdata/R scripts/Sources/MA_plots.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Test for amino acid biases:
Src <- paste0(libPath, "/extdata/R scripts/Sources/AA_biases_test.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Isobaric label purity correction
Src <- paste0(libPath, "/extdata/R scripts/Sources/isobarCorr.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

ev$"Unique State" <- do.call(paste, c(ev[, c("Modified sequence", "Charge")], sep = ""))

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

# DIA-only: MS2-based correction of MS1-based quantitative values
Src <- paste0(libPath, "/extdata/R scripts/Sources/MS2corr2MS1.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Optional - Normalize PSM MS1 intensities, then, if applicable, MS2 reporter (Isobaric labelling) or fragment (DIA) intensities
Src <- paste0(libPath, "/extdata/R scripts/Sources/evNorm.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - ROC analysis for optimizing a threshold to include/exclude peptides mapped to specific GO terms
Src <- paste0(libPath, "/extdata/R scripts/Sources/ROC1.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Create modified peptides table
Src <- paste0(libPath, "/extdata/R scripts/Sources/pepMake.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Peptidoforms-level, calculate quantitation and test for outliers
## (future option, or when executing line-by-line: remove outliers/samples which will not be used)
## First visualize data: are there any clear outliers?
# Calculate single channel intensities and total intensity
#
# To do:
# - If outlier is the only reference sample in a reference group, unfortunately you will have to remove the whole group
# - Add Pearson correlation heatmap amongst visualisations to base decision to remove outliers, it is very good!
#
pep.ref %<o% setNames("Evidence intensities - ", "Original")
if (!"Use" %in% colnames(Exp.map)) { Exp.map$Use <- TRUE } else {
  if (class(Exp.map$Use) == "character") {
    Exp.map$Use[which(Exp.map$Use == "T")] <- "TRUE"
    Exp.map$Use[which(Exp.map$Use == "F")] <- "FALSE"
    Exp.map$Use <- as.logical(Exp.map$Use)
    Exp.map$Use[which(is.na(Exp.map$Use))] <- TRUE
  }
}
source(parSrc, local = FALSE)
exports <- list("smpls", "Exp.map", "pep.ref", "LabelType", "wd")
if (LabelType == "Isobaric") {
  tmp <- ev[, c("MQ.Exp", "Modified sequence",
                paste0(ev.ref[length(ev.ref)], as.character(sort(as.numeric(unique(Exp.map$"Isobaric label"))))))]
  exports <- append(exports, "ev.ref")
}
if (LabelType == "LFQ") {
  tmp <- ev[, c("MQ.Exp", "Modified sequence", ev.col[length(ev.col)])]
  exports <- append(exports, "ev.col")
}
readr::write_rds(tmp, paste0(wd, "/tmp.RDS"))
smpls <- unique(Exp.map$Ref.Sample.Aggregate[which(Exp.map$Use)])
clusterExport(parClust, exports, envir = environment())
invisible(clusterCall(parClust, function(x) {
  library(data.table)
  tmp <<- readr::read_rds(paste0(wd, "/tmp.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp.RDS"))
tmp4 <- setNames(parLapply(parClust, smpls, function(smpl) { #smpl <- smpls[1]
  m <- match(smpl, Exp.map$Ref.Sample.Aggregate)
  mqe <- unlist(Exp.map$MQ.Exp[m])
  w2 <- which(tmp$MQ.Exp %in% mqe)
  tmp2 <- data.frame(mod = NA, Intensity = NA)
  if (length(w2)) {
    if (LabelType == "Isobaric") {
      j <- as.character(sort(as.numeric(Exp.map$"Isobaric label"[m])))
      tmp3 <- tmp[w2, paste0(ev.ref[length(ev.ref)], j), drop = FALSE]
      for (k in j) {
        kk <- paste0(ev.ref[length(ev.ref)], j)
        tmp3[which(!proteoCraft::is.all.good(tmp3[[kk]], 2)), kk] <- NA
      }
      if (length(j) > 1) { tmp3 <- apply(tmp3, 1, sum, na.rm = TRUE) } # Ultra-rare cases where the same parent sample is in different isobaric channels in different fractions
      tmp2 <- data.table(mod = tmp$"Modified sequence"[w2],
                         Intensity = unlist(tmp3))
    }
    if (LabelType == "LFQ") {
      tmp2 <- data.table(mod = tmp$"Modified sequence"[w2],
                         Intensity = tmp[w2, ev.col[length(ev.col)]])
      tmp2$Intensity[which(!proteoCraft::is.all.good(tmp2$Intensity, 2))] <- NA
    }
    tmp2 <- tmp2[, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(mod)]
    tmp2 <- as.data.frame(tmp2)
  }
  return(tmp2)
}), smpls)
for (smpl in smpls) { #smpl <- smpls[1]
  tmp <- tmp4[[smpl]]
  pep[[paste0(pep.ref["Original"], smpl)]] <- 0
  w3 <- which(pep$"Modified sequence" %in% tmp$mod)
  pep[w3, paste0(pep.ref["Original"], smpl)] <- tmp$Intensity[match(pep$"Modified sequence"[w3], tmp$mod)]
}
kol <- paste0(pep.ref["Original"], RSA$values)
kol <- kol[which(kol %in% colnames(pep))]
data <- pep[, c("Modified sequence", kol)]
w <- which(rowSums(data[, kol], na.rm = TRUE) > 0)
data <- data[w,]
pc1 <- prcomp(t(data[, kol]), scale. = TRUE)
dir <- paste0(wd, "/Workflow control/Peptides/PCA plot")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
if (length(pc1$rotation)) {
  scores1 <- as.data.frame(pc1$x)
  if ("PC2" %in% colnames(scores1)) {
    rownames(scores1) <- gsub(topattern(pep.ref["Original"]), "", rownames(scores1))
    scores1[, RSA$names] <- Isapply(strsplit(rownames(scores1), "___"), unlist)
    scores1$Use <- Exp.map$Use[match(rownames(scores1), Exp.map$Ref.Sample.Aggregate)]
    rownames(scores1) <- NULL
    pv1 <- round(100*(pc1$sdev)^2 / sum(pc1$sdev^2), 0)
    pv1 <- pv1[which(pv1 > 0)]
    pv1_ <- paste0("Original: ", paste(vapply(seq_along(pv1), function(x) {
      paste0("PC", x, ": ", pv1[x], "%")
    }, ""), collapse = ", "))
    w <- which(vapply(VPAL$names, function(x) { length(unique(scores1[[x]])) }, 1) > 1)
    w <- w[which(tolower(substr(names(w), 1, 3)) != "rep")]
    scores1$Samples_group <- do.call(paste, c(scores1[, VPAL$names[w], drop = FALSE], sep = " "))
    scores1$Label <- do.call(paste, c(scores1[, RSA$names, drop = FALSE], sep = "-"))
    outlierAnnot_shape %<o% "Replicate"
    outlierAnnot_color %<o% "Samples_group"
    ttl <- "PCA plot - Preliminary - peptide level"
    xLab <- paste0("PC1 = ", pv1[1], "%")
    yLab <- paste0("PC2 = ", pv1[2], "%")
    plot <- ggplot(scores1, aes(x = PC1, y = PC2, colour = .data[[outlierAnnot_color]])) +
      geom_point(aes(shape = .data[[outlierAnnot_shape]])) +
      ggpubr::stat_conf_ellipse(aes(fill = .data[[outlierAnnot_color]]),
                                alpha = 0.2, geom = "polygon", show.legend = FALSE) +
      scale_color_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() +
      xlab(xLab) + ylab(yLab) +
      geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
      ggtitle(ttl#, subtitle = pv1_
              ) +
      geom_text_repel(aes(label = Label), size = 2.5, show.legend = FALSE)
    #poplot(plot)
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    })
    ReportCalls <- AddPlot2Report()
    nReps <- max(as.numeric(Rep))
    Symb <- rep(c("circle", "diamond", "square", "cross", "x"), nReps)[seq_len(nReps)]             
    Symb <- Symb[as.numeric(scores1[[outlierAnnot_shape]])]
    # Custom color scale
    scores1$"Samples group" <- factor(scores1$Samples_group)
    if ("PC3" %in% colnames(scores1)) {
      plot_lyPCA <- plot_ly(scores1, x = ~PC1, y = ~PC2, z = ~PC3,
                            text = ~Label, type = "scatter3d", mode = "markers",
                            color = ~get(outlierAnnot_color), colors = "viridis",
                            symbol = I(Symb))
    } else {
      plot_lyPCA <- plot_ly(scores1, x = ~PC1, y = ~PC2,
                            text = ~Label, type = "scatter", mode = "markers",
                            color = ~`Samples group`, colors = "viridis",
                            symbol = I(Symb))
    }
    plot_lyPCA %<o% layout(plot_lyPCA, title = ttl)
    renderPlotly({ plot_lyPCA <- plot_lyPCA })
    saveWidget(plot_lyPCA, paste0(wd, "/Workflow control/Peptides/PCA plot/", ttl, ".html"),
               selfcontained = TRUE)
    #system(paste0("open \"", wd, "/Workflow control/Peptides/PCA plot/", ttl, ".html"))
  } else {
    stop("There was only one component to the PCA, something must've gone wrong when generating the peptides table!") #(I think this will never happen, the previous check should be identical...?)
  }
} else { stop("There was only one component to the PCA, something must've gone wrong when generating the peptides table!") }

# Peptides heatmap function - useful later and for checking normalisations
pepHtmp %<o% function(intProt = prot.list_pep,
                      Pep = pep,
                      ref = pep.ref[length(pep.ref)],
                      dstDir,
                      ttlRoot = "Peptides log2 heatmap",
                      DB = db,
                      Experiment = Exp,
                      RSA = Ref.Sample.Aggregate,
                      is.log = FALSE,
                      cl,
                      N.clust,
                      N.reserved = 1) {
  TESTING <- FALSE
  #TESTING <- TRUE
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  stopifnot(!is.null(intProt), length(intProt) > 0,
            !is.null(Pep), length(Pep) > 0)
  if ((is.null(ttlRoot))||(!"character" %in% class(ttlRoot))||(!nchar(ttlRoot))) {
    warning("Argument \"ttlRoot\" cannot be NULL!")
    ttlRoot <- "Peptides log2 heatmap"
  }
  if (!dir.exists(dstDir)) { dir.create(dstDir, recursive = TRUE) }
  if (exists("dirlist")) { assign("dirlist", unique(c(dirlist, dstDir)), envir = .GlobalEnv) }
  #
  g <- paste0(ref, RSA$values)
  g <- g[which(g %in% colnames(Pep))]
  g1 <- as.data.frame(t(as.data.frame(strsplit(gsub(topattern(ref), "", g), "___"))))
  colnames(g1) <- RSA$names
  g1$Col <- g
  test <- RSA$names[which(RSA$names != "Replicate")]
  for (i in rev(test)) { g1 <- g1[order(g1[[i]]),] }
  g <- g1$Col
  g1 <- gsub(topattern(ref), "", g)
  g1 <- cleanNms(g1, Experiment = Experiment)
  DB <- DB[match(intProt, DB$"Protein ID"), c("Common Name", "Protein ID", "Sequence")]
  tmpPep <- Pep[grsep2(intProt, Pep$Proteins), c("Proteins", "Sequence", "Modified sequence", g)]
  cleanUp <- FALSE
  if (nrow(tmpPep)) {
    # Create cluster (some steps are slow otherwise)
    if (misFun(cl)) {
      require(parallel)
      dc <- detectCores()
      if (misFun(N.reserved)) { N.reserved <- 1 }
      if (misFun(N.clust)) {
        N.clust <- max(c(dc-N.reserved, 1))
      } else {
        if (N.clust > max(c(dc-N.reserved, 1))) {
          warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
          N.clust <- max(c(dc-N.reserved, 1))
        }
      }
      cl <- makeCluster(N.clust, type = "SOCK")
      cleanUp <- TRUE
    }
    clusterExport(cl, list("DB", "tmpPep", "dir", "g", "g1", "ttlRoot", "is.log", "dstDir"), envir = environment())
    l <- length(intProt)
    tempDat <- data.frame(Protein = rep(intProt, 2),
                          Method = c(rep("Mean", l), rep("ZSc", l)))
    f0 <- function(x) { #x <- tempDat[1,]
      plp <- unlist(x[[1]])
      meth <- x[[2]]
      Plp <- paste(DB[which(DB$"Protein ID" == plp),
                      c("Common Name", "Protein ID")], collapse = " - ")
      grs <- proteoCraft::grsep2(plp, tmpPep$Proteins)
      if (length(grs)) {
        Seq <- DB$Sequence[match(plp, DB$"Protein ID")]
        temp <- tmpPep[grs, c("Sequence", "Modified sequence", g)]
        temp$tst1 <- vapply(temp$Sequence, function(x) { nchar(unlist(strsplit(Seq, x))[1]) }, 1)
        temp$tst2 <- vapply(temp$Sequence, nchar, 1)
        temp$tst3 <- vapply(temp$"Modified sequence", nchar, 1)
        temp <- temp[order(temp$tst1, temp$tst2, temp$tst3),]
        w <- which(rowSums(temp[, g], na.rm = TRUE) == 0)
        if (length(w)) {
          warning("Peptide(s) found with all invalid intensity values!")
          w <- which(rowSums(temp[, g], na.rm = TRUE) > 0)
          temp <- temp[w,]
        }
        if (nrow(temp)) {
          colnames(temp)[match(g, colnames(temp))] <- g1
          # Create heatmap
          temp <- temp[, c("Modified sequence", g1)]
          if (!is.log) {
            warning(" Converting data to log2...")
            temp[, g1] <- suppressWarnings(log2(temp[, g1]))  
          }
          w <- which(!is.finite(as.matrix(temp[, g1])), arr.ind = TRUE)
          temp[, g1][w] <- NA
          M <- rowMeans(temp[, g1], na.rm = TRUE)
          if (meth == "ZSc") { SD <- apply(temp[, g1], 1, sd, na.rm = TRUE) }
          temp[, g1] <- sweep(temp[, g1], 1, M, "-")
          if (meth == "ZSc") { temp[, g1] <- sweep(temp[, g1], 1, SD, "/") }
          temp2 <- proteoCraft::dfMelt(temp[, g1], c("Sample", "value"))
          temp2$"Modified sequence" <- temp$"Modified sequence"
          temp2$Sample <- as.character(temp2$Sample)
          temp2$Xmin <- match(temp2$Sample, colnames(temp))-1
          temp2$Xmax <- temp2$Xmin+1
          temp2$Ymax <- nrow(temp):1
          temp2$Ymin <- temp2$Ymax-1
          hlab <- aggregate(temp2[, c("Xmin", "Xmax")], list(temp2$Sample), unique)
          colnames(hlab)[1] <- "Sample"
          hlab$X <- rowMeans(hlab[, c("Xmin", "Xmax")])
          vlab <- aggregate(temp2[, c("Ymin", "Ymax")], list(temp2$"Modified sequence"), unique)
          colnames(vlab)[1] <- "Modified sequence"
          vlab$Y <- rowMeans(vlab[, c("Ymin", "Ymax")])
          Xscale <- max(temp2$Xmax)
          Yscale <- max(temp2$Ymax)
          # Create heatmap plot
          Xlim <- c(-1, Xscale+15)
          Ylim <- c(-6, Yscale+20)
          temp2a <- temp2[, c("Xmin", "Ymin", "value")]
          Splits <- 20
          XScale2 <- Xscale*0.1/Splits
          temp2b <- data.frame(Xmin = Xscale/2 - XScale2*((-Splits/2):(Splits/2)),
                               Ymin = -5)
          temp2b$Xmax <- temp2b$Xmin + XScale2
          temp2b$value <- min(temp2a$value, na.rm = TRUE) + (Splits:0)*(max(temp2a$value, na.rm = TRUE)-min(temp2a$value, na.rm = TRUE))/Splits
          temp2b$Label <- ""
          temp2b$Label[c(1, Splits/2+1, Splits+1)] <- round(temp2b$value[c(1, Splits/2+1, Splits+1)], 1)
          # Create graph
          ttl <- paste0(ttlRoot, " - ", gsub("/", "-", Plp), c("", " (Z-scored)")[(meth == "ZSc")+1])
          heatmap.plot <- ggplot2::ggplot(temp2a) +
            ggplot2::geom_rect(ggplot2::aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
            ggplot2::geom_rect(data = temp2b, ggplot2::aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymin+1, fill = value)) +
            ggplot2::geom_text(data = hlab, ggplot2::aes(x = X, label = Sample), y = Yscale+1, angle = 60, hjust = 0, size = 2.5) +
            ggplot2::geom_text(data = vlab, ggplot2::aes(y = Y, label = `Modified sequence`), x = Xscale+0.5, hjust = 0, size = 1.7) +
            ggplot2::geom_text(data = temp2b, ggplot2::aes(x = Xmin, label = Label), y = -3.5, hjust = 0.5, size = 2) +
            ggplot2::ggtitle(paste0(ttlRoot, ", ", c("mean-normalized", "Z-scored")[(meth == "ZSc")+1]),
                             subtitle = Plp) +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                           axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), 
                           panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
                           plot.margin = ggplot2::margin(0, 0, 0, 0, "cm")) +
            ggplot2::coord_fixed(0.5) +
            ggplot2::scale_fill_gradient2(low = "red", mid = "black", high = "green", na.value = "lightblue") +
            ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + ggplot2::theme(legend.position = "none") +
            ggplot2::xlim(Xlim[1], Xlim[2]) + ggplot2::ylim(Ylim[1], Ylim[2])
          #proteoCraft::poplot(heatmap.plot)
          suppressMessages({
            ggplot2::ggsave(paste0(dstDir, "/", ttl, ".jpeg"), heatmap.plot,
                            dpi = 600, width = 20, height = 12, units = "in")
            ggplot2::ggsave(paste0(dstDir, "/", ttl, ".pdf"), heatmap.plot,
                            dpi = 600, width = 20, height = 12, units = "in")
          })
          #system(paste0("open \"", dstDir, "/", ttl, ".jpeg", "\""))
          #system(paste0("open \"", dstDir, "/", ttl, ".pdf", "\""))
        }
      }
    }
    environment(f0) <- .GlobalEnv
    parApply(cl, tempDat, 1, f0)
  }
  if (cleanUp) { stopCluster(cl) }
}
# Another useful function for checking peptide normalisations
pepPlotFun %<o% function(df1,
                         df2,
                         ttl,
                         dstDir,
                         save = TRUE,
                         xpMap = Exp.map,
                         VPAL = Volcano.plots.Aggregate.Level) {
  tst1 <- df1
  tst2 <- df2
  colnames(tst1) <- gsub(".* - ", "", colnames(tst1))
  colnames(tst2) <- gsub(".* - ", "", colnames(tst2))
  tst1 <- proteoCraft::dfMelt(tst1)
  tst2 <- proteoCraft::dfMelt(tst2)
  tst1$Norm <- "Original"
  tst2$Norm <- "Re-normalized"
  tst1 <- rbind(tst1, tst2)
  rm(tst2)
  tst1 <- tst1[which(is.finite(tst1$value)),]
  g <- grepl("_REF\\.to\\.REF_", tst1$variable)
  tst1$Type <- c("Samples", "References")[g+1]
  tst1$Group <- xpMap[match(tst1$variable, xpMap$Ref.Sample.Aggregate),
                      VPAL$column]
  tst1$variable <- cleanNms(tst1$variable)
  w <- which(g)
  tst1$variable[w] <- gsub_Rep("_REF\\.to\\.REF_.*", "", tst1$variable[w])
  tst1$Group <- cleanNms(tst1$Group)
  tst1$Group[which(is.na(tst1$Group))] <- "References"
  tst1$Norm <- as.factor(tst1$Norm)
  tst1$Type <- as.factor(tst1$Type)
  plot <- ggplot(tst1) +
    geom_density(stat = "density", alpha = 0.1, aes(x = value, colour = variable, fill = variable)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(Norm~Group) + ggtitle(ttl) + theme_bw()
  if (grepl("ratio", ttl, ignore.case = TRUE)) { ntrcpt <- 0 } else {
    ntrcpt <- median(tst1$value[which((tst1$Group != "References")&
                                        (tst1$Norm == "Original"))])
  }
  plot <- plot + geom_vline(xintercept = ntrcpt, linetype = "dashed")
  print(plot) # This type of QC plot does not need to pop up, the side panel is fine
  if (save) {
    suppressMessages({
      ggsave(paste0(dir[1], "/", ttl, ".jpeg"), plot, dpi = 150)
      ggsave(paste0(dir[1], "/", ttl, ".pdf"), plot, dpi = 150)
    })
  }
}

####################################################
# Optional: choose whether to remove any outliers  #
####################################################
Src <- paste0(libPath, "/extdata/R scripts/Sources/remove_Outliers.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

g <- grep(topattern(pep.ref["Original"]), colnames(pep), value = TRUE)
# View(pep[, g])
test <- rowSums(pep[, g])
l <- length(which(test == 0))
if (l) {
  msg <- paste0("Removing ", l, " peptide", c("", "s")[(l > 1)+1], " with invalid expression values - this is unexpected, investigate!")
  ReportCalls <- AddMsg2Report(Space = FALSE, Warning = TRUE)
  pep <- pep[which(test > 0),]
  w <- which(ev$id %in% unique(as.integer(unlist(strsplit(pep$"Evidence IDs", ";")))))
  ev <- ev[w,]
}
DatAnalysisTxt <- paste0(DatAnalysisTxt, " The long format ", c("evidence.txt", "main report", "psm.tsv")[match(SearchSoft, SearchSoftware)],
                         " table was consolidated into a wide format peptidoforms table, summing up quantitative values where necessary.")

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Impute missing peptide intensities
Src <- paste0(libPath, "/extdata/R scripts/Sources/pep_Impute.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Re-normalize peptide intensities
rfnm <- c("Original", "Imputation")[Impute+1]
Src <- paste0(libPath, "/extdata/R scripts/Sources/pepNorm_VarPlot.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# - Run normalisations
if (Param$Norma.Pep.Intens) {
  Src <- paste0(libPath, "/extdata/R scripts/Sources/pepNorm_Main.R")
  #rstudioapi::documentOpen(Src)
  #rstudioapi::documentOpen(nrmSrc)
  source(Src, local = FALSE)
}
#
rfnm <- names(pep.ref)[length(pep.ref)]
Src <- paste0(libPath, "/extdata/R scripts/Sources/pepNorm_VarPlot.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#View(pep[, grep(topattern(pep.ref[length(pep.ref)]), colnames(pep), value = TRUE)]) # Check final data visually

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)

# Calculate peptide ratios
# Step 1: calculate average references
pep.ratios.root %<o% "Ratio"
pep.ratios.ref %<o% c(Ratios = paste0("log2(", pep.ratios.root, ") - "))
# !!!!!
# !!!!! This code is mirrored in Prot.Quant!!!
# !!!!! When modifying it, either update that too, or create a function for both to ensure consistency!!!
# !!!!! If choosing the second option, careful to remember what is/isn't log scale!!!
# !!!!!
# Please rewrite me! I was written by a beginner (you) and I look very boorish!
kount <- 0
for (i in RRG$values) { #i <- RRG$values[1]
  j <- setNames(unlist(strsplit(i, "___")), RRG$names)
  temp <- lapply(RRG$names, function(x) { which(Exp.map[[x]] == j[[x]]) })
  temp2 <- sort(unique(unlist(temp)))
  test <- vapply(temp2, function(x) { sum(vapply(temp, function(y) { x %in% unlist(y) }, TRUE)) }, 1)
  temp2 <- temp2[which(test == length(temp))]
  temp3 <- Exp.map[temp2,]
  temp3 <- temp3[which(temp3$Reference),]
  b <- temp3$Ref.Sample.Aggregate
  if (length(b)) {
    kol <- paste0(pep.ref[length(pep.ref)], b)
    w <- which(kol %in% colnames(pep))
    if (length(w)) {
      b1 <- apply(pep[, kol[w], drop = FALSE], 1, function(x) {
        2^log_ratio_av(log2(x))
      })
      b2 <- paste0(pep.ref[length(pep.ref)], i, ".REF")
      #print(b2)
      pep[, b2] <- b1
      kount <- kount + 1
      if (kount == 1) {
        ratios.2.ref %<o% data.frame(Name = b2, Source = paste(b, collapse = ";"))
      } else { ratios.2.ref <- rbind(ratios.2.ref, c(b2, paste(b, collapse = ";"))) }
    } else { warning(paste0("Empty group: ", i, ", skipping!")) }
  } else { warning(paste0("There is no reference for level ", i)) }
}
ratios.2.ref$Source <- strsplit(ratios.2.ref$Source, ";")
ratios.2.ref$Used_by <- list(NA)
# Step 2: calculate individual ratios to relevant reference
# If the reference is an average, we will also calculate individual ref ratios to it;
# this will be useful further down the line.
for (i in RRG$values) { #i <- RRG$values[1]
  j <- set_names(unlist(strsplit(i, "___")), RRG$names)
  temp <- lapply(RRG$names, function(x) { which(Exp.map[[x]] == j[[x]]) })
  temp2 <- sort(unique(unlist(temp)))
  test <- vapply(temp2, function(x) { sum(vapply(temp, function(y) { x %in% unlist(y) }, TRUE)) }, 1)
  temp2 <- temp2[which(test == length(temp))]
  temp3 <- Exp.map[temp2,]
  # Get reference
  k <- j[which(names(j) %in% RRG$names)]
  b <- paste0(pep.ref[length(pep.ref)], paste(k, collapse = "___"), ".REF")
  b1 <- pep[[b]]
  a <- temp3$Ref.Sample.Aggregate
  if (length(which(temp3$Reference)) == 1) { # If there is only one ref in the group, remove it as there is no point calculating a ratio to itself
    a <- temp3$Ref.Sample.Aggregate[which(!temp3$Reference)]
  }
  if ((length(a) == 0)||(!sum(paste0(pep.ref[length(pep.ref)], a) %in% colnames(pep)))) {
    warning(paste0("There are no ratios to calculate for level ", i))
  } else {
    pep[, paste0(pep.ratios.root, " - ", a)] <- sweep(pep[, paste0(pep.ref[length(pep.ref)], a), drop = FALSE], 1, b1, "/")
    pep[, paste0(pep.ratios.ref, a)] <- log2(pep[, paste0(pep.ratios.root, " - ", a)])
    a2 <- cleanNms(a)
    cat(paste0("Median log2 ratio: ", paste(paste0("\n - ", a2, " -> ", apply(pep[,paste0(pep.ratios.ref, a), drop = FALSE], 2, function(x) { signif(median(is.all.good(x)), 3) })), collapse = ""), "\n"))
  }
  ratios.2.ref[["Used_by"]][which(ratios.2.ref$Name == b)] <- list(a)
}
ratios.2.ref <- listMelt(ratios.2.ref$Used_by, ratios.2.ref$Name, c("Ratio", "Reference"))
#View(ratios.2.ref)
# View(pep[,grep(topattern("Ratio "), colnames(pep))])

# Visualize:
kol <- paste0(pep.ratios.ref[1], RSA$values)
kol <- kol[which(kol %in% colnames(pep))]
temp <- pep[, kol]
colnames(temp) <- gsub_Rep(topattern(pep.ratios.ref[1]), "", colnames(temp))
temp <- reshape2::melt(temp, measure.vars = colnames(temp))
temp <- temp[which(is.all.good(temp$value, 2)),]
temp[, RSA$names] <- Isapply(strsplit(as.character(temp$variable), "___"), unlist)
temp$variable <- as.factor(cleanNms(as.character(temp$variable)))
temp[, c("Ratios group", "#")] <- Isapply(strsplit(as.character(temp$variable), "_REF\\.to\\.REF_"), unlist)
aggr <- RSA$names
tst <- vapply(aggr, function(a) { length(get(substr(a, 1, 3))) }, 1)
aggr <- aggr[which(tst > 1)]
facets <- (length(aggr) > 0)
if (facets) {
  if (length(aggr) > 1) { form <- paste0(aggr[1], "~", paste(aggr[2:length(aggr)], collapse = "+")) } else {
    form <- paste0(".~", aggr[1])
  }
  form <- as.formula(form)
}
dir <- paste0(wd, "/Workflow control/Peptides/Ratios")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
ttl <- "Peptide ratios distribution"
plot <- ggplot(temp) +
  geom_histogram(aes(x = value, fill = variable), bins = ceiling((max(temp$value)-min(temp$value))/0.05)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d(begin = 0.25) +
  ggtitle(ttl) + theme_bw() + theme(legend.position = "none")
if (facets) { plot <- plot + facet_grid(form) }
#poplot(plot, 12, 20)
print(plot)
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
})
ReportCalls <- AddPlot2Report()

# Create peptide-level Ref-to-Ref ratios (useful for PTMs analysis):
RatConGrps %<o% Param$Ratios.Contaminant.Groups
if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
  pep.Ref.Ratios %<o% NULL
}
if (Param$Ratios.Thresholds == threshMsg) {
  warning("I am duplicate code! Make me a function!!!")
  if (RefRat_Mode == "1") {
    temp <- lapply(RG$values, function(x) { #x <- RG$values[1]
      if (RatConGrps == "Ratio groups") {
        x1 <- unique(Exp.map[which(Exp.map[[RG$column]] == x), RRG$column])
      }
      if (RatConGrps == "Experiments") {
        x1 <- unique(Exp.map$Experiment[which(Exp.map[[RG$column]] == x)])
        x1 <- unique(Exp.map[which(Exp.map$Experiment == x1), RRG$column])
      }
      if (RatConGrps == "Whole dataset") {
        x1 <- unique(Exp.map[[RRG$column]])
      }
      xr1 <- paste0(pep.ref[length(pep.ref)], x1, ".REF")
      xr1 <- xr1[which(xr1 %in% colnames(pep))]
      if (length(xr1) > 1) {
        perm <- gtools::permutations(length(xr1), 2, xr1)
        tmp <- apply(perm, 1, function(y) { log2(pep[[y[[1]]]]/pep[[y[[2]]]]) })
      } else {
        x2 <- unique(Exp.map$Ref.Sample.Aggregate[which((Exp.map[[RG$column]] == x)&(Exp.map$Reference))])
        xr2 <- paste0(pep.ref[length(pep.ref)], x2)
        tmp <- vapply(xr2, function(y) { log2(pep[[y]]/pep[[xr1]]) }, 1)
      }
      tmp <- as.data.frame(tmp)
      colnames(tmp) <- paste0(pep.ratios.ref[1], x, "_REF.to.REF_", 1:ncol(tmp))
      return(tmp)
    })
  }
  if (RefRat_Mode == "2") {
    temp <- lapply(RG$values, function(rtGrp) { #rtGrp <- RG$values[1]
      if (RatConGrps == "Ratio groups") {
        em <- Exp.map[which(Exp.map[[RG$column]] == rtGrp),]
      }
      if (RatConGrps == "Experiments") {
        xp <- unique(Exp.map$Experiment[which(Exp.map[[RG$column]] == rtGrp)])
        em <- Exp.map[which(Exp.map$Experiment == xp),]
      }
      if (RatConGrps == "Whole dataset") {
        em <- Exp.map
      }
      grps <- unique(em[[VPAL$column]])
      x <- lapply(grps, function(grp) {  #grp <- grps[1]
        y <- unique(em$Ref.Sample.Aggregate[which(em[[VPAL$column]] == grp)])
        y <- paste0(pep.ref[length(pep.ref)], y)
        y <- y[which(y %in% colnames(pep))]
        if (length(y)) {
          y <- gtools::permutations(length(y), 2, y) # Important to use permutations and not combinations since we want symmetry!
          y <- apply(y, 1, function(z) {
            suppressWarnings(log2(pep[[z[[1]]]]/pep[[z[[2]]]]))
          })
        } else { y <- NULL }
        return(y)
      })
      x <- as.data.frame(x)
      colnames(x) <- paste0(pep.ratios.ref[1], rtGrp, "_REF.to.REF_", 1:ncol(x))
      return(x)
    })
  }
  L <- length(temp)
  if (L) {
    names(temp) <- NULL
    temp <- data.frame(temp, check.names = FALSE)
    pep.Ref.Ratios %<o% temp
  } else { stop("No Ref to Ref columns were generated, investigate!!!") }
  pep[, colnames(pep.Ref.Ratios)] <- pep.Ref.Ratios
  # Visualize:
  temp <- pep[, colnames(pep.Ref.Ratios)]
  colnames(temp) <- gsub(topattern(pep.ratios.ref[length(pep.ratios.ref)]), "", colnames(temp))
  temp <- suppressMessages(reshape2::melt(temp))
  temp <- temp[which(is.all.good(temp$value, 2)),]
  tmp <- data.frame(Var = as.character(unique(temp$variable)))
  tmp$Var2 <- cleanNms(tmp$Var)
  tmp$Var2 <- gsub("_REF\\.to\\.REF_", " - ", tmp$Var2)
  tmp2 <- Isapply(strsplit(as.character(tmp$Var), "_REF\\.to\\.REF_"), unlist)
  colnames(tmp2) <- c("Ratios group", "#")
  tmp2$"Ratios group" <- cleanNms(tmp2$"Ratios group")
  temp[, colnames(tmp2)] <- tmp2[match(temp$variable, tmp$Var),]
  temp$variable <- tmp$Var2[match(temp$variable, tmp$Var)]
  temp$`#` <- as.integer(temp$`#`)
  dir <- paste0(wd, "/Workflow control/Peptides/Ratios")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  ttl <- "Peptide ratios distribution (reference to reference)"
  plot <- ggplot(temp) +
    geom_histogram(aes(x = value, fill = variable), bins = ceiling((max(temp$value)-min(temp$value))/0.05)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_d(begin = 0.25) +
    ggtitle(ttl) + theme_bw() + theme(legend.position = "none") +
    #facet_grid(`Ratios group` ~ `#`)
    facet_grid(.~`Ratios group`)
  print(plot) # This type of QC plot does not need to pop up, the side panel is fine
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 4*(length(unique(temp$"Ratios group"))+1),
           height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 4*(length(unique(temp$"Ratios group"))+1),
           height = 10, units = "in")
  })
  ReportCalls <- AddPlot2Report()
}

#### Code chunk - Normalize peptide ratios
# Sort of legacy, I really do not think that this is a good idea at this stage if intensities have been well normalized.
# Even if we are dealing with SILAC (not supported yet) where ratios are better measured and should be more precise,
# I would simply re-calculate updated intensities to reflect the ratios at an early stage, then process intensities from there.
if (Param$Norma.Pep.Ratio) {
  # Step 1:
  # Classic normalisation by the median for each column;
  # NB: This does not work the same way as intensities, so here the global scale of ratios need not be conserved:
  # We really want them centred on 0 in log scale.
  Norm.log2.Pep.Ratios <- c()
  NormGrps <- unique(pep$`Normalisation group`)
  g <- grep(topattern(pep.ratios.ref[1]), colnames(pep), value = TRUE)
  pep[, paste0("norm. ", g)] <- NA
  for (nrmgrp in NormGrps) {
    w <- which(pep$`Normalisation group` == nrmgrp)
    for (k in g) {
      m <- median(is.all.good(pep[w, k]))
      #m <- mlv(is.all.good(pep[[k]]), method = "Parzen")[1]
      pep[w, paste0("norm. ", k)] <- pep[w, k]-m
      Norm.log2.Pep.Ratios[paste0(nrmgrp, " - ", gsub(topattern(pep.ratios.ref[1]), "", k))] <- m
    }
  }
  pep.ratios.ref <- unique(c(pep.ratios.ref, paste0("norm. ", pep.ratios.ref[1])))
  # Step 2
  # Optional advanced normalisation:
  if (("Adv.Norma.Pep.Ratio" %in% colnames(Param))&&(Param$Adv.Norma.Pep.Ratio != FALSE)) {
    if (Param$Adv.Norma.Pep.Ratio.Type == "C") {
      k <- Adv.Norma.Pep.Ratio.Type.Group$column
      test <- vapply(Adv.Norma.Pep.Ratio.Type.Group$values, function(i) { #i <- Adv.Norma.Pep.Ratio.Type.Group$values[1]
        i <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[k]] == i)]
        return(length(which(paste0(pep.ratios.ref[length(pep.ratios.ref)], i) %in% colnames(pep))))
      }, 1)
      agg <- Adv.Norma.Pep.Ratio.Type.Group$values[which(test > 1)]
      exports <- list("agg", "Adv.Norma.Pep.Ratio.Type.Group", "Exp.map", "pep.ratios.ref", "pep", "Param")
      clusterExport(parClust, exports, envir = environment())
      invisible(clusterCall(parClust, function() {
        library(proteoCraft)
        return()
      }))
      norm_temp <- parSapply(parClust, 0:length(agg), function(i) { #i <- 1
        if (i == 0) {
          kol <- grep(paste0(topattern(pep.ratios.ref[1]), ".+_REF\\.to\\.REF_"), colnames(pep), value = TRUE)
        } else {
          x <- agg[i]
          j <- unlist(strsplit(x, "___"))
          names(j) <- Adv.Norma.Pep.Ratio.Type.Group$names
          temp <- lapply(Adv.Norma.Pep.Ratio.Type.Group$names, function(x) {
            which(Exp.map[[x]] == j[[x]])
          })
          temp2 <- sort(unique(unlist(temp)))
          test <- vapply(temp2, function(x) { sum(vapply(temp, function(y) { x %in% unlist(y) }, TRUE)) }, 1)
          temp2 <- temp2[which(test == length(temp))]
          temp3 <- Exp.map[temp2,]
          kol <- paste0(pep.ratios.ref[length(pep.ratios.ref)], temp3$Ref.Sample.Aggregate)
          w <- which(kol %in% colnames(pep))
          temp3 <- temp3[w,]
          kol <- kol[w]
          if (length(w) <= 1) {
            warning(paste0("There are no columns to normalize for samples group ", x))
          }
        }
        if (length(kol) > 1) {
          temp2 <- temp <- pep[, c("Modified sequence", kol)]
          colnames(temp2) <- gsub("^AdvNorm\\.norm\\.", "AdvNorm. ", colnames(temp2))
          for (nrmgrp in unique(pep$`Normalisation group`)) {
            w <- which(pep$`Normalisation group` == nrmgrp)
            tmp <- AdvNorm.IL(temp[w,], "Modified sequence", kol, TRUE, 5)
            colnames(tmp) <- gsub("^AdvNorm\\.norm\\.", "AdvNorm. ", colnames(tmp))
            temp2[w, colnames(tmp)] <- tmp[, colnames(tmp)]
          }
          #cat("Advanced peptides ratio normalisation done for samples group", x, "\n")
          temp2$"Modified sequence" <- NULL
          return(list(temp2))
        }
      })
    } else {
      stop("Not implemented yet, I need to update the current AdvNorm function as it takes too long or crashes.")
    }
    for (i in seq_along(norm_temp)) {
      tmp <- norm_temp[[i]]
      pep[, colnames(tmp)] <- tmp
    }
    pep.ratios.ref <- unique(c(pep.ratios.ref, paste0("AdvNorm. ", pep.ratios.ref[1])))
  }
  if (Param$Norma.Pep.Ratio.show) {
    for (i in Ratios.Plot.split$values) { #i <- Ratios.Plot.split$values[1]
      j <- unlist(strsplit(i, "___"))
      names(j) <- unlist(Aggregate.map$Characteristics[which(Aggregate.map$Aggregate.Name == Ratios.Plot.split$aggregate)])
      #k <- lapply(c(seq_along(j)), function(x) { which((Exp.map[[names(j)[x]]] == j[x])&(!Exp.map$Reference)) })
      k <- lapply(c(seq_along(j)), function(x) { which(Exp.map[[names(j)[x]]] == j[x]) })
      l <- sort(unique(unlist(k)))
      test <- vapply(l, function(x) { sum(vapply(k, function(y) { x %in% y }, TRUE)) == length(k) }, 1)
      temp <- Exp.map$Ref.Sample.Aggregate[l[which(test)]]
      a1 <- paste0(pep.ratios.ref[1], temp)
      a2 <- paste0(pep.ratios.ref[length(pep.ratios.ref)], temp)
      a1 <- a1[which(a1 %in% colnames(pep))]
      a2 <- a2[which(a2 %in% colnames(pep))]
      if (length(a1)) {
        temp <- pep[, c("Modified sequence", a1, a2)]
        temp <- reshape2::melt(temp)
        temp$Norm <- grepl(topattern(pep.ratios.ref[length(pep.ratios.ref)]), temp$variable)
        temp$Norm <- c("Original", "Normalised")[temp$Norm + 1]
        temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
        temp$variable <- as.character(temp$variable)
        temp$variable[which(temp$Norm == "Original")] <- gsub_Rep(topattern(pep.ratios.ref[1]), "", temp$variable[which(temp$Norm == "Original")])
        temp$variable[which(temp$Norm == "Normalised")] <- gsub_Rep(topattern(pep.ratios.ref[length(pep.ratios.ref)]), "", temp$variable[which(temp$Norm == "Normalised")])
        temp2 <- Isapply(strsplit(temp$variable, "___"), unlist)
        colnames(temp2) <- unlist(Aggregate.map$Characteristics[which(Aggregate.map$Aggregate.Name == RSA$aggregate)])
        temp[, colnames(temp2)] <- temp2
        if (length(Ratios.Plot.wrap$names) > 1) {
          temp$Wrap <- do.call(paste, c(temp2[, Ratios.Plot.wrap$names], sep = "_"))
        } else { temp$Wrap <- temp2[[Ratios.Plot.wrap$names]] }
        if (length(Ratios.Plot.colour$names) > 1) {
          temp$Colour <- do.call(paste, c(temp2[, Ratios.Plot.colour$names], sep = "_"))
        } else { temp$Colour <- temp2[[Ratios.Plot.colour$names]] }
        temp$X <- do.call(paste, c(temp[, c("Norm", "Colour")], sep = "_"))
        temp$X <- factor(temp$X, levels = unique(unlist(vapply(c("Original", "Normalised"), function(x) {
          paste(x, unique(temp2[[Ratios.Plot.colour$names]]), sep = "_")
        }, ""))))
        temp <- temp[which(is.all.good(temp$value, 2)),]
        dir <- paste0(wd, "/Workflow control/Peptides/Ratios")
        if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
        dirlist <- unique(c(dirlist, dir))
        ttl <- paste0("Peptide ratios distribution_sample group: ", i)
        plot <- ggplot(temp) +
          geom_violin(aes(x = X, y = value, color = Colour, fill = Colour), alpha = 0.25) +
          geom_boxplot(aes(x = X, y = value, color = Colour, fill = Colour), alpha = 0.5) +
          theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          scale_color_viridis_d(begin = 0.25) +
          scale_fill_viridis_d(begin = 0.25) +
          facet_grid(. ~ Norm, scales = "free", space = "free") +
          ggtitle(ttl)
        if (length(unique(temp$Wrap)) > 1) { plot <- plot + facet_wrap(~Wrap) }
        print(plot) # This type of QC plot does not need to pop up, the side panel is fine
        suppressMessages({
          ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".jpeg"), plot,
                 dpi = 300, width = 10, height = 10, units = "in")
          ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".pdf"), plot,
                 dpi = 300, width = 10, height = 10, units = "in")
        })
        ReportCalls <- AddPlot2Report(Title = ttl)
      } else { warning(paste0("Nothing to plot for level ", i)) }
    }
    # Also look at Ref-to-Ref ratios:
    a1 <- grep(paste0(topattern(pep.ratios.ref[1]), ".+_REF.to.REF_[0-9]+$"), colnames(pep), value = TRUE)
    a2 <- grep(paste0(topattern(pep.ratios.ref[length(pep.ratios.ref)]), ".+_REF.to.REF_[0-9]+$"), colnames(pep), value = TRUE)
    if (length(a1)) {
      temp <- pep[, c("Modified sequence", a1, a2)]
      temp <- reshape2::melt(temp)
      temp$Norm <- grepl(topattern(pep.ratios.ref[length(pep.ratios.ref)]), temp$variable)
      temp$Norm <- c("Original", "Normalised")[temp$Norm + 1]
      temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
      temp$variable <- as.character(temp$variable)
      temp$variable[which(temp$Norm == "Original")] <- gsub_Rep(topattern(pep.ratios.ref[1]), "", temp$variable[which(temp$Norm == "Original")])
      temp$variable[which(temp$Norm == "Normalised")] <- gsub_Rep(topattern(pep.ratios.ref[length(pep.ratios.ref)]), "", temp$variable[which(temp$Norm == "Normalised")])
      temp$Ratios.Group <- gsub_Rep("_REF.to.REF_[0-9]+$", "", temp$variable)
      temp <- temp[which(is.all.good(temp$value, 2)),]
      dir <- paste0(wd, "/Workflow control/Peptides/Ratios")
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      ttl <- "Peptide ratios distribution_References-to-References"
      plot <- ggplot(temp) +
        geom_violin(aes(x = variable, y = value, color = variable, fill = variable), alpha = 0.25) +
        geom_boxplot(aes(x = variable, y = value, color = variable, fill = variable), alpha = 0.5) +
        scale_color_viridis_d(begin = 0.25) +
        scale_fill_viridis_d(begin = 0.25) +
        theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        facet_wrap(~Ratios.Group, scales = "free") +
        ggtitle(ttl)
      print(plot) # This type of QC plot does not need to pop up, the side panel is fine
      suppressMessages({
        ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".jpeg"), plot,
               dpi = 300, width = 10, height = 10, units = "in")
        ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".pdf"), plot,
               dpi = 300, width = 10, height = 10, units = "in")
      })
      ReportCalls <- AddPlot2Report(Title = ttl)
      DatAnalysisTxt <- paste0(DatAnalysisTxt, " Peptide ratios were then re-normalized.")
    } else { warning("Nothing to plot for Reference-to-Reference ratios!") }
  }
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
# It makes sense to close/re-create parallel clusters regularly to reduce memory usage + avoid corruption
stopCluster(parClust)
source(parSrc, local = FALSE)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

useSAM %<o% ((names(pvalue.col)[which(pvalue.use)] == "Student")&&(useSAM_thresh)&&)

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
            "character" %in% class(tst3$x),
            length(which(temp$value != temp$Accession)) == 0)
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
    temp[[i]] <- parSapply(parClust, temp[[i]], function(x) { paste(unique(unlist(x)), collapse = ";") }) # Do not use sort here or it will break the correspondance between "GO" and "GO-ID"
  }
  tst1 <- unlist(strsplit(temp$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(temp$GO, ";"))
  tst3 <- data.table(A1 = tst1, A2 = tst2)
  tst3 <- tst3[, list(x = unique(A2)), by = list(Group.1 = A1)]
  tst3 <- as.data.frame(tst3)
  stopifnot(length(tst3$x) == length(unique(tst3$x)), "character" %in% class(tst3$x))
  #
  PG[, kol2] <- temp[match(PG$id, temp$Group.1), kol2]
  #
  #View(tst3[which(vapply(tst3$x, length, 1) > 1),])
  #View(tst3[which(vapply(tst3$x, length, 1) == 0),])
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
  Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_prepare.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

# Define design matrix and contrasts (limma) 
#  From Exp.map to design matrix
Coefficients %<o% Factors[which(!Factors %in% c("Experiment", "Replicate"))]
w <- which(vapply(Coefficients, function(x) { length(unique(FactorsLevels[[x]])) < nrow(Exp.map) }, 1))
Coefficients <- Coefficients[w]
w <- which(vapply(Coefficients, function(x) { length(unique(Exp.map[[x]])) > 1 }, 1))
Coefficients <- Coefficients[w]
EM <- Exp.map
EM <- EM[order(#EM[[RRG$column]], # Do not use RRG here
  EM[[RG$column]], # Cf. below: safe because each RG contains at least one ref samples group
  EM$Reference, # Very important: if any level is dropped from the design matrix, it must be a reference!
  # (Otherwise every gets confusing and  my head starts hurting...)
  EM[[VPAL$column]],
  EM$Replicate),]
Group_ <- do.call(paste, c(EM[, Coefficients, drop = FALSE], sep = "_"))
Group_ <- as.factor(Group_)
EM$Group_ <- Group_
if (Nested) {
  EM$Replicate_ <- Replicate_ <- as.factor(EM$Replicate)
  designMatr %<o% model.matrix(~0 + Replicate_ + Group_)
} else {
  designMatr %<o% model.matrix(~0 + Group_)
}
rownames(designMatr) <- EM$Ref.Sample.Aggregate
#
# Define contrasts
expContrasts %<o% list()
for (ratGrp in RG$values) { #ratGrp <- RG$values[1]
  em <- EM[which(EM[[RG$column]] == ratGrp),]
  grp1 <- unique(em$Group_[which(!em$Reference)])
  vpal1 <- em[match(grp1, em$Group_), VPAL$column]
  grp0 <- unique(em$Group_[which(em$Reference)])
  expContrasts[[ratGrp]] <- plyr::rbind.fill(lapply(grp0, function(g0) {
    data.frame(x1 = paste0("Group_", grp1),
               x0 = paste0("Group_", g0),
               name = vpal1)
  }))
}
expContrasts <- plyr::rbind.fill(expContrasts)
expContrasts$All <- lapply(1:nrow(expContrasts), function(x) { gsub("^Group_", "", expContrasts[x, c("x1", "x0")]) })
expContrasts$Contrasts <- apply(expContrasts[, c("x1", "x0")], 1, function(x) {
  x <- x[which(x %in% colnames(designMatr))]
  paste(x, collapse = " - ")
})
contrCall <- paste0("contrMatr %<o% makeContrasts(",
                    paste(expContrasts$Contrasts, collapse = ", "),
                    ", levels = designMatr)")
#cat(contrCall, "\n")
eval(parse(text = contrCall), envir = .GlobalEnv)
# (NB: Contrasts could be renamed to something shorter, e.g. makeContrasts(Comp1 = A - B, Comp2 = A - C)
#
bhFDRs %<o% sort(BH.FDR, decreasing = FALSE)
samRoots %<o% c(samRoot,
                paste0("SAM regulated-FDR=", paste(100*bhFDRs, collapse = "/"), "% FDR - "))
samSubDir %<o% "Reg. analysis/SAM"
ebamSubDir %<o% "Reg. analysis/EBAM"
ebamRoot %<o% paste0("EBAM regulated-FDR=", paste(100*bhFDRs, collapse = "/"), "% FDR - ")
#
if ("Mirror.Ratios" %in% colnames(Param)) { Mirror.Ratios <- Param$Mirror.Ratios <- as.logical(Param$Mirror.Ratios) }
if (!is.logical(Mirror.Ratios)) {
  warning("I could not make sense of the value of parameter \"Mirror.Ratios\", defaulting to FALSE!")
  Mirror.Ratios <- FALSE
}
if (Mirror.Ratios) {
  stop("The code has changed, this is now deprecated: if you want to re-use this then first check all normalisation steps where ratios are calculated after this stage!!!")
}
#
samDir <- paste0(wd, "/", samSubDir)
ebamDir <- paste0(wd, "/", ebamSubDir)
pvalue.use %<o% (names(pvalue.col) == Param$P.values.type)
AltHyp %<o% c(c("greater", "lower")[Mirror.Ratios+1], "two.sided")[TwoSided+1]
Av_SE_fun %<o% function(vect) {
  res <- proteoCraft::is.all.good(as.numeric(vect))
  if (length(res)) { res <- c(mean(res), sd(res)/sqrt(length(res))) } else {
    res <- unique(vect[which((!is.nan(vect))&(!is.na(vect)))])
    if (length(res) == 1) { return(c(res, NA)) } else { return(c(NA, NA)) }
    return(res)
  }
}
limma.one.sided %<o% function(myFit, lower) {
  se.coef <- sqrt(myFit$s2.post) * myFit$stdev.unscaled
  df.total <- myFit$df.prior + myFit$df.residual
  rs <- pt(myFit$t, df = df.total, lower.tail = lower)
  return(rs[, colnames(myFit$p.value)])
}
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
filter_types[grep("^dat.+2$", filter_types, invert = TRUE)] <- substr(filter_types[which(!grepl("^dat.+2$", filter_types))], 1, 3)
filter_types[grep("^dat.+2$", filter_types)] <- "dat2"
filter_types <- unique(c("con", filter_types))
if ("ref" %in% filter_types) {
  if ((RRG$aggregate != RG$aggregate)||(Nested)) {
    warning("Grouping filter by reference is not feasible if replicates are paired!")
    filter_types <- filter_types[which(filter_types != "ref")]
  } else {
    if (sum(vapply(RG$names, function(x) {! x %in% RSA$names }, TRUE)) > 0) {
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
  expContrasts_F$Contrasts <- tmp <- apply(expContrasts_F[, c("x1", "x0")], 1, function(x) {
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
  if (l > 1) {
    tmp2 <- unlist(vapply(1:(l-1), function(x) {
      paste0("(", tmp[x], ") - (", tmp[(x+1):l], ")")
    }, ""))
    tmp2Nm <- unlist(vapply(1:(l-1), function(x) {
      paste0("(", expContrasts_F$name[x], ") - (", expContrasts_F$name[(x+1):l], ")")
    }, ""))
    tmp2Tbl <- data.frame(Contrasts = tmp2,
                          x1 = gsub("^\\(|\\) - \\(.*", "", tmp2),
                          x0 = gsub(".*\\) - \\(|\\)$", "", tmp2),
                          name = tmp2Nm,
                          Type = "Double")
    tmp2Tbl$All <- lapply(1:nrow(tmp2Tbl), function(x) { gsub("^Group_", "", unlist(strsplit(unlist(tmp2Tbl[x, c("x1", "x0")]), " - "))) })
    tmp <- c(tmp, tmp2)
    expContrasts_F <- rbind(expContrasts_F, tmp2Tbl)
  }
}

#### Code chunk - Modified peptides analysis
Src <- paste0(libPath, "/extdata/R scripts/Sources/Cytoscape_init.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
modPepSrc <- paste0(libPath, "/extdata/R scripts/Sources/modPeptides.R")
#rstudioapi::documentOpen(modPepSrc)
source(modPepSrc, local = FALSE)

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
# It makes sense to close/re-create parallel clusters regularly to reduce memory usage + avoid corruption
stopCluster(parClust)
source(parSrc, local = FALSE)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - Create output tables
## PSMs
dir <- paste0(wd, "/Tables")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
w <- which(vapply(colnames(ev), function(x) { "list" %in% class(ev[[x]]) }, TRUE))
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
MakeRatios <- TRUE # (Used by the sourced, core sub-script)
intNms <- function(nms, topLvl = FALSE, type = "PG") {
  m <- match(type, c("pep", "PG"))
  root <- c("Intensity", "Expression")[m]
  mode <- topLvl+1
  sapply(nms, function(nm) {
    if (nm %in% c("Original", "Intensity", "Expression", "Original int.", "Intensity int.", "Expression int.")) {
      nm <- c(c("int.", "expr.")[m], root)[mode]
    } else {
      if (nm %in% c("ReNorm.", "ReNorm. int.")) { nm <- "re-norm" } else { nm <- substr(nm, 1, min(c(3, nchar(nm)))) }
      nm <- paste0(nm, ". ", c(c("int.", "expr.")[m], root)[mode])
    }
    paste0("log10(", nm, ")")
  })
}
ratNms <- function(nms, topLvl = FALSE) {
  mode <- topLvl+1
  sapply(nms, function(nm) {
    if (nm %in% c("Original", "Ratios", "Original rat.", "Ratios rat.")) {
      nm <- c("rat.", "Ratio")[mode]
    } else {
      if (nm %in% c("ReNorm.", "ReNorm. rat.")) { nm <- "Re-norm" } else { nm <- substr(nm, 1, min(c(3, nchar(nm)))) }
      nm <- paste0(nm, ". ", c("rat.", "ratios")[mode])
    }
    paste0("log2(", nm, ")")
  })
}
for (nm in names(pep.ref)) { #nm <- names(pep.ref[1])
  rpl <- intNms(nm, type = "pep")
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
fl <- system.file("extdata", "Report - column names - with replicates.xlsx", package = "proteoCraft")
styleNms <- openxlsx2::read_xlsx(fl, "tmp", colNames = FALSE)[,1]
WorkBook %<o% wb_load(fl)
repFl <- paste0(wd, "/Tables/Report_", dtstNm, ".xlsx")
WorkBook <- wb_add_data(WorkBook, "Description", dtstNm, wb_dims(2, 5))
WorkBook <- wb_add_data(WorkBook, "Description", format(Sys.Date(), "%d/%m/%Y"), wb_dims(3, 5))
WorkBook <- wb_add_data(WorkBook, "Description", WhoAmI, wb_dims(4, 5))
tmp <- loadedPackages(TRUE)
WorkBook <- wb_add_data(WorkBook, "Description", tmp$Version[grep("proteoCraft", tmp$Name)], wb_dims(5, 5))
WorkBook <- wb_set_base_font(WorkBook, 11, font_name = "Calibri")
cat(" - Writing Excel report...\n")
#
# Function for editing our header
KolEdit <- function(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl) {
  #KolNames <- ColumnsTbl$Col; intTbl = intColsTbl; ratTbl = ratColsTbl
  klnms <- KolNames
  KolNames <- gsub("Peptides?", "Pep.", KolNames)
  KolNames <- gsub("Evidences?", c("PSMs", "Ev.")[(SearchSoft == "MAXQUANT")+1], KolNames)
  KolNames <- gsub("Spectr((al)|(um))", "Spec.", KolNames)
  KolNames <- gsub("Razor", "Raz.", KolNames)
  KolNames <- gsub("Unique", "Uniq.", KolNames)
  KolNames <- gsub("MS/MS", "MS2", KolNames)
  # This would be the place to edit PER sample evidence and MS/MS columns,
  # currently not needed because those do not exist yet
  for (nm in names(intTbl)) { #nm <- names(intTbl)[1] #nm <- names(intTbl)[2]
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
    stopifnot(length(g) == 0)
  }
  #
  # Those names must be unique if the data is to be written as a table!
  # Which is annoying, because this limits how much fat we can cut
  tst <- aggregate(KolNames, list(KolNames), c)
  tst$L <- vapply(tst$x, length, 1)
  tst <- tst[which(tst$Group.1 != ""),]
  stopifnot(max(tst$L) == 1)
  #tst$x[which(tst$L > 1)]
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
  ev$"In list" <- c("", "+")[w+1]
}
QualFilt %<o% c("In list", "Potential contaminant", "Only identified by site")
if ((DiscFilt)&&(DiscFiltMode == DiscFiltModes[3])) { QualFilt <- c(QualFilt, DiscFiltCols) }
II <- setNames(1, "All peptidoforms")
if ((exists("PTMs_pep"))&&(length(PTMs_pep))) {
  stop("Really? There is no PTM-modified class of peptides to analyze? Why did you run this workflow then?")
} else {
  Mod2Write <- names(PTMs_pep)
  II[paste0(Mod2Write, "-mod. pept.")] <- 1+(seq_along(length(Mod2Write)))
}
for (ii in II) { #ii <- II[1] #ii <- II[2]
  tblMode <- tblMode2 <- "pep"
  TbNm <- names(II)[ii]
  tempData <- get(tblMode)
  intRf <- pep.ref
  ratRf <- pep.ratios.ref
  if (ii > 1) {
    Ptm <- Mod2Write[[ii-1]]
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
    if (ii > 1) { CoreCol <- c(CoreCol, "Name", paste0(Ptm, "-site(s)")) }
    CoreCol2 %<o% c("Leading razor proteins",
                    "Protein names", "Gene names", "Protein group IDs", "Razor protein group ID",
                    grep(" (Probabilities|Score Diffs)$", colnames(tempData), value = TRUE),
                    "Normalisation group")
    evcol <- "Evidence IDs"
    spcol <- "MS/MS count"
    gel <- setNames(lapply(intRf, function(rf) {
      x <- c(paste0(rf, RSA$values),
             paste0("Mean ", rf, VPAL$values))
      return(x[which(x %in% colnames(tempData))])
    }), intRf)
    if (ii == 1) {
      # Log transform for normal tables - not necessary for PTM-modified tables as we already transformed
      gel2 <- setNames(lapply(intRf, function(rf) {
        x <- c(paste0(rf, RSA$values),
               paste0("Mean ", rf, VPAL$values))
        w <- which(x %in% colnames(tempData))
        rf2 <- gsub("Evidence intensities - ", "log10(Int.) - ", rf)
        x <- c(paste0(rf2, RSA$values),
               paste0("Mean ", rf, VPAL$values))
        return(x[w])
      }), intRf)
      for (rf in intRf) { #rf <- intRf[1]
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
      intRf <- sapply(intRf, function(rf) {
        gsub("Evidence intensities - ", "log10(Int.) - ", rf)
      })
      gel <- unlist(gel2)
    }
    #
    smpls <- c(VPAL$values, RSA$values)
    if (length(Exp) == 1) { smpls <- gsub(topattern(paste0(Exp, "___")), "", smpls) }
    smpls <- gsub("___", " ", smpls)
    intColsTbl <- setNames(lapply(names(intRf), function(nm) { #nm <- names(intRf)[1]
      res <- data.frame(Log = c(paste0("Mean ", intRf[nm], VPAL$values),
                                paste0(intRf[nm], RSA$values)),
                        Type = c(rep("Average", length(VPAL$values)),
                                 rep("Individual", length(RSA$values))),
                        Sample = smpls)
      w <- which(res$Log %in% colnames(tempData))
      return(res[w,])
    }), names(intRf))
    w <- which(vapply(intColsTbl, nrow, 1) > 0)
    intColsTbl <- intColsTbl[w]; intRf <- intRf[w]
    quantCols <- intCols <- lapply(intColsTbl, function(x) { x$Log })
    ratColsTbl <- setNames(lapply(names(ratRf), function(nm) {
      res <- data.frame(Log = c(paste0("Mean ", ratRf[nm], VPAL$values),
                                paste0(ratRf[nm], RSA$values)),
                        Type = c(rep("Average", length(VPAL$values)),
                                 rep("Individual", length(RSA$values))),
                        Sample = smpls)
      w <- which(res$Log %in% colnames(tempData))
      return(res[w,])
    }), names(ratRf))
    w <- which(vapply(ratColsTbl, nrow, 1) > 0)
    ratColsTbl <- ratColsTbl[w]; ratRf <- ratRf[w]
    ratCols <- lapply(ratColsTbl, function(x) { x$Log })
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
    if (ii > 1) {
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
    if (ii > 1) { kol <- c(kol, "Code") }
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
    intColsTbl <- lapply(intColsTbl, function(x) {
      x$Log <- cleanNms(x$Log, start = FALSE)
      x
    })
    ratColsTbl <- lapply(ratColsTbl, function(x) {
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
    if ((ii > 1)&&(F.test)) {
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
    if (ii > 1) {
      # - Summary Ratios: P-values and significance
      ColumnsTbl[["P-values"]] <- gpl
      # - Significant
      ColumnsTbl[["Significant"]] <- signcol
      # - Regulated
      ColumnsTbl[["Regulated"]] <- regcol
      # F-test
      if ((ii > 1)&&(F.test)) {
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
      for (i in 1:nrow(AnnotTbl)) { ColumnsTbl[[paste0(AnnotTbl$Name[i], " annotations")]] <- AnnotTbl$Columns[[i]] }
    }
    # - PEP
    ColumnsTbl[["PEP"]] <- "PEP"
    # - Filters
    ColumnsTbl[["Filters"]] <- qualFlt
    # Melt
    ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }, 1) > 0)]
    ColumnsTbl <- set_colnames(reshape::melt.list(ColumnsTbl), c("Col", "Grp"))
    #tst <- aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), length); View(tst)
    #tst <- aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), unique); View(tst)
    #tst <- aggregate(1:nrow(ColumnsTbl), list(ColumnsTbl$Col), unique); w <- which(vapply(tst$x, length, 1) > 1); setNames(tst$x[w], tst$Group.1[w])
    stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
    ColumnsTbl$Class <- ""
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Peptides information"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evcol))] <- "Evidence IDs"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(spcol))] <- "Spectral count"
    for (nm in names(intRf)) { #nm <- names(intRf)[1]
      rpl <- intNms(nm, TRUE, type = "pep")
      ColumnsTbl$Class[grep(topattern(intRf[nm]), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[grep(topattern(paste0("Mean ", intRf[nm])), ColumnsTbl$Col)] <- rpl
    }
    for (nm in names(ratRf)) { #nm <- names(ratRf)[1]
      rpl <- ratNms(nm, TRUE)
      ColumnsTbl$Class[grep(topattern(ratRf[nm]), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[grep(topattern(paste0("Mean ", ratRf[nm])), ColumnsTbl$Col)] <- rpl
    }
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "P-values")] <- gsub(" - $", "", pvalue.col[which(pvalue.use)])
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% signcol)] <- "Significant"
    if ((ii > 1)&&(F.test)) {
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
    Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
}
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/Write_Excel_end_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#WorkBook$get_active_sheet()
#xl_open(repFl)

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)

#### Code chunk - Venn diagrams
Src <- paste0(libPath, "/extdata/R scripts/Sources/Venn_diagrams.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Finalize analysis and export results
# Remove empty directories:
#dirlist <- list.dirs()
dirlist <- dirlist[order(nchar(dirlist), decreasing = TRUE)]
for (dir in dirlist) { #d <- dirlist[1]
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
  MatMetCalls$Calls <- append(MatMetCalls$Calls, paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$DatAnalysis[", i,"], prop = WrdFrmt$",
                                                        c("Body", "Template_text")[(MatMetCalls$Texts$DatAnalysis[i] == "TEMPLATE")+1], "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#

# Write SDRF file in case you want to submit to PRIDE
Src <- paste0(libPath, "/extdata/R scripts/Sources/SDRF_4_PRIDE.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Finalize analysis
Src <- paste0(libPath, "/extdata/R scripts/Sources/Finalize_analysis.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

### That's it, done!
#openwd(outDir)
#rm(list = ls())
