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
scrptTypeFull %<o% "withReps_PG_and_PTMs"

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
xplorSrc %<o% paste0(libPath, "/extdata/R scripts/Sources/xplorData.R")

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
                     "tidyr", "ggplotify", "jpeg", "scattermore", "rpanel", "stringi", "lmtest", "ssh", "taxize", "arrow", "unimod",
                     "ggdendro", "colorspace", "factoextra", "NbClust", "BH", "plogr", "iq", "Rtsne"))
bioc_req <- unique(c(bioc_req, "biomaRt", "GO.db", "UniProt.ws", "limma", "sva", "qvalue", "MSnbase", "DEP",
                     "Rgraphviz", "RCy3", "siggenes", "DEqMS", "pRoloc", "pRolocGUI", "rbioapi", "png", "Rhdf5lib"))
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
ReportCalls <- AddTxt2Report("Processing PSMs...")
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
source(parSrc, local = FALSE)

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
source(parSrc, local = FALSE)

#### Code chunk - ROC analysis for optimizing a threshold to include/exclude peptides mapped to specific GO terms
Src <- paste0(libPath, "/extdata/R scripts/Sources/ROC1.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Create modified peptides table
Src <- paste0(libPath, "/extdata/R scripts/Sources/pepMake.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Peptidoforms-level, calculate quantitation and test for outliers
Src <- paste0(libPath, "/extdata/R scripts/Sources/pepQuant.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

####################################################
# Optional: choose whether to remove any outliers  #
####################################################
Src <- paste0(libPath, "/extdata/R scripts/Sources/remove_Outliers.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

if (length(inDirs) > 1) {
  paste0(DatAnalysisTxt, " PSM tables ")
  if (length(unique(SearchSoft)) > 1) {
    DatAnalysisTxt <- paste0(DatAnalysisTxt, "from the different search engines were converted to a similar format and ")
  } else {
    DatAnalysisTxt <- paste0(DatAnalysisTxt, "were ")
  }
  paste0(DatAnalysisTxt, "combined into a single table, then this ")
} else {
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " The long format PSMs table ")
}

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
DatAnalysisTxt <- paste0(DatAnalysisTxt, "was transformed into a wide format peptidoforms table, summing up quantitative values where necessary.")

LocAnalysis2 %<o% FALSE

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
# This code says: "Please rewrite me! I was written by a beginner (you) and I look very boorish!"
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
#poplot(plot, 12, 22)
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
if (("Ratios.Thresholds" %in% colnames(Param))&&(Param$Ratios.Thresholds == threshMsg)) {
  pep.Ref.Ratios %<o% make_RefRat()
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
# Legacy, I really do not think that this is a good idea at this stage if intensities have been well normalized.
# Even if we are dealing with SILAC (not supported yet) where ratios are better measured and should be more precise,
# I would simply re-calculate updated intensities to reflect the ratios at an early stage, then focus on intensities from there.
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
      #DatAnalysisTxt <- paste0(DatAnalysisTxt, " Peptide ratios were then re-normalized.")
    } else { warning("Nothing to plot for Reference-to-Reference ratios!") }
  }
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
# It makes sense to close/re-create parallel clusters regularly to reduce memory usage
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
stopCluster(parClust)
source(parSrc, local = FALSE)

#### Code chunk - Assemble protein groups
ReportCalls <- AddSpace2Report()
ReportCalls <- AddTxt2Report("Starting protein groups assembly...")
if ("Prot.Only.with.at.least" %in% colnames(Param)) {
  NPep <- as.integer(Param$Prot.Only.with.at.least)
  if ((is.na(NPep))||(NPep <= 0)) {
    warning("Invalid `\"Prot.Only.with.at.least\" parameter, defaulting to 1!")
    NPep <- 1
  }
} else { NPep <- 1 }
.obj <- unique(c("NPep", .obj))
#
tm1 <- Sys.time()
source(parSrc, local = FALSE)
Src <- paste0(libPath, "/extdata/R scripts/Sources/PG_assemble.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#loadFun("PG_assembly.RData")
tm2 <- Sys.time()
#
PG %<o% PG_assembly$Protein.groups
pep <- PG_assembly$Peptides
db <- PG_assembly$Database
if ("Evidences" %in% names(PG_assembly)) { ev <- PG_assembly$Evidences }
msg <- paste0(nrow(PG), " protein groups assembled in ", gsub("^Time difference of ", "", capture.output(tm2-tm1)))
ReportCalls <- AddMsg2Report(Space = FALSE, Print = FALSE)



# Here would be a good place to check protein taxonomy and [if necessary/as per parameters] split hybrid groups!
# Should be controlled by a parameter only showing up if taxonomy is present in db and has more than one value!
warning("(TO DO: add 'split-by-taxonomy?' here!)")
# Default = TRUE
# Don't forget to update PG IDs in PG, ev and pep afterwards! Also check potential contaminant column!



# Check those rare proteins IDs which are not in the search DB (should be contaminants, there is a minor inconsistency in how they are )
tst <- unlist(strsplit(pep$Proteins, ";"))
if (length(tst)) {
  msg <- paste0("These protein accessions in peptides are not in the database: ", paste(tst[which(!tst %in% db$`Protein ID`)], collapse = " - "))
  ReportCalls <- AddMsg2Report(Space = FALSE, Print = FALSE)
}

# Basic fix, because I do not like the way I was doing Quality filters up to now
g <- grep("^Quality filter: ", colnames(PG), value = TRUE)
if (length(g)) {
  for (h in g) { #h <- g[1]
    PG[[h]] <- c("no -> dubious!", "")[match(PG[[h]], c("", "Keep"))]
  }
}
#
if (tstOrg) {
  test <- vapply(strsplit(PG$`Protein IDs`, ";"), function(x) {
    paste(sort(unique(c(db[match(x, db$`Protein ID`), dbOrgKol]))), collapse = ";")
  }, "")
  pgOrgKol %<o% c("Organism", "Organism(s)")[(sum(grepl(";", test))>0)+1]
  PG[[pgOrgKol]] <- test
}
DatAnalysisTxt <- paste0(DatAnalysisTxt, " Protein groups were inferred from observed peptides.")

# Some stats on protein groups
tmp <- aggregate(PG$id, list(PG$`Peptides count`), length)
colnames(tmp) <- c("Peptides count", "Protein groups")
tmp$"log10(Protein groups count)" <- log10(tmp$"Protein groups")
pal <- colorRampPalette(c("brown", "yellow"))(max(tmp$"Peptides count")-1)
tmp$Colour <- c("blue", pal)[tmp$`Peptides count`]
tmp2 <- summary(PG$`Peptides count`)
tmp2 <- data.frame(Variable = c(names(tmp2), "", "Protein groups", "Protein groups with 2+ peptidoforms"),
                   Value = c(as.character(signif(as.numeric(tmp2), 3)),
                             "",
                             as.character(c(nrow(PG), sum(PG$"Peptides count" >= 2)))))
tmp2$Txt <- apply(tmp2[, c("Variable", "Value")], 1, function(x) {
  x <- x[which(x != "")]
  if (length(x)) { x <- paste(x, collapse = ": ") } else { x <- "" }
  return(x)
})
tmp2$X <- max(as.numeric(tmp2$Value[match("Max.", tmp2$Variable)]))*0.98
tmp2$Y <- max(tmp$"log10(Protein groups count)")*(0.98-(0:(nrow(tmp2)-1))*0.02)
ttl <- "Peptidoforms per PG"
plot <- ggplot(tmp) + geom_col(aes(x = `Peptides count`, y = `log10(Protein groups count)`, fill = Colour),
                               colour = NA) +
  geom_text(data = tmp2, aes(x = X, y = Y, label = Txt), hjust = 1, size = 3) +
  scale_fill_identity() + theme_bw() + ggtitle(ttl)
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
dir <- paste0(wd, "/Summary plots")
dirlist<- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300)
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)
})

source(parSrc, local = FALSE)
tmp1 <- strsplit(pep$Proteins, ";")
tmp2 <- db[which(db$"Protein ID" %in% unlist(tmp1)), c("Protein ID", "Gene", "Name", "Common Name")]
exports <- list("tmp1", "tmp2")
clusterExport(parClust, "tmp2", envir = environment())
tmp <- parSapply(parClust, tmp1, function(x) {
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
  if (i == "No Isoforms") { j <- i } else { j <- gsub("s$", "", i) }
  if (!j %in% colnames(db)) {
    j <- paste0(gsub("s$", "", i), c(" ID", " IDs", ""))
    w <- which(j %in% colnames(db))
    if (length(w)) { j <- j[w[1]] } else {
      warning(paste0("No near matching column name found for \"", i, "\" in the protein data base table."))
    }
  }
  if (length(j) == 1) {
    tmp4 <- db[[j]]
    exports <- list("i", "tmp4")
    clusterExport(parClust, exports, envir = environment())
    PG[[i]] <- parSapply(parClust, tmp, function(x) {
      x <- unlist(x)
      m1 <- match(x, tmp2)
      if (tstFllID) {
        m1 <- data.frame(m1 = m1, m2 = match(x, tmp3))
        m1 <- apply(m1, 1, function(y) {
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
if (length(w) == 2) {
  temp <- PG[, genkol]
  for (i in genkol) { temp[[i]] <- strsplit(temp[[i]], ";") }
  PG$Genes <- apply(temp, 1, function(x) { paste(sort(unique(unlist(x))), collapse = ";") })
  PG$"Gene names" <- NULL
} else { if (length(w) == 1) { colnames(PG)[which(colnames(PG) %in% genkol)] <- "Genes" } }
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
#
# If Arabidopsis:
if (("ARATH" %in% db$Organism)||(3702 %in% db$TaxID)) {
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
  if (sum(kol %in% colnames(db)) == 2) {
    #View(db[, kol])
    w <- which((nchar(db$TAIR) > 0)&(nchar(db$TAIR_v2) > 0))
    tmp4 <- db[w, kol]
    tmp4$TAIR <- strsplit(tmp4$TAIR, ";")
    tmp4$TAIR_v2 <- strsplit(tmp4$TAIR_v2, ";")
    tmp4 <- apply(tmp4, 1, unique)
    tmp4 <- parSapply(parClust, tmp4, function(x) { paste(unlist(x), collapse = ";") })
    db$TAIR[w] <- tmp4
    w <- which(db$TAIR == "")
    db$TAIR[w] <- db$TAIR_v2[w]
    db$TAIR_v2 <- NULL
    db$TAIR[which(is.na(db$TAIR))] <- ""
  }
  tmp2 <- listMelt(strsplit(PG$`Leading protein IDs`, ";"), PG$id)
  tmp2$TAIR <- db$TAIR[match(tmp2$value, db$`Protein ID`)]
  tmp3 <- listMelt(strsplit(tmp2$TAIR, ";"), tmp2$L1)
  tmp3 <- aggregate(tmp3$value, list(tmp3$L1), function(x) { paste(unique(x), collapse = ";") })
  PG$TAIR <- ""
  w <- which(PG$id %in% tmp3$Group.1)
  PG$TAIR[w] <- tmp3$x[match(PG$id[w], tmp3$Group.1)]
}

# Here, if this is a BioID type experiment, we also want to mark protein groups which have Biotin peptides:
IsBioID2 %<o% FALSE
if (IsBioID) {
  wbiot %<o% grep("biot", Modifs$"Full name", ignore.case = TRUE)
  l <- length(wbiot)
  if (length(wbiot)) {
    if (l == 1) { tmp <- Modifs$"Full name"[wbiot] } else {
      tmp <- paste0(paste(Modifs$"Full name"[wbiot[1:(l-1)]], collapse = "\", \""), "\" and \"", Modifs$"Full name"[wbiot[l]])
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
      temp <- aggregate(temp$value, list(temp$L1), function(x) { paste(sort(x), collapse = ";") })
      PG[wpg, "Biot. peptide IDs"] <- temp$x[match(PG$id[wpg], temp$Group.1)]
      PG[["Biot. peptides count"]] <- vapply(strsplit(PG[["Biot. peptide IDs"]], ";"), length, 1)
      PG[["Biot. peptides [%]"]] <- round(100*PG[["Biot. peptides count"]]/PG$"Peptides count", 1)
      IsBioID2 <- TRUE
    } else { warning("I could not find any biotinylated peptides!") }
  } else { warning("I could not identify any biotinylated PTMs in the modifications table, did you include them in the search?") }
}
# First sequence
m <- match(vapply(strsplit(PG$"Leading protein IDs", ";"), function(x) { x[[1]] }, ""),
           db$"Protein ID")
PG$"Protein ID (1st accession)" <- db$`Protein ID`[m]
PG$"Sequence (1st accession)" <- db$Sequence[m]

# Number of spectra, PSMs and peptides per sample:
source(parSrc, local = FALSE)
invisible(clusterCall(parClust, function() {
  library(proteoCraft)
  library(reshape)
  library(data.table)
  return()
}))
temp_PG <- data.frame(id = PG$id,
                      Accession1 = PG$`Protein ID (1st accession)`)
temp_PG$Pep <- parLapply(parClust, strsplit(PG$"Peptide IDs", ";"), as.integer)
tmp <- pep[, c("id", "Sequence")]
clusterExport(parClust, "tmp", envir = environment())
temp_PG$Pep <- parLapply(parClust, temp_PG$Pep, function(x) { tmp$Sequence[match(x, tmp$id)] })
temp_PG$Seq <- db$Sequence[match(temp_PG$Accession1, db$"Protein ID")]
exports <-
if (!"Sequence coverage [%]" %in% colnames(PG)) {
  exports <- list("Coverage")
  clusterExport(parClust, exports, envir = environment())
  PG$"Sequence coverage [%]" <- round(100*parApply(parClust, temp_PG[, c("Seq", "Pep")], 1, function(x) {
    Coverage(x[[1]], x[[2]])
  }), 1)
}
CreateMSMSKol %<o% (("MS/MS IDs" %in% colnames(ev))&&(class(ev$"MS/MS IDs") %in% c("integer", "character")))
if (CreateMSMSKol) {
  # There appear to be no MSMS IDs for DIA in MaxQuant.
  ev$temp <- parLapply(parClust, strsplit(as.character(ev$"MS/MS IDs"), ";"), as.integer)
  #PG[, paste0("Spectr", c("al count", "um IDs"))]
  temp <- listMelt(lapply(strsplit(PG$`Evidence IDs`, ";"), as.integer), PG$id, c("Ev_id", "PG_id"))
  temp$MSMSIDs <- ev$temp[match(temp$Ev_id, ev$id)]
  temp <- temp[which(vapply(temp$MSMSIDs, length, 1) > 0),] # Remove Match-Between-Runs evidences (no MS/MS)
  temp <- listMelt(temp$MSMSIDs, temp$PG_id, c("MSMSIDs", "PG_id"))
  temp <- do.call(data.frame, aggregate(temp$MSMSIDs, list(temp$PG_id), function(x) {
    x <- unique(x)
    return(c(Count = length(x), List = list(x)))
  }))
  temp$x.Count <- unlist(temp$x.Count)
  temp$Pasted <- vapply(temp$x.List, function(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") }, "")
  PG[, paste0("Spectr", c("al count", "um IDs"))] <- temp[match(PG$id, temp$Group.1), c("x.Count", "Pasted")]
  ev$temp <- NULL
}
temp_ev <- ev[, c("id", "MQ.Exp", "Protein group IDs", "Peptide ID")]
if (CreateMSMSKol) { temp_ev$"MS/MS IDs" <- ev$"MS/MS IDs" }
temp_pep <- pep[, c("id", "Sequence")]
clusterExport(parClust, exports, envir = environment())
Samplez <- list(Sample = c(), Group = c())
MQ_Exp <- list()
for (gr in VPAL$values) { #gr <- VPAL$values[1]
  wh <- which(Exp.map[[VPAL$column]] == gr)
  smplz <- Exp.map$Ref.Sample.Aggregate[wh]
  mqexp <- setNames(lapply(smplz, function(x) { Exp.map$MQ.Exp[match(x, Exp.map$Ref.Sample.Aggregate)] }), smplz)
  mqexp[[gr]] <- Exp.map$MQ.Exp[wh]
  Samplez$Sample <- c(Samplez$Sample, smplz, gr)
  Samplez$Group <- c(Samplez$Group, rep(gr, length(smplz)), gr)
  MQ_Exp[[gr]] <- mqexp
}
Samplez <- data.frame(Sample = Samplez$Sample,
                      Group = Samplez$Group)
exports <- list("Exp.map", "temp_ev", "temp_pep", "temp_PG", "IsBioID2", "MQ_Exp", "Samplez", "Modifs", "CreateMSMSKol")
if (IsBioID2) { exports <- append(exports, "wbiot") }
clusterExport(parClust, exports, envir = environment())
temp <- parApply(parClust, Samplez, 1, function(Smpl) { #Smpl <- unlist(Samplez[1,])
  smpl <- Smpl[[1]]
  gr <- Smpl[[2]]
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
  res[, kolk] <- 0
  res[, koli] <- ""
  res[[paste0("Sequence coverage [%] - ", smpl)]] <- 0
  mqexp <- MQ_Exp[[gr]]
  w <- which(temp_ev$MQ.Exp %in% unique(unlist(mqexp[[smpl]])))
  if (length(w)) {
    e <- temp_ev[w, , drop = FALSE]
    temp1 <- lapply(strsplit(e$"Protein group IDs", ";"), as.integer)
    temp1 <- listMelt(temp1, e$"Peptide ID")
    temp1 <- do.call(data.frame, aggregate(temp1$L1, list(temp1$value), function(x) {
      x <- unique(x)
      return(c(Count = length(x), List = list(x)))
    }))
    temp1$x.Count <- unlist(temp1$x.Count)
    temp1$Pasted <- vapply(temp1$x.List, function(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") }, "")
    tmp <- temp1$x.List
    tmp <- listMelt(tmp)
    tmp$Seq <- temp_pep$Sequence[match(as.numeric(tmp$value), temp_pep$id)]
    tmp <- data.table(Seq = tmp$Seq, Row = tmp$L1)
    tmp <- tmp[, list(Seq = list(unique(Seq))), by = Row]
    tmp <- as.data.frame(tmp)
    temp1$Pepseq <- tmp$Seq[match(tmp$Row, 1:nrow(temp1))]
    temp2 <- lapply(strsplit(e$"Protein group IDs", ";"), as.integer)
    temp2 <- listMelt(temp2, e$id)
    temp2 <- do.call(data.frame, aggregate(temp2$L1, list(temp2$value), function(x) {
      c(Count = length(x), IDs = paste(sort(as.numeric(x)), collapse = ";"))
    }))
    temp2$x.Count <- as.integer(temp2$x.Count)
    if (CreateMSMSKol) {
      w3 <- which(e$"MS/MS IDs" != "")
      temp3 <- lapply(strsplit(e$"Protein group IDs"[w3], ";"), as.integer)
      temp3 <- listMelt(temp3, e$"MS/MS IDs"[w3])
      temp3$L1 <- lapply(strsplit(temp3$L1, ";"), as.integer)
      temp3 <- listMelt(temp3$L1, temp3$value)
      temp3 <- do.call(data.frame, aggregate(temp3$value, list(temp3$L1), function(x) {
        x <- unique(x)
        return(c(Count = length(x), List = list(x)))
      }))
      temp3$x.Count <- unlist(temp3$x.Count)
      temp3$Pasted <- vapply(temp3$x.List, function(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") }, "")
      w <- which(res$id %in% temp3$Group.1)
      m <- match(res$id[w], temp3$Group.1)
      res[w, kols] <- temp3[m, c("x.Count", "Pasted")]
    }
    w <- which(res$id %in% temp1$Group.1)
    m <- match(res$id[w], temp1$Group.1)
    res[w, kolp] <- temp1[m, c("x.Count", "Pasted")]
    temp_PG$Pep <- NA
    temp_PG$Pep[w] <- temp1$Pepseq[m]
    res[w, paste0("Sequence coverage [%] - ", smpl)] <- round(100*apply(temp_PG[w, c("Seq", "Pep")], 1, function(x) {
      Coverage(x[[1]], x[[2]])
    }), 1)
    w <- which(res$id %in% temp2$Group.1)
    m <- match(res$id[w], temp2$Group.1)
    res[w, kole] <- temp2[m, c("x.Count", "x.IDs")]
    if (IsBioID2) {
      kolB <- paste0("Biot. ", paste0(tolower(substr(kol, 1, 1)), substr(kol, 2, nchar(kol))))
      kolBe <- grep("^Biot\\. evidence", kolB, value = TRUE)
      kolBp <- grep("^Biot\\. peptide", kolB, value = TRUE)
      if (CreateMSMSKol) {
        kolBs <- grep("^Biot\\. spectr", kolB, value = TRUE)
        kolB <- kolB[which(!kolB %in% kolBs)]; rm(kolBs) #I don't think we need those columns now... too many is too many
      }
      kolBk <- grep(" count - ", kolB, value = TRUE)
      kolBi <- grep(" IDs - ", kolB, value = TRUE)
      res[, kolBk] <- 0
      res[, kolBi] <- ""
      g <- grep(topattern(Modifs$Mark[wbiot], start = FALSE), e$"Modified sequence")
      if (length(g)) {
        eB <- e[g, , drop = FALSE]
        temp1 <- listMelt(lapply(strsplit(eB$"Protein group IDs", ";"), as.integer), eB$"Peptide ID")
        temp1 <- do.call(data.frame, aggregate(temp1$L1, list(temp1$value), function(x) {
          x <- unique(x)
          return(c(Count = length(x), IDs = paste(sort(as.numeric(x)), collapse = ";")))
        }))
        temp1$x.Count <- as.integer(temp1$x.Count)
        temp1 <- do.call(data.frame, temp1)
        temp2 <- strsplit(eB$"Protein group IDs", ";")
        temp2 <- listMelt(temp2, eB$id)
        temp2 <- do.call(data.frame, aggregate(temp2$L1, list(temp2$value), function(x) {
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
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
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
  g <- grsep2(prot.list, PG$"Protein IDs")
  PG$`In list`[g] <- "+"
  PG$"Potential contaminant"[g] <- ""
  ev$"Potential contaminant"[grsep2(prot.list, ev$Proteins)] <- ""
  pep$"Potential contaminant"[grsep2(prot.list, pep$Proteins)] <- ""
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Calculate protein group-level quantitative values
# You may want to exclude peptides based on how many samples they are found in:
kol <- lapply(VPAL$values, function(x) {
  x <- paste0(pep.ref["Original"], Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)])
  return(x[which(x %in% colnames(pep))])
})
if (!exists("PepFoundInAtLeast")) { PepFoundInAtLeast %<o% 1 }
if ("PepFoundInAtLeast" %in% colnames(Param)) {
  PepFoundInAtLeast %<o% suppressWarnings(as.integer(Param$PepFoundInAtLeast))
  if ((is.na(PepFoundInAtLeast))||(PepFoundInAtLeast < 1)) {
    warning("Invalid \"PepFoundInAtLeast\" parameter, defaulting to 1")
    PepFoundInAtLeast <- 1
  }
}
tst <- apply(pep[, unlist(kol)], 1, function(x) { sum(proteoCraft::is.all.good(x) > 0) })
Pep2Use %<o% which(tst >= PepFoundInAtLeast)
maxAllowed <- max(c(2, length(Rep)-1))
if (!exists("PepFoundInAtLeastGrp")) { PepFoundInAtLeastGrp %<o% maxAllowed }
if ("PepFoundInAtLeastGrp" %in% colnames(Param)) {
  PepFoundInAtLeastGrp %<o% suppressWarnings(as.integer(Param$PepFoundInAtLeastGrp))
  if ((is.na(PepFoundInAtLeastGrp))||(PepFoundInAtLeastGrp < 1)||(PepFoundInAtLeastGrp > length(Rep)-1)) {
    warning(paste0("Invalid \"PepFoundInAtLeastGrp\" parameter, defaulting to ", maxAllowed))
    PepFoundInAtLeastGrp <- maxAllowed
  }
}
tst <- pep[, unlist(kol)]
clusterExport(parClust, list("tst", "kol", "PepFoundInAtLeastGrp"), envir = environment())
tst <- parSapply(parClust, 1:nrow(pep), function(x) {
  x <- max(vapply(kol, function(kl) { sum(proteoCraft::is.all.good(unlist(tst[x, kl])) > 0) }, 1) >= PepFoundInAtLeastGrp)
  return(x)
}) > 0
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Pep2Use %<o% which(tst)
#length(Pep2Use)/nrow(pep)

# Currently 9 methods are implemented, but the last 3 are just for testing and do not create all columns required for the script to complete
# (These could be added relatively easily if necessary though)
QuantData <- setNames(paste0("quant.data", seq_along(QuantMethods)), QuantMethods)
QuantMethods_all <- FALSE
.obj <- unique(c("QuantData", "QuantMethods_all", .obj))
exprsCol <- paste0("log10(Expr.) - ", RSA$values)
# Weights:
# - Higher for peptides with low intra-sample group CV on average
# - Higher for peptides with low PEP
if ((length(inDirs) == 1)&&(QuantUMS)) {
  weightsInsrt <- "mean DiaNN \"Quantity Quality\"" 
  pep$Weights <- pep$"Quantity Quality"
} else {
  weightsInsrt <- "-log(PEP)/CV" 
  source(parSrc, local = FALSE)
  Kols <- paste0(pep.ref[length(pep.ref)], Exp.map$Ref.Sample.Aggregate)
  Kols <- Kols[which(Kols %in% colnames(pep))]
  tmp <- pep[, Kols]
  clusterExport(parClust, list("Exp.map", "VPAL", "pep.ref", "is.all.good", "tmp"), envir = environment())
  CV <- parSapply(parClust, VPAL$values, function(x) { #x <- VPAL$values[1]
    smpls <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)]
    kols <- paste0(pep.ref[length(pep.ref)], smpls)
    kols <- kols[which(kols %in% colnames(tmp))]
    x <- apply(tmp[, kols, drop = FALSE], 1, function(y) {
      y <- is.all.good(log10(unlist(y)))
      if (length(y)) {
        y <- sd(y)/mean(y)
      } else { y <- NA }
      return(y)
    })
    return(x)
  })
  CV <- rowMeans(CV, na.rm = TRUE)
  pep$Weights <- -log10(pep$PEP)/CV
}
summary(pep$Weights)
m <- max(is.all.good(pep$Weights))
pep$Weights <- pep$Weights/m
pep$Weights[which((is.na(pep$Weights))|(pep$Weights < 0.001))] <- 0.001
if ("Prot.Quant.Mod.Excl.is.strict" %in% colnames(Param)) {
  Mod.Excl.is.strict <- Param$Prot.Quant.Mod.Excl.is.strict
  if (!Mod.Excl.is.strict %in% c(1, 0, TRUE, FALSE)) { Mod.Excl.is.strict <- FALSE }
} else { Mod.Excl.is.strict <- FALSE }
Discard.unmod <- Mod.Excl.is.strict+1
if (Discard.unmod == 1) { Discard.unmod <- as.logical(Discard.unmod) }
.obj <- unique(c("Mod.Excl.is.strict", "Discard.unmod", .obj))
if (!grepl("^Prot\\.Quant", Param$QuantMeth)) {
  stop("NB: currently only methods 1 to 6 provide all necessary columns, not just quantitative columns per sample but also reference columns, ratios, etc... Until those are added they cannot be used by this script. Defaulting to method 3!")
  Param$QuantMeth <- "Prot.Quant.Unique"
}
if ((Param$QuantMeth %in% c("Prot.Quant", "Prot.Quant + weights", "Prot.Quant.Unique", "Prot.Quant.Unique + weights"))||(QuantMethods_all)) {
  # Calculates individual protein group expression and ratio values per sample, excluding some modifications and
  # (usually) requiring at least 2 peptides with different modified sequence.
  # This code has on occasion failed for very large datasets with the following error:
  # "Error in serialize(data, node$con) : error writing to connection"
  # Admittedly, when I write code, I tend to favour complexity over speed.
  # To avoid this issue, I added a bit of code to estimate how many threads to use:
  # - Size of clusters should not exceed N-1 threads (N = the computer's number of vCPUs).
  # - Somehow the LFQ function seems to require less than 7x the size of the pep & PG data frames in RAM
  #   per vCPU, regardless of number of threads.
  #   This is true even after I removed useless variable duplicates from the code!
  #   This is at least true for one project (ADAPTED project E), not sure yet about other datasets.
  #
  # This function is adapted from: https://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine
  # The sections for Linux and Mac-OS are untested, guesses based on: https://apple.stackexchange.com/questions/4286/is-there-a-mac-os-x-terminal-version-of-the-free-command-in-linux-systems
  # They are highly likely to need some adjustments before they work, but I have no way to test them!!!
  get_free_ram <- function() {
    sysnm <- tolower(Sys.info()[["sysname"]])
    if (sysnm == "windows") {
      x <- try(system2("wmic", args =  "OS get FreePhysicalMemory /value", stdout = TRUE), silent = TRUE)
      if (!"try-error" %in% class(x)) {
        x <- as.numeric(gsub("^FreePhysicalMemory=|\r$", "", grep("FreePhysicalMemory", x, value = TRUE)))/1000000
      } else {
        x <- try(system2("systeminfo.exe", stdout = TRUE), silent = TRUE)
        if (!"try-error" %in% class(x)) {
          x <- grep("^Total Physical Memory:     ", x, value = TRUE)
          x <- gsub("^Total Physical Memory:     |,", "", x)
          if (grepl(" B$", x)) { x <- as.numeric(gsub(" KB$", "", x))/1024^3 }
          if (grepl(" KB$", x)) { x <- as.numeric(gsub(" KB$", "", x))/1024^2 }
          if (grepl(" MB$", x)) { x <- as.numeric(gsub(" MB$", "", x))/1024 } # The only case I expect
          if (grepl(" GB$", x)) { x <- as.numeric(gsub(" GB$", "", x)) }
          if (grepl(" TB$", x)) { x <- as.numeric(gsub(" TB$", "", x))*1024 }
        } else {
          x <- try(system2("powershell", args = c("Get-WmiObject", "-Class", "Win32_ComputerSystem"),
                           stdout = TRUE), silent = TRUE)
          if (!"try-error" %in% class(x)) {
            x <- grep("^TotalPhysicalMemory : ", x, value = TRUE)
            x <- as.numeric(gsub("^TotalPhysicalMemory : ", "", x))/1024^3
          } else {
            # (If I really cannot detect RAM by any method, I will have to trust in the Omnissiah that I have enough...)
            x <- 10000
          }
        }
      }
    }
    if (sysnm %in% c("mac", "macos")) {
      warning("This code has not been tested and is likely to fail! Check it the first time you run on a Mac-OS system!")
      x <- try(system2("$ vm_stat", stdout = TRUE), silent = TRUE)
      if (!"try-error" %in% class(x)) {
        ps <- as.integer(gsub(".+\\(page size of | bytes\\)", "",x[1]))
        x <- as.numeric(gsub("^Pages free: +|\\.$", "", x[2]))*ps/1000000
      } else {
        x <- 10000 # If I cannot detect RAM, I will have to trust in the Omnissiah that I have enough...
      }
    }
    if (sysnm == "linux") {
      warning("This code has not been tested and is likely to fail! Check it the first time you run on a Linux system!")
      x <- system2("$free", stdout = TRUE)
      x <- as.integer(unlist(strsplit(x[1], " +"))[4])/1000000
    }
    if (!sysnm %in% c("windows", "mac", "macos", "linux")) { stop("Unknown computer system type!") }
    return(x)
  }
  free_ram <- get_free_ram()
  mem_per_thread <- as.numeric(gsub(" Gb$", "", (format(object.size(pep) + object.size(PG), "Gb"))))*3
  N.clust2 <- N.clust
  while (N.clust2*mem_per_thread > free_ram) { N.clust2 <- N.clust2-1 }
  tmpMod2Xclud <- Modifs$`Full name`[match(Mod2Xclud$Mark, Modifs$Mark)]
  if ((Param$QuantMeth == "Prot.Quant")||(QuantMethods_all)) {
    # Classic profiles method, no weights
    DatAnalysisTxt <- gsub("\\.$",
                           paste0(", and quantified using an in-house algorithm which: ",
                                  "i) computes a mean protein-level profile across samples using individual, normalized peptidoform profiles (\"relative quantitation\" step), ",
                                  "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                  "only unique and razor peptidoforms were used; ",
                                  tmpMod2Xclud, " peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                  collapse = " "), DatAnalysisTxt)
    source(parSrc, local = FALSE)
    quant.data1 <- Prot.Quant(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep[Pep2Use,], id = "id",
                              Summary.method = "mean", #Summary.weights = "Weights",
                              Intensity.weights = FALSE,
                              experiments.map = Exp.map, param = Param,
                              Pep.Intens.root = pep.ref[length(pep.ref)],
                              Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                              log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                              Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                              Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                              Discard.unmod = Discard.unmod,
                              Min.N = 1, Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              Refs_Mode = RefRat_Mode,
                              cl = parClust)
    #colnames(quant.data1)
    saveFun(quant.data1, file = "quant.data1.RData")
    #loadFun("quant.data1.RData")
  }
  if ((Param$QuantMeth == "Prot.Quant + weights")||(QuantMethods_all)) {
    # Classic profiles method, weights
    DatAnalysisTxt <- gsub("\\.$",
                           paste0(", and quantified using an in-house algorithm which: ",
                                  "i) computes a protein group-level profile across samples as the weighted mean of individual, normalized peptidoform profiles (\"relative quantitation\" step, weights = ",
                                  weightsInsrt,
                                  ", where PEP and CV = the peptidoform's Posterior Error Probability and Coefficient of Variation, resp.), ",
                                  "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                  "only unique and razor peptidoforms were used; ",
                                  tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                  collapse = " "), DatAnalysisTxt)
    source(parSrc, local = FALSE)
    quant.data2 <- Prot.Quant(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep[Pep2Use,], id = "id",
                              Summary.method = "weighted.mean", Summary.weights = "Weights",
                              Intensity.weights = FALSE,
                              experiments.map = Exp.map, param = Param,
                              Pep.Intens.root = pep.ref[length(pep.ref)],
                              Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                              log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                              Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                              Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                              Discard.unmod = Discard.unmod,
                              Min.N = 1, Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              Refs_Mode = RefRat_Mode,
                              cl = parClust)
    #colnames(quant.data2)
    saveFun(quant.data2, file = "quant.data2.RData")
    #loadFun("quant.data2.RData")
  }
  if ((Param$QuantMeth == "Prot.Quant.Unique")||(QuantMethods_all)) {
    # Profiles method, prefer unique, no weights
    DatAnalysisTxt <- gsub("\\.$",
                           paste0(", and quantified using an in-house algorithm which: ",
                                  "i) computes a mean protein-level profile across samples using individual, normalized peptidoform profiles (\"relative quantitation\" step), ",
                                  "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                  "for protein groups with at least 3 unique peptidoforms, only unique ones were used, otherwise razor peptidoforms were also included; ",
                                  tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                  collapse = " "), DatAnalysisTxt)
    source(parSrc, local = FALSE)
    quant.data3 <- Prot.Quant(Prot = PG, Mode = "PreferUnique",
                              Peptide.IDs = "Razor peptide IDs", Unique.peptide.IDs = "Unique peptide IDs",
                              Pep = pep[Pep2Use,], id = "id",
                              Summary.method = "mean", #Summary.weights = "Weights",
                              Intensity.weights = FALSE,
                              experiments.map = Exp.map, param = Param,
                              Pep.Intens.root = pep.ref[length(pep.ref)],
                              Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                              log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                              Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                              Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                              Discard.unmod = Discard.unmod,
                              Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              Refs_Mode = RefRat_Mode,
                              cl = parClust)
    #colnames(quant.data3)
    saveFun(quant.data3, file = "quant.data3.RData")
    #loadFun("quant.data3.RData")
  }
  if ((Param$QuantMeth == "Prot.Quant.Unique + weights")||(QuantMethods_all)) {
    # Profiles method, prefer unique, weights
    DatAnalysisTxt <- gsub("\\.$",
                           paste0(", and quantified using an in-house algorithm which: ",
                                  "i) computes a protein-level profile across samples as the weighted mean of individual, normalized peptidoform profiles (\"relative quantitation\" step, weights = ",
                                  weightsInsrt,
                                  ", where PEP and CV = the peptidoform's Posterior Error Probability and Coefficient of Variation, resp.), ",
                                  "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                  "for protein groups with at least 3 unique peptidoforms, only unique ones were used, otherwise razor peptidoforms were also included; ",
                                  tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                  collapse = " "), DatAnalysisTxt)
    source(parSrc, local = FALSE)
    quant.data4 <- Prot.Quant(Prot = PG, Mode = "PreferUnique",
                              Peptide.IDs = "Razor peptide IDs", Unique.peptide.IDs = "Unique peptide IDs",
                              Pep = pep[Pep2Use,], id = "id",
                              Summary.method = "weighted.mean", Summary.weights = "Weights",
                              Intensity.weights = FALSE,
                              experiments.map = Exp.map, param = Param,
                              Pep.Intens.root = pep.ref[length(pep.ref)],
                              Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                              log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                              Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                              Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                              Discard.unmod = Discard.unmod,
                              Min.N = 1, Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              Refs_Mode = RefRat_Mode,
                              cl = parClust)
    #colnames(quant.data4)
    saveFun(quant.data4, file = "quant.data4.RData")
    #loadFun("quant.data4.RData")
  }
}
if ((Param$QuantMeth == "Prot.Quant2 + weights")||(QuantMethods_all)) {
  # Averaging method, weights
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using an in-house algorithm which: ",
                                "i) computes protein group expression values as the mean of peptidoform intensity values (\"relative quantitation\" step), then ",
                                # "i) computes protein group expression values as the weighted mean of peptidoform intensity values (\"relative quantitation\" step, weights = ",
                                weightsInsrt,
                                ", where PEP and CV = the peptidoform's Posterior Error Probability and Coefficient of Variation, resp.), then ",
                                "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                "only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, " peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
  quant.data5 <- Prot.Quant2(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep[Pep2Use,], id = "id",
                             Summary.method = "weighted.mean", Summary.weights = "Weights",
                             experiments.map = Exp.map, param = Param,
                             Pep.Intens.root = pep.ref[length(pep.ref)],
                             Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                             log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                             Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                             Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                             Discard.unmod = Discard.unmod,
                             Min.N = 1,
                             Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1])
  #colnames(quant.data5)
  saveFun(quant.data5, file = "quant.data5.RData")
  #loadFun("quant.data5.RData")
}
if ((Param$QuantMeth == "Prot.Quant2")||(QuantMethods_all)) {
  # Averaging method, no weights
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using an in-house algorithm which: ",
                                "i) computes protein group expression values as the mean of peptidoform intensity values (\"relative quantitation\" step), then ",
                                "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                "only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
  quant.data6 <- Prot.Quant2(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep[Pep2Use,], id = "id",
                             Summary.method = "mean", #Summary.weights = "Weights",
                             experiments.map = Exp.map, param = Param,
                             Pep.Intens.root = pep.ref[length(pep.ref)],
                             Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                             log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                             Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                             Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                             Discard.unmod = Discard.unmod,
                             Min.N = 1,
                             Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1]
  )
  #colnames(quant.data6)
  saveFun(quant.data6, file = "quant.data6.RData")
  #loadFun("quant.data6.RData")
}
if ((Param$QuantMeth == "IQ_MaxLFQ")||(QuantMethods_all)) {
  # MaxLFQ
  if (!require(iq, quietly = TRUE)) { pak::pak("iq") }
  library(iq)
  kol <- grep(topattern(pep.ref[length(pep.ref)]), colnames(pep), value = TRUE)
  kol <- grep("\\.REF$", kol, value = TRUE, invert = TRUE)
  Pep2Use2 <- Pep2Use
  if (nrow(Mod2Xclud)) {
    pat <- paste(apply(Mod2Xclud, 1, function(x) {
      paste0(x[[2]], "(", x[[1]], ")", collapse = "|")
    }), collapse = "|")
    if (!Mod.Excl.is.strict) {
      Pep2Use2 <- Pep2Use2[grep(pat, pep$"Modified sequence"[Pep2Use2], invert = TRUE)]
    } else {
      seq <- pep$Sequence[Pep2Use2][grep(pat, pep$"Modified sequence"[Pep2Use2])]
      Pep2Use2 <- Pep2Use2[which(!pep$Sequence[Pep2Use2] %in% seq)]
    }
  }
  # That part rewritten to bypass the very slow iq::create_protein_list() 
  tmpPGlist <- listMelt(strsplit(PG$"Peptide IDs", ";"),
                        1:nrow(PG),
                        c("id", "Row"))
  tmpPGlist$id <- as.integer(tmpPGlist$id)
  tmpPGlist <- tmpPGlist[which(tmpPGlist$id %in% Pep$id[Pep2Use2]),]
  tmpPep2 <- tmpPep[, c("id", kol)]
  tmpPGlist[, kol] <- tmpPep2[match(tmpPGlist$id, tmpPep2$id), kol]
  colnames(tmpPGlist) <- gsub(topattern(pep.ref[length(pep.ref)]), "", colnames(tmpPGlist))
  kol2 <- gsub(topattern(pep.ref[length(pep.ref)]), "", kol)
  tmpPGlist <- as.data.table(tmpPGlist[, c("Row", kol2)])
  tmpPGlist <- split(tmpPGlist[, ..kol2], tmpPGlist$Row)
  names(tmpPGlist) <- Prot$`Leading protein IDs`[as.integer(names(tmpPGlist))]
  #
  # This part below is horrendously SLOWWWWWWWW...
  quant.data7 <- iq::create_protein_table(tmpPGlist)
  
  quant.data7 <- create_protTbl(tmpPGlist, cl = parClust)
  
  
  # One more reason to use in house quant...
  # I need to look into what that function is doing...
  # Also to put a warning in the parameters app if this becomes an option...
  #
  quant.data7 <- quant.data7$estimate/log2(10) # Convert to log10
  quant.data7 <- as.data.frame(quant.data7)
  colnames(quant.data7) <- paste0("log10(Expr.) - ", colnames(quant.data7))
  quant.data7$"Leading protein IDs" <- names(tmpPGlist)
  w <- which(!PG$"Leading protein IDs" %in% quant.data7$"Leading protein IDs")
  temp <- as.data.frame(matrix(rep(NA, length(w)*length(exprsCol)), ncol = length(exprsCol)))
  colnames(temp) <- exprsCol
  temp$"Leading protein IDs" <- PG$"Leading protein IDs"[w]
  quant.data7 <- rbind(quant.data7, temp)
  quant.data7 <- quant.data7[match(PG$"Leading protein IDs", quant.data7$"Leading protein IDs"),]
  #colnames(quant.data7)
  saveFun(quant.data7, file = "quant.data7.RData")
  #loadFun("quant.data7.RData")
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using the MaxLFQ algorithm as implemented in the iq package. ",
                                "Only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
}
if ((Param$QuantMeth == "Top3")||(QuantMethods_all)) {
  # Top3
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using a local implementation of the Top3 algorithm, allowing for a value to be calculated also if only 1 or 2 peptidoforms were present (the resulting systematic bias between proteins with 1, 2 or 3 peptidoforms is then corrected for); ",
                                "only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
  quant.data8 <- TopN(3, PG, "Razor peptide IDs", pep[Pep2Use,],
                      Pep.Intens.Nms = paste0(pep.ref[length(pep.ref)], RSA$values),
                      log.Pep.Intens = FALSE, Mods = Mod4Quant, Out.Norm = FALSE, corr = "global")
  colnames(quant.data8) <- paste0("log10(Expr.) - ", RSA$values)
  saveFun(quant.data8, file = "quant.data8.RData")
  #loadFun("quant.data8.RData")
}
if ((Param$QuantMeth == "Top1")||(QuantMethods_all)) {
  # Top1
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using the highest peptidoform intensity value in each sample. ",
                                "Only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
  quant.data9 <- TopN(1, Prot = PG, Pep = pep[Pep2Use,],
                      Pep.Intens.Nms = paste0(pep.ref[length(pep.ref)], RSA$values),
                      log.Pep.Intens = FALSE, Mods = Mod4Quant, Out.Norm = FALSE, corr = "global")
  colnames(quant.data9) <- paste0("log10(Expr.) - ", RSA$values)
  saveFun(quant.data9, file = "quant.data9.RData")
  #loadFun("quant.data9.RData")
}
if (QuantMethods_all) {
  # Comparison of LFQ profiles and summation approaches:
  cv <- CV <- CV_SD <- list()
  MaxPep <- 5
  # Create average intensity normalisation across methods
  # We want to look at the variability between methods,
  # but need to correct for differences in estimates of absolute abundance,
  # as we are really interested in precision of sample-to-sample variability
  NormFact1 <- sapply(seq_along(QuantMethods), function(i) {
    tmp <- get(QuantData[i])
    tmp$"Leading protein IDs" <- PG$"Leading protein IDs"
    tmp2 <- tmp[, exprsCol]
    tmp2 <- rowMeans(tmp2, na.rm = TRUE)
    w <- which(PG$"Leading protein IDs" %in% tmp$"Leading protein IDs")
    res <- rep(NA, nrow(PG))
    res[w] <- tmp2[match(PG$"Leading protein IDs"[w], tmp$`Leading protein IDs`)]
    return(res)
  })
  NormFact2 <- rowMeans(NormFact1, na.rm = TRUE)
  for (i in seq_along(QuantMethods)) {
    tmp1 <- get(QuantData[i])
    tmp1$"Leading protein IDs" <- PG$"Leading protein IDs"
    tmp1 <- tmp1[, c(exprsCol, "Leading protein IDs")]
    tmp2 <- sapply(VPAL$values, function(x) {
      kol <- paste0("log10(Expr.) - ", Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)])
      x <- apply(tmp1[ , kol, drop = FALSE], 1, function(y) {
        y <- is.all.good(y)
        sd(y)/mean(y)
      })
    })
    tmp2 <- data.frame(CV = rowMeans(tmp2, na.rm = TRUE))
    colnames(tmp1) <- gsub(topattern("log10(Expr.) - "), "", colnames(tmp1))
    tmp1[, c("id", "Razor + unique peptides")] <- PG[match(tmp1$"Leading protein IDs", PG$"Leading protein IDs"),
                                                     c("id", "Peptide counts (razor+unique)")]
    tmp1$"Razor + unique peptides" <- as.integer(gsub(";.*", "", tmp1$"Razor + unique peptides"))
    tmp2$"Razor + unique peptides" <- tmp1$"Razor + unique peptides"
    tmp2$`Razor + unique peptides`[which(tmp2$`Razor + unique peptides` >= MaxPep)] <- paste0(MaxPep, " or more") # No need to look at too rare cases
    tmp2$`Razor + unique peptides` <- as.factor(tmp2$`Razor + unique peptides`)
    cv[[QuantMethods[i]]] <- tmp2[which(is.all.good(tmp2$CV, 2)),]
    CV[[QuantMethods[i]]] <- mean(cv[[QuantMethods[i]]]$CV)
    tmp3 <- tmp1
    m <- match(tmp3$`Leading protein IDs`, PG$`Leading protein IDs`)
    tmp3[, gsub(topattern("log10(Expr.) - "), "", exprsCol)] <- sweep(tmp3[, gsub(topattern("log10(Expr.) - "), "", exprsCol)], 1, NormFact1[m,i], "-")+NormFact2[m]
    tmp1 <- reshape2::melt(tmp1, id.vars = c("id", "Razor + unique peptides", "Leading protein IDs"))
    tmp3 <- reshape2::melt(tmp3, id.vars = c("id", "Razor + unique peptides", "Leading protein IDs"))
    colnames(tmp1) <- c("id", "Razor + unique peptides", "Leading protein IDs", "Sample", QuantMethods[i])
    tmp1[[paste0(QuantMethods[i], " - renorm")]] <- tmp3$value
    tmp1$Sample <- as.character(tmp1$Sample)
    tmp1[, RSA$names] <- Isapply(strsplit(tmp1$Sample, "___"), unlist)
    tmp1$Sample <- cleanNms(tmp1$Sample)
    if (i == 1) { tmp <- tmp1 } else {
      tmp[[QuantMethods[i]]] <- NA
      w <- which(tmp$"Leading protein IDs" %in% tmp1$"Leading protein IDs")
      tmp[w, QuantMethods[i]] <- tmp1[match(tmp$"Leading protein IDs"[w], tmp1$"Leading protein IDs"), QuantMethods[i]] 
    }
  }
  #plot <- ggplot(tmp) + geom_scattermore(aes(x = Prot.Quant, y = Prot.Quant2, colour = Sample), size = 0.1, alpha = 0.01) +
  #  coord_fixed() + theme_bw() + facet_wrap(~Sample)
  comb <- gtools::combinations(length(QuantMethods), 2, QuantMethods)
  dir <- paste0(wd, "/Workflow control/Protein groups/Quantitative methods")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  for (i in 1:nrow(comb)) { #i <- 1
    tmp1 <- tmp
    tmp1$X <- tmp1[[comb[i, 1]]]
    tmp1$Y <- tmp1[[comb[i, 2]]]
    tmp1 <- tmp1[which((is.all.good(tmp1$X, 2))&(is.all.good(tmp1$Y, 2))),]
    tmp1$`Razor + unique peptides`[which(tmp1$`Razor + unique peptides` >= MaxPep)] <- paste0(MaxPep, " or more") # No need to look at too rare cases
    tmp1$`Razor + unique peptides` <- as.factor(tmp1$`Razor + unique peptides`)
    tmp2 <- tmp1
    tmp2$X <- tmp2[[paste0(comb[i, 1], " - renorm")]]
    tmp2$Y <- tmp2[[paste0(comb[i, 2], " - renorm")]]
    tst <- vapply(RSA$names, function(x) { length(get(substr(x, 1, 3))) }, 1)
    tst <- RSA$names[which(tst > 1)]
    form <- as.formula(paste0(tst[1], "~", paste(tst[2:length(tst)], collapse = "+")))
    ttl1 <- paste0("LFQ method comparison: ", comb[i, 1], " VS ", comb[i, 2])
    plot1 <- ggplot(tmp1) +
      geom_scattermore(aes(x = X, y = Y, colour = `Razor + unique peptides`),
                       shape = 16, size = 0.1, alpha = 0.1) +
      scale_color_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() + facet_grid(form) + xlab(comb[i, 1]) + ylab(comb[i, 2]) + ggtitle(ttl1)
    #poplot(plot1, 12, 20)
    ttl1a <- gsub("\\:", " -", ttl1)
    suppressMessages({
      ggsave(paste0(dir, "/", ttl1a, ".jpeg"), plot1, dpi = 300)
      ggsave(paste0(dir, "/", ttl1a, ".pdf"), plot1, dpi = 300)
    })
    ReportCalls <- AddPlot2Report(plot1, Title = ttl1a)
    ttl2 <- paste0("LFQ method comparison (renormalized): ", comb[i, 1], " VS ", comb[i, 2])
    plot2 <- ggplot(tmp1) +
      geom_scattermore(aes(x = X, y = Y, colour = `Razor + unique peptides`),
                       shape = 16, size = 0.1, alpha = 0.1) +
      scale_color_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() + facet_grid(form) + xlab(comb[i, 1]) + ylab(comb[i, 2]) + ggtitle(ttl2)
    print(plot2) # This type of QC plot does not need to pop up, the side panel is fine
    ttl2a <- gsub("\\:", " -", ttl1)
    suppressMessages({
      ggsave(paste0(dir, "/", ttl2a, ".jpeg"), plot2, dpi = 300)
      ggsave(paste0(dir, "/", ttl2a, ".pdf"), plot2, dpi = 300)
    })
    ReportCalls <- AddPlot2Report(plot2, Title = ttl2a)
  }
  for (meth in QuantMethods) {
    cat(paste0(meth, ": ", nrow(cv[[meth]]), " quantified protein groups, mean intra-sample-groups CV = ", round(CV[[meth]]*100, 2), "%\n"))
  }
  tst <- reshape2::melt(cv)
  tst$variable <- NULL
  colnames(tst) <- c("Razor + unique peptides", "CV", "Method")
  ttl <- "Distribution of intra-sample groups coefficients of variation"
  plot <- ggplot(tst) +
    geom_histogram(aes(x = CV, fill = Method), bins = 100) + theme_bw() + ggtitle(ttl) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_d(begin = 0.25) +
    facet_grid(Method~`Razor + unique peptides`) +
    theme(strip.text.y.right = element_text(angle = 0))
  print(plot) # This type of QC plot does not need to pop up, the side panel is fine
  ttla <- gsub(":", " -", ttl)
  suppressMessages({
    ggsave(paste0(dir, "/", ttla, ".jpeg"), plot, dpi = 300)
    ggsave(paste0(dir, "/", ttla, ".pdf"), plot, dpi = 300)
  })
  ReportCalls <- AddPlot2Report(Title = ttla)
}
ReportCalls <- AddTxt2Report(paste0("Protein groups quantitation done using method: ", Param$QuantMeth))
if (!exists(QuantData[match(Param$QuantMeth, QuantMethods)])) { # (for when rerunning script without rerunning quant methods)
  load(paste0(QuantData[match(Param$QuantMeth, QuantMethods)], ".RData"))
}
quant.data <- get(QuantData[match(Param$QuantMeth, QuantMethods)])
Prot.Expr.Root %<o% c(Original = "log10(Expr.) - ")
Prot.Rat.Root %<o% pep.ratios.ref[length(pep.ratios.ref)]
#write.csv(quant.data, file = "Quantitative data.csv", row.names = FALSE)
.obj <- unique(c("quant.data", "Prot.Expr.Root", "Prot.Rat.Root", .obj))
DatAnalysisTxt <- paste0(DatAnalysisTxt, " Estimated expression values were log10-converted...")

### Code chunk
# Sometimes, you want to compare, e.g., two control groups C1 and C2 to one treated sample S.
# In this case, because the workflow expects one reference per group, we would calculate ratios as C1/S and C2/S.
# But since S is the perturbation, people still like to see their ratios the other way, i.e. S/C1 and S/C2.
# This can be set with the "Mirror.Ratios" parameter. The reversion is applied here.
# This is not used currently, check validity before applying!!!
if ("Mirror.Ratios" %in% colnames(Param)) { Mirror.Ratios <- Param$Mirror.Ratios <- as.logical(Param$Mirror.Ratios) }
if (!is.logical(Mirror.Ratios)) {
  warning("I could not make sense of the value of parameter \"Mirror.Ratios\", defaulting to FALSE!")
  Mirror.Ratios <- FALSE
}
if (Mirror.Ratios) {
  stop("The code has changed, this is now deprecated: if you want to re-use this then first check all normalisation steps where ratios are calculated after this stage!!!")
  if (!exists("InvertedOnce")) {
    Prot.Rat.Root.BckUp <- Prot.Rat.Root
    pep.ratios.ref.BckUp <- pep.ratios.ref    
    .obj <- unique(c("Prot.Rat.Root.BckUp", "pep.ratios.ref.BckUp", .obj))
  }
  # Protein groups
  G <- grep(topattern(Prot.Rat.Root.BckUp), colnames(quant.data), value = TRUE)
  for (g in G) {
    quant.data[[g]] <- -quant.data[[g]]# Because this data is already log-transformed
    colnames(quant.data)[which(colnames(quant.data) == g)] <- paste0("-", g)
  }
  if (grepl("-log2", Prot.Rat.Root.BckUp)) { Prot.Rat.Root <- gsub("-log2", "log2", Prot.Rat.Root.BckUp) } else {
    Prot.Rat.Root <- gsub("log2", "-log2", Prot.Rat.Root.BckUp)
  }
  # Peptides
  G <- grep("^Ratio - ", colnames(pep), value = TRUE)
  for (g in G) {
    pep[[g]] <- 1/pep[[g]]# Because this data is not log-transformed...
    colnames(pep)[which(colnames(pep) == g)] <- paste0("1/", g)
  }
  for (prf in pep.ratios.ref.BckUp) { #prf <- pep.ratios.ref.BckUp[1]
    G <- grep(topattern(prf), colnames(pep), value = TRUE)
    for (g in G) {
      pep[[g]] <- -pep[[g]] # ... but this data is!
      colnames(pep)[which(colnames(pep) == g)] <- paste0("-", g)
    }
  }
  if (grepl("-log2", pep.ratios.ref.BckUp)) { pep.ratios.ref <- gsub("-log2", "log2", pep.ratios.ref.BckUp) } else {
    pep.ratios.ref <- gsub("log2", "-log2", pep.ratios.ref.BckUp)
  }
  # No need to do Ref.Ratios, they get created later...
  # ... but pep.Ref.Ratios already exist (and are log-transformed)
  if (Param$Ratios.Thresholds == threshMsg) {
    for (k in colnames(pep.Ref.Ratios)) { pep.Ref.Ratios[[k]] <- -pep.Ref.Ratios[[k]] }
    colnames(pep.Ref.Ratios) <- gsub("log2", "-log2", colnames(pep.Ref.Ratios))
    # To test that the evidence file does not contain ratios:
    stopifnot(length(grep("[Rr]atio[ -\\._$]", colnames(ev), value = TRUE)) == 0)
    # I do not expect this, but this could happen under two scenarii:
    # - I change this workflow significantly
    # - I want to adapt this workflow to SILAC, where MaxQuant directly - and precisely - measures ratios
    InvertedOnce %<o% TRUE
  }
}

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Re-normalize protein group expression values
Src <- paste0(libPath, "/extdata/R scripts/Sources/PG_ReNorm.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Test expression values:
g <- grep(topattern(Prot.Expr.Root), colnames(quant.data), value = TRUE)
g <- grep(": SD$", g, value = TRUE, invert = TRUE)
g <- grep("\\.REF$", g, value = TRUE, invert = TRUE)
test <- quant.data[,g]
colnames(test) <- gsub(topattern(Prot.Expr.Root), "", colnames(test))
test <- test[which(apply(test, 1, function(x) { length(is.all.good(x)) }) > 0),]
test <- suppressMessages(melt.data.frame(test))
test$variable <- as.character(test$variable)
test[, RSA$names] <- ""
#w <- grepl("\\.REF$", test$variable)
w <- rep(FALSE, nrow(test))
test[which(!w), RSA$names] <- Isapply(strsplit(test$variable[which(!w)], "___"), unlist)
#test[which(w), RRG$names] <- a <- Isapply(strsplit(gsub_Rep("\\.REF$", "", test$variable[which(w)]), "___"), unlist)
#test$REF <- c("Individual", "Reference")[1+grepl("\\.REF$", test$variable)]
a <- RSA$names
w <- which(vapply(a, function(x) { length(unique(test[[x]])) }, 1) > 1)
if (length(w)) { a <- a[w] }
test[[a[1]]] <- factor(test[[a[1]]], levels = sort(unique(test[[a[1]]])))
test <- test[which(is.all.good(test$value, 2)),]
test2 <- set_colnames(aggregate(test$value, list(test$variable), median), c("variable", "value"))
test2[, a] <- test[match(test2$variable, test$variable), a]
MinMax <- c(min(test$value), max(test$value))
nbinz <- ceiling((MinMax[2]-MinMax[1])/0.1)
binz <- c(0:nbinz)/nbinz
binz <- binz*(MinMax[2]-MinMax[1])+MinMax[1]
binz[1] <- binz[1]-0.000001
testI <- data.frame(Intensity = (binz[2:(nbinz+1)]+binz[1:nbinz])/2)
for (v in unique(test$variable)) {
  wv <- which(test$variable == v)
  testI[[v]] <- vapply(1:nbinz, function(x) {
    sum((test$value[wv] > binz[x])&(test$value[wv] <= binz[x+1]))
  }, 1)
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
  scale_y_continuous(limits = c(0, max(testI$value)*1.1), expand = c(0, 0))
if (length(a) == 1) { plot <- plot + facet_wrap(as.formula(paste0("~", a))) } else {
  plot <- plot + facet_grid(as.formula(paste0(a[1], "~", paste(a[2:length(a)], collapse = "+"))))
}
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")
})
ReportCalls <- AddPlot2Report()
#
# Test ratio values:
g <- grep(topattern(Prot.Rat.Root), colnames(quant.data), value = TRUE)
g <- grep(": SD$", g, value = TRUE, invert = TRUE)
g <- grep("REF\\.to\\.REF", g, value = TRUE, invert = TRUE)
test <- quant.data[,g]
colnames(test) <- gsub(topattern(Prot.Rat.Root), "", colnames(test))
test <- test[which(apply(test, 1, function(x) { length(is.all.good(x)) }) > 0),]
test <- suppressMessages(melt.data.frame(test))
test$variable <- as.character(test$variable)
test[, RSA$names] <- ""
#w <- grepl("\\.REF$", test$variable)
w <- rep(FALSE, nrow(test))
test[which(!w), RSA$names] <- Isapply(strsplit(test$variable[which(!w)], "___"), unlist)
#test[which(w), RRG$names] <- a <- Isapply(strsplit(gsub_Rep("\\.REF$", "", test$variable[which(w)]), "___"), unlist)
#test$REF <- c("Individual", "Reference")[1+grepl("\\.REF$", test$variable)]
a <- RSA$names
w <- which(vapply(a, function(x) { length(unique(test[[x]])) }, 1) > 1)
if (length(w)) { a <- a[w] }
test[[a[1]]] <- factor(test[[a[1]]], levels = sort(unique(test[[a[1]]])))
test <- test[which(is.all.good(test$value, 2)),]
test2 <- set_colnames(aggregate(test$value, list(test$variable), median), c("variable", "value"))
test2[, a] <- test[match(test2$variable, test$variable), a]
MinMax <- c(min(test$value), max(test$value))
nbinz <- ceiling((MinMax[2]-MinMax[1])/0.1)
binz <- c(0:nbinz)/nbinz
binz <- binz*(MinMax[2]-MinMax[1])+MinMax[1]
binz[1] <- binz[1]-0.000001
testR <- data.frame(Intensity = (binz[2:(nbinz+1)]+binz[1:nbinz])/2)
for (v in unique(test$variable)) {
  wv <- which(test$variable == v)
  testR[[v]] <- vapply(1:nbinz, function(x) {
    sum((test$value[wv] > binz[x])&(test$value[wv] <= binz[x+1]))
  }, 1)
}
testR <- reshape2::melt(testR, id.vars = "Intensity")
testR[, a] <- test[match(testR$variable, test$variable), a]
testR$variable <- cleanNms(testR$variable)
dir <- paste0(wd, "/Workflow control/Protein groups/Expression")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
ttl <- "Protein groups - distribution of Ratios"
plot <- ggplot(testR) +
  geom_area(aes(x = Intensity, y = value, fill = variable, group = variable,
                colour = variable), alpha = 0.25) +
  scale_color_viridis(begin = 0.25, discrete = TRUE, option = "D") +
  scale_fill_viridis(begin = 0.25, discrete = TRUE, option = "D") +
  geom_vline(data = test2, aes(xintercept = value), linetype = "dashed", color = "grey") +
  ggtitle(ttl) + theme_bw() + theme(legend.position = "none", strip.text.y = element_text(angle = 0)) +
  scale_y_continuous(limits = c(0, max(testR$value)*1.1), expand = c(0, 0))
if (length(a) == 1) { plot <- plot + facet_wrap(as.formula(paste0("~", a))) } else {
  plot <- plot + facet_grid(as.formula(paste0(a[1], "~", paste(a[2:length(a)], collapse = "+"))))
}
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")
})
ReportCalls <- AddPlot2Report()

# For Sub-cellular localisation analysis
# (or any other analysis where we want to apply some additional a priori known normalisation factor):
# Re-scale to old sub-cellular fraction relative levels, as measured:
#  - Sample specific evidence intensity range prior to evidence normalisation,
#  - The "Proportion" columns in Exp.map, which is the ratio of loaded to total available sample amount.
# Because we would still want to have normalized data, if either of the following arguments is TRUE,
# the re-scaling is done within the VPAL groups:
# - "Norma.Ev.Intens"
# - "Norma.Pep.Intens"
# - "Norma.Pep.Ratio"
# - "Norma.Prot.Intens"
# - "Norma.Prot.Ratio"
if (LocAnalysis) {
  if (!"Proportion" %in% colnames(Exp.map)) {
    warning("The \"Proportion\" column is absent from the Experiment map, we will assume that the same proportion of each fraction was processed and analyzed.")
    Exp.map$Proportion <- 1 # In case the column was omitted
  } else { Exp.map$Proportion <- as.numeric(Exp.map$Proportion) }
  #
  Prot.Expr.Root2 %<o% setNames(gsub("Expr\\.", "resc. Expr.", Prot.Expr.Root), "SubCell. Profile")
  pep.ref2 %<o% setNames(paste0("resc. ", pep.ref[1]), "SubCell. Profile")
  # Best way to do this:
  # - Compare final values at peptidoforms level
  # - Get original, uncorrected evidences values and aggregate them into peptidoforms level values
  # - Apply the corresponded ratios to protein groups and peptidoforms-level data to recreate original scale
  #
  # Also, take into account the proportion of each fraction which was actually processed
  prevRef <- pep.ref[length(pep.ref)]
  kol1 <- paste0(prevRef, RSA$values)
  WhInColNms <- which(kol1 %in% colnames(pep))
  kol1 <- kol1[WhInColNms]
  tmp1 <- pep[, c("Modified sequence", kol1)]
  kol2 <- kol2a <- c("Modified sequence", "Raw file path", "MQ.Exp", "Experiment", ev.col["Original"])
  if (LabelType == "Isobaric") {
    kol2b <- grep(paste0(topattern(ev.ref["Original"]), "[0-9]+$"), colnames(ev), value = TRUE)
    chan <- gsub(topattern(ev.ref["Original"]), "", kol2b)
    kol2 <- c(kol2a, kol2b)
    tmp2 <- ev[, kol2]
    tst <- apply(tmp2[, kol2b], 1, sum)
    tmp2[, kol2b] <- sweep(tmp2[, kol2b], 1, tmp2[[ev.col["Original"]]]/tst, "*")
    tmp2a <- aggregate(tmp2[, kol2b], list(tmp2$"Modified sequence", tmp2$MQ.Exp), sum)
    tmp2 <- data.frame(`Modified sequence` = pep$"Modified sequence", check.names = FALSE)
    for (mqexp in MQ.Exp) { #mqexp <- MQ.Exp[1]
      em <- Exp.map[which(Exp.map$MQ.Exp == mqexp),]
      m <- match(chan, em$"Isobaric label")
      w <- which(!is.na(m))
      kol2c <- paste0(pep.ref["Original"], em$Ref.Sample.Aggregate[m[w]])
      tmp2[, kol2c] <- 0
      wa <- which(tmp2a$Group.2 == mqexp)
      wb <- which(tmp2$"Modified sequence" %in% tmp2a$Group.1[wa])
      tmp2[wb, kol2c] <- tmp2a[wa[match(tmp2$"Modified sequence"[wb], tmp2a$Group.1[wa])], kol2b[w]]
    }
  } else {
    tmp2a <- as.data.table(ev[, kol2])
    colnames(tmp2a)[which(colnames(tmp2a) == ev.col["Original"])] <- "Int"
    tmp2a <- tmp2a[, list(x = sum(Int)),
                   by = list(`Modified sequence` = `Modified sequence`,
                             `Raw file` = `Raw file path`,
                             MQ.Exp = MQ.Exp,
                             Experiment = Experiment)]
    tmp3 <- listMelt(Exp.map$MQ.Exp, Exp.map$Ref.Sample.Aggregate)
    tmp2a$RSA <- tmp3$L1[match(tmp2a$MQ.Exp, tmp3$value)]
    tmp2 <- data.frame(`Modified sequence` = pep$"Modified sequence", check.names = FALSE)
    for (rsa in tmp3$L1) { #rsa <- tmp3$L1[1]
      kol2c <- paste0(pep.ref["Original"], rsa)
      tmp2[[kol2c]] <- 0
      w <- which(tmp2a$RSA == rsa)
      tmp2b <- copy(tmp2a)
      tmp2b <- tmp2b[w, list(x = sum(x)), by = list(`Modified sequence`)]
      tmp2b <- as.data.frame(tmp2b)
      w <- which(tmp2$"Modified sequence" %in% tmp2b$`Modified sequence`)
      tmp2[w, kol2c] <- tmp2b$x[match(tmp2$"Modified sequence"[w], tmp2b$`Modified sequence`)]
    }
  }
  # re-order
  kol2c <- paste0(pep.ref["Original"], RSA$values)
  kol2c <- kol2c[WhInColNms]
  tmp2 <- tmp2[, c("Modified sequence", kol2c)]
  # Calculate ratio before/after, per samples group
  BefAft <- tmp2[, kol2c]/tmp1[, kol1]
  colnames(BefAft) <- RSA$values[WhInColNms]
  BefAft <- vapply(VPAL$values, function(x) {
    x <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)]
    x <- x[which(x %in% RSA$values[WhInColNms])]
    x <- median(unlist(BefAft[, x]), na.rm = TRUE)
    return(x)
  }, 1)
  # This is the ratio between values at the start of the workflow and final peptidoforms values 
  #
  # Next, we can apply proportion of total material loaded
  # Now...
  # At this stage we should have nicely normalized, i.e., "aligned", data
  # We should thus work the same way within groups of replicates.
  # Thus, if there are different loaded amounts per individual sample,
  # when we apply corrections aimed at restoring original samples' relative intensity scales,
  # we should average at this stage over replicates of the same condition.
  Props <- Exp.map[match(RSA$values[WhInColNms], Exp.map$Ref.Sample.Aggregate),
                   c("Proportion", VPAL$column)]
  Props <- aggregate(Props$Proportion, list(Props[[VPAL$column]]), mean)
  Props <- Props$x[match(VPAL$values, Props$Group.1)]
  BefAft <- BefAft/Props
  # Note: this is done currently at sample group level. It may make sense to do it at subcellular fraction level in the future...
  # but there are also risks, e.g. when the perturbation studied changes cell morphology dramatically.
  # => make it an option?
  #
  # Apply correction:
  # We want to calculate the average change within sample groups (compartments/fractions x treatment group)
  # so as to still profit from the normalisation did on peptides/proteins.
  # (Presumably it is ok to normalize within those groups, as samples should be replicates)
  # For now the covariates aggregate used is VPAL, but we could map it to a custom one using a new Param
  for (grp in VPAL$values) { #grp <- VPAL$values[1]
    smpls <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == grp)]
    PepKol1 <- paste0(prevRef, smpls)
    PepKol2 <- paste0(pep.ref2, smpls) 
    w <- which(PepKol1 %in% colnames(pep))
    pep[, PepKol2[w]] <- pep[, PepKol1[w]]*BefAft[grp]
    PGKol1 <- paste0(Prot.Expr.Root, smpls)
    PGKol2 <- paste0(Prot.Expr.Root2, smpls) 
    w <- which(PGKol1 %in% colnames(quant.data))
    quant.data[, PGKol2[w]] <- quant.data[, PGKol1[w]] + log10(BefAft[grp])
    # NB: Ratios are not re-scaled...
    # and anyway, the ratios should be done within groups which should correspond to the correct proportions 
  }
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " Expression values were re-scaled per fraction to reflect total original protein amount.")
  # Visualize
  dir <- paste0(wd, "/Workflow control/Re-scaling/")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  for (i in 1:2) {
    if (i == 1) {
      rt1 <- prevRef
      rt2 <- pep.ref2
      ttl <- "Peptide intensities - effect of post-normalisation rescaling"
      temp <- pep
    }
    if (i == 2) {
      rt1 <- Prot.Expr.Root
      rt2 <- Prot.Expr.Root2
      ttl <- "PGs expression values - effect of post-normalisation rescaling"
      temp <- quant.data
    }
    kol1 <- paste0(rt1, RSA$values)
    kol2 <- paste0(rt2, RSA$values)
    w <- which((kol1 %in% colnames(temp))&(kol2 %in% colnames(temp)))
    tst <- temp[, c(kol1[w], kol2[w])]
    tst <- reshape2::melt(tst, measure.vars = c(kol1[w], kol2[w]))
    tst$value <- suppressWarnings(log10(tst$value))
    tst <- tst[which(is.all.good(tst$value, 2)),]
    tst$variable <- as.character(tst$variable)
    tst2 <- data.frame(variable = unique(tst$variable))
    tst2[, c("Type", "Sample")] <- as.data.frame(t(sapply(strsplit(tst2$variable, " - "), unlist)))
    tst2[, RSA$names] <- as.data.frame(t(sapply(strsplit(tst2$Sample, "___"), unlist)))
    tst2$Group <- Exp.map[match(tst2$Sample, Exp.map$Ref.Sample.Aggregate), Volcano.plots.Aggregate.Level$aggregate]
    tst2$Sample <- cleanNms(tst2$Sample)
    tst2$Group <- cleanNms(tst2$Group)
    tst2$Type <- factor(tst2$Type, levels = gsub(" - $", "", c(rt1, rt2)))
    tst[, colnames(tst2)] <- tst2[match(tst$variable, tst2$variable), colnames(tst2)]
    plot <- ggplot(tst) +
      geom_violin(aes(x = Sample, y = value, colour = Group, fill = Group), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Group, fill = Group), alpha = 0.5) +
      scale_color_viridis(begin = 0.25, discrete = TRUE, option = "C") +
      scale_fill_viridis(begin = 0.25, discrete = TRUE, option = "C") +
      facet_grid(Type~.) + theme_bw() + ggtitle(ttl) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    if (i == 2) {
      print(plot) # This type of QC plot does not need to pop up, the side panel is fine
    }
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    })
    ReportCalls <- AddPlot2Report()
  }
}

#### Code chunk - Optional - Apply True-Discovery or negative filter
# Let's assume that our data is contaminated (e.g. impure fractions),
# but that we have another source to assess the validity of protein groups discoveries.
# Here we can load a table of known valid protein groups, with one column for each value of RG
# Any protein group not matching that table, and its peptides, will have all quantitative values set to NA from here on.
if (DiscFilt) {
  # Create filter
  DiscFiltFilt %<o% strsplit(PG$"Leading protein IDs", ";")
  DiscFiltFilt <- listMelt(DiscFiltFilt, PG$id)
  colnames(DiscFiltFilt) <- c("Leading protein ID", "PG ID")
  w <- which(DiscFiltFilt$"Leading protein ID" %in% DiscFiltTbl$"Protein ID")
  if (DiscFiltMode %in% DiscFiltModes[1:2]) {
    # Apply TRUE/FALSE from loaded filter
    for (grp in RG$values) { #grp <- RG$values[1]
      DiscFiltFilt[[grp]] <- c(FALSE, TRUE)[match(DiscFiltMode, DiscFiltModes)]
      tmp <- DiscFiltTbl[match(DiscFiltFilt$"Leading protein ID"[w], DiscFiltTbl$`Protein ID`), grp]
      if (DiscFiltMode == DiscFiltModes[2]) { tmp <- !tmp }
      DiscFiltFilt[w, grp] <- tmp
    }
    if (DiscFiltMode == DiscFiltModes[1]) {
      # We only remove a PG if no leading protein is TRUE in the filter
      DiscFiltFilt <- aggregate(DiscFiltFilt[, RG$values], list(DiscFiltFilt$"PG ID"), function(x) { as.logical(max(x)) })
    }
    if (DiscFiltMode == DiscFiltModes[2]) {
      # We remove a PG if any leading protein is FALSE in the filter
      DiscFiltFilt <- aggregate(DiscFiltFilt[, RG$values], list(DiscFiltFilt$"PG ID"), function(x) { as.logical(min(x)) })
    }
    colnames(DiscFiltFilt) <- c("PG ID", RG$values)
    DiscFiltFilt <- DiscFiltFilt[match(PG$id, DiscFiltFilt$"PG ID"),] # Re-order
    # Apply filter to quantitative data
    for (grp in RG$values) { #grp <- RG$values[1]
      w <- which(!DiscFiltFilt[[grp]])
      em <- Exp.map[which(Exp.map[[RG$column]] == grp),]
      kol <- c(paste0(Prot.Expr.Root, em$Ref.Sample.Aggregate),
               grep(topattern(paste0(Prot.Expr.Root, grp, ".REF")), colnames(quant.data), value = TRUE),
               paste0(Prot.Rat.Root, em$Ref.Sample.Aggregate),
               grep(topattern(paste0(Prot.Rat.Root, grp, "_REF.to.REF_")), colnames(quant.data), value = TRUE))
      kol <- kol[which(kol %in% colnames(quant.data))]
      quant.data[w, kol] <- NA
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " Data was filtered to only include proteins identified in the provided true-discovery filter.")
    ReportCalls <- AddTxt2Report("Removing proteins not included in the TRUE-Discovery filter!")
  }
  if (DiscFiltMode == DiscFiltModes[3]) {
    if (length(unique(RG$values)) == 1) {
      DiscFiltCols <- DiscFiltCol
    } else { DiscFiltCols <- paste0(DiscFiltCol, " - ", RG$values) }
    DiscFiltFilt[, RG$values] <- ""
    for (grp in RG$values) { DiscFiltTbl[[grp]] <- c("", "+")[DiscFiltTbl[[grp]]+1] }
    DiscFiltFilt[w, RG$values] <- DiscFiltTbl[match(DiscFiltFilt$`Leading protein ID`[w], DiscFiltTbl$`Protein ID`), RG$values]
    DiscFiltFilt <- aggregate(DiscFiltFilt[, RG$values], list(DiscFiltFilt$`PG ID`), function(x) {
      c("", "+")[("+" %in% unlist(x))+1]
    })
    colnames(DiscFiltFilt) <- c("id", DiscFiltCols)
    PG[, DiscFiltCols] <- ""
    w <- which(PG$id %in% DiscFiltFilt$id)
    PG[w, DiscFiltCols] <- DiscFiltFilt[match(PG$id[w], DiscFiltFilt$id), DiscFiltCols]
  }
}

# Code chunk - Add quant data to PG:
PG <- PG[, which(!colnames(PG) %in% colnames(quant.data))]
PG[, colnames(quant.data)] <- quant.data
if (Param$Prot.Only.with.Quant) {
  colnames(quant.data)
  test1 <- apply(quant.data[,grep(topattern(Prot.Expr.Root), colnames(quant.data), value = TRUE)],
                 1, function(x) {length(is.all.good(x))})
  a <- grep(topattern(Prot.Rat.Root), colnames(quant.data), value = TRUE)
  a <- a[which(!grepl(": SD$|: -log10\\(peptides Pvalue\\)$", a))]
  test2 <- apply(quant.data[,a],
                 1, function(x) {length(is.all.good(x))})
  PG <- PG[which((test1 > 0)|(test2 > 0)),]
}
if (!"Peptides count" %in% colnames(PG)) {
  PG$"Peptides count" <- vapply(strsplit(PG$"Peptide IDs", ";"), length, 1)
}
#if ("Prot.Only.with.at.least" %in% colnames(Param)) {
#  n <- as.numeric(Param$Prot.Only.with.at.least)
#  if ((is.na(n))||(n != round(n, 0))) {
#    warning("Protein groups could not be filtered by number of peptides because the value entered could not be recognized as a number.")
#  } else {
#    w <- which(PG$"Peptides count" >= n)
#    quant.data <- quant.data[w,]
#    PG <- PG[w,]
#  }
#}
#pep_bckp %<o% pep
#pep <- pep[which(pep$id %in% unique(unlist(strsplit(PG$"Peptide IDs", ";")))),]

if (!Param$Plot.labels %in% colnames(PG)) {
  tmp <- gsub("\\.", " ", Param$Plot.labels)
  if (tmp %in% colnames(PG)) {
    warning(paste0("Protein groups table column \"", Param$Plot.labels, "\" not found, the column is called \"", tmp, "\" (check parameter \"Plot.labels\")"))
    Param$Plot.labels <- tmp
  } else {
    tmp <- c("Common Name (short)", "Common Names", "Names", "Protein IDs", "Common.Names.short", "Common.Names", "Protein.IDs")
    w <- which(tmp %in% colnames(PG))
    warning(paste0("Protein groups table column \"", Param$Plot.labels, "\" not found (check parameter \"Plot.labels\"), defaulting to \"", tmp[w[1]], "\""))
    Param$Plot.labels <- tmp[w[1]]
  }
}

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)

#### Code chunk - samples Pearson correlation heatmap
Src <- paste0(libPath, "/extdata/R scripts/Sources/pearsonCorrMap.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Calculate average intensities and ratios, and perform a few statistical tests
# At least, Welch's t-test and moderated t-test; for unpaired replicates a permutations t-test is also performed.
# By default in the "t.test" function var.equal is set to FALSE, which performs a Welch's t.test.
# This is as well, since Welch's t.test apparently never performs worse than Student's.
Av_SE_fun %<o% function(vect) {
  res <- proteoCraft::is.all.good(as.numeric(vect))
  if (length(res)) { res <- c(mean(res), sd(res)/sqrt(length(res))) } else {
    res <- unique(vect[which((!is.nan(vect))&(!is.na(vect)))])
    if (length(res) == 1) { return(c(res, NA)) } else { return(c(NA, NA)) }
    return(res)
  }
}
# 
# Internal function nicked from https://github.com/cran/miRtest/blob/master/R/miRtest.R
# Barely modified
limma.one.sided %<o% function(myFit, lower) {
  se.coef <- sqrt(myFit$s2.post) * myFit$stdev.unscaled
  df.total <- myFit$df.prior + myFit$df.residual
  rs <- pt(myFit$t, df = df.total, lower.tail = lower)
  return(rs[, colnames(myFit$p.value), drop = FALSE])
}
ReportCalls <- AddTxt2Report("Calculating average intensities and ratios and performing statistical tests...")
# Alternative hypothesis for tests:
#AltHyp <- c("two.sided", "greater")[IsPullDown+1]
AltHyp %<o% c(c("greater", "lower")[Mirror.Ratios+1], "two.sided")[TwoSided+1]
# Why use two-sided in most cases?
# For pull-downs, we can still learn something on the left, e.g.:
# - what is enriched specifically in ctrl?
# - see how much variability we see on the left to compare with right
# In addition, it is difficult to implement for the other tests.

# Define design matrix and contrasts (limma) 
#  From Exp.map to design matrix
#Coefficients %<o% Factors[which(!Factors %in% c("Experiment", "Replicate"))]
Coefficients %<o% VPAL$names
w <- which(vapply(Coefficients, function(x) { length(unique(FactorsLevels[[x]])) < nrow(Exp.map) }, TRUE))
Coefficients <- Coefficients[w]
w <- which(vapply(Coefficients, function(x) { length(unique(Exp.map[[x]])) > 1 }, TRUE))
Coefficients <- Coefficients[w]
expMap %<o% Exp.map
expMap <- expMap[order(#expMap[[RRG$column]], # Do not use RRG here
  expMap[[RG$column]], # Cf. below: safe because each RG contains at least one ref samples group
  expMap$Reference, # Very important: if any level is dropped from the design matrix, it must be a reference!
  # (Otherwise every gets confusing and  my head starts hurting...)
  expMap[[VPAL$column]],
  expMap$Replicate),]
#
# Replace hyphens by dots to avoid issues with evaluating contrasts
for (Coeff in Coefficients) {
  l <- length(grep("-", Exp.map[[Coeff]])) # Not expMap in case we are re-running a small chunk
  if (l) {
    nuCoeff <- paste0(Coeff, "___")
    stopifnot(!nuCoeff %in% colnames(Exp.map)) # Not expMap in case we are re-running a small chunk
    Coefficients[which(Coefficients == Coeff)] <- nuCoeff
    expMap[[nuCoeff]] <- gsub("-", ".", expMap[[Coeff]])
  }
}
#
Group_ <- do.call(paste, c(expMap[, Coefficients, drop = FALSE], sep = "_"))
Group_ <- as.factor(Group_)
expMap$Group_ <- Group_
if (Nested) {
  expMap$Replicate_ <- Replicate_ <- as.factor(expMap$Replicate)
  designMatr %<o% model.matrix(~0 + Replicate_ + Group_)
} else {
  designMatr %<o% model.matrix(~0 + Group_)
}
rownames(designMatr) <- expMap$Ref.Sample.Aggregate
#
# Define contrasts
expContrasts %<o% list()
for (ratGrp in RG$values) { #ratGrp <- RG$values[1]
  em <- expMap[which(expMap[[RG$column]] == ratGrp),]
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
PolySTestRoot %<o% paste0("EBAM regulated-FDR=", paste(100*bhFDRs, collapse = "/"), "% FDR - ")
#
samDir <- paste0(wd, "/", samSubDir)
ebamDir <- paste0(wd, "/", ebamSubDir)
#
dataType <- "PG"
Src <- paste0(libPath, "/extdata/R scripts/Sources/Av_and_Stat_tests.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

### Visualize and check P-values
Src <- paste0(libPath, "/extdata/R scripts/Sources/pVal_check.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

useSAM %<o% ((names(pvalue.col)[which(pvalue.use)] == "Student")&&(useSAM_thresh))
if (useSAM) {
  # In this case, we bypass the original decision and base it off SAM even though we plot Student's P-values
  for (i in names(SAM_thresh)) { #i <- names(SAM_thresh)[1]
    dec <- SAM_thresh[[i]]$decision
    mKol <- rev(colnames(dec))[1]
    FCkol <- paste0("Mean ", Prot.Rat.Root, i)
    stopifnot(FCkol %in% names(PG))
    regKol <- paste0("Regulated - ", i)
    PG[[regKol]] <- "non significant"
    fdrs <- as.numeric(gsub("FDR$", "", colnames(dec)[which(colnames(dec) != mKol)]))
    fdrs <- sort(fdrs, decreasing = TRUE)
    for (f in fdrs) { #f <- fdrs[1]
      w <- which(PG[[mKol]] %in% dec[which(dec[[paste0(f, "FDR")]] == "+"), mKol])
      if (length(w)) {
        PG[which(PG[w, FCkol] > 0), regKol] <- paste0("up, FDR = ", f*100, "%")
        if (TwoSided) {
          PG[which(PG[w, FCkol] < 0), regKol] <- paste0("down, FDR = ", f*100, "%")
        }
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
Src <- paste0(libPath, "/extdata/R scripts/Sources/Annotate_me.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Get or generate subcellular localisation markers - for later use
Src <- paste0(libPath, "/extdata/R scripts/Sources/SubCellMark.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)

#### Code chunk - ROC analysis
if (length(ROC_GOterms)) {
  library(ggplot2)
  msg <- "ROC analysis"
  ReportCalls <- AddMsg2Report()
  ROC_GOterms <- unique(unlist(c(ROC_GOterms, GO_terms$Offspring[match(ROC_GOterms, GO_terms$ID)])))
  PG$"ROC - True Positive" <- FALSE
  PG$"ROC - True Positive"[grsep2(ROC_GOterms, PG$`GO-ID`)] <- TRUE
  if (sum(PG$"True Positive") < 10) {
    msg <- "Not enough TRUE positives for ROC analysis (min = 10), skipping!"
    ReportCalls <- AddMsg2Report(Warning = TRUE, Print = FALSE)
  } else {
    dir <- paste0(wd, "/ROC analysis")
    dirlist <- unique(c(dirlist, dir))
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    pack <- "pROC"
    bioc_req <- unique(c(bioc_req, pack))
    biocInstall(pack)
    w1 <- which(paste0(pvalue.col[which(pvalue.use)], VPAL$values) %in% colnames(PG))
    if (!length(w1)) { stop("Where are my P-values?!?!%&$+*!*") } else {
      for (grp in VPAL$values[w1]) { #grp <- VPAL$values[w1[1]]
        grp2 <- cleanNms(grp)
        pkol <- paste0(pvalue.col[which(pvalue.use)], grp)
        w2 <- which(is.all.good(PG[[pkol]], 2))
        PG$"ROC - Predictor" <- 10^(-PG[[pkol]])
        rocobj <- pROC::roc(PG[w2, ], "ROC - True Positive", "ROC - Predictor")
        ttl <- paste0("ROC analysis - ", grp2)
        plot <- pROC::ggroc(rocobj, color = "red") +
          geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dotted") +
          theme_bw() + ggtitle(ttl)
        poplot(plot)
        suppressMessages({
          ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300)
          ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)
        })
        ReportCalls <- AddPlot2Report(Title = gsub(": ?", " - ", ttl))
      }
    }
  }
  PG$"ROC - Predictor" <- NULL
  PG$"ROC - True Positive" <- NULL
}

#### Code chunk - Estimate P-value significance for a set of accepted FDRs
## NB: For graphical reasons (volcano plots), there is only support for 4 different FDR values. This should suffice anyway.
a <- sapply(strsplit(Param$Plot.metrics, ";"), function(x) { strsplit(x, ":") })
a[[2]][2] <- gsub("\\.$", "", pvalue.col[which(pvalue.use)])
Param$Plot.metrics <- paste(vapply(a, paste, "", collapse = ":"), collapse = ";")
FDR.thresholds %<o% c()

A <- VPAL$values
test <- vapply(A, function(x) { #x <- A[6]
  x <- paste0(pvalue.col[which(pvalue.use)], x)
  r <- x %in% colnames(PG)
  if (r) { r <- length(is.all.good(as.numeric(PG[[x]]))) > 0 }
  return(r)
}, TRUE)
A <- A[which(test)]
PG <- PG[, grep("^Significant-FDR=", colnames(PG), invert = TRUE)]
for (a in A) { #a <- A[1]
  temp <- FDR(data = PG,
              aggregate = a,
              pvalue_root = pvalue.col[which(pvalue.use)],
              fdr = BH.FDR, returns = c(TRUE, TRUE), method = "BH")
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
  if (length(pkol)) {
    msg <- "Adjusting P-values..."
    ReportCalls <- AddMsg2Report(Space = FALSE, Print = FALSE)
    for (pk in pkol) { #pk <- pkol[1]
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
    #                              slope = rep(0, 2),
    #                              xintercept = rep(NA, 2),
    #                              colour = colorRampPalette(c("orange", "red"))(2),
    #                              label = paste0(c(0.01, 0.05)*100, "% P-value"))
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " P-values were adjusted using the Benjamini-Hochberg (FDR) method.")
  } else { stop("There should be P-value columns in the protein groups table at this stage!") }
}

# Create list of control ratio values for the purpose of identifying vertical thresholds for plots:
## This should be, for each volcano plot, a list of ratios with the same name.
# Note: It's taken me forever to get this right, but I think this should be ok now!!!
if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
  plotMetr <- as.data.frame(strsplit(unlist(strsplit(Param$Plot.threshold.metrics, ";")), ":"))
  plotMetr <- as.data.frame(t(plotMetr)) 
  rownames(plotMetr) <- NULL
  colnames(plotMetr) <- c("Levels", "Axis")
  a2 <- set_colnames(as.data.frame(t(sapply(strsplit(unlist(strsplit(Param$Plot.threshold.values, split = "; *")), split = ": *"), unlist))),
                     c("Direction", "Text.value"))
  plotMetr$Text.value <- a2$Text.value[match(plotMetr$Levels, a2$Direction)]
  w <- which(plotMetr$Axis == "X")
  m <- w[match(c("down", "up"), plotMetr$Levels[w])]
  plotMetr$Text.value[w] <- as.character(c(-Param$Ratios.Contamination.Rates, Param$Ratios.Contamination.Rates))
  Param$Plot.threshold.values <- do.call(paste, c(plotMetr[, c("Levels", "Text.value")], sep = ": ", collapse = ";"))
  Ref.Ratios %<o% NULL
}
if (Param$Ratios.Thresholds == threshMsg) {
  Ref.Ratios %<o% setNames(lapply(VPAL$values, function(x) { #x <- VPAL$values[1]
    if (RatConGrps == "Ratio groups") {
      x1 <- unique(Exp.map[which(Exp.map[[VPAL$column]] == x), RG$column])
    }
    if (RatConGrps == "Experiments") {
      x1 <- unique(Exp.map$Experiment[which(Exp.map[[VPAL$column]] == x)])
      x1 <- unique(Exp.map[which(Exp.map$Experiment == x1), RG$column])
    }
    if (RatConGrps == "Whole dataset") {
      x1 <- unique(Exp.map[[RG$column]])
    }
    x <- unique(Exp.map[which(Exp.map[[VPAL$column]] == x), RG$column])
    x <- grep(paste0(topattern(paste0(Prot.Rat.Root, x1, "_REF.to.REF_")), "[0-9]+"), colnames(quant.data), value = TRUE)
    if (length(x)) { x <- is.all.good(as.numeric(unlist(quant.data[, x]))) } else { x <- NULL }
    return(x)
  }), VPAL$values)
}

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - t-test volcano plot(s)
#
PG$"1-PEP" <- 1 - PG$PEP
PG$"log10(1-PEP)" <- log10(PG$"1-PEP")
PG$"log10(Peptides count)" <- log10(PG$"Peptides count")
a <- grep(topattern(Prot.Expr.Root), colnames(PG), value = TRUE)
a <- grep("\\.REF$", a, value = TRUE, invert = TRUE)
PG$"Av. log10 abundance" <- apply(PG[, a], 1, function(x) { mean(is.all.good(unlist(x))) })
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
tempVP <- try(Volcano.plot(Prot = PG,
                           mode = "custom",
                           experiments.map = Exp.map,
                           X.root = paste0("Mean ", Prot.Rat.Root),
                           Y.root = pvalue.col[which(pvalue.use)],
                           aggregate.map = Aggregate.map,
                           aggregate.name = VPAL$aggregate,
                           aggregate.list = Aggregate.list, parameters = Param,
                           save = c("jpeg", "pdf"), labels = c("FDR", "both")[useSAM+1],
                           Ref.Ratio.values = Ref.Ratios,
                           Ref.Ratio.method = paste0("obs", RefRat_Mode),
                           ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                           FDR.thresh = FDR.thresholds,
                           arbitrary.lines = arbitrary.thr,
                           proteins = prot.list, proteins_split = protsplit,
                           return = TRUE, return.plot = TRUE,
                           title = "Volcano plot ",
                           subfolder = subDr,
                           subfolderpertype = FALSE, Symmetrical = TwoSided,
                           Alpha = "Rel. log10(Peptides count)",
                           Size = "Rel. av. log10 abundance", Size.max = 2,
                           plotly = create_plotly, plotly_local = create_plotly_local,
                           plotly_labels = PrLabKol,
                           cl = parClust,
                           SAM = useSAM, curved_Thresh = SAM_thresh, saveData = TRUE
                           ))
if (!class(tempVP) %in% c("try-error", "character")) {
  #
  # Save plotly plots
  dr <- paste0(wd, "/", subDr)
  myPlotLys <- tempVP$"Plotly plots"
  Src <- paste0(libPath, "/extdata/R scripts/Sources/save_Plotlys.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
  #
  VP_list <- tempVP
  insrt <- ""
  Src <- paste0(libPath, "/extdata/R scripts/Sources/thresholds_Excel.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
  #
  thresh <- lapply(names(tempVP$Thresholds$Absolute), function(x) {
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
  dms <- wb_dims(2, 1)
  wb <- wb_add_data(wb, "Thresholds", "Absolute thresholds", dms)
  wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                    italic = "true", underline = "single")
  dms <- wb_dims(3, 2)
  wb <- wb_add_data_table(wb, "Thresholds", thresh, dms,
                          col_names = TRUE, table_style = "TableStyleMedium2",
                          banded_rows = TRUE, banded_cols = FALSE)
  dms <- wb_dims(3+nrow(thresh)+3, 1)
  wb <- wb_add_data(wb, "Thresholds", "FDR thresholds", dms)
  wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                    italic = "true", underline = "single")
  dms <- wb_dims(3+nrow(thresh)+4, 2)
  wb <- wb_add_data_table(wb, "Thresholds", fdrThresh, dms,
                          col_names = TRUE, table_style = "TableStyleMedium2",
                          banded_rows = TRUE, banded_cols = FALSE)
  wb <- wb_set_col_widths(wb, "Thresholds", 1, 3)
  tmp1 <- rbind(colnames(thresh), thresh)
  colnames(tmp1) <- paste0("V", 1:ncol(tmp1))
  tmp2 <- rbind(colnames(fdrThresh), fdrThresh)
  colnames(tmp2) <- paste0("V", 1:ncol(tmp2))
  tst <- plyr::rbind.fill(tmp1, tmp2)
  tst <- setNames(apply(tst, 2, function(x) { max(nchar(x), na.rm = TRUE) }), NULL)
  wb <- wb_set_col_widths(wb, "Thresholds", 1:(length(tst)+1), c(3, tst))
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
} else { stop("MAJOR ERROR: No volcano plots were created, investigate!") }
# Also calculate Q-values - for now, the plot is created but not saved!
if (("Q.values" %in% colnames(Param))&&(is.logical(Param$Q.values))&&(Param$Q.values)) {
  require(qvalue)
  pkol <- grep(topattern(pvalue.col[which(pvalue.use)]), colnames(PG), value = TRUE)
  if (length(pkol)) {
    msg <- "Calculating Q-values..."
    ReportCalls <- AddMsg2Report(Space = FALSE, Print = FALSE)
    for (pk in pkol) { #pk <- pkol[2]
      temp <- 10^(-PG[[pk]])
      wag <- which(is.all.good(temp, 2))
      pi0 <- qvalue::pi0est(temp)$pi0[1]
      temp <- try(qvalue::qvalue(temp[wag]), silent = TRUE) # For now we do not explicitly set pi0
      if ("try-error" %in% class(temp)) {
        temp <- try(qvalue::qvalue(temp[wag], pi0 = pi0), silent = TRUE)
        while (("try-error" %in% class(temp))&&(pi0 <= 1)) {
          pi0 <- pi0 + 0.05
          temp <- try(qvalue::qvalue(temp[wag], pi0 = pi0), silent = TRUE)
        }
      }
      if (!"try-error" %in% class(temp)) {
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
    # Probably quite deprecated... check arguments before running
    subDr <- "Reg. analysis/t-tests"
    tempVP2 <- Volcano.plot(Prot = PG, mode = "standard-ratios", experiments.map = Exp.map,
                            X.root = paste0("Mean ", Prot.Rat.Root),
                            Y.root = "-log10(Qvalue) - ",
                            aggregate.map = Aggregate.map,
                            aggregate.name = VPAL$aggregate,
                            aggregate.list = Aggregate.list, parameters = Param,
                            save = FALSE, labels = "both",
                            Ref.Ratio.values = Ref.Ratios,
                            Ref.Ratio.method = paste0("obs", RefRat_Mode),
                            ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                            arbitrary.thresh = qval.thresh,
                            proteins = prot.list, proteins_split = protsplit,
                            return = TRUE, return.plot = TRUE,
                            title = "Q-values volcano plot ", subfolder = subDr,
                            subfolderpertype = FALSE, Symmetrical = TwoSided,
                            Alpha = "Rel. log10(Peptides count)",
                            Size = "Rel. av. log10 abundance", Size.max = 2,
                            plotly = create_plotly, plotly_local = create_plotly_local,
                            cl = parClust)
    #
    # Save plotly plots
    dr <- paste0(wd, "/", subDr)
    myPlotLys <- tempVP2$"Plotly plots"
    Src <- paste0(libPath, "/extdata/R scripts/Sources/save_Plotlys.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    VP_list <- tempVP2
    insrt <- "_Qvalues"
    Src <- paste0(libPath, "/extdata/R scripts/Sources/thresholds_Excel.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " Q-values were computed using package qvalue.")
  }
}

# Specificity mark for untested proteins
# We can assume that proteins with PSMs only in the specific pull-down samples are actually specifically enriched!
# -> Label them as such!
# NB: This used to be for pull-downs only but I am extending it to the whole workflow for now.
#     This could also be parameter-controlled.
# NB: For the F-test, this is done within the function
#if (IsPullDown) {
MinPep4Spec %<o% 2
if ("Min.pep.per.sample.for.spec" %in% colnames(Param)) {
  MinPep4Spec <- Param$Min.pep.per.sample.for.spec
  if (!is.numeric(MinPep4Spec)) {
    warning("Invalid value for parameter \"Min.pep.per.sample.for.spec\", defaulting to 2!")
    MinPep4Spec <- 2
  }
}
for (a in RG$values) { #a <- RG$values[1]
  e <- Exp.map[which(Exp.map[[RG$column]] == a),]
  grp1 <- unique(e[which(!e$Reference %in% c(TRUE, "TRUE")), VPAL$column])
  kolTR <- paste0("Regulated - ", grp1)
  w <- which(kolTR %in% colnames(PG))
  if (length(w)) {
    grp1 <- grp1[w]
    kolTR <- kolTR[w]
    grp0 <- unique(e[which(e$Reference %in% c(TRUE, "TRUE")), VPAL$column])
    e0 <- unique(e$Ref.Sample.Aggregate[which(e[[VPAL$column]] == grp0)])
    kole0 <- paste0("Evidences count - ", e0)
    kolp0 <- paste0("Peptides count - ", e0)
    for (i in seq_along(grp1)) { #i <- 1
      e1i <- unique(e$Ref.Sample.Aggregate[which(e[[VPAL$column]] == grp1[i])])
      kole1i <- paste0("Evidences count - ", e1i)
      kolp1i <- paste0("Peptides count - ", e1i)
      tst1i <- rowSums(PG[, kolp1i] >= MinPep4Spec, na.rm = TRUE) == length(kolp1i) # Are there at least MinPep4Spec peptidoforms in grp1i?...
      tst0 <- rowSums(PG[, kolp0] == 0, na.rm = TRUE) == length(kolp0) #... and concurrently none in grp0?
      w1i <- which(tst1i & tst0)
      if (length(w1i)) {
        pepmin <- apply(PG[w1i, kolp1i, drop = FALSE], 1, min)
        evcount <- apply(PG[w1i, kole1i, drop = FALSE], 1, sum)
        txtup <- paste0("Specific: at least ", pepmin, " pep./sample (", evcount, " PSMs tot.)")
        PG[w1i, kolTR[i]] <- txtup
      }
    }
  } else { warning("I would expect \"Regulated ...\" columns in the PG table by this stage!") }
}
#}

# Create t-test filters:
## These can then be used for further steps down the line, such as volcano plots, etc...
Reg_filters %<o% list()
filter_types %<o% tolower(unlist(strsplit(Param$Filters.type, ";")))
filter_types[grep("^dat.+2$", filter_types, invert = TRUE)] <- substr(filter_types[which(!grepl("^dat.+2$", filter_types))], 1, 3)
filter_types[grep("^dat.+2$", filter_types)] <- "dat2"
filter_types <- unique(c("con", filter_types))
if ("ref" %in% filter_types) {
  if ((RRG$aggregate != RG$aggregate)||(Nested)) {
    warning("Grouping filter by reference is not feasible if replicates are paired!")
    filter_types <- filter_types[which(filter_types != "ref")]
  } else {
    if (sum(vapply(RG$names, function(x) { !x %in% RSA$names }, TRUE)) > 0) {
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
  g2$Ref <- vapply(tmp, function(x) { unique(tst[which((Exp.map$Reference)&(tst == x))]) }, "")
  for (i in unique(g2$Ref)) {
    w <- which(g2$Ref == i)
    u <- grep("^up|^Specific", unique(as.character(PG[, g[w]])), value = TRUE)
    d <- grep("^down", unique(as.character(PG[, g[w]])), value = TRUE)
    Reg_filters$"t-tests"$"By reference"[[i]] <- list(Columns = g[w],
                                                      Filter_up = sort(which(apply(PG[, g[w], drop = FALSE], 1, function(x) {
                                                        length(which(x %in% up))
                                                      }) > 0)),
                                                      Filter_down = sort(which(apply(PG[, g[w], drop = FALSE], 1, function(x) {
                                                        length(which(x %in% down))
                                                      }) > 0)),
                                                      Filter = sort(which(apply(PG[, g[w], drop = FALSE], 1, function(x) {
                                                        length(which(x %in% c(up, down)))
                                                      }) > 0)))
  }
}
if (sum(c("dat", "dat2") %in% filter_types)) {
  Reg_filters$"t-tests"$"Whole dataset" <- list(Columns = g,
                                                Filter_up = sort(which(apply(PG[, g, drop = FALSE], 1, function(x) {
                                                  length(which(x %in% up))
                                                }) > 0)),
                                                Filter_down = sort(which(apply(PG[, g, drop = FALSE], 1, function(x) {
                                                  length(which(x %in% down))
                                                }) > 0)),
                                                Filter = sort(which(apply(PG[, g, drop = FALSE], 1, function(x) {
                                                  length(which(x %in% c(up, down)))
                                                }) > 0)))
}
#
# Z-scored clustering heatmaps of regulated proteins
clustMode <- "t-tests"
clstSrc <- paste0(libPath, "/extdata/R scripts/Sources/cluster_Heatmap_Main.R")
#rstudioapi::documentOpen(clstSrc)
source(clstSrc, local = FALSE)
#

rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)
gc()
try({ stopCluster(parClust) }, silent = TRUE)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - F-test
#Param <- Param.load()
F.test %<o% FALSE
if (("F.test" %in% colnames(Param))&&(Param$F.test)) {
  if ((length(VPAL$values) == 2)&&(pvalue.col[pvalue.use] == "Moderated t-test -log10(Pvalue) - ")) {
    F.test <- FALSE
  } else {
    F.test <- TRUE
  }
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
  # tmp1 <- unlist(strsplit(Param$F.test_factors, "_\\|_"))
  # tmp2 <- unlist(strsplit(Param$F.test_factors_ref, "_\\|_"))
  # stopifnot(length(tmp2) == length(tmp2))
  # tmp1 <- strsplit(tmp1, "_;_")
  # tmp2 <- strsplit(tmp2, "_;_")
  # stopifnot(sum(!vapply(tmp1, function(x) { length(x) %in% c(2:3)[Nested+1] }, TRUE)) == 0,
  #           sum(!vapply(tmp2, function(x) { length(x) == 2 }, TRUE)) == 0) # If this breaks, this will mean that these parameters are not built the way I remember them to be.
  # # Anyway, this code is in dire need of a refresher! It would be a welcome occasion to fix it.
  # #
  # tmp1 <- lapply(tmp1, function(x) { Factors[x] })
  expMap_F %<o% Exp.map[match(rownames(designMatr), Exp.map$Ref.Sample.Aggregate),]
  # Replace hyphens by dots to avoid issues with evaluating contrasts
  for (nuCoeff in Coefficients) {
    Coeff <- gsub("___$", "", nuCoeff)
    stopifnot(Coeff %in% colnames(Exp.map))
    l <- length(grep("-", Exp.map[[Coeff]])) # Not expMap_F in case we are re-running a small chunk
    if (l) {
      stopifnot(!nuCoeff %in% colnames(Exp.map)) # Not expMap_F in case we are re-running a small chunk
      expMap_F[[nuCoeff]] <- gsub("-", ".", expMap_F[[Coeff]])
    }
  }
  Group_ <- do.call(paste, c(expMap_F[, Coefficients, drop = FALSE], sep = "_"))
  Group_ <- as.factor(Group_)
  expMap_F$Group_ <- Group_
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
    tmp2 <- unlist(sapply(1:(l-1), function(x) {
      paste0("(", tmp[x], ") - (", tmp[(x+1):l], ")")
    }))
    tmp2Nm <- unlist(sapply(1:(l-1), function(x) {
      paste0("(", expContrasts_F$name[x], ") - (", expContrasts_F$name[(x+1):l], ")")
    }))
    tmp2Tbl <- data.frame(Contrasts = tmp2,
                          x1 = gsub("^\\(|\\) - \\(.*", "", tmp2),
                          x0 = gsub(".*\\) - \\(|\\)$", "", tmp2),
                          name = tmp2Nm,
                          Type = "Double")
    tmp2Tbl$All <- lapply(1:nrow(tmp2Tbl), function(x) { gsub("^Group_", "", unlist(strsplit(unlist(tmp2Tbl[x, c("x1", "x0")]), " - "))) })
    tmp <- c(tmp, tmp2)
    expContrasts_F <- rbind(expContrasts_F, tmp2Tbl)
  }
  #
  F_Root %<o% "mod. F-test -log10(Pvalue)"
  dataType <- "PG"
  #
  FSrc %<o% paste0(libPath, "/extdata/R scripts/Sources/run_F_test.R")
  #rstudioapi::documentOpen(FSrc)
  tstFtst <- try(source(FSrc, local = FALSE), silent = TRUE)
  #
  if (!"try-error" %in% class(tstFtst)) {
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
    # Also a posteriori F-test P-values histogram:
    nbin <- 20
    bd <- (0:nbin)/nbin
    if (F_Root %in% colnames(PG)) {
      temp <- data.frame(value = is.all.good(10^(-PG[[F_Root]])))
      ttl <- "Histogram: F-test moderated Pvalue"
      plot <- ggplot(temp, aes(x = value)) +
          geom_histogram(bins = nbin, colour = "black", alpha = 0.25, fill = "green") +
          guides(fill = "none") + theme_bw() + ggtitle(ttl)
      poplot(plot)
      dir <- paste0(wd, "/Workflow control/Protein groups/P-values")
      dirlist <- unique(c(dirlist, dir))
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      ttla <- gsub(": ?", " - ", ttl)
      suppressMessages({
        ggsave(paste0(dir, "/", ttla, ".jpeg"), plot, dpi = 300)
        ggsave(paste0(dir, "/", ttla, ".pdf"), plot, dpi = 300)
      })
      ReportCalls <- AddPlot2Report(Title = ttla)
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
                                                    Filter_up = sort(which(apply(F_test_data[, g, drop = FALSE], 1, function(x) {
                                                      length(which(x %in% up))
                                                    }) > 0)),
                                                    Filter_down = sort(which(apply(F_test_data[, g, drop = FALSE], 1, function(x) {
                                                      length(which(x %in% down))
                                                    }) > 0)),
                                                    Filter = sort(which(apply(F_test_data[, g, drop = FALSE], 1, function(x) {
                                                      length(which(x %in% c(up, down)))
                                                    }) > 0)))
    }
    if (("Q.values" %in% colnames(Param))&&(Param$Q.values)) {
      require(qvalue)
      if (F_Root %in% colnames(F_test_data)) {
        temp <- 10^(-F_test_data[[F_Root]])
        wag <- which(is.all.good(temp, 2))
        pi0 <- qvalue::pi0est(temp)$pi0[1]
        temp <- try(qvalue::qvalue(temp[wag]), silent = TRUE) # For now we do not explicitly set pi0
        if ("try-error" %in% class(temp)) {
          temp <- try(qvalue::qvalue(temp[wag], pi0 = pi0), silent = TRUE)
          while (("try-error" %in% class(temp))&&(pi0 <= 1)) {
            pi0 <- pi0 + 0.05
            temp <- try(qvalue(temp[wag], pi0 = pi0), silent = TRUE)
          }
        }
        if (!"try-error" %in% class(temp)) {
          F_test_data[["-log10(Qvalue)"]] <- NA
          F_test_data[["local FDR"]] <- NA
          F_test_data[wag, "-log10(Qvalue)"] <- -log10(temp$qvalues)
          F_test_data[wag, "local FDR"] <- temp$lfdr
        } else { warning(paste0("F-test: Q-values calculation failed for ", F_Root, "! No q-values column will be created.")) }
      }
    }
    #
    # Z-scored clustering heatmaps of regulated proteins
    clustMode <- "F-tests"
    Src <- paste0(libPath, "/extdata/R scripts/Sources/cluster_Heatmap_Main.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    rm(list = ls()[which(!ls() %in% .obj)])
    Script <- readLines(ScriptPath)
    gc()
    invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
    saveImgFun(BckUpFl)
    #loadFun(BckUpFl)
    source(parSrc, local = FALSE)
  } else { warning("F-test analysis failed, check your parameters!")}
}
# Mat-meth text
tmp <- BH.FDR*100
l <- length(tmp)
if (l > 1) { tmp <- paste0(paste(tmp[1:(l-1)], collapse = "%, "), " and ", tmp[l], "%") }
tmp2 <- Param$Ratios.Contamination.Rates
tmpPVal <- gsub(" -log10\\(pvalue\\) - ", "", gsub("welch", "Welch", tolower(pvalue.col[which(pvalue.use)])))
if (grepl("^moderated", tmpPVal, ignore.case = TRUE)) { tmpPVal <- paste0(tmpPVal, " (limma)") }
if (grepl("^((EBA)|(S))AM ", tmpPVal, ignore.case = TRUE)) { tmpPVal <- paste0(tmpPVal, " (siggenes)") }
if (grepl("^permutations ", tmpPVal, ignore.case = TRUE)) { tmpPVal <- paste0(tmpPVal, " (coin)") }
if (grepl("^(ODP)|(LRT) ", tmpPVal, ignore.case = TRUE)) { tmpPVal <- paste0(tmpPVal, " (edge)") }
DatAnalysisTxt <- paste0(DatAnalysisTxt, " Average log10 expression values were tested for significance using a ",
                         c("two", "one")[match(AltHyp, c("two.sided", "greater", "lower"))], "-sided ",
                         tmpPVal, " per samples group",
                         c("", " and a moderated ANOVA (limma, F-test run using voomaLmFit for heteroskedasticity correction with individual moderated t-tests as post-hoc tests)")[F.test+1],
                         ". Significance thresholds were calculated using the Benjamini-Hochberg procedure for False Discovery Rate (FDR) values of ", tmp,
                         ". For all tests, differentially expressed protein groups were defined as those with a significant P-value and a",
                         c("n absolute", "")[IsPullDown+1],
                         " log2 average ratio greater than ",
                         c(paste0(Param$Ratios.Contamination.Rates*100, "% of ",
                                  c(paste0("control-to", c("-average", "")[Nested + 1],
                                           "-control"),
                                    "intra-sample groups")[match(RefRat_Mode, c("1", "2"))],
                                  " ratios"),
                           Param$Ratios.Contamination.Rates)[match(Param$Ratios.Thresholds,
                                                                   threshOpt)], ".")

#### Code chunk - SAINTexpress
Src <- paste0(libPath, "/extdata/R scripts/Sources/SAINTexpress.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Now let's create a table of regulated protein groups per test made:
g <- grep("^Regulated - ", colnames(PG), value = TRUE)
regPG_TTest <- data.frame(Test = gsub("^Regulated - ", "", g))
if (TwoSided) { dir <- c("up", "down") } else { dir <- "up" }
kolstms <- c("count", "PG IDs", "Leading Protein IDs", "Genes")
for (d in dir) { #d <- "up"
  for (i in seq_along(BH.FDR)) { #i <- 1
    tmp <- paste0(d, ", FDR = ", BH.FDR[1:i]*100, "%")
    kolnms <- paste0(tmp[i], " - ", kolstms)
    tmp2 <- set_colnames(Isapply(g, function(x) {
      w <- which(PG[[x]] %in% tmp)
      x1 <- length(w)
      x2 <- paste0(PG$id[w], collapse = ", ")
      x3 <- paste0(PG$"Leading protein IDs"[w], collapse = ", ")
      x4 <- paste0(PG$Genes[w], collapse = ", ")
      return(c(x1, x2, x3, x4))
    }), kolnms)
    tmp2[[kolnms[1]]] <- as.numeric(tmp2[[kolnms[1]]])
    tst <- unique(tmp2[[kolnms[1]]])
    tst <- tst[which(tst > 0)]
    if (length(tst)) { regPG_TTest[, kolnms] <- tmp2 }
  }
}
if (IsPullDown) {
  tmp <- grep("^Specific: ", unique(unlist(PG[,g])), value = TRUE)
  if (length(tmp)) {
    kolnms <- paste0("Specific - ", kolstms)
    tmp2 <- set_colnames(Isapply(g, function(x) {
      w <- which(PG[[x]] %in% tmp)
      x1 <- length(w)
      x2 <- paste0(PG$id[w], collapse = ", ")
      x3 <- paste0(PG$"Leading protein IDs"[w], collapse = ", ")
      x4 <- paste0(PG$Genes[w], collapse = ", ")
      return(c(x1, x2, x3, x4))
    }), kolnms)
    tmp2[[kolnms[1]]] <- as.numeric(tmp2[[kolnms[1]]])
    tst <- unique(tmp2[[kolnms[1]]])
    tst <- tst[which(tst > 0)]
    if (length(tst)) { regPG_TTest[, kolnms] <- tmp2 }
  }
}
if (ncol(regPG_TTest) > 1) {
  dir <- paste0(wd, "/Tables")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  regPG_TTest <- regPG_TTest[, c("Test", as.character(sapply(kolstms, function(x) {
    grep(paste0(" - ", x, "$"), colnames(regPG_TTest), value = TRUE)
  })))]
  write.csv(regPG_TTest, file = paste0(dir, "/Reg. PGs - t-test.csv"), row.names = FALSE)
} else { warning("The t-test(s) did not identify any regulated protein groups!") }
if (F.test) {
  g <- grep("^mod\\. F-test Regulated - ", colnames(F_test_data), value = TRUE)
  regPG_FTest <- data.frame(Test = gsub("^Regulated - ", "", g))
  # For the F-test we always want to include both directions because there can be up and down for a pull-down if doing a secondary comparison
  for (d in c("down", "up")) { #d <- "down"
    for (i in seq_along(BH.FDR)) { #i <- 1
      tmp <- paste0(d, ", FDR = ", BH.FDR[1:i]*100, "%")
      kolnms <- paste0(tmp[i], " - ", kolstms)
      tmp2 <- set_colnames(Isapply(g, function(x) { #x <- g[1]
        w <- which(F_test_data[[x]] %in% tmp)
        x1 <- length(w)
        x2 <- paste0(PG$id[w], collapse = ", ")
        x3 <- paste0(PG$"Leading protein IDs"[w], collapse = ", ")
        x4 <- paste0(PG$Genes[w], collapse = ", ")
        return(c(x1, x2, x3, x4))
      }), kolnms)
      tmp2[[kolnms[1]]] <- as.numeric(tmp2[[kolnms[1]]])
      tst <- unique(tmp2[[kolnms[1]]])
      tst <- tst[which(tst > 0)]
      if (length(tst)) { regPG_FTest[, kolnms] <- tmp2 }
    }
  }
  if (IsPullDown) {
    tmp <- grep("^Specific: ", unique(unlist(F_test_data[,g])), value = TRUE)
    if (length(tmp)) {
      kolnms <- paste0("Specific - ", kolstms)
      tmp2 <- set_colnames(Isapply(g, function(x) {
        w <- which(F_test_data[[x]] %in% tmp)
        x1 <- length(w)
        x2 <- paste0(PG$id[w], collapse = ", ")
        x3 <- paste0(PG$"Leading protein IDs"[w], collapse = ", ")
        x4 <- paste0(PG$Genes[w], collapse = ", ")
        return(c(x1, x2, x3, x4))
      }), kolnms)
      tmp2[[kolnms[1]]] <- as.numeric(tmp2[[kolnms[1]]])
      tst <- unique(tmp2[[kolnms[1]]])
      tst <- tst[which(tst > 0)]
      if (length(tst)) { regPG_FTest[, kolnms] <- tmp2 }
    }
  }
  if (ncol(regPG_FTest) > 1) {
    dir <- paste0(wd, "/Tables")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    regPG_FTest <- regPG_FTest[, c("Test", as.character(sapply(kolstms, function(x) {
      grep(paste0(" - ", x, "$"), colnames(regPG_FTest), value = TRUE)
    })))]
    write.csv(regPG_FTest, file = paste0(dir, "/Reg. PGs - F-test.csv"), row.names = FALSE)
  } else { warning("The F-test did not identify any regulated protein groups!") }
}
# To do: SAINTexpress!

# Summary table and heatmap of number of regulated protein groups
Tsts <- c("t-tests", "F-tests")
WhTsts <- which(Tsts %in% names(Reg_filters))
for (tt in WhTsts) { #tt <- WhTsts[1]
  tstrt <- Tsts[tt]
  stopifnot(!is.na(tstrt))
  dir <- paste0(wd, "/Reg. analysis/", tstrt)
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  filt <- Reg_filters[[tstrt]]
  By <- c("By condition", "By reference", "By analysis")
  By <- By[which(By %in% names(filt))]
  if (length(By)) {
    for (bee in By) { #bee <- By[1]
      flt <- filt[[bee]]
      if (length(flt) >= 2) {
        flt <- flt[order(names(flt))]
        N <- length(flt)
        temp <- as.data.frame(matrix(rep("", (N+1)^2), ncol = N+1))
        temp[2:(N+1), 1] <- temp[1, 2:(N+1)] <- names(flt)
        for (i in 2:(N+1)) {
          x <- flt[[temp[i, 1]]]$Filter
          temp[i, 2:(N+1)] <- vapply(temp[1, 2:(N+1)], function(y) { sum(x %in% flt[[y]]$Filter) }, 1)
        }
        nms <- cleanNms(names(flt))
        tst <- vapply(strsplit(nms, " - "), length, 1)
        tst <- (min(tst) > 1)&(length(unique(tst)) == 1)
        if (tst) {
          tst <- as.data.frame(t(sapply(strsplit(nms, " - "), unlist)))
          l <- apply(tst, 2, function(x) { length(unique(x)) })
          tst <- tst[, which(l > 1), drop = FALSE]
          nms <- do.call(paste, c(tst, sep = " - "))
        }
        temp[2:(N+1), 1] <- temp[1, 2:(N+1)] <- nms
        nm <- paste0("N. of co-regulated PGs\n", tstrt, "\n(", tolower(bee), ")")
        write.csv(temp, file = paste0(dir, "/", gsub("\n", " - ", gsub("\n\\(", " (", nm)), ".csv"), row.names = FALSE)
        temp2 <- temp[2:(N+1), 2:(N+1)]
        colnames(temp2) <- temp[1, 2:(N+1)]
        rownames(temp2) <-  temp[2:(N+1), 1]
        for (i in 1:nrow(temp2)) { temp2[[i]] <- as.numeric(temp2[[i]]) }
        if (max(is.all.good(unlist(temp2))) > 0) {
          temp2 <- as.matrix(temp2)
          basic.heatmap(temp2,
                        "N. of co-regulated PGs",
                        paste0(tstrt, "\n(", tolower(bee), ")"),
                        save = c("pdf", "jpeg"),
                        folder = dir)
        } else { warning(paste0("Not a single regulated protein group in any of the ", tstrt, " performed, skipping.")) }
      } else {
        bb <- gsub("By ", "", bee)
        if (bb == "condition") { msg <- paste0(tstrt, " analysis: only one ", bb, " tested, skipping.") } else {
          msg <- paste0(tstrt, ": only one ", bb, " ", c("performed", "used")[(bb == "reference")+1], " -> skipping.")
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
  Src <- paste0(libPath, "/extdata/R scripts/Sources/GSEA.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

# Weighted Gene Correlation Networks Analysis (WGCNA)
if (!exists("runWGCNA")) { runWGCNA %<o% FALSE }
if (runWGCNA) {
  Src <- paste0(libPath, "/extdata/R scripts/Sources/WGCNA.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

#### Code chunk - Heatmaps with clustering at samples and protein groups level, highlighting proteins of interest
clustMode <- "standard"
Src <- paste0(libPath, "/extdata/R scripts/Sources/cluster_Heatmap_Main.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

try({ stopCluster(parClust) }, silent = TRUE)
source(parSrc, local = FALSE)
rm(list = ls()[which(!ls() %in% .obj)])
gc()
Script <- readLines(ScriptPath)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
Script <- readLines(ScriptPath)

#### Code chunk - Dimensionality reduction plots
dmrdSrc <- paste0(libPath, "/extdata/R scripts/Sources/dimRed_plots.R")
#rstudioapi::documentOpen(dmrdSrc)
source(dmrdSrc, local = FALSE)

#### Code chunk - Protein group profile plots and ranked abundance plots
Src <- paste0(libPath, "/extdata/R scripts/Sources/profile_and_rankedAbund_plots.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Visualize results
if (!exists("xplorSrc")) {
  xplorSrc <- paste0(libPath, "/extdata/R scripts/Sources/xplorData.R") # Backwards compatibility
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
  for (i in A) { #i <- A[1]
    i1 <- unlist(strsplit(i, "___"))
    e <- lapply(seq_along(o), function(x) { which(Exp.map[[o[x]]] == i1[x]) })
    l <- unique(unlist(e))
    t <- vapply(l, function(x) { length(which(unlist(e) == x)) == length(o) }, TRUE)
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
      if (length(t1) <= 1) {
        if (!length(t1)) { cat("   There is no valid data for aggregate", i, "\n")
        } else { cat("   There is only a single time point for aggregate", i, "\n") }
      } else {
        test <- apply(PG[,c(t1, t2)], 1, function(x) {length(is.all.good(x)) == length(tp)*2})
        col <- c("Protein IDs", "Names", "ID")
        col <- col[which(col %in% colnames(PG))]
        temp1 <- PG[which(test), c(col, Param$Plot.labels, t1)]
        temp1$IDs <- as.character(c(1:nrow(temp1)))
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
    kount <- 0
    for (tp in A) {
      if (tp %in% names(temp)) {
        kount <- kount + 1
        if (kount == 1) {
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
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 600, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 600, width = 10, height = 10, units = "in")
    })
    ReportCalls <- AddPlot2Report()
    if (create_plotly) {
      #test <- aggregate(tmp$`log2(Ratio)`, list(tmp$IDs), function(x) {length(is.all.good(x)) == length(Tim)-1})
      tmp2 <- aggregate(tmp$`log2(Ratio)`, list(tmp$IDs), function(x) { max(abs(is.all.good(x)))})
      tmp2 <- tmp2[order(tmp2$x, decreasing = TRUE),]
      tmp2 <- tmp2$Group.1[1:min(1000, nrow(tmp2))]
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
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)

#### Code chunk - Sub-Cellular localisation analysis
# Includes:
# - pRoloc-based prediction of localisation
# - Re-localisation analysis based on testing the Sums of Squared Differences of Profiles
Src <- paste0(libPath, "/extdata/R scripts/Sources/SubCellLoc.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

### Check that Cytoscape is installed and can run, then launch it.
Src <- paste0(libPath, "/extdata/R scripts/Sources/Cytoscape_init.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
runClueGO <- runClueGO&CytoScape

#### Code chunk - Gene Ontology terms enrichment analysis
if (enrichGO||globalGO) {
  msg <- "Gene Ontology terms enrichment analysis"
  ReportCalls <- AddMsg2Report(Space = FALSE)
  packs <- c("GO.db", "topGO")
  for (pack in packs) {
    bioc_req <- unique(c(bioc_req, pack))
    biocInstall(pack)
  }
  GO_enrich.dat %<o% list()
  GO_enrich.FCRt %<o% list()
  GO_enrich.tbl %<o% list()
  GO_Plots %<o% list()
  Reg_GO_terms %<o% list()
  GO.enrich.MultiRefs %<o% (("GO.enrichment.Ref.Aggr" %in% colnames(Param))&&(!Param$GO.enrichment.Ref.Aggr %in% c("", "NA", NA)))
  if (GO.enrich.MultiRefs) { parse.Param.aggreg.2("GO.enrichment.Ref.Aggr") }
  #
  if (runClueGO) {
    # Initialize ClueGO
    Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_init.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
  #
  if (enrichGO) {
    Tsts <- c("t-tests", "F-tests", "Localisation", "SAINTexpress")
    WhTsts <- which(Tsts %in% names(Reg_filters))
    if (!exists("SSD.Root")) { SSD.Root <- "" }
    for (tt in WhTsts) { #tt <- WhTsts[1] #tt <- WhTsts[2] #tt <- WhTsts[3] #tt <- WhTsts[4]
      tstrt <- Tsts[tt]
      stopifnot(!is.na(tstrt))
      ReportCalls <- AddMsg2Report(Msg = paste0("\n - ", tstrt), Space = FALSE)
      dir <- paste0(wd, "/Reg. analysis/GO enrich/", tstrt)
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      filt <- Reg_filters[[tstrt]]
      #By <- c("By condition", "By reference", "By analysis", "Whole dataset")
      By <- "By condition"
      By <- By[which(By %in% names(filt))]
      if (length(By)) {
        for (bee in By) { #bee <- By[1]
          flt <- filt[[bee]]
          if (bee == "Whole dataset") { flt <- list("Whole dataset" = flt) }
          tstbee <- paste0(tstrt, "_", tolower(bee))
          if (length(flt)) {
            flt <- flt[order(names(flt))]
            reg <- setNames(lapply(flt, function(x) { list(x$Columns) }), names(flt))
            reg <- set_colnames(reshape::melt(reg), c("Name", "Bleh", "For"))
            reg$Bleh <- NULL
            fcr <- c(paste0("Mean ", c(Prot.Rat.Root, "log2(Ratio) - ", SSD.Root)), "log2(FC) - ")[tt]
            reg$ParentFC <- gsub(paste0(".*Re", c(rep("gulated", 2), "-localized", "gulated")[tt], " - "), fcr, reg$Name)
            reg$FCname <- paste0(fcr, reg$For)
            tmpdat <- get(c("PG", "F_test_data", "PG", "allSAINTs")[tt])
            UF <- unique(reg$For)
            temp <- as.data.frame(sapply(UF, function(x) { #x <- UF[1]
              x <- reg$ParentFC[which(reg$For == x)]
              if (length(x) > 1) {
                x <- apply(tmpdat[, x], 1, log_ratio_av)
              } else { x <- tmpdat[[x]] }
              return(x)
            }))
            colnames(temp) <- paste0(fcr, UF)
            for (i in 1:nrow(reg)) {
              Reg_filters[[tstrt]][[bee]][[reg$For[i]]]$Ratios <- temp[[paste0(fcr, reg$For[i])]]
            }
            if (tt %in% 1:3) {
              Kol2Add <- c("Leading protein IDs", "Protein IDs", "id", "Protein names", "No Isoforms", "Names", "Genes",
                           "Common Names", Param$Plot.labels, "GO", "GO-ID")
              Kol2Add <- Kol2Add[which(Kol2Add %in% colnames(PG))]
              temp[, Kol2Add] <- PG[, Kol2Add]
            } else {
              Kol2Add <- c("No Isoforms", "Gene", "Common Name", "GO", "GO-ID")
              temp$Protein <- allSAINTs$Protein
              temp$PG_id <- allSAINTs$PG_id
              temp[, Kol2Add] <- db[match(allSAINTs$Protein, db$`Protein ID`), Kol2Add]
            }
            GO_enrich.dat[[tstbee]] <- temp
            GO_enrich.FCRt[[tstbee]] <- fcr
            GO_enrich.tbl[[tstbee]] <- reg
            flt <- setNames(lapply(UF, function(x) { flt[[x]]$Filter }), UF)
            #flt <- setNames(lapply(UF, function(x) { flt[[x]]$Filter_down }), UF)
            #flt <- setNames(lapply(UF, function(x) { flt[[x]]$Filter_up }), UF)
            # see function code for defaults)
            ttr <- btr <- ""
            if (length(By) > 1) { ttr <- btr <- paste0(tolower(bee), "_") }
            Ref.Filt <- tmpFilt <- setNames(lapply(names(flt), function(x) { 1:nrow(GO_enrich.dat[[tstbee]]) }), names(flt))
            if (tt %in% c(1, 4)) {
              if (GO.enrich.MultiRefs) {
                Ref.Filt <- try(setNames(lapply(names(flt), function(x) { #x <- names(flt)[1]
                  m <- Exp.map[which(Exp.map[[VPAL$column]] == x), , drop = FALSE]
                  y <- unique(m[[GO.enrichment.Ref.Aggr$column]])
                  stopifnot(length(y) == 1)
                  w <- which(Exp.map[[GO.enrichment.Ref.Aggr$column]] == y)
                  w1 <- which(apply(PG[, paste0(Prot.Expr.Root, Exp.map$Ref.Sample.Aggregate[w])], 1, function(x) {
                    sum(is.all.good(x, 2))
                  }) > 0)
                  if (tt == 4) {
                    w1 <- which(temp$PG_id %in% PG$id[w1])
                  }
                  w2 <- flt[[x]] # Required for if we are imputing missing values
                  return(sort(unique(c(w1, w2))))
                }), names(flt)), silent = TRUE)
                if ("try-error" %in% class(Ref.Filt)) {
                  warning("Invalid \"GO.enrichment.Ref.Aggr\" argument: multiple references for GO enrichment are only feasible if each enrichment filter maps to a single reference! Skipping...")
                  GO.enrich.MultiRefs <- FALSE
                  Ref.Filt <- tmpFilt
                }
              }
            }
            if (tt == 2) {
              expMap <- Exp.map
              # Replace hyphens by dots to avoid issues with evaluating contrasts
              for (nuCoeff in Coefficients) {
                Coeff <- gsub("___$", "", nuCoeff)
                stopifnot(Coeff %in% colnames(Exp.map))
                l <- length(grep("-", Exp.map[[Coeff]])) # Not expMap in case we are re-running a small chunk
                if (l) {
                  stopifnot(!nuCoeff %in% colnames(Exp.map)) # Not expMap in case we are re-running a small chunk
                  expMap[[nuCoeff]] <- gsub("-", ".", expMap[[Coeff]])
                }
              }
              Group_ <- do.call(paste, c(expMap[, Coefficients, drop = FALSE], sep = "_"))
              Group_ <- as.factor(Group_)
              expMap$Group_ <- Group_
              cM <- as.data.frame(contrMatr_F)
              Ref.Filt <- setNames(lapply(names(flt), function(x) { #x <- names(flt)[1]
                x <- expMap$Ref.Sample.Aggregate[which(expMap$Group_ %in% unlist(expContrasts_F$All[match(x, expContrasts_F$name)]))]
                kol <- paste0(Prot.Expr.Root, x)
                which(parApply(parClust, PG[, kol], 1, function(x) { length(proteoCraft::is.all.good(x)) }) > 0)
              }), names(flt))
            }
            if ((length(Ref.Filt) > 1)||(!is.na(Ref.Filt))) {
              flt <- setNames(lapply(names(flt), function(x) { flt[[x]][which(flt[[x]] %in% Ref.Filt[[x]])] }),
                              names(flt))
            }
            # Also save the reference filters 
            nms <- names(Reg_filters[[tstrt]][[bee]])
            for (nm in nms) {
              Reg_filters[[tstrt]][[bee]][[nm]]$Background_filter <- Ref.Filt[[nm]]
            }
            #
            Mode <- "regulated"
            if (tt %in% c(1, 3)) { dataType <- "PG" }
            if (tt == 4) { dataType <- "Prot" }
            if ((!exists("GO_mappings"))&&(file.exists("GO_mappings.RData"))) {
              loadFun("GO_mappings.RData")
            }
            if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) { loadFun("GO_terms.RData") }
            #
            Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_enrich.R")
            #rstudioapi::documentOpen(Src)
            source(Src, local = FALSE)
            #
            if (runClueGO) {
              clueGO_outDir <- dir
              clueGO_type <- "Enrichment (Right-sided hypergeometric test)"
              Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_enrich.R")
              #rstudioapi::documentOpen(Src)
              source(Src, local = FALSE)
            }
            #
            # Cleanup - do it now, not within sources!
            try(rm(list = allArgs), silent = TRUE)
            #
            GO_Plots[[tstbee]] <- goRES
            #
            if (!is.null(names(GO_Plots[[tstbee]]))) {
              if ("All_GO_terms" %in% names(GO_Plots[[tstbee]])) {
                GO_terms <- GO_Plots[[tstbee]]$All_GO_terms
                GO_Plots[[tstbee]]$All_GO_terms <- NULL
              }
              n2 <- names(GO_Plots[[tstbee]]$GO_plots)
              dir2 <- paste0(wd, "/", dir)
              for (ttl in n2) { #ttl <- n2[1]
                plot <- GO_Plots[[tstbee]]$GO_plots[[ttl]]
                ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
              }
              if ((create_plotly)&&(!create_plotly_local)) { plot_ly[[paste0("GO plots - Regulated vs Observed - ", tstbee)]] <- GO_Plots[[tstbee]]$GO_plot_ly }
              temp <- GO_Plots[[tstbee]]$GO_terms
              temp$Mapping <- NULL
              if ("Offspring" %in% colnames(temp)) {
                temp$Offspring <- vapply(temp$Offspring, paste, "", collapse = ";")
              }
              temp$`Protein table row(s)` <- NULL
              colnames(temp) <- cleanNms(colnames(temp))
              gn <- grep("^Genes", colnames(temp), value = TRUE)
              pr <- grep("^Proteins", colnames(temp), value = TRUE)
              pg <- grep("^PG IDs", colnames(temp), value = TRUE)
              kn <- grep("^Count", colnames(temp), value = TRUE)
              pv <- grep("^Pvalue", colnames(temp), value = TRUE)
              zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
              lf <- grep("^logFC", colnames(temp), value = TRUE)
              si <- grep("^Significance", colnames(temp), value = TRUE)
              #lp <- grep("^Leading protein IDs", colnames(temp), value = TRUE)
              #kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, lp))]
              kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, pr))]
              #temp <- temp[, c(kl, si, gn, pg, lp, kn, pv, zs, lf)]
              temp <- temp[, c(kl, si, gn, pg, pr, kn, pv, zs, lf)]
              w <- apply(temp[, pv, drop = FALSE], 1, function(x) { sum(!is.na(x)) }) > 0
              temp <- temp[w,]
              tst <- apply(temp[, kn, drop = FALSE], 1, function(x) { sum(x[which(!is.na(x))]) })
              temp <- temp[order(tst, decreasing = TRUE),]
              temp <- temp[order(temp$Ontology, decreasing = FALSE),]
              Reg_GO_terms[[tstbee]] <- temp
              #temp <- Reg_GO_terms[[tstbee]]
              write.csv(temp, file = paste0(dir, "/GO terms - ", tstbee, ".csv"), row.names = FALSE)
              w <- which(vapply(colnames(temp), function(x) { "character" %in% class(temp[[x]]) }, TRUE))
              if (length(w)) {
                for (i in w) { #i <- w[1]
                  w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
                  if (length(w1)) {
                    temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1, ExcelMax-3), "...")
                  }
                }
              }
              require(openxlsx)
              HdrStl <- createStyle(textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE,
                                    numFmt = "TEXT", fontSize = 12)
              wb <- createWorkbook()
              kount <- 0
              for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1]
                w <- which(temp$Ontology == Ontologies[ont])
                if (length(w)) {
                  kount <- kount + 1
                  addWorksheet(wb, ont)
                  writeData(wb, ont, temp[w,])
                  setRowHeights(wb, ont, 1, 60)
                  setColWidths(wb, ont, 1:ncol(temp), 12)
                  setColWidths(wb, ont, which(colnames(temp) == "Term"), 45)
                  setColWidths(wb, ont, which(colnames(temp) %in% gn), 20)
                  setColWidths(wb, ont, which(colnames(temp) %in% pr), 20)
                  setColWidths(wb, ont, which(colnames(temp) %in% pg), 20)
                  addStyle(wb, ont, HdrStl, 1, 1:ncol(temp))
                  addStyle(wb, ont, createStyle(numFmt = "0"), 2:(length(w)+1),
                           which(colnames(temp) %in% kn), gridExpand = TRUE)
                  addStyle(wb, ont, createStyle(numFmt = "0.000"), 2:(length(w)+1),
                           which(colnames(temp) %in% c(zs, lf, pv)), gridExpand = TRUE)
                }
              }
              if (kount) {
                saveWorkbook(wb, paste0(dir, "/GO terms - ", tstbee, ".xlsx"), overwrite = TRUE)
                #openXL(paste0(dir, "/GO terms - ", tstbee, ".xlsx"))
                Kol2 <- grep(paste0("^Significance - .+", " ", max(BH.FDR)*100, "%"),
                             colnames(Reg_GO_terms[[tstbee]]), value = TRUE)
                if (length(Kol2)) {
                  if (length(Kol2) > 1) {
                    w <- which(apply(Reg_GO_terms[[tstbee]][,Kol2], 1, function(x) {"+" %in% x}))
                  } else { w <- which(vapply(Reg_GO_terms[[tstbee]][,Kol2], function(x) { "+" %in% x }, TRUE)) }
                  write.csv(Reg_GO_terms[[tstbee]][w,], file = paste0(dir, "/Regulated GO terms - ", tstbee, ".csv"),
                            row.names = FALSE)
                }
                # Summary table and heatmap of number of regulated GO terms
                Kol3 <- grep(paste0("^Significance - [^ ]+ [1-9][0-9]*\\.*[0-9]*%$"),
                             colnames(GO_Plots[[tstbee]]$GO_terms), value = TRUE)
                tst <- as.numeric(gsub(paste0("^Significance - [^ ]+ |%$"), "", Kol3))
                Kol3 <- Kol3[which(tst == max(tst))]
                N <- length(Kol3)
                if (N > 1) {
                  temp <- as.data.frame(matrix(rep("", (N+1)^2), ncol = N+1))
                  W <- lapply(Kol3, function(x) { which(GO_Plots[[tstbee]]$GO_terms[[x]] == "+") })
                  names(W) <- gsub(paste0("^Significance - | ", max(BH.FDR)*100, "%$"), "", Kol3)
                  temp[2:(N+1), 1] <- temp[1, 2:(N+1)] <- names(W)
                  for (i in 2:(N+1)) { #i <- 2
                    temp[i, 2:(N+1)] <- vapply(Kol3, function(x) {
                      sum((GO_Plots[[tstbee]]$GO_terms[[x]] == "+")&(GO_Plots[[tstbee]]$GO_terms[[Kol3[i]]] == "+"),
                          na.rm = TRUE)
                    }, 1)
                  }
                  names(W) <- cleanNms(gsub(" [0-9]+%$", "", names(W)))
                  tst <- vapply(strsplit(names(W), " - "), length, 1)
                  tst <- (min(tst) > 1)&(length(unique(tst)) == 1)
                  if (tst) {
                    tst <- as.data.frame(t(sapply(strsplit(names(W), " - "), unlist)))
                    l <- apply(tst, 2, function(x) { length(unique(x)) })
                    tst <- tst[, which(l > 1)]
                    names(W) <- do.call(paste, c(tst, sep = " - "))
                  }
                  temp[2:(N+1), 1] <- temp[1, 2:(N+1)] <- names(W)
                  nm <- paste0("N. of co-regulated GO terms\n", tstrt, "\n(", tolower(bee), ")")
                  write.csv(temp, file = paste0(dir, "/", gsub("\n", " - ", gsub("\n\\(", " (", nm)), ".csv"), row.names = FALSE)
                  temp2 <- temp[2:(N+1), 2:(N+1)]
                  colnames(temp2) <- temp[1, 2:(N+1)]
                  rownames(temp2) <-  temp[2:(N+1), 1]
                  for (i in 1:nrow(temp2)) { temp2[[i]] <- as.numeric(temp2[[i]]) }
                  if (max(is.all.good(unlist(temp2))) > 0) {
                    temp2 <- as.matrix(temp2)
                    basic.heatmap(temp2, "N. of co-regulated GO terms", paste0(tstrt, "\n(", tolower(bee), ")"),
                                  save = c("pdf", "jpeg"), folder = dir)
                  }
                } else { warning(paste0(tstrt, " performed for only one condition, skipping.")) }
              }
            }
          } else { warning(paste0("Filter ", tstbee, " has length 0, skipping.")) }
        }
      } else { warning(paste0("No filters available for ", tstrt, ", skipping.")) }
    }
  }
  if (globalGO) {
    msg <- " - Dataset"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    #
    dir <- paste0(wd, "/Reg. analysis/GO enrich/Dataset")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    #
    Mode <- "dataset"
    dataType <- "PG"
    if ((!exists("GO_mappings"))&&(file.exists("GO_mappings.RData"))) {
      loadFun("GO_mappings.RData")
    }
    if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) {
      loadFun("GO_terms.RData")
    }
    #
    #try(rm(list = allArgs), silent = TRUE)
    #
    Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_enrich.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    # Off because it was taking too long and failing sometimes
    # ClueGO really cannot handle too large gene lists, which are the norm here for dataset analysis
    #if (runClueGO) {
    #  clueGO_outDir <- dir
    #  clueGO_type <- "Enrichment/Depletion (Two-sided hypergeometric test)"
    #  Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_enrich.R")
    #  rstudioapi::documentOpen(Src)
    #  source(Src, local = FALSE)
    #}
    #
    # Cleanup - do it now, not within sources!
    try(rm(list = allArgs), silent = TRUE)
    #
    # Quick Fisher exact test on GO terms of interest - looking only at observed data per group
    Src <- paste0(libPath, "/extdata/R scripts/Sources/interestGO_Fisher.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    GO_Plots_2 %<o% goRES
    #
    if ((!is.null(GO_Plots_2))&&("All_GO_terms" %in% names(GO_Plots_2))) {
      GO_terms <- GO_Plots_2$All_GO_terms
      GO_Plots_2$All_GO_terms <- NULL
    }
    if ((create_plotly)&&(!create_plotly_local)) { plot_ly$"GO plots - Observed dataset vs Theoretical proteome" <- GO_Plots_2$GO_plot_ly }
  }
  DatAnalysisTxt <- paste0(DatAnalysisTxt,
                           " GO terms enrichment analysis was performed, comparing for each test regulated against observed protein groups, using topGO",
                           c("", "and ClueGO")[runClueGO+1], ".")
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

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
# It makes sense to close/re-create parallel clusters regularly to reduce memory usage
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Proteomic ruler
if (protrul) {
  ProtRulRoot %<o% "log10(est. copies/cell) - "
  tempPG <- PG
  exprsRt <- paste0("Mean ", prtRfRoot)
  if (LocAnalysis) {
    tempPG <- tempPG[, grep(topattern(exprsRt), colnames(tempPG), invert = TRUE)]
    for (grp2 in SubCellFracAggr2$values) { #grp2 <- SubCellFracAggr2$values[1]
      em2 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == grp2),]
      for (grp in unique(em2[[SubCellFracAggr$column]])) { #grp <- unique(em2[[SubCellFracAggr$column]])[1]
        em <- em2[which(em2[[SubCellFracAggr$column]] == grp),]
        kol <- paste0(prtRfRoot, em$Ref.Sample.Aggregate)
        tempPG[[paste0(prtRfRoot, grp)]] <- apply(10^tempPG[, kol], 1, function(x) {
          log10(sum(is.all.good(x)))
        })
      }
      tempPG[[paste0(exprsRt, grp2)]] <- apply(tempPG[, paste0(prtRfRoot, unique(em2[[SubCellFracAggr$column]]))], 1, function(x) {
        mean(is.all.good(x))
      })
    }
  }
  temp <- try(Prot.Ruler(tempPG, db, exprsRt, NuclL = ProtRulNuclL), silent = TRUE)
  if (class(temp) == "list") {
    db <- temp$Database
    temp <- temp$Protein.groups
    kol <- c(grep(topattern(exprsRt), colnames(temp), value = TRUE), grep(topattern(ProtRulRoot), colnames(temp), value = TRUE))
    PG[, kol] <- temp[, kol]
    protrul <- TRUE
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " Protein group copy numbers per cell were estimated using a variant of the proteome ruler logic, normalizing to scaled values of all identified histones.")
  } else {
    warning("Failed to run Prot.Ruler function; is the remote NCBI server available?")
    protrul <- FALSE
  }
  rm(temp)
}

#### Code chunk - Summary table and QC plots
if ((create_plotly)&&(!create_plotly_local)) { # This code is so old it's auld!
  plot_ly_addresses %<o% data.frame(`Plot type` = NA, `Name` = NA, `Address` = NA)
  kount <- 1
  for (n in seq_along(plot_ly)) { #n <- 1
    i <- plot_ly[[n]]
    if (length(i)) {
      plot_ly_addresses[kount,] <- c(names(plot_ly)[n], "", "")
      for (j in seq_along(i)) { #j <- 1
        kount <- kount + 1
        k <- i[[j]]
        plot_ly_addresses[kount,] <- c("",
                                       names(i)[j],
                                       paste0(gsub("\\.embed$", "", k$embed_url), "/?share_key=", k$share_key))
      }
      kount <- kount + 1
    }
  }
  class(plot_ly_addresses$Address) <- "hyperlink"
  wb <- createWorkbook()
  sheet  <- addWorksheet(wb, "Plotly plot addresses")
  writeData(wb, sheet = "Plotly plot addresses", x = plot_ly_addresses)
  style1 <- createStyle(textDecoration = "bold")
  addStyle(wb, sheet = "Plotly plot addresses", style = style1, rows = 2:nrow(plot_ly_addresses)+1, cols = 1, gridExpand = FALSE, stack = FALSE)
  style2 <- createStyle()
  addStyle(wb, sheet = "Plotly plot addresses", style = style2, rows = 2:nrow(plot_ly_addresses)+1, cols = 2, gridExpand = FALSE, stack = FALSE)
  style3 <- createStyle(textDecoration = "italic")
  addStyle(wb, sheet = "Plotly plot addresses", style = style3, rows = 2:nrow(plot_ly_addresses)+1, cols = 3, gridExpand = FALSE, stack = FALSE)
  setColWidths(wb, sheet = "Plotly plot addresses", cols = 1, widths = 25)
  setColWidths(wb, sheet = "Plotly plot addresses", cols = 2, widths = 50)
  setColWidths(wb, sheet = "Plotly plot addresses", cols = 3, widths = 60)
  freezePane(wb, 1, firstRow = TRUE)
  saveWorkbook(wb, file = "Tables/Plotly plot addresses.xlsx", overwrite = TRUE)
}
mods <- setNames(Modifs$Mark[which(Modifs$Type == "Variable")],
                 nm = Modifs$"Full name"[which(Modifs$Type == "Variable")])
tmp <- aggregate(Frac.map$"Raw file", list(Frac.map$MQ.Exp), length)
tmp <- round(mean(tmp$x)) # Size of a full fraction set, rounding for cases where we removed some fractions
defSc <- 60 # (Non-strict) default maximum number of files to look at per plot
if (tmp > defSc) {
  # If one set of fractions is larger than defSc
  sc <- tmp
} else {
  # What is closest to default: n or n+1 full sets of fractions?
  tst <- (defSc %% tmp) >= defSc/2
  # Identify N = fixed number of files from full fraction sets we can fit in one plot
  N <- c(floor(defSc/tmp), ceiling(defSc/tmp))[tst+1]*tmp
  # If we divide the total number of files by that number, how many plots do we have?
  Nplts <- ceiling(nrow(Frac.map)/N)
  # So that makes how many files per plot:
  sc <- ceiling(nrow(Frac.map)/Nplts)
}
sc <- max(c(sc, 1))
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
Exp_summary$"Biological sample"[1] <- "All samples"
Exp_summary <- Exp_summary[, c("Sample", "Biological sample",
                               colnames(Exp_summary)[which(!colnames(Exp_summary) %in% c("Sample", "Biological sample"))])]
write.csv(Exp_summary, paste0(wd, "/Workflow control/Summary.csv"), row.names = FALSE)
#Exp_summary <- read.csv(paste0(wd, "/Workflow control/Summary.csv"), check.names = FALSE)

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)

#### Code chunk - XML coverage columns
Src <- paste0(libPath, "/extdata/R scripts/Sources/xml_Coverage_columns.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
# To do: also PTMs in a different color (one for all)

#### Code chunk - GO term columns
GO_PG_col %<o% unique(unlist(strsplit(Param$GO.tabs, ";")))
GO_filt %<o% length(GO_PG_col) > 0
if (GO_filt) {
  if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) { loadFun("GO_terms.RData") }
  GO_PG_col <- GO_PG_col[which(GO_PG_col %in% GO_terms$ID)]
  GO_filt <- length(GO_PG_col) > 0
}
if (GO_filt) {
  tmp <- listMelt(strsplit(PG$`GO-ID`, ";"), 1:nrow(PG), c("Term", "Row"))
  Offspring <- setNames(lapply(GO_PG_col, function(x) { #x <- "GO:0009725"
    ont <- Ontology(x)
    x <- c(x, get(paste0("GO", ont, "OFFSPRING"))[[x]])
    x <- x[which(!is.na(x))]
    return(x)
  }), GO_PG_col)
  tmp <- tmp[which(tmp$Term %in% unlist(Offspring)),]
  tmp <- aggregate(tmp$Row, list(tmp$Term), c)
  colnames(tmp) <- c("Term", "Rows")
  w <- which(vapply(GO_PG_col, function(x) { sum(Offspring[[x]] %in% tmp$Term) }, 1) == 0)
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
    for (go in GO_PG_col) { #go <- GO_PG_col[1]
      w <- which(tmp$Term %in% Offspring[[go]])
      w2 <- unique(unlist(tmp$Rows[w]))
      PG[w2, GO_PG_col2[go]] <- "+"
    }
    #View(PG[, GO_PG_col2])
    tst <- setNames(vapply(GO_PG_col2, function(x) {
      sum(PG[[x]] == "+")
    }, 1), GO_PG_col2)
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
w <- which(vapply(colnames(ev), function(x) { "list" %in% class(ev[[x]]) }, TRUE))
if (length(w)) { for (i in w) { ev[[i]] <- parSapply(parClust, ev[[i]], paste, collapse = ";") } }
data.table::fwrite(ev, paste0(dir, "/evidence.tsv"), sep = "\t", row.names = FALSE, na = "NA")
#
## Main peptidoforms- and protein groups-level, multi-tabs report
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
for (nm in names(Prot.Expr.Root)) { #nm <- names(pep.ref[1])
  rpl <- intNms(nm)
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
for (nm in unique(c(names(pep.ratios.ref), names(Prot.Rat.Root)))) { #nm <- names(pep.ratios.ref[1])
  rpl <- ratNms(nm)
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Ratios"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Ratios"
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
  #KolNames <- xlTabs[[sheetnm]]
  klnms <- KolNames
  KolNames <- gsub("Peptides?", "Pep.", KolNames)
  KolNames <- gsub("Evidences?", "PSMs", KolNames)
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
QualFilt %<o% c("In list", "Potential contaminant", "Only identified by site",
                grep("^Quality filter: ", colnames(PG), value = TRUE),
                grep("^Quantity Quality$", colnames(pep), value = TRUE))
if ((DiscFilt)&&(DiscFiltMode == DiscFiltModes[3])) { QualFilt <- c(QualFilt, DiscFiltCols) }
II <- setNames(1, "All peptidoforms")
if ((exists("PTMs_pep"))&&(length(PTMs_pep))) {
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
    CoreCol2 %<o% c("Leading proteins", "Leading razor proteins",
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
    kol <- c(CoreCol, "In list", CoreCol2, evcol, spcol, "PEP", quantcol, signcol, regcol, qualFlt[which(qualFlt != "In list")], aacol)
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
      #mnratcolF %<o% grep("Mean log2\\(Ratio\\) - ", colnames(tempPepF), value = TRUE)
      #m <- match(mnratcolF, colnames(tempPepF))
      #colnames(tempPepF)[m] <- paste0("mod. F-test ", mnratcolF)
      #mnratcolF <- paste0("mod. F-test ", mnratcolF)
      pvalcolF %<o% F_Root
      signcolF %<o% grep("^mod\\. F-test Significant", colnames(tempPepF), value = TRUE)
      regcolF %<o% grep("^mod\\. F-test Regulated", colnames(tempPepF), value = TRUE)
      Fkol %<o% c(regcolF, #mnratcolF,
                  pvalcolF, signcolF)
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
        #ColumnsTbl[["F-test summary Ratios"]] <- mnratcolF
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
    #aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), unique)
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
    #saveFun(WorkBook, file = "WorkBook_bckp.RData")
    #wb_save(WorkBook, paste0(wd, "/tst.xlsx")); xl_open(paste0(wd, "/tst.xlsx"))
    #loadFun("WorkBook_bckp.RData")
  }
}
TbNm <- "Protein groups"
tblMode <- tblMode2 <- "PG"
# Function for editing the header
KolEdit <- function(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl, locTbl = SubCellLocTbl) {
  #KolNames <- xlTabs[[sheetnm]]
  klnms <- KolNames
  KolNames <- gsub("Peptides?", "Pep.", KolNames)
  KolNames <- gsub("Evidences?", "PSMs", KolNames)
  KolNames <- gsub("Spectr((al)|(um))", "Spec.", KolNames)
  KolNames <- gsub("Razor", "Raz.", KolNames)
  KolNames <- gsub("Unique", "Uniq.", KolNames)
  KolNames <- gsub("MS/MS", "MS2", KolNames)
  for (nm in names(intTbl)) { #nm <- names(intTbl)[1] #nm <- names(intTbl)[2]
    m <- match(intTbl[[nm]]$Log, KolNames)
    w <- which(!is.na(m))
    if (length(w)) {
      rpl <- intNms(nm)
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
  if (!missing("locTbl")) {
    for (nm in names(locTbl)) {
      m <- match(locTbl[[nm]]$Columns, KolNames)
      w <- which(!is.na(m))
      if (length(w)) {
        KolNames[m[w]] <- paste0(nm, locTbl[[nm]]$Sample[w])
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
  KolNames <- gsub(".*Pvalue\\)( - )?", "-log10 pval. ", KolNames)
  KolNames <- gsub(".*Significant-", "signif. ", KolNames)
  KolNames <- gsub(".*Regulated - ", "reg. ", KolNames)
  KolNames <- gsub("log10\\(est\\. copies/cell\\) (- )?", "ProtRul. ", KolNames)
  KolNames <- gsub(paste0(topattern("Sequence coverage [%] ", start = FALSE), "(- )?"), "Cov. ", KolNames)
  KolNames[which(KolNames == "Max. theoretical sequence coverage [%]")] <- "Theoretical max."
  KolNames[which(KolNames == "Sequence coverage [%]")] <- "All peptides"
  KolNames[which(KolNames == "Uniq. + razor sequence coverage [%]")] <- "Unique + razor"
  KolNames[which(KolNames == "Uniq. sequence coverage [%]")] <- "Unique"
  KolNames <- gsub("^Cluster \\([^\\)]+\\) - ", "Clust. ", KolNames)
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
  return(KolNames) #View(KolNames)
}
#
tempData <- get(tblMode)
CoreCol <- c("Leading protein IDs", "Common Names", "Genes")
CoreCol2 <- c("Protein IDs", "Names", "id", "Mol. weight [kDa]")
pepevspeccol <- c("Peptides count",
                  grep("^Peptides count - ", colnames(tempData), value = TRUE),
                  "Peptide IDs", "Razor peptide IDs", "Unique peptide IDs",
                  grep("^Peptide IDs - ", colnames(tempData), value = TRUE),
                  "psmS count",
                  grep("^Evidences count - ", colnames(tempData), value = TRUE),
                  "Evidence IDs",
                  grep("^Evidence IDs - ", colnames(tempData), value = TRUE),
                  "Spectral count",
                  grep("^Spectral count - ", colnames(tempData), value = TRUE),
                  "Spectrum IDs",
                  grep("^Spectrum IDs - ", colnames(tempData), value = TRUE),
                  "Biot. peptides count",
                  grep("^Biot\\. peptides count - ", colnames(tempData), value = TRUE),
                  "Biot. peptide IDs",
                  grep("^Biot\\. peptides IDs - ", colnames(tempData), value = TRUE),
                  "Biot. peptides [%]",
                  "Biot. evidences count",
                  grep("^Biot\\. evidences count - ", colnames(tempData), value = TRUE),
                  "Biot. evidence IDs",
                  grep("^Biot\\. evidence IDs - ", colnames(tempData), value = TRUE),
                  "Biot. spectral count",
                  grep("^Biot\\. spectral count - ", colnames(tempData), value = TRUE),
                  "Biot. spectrum IDs",
                  grep("^Biot\\. spectrum IDs - ", colnames(tempData), value = TRUE))
pepevspeccol <- pepevspeccol[which(pepevspeccol %in% colnames(tempData))]
pepcountcol1 <- grep("[Pp]eptides count$", pepevspeccol, value = TRUE)
pepcountcol2 <- grep("[Pp]eptides count - ", pepevspeccol, value = TRUE)
pepidcol1 <- grep("[Pp]eptide IDs$", pepevspeccol, value = TRUE)
pepidcol2 <- grep("[Pp]eptide IDs - ", pepevspeccol, value = TRUE)
evcountcol1 <- grep("[Ev]vidences count$", pepevspeccol, value = TRUE)
evcountcol2 <- grep("[Ev]vidences count - ", pepevspeccol, value = TRUE)
evidcol1 <- grep("[Ev]vidence IDs$", pepevspeccol, value = TRUE)
evidcol2 <- grep("[Ev]vidence IDs - ", pepevspeccol, value = TRUE)
speccountcol1 <- grep("[Ss]pectral count$", pepevspeccol, value = TRUE)
speccountcol2 <- grep("[Ss]pectral count - ", pepevspeccol, value = TRUE)
specidcol1 <- grep("[Ss]pectrum IDs$", pepevspeccol, value = TRUE)
specidcol2 <- grep("[Ss]pectrum IDs - ", pepevspeccol, value = TRUE)
pepevspeccol <- c(pepcountcol1, pepcountcol2, pepidcol1, pepidcol2,
                  evcountcol1, evcountcol2, evidcol1, evidcol2,
                  speccountcol1, speccountcol2, specidcol1, specidcol2)
pepevspeccola <- grep("((([Ss]pectral|[Pp]eptides|[Ee]vidences) count)|(([Ss]pectrum|[Pp]eptide|[Ee]vidence) IDs))$", pepevspeccol, value = TRUE)
pepevspeccolb <- grep("((([Ss]pectral|[Pp]eptides|[Ee]vidences) count)|(([Ss]pectrum|[Pp]eptide|[Ee]vidence) IDs)) - ", pepevspeccol, value = TRUE)
kol <- c(CoreCol, "In list", CoreCol2, pepevspeccol)
intRf <- Prot.Expr.Root
names(intRf) <- paste0(names(intRf), " int.")
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
gel <- unlist(intCols)
for (gl in gel) {
  w <- which(is.infinite(tempData[[gl]]))
  tempData[w, gl] <- NA
}
quantcol <- gel
ratRf <- Prot.Rat.Root
names(ratRf) <- paste0(names(ratRf), " rat.")
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
quantcol <- c(quantcol, grl)
quantCols[names(ratRf)] <- ratCols
if (protrul) {
  tmp <- grep(topattern("log10(est. copies/cell) - ", start = FALSE), colnames(tempData), value = TRUE)
  quantcol <- c(quantcol, tmp)
  quantCols[["Proteomic ruler"]] <- tmp
}
pvalcol <- grep(topattern(pvalue.col[pvalue.use]), colnames(tempData), value = TRUE)
regcol <- grep("^((Enriched)|(Regulated)) - ", colnames(tempData), value = TRUE)
signcol <- grep("^Significant-FDR=[1-9][0-9]*\\.*[0-9]*% - ", colnames(tempData), value = TRUE)
signcol <- grep(" - Analysis_[0-9]+", signcol, invert = TRUE, value = TRUE)
covcol <- c(xmlCovCol,
            c("Sequence coverage [%]",
              "Unique + razor sequence coverage [%]",
              "Unique sequence coverage [%]")[1:c(1, 3)[isEukaLike+1]],
            grep(topattern("Sequence coverage [%] - "), colnames(tempData), value = TRUE)) # The complicated way, but ensures the order is correct
kol <- c(kol, "Mol. weight [kDa]", covcol, "PEP", covcol, quantcol, pvalcol, regcol, signcol)
if ((exists("KlustKols"))&&(length(KlustKols))) { kol <- c(kol, KlustKols) }
qualFlt <- QualFilt
if (length(GO_PG_col)) { qualFlt <- c(qualFlt, GO_PG_col2) }
kol <- c(kol, qualFlt[which(qualFlt != "In list")])
if (Annotate) { kol <- c(kol, annot.col) }
if (Annotate&&LocAnalysis) {
  PG$Marker <- ""
  w <- which(PG$Label %in% names(SubCellMark2))
  PG$Marker[w] <- SubCellMark2[match(PG$Label[w], names(SubCellMark2))]
  SubCellLocTbl <- list()
  lokol1 <- paste0("Localisation - ", SubCellFracAggr$values)
  w <- which(lokol1 %in% colnames(PG))
  lokol1 <- lokol1[w]
  SubCellLocTbl$Localisation <- list(Columns = lokol1,
                                     Sample = SubCellFracAggr$values[w])
  kol <- c(kol, "Marker", lokol1)
  if (LocAnalysis2) {
    rt2 <- SSD.Root
    lokol2 <- grep(topattern(rt2), colnames(PG), value = TRUE)
    rt3 <- paste0("Mean ", SSD.Root)
    lokol3 <- grep(topattern(rt3), colnames(PG), value = TRUE)
    SubCellLocTbl$SSDs <- list(Columns = c(lokol3, lokol2),
                               Sample = c(gsub(topattern(rt3), "", lokol3),
                                          gsub(topattern(rt2), "", lokol2)))
    rt <- SSD.Pval.Root
    lokol4 <- grep(topattern(rt), colnames(PG), value = TRUE)
    SubCellLocTbl$"SSD P-vals." <- list(Columns = lokol4,
                                        Sample = gsub(topattern(rt), "", lokol4))
    rt <- "Signif. SSDs-FDR="
    lokol5 <- grep(topattern(rt), colnames(PG), value = TRUE)
    SubCellLocTbl$"Signif. SSDs" <- list(Columns = lokol5,
                                         Sample = gsub(topattern(rt), "", lokol5))
    rt <- "Re-localized"
    lokol6 <- grep(topattern(rt), colnames(PG), value = TRUE)
    SubCellLocTbl$"Re-Loc." <- list(Columns = lokol6,
                                    Sample = gsub(topattern(rt), "", lokol6))
    kol <- c(kol, lokol6, lokol2, lokol3, lokol4, lokol5)
  }
} else { SubCellLocTbl <- NULL}
kol <- unique(kol[which(kol %in% colnames(tempData))])
tempData <- tempData[, kol]
#
if (F.test) {
  tmpPGf <- F_test_data
  tmpPGf <- tmpPGf[, which(!colnames(tmpPGf) %in% c(Param$Plot.labels, "Rel. log10(Peptides count)", "Av. log10 abundance"))]
  colnames(tmpPGf) <- cleanNms(colnames(tmpPGf), start = FALSE)
  #mnratcolF %<o% grep("Mean log2\\(Ratio\\) - ", colnames(tmpPGf), value = TRUE)
  #m <- match(mnratcolF, colnames(tmpPGf))
  #colnames(tmpPGf)[m] <- paste0("mod. F-test ", mnratcolF)
  #mnratcolF <- paste0("mod. F-test ", mnratcolF)
  pvalcolF %<o% F_Root
  signcolF %<o% grep("^mod\\. F-test Significant", colnames(tmpPGf), value = TRUE)
  regcolF %<o% grep("^mod\\. F-test Regulated", colnames(tmpPGf), value = TRUE)
  for (k in regcolF) {
    tmpPGf[which(tmpPGf[[k]] == "non significant"), k] <- "n.s."
    tmpPGf[which(tmpPGf[[k]] == ""), k] <- "n.t."
  }
  Fkol %<o% c(regcolF, #mnratcolF,
              pvalcolF, signcolF)
  tempData[, Fkol] <- tmpPGf[match(tempData$`Protein IDs`, tmpPGf$`Protein IDs`), Fkol]
}
if (Annotate&&LocAnalysis) {
  if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) { loadFun("GO_terms.RData") }
  GOCC <- GO_terms$ID[which(GO_terms$Ontology == "CC")]
  tempData$"GO-ID (CC)" <- lapply(strsplit(tempData$`GO-ID`, ";"), function(x) { x[which(x %in% GOCC)] })
  w <- which(vapply(tempData$"GO-ID (CC)", length, 1) > 0)
  tempData$"GO (CC)" <- ""
  tempData$"GO (CC)"[w] <- vapply(tempData$"GO-ID (CC)"[w], function(x) {
    paste(GO_terms$Term[match(x, GO_terms$ID)], collapse = ";")
  }, "")
  if (LocAnalysis2) {
    for (k in lokol6) {
      tempData[which(tempData[[k]] == "non significant"), k] <- "n.s."
      tempData[which(tempData[[k]] == ""), k] <- "n.t."
      tempData[[k]] <- gsub("^up, FDR = ", "re-loc., FDR = ", tempData[[k]])
    }
  }
}
#
# Which columns are affected by each style
# - IDs
ColumnsTbl <- list(IDs = c(CoreCol, CoreCol2))
# - Peptide and evidence counts and IDs
ColumnsTbl[["Global Pep. IDs"]] <- pepidcol1
ColumnsTbl[["Global Pep. counts"]] <- pepcountcol1
ColumnsTbl[["Pep. IDs"]] <- pepidcol2
ColumnsTbl[["Pep. counts"]] <- pepcountcol2
ColumnsTbl[["Global Ev. IDs"]] <- evidcol1
ColumnsTbl[["Global Ev. counts"]] <- evcountcol1
ColumnsTbl[["Ev. IDs"]] <- evidcol2
ColumnsTbl[["Ev. counts"]] <- evcountcol2
ColumnsTbl[["Global Spec. IDs"]] <- specidcol1
ColumnsTbl[["Global Spec. counts"]] <- speccountcol1
ColumnsTbl[["Spec. IDs"]] <- specidcol2
ColumnsTbl[["Spec. counts"]] <- speccountcol2
if (IsBioID2) {
  # Some of these may currently be empty - we don't want to overload the files with columns
  ColumnsTbl[["Global Biot. Pep. IDs"]] <- biotpepidcol1
  ColumnsTbl[["Global Biot. Pep. counts"]] <- biotpepcountcol1
  ColumnsTbl[["Biot. Pep. IDs"]] <- biotpepidcol2
  ColumnsTbl[["Biot. Pep. counts"]] <- biotpepcountcol2
  ColumnsTbl[["Global Biot. Ev. IDs"]] <- biotevidcol1
  ColumnsTbl[["Global Biot. Ev. counts"]] <- biotevcountcol1
  ColumnsTbl[["Biot. Ev. IDs"]] <- biotevidcol2
  ColumnsTbl[["Biot. Ev. counts"]] <- biotevcountcol2
  ColumnsTbl[["Biot. Pep. %"]] <- "Biot. peptides [%]"
}
# Quantitation
# - Expression values
for (nm in names(intRf)) { #nm <- names(intRf[1])
  rpl <- intNms(nm)
  ColumnsTbl[[paste0(rpl, ", avg.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Average")]
  ColumnsTbl[[paste0(rpl, ", indiv.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Individual")]
}
# - Ratios
for (nm in names(ratRf)) { #nm <- names(ratRf[1])
  rpl <- ratNms(nm)
  ColumnsTbl[[paste0(rpl, ", avg.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Average")]
  ColumnsTbl[[paste0(rpl, ", indiv.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Individual")]
}
# - P-values
ColumnsTbl[["P-values"]] <- pvalcol
# - Significant
ColumnsTbl[["Significant"]] <- signcol
# - Regulated
ColumnsTbl[["Regulated"]] <- regcol
# - Proteome Ruler
if (protrul) {
  ColumnsTbl[["Proteome Ruler"]] <- grep(topattern("log10(est. copies/cell) - ", start = FALSE), colnames(tempData), value = TRUE)
}
#
# F-test
if (F.test) {
  #ColumnsTbl[["F-test summary Ratios"]] <- mnratcolF
  ColumnsTbl[["F-test P-values"]] <- pvalcolF
  ColumnsTbl[["F-test significant"]] <- signcolF
  ColumnsTbl[["F-test regulated"]] <- regcolF
}
# - Annotations
if (Annotate) {
  annot.col2 <- gsub("_names$", " names", annot.col)
  AnnotTbl$Columns <- list(c("GO", "GO-ID"), c("Taxonomy", "TaxID"), NA, NA, NA, NA, "EMBL", NA)
  annot <- c("InterPro", "Pfam", "PIRSF", "PROSITE")
  for (i in annot) { AnnotTbl$Columns[match(i, AnnotTbl$Name)] <- list(c(i, paste0(i, " names"))) }
  AnnotTbl$Columns[match("Other", AnnotTbl$Name)] <- list(annot.col2[which(!annot.col2 %in% unlist(AnnotTbl$Columns))])
  for (i in 1:nrow(AnnotTbl)) { ColumnsTbl[[paste0(AnnotTbl$Name[i], " annotations")]] <- AnnotTbl$Columns[[i]] }
  if (LocAnalysis) {
    ColumnsTbl[["GO annotations"]] <- c(ColumnsTbl[["GO"]], "GO (CC)")
    ColumnsTbl[["GO-ID annotations"]] <- c(ColumnsTbl[["GO-ID"]], "GO-ID (CC)")
    ColumnsTbl[["Localisation"]] <- lokol1
    ColumnsTbl[["Marker"]] <- "Marker"
    if (LocAnalysis2) {
      ColumnsTbl[["SSDs"]] <- lokol2
      ColumnsTbl[["Mean SSDs"]] <- lokol3
      ColumnsTbl[["SSDs P-values"]] <- lokol4
      ColumnsTbl[["SSDs significant"]] <- lokol5
      ColumnsTbl[["Re-localized"]] <- lokol6
    }
  }
}
# - PEP
ColumnsTbl[["PEP"]] <- "PEP"
# - Filters
ColumnsTbl[["Filters"]] <- qualFlt
# - Clusters
if ((exists("KlustKols"))&&(length(KlustKols))) { ColumnsTbl[["Cluster"]] <- KlustKols }
# - Coverage
ColumnsTbl[["Coverage"]] <- covcol
# Melt
ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }, 1) > 0)]
ColumnsTbl <- listMelt(ColumnsTbl, names(ColumnsTbl), c("Col", "Grp"))
stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
#tst <- aggregate(ColumnsTbl$Col, list(ColumnsTbl$Col), length)
#tst[which(tst$x > 1),]
ColumnsTbl$Class <- ""
ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Protein Group information"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(pepcountcol1, pepcountcol2))] <- "Peptides count"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(pepidcol1, pepidcol2))] <- "Peptide IDs"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evcountcol1, evcountcol2))] <- "Evidences count"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evidcol1, evidcol2))] <- "Evidence IDs"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(speccountcol1, speccountcol2))] <- "Spectral count"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(specidcol1, specidcol2))] <- "Spectrum IDs"
if (IsBioID2) {
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(biotpepcountcol1, biotpepcountcol2))] <- "Biotin peptides count"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(biotpepidcol1, biotpepidcol2))] <- "Biotin peptide IDs"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(biotevcountcol1, biotevcountcol2))] <- "Biotin evidences count"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(biotevidcol1, biotevidcol2))] <- "Biotin evidence IDs"
}
for (nm in names(intRf)) { #nm <- names(intRf)[1]
  rpl <- intNms(nm, TRUE)
  kl <- c(paste0("Mean ", intRf[[nm]], VPAL$values), paste0(intRf[[nm]], RSA$values))
  kl <- kl[which(kl %in% ColumnsTbl$Col)]
  ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
}
ColumnsTbl$Class[grep("Ruler", ColumnsTbl$Grp)] <- "log10(est. copies/cell)"
for (nm in names(ratRf)) { #nm <- names(ratRf)[1]
  rpl <- ratNms(nm, TRUE)
  kl <- c(paste0("Mean ", ratRf[[nm]], VPAL$values), paste0(ratRf[[nm]], RSA$values))
  kl <- kl[which(kl %in% ColumnsTbl$Col)]
  ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
}
ColumnsTbl$Class[which(ColumnsTbl$Grp == "Proteome Ruler")] <- "log10(est. copies/cell)"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% pvalcol)] <- "P-value"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% signcol)] <- "Significant"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
if (F.test) {
  #ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test summary Ratios")] <- "F-test"
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test P-values")] <- "F-test"
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test significant")] <- "F-test"
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test regulated")] <- "F-test"
}
if (Annotate) {
  ColumnsTbl$Class[grep("[Aa]nnotations", ColumnsTbl$Grp)] <- "Annotations"
  if (LocAnalysis) {
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "Localisation")] <- "Localisation"
    if (LocAnalysis2) {
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "SSDs P-values")] <- gsub(" - $", "", SSD.Pval.Root)
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "SSDs significant")] <- paste0(gsub(" - $", "", SSD.Pval.Root))
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "Re-localized")] <- "Re-localized"
    }
  }
}
ColumnsTbl$Class[grep("[Ss]equence coverage \\[%\\]", ColumnsTbl$Col)] <- "Sequence coverage [%]"
ColumnsTbl$Class[grep("^1st ID cov\\.", ColumnsTbl$Col)] <- "1st accession sequence coverage (peptides)"
if ((exists("KlustKols"))&&(length(KlustKols))) {
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "Cluster")] <- paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ")")
}
ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters", "Negative filter"))] <- "QC filters"
stopifnot(min(nchar(ColumnsTbl$Class)) > 0)
#View(ColumnsTbl[which(nchar(ColumnsTbl$Class) == 0),])
w <- c(which(ColumnsTbl$Class == "General Protein Group information"),
       as.integer(unlist(lapply(names(intCols), function(nm) {
         rpl <- intNms(nm, TRUE)
         which(ColumnsTbl$Class == rpl)
       }))),
       which(ColumnsTbl$Class == "log10(est. copies/cell)"),
       as.integer(unlist(lapply(names(ratCols), function(nm) {
         rpl <- ratNms(nm, TRUE)
         which(ColumnsTbl$Class == rpl)
       }))),
       which(ColumnsTbl$Class == "P-value"),
       which(ColumnsTbl$Class == "Significant"),
       which(ColumnsTbl$Class == "Regulated"),
       #which(ColumnsTbl$Class == "F-test: Mean log2(Ratio)"),
       #which(ColumnsTbl$Class == "F-test: -log10(P-values)"),
       #which(ColumnsTbl$Class == "F-test: Significant"),
       #which(ColumnsTbl$Class == "F-test: Regulated"),
       which(ColumnsTbl$Class == "F-test"))
if (Annotate&&LocAnalysis) {
  w <- c(w,
         which(ColumnsTbl$Class == gsub(" - $", "", SSD.Pval.Root)),
         which(ColumnsTbl$Class == paste0(gsub(" - $", "", SSD.Pval.Root))),
         which(ColumnsTbl$Class == "Re-localized"))
}
w <- c(w,
       which(ColumnsTbl$Class == "QC filters"),
       which(ColumnsTbl$Class == "Sequence coverage [%]"),
       which(ColumnsTbl$Class == "1st accession sequence coverage (peptides)"),
       which(ColumnsTbl$Class == paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ")")),
       which(ColumnsTbl$Class == "Peptides count"),
       which(ColumnsTbl$Class == "Peptide IDs"),
       which(ColumnsTbl$Class == "Evidences count"),
       which(ColumnsTbl$Class == "Evidence IDs"),
       which(ColumnsTbl$Class == "Spectral count"),
       which(ColumnsTbl$Class == "Spectrum IDs"),
       which(ColumnsTbl$Class == "Biotin peptides count"),
       which(ColumnsTbl$Class == "Biotin peptide IDs"),
       which(ColumnsTbl$Class == "Biotin evidences count"),
       which(ColumnsTbl$Class == "Biotin evidence IDs"),
       which(ColumnsTbl$Class == "Annotations"))
stopifnot(length(w) == nrow(ColumnsTbl))
#View(ColumnsTbl[which(!1:nrow(ColumnsTbl) %in% w),])
ColumnsTbl <- ColumnsTbl[w,]
ColumnsTbl$Hide <- ColumnsTbl$Class %in% c("Peptide IDs", "Peptides count", "Evidence IDs", "Evidences count", "Spectral count", "Spectrum IDs",
                                           "Biotin peptides count", "Biotin peptide IDs", "Biotin evidences count", "Biotin evidence IDs",
                                           "Annotations",
                                           unique(grep("^F-test: ", ColumnsTbl$Class, value = TRUE)) # Let's not show too much stuff
                                           )
if (length(intCols) > 1) {
  for (nm in names(intCols)[1:(length(intCols) - 1)]) {
    rpl <- intNms(nm)
    ColumnsTbl$Hide[which(ColumnsTbl$Class == rpl)] <- TRUE
  }
}
if (length(ratCols) > 1) {
  for (nm in names(ratCols)[1:(length(ratCols) - 1)]) {
    rpl <- ratNms(nm)
    ColumnsTbl$Hide[which(ColumnsTbl$Class == rpl)] <- TRUE
  }
}
#
if (MakeRatios) { a <- KolEdit(ColumnsTbl$Col, intColsTbl, ratColsTbl) } else { a <- KolEdit(ColumnsTbl$Col, intColsTbl) }
ColumnsTbl$edit_Col <- unlist(a)
#wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#saveFun(WorkBook, file = "WorkBook_bckp.RData")
#wb_save(WorkBook, paste0(wd, "/tst.xlsx")); xl_open(paste0(wd, "/tst.xlsx"))
#loadFun("WorkBook_bckp.RData")
if (saintExprs) {
  tblMode <- tblMode2 <- TbNm <- "SAINTexpress"
  # Function for editing the header
  KolEdit <- function(KolNames, #intTbl = intColsTbl,
                      ratTbl = ratColsTbl) {
    #KolNames <- xlTabs[[sheetnm]]
    klnms <- KolNames
    KolNames <- gsub("^log2\\(FC\\) - ", "log2(FC) ", KolNames)
    KolNames <- gsub("^AvgP - ", "AvgP ", KolNames)
    KolNames <- gsub("^MaxP - ", "MaxP ", KolNames)
    KolNames <- gsub("^TopoAvgP - ", "TopoAvgP ", KolNames)
    KolNames <- gsub("^TopoMaxP - ", "TopoMaxP ", KolNames)
    KolNames <- gsub("^SaintScore - ", "SaintScore ", KolNames)
    KolNames <- gsub("^OddsScore - ", "OddsScore ", KolNames)
    KolNames <- gsub("^boosted_by - ", "Boosted by... ", KolNames)
    KolNames <- gsub("^BFDR - ", "BFDR ", KolNames)
    KolNames <- cleanNms(KolNames, start = FALSE)
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
    return(KolNames) #View(KolNames)
  }
  #
  tempData <- allSAINTs
  CoreCol <- c("Protein", "Gene", "Common Name")
  quantCols <- quantcol <- "Av. log10 abundance"
  ratRf <- setNames("log2(FC) - ", "log2(rat.), avg.")
  ratColsTbl <- setNames(lapply(names(ratRf), function(nm) {
    res <- data.frame(Log = paste0(ratRf[nm], VPAL$values),
                      Type = "Average",
                      Sample = VPAL$values)
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
  quantcol <- c(quantcol, grl)
  quantCols[names(ratRf)] <- ratCols
  fdrcol <- grep("^BFDR - ", colnames(tempData), value = TRUE)
  avgPcol <- grep("^AvgP - ", colnames(tempData), value = TRUE)
  maxPcol <- grep("^MaxP - ", colnames(tempData), value = TRUE)
  topoAvgPcol <- grep("^TopoAvgP - ", colnames(tempData), value = TRUE)
  topoMaxPcol <- grep("^TopoMaxP - ", colnames(tempData), value = TRUE)
  saintcol <- grep("^SaintScore - ", colnames(tempData), value = TRUE)
  oddscol <- grep("^OddsScore - ", colnames(tempData), value = TRUE)
  boostcol <- grep("^boosted_by - ", colnames(tempData), value = TRUE)
  qualFlt <- "Potential contaminant"
  kol <- c(CoreCol, "In list", quantcol, fdrcol, avgPcol, maxPcol, topoAvgPcol, topoMaxPcol,
           #saintcol, oddscol, # Those 2 columns look a bit useless... 
           boostcol, qualFlt[which(qualFlt != "In list")])
  kol <- unique(kol)
  kol <- kol[which(kol %in% colnames(tempData))]
  tempData <- tempData[, kol]
  #
  # Re-order
  for (k in rev(fdrcol)) {
    tempData <- tempData[order(tempData[[k]], decreasing = FALSE),]
  }
  m <- match(prot.list, tempData$Protein)
  m <- m[which(!is.na(m))]
  w <- c(m, which(!tempData$Protein %in% prot.list))
  tempData <- tempData[w,]
  #
  # Which columns are affected by each style
  # - IDs
  ColumnsTbl <- list(IDs = CoreCol)
  # Quantitation
  # - Ratios
  for (nm in names(ratRf)) { #nm <- names(ratRf[1])
    ColumnsTbl[[names(ratRf)]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Average")]
  }
  # SAINTexpress stats
  ColumnsTbl$"P-values" <- fdrcol
  ColumnsTbl$AvgP <- avgPcol
  ColumnsTbl$MaxP <- maxPcol
  ColumnsTbl$topoAvgP <- topoAvgPcol
  ColumnsTbl$topoMaxP <- topoMaxPcol
  ColumnsTbl$SaintScore <- saintcol
  ColumnsTbl$OddsScore <- oddscol
  ColumnsTbl$boosted_by <- boostcol
  # - Filters
  ColumnsTbl$Filters <- qualFlt
  # Melt
  ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }, 1) > 0)]
  ColumnsTbl <- listMelt(ColumnsTbl, names(ColumnsTbl), c("Col", "Grp"))
  stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
  #tst <- aggregate(ColumnsTbl$Col, list(ColumnsTbl$Col), length)
  #tst[which(tst$x > 1),]
  ColumnsTbl$Class <- ""
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Protein Group information"
  for (nm in names(ratRf)) { #nm <- names(ratRf)[1]
    rpl <- "log2(rat.), avg."
    kl <- paste0(ratRf[[nm]], VPAL$values)
    kl <- kl[which(kl %in% ColumnsTbl$Col)]
    ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
  }
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% fdrcol)] <- "P-values"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(avgPcol, maxPcol, topoAvgPcol, topoMaxPcol))] <- "SAINTexpress P"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% saintcol)] <- "SAINTexpress score"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% oddscol)] <- "SAINTexpress odds score"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% boostcol)] <- "SAINTexpress boost"
  ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters", "Negative filter"))] <- "QC filters"
  stopifnot(min(nchar(ColumnsTbl$Class)) > 0)
  #View(ColumnsTbl[which(nchar(ColumnsTbl$Class) == 0),])
  w <- c(which(ColumnsTbl$Class == "General Protein Group information"),
         as.integer(unlist(lapply(names(ratCols), function(nm) {
           which(ColumnsTbl$Class == "log2(rat.), avg.")
         }))),
         which(ColumnsTbl$Grp == "P-values"),
         which(ColumnsTbl$Grp == "AvgP"),
         which(ColumnsTbl$Grp == "MaxP"),
         which(ColumnsTbl$Grp == "topoAvgP"),
         which(ColumnsTbl$Grp == "topoMaxP"),
         which(ColumnsTbl$Grp == "SaintScore"),
         which(ColumnsTbl$Grp == "OddsScore"),
         which(ColumnsTbl$Grp == "boosted_by"),
         which(ColumnsTbl$Class == "QC filters"))
  stopifnot(length(w) == nrow(ColumnsTbl))
  #View(ColumnsTbl[which(!1:nrow(ColumnsTbl) %in% w),])
  ColumnsTbl <- ColumnsTbl[w,]
  ColumnsTbl$Hide <- ColumnsTbl$Class %in% "SAINTexpress P" # Let's not show too much stuff
  #
  a <- KolEdit(ColumnsTbl$Col, ratColsTbl)
  ColumnsTbl$edit_Col <- unlist(a)
  Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
  #saveFun(WorkBook, file = "WorkBook_bckp.RData")
  #wb_save(WorkBook, paste0(wd, "/tst.xlsx")); xl_open(paste0(wd, "/tst.xlsx"))
  #loadFun("WorkBook_bckp.RData")
}
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/Write_Excel_end_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#WorkBook$get_active_sheet()
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
  tmp$Razor <- lapply(strsplit(tmp$Razor, ";"), function(x) {
    as.logical(toupper(x))
  })
  tmp$IDs <- lapply(strsplit(tmp$IDs, ";"), as.numeric)
  tmp$RazorIDs <- apply(tmp[, c("IDs", "Razor")], 1, function(x) {
    x[[1]][which(x[[2]])]
  })
  for (i in Exp.map$Ref.Sample.Aggregate[which(as.logical(Exp.map$Use))]) { #i <- Exp.map$Ref.Sample.Aggregate[which(as.logical(Exp.map$Use))][1]
    i2 <- cleanNms(i, rep = ".")
    kol <- paste0("LFQIntensity_", i2)
    AmicTbl[[kol]] <- PG[[paste0(prtRfRoot, i)]]/log10(2)
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2)), kol] <- NaN
    kol <- paste0("razorUniqueCount_", i2)
    tmp$Tmp <- strsplit(PG[[paste0("Peptide IDs - ", i)]], ";")
    AmicTbl[[kol]] <- apply(tmp[, c("IDs", "Tmp")], 1, function(x) {
      sum(x[[2]] %in% x[[1]])
    })
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2)), kol] <- 0
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
  for (g in grps) { #g <- grps[1]
    gEd <- cleanNms(g, rep = ".")
    m <- Exp.map[which(Exp.map[[VPAL$column]] == g),]
    ratgrps <- unique(m[[RG$column]])
    g0 <- unique(Exp.map[which((Exp.map[[RG$column]] == ratgrps)&(Exp.map$Reference)),
                         VPAL$column])
    gEd0 <- cleanNms(g0, rep = ".")
    kol <- paste0("P.Value_", gEd, "__vs__", paste(gEd0, collapse = "&"))
    AmicTbl[[kol]] <- 10^(-PG[[paste0(pvalue.col[which(pvalue.use)], g)]])
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2)), kol] <- NaN
    kol <- paste0("adj.P.Val_", gEd, "__vs__", paste(gEd0, collapse = "&"))
    PVkol <- paste0(pvalue.col[which(pvalue.use)], g)
    AmicTbl[[kol]] <- p.adjust(10^(-PG[[PVkol]]), method = "BH")
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2)), kol] <- NaN
    kol <- paste0("logFC_", gEd, "__vs__", paste(gEd0, collapse = "&"))
    AmicTbl[[kol]] <- PG[[paste0("Mean ", Prot.Rat.Root, g)]]
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2)), kol] <- NaN
    kol <- paste0("AveExpr_", gEd, "__vs__", paste(gEd0, collapse = "&"))
    AmicTbl[[kol]] <- PG[[paste0("Mean ", prtRfRoot, g)]]/log10(2)
    AmicTbl[which(!is.all.good(AmicTbl[[kol]], 2)), kol] <- NaN
  }
  tst <- apply(AmicTbl[, grep("^AveExpr_", colnames(AmicTbl), value = TRUE), drop = FALSE], 1, function(x) {
    length(is.all.good(x))
  }) > 0
  AmicTbl$quantified <- c("", "+")[((AmicTbl$razorUniqueCount >= 2)&(tst))+1]
  AmicTbl <- AmicTbl[which(AmicTbl$quantified == "+"),]
  dir <- paste0(wd, "/Amica")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  data.table::fwrite(AmicTbl, paste0(wd, "/Amica/Amica_file.csv"), row.names = FALSE, na = "NaN", sep = "\t", quote = FALSE)
  data.table::fwrite(AmicaDesign, paste0(wd, "/Amica/Experimental_design.csv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")
  #data.table::fwrite(AmicTbl[1:500,], paste0(wd, "/Amica/Amica_file_short.csv"), row.names = FALSE, na = "NaN", quote = FALSE)
}

saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Venn diagrams
Src <- paste0(libPath, "/extdata/R scripts/Sources/Venn_diagrams.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Coverage maps, XICs and heatmaps for proteins of interest
protlspep <- prot.list_pep
if (length(protlspep)) {
  test <- vapply(protlspep, function(i) { length(grsep2(i, PG$"Leading protein IDs")) }, 1)
  if (0 %in% test) {
    w <- which(test == 0)
    for (w1 in w) {
      m <- match(protlspep[w1], db$"Protein ID")
      nm <- paste0(db$"Protein ID"[m], " - ", db$"Common Name"[m])
      warning(paste0("Protein of interest ",nm, " was not found in the dataset!"))
    }
    protlspep <- protlspep[which(test > 0)]
  }
}
if (length(protlspep)) { # Coverage
  setwd(wd) # To make sure we are in the working directory
  dir <- paste0(wd, "/Coverage")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  xKol <- paste0(pep.ref[length(pep.ref)], Exp.map$Ref.Sample.Aggregate) 
  tmpDB <- db[match(protlspep, db$"Protein ID"), c("Common Name", "Protein ID", "Sequence")]
  tst <- lapply(protlspep, function(x) { grsep2(x, pep$Proteins) })
  w <- which(vapply(tst, length, 1) > 0)
  prots <- protlspep[w]
  tst <- unique(unlist(tst))
  tmpPep <- pep[tst, c("Proteins", "Modified sequence", xKol)]
  source(parSrc, local = FALSE)
  clusterExport(parClust, list("tmpDB", "tmpPep", "pep.ref", "xKol", "wd", "VPAL", "Exp.map", "Exp"), envir = environment())
  lst <- parLapply(parClust, prots, function(i) { #i <- prots[1]
    nm <- tmpDB$"Common Name"[match(i, tmpDB$"Protein ID")]
    nm <- gsub("[<>:\"/\\\\\\|\\?]", "-", nm)
    if (nchar(nm) > 20) { nm <- paste0(gsub(" $", "", substr(nm, 1, 17)), "...") }
    seq <- setNames(tmpDB$Sequence[which(tmpDB$"Protein ID" == i)],
                    paste(tmpDB[which(tmpDB$"Protein ID" == i), c("Protein ID", "Common Name")], collapse = " - "))
    grs <- proteoCraft::grsep2(i, tmpPep$Proteins)
    drLst <- dir <- paste0(wd, "/Coverage/", i)
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    p <- tmpPep[grs,]
    m <- apply(p[, xKol], 1, function(x) {
      10^mean(proteoCraft::is.all.good(log10(x)))
    })
    w <- which(m == 0)
    if (length(w)) { stop("I didn't expect this, investigate!")}
    p1 <- data.frame(Sequence = p$"Modified sequence",
                     Intensity = m)
    print(proteoCraft::Coverage(seq, p1$Sequence))
    ttl <- gsub(":|/", "-", names(seq))
    setwd(dir) # To control precisely where it is saved
    proteoCraft::Coverage(seq, p1$Sequence, Mode = "Align2", title = paste0("Coverage map - ", nm), save = c("jpeg", "pdf"),
                    intensities = p1$Intensity, display = FALSE)
    setwd(wd) # To make sure I return to the working directory
    for (j in VPAL$values) { #j <- VPAL$values[1]
      sm <- Exp.map[which(Exp.map[[VPAL$column]] == j),]
      m <- apply(p[, paste0(pep.ref[length(pep.ref)], sm$Ref.Sample.Aggregate)], 1, function(x) {
        10^mean(proteoCraft::is.all.good(log10(x)))
      })
      w <- which(m > 0)
      if (length(w)) {
        p1 <- data.frame(Sequence = p$"Modified sequence"[w],
                         Intensity = m[w])
        dir2 <- paste0(dir, "/", gsub(":|\\*|\\?|<|>|\\|", "-", proteoCraft::cleanNms(j, rep = "_")))
        if (!dir.exists(dir2)) { dir.create(dir2, recursive = TRUE) }
        drLst <- unique(c(drLst, dir2))
        ttl <- paste0(gsub(":|/", "-", names(seq)), " - ", proteoCraft::cleanNms(j, rep = "_"))
        setwd(dir2) # To control precisely where it is saved
        proteoCraft::Coverage(seq, p1$Sequence, Mode = "Align2", save = c("jpeg", "pdf"), title = paste0("Coverage map - ", nm),
                        intensities = p1$Intensity, display = FALSE)
        setwd(wd) # To make sure I return to the working directory
        for (k in sm$Ref.Sample.Aggregate) { #k <- sm$Ref.Sample.Aggregate[1]
          m <- p[[paste0(pep.ref[length(pep.ref)], k)]]
          w <- which(m > 0)
          if (length(w)) {
            p1 <- data.frame(Sequence = p$"Modified sequence"[w],
                             Intensity = m[w])
            ttl <- paste0(gsub(":|/", "-", names(seq)), " - ", proteoCraft::cleanNms(k, rep = "_"))
            setwd(dir2) # To control precisely where it is saved
            proteoCraft::Coverage(seq, p1$Sequence, Mode = "Align2", save = c("jpeg", "pdf"),
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
        XICs <- parLapply(parClust, XIC_fls, function(x) {
          res <- arrow::read_parquet(x)
          res$"Mod. seq." <- proteoCraft::gsub_Rep("[0-9]+$", "", res$pr)
          res <- res[which(res$"Mod. seq." %in% u),]
          nm <- gsub(".*/|\\.xic\\.parquet$", "", x)
          res$File <- nm
          res$"Seq_Run" <- do.call(paste, c(res[, c("pr", "File")], sep = ">>>"))
          res$File <- factor(res$File, levels = tmp)
          res <- res[which(res$feature != "index"),]
          return(res)
        })
        #View(XICs[[1]])
        XICs <- plyr::rbind.fill(XICs)
        #View(XICs[1:100,])
        #
        m <- match(XICs$"Mod. seq.", ev$"Mod. seq. (DiaNN format)")
        XICs[, c("Proteins", "Sequence", "ModSeq", "PEP", "Quantity Quality")] <- ev[m,
                                                                                     c("Proteins", "Sequence", "Modified sequence", "PEP", "Quantity Quality")]
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
        for (pr in protlspep) { #pr <- protlspep[1]
          xicDir2 <- paste0(wd, "/XIC/", pr)
          if (!dir.exists(xicDir2)) { dir.create(xicDir2, recursive = TRUE) }
          dirlist <- unique(c(dirlist, xicDir2))
          g <- grsep2(pr, XICs$Proteins)
          if (length(g)) {
            XIC <- XICs[g,]
            pkBnds <- Boundaries[which(Boundaries$Seq_Run %in% XIC$Seq_Run),]
            u <- unique(XIC$ModSeq)
            clusterExport(parClust, list("XIC", "pkBnds", "xicDir2", "pr"), envir = environment())
            invisible(parLapply(parClust, u, function(sq) { #sq <- u[1] #sq <- u[2]
              sq2 <- gsub("^_|_$", "", sq)
              ppXIC <- XIC[which(XIC$ModSeq == sq),]
              yMax <- aggregate(ppXIC$value, list(ppXIC$File), max)
              xMin <- min(ppXIC$rt)
              bnds <- pkBnds[which(pkBnds$Seq_Run %in% ppXIC$Seq_Run),]
              bnds$yMax <- yMax$x[match(bnds$File, yMax$Group.1)]
              wMS1 <- which(ppXIC$feature == "ms1")
              wMS2 <- which(ppXIC$feature != "ms1")
              aNNOt <- aggregate(ppXIC[, c("PEP", "Quantity Quality")], list(ppXIC$File), function(x) { signif(mean(x, na.rm = TRUE), 3) })
              colnames(aNNOt)[1] <- "File"
              aNNOt$PEP <- paste0("PEP = ", aNNOt$PEP)
              aNNOt$"Quantity Quality" <- paste0("Quantity Quality = ", aNNOt$"Quantity Quality")
              aNNOt$Text <- do.call(paste, c(aNNOt[, c("PEP", "Quantity Quality")], sep = "\n"))
              aNNOt$y <- yMax$x[match(aNNOt$File, yMax$Group.1)]
              plot <- ggplot2::ggplot() + ggplot2::scale_y_continuous(expand = c(0, 10))
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
              #proteoCraft::poplot(plot, 12, 22)
              #
              suppressMessages({
                ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".jpeg"), plot, dpi = 450, height = 10, width = 10)
                ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".pdf"), plot, height = 10, width = 10)
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
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))

#### Code chunk - peptide tables for visualizing the coverage of proteins of interest in 3D using SCV
# Edit me to also just use the newer cov3D function!!!
if (length(protlspep)) {
  # From https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r
  valid_url <- function(url_in, t = 2){
    con <- url(url_in)
    check <- suppressWarnings(try(open.connection(con, open = "rt", timeout = t), silent = TRUE)[1])
    suppressWarnings(try(close.connection(con), silent = TRUE))
    ifelse(is.null(check), TRUE, FALSE)
  }
  # Mods and their mass shifts
  SCV_PTMs <- TRUE
  if (!"Mass shift" %in% colnames(Modifs)) {
    if ("UniMod" %in% colnames(Modifs)) {
      if (!require("unimod", quietly = TRUE)) { suppressMessages(devtools::install_github("rformassspectrometry/unimod")) }
      require(unimod)
      UniMod <- unimod::modifications
      Modifs$"Mass shift" <- UniMod$MonoMass[match(Modifs$UniMod, UniMod$UnimodId)]
    } else {
      if ("MAXQUANT" %in% SearchSoft) {
        if ((length(MQFold) == 1)&&(dir.exists(MQFold))) {
          modFls <- paste0(MQFold, "/bin/conf/modifications", c("", ".local"), ".xml")
          modFls <- modFls[which(file.exists(modFls))]
        } else {
          if ((exists("mqFld"))&&(dir.exists(mqFld))) { dflt <- mqFld } else { dflt <- "C:" }
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
          modFls <- lapply(modFls, function(modFl) { #modFl <- modFls[1]
            xml_lst <- as_list(read_xml(modFl))
            xml_lst <- xml_lst[[1]]
            xml_lst <- as.data.frame(t(sapply(xml_lst, function(x) {
              #x <- xml_lst[[1]]
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
            w <- which(vapply(colnames(Modifs), function(x) { "list" %in% class(Modifs[[x]]) }, 1))
            for (i in w) { temp[[i]] <- vapply(temp[[i]], paste, 1, collapse = ", ") }
            write.csv(temp, "Workflow control/Modifications.csv", row.names = FALSE)
          }
          Modifs$"Mass shift" <- vapply(strsplit(Modifs$Composition, " "), function(x) {
            #x <- strsplit(Modifs$Composition, " ")[1]
            x <- unlist(x)
            x <- as.data.frame(t(sapply(strsplit(gsub("\\)$", "", x), "\\("), function(y) {
              if (length(y) == 1) { y <- c(y, 1) }
              return(y)
            })))
            m <- match(x[[1]], IsotopeProbs$Atom)
            stopifnot(sum(is.na(m)) == 0)
            # For now the code above throws an error if an elements is missing from the table
            # If it ever does, I should expand the table to add isotopic probabilities for more elements!!!
            x <- sum(as.numeric(gsub("_.+", "", IsotopeProbs$Monoisotopic[m]))*as.integer(x[[2]]))
            return(x)
          }, 1)
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
  #for (plp in protlspep) { #plp <- protlspep[1]
  Tst <- parSapply(parClust, protlspep, function(plp) { #plp <- protlspep[1]
    grs <- proteoCraft::grsep2(plp, prVect)
    OutCome <- FALSE
    if (length(grs)) {
      dir <- paste0(wd, "/Coverage/", plp)
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      # For each protein we want to download:
      # - all models for all fragments.
      # - latest version only!
      # PDB: we get PDB IDs from parsing the txt file
      kPBD <- 0
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
      kAlpha <- 0
      tstF <- TRUE # Continue looking for the next fragment?
      while (tstF) {
        kAlpha <- kAlpha+1
        kV <- 0
        while ((kV == 0)||(tstV)) {
          kV <- kV + 1
          mdlNm <- paste0("AF-", plp,"-F", kAlpha, "-model_v", kV, ".pdb")
          url <- paste0("https://alphafold.ebi.ac.uk/files/", mdlNm)
          tstV <- valid_url(url) # We want to find out which is the latest v version of a model for that protein
        }
        kV <- kV - 1 # The last is always a failure
        if (kV) { # Did we find a valid url?
          mdlNm <- paste0("AF-", plp,"-F", kAlpha, "-model_v", kV, ".pdb")
          url <- paste0("https://alphafold.ebi.ac.uk/files/", mdlNm)
          download.file(url, paste0(dir, "/", mdlNm))
        }
        tstF <- kV > 0
      }
      kAlpha <- kAlpha-1 # The last is always a failure
      if (kPBD||kAlpha) {
        # We have found at least one model model which can be used to visualize coverage
        # Let's write peptidoforms
        tmp <- modSq[grs]
        if (SCV_PTMs) {
          tmp <- gsub("_", "", tmp)
          tmp <- gsub("\\)", "]_",gsub("\\(", "_[", tmp))
          tmp <- strsplit(tmp, "_")
          tmp <- vapply(tmp, function(x) { #x <- tmp[1]
            x <- unlist(x)
            w <- grep("\\[.+\\]", x)
            x[w] <- paste0("[", round(Modifs$"Mass shift"[match(x[w], paste0("[", Modifs$Mark, "]"))], 0), "]")
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
strngSrc <- paste0(libPath, "/extdata/R scripts/Sources/STRINGdb.R")
#rstudioapi::documentOpen(strngSrc)
source(strngSrc, local = FALSE)

#### Code chunk - For pull-downs: create table summarizing types of evidence for all proteins of interest
# if (IsPullDown) {
#   g <- paste0("Regulated - ", unique(Exp.map[which(!Exp.map$Reference), VPAL$column]))
#   test <- apply(PG[, g, drop = FALSE], 1, function(x) {
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
#       temp <- sapply(unique(ev$Type), function(x) {
#         vapply(unique(ev$"Raw file"), function(y) { length(which((temp$Type == x)&(temp$"Raw file" == y))) }, 1)
#       })
#       write.csv(temp, paste0(dir, "/Ev table - ", i, ".csv"))
#     }
#   }
# }

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
# Save decisions - should this go?
save(AllAnsw, file = "All_decisions.RData")

# Finalize reports
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

# Visualize results - already run earlier, could be expanded here to add more visualizations at further stages
#rstudioapi::documentOpen(xplorSrc)
#source(xplorSrc, local = FALSE)
