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
#ReUseAnsw %<o% FALSE
scrptType %<o% "noReps"
scrptTypeFull %<o% "noReps_PG_and_PTMs"
ExcelMax %<o% 32767L


# Parameters used by the start analysis script:
###-|-### Workflows: setNames(c("Discovery -> no comparisons, complex sample", "Band ID -> no comparisons, focus on coverage and the top proteins", "Regulation -> ratio analysis (up and down)", "Pull-down (incl. BioID) -> ratio analysis (choice between up only or up and down)"), c("Discovery", "Band ID", "Regulation", "Pull-down"))
###-|-### Replicates? FALSE
###-|-### External dependencies: Excel (loose); ScanHeadsman (loose)

### Packages
## For convenience all (or most) of the packages used are installed here:
## CRAN packages:
if(!exists("cran_req")) { cran_req <- "pak" }
cran_req %<o% cran_req
if(!exists("bioc_req")) { bioc_req <- c() } 
bioc_req %<o% bioc_req
cran_req <- unique(c(cran_req, "pak", "fs", "shiny", "renv", "R.utils", "data.table", "devtools", "qs2", "shinyWidgets", "DT", "shinyBS", "stringr",
                     "gplots", "ggplot2", "ggpubr", "reshape", "reshape2", "compiler", "stats", "rgl", "ggrepel", "rstudioapi", "gtools", "minpack.lm",
                     "parallel", "openxlsx", "openxlsx2", "openssl", "plotly", "Peptides", "venn", "ggdendro", "ggpubr", "colorspace", "ggnewscale",
                     "viridis", "factoextra", "NbClust", "gridExtra", "svDialogs", "htmlwidgets", "magrittr", "tibble", "fs", "officer", "snow",
                     "imputeLCMD", "ggplotify", "cowplot", "plyr", "shinyjs", "shinyFiles", "TeachingDemos", "shinycssloaders", "jpeg", "stringi",
                     "readr", "ssh", "taxize", "arrow", "iq", "Rtsne"))
bioc_req <- unique(c(bioc_req, "UniProt.ws", "pcaMethods", "impute", "GO.db", "topGO", "pcaMethods",
                     "limpa", "QFeatures"))
inst <- as.data.frame(installed.packages())
for (pack in cran_req) {
  if (!pack %in% inst$Package) {
    if (pack %in% c("pak", #"shiny",
                    "uchardet", #"openxlsx2",
                    "taxize")) {
      # Exceptions where for now we want a specific version to be installed,
      # or have to help the installer so it finds the right location
      if (pack == "pak") {
        install.packages("pak")
      }
      # if (pack == "shiny") { # Should be fixed now
      #   install.packages("https://cran.r-project.org/src/contrib/Archive/shiny/shiny_1.7.5.tar.gz")
      # }
      if (pack == "uchardet") {
        url <- "https://cran.r-project.org/src/contrib/Archive/uchardet/uchardet_1.1.1.tar.gz"
        destfile <- "uchardet_1.1.1.tar.gz"
        tst <- try(download.file(url, destfile, "curl"), silent = TRUE)
        if (inherits(tst, "try-error")) { try(download.file(url, destfile, "wget"), silent = TRUE) }
        install.packages(destfile)
        unlink(destfile)
      }
      # if (pack == "openxlsx2") {
      #   pak::pkg_install("JanMarvin/openxlsx2@v1.10", ask = FALSE) # ... until I can figure out what is happening...
      # }
      # if (pack == "myTAI") {
      #   pak::pkg_install("drostlab/myTAI@v0.9.3", ask = FALSE)
      # }
      if (pack == "taxize") {
        pak::pkg_install("ropensci/bold", ask = FALSE)
        pak::pkg_install("ropensci/taxize", ask = FALSE)
      }
    } else {
      pak::pkg_install(pack, ask = FALSE)
      #install.packages(pack)
    }
    inst <- as.data.frame(installed.packages())
  }
}
## Bioconductor packages:
biocInstall %<o% function(pack, load = TRUE) {
  inst <- as.data.frame(installed.packages())
  if (!pack %in% inst$Package) {
    pak::pkg_install(pack, ask = FALSE)
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

# shiny used to cause issues with View() for data.frames, possibly by importing jsonlite,
# which I know also caused the same problem (and is also imported by DT)
# The problem seems fixed, but if it re-appears, do not load the full packages,
# instead call functions using packageName::function()
ScriptPath %<o% normalizePath(gtools::script_file(), winslash = "/")
RunByMaster %<o% grepl(" - master script\\.R$", ScriptPath)
if (RunByMaster) { ScriptPath <- BehindTheScenes$ScriptFile }
Script %<o% readLines(ScriptPath)

# Update the proteoCraft package?
# msg <- "Should we update the proteoCraft package?"
# updt_proteoCraft %<o% c(TRUE, FALSE)[match(svDialogs::dlg_message(msg, "yesno")$res, c("yes", "no"))]
updt_proteoCraft %<o% FALSE

# Define input, output, project folder etc...
Src <- paste0(libPath, "/extdata/Sources/Start_analysis.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

AnalysisParam %<o% list("Input folder" = inDirs,
                        "Temp. folder" = wd,
                        "N. of threads" = N.clust)

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

#### Code chunk - MS files map
Src <- paste0(libPath, "/extdata/Sources/noRep_Fractions_Map_editor.R")
#rstudioapi::documentOpen(Src)
tstFrMp <- FALSE
while (!tstFrMp) {
  source(Src, local = FALSE)
}

#FracMap <- read.csv(FracMapPath, check.names = FALSE)
# No need to reload the local copy, values are updated in environment
m <- match(FracMap$`Raw files name`[which(FracMap$Use)], rawFiles2)
m <- m[which(!is.na(m))]
rawFiles <- rawFiles[m]
rawFiles2 <- rawFiles2[m]
FracMap <- FracMap[which(FracMap$Use),]
ev <- ev[which(ev$`Raw file path` %in% FracMap$`Raw file`),]
ev$Experiment <- FracMap$`Parent sample`[match(ev$`Raw file path`, FracMap$`Raw file`)]
exp <- unique(FracMap$"Parent sample")
if (length(rawFiles) < nrow(FracMap)) {
  warning("Some samples included in the analysis did not result in any identication!")
}
#
AnalysisParam$Type <- WorkFlow
MakeRatios %<o% FALSE
RatiosThresh_2sided %<o% TRUE
# No need to reload the local copy, values are updated in environment
if (WorkFlow %in% c("Regulation", "Pull-down")) {
  if (length(unique(FracMap$`Parent sample`)) == 1L) {
    warning("Only one sample, skipping ratios analysis!")
    WorkFlow <- "Discovery"
  } else { MakeRatios <- TRUE }
}

#### Code chunk - Samples/Experiment map
Src <- paste0(libPath, "/extdata/Sources/noRep_Experiment_Map_editor.R")
#rstudioapi::documentOpen(Src)
tstXpMp <- FALSE
while (!tstXpMp) {
  source(Src, local = FALSE)
}

# Labeling
AnalysisParam$"Label type" <- LabelType
int.cols %<o%  c()
#if (LabelType == "Isobaric") { int.cols["Original"] <- int.col <- paste0("Reporter intensity ", Labels) }
if (LabelType == "LFQ") {
  int.col %<o% "Intensity"
  int.cols["Original"] <- int.col
}

if (WorkFlow == "Pull-down") {
  IsBioID %<o% c(TRUE, FALSE)[match(dlg_message("Is this a BioID (and variants) experiment?", "yesno")$res, c("yes", "no"))]
  if (IsBioID) { AnalysisParam$"Type - advanced" <- "BioID" }
} else { IsBioID %<o% FALSE }
IsBioID2 %<o% FALSE

#### Code chunk - Load and process search database(s)
Src <- paste0(libPath, "/extdata/Sources/Process_Fasta_DBs.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

evNm %<o% "PSM"
#evNm %<o% c("PSM", "Evidence")[(SearchSoft == "MAXQUANT")+1L]

#### Code chunk - Load and process annotations
## This includes a QC step in case the database differs slightly from the one used by MQ, or if somehow some IDs have not been properly parsed.
Src <- paste0(libPath, "/extdata/Sources/Load_Annotations.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
source(parSrc, local = FALSE)
Src <- paste0(libPath, "/extdata/Sources/GO_prepare.R") # Doing this earlier but do keep latter instance for now
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
AnalysisParam$Annotations <- Annotate

#### Code chunk - Define analysis parameters
paramSrc <- paste0(libPath, "/extdata/Sources/noRep_Parameters_editor_Main.R")
#rstudioapi::documentOpen(paramSrc)
source(paramSrc, local = FALSE)

ev$"Protein group IDs" <- NULL
ev$"Peptide ID" <- NULL

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

# Create Materials and Methods template
Src <- paste0(libPath, "/extdata/Sources/autoMatMet.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Update peptide-to-protein mappings
Src <- paste0(libPath, "/extdata/Sources/checkPep2Prot.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
w <- which(ev$Proteins == "")
if (length(w)) {
  warning(paste0("Removing ", length(w), " peptide evidences with no matches to the search database."))
  ev <- ev[which(ev$Proteins != ""),]
}

# Gene names
# Simplify Gene columns
genkol %<o% c("Gene", "Gene name")
w <- which(genkol %in% colnames(db))
temp <- db[, genkol[w]]
if (length(w) == 2L) {
  for (i in genkol) { temp[[i]] <- strsplit(temp[[i]], ";") }
  temp <- apply(temp, 1L, \(x) { paste(sort(unique(unlist(x))), collapse = ";") })
}
temp2 <- as.data.table(listMelt(strsplit(ev$Proteins, ";"), 1L:nrow(ev)))
temp2$Gene <- temp[match(temp2$value, db$`Protein ID`)]
temp2 <- temp2[, list(Genes = paste(Gene, collapse = ";")), keyby = list(id = L1)]
temp2 <- as.data.frame(temp2)
ev$"Gene names" <- temp2$Genes[match(ev$id, temp2$id)]

# Deal with PTM-enriched data
FracMap$`PTM-enriched`[which(FracMap$`PTM-enriched` == "NA")] <- NA
EnrichedPTMs %<o% unique(FracMap$`PTM-enriched`[which(!is.na(FracMap$`PTM-enriched`))])
PTMriched %<o% (length(EnrichedPTMs) > 0L)
PTMev2Remov %<o% list()
if (PTMriched) {
  for (ptm in EnrichedPTMs) {
    dir <- paste0(wd, "/Workflow control")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    ev[[ptm]] <- grepl(paste0("\\(", Modifs$Mark[match(ptm, Modifs$`Full name`)], "\\)"),
                       ev$`Modified sequence`)
    ttl <- paste0(ptm, "-enrichment efficiency")
    w <- which(ev$`Raw file path` %in% FracMap$`Raw file`[which(FracMap$`PTM-enriched` == ptm)])
    tmp <- ev[w, c("Raw file", ptm)]
    tmp <- aggregate(tmp[[ptm]], list(tmp$`Raw file`, tmp[[ptm]]), length)
    colnames(tmp) <- c("Raw file", ptm, "Count")
    tmp[[ptm]] <- c("-", "+")[tmp[[ptm]] + 1L]
    tmp[[ptm]] <- factor(tmp[[ptm]], levels = c("-", "+"))
    tmp$"Raw file" <- gsub(".*/|\\.[^\\.]+$", "", tmp$"Raw file")
    plot <- ggplot(tmp) + geom_bar(aes(x = `Raw file`, y = Count, fill = .data[[ptm]]), stat = "identity") +
      scale_fill_viridis(discrete = TRUE, option = "H", begin = 0.25, end = 0.8) + ggtitle(ttl) + 
      theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
    poplot(plot)
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150L)
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150L)
    })
    #
    PTMev2Remov[[ptm]] <- lapply(Exp, \(exp) { #exp <- Exp[1L]
      res <- c()
      fm <- FracMap[which(FracMap$`Parent sample` == exp),]
      tst <- unique(fm$`PTM-enriched`[which(fm$`Parent sample` == exp)])
      if ((ptm %in% tst)&&(length(tst) > 1L)) {
        # We need to filter out peptides with the mark from all non-enriched samples
        smpls1 <- fm$`Raw file`[which((is.na(fm$`PTM-enriched`))|(fm$`PTM-enriched` != ptm))] # Non-enriched/flow through samples
        smpls2 <- fm$`Raw file`[which(fm$`PTM-enriched` == ptm)] # Enriched samples
        if (!smpls1 %in% unique(ev$`Raw file path`)) { # Preempt bugs
          stop("Debug me pleeeeeeease!!!")
        }
        # From smpls1 we remove peptides bearing the PTM
        res <- ev$id[which((ev[[ptm]])&(ev$`Raw file path` %in% smpls1))]
      }
      # And from smpls2 (enriched samples) we remove peptides not bearing the PTM
      res <- c(res, ev$id[which((!ev[[ptm]])&(ev$`Raw file path` %in% smpls2))])
      return(res)
    })
  }
}
PTMev2Remov <- unlist(PTMev2Remov)
l <- length(PTMev2Remov)
if (l) {
  msg <- paste0("Removing ", l, " PSMs (", signif(100*l/nrow(ev), 2L), "%) with enriched PTMs from non-enriched samples.")
  ev <- ev[which(!ev$id %in% PTMev2Remov),]
}

# Negative filter
if (!"Negative Filter" %in% colnames(SamplesMap)) {
  SamplesMap$"Negative Filter" <- FALSE
}
SamplesMap$"Negative Filter"[which(is.na(SamplesMap$"Negative Filter"))] <- FALSE
NegFilt %<o% (sum(SamplesMap$`Negative Filter`) > 0L)

# Filter PSMs
# Optional (not used): remove "wrongly assigned" charge 1 evidences
RemovZ1 %<o% FALSE
w1 <- which(ev$Charge == 1L)
wHt1 <- which(ev$Charge > 1L)
if ((RemovZ1)&&(length(w1))) {
  cat("Removing the following presumably bogus identifications with Z=1:\n",
      paste(unique(ev$`Modified sequence`[w1]), collapse = "\n"))
  ev <- ev[wHt1,]
}
#
# Remove evidences with null intensity values
w1 <- which((is.all.good(ev$Intensity, 2L))&(ev$Intensity > 0))
w2 <- which(!1L:nrow(ev) %in% w1)
l2 <- length(w2)
if (l2) {
  warning(paste0("Removing ", l2, " (", round(100*l2/nrow(ev), 2L), "%) peptide evidences with invalid or null intensity values!"))
  nullEv <- ev[w2,]
  ev <- ev[w1,]
}
w1 <- which((is.na(ev$Reverse))|(ev$Reverse != "+"))
w2 <- which(!1L:nrow(ev) %in% w1)
l2 <- length(w2)
if (l2) {
  warning(paste0("Removing ", l2, " (", round(100*l2/nrow(ev), 2L), "%) reverse peptide evidences!"))
  revEv <- ev[w2,]
  ev <- ev[w1, ]
}
w1 <- which((is.na(ev$"Potential contaminant"))|(ev$"Potential contaminant" != "+"))
w2 <- which(!1L:nrow(ev) %in% w1)
l2 <- length(w2)
if (l2) {
  message(paste0(round(100*l2/nrow(ev), 2L), "% of valid identifications are potential contaminants."))
  #  contEv <- ev[-w, ]
  #  ev <- ev[w, ]
}

Exp <- expOrder[which(expOrder %in% ev$Experiment)] # Update experiments

# DIA-only: MS2-based correction of MS1-based quantitative values
Src <- paste0(libPath, "/extdata/Sources/MS2corr2MS1.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Pepper correction
# (Maybe this should be done at peptides level?)
Src <- paste0(libPath, "/extdata/Sources/run_Pepper.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Create modified peptides table
Src <- paste0(libPath, "/extdata/Sources/pepMake.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
for (e in Exp) { #e <- Exp[1L]
  temp <- ev[which(ev$Experiment == e),]
  temp <-  set_colnames(aggregate(temp[[int.col]], list(temp$"Modified sequence"), \(x) { sum(is.all.good(x)) }),
                        c("Modified sequence", int.col))
  pep[, paste0(int.col, " - ", e)] <- 0
  w <- which(pep$"Modified sequence" %in% temp$"Modified sequence")
  pep[w, paste0(int.col, " - ", e)] <- temp[match(pep$"Modified sequence"[w], temp$"Modified sequence"), int.col]
}
ev_to_pep <- match(pep$`Modified sequence`, ev$`Modified sequence`)
if (PTMriched) {
  pep[, EnrichedPTMs] <- ev[ev_to_pep, EnrichedPTMs]
}

#### Code chunk - optionally impute missing expression values
if ((length(Exp) > 1L)&&(Impute)) {
  kol <- grep(topattern(int.col), colnames(pep), value = TRUE)
  tst <- length(which(!is.all.good(log10(unlist(pep[, kol])), 2L)))
  if (length(tst)) {
    cat("Incomplete peptides-level expression values matrix.\nImputing missing values with random draws from a gaussian distribution centered on the lowest observed value and with SD = 1/5 that of the data.\n")
    g <- grep(topattern(paste0(int.col, " - ")), colnames(pep), value = TRUE)
    temp <- log10(as.matrix(pep[, g]))
    w1 <- which((is.na(temp) | is.infinite(temp)), arr.ind = TRUE)
    w2 <- which(!(is.na(temp) | is.infinite(temp)), arr.ind = TRUE)
    Min <- min(unlist(temp[w2]))
    SD <- sd(unlist(temp[w2]))
    temp[w1] <- rnorm(nrow(w1), Min, SD/5)
    kol2 <- gsub(topattern(int.col), paste0("Imput. ", int.cols["Original"]), g)
    colnames(temp) <- kol2
    temp <- 10L^temp
    temp[w2] <- pep[, g][w2]
    temp <- as.data.frame(temp)
    temp[[paste0("Imput. ", int.cols["Original"])]] <- rowSums(temp[, kol2])
    pep[, colnames(temp)] <- temp
    int.cols["Imputed"] <- int.col <- paste0("Imput. ", int.cols["Original"])
  }
}

if (prot.list.Cond) {
  pep$`In list` <- ""
  g <- grsep2(prot.list, pep$"Proteins")
  pep$`In list`[g] <- "+"
  if ("Potential contaminant" %in% colnames(pep)) {
    pep$"Potential contaminant"[g] <- ""
  }
  ev$`In list` <- ""
  g <- which(ev$`Modified sequence` %in% pep$`Modified sequence`[g])
  ev$`In list`[g] <- "+"
  if ("Potential contaminant" %in% colnames(ev)) {
    ev$"Potential contaminant"[g] <- ""
  }
}

## Peptides level:
# Intensity distribution:
kol <- c("Modified sequence", paste0(int.col, " - ", Exp))
kol2 <- "Modified sequence"
form <- ".~Sample"
if (PTMriched) {
  kol <- c(kol, EnrichedPTMs)
  kol2 <- c(kol2, EnrichedPTMs)
  form <- gsub("^\\.", "PTMs", form)
}
form <- as.formula(form)
temp <- pep[, kol]
temp <- dfMelt(temp, id.vars = kol2)
temp$Sample <- gsub_Rep(topattern(paste0(int.col, " - ")), "", temp$variable)
temp$variable <- NULL
temp$Sample <- factor(temp$Sample, levels = Exp)
temp$value <- log10(temp$value)
temp <- temp[which(is.all.good(temp$value, 2L)),]
if (PTMriched) {
  for (ptm in EnrichedPTMs) { temp[[ptm]] <- c("", ptm)[temp[[ptm]]+1L] }
  if (length(EnrichedPTMs) > 1L) {
    temp$PTMs <- do.call(paste, c(temp[, EnrichedPTMs], sep = "-"))
  } else {
    temp$PTMs <- temp[[EnrichedPTMs]]
  }
}
temp$Sample <- factor(temp$Sample, levels = Exp)
ttl <- "Density plot - Peptides level"
if (prot.list.Cond) {
  temp$"In list" <- pep$"In list"[match(temp$`Modified sequence`, pep$`Modified sequence`)]
  plot <- ggplot(temp) + geom_histogram(aes(x = value, fill = `In list`), bins = 100L)
} else {
  plot <- ggplot(temp) + geom_histogram(aes(x = value, fill = Sample), bins = 100L)
}
plot <- plot +
  facet_grid(form) +
  theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
  scale_y_continuous(expand = c(0L, 0L)) +
  scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
  ggtitle(ttl) + xlab("log10(Peptides Intensity)")
poplot(plot, 12L, 22L)
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
})

# Correlation:
if (length(Exp) > 1L) {
  temp <- pep[, c("Modified sequence", paste0(int.col, " - ", Exp))]
  for (e in Exp) {
    temp[[paste0("log10(Intensity) - ", e)]] <- log10(temp[[paste0(int.col, " - ", e)]])
  }
  comb <- as.data.frame(gtools::combinations(length(Exp), 2L, Exp))
  temp2 <- temp[, grep(topattern("log10(Intensity) - "), colnames(temp), value = TRUE)]
  source(parSrc, local = FALSE)
  tmpFl <- tempfile(fileext = ".rds")
  clusterExport(parClust, "tmpFl", envir = environment())
  readr::write_rds(temp2, tmpFl)
  invisible(clusterCall(parClust, \(x) {
    assign("temp2", readr::read_rds(tmpFl), envir = .GlobalEnv)
    return()
  }))
  unlink(tmpFl)
  temp2 <- parApply(parClust, comb, 1L, \(x) {
    temp3 <- temp2[, paste0("log10(Intensity) - ", unlist(x))]
    temp3$X <- x[[1L]]
    temp3$Y <- x[[2L]]
    temp3$Comparison <- paste0(x[[1L]], " (X) vs ", x[[2L]], " (Y)")
    colnames(temp3)[1L:2L] <- c("log10(X intensity)", "log10(Y intensity)")
    return(temp3)
  })
  temp2 <- plyr::rbind.fill(temp2)
  temp2$"Modified sequence" <- temp$"Modified sequence"
  clusterExport(parClust, "is.all.good", envir = environment())
  test <- parApply(parClust, temp2[, c("log10(X intensity)", "log10(Y intensity)")], 1L, \(x) {
    length(is.all.good(x))
  }) == 2L
  temp2 <- temp2[which(test),]
  temp2$X <- factor(temp2$X, levels = Exp)
  temp2$Y <- factor(temp2$Y, levels = Exp)
  temp3 <- as.data.frame(t(sapply(unique(temp2$Comparison), \(x) { #x <- unique(temp2$Comparison)[1L]
    x1 <- temp2[which(temp2$Comparison == unlist(x)), c("log10(X intensity)", "log10(Y intensity)")]
    x1 <- x1$"log10(Y intensity)"-x1$"log10(X intensity)"
    return(setNames(c(x, paste0("Median = ", round(median(x1), 3L)), paste0("S.D. = ", round(sd(x1), 3L))),
                    c("Comparison", "Median", "SD")))
  })))
  temp3$X <- gsub(" \\(X\\).+", "", temp3$Comparison)
  temp3$Y <- gsub(".+\\(X\\) vs ", "", gsub(" \\(Y\\)$", "", temp3$Comparison))
  temp3$X <- factor(temp3$X, levels = Exp)
  temp3$Y <- factor(temp3$Y, levels = Exp)
  temp3$R <- vapply(strsplit(gsub(" \\([XY]\\)", "", temp3$Comparison), " vs "), \(x) {
    return(paste0("R = ", round(cor(pep[[paste0(int.col, " - ", x[[1L]])]],
                                    pep[[paste0(int.col, " - ", x[[2L]])]]), 3L)))
  }, "")
  x_min <- min(temp2$`log10(X intensity)`)
  y_min <- min(temp2$`log10(Y intensity)`)
  y_max <- max(temp2$`log10(Y intensity)`)
  # From: https://slowkow.com/notes/ggplot2-color-by-density/
  get_density %<o% function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  tmp2 <- temp2[, c("log10(X intensity)", "log10(Y intensity)", "Comparison", "Modified sequence")]
  tmpFl <- tempfile(fileext = ".rds")
  clusterExport(parClust, list("tmpFl", "get_density"), envir = environment())
  readr::write_rds(tmp2, tmpFl)
  invisible(clusterCall(parClust, \(x) {
    assign("tmp2", readr::read_rds(tmpFl), envir = .GlobalEnv)
    return()
  }))
  unlink(tmpFl)
  comps <- unique(temp2$Comparison)
  tmp2D <- setNames(parLapply(parClust, comps, \(cmp) { #cmp <- comps[1L]
    w <- which(tmp2$Comparison == cmp)
    x <- try({
      y <- get_density(tmp2$`log10(X intensity)`[w], tmp2$`log10(Y intensity)`[w], n = 500L)
      list(Success = TRUE, Density = data.frame(Density = y, Seq = tmp2$"Modified sequence"[w]))
    }, silent = TRUE)
    if (inherits(x, "try-error")) { x <- list(Success = FALSE) }
    return(x)
  }), comps)
  tmp2D <- tmp2D[which(vapply(tmp2D, \(x) { x$Success }, TRUE))]
  tmp2D <- lapply(tmp2D, \(x) { x$Density })
  if (length(tmp2D)) {
    comps <- names(tmp2D)
    temp2 <- temp2[which(temp2$Comparison %in% comps),]
    temp2$Density <- 0
    for (cmp in comps) { #cmp <- comps[1L]
      w <- which(temp2$Comparison == cmp)
      tmp2Dw <- tmp2D[[cmp]]
      temp2$Density[w] <- tmp2Dw$Density[match(temp2$`Modified sequence`[w], tmp2Dw$Seq)]
    }
    ttl <- "Correlation plot - Peptides level"
    plot <- ggplot(temp2) +
      geom_point(aes(x = `log10(X intensity)`, y = `log10(Y intensity)`,
                     colour = Density), size = 0.1) +
      geom_abline(intercept = 0, slope = 1, colour = "red") +
      geom_text(data = temp3, x = x_min, y = y_max - 0.01*(y_max - y_min),
                aes(label = R), hjust = 0, cex = 2.5) +
      geom_text(data = temp3, x = x_min, y = y_max - 0.06*(y_max - y_min),
                aes(label = Median), hjust = 0, cex = 2.5) +
      geom_text(data = temp3, x = x_min, y = y_max - 0.11*(y_max - y_min),
                aes(label = SD), hjust = 0, cex = 2.5) +
      facet_grid(Y~X, drop = TRUE) + coord_fixed(1L) +
      scale_color_viridis() +
      theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
                         strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) + ggtitle(ttl)
    poplot(plot, 12L, 22L)
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    })
  } else { warning("Do we really have enough data to continue? Investigate...") }
}
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))

# Calculate peptide ratios
rat.col %<o% "log2(Ratio)"
rat.cols %<o% c()
if (MakeRatios) {
  rat.cols["Original"] <- rat.col
  rat.grps %<o% unique(SamplesMap$`Ratios group`)
  rat.grps <- rat.grps[which(!is.na(rat.grps))]
  for (gr in rat.grps) { #gr <- rat.grps[1L]
    m <- SamplesMap[which(SamplesMap$`Ratios group` == gr),]
    if (sum(c(TRUE, FALSE) %in% m$Reference) < 2L) {
      warning(paste0("Ratios group ", gr, " - there should be reference (control) and non-reference samples!"))
    } else {
      ref <- apply(pep[, paste0(int.col, " - ", m$Experiment[which(m$Reference)]), drop = FALSE], 1L, \(x) {
        x <- is.all.good(x)
        l <- length(x)
        if (!l) { x <- NA } else { x <- prod(x)^(1/l) }
        return(x)
      })
      w <- which(!m$Reference)
      for (x in w) { #x <- w[1L]
        pep[[paste0(rat.col, " - ", m$Experiment[x])]] <- log2(pep[[paste0(int.col, " - ", m$Experiment[x])]]/ref)
      }
    }
  }
  # Peptide ratios distribution:
  dir <- paste0(wd, "/Workflow control")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  temp <- pep[, c("Modified sequence", grep(topattern(paste0(rat.col, " - ")), colnames(pep), value = TRUE))]
  temp <- dfMelt(temp, id.vars = "Modified sequence")
  temp$Sample <- gsub(topattern(paste0(rat.col, " - ")) , "", temp$variable)
  temp$Sample <- factor(temp$Sample, levels = Exp)
  temp <- temp[which(is.all.good(temp$value, 2L)),]
  ttl <- "Ratios density plot - Peptides level"
  plot <- ggplot(temp) + geom_histogram(aes(x = value, fill = Sample), bins = 100L) +
    geom_vline(xintercept = RatiosThresh, colour = "red") +
    ggtitle(ttl) + facet_grid(.~Sample) +
    scale_y_continuous(expand = c(0L, 0L)) +
    scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
    theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
    xlab("log2(Ratio)")
  if (RatiosThresh_2sided) { plot <- plot + geom_vline(xintercept = -RatiosThresh, colour = "red") }
  poplot(plot, 12L, 22L)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
  })
}

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Assemble protein groups
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

# Basic fix, because I do not like the way I was doing Quality filters up to now
g <- grep("^Quality filter: ", colnames(PG), value = TRUE)
if (length(g)) {
  for (h in g) { #h <- g[1L]
    PG[[h]] <- c("no -> dubious!", "")[match(PG[[h]], c("", "Keep"))]
  }
}
#
if (prot.list.Cond) {
  PG$"In list" <- ""
  g <- grsep2(prot.list, PG$"Protein IDs")
  PG$`In list`[g] <- "+"
  PG$"Potential contaminant"[g] <- ""
  ev$"Potential contaminant"[grsep2(prot.list, ev$Proteins)] <- ""
}
#
if (!"Protein group IDs" %in% colnames(ev)) {
  temp <- listMelt(strsplit(PG$`Evidence IDs`, ";"), PG$id)
  temp <- temp[order(temp$L1, decreasing = FALSE),]
  temp <- data.table(value = temp$value, L1 = temp$L1)
  temp <- temp[, list(x = paste(L1, collapse = ";")), by = list(Group.1 = value)]
  temp <- as.data.frame(temp)
  temp$Group.1 <- as.integer(temp$Group.1)
  ev$"Protein group IDs" <- temp$x[match(ev$id, temp$Group.1)]
}
if (tstOrg) {
  test <- listMelt(strsplit(PG$`Protein IDs`, ";"), PG$id)
  test$Org <- db[match(test$value, db$`Protein ID`), dbOrgKol]
  test <- test[order(test$Org, decreasing = FALSE),]
  test <- aggregate(test$Org, list(test$L1), \(x) { paste(unique(x), collapse = ";") })
  pgOrgKol %<o% c("Organism", "Organism(s)")[(sum(grepl(";", test))>0L)+1L]
  PG[[pgOrgKol]] <- test$x[match(PG$id, test$Group.1)]
}

# Some stats on protein groups
tmp <- aggregate(PG$id, list(PG$`Peptides count`), length) # Faster than data.table here
colnames(tmp) <- c("Peptides count", "Protein groups")
tmp$"log10(Protein groups count)" <- log10(tmp$"Protein groups")
pal <- colorRampPalette(c("brown", "yellow"))(max(tmp$"Peptides count")-1L)
tmp$Colour <- c("blue", pal)[tmp$`Peptides count`]
tmp2 <- summary(PG$`Peptides count`)
tmp2 <- data.frame(Variable = c(names(tmp2), "", "Protein groups", "Protein groups with 2+ peptidoforms"),
                   Value = c(as.character(signif(as.numeric(tmp2), 3L)),
                             "",
                             as.character(c(nrow(PG), sum(PG$"Peptides count" >= 2L)))))
tmp2$Txt <- apply(tmp2[, c("Variable", "Value")], 1L, \(x) {
  x <- x[which(x != "")]
  if (length(x)) { x <- paste(x, collapse = ": ") } else { x <- "" }
  return(x)
})
tmp2$X <- max(as.numeric(tmp2$Value[match("Max.", tmp2$Variable)]))*0.98
tmp2$Y <- max(tmp$"log10(Protein groups count)")*(0.98-(0L:(nrow(tmp2)-1L))*0.02)
ttl <- "Peptidoforms per PG"
plot <- ggplot(tmp) + geom_col(aes(x = `Peptides count`, y = `log10(Protein groups count)`, fill = Colour),
                               colour = NA) +
  geom_text(data = tmp2, aes(x = X, y = Y, label = Txt), hjust = 1, size = 3L) +
  scale_fill_identity() + theme_bw() + ggtitle(ttl)
poplot(plot)
dir <- paste0(wd, "/Summary plots")
#dirlist<- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300L)
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L)
})

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
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
# Simplify Gene columns
genkol <- c("Genes", "Gene names")
w <- which(genkol %in% colnames(PG))
if (length(w) == 2L) {
  temp <- PG[, genkol]
  for (i in genkol) { temp[[i]] <- strsplit(temp[[i]], ";") }
  PG$Genes <- apply(temp, 1L, \(x) { paste(sort(unique(unlist(x))), collapse = ";") })
  PG$"Gene names" <- NULL
} else { if (length(w) == 1L) { colnames(PG)[which(colnames(PG) %in% genkol)] <- "Genes" } }

# Here, if this is a BioID type experiment, we also want to mark protein groups which have Biotin peptides:
if (IsBioID) {
  wbiot <- grep("biot", Modifs$"Full name", ignore.case = TRUE)
  l <- length(wbiot)
  if (length(wbiot)) {
    tmp <- if (l == 1L) { Modifs$"Full name"[wbiot] } else {
      paste0(paste(Modifs$"Full name"[wbiot[1L:(l-1L)]], collapse = "\", \""), "\" and \"", Modifs$"Full name"[wbiot[l]])
    }
    warning(paste0("Modifications \"", tmp, "\" were detected as biotinylations, check that this is correct!"))
    AnalysisParam$"Type - advanced (mod)" <- list(Modifs$"Full name"[wbiot])
  } else {
    warning("I could not identify any biotinylated PTMs in the modifications table, did you include them in the search?")
    IsBioID <- FALSE
    AnalysisParam$"Type - advanced" <- "BioID... - except I could not find biotinylations among the PTMs!"
  }
}
if (IsBioID) {
  g <- grep(topattern(Modifs$Mark[wbiot], start = FALSE), pep$"Modified sequence")
  if (length(g)) {
    wpg <- unique(unlist(strsplit(pep$"Protein group IDs"[g], ";")))
    wpg <- which(PG$id %in% wpg)
    PG[["Biot. peptide IDs"]] <- ""
    temp <- listMelt(strsplit(PG$"Peptide IDs", ";"), PG$id)
    temp <- temp[which(temp$value %in% pep$id[g]),]
    temp <- aggregate(temp$value, list(temp$L1), \(x) { paste(sort(x), collapse = ";") })
    PG[wpg, "Biot. peptide IDs"] <- temp$x[match(PG$id[wpg], temp$Group.1)]
    PG[["Biot. peptides count"]] <- lengths(strsplit(PG[["Biot. peptide IDs"]], ";"))
    PG[["Biot. peptides [%]"]] <- round(100*PG[["Biot. peptides count"]]/PG$"Peptides count", 1L)
    IsBioID2 <- TRUE
  } else {
    warning("I could not find any biotinylated peptides!")
    IsBioID <- FALSE
    AnalysisParam$"Type - advanced" <- "BioID... - except that no biotinylated peptides were detected!"
  }
}

# Number of spectra, evidences and peptides per sample:
source(parSrc, local = FALSE)
invisible(clusterCall(parClust, \() {
  library(proteoCraft)
  library(reshape)
  library(data.table)
  return()
}))
temp_PG <- data.frame(id = PG$id,
                      Accession1 = vapply(strsplit(PG$"Leading protein IDs", ";"), \(x) { unlist(x)[1L] }, ""))
temp_PG$Pep <- parLapply(parClust, strsplit(PG$"Peptide IDs", ";"), as.integer)
tmp <- pep[, c("id", "Sequence")]
clusterExport(parClust, "tmp", envir = environment())
temp_PG$Pep <- parLapply(parClust, temp_PG$Pep, \(x) { tmp$Sequence[match(x, tmp$id)] })
temp_PG$Seq <- db$Sequence[match(temp_PG$Accession1, db$"Protein ID")]
if (!"Sequence coverage [%]" %in% colnames(PG)) {
  exports <- list("Coverage")
  clusterExport(parClust, exports, envir = environment())
  PG$"Sequence coverage [%]" <- round(100*parApply(parClust, temp_PG[, c("Seq", "Pep")], 1L, \(x) {
    Coverage(x[[1L]], x[[2L]])
  }), 1L)
}
CreateMSMSKol %<o% (("MS/MS IDs" %in% colnames(ev))&&(class(ev$"MS/MS IDs") %in% c("integer", "character")))
if (CreateMSMSKol) {
  # Somehow there are no MSMS IDs for DIA, but this may be a late tables writing bug.
  ev$temp <- parLapply(parClust, strsplit(as.character(ev$"MS/MS IDs"), ";"), as.integer)
  #PG[, paste0("Spectr", c("al count", "um IDs"))]
  temp <- reshape::melt(setNames(lapply(strsplit(PG$`Evidence IDs`, ";"), as.integer), PG$id))
  temp$MSMSIDs <- ev$temp[match(temp$value, ev$id)]
  temp <- temp[which(lengths(temp$MSMSIDs) > 0L),] # Remove Match-Between-Runs evidences (no MS/MS)
  temp <- reshape::melt(setNames(temp$MSMSIDs, temp$L1))
  temp <- do.call(data.frame, aggregate(temp$value, list(temp$L1), \(x) {
    x <- unique(x)
    return(c(Count = length(x), List = list(x)))
  }))
  temp$x.Count <- unlist(temp$x.Count)
  temp$Pasted <- vapply(temp$x.List, \(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") }, "")
  PG[, paste0("Spectr", c("al count", "um IDs"))] <- temp[match(PG$id, temp$Group.1), c("x.Count", "Pasted")]
  ev$temp <- NULL
}
temp_ev <- ev[, c("id", "Experiment", "Protein group IDs", "Peptide ID")]
if (CreateMSMSKol) { temp_ev$"MS/MS IDs" <- ev$"MS/MS IDs" }
temp_pep <- pep[, c("id", "Sequence")]
source(parSrc)
tmpFl1 <- tempfile(fileext = ".rds")
tmpFl2 <- tempfile(fileext = ".rds")
tmpFl3 <- tempfile(fileext = ".rds")
exports <- list("wd", "IsBioID2", "Exp", "Modifs", "CreateMSMSKol", "tmpFl1", "tmpFl2", "tmpFl3")
if (IsBioID2) { exports <- append(exports, "wbiot") }
clusterExport(parClust, exports, envir = environment())
readr::write_rds(temp_ev, tmpFl1)
readr::write_rds(temp_pep, tmpFl2)
readr::write_rds(temp_PG, tmpFl3)
invisible(clusterCall(parClust, \(x) {
  library(data.table)
  assign("temp_ev", readr::read_rds(tmpFl1), envir = .GlobalEnv)
  assign("temp_pep", readr::read_rds(tmpFl2), envir = .GlobalEnv)
  assign("temp_PG", readr::read_rds(tmpFl3), envir = .GlobalEnv)
  return()
}))
unlink(tmpFl1)
unlink(tmpFl2)
unlink(tmpFl3)
clusterExport(parClust, list("listMelt", "Coverage", "topattern"), envir = environment())
temp <- parLapply(parClust, Exp, \(exp) { #exp <- Exp[1L]
  res <- temp_PG[, "id", drop = FALSE]
  kol <- c()
  if (CreateMSMSKol) {
    kols <- paste0("Spectr", c("al count", "um IDs"), " - ", exp)
    kol <- c(kol, kols)
  }
  kole <- paste0("Evidence", c("s count", " IDs"), " - ", exp)
  kolp <- paste0("Peptide", c("s count", " IDs"), " - ", exp)
  kol <- c(kol, kole, kolp)
  kolk <- grep(" count - ", kol, value = TRUE)
  koli <- grep(" IDs - ", kol, value = TRUE)
  res[, kolk] <- 0L
  res[, koli] <- ""
  res[[paste0("Sequence coverage [%] - ", exp)]] <- 0
  w <- which(temp_ev$Experiment == exp)
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
    res[w, paste0("Sequence coverage [%] - ", exp)] <- round(100*apply(temp_PG[w, c("Seq", "Pep")], 1L, \(x) {
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
        temp1 <- reshape::melt(setNames(lapply(strsplit(eB$"Protein group IDs", ";"), as.integer), eB$"Peptide ID"))
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
#View(PG[, grep("^Spectr|^Peptide|^Evidence", colnames(PG))])

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

# Coverage columns
PG$"1st accession" <- vapply(strsplit(PG$`Leading protein IDs`, ";"), \(x) { unlist(x)[1L] }, "")
PG$"Sequence (1st accession)" <- db$Sequence[match(PG$`1st accession`, db$`Protein ID`)]
source(parSrc, local = FALSE)
tmpLst <- tmpRws <- c()
for (exp in Exp) { #exp <- Exp[1L]
  temp <- PG[, c("1st accession", "Sequence (1st accession)", paste0("Peptide IDs - ", exp))]
  w <- which(temp[[paste0("Peptide IDs - ", exp)]] != "")
  temp <- temp[w,]
  tmpRws[[exp]] <- w
  tmp <- listMelt(strsplit(temp[[paste0("Peptide IDs - ", exp)]], ";"), temp$`1st accession`)
  tmp$ModSeq <- pep$`Modified sequence`[match(as.integer(tmp$value), pep$id)]
  tmp <- aggregate(tmp$ModSeq, list(tmp$L1), c)
  temp$Peptides <- tmp$x[match(temp$`1st accession`, tmp$Group.1)]
  names(temp$"Sequence (1st accession)") <- temp$"1st accession"
  tmpLst[[exp]] <- temp[, c("Sequence (1st accession)", "Peptides")]
}
tmp <- do.call(rbind, tmpLst)
tmpFl <- tempfile(fileext = ".rds")
clusterExport(parClust, "tmpFl", envir = environment())
readr::write_rds(tmp, tmpFl)
invisible(clusterCall(parClust, \(x) {
  assign("tmp", readr::read_rds(tmpFl), envir = .GlobalEnv)
  return()
}))
unlink(tmpFl)
clusterExport(parClust, "Coverage", envir = environment())
tmp$Coverage <- parSapply(parClust, seq_len(nrow(tmp)), \(x) {
  round(100*Coverage(tmp$"Sequence (1st accession)"[x], tmp$"Peptides"[[x]]), 1L)
})
nRws <- lengths(tmpRws)
stopifnot(nrow(tmp) == sum(nRws))
for (exp in Exp) {
  rg <- seq_len(nRws[exp])
  m <- match(exp, Exp)
  if (m > 1L) {
    rg <- rg + sum(nRws[1L:(m-1L)])
  }
  PG[[paste0("Sequence coverage [%] - ", exp)]] <- 0
  PG[tmpRws[[exp]], paste0("Sequence coverage [%] - ", exp)] <- tmp$Coverage[rg]
}

#### Code chunk - Calculate protein group-level quantitative values
quntSrc <- paste0(libPath, "/extdata/Sources/PG_Quant.R")
#rstudioapi::documentOpen(quntSrc)
source(quntSrc, local = FALSE)

#### Code chunk - Re-normalize protein group expression values
# Normalize (Levenberg-Marquardt)
if ((length(Exp) > 1L)&&(NormalizePG)) {
  g <- grep(topattern(PG.int.col), colnames(PG), value = TRUE)
  temp <- PG[, c("id", g)]
  m <- apply(temp[, g], 2L, \(x) { median(is.all.good(x)) })
  M <- median(is.all.good(unlist(temp[, g])))
  temp[, g] <- sweep(temp[, g], 2L, m, "-") + M
  temp <- AdvNorm.IL(temp[, c("id", g)], "id", exprs.col = g, exprs.log = TRUE)
  PG[, gsub(topattern(PG.int.col), paste0("Norm. ", PG.int.cols["Original"]), g)] <- temp[, paste0("AdvNorm.", g)]
  PG.int.cols["Normalized"] <- PG.int.col <- paste0("Norm. ", PG.int.cols["Original"])
  if (MakeRatios) {
    for (gr in unique(SamplesMap$`Ratios group`)) { #gr <- unique(SamplesMap$`Ratios group`)[1L]
      m <- SamplesMap[which(SamplesMap$`Ratios group` == gr),]
      rf <- m$Experiment[which(m$Reference)]
      for (i in m$Experiment[which(!m$Reference)]) { #i <- m$Experiment[which(!m$Reference)][1L]
        PG[[paste0("Norm. ", rat.cols["Original"], " - ", i)]] <- (PG[[paste0(PG.int.col, i)]] - PG[[paste0(PG.int.col, rf)]])/log10(2L)
      }
    }
    PG.rat.cols["Normalized"] <- PG.rat.col <- paste0("Norm. ", PG.rat.cols["Original"])
  }
}

# Average expression columns:
if (length(Exp) > 1L) {
  for (i in PG.int.cols) {
    PG[[paste0("Mean ", gsub(" - $", "", i))]] <- apply(PG[, grep(topattern(i), colnames(PG), value = TRUE)],
                                                        1L, mean, na.rm = TRUE)
  }
}

#### Code chunk - Prepare Annotations and (if applicable) GO terms
# NB: I used to get functional annotations for all proteins in the protein group.
#     However we are now - and I think with reason - only using annotations from the leading protein(s)!
setwd(wd)
Src <- paste0(libPath, "/extdata/Sources/Annotate_me.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# GO filters
if (GO_filt) {
  for (goID in GO_filter) { #goID <- GO_filter[1L]
    # Get children terms
    gofilter <- unique(unlist(c(goID,
                                GOBPOFFSPRING[[goID]],
                                GOCCOFFSPRING[[goID]],
                                GOMFOFFSPRING[[goID]])))
    gofilter <- gofilter[which(!is.na(gofilter))]
    if (sum(gofilter %in% AllTerms)) {
      PG[[goID]] <- ""
      PG[grsep(gofilter, x = PG$`GO-ID`), goID] <- "+"
      #which(vapply(strsplit(PG$`GO-ID`[-wtst],";"), \(x) { sum(x %in% gofilter) }, 1L) != 0L)
      #which(vapply(strsplit(PG$`GO-ID`[wtst],";"), \(x) { sum(x %in% gofilter) }, 1L == 0L)
    }
  }
  Src <- paste0(libPath, "/extdata/Sources/interestGO_Fisher.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

#### Code chunk - Correlation and distribution plots
# LFQ correlation plots:
# Idea for this plot: map color to protein type (list, organism or contaminant) as in profile plots
dir <- paste0(wd, "/Workflow control")
if (length(Exp) > 1L) {
  comb <- as.data.frame(gtools::combinations(length(Exp), 2L, Exp))
  temp <- setNames(lapply(names(PG.int.cols), \(klnm) { #klnm <- names(PG.int.cols)[1L]
    kol <- PG.int.cols[klnm]
    kolZ <- grep(topattern(kol), colnames(PG), value = TRUE)
    if (!length(kolZ)) { return() }
    temp2 <- PG[, kolZ] 
    source(parSrc, local = FALSE)
    tmpFl <- tempfile(fileext = ".rds")
    clusterExport(parClust, list("kol", "klnm", "tmpFl"), envir = environment())
    readr::write_rds(temp2, tmpFl)
    clusterCall(parClust, \(x) {
      assign("temp2", readr::read_rds(tmpFl), envir = .GlobalEnv)
      return()
    })
    unlink(tmpFl)
    temp2 <- parApply(parClust, comb, 1L, \(x) {
      temp3 <- temp2[, paste0(kol, unlist(x))]
      temp3$X <- x[[1L]]
      temp3$Y <- x[[2L]]
      temp3$Comparison <- paste0(x[[1L]], " (X) vs ", x[[2L]], " (Y)")
      colnames(temp3)[1L:2L] <- c("log10(X LFQ)", "log10(Y LFQ)")
      temp3$Type <- klnm
      return(temp3)
    })
    temp2 <- plyr::rbind.fill(temp2)
    temp2$"Common Name (short)" <- PG$"Common Name (short)"
    temp2$Type <- klnm
    return(temp2)
  }), names(PG.int.cols))
  temp <- plyr::rbind.fill(temp)
  clusterExport(parClust, "is.all.good", envir = environment())
  test <- parApply(parClust, temp[, c("log10(X LFQ)", "log10(Y LFQ)")], 1L, \(x) {
    length(is.all.good(x))
  }) == 2L
  temp <- temp[which(test),]
  w <- aggregate(1L:nrow(temp), list(temp$Type), list)
  temp2 <- data.frame(Comparison = rep(unique(temp$Comparison), length(PG.int.cols)))
  temp2$Type <- as.character(sapply(names(PG.int.cols), \(x) { rep(x, length(unique(temp$Comparison))) }))
  temp2[, c("Median", "SD", "R")] <- as.data.frame(t(apply(temp2[, c("Type", "Comparison")], 1L, \(x) {
    #x <- temp2[1L, c("Type", "Comparison")]
    x1 <- temp[which((temp$Type == x[[1L]])&(temp$Comparison == x[[2L]])), c("log10(X LFQ)", "log10(Y LFQ)")]
    x2 <- x1$"log10(Y LFQ)"-x1$"log10(X LFQ)"
    return(c(paste0("Median = ", round(median(x2), 3L)),
             paste0("S.D. = ", round(sd(x2), 3L)),
             paste0("R = ", round(cor(x1$"log10(X LFQ)", x1$"log10(Y LFQ)"), 3L))))
  })))
  temp$Comparison <- gsub(" vs ", "\nvs\n", temp$Comparison)
  temp2$Comparison <- gsub(" vs ", "\nvs\n", temp2$Comparison)
  temp$Type <- factor(temp$Type, levels = names(PG.int.cols))
  temp2$Type <- factor(temp2$Type, levels = names(PG.int.cols))
  temp$X <- factor(temp$X, levels = Exp)
  temp$Y <- factor(temp$Y, levels = Exp)
  temp2[, c("X", "Y")] <- temp[match(temp2$Comparison, temp$Comparison), c("X", "Y")]
  x_min <- min(temp$`log10(X LFQ)`)
  y_min <- min(temp$`log10(Y LFQ)`)
  y_max <- max(temp$`log10(Y LFQ)`)
  tmp <- temp[, c("log10(X LFQ)", "log10(Y LFQ)", "Comparison", "Common Name (short)")]
  tmpFl <- tempfile(fileext = ".rds")
  clusterExport(parClust, list("get_density", "tmpFl"), envir = environment())
  readr::write_rds(tmp, tmpFl)
  clusterCall(parClust, \(x) {
    assign("tmp", readr::read_rds(tmpFl), envir = .GlobalEnv)
    return()
  })
  unlink(tmpFl)
  comps <- unique(temp$Comparison)
  tmpD <- setNames(parLapply(parClust, comps, \(cmp) { #cmp <- comps[1L]
    w <- which(tmp$Comparison == cmp)
    x <- get_density(tmp$`log10(X LFQ)`[w], tmp$`log10(Y LFQ)`[w], n = 500L)
    x <- data.frame(Density = x,
                    Nm = tmp$"Common Name (short)"[w])
    return(x)
  }), comps)
  temp$Intensity <- 0
  for (cmp in comps) { #cmp <- comps[1L]
    w <- which(temp$Comparison == cmp)
    tmpDw <- tmpD[[cmp]]
    temp$Intensity[w] <- tmpDw$Density[match(temp$"Common Name (short)"[w], tmpDw$Nm)]
  }
  ttl <- "Correlation plot - PGs level"
  if (prot.list.Cond) {
    temp$"In list" <- PG$`In list`[match(temp$`Common Name (short)`, PG$`Common Name (short)`)]
  }
  plot <- ggplot(temp) + geom_point(aes(x = `log10(X LFQ)`, y = `log10(Y LFQ)`, color = Intensity), size = 0.1) +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_text(data = temp2, x = x_min, y = y_max - 0.01*(y_max - y_min), aes(label = R), hjust = 0, cex = 2.5) +
    geom_text(data = temp2, x = x_min, y = y_max - 0.06*(y_max - y_min), aes(label = Median), hjust = 0, cex = 2.5) +
    geom_text(data = temp2, x = x_min, y = y_max - 0.11*(y_max - y_min), aes(label = SD), hjust = 0, cex = 2.5) +
    facet_grid(Y~Type+X) + ggtitle(ttl) + coord_fixed(1) +
    scale_color_viridis() +
    theme_bw() + theme(strip.text.y = element_text(angle = -90, vjust = 0, hjust = 0.5),
                       strip.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5))
  if (prot.list.Cond) {
    plot <- plot + geom_point(data = temp[which(temp$`In list` == "+"),], aes(x = `log10(X LFQ)`, y = `log10(Y LFQ)`), size = 1L, shape = 1L, color = "red") +
      geom_text(data = temp[which(temp$`In list` == "+"),], aes(label = `Common Name (short)`, x = `log10(X LFQ)`, y = `log10(Y LFQ)`),
                size = 2L, hjust = 0, vjust = 0, color = "red")
  }
  poplot(plot, 12L, 22L)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
  })
}
g <- as.character(sapply(PG.int.cols, \(x) { grep(topattern(x), colnames(PG), value = TRUE) }))
test <- vapply(g, \(x) { length(is.all.good(PG[[x]])) }, 1L)
print(test)
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))

# Intensities distribution:
long.dat <- list()
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
temp <- PG[, c("Common Name (short)", as.character(sapply(PG.int.cols, \(x) { paste0(x, Exp) })))]
temp <- reshape::melt(temp, id.vars = "Common Name (short)")
colnames(temp)[which(colnames(temp) == "value")] <- "log10(Intensity)"
if (length(Exp) > 1L) {
  temp2 <- PG[, c("Common Name (short)", paste0("Mean ", gsub(" - $", "", PG.int.cols)))]
  temp2 <- dfMelt(temp2, id.vars = "Common Name (short)")
  temp$"Mean log10(Intensity)" <- temp2$value
}
temp$Experiment <- gsub(paste0(".*", topattern(PG.int.cols["Original"], start = FALSE)), "", temp$variable)
temp$Experiment <- factor(temp$Experiment, levels = SamplesMap$Experiment)
temp$Type <- gsub(paste0(" ?", topattern(PG.int.cols["Original"], start = FALSE), ".*$"), "", temp$variable)
temp$Type[which(temp$Type == "")] <- "Orig."
temp$Type <- paste0("log10(", tolower(temp$Type), " LFQ)")
temp$Type <- factor(temp$Type, levels = paste0("log10(", c("orig", "imput", "norm"), ". LFQ)"))
long.dat$intens <- temp
temp <- temp[which(is.all.good(temp$"log10(Intensity)", 2L)),]
ttl <- "LFQ density plot - PGs level"
if (prot.list.Cond) {
  temp$"In list" <- PG$"In list"[match(temp$`Common Name (short)`, PG$`Common Name (short)`)]
  plot <- ggplot(temp) + geom_histogram(aes(x = `log10(Intensity)`, fill = `In list`), bins = 100L)
} else {
  plot <- ggplot(temp) + geom_histogram(aes(x = `log10(Intensity)`, fill = Type), bins = 100L)
}
plot <- plot +
  theme_bw() + theme(#axis.line = element_blank(),
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    #axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    #legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank()) +
  ggtitle(ttl) + facet_grid(Type~Experiment) +
  scale_y_continuous(expand = c(0L, 0L)) +
  scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
  xlab("log10(Protein Groups LFQ)")
#poplot(plot)
M <- c()
if (GO_filt) {
  CompGOTerms %<o%setNames(c("GO:0005634", "GO:0005654", "GO:0000785", "GO:0005730", "GO:0005635", "GO:0005737", "GO:0005829", "GO:0005783", "GO:0005794", "GO:0031988", "GO:0005739", "GO:0009536", "GO:0005886", "GO:0031012", "GO:1903561"),
                           c("Nucleus", "Nucleoplasm", "Chromatin", "Nucleolus", "Nuclear envelope", "Cytoplasm", "Cytosol", "ER", "Golgi", "Vesicle", "Mitochondrion", "Plastid", "Plasma membrane", "Extracellular matrix", "Extracellular vesicle"))
  GO_filter1 %<o% unique(c(GO_filter, CompGOTerms))
  for (goID in GO_filter1) { #goID <- GO_filter1[1L]
    # Get children terms
    gofilter <- unique(unlist(c(goID,
                                GOBPOFFSPRING[[goID]],
                                GOCCOFFSPRING[[goID]],
                                GOMFOFFSPRING[[goID]])))
    gofilter <- gofilter[which(!is.na(gofilter))]
    if (sum(gofilter %in% AllTerms)) {
      temp[[goID]] <- 0L
      wtst <- grsep2(gofilter, PG$`GO-ID`)
      #which(vapply(strsplit(PG$`GO-ID`[-wtst],";"), \(x) { sum(x %in% gofilter) }, 1L) != 0L)
      #which(vapply(strsplit(PG$`GO-ID`[wtst],";"), \(x) { sum(x %in% gofilter) }, 1L) == 0L)
      wtst2 <- which(temp$"Common Name (short)" %in% PG$"Common Name (short)"[wtst])
      temp[wtst2, goID] <- 1L
      w <- which(temp$Type == "log10(orig. LFQ)")
      wtst3 <- which(temp$"Common Name (short)"[w] %in% PG$"Common Name (short)"[wtst])
      if (length(wtst3) > 1L) {
        d <- density(temp[wtst3, "log10(Intensity)"])
        M[goID] <- d$x[which.max(d$y)]
      }
    }
  }
}
if (length(M)) {
  leg2 <- get_legend(plot)
  plot2 <- plot + theme(legend.position = "none")
  #
  temp2 <- temp[which(temp$Type == "log10(orig. LFQ)"),
                c("log10(Intensity)", "Experiment", names(M))]
  temp2 <- dfMelt(temp2, id.vars = c("log10(Intensity)", "Experiment"))
  colnames(temp2) <- c("log10(Intensity)", "Experiment", "GO term", "+")
  temp2 <- temp2[which(temp2$"+" == 1L),]
  if (globalGO) {
    temp2$`GO term` <- GO_terms$Term[match(temp2$`GO term`, GO_terms$ID)]
    temp2$`GO term` <- gsub("\\]$", "", gsub(" \\[", "\n",  temp2$`GO term`))
  }
  ttl2 <- "LFQ density plot - PGs level, GO terms"
  plot3 <- ggplot(temp2) + geom_density(stat = "density", aes(x = `log10(Intensity)`, fill = `GO term`),
                                        colour = "black", alpha = 0.5) +
    geom_density(data = temp, stat = "density", aes(x = `log10(Intensity)`), colour = "red",
                 linewidth = 1L, linetype = "dotted", show.legend = FALSE) +
    scale_y_continuous(expand = c(0L, 0L)) +
    scale_fill_viridis(option = "C", discrete = TRUE) + theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.background = element_blank(),
          strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.y = element_blank(),
          axis.title.y = element_blank(),
          plot.subtitle = element_text(face = "italic", color = "red"))
  if (globalGO) { plot3 <- plot3 + facet_grid(`GO term`~Experiment)  } else { plot2 <- plot2 + facet_wrap(~Experiment) }
  plot3a <- plot3 + ggtitle(ttl2, subtitle = "(dotted line = sample's reference density)")
  plot3 <- plot3 + ggtitle("", subtitle = "(dotted line = sample's reference density)")
  #poplot(plot2)
  #poplot(plot3)
  #poplot(plot3a)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl2, ".jpeg"), plot3a, dpi = 300L, width = 10L, height = 10L, units = "in")
    ggsave(paste0(dir, "/", ttl2, ".pdf"), plot3a, dpi = 300L, width = 10L, height = 10L, units = "in")
  })
  leg3 <- get_legend(plot3)
  plot3 <- plot3 + theme(legend.position = "none")
  g2 <- ggplotGrob(plot2)
  g3 <- ggplotGrob(plot3)
  maxWidth <- grid::unit.pmax(g2$widths[2L:5L], g3$widths[2L:5L])
  g2$widths[2L:5L] <- as.list(maxWidth)
  g3$widths[2L:5L] <- as.list(maxWidth)
  plot4 <- arrangeGrob(grobs = list(g2, leg2, g3, leg3),
                       widths = c(4L, 1L),
                       heights = c(length(unique(temp$Type)), length(M)/5),
                       layout_matrix = rbind(c(1L, 2L), c(3L, 4L)),
                       padding = 5L)
  plot4 <- as.ggplot(plot4)
  poplot(plot4)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot4, dpi = 300L, width = 10L, height = 10L, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot4, dpi = 300L, width = 10L, height = 10L, units = "in")
  })
  #poplot(plot3)
} else {
  poplot(plot, 12L, 22L)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
  })
}

# Gene-Set Enrichment Analysis (GSEA)
if (runGSEA) {
  dataType <- "PG"
  GSEAmode <- "standard"
  Src <- paste0(libPath, "/extdata/Sources/GSEA.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}
#

if (MakeRatios) {
  FC_filt %<o% c()
  FC_Smpls %<o% list()
  # Fold change filters:
  ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1L]
  rat.grps <- unique(SamplesMap$`Ratios group`)
  rat.grps <- rat.grps[which(!is.na(rat.grps))]
  for (grp in rat.grps) { #grp <- rat.grps[1L]
    SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
    smpl0 <- SmplMp$Experiment[which(SmplMp$Reference)]
    stopifnot(length(smpl0) == 1L)
    smpl1 <- SmplMp$Experiment[which(!SmplMp$Reference)]
    e0 <- PG[[paste0(ref, smpl0)]]
    if (NegFilt) {
      nf <- SmplMp$Experiment[which(SmplMp$"Negative Filter")]
      if (length(nf)) {
        smpl1 <- smpl1[which(!smpl1 %in% nf)]
        nftst <- apply(PG[, paste0(ref, nf), drop = FALSE], 1L, \(x) { length(is.all.good(x)) }) > 0L
      }
    }
    FC_filt <- append(FC_filt, setNames(lapply(smpl1, \(x) {
      e1 <- PG[[paste0(ref, x)]]
      r1 <- PG[[paste0(PG.rat.col, x)]]
      if (RatiosThresh_2sided) { r1 <- abs(r1) }
      w <- which(((r1 >= RatiosThresh)&(is.all.good(e1, 2L)))|(is.all.good(e1, 2L)&(!is.all.good(e0, 2L))))
      if ((NegFilt)&&(length(nf))) { w <- w[which(!nftst[w])] }
      return(w)
    }), smpl1))
    FC_Smpls[[grp]] <- list(Numerator = smpl1, Denominator = smpl0)
  }
  #g <- grep(topattern(PG.rat.col), colnames(PG), value = TRUE)
  #test <- setNames(vapply(g, \(x) { length(is.all.good(PG[[x]])) }, 1L), gsub(topattern(PG.rat.col), "", g))
  #print(test)
  #
  # Ratios distribution:
  dir <- paste0(wd, "/Workflow control")
  temp <- PG[, c("Common Name (short)", as.character(unlist(sapply(PG.rat.cols, \(x) {
    grep(topattern(x), colnames(PG), value = TRUE)
  }))))]
  temp <- reshape::melt(temp, id.vars = "Common Name (short)")
  colnames(temp)[which(colnames(temp) == "value")] <- "log2(Ratio)"
  temp$Experiment <- gsub(paste0(".*", topattern(PG.rat.cols["Original"], start = FALSE)), "", temp$variable)
  temp$Experiment <- factor(temp$Experiment, levels = SamplesMap$Experiment)
  temp$Type <- gsub(paste0(" ?", topattern(PG.rat.cols["Original"], start = FALSE), ".*$"), "", temp$variable)
  temp$Type[which(temp$Type == "")] <- "Orig."
  temp$Type <- paste0("log2(", tolower(temp$Type), " ratio)")
  temp$Type <- factor(temp$Type, levels = paste0("log2(", c("orig", "norm", "imput"), ". ratio)"))
  long.dat$ratios <- temp
  temp <- temp[which(is.all.good(temp$"log2(Ratio)", 2L)),]
  ttl <- "Ratios density plot - PGs level"
  plot <- ggplot(temp) + geom_histogram(aes(x = `log2(Ratio)`, fill = Type), bins = 100L) +
    geom_vline(xintercept = RatiosThresh, colour = "red") +
    theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
    scale_y_continuous(expand = c(0L, 0L)) +
    scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
    ggtitle(ttl) + facet_grid(Type~Experiment) +
    xlab("log2(Ratio)")
  if (RatiosThresh_2sided) { plot <- plot + geom_vline(xintercept = -RatiosThresh, colour = "red") }
  poplot(plot)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
  })
  # MA plots:
  temp <- long.dat$intens
  temp <- temp[which(temp$Experiment %in% long.dat$ratios$Experiment),]
  tst1 <- do.call(paste, c(temp[, c("Common Name (short)", "Experiment")], sep = "___"))
  tst2 <- do.call(paste, c(long.dat$ratios[, c("Common Name (short)", "Experiment")], sep = "___"))
  temp$"log2(Ratio)" <- long.dat$ratios$`log2(Ratio)`[match(tst1, tst2)]
  temp <- temp[which(is.all.good(temp$`log2(Ratio)`, 2L)),]
  ttl <- "MA plots - PGs level"
  plot <- ggplot(temp) + geom_point(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`, colour = Type), size = 0.1) +
    geom_hline(yintercept = 0, linewidth = 0.8, linetype = "dashed") +
    ggtitle(ttl) + coord_fixed(log10(2L)/log2(10L)) + theme_bw() + facet_grid(Type~Experiment) +
    scale_color_viridis(option = "D", discrete = TRUE, begin = 0.25) +
    xlab("A = mean log10(Intensity)") + ylab("M = sample log2(Ratio)") +
    theme(strip.text.y = element_text(angle = 0))
  poplot(plot)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
  })
  # "Regulated/Enriched" columns
  ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1L]
  for (grp in rat.grps) { #grp <- rat.grps[1L]
    SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
    if (sum(c(TRUE, FALSE) %in% SmplMp$Reference) < 2L) {
      stop("The reference column should include TRUE and FALSE values!")
    } else {
      w0 <- which(SmplMp$Reference)
      stopifnot(length(w0) == 1L)
      w1 <- which(!SmplMp$Reference)
      xp0 <- SmplMp$Experiment[w0]
      for (x in w1) { #x <- w1[1L]
        xp1 <- SmplMp$Experiment[x]
        wonly <- which(is.all.good(PG[[paste0(ref, xp1)]], 2L)&(!is.all.good(PG[[paste0(ref, xp0)]], 2L)))
        PG[[paste0("Regulated - ", xp1)]] <- ""
        if (RatiosThresh_2sided) {
          wup <- which(PG[[paste0(PG.rat.col, xp1)]] >= RatiosThresh)
          wdwn <- which(PG[[paste0(PG.rat.col, xp1)]] <= -RatiosThresh)
          PG[wup, paste0("Regulated - ", xp1)] <- "up"
          PG[wdwn, paste0("Regulated - ", xp1)] <- "down"
        } else {
          PG[[paste0("Enriched - ", xp1)]] <- c("", "+")[(PG[[paste0(PG.rat.col, xp1)]] >= RatiosThresh)+1L]
        }
        PG[wonly, paste0("Regulated - ", xp1)] <- "specific"
      }
    }
  }
  # Modified peptides
  if (PTMriched) {
    PTMs_FC_filt %<o% list()
    PTMs_FC_Smpls %<o% list()
    PTM_normalize %<o% list()
    PTMs_pep %<o% list()
    PTMs_intRf %<o% rev(int.cols[which(int.cols != paste0("Imput. ", int.cols["Original"]))])[1L]
    PTMs_ratRf %<o% rev(rat.cols[which(rat.cols != paste0("Imput. ", rat.cols["Original"]))])[1L]
    PTMs_intNm0 <- names(PTMs_intRf)
    PTMs_ratNm0 <- names(PTMs_ratRf)
    for (ptm in EnrichedPTMs) { #ptm <- EnrichedPTMs[1L]
      p <- Modifs$Mark[match(ptm, Modifs$`Full name`)]
      ppat <- paste0("\\(", p, "\\)|\\(", p, ",|,", p, "\\)|,", p, ",") # Pattern to catch all instances of the mod
      PTMs_FC_filt[[ptm]] <- c()
      PTMs_FC_Smpls[[ptm]] <- list()
      PTM_normalize[[ptm]] <- TRUE
      ptmpep <- pep[which(pep[[ptm]]),]
      a <- unlist(strsplit(gsub("\\)$", "", ptm), "\\("))
      if (length(a) > 1L) {
        Ptm <- paste0(toupper(substr(a[1L], 1L, 1L)), substr(a[1L], 2L, nchar(ptm)), "(", a[2L], ")")
      } else {
        Ptm <- paste0(toupper(substr(a[1L], 1L, 1L)), substr(a[1L], 2L, nchar(ptm)))
      }
      p <- Modifs$Mark[match(ptm, Modifs$`Full name`)]
      ptmsh <- substr(p, 1L, 1L)
      temp <- ptmpep[, c("Modified sequence", "Leading proteins")]
      temp$"Leading proteins" <- strsplit(temp$"Leading proteins", ";")
      temp$"Modified sequence" <- gsub(paste0("[^A-Z", ptmsh, "]"), "",
                                       gsub(ppat, ptmsh, temp$"Modified sequence"))
      ptmpep[, c("Match(es)", paste0(Ptm, "-site(s)"))] <- ""
      dbsmall <- db[which(db$"Protein ID" %in% unique(unlist(temp$"Leading proteins"))), c("Protein ID", "Sequence")]
      # On I/L ambiguity remaining even with newer DIA methods taking into account RT, IM and fragments intensity,
      # see https://github.com/vdemichev/DiaNN/discussions/1631
      dbsmall$"Seq*" <- gsub("I", "L", dbsmall$Sequence)
      temp$"ModSeq*" <- gsub("I", "L", temp$`Modified sequence`)
      #
      kol <- c("Leading proteins", "ModSeq*")
      temp$`ModSeq*` <- strsplit(temp$`ModSeq*`, "")
      ptmpep[, c("Match(es)", paste0(Ptm, "-site(s)"))] <- as.data.frame(t(apply(temp[, kol], 1L, \(x) {
        m <- unlist(x[[2L]])
        m <- data.frame(Seq = m, Mod.seq = m, Test = FALSE)
        w1 <- which(m$Seq == ptmsh)
        w2 <- which(m$Seq != ptmsh)
        m$Mod.seq[w1-1L] <- paste0(ptmsh, m$Mod.seq[w1-1L])
        m$Test[w1-1L] <- TRUE
        m <- m[w2,]
        l <- nrow(m)
        m$Offset <- 0L:(l-1L)
        q <- unlist(x[[1L]])
        mtch <- match(q, dbsmall$"Protein ID")
        wN <- which(!is.na(mtch))
        mtch <- mtch[wN]
        q <- dbsmall$"Protein ID"[mtch[wN]]
        if (length(mtch)) {
          seq <- strsplit(dbsmall$"Seq*"[mtch], "")
          matches <- lapply(seq, \(S) { #S <- seq[1L]
            S <- unlist(S)
            lS <- length(S)
            m1 <- m
            m1$Match <- apply(m1[, c("Seq", "Offset")], 1L, \(y) { which(S == y[1L]) - as.numeric(y[2L]) })
            M <- unlist(m1$Match)
            M <- M[which(M > 0L)]
            M <- aggregate(M, list(M), length)
            M <- M[order(-M$x),]
            M <- M$Group.1[which(M$x == l)]
            # Check that peptides are tryptic:
            #test <- sapply(M, \(y) {
            #  # r1: on the N-terminal end, is the peptide preceded by K, R or (if starting at position 2, M)?
            #  if (y > 1L) {
            #    if (y == 2) { r1 <- S[y-1] %in% c("K", "R", "M") } else { r1 <- S[y-1] %in% c("K", "R") }
            #  } else { r1 <- TRUE }
            #  # r2: on the C-terminal end, is this a tryptic peptide or the last peptide in the protein?
            #  r2 <- (m1$Seq[l] %in% c("K", "R"))|(y+l-1 == lS)
            #  return(r1+r2 == 2)
            #})
            #M <- M[which(test)]
            return(M)
          })
          names(matches) <- q
          matches <- matches[which(lengths(matches) > 0L)]
          if (length(matches)) {
            matches <- set_colnames(reshape::melt(matches), c("Match", "Protein"))
            matches <- aggregate(matches$Protein, list(matches$Match), paste, collapse = ";")
            colnames(matches) <- c("Match", "Proteins")
            w <- which(m$Test)
            matches$Sites <- sapply(matches$Match, \(y) {
              y <- paste(sapply(w, \(z) { paste0(m$Mod.seq[z], y+m$Offset[z]) }), collapse = "-")
            })
            matches$Match <- apply(matches[, c("Match", "Proteins")], 1L, paste, collapse = " ")
            matches$Sites <- apply(matches[, c("Sites", "Proteins")], 1L, paste, collapse = " ")
            matches <- apply(matches[, c("Match", "Sites")], 2L, paste, collapse = "/")
          } else { matches <- c(NA, NA) }
        } else { matches <- c(NA, NA) }
        return(matches)
      })))
      ptmpep[[paste0(Ptm, "-site")]] <- gsub(" .+", "", ptmpep[[paste0(Ptm, "-site(s)")]])
      ptmpep <- ptmpep[which(!is.na(ptmpep$`Match(es)`)),]
      ptmpep$tmp1 <- gsub("^_|_$", "", ptmpep$`Modified sequence`)
      ptmpep$tmp2 <- ptmpep[[paste0(Ptm, "-site(s)")]]
      nc <- nchar(ptmpep$tmp2)
      w <- which(nc > 25L)
      ptmpep$tmp2[w] <- paste0(substr(ptmpep$tmp2[w], 1L, 22L), "...")
      ptmpep$Code <- apply(ptmpep[, paste0("tmp", as.character(1L:2L))], 1L, paste, collapse = "\n")
      ptmpep$tmp1 <- NULL
      ptmpep$tmp2 <- NULL
      ptmpep$Name <- ""
      w <- which(ptmpep$"Leading proteins" != "")
      ptmpep$Name[w] <- vapply(strsplit(gsub("[/,;].+$", "", ptmpep[w, paste0(Ptm, "-site(s)")]), " "), \(x) {
        paste0(x[[1L]], " ", db$"Common Name"[match(x[[2L]], db$"Protein ID")])
      }, "")
      ptmpep$Name[which(ptmpep$Name == "")] <- paste0("Unknown ", ptm, "-modified peptide #", seq_along(which(ptmpep$Name == "")))
      #View(ptmpep[,c("Match(es)", "Modified sequence", "Code", paste0(Ptm, "-site(s)"))])
      if (grepl("^[P,p]hospho( \\([A-Z]+\\))?$", ptm)) {
        p_col <- paste0(gsub(" |\\(|\\)", ".", ptm), ".Probabilities")
        scd_col <- paste0(gsub(" |\\(|\\)", ".", ptm), ".Score.Diffs")
        if (sum(c(p_col, scd_col) %in% colnames(ptmpep)) == 2L) {
          temp <- try(phos_QC(ptmpep, P_col = p_col, ScD_col = scd_col), silent = TRUE)
          if (inherits(temp, "try-error")) {
            warning("No phospho QC performed, check colnames or the phos_QC function!")
          } else { ptmpep$High_Quality_ptmpep <- temp }
        }
      }
      #
      # Optional: normalize to parent protein group
      if (PTM_normalize[[ptm]]) {
        # Step 1: normalize ratios:
        a <- grep(topattern(PTMs_ratRf[PTMs_ratNm0]), colnames(ptmpep), value = TRUE)
        a1 <- gsub(topattern(paste0(PTMs_ratRf[PTMs_ratNm0], " - ")), PG.rat.col, a)
        stopifnot(sum(!a1 %in% colnames(PG)) == 0L)
        temp <- Isapply(strsplit(ptmpep$"Protein group IDs", ";"), \(x) { #x <- strsplit(ptmpep$"Protein group IDs", ";")[1L]
          x <- unlist(x)
          y <- PG[match(x, PG$id), a1, drop = FALSE]
          if (length(x) > 1L) { y <- apply(y, 2L, \(x) { mean(is.all.good(x)) }) }
          return(unlist(y))
        })
        ptmpep[, paste0("ReNorm. ", a)] <- ptmpep[, a] - temp # It's log data so "-", not "/"
        PTMs_ratRf["Re-normalized"] <- paste0("ReNorm. ", PTMs_ratRf[PTMs_ratNm0])
        # Step 2: adjust expression values to reflect the re-normalized ratios:
        # Intensities are not logged (for peptides) so we need to delog ratios and work in multiplicative, not additive, mode 
        PTMs_intRf["Re-normalized"] <- "ReNorm. int."
        temp <- set_colnames(data.frame(matrix(rep(0L, nrow(ptmpep)*length(Exp)),
                                               ncol = length(Exp))),
                             paste0(PTMs_intRf["Re-normalized"], " - ", Exp))
        for (grp in rat.grps) { #grp <- rat.grps[1L]
          e <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
          smpls0 <- e$Experiment[which(e$Reference %in% c("TRUE", TRUE))]
          stopifnot(length(smpls0) == 1L) # This is a script without replicates! Only one reference is allowed per sample group!
          smpls1 <- e$Experiment[which(!e$Reference %in% c("TRUE", TRUE))]
          # Pre-normalisation values
          cole0 <- paste0(PTMs_intRf[PTMs_intNm0], " - ", smpls0)
          cole1 <- paste0(PTMs_intRf[PTMs_intNm0], " - ", smpls1)
          # Re-normalized intensities (what we are calculating) 
          cols0 <- paste0(PTMs_intRf["Re-normalized"], " - ", smpls0)
          cols1 <- paste0(PTMs_intRf["Re-normalized"], " - ", smpls1)
          # The re-normalised ratios we use for that // these are log2-transformed!
          colr0 <- paste0(PTMs_ratRf["Re-normalized"], " - ", smpls0)
          colr0 <- colr0[which(colr0 %in% colnames(ptmpep))]
          colr1 <- paste0(PTMs_ratRf["Re-normalized"], " - ", smpls1)
          #
          totB <- apply(ptmpep[, c(cole0, cole1)], 1L, sum, na.rm = TRUE)
          if (length(colr0)) {
            av <- apply(ptmpep[, cole0, drop = FALSE], 1L, mean, na.rm = TRUE)
            temp[, c(cols0, cols1)] <- sweep(2L^ptmpep[, c(colr0, colr1)], 1L, av, "*")
            #tst <- log10(temp[, cols1]/temp[, cols0])
            #tst2 <- apply(tst, 2L, \(x) { summary(is.all.good(x)) })
          } else {
            temp[, cols0] <- ptmpep[, cole0] # This stays the same as before
            temp[, cols1] <- ptmpep[, cole0]*(2L^ptmpep[, colr1])
          }
          #tst <- log10(temp[, cols1]/temp[, cols0])
          #tst2 <- apply(tst, 2L, \(x) { summary(is.all.good(x)) })
          # Note: the price of normalisation is that often there is no valid parent protein so a lot of NAs are introduced
          #
          totA <- apply(temp[, c(cols0, cols1)], 1L, sum, na.rm = TRUE)
          temp[, c(cols0, cols1)] <- temp[, c(cols0, cols1)]*totB/totA # Re-apply correct scale
        }
        ptmpep[, colnames(temp)] <- temp
        #kol <- grep("ReNorm. log2", colnames(ptmpep), value = TRUE)
        #tst <- apply(ptmpep[, kol], 2L, \(x) { summary(is.all.good(x)) })
      }
      #
      # Fold change filters:
      Int <- PTMs_intRf[length(PTMs_intRf)]
      Rat <- PTMs_ratRf[length(PTMs_ratRf)]
      for (grp in rat.grps) { #grp <- rat.grps[1L]
        SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
        smpl0 <- SmplMp$Experiment[which(SmplMp$Reference)]
        stopifnot(length(smpl0) == 1L)
        smpl1 <- SmplMp$Experiment[which(!SmplMp$Reference)]
        e0 <- ptmpep[[paste0(Int, " - ", smpl0)]]
        if (NegFilt) {
          nf <- SmplMp$Experiment[which(SmplMp$"Negative Filter")]
          if (length(nf)) {
            smpl1 <- smpl1[which(!smpl1 %in% nf)]
            nftst <- apply(ptmpep[, paste0(Int, " - ", nf), drop = FALSE], 1L, \(x) {
              length(is.all.good(x))
            }) > 0L
          }
        }
        PTMs_FC_filt[[ptm]] <- append(PTMs_FC_filt[[ptm]], setNames(lapply(smpl1, \(x) {
          e1 <- ptmpep[[paste0(Int, " - ", x)]]
          r1 <- ptmpep[[paste0(Rat, " - ", x)]]
          if (RatiosThresh_2sided) { r1 <- abs(r1) }
          w <- which(((r1 >= RatiosThresh)&(is.all.good(e1, 2L)))|(is.all.good(e1, 2L)&(!is.all.good(e0, 2L))))
          if ((NegFilt)&&(length(nf))) { w <- w[which(!nftst[w])] }
          return(w)
        }), smpl1))
        PTMs_FC_Smpls[[ptm]][[grp]] <- list(Numerator = smpl1, Denominator = smpl0)
      }
      #
      # "Regulated/Enriched" columns
      Int <- PTMs_intRf[length(PTMs_intRf)]
      Rat <- PTMs_ratRf[length(PTMs_ratRf)]
      for (grp in rat.grps) { #grp <- rat.grps[1L]
        SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
        if (sum(c(TRUE, FALSE) %in% SmplMp$Reference) < 2L) {
          stop("The reference column should include TRUE and FALSE values!")
        } else {
          w0 <- which(SmplMp$Reference)
          stopifnot(length(w0) == 1L)
          w1 <- which(!SmplMp$Reference)
          xp0 <- SmplMp$Experiment[w0]
          for (x in w1) { #x <- w1[1L]
            xp1 <- SmplMp$Experiment[x]
            wonly <- which(is.all.good(ptmpep[[paste0(Int, " - ", xp1)]], 2L)&(!is.all.good(ptmpep[[paste0(Int, " - ", xp0)]], 2L)))
            PG[[paste0("Regulated - ", xp1)]] <- ""
            if (RatiosThresh_2sided) {
              wup <- which(ptmpep[[paste0(Rat, " - ", xp1)]] >= RatiosThresh)
              wdwn <- which(ptmpep[[paste0(Rat, " - ", xp1)]] <= -RatiosThresh)
              ptmpep[wup, paste0("Regulated - ", xp1)] <- "up"
              ptmpep[wdwn, paste0("Regulated - ", xp1)] <- "down"
            } else {
              ptmpep[[paste0("Enriched - ", xp1)]] <- c("", "+")[(ptmpep[[paste0(Rat, " - ", xp1)]] >= RatiosThresh)+1L]
            }
            ptmpep[wonly, paste0("Regulated - ", xp1)] <- "specific"
          }
        }
      }
      # Ratios distribution:
      dir <- paste0(wd, "/Workflow control")
      temp <- ptmpep[, c("Code", as.character(sapply(paste0(PTMs_ratRf, " - "), \(x) {
        grep(topattern(x), colnames(ptmpep), value = TRUE)
      })))]
      temp <- dfMelt(temp, id.vars = "Code")
      colnames(temp)[which(colnames(temp) == "value")] <- "log2(Ratio)"
      temp$Experiment <- gsub(paste0(".*", topattern(paste0(PTMs_ratRf["Original"], " - "), start = FALSE)), "", temp$variable)
      temp$Experiment <- factor(temp$Experiment, levels = SamplesMap$Experiment)
      temp$Type <- gsub(paste0(" ?", topattern(PTMs_ratRf["Original"], start = FALSE), ".*$"), "", temp$variable)
      temp$Type[which(temp$Type == "")] <- "Orig."
      temp$Type <- paste0("log2(", tolower(temp$Type), " ratio)")
      temp$Type <- factor(temp$Type, levels = paste0("log2(", c("orig", "norm", "imput", "renorm"), ". ratio)"))
      long.dat$ratios <- temp
      temp <- temp[which(is.all.good(temp$"log2(Ratio)", 2L)),]
      ttl <- paste0("Ratios density plot - ", ptm, "-modified peptides")
      plot <- ggplot(temp) + geom_histogram(aes(x = `log2(Ratio)`, fill = Type), bins = 100L) +
        geom_vline(xintercept = RatiosThresh, colour = "red") +
        theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
        scale_y_continuous(expand = c(0L, 0L)) +
        scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
        ggtitle(ttl) + facet_grid(Type~Experiment) +
        xlab("log2(Ratio)")
      if (RatiosThresh_2sided) { plot <- plot + geom_vline(xintercept = -RatiosThresh, colour = "red") }
      poplot(plot)
      suppressMessages({
        ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
        ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
      })
      #
      # Gene-Set Enrichment Analysis (GSEA)
      if (runGSEA) {
        dataType <- "modPeptides"
        Src <- paste0(libPath, "/extdata/Sources/GSEA.R")
        #rstudioapi::documentOpen(Src)
        source(Src, local = FALSE)
      }
      #
      PTMs_pep[[ptm]] <- ptmpep
    }
  }
}

# Biotinylated peptides
if (IsBioID2) {
  for (exp in Exp) { #exp <- Exp[1L]
    PG[[paste0("Biot. peptides count - ", exp)]] <- 0L
    e <- ev[which(ev$Experiment == exp),]
    g <- grep(topattern(Modifs$Mark[wbiot], start = FALSE), e$"Modified sequence")
    if (length(g)) {
      e <- e[g,]
      temp2 <- listMelt(strsplit(e$"Protein group IDs", ";"), e$"Modified sequence")
      temp2 <- aggregate(temp2$L1, list(temp2$value), \(x) { length(unique(x)) })
      w2 <- which(PG$id %in% temp2$Group.1)
      PG[w2, paste0("Biot. peptides count - ", exp)] <- temp2$x[match(PG$id[w2], temp2$Group.1)]
    }
  }
}

#### Code chunk - Proteomic ruler
if (protrul) {
  ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1L]
  if (length(Exp) > 1L) { ref <- c(ref, paste0("Mean ", gsub(" - $", "", ref))) }
  temp <- try(Prot.Ruler(PG, db, ref, NuclL = ProtRulNuclL), silent = TRUE)
  if ((!inherits(temp, "try-error"))&&(!is.logical(temp))) {
    PG <- temp$Protein.groups
    db <- temp$Database
  } else { protrul <- FALSE }
}

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Summary table and QC plots
source(parSrc, local = FALSE)
Exp_summary %<o% MQ.summary(wd = wd, ev = ev, pg = PG, mods = setNames(Modifs$Mark, Modifs$"Full name"),
                            raw.files = rawFiles, sc = max(c(20L, round(length(rawFiles2)/length(Exp)))),
                            save = c("jpeg", "pdf"), cl = parClust,
                            MQtxt = inDirs[which(SearchSoft == "MAXQUANT")])
write.csv(Exp_summary, paste0(wd, "/Workflow control/Summary.csv"), row.names = FALSE)
#Exp_summary <- read.csv(paste0(wd, "/Workflow control/Summary.csv"), check.names = FALSE)

# Plot of contamination levels per sample
tmp <- ev[, c("Proteins", "Intensity", "Experiment")]
tmp2 <- listMelt(strsplit(tmp$Proteins, ";"), 1L:nrow(tmp), c("Protein", "Row"))
m <- match(gsub("^CON__", "", tmp2$Protein), gsub("^CON__", "", db$`Protein ID`))
w <- which(!is.na(m))
tmp2 <- tmp2[w,]; m <- m[w]
# For our purpose here we must match contaminant proteins.
tmp2$Cont <- db$`Potential contaminant`[m]
if (tstOrg) {
  tmp2$Organism <- db[m, dbOrgKol]
  tmp2 <- as.data.table(tmp2)
  f0 <- function(x) { c("Target", "Contaminant")[("Contaminant" %in% x)+1L] }
  tmp2 <- tmp2[, list(x = f0(Organism)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
} else {
  tmp2 <- as.data.table(tmp2)
  f0 <- function(x) { c("Target", "Contaminant")[("+" %in% x)+1L] }
  tmp2 <- tmp2[, list(x = f0(Cont)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
}
tmp$Organism <- tmp2$x[match(1L:nrow(tmp), tmp2$Group.1)]
tmp <- tmp[which(!is.na(tmp$Organism)),]
tmp$Organism <- factor(tmp$Organism, levels = c("Contaminant", "Target"))
tmp <- aggregate(tmp$Intensity, list(tmp$Experiment, tmp$Organism), sum)
colnames(tmp) <- c("Experiment", "Organism", "Total intensity")
tmp$Experiment <- factor(tmp$Experiment, levels = SamplesMap$Experiment)
ttl <- "Contributions to TIC"
plot <- ggplot(tmp) + geom_bar(stat = "identity", aes(x = Experiment, y = `Total intensity`, fill = Organism)) +
  theme_bw() + scale_fill_viridis(discrete = TRUE, begin = 0.8, end = 0.2) +
  ggtitle(ttl, subtitle = "Summed TIC for each class of identified peptides")
poplot(plot)
suppressMessages({
  ggsave(paste0(wd, "/Summary plots/", ttl, ".jpeg"), plot, dpi = 150L, width = 10L, height = 10L, units = "in")
  ggsave(paste0(wd, "/Summary plots/", ttl, ".pdf"), plot, dpi = 150L, width = 10L, height = 10L, units = "in")
})

# Test for amino acid biases:
Src <- paste0(libPath, "/extdata/Sources/AA_biases_test.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - Coverage maps for proteins of interest
Src <- paste0(libPath, "/extdata/Sources/noRep_protPlots.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

### Check that CytoScape is installed and can run, then launch it.
#CytoScape <- TRUE #You may need to reset this
Src <- paste0(libPath, "/extdata/Sources/Cytoscape_init.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Initialize ClueGO
if (runClueGO) {
  Src <- paste0(libPath, "/extdata/Sources/ClueGO_init.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

#### Code chunk - Gene Ontology terms enrichment analysis
source(parSrc, local = FALSE)
create_plotly %<o% TRUE
goSrc <- paste0(libPath, "/extdata/Sources/noRep_GO.R")
#rstudioapi::documentOpen(goSrc)
source(goSrc, local = FALSE)

#### Code chunk - Venn diagrams
setwd(wd)
packs <- c("ggplot2", "ggpolypath", "venn")
for (pck in packs) {
  if (!require(pck, character.only = TRUE)) {
    pak::pak(pck, upgrade = FALSE, ask = FALSE)
  }
}
for (pck in packs) {
  library(pck, character.only = TRUE)
}
HdrStlVenn <- createStyle(textDecoration = "bold", halign = "left", valign = "bottom", wrapText = TRUE,
                          numFmt = "TEXT", fontSize = 12L, textRotation = 60)
wb <- createWorkbook()
wbKount <- 0L
VennMx <- 7L
if (Venn_Obs) {
  dir <- paste0(wd, "/Venn diagrams")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  ref <- PG.int.cols["Original"]
  Grps <- unique(SamplesMap$`Ratios group`)
  for (grp in Grps) {
    sm <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
    Xp <- sm$Experiment
    test <- setNames(lapply(Xp, \(exp) {
      x <- PG[[paste0(ref, exp)]]
      which(is.all.good(x, 2L))
    }), Xp)
    w <- which(lengths(test) > 0L)
    VennExp <- Xp[w]
    OK <- length(w) > 1L
    if (OK) {
      if (length(w) > VennMx) {
        msg <- paste0("Too many samples, select at least 2 and up to ", VennMx,
                      " to include in the Venn diagram (comma-separated):")
        if (length(Grps) > 1L) { msg <- paste0("Ratios group ", grp, ": ", msg) }
        tst <- !sm$Reference[match(VennExp, sm$Experiment)]
        tmp <- setNames(vapply(VennExp, \(x) { paste(c(x, rep(" ", 200L-nchar(x))), collapse = "") }, ""), VennExp)
        VennExp <- names(tmp)[match(dlg_list(tmp, tmp[which(tst)[1L:min(c(sum(tst, VennMx)))]], TRUE, msg)$res, tmp)]
        if (length(VennExp) < 2L) {
          msg <- "Skipping per-sample observations Venn diagrams"
          if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected at least 2 samples!")
          warning(msg)
          OK <- FALSE
        }
        if (length(VennExp) > VennMx) {
          msg <- "Skipping per-sample observations Venn diagrams"
          if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected ", VennMx,
                        " samples at most!")
          warning(msg)
          OK <- FALSE
        }
      }
    }
    if (OK) {
      ttl <- "Observations_LFQ_Venn_diagram_-_global"
      SheetNm <- "Samples composition"
      if (length(Grps) > 1L) {
        ttl <- paste0(ttl, " - group ", grp)
        SheetNm <- paste0(SheetNm, "_grp", grp)
      }
      test <- test[VennExp]
      AnalysisParam$"Venn diagram (LFC) - samples" <- list(VennExp)
      plot <- venn(test, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
      plot <- plot + ggtitle("Venn diagram", subtitle = "Global, LFQ") +
        theme(plot.title = element_text(size = 15L), plot.subtitle = element_text(size = 10L))
      poplot(plot)
      suppressMessages({
        ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 150L)
        ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150L)
      })
      #system(paste0("open \"", dir, "/", ttl, ".jpg", "\""))
      wbKount <- wbKount+1L
      if (SheetNm %in% names(wb)) { removeWorksheet(wb, SheetNm) } 
      addWorksheet(wb, SheetNm)
      writeData(wb, SheetNm, PG[, c("id", "Leading protein IDs", "Genes")], 1L, 1L)
      l <- length(test)
      tmp <- sapply(names(test), \(smpl) {
        res <- rep("", nrow(PG))
        res[test[[smpl]]] <- "+"
        return(res)
      })
      writeData(wb, SheetNm, tmp, 4L, 1L)
      setRowHeights(wb, SheetNm, 1L, 120L)
      addStyle(wb, SheetNm, HdrStlVenn, 1L, 1L:(l+3L))
    } else { warning("Skipping LFC Venn diagrams: not enough valid samples!") }
  }
}
if (Venn_Ratios) {
  VennTypes <- ""
  if (RatiosThresh_2sided) { VennTypes <- c(VennTypes, "up", "down") }
  Grps <- unique(SamplesMap$`Ratios group`)
  for (grp in Grps) {
    sm <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
    VennExp <- sm$Experiment[which(!sm$Reference)]
    fc_filt <- FC_filt[VennExp]
    tst <- lengths(fc_filt)
    w <- which(tst == 0L)
    if (0L %in% tst) {
      l <- length(w)
      tmp <- names(fc_filt)[which(tst == 0L)]
      if (l > 1L) { tmp <- paste0(paste(tmp[1L:(l-1L)], collapse = ", "), " and ", tmp[l]) }
      warning(paste0("No filtered proteins for samples ", tmp, ", they will be skipped!"))
    }
    fc_filt <- fc_filt[which(tst > 0L)]
    VennExp <- names(fc_filt)
    OK <- length(VennExp) > 1L
    if (OK) {
      if (length(fc_filt) > VennMx) {
        msg <- paste0("Too many samples, select at least 2 and up to ", VennMx,
                      " to include in the Venn diagram (comma-separated):")
        if (length(Grps) > 1L) { msg <- paste0("Ratios group ", grp, ": ", msg) }
        VennExp <- dlg_list(VennExp, VennExp[1L:VennMx], TRUE, msg)$res
        if (length(VennExp) < 2L) {
          msg <- "Skipping ratios Venn diagrams"
          if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected at least 2 samples!")
          warning(msg)
          OK <- FALSE
        }
        if (length(VennExp) > VennMx) {
          msg <- "Skipping ratios Venn diagrams"
          if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected ", VennMx,
                        " samples at most!")
          warning(msg)
          OK <- FALSE
        }
      }
      if (OK) {
        fc_filt <- fc_filt[VennExp]
        AnalysisParam$"Venn diagram (ratios) - samples" <- list(VennExp)
        if (RatiosThresh_2sided) {
          updowntst <- setNames(lapply(VennExp, \(x) {
            #x <- VennExp[1L]
            rs <- sign(PG[fc_filt[[x]], paste0(PG.rat.cols["Original"], x)])
            w <- which(is.na(rs))
            rs[w] <- c(1L, -1L)[is.na(PG[fc_filt[[x]][w], paste0(PG.int.cols["Original"], x)])+1L]
            return(rs)
          }), VennExp)
        }
        setwd(wd); suppressWarnings(dir.create("Venn diagrams"))
        for (vt in VennTypes) { #vt <- ""
          ttl <- gsub("_\\(\\)$", "", paste0("Ratios_Venn_diagram_-_global", "_(", vt, ")"))
          SheetNm <- paste0(c("Up/down", "Up", "Down")[match(vt, VennTypes)], "-reg. PGs")
          if (length(Grps) > 1L) {
            msg <- paste0(msg, " for ratios group ", grp)
            ttl <- paste0(ttl, " - group ", grp)
            SheetNm <- paste0(SheetNm, "_grp", grp)
          }
          sbttl <- rat.col
          flt <- fc_filt
          if (vt %in% c("up", "down")) {
            sbttl <- paste0(sbttl, ", ", vt)
            flt <- setNames(lapply(VennExp, \(x) {
              flt[[x]][which(updowntst[[x]] == c(1L, -1L)[match(vt, c("up", "down"))])]
            }), VennExp)
          }
          w <- which(lengths(flt) > 0L)
          if (length(w) > 1L) {
            flt <- flt[w]
            plot <- venn(flt, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
            plot <- plot + ggtitle("Venn diagram", subtitle = sbttl) +
              theme(plot.title = element_text(size = 15L), plot.subtitle = element_text(size = 10L))
            poplot(plot)
            suppressMessages({
              ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 150L)
              ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150L)
            })
            #system(paste0("open \"", dir, "/", ttl, ".jpg", "\""))
            wbKount <- wbKount+1L
            if (SheetNm %in% names(wb)) { removeWorksheet(wb, SheetNm) } 
            addWorksheet(wb, SheetNm)
            writeData(wb, SheetNm, PG[, c("id", "Leading protein IDs", "Genes")], 1L, 1L)
            l <- length(flt)
            tmp <- sapply(names(flt), \(smpl) {
              res <- rep("", nrow(PG))
              res[flt[[smpl]]] <- "+"
              return(res)
            })
            writeData(wb, SheetNm, tmp, 4L, 1L)
            setRowHeights(wb, SheetNm, 1L, 120L)
            addStyle(wb, SheetNm, HdrStlVenn, 1L, 1L:(l+3L))
          } else { message(gsub(" \\(\\)$", "", paste0("No overlaps to plot (", vt, ")"))) }
        }
      }
    } else {
      msg <- "Skipping ratios Venn diagrams"
      if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
      msg <- paste0(msg, ": not enough valid samples!")
      warning(msg)
    }
  }
}
if (wbKount) { saveWorkbook(wb, paste0(wd, "/Venn diagrams/Venn diagrams.xlsx"), overwrite = TRUE) }
setwd(wd)

#### Code chunk - Dimensionality reduction plots
## Process data
if (length(Exp) > 2L) {
  setwd(wd); suppressWarnings(dir.create("PCA plots"))
  temp <- PG[, grep(topattern(PG.int.col), colnames(PG), value = TRUE)]
  colnames(temp) <- gsub(topattern(PG.int.col), "", colnames(temp))
  ## Impute data
  ## Here we have no way to decide between MAR/MCAR/MNAR,
  ## so we will instead replace every missing value with a random value drawn from a gaussian distribution of reduced m and sd
  m <- median(is.all.good(unlist(temp)))
  sd <- sd(is.all.good(unlist(temp)))
  for (i in colnames(temp)) {
    w <- which(!is.all.good(temp[[i]], 2L))
    temp2 <- rnorm(length(w), m-3, sd/2)
    temp[w, i] <- temp2
  }
  nrm <- PG[[paste0("Mean ", gsub(" - $", "", PG.int.col))]]
  w <- which(is.all.good(nrm, 2L))
  temp <- sweep(temp[w,], 1L, nrm[w], "-")
  rownames(temp) <- PG$"Leading protein IDs"[w]
  temp <- temp + rnorm(nrow(temp)*ncol(temp), 0, 10L^-24L) # Add small random value in case we (very rarely) get non-unique values per row
  ## 1/ Samples level:
  pc <- prcomp(t(temp), scale. = TRUE)
  scores <- as.data.frame(pc$x)
  pv <- round(100*(pc$sdev)^2L / sum(pc$sdev^2L), 0L)
  pv <- pv[which(pv > 0)]
  pv <- paste0("Components: ", paste(vapply(seq_along(pv), \(x) {
    paste0("PC", x, ": ", pv[x], "%")
  }, ""), collapse = ", "))
  scores$Sample <- rownames(scores)
  rownames(scores) <- NULL
  tst <- apply(scores[, c("PC1", "PC2")], 2L, \(x) { length(unique(x)) })
  if (min(tst) > 1L) {
    ttl <- "PCA plot - Samples-level"
    plot <- ggplot(scores) + geom_point(aes(x = PC1, y = PC2, colour = Sample)) +
      coord_fixed() + ggtitle(ttl, subtitle = pv) +
      geom_text_repel(aes(x = PC1, y = PC2, label = Sample, colour = Sample), size = 2.5) + 
      theme_bw() + theme(legend.position = "none")
    #poplot(plot)
    if ("PC3" %in% colnames(scores)) {
      plot_lyPCA <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Sample, text = ~Sample, type = "scatter3d",
                            mode = "markers")
      plot_lyPCA <- add_trace(plot_lyPCA, scores, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "text",
                              showlegend = FALSE)
      plot_lyPCA <- layout(plot_lyPCA, title = ttl)
      tst <- try(saveWidget(partial_bundle(plot_lyPCA), paste0(wd, "/PCA plots/", ttl, ".html")), silent = TRUE)
      if (inherits(tst, "try-error")) { tst <- try(saveWidget(plot_lyPCA, paste0(wd, "/PCA plots/", ttl, ".html")), silent = TRUE) }
      if (!inherits(tst, "try-error")) { system(paste0("open \"", wd, "/PCA plots/", ttl, ".html")) }
    } else { poplot(plot, width = 18L) }
    suppressMessages({
      ggsave(paste0(wd, "/PCA plots/", ttl, ".jpeg"), plot, dpi = 150L)
      ggsave(paste0(wd, "/PCA plots/", ttl, ".pdf"), plot, dpi = 150L)
    })
    ## 2/ Protein groups level:
    pc <- prcomp(temp, scale. = TRUE)
    scores <- as.data.frame(pc$x)
    pv <- round(100*(pc$sdev)^2L / sum(pc$sdev^2L), 0L)
    pv <- pv[which(pv > 0)]
    pv <- paste0("Components: ", paste(vapply(seq_along(pv), \(x) {
      paste0("PC", x, ": ", pv[x], "%")
    }, ""), collapse = ", "))
    scores$"Leading protein IDs" <- rownames(scores)
    rownames(scores) <- NULL
    scores[, c("Protein IDs", "Protein group")] <- PG[match(scores$"Leading protein IDs", PG$"Leading protein IDs"),
                                                      c("Protein IDs", "Label")]
    scores$Alpha <- (scores$PC1^2L + scores$PC2^2L)
    scores$Direction <- apply(temp[w,], 1L, \(x) {
      wh <- which(is.all.good(10L^x, 2L))
      return(weighted.mean(c(seq_along(colnames(temp)))[wh], 10L^x[wh]))
    })
    scores$Class <- ""
    breaks <- seq_along(Exp)
    labels <- Exp
    if (prot.list.Cond) {
      g1 <- grsep2(prot.list, scores$"Protein IDs")
      if (length(g1)) {
        g2 <- grsep2(prot.list, scores$"Protein IDs", invert = TRUE)
        scores2 <- scores[g1,]
        scores <- scores[g2,]
        scores2$Alpha <- 1L
        scores2$Class <- length(Exp)+1L
      }
    }
    ttl <- "PCA plot - Protein groups (PG-level)"
    plot <- ggplot(scores) + geom_point(aes(x = PC1, y = PC2, alpha = Alpha, colour = Direction, text = `Protein group`), shape = 1) +
      coord_fixed() + ggtitle(ttl, subtitle = pv) +
      theme_bw() + guides(alpha = "none", shape = "none") +
      scale_color_gradientn(colors = hcl.colors(length(Exp)), breaks = breaks, labels = labels, guide = "legend")
    if ((prot.list.Cond)&&(length(g1)))
      plot <- plot + geom_point(data = scores2, colour = "red", shape = 2L, aes(x = PC1, y = PC2, text = `Protein group`))
  }
  suppressMessages({
    ggsave(paste0(wd, "/PCA plots/", ttl, ".jpeg"), plot, dpi = 150L)
    ggsave(paste0(wd, "/PCA plots/", ttl, ".pdf"), plot, dpi = 150L)
  })
  if ("PC3" %in% colnames(scores)) {
    plot_lyPCAProt <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Direction, text = ~`Protein group`,
                              type = "scatter3d", mode = "markers", showlegend = FALSE, marker = list(size = 1L))
    plot_lyPCAProt <- layout(plot_lyPCAProt, title = ttl)
    if ((prot.list.Cond)&&(length(g1))) {
      scores3 <- scores2
      scores3$"Protein group" <- gsub(" - | ?, ?", "<br>", scores3$"Protein group")
      plot_lyPCAProt <- add_trace(plot_lyPCAProt, data = scores3, x = ~PC1, y = ~PC2, z = ~PC3,
                                  type = "scatter3d", mode = "markers+text", color = I("red"), marker = list(size = 5L),
                                  showlegend = FALSE, textposition = "bottom right")
    }
  } else { plot_lyPCAProt <- ggplotly(plot, tooltip = "text") }
  tst <- try(saveWidget(partial_bundle(plot_lyPCAProt), paste0(wd, "/PCA plots/", ttl, ".html")), silent = TRUE)
  if (inherits(tst, "try-error")) { tst <- try(saveWidget(plot_lyPCAProt, paste0(wd, "/PCA plots/", ttl, ".html")), silent = TRUE) }
  if (!inherits(tst, "try-error")) { system(paste0("open \"", wd, "/PCA plots/", ttl, ".html")) }
  #plot <- plot + geom_text_repel(data = scores, aes(x = PC1, y = PC2, label = `Protein group`, alpha = Alpha), size = 2.5)
  if ((prot.list.Cond)&&(length(g1))) {
    plot <- plot + geom_text_repel(data = scores2, colour = "red", size = 2.5,
                                   aes(x = PC1, y = PC2, label = `Protein group`))
  }
  #poplot(plot, width = 18L)
  suppressMessages({
    ggsave(paste0(wd, "/PCA plots/", ttl, " (labels).jpeg"), plot, dpi = 150L)
    ggsave(paste0(wd, "/PCA plots/", ttl, " (labels).pdf"), plot, dpi = 150L)
  })
} else { warning("No PCA plots drawn: samples are too similar!") }

# Prepare data for clustering and dimensionality reduction plots
Src <- paste0(libPath, "/extdata/Sources/cluster_Heatmap_Prep.R")
#rstudioapi::documentOpen(Src)
dataType <- "PG"
source(Src, local = FALSE)
dataType <- "peptides"
source(Src, local = FALSE)

#### Code chunk - Heatmaps with clustering at samples and protein groups level, highlighting proteins of interest
clustMode <- "standard"
dataType <- "PG"
clstSrc <- paste0(libPath, "/extdata/Sources/cluster_Heatmap_Main.R")
#rstudioapi::documentOpen(clstSrc)
source(clstSrc, local = FALSE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - profile plots and ranked abundance plots
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

# Negative filter
if (NegFilt) {
  e <- ev[grep("-MATCH$", ev$Type, invert = TRUE),]
  nf <- SamplesMap$Experiment[which(SamplesMap$"Negative Filter")]
  e <- e[which(e$Experiment %in% nf),]
  PG$"Direct identification in negative filter sample(s)" <- vapply(strsplit(PG$`Evidence IDs`, ";"), \(x) {
    sum(unlist(x) %in% e$id)
  }, 1L) > 0L
  PG$"Direct identification in negative filter sample(s)" <- c("", "+")[match(PG$"Direct identification in negative filter sample(s)", c(FALSE, TRUE))]
  if (MakeRatios) {
    exp <- SamplesMap$Experiment[which(!SamplesMap$Reference)]
    PG[which(PG$"Direct identification in negative filter sample(s)" == "+"),
       paste0(c("Enriched", "Regulated")[RatiosThresh_2sided+1L], " - ", exp)] <- ""
  }
}

# Backup data/update cluster
stopClust <- TRUE
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

#### Code chunk - XML coverage columns
Src <- paste0(libPath, "/extdata/Sources/xml_Coverage_columns.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
# Calculate maximum expected coverage per group
source(parSrc, local = FALSE)
if (WorkFlow == "Band ID") {
  m <- match(frstProt, db$`Protein ID`)
  Dig <- data.frame(ID = frstProt, Seq = db$Sequence[m])
  invisible(clusterCall(parClust, \(x) {
    library(proteoCraft)
    return()
  }))
  Dig$Digest <- Digest(Dig$Seq, min = MinPepSz, missed = Missed, cl = parClust)
  Dig$MaxCov <- parApply(parClust, Dig[, c("Seq", "Digest")], 1L, \(x) {
    Coverage(x[[1L]], x[[2L]])
  })
  PG$`Max. theoretical sequence coverage [%]` <- round(100*Dig$MaxCov, 1L)
}

#### Code chunk - Create output tables
## PSMs
dir <- paste0(wd, "/Tables")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
w <- which(vapply(colnames(ev), \(x) { inherits(ev[[x]], "list") }, TRUE))
if (length(w)) { for (i in w) { ev[[i]] <- parSapply(parClust, ev[[i]], paste, collapse = ";") } }
data.table::fwrite(ev, paste0(dir, "/evidence.tsv"), sep = "\t", row.names = FALSE, na = "NA")
#
## Main peptidoforms- and protein groups-level, multi-tabs report
xlSrc <- paste0(libPath, "/extdata/Sources/noRep_Write_Excel.R")
#rstudioapi::documentOpen(xlSrc)
source(xlSrc, local = FALSE)
#xl_open(repFl)

# Save special quantitative table for proteins of interest
if ((length(Exp) > 1L)&&(prot.list.Cond)) {
  dir <- paste0(wd, "/Proteins of interest")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  w <- grsep2(prot.list, PG$"Leading protein IDs")
  if (length(w)) {
    kols <- grep(topattern(PG.int.col), colnames(PG), value = TRUE)
    temp <- PG[w, c("Label", "Peptide IDs", kols)]
    for (kol in kols) {
      kol2 <- gsub(" log10\\(", " ", gsub("\\) - ", " - ", kol))
      temp[[kol2]] <- suppressWarnings(10L^temp[[kol]])
      temp[which(!is.all.good(temp[[kol2]], 2L)), kol] <- NA
      temp[[kol]] <- NULL
    }
    data.table::fwrite(temp, paste0(dir, "/Protein of interest profiles.tsv"),
                       sep = "\t", row.names = FALSE, na = "NA")
  }
}

#### Code chunk - Heatmap(s) of peptides mapping to proteins of interest
if ((length(Exp) > 1L)&&(!is.null(prot.list))&&(length(prot.list))) {
  dir <- paste0(wd, "/Heatmaps")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  ref <- int.cols[which(names(int.cols) == "Imputed")-1L]
  kol <- paste0(ref, " - ", Exp)
  kol <- kol[which(kol %in% colnames(pep))]
  StdWdth <- 6L
  for (plp in prot.list) { #plp <- prot.list[1L]
    Plp <- paste(db[which(db$"Protein ID" == plp), c("Common Name", "Protein ID")], collapse = " - ")
    grs <- grsep2(plp, pep$Proteins)
    if (length(grs)) {
      temp <- pep[grs, c("Sequence", "Modified sequence", kol)]
      Seq <- db$Sequence[match(plp, db$`Protein ID`)]
      temp$tst1 <- vapply(temp$Sequence, \(x) { nchar(unlist(strsplit(Seq, x))[1L]) }, 1L)
      temp$tst2 <- vapply(temp$Sequence, nchar, 1L)
      temp$tst3 <- vapply(temp$"Modified sequence", nchar, 1L)
      temp <- temp[order(temp$tst1, temp$tst2, temp$tst3),]
      w <- which(rowSums(temp[, kol, drop = FALSE], na.rm = TRUE) == 0)
      if (length(w)) {
        warning("Peptide(s) found with all invalid intensity values!")
        w <- which(rowSums(temp[, kol, drop = FALSE], na.rm = TRUE) > 0)
        temp <- temp[w,]
      }
      if (nrow(temp)) {
        colnames(temp) <- gsub(topattern(paste0(ref, " - ")), "", colnames(temp))
        # Create heatmap
        temp <- temp[, c("Modified sequence", Exp)]
        tst <- apply(temp[, Exp, drop = FALSE], 1L, \(x) { mean(x[which(x > 0)]) })
        temp[, Exp] <- sweep(temp[, Exp, drop = FALSE], 1L, tst, "/")
        temp2 <- set_colnames(dfMelt(temp, id.vars = "Modified sequence"),
                              c("Modified sequence", "Sample", "value"))
        temp2$"Modified sequence" <- gsub("^_|_$", "", temp2$"Modified sequence")
        temp2$Sample <- as.character(temp2$Sample)
        temp2$value <- suppressWarnings(log2(temp2$value))+StdWdth/2
        w <- which(!is.all.good(temp2$value, 2L))
        temp2$value[w] <- NA
        temp2$value[which(temp2$value < -StdWdth/2)] <- -StdWdth/2
        temp2$Xmin <- match(temp2$Sample, colnames(temp))-1L
        temp2$Xmax <- temp2$Xmin+1L
        temp2$Ymax <- nrow(temp):1L
        temp2$Ymin <- temp2$Ymax-1L
        hlab <- aggregate(temp2[, c("Xmin", "Xmax")], list(temp2$Sample), unique)
        colnames(hlab)[1L] <- "Sample"
        hlab$X <- rowMeans(hlab[, c("Xmin", "Xmax")])
        vlab <- aggregate(temp2[, c("Ymin", "Ymax")], list(temp2$"Modified sequence"), unique)
        colnames(vlab)[1L] <- "Modified sequence"
        vlab$Y <- rowMeans(vlab[, c("Ymin", "Ymax")])
        Xscale <- max(temp2$Xmax)
        Yscale <- max(temp2$Ymax)
        # Create heatmap plot
        Xlim <- c(-1L, Xscale+15L)
        Ylim <- c(-6L, Yscale+20L)
        temp2a <- temp2[, c("Xmin", "Ymin", "value")]
        Splits <- 20L
        XScale2 <- Xscale*0.1/Splits
        temp2b <- data.frame(Xmin = c(Xscale/2 - 8*XScale2*((-Splits/2L):(Splits/2L)), 4),
                             Ymin = Yscale+10L)
        temp2b$Xmax <- temp2b$Xmin + XScale2*8
        temp2b$Xmax[nrow(temp2b)] <- temp2b$Xmin[nrow(temp2b)]+0.5
        temp2b$value <- c((Splits:0L)*StdWdth/Splits - StdWdth/2, NA)
        temp2b$Label <- ""
        temp2b$Label[c(1L, Splits/2L+1L, Splits+1L, Splits+2L)] <- c(paste0("log10(max)", c("", -StdWdth/2, -StdWdth)), "0 or NA")
        # Create graph
        ttl <- paste0("Peptides log2 heatmap - ", gsub("/", "-", Plp))
        heatmap.plot <- ggplot(temp2a) +
          geom_rect(aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
          geom_rect(data = temp2b, aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymin+1, fill = value)) +
          geom_text(data = hlab, aes(x = X, label = Sample), y = Yscale+1, angle = 60, hjust = 0, size = 2.5) +
          geom_text(data = vlab, aes(y = Y, label = `Modified sequence`), x = Xscale+0.5, hjust = 0, size = 1.7) +
          geom_text(data = temp2b, aes(x = Xmin, label = Label), y = Yscale+11.25, angle = 60, hjust = 0, vjust = 0.5, size = 2.5) +
          ggtitle("Peptides heatmap (log2 scale, normalised)", subtitle = Plp) +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                panel.background = element_rect(fill = "transparent", color = NA),
                plot.margin = margin(0L, 0L, 0L, 0L, "cm")) +
          coord_fixed(0.5) +
          scale_fill_gradient2(low = "#D55E00", mid = "black", high = "#009E73", na.value = "lightblue", breaks = c(-StdWdth/2, 0, StdWdth/2)) +
          xlab(NULL) + ylab(NULL) + theme(legend.position = "none") +
          xlim(Xlim[1L], Xlim[2L]) + ylim(Ylim[1L], Ylim[2L])
        poplot(heatmap.plot)
        suppressMessages({
          ggsave(paste0(dir, "/", ttl, ".jpeg"), heatmap.plot,
                 dpi = 600L, width = 20L, height = 12L, units = "in")
          ggsave(paste0(dir, "/", ttl, ".pdf"), heatmap.plot,
                 dpi = 600L, width = 20L, height = 12L, units = "in")
        })
        #system(paste0("open \"", dir, "/", ttl, ".jpeg", "\""))
        #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
      }
    }
  }
  setwd(wd)
}

#### Code chunk - peptide tables for visualizing the coverage of proteins of interest in 3D
Src <- paste0(libPath, "/extdata/Sources/3d_Cov.R")
#rstudioapi::documentOpen(xlSrc)
source(Src)

# Save parameters
temp <- data.frame(Parameter = names(AnalysisParam),
                   Value = vapply(AnalysisParam, \(x) { paste(unlist(x), collapse = ";") }, ""),
                   row.names = NULL)
write.csv(temp, "Analysis parameters.csv", row.names = FALSE)
saveFun(AnalysisParam, "AnalysisParam.RData")
#
# Save copy of this script to local work directory
fs::file_copy(ScriptPath, wd, overwrite = TRUE)

# Backup data/update cluster
#rstudioapi::documentOpen(bckpSrc)
source(bckpSrc, local = FALSE)
#loadFun(BckUpFl)

# Write PTMs table
temp <- Modifs
w <- which(vapply(colnames(Modifs), \(x) { inherits(Modifs[[x]], "list") }, TRUE))
for (i in w) { temp[[i]] <- vapply(temp[[i]],  paste, "", collapse = ", ") }
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
write.csv(temp, paste0(dir, "/Modifications.csv"), row.names = FALSE)

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
#loadFun(BckUpFl)
source(Src, local = FALSE)

### That's it, done!
#openwd(outDir)
#rm(list = ls())
