#### Code chunk - Initialization
if (!interactive()) { stop("This script should only be run within an interactive R session!") }
options(stringsAsFactors = FALSE)
options(install.packages.compile.from.source = "never")
#rm(list = ls()[which(!ls() %in% c("dtstNm", "wd", "indir", "outdir"))])

## The proteoCraft package can be re-installed at any time in the workflow (there is a specific script for this in the package's library folder),
## or just load it here:
if (exists(".obj")) { rm(".obj") }
myPackNm %<o% "proteoCraft"
library(myPackNm, character.only = FALSE)
dirlist %<o% c() # This should go!!!
ReUseAnsw %<o% FALSE
ReLoadPSMsBckp %<o% FALSE

RPath %<o% as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match(myPackNm, RPath$Package)], winslash = "/")
libPath %<o% paste0(RPath, "/", myPackNm)
homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/", myPackNm)
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1, 1)) }
parSrc %<o% paste0(libPath, "/extdata/R scripts/Sources/make_check_Cluster.R")
fls <- paste0(homePath, "/", c("Regulation analysis - master script.R",
                               "Regulation analysis - detailed script.R",
                               "Regulation analysis - detailed script_pepOnly.R",
                               "No replicates analysis - detailed script.R",
                               "Reload_renv_from_lock_file.R",
                               "Default_locations.xlsx",
                               "LC_columns.xlsx"))
tst <- sum(!file.exists(fls))
if (tst) { eval(parse(text = paste0(myPackNm, "::Configure()"))) }
scrptType %<o% "withReps"
scrptTypeFull %<o% "withReps_PG_and_PTMs"

# Parameters used by the master script:
###-|-### Workflows: setNames(c("Differential Protein Expression analysis", "Pull-Down (e.g. co-IP)", "Biotin-based Pull-Down (BioID, TurboID, APEX...)", "Time Course","SubCellular Localisation analysis"), c("REGULATION", "PULLDOWN", "BIOID", "TIMECOURSE", "LOCALISATION"))
###-|-### Replicates? TRUE
###-|-### External dependencies: Excel (loose); ScanHeadsman (loose); Cytoscape (loose); saintExpress (auto)

### Packages
## For convenience all (or most) of the packages used are loaded or installed here:
## CRAN packages:
if(!exists("cran_req")) { cran_req %<o% "pak" } else { cran_req %<o% cran_req }
if(!exists("bioc_req")) { bioc_req %<o% c() } else { bioc_req %<o% bioc_req }
cran_req <- unique(c(cran_req, "pak", "fs", "shiny", "renv", "R.utils", "data.table", "devtools", "iq", "Rtsne", #"uchardet", # Should not be necessary anymore since Rcy3 replaced it with stringi in version 2.24.0
                     "qs2", "shinyWidgets", "DT",
                     "shinyBS", "stringr", "gplots", "ggplot2", "ggpubr", "gtools", "reshape", "reshape2", "compiler", "stats", "rgl", "ggrepel", "rstudioapi",
                     "modeest", "minpack.lm", "snow", "viridis", "pcaMethods", "impute", "imputeLCMD", "parallel", "coin", "openxlsx", "openxlsx2", "plotly",
                     "Peptides", "xml2", "pdftools", "statmod", "ggpolypath", "venn", "gridExtra", "svDialogs", "htmlwidgets", "magrittr", "tibble", "officer",
                     "hexbin", "igraph", "matlib", "umap", "plyr", "ggnewscale", "shinyjs", "shinyFiles", "TeachingDemos", "shinycssloaders", "tidyr",
                     "ggplotify", "jpeg", "scattermore", "rpanel", "stringi", "lmtest", "ssh", "taxize", "ggdendro", "colorspace", "factoextra", "NbClust",
                     "BH", "plogr", "unimod"))
bioc_req <- unique(c(bioc_req, "biomaRt", "GO.db", "UniProt.ws", "limma", "sva", "qvalue", "MSnbase", "DEP",
                     "Rgraphviz", "RCy3", "siggenes", "pRoloc", "pRolocGUI", "DEqMS", "rawrr", "rbioapi", "png", "Rhdf5lib"))
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
#devtools::install_github("cpanse/rawrr")
#rawrr::installRawFileReaderDLLs() # Deprecated
tst <- try(normalizePath(rawrr:::.rawrrAssembly(), winslash = "/"), silent = TRUE)
if (("try-error" %in% class(tst))||(!file.exists(tst))) {
  rawrr::installRawrrExe()
}

# Fast save and load functions
Src <- paste0(libPath, "/extdata/R scripts/Sources/Save_Load_fun.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

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
LocAnalysis %<o% (WorkFlow %in% c("LOCALISATION", "LOCALIZATION"))

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

# MS raw files map
tstFrMp <- FALSE
while(!tstFrMp) {
  Src <- paste0(libPath, "/extdata/R scripts/Sources/Fractions_Map_editor.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

# Create Experimental Factors
minFact <- c("Experiment", "Replicate")
minFactDesc <- setNames(c("grouping of sample groups to compare", "maximum number of replicates per sample group"),
                        minFact)
if (file.exists("Factors.RData")) {
  load("Factors.RData")
  Factors <- setNames(tmp$Factors, substr(tmp$Factors, 1, 3))
  Factors <- Factors[which(!is.na(Factors))]
  if (length(Factors)) { FactorsLevels %<o% tmp$Levels[Factors] } else { rm(Factors) }
}
if (WorkFlow == "PULLDOWN") {
  minFact <- c(minFact, "Target")
  minFactDesc["Target"] <- "ID (e.g. UniProtKB accession) in the search database(s) of the bait protein"
}
if (WorkFlow == "TIMECOURSE") {
  minFact <- c(minFact, "Time.point")
  minFactDesc["Time.point"] <- "number without unit"
}
if (WorkFlow == "LOCALISATION") {
  minFact <- c(minFact, "Compartment")
  minFactDesc["Compartment"] <- "compartment/subcellular fractionation"
}
if (LabelType == "Isobaric") {
  minFact <- c(minFact, "Isobaric.set")
  minFactDesc["Isobaric.set"] <- "grouping of individual isobarically-labelled samples pooled together (integer)"
}
if (!exists("Factors")) {
  Factors <- minFact
  FactorsLevels %<o% setNames(lapply(Factors, function(x) { c("") }), Factors)
} else { Factors %<o% unique(c(Factors, minFact)) }
if (WorkFlow == "BIOID") {
  Factors <- unique(c(Factors, "Target")) # In this case target is not obligatory
  minFactDesc["Target"] <- "ID (e.g. UniProtKB accession) in the search database(s) of the bait protein"
}
if (!exists("FactorsLevels")) {
  FactorsLevels %<o% setNames(lapply(Factors, function(x) { c("") }), Factors)
}
w <- which(!Factors %in% names(FactorsLevels))
if (length(w)) {
  FactorsLevels[Factors[w]] <- c()
}
Factors %<o% Factors[which(!is.na(Factors))]
FactorsLevels %<o% FactorsLevels[Factors]
#rm(Factors, FactorsLevels)
# Do not use my usual meta-coding approach!!! It is flexible... but not reactive, and almost as bad as shiny to debug!
appNm <- paste0(dtstNm, " - Exp. Factors")
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  titlePanel(tag("u", "Experimental Factors editor"),
             appNm),
  h2(dtstNm), 
  h3("Enter the names and levels of all the Factors describing the samples."),
  br(),
  h4("The following Factors are included by defaults and cannot be removed:"),
  uiOutput("minFact"),
  br(),
  sidebarLayout(
    sidebarPanel(
      textInput("nuFact", "Experimental Factor name(s)", ""),
      actionButton("addFactor", "Create Factor(s)"),
      actionButton("rmvFactor", "Remove Factor(s)"),
      h5(em("Factor names must be:")),
      h5(em(" - at least 3 characters long")),
      h5(em(" - start with a capital letter,")),
      h5(em(" - followed by only lower case letters, numbers or dots.")),
      br(),
      br(),
    ),
    mainPanel(
      h5(em("Add Factor levels:")),
      h5(em(" - use underscores, no spaces or hyphens within level names; sequences of 3 or more consecutive underscores are forbidden!\"\"")),
      h5(em(" - use spaces to add multiple levels at a time")),
      br(),
      span(uiOutput("Message"), style = "color:red"),
      h3("Factor levels:"),
      uiOutput("Factors"),
      br(),
      br(),
      actionButton("saveBtn", "Save")
    )
  )
)
server <- function(input, output, session) {
  # Initialize reactive variables
  FACT <- reactiveVal(Factors)
  FACTLevels <- reactiveVal(FactorsLevels)
  intFact <- c("Replicate", "Isobaric.set")
  dfltInt <- c(2, 1)
  # Create function to update output$Factors for UI
  updtFactUI <- function(reactive = TRUE) {
    if (reactive) {
      FAKT <- FACT()
      FAKTLevels <- FACTLevels()
    } else {
      FAKT <- Factors
      FAKTLevels <- FactorsLevels
    }
    L <- length(FAKT)
    # Update UI
    return(renderUI({
      lst <- vector("list", L*3)
      for (i in 1:L) {
        j <- i*3-2
        Fact <- FAKT[i]
        if (Fact %in% intFact) {
          miN <- dfltInt[match(Fact, intFact)]
          dflt <- max(c(miN, FAKTLevels[[Fact]]))
          txt <- as.character(dflt)
          if (!length(txt)) { txt <- "" }
          lst[[j]] <- list(numericInput(paste0(Fact, "_lev"),
                                        paste0("N. of ", gsub("\\.", " ", Fact), "s = ", txt), dflt, miN, step = 1))
        }
        if (Fact == "Time.point") {
          lst[[j]] <- list(textInput(paste0(Fact, "_lev"),
                                     paste0(Fact, ", levels (numeric(s) only!) = ", paste(FAKTLevels[[Fact]], collapse = " / ")), ""))
        }
        if (!Fact %in% c(intFact, "Time.point")) {
          lst[[j]] <- list(textInput(paste0(Fact, "_lev"),
                                     paste0(Fact, ", levels = ", paste(FAKTLevels[[Fact]], collapse = " / ")), ""))
        }
        lst[[j+1]] <- list(actionButton(paste0(Fact, "_levAdd"), "Add level(s)"))
        lst[[j+2]] <- list(actionButton(paste0(Fact, "_levRmv"), "Remove level(s)"))
      }
      return(lst)
    }))
  }
  #
  output$minFact <- renderUI({ HTML(paste0(" - ", minFact, ": ", minFactDesc, collapse = "<br>")) })
  output$Message <- renderUI({ em(" ") })
  # Initialize
  output$Factors <- updtFactUI(reactive = FALSE)
  #
  # Observers for already extent factors
  #  - Add new Factor level
  sapply(Factors, function(Fact) {
    observeEvent(input[[paste0(Fact, "_levAdd")]], {
      vals <- input[[paste0(Fact, "_lev")]]
      vals <- vals[which(!is.na(vals))]
      if (length(vals)) {
        if ((is.character(vals))&&(Fact != "Target")) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
        tmp2 <- FACTLevels()
        if (Fact %in% intFact) {
          tmp2[[Fact]] <- 1:as.integer(max(c(vals, dfltInt[match(Fact, intFact)])))
        } else {
          if (Fact == "Time.point") {
            vals <- suppressWarnings(as.numeric(vals))
            vals <- vals[which(!is.na(vals))]
          }
          tmp <- unique(c(FACTLevels()[[Fact]], vals))
          tmp <- tmp[which(tmp != "")]
          tmp2[[Fact]] <- tmp
        }
        FACTLevels(tmp2)
        output$Factors <- updtFactUI()
      }
    })
  })
  #  - Remove Factor level
  sapply(Factors, function(Fact) {
    observeEvent(input[[paste0(Fact, "_levRmv")]], {
      vals <- input[[paste0(Fact, "_lev")]]
      vals <- vals[which(!is.na(vals))]
      if (length(vals)) {
        if ((is.character(vals))&&(Fact != "Target")) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
        tmp2 <- FACTLevels()
        if (Fact %in% intFact) {
          tmp2[[Fact]] <- 1:(max(c(as.integer(vals)-1, dfltInt[match(Fact, intFact)])))
        } else {
          if (Fact == "Time.point") {
            vals <- suppressWarnings(as.numeric(vals))
            vals <- vals[which(!is.na(vals))]
          }
          tmp2[[Fact]] <- tmp2[[Fact]][which(!tmp2[[Fact]] %in% vals)]
        }
        FACTLevels(tmp2)
        output$Factors <- updtFactUI()
      }
    })
  })
  #
  # Create new Factor(s)
  observeEvent(input$addFactor, {
    msg <- " "
    Facts <- gsub("[^A-Z,a-z,0-9]", "\\.", unlist(strsplit(input$nuFact, " +")))
    sapply(Facts, function(Fact) {
      if ((nchar(Fact) < 3)||(grepl("^[0-9]", Fact))||(substr(Fact, 1, 3) %in% substr(FACT(), 1, 3))) {
        msg <- "Invalid Factor name! Must be at least 3 characters long and start with a capital letter! The first 3 characters must be unique to this Factor!"
      } else {
        Fact <- paste0(toupper(substr(Fact, 1, 1)), tolower(substr(Fact, 2, nchar(Fact))))
        if (!Fact %in% FACT()) {
          FACT(c(FACT(), Fact))
          tmp <- FACTLevels()
          tmp[Fact] <- list(NULL)
          FACTLevels(tmp)
          # Also, ABSOLUTELY crucial: create new observers for level addition/removal!!!
          observeEvent(input[[paste0(Fact, "_levAdd")]], {
            vals <- input[[paste0(Fact, "_lev")]]
            vals <- vals[which(!is.na(vals))]
            if (length(vals)) {
              if ((is.character(vals))&&(Fact != "Target")) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
              tmp2 <- FACTLevels()
              if (Fact %in% intFact) {
                tmp2[[Fact]] <- 1:as.integer(max(c(vals, dfltInt[match(Fact, intFact)])))
              } else {
                if (Fact == "Time.point") {
                  vals <- suppressWarnings(as.numeric(vals))
                  vals <- vals[which(!is.na(vals))]
                }
                tmp <- unique(c(FACTLevels()[[Fact]], vals))
                tmp <- tmp[which(tmp != "")]
                tmp2[[Fact]] <- tmp
              }
              FACTLevels(tmp2)
              output$Factors <- updtFactUI()
            }
          })
          observeEvent(input[[paste0(Fact, "_levRmv")]], {
            vals <- input[[paste0(Fact, "_lev")]]
            vals <- vals[which(!is.na(vals))]
            if (length(vals)) {
              if ((is.character(vals))&&(Fact != "Target")) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
              tmp2 <- FACTLevels()
              if (Fact %in% intFact) {
                tmp2[[Fact]] <- 1:(max(c(as.integer(vals)-1, dfltInt[match(Fact, intFact)])))
              } else {
                if (Fact == "Time.point") {
                  vals <- suppressWarnings(as.numeric(vals))
                  vals <- vals[which(!is.na(vals))]
                }
                tmp <- unique(c(FACTLevels()[[Fact]], vals))
                tmp2[[Fact]] <- tmp2[[Fact]][which(!tmp2[[Fact]] %in% vals)]
              }
              FACTLevels(tmp2)
              output$Factors <- updtFactUI()
            }
          })
        }
      }
    })
    output$Message <- renderUI({ em(msg) })
    output$Factors <- updtFactUI()
  })
  # Remove extant Factor
  observeEvent(input$rmvFactor, {
    msg <- " "
    tmp <- gsub("[^A-Z,a-z,0-9]", "\\.", input$nuFact)
    if ((nchar(tmp) < 3)||(grepl("^[0-9]", tmp))) {
      msg <- "Invalid Factor name! Must be at least 3 characters long and start with a capital letter!"
    } else {
      tmp <- paste0(toupper(substr(tmp, 1, 1)), tolower(substr(tmp, 2, nchar(tmp))))
      if (!tmp %in% FACT()) {
        msg <- "Cannot remove a non-existent Factor!"
      } else {
        if (tmp %in% minFact) {
          msg <- paste0("Factor ", tmp, " is included by default in this workflow and cannot be removed!")
        } else {
          tmp2 <- FACT()
          tmp2 <- tmp2[which(tmp2 != tmp)]
          FACT(tmp2)
          FACTLevels(FACTLevels()[FACT()])
        }
      }
    }
    output$Message <- renderUI({ em(msg) })
    output$Factors <- updtFactUI()
  })
  #
  observeEvent(input$saveBtn, {
    assign("Factors", FACT(), envir = .GlobalEnv)
    assign("FactorsLevels", FACTLevels(), envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
eval(parse(text = runApp), envir = .GlobalEnv)
#
if ("Target" %in% Factors) {
  FactorsLevels$Target <- FactorsLevels$Target[which(!is.na(FactorsLevels$Target))]
  if (length(FactorsLevels$Target) == 1) {
    FactorsLevels$Target <- c(FactorsLevels$Target, "Control")
    if (length(FactorsLevels$Target) == 1) {
      FactorsLevels$Target <- c(FactorsLevels$Target, "Ctrl")
    }
  }
}
FactorsLevels <- setNames(lapply(Factors, function(fct) {
  x <- FactorsLevels[[fct]]
  x[which(!is.na(x))]
}), Factors)
Factors <- Factors[which(sapply(FactorsLevels[Factors], length) > 0)]
names(Factors) <- substr(Factors, 1, 3)
Factors <- Factors[c("Exp", names(Factors)[which(!names(Factors) %in% c("Exp", "Rep"))], "Rep")]
FactorsLevels <- FactorsLevels[Factors]
tmp <- list(Factors = Factors, Levels = FactorsLevels)
save(tmp, file = "Factors.RData")

#### Code chunk - Edit Experiment map
tstXpMp <- FALSE
while(!tstXpMp) {
  Src <- paste0(libPath, "/extdata/R scripts/Sources/Experiment_Map_editor.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}
#

#### Code chunk - Load and process search database(s)
Src <- paste0(libPath, "/extdata/R scripts/Sources/Process_Fasta_DBs.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

evNm %<o% c("PSM", "Evidence")[(SearchSoft == "MAXQUANT")+1]

#### Code chunk - Load and process annotations
## This includes a QC step in case the database differs slightly from the one used by MQ, or if somehow some IDs have not been properly parsed.
FndAnnotFl <- "Parsed_annotations.RData" %in% list.files(wd)
if (FndAnnotFl) {
  msg <- "Functional annotations backup detected - reuse it?"
  tst <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  if (tst) { msg <- "" } else { FndAnnotFl <- FALSE }
}
GO.col %<o% c("GO", "GO-ID")
ObjNm <- "Annotate"
if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
  if (FndAnnotFl) { tmp <- TRUE } else {
    msg <- "Can you provide functional annotations? (required for GO analysis)"
    tmp <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  }
  ObjNm %<c% tmp
  AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
  tmp <- AllAnsw[1,]
  tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
  tmp$Value <- list(get(ObjNm))
  m <- match(ObjNm, AllAnsw$Parameter)
  if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
}
if (Annotate) {
  if (FndAnnotFl) { loadFun("Parsed_annotations.RData") }
  if (exists("Parsed_annotations")) { Parsed_annotations <- Parsed_annotations } else {
    AnnotFls <- sapply(gsub("\\.fa((s(ta(\\.fas)?)?)|a?)?$", ".txt", fastasTbl$Full), function(x) {
      x2 <- gsub(".+/", "D:/Fasta_databases/", x)
      if (!file.exists(x)) {
        if (file.exists(x2)) {
          fs::file_copy(x2, wd)
          x <- x2
        } else { x <- NA }
      }
      return(x)
    })
    AnnotFls %<o% AnnotFls[which(!is.na(AnnotFls))]
    if (!length(AnnotFls)) {
      moar <- TRUE
      kount <- 0
      while (moar) {
        if (kount == 0) {
          msg <- "No functional annotation file detected. Select one (or more)?"
        } else { msg <- "Select more?" }
        moar <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
        if (moar) {
          msg <- "Choose annotation file(s):"
          filt <- matrix(c("Annotations txt file", "*.txt"), ncol = 2)
          AnnotFls <- c(AnnotFls, choose.files("D:/Fasta_databases/*.txt", msg, TRUE, filt))
          AnnotFls <- normalizePath(AnnotFls, winslash = "/")
          AnnotFls <- AnnotFls[which(!is.na(AnnotFls))]
          kount <- kount + 1
        }
      }
    }
    if (!length(AnnotFls)) {
      warning("No annotations file(s) provided, skipping annotations!")
      Annotate <- FALSE
    } else {
      source(parSrc, local = FALSE)
      Parsed_annotations <- lapply(AnnotFls, function(x) { #x <- AnnotFls[1]
        # If the annotations is not present locally, make a local copy
        if (!file.exists(basename(x))) { fs::file_copy(x, wd) }
        # Parse it
        return(Format.DB_txt(x, usePar = TRUE, cl = parClust))
      })
      Parsed_annotations <- dplyr::bind_rows(Parsed_annotations)
      tst1 <- unlist(strsplit(Parsed_annotations$`GO-ID`, ";"))
      tst2 <- unlist(strsplit(Parsed_annotations$GO, ";"))
      tst3 <- data.table(A1 = tst1, A2 = tst2)
      tst3 <- tst3[, list(x = unique(A2)), by = list(Group.1 = A1)]
      tst3 <- as.data.frame(tst3)
      stopifnot(length(tst3$x) == length(unique(tst3$x)), "character" %in% class(tst3$x))
      #View(tst3[which(sapply(tst3$x, length) > 1),])
      #View(tst3[which(sapply(tst3$x, length) == 0),])
      saveFun(Parsed_annotations, file = "Parsed_annotations.RData")
      #loadFun("Parsed_annotations.RData")
      #View(Parsed_annotations[, c("GO", "GO-ID")])
    }
  }
  # Check GO for name degeneracies
  # (the same term has had different names in different files, and we allow for multiple files)
  w <- which(nchar(Parsed_annotations$GO) > 0)
  # rows-to-names
  tst1 <- listMelt(strsplit(Parsed_annotations$GO[w], ";"), w, c("name", "row"))
  # names-to-rows-to-IDs
  tst2 <- set_colnames(aggregate(tst1$row, list(tst1$name), list), c("name", "rows"))
  tst2$ID <- gsub(".*\\[|\\]$", "", tst2$name)
  tst1$ID <- tst2$ID[match(tst1$name, tst2$name)]
  tst3 <- set_colnames(aggregate(tst2$name, list(tst2$ID), unique), c("ID", "names"))
  tst3$L <- sapply(tst3$names, length)
  if (max(tst3$L) > 1) { # this would indicate degeneracy and the need to fix names
    tst3 <- tst3[which(tst3$L > 1),]
    tst3$name <- sapply(tst3$names, function(x) { x[[1]] })
    rws <- unique(unlist(tst2$rows[which(tst2$ID %in% tst3$ID)]))
    tst1 <- tst1[which(tst1$row %in% rws),]
    tst2 <- tst2[which(tst2$ID %in% tst3$ID),]
    tst2$name <- tst3$name[match(tst2$ID, tst3$ID)]
    w2 <- which(tst1$ID %in% tst2$ID)
    tst1$name[w2] <- tst2$name[match(tst1$ID[w2], tst2$ID)]
    tst1 <- as.data.table(tst1)
    tst1 <- tst1[, list(GO = paste(name, collapse = ";"),
                        ID = paste(ID, collapse = ";")),
                 by = list(tst1$row)]
    tst1 <- as.data.frame(tst1)
    Parsed_annotations[w, c("GO", "GO-ID")] <- tst1[match(w, tst1$tst1), c("GO", "ID")]
    saveFun(Parsed_annotations, file = "Parsed_annotations.RData")
  }
}
Annotate <- exists("Parsed_annotations") # Update condition
if (Annotate) {
  #db <- db[, which(!colnames(db) %in% annot.col)]
  kol <- colnames(Parsed_annotations)
  annot.col %<o% kol[which(!kol %in% c("Accession", "id", "ID", "Names", "Sequence", "MW (Da)"))]
  annot.col2 %<o% annot.col[which(!annot.col %in% colnames(db))]
  annot.col3 %<o% annot.col[which(annot.col %in% colnames(db))]
  if (length(annot.col2)) { db[, annot.col2] <- NA }
  mtch <- match(db$`Protein ID`, Parsed_annotations$Accession)
  db[, annot.col2] <- Parsed_annotations[mtch, annot.col2]
  if (length(annot.col3)) {
    for (kol in annot.col3) { #kol <- annot.col3[1]
      w <- which(!is.na(mtch)) # check that there is a valid match...
      w <- w[which((is.na(db[w, kol]))|(db[w, kol] %in% c("", "NA", "NaN")))] #... and that it is useful!
      db[w, kol] <- Parsed_annotations[mtch[w], kol]
      w <- which((is.na(db[[kol]]))|(db[[kol]] %in% c("", "NA", "NaN"))) #... and that it is useful!
      db[w, kol] <- ""
    }
  }
  tst1 <- unlist(strsplit(db$`GO-ID`, ";"))
  tst2 <- unlist(strsplit(db$GO, ";"))
  tst3 <- data.table(A1 = tst1, A2 = tst2)
  tst3 <- tst3[, list(x = unique(A2)), by = list(Group.1 = A1)]
  tst3 <- as.data.frame(tst3)
  #View(tst3[which(sapply(tst3$x, length) > 1),])
  #View(tst3[which(sapply(tst3$x, length) == 0),])
  stopifnot(length(tst3$x) == length(unique(tst3$x)), "character" %in% class(tst3$x))
  db$Ontology <- NULL # Temporary fix for now, this column is broken
  #
  write.csv(db, "Parsed, annotated search db.csv", row.names = FALSE)
  #db <- read.csv("Parsed, annotated search db.csv", check.names = FALSE, colClasses = "character")
}

# Create experiment Factors shortcuts
Aggregates %<o% Factors
a <- substr(Aggregates, 1, 3)
if (length(unique(a)) != length(a)) {
  stop("Factors must start with a unique 3 letter tag! Edit the experiment map columns and restart.")
} else { names(Aggregates) <- a }
for (i in Aggregates) {
  temp <- sort(unique(Exp.map[[i]]))
  tst <- sum(as.character(suppressWarnings(as.numeric(temp))) != temp) # Are they numerics?
  if ((!is.na(tst))&&(!tst)) { temp <- as.numeric(temp) }
  tst <- sum(as.character(suppressWarnings(as.integer(temp))) != temp) # Are they integers?
  if ((!is.na(tst))&&(!tst)) { temp <- as.integer(temp) }
  substr(i, 1, 3) %<c% temp
}
if (exists("Rep")) { names(Rep) <- gsub("rep", "", Rep) }
if (LabelType == "Isobaric") {
  if (!exists("Iso")) {
    stop("See older versions of this script - but in the newer version I would always expect Iso to exist here.")
  } else { test.iso.set %<o% TRUE }
}
MQ.Frac %<o% sort(suppressWarnings(unique(as.integer(unlist(strsplit(as.character(Exp.map$Fractions), ";"))))), na.last = TRUE)

# Test
if ((SearchSoft == "MAXQUANT")&&(LabelType == "Isobaric")) {
  tst <- sum(!paste0("Reporter intensity ", Exp.map$Isobaric.label) %in% colnames(ev))
  stopifnot(tst == 0)
  # If this is not correct, then we should consider re-introducing older "Code chunk - Isobarically-labelled samples only! Check whether 1 must be subtracted from MaxQuant's isobaric labels index"
}

saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

# Fractions map
Exp.map$Use <- as.logical(Exp.map$Use)
MQ.Exp <- MQ.Exp[which(MQ.Exp %in% unique(unlist(Exp.map$MQ.Exp[which(Exp.map$Use)])))]
if (file.exists(FracMapPath)) {
  Frac.map %<o% read.csv(FracMapPath, check.names = FALSE)
  if (("Parent sample" %in% colnames(Frac.map))&&(!"MQ.Exp" %in% colnames(Frac.map))) {
    Frac.map$MQ.Exp <- Frac.map$"Parent sample"
  }
  Frac.map <- Frac.map[which(Frac.map$Use),]
  MQ.Exp <- MQ.Exp[which(MQ.Exp %in% unique(unlist(Frac.map$MQ.Exp[which(Frac.map$Use)])))]
  Exp.map <- Exp.map[which(sapply(Exp.map$MQ.Exp, function(x) { sum(x %in% MQ.Exp) }) > 0),]
  Frac.map <- Frac.map[which(Frac.map$MQ.Exp %in% MQ.Exp),]
  #Frac.map$"Raw file" <- gsub("\\\\", "/", Frac.map$"Raw file") # Redundant
  #ev$"Raw file path" <- gsub_Rep("\\\\", "/", ev$"Raw file path")
  m <- match(ev$`Raw file`, rawFiles2)
  if (!sum(!is.na(m))) {
    stop()
    #ev$`Raw file` <- gsub_Rep("\\.[^\\.]+$", "", ev$`Raw file`)
    #m <- match(ev$`Raw file`, rawFiles2)
  }
  #rawFiles <- gsub_Rep("\\\\", "/", rawFiles)
  # if (!"Raw file path" %in% colnames(ev)) { ev$"Raw file path" <- rawFiles[match(ev$`Raw file`, rawFiles2)] } else {
  #   ev$"Raw file path" <- gsub_Rep("\\\\", "/", ev$"Raw file path")
  # }
  # if (SearchSoft != "DIANN") {
  #   Frac.map$"Raw file_Full" <- Frac.map$"Raw file"
  #   if (SearchSoft == "MAXQUANT") {
  #     Frac.map$"Raw file" <- gsub(".*/|\\.((raw)|(mzX?ML)|(d)|(dia))$", "", Frac.map$"Raw file", ignore.case = TRUE)
  #   } else {
  #     Frac.map$"Raw file" <- gsub("\\.((raw)|(mzX?ML)|(d)|(dia))$", "", Frac.map$"Raw file", ignore.case = TRUE)
  #   }
  # }
  if ("MQ.Exp" %in% colnames(ev)) {
    tst <- sum(is.na(ev$MQ.Exp)) == nrow(ev)
    if (tst) { ev$MQ.Exp <- NULL }
  }
  test <- sort(unique(Frac.map$MQ.Exp))
  if (sum(test != sort(MQ.Exp))) { stop("Column \"Experiment\" not defined properly in Fractions map!") }
  Frac.map <- Frac.map[which(Frac.map$MQ.Exp %in% MQ.Exp),]
  Frac.map$Experiment <- sapply(Frac.map$MQ.Exp, function(x) { #x <- Frac.map$MQ.Exp[1]
    x1 <- unlist(unique(Exp.map$Experiment[which(sapply(Exp.map$MQ.Exp, function(y) { x %in% unlist(y) }))]))
    if (length(x1) == 1) { return(x1) } else {
      if (length(x1) > 1) { stop(paste0(x, " - each raw file must be mapped to exactly one Experiment!")) } else {
        return(NA)
      }
    }
  })
  Frac.map <- Frac.map[which(!is.na(Frac.map$Experiment)),]
  if (("Replicate" %in% colnames(Frac.map))&&(sum(sort(unique(Frac.map$Replicate)) != sort(Rep)) != 0)) {
    stop("Replicates from Fractions map and Experiment map do not match!")
  }
  # if (SearchSoft == "DIANN") {
  #   if (length(which(!unique(ev$"Raw file path") %in% unique(Frac.map$"Raw file")))) {
  #     warning("Some evidences do not have corresponding raw files in Fractions map and shall thus be removed. Check that this is ok!")
  #     ev <- ev[which(ev$"Raw file path" %in% unique(Frac.map$"Raw file")),]
  #   }
  #   ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
  #   if (length(which(!unique(Frac.map$"Raw file") %in% unique(ev$"Raw file path")))) {
  #     warning("There are some raw files for which not a single peptide evidence is present in the data!")
  #     Frac.map <- Frac.map[which(Frac.map$"Raw file" %in% unique(ev$"Raw file path")),]
  #   }
  # } else {
  #   ev$"Raw file" <- gsub_Rep("\\.((raw)|(mzX?ML)|(d)|(dia))$", "", ev$"Raw file", ignore.case = TRUE)
  #   Frac.map$"Raw file" <- gsub("\\.((raw)|(mzX?ML)|(d)|(dia))$", "", Frac.map$"Raw file", ignore.case = TRUE)
  #   if (length(which(!unique(ev$"Raw file") %in% unique(Frac.map$"Raw file")))) {
  #     warning("Some evidences do not have corresponding raw files in Fractions map and shall thus be removed. Check that this is ok!")
  #     ev <- ev[which(ev$"Raw file" %in% unique(Frac.map$"Raw file")),]
  #   }
  #   ev2fr %<o% match(ev$"Raw file", Frac.map$"Raw file")
  #   if (length(which(!unique(Frac.map$"Raw file") %in% unique(ev$"Raw file")))) {
  #     warning("There are some raw files for which not a single peptide evidence is present in the data!")
  #     Frac.map <- Frac.map[which(Frac.map$"Raw file" %in% unique(ev$"Raw file")),]
  #   }
  # }
  stopifnot("Raw file" %in% colnames(Frac.map),
            "Raw file path" %in% colnames(ev))
  ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
  tst <- is.na(ev2fr)
  if (sum(tst)) {
    w <- which(tst)
    for (k in c("Raw file", "Raw files name")) { #k <- "Raw file" #k <-  "Raw files name"
      w2 <- which(Frac.map[[k]] %in% ev$"Raw file"[w])
      if (length(w2)) {
        m <- match(Frac.map[w2, k], ev$"Raw file")
        tst2 <- gsub(".*/|\\.[^\\.]+$", "", Frac.map$"Raw file"[w2]) == gsub(".*/|\\.[^\\.]+$", "", ev$"Raw file path"[m])
        stopifnot(sum(is.na(tst2)) == 0, sum(!tst2) == 0)
        Frac.map$"Raw file"[w2] <- ev$"Raw file path"[m]
        ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
      }
    }
  }
  tst <- sum(is.na(ev2fr))
  if (tst) {
    warning(paste0("Removing ", tst, " PSMs not matching selected raw files..."))
    ev <- ev[which(!is.na(ev2fr)),]
    ev2fr <- ev2fr[which(!is.na(ev2fr))]
  }
  if (!"MQ.Exp" %in% colnames(ev)) {
    ev$MQ.Exp <- Frac.map$MQ.Exp[ev2fr]
    w <- which(is.na(ev$MQ.Exp))
    if (length(w)) { warning("Some PSMs do not have a corresponding Samples, is this expected?\n(This is normal if you decided to not use all samples...)") }
    if (LabelType == "Isobaric") { Iso <- sort(unique(Exp.map$Isobaric.set)) }
  }
  # if (exists("MQ.Exp2discrd")) {
  #   ev <- ev[which(!ev$MQ.Exp %in% MQ.Exp2discrd),]
  #   ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
  #   Frac.map <- Frac.map[which(!Frac.map$MQ.Exp %in% MQ.Exp2discrd),]
  # }
  if ("Replicate" %in% colnames(Frac.map)) {
    Frac.map$Unique.Frac.ID <- apply(Frac.map[, c("Replicate", "MQ.Exp", "Fraction")], 1, function(x) {
      paste0("Rep", paste(x, collapse = "_"))
    })
    ev$Replicate <- Frac.map$Replicate[ev2fr]
  } else { Frac.map$Unique.Frac.ID <- apply(Frac.map[, c("MQ.Exp", "Fraction")], 1, paste, collapse = "_") }
  if (LabelType == "Isobaric") {
    if (!"Isobaric.set" %in% colnames(Frac.map)) {
      if (test.iso.set) { stop("The fractions map does not include an isobaric set column!") } else {
        Frac.map$Isobaric.set <- 1
        ev$Isobaric.set <- 1
      }
    } else { ev$Isobaric.set <- Frac.map$Isobaric.set[ev2fr] }
  }
  # Filter for ones to keep
  MQ.Exp <- MQ.Exp[which(MQ.Exp %in% Frac.map$MQ.Exp)]
  Exp.map <- Exp.map[which(sapply(Exp.map$MQ.Exp, function(x) { sum(x %in% MQ.Exp) > 0 })),]
  Frac.map$Fraction <- as.numeric(gsub("^ | $", "", Frac.map$Fraction))
  ev <- ev[which(ev$MQ.Exp %in% MQ.Exp),]
  ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
  # Final test
  test <- data.frame(Ev = sort(unique(ev$MQ.Exp)),
                     Fraction.map = sort(unique(Frac.map$MQ.Exp)),
                     Exp.map = sort(unique(unlist(Exp.map$MQ.Exp))),
                     MQ.Exp = sort(MQ.Exp))
  if (max(apply(test, 1, function(x) { length(unique(x)) })) > 1) {
    stop(paste0("The MQ.Exp object, the Fractions and Experiment map, and the ", evNm, " file do not match!!!"))
    #apply(test, 1, unique)
  }
  if ("Norma.groups" %in% colnames(Frac.map)) {
    warning("Column \"Norma.groups\" in the Fractions Map is deprecated, use \"PTM-enriched\" instead!")
  }
  Unique.Frac %<o% data.frame(Unique.Frac.ID = unique(Frac.map$Unique.Frac.ID))
  Unique.Frac$Raw.files <- sapply(Unique.Frac$Unique.Frac.ID, function(x) {
    list(Frac.map$"Raw file"[which(Frac.map$Unique.Frac.ID == x)])
  })
  ev$"Unique Frac" <- Frac.map$Unique.Frac.ID[ev2fr]
} else {
  kol <- "MQ.Exp"
  if (length(Exp) == 1) {
    ev$Experiment <- Exp
    kol <- c(kol, "Experiment")
  } else { stop("There are several Experiments in Exp.map, yet Fractions map is missing!") }
  if (LabelType == "Isobaric") {
    if (length(Iso) == 1) {
      ev$Isobaric.set <- Iso
      kol <- c(kol, "Isobaric.set")
    } else { stop("There are several Isobaric sets in Exp.map, yet Fractions map is missing!") }
  }
  Frac.map %<o% set_colnames(aggregate(ev[, kol], list(ev$"Raw file path"), unique), c("Raw file", kol))
  # if (exists("MQ.Exp2discrd")) {
  #   ev <- ev[which(!ev$MQ.Exp %in% MQ.Exp2discrd),]
  #   Frac.map <- Frac.map[which(!Frac.map$MQ.Exp %in% MQ.Exp2discrd),]
  # }
}
# Filter
Exp.map$Use <- as.logical(Exp.map$Use)
Frac.map$Use <- as.logical(Frac.map$Use)
Exp.map <- Exp.map[which(Exp.map$Use),]
Frac.map <- Frac.map[which(Frac.map$Use),]
stopifnot(sum(sort(unique(Exp.map$Experiment)) != sort(unique(Frac.map$Experiment))) == 0)
stopifnot(sum(sort(unique(unlist(Exp.map$MQ.Exp))) != sort(unique(unlist(Frac.map$MQ.Exp)))) == 0)
ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file") # Update it
tst <- listMelt(lapply(1:nrow(Exp.map), function(x) { Exp.map$MQ.Exp[x] }), 1:nrow(Exp.map))
tst$L1 <- as.integer(tst$L1)
tst <- tst[order(tst$L1, tst$value),]
MQ.Exp <- sort(unique(tst$value))
if (LabelType == "LFQ") { stopifnot(length(MQ.Exp) == length(unique(MQ.Exp)))}
ev <- ev[which(ev$MQ.Exp %in% MQ.Exp),]
ev <- ev[which(ev$"Raw file path" %in% Frac.map$"Raw file"),]
# Update values
ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
rawFiles %<o% unique(ev$"Raw file path")
rawFiles2 %<o% unique(ev$`Raw file`)
#
if (!"MQ.Exp" %in% colnames(ev)) { ev$MQ.Exp <- Frac.map$MQ.Exp[ev2fr] }
if (!"Parent sample" %in% colnames(ev)) { ev$"Parent sample" <- ev$MQ.Exp } # Synonym for now
if (!"Experiment" %in% colnames(ev)) { ev$Experiment <- Frac.map$Experiment[ev2fr] }
if (!"Fraction" %in% colnames(ev)) { ev$Fraction <- Frac.map$Fraction[ev2fr] }
if (!"Replicate" %in% colnames(ev)) {
  tmp <- listMelt(Exp.map$MQ.Exp, Exp.map$Replicate)
  ev$Replicate <- tmp$L1[match(ev$MQ.Exp, tmp$value)]
}
tmp <- listMelt(Exp.map$MQ.Exp, Exp.map$Experiment)
ev$Experiment <- tmp$L1[match(ev$MQ.Exp, tmp$value)]

# Intensity columns at PSM level
ev.col %<o% c(Original = "Intensity")
if (LabelType == "Isobaric") {
  ev.ref %<o% c(Original = "Reporter intensity ")
  korkol <- grep(topattern(paste0(ev.ref["Original"], "corrected ")), colnames(ev), value = TRUE)
  ev <- ev[, which(!colnames(ev) %in% korkol)]
}

# Update Factors
FactorsLevels <- setNames(lapply(Factors, function(Fact) {
  unique(Exp.map[[Fact]])
}), Factors)
w <- which(sapply(FactorsLevels, length) > 0)
Factors <- Factors[w]
FactorsLevels <- FactorsLevels[Factors]

# Protein headers for shiny
Src <- paste0(libPath, "/extdata/R scripts/Sources/protHeaders_for_shiny.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# Proteins of interest
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# Targets
# Sometimes the user does not fill the Target factor with valid protein IDs... but this is what we would actually need.
# Here, if necessary, we will remap those to valid IDs:
Src <- paste0(libPath, "/extdata/R scripts/Sources/Targets.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# KnockOut, KnockIn or KnockDown
tst <- tolower(gsub("[- _]", "", Factors))
if (sum(c("knockout", "knockin", "knockdown") %in% tst)) {
  w <- which(c("knockout", "knockin", "knockdown") %in% tst)
  # There should be only one for now, because all three share the same 3-characters root = "Kno"
  # This should evolve but will be difficult, knowing how complex this script is now.
  prot.list %<o% unique(c(unique(Exp.map[[Factors["Kno"]]]), prot.list))
  prot.list_pep %<o% unique(c(unique(Exp.map[[Factors["Kno"]]]), prot.list_pep))
}
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList2.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# Protein headers for shiny (update)
Src <- paste0(libPath, "/extdata/R scripts/Sources/protHeaders_for_shiny.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Define analysis parameters
Src <- paste0(libPath, "/extdata/R scripts/Sources/rep_Parameters_editor_Main.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
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

# Start of processing of evidences table
ReportCalls <- AddSpace2Report()
ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_fpar(Report, fpar(ftext(\"Processing of the ",
                                                      c("evidence", "PSMs", "PSMs", "identifications")[mSft],
                                                      " table\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))"))
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
}
ev <- ev[wHt1,]

w <- grep("CONTAMINANT", colnames(ev), ignore.case = TRUE)
if (length(w) > 1) { warning("Hmmm..., you might wanna check what is happening here...") } else {
  colnames(ev)[w] <- "Potential contaminant"
}
for (i in c("Potential contaminant", "Reverse")) {
  w <- which(is.na(ev[[i]]))
  ev[w, i] <- ""
}

#### Code chunk - Update evidences-to-protein mappings
# FragPipe only reports one protein per PSM, and among other search software, MaxQuant at least is also not fully exhaustive.
# While this may be correct from their point of view, I think I should report all matches.
setwd(wd)
#ev$Proteins <- gsub(";CON_", ";", gsub("^CON_", "", gsub(";CON__", ";", gsub("^CON__", "", ev$Proteins))))
if (Update_Prot_matches) {
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              paste0("body_add_fpar(Report, fpar(ftext(\" - Checking ", names(SearchSoft), "'s peptide sequence to protein assignments:\", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$just))"))
  Reuse_Prot_matches <- FALSE
  ObjNm <- "Reuse_Prot_matches"
  if ("evmatch.RData" %in% list.files(wd)) {
    if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
      msg <- "Some backed-up protein matches are available in the folder, do you want to reuse them?\n"
      tmp <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
      if (is.na(tmp)) { tmp <- FALSE }
      ObjNm %<c% tmp
      tmp <- AllAnsw[1,]
      tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
      tmp$Value <- list(get(ObjNm))
      m <- match(ObjNm, AllAnsw$Parameter)
      if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
    }
  }
  if (Reuse_Prot_matches) { loadFun(paste0(wd, "/evmatch.RData")) }
  if (exists("evmatch")) {
    temp <- evmatch
  } else {
    source(parSrc, local = FALSE)
    temp <- ProtMatch2(unique(ev$Sequence), db,
                       cl = parClust) # (ignore the warning for now until we remove contaminant evidences)
  }
  wh1 <- which(ev$Sequence %in% temp$Sequence)
  wh2 <- which(temp$Sequence %in% ev$Sequence)
  mtch1 <- match(ev$Sequence[wh1], temp$Sequence)
  if (!"Proteins" %in% colnames(ev)) { ev$Proteins <- "" } else {
    source(parSrc, local = FALSE)
    tmpPs <- unique(temp$Sequence[wh2])
    tmpE <- ev$Proteins[match(tmpPs, ev$Sequence)]
    tmpP <- temp$Proteins[match(tmpPs, temp$Sequence)]
    tmpE <- strsplit(tmpE, ";")
    tmpP <- strsplit(tmpP, ";")
    f0 <- function(x) { paste(sort(x), collapse = ";") }
    tst1 <- parSapply(parClust, tmpE, f0)
    tst2 <- parSapply(parClust, tmpP, f0)
    wN <- which(tst1 != tst2)
    ViewTst <- FALSE
    if ((length(wN))&&(ViewTst)) {
      # Below some code to help with investigations...
      tst <- data.frame(Seq = tmpPs[wN],
                        Orig = tst1[wN],
                        Corr = tst2[wN])
      tst$Orig <- strsplit(tst$Orig, ";")
      #View(tst)
      tst$Corr <- strsplit(tst$Corr, ";")
      f0 <- function(x, y) { x[which(!x %in% y)] }
      tst$In_Orig_only <- lapply(1:nrow(tst), function(x) { f0(unlist(tst$Orig[[x]]),
                                                               unlist(tst$Corr[[x]])) })
      tst$In_Corr_only <- lapply(1:nrow(tst), function(x) { f0(unlist(tst$Corr[[x]]),
                                                               unlist(tst$Orig[[x]])) })
      tst$In_Orig_only_in_db <- sapply(tst$In_Orig_only, function(x) { sum(x %in% db$`Protein ID`) })
      sum(unlist(tst$In_Orig_only_in_db))
      sum(!unlist(tst$In_Orig_only_in_db))
      w <- which(tst$In_Orig_only_in_db > 0)
      View(tst[w,])
      # lapply(1:nrow(tst), function(x) { f0(unlist(tst$Orig[[x]]),
      #                                      unlist(tst$Corr[[x]])) })
    }
    tst <- sum(tst1 != tst2)
    if (tst) {
      msg <- paste0("Corrected ", tst, " out of ", nrow(temp), " assignments (~", round(tst/nrow(temp), 2), "%)!
Prior investigations have identified 3 cases:
 1) Protein present in original column but not corrected results:
   a) If redundant entries were present in the search fasta, then they would have usually been filtered by this workflow when it loads and parses the fasta database.
   b) More rarely, the accession is present in the db... but for every peptide we have checked so far the protein sequence was not compatible with it (even allowing for I/L ambiguity)!!!
 2) Present in corrected but not original. We have multiple times verified that some search software (at least MaxQuant) miss some proteins despite a peptide being a canonical digestion product for it!
All this supports the need to stringently check a search engine's assignments!
Of course, in some cases (e.g. all assignments corrected) discrepancies can also be due to differences in how fasta headers are processed to extract protein accessions!
")
      if (SearchSoft == "FRAGPIPE") {
        msg <- paste0(msg, " In addition, FragPipe only reports one protein per PSM for some reason...
Hence it is to be expected that a very high percentage of assignments should be corrected there.
")
      }
    } else { msg <- "All assignments validated.\n" }
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
  }
  ev$Proteins[wh1] <- temp$Proteins[mtch1]
  if (!Reuse_Prot_matches) {
    # Save in case you need to reload:
    evmatch <- ev[, c("Sequence", "Proteins")]
    saveFun(evmatch, file = paste0(wd, "/evmatch.RData"))
  }
} else {
  kol <- c("Leading proteins", "Proteins")
  kol <- kol[which(kol %in% colnames(ev))]
  tmp <- ev[, kol, drop = FALSE]
  for (k in kol) { tmp[[k]] <- strsplit(tmp[[k]], ";") }
  ev$Proteins <- parApply(parClust, tmp, 1, function(x) {
    paste(unique(unlist(x)), collapse = ";")
  })
}
tst <- unique(unlist(strsplit(ev$Proteins, ";")))
if ("NA" %in% tst) { stop("\"NA\" is not an accepted protein accession!") }
#View(ev[, c("Sequence", "Proteins")])

#test <- c()
#if (!is.null(prot.list)) {
#  a <- paste0(";", ev$Proteins, ";")
#  b <- prot.list
#  a1 <- paste0(";", b, ";")
#  test <- unique(unlist(sapply(a1, function(x) {ev$id[grep(x, a)]})))
#}
#kol <- which(toupper(colnames(ev)) %in% c("CONTAMINANT", "POTENTIAL CONTAMINANT"))
#ev <- ev[which((is.na(ev[[kol]]))|(ev[[kol]] == "")|(ev$id %in% test)),]
# Test if there are still any evidences without any matching proteins from the database
ev$"Tryptic peptide?" <- TRUE
w <- which(ev$Proteins == "")
l <- length(w)
if (l) {
  tst <- (l>1)+1
  msg <- paste0("There ", c("is", "are")[tst], " ", length(w), " PSM", c("", "s")[tst],
                " (", round(100*l/nrow(ev)), "%) without any matching protein from the database, probably from ",
                c("a ", "")[tst], "non-tryptic peptide", c("", "s")[tst], ".")
  ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
  ev$"Tryptic peptide?"[w] <- FALSE
  temp <- data.frame(Seq = unique(ev$Sequence[w]))
  temp$Proteins <- sapply(temp$Seq, function(x) {
    paste(db$"Protein ID"[grep(x, db$Sequence)], collapse = ";")
  })
  ev$Proteins[w] <- temp$Proteins[match(ev$Sequence[w], temp$Seq)]
  ev <- ev[which(ev$Proteins != ""),]
}
# Also remove those protein columns we will re-create later
ev$"Leading proteins" <- NULL
ev$"Leading razor protein" <- NULL

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

# Filter to keep only evidences with valid quantitative values:
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
  u <- sort(as.numeric(unique(Exp.map$Isobaric.label)))
  w1 <- which(tmpIso %in% u)
  w2 <- which(!tmpIso %in% u)
  if (length(w2)) {
    ev <- ev[, which(!colnames(ev) %in% unlist(sapply(tmpIso[w2], function(x) {
      paste0(c("Reporter intensity corrected ", "Reporter intensity ", "Reporter intensity count "), x)
    })))]
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
###########################################################################################################
# MA plots                                                                                                #
###########################################################################################################
# This will create MA plots to check whether any shape correction is required on the data (loess or vsn)
# If you think shape correction is required, set parameter "Norma.Pep.Intens.Shape" to loess or vsn,
# then reload parameters.
#
# What are the options here?
# We can have:
# - label-free, one experiment => skip
# - label-free, multiple MQ experiments (+/- fractions)
# - isobaric (+/- MQ experiments/fractions)
# Update:
ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file") # Update again
#
if ((length(MQ.Exp) > 1)||(LabelType == "Isobaric")) { # Should be always TRUE
  source(parSrc, local = FALSE)
  data <- ev
  colnames(data)[which(colnames(data) == "MQ.Exp")] <- "Parent sample"
  if (("PTM-enriched" %in% colnames(Frac.map))&&(sum(Modifs$"Full name" %in% Frac.map$"PTM-enriched"))) {
    data$"PTM-enrich." <- Frac.map$"PTM-enriched"[ev2fr]
    data$"PTM-enrich."[which(is.na(data$"PTM-enrich."))] <- ""
  }
  data <- data[which(data$Reverse != "+"),]
  data <- data[which((is.na(data$"Potential contaminant"))|(data$"Potential contaminant" != "+")),]
  if (LabelType == "Isobaric") {
    kol <- grep(paste0(topattern(ev.ref["Original"]), "[0-9]+$"), colnames(data), value = TRUE)
  } else { kol <- ev.col["Original"] }
  w <- which(rowSums(data[, kol, drop = FALSE], na.rm = TRUE) > 0)
  data <- data[w,]
  if (!"Fraction" %in% colnames(data)) { data$Fraction <- 1 }
  Fraction <- sort(unique(data$Fraction), decreasing = FALSE)
  Experiment <- Exp
  kols <- c("Parent sample", "Fraction", "Experiment")
  if (LabelType == "Isobaric") {
    X <- "Label"
    kols <- c("Fraction", "Parent sample", "Experiment") # The order matters!
    tst <- sapply(kols, function(x) { length(unique(data[[x]])) })
    w1 <- which(tst > 1)
    w2 <- which(tst >= 1) 
    if (length(w1)) { Y <- kols[w1[1]] } else { Y <- kols[w2[1]] }
  }
  if (LabelType == "LFQ") {
    kols <- c("Parent sample", "Fraction", "Experiment") # The order matters!
    tst <- sapply(kols, function(x) { length(unique(data[[x]])) })
    w1 <- which(tst > 1)
    w2 <- which(tst >= 1) 
    X <- kols[w1[1]]
    if (length(w1) > 1) { Y <- kols[w1[2]] } else { Y <- kols[w2[2]] }
  }
  kols <- kols[which(!kols %in% c(X, Y))]
  if (length(kols) > 1) {
    data$Group <- apply(data[, kols], 1, paste, collapse = " ")
    Grpkol <- "Group"
  } else { Grpkol <- kols }
  grps <- sort(unique(data[[Grpkol]]))
  ReportCalls <- AddSpace2Report()
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              paste0("body_add_fpar(Report, fpar(ftext(\"MA plot", c("", "s")[(length(grps) > 1)+1], ":\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))"))
  ReportCalls$Objects$MA_groups <- c()
  ReportCalls$Plots$MA_plots <- list()
  ReportCalls$Calls <- append(ReportCalls$Calls, list())
  LRepCalls <- length(ReportCalls$Calls)
  dir <- paste0(wd, "/Workflow control/MA plots")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  #
  for (grp in grps) { #grp <- grps[1]
    data2 <- data[which(data[[Grpkol]] == grp),]
    lsKl <- c("Modified sequence", Y)
    if (LabelType == "LFQ") { lsKl <- c(lsKl, X) }
    if ("PTM-enrich." %in% colnames(data2)) { lsKl <- c(lsKl, "PTM-enrich.") }
    ls <- lapply(lsKl, function(kl) { data2[[kl]] })
    tmp <- do.call(paste, c(data2[, lsKl], sep = "---"))
    if (LabelType == "Isobaric") {
      kol2 <- gsub(topattern(ev.ref["Original"]), "", kol)
      data2a <- as.data.table(data2[, kol])
      colnames(data2a) <- kol
      data2a$Group <- tmp
      data2a <- data2a[, lapply(.SD, sum, na.rm = TRUE), keyby = Group]
      colnames(data2a) <- c("Group", kol)
      data2a <- as.data.frame(data2a)
      data2a[, kol2] <- log10(data2a[, kol])
      data2a <- data2a[, which(!colnames(data2a) %in% kol)]
      data2a[, lsKl] <- data2[match(data2a$Group, tmp), lsKl]
      data2a$Group <- NULL
    }
    if (LabelType == "LFQ") {
      if (X == "Parent sample") { kol2 <- get("MQ.Exp") } else { kol2 <- get(X) }
      data2a <- data.table(Intensity = data2[[ev.col["Original"]]], Group = tmp)
      data2a <- data2a[, list(`log10(Intensity)` = sum(Intensity, na.rm = TRUE)),
                       keyby = Group]
      data2a$`log10(Intensity)` <- log10(data2a$`log10(Intensity)`)
      data2a <- as.data.frame(data2a)
      data2a[, lsKl] <- data2[match(data2a$Group, tmp), lsKl]
      data2a$Group <- NULL
      data2 <- spread(data2a, X, "log10(Intensity)")
    }
    tmp <- data2[, kol2]
    avkol <- rowMeans(data2[, kol2], na.rm = TRUE)
    data2[, kol2] <- sweep(data2[, kol2], 1, avkol, "-")/log10(2)
    kol <- c("Modified sequence", Y)
    if ("PTM-enrich." %in% colnames(data2)) { kol <- c(kol, "PTM-enrich.") }
    #data2 <- as.data.table(data2)
    data2 <- set_colnames(reshape2::melt(data2, id.vars = kol),
                          c(kol, X, "log2(Ratio)"))
    data2 <- as.data.frame(data2)
    data2$"Mean log10(Intensity)" <- avkol
    w1 <- is.all.good(data2$"Mean log10(Intensity)", 2)
    w2 <- is.all.good(data2$"log2(Ratio)", 2)
    w <- which(w1 & w2)
    if (length(w)) {
      data2 <- data2[w,]
      tst <- sapply(c(Y, X), function(x) { length(unique(data2[[x]])) })
      wrpKl <- c(Y, X)
      if (1 %in% tst) {
        wrpKl <- c(Y, X)[which(tst > 1)]
        #if ((length(wrpKl) == 1)&&(grepl(" ", wrpKl))) { wrpKl <- paste0("`", wrpKl, "`") }
      } else {
        if ((tst[1] >= tst[2]*3)||(tst[1] <= tst[2]/3)) {
          wrpKl <- paste0(paste0("`", X, "`"), "+", paste0("`", Y, "`"))
          tmp <- data2[, c(X, Y)]
          source(parSrc, local = FALSE)
          clusterExport(parClust, "tmp", envir = environment())
          data2[[wrpKl]] <- as.factor(parApply(parClust, data2[, c(X, Y)], 1, function(x) { paste(x, collapse = " ") }))
        } else {
          wrpKl <- c(Y, X)
        }
      }
      annot <- aggregate(data2$"log2(Ratio)", list(data2[[wrpKl]]), function(x) {
        x <- is.all.good(x)
        c(paste0("Median: ", signif(median(x), 3)), paste0("IQR: ", signif(IQR(x), 3)))
      })
      annot[, c("Median", "IQR")] <- do.call(as.data.frame, annot)
      annot$x <- NULL
      colnames(annot)[1:length(wrpKl)] <- wrpKl
      annot$Amax <- max(is.all.good(data2$"Mean log10(Intensity)"))*1.1
      annot$Amin <- min(is.all.good(data2$"Mean log10(Intensity)"))*1.1
      annot$Mmax <- max(is.all.good(data2$"log2(Ratio)"))*1.1
      annot$Mmin <- min(is.all.good(data2$"log2(Ratio)"))*1.1
      annot2 <- annot[, c(wrpKl, "Amax", "Mmin", "Mmax")] 
      annot2 <- rbind(annot2, annot2)
      annot2$Tag <- c(annot$Median, annot$IQR)
      data2[[Y]] <- factor(data2[[Y]], levels = get(Y))
      ttl <- "MA plot - PSMs"
      if (length(grps) > 1) { ttl <- paste0(ttl, " - ", grp) }
      ylim <- max(c(abs(c(annot$Mmax, annot$Mmin, (annot$Amax-annot$Amin)/4))))
      annot2$Y <- ylim*0.9
      w <- grep("^IQR: ", annot2$Tag)
      annot2$Y[w] <- -ylim*0.9
      if ("PTM-enrich." %in% colnames(data2)) {
        data2 <- data2[order(data2$"PTM-enrich."),]
        data2$"PTM-enrich." <- factor(data2$"PTM-enrich.", levels = c("", Modifs$`Full name`))
        data2$Alpha <- c(0.1, 1)[(data2$"PTM-enrich." != "")+1]
        myColors <- setNames(c("darkblue", "red"), c("Not-enriched/Flow-through", "Enriched"))
        plot <- ggplot(data2) +
          geom_scattermore(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`, colour = `PTM-enrich.`#, alpha = Alpha
                               ), shape = 16, size = 1#, alpha = 0.1
                           )
        if (length(unique(data$`PTM-enrich.`)) == 1) { plot <- plot + guides(colour = "none") }
      } else {
        plot <- ggplot(data2) + geom_scattermore(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`,
                                                     colour = `Parent sample`), size = 1#, alpha = 0.1
                                                 )
      }
      plot <- plot + scale_color_viridis(begin = 0.25, discrete = TRUE, option = "D") +
        geom_hline(yintercept = 0, colour = "grey") + xlab("A = mean log10(Intensity)") + ylab("M = sample log2(Ratio)") +
        geom_smooth(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`), color = "purple", formula = "y ~ s(x, bs = 'cs')", # Note that this fails sometimes, may be a package version issue
                    linewidth = 0.1, linetype = "dashed", method = 'gam') +
        geom_text(data = annot2, aes(x = Amax, y = Y, label = Tag), hjust = 1, cex = 2) +
        coord_fixed(log10(2)/log2(10)) + theme_bw() + ggtitle(ttl) + ylim(-ylim, ylim) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"), axis.text = element_text(size = 3),
              strip.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 7),
              strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 7))
      if (length(wrpKl) == 1) {
        if (grepl(" ", wrpKl)) { wrpKl <- paste0("`", wrpKl, "`") }
        plot <- plot + facet_wrap(as.formula(paste0(" ~ ", wrpKl)), drop = TRUE)
      } else {
        plot <- plot + facet_grid(as.formula(paste0("`", Y, "` ~ `", X, "`")), drop = TRUE)
      }
      #poplot(plot, 12, 22)
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, width = 10, height = 10, units = "in")
      ReportCalls <- AddPlot2Report(Space = FALSE)
    } else {
      msg <- paste0("Not enough valid data", c("", paste0(" for group ", grp))[(length(grps) > 1) + 1], " to draw an MA plot!")
      ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
    }
  }
  ReportCalls$Calls[[LRepCalls]] <- append(ReportCalls$Calls[[LRepCalls]],
                                           "body_add_par(Report, \"\", style = \"Normal\")")
} else {
  stop("Uh, I think you have the wrong analysis pipeline here...\nwhere are my sample groups and replicates!?!?!")
}

# Test for amino acid biases:
AA_biases %<o% AA_bias(Ev = ev, DB = db)
View(AA_biases)
write.csv(AA_biases, paste0(wd, "/Workflow control/AA_biases.csv"), row.names = FALSE)
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\"Amino acids frequency biases:\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_table(Report, AA_biases)")
ReportCalls <- AddSpace2Report()

if ((LabelType == "Isobaric")&&("Label.Purities.file" %in% colnames(Param))&&(!Param$Label.Purities.file %in% c("", " ", "NA", NA))&&(file.exists(Param$Label.Purities.file))) {
  if (!require("matlib", quietly = TRUE)) { install.packages("matlib") }
  cran_req <- c(cran_req, "matlib")
  Iso.purity <- read.csv(Param$Label.Purities.file)
  msg <- "Correct evidence reporter intensities for labels purity."
  ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " Reporter intensities were corrected for ", IsobarLab, " label purity factors.")
  if (colnames(Iso.purity)[1] == "New.format") {
    Iso.purity <- read.csv(Param$Label.Purities.file, check.names = FALSE)
    Iso.purity <- Iso.purity[, 2:ncol(Iso.purity)]
    Iso.purity$"Isobaric set" <- lapply(strsplit(as.character(Iso.purity$"Isobaric set"), ";"), as.integer)
    test <- sort(unique(unlist(Iso.purity$"Isobaric set")))
    if (exists("Iso")) {
      if (!sum(!Iso %in% test)) {
        test <- data.frame("Isobaric set" = test, check.names = FALSE)
        test$"Purity table" <- sapply(test$"Isobaric set", function(x) {
          list(Iso.purity[which(sapply(Iso.purity$"Isobaric set", function(y) { x %in% y })),
                          which(colnames(Iso.purity) != "Isobaric set")])
        })
        kol <- paste0("corr. ", ev.ref["Original"], get(IsobarLab))
        kol <- kol[which(paste0(ev.ref["Original"], get(IsobarLab)) %in% colnames(ev))]
        ev[, kol] <- NA
        for (i in Iso) { #i <- Iso[1]
          wI <- which(ev$"Raw file path" %in% Frac.map$"Raw file"[which(Frac.map$Isobaric.set == i)])
          e <- ev[wI,]
          pt <- test$"Purity table"[[which(test$"Isobaric set" == i)]]
          lb <- Exp.map$Isobaric.label[which(Exp.map$Isobaric.set == i)]
          lb2 <- Exp.map$Isobaric.label.details[which(Exp.map$Isobaric.set == i)]
          kol <- paste0(ev.ref[length(ev.ref)], lb)
          w <- which(kol %in% colnames(e))
          kol <- kol[w] ; lb <- lb[w] ; lb2 <- lb2[w]
          o <- order(as.numeric(lb))
          kol <- kol[o] ; lb <- lb[o] ; lb2 <- lb2[o]
          w <- match(lb2, pt$"Isobaric label details")
          pt <- pt[w,]
          w <- which(colnames(pt) != "Isobaric label details")
          A <- matrix(rep(0, length(lb)*(length(lb)+length(w)-1)), ncol = length(lb))
          for (j in 1:length(lb2)) {
            tmp <- as.numeric(pt[j, w])
            tmp[which(is.na(tmp))] <- 0
            A[j:(j+length(w)-1), j] <- tmp
          }
          w <- which(sapply(1:length(w), function(x) {
            sum(sapply(1:length(lb), function(y) { A[y+x-1, y] == 100 }))
          }) == length(lb))
          A <- A[(1:length(lb))+w-1,]
          source(parSrc, local = FALSE)
          exports <- list("A", "e", "kol")
          clusterExport(parClust, exports, envir = environment())
          clusterCall(parClust, function() library(matlib))
          clusterCall(parClust, function() library(proteoCraft))
          temp <- as.data.frame(t(parApply(parClust, e[,kol], 1, function(x) {
            b <- as.numeric(x)
            b[which(!is.all.good(b, 2))] <- 0
            sb <- sum(b)
            if (sb > 0) {
              #showEqn(round(A, 3), b)
              #res <- as.numeric(gsub(".+= +", "", Solve(A, b))) # I stopped using this since it seems to return approximations sometimes
              res <- solve(A, b)
              # There are cases where we will get negative values which we will have to truncate:
              res[which(res < 0)] <- 0
              res <- #round(
                res*sb/sum(res)
              #, 0)# I used to round here, but do not think it's necessary
              # print(res/b)
            } else { res <- rep(NA, length(kol)) }
            res <- setNames(res, kol)
            return(res)
          })))
          e[, gsub(topattern(ev.ref[length(ev.ref)]), paste0("corr. ", ev.ref["Original"]), kol)] <- temp[, kol]
          ev[wI,] <- e
        }
        ev.ref["Corrected"] <- paste0("corr. ", ev.ref["Original"])
        cat("Reporter Intensity correction done!\n")
      } else { warning("Some Isobaric sets are not defined in the labels purity table, skipping!") }
    } else { stop("There should be an \"Iso\" (Isobaric sets) object by now!") }
  } else {
    Iso.purity$Isobaric.set <- as.character(Iso.purity$Isobaric.set)
    test <- sort(unique(as.integer(unlist(strsplit(Iso.purity$Isobaric.set, ";")))))
    if (exists("Iso")) {
      if (!sum(!Iso %in% test)) {
        test <- test[which(test %in% Iso)]
        test <- data.frame(Isobaric.set = test)
        test$Purity.table <- sapply(test$Isobaric.set, function(x) {
          list(Iso.purity[which(sapply(lapply(strsplit(Iso.purity$Isobaric.set, ";"), as.integer), function(y) { x %in% unlist(y) })),
                          which(colnames(Iso.purity) != "Isobaric.set")])
        })
        kol <- paste0("corr. ", ev.ref["Original"], get(IsobarLab))
        kol <- kol[which(paste0(ev.ref["Original"], get(IsobarLab)) %in% colnames(ev))]
        ev[, kol] <- NA
        for (i in Iso) { #i <- Iso[1]
          wI <- which(ev$"Raw file path" %in% Frac.map$"Raw file"[which(Frac.map$Isobaric.set == i)])
          e <- ev[wI,]
          pt <- test$Purity.table[[which(test$Isobaric.set == i)]]
          lb <- Exp.map$Isobaric.label[which(Exp.map$Isobaric.set == i)]
          lb2 <- Exp.map$Isobaric.label.details[which(Exp.map$Isobaric.set == i)]
          kol <- paste0(ev.ref[length(ev.ref)], lb)
          w <- which(kol %in% colnames(e))
          kol <- kol[w] ; lb <- lb[w] ; lb2 <- lb2[w]
          o <- order(as.numeric(lb))
          kol <- kol[o] ; lb <- lb[o] ; lb2 <- lb2[o]
          pt <- pt[which(pt$Isobaric.label.details %in% lb2),]
          kol2 <- c("MI_minus.2", "MI_minus.1", "Monoisotopic", "MI_plus.1", "MI_plus.2")
          kol3 <- c("Who_minus.2", "Who_minus.1", "Who_plus.1", "Who_plus.2")
          pt[,kol2] <- sweep(pt[,kol2], 1, rowSums(pt[,kol2]), "/")
          A <- as.data.frame(matrix(rep(0, length(lb)*(length(lb)+4)), nrow = length(lb)))
          colnames(A) <- c("-2", "-1", lb2, "+1", "+2")
          kount <- 0
          for (l in lb2) { #l <- lb2[2]
            kount <- kount+1
            tmp <- unlist(pt[which(pt$Isobaric.label.details == l), kol3])
            wc <- wc2 <- which(tmp != "")
            wc2[which(wc2 > 2)] <- wc2[which(wc2 > 2)]+1
            fact <- pt[which(pt$Isobaric.label.details == l), kol2]
            fact <- fact[sort(c(3, wc2))]
            names(fact) <- c(l, tmp[wc])[order(c(3, wc2))]
            fact <- fact[which(names(fact) %in% colnames(A))]
            if (length(fact)) {
              mcola <- sapply(names(fact), function(x) {which(colnames(A) == x)})
              A[kount,mcola] <- fact
            }
          }
          A <- as.matrix(A[,3:(length(lb)+2)])
          exports <- list("A", "e", "kol")
          clusterExport(parClust, exports, envir = environment())
          clusterCall(parClust, function() library(matlib))
          clusterCall(parClust, function() library(proteoCraft))
          temp <- as.data.frame(t(parApply(parClust, e[,kol], 1, function(x) {
            b <- as.numeric(x)
            b[which(!is.all.good(b, 2))] <- 0
            sb <- sum(b)
            if (sb > 0) {
              #showEqn(round(A, 3), b)
              #res <- as.numeric(gsub(".+= +", "", Solve(A, b))) # I stopped using this since it seems to return approximations sometimes
              res <- solve(A, b)
              # There are cases where we will get negative values which we will have to truncate:
              res[which(res < 0)] <- 0
              res <- #round(
                res*sb/sum(res)
              #, 0)# I used to round here, but do not think it's necessary
              # print(res/b)
            } else { res <- rep(NA, length(kol)) }
            return(res)
          })))
          e[, gsub(topattern(ev.ref[length(ev.ref)]), paste0("corr. ", ev.ref["Original"]), kol)] <- temp
          ev[wI,] <- e
        }
        ev.ref["Corrected"] <- paste0("corr. ", ev.ref["Original"])
        cat("Reporter Intensity correction done!\n")
      } else { warning("Some Isobaric sets are not defined in the labels purity table, skipping!") }
    } else { stop("There should be an \"Iso\" (Isobaric sets) object by now!") }
  }
}
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

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Optional - Normalize evidence MS1 intensities, then, if applicable, MS2 reporter (Isobaric labelling) or fragment (DIA) intensities
# Step 0 for DIA measurements
# For DIA we have more accurate estimates of the ratio between samples from measurements of fragments (MS2-based quant)
# We will re-scale total precursor intensities across samples based on the relative amounts of each sample
#Param$Norma.Ev.Intens <- FALSE
if (Param$Norma.Ev.Intens) {
  msg <- paste0(evNm, "s-level normalisations:\n------------------------------------\n")
  ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
  # Define groups - this will ensure that, if phospho (or other) -enrichment took place, these peptides will be normalized separately
  ev$"Normalisation group" <- "Standard"
  if ("PTM-enriched" %in% colnames(Frac.map)) {
    # (In column "PTM enriched", use NA to indicate no enrichment!!!)
    if (!"PTM-enriched" %in% colnames(Frac.map)) { Frac.map$"PTM-enriched" <- NA }
    ptmChck <- unique(Frac.map$"PTM-enriched")
    ptmChck <- ptmChck[which(!is.na(ptmChck))]
    if (sum(!ptmChck %in% Modifs$`Full name`)) { stop("Some of the modifications in column \"PTM-enriched\" of Fractions map are invalid!") }
    # Here is what we want to do for those modifications:
    # - For enriched samples, keep only peptides with the target modification
    # - For other samples, remove all peptides with the modification which were also found in enriched samples
    if (length(ptmChck)) {
      for (ptm in ptmChck) { #ptm <- ptmChck[1]
        # Below "modified" means "modified with ptm" and "enriched" means "enriched for ptm"
        mrk <- Modifs$Mark[match(ptm, Modifs$`Full name`)]
        rw1 <- Frac.map$"Raw file"[which(Frac.map$"PTM-enriched" == ptm)]
        rw0 <- Frac.map$"Raw file"[which((is.na(Frac.map$"PTM-enriched"))|(Frac.map$"PTM-enriched" != ptm))]
        w1 <- ev$id[which(ev$"Raw file path" %in% rw1)] # PSMs from enriched runs
        ev$"Normalisation group"[match(w1, ev$id)] <- ptm
        w0 <- ev$id[which(ev$"Raw file path" %in% rw0)] # PSMs from non-enriched runs
        w2 <- ev$id[which(!ev$"Raw file path" %in% c(rw0, rw1))] # Any others
        w1a <- w1[grep(mrk, ev$"Modified sequence"[match(w1, ev$id)])] # Modified PSMs from enriched samples (i.e. what we were trying to enrich!)
        if (length(w1a)) {
          w0a <- w0[which(!ev$"Modified sequence"[match(w0, ev$id)] %in% unique(ev$"Modified sequence"[match(w1a, ev$id)]))] # Un-modified PSMs from non-enriched runs
          l1 <- length(w1)-length(w1a) # This is the number of un-modified PSMs we are removing from enriched runs
          l0 <- length(w0)-length(w0a) # This is the number of modified PSMs we are removing from non-enriched runs
          if (l1) {
            msg <- paste0("Removing ", l1, " peptide evidences without the ", ptm, " modification from ", ptm,
                          "-enriched raw files (", round(100*l1/length(w1), 2), "%)!")
            ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
          }
          if (l0) {
            msg <- paste0("Removing ", l0, " ", ptm, "-modified peptide evidences from un-enriched samples!")
            ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
          }
          ev <- ev[which(ev$id %in% c(w0a, w1a, w2)),]
        } else {
          msg <- paste0("Not a single ", ptm, "-modified evidence found in ", ptm, "-enriched raw files, investigate!")
          ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
        }
      }
    }
  }
  nms <- Norm.Groups$names
  w <- which(!nms %in% colnames(ev))
  if (length(w)) { ev[, nms[w]] <- Exp.map[match(ev$MQ.Exp, Exp.map$MQ.Exp), nms[w]] }
  ev$"Normalisation group" <- do.call(paste, c(ev[, c(nms, "Normalisation group")], sep = "_"))
  # Now, those groups are fine, but we also want to normalize per fraction at this stage so that each series of
  # equivalent fractions from different samples get normalized to each other.
  mrmgrps <- unique(ev$"Normalisation group")
  fr <- unique(ev$Fraction)
  tmp <- data.frame(Grp = as.character(sapply(mrmgrps, function(x) { rep(x, length(fr)) })),
                    Frac = as.character(rep(fr, length(mrmgrps))))
  tmp$Nm <- apply(tmp, 1, paste, collapse = "_")
  ev$"Normalisation group + Fraction" <- NA
  for (i in 1:nrow(tmp)) {
    w <- which((ev$"Normalisation group" == tmp$Grp[i])&(ev$Fraction == tmp$Frac[i]))
    ev$"Normalisation group + Fraction"[w] <- tmp$Nm[i]
  }
  Norma.Ev.Intens.Groups %<o% set_colnames(aggregate(ev$"Normalisation group + Fraction", list(ev$"Raw file path"), unique),
                                           c("Raw file", "Groups"))
  tst <- aggregate(Norma.Ev.Intens.Groups$"Raw file", list(Norma.Ev.Intens.Groups$Groups), length)
  w <- which(tst$x > 1)
  #tst <- aggregate(Norma.Ev.Intens.Groups$"Raw file", list(Norma.Ev.Intens.Groups$Groups), list)
  if (length(w)) {
    msg <- " - Normalizing MS1-level PSM intensities"
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
    if (length(w) > 1) { cat(" (per fraction/PTM enrichment group)\n") }
    msg <- "   - Classic normalisation to the median"
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
    Norma.Ev.Intens.Groups <- Norma.Ev.Intens.Groups[which(Norma.Ev.Intens.Groups$Groups %in% tst$Group.1[w]),]
    # (Per fractions X PTM enrichment group)
    # Step 1a:
    ev.col["Normalisation"] <- paste0("norm. ", ev.col["Original"])
    ev[[ev.col["Normalisation"]]] <- NA
    Norm.Ev %<o% data.frame(Group = unique(Norma.Ev.Intens.Groups$Groups))
    Grps2 <- MQ.Exp
    Grps2Kol <- "MQ.Exp"
    if (LabelType == "Isobaric") {
      Grps2 <- Iso
      Grps2Kol <- "Isobaric.set"
    }
    for (grp2 in Grps2) { Norm.Ev[[paste0("Grp", grp2)]] <- 1 }
    for (grp in Norm.Ev$Group) { #grp <- Norm.Ev$Group[1]
      r <- Norma.Ev.Intens.Groups$"Raw file"[which(Norma.Ev.Intens.Groups$Groups == grp)]
      wg <- which(ev$"Raw file path" %in% r)
      M <- 10^median(is.all.good(log10(unlist(ev[wg, ev.col["Original"]])))) # For preserving original scale
      #M <- 10^mlv(is.all.good(log10(unlist(w[wg, ev.col["Original"]]))), method = "Parzen")[1]
      for (grp2 in Grps2) { #grp2 <- Grps2[1]
        w2 <- which(ev[wg, Grps2Kol] == grp2)
        if (length(w2)) {
          m <- 10^median(is.all.good(log10(ev[wg[w2], ev.col["Original"]])))
          #m <- 10^mlv(is.all.good(log10(ev[wg[w2], ev.col["Original"]])), method = "Parzen")[1]
          ev[wg[w2], ev.col["Normalisation"]] <- ev[wg[w2], ev.col["Original"]]*M/m
          Norm.Ev[match(grp, Norm.Ev$Group), paste0("Grp", grp2)] <- m/M
        }
      }
    }
    if (("Adv.Norma.Ev.Intens" %in% colnames(Param))&&(Param$Adv.Norma.Ev.Intens != FALSE)) {
      msg <- c("   - Levenberg-Marquardt normalisation to finely align values..."#, 
               #"     (this can take some time)",
               #"     ..."
               )
      ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Print = FALSE)
      cat(paste0(msg, "\n", collapse = "\n"))
      txtAdv <- "Evidence MS1 intensities"
      ev.col["Advanced normalisation"] <- paste0("AdvNorm. ", ev.col["Original"])
      ev[[ev.col["Advanced normalisation"]]] <- NA
      AdvNorm.Ev %<o% data.frame(Group = unique(Norma.Ev.Intens.Groups$Groups))
      for (grp2 in Grps2) { AdvNorm.Ev[[paste0("Grp", grp2)]] <- 1 }
      for (grp in Norm.Ev$Group) { #grp <- Norm.Ev$Group[1]
        w <- which(Norma.Ev.Intens.Groups$Groups == grp)
        r <- Norma.Ev.Intens.Groups$"Raw file"[which(Norma.Ev.Intens.Groups$Groups == grp)]
        wg <- which(ev$"Raw file path" %in% r)
        tmp <- data.table(Uniq = ev$`Unique State`[wg], Grp = ev[wg, Grps2Kol], Int = ev[wg, ev.col["Normalisation"]])
        tmp <- tmp[, list(x = sum(Int, na.rm = TRUE)), keyby = list(Group.1 = Uniq, Group.2 = Grp)]
        tmp <- as.data.frame(tmp)
        tmp <- suppressMessages(reshape2::dcast(tmp, Group.1~Group.2))
        colnames(tmp) <- c("Unique State", paste0("Grp", colnames(tmp)[2:ncol(tmp)]))
        kol <- paste0("Grp", Grps2)
        w <- which(kol %in% colnames(tmp))
        kol <- kol[w] 
        tmp2 <- AdvNorm.IL(tmp, "Unique State", kol, FALSE, 5)
        m2 <- setNames(sapply(Grps2[w], function(x) {
          mean(tmp[[paste0("Grp", x)]]/tmp2[[paste0("AdvNorm.Grp", x)]], na.rm = TRUE)
        }), kol)
        for (grp2 in Grps2) {
          w2 <- which(ev[wg, Grps2Kol] == grp2)
          ev[wg[w2], ev.col["Advanced normalisation"]] <- ev[wg[w2], ev.col["Normalisation"]]/m2[paste0("Grp", grp2)]
        }
        AdvNorm.Ev[match(grp, AdvNorm.Ev$Group), kol] <- 1/m2
      }
      cat("     Done!\n")
    }
    cat("\n")
  } else {
    cat(" - Skipping MS1-level PSM intensity normalisation since all groups contain exactly 1 MS file.\n\n")
  }
  if (LabelType == "Isobaric") {
    msg <- " - Normalizing Reporter intensities"
    if (length(Iso) > 1) { msg <- c(msg, paste0("   (per ", IsobarLab, " sample)\n")) }
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
    msg <- " - Classic normalisation to the median"
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
    # Per combined sample
    er1 <- ev.ref["Normalisation"] <- paste0("Norm. ", ev.ref["Original"])
    er0 <- ev.ref[match("Normalisation", names(ev.ref))-1]
    k0 <- paste0(er0, get(IsobarLab))
    k1 <- paste0(er1, get(IsobarLab))
    ev[, k1] <- NA
    Norm.Ev.RepIntens %<o% data.frame(Group = Iso)
    for (ch in get(IsobarLab)) { Norm.Ev.RepIntens[[paste0("Channel_", ch)]] <- 1 }
    for (i in Iso) { #i <- Iso[1]
      wg <- which(ev$Isobaric.set == i)
      M3 <- 10^median(is.all.good(log10(unlist(ev[wg, k0])))) # For preserving original scale
      #M <- 10^mlv(is.all.good(log10(unlist(w[wg, ev.col["Original"]]))), method = "Parzen")[1]
      m3 <- sapply(get(IsobarLab), function(ch) { 10^median(is.all.good(log10(ev[wg, paste0(er0, ch)]))) })
      ev[wg, k1] <- sweep(ev[wg, k0], 2, M3/m3, "*")
      Norm.Ev.RepIntens[match(i, Norm.Ev.RepIntens$Group), paste0("Channel_", get(IsobarLab))] <- m3/M3
    }
    tstAdvNrm <- FALSE
    if (("Adv.Norma.Ev.Intens" %in% colnames(Param))&&(Param$Adv.Norma.Ev.Intens != FALSE)) {
      if (Param$Adv.Norma.Ev.Intens.Type == "C") {
        msg <- c("   - Levenberg-Marquardt normalisation...", "     (please wait)", "     ...")
        ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Print = FALSE)
        cat(paste0(msg, "\n", collapse = "\n"))
        # Per combined sample
        # Because this is computationally expensive, we are doing it per Fraction then averaging:
        # the values should be the same across fractions
        er1 <- ev.ref["Advanced normalisation"] <- paste0("AdvNorm. ", ev.ref["Original"])
        er0 <- ev.ref["Normalisation"]
        k0 <- paste0(er0, get(IsobarLab))
        k1 <- paste0(er1, get(IsobarLab))
        ev[, k1] <- NA
        AdvNorm.Ev.RepIntens %<o% data.frame(Group = Iso)
        tmpEv <- ev[, c("Isobaric.set", "Unique State", k0)]
        if ("Fraction" %in% colnames(ev)) { tmpEv$Fraction <- ev$Fraction } else { tmpEV$Fraction <- 1 }
        clusterCall(parClust, function() library(proteoCraft))
        m4 <- tstRI <- list()
        msg <- "     Estimating normalisation factors within..."
        ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Print = FALSE)
        for (i in Iso) {#i <- Iso[1]
          msg <- paste0("     ...", IsobarLab, " sample ", i, "...")
          ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Print = FALSE)
          wi <- which(tmpEv$Isobaric.set == i)
          tmp <- aggregate(tmpEv[wi, k0], list(tmpEv$`Unique State`[wi], tmpEv$Fraction[wi]), sum, na.rm = TRUE)
          colnames(tmp)[which(colnames(tmp) == "Group.1")] <- "Unique State"
          colnames(tmp)[which(colnames(tmp) == "Group.2")] <- "Fraction"
          Fr <- unique(tmp$Fraction)
          clusterExport(parClust, list("tmp", "k0", "Fr"), envir = environment())
          # Create normalized data
          tmp2 <- setNames(parLapply(parClust, Fr, function(fr) {
            AdvNorm.IL(tmp[which(tmp$Fraction == fr),], "Unique State", k0, FALSE, 5)
          }), paste0("Fr. ", Fr))
          # Compute normalisation factors
          tmp2F <- as.data.frame(sapply(Fr, function(fr) {
            sapply(1:length(get(IsobarLab)), function(x) {
              setNames(mean(tmp[which(tmp$Fraction == fr), k0[x]]/(tmp2[[paste0("Fr. ", fr)]][, paste0("AdvNorm.", k0[x])]),
                            na.rm = TRUE), paste0("Ch. ", get(IsobarLab)[x]))
            })
          }))
          colnames(tmp2F) <- paste0("Fr. ", Fr)
          stopifnot(!sum(rownames(tmp2F) != paste0("Ch. ", get(IsobarLab))))
          tstRI[[i]] <- tmp2F
          tmp2F$Label <- get(IsobarLab)
          # Number of valid values
          tmp2K <- as.data.frame(sapply(Fr, function(fr) {
            sapply(1:length(get(IsobarLab)), function(x) {
              length(is.all.good(log10(tmp[which(tmp$Fraction == fr), k0[x]])))
            })
          }))
          colnames(tmp2K) <- paste0("Fr. ", Fr)
          rownames(tmp2F) <- paste0("Ch. ", get(IsobarLab))
          tmp2K$Label <- get(IsobarLab)
          # Calculate weighted mean
          tmp2NormFact <- setNames(sapply(get(IsobarLab), function(x) {
            weighted.mean(tmp2F[match(x, tmp2F$Label), paste0("Fr. ", Fr)],
                          tmp2K[match(x, tmp2K$Label), paste0("Fr. ", Fr)])
          }), paste0("Ch ", get(IsobarLab)))
          tmp2NormFact[which(is.na(tmp2NormFact))] <- 1 # Better not normalize than corrupt data!
          m4[[i]] <- tmp2NormFact
        }
        # Apply results
        m4 <- as.data.frame(t(sapply(m4, unlist)))
        for (i in Iso) {
          wi <- which(tmpEv$Isobaric.set == i)
          ev[wi, k1] <- sweep(ev[wi, k0], 2, unlist(m4[which(rownames(m4) == i), paste0("Ch ", get(IsobarLab))]), "/")
        }
        AdvNorm.Ev.RepIntens[, paste0("Channel_", get(IsobarLab))] <- 1/m4[match(AdvNorm.Ev.RepIntens$Group, rownames(m4)), paste0("Ch ", get(IsobarLab))]
        tstAdvNrm <- TRUE
        cat("     Done!\n")
      } else {
        stop("Not implemented yet, I need to update the current AdvNorm function as it takes too long or crashes.")
      }
    }
    # Combine results of all reporter intensity normalisation steps
    Norm.Ev.RepIntens.All %<o% Norm.Ev.RepIntens
    rownames(Norm.Ev.RepIntens.All) <- Norm.Ev.RepIntens.All$Group; Norm.Ev.RepIntens.All$Group <- NULL
    Norm.Ev.RepIntens.All <- suppressMessages(reshape2::melt(Norm.Ev.RepIntens.All))
    Norm.Ev.RepIntens.All$Iso <- Iso
    Norm.Ev.RepIntens.All$Channel <- as.numeric(gsub("^Channel_", "", Norm.Ev.RepIntens.All$variable))
    Norm.Ev.RepIntens.All$Fraction <- "All"
    if (tstAdvNrm) {
      tmp2tst <- suppressMessages(melt(tstRI))
      colnames(tmp2tst)[which(colnames(tmp2tst) == "L1")] <- "Iso"
      tmp2tst$Channel <- get(IsobarLab)
      tmp2tst$Fraction <- as.numeric(gsub("^Fr\\. ", "", tmp2tst$variable))
      for (i in Iso) {
        for (ch in get(IsobarLab)) {
          w1 <- which((Norm.Ev.RepIntens.All$Iso == i)&(Norm.Ev.RepIntens.All$Channel == ch))
          w2 <- which((tmp2tst$Iso == i)&(tmp2tst$Channel == ch))
          tmp2tst$value[w2] <- tmp2tst$value[w2]*Norm.Ev.RepIntens.All$value[w1]
        }
      }
      tmp2tst$Fraction <- factor(tmp2tst$Fraction, levels = 1:max(tmp2tst$Fraction))
      Norm.Ev.RepIntens.All <- tmp2tst
    }
    Norm.Ev.RepIntens.All$Channel <- factor(Norm.Ev.RepIntens.All$Channel, levels = 0:max(Norm.Ev.RepIntens.All$Channel))
    Norm.Ev.RepIntens.All$Iso <- factor(paste0(IsobarLab, " sample ", Norm.Ev.RepIntens.All$Iso), levels = paste0(IsobarLab, " sample ", Iso))
    w <- (("Adv.Norma.Ev.Intens" %in% colnames(Param))&(Param$Adv.Norma.Ev.Intens != FALSE))+1
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " MS1 intensities and ", IsobarLab, " reporter intensities were re-normalized ",
                             c("to the median", "using the Levenberg-Marquardt procedure to minimize sample-to-sample differences")[w], ".")
    # Visualize results
    dir <- paste0(wd, "/Workflow control/", evNm, "s/Reporter intensities")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    ttl <- "Reporter intensities trend VS fractions"
    plot <- ggplot(Norm.Ev.RepIntens.All) +
      geom_tile(aes(x = Fraction, y = Channel, fill = value), linewidth = 1, width = 1) +
      scale_fill_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() + ylab(paste0(IsobarLab, " channel")) + ggtitle(ttl)
    if (length(Iso) > 1) { plot <- plot + facet_wrap(~Iso) }
    print(plot) # This type of QC plot does not need to pop up, the side panel is fine
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
  } else {
    w <- (("Adv.Norma.Ev.Intens" %in% colnames(Param))&(Param$Adv.Norma.Ev.Intens != FALSE))+1
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " MS1 intensities were re-normalized ",
                             c("to the median", "using the Levenberg-Marquardt procedure to minimize sample-to-sample differences")[w], ".")
  }
  if ((LabelType == "LFQ")&&(Param$Label == "DIA")&&("MS2_intensities" %in% colnames(ev))) {
    # If isobaric, re-scale MS2 intensities to apply MS1 normalisation factors
    # Important: this bit must remain after the normalisation of MS1 intensities
    DatAnalysisTxt <- gsub("\\.$", ", then individual MS2 intensities were re-scaled using MS1 intensities.", DatAnalysisTxt)
    kol <- gsub("Intensity", "MS2 intensities", ev.col[length(ev.col)])
    stopifnot(grepl("MS2 intensities", kol))
    ev[[kol]] <- ev$MS2_intensities # (Let's keep this as a numeric list)
    for (smpl in RSA$values) {
      mqe <- unlist(Exp.map$MQ.Exp[match(smpl, Exp.map$Ref.Sample.Aggregate)])
      w <- which(ev$MQ.Exp %in% mqe)
      mRt <- median(is.all.good(ev[w, ev.col[length(ev.col)]]/ev[w, ev.col["Original"]]))
      clusterExport(parClust, "mRt", envir = environment())
      ev[[kol]][w] <- parLapply(parClust, ev[[kol]][w], function(x) { x*mRt })
    }
    if (Param$Norma.Ev.Intens&&Param$Norma.Ev.Intens.show) {
      kol2 <- unique(c("id", "MQ.Exp", "MS2_intensities", kol))
      tst <- set_colnames(ev[, kol2], kol2)
      # Here it is easier to sum per row (otherwise this makes for very slow processing, creates a very huge table and plot, with little added value)
      kolz <- colnames(tst)[which(!colnames(tst) %in% c("id", "MQ.Exp"))]
      for (kl in kolz) {
        if (class(tst[[kl]]) != "list") { tst[[kl]] <- sapply(strsplit(tst[[kl]], ";"), as.numeric) }
        tst[[kl]] <- parSapply(parClust, tst[[kl]], sum)
      }
      tst <- reshape2::melt(tst, id.vars = c("id", "MQ.Exp"))
      tst$value <- log10(tst$value)
      tst$variable <- as.character(tst$variable)
      tst$Norm <- "Original"
      tst$Norm[which(tst$variable == kol)] <- "Normalised"
      tst$Norm <- factor(tst$Norm, levels = c("Original", "Normalised"))
      ttl <- "Evidence intensity normalisation"
      dir <- paste0(wd, "/Workflow control/", evNm, "s/Normalisation")
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      plot <- ggplot(tst) +
        geom_violin(aes(x = MQ.Exp, y = value, color = Norm, fill = Norm), alpha = 0.25) +
        geom_boxplot(aes(x = MQ.Exp, y = value, color = Norm, fill = Norm), alpha = 0.5) +
        scale_color_viridis_d(begin = 0.25) +
        scale_fill_viridis_d(begin = 0.25) +
        facet_wrap(~Norm, scales = "free") + theme_bw() + ggtitle(ttl) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      print(plot) # This type of QC plot does not need to pop up, the side panel is fine
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ReportCalls <- AddPlot2Report()
    }
  }
}
# If isobaric, re-scale reporter intensities to total evidence intensities:
# Important: this bit must remain after the normalisation of MS1 intensities
if (LabelType == "Isobaric") {
  DatAnalysisTxt <- gsub("\\.$", ", then reporter intensities were re-scaled using MS1 intensities.", DatAnalysisTxt)
  ev.ref["Adjusted"] <- paste0("adj. ", ev.ref["Original"])
  k0 <- paste0(ev.ref[match("Adjusted", names(ev.ref))-1], get(IsobarLab))
  k1 <- paste0(ev.ref["Adjusted"], get(IsobarLab))
  temp <- rowSums(ev[, k0], na.rm = TRUE)
  ev[, k1] <- sweep(ev[,k0], 1, ev[, ev.col[length(ev.col)]]/temp, "*")
  if (Param$Norma.Ev.Intens&&Param$Norma.Ev.Intens.show) {
    er0 <- ev.ref[match("Normalisation", names(ev.ref))-1]
    er1 <- ev.ref["Adjusted"]
    a0 <- paste0(er0, get(IsobarLab))
    a1 <- paste0(er1, get(IsobarLab))
    test <- ev[,c("MQ.Exp", a0, a1)]
    test <- reshape2::melt(test, id.vars = "MQ.Exp")
    test$Norm <- NA
    test$Norm[which(test$variable %in% a0)] <- "Original"
    test$Norm[which(test$variable %in% a1)] <- "Normalised"
    test$Norm <- factor(test$Norm, levels = c("Original", "Normalised"))
    test$value <- log10(test$value)
    ttl <- "Evidence intensity normalisation"
    dir <- paste0(wd, "/Workflow control/", evNm, "s/Normalisation")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    test$Channel <- as.numeric(gsub(topattern(c(er0, er1)), "", as.character(test$variable)))
    test$Channel <- factor(test$Channel, levels = sort(unique(test$Channel)))
    plot <- ggplot(test) +
      geom_violin(aes(x = Channel, y = value, color = Channel, fill= Channel), alpha = 0.25) +
      geom_boxplot(aes(x = Channel, y = value, color = Channel, fill = Channel), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(MQ.Exp ~ Norm, scales = "free", space = "free") + theme_bw() + ggtitle(ttl) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(plot) # This type of QC plot does not need to pop up, the side panel is fine
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
  }
  #Isobaric data: valid values
  kol <- grep(topattern(ev.ref["Original"]), colnames(ev), value = TRUE)
  kol <- grep(" count ", kol, value = TRUE, invert = TRUE)
  tst1 <- sapply(MQ.Exp, function(x) { sapply(gsub(topattern(ev.ref["Original"]), "", kol), function(y) {
    w <- which((Exp.map$MQ.Exp == x)&(Exp.map$Isobaric.label == y))
    if (length(w)) { res <- cleanNms(Exp.map$Ref.Sample.Aggregate[w]) } else { res <- "" }
    return(res)
  })})
  tst2 <- set_rownames(sapply(MQ.Exp, function(x) { sapply(kol, function(y) {
    sum(ev[which(ev$MQ.Exp == x), y] > 0)
  })}), rownames(tst1))
  tst3 <- set_rownames(sapply(MQ.Exp, function(x) { sapply(kol, function(y) {
    round(median(is.all.good(log10(ev[which(ev$MQ.Exp == x), y]))), 2)
  })}), rownames(tst1))
  tst <- rbind(rep("", length(MQ.Exp)),
               c("Samples", rep("", length(MQ.Exp)-1)),
               colnames(tst1),
               tst1,
               rep("", length(MQ.Exp)),
               c("Number of valid values", rep("", length(MQ.Exp)-1)), 
               colnames(tst2),
               tst2,
               rep("", length(MQ.Exp)),
               c("Median log10", rep("", length(MQ.Exp)-1)), 
               colnames(tst3),
               tst3,
               rep("", length(MQ.Exp)))
  colnames(tst) <- NULL
  data.table::fwrite(tst, paste0(wd, "/Workflow control/Valid values test.csv"), quote = FALSE, sep = ",", col.names = FALSE, na = "NA")
  #system(paste0("open \"", wd, "/Workflow control/Valid values test.csv\""))
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              "body_add_fpar(Report, fpar(ftext(\"Number of valid values per sample and channel:\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
  ReportCalls$Objects$Valid_values <- as.data.frame(tst)
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_table(Report, ReportCalls$Objects$Valid_values)")
  ReportCalls <- AddSpace2Report()
  ## To do here:
  ## Format table
  ## Add fraction-specificity
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Create modified peptides table
tmp <- as.character(ev$id)
tmp2 <- data.table(id = tmp, mqxp = ev$MQ.Exp, mod = ev$"Modified sequence")
tmp2 <- tmp2[order(ev$id, decreasing = FALSE),]
source(parSrc, local = FALSE)
exports <- list("tmp2", "Exp.map")
clusterExport(parClust, exports, envir = environment())
clusterCall(parClust, function() library(data.table))
tmp4 <- setNames(parLapply(parClust, c("ALLMYSAMPLESTUDUDUDUMMMDADA", RSA$values), function(i) {
  tmp3 <- copy(tmp2) # Because of how data.tables work! No idea how that plays with clusters...
  if (i != "ALLMYSAMPLESTUDUDUDUMMMDADA") {
    mqe <- unlist(Exp.map$MQ.Exp[which(Exp.map$Ref.Sample.Aggregate == i)])
    tmp3 <- tmp3[which(tmp3$mqxp %in% mqe), c("id", "mod")]
  }
  tmp3 <- as.data.frame(tmp3[, list(IDs = paste(id, collapse = ";")),
                             keyby = list(ModSeq = mod)])
  return(tmp3)
}), c("ALLMYSAMPLESTUDUDUDUMMMDADA", RSA$values))
pep %<o% set_colnames(tmp4[["ALLMYSAMPLESTUDUDUDUMMMDADA"]], c("Modified sequence", "Evidence IDs"))
for (i in RSA$values) {
  tmp <- tmp4[[i]]
  ki <- paste0("Evidence IDs - ", i)
  pep[[ki]] <- tmp$IDs[match(pep$"Modified sequence", tmp$ModSeq)]
  pep[which(is.na(pep[[ki]])), ki] <- ""
}
pep$id <- c(1:nrow(pep))
rvmtch2 <- match(pep$"Modified sequence", ev$"Modified sequence")
if ("Modified sequence_verbose" %in% colnames(ev)) {
  pep$"Modified sequence_verbose" <- ev$"Modified sequence_verbose"[rvmtch2]
}
pep$Sequence <- ev$Sequence[rvmtch2]
pep$Proteins <- ev$Proteins[rvmtch2]
tmp <- data.table(mod = ev$"Modified sequence", PEP = ev$PEP)
tmp <- as.data.frame(tmp[, list(PEP = min(PEP, na.rm = TRUE)), by = list(mod)])
pep$PEP <- tmp$PEP[match(pep$"Modified sequence", tmp$mod)]
mtch <- match(ev$Sequence, pep$Sequence)
mtch2 <- match(ev$"Modified sequence", pep$"Modified sequence")
# Amino Acid counts
sq <- pep$Sequence
clusterExport(parClust, "sq", envir = environment())
tmp <- parSapply(parClust, proteoCraft::AA, function(aa) {
  nchar(sq) - nchar(gsub(aa, "", sq))
})
colnames(tmp) <- paste0(colnames(tmp), " Count")
pep[, colnames(tmp)] <- tmp
ev[, paste0(AA, " Count")] <- pep[mtch, paste0(AA, " Count")]
for (aa in c("O", "U")) { # Only keep the selenocysteine and pyrrolysine amino acid columns if they are non-empty (NB: pyrrolysine should really only be in some bacteria)
  if (sum(ev[[paste0(aa, " Count")]]) == 0) {
    ev[[paste0(aa, " Count")]] <- NULL
    pep[[paste0(aa, " Count")]] <- NULL
  }
}
#
pep$Length <- nchar(pep$Sequence)
ev$Length <- pep$Length[mtch]
tmp <- data.table(mod = ev$"Modified sequence", Intensity = ev$Intensity)
tmp$Intensity[which(!is.all.good(tmp$Intensity, 2))] <- NA
w2 <- which(ev$MQ.Exp %in% unique(unlist(Exp.map$MQ.Exp[which(Exp.map$Use)])))
tmp2 <- copy(tmp)
tmp2 <- tmp2[w2, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(mod)]
pep$Intensity <- tmp2$Intensity[match(pep$"Modified sequence", tmp2$mod)]
if ((sum(!Exp.map$Use))||(length(unique(ev$MQ.Exp)) > length(unique(unlist(Exp.map$MQ.Exp))))) {
  tmp3 <- copy(tmp)
  tmp3 <- tmp3[, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(mod)]
  pep$"Intensity (all samples including unused)" <- tmp3$Intensity[match(pep$"Modified sequence", tmp3$mod)]
}
ev$"Peptide ID" <- pep$id[mtch2]
#
# Peptide name
tmp <- pep[, c("Proteins", "Modified sequence")]
w <- which(nchar(tmp$Proteins) > 13)
tmp$Proteins[w] <- paste0(substr(tmp$Proteins[w], 1, 12), "...")
tmp$"Modified sequence" <- gsub("_", "", tmp$"Modified sequence")
pep$"Peptide name" <- do.call(paste, c(tmp, sep = "\n"))
ev$"Peptide name" <- pep$"Peptide name"[mtch2]
#
# Modification site probabilities, score differences, site IDs, etc...
wPr <- which(paste0(Modifs$"Full name", " Probabilities") %in% colnames(ev))
if (length(wPr)) {
  for (i in wPr) { #i <- wPr[1]
    for (j in c(" Probabilities", " Score Diffs")) { #j <- " Probabilities"
      j1 <- paste0(Modifs$"Full name"[i], j)
      temp <- ev[, c("Modified sequence", j1, "PEP")]
      temp <- temp[which(temp[[j1]] != ""),]
      a0 <- unique(temp$"Modified sequence")
      clusterExport(parClust, list("temp", "j1"), envir = environment())
      b1 <- parLapply(parClust, a0, function(x) { temp[which(temp$"Modified sequence" == x), j1] })
      b2 <- parLapply(parClust, a0, function(x) { temp$PEP[which(temp$"Modified sequence" == x)] })
      l <- sapply(b1, function(x) { length(unique(x)) })
      wb <- which(l == 1)
      if (length(wb)) { b1[wb] <- sapply(b1[wb], function(x) {unique(x)}) } 
      wb <- which(l > 1)
      if (length(wb)) {
        temp1 <- cbind(b1[wb], b2[wb])
        clusterExport(parClust, "temp1", envir = environment())
        b1[wb] <- parApply(parClust, temp1, 1, function(x) {
          p <- 1-x[[2]]
          p <- p/sum(p)
          x <- x[[1]]
          s <- unlist(strsplit(unique(gsub("\\([^\\)]+\\)", "", x)), ""))
          an <- proteoCraft::Isapply(x, function(y) {
            as.numeric(proteoCraft::annot_to_tabl(y, Nterm = FALSE, Cterm = FALSE, numeric_data = TRUE)[[1]]$Annotations)
          })
          an <- sweep(an, 1, p, "*")
          n <- which(apply(an, 2, function(y) { length(proteoCraft::is.all.good(y)) }) == 0)
          an <- colSums(an, na.rm = TRUE)
          an[n] <- NA
          n <- which(!is.na(an))
          an <- paste0("(", round(an[n], 3), ")")
          x <- paste(proteoCraft::insertElems(s, n-1, an), collapse = "")
          return(x)
        })
      }
      b2 <- unlist(b1)
      if (length(b2) != length(b1)) { stop("Something went awry!") }
      pep[[j1]] <- ""
      wp <- which(pep$"Modified sequence" %in% a)
      pep[wp,j1] <- b2[match(pep$"Modified sequence"[wp], a)]
    }
    j1 <- paste0(j, " site IDs")
    pep[[j1]] <- NA
    w <- which(!is.na(ev[[j1]]))
    if (length(w)) {
      e <- ev[w,]
      temp <- aggregate(as.character(e[[j1]]), list(e$"Modified sequence"), function(x) {
        paste(sort(unlist(strsplit(x, ";"))), collapse = ";")
      })
      w1 <- which(pep$"Modified sequence" %in% e$"Modified sequence")
      pep[w1, j1] <- temp$x[match(pep$"Modified sequence"[w1], temp$Group.1)]
    }
  }
}

pep[, c("Reverse", "Potential contaminant")] <- ev[rvmtch2, c("Reverse", "Potential contaminant")]
if ((Param$Norma.Pep.Intens)||(Param$Norma.Pep.Ratio)) {
  if (Param$Norma.Ev.Intens) {
    tst <- aggregate(ev$"Normalisation group", list(ev$"Modified sequence"), unique)
    if ("character" %in% class(tst$x)) { # It could be that the same sequence is in different normalization groups - if so we do not want to re-use PSM-level normalisation groups!
      pep$"Normalisation group" <- ev$"Normalisation group"[rvmtch2]
    } else { pep$"Normalisation group" <- "Standard" }
  } else { pep$"Normalisation group" <- "Standard" }
  pep$MQ.Exp <- ev$MQ.Exp[match(pep$"Modified sequence", ev$"Modified sequence")]
  nms <- Norm.Groups$names
  w <- which(!nms %in% colnames(pep))
  if (length(w)) { pep[, nms[w]] <- Exp.map[match(pep$MQ.Exp, Exp.map$MQ.Exp), nms[w]] }
  pep$"Normalisation group" <- do.call(paste, c(pep[, c(nms, "Normalisation group")], sep = "_"))
}
if ("Quantity Quality" %in% colnames(ev)) { # Dia-NN specific:
  tmp <- data.table(Qual = ev$"Quantity Quality", ModSeq = ev$"Modified sequence")
  tmp <- tmp[, list(Qual= mean(Qual)), by = list(ModSeq)]
  tmp <- as.data.frame(tmp)
  pep$"Quantity Quality" <- tmp$Qual[match(pep$"Modified sequence", tmp$ModSeq)]
}
# Interesting note by Vadim on the balance between using more noisy peptides for quant vs fewer high quality ones:
#   https://github.com/vdemichev/DiaNN/issues/1102
# As I expected, they find out that using more is better than filtering.

#### Code chunk - Peptidoforms-level, calculate quantitation and test for outliers
## (future option, or when executing line-by-line: remove outliers/samples which will not be used)
## First visualize data: are there any clear outliers?
# Calculate single channel intensities and total intensity
#
# To do:
# - If outlier is the only reference sample in a reference group, unfortunately you will have to remove the whole group
# - Add Pearson correlation heatmap amongst visualizations to base decision to remove outliers, it is very good!
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
exports <- list("smpls", "Exp.map", "tmp", "pep.ref", "LabelType", "is.all.good")
if (LabelType == "Isobaric") {
  tmp <- ev[, c("MQ.Exp", "Modified sequence",
                paste0(ev.ref[length(ev.ref)], as.character(sort(as.numeric(unique(Exp.map$Isobaric.label))))))]
  exports <- append(exports, "ev.ref")
}
if (LabelType == "LFQ") {
  tmp <- ev[, c("MQ.Exp", "Modified sequence", ev.col[length(ev.col)])]
  exports <- append(exports, "ev.col")
}
smpls <- unique(Exp.map$Ref.Sample.Aggregate[which(Exp.map$Use)])
clusterExport(parClust, exports, envir = environment())
clusterCall(parClust, function() library(data.table))
tmp4 <- setNames(parLapply(parClust, smpls, function(smpl) { #smpl <- smpls[1]
  m <- match(smpl, Exp.map$Ref.Sample.Aggregate)
  mqe <- unlist(Exp.map$MQ.Exp[m])
  w2 <- which(tmp$MQ.Exp %in% mqe)
  tmp2 <- data.frame(mod = NA, Intensity = NA)
  if (length(w2)) {
    if (LabelType == "Isobaric") {
      # Not tested! If buggy, look at older versions from before the Feb. 2024 rewrite!
      j <- as.character(sort(as.numeric(Exp.map$Isobaric.label[m])))
      tmp3 <- tmp[, paste0(ev.ref[length(ev.ref)], j)]
      for (k in j) {
        kk <- paste0(ev.ref[length(ev.ref)], j)
        tmp3[which(!is.all.good(tmp3[[kk]], 2)), kk] <- NA
      }
      if (length(j) > 1) { tmp3 <- apply(tmp3, 1, sum, na.rm = TRUE) } # Ultra-rare cases where the same parent sample is in different isobaric channels in different fractions
      tmp2 <- data.table(mod = tmp[w2, "Modified sequence"],
                         Intensity = tmp3)
    }
    if (LabelType == "LFQ") {
      tmp2 <- data.table(mod = tmp[w2, "Modified sequence"],
                         Intensity = tmp[w2, ev.col[length(ev.col)]])
      tmp2$Intensity[which(!is.all.good(tmp2$Intensity, 2))] <- NA
    }
    tmp2 <- tmp2[, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(mod)]
    tmp2 <- as.data.frame(tmp2)
  }
  return(tmp2)
}), smpls)
for (smpl in smpls) {
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
    pv1 <- paste0("Original: ", paste(sapply(1:length(pv1), function(x) {
      paste0("PC", x, ": ", pv1[x], "%")
    }), collapse = ", "))
    w <- which(sapply(VPAL$names, function(x) { length(unique(scores1[[x]])) }) > 1)
    w <- w[which(tolower(substr(names(w), 1, 3)) != "rep")]
    scores1$Samples_group <- do.call(paste, c(scores1[, VPAL$names[w], drop = FALSE], sep = " "))
    scores1$Label <- apply(scores1[, RSA$names, drop = FALSE], 1, paste, collapse = "-")
    ttl <- "PCA plot - Preliminary - peptide level"
    plot <- ggplot(scores1) +
      geom_point(aes(x = PC1, y = PC2, colour = Samples_group, shape = Use)) +
      scale_color_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() +
      geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
      ggtitle(ttl, subtitle = pv1) +
      geom_text_repel(aes(x = PC1, y = PC2, label = Label, colour = Samples_group),
                      size = 2.5, show.legend = FALSE)
    #poplot(plot)
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
    Symb <- rep(c("circle", "diamond", "square", "cross", "x"), max(as.numeric(Rep)))[1:max(as.numeric(Rep))]             
    Symb <- Symb[as.numeric(scores1$Replicate)]
    # Custom color scale
    scores1$"Samples group" <- factor(scores1$Samples_group)
    if ("PC3" %in% colnames(scores1)) {
      plot_lyPCA <- plot_ly(scores1, x = ~PC1, y = ~PC2, z = ~PC3,
                            text = ~Label, type = "scatter3d", mode = "markers",
                            color = ~`Samples group`, colors = "viridis",
                            symbol = I(Symb))
    } else {
      plot_lyPCA <- plot_ly(scores1, x = ~PC1, y = ~PC2,
                            text = ~Label, type = "scatter", mode = "markers",
                            color = ~`Samples group`, colors = "viridis",
                            symbol = I(Symb))
    }
    plot_lyPCA %<o% layout(plot_lyPCA, title = ttl)
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
        temp$tst1 <- sapply(temp$Sequence, function(x) { nchar(unlist(strsplit(Seq, x))[1]) })
        temp$tst2 <- sapply(temp$Sequence, nchar)
        temp$tst3 <- sapply(temp$"Modified sequence", nchar)
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
            ggplot2::ggtitle(paste0(ttlRoot, ", ", c("mean-normalised", "Z-scored")[(meth == "ZSc")+1]),
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
          ggplot2::ggsave(paste0(dstDir, "/", ttl, ".jpeg"), heatmap.plot,
                          dpi = 600, width = 20, height = 12, units = "in")
          ggplot2::ggsave(paste0(dstDir, "/", ttl, ".pdf"), heatmap.plot,
                          dpi = 600, width = 20, height = 12, units = "in")
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
    ggsave(paste0(dir[1], "/", ttl, ".jpeg"), plot, dpi = 150)
    ggsave(paste0(dir[1], "/", ttl, ".pdf"), plot, dpi = 150)
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
if (Impute) {
  DatAnalysisTxt <- paste0(DatAnalysisTxt,
                           " Missing values were imputed using two different strategies: i) the KNN (K-Nearest Neighbours) method for Missing-At-Random values within sample groups, and ii) the QRICL (Quantile Regression Imputation of Left-Censored data) method for Missing-Not-At-Random values.")
  msg <- "Imputing missing values..."
  ReportCalls <- AddMsg2Report(Space = FALSE)
  ref <- pep.ref[length(pep.ref)]
  if ("Imputation" %in% names(pep.ref)) { ref <- pep.ref[match("Imputation", names(pep.ref))-1] }
  kol <- grep(topattern(ref), colnames(pep), value = TRUE)
  groups <- Exp.map[match(gsub(topattern(ref), "", kol), Exp.map$Ref.Sample.Aggregate), VPAL$column]
  temp <- Data_Impute2(pep[,kol], groups, is.log = FALSE)
  temp2 <- temp$Imputed_data
  colnames(temp2) <- gsub(topattern(ref), paste0("imput. ", pep.ref["Original"]), colnames(temp2))
  pep[, colnames(temp2)] <- temp2
  pep.ref["Imputation"] <- paste0("imput. ", pep.ref["Original"])
  rm(list = ls()[which(!ls() %in% .obj)])
  Script <- readLines(ScriptPath)
}

#### Code chunk - Re-normalize peptide intensities
# Check variance/intensity dependency before normalisation
rfnm <- "Original"
ttl <- "Peptides Variance vs Intensity dependency"
kol <- grep(topattern(pep.ref[rfnm]), colnames(pep), value = TRUE)
temp <- pep[, kol]
w <- which(apply(temp, 1, function(x) { length(is.all.good(x)) }) > 0)
temp <- temp[w,]
Aggr <- VPAL
LocAnalysis2 %<o% FALSE
if (LocAnalysis) { Aggr <- parse.Param.aggreg("Exp;Com") }
tst <- sapply(Aggr$values, function(x) {
  em <- Exp.map[which(Exp.map[[Aggr$column]] == x),]
  kl <- paste0(pep.ref[rfnm], em$Ref.Sample.Aggregate)
  return(rowMeans(temp[, kl, drop = FALSE], na.rm = TRUE))
})
colnames(tst) <- cleanNms(colnames(tst))
tst <- apply(tst, 1, function(x) { colnames(tst)[which(x == max(x))][1] })
tst2 <- aggregate(ev$Charge, list(ev$`Modified sequence`), function(x) { round(mean(x)) })
temp2 <- data.frame("log10(Mean)" = log10(rowMeans(temp, na.rm = TRUE)),
                    "log10(Variance)" = log10(apply(temp, 1, var, na.rm = TRUE)),
                    "Strongest in..." = tst,
                    check.names = FALSE)
temp2$"Main charge" <- paste0("Z = ", tst2$x[match(pep$`Modified sequence`[w], tst2$Group.1)])
temp2$"Main charge" <- factor(temp2$"Main charge", levels = paste0("Z = ", as.character(1:8)))
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
plot <- ggplot(temp2) +
  geom_scattermore(aes(x = `log10(Mean)`, y = `log10(Variance)`, colour = `Strongest in...`),
                   size = 1, alpha = 1) + ggtitle(ttl, subtitle = rfnm) + coord_fixed(0.3) + theme_bw() +
  facet_grid(`Strongest in...` ~ `Main charge`) + theme(strip.text.y.right = element_text(angle = 0)) +
  scale_colour_viridis_d(begin = 0.25)
#poplot(plot, 12, 22)
ggsave(paste0(dir, "/", ttl, " - ", rfnm, ".jpeg"), plot, width = 10, height = 10, units = "in", dpi = 300)
ggsave(paste0(dir, "/", ttl, " - ", rfnm, ".pdf"), plot, width = 10, height = 10, units = "in", dpi = 300)
ReportCalls <- AddPlot2Report(Title = paste0(ttl, " - ", rfnm))

if (Param$Norma.Pep.Intens) {
  if ("Batch.correction_PTM.subgroups" %in% colnames(Param)) {
    warning("Argument \"Batch.correction_PTM.subgroups\" is deprecated!")
  }
  NormGrps %<o% setNames(aggregate(pep$id, list(pep$"Normalisation group"), list), c("Group", "IDs"))
  nrmDr <- paste0(wd, "/Workflow control/Peptides/Intensities")
  pep_intens_norm <- list()
  shapenorm %<o% FALSE
  irsnorm %<o% FALSE
  irsok %<o% FALSE
  advnorm %<o% FALSE
  combatnorm %<o% FALSE
  pep.ref["Classic normalisation"] <- paste0("norm. ", pep.ref["Original"])
  if (("Norma.Pep.Intens.Shape" %in% colnames(Param))&&(toupper(Param$Norma.Pep.Intens.Shape) %in% c("VSN", "LOESS"))) {
    shapenorm %<o% TRUE
    shapenormtype %<o% toupper(Param$Norma.Pep.Intens.Shape)
    pep.ref[paste0(shapenormtype, " normalisation")]  <- paste0(tolower(shapenormtype), ". ", pep.ref["Original"])
  }
  if ((LabelType == "Isobaric")&&("Norma.Pep.Intens.IRS" %in% colnames(Param))&&(Param$Norma.Pep.Intens.IRS == TRUE)) {
    irsnorm %<o% TRUE
    irsok %<o% FALSE
    pep.ref["IRS normalisation"] <- paste0("irs. ", pep.ref["Original"])
  }
  if (("Adv.Norma.Pep.Intens" %in% colnames(Param))&&(Param$Adv.Norma.Pep.Intens == TRUE)) {
    advnorm %<o% TRUE
    pep.ref["Advanced (Levenberg-Marquardt) normalisation"] <- paste0("AdvNorm. ", pep.ref["Original"])
  }
  if (("Batch.correction" %in% colnames(Param))&&(!as.character(Param$Batch.correction) %in% c("", "F", "FALSE"))) {
    combatnorm %<o% TRUE
    pep.ref["SVA batch correction"] <- paste0("sva. ", pep.ref["Original"])
  }
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " Peptidoform intensity values were re-normalized as follows:")
  # Step 0:
  # Let's create the data frame which we will be normalizing;
  # For ease of processing, we will be converting it to log10:
  # (no need to do this by group)
  cat(" -> log10 transformation\n")
  currstep <- names(pep.ref)[match("Classic normalisation", names(pep.ref))-1]
  currref <- pep.ref[currstep]
  kol <- paste0(currref, RSA$values)
  kol <- kol[which(kol %in% colnames(pep))]
  temp_norm <- pep[,kol]
  for (i in kol) { temp_norm[[i]] <- log10(as.numeric(temp_norm[[i]])) }
  pepAllGood <- set_colnames(as.data.frame(apply(temp_norm[,kol], 1, function(x) { length(is.all.good(x)) }) > 0),
                             currstep)
  if (sum(!pepAllGood[[currstep]]) > 0) { warning("Quantitative peptide loss!") }
  pep_intens_norm[[currstep]] <- temp_norm
  # View(pep_intens_norm[[currstep]])
  # Step 1:
  # Classic normalisation by the median for each MQ.Exp:
  # This is fast and usually a good start:
  TxtSteps <- c("median-normalized")
  currstep <- "Classic normalisation"
  laststep <- names(pep.ref)[match(currstep, names(pep.ref))-1]
  wAG <- which(pepAllGood[[laststep]])
  temp_norm <- pep_intens_norm[[laststep]]
  temp_norm2 <- temp_norm*NA
  colnames(temp_norm2) <- gsub(topattern(pep.ref[laststep]), pep.ref[currstep], colnames(temp_norm2))
  for (grp in NormGrps$Group) { #grp <- NormGrps$Group[1]
    grpMtch <- match(NormGrps$IDs[[match(grp, NormGrps$Group)]], pep$id[wAG])
    grpMtch <- grpMtch[which(!is.na(grpMtch))]
    M <- median(is.all.good(unlist(temp_norm[grpMtch,])))
    #M <- mlv(unlist(temp_norm[grpMtch,]), method = "Parzen")[1]
    for (i in RSA$values) { #i <- RSA$values[1]
      k1 <- paste0(pep.ref[laststep], i)
      if (k1 %in% colnames(temp_norm)) {
        k2 <- paste0(pep.ref[currstep], i)
        m <- median(is.all.good(temp_norm[grpMtch, k1]))
        #m <- mlv(is.all.good(pep[[k1]]), method = "Parzen")[1]
        temp_norm2[grpMtch, k2] <- temp_norm[grpMtch, k1]+M-m
      }
    }
  }
  pepAllGood[[currstep]] <- apply(temp_norm2, 1, function(x) { length(is.all.good(x)) }) > 0
  tst1 <- sum(pepAllGood[[laststep]])
  tst2 <- sum(pepAllGood[[currstep]])
  if (tst2 < tst1) { warning("Quantitative peptide loss!") }
  pep_intens_norm[[currstep]] <- temp_norm2
  cat(" ->", currstep, "done!\n")
  # View(pep_intens_norm[[currstep]])
  # Step 2:
  # Some normalisations are meant to correct for some systematic errors, e.g. LOESS or VSN.
  # (I call them "shape" normalisations because of how they affect the shape of an MA plot.)
  if (shapenorm) {
    TxtSteps <- unique(c(TxtSteps, paste0("corrected for intensity-related variance biases using ", shapenormtype)))
    currstep <- paste0(shapenormtype, " normalisation")
    laststep <- names(pep.ref)[match(currstep, names(pep.ref))-1]
    wAG <- which(pepAllGood[[laststep]])
    cran_req <- unique(c(cran_req, "hexbin"))
    if (!require("hexbin", character.only = TRUE, quietly = TRUE)) { install.packages("hexbin") }
    require("hexbin", character.only = TRUE)
    if (shapenormtype == "LOESS") { pac <- "affy" }
    if (shapenormtype == "VSN") { pac <- "vsn" }
    bioc_req <- unique(c(bioc_req, pack))
    biocInstall(pack)
    temp_norm2 <- pep_intens_norm[[laststep]][wAG,]
    Qkol <- colnames(temp_norm2)
    temp_norm2$id <- pep$id[wAG]
    temp_norm2$Group <- pep$"Normalisation group"[wAG]
    w <- which(apply(temp_norm2[,Qkol], 1, function(x) { length(is.all.good(x)) == length(Qkol) }))
    A <- apply(temp_norm2[,Qkol], 1, function(x) { mean(is.all.good(x)) })
    # Visualize
    temp_plot <- reshape2::melt(temp_norm2, id.vars = c("id", "Group"))
    temp_plot$A <- rep(A, length(Qkol))
    temp_plot <- temp_plot[which(is.all.good(temp_plot$value, 2)),]
    temp_plot$M <- temp_plot$value - temp_plot$A
    temp_plot$variable <- gsub_Rep(topattern(pep.ref[laststep]), "", temp_plot$variable)
    u <- unique(cleanNms(temp_plot$variable))
    nu <- cleanNms(u)
    temp_plot$variable <- nu[match(x, u)]
    annot <- data.frame(variable = unique(temp_plot$variable))
    annot$Median <- sapply(annot$variable, function(x) {
      paste0("Median: ", round(median(is.all.good(temp_plot$M[which(temp_plot$variable == x)])), 3))
    })
    annot$IQR <- sapply(annot$variable, function(x) {
      paste0("IQR: ", round(IQR(is.all.good(temp_plot$M[which(temp_plot$variable == x)])), 3))
    })
    annot$Amax <- max(is.all.good(temp_plot$A))*1.1
    annot$Amin <- min(is.all.good(temp_plot$A))*1.1
    annot$Mmax <- max(is.all.good(temp_plot$M))*1.1
    annot$Mmin <- min(is.all.good(temp_plot$M))*1.1
    annot2 <- annot[,c("variable", "Amax", "Mmin", "Mmax")] 
    annot2 <- rbind(annot2, annot2)
    annot2$Label <- c(annot$Median, annot$IQR)
    ttl <- paste0("MA plot - ", shapenormtype, " normalisation_before")
    ylim <- max(c(abs(c(annot$Mmax, annot$Mmin, (annot$Amax-annot$Amin)/4))))
    annot2$Y <- ylim*0.9
    w <- grep("^IQR: ", annot2$Label)
    annot2$Y[w] <- -ylim*0.9
    l1 <- length(unique(temp_plot$variable))
    l2 <- length(unique(temp_plot$Group))
    nkol <- max(c(1, round(sqrt(l1*l2))))
    if ((l2 > 1)&&((nkol %% l2) != 0)) { nkol <- ceiling(nkol/l2)*l2 }
    while (nkol > l1*l2) { nkol <- nkol-1 }
    plot <- ggplot(temp_plot) +
      geom_scattermore(aes(x = A, y = M, colour = Group), size = 1, alpha = 1) +
      geom_hline(yintercept = 0, colour = "grey") + geom_smooth(aes(x = A, y = M), color = "red", linewidth = 0.8, linetype = "dashed") +
      geom_text(data = annot2, aes(x = Amax, y = Y, label = Label), hjust = 1, cex = 2) +
      scale_color_viridis_d(begin = 0.25) +
      facet_wrap(~variable+Group, ncol = nkol) + coord_fixed(log10(2)) + theme_bw() + ggtitle(ttl)
    print(plot) # This type of QC plot does not need to pop up, the side panel is fine
    dir <- paste0(nrmDr, "/", currstep)
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
    ImpGrps <- Exp.map[match(gsub(topattern(pep.ref[laststep]), "", Qkol), Exp.map$Ref.Sample.Aggregate),
                       VPAL$column]
    temp_norm3 <- temp_norm2[, Qkol]*NA
    Pos <- matrix(rep(FALSE, nrow(temp_norm2)*length(Qkol)), ncol = length(Qkol))
    for (grp in NormGrps$Group) {
      w <- which(temp_norm2$Group == grp)
      # Imputation to allow normalisation
      # We will fill gaps with normal distribution-based simulated data;
      # For each row, the mean of the distribution will be the mean of the row;
      # and the sd will be based on a LOESS regression estimate (nearest neighbours if missing):
      test <- Data_Impute2(temp_norm2[w, Qkol], ImpGrps)
      temp2 <- test$Imputed_data
      # Normalisation proper:
      if (shapenormtype == "LOESS") { tmp_nrm3 <- normalizeCyclicLoess(as.matrix(temp2)) }
      if (shapenormtype == "VSN") { tmp_nrm3 <- justvsn(as.matrix(10^temp2))/log2(10) }
      temp_norm3[w,] <- tmp_nrm3
      Pos[w,] <- test$Positions_Imputed
    }
    sd1 <- meanSdPlot(as.matrix(temp2))
    sd2 <- meanSdPlot(as.matrix(temp_norm3))
    ttl1 <- paste0("MA plot - ", shapenormtype, " mean SD plot_before")
    ttl2 <- paste0("MA plot - ", shapenormtype, " mean SD plot_after")
    plot1 <- sd1$gg + theme_bw() + ggtitle(ttl1)
    plot2 <- sd2$gg + theme_bw() + ggtitle(ttl2)
    #poplot(plot1, 12, 20)
    print(plot2) # This type of QC plot does not need to pop up, the side panel is fine
    ggsave(paste0(dir, "/", ttl1, ".jpeg"), plot1, dpi = 150, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl2, ".jpeg"), plot2, dpi = 150, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl1, ".pdf"), plot1, dpi = 150, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl2, ".pdf"), plot2, dpi = 150, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report(Plot = plot1, Title = ttl1, Space = FALSE)
    ReportCalls <- AddPlot2Report(Plot = plot2, Title = ttl2)
    # Remove imputed data:
    w <- which(Pos, arr.ind = TRUE)
    if (nrow(w)) { temp_norm3[w] <- temp_norm2[,Qkol][w] }
    colnames(temp_norm3) <- gsub(topattern(pep.ref[laststep]), pep.ref[currstep], Qkol)
    An <- apply(temp_norm3, 1, function(x) { mean(is.all.good(x)) })
    # Visualize:
    temp_plot <- temp_norm3
    temp_plot[, c("id", "Group")] <- temp_norm2[, c("id", "Group")]
    temp_plot <- reshape2::melt(temp_plot, id.vars = c("id", "Group"))
    temp_plot$A <- rep(An, length(Qkol))
    temp_plot <- temp_plot[which(is.all.good(temp_plot$value, 2)),]
    temp_plot$M <- temp_plot$value - temp_plot$A
    temp_plot$variable <- gsub_Rep(topattern(pep.ref[currstep]), "", temp_plot$variable)
    temp_plot$variable <- cleanNms(temp_plot$variable)
    annot <- data.frame(variable = unique(temp_plot$variable))
    annot$Median <- sapply(annot$variable, function(x) {
      paste0("Median: ", round(median(is.all.good(temp_plot$M[which(temp_plot$variable == x)])), 3))
    })
    annot$IQR <- sapply(annot$variable, function(x) {
      paste0("IQR: ", round(IQR(is.all.good(temp_plot$M[which(temp_plot$variable == x)])), 3))
    })
    annot$Amax <- max(is.all.good(temp_plot$A))*1.1
    annot$Amin <- min(is.all.good(temp_plot$A))*1.1
    annot$Mmax <- max(is.all.good(temp_plot$M))*1.1
    annot$Mmin <- min(is.all.good(temp_plot$M))*1.1
    annot2 <- annot[,c("variable", "Amax", "Mmin", "Mmax")] 
    annot2 <- rbind(annot2, annot2)
    annot2$Label <- c(annot$Median, annot$IQR)
    ttl <- paste0("MA plot - ", shapenormtype, " normalisation_after")
    ylim <- max(c(abs(c(annot$Mmax, annot$Mmin, (annot$Amax-annot$Amin)/4))))
    annot2$Y <- ylim*0.9
    w <- grep("^IQR: ", annot2$Label)
    annot2$Y[w] <- --ylim*0.9
    plot <- ggplot(temp_plot) +
      geom_scattermore(aes(x = A, y = M, colour = Group), size = 1, alpha = 1) +
      geom_hline(yintercept = 0, colour = "grey") + geom_smooth(aes(x = A, y = M), color = "red", linewidth = 0.8, linetype = "dashed") +
      geom_text(data = annot2, aes(x = Amax, y = Y, label = Label), hjust = 1, cex = 2) +
      scale_color_viridis_d(begin = 0.25) +
      facet_wrap(~variable+Group, ncol = nkol) + coord_fixed(log10(2)) + theme_bw() + ggtitle(ttl)
    print(plot) # This type of QC plot does not need to pop up, the side panel is fine
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
    temp_norm2 <- pep_intens_norm[[laststep]]*NA
    temp_norm2[wAG,] <- temp_norm3
    colnames(temp_norm2) <- colnames(temp_norm3)
    pepAllGood[[currstep]] <- apply(temp_norm2, 1, function(x) { length(is.all.good(x)) }) > 0
    tst1 <- sum(pepAllGood[[laststep]])
    tst2 <- sum(pepAllGood[[currstep]])
    if (tst2 < tst1) { warning("Quantitative peptide loss!") }
    pep_intens_norm[[currstep]] <- temp_norm2[, grep(topattern(pep.ref[currstep]), colnames(temp_norm2))]
    # View(pep_intens_norm[[currstep]])
    cat(" ->", currstep, "done!\n")
  }
  # Step 3:
  # Adapted from: https://pwilmart.github.io/IRS_normalisation/understanding_IRS.html
  # NB: Since we have log10 data here already, this has been rewritten for log-transformed data.
  # IRS is done per row so I do not think that the "Normalisation group" column are relevant here.
  if (irsnorm) {
    if (length(Iso) <= 1) {
      if (length(Iso) == 1) {
        warning("Skipping IRS normalisation: there is only 1 Isobarically labelled sample (aka \"Isobaric Set\")")
        pep.ref <- pep.ref[which(names(pep.ref) != "IRS normalisation")]
      } else { stop("\"LabelType\" is \"Isobaric\" but \"Iso\" has invalid length, investigate!") }
    } else {
      currstep <- "IRS normalisation"
      laststep <- names(pep.ref)[match(currstep, names(pep.ref))-1]
      wAG <- which(pepAllGood[[laststep]])
      dir <- paste0(nrmDr, "/IRS normalisation")
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      temp_norm2 <- pep_intens_norm[[laststep]]
      temp_norm3 <- temp_norm2[wAG,]
      design <- data.frame(Columns = colnames(temp_norm3))
      design$Iso <- sapply(gsub(topattern(pep.ref[laststep]), "", design$Columns), function(x) {
        Exp.map$Isobaric.set[which(Exp.map$Ref.Sample.Aggregate == x)]
      })
      kolmn <- VPAL$column
      design[[kolmn]] <- sapply(gsub(topattern(pep.ref[laststep]), "", design$Columns), function(x) {
        unique(Exp.map[which(Exp.map$Ref.Sample.Aggregate == x), kolmn])
      })
      design$Original.order <- 1:nrow(design)
      design <- design[order(design[[kolmn]]),]
      design <- design[order(design$Iso),]
      design$New.order <- 1:nrow(design)
      #tst <- aggregate(design$Iso, list(design[[kolmn]]), length); stopifnot(min(tst$x) == length(Iso))
      #tst <- aggregate(design[[kolmn]], list(design$Iso), length); stopifnot(min(tst$x) == length(A$values))
      temp_norm3 <- temp_norm3[, design$Columns]
      kol <- colnames(temp_norm3)
      Columns <- sapply(unique(design$Iso), function(x) { #x <- unique(design$Iso)[1]
        list(lapply(unique(design[[kolmn]]), function(y) {
          design$Columns[which((design$Iso == x)&(design[[kolmn]] == y))]
        }))
      })
      irsok1 <- FALSE
      if (!Param$Norma.Pep.Intens.IRS_Ref_channels %in% c(NA, "NA", FALSE, "FALSE")) {
        irstmp <- as.integer(unlist(strsplit(Param$Norma.Pep.Intens.IRS_Ref_channels, ";")))
        irstmp <- irstmp[sort(as.numeric(unique(Exp.map$Isobaric.set)))]
        if (length(irstmp) == length(Iso)) {
          irstmp <- data.frame(Channel = irstmp, Iso = as.integer(Iso))
          irstmp$Ref.Sample.Aggregate <- apply(irstmp[, c("Channel", "Iso")], 1, function(x) {
            #x <- irstmp[1, c("Channel", "Iso")]
            Exp.map$Ref.Sample.Aggregate[which((as.integer(Exp.map$Isobaric.label) == x[[1]])&(as.integer(Exp.map$Isobaric.set) == x[[2]]))]
          })
          irs1 <- as.data.frame(sapply(irstmp$Ref.Sample.Aggregate, function(x) {
            temp_norm3[[paste0(pep.ref[laststep], x)]]
          }))
          irsok1 <- TRUE
        }
        Exp.map$Use[match(irstmp$Ref.Sample.Aggregate, Exp.map$Ref.Sample.Aggregate)] <- FALSE
      }
      irsok2 <- FALSE
      ag <- parse.Param.aggreg(Param_filter(Param$Ratios.Groups.Ref.Aggregate.Level, "Rep"))
      ag <- ag$column
      test <- sapply(Iso, function(x) { list(Exp.map[which(Exp.map$Isobaric.set == x), ag]) })
      test2 <- sapply(unique(unlist(test)), function(x) { sum(sapply(test, function(y) { x %in% y })) })
      if (min(test2) >= length(Iso)) {
        if (max(test2) <= length(Iso)) {# Here in absence of an internal reference channel we are creating our average.
          if (!irsok1) {
            cat("Since no reference channels were provided, we will attempt to generate some dummy ones by averaging within TMT samples.\n")
          }
          irs2 <- as.data.frame(sapply(Columns, function(x) {
            apply(temp_norm3[,unlist(x)], 1, function(y) { log10(sum(10^is.all.good(y))) })
          }))
          irsok2 <- TRUE
        } else { stop("There seems to be an issue with the isobaric set column of the experiment map.") }
      } else { if (!irsok1) {
        warning("Not all conditions are contained in each Isobaric sample, IRS normalisation cannot be applied!")
      } }
      if (irsok1 + irsok2 == 2) {
        # Try to extend the validity of irs by replacing missing values with estimates from the average of values per Isobaric sample
        # w <- which(apply(irs1, 1, function(x) { length(is.all.good(x)) }) == 0 #( This would be the safest way, possibly)
        w <- which(apply(irs1, 1, function(x) { length(is.all.good(x)) }) < length(Iso))
        irs1[w,] <- irs2[w,]
      }
      if (irsok1) { irs <- irs1 } else { if (irsok2) { irs <- irs2 } }
      irsok <- (irsok1 + irsok2 > 0)
      if (irsok) {
        TxtSteps <- unique(c(TxtSteps, paste0("corrected against the ", IsobarLab, " batch effect using Internal-Reference Scaling")))
        dir <- paste0(nrmDr, "/IRS normalisation")
        if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
        dirlist <- unique(c(dirlist, dir))
        colnames(irs) <- paste0("sum", Iso)
        irs$Average <- apply(irs[, paste0("sum", Iso)], 1, function(x) { mean(is.all.good(x)) })
        # Scaling factor vectors
        for (i in Iso) { irs[[paste0("fact", i)]] <- irs$Average-irs[[paste0("sum", i)]] }
        # make new data frame with IRS normalized data
        data_irs <- temp_norm3[, unlist(Columns[[Iso[1]]])]+irs[[paste0("fact", Iso[1])]]
        for (i in Iso[2:length(Iso)]) {
          data_irs <- cbind(data_irs, temp_norm3[, unlist(Columns[[i]])]+irs[[paste0("fact", i)]])
        }
        design <- design[order(design$Original.order),]
        data_irs <- data_irs[, design$Columns]
        temp_norm3 <- temp_norm3[, design$Columns]
        kol <- colnames(data_irs)
        temp1 <- unlist(data_irs)
        temp2 <- unlist(temp_norm3)
        w <- which(!is.all.good(temp2, 2))
        temp1[w] <- temp2[w]
        data_irs <- as.data.frame(matrix(temp1, ncol = length(kol)))
        colnames(data_irs) <- gsub(topattern(pep.ref[laststep]), pep.ref[currstep], kol)
        # See what the IRS data look like:
        temp1 <- cbind(temp_norm2, data_irs)
        temp1$"Modified sequence" <- pep$"Modified sequence"[wAG]
        temp1 <- reshape2::melt(temp1, id.vars = "Modified sequence")
        temp1$variable <- as.character(temp1$variable)
        temp1 <- temp1[which(is.all.good(temp1$value, 2)),]
        temp1$Norm <- c("Original", "IRS")[grepl("^irs\\.", temp1$variable)+1]
        temp1$Norm <- factor(temp1$Norm, levels = c("Original", "IRS"))
        temp1$variable <- gsub_Rep(paste(c(topattern(pep.ref[laststep]), topattern(pep.ref[currstep])), collapse = "|"),
                                    "", temp1$variable)
        tmp1 <- unique(temp1$variable) 
        tmp2 <- Isapply(strsplit(tmp1, "___"), unlist)
        temp1[, RSA$names] <- tmp2[match(temp1$variable, tmp1),]
        temp1$Channel <- Exp.map$Isobaric.label[match(temp1$variable, Exp.map$Ref.Sample.Aggregate)]
        C <- cleanNms(gsub(topattern(pep.ref[laststep]), "", unlist(Columns)))
        temp1$variable <- factor(cleanNms(temp1$variable), levels = C)
        ttl1 <- "Peptides IRS intensity normalisation"
        ttl2 <- "Peptides IRS intensity normalisation - density distribution"
        plot1 <- ggplot(temp1) +
          geom_violin(aes(x = variable, y = value, color = Channel, fill = Channel), alpha = 0.25) +
          geom_boxplot(aes(x = variable, y = value, color = Channel, fill = Channel), alpha = 0.5) +
          theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          scale_y_continuous(expand = c(0, 0)) +
          scale_color_viridis_d(begin = 0.25) +
          scale_fill_viridis_d(begin = 0.25) +
          facet_grid(Norm~.) + ggtitle(ttl1)
        print(plot1) # This type of QC plot does not need to pop up, the side panel is fine
        ggsave(paste0(dir, "/", ttl1, ".jpeg"), plot1, dpi = 300, width = 10, height = 10, units = "in")
        ggsave(paste0(dir, "/", ttl1, ".pdf"), plot1, dpi = 300, width = 10, height = 10, units = "in")
        ReportCalls <- AddPlot2Report(Plot = plot1, Title = ttl1)
        plot2 <- ggplot(temp1) +
          geom_density(stat = "density", aes(x = value, color = Channel)) +
          scale_y_continuous(expand = c(0, 0)) +
          scale_color_viridis_d(begin = 0.25) +
          theme_bw() + facet_wrap(~Norm) + ggtitle(ttl2)
        print(plot2) # This type of QC plot does not need to pop up, the side panel is fine
        ggsave(paste0(dir, "/", ttl2, ".jpeg"), plot2, dpi = 300, width = 10, height = 10, units = "in")
        ggsave(paste0(dir, "/", ttl2, ".pdf"), plot2, dpi = 300, width = 10, height = 10, units = "in")
        ReportCalls <- AddPlot2Report(Plot = plot2, Title = ttl2)
        # PCA plot of IRS corrected data
        temp1 <- data_irs
        test <- apply(temp1, 1, function(x) { length(is.all.good(x)) })
        temp1 <- temp1[which(test == length(kol)),]
        colnames(temp1) <- gsub(topattern(pep.ref[currstep]), "", colnames(temp1))
        pc1 <- prcomp(t(temp1), scale. = TRUE)
        scores1 <- as.data.frame(pc1$x)
        if ("PC2" %in% colnames(scores1)) {
          scores1$RSA <- gsub(topattern(pep.ref[currstep]), "", rownames(scores1))
          scores1$Label <- rownames(scores1) <- cleanNms(rownames(scores1))
          scores1$Iso <- sapply(scores1$RSA, function(x) {
            Exp.map$Isobaric.set[match(x, Exp.map$Ref.Sample.Aggregate)]
          })
          scores1[, RSA$names] <- Isapply(strsplit(scores1$RSA, "___"), unlist)
          scores1$Iso <- factor(scores1$Iso, levels = as.numeric(Iso))
          pv1 <- round(100*(pc1$sdev)^2 / sum(pc1$sdev^2), 0)
          pv1 <- pv1[which(pv1 > 0)]
          pv1 <- paste0("Original: ", paste(sapply(1:length(pv1), function(x) {
            paste0("PC", x, ": ", pv1[x], "%")
          }), collapse = ", "))
          tst <- sapply(RSA$names, function(x) {length(unique(scores1[[x]]))})
          tst <- names(tst)[which(tst > 1)[1]]
          scores1[[tst]] <- as.factor(scores1[[tst]])
          nm1 <- "PCA plot - IRS normalisation result"
          plot <- ggplot(scores1) +
            geom_point(aes(x = PC1, y = PC2, colour = Iso, shape = .data[[tst]])) +
            scale_shape_manual(values = 1:nlevels(scores1[[tst]])) +
            scale_color_viridis_c(begin = 0.25) +
            coord_fixed() + theme_bw() +
            geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
            ggtitle("Data after IRS normalisation", subtitle = pv1) +
            geom_text_repel(aes(x = PC1, y = PC2, label = Label),
                            size = 2.5, show.legend = FALSE)
          ggsave(paste0(dir, "/", nm1, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
          ggsave(paste0(dir, "/", nm1, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
          ReportCalls <- AddPlot2Report(Title = nm1)
          Symb <- rep(c("circle", "diamond", "square", "cross", "x"), max(as.numeric(Rep)))[1:max(as.numeric(Rep))]             
          Symb <- Symb[as.numeric(scores1$Replicate)]
          #
          scores1$"Isobaric set" <- factor(scores1$Iso)
          if ("PC3" %in% colnames(scores1)) {
            plot_lyPCA1 <- plot_ly(scores1, x = ~PC1, y = ~PC2, z = ~PC3,
                                   text = ~Label, color = ~`Isobaric set`, colors = "viridis",
                                   type = "scatter3d", mode = "markers", symbol = I(Symb))
            # plot_lyPCA1 <- add_trace(plot_lyPCA1, scores1, x = ~PC1, y = ~PC2, z = ~PC3,
            #                          type = "scatter3d", mode = "text", showlegend = FALSE)
          } else {
            plot_lyPCA1 <- plot_ly(scores1, x = ~PC1, y = ~PC2,
                                   text = ~Label, color = ~`Isobaric set`, colors = "viridis",
                                   type = "scatter", mode = "markers", symbol = I(Symb))
            # plot_lyPCA1 <- add_trace(plot_lyPCA1, scores1, x = ~PC1, y = ~PC2,
            #                          type = "scatter", mode = "text", showlegend = FALSE)
          }
          plot_lyPCA1 %<o% layout(plot_lyPCA1, title = ttl1)
          saveWidget(plot_lyPCA1, paste0(dir, "/", nm1, ".html"), selfcontained = TRUE)
          #system(paste0("open \"", dir, "/", nm1, ".html"))
        }
        temp_norm2 <- temp_norm2*NA
        kol <- colnames(temp_norm2) <- colnames(data_irs)
        temp_norm2[wAG,] <- data_irs
        pepAllGood[[currstep]] <- apply(temp_norm2, 1, function(x) { length(is.all.good(x)) }) > 0
        tst1 <- sum(pepAllGood[[laststep]])
        tst2 <- sum(pepAllGood[[currstep]])
        if (tst2 < tst1) { warning("Quantitative peptide loss!") }
        pep_intens_norm[[currstep]] <- temp_norm2
        # View(pep_intens_norm[[currstep]])
        cat(" ->", currstep, "done!\n")
      } else {
        pep.ref <- pep.ref[which(names(pep.ref) != "IRS normalisation")]
        currstep <- laststep
      }
    }
  }
  # Step 4: Use the Levenberg-Marquardt method to align the different datasets:
  # Should be MQ.Exp-specific!
  if (advnorm) {
    TxtSteps <- unique(c(TxtSteps, "re-normalized using the Levenberg-Marquardt procedure"))
    currstep <- "Advanced (Levenberg-Marquardt) normalisation"
    laststep <- names(pep.ref)[match(currstep, names(pep.ref))-1]
    wAG <- which(pepAllGood[[laststep]])
    temp_norm2 <- pep_intens_norm[[laststep]]
    temp_norm3 <- temp_norm2[wAG,]
    kol <- colnames(temp_norm3)
    RefGrp <- RSA$value
    RefGrp <- RefGrp[which(paste0(pep.ref[laststep], RefGrp) %in% kol)]
    RefGrp <- setNames(sapply(RefGrp, function(x) {
      Exp.map[match(x, Exp.map$Ref.Sample.Aggregate), RRG$column]
    }), RefGrp)
    tst <- aggregate(RefGrp, list(RefGrp), length)
    tst <- tst$x[match(RefGrp, tst$Group.1)]
    w1 <- which(tst > 1)
    w2 <- which(tst <= 1)
    if (length(w2)) {
      warning("Removing some invalid (normalisation?) columns!")
      RefGrp <- RefGrp[w1]
      kol <- kol[w1]
      temp_norm3 <- temp_norm3[, kol]
    }
    temp_norm4 <- temp_norm3*NA
    colnames(temp_norm4) <- gsub(topattern(pep.ref[laststep]), pep.ref[currstep], colnames(temp_norm4))
    temp_norm3$"Modified sequence" <- pep$"Modified sequence"[wAG]
    NormGrps2 <- pep$`Normalisation group`[wAG]
    #temp_norm3$"Modified sequence" <- NULL
    exports <- list("MQ.Exp", "temp_norm3", "pep.ref", "laststep", "currstep", "Param", "RSA",
                    "RefGrp", "NormGrps2")
    source(parSrc, local = FALSE)
    clusterExport(parClust, exports, envir = environment())
    clusterCall(parClust, function() library(proteoCraft))
    #if (Param$Adv.Norma.Pep.Intens.Type == "C") { # "C" here means by columns
      refgrp <- unique(RefGrp)
      norm_temp <- parSapply(parClust, refgrp, function(i) { #i <- refgrp[1]
        clmnsa <- paste0(pep.ref[laststep], names(RefGrp)[which(RefGrp == i)])
        temp_norm3b <- temp_norm3a <- temp_norm3[, c("Modified sequence", clmnsa)]
        temp_norm3b[,clmnsa] <- temp_norm3b[,clmnsa]*NA
        for (nrmgrp in unique(NormGrps2)) { #nrmgrp <- unique(NormGrps2)[1]
          wNG <- which(NormGrps2 == nrmgrp)
          tmp_nrm <- AdvNorm.IL(temp_norm3a[wNG,], "Modified sequence", clmnsa, TRUE, 5)
          temp_norm3b[wNG, clmnsa] <- tmp_nrm[, paste0("AdvNorm.", clmnsa)]
        }
        colnames(temp_norm3b) <- gsub(topattern(pep.ref[laststep]), pep.ref[currstep], colnames(temp_norm3b))
        clmnsb <- grep(topattern(pep.ref[currstep]), colnames(temp_norm3b), value = TRUE)
        #cat("Advanced peptides intensity normalisation done for Search-Engine-Experiment", i, "\n")
        return(list(temp_norm3b[,clmnsb]))
      }, USE.NAMES = TRUE)
    # } else {
    #   stop("Not implemented yet, I need to update the current AdvNorm function as it takes too long or crashes.")
    # }
    for (i in refgrp) { #i <- refgrp[1]
      tmp <- norm_temp[[i]]
      temp_norm4[, colnames(tmp)] <- tmp
    }
    colnames(temp_norm2) <- gsub(topattern(pep.ref[laststep]), pep.ref[currstep], colnames(temp_norm2))
    temp_norm2 <- temp_norm2[, which(colnames(temp_norm2) %in% colnames(temp_norm4))]
    temp_norm2[wAG,] <- temp_norm4
    kol <- colnames(temp_norm2)
    pepAllGood[[currstep]] <- apply(temp_norm2, 1, function(x) { length(is.all.good(x)) }) > 0
    tst1 <- sum(pepAllGood[[laststep]])
    tst2 <- sum(pepAllGood[[currstep]])
    if (tst2 < tst1) { warning("Quantitative peptide loss!") }
    pep_intens_norm[[currstep]] <- temp_norm2
    # View(pep_intens_norm[[currstep]])
    cat(" ->", currstep, "done!\n")
  }
  # Step 5: Batch correction:
  if (combatnorm) {
    test <- sum(!unlist(strsplit(as.character(Param$Batch.correction), ";")) %in% names(Aggregates))
    if (sum(test)) {
      warning("Skipping Batch correction as the supplied aggregate is invalid!")
      combatnorm <- FALSE
    } else {
      if (length(Batch.correction$values) == 1) {
        warning("Skipping Batch correction as there is only one batch!")
        combatnorm <- FALSE
      } else {
        bioc_req <- unique(c(bioc_req, "sva"))
        biocInstall("sva")
        tmp <- tolower(paste(Batch.correction$names, collapse = "/"))
        TxtSteps <- unique(c(TxtSteps, paste0("corrected", c("", " once again")[(sum(grepl("Internal-Reference Scaling", TxtSteps))>0)+1],
                                              " using the ComBat function from the sva package to remove the ",
                                              c(paste0(tmp, "-related "), "")[(tmp == "batch")+1], "batch effect")))
        currstep <- "SVA batch correction"
        laststep <- names(pep.ref)[match(currstep, names(pep.ref))-1]
        wAG <- which(pepAllGood[[laststep]])
        dir <- paste0(nrmDr, "/Batch correction")
        if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
        dirlist <- unique(c(dirlist, dir))
        batch_aggr <- Batch.correction$column
        # First let's look at the pre-batch correction data with a PCA plot:
        temp2 <- pep_intens_norm[[laststep]][wAG,]
        kol <- colnames(temp2)
        test <- apply(temp2, 1, function(x) {length(is.all.good(x))})
        temp2 <- temp2[which(test == length(kol)),]
        colnames(temp2) <- gsub(topattern(pep.ref[laststep]), "", colnames(temp2))
        pc2 <- prcomp(t(temp2), scale. = TRUE)
        scores2 <- as.data.frame(pc2$x)
        colnames(scores2) <- paste0("Original_", colnames(scores2))
        scores2$Ref.Sample.Aggregate <- rownames(scores2) <- gsub_Rep(topattern(pep.ref[laststep]), "",
                                                                       rownames(scores2))
        scores2[, RSA$names] <- Isapply(strsplit(scores2$Ref.Sample.Aggregate, "___"), unlist)
        scores2$Label <- cleanNms(rownames(scores2))
        scores2$Batch <- as.factor(Exp.map[match(scores2$Ref.Sample.Aggregate, Exp.map$Ref.Sample.Aggregate),
                                           batch_aggr])
        pv2 <- round(100*(pc2$sdev)^2 / sum(pc2$sdev^2), 0)
        pv2 <- pv2[which(pv2 > 0)]
        pv2 <- paste0("Original: ", paste(sapply(1:length(pv2), function(x) {
          paste0("PC", x, ": ", pv2[x], "%")
        }), collapse = ", "))
        print(pv2)
        scores2$Corrected <- "Original"
        tst <- sapply(RSA$names, function(x) { length(unique(scores2[[x]])) })
        tst <- names(tst)[which(tst > 1)[1]]
        if ("Original_PC2" %in% colnames(scores2)) {
          nm2 <- "PCA plot - SVA batch correction_before"
          plot <- ggplot(scores2) +
            geom_point(aes(x = Original_PC1, y = Original_PC2, color = .data[[tst]], shape = Batch)) +
            scale_shape_manual(values = (14+1:nlevels(scores2$Batch))%%25) +
            scale_color_viridis(begin = 0.25, discrete = TRUE, option = "D") +
            coord_fixed() + theme_bw() +
            geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
            ggtitle(nm2, subtitle = pv2) +
            geom_text_repel(aes(x = Original_PC1, y = Original_PC2, label = Label),
                            size = 2.5, show.legend = FALSE)
          #poplot(plot)
          ggsave(paste0(dir, "/", nm2, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
          ggsave(paste0(dir, "/", nm2, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
          ReportCalls <- AddPlot2Report(Title = nm2)
          Symb <- rep(c("circle", "diamond", "square", "cross", "x"), max(as.numeric(Rep)))[1:max(as.numeric(Rep))]
          Symb <- Symb[as.numeric(scores2$Replicate)]
          #
          scores2$"Samples group" <- factor(scores2$Batch)
          if ("Original_PC3" %in% colnames(scores2)) {
            plot_lyPCA2 <- plot_ly(scores2, x = ~Original_PC1, y = ~Original_PC2, z = ~Original_PC3,
                                   color = ~`Samples group`, colors = "viridis",
                                   text = ~Label, type = "scatter3d", mode = "markers", symbol = I(Symb))
            plot_lyPCA2 <- add_trace(plot_lyPCA2, scores2, x = ~Original_PC1, y = ~Original_PC2, z = ~Original_PC3,
                                     type = "scatter3d", mode = "text", showlegend = FALSE)
          } else {
            plot_lyPCA2 <- plot_ly(scores2, x = ~Original_PC1, y = ~Original_PC2,
                                   color = ~`Samples group`, colors = "viridis",
                                   text = ~Label, type = "scatter", mode = "markers", symbol = I(Symb))
            plot_lyPCA2 <- add_trace(plot_lyPCA2, scores2, x = ~Original_PC1, y = ~Original_PC2,
                                     type = "scatter", mode = "text", showlegend = FALSE)
          }
          plot_lyPCA2 %<o% layout(plot_lyPCA2, title = nm2)
          saveWidget(plot_lyPCA2, paste0(dir, "/", nm2, ".html"))
          #system(paste0("open \"", dir, "/", nm2, ".html"))
        } else { warning("PCA failed, investigate!") }
        # Now let's do some imputation:
        # (we will have to remember to remove those afterwards!)
        #grps <- Exp.map$Experiment[match(gsub(topattern(pep.ref[laststep]), "", kol), Exp.map$Ref.Sample.Aggregate)]
        grps <- Exp.map[match(gsub(topattern(pep.ref[laststep]), "", kol), Exp.map$Ref.Sample.Aggregate), VPAL$column]
        # Treat Mixed as one group
        grps <- gsub(".+___mixed___.+|^mixed___.+|.+___mixed$", "Mixed", grps, ignore.case = TRUE)
        temp <- pep_intens_norm[[laststep]][wAG,]
        test <- Data_Impute2(temp, grps)
        temp2 <- test$Imputed_data
        Pos <- test$Positions_Imputed
        #test <- t(apply(temp2, 1, function(x) { aggregate(x, list(grps), function(y) { length(is.all.good(y)) })[[2]] }))
        #test2 <- apply(test, 1, function(x) { length(which(x == 0)) })
        #WhTest <- which(test2 == 0)
        WhTest <- which(apply(temp, 1, function(x) { length(is.all.good(x)) }) > 0)
        # Now for the actual batch correction:
        col2 <- gsub(topattern(pep.ref[laststep]), "", kol)
        e <- Exp.map[match(col2, Exp.map$Ref.Sample.Aggregate), Aggregates]
        batches <- "Batches"
        e[[batches]] <- e[[batch_aggr]]
        modcombat <- model.matrix(~1, data = e)
        tstuniq <- apply(temp2[WhTest,], 1, function(x) { length(unique(x)) }) > 1
        Whuniq <- WhTest[which(tstuniq)]
        tempmatr <- as.matrix(temp2[Whuniq,])
        # (I checked and it's considered best to apply ComBat to log-transformed data.)
        NormGrps2 <- pep$`Normalisation group`[wAG[Whuniq]]
        sva.exprs_tmp <- setNames(lapply(unique(NormGrps2), function(nrmgrp) { #nrmgrp <- unique(NormGrps2)[1]
          w <- which(NormGrps2 == nrmgrp)
          ComBat(dat = tempmatr[w,], batch = e[[batches]], mod = modcombat)
        }), unique(NormGrps2))
        sva.exprs <- pep_intens_norm[[laststep]]*NA
        for (nrmgrp in unique(NormGrps2)) {
          w1 <- wAG[Whuniq]
          w2 <- which(NormGrps2 == nrmgrp)
          sva.exprs[w1[w2],] <- sva.exprs_tmp[[nrmgrp]]
        }
        # Remove imputed data:
        w <- which(Pos, arr.ind = TRUE)
        if (nrow(w)) { sva.exprs[wAG,][w] <- temp[w] }
        col3 <- paste0(pep.ref[currstep], col2)
        colnames(sva.exprs) <- col3 
        # PCA plot of results
        temp3 <- sva.exprs
        test <- apply(temp3[,col3], 1, function(x) { length(is.all.good(x)) })
        temp3 <- temp3[which(test == length(col3)),]
        colnames(temp3) <- gsub(topattern(pep.ref[currstep]), "", colnames(temp3))
        pc3 <- prcomp(t(temp3), scale. = TRUE)
        scores3 <- as.data.frame(pc3$x)
        colnames(scores3) <- paste0("Corrected_", colnames(scores3))
        scores3$Ref.Sample.Aggregate <- rownames(scores3)
        scores3[, RSA$names] <- Isapply(strsplit(rownames(scores3), "___"), unlist)
        scores3$Label <- cleanNms(rownames(scores3))
        scores3$Batch <- as.factor(Exp.map[match(scores3$Ref.Sample.Aggregate, Exp.map$Ref.Sample.Aggregate),
                                           batch_aggr])
        pv3 <- round(100*(pc3$sdev)^2 / sum(pc3$sdev^2), 0)
        pv3 <- pv3[which(pv3 > 0)]
        pv3 <- paste0("Corrected: ", paste(sapply(1:length(pv3), function(x) {
          paste0("PC", x, ": ", pv3[x], "%")
        }), collapse = ", "))
        print(pv3)
        scores3$Corrected <- "Corrected"
        temp2 <- set_colnames(scores2, gsub("^Original_", "", colnames(scores2)))
        temp3 <- set_colnames(scores3, gsub("^Corrected_", "", colnames(scores3)))
        colu <- unique(c(colnames(temp2), colnames(temp3)))
        colu <- colu[which((colu %in% colnames(temp2))&(colu %in% colnames(temp3)))]
        scores <- rbind(temp2[,colu], temp3[,colu])
        scores$Corrected <- factor(scores$Corrected, levels = c("Original", "Corrected"))
        tst <- sapply(RSA$names, function(x) { length(unique(scores2[[x]])) })
        tst <- names(tst)[which(tst > 1)[1]]
        if ("PC2" %in% colnames(scores)) {
          nm3 <- "PCA plot - SVA batch correction_after"
          nm0 <- "PCA plot - SVA batch correction"
          plot <- ggplot(scores3) +
            geom_point(aes(x = Corrected_PC1, y = Corrected_PC2, color = .data[[tst]], shape = Batch)) +
            scale_shape_manual(values = (14+1:nlevels(scores3$Batch))%%25) +
            scale_color_viridis(begin = 0.25, discrete = TRUE, option = "D") +
            coord_fixed() + theme_bw() +
            geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
            ggtitle(nm3, subtitle = pv3) +
            geom_text_repel(aes(x = Corrected_PC1, y = Corrected_PC2, label = Label),
                            size = 2.5, show.legend = FALSE)
          #poplot(plot)
          ggsave(paste0(dir, "/", nm3, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
          ggsave(paste0(dir, "/", nm3, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
          ReportCalls <- AddPlot2Report(Title = nm3)
          Symb <- rep(c("circle", "diamond", "square", "cross", "x"), max(as.numeric(Rep)))[1:max(as.numeric(Rep))]             
          Symb <- Symb[as.numeric(scores3$Replicate)]
          #
          scores3$"Samples group" <- factor(scores3$Batch)
          if ("PC3" %in% colnames(scores)) {
            plot_lyPCA3 <- plot_ly(scores3, x = ~Corrected_PC1, y = ~Corrected_PC2, z = ~Corrected_PC3,
                                   color = ~`Samples group`, colors = "viridis",
                                   text = ~Label, type = "scatter3d", mode = "markers", symbol = I(Symb))
            plot_lyPCA3 <- add_trace(plot_lyPCA3, scores3, x = ~Corrected_PC1, y = ~Corrected_PC2, z = ~Corrected_PC3,
                                     type = "scatter3d", mode = "text", showlegend = FALSE)
          } else {
            plot_lyPCA3 <- plot_ly(scores3, x = ~Corrected_PC1, y = ~Corrected_PC2,
                                   color = ~`Samples group`, colors = "viridis",
                                   text = ~Label, type = "scatter", mode = "markers", symbol = I(Symb))
            plot_lyPCA3 <- add_trace(plot_lyPCA3, scores3, x = ~Corrected_PC1, y = ~Corrected_PC2,
                                     type = "scatter", mode = "text", showlegend = FALSE)
          }
          plot_lyPCA3 %<o% layout(plot_lyPCA3, title = nm3)
          saveWidget(plot_lyPCA3, paste0(dir, "/", nm3, ".html"))
          #system(paste0("open \"", dir, "/", nm3, ".html"))
          # NB: There is currently no way to create a 3D, faceted plot in plotly for R that I know of)
          # Ok, there is, embedding several plots in a shiny app!
          plot <- ggplot(scores) +
            geom_point(aes(x = PC1, y = PC2, color = .data[[tst]], shape = Batch)) +
            scale_shape_manual(values = (14+1:nlevels(scores$Batch))%%25) +
            scale_color_viridis(begin = 0.25, discrete = TRUE, option = "D") +
            coord_fixed() + theme_bw() +
            geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
            ggtitle(nm0,
                    subtitle = paste(pv2, pv3, sep = "\n")) +
            geom_text_repel(aes(x = PC1, y = PC2, label = Label),
                            size = 2.5, show.legend = FALSE) + 
            facet_wrap(~ Corrected)
          ggsave(paste0(dir, "/", nm0, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
          ggsave(paste0(dir, "/", nm0, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
          ReportCalls <- AddPlot2Report(Title = nm0)
        } else { warning("PCA failed, investigate!") }
        KeepComBatRes %<o% TRUE
        PCs <- data.frame("Component" = paste0("PC", as.character(1:length(pc2$sdev))),
                          "Before (%)" = round(100*(pc2$sdev)^2 / sum(pc2$sdev^2), 0),
                          "After (%)" = round(100*(pc3$sdev)^2 / sum(pc2$sdev^2), 0))
        msg <- "Keep results from ComBat batch correction? (untick to cancel correction)"
        appNm <- paste0(dtstNm, " - ComBat batch correction")
        ui <- fluidPage(
          useShinyjs(),
          extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
          tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
          titlePanel(tag("u", "ComBat batch correction"),
                     appNm),
          br(),
          fluidRow(column(5,
                          checkboxInput("KeepResults", msg, KeepComBatRes),
                          actionButton("saveBtn", "Save"),
                          h4("Recommended criteria:"),
                          h5(HTML("&nbsp;Does the original grouping follow known batches?")),
                          h5(HTML("&nbsp;&nbsp;-> If no: only accept the correction if it improves the apparent grouping of expected sample groups.")),
                          h5(HTML("&nbsp;&nbsp;-> If yes: accept the correction if...")),
                          h5(HTML("&nbsp;&nbsp;&nbsp;- ... it removes the original grouping by batches...")),
                          h5(HTML("&nbsp;&nbsp;&nbsp;- ... or it improves the apparent grouping of expected sample groups.")),
                          withSpinner(DTOutput("PCs")),
                          br(),
                          br(),
                          br()),
                   column(7,
                          withSpinner(plotlyOutput("Before", height = "550px")),
                          withSpinner(plotlyOutput("After", height = "550px")))),
          br(),
          br()
        )
        server <- function(input, output, session) {
          output$Before <- renderPlotly(plot_lyPCA2)
          output$After <- renderPlotly(plot_lyPCA3)
          output$PCs <- renderDT({ PCs },
                                 FALSE,
                                 escape = FALSE,
                                 selection = "none",
                                 editable = FALSE,
                                 rownames = FALSE,
                                 options = list(
                                   dom = 't',
                                   paging = FALSE,
                                   ordering = FALSE
                                 ),
                                 callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
          observeEvent(input[["KeepResults"]], {
            assign("KeepComBatRes", as.logical(input[["KeepResults"]]), envir = .GlobalEnv)
          })
          observeEvent(input$saveBtn, { stopApp() })
          #observeEvent(input$cancel, { stopApp() })
          session$onSessionEnded(function() { stopApp() })
        }
        eval(parse(text = runApp), envir = .GlobalEnv)
        #
        tmp <- AllAnsw[1,]
        tmp[, c("Parameter", "Message")] <- c("KeepComBatRes", msg)
        tmp$Value <- list(KeepComBatRes)
        m <- match("KeepComBatRes", AllAnsw$Parameter)
        if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
        if (KeepComBatRes) {
          write(c("Decision:", "---------", "", "The ComBat batch correction results were accepted."),
                file = paste0(wd, "/Workflow control/Peptides/Intensities/Batch correction/Decision.txt"))
          pepAllGood[[currstep]] <- apply(sva.exprs, 1, function(x) { length(is.all.good(x)) }) > 0
          tst1 <- sum(pepAllGood[[laststep]])
          tst2 <- sum(pepAllGood[[currstep]])
          if (tst2 < tst1) { warning("Quantitative peptide loss!") }
          pep_intens_norm[[currstep]] <- sva.exprs
          cat(" ->", currstep, "done!\n")
        } else {
          write(c("Decision:", "---------", "", "The ComBat batch correction results were rejected."),
                file = "Workflow control/Peptides/Intensities/Batch correction/Decision.txt")
          pep.ref <- pep.ref[1:(which(pep.ref == paste0("sva. ", pep.ref["Original"]))-1)]
          TxtSteps <- grep("ComBat", TxtSteps, invert = TRUE)
          cat(" ->", currstep, "normalisation rejected by user, reverting...\n")
          currstep <- laststep
        }
        # View(pep_intens_norm[[currstep]])
      }
    }
  }
  # De-log:
  for (i in names(pep_intens_norm)) {
    for (j in colnames(pep_intens_norm[[i]])) { pep_intens_norm[[i]][[j]] <- 10^pep_intens_norm[[i]][[j]] }
  }
  # Sometimes, we skip a step because some preliminary test concluded it did not make sense
  # We should then remove any pep.ref created which we do not need anymore
  pep.ref <- pep.ref[which(names(pep.ref) %in% c(names(pep.ref)[1:(match("Classic normalisation", names(pep.ref))-1)],
                                                 names(pep_intens_norm)))]
  # Visualize:
  if (Param$Norma.Pep.Intens.show) {
    a1 <- paste0(pep.ref["Original"], RSA$values)
    a1 <- a1[which(a1 %in% colnames(pep))]
    kol <- c("Modified sequence", "Normalisation group")
    temp <- pep[, c(kol, a1)]
    temp <- reshape2::melt(temp, id.vars = kol)
    temp$Norm <- "Original"
    temp$Sample <- gsub_Rep(topattern(pep.ref["Original"]), "", temp$variable)
    for (i in names(pep_intens_norm)[2:length(names(pep_intens_norm))]) { #i <- names(pep_intens_norm)[2]
      temp2 <- pep_intens_norm[[i]]
      temp2[, kol] <- pep[, kol]
      temp2 <- reshape2::melt(temp2, id.vars = kol)
      temp2$Norm <- i
      temp2$Sample <- gsub_Rep(topattern(pep.ref[i]), "", temp2$variable)
      temp <- rbind(temp, temp2)
    }
    colnames(temp)[which(colnames(temp) == "Normalisation group")] <- "Normalisation_group"
    temp$Norm <- factor(temp$Norm, levels = unique(temp$Norm))
    temp$value <- log10(temp$value)
    temp[, RSA$names] <- Exp.map[match(temp$Sample, Exp.map$Ref.Sample.Aggregate), RSA$names]
    temp$Sample <- cleanNms(temp$Sample)
    lev <- cleanNms(RSA$values)
    temp$Sample <- factor(temp$Sample, levels = lev)
    temp <- temp[which(is.all.good(temp$value, 2)),]
    w <- which(sapply(VPAL$names, function(x) { length(unique(temp[[x]])) }) > 1)
    w <- w[which(tolower(substr(names(w), 1, 3)) != "rep")]
    temp$Samples_group <- do.call(paste, c(temp[, VPAL$names[w], drop = FALSE], sep = " "))
    temp$Replicate <- as.factor(temp$Replicate)
    ttl <- "Peptides intensity normalisation"
    dir <- paste0(nrmDr, "/Normalisation - summary")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    kolz <- "."
    if (length(unique(temp$Normalisation_group)) > 1) { kolz <- c(kolz, "Normalisation_group") } else {
      if (length(tst) > 1) { kolz <- c(kolz, tst[2:length(tst)]) }
    }
    if (length(kolz) > 2) { kolz <- kolz[3:length(kolz)] } else { kolz <- kolz[1]}
    form <- as.formula(paste0("Norm~", paste0(kolz, collapse = "+")))
    plot <- ggplot(temp) +
      geom_violin(aes(x = Sample, y = value, color = Samples_group, fill = Samples_group), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, color = Samples_group, fill = Samples_group), alpha = 0.5) +
      scale_color_viridis(begin = 0.25, discrete = TRUE, option = "C") +
      scale_fill_viridis(begin = 0.25, discrete = TRUE, option = "C") +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      geom_vline(xintercept = 0) + ggtitle(ttl) + facet_grid(form) +
      theme(strip.text.y = element_text(angle = 0))
    print(plot) # This type of QC plot does not need to pop up, the side panel is fine
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
  }
  #
}
# Here we want to make sure that we only have as valid experiment factors (aggregates) values the final ones (e.g. removing ones needed only for normalisation, such as a mixed channel for IRS)
if ((Param$Norma.Pep.Intens)&&(irsok)&&(!Param$Norma.Pep.Intens.IRS_Ref_channels %in% c("", NA, "NA", FALSE, "FALSE"))) {
  warning("Removing mixed channels used for IRS normalisation")
}
Exp.map <- Exp.map[which(as.character(Exp.map$Use) %in% c("T", "TRUE")),]
if (Param$Norma.Pep.Intens) {
  for (i in names(pep.ref)) { #i <- "Original"
    tmp <- pep_intens_norm[[i]]
    tst <- gsub(topattern(pep.ref[[i]]), "", colnames(tmp)) %in% Exp.map$Ref.Sample.Aggregate
    pep_intens_norm[[i]] <- tmp[, which(tst)]
  }
  temp <- pep_intens_norm[[names(pep.ref)[length(pep.ref)]]]
  pep[, colnames(temp)] <- temp
  saveFun(pep_intens_norm, paste0(nrmDr, "/pep_intens_norm.RData"))
}
for (i in Param.aggreg) { #i <- Param.aggreg[1]
  tmp <- get(i)
  tmp$values <- unique(Exp.map[[tmp$column]])
  assign(i, tmp)
}
l <- sapply(Aggregate.map$Characteristics, length)
i1 <- which(l == 1)
i2 <- which(l > 1)
for (i in i1) { #i <- i1[1]
  Aggregate.list[[i]] <- sort(unique(Exp.map[[Aggregate.map$Characteristics[[i]]]]))
}
for (i in i2) { #i <- i2[1]
  Aggregate.list[[i]] <- sort(unique(Exp.map[[Aggregate.map$Aggregate.Name[[i]]]]))
}
for (i in names(Aggregate.list)) {
  assign(i, Aggregate.list[[i]])
}
# Also update Factors
FactorsLevels <- setNames(lapply(Factors, function(Fact) {
  unique(Exp.map[[Fact]])
}), Factors)
w <- which(sapply(FactorsLevels, length) > 0)
Factors <- Factors[w]
FactorsLevels <- FactorsLevels[Factors]
# Update abbreviations (RSA, VPAL, RRG and RG) in case we removed Mixed channels
RSA <- Ref.Sample.Aggregate
VPAL <- Volcano.plots.Aggregate.Level
RRG <- Ratios.Ref.Groups
RG <- Ratios.Groups
#
l <- length(TxtSteps)
if (l > 1) { TxtSteps <- paste0(paste(TxtSteps[1:(l-1)], collapse = ", "), " and ", TxtSteps[l]) }
DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", TxtSteps, ". Peptidoform-level ratios were then calculated.")

# Check variance/intensity dependency after normalisation
rfnm <- names(pep.ref)[length(pep.ref)]
ttl <- "Peptides Variance vs Intensity dependency"
kol <- grep(topattern(pep.ref[rfnm]), colnames(pep), value = TRUE)
temp <- pep[, kol]
w <- which(apply(temp, 1, function(x) { length(is.all.good(x)) }) > 0)
temp <- temp[w,]
Aggr <- VPAL
if (LocAnalysis) { Aggr <- parse.Param.aggreg("Exp;Com") }
tst <- sapply(Aggr$values, function(x) {
  em <- Exp.map[which(Exp.map[[Aggr$column]] == x),]
  kl <- paste0(pep.ref[rfnm], em$Ref.Sample.Aggregate)
  return(rowMeans(temp[, kl, drop = FALSE], na.rm = TRUE))
})
kol <- colnames(tst) <- cleanNms(colnames(tst))
tst <- apply(tst, 1, function(x) { kol[which(x == max(x))][1] })
tst2 <- aggregate(ev$Charge, list(ev$`Modified sequence`), function(x) { round(mean(x)) })
temp2 <- data.frame("log10(Mean)" = log10(rowMeans(temp, na.rm = TRUE)),
                    "log10(Variance)" = log10(apply(temp, 1, var, na.rm = TRUE)),
                    "Strongest in..." = tst,
                    check.names = FALSE)
w <- which(!is.na(temp2$`Strongest in...`))
temp2 <- temp2[w,]
temp2$"Main charge" <- paste0("Z = ", tst2$x[match(pep$`Modified sequence`[w], tst2$Group.1)])
temp2$"Main charge" <- factor(temp2$"Main charge", levels = paste0("Z = ", as.character(1:8)))
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
plot <- ggplot(temp2) +
  geom_scattermore(aes(x = `log10(Mean)`, y = `log10(Variance)`, colour = `Strongest in...`),
                   size = 1, alpha = 1) + ggtitle(ttl, subtitle = rfnm) + coord_fixed(0.3) + theme_bw() +
  scale_color_viridis_d(begin = 0.25) +
  facet_grid(`Strongest in...` ~ `Main charge`) + theme(strip.text.y.right = element_text(angle = 0))
poplot(plot, 12, 22)
ggsave(paste0(dir, "/", ttl, " - ", rfnm, ".jpeg"), plot, width = 10, height = 10, units = "in", dpi = 300)
ggsave(paste0(dir, "/", ttl, " - ", rfnm, ".pdf"), plot, width = 10, height = 10, units = "in", dpi = 300)
ReportCalls <- AddPlot2Report(Title = paste0(ttl, " - ", rfnm))

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
  temp <- sapply(RRG$names, function(x) { list(which(Exp.map[[x]] == j[[x]])) })
  temp2 <- sort(unique(unlist(temp)))
  test <- sapply(temp2, function(x) { sum(sapply(temp, function(y) {x %in% unlist(y)})) })
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
  temp <- sapply(RRG$names, function(x) { list(which(Exp.map[[x]] == j[[x]])) })
  temp2 <- sort(unique(unlist(temp)))
  test <- sapply(temp2, function(x) { sum(sapply(temp, function(y) {x %in% unlist(y)})) })
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
tst <- sapply(aggr, function(a) { length(get(substr(a, 1, 3))) })
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
ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
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
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 4*(length(unique(temp$"Ratios group"))+1),
         height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 4*(length(unique(temp$"Ratios group"))+1),
         height = 10, units = "in")
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
      test <- sapply(Adv.Norma.Pep.Ratio.Type.Group$values, function(i) { #i <- Adv.Norma.Pep.Ratio.Type.Group$values[1]
        i <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[k]] == i)]
        return(length(which(paste0(pep.ratios.ref[length(pep.ratios.ref)], i) %in% colnames(pep))))
      })
      agg <- Adv.Norma.Pep.Ratio.Type.Group$values[which(test > 1)]
      exports <- list("agg", "Adv.Norma.Pep.Ratio.Type.Group", "Exp.map", "pep.ratios.ref", "pep", "Param")
      clusterExport(parClust, exports, envir = environment())
      clusterCall(parClust, function() library(proteoCraft))
      norm_temp <- parSapply(parClust, 0:length(agg), function(i) { #i <- 1
        if (i == 0) {
          kol <- grep(paste0(topattern(pep.ratios.ref[1]), ".+_REF\\.to\\.REF_"), colnames(pep), value = TRUE)
        } else {
          x <- agg[i]
          j <- unlist(strsplit(x, "___"))
          names(j) <- Adv.Norma.Pep.Ratio.Type.Group$names
          temp <- sapply(Adv.Norma.Pep.Ratio.Type.Group$names, function(x) {
            list(which(Exp.map[[x]] == j[[x]]))
          })
          temp2 <- sort(unique(unlist(temp)))
          test <- sapply(temp2, function(x) { sum(sapply(temp, function(y) {x %in% unlist(y)})) })
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
          return(list(temp2))
          temp2$"Modified sequence" <- NULL
        }
      })
    } else {
      stop("Not implemented yet, I need to update the current AdvNorm function as it takes too long or crashes.")
    }
    for (i in 1:length(norm_temp)) {
      tmp <- norm_temp[[i]]
      pep[, colnames(tmp)] <- tmp
    }
    pep.ratios.ref <- unique(c(pep.ratios.ref, paste0("AdvNorm. ", pep.ratios.ref[1])))
  }
  if (Param$Norma.Pep.Ratio.show) {
    for (i in Ratios.Plot.split$values) { #i <- Ratios.Plot.split$values[1]
      j <- unlist(strsplit(i, "___"))
      names(j) <- unlist(Aggregate.map$Characteristics[which(Aggregate.map$Aggregate.Name == Ratios.Plot.split$aggregate)])
      #k <- sapply(c(1:length(j)), function(x) {list(which((Exp.map[[names(j)[x]]] == j[x])&(!Exp.map$Reference)))})
      k <- sapply(c(1:length(j)), function(x) {list(which(Exp.map[[names(j)[x]]] == j[x]))})
      l <- sort(unique(unlist(k)))
      test <- sapply(l, function(x) { sum(sapply(k, function(y) {x %in% y})) == length(k) })
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
        temp$X <- factor(temp$X, levels = unlist(sapply(c("Original", "Normalised"), function(x) {paste(x, unique(temp2[[Ratios.Plot.colour$names]]), sep = "_")})))
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
        ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
        ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
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
      ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ReportCalls <- AddPlot2Report(Title = ttl)
      DatAnalysisTxt <- paste0(DatAnalysisTxt, " Peptide ratios were then re-normalized.")
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
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\"Starting protein groups assembly:\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
if ("Prot.Only.with.at.least" %in% colnames(Param)) {
  NPep <- as.integer(Param$Prot.Only.with.at.least)
  if ((is.na(NPep))||(NPep <= 0)) {
    warning("Invalid `\"Prot.Only.with.at.least\" parameter, defaulting to 1!")
    NPep <- 1
  }
} else { NPep <- 1 }
.obj <- unique(c(.obj, "NPep"))
tm1 <- Sys.time()
source(parSrc, local = FALSE)
PG_assembly <- PG_assemble(pep, "id", "Proteins", db, "Evidence IDs", ev, Custom_PGs = custPGs, Npep = NPep,
                           cl = parClust)
saveFun(PG_assembly, file = "PG_assembly.RData")
#loadFun("PG_assembly.RData")
tm2 <- Sys.time()
PG %<o% PG_assembly$Protein.groups
pep <- PG_assembly$Peptides
db <- PG_assembly$Database
if ("Evidences" %in% names(PG_assembly)) { ev <- PG_assembly$Evidences }
msg <- paste0(nrow(PG), " protein groups assembled in ", gsub("^Time difference of ", "", capture.output(tm2-tm1)))
ReportCalls <- AddMsg2Report(Space = FALSE, Print = FALSE)

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
if (tstorg) {
  test <- sapply(strsplit(PG$`Protein IDs`, ";"), function(x) {
    paste(sort(unique(c(db[match(x, db$`Protein ID`), dbOrgKol]))), collapse = ";")
  })
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
ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300)
ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)

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
      return(paste(x, collapse = ";"))
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
IsBioID %<o% (gsub(" |_|-|\\.", "", toupper(Param$Type)) == "BIOID")
IsBioID2 %<o% FALSE
if (IsBioID) {
  wbiot <- grep("biot", Modifs$"Full name", ignore.case = TRUE)
  .obj <- unique(c(.obj, "wbiot"))
  l <- length(wbiot)
  if (length(wbiot)) {
    if (l == 1) { tmp <- Modifs$"Full name"[wbiot] } else { tmp <- paste0(paste(Modifs$"Full name"[wbiot[1:(l-1)]], collapse = "\", \""), "\" and \"", Modifs$"Full name"[wbiot[l]]) }
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
      PG[["Biot. peptides count"]] <- sapply(strsplit(PG[["Biot. peptide IDs"]], ";"), length)
      PG[["Biot. peptides [%]"]] <- round(100*PG[["Biot. peptides count"]]/PG$"Peptides count", 1)
      IsBioID2 <- TRUE
    } else { warning("I could not find any biotinylated peptides!") }
  } else { warning("I could not identify any biotinylated PTMs in the modifications table, did you include them in the search?") }
}
# First sequence
PG$"Sequence (1st accession)" <- db$Sequence[match(sapply(strsplit(PG$"Leading protein IDs", ";"), function(x) { x[[1]] }),
                                                   db$"Protein ID")]

# Number of spectra, evidences and peptides per sample:
source(parSrc, local = FALSE)
clusterCall(parClust, function() library(proteoCraft))
clusterCall(parClust, function() library(reshape))
clusterCall(parClust, function() library(data.table))
temp_PG <- data.frame(id = PG$id, Accession1 = sapply(strsplit(PG$"Leading protein IDs", ";"), function(x) { unlist(x)[1] }))
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
  temp <- temp[which(sapply(temp$MSMSIDs, length) > 0),] # Remove Match-Between-Runs evidences (no MS/MS)
  temp <- listMelt(temp$MSMSIDs, temp$PG_id, c("MSMSIDs", "PG_id"))
  temp <- do.call(data.frame, aggregate(temp$MSMSIDs, list(temp$PG_id), function(x) {
    x <- unique(x)
    return(c(Count = length(x), List = list(x)))
  }))
  temp$x.Count <- unlist(temp$x.Count)
  temp$Pasted <- sapply(temp$x.List, function(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") })
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
    temp1$Pasted <- sapply(temp1$x.List, function(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") })
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
      temp3$Pasted <- sapply(temp3$x.List, function(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") })
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
  max(sapply(kol, function(y) { sum(proteoCraft::is.all.good(unlist(tst[x, y])) > 0) }) >= PepFoundInAtLeastGrp)
}) > 0
Pep2Use %<o% which(tst)
#length(Pep2Use)/nrow(pep)

# Currently 9 methods are implemented, but the last 3 are just for testing and do not create all columns required for the script to complete
# (These could be added relatively easily if necessary though)
QuantData <- setNames(paste0("quant.data", 1:length(QuantMethods)), QuantMethods)
QuantMethods_all <- FALSE
.obj <- unique(c(.obj, "QuantData", "QuantMethods_all"))
expscol <- paste0("log10(Expr.) - ", RSA$values)
# Weights:
# - Higher for peptides with low intra-sample group CV on average
# - Higher for peptides with low PEP
if (QuantUMS) {
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
.obj <- unique(c(.obj, "Mod.Excl.is.strict", "Discard.unmod"))
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
    quant.data3 <- Prot.Quant(Prot = PG,  Mode = "PreferUnique",
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
    quant.data4 <- Prot.Quant(Prot = PG,  Mode = "PreferUnique",
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
  if (!require("iq", quietly = TRUE)) { install.packages("iq") }
  library("iq")
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
  temp <- lapply(strsplit(pep$`Protein group IDs`[Pep2Use2], ";"), as.integer)
  temp <- listMelt(temp, pep$id[Pep2Use2])
  temp[, c("Modified sequence", kol)] <- pep[match(temp$id, pep$id), c("Modified sequence", kol)]
  colnames(temp) <- gsub(topattern(pep.ref[length(pep.ref)]), "", colnames(temp))
  temp <- reshape2::melt(temp, id.vars = c("PG", "id", "Modified sequence"))
  temp <- temp[which((is.all.good(temp$value, 2))&(temp$value > 0)),]
  temp$protein_list <- PG$"Leading protein IDs"[match(temp$PG, PG$id)]
  colnames(temp)[which(colnames(temp) == "variable")] <- "sample_list"
  temp$id <- apply(temp[, c("id", "Modified sequence")], 1, function(x) { paste(as.character(x), collapse = "") })
  temp$quant <- log2(temp$value)
  temp$PG <- NULL
  temp$"Modified sequence" <- NULL
  temp$value <- NULL
  temp2 <- iq::create_protein_list(temp)
  quant.data7 <- iq::create_protein_table(temp2)
  quant.data7 <- quant.data7$estimate/log2(10)
  quant.data7 <- as.data.frame(quant.data7)
  colnames(quant.data7) <- paste0("log10(Expr.) - ", colnames(quant.data7))
  quant.data7$"Leading protein IDs" <- names(temp2)
  w <- which(!PG$"Leading protein IDs" %in% quant.data7$"Leading protein IDs")
  temp <- as.data.frame(matrix(rep(NA, length(w)*length(expscol)), ncol = length(expscol)))
  colnames(temp) <- expscol
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
  NormFact1 <- sapply(1:length(QuantMethods), function(i) {
    tmp <- get(QuantData[i])
    tmp$"Leading protein IDs" <- PG$"Leading protein IDs"
    tmp2 <- tmp[,expscol]
    tmp2 <- rowMeans(tmp2, na.rm = TRUE)
    w <- which(PG$"Leading protein IDs" %in% tmp$"Leading protein IDs")
    res <- rep(NA, nrow(PG))
    res[w] <- tmp2[match(PG$"Leading protein IDs"[w], tmp$`Leading protein IDs`)]
    return(res)
  })
  NormFact2 <- rowMeans(NormFact1, na.rm = TRUE)
  for (i in 1:length(QuantMethods)) {
    tmp1 <- get(QuantData[i])
    tmp1$"Leading protein IDs" <- PG$"Leading protein IDs"
    tmp1 <- tmp1[, c(expscol, "Leading protein IDs")]
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
    tmp3[, gsub(topattern("log10(Expr.) - "), "", expscol)] <- sweep(tmp3[, gsub(topattern("log10(Expr.) - "), "", expscol)], 1, NormFact1[m,i], "-")+NormFact2[m]
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
    tst <- sapply(RSA$names, function(x) { length(get(substr(x, 1, 3))) })
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
    ggsave(paste0(dir, "/", ttl1a, ".jpeg"), plot1, dpi = 300)
    ggsave(paste0(dir, "/", ttl1a, ".pdf"), plot1, dpi = 300)
    ReportCalls <- AddPlot2Report(Plot = plot1, Title = ttl1a)
    ttl2 <- paste0("LFQ method comparison (renormalized): ", comb[i, 1], " VS ", comb[i, 2])
    plot2 <- ggplot(tmp1) +
      geom_scattermore(aes(x = X, y = Y, colour = `Razor + unique peptides`),
                       shape = 16, size = 0.1, alpha = 0.1) +
      scale_color_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() + facet_grid(form) + xlab(comb[i, 1]) + ylab(comb[i, 2]) + ggtitle(ttl2)
    print(plot2) # This type of QC plot does not need to pop up, the side panel is fine
    ttl2a <- gsub("\\:", " -", ttl1)
    ggsave(paste0(dir, "/", ttl2a, ".jpeg"), plot2, dpi = 300)
    ggsave(paste0(dir, "/", ttl2a, ".pdf"), plot2, dpi = 300)
    ReportCalls <- AddPlot2Report(Plot = plot2, Title = ttl2a)
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
  ggsave(paste0(dir, "/", ttla, ".jpeg"), plot, dpi = 300)
  ggsave(paste0(dir, "/", ttla, ".pdf"), plot, dpi = 300)
  ReportCalls <- AddPlot2Report(Title = ttla)
}
ReportCalls$Calls <- append(ReportCalls$Calls,
                            paste0("body_add_fpar(Report, fpar(ftext(\"Protein groups quantitation done using method: ", Param$QuantMeth, "\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))"))
if (!exists(QuantData[match(Param$QuantMeth, QuantMethods)])) { # (for when rerunning script without rerunning quant methods)
  load(paste0(QuantData[match(Param$QuantMeth, QuantMethods)], ".RData"))
}
quant.data <- get(QuantData[match(Param$QuantMeth, QuantMethods)])
Prot.Expr.Root %<o% c(Original = "log10(Expr.) - ")
Prot.Rat.Root %<o% pep.ratios.ref[length(pep.ratios.ref)]
#write.csv(quant.data, file = "Quantitative data.csv", row.names = FALSE)
.obj <- unique(c(.obj, "quant.data", "Prot.Expr.Root", "Prot.Rat.Root"))
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
    .obj <- unique(c(.obj, "Prot.Rat.Root.BckUp", "pep.ratios.ref.BckUp"))
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
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Re-normalize protein group expression values
Norma.Prot.Ratio.classic %<o% FALSE
Norma.Prot.Ratio.to.Biot %<o% FALSE
Norma.Prot.Ratio.to.proteins %<o% FALSE
if (("Norma.Prot.Ratio" %in% colnames(Param))&&(Param$Norma.Prot.Ratio)) {
  quant.data.norm <- quant.data
  ExpKol <- paste0(Prot.Expr.Root, Exp.map$Ref.Sample.Aggregate)
  RatKol <- paste0(Prot.Rat.Root, Exp.map$Ref.Sample.Aggregate)
  whichRat <- which(RatKol %in% colnames(quant.data.norm))
  RatKol <- RatKol[whichRat]
  RefKol <- paste0(Prot.Expr.Root, RRG$values, ".REF")
  Norma.Prot.Ratio.classic <- TRUE # The default
  if ((IsBioID2)&&("Norma.Prot.Ratio.to.Biot" %in% colnames(Param))&&(is.logical(Param$Norma.Prot.Ratio.to.Biot))&&(Param$Norma.Prot.Ratio.to.Biot)) {
    nrmFlt <- which(PG$"Biot. peptides count" > 0)
    if (length(w)) {
      Norma.Prot.Ratio.classic <- Norma.Prot.Ratio.to.proteins <- FALSE
      Norma.Prot.Ratio.to.Biot <- TRUE
      txTmp <- "biotinylated protein groups"
    }
  }
  if (("Norma.Prot.Ratio.to.proteins" %in% colnames(Param))&&
      (is.character(Param$Norma.Prot.Ratio.to.proteins))&&
      (!Param$Norma.Prot.Ratio.to.proteins %in% c("", "NA"))&&
      (!Norma.Prot.Ratio.to.Biot)) {
    Prot.Ratio.ref.Acc %<o% unique(unlist(strsplit(Param$Norma.Prot.Ratio.to.proteins, ";")))
    Prot.Ratio.ref.Acc <- Prot.Ratio.ref.Acc[which(Prot.Ratio.ref.Acc %in% db$"Protein ID")]
    if (length(Prot.Ratio.ref.Acc)) {
      Prot.Ratio.ref.Acc <- Prot.Ratio.ref.Acc[which(Prot.Ratio.ref.Acc %in% unique(unlist(strsplit(PG$"Protein IDs", ";"))))] # Having "Leading protein IDs" here was too stringent
      Norma.Prot.Ratio.classic <- FALSE
      l <- length(Prot.Ratio.ref.Acc)
      if (l) {
        temp <- strsplit(PG$"Protein IDs", ";")
        temp <- listMelt(temp, PG$id)
        temp$Norm <- temp$value %in% Prot.Ratio.ref.Acc
        temp <- temp$L1[which(temp$Norm)]
        nrmFlt <- which(PG$id %in% temp)
        if (length(nrmFlt)) {
          Norma.Prot.Ratio.classic <- Norma.Prot.Ratio.to.Biot <- FALSE
          Norma.Prot.Ratio.to.proteins <- TRUE
          txTmp <- paste0("the following proteins: ",
                          paste(db$`Full Name`[match(Prot.Ratio.ref.Acc, db$`Protein ID`)],
                                collapse = "-"))
        }
      }
    }
  }
  if (Norma.Prot.Ratio.classic) {
    nrmFlt <- 1:nrow(quant.data)
    Norma.Prot.Ratio.to.Biot <- Norma.Prot.Ratio.to.proteins <- FALSE
    txTmp <- "all protein groups"
  }
  #
  # nrmFlt is the rows filter we will use to ensure that we normalize to the correct proteins:
  # - Norma.Prot.Ratio.classic: the filter includes all rows
  # - Norma.Prot.Ratio.to.Biot: the filter includes biotinylated proteins
  # - Norma.Prot.Ratio.to.proteins: the filter includes specific user-defined protein(s)
  #
  tst <- sum(c(Norma.Prot.Ratio.to.Biot, Norma.Prot.Ratio.to.proteins, Norma.Prot.Ratio.classic))
  if (tst) {
    stopifnot(tst == 1)
    globMed1 <- median(is.all.good(unlist(quant.data.norm[#nrmFlt # For global scales, we use all
      , c(ExpKol, RefKol)])))
    # First normalize within groups
    for (grp in RG$values) { #grp <- RG$values[1]
      smpls <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[RG$column]] == grp)]
      xpKol <- paste0(Prot.Expr.Root, smpls)
      rtKol <- paste0(Prot.Rat.Root, smpls)
      wRat <- which(rtKol %in% colnames(quant.data.norm))
      rtKol <- rtKol[wRat]
      rfgrps <- unique(Exp.map[which(Exp.map[[RG$column]] == grp), RRG$column])
      rfKol <- paste0(Prot.Expr.Root, rfgrps, ".REF")
      rtMed <- apply(quant.data.norm[nrmFlt, rtKol], 2, function(x) { median(is.all.good(x)) })
      xpMed1 <- median(is.all.good(unlist(quant.data.norm[nrmFlt, xpKol]))) # Before, for re-scaling
      # 1a - apply normalization at ratios level...
      quant.data.norm[, rtKol] <- sweep(quant.data.norm[, rtKol, drop = FALSE],
                                        2,
                                        rtMed,
                                        "-")
      # 1b - ... and probagate it to the corresponding columns
      quant.data.norm[, xpKol[wRat]] <- sweep(quant.data.norm[, xpKol[wRat], drop = FALSE],
                                              2,
                                              rtMed/log2(10), # Convert from base 2 to 10 
                                              "-")
      # 2 - restore original scaling
      xpMed2 <- median(is.all.good(unlist(quant.data.norm[nrmFlt, xpKol]))) # After, for re-scaling
      quant.data.norm[, c(xpKol, rfKol)] <- sweep(quant.data.norm[, c(xpKol, rfKol)],
                                                  2,
                                                  xpMed1-xpMed2,
                                                  "+")
    }
    # ... then re-normalize globally.
    globMed2 <- median(is.all.good(unlist(quant.data.norm[#nrmFlt
      , c(ExpKol, RefKol)])))
    quant.data.norm[, c(ExpKol, RefKol)] <- sweep(quant.data.norm[, c(ExpKol, RefKol)],
                                                  2,
                                                  globMed1-globMed2,
                                                  "+")
    # Sometimes, if you are recycling to e.g. a single protein, you may get randomly singular values in one row
    # This breaks stats so we need to fix that.
    #  - Expression:
    tst <- apply(quant.data.norm[, c(ExpKol, RefKol)], 1, function(x) {
      x <- x[which(!is.na(x))]
      length(unique(x))
    })
    if (1 %in% tst) {
      #aggregate(tst, list(tst), length)
      # Add a small value
      w <- which(tst > 0)
      quant.data.norm[w, c(ExpKol, RefKol)] <- quant.data.norm[w, c(ExpKol, RefKol)] + runif(length(w)*length(c(ExpKol, RefKol)), min = 0, max = 1e-10)
    }
    #  - Ratios:
    tst <- apply(quant.data.norm[, RatKol], 1, function(x) {
      x <- x[which(!is.na(x))]
      length(unique(x))
    })
    if (1 %in% tst) {
      #aggregate(tst, list(tst), length)
      # Add a small value
      w <- which(tst > 0)
      quant.data.norm[w, RatKol] <- quant.data.norm[w, RatKol] + runif(length(w)*length(RatKol), min = 0, max = 1e-10)
    }
    #
    ttl <- paste0("re-normalisation (median of ",
                  c("Biot. PGs", "specific proteins", "all PGs")[which(c(Norma.Prot.Ratio.to.Biot,
                                                                         Norma.Prot.Ratio.to.proteins,
                                                                         Norma.Prot.Ratio.classic))],
                  ")")
    DatAnalysisTxt <- gsub("\\.\\.\\.$", paste0(" and re-normalized to the median value of ", txTmp, "."), DatAnalysisTxt)
    ReportCalls$Calls <- append(ReportCalls$Calls,
                                paste0("body_add_fpar(Report, fpar(ftext(\"Re-normalizing Protein groups-level expression values  the median value of ",
                                       txTmp, "\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))"))
  }
  if (Norma.Prot.Ratio.classic) {
    if (!"Adv.Norma.Prot.Intens" %in% colnames(Param)) {
      ObjNm <- "Norma.Prot.Ratio.Adv"
      if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
        msg <- "Do you want to apply Levenberg-Marquardt (slow) re-normalisation at protein groups level?"
        tmp <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
        if (is.na(tmp)) { tmp <- FALSE }
        ObjNm %<c% tmp
        tmp <- AllAnsw[1,]
        tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
        tmp$Value <- list(get(ObjNm))
        m <- match(ObjNm, AllAnsw$Parameter)
        if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
      }
    } else { Norma.Prot.Ratio.Adv <- Param$Adv.Norma.Prot.Intens }
    if (Norma.Prot.Ratio.Adv) {
      ReportCalls$Calls <- append(ReportCalls$Calls,
                                  "body_add_fpar(Report, fpar(ftext(\"Applying Levenberg-Marquardt normalisation at Protein groups-level...\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
      ids <- PG$id[nrmFlt]
      # Let's first normalize per ratio group...
      clusterExport(parClust, list("quant.data.norm", "Exp.map", "RG", "RRG", "nrmFlt", "is.all.good",
                                   "Prot.Expr.Root", "Prot.Rat.Root", "ids"), envir = environment())
      tmpNrm <- parLapply(parClust, RG$values, function(grp) { #grp <- RG$values[1]
        smpls <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[RG$column]] == grp)]
        xpKol <- paste0(Prot.Expr.Root, smpls)
        rfgrps <- unique(Exp.map[which(Exp.map[[RG$column]] == grp), RRG$column])
        rfKol <- paste0(Prot.Expr.Root, rfgrps, ".REF")
        grpMed1 <- median(is.all.good(unlist(quant.data.norm[#nrmFlt
          , c(xpKol, rfKol)])))
        tmp <- quant.data.norm[nrmFlt, xpKol]
        tmp$id <- ids
        xpMed1 <- apply(quant.data.norm[nrmFlt, xpKol], 2, median, na.rm = TRUE)
        tmp2 <- proteoCraft::AdvNorm.IL(tmp, "id", xpKol, TRUE, 5)
        tmp2$id <- NULL
        colnames(tmp2) <- xpKol
        xpMed2 <- apply(tmp2, 2, median, na.rm = TRUE)
        tmp2 <- sweep(quant.data.norm[, xpKol],
                      2, 
                      xpMed2-xpMed1,
                      "+")
        # References
        rrg <- RRG$values[which(sapply(RRG$values, function(x) {
          unique(Exp.map[which(Exp.map[[RRG$column]] == x), RG$column] == grp)
        }))]
        #xpMed2-xpMed1
        rfMed1 <- apply(quant.data.norm[nrmFlt, rfKol, drop = FALSE], 2, median, na.rm = TRUE)
        rfMed2 <- sapply(rrg, function(x) { #x <- rrg[1]
          x <- Exp.map$Ref.Sample.Aggregate[which((Exp.map[[RRG$column]] == x)&(Exp.map$Reference))]
          median(rowMeans(tmp2[, paste0(Prot.Expr.Root, x), drop = FALSE], na.rm = TRUE), na.rm = TRUE)
        })
        #rfMed2-rfMed1
        tmp3 <- sweep(quant.data.norm[, rfKol, drop = FALSE],
                      2, 
                      rfMed2-rfMed1,
                      "+")
        resXp <- cbind(tmp2, tmp3)
        # Re-apply original scaling
        grpMed2 <- median(is.all.good(unlist(resXp[#nrmFlt
          , c(xpKol, rfKol)])))
        resXp <- sweep(resXp,
                       2, 
                       grpMed1-grpMed2,
                       "+")
        return(resXp)
      })
      tmpNrm2 <- tmpNrm
      tmpNrm <- do.call(cbind, tmpNrm)
      quant.data.norm[, colnames(tmpNrm)] <- tmpNrm
      if (length(RG$values) > 1) {
        # ... then between ratio groups...
        tmpNrm2 <- as.data.frame(sapply(tmpNrm2, function(x) { #x <- tmpNrm[[1]]
          k <- grep("\\.REF$", colnames(x), invert = TRUE, value = TRUE)
          rowMeans(x[nrmFlt, k], na.rm = TRUE)
        }))
        colnames(tmpNrm2) <- RG$values
        globMed1 <- median(is.all.good(unlist(tmpNrm2)))
        xpMed1 <- apply(tmpNrm2, 2 , median, na.rm = TRUE)
        tmpNrm2$id <- ids
        tmp2 <- proteoCraft::AdvNorm.IL(tmpNrm2, "id", RG$values, TRUE, 5)
        tmp2$id <- NULL
        colnames(tmp2) <- RG$values
        xpMed2 <- apply(tmp2, 2, median, na.rm = TRUE)
        globMed2 <- median(is.all.good(unlist(tmp2)))
        xpMed1 <- c(sapply(Exp.map[[RG$column]], function(x) { xpMed1[x] }),
                    sapply(RRG$values, function(x) {
                      xpMed1[Exp.map[match(x, Exp.map[[RRG$column]]) # Only works because RRGs are always contained in RGs
                                     , RG$column]]
                    }))
        xpMed2 <- c(sapply(Exp.map[[RG$column]], function(x) { xpMed2[x] }),
                    sapply(RRG$values, function(x) {
                      xpMed2[Exp.map[match(x, Exp.map[[RRG$column]]) # Only works because RRGs are always contained in RGs
                                     , RG$column]]
                    }))
        quant.data.norm[, c(ExpKol, RefKol)] <- sweep(quant.data.norm[, c(ExpKol, RefKol)],
                                                      2, xpMed2+globMed1-xpMed1-globMed2,
                                                      "+")
      }
      # Re-calculate ratios
      Rat <- as.data.frame(sapply(Exp.map$Ref.Sample.Aggregate[whichRat], function(x) {
        #x <- Exp.map$Ref.Sample.Aggregate[whichRat][1]
        xp <- paste0(Prot.Expr.Root, x)
        rf <- paste0(Prot.Expr.Root, Exp.map[match(x, Exp.map$Ref.Sample.Aggregate), RRG$column], ".REF")
        x1 <- quant.data.norm[[xp]]
        x0 <- quant.data.norm[[rf]]
        return((x1-x0)/log10(2))
      }))
      colnames(Rat) <- RatKol
      quant.data.norm[, RatKol] <- Rat
      #
      ttl <- "re-normalisation (Levenberg-Marquardt)"
      DatAnalysisTxt <- gsub(" and re-normalized to the median value of .+",
                             paste0(" and re-normalized using the Levenberg-Marquardt procedure."), DatAnalysisTxt)
      ReportCalls$Calls <- append(ReportCalls$Calls,
                                  "body_add_fpar(Report, fpar(ftext(\"Done!\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
    }
  }
  if (Norma.Prot.Ratio.classic||Norma.Prot.Ratio.to.Biot||Norma.Prot.Ratio.to.proteins) {
    # Once this is done, update ref-to-ref ratios
    res <- make_RefRat(data = quant.data.norm,
                       int.root = Prot.Expr.Root,
                       rat.root = Prot.Rat.Root,
                       logInt = 10)
    quant.data.norm[, colnames(res)] <- res
    #... ok, it's done, but what about peptides?
    # If we want to normalize, e.g., peptide fold changes to parent protein fold changes,
    # then we must first propagate the normalisations to the peptides
    kPpE <- c(paste0(pep.ref[length(pep.ref)], Exp.map$Ref.Sample.Aggregate),
              paste0(pep.ref[length(pep.ref)], unique(Exp.map[[RRG$column]]), ".REF"))
    kPpR <- paste0(pep.ratios.ref[length(pep.ratios.ref)], Exp.map$Ref.Sample.Aggregate)
    kPrE <- c(paste0(Prot.Expr.Root, Exp.map$Ref.Sample.Aggregate),
              paste0(Prot.Expr.Root, unique(Exp.map[[RRG$column]]), ".REF"))
    kPrR <-  paste0(Prot.Rat.Root, Exp.map$Ref.Sample.Aggregate)
    kPpR <- kPpR[which(kPpR %in% colnames(pep))]
    kPrR <- kPrR[which(kPrR %in% colnames(quant.data))]
    pep.quant.data.norm <- pep[, c(kPpE, kPpR)]
    diffE <- 10^sapply(kPrE, function(x) {
      mean(is.all.good(unlist(quant.data.norm[x] - quant.data[x])))
    })
    diffR <- sapply(kPrR, function(x) {
      mean(is.all.good(unlist(quant.data.norm[x] - quant.data[x])))
    })
    pep.quant.data.norm[, kPpE] <- sweep(pep.quant.data.norm[, kPpE], 2, diffE, "*")
    pep.quant.data.norm[, kPpR] <- sweep(pep.quant.data.norm[, kPpR], 2, diffR, "+")
    # Also pep.ref.ratios
    if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
      pep.Ref.Ratios.norm %<o% NULL
    }
    if (Param$Ratios.Thresholds == threshMsg) {
      pep.Ref.Ratios.norm %<o% make_RefRat(data = pep.quant.data.norm)
    }
    #
    # And now check and visualize
    dirPG <- paste0(wd, "/Workflow control/Protein groups")
    if (!dir.exists(dirPG)) { dir.create(dirPG, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dirPG))
    dirPep <- paste0(wd, "/Workflow control/Peptides")
    if (!dir.exists(dirPep)) { dir.create(dirPep, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dirPep))
    nrmPlots <- list()
    # - PG expression values
    Samples <- cleanNms(gsub(topattern(Prot.Expr.Root), "", ExpKol))
    temp1 <- data.table(quant.data[, ExpKol])
    colnames(temp1) <- Samples
    temp1 <- melt(temp1, measure.vars = Samples)
    temp2 <- data.table(quant.data.norm[, ExpKol])
    colnames(temp2) <- Samples
    temp2 <- melt(temp2, measure.vars = Samples)
    temp1$Norm <- "Original"
    temp2$Norm <- "Normalised"
    temp <- rbind(temp1, temp2)
    rm(temp1, temp2)
    temp$Sample <- factor(temp$variable, levels = Samples)
    temp$variable <- NULL
    temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
    temp <- temp[which(is.all.good(temp$value, 2)),]
    ttlI1 <- paste0("Protein groups ", ttl, ", expression")
    intPlot1 <- ggplot(temp) +
      geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(Norm~.) + ggtitle(ttlI1) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    #print(intPlot1)
    nrmPlots[["PG_int"]] <- list(Path = paste0(dirPG, "/", ttlI1),
                                 Plot = plotEval(intPlot1),
                                 Ext = "jpeg")
    #ggsave(paste0(dirPG, "/", ttlI1, ".jpeg"), intPlot1, dpi = 300)
    #ggsave(paste0(dirPG, "/", ttlI1, ".pdf"), intPlot1, dpi = 300)
    #ReportCalls <- AddPlot2Report(Plot = intPlot1, Title = ttlI1)
    #
    # - PG ratios
    Samples <- cleanNms(gsub(topattern(Prot.Rat.Root), "", RatKol))
    temp1 <- data.table(quant.data[, RatKol])
    colnames(temp1) <- Samples
    temp1 <- melt(temp1, measure.vars = Samples)
    temp2 <- data.table(quant.data.norm[, RatKol])
    colnames(temp2) <- Samples
    temp2 <- melt(temp2, measure.vars = Samples)
    temp1$Norm <- "Original"
    temp2$Norm <- "Normalised"
    temp <- rbind(temp1, temp2)
    rm(temp1, temp2)
    temp$Sample <- factor(temp$variable, levels = Samples)
    temp$variable <- NULL
    temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
    temp <- temp[which(is.all.good(temp$value, 2)),]
    ttlR1 <- paste0("Protein groups ", ttl, ", ratios")
    ratPlot1 <- ggplot(temp) +
      geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(Norm~.) + ggtitle(ttlR1) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    #print(ratPlot1)
    nrmPlots[["PG_rat"]] <- list(Path = paste0(dirPG, "/", ttlR1),
                                 Plot = plotEval(ratPlot1),
                                 Ext = "jpeg")
    #ggsave(paste0(dirPG, "/", ttlR1, ".jpeg"), ratPlot1, dpi = 300)
    #ggsave(paste0(dirPG, "/", ttlR1, ".pdf"), ratPlot1, dpi = 300)
    #ReportCalls <- AddPlot2Report(Plot = ratPlot1, Title = ttlR1)
    #
    # - Peptide intensities
    Samples <- cleanNms(gsub(topattern(pep.ref[length(pep.ref)]), "", kPpE))
    temp1 <- data.table(pep[, kPpE])
    colnames(temp1) <- Samples
    temp1 <- melt(temp1, measure.vars = Samples)
    temp2 <- data.table(pep.quant.data.norm[, kPpE])
    colnames(temp2) <- Samples
    temp2 <- melt(temp2, measure.vars = Samples)
    temp1$Norm <- "Original"
    temp2$Norm <- "Normalised"
    temp <- rbind(temp1, temp2)
    rm(temp1, temp2)
    temp$Sample <- factor(temp$variable, levels = Samples)
    temp$variable <- NULL
    temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
    temp$value <- suppressWarnings(log10(temp$value)) # Peptide intensities are not log-transformed!
    temp <- temp[which(is.all.good(temp$value, 2)),]
    ttlI2 <- paste0("Peptides ", ttl, " from PGs, intensity")
    intPlot2 <- ggplot(temp) +
      geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(Norm~.) + ggtitle(ttlI2) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    #print(intPlot2)
    nrmPlots[["Pep_int"]] <- list(Path = paste0(dirPep, "/", ttlI2),
                                  Plot = plotEval(intPlot2),
                                  Ext = "jpeg")
    #ggsave(paste0(dirPep, "/", ttlI2, ".jpeg"), intPlot2, dpi = 300)
    #ggsave(paste0(dirPep, "/", ttlI2, ".pdf"), intPlot2, dpi = 300)
    #ReportCalls <- AddPlot2Report(Plot = intPlot2, Title = ttlI2)
    #
    # - Peptides ratios
    Samples <- cleanNms(gsub(topattern(pep.ratios.ref[length(pep.ratios.ref)]), "", kPpR))
    temp1 <- data.table(pep[, kPpR])
    colnames(temp1) <- Samples
    temp1 <- melt(temp1, measure.vars = Samples)
    temp2 <- data.table(pep.quant.data.norm[, kPpR])
    colnames(temp2) <- Samples
    temp2 <- melt(temp2, measure.vars = Samples)
    temp1$Norm <- "Original"
    temp2$Norm <- "Normalised"
    temp <- rbind(temp1, temp2)
    rm(temp1, temp2)
    temp$Sample <- factor(temp$variable, levels = Samples)
    temp$variable <- NULL
    temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
    temp <- temp[which(is.all.good(temp$value, 2)),]
    ttlR2 <- paste0("Peptides ", ttl, " from PGs, ratios")
    ratPlot2 <- ggplot(temp) +
      geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(Norm~.) + ggtitle(ttlR2) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    #print(ratPlot2)
    nrmPlots[["Pep_rat"]] <- list(Path = paste0(dirPep, "/", ttlR2),
                                  Plot = plotEval(ratPlot2),
                                  Ext = "jpeg")
    #ggsave(paste0(dirPep, "/", ttlR2, ".jpeg"), ratPlot2, dpi = 300)
    #ggsave(paste0(dirPep, "/", ttlR2, ".pdf"), ratPlot2, dpi = 300)
    #ReportCalls <- AddPlot2Report(Plot = ratPlot2, Title = ttlR2)
    #
    nrmPlots2 <- lapply(nrmPlots, function(x) {
      x$Ext <- "pdf"
      return(x)
    })
    nrmPlots <- c(nrmPlots, nrmPlots2)
    # Save plots
    tst <- parLapply(parClust, nrmPlots, function(x) { ggplot2::ggsave(paste0(x$Path, ".", x$Ext), x$Plot, dpi = 300) })
    #
    if (exists("acceptPGNorm")) {
      acceptPGNorm <- as.logical(acceptPGNorm)
      if ((length(acceptPGNorm) != 1)||(is.na(acceptPGNorm))) { acceptPGNorm <- TRUE }
    } else { acceptPGNorm %<o% TRUE }
    msg <- "Accept re-normalisation results?"
    appNm <- paste0(dtstNm, " - PG re-normalisation")
    ui <- fluidPage(
      useShinyjs(),
      extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
      titlePanel(tag("u", "Protein Groups (and peptides) re-normalisation"),
                 appNm),
      br(),
      checkboxInput("acceptPGNorm", msg, acceptPGNorm),
      actionButton("saveBtn", "Save"),
      # fluidRow(column(4, withSpinner(plotOutput("pgIntensities"))),
      #          column(4, withSpinner(plotOutput("pgRatios")))),
      #br(),
      # fluidRow(column(4, withSpinner(plotOutput("pepIntensities"))),
      #          column(4, withSpinner(plotOutput("pepRatios")))),
      fluidRow(column(6, withSpinner(imageOutput("pgIntensities", width = "100%", height = "100%"))),
               column(6, withSpinner(imageOutput("pgRatios", width = "100%", height = "100%")))),
      br(),
      br(),
      br(),
      fluidRow(column(6, withSpinner(imageOutput("pepIntensities", width = "100%", height = "100%"))),
               column(6, withSpinner(imageOutput("pepRatios", width = "100%", height = "100%")))),
      br(),
      br()
    )
    IMGs <- c(paste0(dirPG, "/", c(ttlI1, ttlR1), ".jpeg"),
              paste0(dirPep, "/", c(ttlI2, ttlR2), ".jpeg"))
    IMGsDims <- as.data.frame(t(parSapply(parClust, IMGs, function(x) { #x <- IMGs[1]
      a <- jpeg::readJPEG(x)
      setNames(dim(a)[1:2], c("height", "width"))
    })))
    IMGsDims$height <- screenRes$width*0.35*IMGsDims$height/max(IMGsDims$height)
    IMGsDims$width <- screenRes$width*0.35*IMGsDims$width/max(IMGsDims$width)
    server <- function(input, output, session) {
      # output$pgIntensities <- renderPlot(intPlot1)
      # output$pgRatios <- renderPlot(ratPlot1)
      # output$pepIntensities <- renderPlot(intPlot2)
      # output$pepRatios <- renderPlot(ratPlot2)
      output$pgIntensities <- renderImage({
        list(src = IMGs[1], height = IMGsDims$height[1], width = IMGsDims$width[1])
      }, deleteFile = FALSE)
      output$pgRatios <- renderImage({
        list(src = IMGs[2], height = IMGsDims$height[2], width = IMGsDims$width[2])
      }, deleteFile = FALSE)
      output$pepIntensities <- renderImage({
        list(src = IMGs[3], height = IMGsDims$height[3], width = IMGsDims$width[3])
      }, deleteFile = FALSE)
      output$pepRatios <- renderImage({
        list(src = IMGs[4], height = IMGsDims$height[4], width = IMGsDims$width[4])
      }, deleteFile = FALSE)
      #
      observeEvent(input$acceptPGNorm, { assign("acceptPGNorm", as.logical(input[["acceptPGNorm"]]), envir = .GlobalEnv) })
      observeEvent(input$saveBtn, { stopApp() })
      #observeEvent(input$cancel, { stopApp() })
      session$onSessionEnded(function() { stopApp() })
    }
    eval(parse(text = runApp), envir = .GlobalEnv)
    #
    if (acceptPGNorm) {
      quant.data <- quant.data.norm
      pep[, colnames(pep.quant.data.norm)] <- pep.quant.data.norm
      if (Param$Ratios.Thresholds == threshMsg) {
        pep.Ref.Ratios <- pep.Ref.Ratios.norm
        pep[, colnames(pep.Ref.Ratios)] <- pep.Ref.Ratios
      }
    }
  }
  ReportCalls <- AddSpace2Report()
}

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
w <- which(sapply(a, function(x) {length(unique(test[[x]]))}) > 1)
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
  testI[[v]] <- sapply(1:nbinz, function(x) {
    sum((test$value[wv] > binz[x])&(test$value[wv] <= binz[x+1]))
  })
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
ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")
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
w <- which(sapply(a, function(x) {length(unique(test[[x]]))}) > 1)
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
  testR[[v]] <- sapply(1:nbinz, function(x) {
    sum((test$value[wv] > binz[x])&(test$value[wv] <= binz[x+1]))
  })
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
ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")
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
      m <- match(chan, em$Isobaric.label)
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
  BefAft <- sapply(VPAL$values, function(x) {
    x <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)]
    x <- x[which(x %in% RSA$values[WhInColNms])]
    x <- median(unlist(BefAft[, x]), na.rm = TRUE)
    return(x)
  })
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
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
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
    ReportCalls$Calls <- append(ReportCalls$Calls,
                                "body_add_fpar(Report, fpar(ftext(\"Removing proteins not identified in the provided true-discovery filter.\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
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
  PG$"Peptides count" <- sapply(strsplit(PG$"Peptide IDs", ";"), length)
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
Script <- readLines(ScriptPath)

#### Code chunk - samples Pearson correlation heatmap
if (LocAnalysis) { prtRfRoot %<o% Prot.Expr.Root2 } else { prtRfRoot %<o% Prot.Expr.Root }
dir <- paste0(wd, "/Pearson correlation map")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
g <- paste0(prtRfRoot, RSA$values)
w <- which(g %in% colnames(PG))
g <- g[w]
smpls <- cleanNms(RSA$values[w])
corMap <- sapply(1:length(g), function(x) {
  sapply(1:length(g), function(y) {
    if (x > y) {
      x <- PG[[g[x]]]
      y <- PG[[g[y]]]
      w <- which((is.all.good(x, 2))&(is.all.good(y, 2)))
      return(cor(x[w], y[w], method = "pearson"))
    } else { return(NA) }
  })
})
rownames(corMap) <- colnames(corMap) <- smpls
corMap <- reshape2::melt(corMap, measure.vars = smpls)
colnames(corMap) <- c("Sample 2", "Sample 1", "Pearson corr.")
corMap <- corMap[which(!is.na(corMap$"Pearson corr.")),]
tmp <- PG[, g[2:length(g)]]
colnames(tmp) <- smpls[2:length(g)]
tmp <- reshape2::melt(tmp, measure.vars = smpls[2:length(g)])
scattrMap <- data.frame("Sample 1" = smpls[1],
                        "Sample 2" = tmp$variable,
                        "X" = PG[[g[1]]],
                        "Y" = tmp$value,
                        check.names = FALSE)
if (length(smpls) > 2) {
  for (i in 2:(length(smpls)-1)) {
    tmp <- tmp[(which(tmp$variable == smpls[i+1])[1]):nrow(tmp),]
    tmp2 <- data.frame("Sample 1" = smpls[i],
                       "Sample 2" = tmp$variable,
                       "X" = PG[[g[i]]],
                       "Y" = tmp$value,
                       check.names = FALSE)
    scattrMap <- rbind(scattrMap, tmp2)
  }
}
scattrMap <- scattrMap[which(!is.na(scattrMap$X)),]
scattrMap <- scattrMap[which(!is.na(scattrMap$Y)),]
Mn <- min(scattrMap$X)
scattrMap$X <- scattrMap$X - Mn
scattrMap$Y <- scattrMap$Y - Mn
Mx <- max(scattrMap$X)
scattrMap$X <- scattrMap$X / Mx
scattrMap$Y <- scattrMap$Y / Mx
scattrMap$`Sample 1` <- factor(scattrMap$`Sample 1`, levels = smpls)
scattrMap$`Sample 2` <- factor(scattrMap$`Sample 2`, levels = smpls)
corMap$`Sample 1` <- factor(corMap$`Sample 1`, levels = smpls)
corMap$`Sample 2` <- factor(corMap$`Sample 2`, levels = smpls)
ttl <- "Samples Pearson correlation map"
plot <- ggplot(scattrMap) +
  geom_bin2d(aes(x = X, y = Y), bins = max(c(1, round(nrow(PG)/200)))) +
  scale_fill_viridis() +
  new_scale("fill") +
  geom_tile(data = corMap, aes(fill = `Pearson corr.`, x = 0.5, y = 0.5), width = 1, height = 1) +
  scale_fill_viridis(option = "B") +
  facet_grid(`Sample 2`~`Sample 1`) + ggtitle(ttl) +
  theme_minimal() + theme(axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          strip.text.y.right = element_text(angle = 0, size = 5),
                          strip.text.x.top = element_text(angle = 90, size = 5))
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 600)
ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 600)

rm(list = ls()[which(!ls() %in% .obj)])
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
ReportCalls$Calls <- append(ReportCalls$Calls,
                            "body_add_fpar(Report, fpar(ftext(\"Calculate average intensities and ratios and performing statistical tests.\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
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
w <- which(sapply(Coefficients, function(x) { length(unique(FactorsLevels[[x]])) < nrow(Exp.map) }))
Coefficients <- Coefficients[w]
w <- which(sapply(Coefficients, function(x) { length(unique(Exp.map[[x]])) > 1 }))
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

### Visualize P-values
#
# Scatter plots:
pkol2 <- gsub(" - $", "", pvalue.col)
temp <- lapply(VPAL$values, function(vpal) { #vpal <- VPAL$values[1]
  j <- paste0(pvalue.col, vpal)
  w <- which(j %in% colnames(PG))
  j <- j[w]
  if (length(j)) {
    x <- set_colnames(PG[, j, drop = FALSE], pkol2[w])
    x$Group <- cleanNms(vpal)
  } else { x <- NA }
  return(x)
})
temp <- temp[which(sapply(temp, function(x) { "data.frame" %in% class(x) }))]
temp <- plyr::rbind.fill(temp)
kol <- colnames(temp)
kol <- kol[which(kol != "Group")]
temp <- temp[which(rowSums(temp[, kol]) > 0), ]
Comb <- gtools::combinations(length(kol), 2, kol)
pvalDir <- paste0(wd, "/Workflow control/Protein groups/P-values")
if (!dir.exists(pvalDir)) { dir.create(pvalDir, recursive = TRUE) }
dirlist <- unique(c(dirlist, pvalDir))
clusterExport(parClust, list("Comb", "temp", "pvalDir", "pvalue.col", "plotEval"), envir = environment())
temp2 <- parLapply(parClust, 1:nrow(Comb), function(i) { #i <- 1
  X <- Comb[i, 1]
  Y <- Comb[i, 2]
  X2 <- gsub(" -log10\\(Pvalue\\)$", "", X)
  Y2 <- gsub(" -log10\\(Pvalue\\)$", "", Y)
  x <- temp[, c(X, Y, "Group")]
  colnames(x) <- c("X", "Y", "Group")
  x$"P-value, X axis" <- proteoCraft::gsub_Rep(" -log10\\(Pvalue\\)$", "", Comb[i, 1])
  x$"P-value, Y axis" <- proteoCraft::gsub_Rep(" -log10\\(Pvalue\\)$", "", Comb[i, 2])
  ttl1 <- paste0("P-values scatter plot - ", X2, " VS ", Y2)
  ttl1a <- paste0(ttl1, " (-log10)")
  Mx <- max(proteoCraft::is.all.good(c(x$X, x$Y)))
  uX <- unique(x$`P-value, X axis`)
  uY <- unique(x$`P-value, Y axis`)
  plot1 <- ggplot2::ggplot(x) +
    scattermore::geom_scattermore(ggplot2::aes(x = X, y = Y, colour = Group),
                                  pixels = c(1024, 1024), pointsize = 3.2) +
    viridis::scale_color_viridis(begin = 0.4, end = 0.8, discrete = TRUE, option = "F") +
    #ggplot2::facet_grid(`P-value, Y axis`~Group+`P-value, X axis`) + # This was for when we saved all as one plot outside the parLapply
    ggplot2::facet_grid(Group~`P-value, Y axis`) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::coord_fixed() + ggplot2::ggtitle(ttl1a) + ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) +
    ggplot2::xlim(0, Mx) + ggplot2::ylim(0, Mx) +
    ggplot2::xlab(uX) + ggplot2::ylab(uY) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  #proteoCraft::poplot(plot1, 12, 20)
  Img1 <- paste0(pvalDir, "/", gsub(":", " - ", ttl1))
  # This was for when we saved all as one plot outside the parLapply:
  #h1 <- (length(pvalue.col)-0.8)*2
  #w1 <- ((length(pvalue.col)-1)*length(unique(x$Group))+1)*2
  w1 <- 5
  h1 <- w1*(length(unique(x$Group)) + 0.5)/3
  ggplot2::ggsave(paste0(Img1, ".jpeg"), plot1, dpi = 150, width = w1,# height = h1,
                  units = "in")
  ggplot2::ggsave(paste0(Img1, ".pdf"), plot1, dpi = 150, width = w1,# height = h1,
                  units = "in")
  plot1 <- plotEval(plot1)
  #system(paste0("open \"", Img1, ".jpeg\""))
  return(list(Plot = plot1, Title = ttl1))
})
for (nm in 1:length(temp2)) {
  ReportCalls <- AddPlot2Report(Plot = temp2[[i]]$Plot, Title = temp2[[i]]$Title, Dir = pvalDir)
}
#
# P-values histogram:
nbin <- 20
bd <- (0:nbin)/nbin
w <- which(sapply(pvalue.col, function(type) { #type <- pvalue.col[1]
  length(grep(topattern(type), colnames(PG))) > 0
}))
temp <- lapply(pvalue.col[w], function(type) { #type <- pvalue.col[w][1]
  kol <- grep(topattern(type), colnames(PG), value = TRUE)
  if (!length(kol)) { stop(type) }
  temp <- PG[, kol, drop = FALSE]
  colnames(temp) <- cleanNms(gsub(topattern(type), "", colnames(temp)))
  temp <- reshape::melt(temp, measure.vars = colnames(temp))
  temp$value <- 10^(-temp$value)
  temp <- temp[which(is.all.good(temp$value, 2)),]
  temp$Bin <- sapply(temp$value, function(x) { min(which(bd >= x))-1 })
  res <- aggregate(temp$Bin, list(temp$variable, temp$Bin), length)
  colnames(res) <- c("Group", "Bin", "Count")
  grps <- unique(res$Group)
  res$Frequency <- NA
  for (grp in grps) {
    w <- which(res$Group == grp)
    res$Frequency[w] <- res$Count[w]/sum(res$Count[w])
  }
  res$"P-value type" <- type
  res[, RSA$names] <- Exp.map[match(res$Variable, cleanNms(Exp.map[[VPAL$column]])),
                              RSA$names]
  res$Low <- 0
  return(res)
})
temp <- plyr::rbind.fill(temp)
temp$"P-value type" <- gsub_Rep(" -log10\\(Pvalue\\) - $", "", temp$"P-value type")
ttl2 <- "P-values histogram"
plot2 <- ggplot(temp) +
  geom_rect(position = "identity",
            aes(xmin = (Bin-1)/nbin, ymin = Low, xmax = Bin/nbin, ymax = Frequency, fill = `P-value type`),
            colour = "black", alpha = 0.5) + ylim(0, 1) +
  scale_fill_viridis(begin = 0.4, end = 0.8, discrete = TRUE, option = "G") +
  coord_fixed(ratio = 0.75) + facet_grid(`P-value type`~Group) + ggtitle(ttl2) + theme_bw() +
  theme(strip.text.y = element_text(angle=0))
#poplot(plot2, 12, 20)
Img2 <- paste0(pvalDir, "/", ttl2)
w2 <- ((length(VPAL$values)+1)*1.25)*2
h2 <- ((length(pvalue.col)+0.2)*1.25)*2
ggsave(paste0(Img2, ".jpeg"), plot2, dpi = 300, width = w2, height = h2, units = "in")
ggsave(paste0(Img2, ".pdf"), plot2, dpi = 300, width = w2, height = h2, units = "in")
#system(paste0("open \"", Img2, ".jpeg\""))
ReportCalls <- AddPlot2Report(Title = ttl2, Dir = pvalDir)
#
# Which type of P-values do we want to use?
Imgs1 <- list.files(pvalDir, "^P-values scatter plot - .*\\.jpeg$", full.names = TRUE)
Img2 <- gsub("(\\.jpeg)+$", ".jpeg", paste0(Img2, ".jpeg"))
Imgs1 <- Imgs1[which(Imgs1 != Img2)]
IMGS <- c(Img2, Imgs1)
msg <- "Confirm which type of P-values to use for t-test volcano plots"
pvalue.use %<o% (names(pvalue.col) == Param$P.values.type)
appNm <- paste0(dtstNm, " - t-test P-values")
Imgs1Nms <- gsub(" t-test|'s", "", gsub("Moderated", "Mod.", gsub("Permutations", "Perm.", gsub(".*/P-values scatter plot - |\\.jpeg", "", Imgs1))))
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  titlePanel(tag("u", "Type of t-test P-values"),
             appNm),
  br(),
  #
  fluidRow(column(2, selectInput("PVal", msg, names(pvalue.col), Param$P.values.type)),
           column(2, actionButton("saveBtn", "Save"))),
  br(),
  br(),
  fluidRow(column(6,
                  withSpinner(imageOutput("Img1", inline = TRUE)),
                  selectInput("XY", "Comparison", Imgs1Nms, Imgs1Nms[1])),
           column(6, withSpinner(imageOutput("Img2", inline = TRUE)))),
  br(),
  br()
)
#h0 <- paste0(round(screenRes$height*0.35), "px")
source(parSrc, local = FALSE)
IMGsDims <- as.data.frame(t(parSapply(parClust, IMGS, function(x) { #x <- IMGs[1]
  a <- jpeg::readJPEG(x)
  setNames(dim(a)[1:2], c("height", "width"))
})))
IMGsDims$height <- screenRes$width*IMGsDims$height/max(IMGsDims$height)*0.7
IMGsDims$width <- screenRes$width*IMGsDims$width/max(IMGsDims$width)*0.7
IMGsDims$height[1] <- IMGsDims$height[1]*0.8
IMGsDims$width[1] <- IMGsDims$width[1]*0.8
IMGsDims$width <- IMGsDims$width*0.8 # Correction for skew
fct <- 3
server <- function(input, output, session) {
  myIMG1 <- reactiveVal(1)
  updtIMG1 <- function(reactive = TRUE) {
    if (!reactive) {
      renderImage({
        list(src = Imgs1[1],
             height = IMGsDims$height[1+1]*fct,
             width = IMGsDims$width[1+1]*fct)
      }, deleteFile = FALSE)
    } else {
      renderImage({
        list(src = Imgs1[myIMG1()],
             height = IMGsDims$height[myIMG1()+1]*fct,
             width = IMGsDims$width[myIMG1()+1]*fct)
      }, deleteFile = FALSE)
    }
  }
  output$Img1 <- updtIMG1(FALSE)
  output$Img2 <- renderImage({
    list(src = IMGS[1], height = IMGsDims$height[1], width = IMGsDims$width[1])
  }, deleteFile = FALSE)
  observeEvent(input$PVal, {
    assign("pvalue.use", names(pvalue.col) == input$PVal, envir = .GlobalEnv)
  }, ignoreInit = FALSE)
  observeEvent(input$XY, {
    myIMG1(match(input$XY, Imgs1Nms))
    output$Img1 <- updtIMG1()
  }, ignoreInit = FALSE)
  observeEvent(input$saveBtn, { stopApp() })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
eval(parse(text = runApp), envir = .GlobalEnv)

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
  #View(tst3[which(sapply(tst3$x, length) > 1),])
  #View(tst3[which(sapply(tst3$x, length) == 0),])
  #
  # Also peptides (minor approximation: use first protein group)
  pep[, kol] <- PG[match(as.integer(gsub(";.*", "", pep$`Protein group ID`)), PG$id), kol]
  #
  PG$Ontology <- NULL # Temporary fix for now, this column is broken
  #
  stopCluster(parClust)
  source(parSrc, local = FALSE)
  Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_prepare.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)

#### Code chunk - ROC analysis
warning("Parameter \"ROC.GO.terms\" is currently not supported in the Parameters shiny App! Add it!")
if (("ROC.GO.terms" %in% colnames(Param))&&(!as.character(Param$ROC.GO.terms) %in% c("", "NA", NA, " "))) {
  ROC_GO.terms %<o% unique(unlist(strsplit(Param$ROC.GO.terms, ";")))
  ROC_GO.terms <- ROC_GO.terms[which(ROC_GO.terms %in% GO_terms$ID)]
  if (length(ROC_GO.terms)) {
    msg <- "ROC analysis"
    ReportCalls <- AddMsg2Report()
    ROC_GO.terms <- unique(unlist(c(ROC_GO.terms, GO_terms$Offspring[match(ROC_GO.terms, GO_terms$ID)])))
    PG$"True Positive" <- FALSE
    PG$"True Positive"[grsep2(ROC_GO.terms, PG$`GO-ID`)] <- TRUE
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
          PG$Predictor <- 10^(-PG[[pkol]])
          rocobj <- roc(PG[w2, ], "True Positive", "Predictor")
          ttl <- paste0("ROC analysis - ", grp2)
          plot <- geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dotted") +
            ggroc(rocobj) + theme_bw() + ggtitle(ttl)
          poplot(plot)
          ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300)
          ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)
          ReportCalls <- AddPlot2Report(Title = gsub(": ?", " - ", ttl))
          PG$Predictor <- NULL
        }
      }
    }
  }
}

#### Code chunk - Estimate P-value significance for a set of accepted FDRs
## NB: For graphical reasons (volcano plots), there is only support for 4 different FDR values. This should suffice anyway.
a <- sapply(strsplit(Param$Plot.metrics, ";"), function(x) {strsplit(x, ":")})
a[[2]][2] <- gsub("\\.$", "", pvalue.col[which(pvalue.use)])
Param$Plot.metrics <- paste(sapply(a, function(x) {paste(x, collapse = ":")}), collapse = ";")
FDR.thresholds %<o% c()

A <- VPAL$values
test <- sapply(A, function(x) { #x <- A[6]
  x <- paste0(pvalue.col[which(pvalue.use)], x)
  r <- x %in% colnames(PG)
  if (r) { r <- length(is.all.good(as.numeric(PG[[x]]))) > 0 }
  return(r)
})
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
                           Size = "Av. log10 abundance", Size.max = 2,
                           plotly = create_plotly, plotly_local = create_plotly_local,
                           plotly_labels = PrLabKol,
                           cl = parClust,
                           SAM = useSAM, curved_Thresh = SAM_thresh))
if (!class(tempVP) %in% c("try-error", "character")) {
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
  # Legacy code for web-hosted plotly plots:
  if ((create_plotly)&&(!create_plotly_local)) { plot_ly$"t-tests" <- tempVP$"Plotly plots" }
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
                            Size = "Av. log10 abundance", Size.max = 2,
                            plotly = create_plotly, plotly_local = create_plotly_local,
                            cl = parClust)
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
    for (i in 1:length(grp1)) { #i <- 1
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
    if (sum(sapply(RG$names, function(x) {! x %in% RSA$names })) > 0) {
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
  for (i in 1:length(g)) {
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
  tst <- apply(Exp.map[, RRG$names, drop = FALSE], 1, paste, collapse = "___")
  tmp <- apply(g2[,RRG$names, drop = FALSE], 1, paste, collapse = "___")
  g2$Ref <- sapply(tmp, function(x) { y <- unique(tst[which((Exp.map$Reference)&(tst == x))]) })
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

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
parLapply(parClust, 1:N.clust, function(x) {
  rm(list = ls())
  gc()
})
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - F-test
#Param <- Param.load()
F.test %<o% FALSE
if (("F.test" %in% colnames(Param))&&(Param$F.test)) {
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
  # stopifnot(sum(!sapply(tmp1, function(x) { length(x) %in% c(2:3)[Nested+1] })) == 0,
  #           sum(!sapply(tmp2, function(x) { length(x) == 2 })) == 0) # If this breaks, this will mean that these parameters are not built the way I remember them to be.
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
  Src <- paste0(libPath, "/extdata/R scripts/Sources/run_F_test.R")
  #rstudioapi::documentOpen(Src)
  tstFtst <- try(source(Src, local = FALSE), silent = TRUE)
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
    # Legacy code for web-hosted plotly plots:
    if ((create_plotly)&&(!create_plotly_local)) { plot_ly$"F-tests" <- F_volc$"Plotly plots" }
    # Also a posteriori F-test P- or Q- values histogram:
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
      ggsave(paste0(dir, "/", ttla, ".jpeg"), plot, dpi = 300)
      ggsave(paste0(dir, "/", ttla, ".pdf"), plot, dpi = 300)
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
      for (i in 1:length(g)) { #i <- 1
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
    rm(list = ls()[which(!ls() %in% .obj)])
    Script <- readLines(ScriptPath)
    gc()
    parLapply(parClust, 1:N.clust, function(x) {
      rm(list = ls())
      gc()
    })
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
DatAnalysisTxt <- paste0(DatAnalysisTxt, " Average log10 expression values were tested for significance using a ",
                         c("two", "one")[match(AltHyp, c("two.sided", "greater", "lower"))], "-sided ",
                         tmpPVal, " per samples group",
                         c("", " and a global F-test")[F.test+1],
                         " (limma). Significance thresholds were calculated using the Benjamini-Hochberg procedure for False Discovery Rate (FDR) values of ", tmp,
                         ". For all tests, regulated protein groups were defined as those with a significant P-value and a",
                         c("n absolute", "")[IsPullDown+1],
                         " log2 ratio greater than ",
                         c(paste0(Param$Ratios.Contamination.Rates*100, "% of ",
                                  c(paste0("control-to", c("-average", "")[Nested + 1],
                                           "-control"),
                                    "intra-sample groups")[match(RefRat_Mode, c("1", "2"))],
                                  " ratios"),
                           Param$Ratios.Contamination.Rates)[match(Param$Ratios.Thresholds,
                                                                   threshOpt)], ".")

#### Code chunk - SAINTexpress
if (saintExprs) { saintExprs <- "Target" %in% colnames(Exp.map) }
if (saintExprs) {
  SaintRoot <- "C:/SAINTexpress"
  SaintDir <- paste0(SaintRoot, "/SAINTexpress_v3.6.3__2018-03-09")
  if (!dir.exists(SaintDir)) { dir.create(SaintDir, recursive = TRUE) }
  SaintEx <- paste0(SaintDir, "/Precompiled_binaries/Windows64/SAINTexpress-int.exe")
  if (!file.exists(SaintEx)) {
    url <- "https://download.sourceforge.net/saint-apms/SAINTexpress_v3.6.3__2018-03-09.tar.gz"
    packs <- c("curl")
    for (pack in packs) {
      if (!suppressMessages(require(pack, character.only = TRUE))) { install.packages(pack, update = FALSE) }
      require(pack, character.only = TRUE)
    }
    cran_req <- unique(c(cran_req, packs))
    destFl <- paste0(SaintRoot, "/SAINTexpress_v3.6.3__2018-03-09.tar.gz")
    kount <- 0
    while ((!kount)||((kount < 5)&&("try-error" %in% class(tst)))) {
      tst <- try(download.file(url, destFl), silent = TRUE)
      kount <- kount+1
    }
    if ("try-error" %in% class(tst)) { saintExprs <- FALSE } else {
      gunzip(destFl)
      utils::untar(gsub("\\.gz$", "", destFl), exdir = SaintRoot)
      unlink(destFl)
      unlink(gsub("\\.gz$", "", destFl))
      saintExprs <- file.exists(SaintEx)
    }
  }
}
if (saintExprs) {
  #
  msg <- "Running SAINTexpress analysis..."
  ReportCalls <- AddMsg2Report(Space = FALSE)
  #
  subDr <- "Reg. analysis/SAINTexpress"
  saintDir <- gsub("/+", "/", paste0(wd, "/", subDr))
  if (!dir.exists(saintDir)) { dir.create(saintDir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, saintDir))
  fcRt <- "log2(FC) - "
  fdrRt <- "-log10(BFDR) - "
  #
  # Copy the SAINTexpress manual
  #url <- "https://raw.githubusercontent.com/bornea/APOSTL/master/wk_images/SAINTexpress-manual.pdf"
  #destfile <- paste0(saintDir, "/SAINT express manual.pdf")
  #tst <- try(download.file(url, destfile, "curl"), silent = TRUE)
  #if ("try-error" %in% class(tst)) { try(download.file(url, destfile, "wget"), silent = TRUE) }
  # Description:
  # http://saint-apms.sourceforge.net/Main.html
  fl <- system.file("extdata", "SAINTexpress-manual.pdf", package = "proteoCraft")
  try(file.copy(fl, saintDir), silent = TRUE)
  #
  Indiq <- c("T", "C")
  if (Mirror.Ratios) { Indiq <- rev(Indiq) } # Deprecate me please!!!
  # Prey and GO tables are not Target specific
  mtch <- listMelt(strsplit(PG$"Protein IDs", ";"), PG$id, c("Protein", "PG id"))
  mtch <- set_colnames(aggregate(mtch$`PG id`, list(mtch$Protein), unique), c("Protein", "PG ids"))
  Interact <- Prey <- data.frame(Protein = mtch$Protein)
  Prey[, c("Sequence", "Gene")] <- db[match(gsub("^CON__", "", Prey$Protein), gsub("^CON__", "", db$`Protein ID`)), c("Sequence", "Gene")]
  Prey$Length <- nchar(Prey$Sequence)
  w <- which(is.na(Prey$Length))
  if (length(w)) {
    Prey$Length[w] <- median(Prey$Length, na.rm = TRUE)
  }
  Prey$Sequence <- NULL
  Prey <- Prey[, c("Protein", "Length", "Gene")]
  w <- which((is.na(Prey$Gene))|(Prey$Gene == ""))
  l <- length(w)
  if (l) { Prey$Gene[w] <- paste0("Dummy_Gene", 1:l) }
  if (Annotate) {
    GO <- listMelt(strsplit(PG$`GO-ID`, ";"), PG$id)
    GO$Proteins <- PG$`Protein IDs`[match(GO$L1, PG$id)]
    GO <- listMelt(strsplit(GO$Proteins, ";"), GO$value)
  }
  Bait <- data.frame(IP_name = gsub("\\.", "", cleanNms(Exp.map$Ref.Sample.Aggregate, rep = "")),
                     Bait = Exp.map$Target,
                     Indicator = Indiq[Exp.map$Reference+1])
  #Bait <- Bait[which((!is.na(Bait$Bait))&(Bait$Bait %in% Prey$Protein)),]
  Bait$Bait[which((is.na(Bait$Bait))|(!Bait$Bait %in% Prey$Protein))] <- "CONTROL"
  # if (("CONTROL" %in% Bait$Bait)&&(!"CONTROL" %in% Prey$Protein)) {
  #   # Add a dummy CONTROL protein if any co-IPs do not have a bait (typically IP- isotype controls)
  #   Prey <- rbind(Prey,
  #                 data.frame(Protein = "CONTROL",
  #                            Length = 10000, # Cannot be NA
  #                            Gene = "dummyGene"))
  # }
  kol <- paste0(Prot.Expr.Root, Exp.map$Ref.Sample.Aggregate)
  klnms <- gsub("\\.", "", cleanNms(Exp.map$Ref.Sample.Aggregate, rep = ""))
  tmp <- PG[, c("id", kol)]
  clusterExport(parClust, list("tmp", "kol"), envir = environment())
  Interact[, klnms] <- as.data.frame(t(parSapply(parClust, mtch$`PG ids`, function(x) {
    m <- match(unlist(x), tmp$id)
    if (!length(m)) { stop(x) }
    x <- tmp[m, kol, drop = FALSE]
    return(as.numeric(apply(x, 2, function(x) { mean(proteoCraft::is.all.good(x)) })))
  })))
  Interact <- reshape::melt(Interact, id.vars = "Protein")
  colnames(Interact) <- c("Protein", "IP_name", "Intensity")
  Interact$IP_name <- as.character(Interact$IP_name)
  Interact$Intensity <- 10^Interact$Intensity
  Interact <- Interact[which(is.all.good(Interact$Intensity, 2)),]
  Interact <- Interact[which(Interact$Intensity > 0),]
  #Interact <- Interact[which(Interact$IP_name %in% Bait$IP_name),]
  Interact$Bait <- Bait$Bait[match(Interact$IP_name, Bait$IP_name)]
  Interact <- Interact[, c("IP_name", "Bait", "Protein", "Intensity")]
  #data.table::fwrite(Bait, baitFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
  #data.table::fwrite(Prey, preyFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
  #data.table::fwrite(Interact, interFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
  Grps <- unique(Exp.map[which(!Exp.map$Reference), VPAL$column])
  Parma <- Param
  Parma$Plot.labels <- "Common Name"
  Parma$Plot.threshold.metrics <- Parma$Plot.threshold.values <- Parma$Plot.threshold.tests <- Parma$Plot.threshold.colours <- ""
  saveRDS(Interact, paste0(saintDir, "/Interact.RDS"))
  saveRDS(Prey, paste0(saintDir, "/Prey.RDS"))
  if (Annotate) { saveRDS(GO, paste0(saintDir, "/GO.RDS")) } 
  kol1 <- c("AvgP", "MaxP", "TopoAvgP", "TopoMaxP", "SaintScore")
  kol2 <- c("OddsScore", "BFDR", "boosted_by")
  clusterExport(parClust,
                list("Exp.map", "Exp", "VPAL", "RG", "Bait", "Annotate", "SaintEx", "saintDir", "wd", "fcRt", "kol1", "kol2"),
                envir = environment())
  saintst <- setNames(parLapply(parClust, Grps, function(grp) { #grp <- Grps[1]
    grp2 <- proteoCraft::cleanNms(grp)
    dr <- paste0(saintDir, "/", grp2)
    if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
    setwd(dr)
    ratgrp <- unique(Exp.map[match(grp, Exp.map[[VPAL$column]]), RG$column])
    m <- Exp.map[which(Exp.map[[RG$column]] == ratgrp),]
    m <- m[which((m[[VPAL$column]] == grp)|(m$Reference)),]
    Bait2 <- Bait[which(Bait$IP_name %in% gsub("\\.", "", proteoCraft::cleanNms(m$Ref.Sample.Aggregate, rep = ""))),]
    Interact <- readRDS(paste0(saintDir, "/Interact.RDS"))
    Prey <- readRDS(paste0(saintDir, "/Prey.RDS"))
    Interact2 <- Interact[which(Interact$IP_name %in% Bait2$IP_name),]
    Prey2 <- Prey[which(Prey$Protein %in% c(Bait2$Bait, Interact2$Protein)),]
    baitFl <- paste0(dr, "/tempBait.txt")
    preyFl <- paste0(dr, "/tempPrey.txt")
    interFl <- paste0(dr, "/tempInteract.txt")
    lstFl <- paste0(dr, "/list.txt") # Important: it saves the list in the local directory!!!
    data.table::fwrite(Bait2, baitFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
    data.table::fwrite(Prey2, preyFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
    data.table::fwrite(Interact2, interFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
    if (Annotate) {
      goFl <- paste0(dr, "/tempGO.txt")
      GO2 <- readRDS(paste0(saintDir, "/GO.RDS"))
      GO2 <- GO2[which(GO2$value %in% Prey$Protein),]
      GO2 <- aggregate(GO2$value, list(GO2$L1), function(x) { paste(unique(x), collapse = " ") })
      data.table::fwrite(GO2, goFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
    }
    fls <- c(SaintEx, c(interFl, preyFl, baitFl, goFl)[1:(3+Annotate)])
    w <- which(!file.exists(fls))
    stopifnot(length(w) == 0)
    #fls[w]
    fls2 <- paste0("\"", fls, "\"")
    # cmd <- paste0(c(fls2[1],
    #                 paste0("-L", sum(Exp.map$Reference)),
    #                 fls2[2:length(fls2)]), collapse = " ")
    if (file.exists(lstFl)) { unlink(lstFl) }
    cmd <- paste0(c(fls2[1],
                    paste0("-L", sum(m$Reference)),
                    fls2[2:length(fls2)]), collapse = " ")
    #cat(cmd)
    #writeClipboard(cmd)
    #cat("   ", grp2, "\n")
    system(cmd)
    #cat("    -> done\n")
    rs <- list(Outcome = file.exists(lstFl))
    if (rs$Outcome) {
      # Read and process results
      tmpDF <- data.table::fread(lstFl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
      tmpDF <- tmpDF[order(tmpDF$SaintScore, decreasing = TRUE),]
      for (k in kol1) {
        w <- which(tmpDF[[k]] == ".")
        if (length(w)) { tmpDF[w, k] <- NA }
        tmpDF[[k]] <- as.numeric(tmpDF[[k]])
      }
      tmpDF[, paste0(kol1, " - ", grp)] <- tmpDF[, kol1]
      for (k in kol2) {
        tmpDF[[paste0(k, " - ", grp)]] <- tmpDF[[k]]
      }
      tmpDF[[paste0(fcRt, grp)]] <- log2(tmpDF$FoldChange)
      tmpDF <- tmpDF[, which(!colnames(tmpDF) %in% c(kol1, kol2))]
      tmpDF$"Av. log10 abundance" <- log10(tmpDF$AvgIntensity)
      rs$Table <- tmpDF
    } else { warning(paste0("SAINTexpress analysis failed for group ", grp2)) }
    setwd(wd) # Important to release the subfolder
    return(rs)
  }), Grps)
  setwd(wd)
  saintst <- saintst[which(sapply(saintst, function(x) { x$Outcome} ))]
  l <- length(saintst)
  if (l) {
    nms <- names(saintst)
    saintst <- setNames(lapply(saintst, function(x) { x$Table }), nms)
    msg <- "   -> Reading results\n"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    allSAINTs %<o% data.frame(Protein = unique(unlist(lapply(saintst, function(x) { x$Prey }))))
    for (i in 1:l) {
      kol <- paste0(c("log2(FC)", kol1, kol2), " - ", nms[i])
      tmp <- saintst[[nms[i]]][, c("Prey", kol)]
      allSAINTs[, kol] <- NA
      w <- which(allSAINTs$Protein %in% tmp$Prey)
      allSAINTs[w, kol] <- tmp[match(allSAINTs$Protein[w], tmp$Prey), kol]
    }
    m <- match(allSAINTs$Protein, db$`Protein ID`)
    kol <- c("Potential contaminant", "Common Name", "Gene")
    allSAINTs[, kol] <- db[m, kol]
    tmp <- listMelt(strsplit(PG$`Protein IDs`, ";"), PG$id)
    tmp <- tmp[which(tmp$value %in% allSAINTs$Protein),]
    tmp$"Av. log10 abundance" <- PG$`Av. log10 abundance`[match(as.numeric(tmp$L1), PG$id)] # Not exactly, but good enough for now (this is using PG-level data for a Protein-level table)
    allSAINTs$"Av. log10 abundance" <- tmp$"Av. log10 abundance"[match(allSAINTs$Protein, tmp$value)]
    # Since here we will be plotting FDR, not P-values,
    # we can use our parameter BH FDR thresholds as they are
    ArbThr <- data.frame(yintercept = -log10(BH.FDR),
                         slope = 0,
                         xintercept = NA,
                         colour = colorRampPalette(c("orange", "red"))(length(BH.FDR)),
                         label = paste0(BH.FDR*100, "% FDR"))
    #
    k1 <- grep("^BFDR - ", colnames(allSAINTs), value = TRUE)
    k2 <- gsub("^BFDR - ", fdrRt, k1)
    tmp1 <- as.matrix(allSAINTs[, k1])
    tmp2 <- -log10(tmp1)
    w <- which(tmp1 == 0, arr.ind = TRUE)
    nr <- nrow(w)
    if (nr) {
      mx <- max(is.all.good(as.numeric(tmp2)))
      ArbThr <- rbind(ArbThr,
                      data.frame(yintercept = mx+0.1,
                                 slope = 0,
                                 xintercept = NA,
                                 colour = "#F000FF",
                                 label = "^ Infinite value (BFDR = 0)"))
      rp <- mx+0.2
      msg <- paste0("      Replacing ", nr, " infinite -log10(BFDR) values with a high but finite value of ", rp, " for the purpose of plotting them!")
      ReportCalls <- AddMsg2Report(Space = FALSE)
      tmp2[w] <- rp
    }
    allSAINTs[, k2] <- tmp2
    #
    # Match to PGs
    tmp1 <- listMelt(setNames(strsplit(PG$`Leading protein IDs`, ";"), PG$id), ColNames = c("Prot", "PG"))
    tmp2 <- listMelt(setNames(strsplit(PG$`Protein IDs`, ";"), PG$id), ColNames = c("Prot", "PG"))
    allSAINTs$PG_id <- tmp1$PG[match(allSAINTs$Protein, tmp1$Prot)]
    w <- which(is.na(allSAINTs$PG_id))
    if (length(w)) {
      allSAINTs$PG_id[w] <- tmp2$PG[match(allSAINTs$Protein[w], tmp2$Prot)]
    }
    allSAINTs$PG_id <- as.integer(allSAINTs$PG_id)
    #
    data.table::fwrite(allSAINTs, paste0(saintDir, "/SAINTexpress interactions list.tsv"),
                       row.names = FALSE, sep = "\t", na = "NA")
    #
    # Volcano plot
    labKol <- setNames(c("Common Name", "Protein", "Gene"),
                       c("Protein name", "Protein ID", "Gene"))
    tempVPip <- Volcano.plot(Prot = allSAINTs, Proteins.col = "Protein",
                             mode = "custom",
                             experiments.map = Exp.map,
                             X.root = fcRt,
                             Y.root = fdrRt,
                             aggregate.map = Aggregate.map,
                             aggregate.name = VPAL$aggregate,
                             aggregate.list = Aggregate.list, parameters = Parma,
                             save = c("jpeg", "pdf"), labels = "thresholds",
                             Ref.Ratio.values = Ref.Ratios,
                             Ref.Ratio.method = paste0("obs", RefRat_Mode),
                             ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                             arbitrary.lines = ArbThr,
                             proteins = prot.list, proteins_split = protsplit,
                             return = TRUE, return.plot = TRUE,
                             title = "SAINTexpress volcano plot ",
                             subfolder = subDr,
                             subfolderpertype = FALSE, Symmetrical = TwoSided,
                             Alpha = "Av. log10 abundance",
                             Size = "Av. log10 abundance", Size.max = 2,
                             plotly = create_plotly, plotly_local = create_plotly_local,
                             plotly_labels = labKol,
                             cl = parClust)
    #
    VP_list <- tempVPip
    insrt <- ""
    Src <- paste0(libPath, "/extdata/R scripts/Sources/thresholds_Excel.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    # Folder cleanup
    unlink(paste0(saintDir, "/", c("Interact", "Prey"), ".RDS"))
    if (Annotate) { unlink(paste0(saintDir, "/GO.RDS")) }
    unlink(paste0(saintDir, "/", cleanNms(nms)))
    #
    # Regulated columns
    # Here we do not use the ones from Volcano.plot because the logic is different than usual
    regkol <- paste0("Regulated - ", nms)
    Reg_filters$"SAINTexpress" <- list()
    if ("con" %in% filter_types) {
      Reg_filters$"SAINTexpress"$"By condition" <- list()
    }
    allSAINTs[, regkol] <- ""
    for (nm in nms) { #nm <- nms[1]
      thresh <- tempVPip$Thresholds$Absolute[[nm]]
      up <- thresh$Value[match("up", thresh$Levels)]
      if (TwoSided) { down <- thresh$Value[match("down", thresh$Levels)] }
      k1 <- paste0("BFDR - ", nm)
      k2 <- paste0("log2(FC) - ", nm)
      k3 <- paste0("Regulated - ", nm)
      allSAINTs[which(allSAINTs[[k1]] <= max(BH.FDR)), k3] <- "too small FC" # base level
      for (f in rev(BH.FDR)) {
        allSAINTs[which((allSAINTs[[k1]] <= f)&(allSAINTs[[k2]] >= up)), k3] <- paste0("up, FDR = ", f*100,"%")
        if (TwoSided) {
          allSAINTs[which((allSAINTs[[k1]] <= f)&(allSAINTs[[k2]] <= down)), k3] <- paste0("down, FDR = ", f*100,"%")
        }
      }
      # Create SAINTexpress-based filters
      if ("con" %in% filter_types) {
        up <- grep("^up|^Specific", unique(unlist(allSAINTs[, regkol])), value = TRUE)
        down <- grep("^down|^Anti-specific", unique(unlist(allSAINTs[, regkol])), value = TRUE) # The "anti-specific" part will only become relevant if in future I disconnect symmetry and/or adding specific tags from pull-down experiments
        Reg_filters$"SAINTexpress"$"By condition"[[nm]] <- list(Columns = k3,
                                                                Filter_up = sort(which(allSAINTs[[k3]] %in% up)),
                                                                Filter_down = sort(which(allSAINTs[[k3]] %in% down)),
                                                                Filter = sort(which(allSAINTs[[k3]] %in% c(up, down))))
      }
    }
    # Mat-meth text
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " Data was also tested with the SAINTexpress algorithm - which directly outputs FDR values - using GO annotations as boosting information.")
    #
    
  }
  dlg_message("TO DO: add SAINTq for PSM-level analysis of DIA data!", "ok")
}
setwd(wd)

# Now let's create a table of regulated protein groups per test made:
g <- grep("^Regulated - ", colnames(PG), value = TRUE)
regPG_TTest <- data.frame(Test = gsub("^Regulated - ", "", g))
if (TwoSided) { dir <- c("up", "down") } else { dir <- "up" }
kolstms <- c("count", "PG IDs", "Leading Protein IDs", "Genes")
for (d in dir) { #d <- "up"
  for (i in 1:length(BH.FDR)) { #i <- 1
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
    for (i in 1:length(BH.FDR)) { #i <- 1
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
          temp[i, 2:(N+1)] <- sapply(temp[1, 2:(N+1)], function(y) { sum(x %in% flt[[y]]$Filter) })
        }
        nms <- cleanNms(names(flt))
        tst <- sapply(strsplit(nms, " - "), length)
        tst <- (min(tst) > 1)&(length(unique(tst)) == 1)
        if (tst) {
          tst <- as.data.frame(t(sapply(strsplit(nms, " - "), unlist)))
          l <- apply(tst, 2, function(x) { length(unique(x)) })
          tst <- tst[, which(l > 1), drop = FALSE]
          nms <- apply(tst, 1, paste, collapse = " - ")
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
dataType <- "PG"
Src <- paste0(libPath, "/extdata/R scripts/Sources/GSEA.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Heatmaps with clustering at samples and protein groups level, highlighting proteins of interest
dir <- paste0(wd, "/Clustering")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
lblKol <- RSA$names
if ((length(Exp) == 1)&&(length(lblKol) > 1))  { lblKol <- lblKol[which(lblKol != "Experiment")] }
wSmpls <- which(paste0(Prot.Expr.Root, Exp.map$Ref.Sample.Aggregate) %in% colnames(PG))
xMap <- Exp.map[wSmpls,]
Samples <- apply(xMap[, lblKol, drop = FALSE], 1, paste, collapse = " ")
if (("GO.enrichment.Ref.Aggr" %in% colnames(Param))&&(Param$GO.enrichment.Ref.Aggr %in% Aggregate.map$Aggregate.Name)) {
  ClustGrp <- Param$GO.enrichment.Ref.Aggr
} else {
  ClustGrp <- Param$GO.enrichment.Ref.Aggr <- "Exp"
}
if ((length(Exp) == 1)&&(nchar(ClustGrp) %% 3 > 0)) { ClustGrp <- Param_filter(ClustGrp, "Exp") }
val <- Aggregate.list[[ClustGrp]]
nms <- unlist(Aggregate.map$Characteristics[which(Aggregate.map$Aggregate.Name == ClustGrp)])
if (length(nms) == 1) { kol <- nms } else { kol = ClustGrp }
ClustGrp <- list(aggregate = ClustGrp, values = val, names = nms, column = kol)
I <- list(Global = Samples[which(xMap[[ClustGrp$column]] %in% ClustGrp$values)])
if (length(ClustGrp$values) > 1) { for (i in ClustGrp$values) {
  iNm <- cleanNms(i)
  I[[iNm]] <- Samples[which(xMap[[ClustGrp$column]] %in% i)]
} }
I <- I[which(sapply(I, length) > 1)]
Heatmaps %<o% list()
NHClust %<o% list()
NVClust %<o% list()
KlustKols %<o% c()
ImputeKlust %<o% TRUE
MaxHClust %<o% min(c(floor(nrow(PG)/2), 100)) # We want at most 20 clusters
MaxVClust %<o% length(VPAL$values)
VClustScl <- setNames(1:MaxVClust, paste0("Cluster", 1:MaxVClust))
HClustScl <- setNames(1:MaxHClust, paste0("Cluster", 1:MaxHClust))
VClustScl <- setNames(rainbow(MaxVClust), VClustScl)
HClustScl <- setNames(rainbow(MaxHClust), HClustScl)
VClusters %<o% list()
# Different options for which proteins to use for vertical clustering (samples level)
VClustUse <- "All" # Use all
# VClustUse <- 500 # Max number of proteins to use (starting from ones with highest CV)
# VClustUse <- 25% # Percentage of proteins to use (starting from ones with highest CV)
# VClustUse <- "DEP" # Use only Differentially Expressed Proteins (from Reg_filters)
VClustUse <- toupper(VClustUse)
KlustRoot %<o% paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ") - ")
normTypes <- c("Norm. by row", "Z-scored")
plotLeatMaps %<o% list()
for (i in names(I)) { #i <- names(I)[1] #i <- names(I)[2]
  smpls <- I[[i]]
  smplsMtch <- match(smpls, Samples)
  tstlblKol <- sapply(lblKol, function(x) { length(unique(xMap[smplsMtch, x])) })
  smpls <- apply(xMap[match(smpls, Samples), lblKol[which(tstlblKol > 1)], drop = FALSE], 1, paste, collapse = " ")
  xMap$Samples <- smpls
  kol <- paste0(prtRfRoot, xMap$Ref.Sample.Aggregate[smplsMtch])
  msg <- paste0("Creating ", c("", paste0(i, " "))[(i == "Global")+1], "heatmap",
                c(paste0(" for ", i), "")[(i == "Global")+1], ".")
  ReportCalls <- AddMsg2Report(Space = FALSE)
  for (normType in normTypes) { #normType <- normTypes[1] #normType <- normTypes[2]
    normTypeInsrt <- c("", paste0(" (", normType, ")"))[match(normType, normTypes)]
    temp <- set_rownames(set_colnames(PG[, kol], smpls), PG$Label)
    Gr <- xMap[smplsMtch, RG$column]
    Gr <- setNames(match(Gr, unique(Gr)), Gr)
    tst <- apply(temp, 1, function(x) { length(is.all.good(x)) })
    filt <- rep(FALSE, nrow(temp))
    temp <- temp[which((tst > 0)&(PG$`Potential contaminant` != "+")),]
    if (ImputeKlust) {
      temp2 <- Data_Impute2(temp, Gr)
      imputed <- temp2$Positions_Imputed
      temp <- temp2$Imputed_data
      rownames(imputed) <- rownames(temp)
      colnames(imputed) <- colnames(temp)
    }
    w <- which(apply(temp, 1, function(x) { length(is.all.good(x)) }) == length(kol))
    filt[which(tst > 0)[w]] <- TRUE
    temp <- temp[w,]
    imputed <- imputed[w,]
    rwMns <- rowMeans(temp)
    if (normType == "rowNorm") {
      temp <- sweep(temp, 1, rwMns, "-")
    }
    if (normType == "Z-scored") {
      SDs <- apply(temp, 1, function(x) { sd(is.all.good(x)) })
      temp <- sweep(sweep(temp, 1, rwMns, "-"), 1, SDs, "/")
    }
    temp2 <- as.matrix(temp)
    temp3 <- temp2 <- temp2 + runif(length(temp2), min = 0, max = 1e-10) # Small error added to avoid duplicate rows where this breaks
    # Data is now normalized and either imputed or filtered, so ready for clustering
    # 1/ At samples level
    # We perform hierarchical clustering in all cases, because we want to see the dendrogram.
    # But the clusters displayed using colours may be generated using either k-means or hierarchical clustering (default).
    temp2 <- t(temp2)
    tst2 <- apply(temp2, 2, function(x) {
      #x <- is.all.good(x) # it's been imputed, there are no missing values
      sd(x)/mean(x)
    })
    temp2 <- temp2[, order(tst2, decreasing = TRUE)]
    if (VClustUse != "ALL") {
      if (!is.na(suppressWarnings(as.integer(round(as.numeric(VClustUse)))))) {
        VClustUse <- as.integer(round(as.numeric(VClustUse)))
        if (VClustUse > 0) { temp2 <- temp2[,1:VClustUse] }
      }
      if (grepl("^[0-9]+\\.?[0-9]*%$", VClustUse)) {
        VClustUse <- as.numeric(gsub("%$", "", VClustUse))/100
        if (VClustUse > 0) { temp2 <- temp2[,1:max(c(1, round(ncol(temp2)*VClustUse)))] }
      }
      if (VClustUse == "DEP") {
        PGflt <- sapply(names(Reg_filters), function(x) { #x <- names(Reg_filters)[1]
          x <- Reg_filters[[x]]$"By condition"
          sapply(names(x), function(y) { x[[y]]$Filter })
        })
        PGflt <- unique(unlist(PGflt))
        temp2 <- temp2[, which(colnames(temp2) %in% PG$Label[PGflt])]
      }
    }
    vcluster <- hclust(dist(temp2))
    vdendro <- as.dendrogram(vcluster)
    # Estimate ideal number of clusters... but ignore it!
    # In fact we know the number of sample groups, so would like to see 1 cluster per group.
    # How well groups and clusters overlap would tell us how well the clustering works, i.e. how different clusters are.
    #
    # Here we use kmeans, but findings apply to any method
    #cat("Estimating optimal number of samples-level clusters...\n")
    vnm <- paste0("Samples-level clusters number analysis - ", i)
    NVClust[[i]] <- NGr <- length(unique(Gr))
    tst <- cluster::clusGap(temp2, stats::kmeans, min(c(nrow(temp2), NGr+1)))
    tst2 <- as.data.frame(tst$Tab)  
    yScl <- max(tst2$gap)
    tst2 <- sapply(1:NGr, function(x) { tst2$gap[x] >= tst2$gap[x+1] - tst2$SE.sim[x+1] })
    tst2 <- which(tst2)
    # I like to do one more, often these methods feel too conservative
    #if (length(tst2)) { NVClust[[i]] <- tst2[1]+1 } # Not used for now: we use NGr instead
    vplot <- fviz_gap_stat(tst)
    tstLy <- capture.output(print(vplot$layers))
    g1 <- grep("geom_vline", tstLy)
    g2 <- grep("\\[\\[[0-9]+\\]\\]", tstLy)
    g2 <- as.numeric(gsub("\\[|\\]", "", tstLy[max(g2[which(g2 < g1)])]))
    vplot$layers[[g2]] <- NULL
    vplot <- vplot +
      geom_vline(xintercept = tst2[1]+1, colour = "red", linetype = "dashed") +
      geom_vline(xintercept = NGr, colour = "orange", linetype = "dashed") +
      geom_text(label = "Optimal number of sample clusters", x = tst2[1]+1-0.2, y = yScl*0.9,
                colour = "red", angle = 90, hjust = 1) +
      geom_text(label = "Number of sample groups", x = NGr+0.2, y = yScl*0.9,
                colour = "orange", angle = 90, hjust = 1) +
      theme_bw() + ggtitle(vnm)
    #poplot(vplot)
    ggsave(paste0(dir, "/", vnm, normTypeInsrt, ".jpeg"), vplot, dpi = 150)
    ggsave(paste0(dir, "/", vnm, normTypeInsrt, ".pdf"), vplot, dpi = 150)
    NVClust[[i]] <- max(c(NGr, 2))
    # 2/ At protein groups level
    # As above, we always draw a dendrogram, but colours will be defined by the clustering approach.
    hcluster <- hclust(dist(temp3))
    hdendro <- as.dendrogram(hcluster)
    hnm <- paste0("Protein groups-level clusters number analysis - ", i)
    NHClust[[i]] <- MaxHClust
    Straps <- 10
    # Here we really want to optimize the number of clusters
    # Apply the same method for optimization for any clustering method
    # Number of cluster should not depend on method
    source(parSrc, local = FALSE)
    clusterExport(parClust, list("temp3", "Straps"), envir = environment())
    tst <- setNames(parLapply(parClust, 2:MaxHClust, function(kl) { #kl <- 2
      try(kmeans(temp3, kl, nstart = Straps)$tot.withinss, silent = TRUE)
    }), 2:MaxHClust)
    tst <- tst[which(sapply(tst, function(x) { !"try-error" %in% class(x) }))]
    tst <- setNames(as.numeric(tst)/kmeans(temp3, 1, nstart = 1)$tot.withinss, names(tst))
    yScl2 <- max(tst)
    tst2 <- data.frame("Number of clusters k" = as.integer(names(tst)),
                       "[tot WSS (k)]/[tot WSS (1)]" = tst,
                       check.names = FALSE)
    
    tst <- parSapply(parClust, 2:MaxHClust, function(kl) {
      kmeans(temp3, kl, nstart = Straps)$tot.withinss
    })/kmeans(temp3, 1, nstart = 1)$tot.withinss
    yScl2 <- max(tst)
    tst2 <- data.frame("Number of clusters k" = 2:MaxHClust,
                       "[tot WSS (k)]/[tot WSS (1)]" = tst,
                       check.names = FALSE)
    # For the elbow detection method, we need to normalize to a 1x1 graph which we can rotate by 45 degrees
    tst2$X1 <- tst2$`Number of clusters k`/MaxHClust # divide by max theoretical number of clusters 
    tst2$Y1 <- tst2$"[tot WSS (k)]/[tot WSS (1)]" # Here no need to normalize, this is a ratio
    Angle <- -pi/4
    meanX <- mean(tst2$X1)
    meanY <- mean(tst2$Y1)
    tst2$X2 <- tst2$X1 - meanX
    tst2$Y2 <- tst2$Y1 - meanY
    tst2$X2 <- tst2$X2*cos(Angle)+tst2$Y2*sin(Angle)
    tst2$Y2 <- -tst2$X2*sin(Angle)+tst2$Y2*cos(Angle)
    tst2$X2 <- tst2$X2 - mean(tst2$X2) + meanX
    tst2$Y2 <- tst2$Y2 - mean(tst2$Y2) + meanY
    w <- rev(which(tst2$Y2 == min(tst2$Y2)))[1] # Again, prefer more rather than fewer clusters
    NHClust[[i]] <- tst2$`Number of clusters k`[w]
    tst2$Size <- 1
    tst2$Size[w] <- 2
    xMin <- min(c(tst2$X1, tst2$X2))
    xMax <- max(c(tst2$X1, tst2$X2))
    xScl <- xMax-xMin
    yMin <- min(c(tst2$Y1, tst2$Y2))
    yMax <- max(c(tst2$Y1, tst2$Y2))
    yScl <- yMax-yMin
    hplot <- ggplot(tst2) +
      geom_segment(x = tst2$X2[w], y = tst2$Y2[w],
                   xend = tst2$X1[w], yend = tst2$Y1[w], color = "grey", linetype = "dotted") +
      geom_point(aes(x = X1, y = Y1, size = Size), color = "blue") +
      geom_point(aes(x = X2, y = Y2, size = Size), color = "red") +
      geom_text(x = xScl*0.98+xMin, y = yMin+yScl*0.98, color = "blue", label = "ratio of tot. WSS", hjust = 1) +
      geom_text(x = xScl*0.98+xMin, y = yMin+yScl*0.962, color = "red", label = "ratio of tot. WSS, -pi/4 rotation", hjust = 1) +
      geom_hline(yintercept = tst2$Y2[w], color = "red", linetype = "dashed") +
      geom_vline(xintercept = tst2$X1[w], color = "deepskyblue", linetype = "dashed") +
      geom_text(x = xMin+0.01*xScl, y = tst2$Y2[w]+0.02*yScl, color = "red", label = "Elbow", hjust = 0) +
      geom_text(x = tst2$X1[w]-0.02*xScl, y = yScl*0.9, angle = 90, color = "deepskyblue",
                label = paste0("Optimal = ", tst2$`Number of clusters k`[w], " clusters"), hjust = 1) +
      ggtitle(hnm) + scale_size_identity() + theme_bw() +
      theme(legend.position = "none") + ylab("Normalised total Within-clusters vs Total Sum of Squares")
    #poplot(hplot)
    ggsave(paste0(dir, "/", hnm, normTypeInsrt, ".jpeg"), hplot, dpi = 150)
    ggsave(paste0(dir, "/", hnm, normTypeInsrt, ".pdf"), hplot, dpi = 150)
    # Apply cutoffs
    if (KlustMeth == 1) {
      VClusters[[i]] <- kmeans(t(temp3), NVClust[[i]], 100)$cluster
      tempClust <- kmeans(temp3, NHClust[[i]], 100)$cluster
    }
    if (KlustMeth == 2) {
      VClusters[[i]] <- cutree(vcluster, NVClust[[i]])
      tempClust <- cutree(hcluster, NHClust[[i]])
    }
    KlKol <- paste0(KlustRoot, i)
    KlustKols <- unique(c(KlustKols, KlKol))
    PG[[KlKol]] <- tempClust[match(PG$Label, names(tempClust))]
    #
    Width <- nrow(temp)
    Height <- ncol(temp)
    #
    whImps <- which(imputed, arr.ind = TRUE)
    temp[whImps] <- NA
    #
    # Get dendrograms
    vddata <- dendro_data(vdendro)
    hddata <- dendro_data(hdendro)
    # vdendro.plot <- ggdendrogram(data = vdendro) +
    #   theme(axis.text.y = element_text(size = 0.1), plot.margin = margin(0, 0, 0, 0, "cm"))
    # poplot(vdendro.plot)
    # Modify dendrograms
    # - Vertical
    vlabs <- ggdendro::label(vddata)
    vlabs$Cluster <- as.factor(VClusters[[i]][match(vlabs$label, smpls)])
    vSeg <- vddata$segments
    # - Horizontal
    hlabs <- ggdendro::label(hddata)
    hlabs$Cluster <- as.factor(tempClust[match(hlabs$label, names(tempClust))])
    hSeg <- hddata$segments
    # Rotate samples-level dendrogram: x -> y, y -> x (for ease of calculation, gets inverted later)
    vlabs$y <- vlabs$x
    vlabs$x <- -0.3
    x <- vSeg$x
    y <- vSeg$y
    vSeg$y <- x
    vSeg$x <- y
    xend <- vSeg$xend
    yend <- vSeg$yend
    vSeg$yend <- xend
    vSeg$xend <- yend
    # Adjust width/height
    # - Vertical dendrogram
    vnch <- Width*0.07*max(nchar(vlabs$label))/12
    xmn <- min(c(vSeg$x, vSeg$xend))
    xmx <- max(c(vSeg$x, vSeg$xend))
    vSeg$x <- -(vSeg$x - xmn)*Width*0.07/xmx - 1.5 - vnch
    vSeg$xend <- -(vSeg$xend - xmn)*Width*0.07/xmx - 1.5 - vnch
    # - Horizontal dendrogram
    hnch <- Height*0.07*max(nchar(hlabs$label))/800
    ymn <- min(c(hSeg$y, hSeg$yend))
    ymx <- max(c(hSeg$y, hSeg$yend))
    hSeg$y <- Height + hnch + 1.5 + (hSeg$y - ymn)*Height*0.07/ymx
    hSeg$yend <- Height + hnch + 1.5 + (hSeg$yend - ymn)*Height*0.07/ymx
    # Order labels by order of appearance
    hlabs <- hlabs[order(hlabs$x, decreasing = FALSE),]
    vlabs <- vlabs[order(vlabs$y, decreasing = FALSE),]
    # Re-order our matrix based on extracted dendrogram labels
    temp <- temp[, match(vlabs$label, colnames(temp))]
    temp <- temp[match(hlabs$label, rownames(temp)),]
    imputed <- imputed[, match(vlabs$label, colnames(imputed))]
    imputed <- imputed[match(hlabs$label, rownames(imputed)),]
    # Re-introduce missing values
    MaxChar <- 13
    hlabs$label2 <- hlabs$label
    w <- which(nchar(hlabs$label2) > MaxChar)
    hlabs$label2[w] <- paste0(substr(hlabs$label2[w], 1, MaxChar-3), "...")
    # Create heatmap
    temp$Rowname <- row.names(temp)
    temp2 <- set_colnames(reshape2::melt(temp, id.vars = "Rowname"), c("Label", "Sample", "value"))
    temp2$Label <- as.character(temp2$Label)
    temp2$Sample <- as.character(temp2$Sample)
    temp2$"Leading protein IDs" <- PG$"Leading protein IDs"[match(temp2$Label, PG$Label)]
    temp2$Xmax <- match(temp2$Label, hlabs$label) # Explicitly should be the case now!
    temp2$label2 <- hlabs$label2[temp2$Xmax]
    temp2$Xmin <- temp2$Xmax-1
    temp2$Ymax <- vlabs$y[match(temp2$Sample, vlabs$label)]
    temp2$Ymin <- temp2$Ymax-1
    w1 <- which(temp2$Colour == "green")
    w2 <- which((temp2$Ymin == max(temp2$Ymin))&(temp2$Colour == "green"))
    # Color and fill scales
    wV <- round(c(1:NVClust[[i]])*MaxVClust/NVClust[[i]])
    wH <- round(c(1:NHClust[[i]])*MaxHClust/NHClust[[i]])
    vClScl <- setNames(VClustScl[wV], 1:length(wV))
    hClScl <- setNames(HClustScl[wH], 1:length(wH))
    VcolScale <- scale_color_manual(name = "Samples cluster", values = vClScl)
    VfillScale <- scale_fill_manual(name = "Samples cluster", values = vClScl)
    HcolScale <- scale_color_manual(name = "Protein groups cluster", values = hClScl)
    HfillScale <- scale_fill_manual(name = "Protein groups cluster", values = hClScl)
    # Create heatmap plot
    Xlim <- c(NA, Width)
    Ylim <- c(-10, max(c(max(hSeg$y) + Height*0.6), 20))
    #
    prot.list.marks <- FALSE
    if (prot.list.Cond) {
      w0 <- which(temp2$Ymin == 0)
      g <- grsep2(prot.list, temp2$"Leading protein IDs"[w0])
      #tst <- unlist(strsplit(temp2$"Leading protein IDs"[w0], ";"))[1]
      #g <- grsep2(tst, temp2$"Leading protein IDs"[w0])
      if (length(g)) {
        prot.list.marks <- TRUE
        Ylim[1] <- -20
        temp2c <- temp2[w0[g], , drop = FALSE]
      }
    }
    # Main data
    temp2a <- temp2[, c("Xmin", "Ymin", "value", "Label", "Sample")]
    # Colour scale
    temp2b <- data.frame(Xmin = 0:round(Width*0.1),
                         Ymin = Ylim[2]*0.8)
    Mn <- min(temp2a$value, na.rm = TRUE)
    Mx <- max(temp2a$value, na.rm = TRUE)
    temp2b$value <- Mn + temp2b$Xmin*(Mx-Mn)/max(temp2b$Xmin)
    temp2b$Xmin <- temp2b$Xmin-Width*0.15
    temp2b$Label <- temp2b$Sample <- NA
    w2a <- 1:nrow(temp2a)
    w2b <- 1:nrow(temp2b) + max(w2a)
    temp2a <- rbind(temp2a, temp2b)
    # Samples-level: how well do clusters fit expectations
    SamplesClust <- vlabs[, c("label", "Cluster")]
    SamplesClust$Group <- as.factor(xMap[match(SamplesClust$label, xMap$Samples), VPAL$column])
    SamplesClust$Cluster <- as.factor(paste0("Cluster", as.character(SamplesClust$Cluster)))
    tstSmplClust <- table(SamplesClust$Group, SamplesClust$Cluster)
    ClustChiSqTst <- suppressWarnings(chisq.test(tstSmplClust))
    #
    xCntr <- Width*0.6
    yScale <- Height*2
    nm <- paste0("Clust. heatmap - ", i)
    GrLab <- data.frame(label = c(nm,
                                  normType,
                                  "Protein group",
                                  "Sample",
                                  paste0("Chi-squared contingency test P-value: ", round(ClustChiSqTst$p.value, 5)),
                                  "Min",
                                  "Max"),
                        x = c(xCntr,
                              xCntr,
                              xCntr,
                              min(vSeg$x)-Width*0.11,
                              min(vSeg$x)-Width*0.095,
                              min(temp2b$Xmin),
                              max(temp2b$Xmin)),
                        y = c(max(c(hSeg$y, hSeg$yend)) + yScale*0.15,
                              max(c(hSeg$y, hSeg$yend)) + yScale*0.125,
                              max(c(hSeg$y, hSeg$yend)) + yScale*0.1,
                              Height*0.5,
                              Height*0.5,
                              temp2b$Ymin[1]-1,
                              temp2b$Ymin[1]-1),
                        angle = c(0, 0, 0, 90, 90, 0, 0),
                        size = c(5, 3.5, 4, 4, 3, 3, 3),
                        fontface = c("bold", "italic", "plain", "plain", "italic", "plain", "plain"))
    Xlim[1] <- min(c(vSeg$x, vSeg$xend))-Width*0.15
    Ylim[1] <- min(c(Ylim[1], min(temp2a$Ymin)))
    Xlim[2] <- max(c(Xlim[2], max(temp2a$Xmin)+1))
    Ylim[2] <- max(c(Ylim[2], max(temp2a$Ymin)+1))
    xCntr <- mean(Xlim)
    yCntr <- mean(Ylim)
    yScale <- Ylim[2]-Ylim[1]
    #
    # Create graph
    heatmap.plot <- ggplot(temp2a[w2a,]) +
      geom_rect(aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value, text = Label)) +
      geom_rect(data = temp2a[w2b,], aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.margin = margin(0, 0, 0, 0, "cm")) +
      #scale_fill_gradient2(low = "darkblue", mid = "lightgrey", high = "darkred") +
      scale_fill_viridis(option = "D") +
      xlab(NULL) + ylab(NULL) + theme(legend.position = "none") +
      xlim(Xlim[1], Xlim[2]) + ylim(Ylim[1], Ylim[2])
    #poplot(heatmap.plot, 12, 20)
    # heatmap.plot <- heatmap.plot + 
    #   geom_text(data = temp2a[which(temp2a$Xmin == min(temp2$Xmin)),],
    #             aes(x = Xmin, y = Ymin, label = Sample))
    # Title, axis labels, colour scale annotations
    #
    # Add labels and dendrograms:
    # - Title, scale, chi-squared test
    heatmap.plot <- heatmap.plot +
      new_scale("size") + scale_size_identity() +
      geom_text(data = GrLab, aes(label = label, x = x, y = y, angle = angle, size = size, fontface = fontface),
                color = "black", hjust = 0.5)
    #poplot(heatmap.plot, 12, 20)
    # - Cluster boxes
    heatmap.plot <- heatmap.plot + new_scale_color() + HcolScale
    if (KlustMeth == 2) {
      Clutst <- aggregate(hlabs$x, list(hlabs$Cluster), function(x) { min(x)-1 })
      colnames(Clutst)[1] <- "Cluster"
      Clutst$xend <- aggregate(hlabs$x, list(hlabs$Cluster), max)$x
      Clutst$mid <- (Clutst$xend+Clutst$x)/2
      heatmap.plot <- heatmap.plot +
        geom_rect(data = Clutst, aes(xmin = x, xmax = xend, colour = Cluster),
                  ymin = 0, ymax = Height, fill = NA) +
        geom_text(data = Clutst, aes(x = mid, label = Cluster, colour = Cluster),
                  y = Height/2, hjust = 0.5, vjust = 0.5, cex = 5)
      #poplot(heatmap.plot2, 12, 20)
    }
    #poplot(heatmap.plot, 12, 20)
    # - Horizontal dendrogram and protein  names
    heatmap.plot <- heatmap.plot +
      geom_segment(data = hSeg, linewidth = 0.5,
                   aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_text(data = hlabs, aes(x = x-0.5, label = label2, colour = Cluster),
                y = Height+0.05, angle = 90, cex = 0.5, hjust = 0, vjust = 0)
    # - Vertical dendrogram and sample names
    heatmap.plot <- heatmap.plot +
      geom_segment(data = vSeg, linewidth = 0.5,
                   aes(x = x, y = y - 0.5, xend = xend, yend = yend - 0.5)) +
      new_scale_color() + VcolScale +
      geom_text(data = vlabs, aes(y = y - 0.5, label = label, colour = Cluster),
                x = -0.5, hjust = 1, vjust = 0.5, cex = 2.5)
    # - Proteins of interest
    if (prot.list.marks) {
      heatmap.plot <- heatmap.plot +
        geom_point(data = temp2c, aes(x = Xmin+0.5), y = -0.5, colour = "red", fill = "red", shape = 17) +
        geom_text(data = temp2c, aes(x = Xmin+0.5, label = label2),
                  y = -1, colour = "red", angle = -60, hjust = 0, cex = 2)
    }
    #poplot(heatmap.plot, 12, 20)
    ggsave(paste0(dir, "/", nm, normTypeInsrt, ".jpeg"), heatmap.plot)
    ggsave(paste0(dir, "/", nm, normTypeInsrt, ".pdf"), heatmap.plot)
    #
    # Plotly version
    tempLy <- temp2a[w2a,]
    tempLy$Sample <- factor(temp2$Sample, levels = smpls)
    ##tempLy$Ymin <- tempLy$Ymin+0.5
    plotleatmap <- plot_ly(data = tempLy, x = ~Xmin, y = ~Ymin, z = ~value, type = "heatmap", hovertext = tempLy$Label)
    # I cannot find a way to remove tick marks!!!!!
    plLyV <- vSeg
    ##plLyV$x <- -Width*0.2-plLyV$x 
    ##plLyV$xend <- -Width*0.2-plLyV$xend
    plLyV$y <- plLyV$y-1
    plLyV$yend <- plLyV$yend-1 
    plotleatmap <- add_segments(plotleatmap, x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                                data = plLyV, inherit = FALSE, color = I("black"),
                                showlegend = FALSE)
    plLyH <- hSeg
    ##plLyH$x <- plLyH$x-0.5
    ##plLyH$xend <- plLyH$xend-0.5
    plLyH$y <- plLyH$y - hnch - 1.5
    plLyH$yend <- plLyH$yend - hnch - 1.5
    plotleatmap <- add_segments(plotleatmap, x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                                data = plLyH, inherit = FALSE, color = I("black"),
                                showlegend = FALSE)
    if (KlustMeth == 2) { # Cluster shapes do not appear to work
      plotleatmap <- add_segments(plotleatmap, x = ~ x, xend = ~ xend,
                                  y = I(-0.2*as.numeric(Clutst$Cluster))-1,
                                  yend = I(-0.2*as.numeric(Clutst$Cluster))-1, data = Clutst,
                                  inherit = FALSE, color = ~ Cluster, showlegend = FALSE)
    }
    vlabs2 <- vlabs
    vlabs2$x <- -vnch/2
    vlabs2$y <- vlabs2$y-1
    plotleatmap <- add_trace(plotleatmap, data = vlabs2, y = ~y, x = ~x, text = ~label,
                             color = I("black"), inherit = FALSE, type = "scatter",
                             mode = "text", showlegend = FALSE)
    plotleatmap <- layout(plotleatmap, title = nm,
                          xaxis = list(title = "Protein groups",
                                       tickmode = "array", tickvals = NULL, showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
                          yaxis = list(title = "Samples",
                                       tickmode = "array", tickvals = NULL, showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE))
    if (prot.list.marks) {
      temp2c$X <- (temp2c$Xmin+temp2c$Xmax)/2
      plotleatmap <- add_trace(plotleatmap, data = temp2c, y = I(-0.501), x = ~X,
                               text = ~Label, color = ~Label, inherit = FALSE,
                               type = "scatter", mode = "markers", showlegend = FALSE)
    }
    saveWidget(plotleatmap, paste0(dir, "/", nm, normTypeInsrt, ".html"))
    #system(paste0("open \"", dir, "/", nm, normTypeInsrt, ".html\""))
    plotLeatMaps[[i]][[normType]] <- plotleatmap
  }
}
saveFun(plotLeatMaps, file = paste0(dir, "/HeatMaps.RData"))
#
kol <- paste0(Prot.Expr.Root, Exp.map$Ref.Sample.Aggregate)
kol <- kol[which(kol %in% colnames(PG))]
temp <- PG[, c("Leading protein IDs", "Protein names", "Genes", "Mol. weight [kDa]",
               kol, KlustKols)]
write.csv(temp, paste0(dir, "/Protein Groups and Clusters.csv"), row.names = FALSE)
#

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)

#### Code chunk - Dimensionality reduction plots
# PCA plots, by sample
dimRedPlotLy %<o% list()
msg <- "PCA plots, by sample"
ReportCalls <- AddMsg2Report(Space = FALSE)
dir <- paste0(wd, "/Dimensionality red. plots/PCA")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
agg <- parse.Param.aggreg(Param_filter(Param$Ratios.Groups.Ref.Aggregate.Level, "Exp"))
kol <- paste0(prtRfRoot, unique(Exp.map$Ref.Sample.Aggregate))
kol <- kol[which(kol %in% colnames(PG))]
temp <- PG[, kol]
colnames(temp) <- gsub(topattern(prtRfRoot), "", colnames(temp))
w <- which((apply(temp, 1, function(x) {length(is.all.good(x))}) == ncol(temp))
           &(PG$`Potential contaminant` != "+"))
temp <- temp[w,]
# We are now (line below) normalizing by average abundance, just in case, however the effect seems minimal:
temp <- sweep(temp, 1, PG$"Av. log10 abundance"[w], "-")
temp <- temp + rnorm(length(unlist(temp)), 0, 10^-9) # To avoid constant/zero columns, add a small random error
pc <- prcomp(t(temp), scale. = TRUE)
scores <- as.data.frame(pc$x)
pv <- round(100*(pc$sdev)^2 / sum(pc$sdev^2), 0)
pv <- pv[which(pv > 0)]
pv <- paste0("Components: ", paste(sapply(1:length(pv), function(x) {
  paste0("PC", x, ": ", pv[x], "%")
}), collapse = ", "))
scores$Sample <- rownames(scores)
scores$Group <- cleanNms(sapply(scores$Sample, function(x) {
  Exp.map[which(Exp.map[[RSA$aggregate]] == x), VPAL$column]
}))
scores[, RSA$names] <- Isapply(strsplit(scores$Sample, "___"), unlist)
if (exists("Tim")) {
  scores$"Time.point" <- as.numeric(scores[[Aggregates[[which(names(Aggregates) == "Tim")]]]])
}
a <- RG$column
scores$Ratios.Groups <- sapply(scores$Sample, function(x) {
  Exp.map[which(Exp.map$Ref.Sample.Aggregate == x),a]
})
scores$Sample <- cleanNms(scores$Sample)
if (exists("Tim")) {
  if (length(Exp) > 1) { form <- "Experiment ~ Time.point" } else { form <- "~Time.point" }
} else { form <- "~Experiment" }
ttl <- "PCA plot - Samples (PG-level)"
if ("PC2" %in% colnames(scores)) {
  plot <- ggplot(scores) +
    geom_point(aes(x = PC1, y = PC2, color = Group)) +
    scale_color_viridis_d(begin = 0.25) +
    coord_fixed() +
    geom_hline(yintercept = 0, colour = "black", alpha = 0.5) +
    geom_vline(xintercept = 0, colour = "black", alpha = 0.5) +
    ggtitle(ttl, subtitle = pv) +
    geom_text_repel(aes(x = PC1, y = PC2, label = Sample),
                    size = 2.5, show.legend = FALSE) + 
    theme_bw()
  if (substr(form, 1, 1) == "~") { plot <- plot + facet_wrap(form) } else { plot <- plot + facet_grid(form) }
  #
  scores$"Samples group" <- factor(scores$Group)
  if ("PC3" %in% colnames(scores)) {
    Symb <- rep(c("circle", "diamond", "square", "cross", "x"), max(as.numeric(Rep)))[1:max(as.numeric(Rep))]             
    Symb <- Symb[as.numeric(scores$Replicate)]
    plot_lyPCAProt <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3,
                              color = ~`Samples group`, colors = "viridis",
                              text = ~Sample, type = "scatter3d", mode = "markers",
                              symbol = I(Symb))
    plot_lyPCAProt <- add_trace(plot_lyPCAProt, scores, x = ~PC1, y = ~PC2, z = ~PC3,
                                type = "scatter3d", mode = "text", showlegend = FALSE)
    plot_lyPCAProt <- layout(plot_lyPCAProt, title = ttl)
    saveWidget(plot_lyPCAProt, paste0(dir, "/", ttl, ".html"))
    dimRedPlotLy[["Samples PCA"]] <- plot_lyPCAProt
    #system(paste0("open \"", dir, "/", ttl, ".html"))
    # NB: There is currently no way to create a 3D, faceted plot in plotly for R that I know of) 
  } else { poplot(plot, width = 18) }
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ReportCalls <- AddPlot2Report()
} else { warning("PCA failed, investigate!") }

# PCA, t-SNE and UMAP plots, by protein group
msg <- "PCA plots, by protein group"
ReportCalls <- AddMsg2Report(Space = FALSE)
dir <- paste0(wd, "/Dimensionality red. plots/PCA")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
agg <- parse.Param.aggreg(Param_filter(Param$Ratios.Groups.Ref.Aggregate.Level, "Exp"))
kol <- paste0(prtRfRoot, unique(Exp.map$Ref.Sample.Aggregate))
kol <- kol[which(kol %in% colnames(PG))]
temp <- PG[, kol]
colnames(temp) <- gsub(topattern(prtRfRoot), "", colnames(temp))
filt <- which((apply(temp, 1, function(x) { length(is.all.good(x)) }) == ncol(temp))
              &(PG$`Potential contaminant` != "+"))
if (length(filt) > 2) {
  temp <- temp[filt,]
  rownames(temp) <- PG$Label[filt]
  # For this type of PCA normalizing properly seems crucial:
  temp <- sweep(temp, 1, PG$"Av. log10 abundance"[filt], "-")
  temp <- temp + rnorm(length(unlist(temp)), 0, 10^-9) # To avoid constant/zero columns, add a small random error
  pc <- prcomp(temp, scale. = TRUE)
  scores <- as.data.frame(pc$x)
  pv <- round(100*(pc$sdev)^2 / sum(pc$sdev^2), 0)
  pv <- pv[which(pv > 0)]
  pv <- paste0("Components: ", paste(sapply(1:length(pv), function(x) {
    paste0("PC", x, ": ", pv[x], "%")
  }), collapse = ", "))
  m <- match(rownames(scores), PG$Label)
  scores[, c("Protein group", "Av. log10 abundance")] <- PG[m, c(Param$Plot.labels, "Av. log10 abundance")]
  scores$Range <- PG$"Rel. av. log10 abundance"[m] # Useful to check if distribution correlates with expression.
  if (Annotate&&LocAnalysis) {
    # Define "Markers"
    CompGOTerms %<o% setNames(c("GO:0005634", "GO:0005654", "GO:0000785", "GO:0005730", "GO:0005635", "GO:0005737", "GO:0005829", "GO:0005783", "GO:0005794", "GO:0031988", "GO:0005739", "GO:0009536", "GO:0005886"),
                              c("Nucleus", "Nucleoplasm", "Chromatin", "Nucleolus", "Nuclear envelope", "Cytoplasm", "Cytosol", "ER", "Golgi", "Vesicle", "Mitochondrion", "Plastid", "Plasma membrane"))
    tstSp <- unlist(strsplit(Param$Search.DB.species, ";"))[1]
    tstSp <- tolower(unlist(strsplit(tstSp, " ")))
    tstSp <- paste0(substr(tstSp[1], 1, 1), substr(tstSp[2], 1, 3))
    tst2 <- capture.output(pRolocmarkers())
    tst2 <- tst2[2:length(tst2)]
    l2 <- length(tst2)
    tst2 <- data.frame(Species = tst2[c(1:(l2/2))*2-1],
                       IDsType = tst2[c(1:(l2/2))*2])
    tst2$SpTag <- gsub("^[^\\[]+\\[|(\\]|_).*:$", "", tst2$Species)
    if (tstSp %in% tst2$SpTag) {
      SubCellMark <- pRolocmarkers(tstSp)
      tst3 <- aggregate(SubCellMark, list(SubCellMark), length)
      tmp <- strsplit(PG$`Leading protein IDs`, ";")
      tmp <- listMelt(tmp)
      tmp <- tmp[which(tmp$value %in% names(SubCellMark)),]
      tmp$Comp <- SubCellMark[match(tmp$value, names(SubCellMark))]
      tst <- aggregate(tmp$L1, list(tmp$Comp), function(x) { list(x) })
      tst <- setNames(tst$x, tst$Group.1)
      ClassNm <- paste(tst$Group.1, collapse = " / ")
      scores$Classifier <- rep("", nrow(scores))
      for (comp in names(tst)) {
        w1 <- which(rownames(scores) %in% PG$Label[unlist(tst[[comp]])])
        scores$Classifier[w1] <- comp
      }
      SubCellMark2 %<o% listMelt(strsplit(PG$`Leading protein IDs`, ";"), PG$id)
      SubCellMark2 <- SubCellMark2[which(SubCellMark2$value %in% names(SubCellMark)),]
      SubCellMark2$Comp <- SubCellMark[match(SubCellMark2$value, names(SubCellMark))]
      SubCellMark2$L1 <- PG$Label[match(SubCellMark2$L1, PG$id)]
      SubCellMark2 <- aggregate(SubCellMark2$Comp, list(SubCellMark2$L1), unique)
      SubCellMark2 <- SubCellMark2[which(sapply(SubCellMark2$x, length) == 1),]
      SubCellMark2$x <- unlist(SubCellMark2$x)
      SubCellMark2 <- setNames(SubCellMark2$x, SubCellMark2$Group.1)
    } else {
      tst <- setNames(lapply(CompGOTerms, function(x) { grep(x, PG$"GO-ID") }), names(CompGOTerms))
      tst2 <- unique(unlist(tst))
      tst2 <- tst2[which(sapply(tst2, function(x) { sum(sapply(tst, function(y) { x %in% y })) }) == 1)]
      tst <- lapply(tst, function(x) { x[which(x %in% tst2)] })
      tst <- tst[which(sapply(tst, length) > 0)]
      SubCellMark2 %<o% listMelt(tst)
      SubCellMark2$value <- PG$Label[SubCellMark2$value]
      SubCellMark2 <- setNames(SubCellMark2$L1, SubCellMark2$value)
    }
    ClassNm <- "Compartment"
    scores$Classifier <- rep("", nrow(scores))
    for (comp in names(tst)) {
      w1 <- which(rownames(scores) %in% PG$Label[unlist(tst[[comp]])])
      scores$Classifier[w1] <- comp
    }
    val <- c("", unique(scores$Classifier))
  } else {
    g <- grep("^Regulated - ", colnames(PG), value = TRUE)
    if (length(g) <= 6) {
      ClassNm <- cleanNms(gsub("^Regulated - ", "", g))
      tmp <- do.call(cbind, (lapply(g, function(x) {
        gsub_Rep("^Anti-specific: .*", "down", gsub_Rep("^Specific: .*", "up", gsub_Rep(", FDR = .+|^(non significant)|(too small FC)$", "", PG[w, x])))
      })))
      tmp <- set_colnames(as.data.frame(tmp), g)
      ClassNm <- paste(ClassNm, collapse = " / ")
      scores$Classifier <- do.call(paste, c(tmp, sep = " / "))
      val <- unique(c(paste(rep(" ", length(g)), collapse = "/"), scores$Classifier))
    } else {
      ClassNm <- "Regulated"
      tmp <- do.call(cbind, (lapply(g, function(x) {
        gsub_Rep("^Anti-specific: .*", "down", gsub_Rep("^Specific: .*", "up", gsub_Rep(", FDR = .+|^(non significant)|(too small FC)$", "", PG[w, x])))
      })))
      scores$Classifier <- apply(tmp, 1, function(x) {
        x <- x[which(x != "")]
        if (length(x)) { x <- "+" } else { x <- "" }
      })
      val <- unique(scores$Classifier)
    }
  }
  ttl <- "PCA plot - Protein groups (PG-level)"
  myColors <- c("black", rainbow(length(val)-1))
  names(myColors) <- val
  colScale <- scale_colour_manual(name = "colour", values = myColors)
  scores <- rbind(scores[which(scores$Classifier == val[1]),],
                  scores[which(scores$Classifier != val[1]),])
  scores$Alpha <- ((scores$Classifier != val[1])/2)+0.1
  wreg <- which(scores$Classifier != val[1])
  if (Mirror.Ratios) { ClassNm <- paste0("...vs: ", ClassNm) }
  plot <- ggplot(scores) +
    geom_scattermore(aes(x = PC1, y = PC2, colour = Classifier, 
                         alpha = Alpha, size = `Av. log10 abundance`),
                     shape = 16) +
    coord_fixed() + colScale +
    geom_hline(yintercept = 0, colour = "black", alpha = 0.5) +
    geom_vline(xintercept = 0, colour = "black", alpha = 0.5) +
    scale_alpha_identity() + scale_size_continuous(range = c(0.01, 3)) +
    ggtitle(ttl, subtitle = pv) + theme_bw() +
    guides(alpha = "none", colour = guide_legend(title = gsub("/", "/\n", ClassNm)))
  #poplot(plot, 12, 22)
  if (length(wreg)) {
    plot <- plot +
      geom_text_repel(data = scores[wreg[1:min(c(length(wreg), 50))],],
                      aes(x = PC1, y = PC2, label = `Protein group`),
                      size = 2.5, show.legend = FALSE)
  }
  #
  scores$Classifier <- factor(scores$Classifier)
  if ("PC3" %in% colnames(scores)) {
    col <- "#808080"
    if (length(wreg)) { col <- c(col, rainbow(length(val)-1)) }
    plot_lyPCAProt2 <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, #alpha = ~Alpha, (fails)
                               color = ~`Classifier`, colors = col, size = ~`Av. log10 abundance`,
                               text = ~paste0("<b>Protein group:</b> ", `Protein group`,
                                              "\n<b>", c("Class", "Compartment")[LocAnalysis+1],
                                              ":</b> ", Classifier),
                               type = "scatter3d", mode = "markers")
    plot_lyPCAProt2 <- layout(plot_lyPCAProt2, title = ttl,
                              legend = list(y = 0.5, title = list(text = paste0("<b>", ClassNm, "</b>"))))
    dimRedPlotLy[["PCA"]] <- plot_lyPCAProt2
    saveWidget(plot_lyPCAProt2, paste0(dir, "/", ttl, ".html"))
    #system(paste0("open \"", dir, "/", ttl, ".html"))
    # NB: There is currently no way to create a 3D, faceted plot in plotly for R that I know of) 
  } else { poplot(plot, 12, 22) }
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ReportCalls <- AddPlot2Report()
  #
  msg <- "t-SNE plots, by protein group"
  ReportCalls <- AddMsg2Report(Space = FALSE)
  cran_req <- unique(c(cran_req, "Rtsne"))
  if (!require("Rtsne", quietly = TRUE)) { install.packages("Rtsne") }
  require(Rtsne)
  tsne <- try(Rtsne(temp, dims = 3, perplexity = 30, verbose = TRUE, max_iter = 500), silent = TRUE)
  dir <- paste0(wd, "/Dimensionality red. plots/t-SNE")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  if (!"try-error" %in% class(tsne)) {
    temp2 <- as.data.frame(tsne$Y); colnames(temp2) <- c("t-SNE Y1", "t-SNE Y2", "t-SNE Y3")
    kol <- c("Protein group", "Av. log10 abundance", "Classifier", "Alpha")
    temp2[, kol] <- scores[match(rownames(temp), rownames(scores)), kol]
    ttl <- "t-SNE plot - Protein groups (PG-level)"
    plot <- ggplot(temp2) +
      geom_scattermore(aes(x = `t-SNE Y1`, y = `t-SNE Y2`, colour = Classifier, alpha = Alpha,
                           size = `Av. log10 abundance`), shape = 16) +
      colScale +
      geom_hline(yintercept = 0, colour = "black", alpha = 0.5) +
      geom_vline(xintercept = 0, colour = "black", alpha = 0.5) +
      scale_alpha_identity() + scale_size_continuous(range = c(0.01, 3)) +
      ggtitle(ttl) + theme_bw() +
      guides(alpha = "none", colour = guide_legend(title = gsub("/", "/\n", ClassNm)))
    if (length(wreg)) {
      plot <- plot +
        geom_text_repel(data = temp2[wreg[1:min(c(length(wreg), 50))],],
                        aes(x = `t-SNE Y1`, y = `t-SNE Y2`, colour = Classifier, label = `Protein group`),
                        size = 2.5, show.legend = FALSE)
    }
    #poplot(plot, width = 18)
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
    plot_lytSNE <- plot_ly(temp2, x = ~`t-SNE Y1`, y = ~`t-SNE Y2`, z = ~`t-SNE Y3`, color = ~Classifier, #alpha = ~Alpha, (fails)
                           colors = col, size = ~`Av. log10 abundance`,
                           text = ~paste0("<b>Protein group:</b> ", `Protein group`,
                                          "\n<b>", c("Class", "Compartment")[LocAnalysis+1],
                                          ":</b> ", Classifier),
                           type = "scatter3d",
                           mode = "markers")
    plot_lytSNE <- layout(plot_lytSNE, title = ttl,
                          legend = list(y = 0.5, title = list(text = paste0("<b>", ClassNm, "</b>"))))
    dimRedPlotLy[["t-SNE"]] <- plot_lytSNE
    saveWidget(plot_lytSNE, paste0(dir, "/", ttl, ".html"))
    #system(paste0("open \"", dir, "/", ttl, ".html"))
  } else { warning(tsne) }
  #
  msg <- "UMAP plots, by protein group"
  ReportCalls <- AddMsg2Report(Space = FALSE)
  cran_req <- unique(c(cran_req, "umap"))
  if (!require("umap", quietly = TRUE)) { install.packages("umap") }
  require(umap)
  UMAP <- try(umap(temp, n_components = 3), silent = TRUE)
  dir <- paste0(wd, "/Dimensionality red. plots/UMAP")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  if (!"try-error" %in% class(UMAP)) {
    UMAPlayout <- data.frame(UMAP$layout)
    kol <- c("Protein group", "Av. log10 abundance", "Classifier", "Alpha")
    UMAPlayout[, kol] <- scores[match(row.names(UMAPlayout), row.names(scores)), kol]
    ttl <- "UMAP plot - Protein groups (PG-level)"
    plot <- ggplot(UMAPlayout) +
      geom_scattermore(aes(x = X1, y = X2, colour = Classifier, alpha = Alpha, size = `Av. log10 abundance`),
                       shape = 16) +
      colScale +
      geom_hline(yintercept = 0, colour = "black", alpha = 0.5) +
      geom_vline(xintercept = 0, colour = "black", alpha = 0.5) +
      scale_alpha_identity() + scale_size_continuous(range = c(0.01, 3)) +
      ggtitle(ttl) + theme_bw() +
      guides(alpha = "none", colour = guide_legend(title = gsub("/", "/\n", ClassNm)))
    if (length(wreg)) {
      plot <- plot +
        geom_text_repel(data = UMAPlayout[wreg[1:min(c(length(wreg), 50))],],
                        aes(x = X1, y = X2, colour = Classifier, label = `Protein group`),
                        size = 2.5, show.legend = FALSE)
    }
    #poplot(plot, width = 18)
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
    plot_lyUMAP <- plot_ly(UMAPlayout, x = ~X1, y = ~X2, z = ~X3, color = ~Classifier, #alpha = ~Alpha, (fails)
                           colors = col, size = ~`Av. log10 abundance`,
                           text = ~paste0("<b>Protein group:</b> ", `Protein group`,
                                          "\n<b>", c("Class", "Compartment")[LocAnalysis+1],
                                          ":</b> ", Classifier),
                           type = "scatter3d",
                           mode = "markers")
    plot_lyUMAP <- layout(plot_lyUMAP, title = ttl,
                          legend = list(y = 0.5, title = list(text = paste0("<b>", ClassNm, "</b>"))))
    dimRedPlotLy[["UMAP"]] <- plot_lyUMAP
    saveWidget(plot_lyUMAP, paste0(dir, "/", ttl, ".html"))
    #system(paste0("open \"", dir, "/", ttl, ".html"))
  } else { warning(umap) }
} else { warning("Not enough observations to create Protein Groups-level dimensionality reduction plots!") }
ReportCalls <- AddSpace2Report()
saveFun(dimRedPlotLy, file = paste0(dir, "/DimRedPlots.RData"))

#### Code chunk - Protein group profile plots and sorting plots
rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
require(RColorBrewer)
require(colorspace)
QuantTypes %<o% c("Expression", "Coverage")
if (CreateMSMSKol) { QuantTypes <- c(QuantTypes, "Spectra") }
Qtstorg2 <- c()
if (tstorg) {
  PG$temp <- PG[[pgOrgKol]]
  PG$temp[which(PG$`Potential contaminant` == "+")] <- "Contaminant"
  tstorg2 <- aggregate(PG$temp, list(PG$temp), length)
  tstorg2 <- tstorg2[order(tstorg2$x, decreasing = TRUE),]
  tstorg2 <- tstorg2$Group.1[which(tstorg2$x > 1)]
  tstorg2 <- tstorg2[which(tstorg2 != "Contaminant")]
}
MainDir <- paste0(wd, "/Sorting plots")
MainDir2 <- paste0(wd, "/Profile plots")
dirlist <- unique(c(dirlist, MainDir, MainDir2))
ggQuant %<o% list()
ggProf %<o% list()
QuantLy %<o% list()
ProfLy %<o% list()
clusterCall(parClust, function() library(ggplot2))
cat("Drawing protein group profile plots\n")
for (QuantType in QuantTypes) { #QuantType <- "Expression"
  cat(" ->", QuantType, "\n")
  QuantLy[[QuantType]] <- list()
  myColors <- setNames("black", "-")
  if (length(tstorg2)) {
    if (length(tstorg2) >= 3) {
      if (length(tstorg2) <= 12) {
        myColors <- c(myColors, setNames(brewer.pal(min(c(12, length(tstorg2))), "Set3"), tstorg2[min(c(12, length(tstorg2)))]))
      }
    } else {
      myColors <- c(myColors, setNames(c("#8DD3C7", "#FFFFB3")[1:length(tstorg2)], tstorg2)) # We never expect more than a handful of organisms  
    }
    myColors[["Contaminant"]] <- "deepskyblue"
  }
  myColors[[c("+", "In list")[(length(tstorg2) > 0)+1]]] <- "red"
  myColors["Specific"] <- "purple"
  SubDir <- paste0(MainDir, "/", QuantType)
  dirlist <- unique(c(dirlist, SubDir))
  if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
  kolnm <- c("log10 LFQ", "Coverage [%]", "Spectral count")[match(QuantType, QuantTypes)]
  ref <- c(paste0("Mean ", prtRfRoot), #prtRfRoot,
           "Sequence coverage [%] - ",
           "Spectral count - ")[match(QuantType, QuantTypes)]
  #ref <- c(prtRfRoot, "Sequence coverage [%] - ", "Spectral count - ")[match(QuantType, QuantTypes)]
  Agg <- get(c("VPAL",# "RSA",
               "RSA",
               "RSA")[match(QuantType, QuantTypes)])
  #Agg <- RSA
  g <- paste0(ref, Agg$values)
  Wh <- which(g %in% colnames(PG))
  g <- g[Wh]
  if (length(g)) {
    varkol <- c("Leading protein IDs", "Protein IDs", "id", "Label")
    temp <- PG[, c(varkol, g)]
    colnames(temp) <- gsub(topattern(ref), "", colnames(temp))
    test <- aggregate(temp$"Label", list(temp$"Label"), length)
    w <- which(test$x > 1)
    if (length(w)) {
      test <- test[w,]
      for (i in test$Group.1) {
        w <- which(temp$"Label" == i)
        temp$"Label"[w[2:length(w)]] <- paste0(temp$"Label"[w[2:length(w)]], "_", 2:length(w))
      }
    }
    temp <- set_colnames(reshape::melt.data.frame(temp, id.vars = varkol), # Checked
                         c(varkol[1:2], "id", "Protein Group", "Sample", "Y"))
    temp <- temp[which(is.all.good(temp$Y, 2)),]
    if (QuantType %in% c("Coverage", "Spectra")) { temp <- temp[which(temp$Y > 0),] }
    mID <- match(temp$id, PG$id)
    if (nrow(temp)) {
      temp$"In list" <- "-"
      if (length(tstorg2)) { temp$Category <- PG$temp[mID] } else {
        temp$Category <- "-"
      }
      if (length(prot.list)) {
        gA1 <- grsep2(prot.list, PG$"Protein IDs") # Better to use "Protein IDs", not "Leading ...", for proteins of interest
        if (length(gA1)) {
          wA1 <- which(temp$id %in% PG$id[gA1])
          temp$"In list"[wA1] <- "+"
          temp$Category[wA1] <- c("+", "In list")[(length(tstorg2) > 0)+1]
        }
      }
      lev <- c("In list", "+", tstorg2, "-", "Contaminant")
      lev <- lev[which(lev %in% unique(temp$Category))]
      temp$Category <- factor(temp$Category, levels = lev)
      temp2 <- data.table(Category = temp$Category, PG = temp$`Protein Group`)
      temp2 <- temp2[, list(x = unique(Category)), by = list(Group.1 = PG)]
      temp2 <- as.data.frame(temp2)
      temp2$Protein_group <- "darkgrey"
      w <- which(temp2$x == "Contaminant")
      suppressWarnings(temp2$Protein_group[w] <- brewer.pal(min(c(length(w), 12)), "Blues"))
      w <- which(temp2$x %in% tstorg2)
      temp2$Protein_group[w] <- rainbow_hcl(length(w))
      w <- which(temp2$x %in% c("In list",  "+"))
      temp2$Protein_group[w] <- "red"
      temp$Protein_group <- temp2$Protein_group[match(temp$`Protein Group`, temp2$Group.1)]
      exports <- list("Exp.map", "Agg", "myColors", "QuantType", "kolnm", "temp", "Exp", "SubDir", "plotEval", "cleanNms", "wd")
      if (QuantType == "Expression") {
        rgKol <- paste0("Regulated - ", Agg$values)
        myFilt <- which(rgKol %in% colnames(PG))
        rgKol <- rgKol[myFilt]
        tmpPG <- PG[, c("id", rgKol)]
        exports <- append(exports, "tmpPG")
      } else { myFilt <- 1:length(Agg$values) }
      clusterExport(parClust, exports, envir = environment())
      f0 <- function(val) {
      #for (val in Agg$values) {
        #val <- Agg$values[1]
        e <- Exp.map[which(Exp.map[[Agg$column]] %in% val),]
        temp2 <- temp[which(temp$Sample %in% e[[Agg$column]]),]
        temp2 <- temp2[order(temp2$Y, decreasing = TRUE),]
        temp2$"Protein Group" <- factor(temp2$"Protein Group", levels = temp2$"Protein Group")
        regtst <- (QuantType == "Expression")
        if (regtst) {
          rg <- paste0("Regulated - ", val)
          temp2$Regulated <- tmpPG[match(temp2$id, tmpPG$id), rg]
          spc <- grep("^Specific", temp2$Regulated)
          if (length(spc)) {
            temp3 <- temp2[spc,]
            temp3$Category <- "Specific"
          } else { regtst <- FALSE }
        }
        intmin <- floor(min(temp2$Y))
        intmax <- ceiling(max(temp2$Y))
        intscale <- intmax-intmin
        #xmax <- max(c(max(c(0, which(temp2$Category == "+")))+round(nrow(temp2)/15), nrow(temp2)))
        xmax <- nrow(temp2)*18/15
        ttl <- paste0("PGs by ", QuantType, " - ", proteoCraft::cleanNms(val))
        plot <- ggplot2::ggplot(temp2) +
          ggplot2::geom_bar(stat = "identity", ggplot2::aes(`Protein Group`, Y,  fill = Category)) +
          ggplot2::scale_fill_viridis_d(begin = 0.25) +
          ggplot2::annotate("text", nrow(temp2)/2, intmax*1.2+intscale*0.04,
                            label = paste0(nrow(temp2), " Protein Groups"),
                            hjust = 0.5, ) +
          ggplot2::theme_bw() + ggplot2::ggtitle(ttl) + ggplot2::ylab(kolnm) + 
          ggplot2::coord_cartesian(xlim = c(1, xmax),
                                   ylim = c(intmin, intmax*1.2+intscale*0.05)) +
          ggplot2::theme(panel.border = ggplot2::element_blank(),
                         panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         axis.line = ggplot2::element_line(colour = "black"),
                         axis.text.x = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank(),
                         plot.margin = ggplot2::margin(r = 100))
        if (regtst) {
          plot <- plot +
            ggplot2::geom_point(data = temp3,
                                ggplot2::aes(`Protein Group`, Y+0.05, fill = Category),
                                shape = 21)
        }
        #proteoCraft::poplot(plot, 12, 22)
        plot_ly <- plotly::ggplotly(plot, tooltip = "Protein Group")
        setwd(SubDir) # For some reason, unless I do this the default selfcontained = TRUE argument gets ignored and
        # a folder with external resources is created for each html plot!
        pth <- paste0(SubDir, "/", ttl)
        plPath <- paste0(pth, ".html")
        tstPL <- try(htmlwidgets::saveWidget(plot_ly, plPath, selfcontained = TRUE),
                     silent = TRUE)
        tstPL <- !"try-error" %in% class(tstPL)
        #system(paste0("open \"", SubDir, "/", ttl, ".html"))
        w <- which(temp2$Category == "In list")
        if (length(w)) {
          plot <- plot +
            ggplot2::geom_text(data = temp2[w,],
                               ggplot2::aes(`Protein Group`, Y + intmax*0.01, colour = Category,
                                            label = `Protein Group`),
                               angle = 45, hjust = 0, cex = 3.5)
        }
        evPlot <- plotEval(plot)
        RES <- list()
        ggplot2::ggsave(paste0(pth, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
        ggplot2::ggsave(paste0(pth, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")
        setwd(wd)
        # RES$JPEG <- list(Plot = evPlot,
        #                  Path = paste0(SubDir, "/", ttl),
        #                  Saved = FALSE,
        #                  Save = "jpeg",
        #                  Param = c(dpi = 150, width = 10, height = 10, units = "in"))
        # RES$PDF <- list(Plot = evPlot,
        #                 Path = paste0(SubDir, "/", ttl),
        #                 Saved = FALSE,
        #                 Save = "pdf")
      #}
        return(list(plotly = plot_ly,
                    plotly_saved = tstPL,
                    plotly_path = plPath,
                    ggplot = list(title = ttl, path = pth, plot = evPlot, dpi = 150, width = 10, height = 10)))
      }
      tmp <- try(parLapply(parClust, Agg$values[myFilt], f0), silent = TRUE)
      if ("try-error" %in% class(tmp)) {
        tmp <- try(lapply(Agg$values[myFilt], f0), silent = TRUE)  
      }
      if ("try-error" %in% class(tmp)) {
        stop("Raaaaaahhhh!!!! This error is caused by how pandoc (called by htmlwidgets above) handles memory. I hoped that the non-parallel lapply above could avoid the issue encountered with parLapply, but this error proves me wrong. If you encounter this, you probably need more RAM.")
      } else {
        wN <- sapply(tmp, function(x) { x$plotly_saved })
        if (length(wN)){
          # The hope here is maybe some succeeded and we can get away with 
          f0 <- function(x) {
            htmlwidgets::saveWidget(x$plotly, x$plotly_path, selfcontained = TRUE)
          }
          tst <- try(parLapply(parClust, tmp[wN], f0), silent = TRUE)
          if ("try-error" %in% class(tmp)) {
            tst <- try(lapply(tmp[wN], f0), silent = TRUE)  
          }
          if ("try-error" %in% class(tst)) {
            stop("Raaaaaahhhh!!!! This error is caused by how pandoc (called by htmlwidgets above) handles memory. I hoped that the non-parallel lapply above could avoid the issue encountered with parLapply, but this error proves me wrong. If you encounter this, you probably need more RAM.")
          }
        }
      }
      tmp <- setNames(tmp, Agg$values[myFilt])
      QuantLy[[QuantType]] <- setNames(lapply(Agg$values[myFilt], function(x) {
        tmp[[x]]$plotly
      }), Agg$values[myFilt])
      ggQuant[[QuantType]] <- setNames(lapply(Agg$values[myFilt], function(x) {
        tmp[[x]]$ggplot
      }), Agg$values[myFilt])
      if (length(g) > 1) {
        SubDir <- paste0(MainDir2, "/", QuantType)
        if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
        temp2 <- temp
        temp2 <- temp2[order(temp2$`Leading protein IDs`, temp2$Sample),]
        m <- match(temp2$Category, c("-", "Contaminant", tstorg2, "+", "In list"))
        temp2$LineType <- c(rep("dotted", 2), rep("dashed", length(tstorg2)), rep("solid", 2))[m]
        #temp2$DotSize <- c(rep(0.2, 2), rep(0.3, length(tstorg2)), rep(0.5, 2))[m]
        temp2$DotSize <- 1
        #temp2$Alpha <- c(rep(0.25, 2), rep(0.5, length(tstorg2)), rep(1, 2))[m]
        temp2$Alpha <- 1
        temp2$Sample <- factor(cleanNms(as.character(temp2$Sample)),
                               levels = cleanNms(Agg$values))
        WhLst <- which(temp2$Sample == rev(levels(temp2$Sample))[1])
        colScale <- scale_colour_manual(name = "Category", values = myColors)
        fillScale <- scale_fill_manual(name = "Category", values = myColors)
        ttl <- paste0("Protein group ",
                      c("Expression", "Coverage", "Spectral Count")[match(QuantType, QuantTypes)],
                      " profiles")
        suppressWarnings({
          plot <- ggplot(temp2, aes(x = Sample, y = Y, group = `Leading protein IDs`, colour = Protein_group,
                                    text = `Protein Group`)) +
            geom_line(aes(linetype = LineType), na.rm = TRUE) + geom_point(aes(size = DotSize), na.rm = TRUE) +
            ggtitle(ttl) + ylab(kolnm) + theme_bw() +
            facet_grid(~Category) + #colScale +
            scale_size_identity(guide = "none") + scale_linetype_identity(guide = "none") +
            theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
            scale_colour_identity()
        })
        #poplot(plot, 12, 22)
        setwd(SubDir)
        plotlyProfiles <- ggplotly(plot, tooltip = "text")
        # Grey plot for shiny app
        plotlyProfiles2 <- plot_ly(temp2)
        plotlyProfiles2 <- add_trace(plotlyProfiles2, x = ~Sample, y = ~Y,
                                     #split = ~`Protein Group`, # This would be correct but makes it slow, so this small approximation seems ok for now
                                     color = I("lightgrey"), type = "scatter",
                                     mode = "lines+markers", text = ~`Protein Group`, connectgaps = FALSE,
                                     name = "", showlegend = FALSE)
        ProfLy[[QuantType]] <- list(data = temp2,
                                    coloured = plotlyProfiles,
                                    grey = plotlyProfiles2)
        saveWidget(plotlyProfiles, paste0(SubDir, "/", ttl, ".html"))
        #system(paste0("open \"", SubDir, "/", ttl, ".html"))
        setwd(wd)
        plot <- plot + geom_text(data = temp2[WhLst,], aes(label = `Protein Group`, alpha = Alpha),
                                 hjust = 0, cex = 2.5) + scale_alpha_identity(guide = "none")
        evPlot <- plotEval(plot)
        #poplot(plot, 12, 22)
        pth <- paste0(SubDir, "/", ttl)
        ggsave(paste0(pth, ".jpeg"), plot, dpi = 150)
        ggsave(paste0(pth, ".pdf"), plot, dpi = 150)
        ggProf[[QuantType]] <- list(title = ttl, path = pth, plot = evPlot, dpi = 150, width = 10, height = 10)
      }
      ReportCalls <- AddSpace2Report()
    }
  }
}
ReportCalls <- AddSpace2Report()
PG$temp <- NULL
saveFun(ggQuant, file = paste0(MainDir, "/ggQuantPlots.RData"))
saveFun(ggProf, file = paste0(MainDir2, "/ggProfilePlots.RData"))
saveFun(QuantLy, file = paste0(MainDir, "/QuantPlots.RData"))
saveFun(ProfLy, file = paste0(MainDir2, "/ProfilePlots.RData"))
for (QuantType in QuantTypes) {
  tst <- parSapply(parClust, ggQuant[[QuantType]], function(x) {
    if (!dir.exists(x$path)) { dir.create(x$path, recursive = TRUE) }
    ggplot2::ggsave(paste0(x$path, ".pdf"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
    ggplot2::ggsave(paste0(x$path, ".jpeg"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
  })
}
tst <- parSapply(parClust, ggProf, function(x) {
  if (!dir.exists(x$path)) { dir.create(x$path, recursive = TRUE) }
  ggplot2::ggsave(paste0(x$path, ".pdf"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
  ggplot2::ggsave(paste0(x$path, ".jpeg"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
})
for (nm in names(ggProf)) {
  ReportCalls <- AddPlot2Report(Plot = ggProf[[nm]]$plot,
                                Dir = dirname(ggProf[[nm]]$path),
                                Title = ggProf[[nm]]$title)
}

# Visualize results
HEIGHT <- "500px"
allProt <- PG$Label
dfltPrt <- PG$Label[grsep2(prot.list[1], PG$"Protein IDs")]
appNm <- paste0(dtstNm, " - Data exploration")
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", "Data exploration"),
             appNm),
  br(),
  em("Here you can explore you data."), br(),
  em("Click "), actionButton("exitBtn", "exit"), em(" to continue along the workflow."), br(),
  br(),
  #
  h4("Clustering heatmap"),
  selectInput("HeatMap", "Select heatmap", names(plotLeatMaps), "Global"),
  radioButtons("NormType", "Normalisation method:", c("Norm. by row", "Z-scored"), "Norm. by row"),
  fluidRow(withSpinner(plotlyOutput("MYplotLeatMap", height = HEIGHT))),
  #
  h4("Dimensionality reduction plots"),
  selectInput("DimRedPlot", "Select method",
              names(dimRedPlotLy)[which(names(dimRedPlotLy) != "Samples PCA")],
              "PCA"),
  fluidRow(
    column(8, withSpinner(plotlyOutput("PGDimRed", height = HEIGHT))),
    column(4, withSpinner(plotlyOutput("SmplsDimRed", height = HEIGHT)))
  ),
  #
  h4("Sorted Quantification plots"),
  fluidRow(column(2, selectInput("QuantType", "Select Quant method", QuantTypes, "Expression")),
           column(4, pickerInput("ProfPlotProt", "Select protein(s) to display", allProt, dfltPrt, TRUE,
                                 pickerOptions(title = "Search me",
                                               `live-search` = TRUE,
                                               actionsBox = TRUE,
                                               deselectAllText = "Clear search",
                                               showTick = TRUE)))),
  withSpinner(plotlyOutput("ProfPlot", height = HEIGHT)),
  # fluidRow(withSpinner(uiOutput("SortedPGPlot", height = HEIGHT))),
  br(),
  br()
)
server <- function(input, output, session) {
  if (!exists("plotLeatMaps")) { loadFun(paste0(wd, "/Clustering/HeatMaps.RData")) }
  if (!exists("dimRedPlotLy")) { loadFun(paste0(wd, "/Dimensionality red. plots/DimRedPlots.RData")) }
  if (!exists("QuantLy")) { loadFun(paste0(wd, "/Sorting plots/QuantPlots.RData")) }
  if (!exists("ProfLy")) { loadFun(paste0(wd, "/Profile plots/ProfilePlots.RData")) }
  #
  #
  HEATMAP <- reactiveVal("Global")
  NORMTYPE <- reactiveVal("Norm. by row")
  DIMRED <- reactiveVal("PCA")
  QUANTTYPE <- reactiveVal("Expression")
  PROFPROT <- reactiveVal(prot.list[1])
  # SORTDSMPL <- reactiveVal(names(QuantLy[["Expression"]])[1])
  #
  updtHtMp <- function(reactive = TRUE) {
    if (reactive) {
      renderPlotly(plotLeatMaps[[HEATMAP()]][[NORMTYPE()]])
    } else {
      renderPlotly(plotLeatMaps[["Global"]][["Norm. by row"]])
    }
  }
  updtDimRed <- function(reactive = TRUE) {
    if (reactive) {
      renderPlotly(dimRedPlotLy[[DIMRED()]])
    } else {
      renderPlotly(dimRedPlotLy[["PCA"]])
    }
  }
  updtProfPlot <- function(reactive = TRUE) {
    if (reactive) {
      qunt <- QUANTTYPE()
      prt <- PROFPROT()
    } else {
      qunt <- "Expression"
      prt <- dfltPrt
    }
    pltly <- ProfLy[[qunt]]$grey
    if (length(prt)) {
      dat <- ProfLy[[qunt]]$data
      dat <- dat[which(dat$`Protein Group` %in% prt),]
      if (nrow(dat)) {
        #if (reactive) { print(PROFPROT()) }
        pltly <- add_trace(pltly, data = dat, x = ~Sample, y = ~Y,
                           split = ~`Protein Group`,
                           color = ~`Protein Group`, type = "scatter",
                           mode = "lines+markers", text = ~`Protein Group`, connectgaps = FALSE,
                           name = "")
      }
    }
    renderPlotly(pltly)
  }
  # updtSortPlot <- function(reactive = TRUE) {
  #   if (reactive) {
  #     ls <- list(list(
  #       selectInput("SortedSample", "Select sample", names(QuantLy[[QUANTTYPE()]]), SORTDSMPL()),
  #       fluidRow(withSpinner(plotlyOutput(renderPlotly(QuantLy[[QUANTTYPE()]][[SORTDSMPL()]]), height = HEIGHT)))
  #     ))
  #   } else {
  #     s <- list(list(
  #       selectInput("QuantPlot", "Select sample", names(QuantLy[["Expression"]]), names(QuantLy[["Expression"]])[1]),
  #       fluidRow(withSpinner(plotlyOutput(renderPlotly(QuantLy[["Expression"]][[1]]), height = HEIGHT)))
  #     ))
  #   }
  #   return(renderUI(ls))
  # }
  output$MYplotLeatMap <- updtHtMp(FALSE)
  output$PGDimRed <- updtDimRed(FALSE)
  output$SmplsDimRed <- renderPlotly(dimRedPlotLy[["Samples PCA"]])
  output$ProfPlot <- updtProfPlot(FALSE)
  # output$SortedPGPlot <- updtSortPlot(FALSE)
  observeEvent(input$HeatMap, {
    HEATMAP(input$HeatMap)
    output$MYplotLeatMap <- updtHtMp()
  })
  observeEvent(input$NormType, {
    NORMTYPE(input$NormType)
    output$MYplotLeatMap <- updtHtMp()
  })
  observeEvent(input$DimRedPlot, {
    DIMRED(input$DimRedPlot)
    output$PGDimRed <- updtDimRed()
  })
  observeEvent(input$QuantType, {
    QUANTTYPE(input$QuantType)
    output$ProfPlot <- updtProfPlot()
  })
  observeEvent(input$ProfPlotProt, {
    PROFPROT(input$ProfPlotProt)
    output$ProfPlot <- updtProfPlot()
  })
  # observeEvent(input$SortedSample, {
  #   SORTDSMPL(input$SortedSample)
  #   output$SortedPGPlot <- updtSortPlot()
  # })
  # Save
  observeEvent(input$exitBtn, { stopApp() })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
eval(parse(text = runApp), envir = .GlobalEnv)
#

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
    e <- sapply(c(1:length(o)), function(x) { which(Exp.map[[o[x]]] == i1[x]) })
    l <- unique(unlist(e))
    t <- sapply(l, function(x) {length(which(unlist(e) == x)) == length(o)})
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
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 600, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 600, width = 10, height = 10, units = "in")
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
Script <- readLines(ScriptPath)

#### Code chunk - Sub-Cellular localisation analysis: pRoloc-based prediction of localisation + analysis of the Sums of Squared Differences of Profiles
if (Annotate&&LocAnalysis) {
  msg <- "Sub-Cellular localisation analysis:"
  ReportCalls <- AddMsg2Report(Space = FALSE)
  dir <- paste0(wd, "/pRoloc")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  packs <- c("pRoloc", "pRolocGUI") # We will not use pRolocGUI as it seems to be using an outdated version of shinydashboardPlus (specifically a retired function)
  for (pack in packs) { #pack <- packs[1]
    bioc_req <- unique(c(bioc_req, pack))
    biocInstall(pack)
  }
  SubCellFracAggr %<o% parse.Param.aggreg(Param_filter(Param$Ratios.Groups.Ref.Aggregate.Level, "Com"))
  # Impute missing values
  AllKol <- paste0(prtRfRoot, RSA$values)
  w <- which(AllKol %in% colnames(PG))
  grps <- Exp.map[match(RSA$values[w], Exp.map$Ref.Sample.Aggregate), VPAL$column] 
  AllKol <- AllKol[w]
  wNC <- which(PG$`Potential contaminant` != "+")
  tempDat <- PG[wNC, AllKol]
  tempDat <- Data_Impute2(tempDat, grps)$Imputed_data
  #
  pRolocData <- SVMparams <- list()
  # For sub-cellular localisation analysis, we need to define a series of compartment markers to predict protein location
  # We can use the very granular markers already defined above (pRoloc or built-in), or define new ones here
  ObjNm <- "CompGOTerms2"
  .obj <- unique(c(.obj, ObjNm))
  if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
    msg <- "Enter a list of GO Cell Compartment (GO CC) terms for compartments of interest (semicolon-separated).
You may also include:
 - 1 for a low-resolution markers list (nucleoplasm - chromatin - cytoplasm)
 - 2 for a more exhaustive list.

Example: \"GO:0031012;2\"
"
    tmp <- unlist(strsplit(dlg_input(msg, 2)$res, "[;,] ?"))
    if ("1" %in% tmp) { tmp <- unique(c(tmp, "GO:0005654", "GO:0000785", "GO:0005737")) }
    if ("2" %in% tmp) { tmp <- unique(c(tmp, CompGOTerms)) }
    tmp <- tmp[which(tmp %in% GO_terms$ID[which(GO_terms$Ontology == "CC")])] # (Also neatly removes "1" and "2"...)
    ObjNm %<c% tmp
    AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
    tmp <- AllAnsw[1,]
    tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
    tmp$Value <- list(get(ObjNm))
    m <- match(ObjNm, AllAnsw$Parameter)
    if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
  }
  if (length(CompGOTerms2)) {
    CompGOTerms2 <- data.frame(Term = CompGOTerms2)
    CompGOTerms2$Offspring <- sapply(CompGOTerms2$Term, function(x) { GOCCOFFSPRING[[x]] })
    allOffspr <- unlist(CompGOTerms2$Offspring)
    allOffspr <- aggregate(allOffspr, list(allOffspr), length)
    CompGOTerms2$Offspring <- lapply(CompGOTerms2$Offspring, function(x) {
      x[which(!x %in% c(CompGOTerms2$Term, allOffspr$Group.1[which(allOffspr$x > 1)]))]
    })
    CompGOTerms2$Name <- gsub(" \\[GO:[0-9]{7}\\]$", "", GO_terms$Term[match(CompGOTerms2$Term, GO_terms$ID)])
    CompGOTerms2$All <- apply(CompGOTerms2[, c("Term", "Offspring")], 1, function(x) { unique(unlist(x)) })
    tst <- setNames(lapply(CompGOTerms2$All, function(x) { grep(paste(x, collapse = "|"), PG$"GO-ID") }),
                    CompGOTerms2$Name)
    tst2 <- unique(unlist(tst))
    tst2 <- tst2[which(sapply(tst2, function(x) { sum(sapply(tst, function(y) { x %in% y })) }) == 1)]
    tst <- lapply(tst, function(x) { x[which(x %in% tst2)] })
    tst <- tst[which(sapply(tst, length) > 0)]
    SubCellMark2 %<o% listMelt(tst) # Overwrite former value
    SubCellMark2$value <- PG$Label[SubCellMark2$value]
    SubCellMark2 <- setNames(SubCellMark2$L1, SubCellMark2$value)
    if (length(SubCellMark2)) {
      msg <- " - Performing per-sample pRoloc analysis"
      ReportCalls <- AddMsg2Report(Space = FALSE)
      ttls <- c()
      pRolocVisMeth <- "t-SNE"
      tmpPG <- PG[, c("id", "Leading protein IDs", "Protein IDs", "Protein names", "Genes", "Label")]
      exports <- list("SubCellMark2", "Exp.map", "SubCellFracAggr", "prtRfRoot", "tempDat", "VPAL", "tmpPG", "Aggregates",
                      "prot.list", "pRolocVisMeth", "Exp", "dir", "wNC")
      clusterExport(parClust, exports, envir = environment())
      clusterCall(parClust, function() {
        library(graphics)
        library(Biobase)
        library(MSnbase)
        library(pRoloc)
        library(proteoCraft)
      })
      f0 <- function(grp) { #grp <- SubCellFracAggr$values[1]
        ttl_s <- c()
        grp1 <- cleanNms(grp)
        w <- which(Exp.map[[SubCellFracAggr$column]] == grp)
        RES <- list(Success = FALSE)
        if (length(w)) {
          em <- Exp.map[w,]
          kol <- paste0(prtRfRoot, unique(em$Ref.Sample.Aggregate))
          kol <- kol[which(kol %in% colnames(tempDat))]
          temp1 <- tempDat[, kol]
          # Add a small noise to avoid non-unicity issue
          sd <- sd(unlist(temp1))
          noise <- rnorm(length(unlist(temp1)), 0, sd*10^-6)
          temp1 <- temp1 + noise
          #
          temp1 <- 10^temp1
          tst <- rowSums(temp1)
          wAG <- which((!is.na(tst))&(is.finite(tst))&(tst > 0))
          if (length(wAG)) {
            #cat(paste0(grp1, ":\n", paste(rep("-", nchar(grp1)+1), collapse = ""), "\n"))
            m <- match(unique(em[[VPAL$column]]), em[[VPAL$column]])
            temp1 <- sweep(temp1[wAG,], 1, tst[wAG], "/")
            rownames(temp1) <- tmpPG$Label[wAG]
            colnames(temp1) <- cleanNms(unique(em[[VPAL$column]]))
            temp1 <- as.matrix(temp1)
            temp2 <- em[m, Aggregates]
            rownames(temp2) <- cleanNms(unique(em[[VPAL$column]]))
            temp2 <- AnnotatedDataFrame(temp2)
            temp3 <- AnnotatedDataFrame(tmpPG[wAG, c("id", "Leading protein IDs", "Protein IDs", "Protein names", "Genes")])
            rownames(temp3) <- tmpPG$Label[wAG]
            MSnData <- MSnSet(exprs = temp1,
                              pData = temp2,
                              fData = temp3)
            MSnData <- addMarkers(MSnData, SubCellMark2)
            #pRolocData[[grp]] <<- MSnData
            #getMarkers(MSnData)
            # Define features of interest
            plotPLmtchs <- FALSE
            if (length(prot.list)) {
              PLmtchs <- grsep2(prot.list, tmpPG$"Leading protein IDs"[wAG])
              if (length(PLmtchs)) {
                plotPLmtchs <- TRUE
                foi1 <- FeaturesOfInterest(description = "Proteins of interest",
                                           fnames = featureNames(MSnData)[PLmtchs])
                description(foi1)
                foi(foi1)
              }
            }
            # the plotDist function is unsatisfactory
            #plotDist(MSnData, featureNames(MSnData))
            #
            # Markers hierarchical dendrogram
            ttl <- paste0("Markers hier. clust. - ", grp1)
            ttl2 <- paste0("Markers hierarchical clustering\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            mrkHClust(MSnData, main = ttl2)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #
            # Average markers class profile plots
            hc <- mrkHClust(MSnData, plot = FALSE)
            ## order of markers according to histogram
            mm <- getMarkerClasses(MSnData)
            m_order <- levels(factor(mm))[order.dendrogram(hc)]
            ## average marker profile
            fmat <- mrkConsProfiles(MSnData)
            ttl <- paste0("Marker classes av. profiles - ", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plotConsProfiles(fmat, order = m_order)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #
            # Dimensionality reduction
            pRolocVisMeth2 <- pRolocVisMeth
            if (pRolocVisMeth == "sLDA") { pRolocVisMeth2 <- "lda" }
            ttl <- paste0(pRolocVisMeth, " - ", grp1)
            ttl2 <- gsub(" - ", "\n", ttl)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plot2D(MSnData, fcol = "markers", main = ttl2, method = pRolocVisMeth2)
            addLegend(MSnData, fcol = "markers", cex = 0.7, where = "bottomright", ncol = 2)
            if (plotPLmtchs) {
              highlightOnPlot(MSnData, foi1, col = "black", lwd = 2)
              legend("topright", c("Proteins of interest"),
                     bty = "n", col = c("purple"),
                     pch = 1)
            }
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #plot3D(MSnData, fcol = "markers", main = ttl2) # Much worse - but also much older - than plotly
            #
            #pRolocVis(MSnData) # Not used: buggy
            #
            # Assess the resolution of the fractionation for the different compartments
            hlq <- QSep(MSnData)
            ttl <- paste0("Subcell. res. heatmap - ", grp1)
            ttl2 <- paste0("Sub-cellular resolution heatmap\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            levelPlot(hlq, main = ttl2)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            ttl <- paste0("Subcell. res. boxplot - ", grp1)
            ttl2 <- paste0("Sub-cellular resolution boxplot\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plot(hlq, main = ttl2)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #
            # Prediction of compartment assignment
            ## ... Unsupervised (kept as an example here, but doesn't look very useful)
            #kcl <- MLInterfaces::MLearn( ~ ., MSnData,  kmeansI, centers = length(Com))
            #plot(kcl, exprs(MSnData))
            #hcl <- MLInterfaces::MLearn( ~ ., MSnData,
            #                             MLInterfaces::hclustI(distFun = dist,
            #                       cutParm = list(k = length(Com))))
            #plot(hcl, labels = FALSE)
            #pcl <- MLearn( ~ ., MSnData,  pamI(dist), k = length(Com))
            #plot(pcl, data = exprs(MSnData))
            #
            ## ... Supervised
            # This is really, reaaally slow!
            cat(" Predicting compartment localisation using SVM (100 iterations), please wait...\n")
            wghts <- classWeights(MSnData, fcol = "markers")
            wghts <- wghts[which(wghts < 1)]
            params <- suppressMessages(suppressWarnings(svmOptimisation(MSnData, "markers", times = 10, xval = 5,
                                                                        class.weights = wghts, verbose = TRUE)))
            # The line above takes a millenium and a half
            # Unfortunately it cannot be parallelized or optimized easily without rewriting pRoloc
            # But this loop could be parallelized at least
            #
            #SVMparams[[grp]] <- params
            ttl <- paste0("SVM opt. boxplot - ", grp1)
            ttl2 <- paste0("SVM optimisation boxplot\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plot(params, main = ttl2)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            ttl <- paste0("SVM opt. heatmap - ", grp1)
            ttl2 <- paste0("SVM optimisation heatmap\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            levelPlot(params)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #f1Count(params)
            params2 <- getParams(params)
            MSnData <- svmClassification(MSnData, params,
                                         fcol = "markers", class.weights = wghts)
            #processingData(MSnData)
            p1 <- getPredictions(MSnData, fcol = "svm")
            p1 <- fData(p1)$svm.pred
            minprob <- median(fData(MSnData)$svm.scores)
            p2 <- getPredictions(MSnData, fcol = "svm", t = minprob)
            p2 <- fData(p2)$svm.pred
            #table(p1, p2)
            ptsze <- exp(fData(MSnData)$svm.scores) - 1
            ttl <- paste0(pRolocVisMeth, " with predictions - ", grp1)
            ttl2 <- gsub(" - ", "\n", ttl)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plot2D(MSnData, fcol = "markers", main = ttl2, method = pRolocVisMeth2)
            addLegend(MSnData, fcol = "markers", cex = 0.7, where = "bottomright", ncol = 2)
            if (plotPLmtchs) {
              highlightOnPlot(MSnData, foi1, col = "black", lwd = 2)
              legend("topright", c("Proteins of interest"),
                     bty = "n", col = c("purple"),
                     pch = 1)
              dev.off()
            }
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            lokol <- paste0("Localisation - ", grp)
            #PG[[lokol]] <<- ""
            w <- which(p2 != "unknown")
            #PG[wNC[wAG[w]], lokol] <<- p2[w]
            #pRolocData[[grp]] <<- MSnData
            #cat("\n")
            RES <- list(Success = TRUE,
                        Titles = ttl_s,
                        Which = wNC[wAG[w]],
                        Column = lokol,
                        Values = p2[w],
                        pRolocData = MSnData,
                        SVMparams = params)
          }
        }
        return(RES)
      }
      tmpClass <- parLapply(parClust, SubCellFracAggr$values, f0)
      names(tmpClass) <- SubCellFracAggr$values
      grps <- names(tmpClass)[which(sapply(tmpClass, function(x) { x$Success }))]
      for (grp in grps) {
        ttls <- c(ttls, tmpClass[[grp]]$Titles)
        lokol <- tmpClass[[grp]]$Column
        w <- tmpClass[[grp]]$Which
        PG[[lokol]] <- ""
        PG[w, lokol] <- tmpClass[[grp]]$Values
        pRolocData[[grp]] <- tmpClass[[grp]]$MSnData
        SVMparams[[grp]] <- tmpClass[[grp]]$SVMparams
      }
      ReportCalls <- AddSpace2Report()
      if (length(ttls)) {
        SilentPDF2JPEG <- function(ttl) {
          suppressMessages(
            suppressWarnings(
              pdf_convert(paste0(ttl, ".pdf"), "jpeg", filenames = paste0(ttl, ".jpeg"), dpi = 600)
            )
          )
        }
        clusterExport(parClust, list("dir", "SilentPDF2JPEG", "pdf_convert"), envir = environment())
        parSapply(parClust, ttls, function(ttl) { try(SilentPDF2JPEG(paste0(dir, "/", ttl)), silent = TRUE) })
      }
      w <- which(!PG$Label %in% names(SubCellMark2))
      #View(PG[w, grep("^Localisation - ", colnames(PG), value = TRUE)])
    }
  }
  #
  # Re-localisation analysis
  SubCellFracAggr2 %<o% parse.Param.aggreg(Param_filter(Param$Volcano.plots.Aggregate.Level, "Com"))
  SSD.Root %<o% "log10(SSD) - "
  SSD.Pval.Root %<o% "Welch's t-test on SSDs -log10(Pvalue) - "
  WhRef <- sapply(SubCellFracAggr2$values, function(x) { unique(Exp.map$Reference[which(Exp.map[[SubCellFracAggr2$column]] == x)]) })
  tst <- sapply(WhRef, length)
  if (max(tst) == 1) {
    tempDat2 <- 10^tempDat
    wh0 <- which(WhRef)
    wh1 <- which(!WhRef)
    if ((length(wh0) == 1)&&(length(wh1))) {
      msg <- " - Performing re-localisation analysis"
      ReportCalls <- AddMsg2Report(Space = FALSE)
      LocAnalysis2 <- TRUE
      EM0 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == SubCellFracAggr2$values[wh0]),]
      EM0$Replicate <- as.numeric(EM0$Replicate)
      EM0 <- EM0[order(EM0$Replicate, EM0$Compartment),]
      grp0 <- SubCellFracAggr2$values[wh0]
      # Currently, this only supports one type of grouping defined by SubCellFracAggr2,
      # itself equivalent to removing "Compartment" from Sample groups ("VPAL")
      comb <- gtools::combinations(max(EM0$Replicate), 2, unique(EM0$Replicate))
      # SSDs should normalize for expression level, otherwise we will pick up proteins whose expression changes!
      # Thus NormSSDs is always TRUE for now.
      NormSSDs <- TRUE
      if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
        RefSSDs %<o% NULL
      }
      if (Param$Ratios.Thresholds == threshMsg) {
        RefSSDs %<o% as.numeric(apply(comb, 1, function(i) { #i <- comb[1,]
          # NB: differs now from the way Ref.Ratios is written - but for a good reason.
          # We operate using different groupings for subcellular re-localisation analysis.
          A <- tempDat2[, paste0(prtRfRoot, EM0$Ref.Sample.Aggregate[which(EM0$Replicate == i[[1]])])]
          B <- tempDat2[, paste0(prtRfRoot, EM0$Ref.Sample.Aggregate[which(EM0$Replicate == i[[2]])])]
          if (NormSSDs) {
            A <- sweep(A, 1, rowSums(A), "/") 
            B <- sweep(B, 1, rowSums(B), "/") 
          }
          res <- log10(rowSums((A-B)^2))
          return(res)
        }))
        RefSSDs %<o% setNames(lapply(SubCellFracAggr2$values, function(x) { RefSSDs }), SubCellFracAggr2$values)
      }
      SSDs %<o% list()
      SSD.FDR.thresh %<o% c()
      for (wh in wh1) { #wh <- wh1[1]
        grp <- SubCellFracAggr2$values[wh]
        EM1 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == grp),]
        EM1$Replicate <- as.numeric(EM1$Replicate)
        EM1 <- EM1[order(EM1$Replicate, EM1$Compartment),]
        # Calculate sum of squared differences
        if (grepl("Rep", SubCellFracAggr$aggregate)) {
          rps <- unique(EM1$Replicate)
          SSDs[[grp]] <- as.data.frame(sapply(rps, function(rp) {
            P0 <- tempDat2[, paste0(prtRfRoot, EM0$Ref.Sample.Aggregate[which(EM0$Replicate == rp)])]
            P1 <- tempDat2[, paste0(prtRfRoot, EM1$Ref.Sample.Aggregate[which(EM1$Replicate == rp)])]
            if (NormSSDs) {
              P0 <- sweep(P0, 1, rowSums(P0), "/") 
              P1 <- sweep(P1, 1, rowSums(P1), "/") 
            }
            res <- P0-P1
            res <- log10(rowSums((P0-P1)^2))
            return(res)
          }))
          colnames(SSDs[[grp]]) <- SSDkols <- paste0(SSD.Root, EM1[match(rps, EM1$Replicate), SubCellFracAggr$column])
        } else {
          comb <- gtools::permutations(max(as.numeric(Rep)), 2, as.numeric(Rep), repeats.allowed = TRUE)
          SSDs[[grp]] <- as.data.frame(apply(comb, 1, function(i) {
            P0 <- tempDat2[, paste0(prtRfRoot, EM0$Ref.Sample.Aggregate[which(EM0$Replicate == i[[1]])])]
            P1 <- tempDat2[, paste0(prtRfRoot, EM1$Ref.Sample.Aggregate[which(EM1$Replicate == i[[2]])])]
            if (NormSSDs) {
              P0 <- sweep(P0, 1, rowSums(P0), "/") 
              P1 <- sweep(P1, 1, rowSums(P1), "/") 
            }
            res <- log10(rowSums((P0-P1)^2))
            return(res)
          }))
          colnames(SSDs[[grp]]) <- SSDkols <- paste0(SSD.Root, apply(comb, 1, function(i) {
            paste0(EM1[match(i[[2]], EM1$Replicate), SubCellFracAggr$column], "_vs_", EM0[match(i[[1]], EM0$Replicate), SubCellFracAggr$column])
          }))
        }
        SSDs[[grp]][[paste0("Mean ", SSD.Root, grp)]] <- apply(SSDs[[grp]][, SSDkols, drop = FALSE], 1, function(x) { log10(mean(10^x)) })
        PG[, colnames(SSDs[[grp]])] <- NA
        PG[wNC, colnames(SSDs[[grp]])] <- SSDs[[grp]]
      }
      # Check distribution
      # Test expression values:
      g <- c(grep(topattern("log10(SSD) - "), colnames(PG), value = TRUE),
             grep(topattern("Mean log10(SSD) - "), colnames(PG), value = TRUE))
      test <- PG[, g]
      colnames(test) <- gsub(topattern("log10(SSD) - ", start = FALSE), "", colnames(test))
      w <- grep("^Mean ", colnames(test))
      colnames(test)[w] <- paste0(gsub("^Mean ", "", colnames(test)[w]), "___Mean")
      test <- test[which(apply(test, 1, function(x) { length(is.all.good(x)) }) > 0),]
      test <- suppressMessages(melt.data.frame(test))
      test$variable <- as.character(test$variable)
      test[, SubCellFracAggr$names] <- ""
      w <- rep(FALSE, nrow(test))
      test[which(!w), SubCellFracAggr$names] <- Isapply(strsplit(test$variable[which(!w)], "___"), unlist)
      a <- SubCellFracAggr$names
      w <- which(sapply(a, function(x) { length(unique(test[[x]])) }) > 1)
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
        testI[[v]] <- sapply(1:nbinz, function(x) {
          sum((test$value[wv] > binz[x])&(test$value[wv] <= binz[x+1]))
        })
      }
      testI <- reshape2::melt(testI, id.vars = "Intensity")
      testI[, a] <- test[match(testI$variable, test$variable), a]
      testI$variable <- cleanNms(testI$variable)
      dir <- paste0(wd, "/Reg. analysis/Localisation")
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      ttl <- "Distribution of log10(SSD) values"
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
      poplot(plot, 12, 22)
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")
      ReportCalls <- AddPlot2Report()
      # Statistical test
      M <- median(unlist(PG[wNC, grep(topattern(SSD.Root), colnames(PG), value = TRUE)]))
      for (wh in wh1) { #wh <- wh1[1]
        grp <- SubCellFracAggr2$values[wh]
        SSDs[[grp]][[paste0(SSD.Pval.Root, grp)]] <- apply(SSDs[[grp]][, SSDkols], 1, function(x) {
          -log10(t.test(x, alternative = "greater", mu = M)$p.value)
        })
        temp <- FDR(SSDs[[grp]], grp, SSD.Pval.Root, fdr = BH.FDR, returns = c(TRUE, TRUE), method = "BH")
        SSDs[[grp]][, gsub("^Significant-", "Signif. SSDs-", colnames(temp$`Significance vector`))] <- temp$`Significance vector`
        SSD.FDR.thresh <- c(SSD.FDR.thresh, temp$Thresholds)
        PG[, colnames(SSDs[[grp]])] <- NA
        PG[wNC, colnames(SSDs[[grp]])] <- SSDs[[grp]]
      }
      temp <- PG[, grep("^Regulated - ", colnames(PG), value = TRUE, invert = TRUE)]
      subDr <- "Reg. analysis/Localisation"
      tempVP3 <- Volcano.plot(Prot = temp,
                              mode = "custom",
                              experiments.map = Exp.map,
                              X.root = paste0("Mean ", SSD.Root),
                              Y.root = SSD.Pval.Root,
                              aggregate.map = Aggregate.map,
                              aggregate.name = SubCellFracAggr2$aggregate,
                              aggregate.list = Aggregate.list, parameters = Param,
                              save = c("jpeg", "pdf"), labels = "FDR",
                              Ref.Ratio.values = RefSSDs,
                              Ref.Ratio.method = paste0("obs", RefRat_Mode),
                              ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                              FDR.thresh = SSD.FDR.thresh, FDR.root = "Signif. SSDs-FDR=",
                              arbitrary.lines = arbitrary.thr,
                              proteins = prot.list, proteins_split = protsplit,
                              return = TRUE, return.plot = TRUE,
                              title = "SSDs volcano plot ",
                              subfolder = subDr,
                              subfolderpertype = FALSE, Symmetrical = FALSE,
                              Alpha = "Rel. log10(Peptides count)",
                              Size = "Av. log10 abundance", Size.max = 2,
                              plotly = create_plotly, plotly_local = create_plotly_local,
                              plotly_labels = PrLabKol,
                              X.normalized = FALSE,
                              cl = parClust)
      if (!class(tempVP3) %in% c("try-error", "character")) {
        #
        VP_list <- tempVP3
        insrt <- ""
        Src <- paste0(libPath, "/extdata/R scripts/Sources/thresholds_Excel.R")
        #rstudioapi::documentOpen(Src)
        source(Src, local = FALSE)
        #
        g <- grep("Regulated - ", colnames(tempVP3$Protein_groups_file), value = TRUE)
        PG[, gsub("^Regulated - ", "Re-localized - ", g)] <- tempVP3$Protein_groups_file[,g]
        volcano.plots$Localisation_Unlabelled <- tempVP3$Plots$Unlabelled
        volcano.plots$Localisation_Labelled <- tempVP3$Plots$Labelled
        n2 <- names(volcano.plots$Localisation_Labelled)
        dir <- paste0(wd, "/Reg. analysis/Localisation")
        for (ttl in n2) {
          plot <- volcano.plots$Localisation_Labelled[[ttl]]
          ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
        }
        if ((create_plotly)&&(!create_plotly_local)) { plot_ly$"Localisation" <- tempVP3$"Plotly plots" }
        # Edit wording + create filters
        g <- grep("^Re-localized - ", colnames(PG), value = TRUE)
        for (gi in g) { #gi <- g[1]
          PG[which(PG[[gi]] == "too small FC"), gi] <- "unchanged distr."
          wUp <- grep("^up, FDR = ", PG[[gi]])
          if (length(wUp)) {
            PrtWidth <- 5
            dir <- paste0(wd, "/Reg. analysis/Localisation/Relocalised proteins")
            if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
            dirlist <- unique(c(dirlist, dir))
            plotNorm <- FALSE
            grp <- gsub("^Re-localized - ", "", gi)
            EM1 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == grp),]
            EM1$Replicate <- as.numeric(EM1$Replicate)
            EM1 <- EM1[order(EM1$Replicate, EM1$Compartment),]
            for (wup in wUp) { #wup <- wUp[1]
              # Protein heatmap
              Seq <- db$Sequence[match(unlist(strsplit(PG$`Leading protein IDs`[wup], ";"))[1], db$`Protein ID`)]
              Prfl <- data.frame(Sample = c(EM0$Ref.Sample.Aggregate, EM1$Ref.Sample.Aggregate))
              if (plotNorm) {
                PGRoot <- Prot.Expr.Root
                PepRoot <- pep.ref[length(pep.ref)]
              } else {
                PGRoot <- Prot.Expr.Root2
                PepRoot <- pep.ref2
              }
              Prfl$value <- 10^as.numeric(PG[wup, paste0(PGRoot, Prfl$Sample)])
              Prfl[, RSA$names] <- Isapply(strsplit(Prfl$Sample, "___"), unlist)
              Prfl$Replicate <- as.integer(Prfl$Replicate)
              Prfl[, c(SubCellFracAggr$aggregate, SubCellFracAggr2$aggregate)] <- Exp.map[match(Prfl$Sample, Exp.map$Ref.Sample.Aggregate),
                                                                                          c(SubCellFracAggr$column, SubCellFracAggr2$column)]
              Prfl$x <- 0
              Prfl$xend <- Prfl$x + Prfl$value/max(Prfl$value)*PrtWidth
              Prfl$y <- 0:(nrow(Prfl)-1)
              Prfl$"SubCell. Frac." <- factor(Prfl$Compartment, levels = Com)
              Prfl$Entity <- "Protein"
              Prfl$Angle <- 0
              Prfl$Size <- 3
              kol <- paste0(pep.ref[length(pep.ref)], Prfl$Sample)
              kol <- paste0(PepRoot, Prfl$Sample)
              pepPrfl <- pep[match(as.integer(unlist(strsplit(PG$`Peptide IDs`[wup], ";"))), pep$id),
                             c("Modified sequence", "Sequence", kol)]
              Match <- sapply(pepPrfl$Sequence, function(x) { nchar(unlist(strsplit(Seq, x))[1]) })
              pepPrfl <- pepPrfl[order(Match),]
              pepPrfl[, kol] <- sweep(pepPrfl[, kol], 1, rowMax(as.matrix(pepPrfl[, kol])), "/")
              mSeq <- unique(pepPrfl$"Modified sequence") 
              colnames(pepPrfl) <- gsub(topattern(PepRoot), "", colnames(pepPrfl))
              pepPrfl$Sequence <- NULL
              pepPrfl <- reshape2::melt(pepPrfl, id.vars = "Modified sequence")
              colnames(pepPrfl) <- c("Entity", "Sample", "value")
              pepPrfl$Sample <- as.character(pepPrfl$Sample)              
              pepPrfl[, RSA$names] <- Isapply(strsplit(pepPrfl$Sample, "___"), unlist)
              pepPrfl$Replicate <- as.integer(pepPrfl$Replicate)
              pepPrfl[, c(SubCellFracAggr$aggregate, SubCellFracAggr2$aggregate)] <- Exp.map[match(pepPrfl$Sample, Exp.map$Ref.Sample.Aggregate),
                                                                                             c(SubCellFracAggr$column, SubCellFracAggr2$column)]
              pepPrfl$x <- PrtWidth + match(pepPrfl$Entity, mSeq)
              pepPrfl$xend <- pepPrfl$x + pepPrfl$value
              pepPrfl$y <- Prfl$y[match(pepPrfl$Sample, Prfl$Sample)]
              pepPrfl$"SubCell. Frac." <- factor(pepPrfl$Compartment, levels = Com)
              pepPrfl$Angle <- 60
              pepPrfl$Size <- 2
              pepPrfl <- pepPrfl[which(pepPrfl$value > 0),]
              temp <- rbind(Prfl, pepPrfl)
              temp2 <- aggregate(temp[, c("x", "Angle", "Size")], list(temp$Entity), unique)
              colnames(temp2) <- c("Entity", "X", "Angle", "Size")
              temp3 <- aggregate(temp$y, list(temp[[SubCellFracAggr$aggregate]]), function(x) { list(Min = min(x),
                                                                                                     Max = max(x)+1) })
              colnames(temp3) <- c("Samples group", "x")
              temp3[, c("Min", "Max")] <- apply(temp3$x, 2, as.numeric)
              temp3$Mean <- apply(temp3[, c("Min", "Max")], 1, mean)
              temp3[, SubCellFracAggr$names] <- Isapply(strsplit(as.character(temp3$`Samples group`), "___"), unlist)
              temp3$`Sample group` <- cleanNms(temp3$`Samples group`)
              temp4 <- aggregate(temp$y, list(temp[[SubCellFracAggr2$aggregate]]), function(x) { list(Min = min(x),
                                                                                                      Max = max(x)+1) })
              colnames(temp4) <- c("Samples group 2", "x")
              temp4[, c("Min", "Max")] <- apply(temp4$x, 2, as.numeric)
              temp4$Mean <- apply(temp4[, c("Min", "Max")], 1, mean)
              temp4$"Samples group 2" <- cleanNms(temp4$"Samples group 2")
              temp5 <- aggregate(temp$y, list(temp$Sample), function(x) { list(Min = min(x), Max = max(x)+1) })
              colnames(temp5) <- c("Sample", "x")
              temp5[, c("Min", "Max")] <- apply(temp5$x, 2, as.numeric)
              temp5$Mean <- apply(temp5[, c("Min", "Max")], 1, mean)
              temp5$"SubCell. frac." <- Exp.map$Compartment[match(temp5$Sample, Exp.map$Ref.Sample.Aggregate)]
              NCSc <- max(nchar(temp5$`SubCell. frac.`))
              ttl <- paste0("Subcell. fract. distr. - ", PG$`Leading protein IDs`[wup])
              plot <- ggplot(temp) +
                geom_rect(aes(xmin = x, xmax = xend, ymin = y, ymax = y+1, fill = `SubCell. Frac.`)) +
                scale_fill_viridis_d(begin = 0.25) +
                geom_segment(aes(x = x, xend = x), y = -2, yend = max(temp$y + 1.5), colour = "lightgrey") +
                geom_segment(x = max(temp$x)+1, xend =  max(temp$x)+1, y = -2, yend = max(temp$y + 1.5), colour = "lightgrey") +
                geom_segment(x = PrtWidth + 0.5, xend = PrtWidth + 0.5, y = -2, yend = max(temp$y + 1.5)) +
                geom_hline(yintercept = max(temp$y + 1.5)) +
                geom_text(data = temp2, aes(x = X, label = Entity, angle = Angle, size = Size), y = max(temp$y + 2.5), hjust = 0) +
                geom_text(data = temp4, aes(y = Mean, label = `Samples group 2`), x = -3-NCSc, hjust = 1, vjust = 0.5, size = 2.5) +
                geom_segment(data = temp4, aes(y = Min+0.1, yend = Max-0.1), x = -2.8-NCSc, xend = -2.8-NCSc) +
                geom_text(data = temp3, aes(y = Mean, label = Replicate), x = -2-NCSc, hjust = 1, vjust = 0.5, size = 1.8) +
                geom_segment(data = temp3, aes(y = Min+0.1, yend = Max-0.1), x = -1.8-NCSc, xend = -1.8-NCSc) +
                geom_text(data = temp5, aes(y = Mean, label = `SubCell. frac.`), x = -1, hjust = 1, vjust = 0.5, size = 1.8) +
                coord_fixed() + theme_bw() + scale_size_identity() +
                ggtitle("Distribution across subcellular fractions", subtitle = PG$Label[wup]) +
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(),  axis.line = element_line(colour = "black"),
                      plot.margin = margin(0, 0, 0, 0, "cm"), legend.position = "none") +
                xlab(NULL) + ylab(NULL) + #theme(legend.position = "none") +
                xlim(-5-NCSc, max(temp$xend+1)) + ylim(0, max(temp$y+20))
              #poplot(plot, 12, 22)
              ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 20, height = 12, units = "in")
              ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 20, height = 12, units = "in")
            }
          }
        }
        g1 <- gsub("^Re-localized - ", "", g)
        up <- grep("^reloc\\., FDR = ", unique(unlist(PG[, g])), value = TRUE)
        Reg_filters$Localisation <- list()
        if ("con" %in% filter_types) {
          Reg_filters$Localisation$"By condition" <- list()
          rat <- paste0("Mean ", SSD.Root, g1)
          for (i in 1:length(g)) { #i <- 1
            Reg_filters$Localisation$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                                     Filter_up = sort(which(PG[[g[i]]] %in% up)),
                                                                     Filter_down = c(),
                                                                     Filter = sort(which(PG[[g[i]]] %in% up)),
                                                                     Ratios = PG[[rat[i]]],
                                                                     Background_filter = 1:nrow(PG))
            # To do here: check that Ratios values above are correct
          }
        }
      } else { warning("No localisation volcano plots created, investigate!") }
    } else {
      if (!length(wh0)) {
        warning("Skipping re-localisation analysis: no reference group...")
      }
      if (length(wh0) > 0) {
        warning("Skipping re-localisation analysis: too many reference groups, only 1 allowed...")
      }
      if (!length(wh1)) {
        warning("Skipping re-localisation analysis: no non-reference group...")
      }
    }
  } else {
    warning("Skipping re-localisation analysis: some subcellular fraction groups contain both reference and non-reference samples...")
  }
  ReportCalls <- AddSpace2Report()
  rm(list = ls()[which(!ls() %in% .obj)])
  Script <- readLines(ScriptPath)
  gc()
  parLapply(parClust, 1:N.clust, function(x) {
    rm(list = ls())
    gc()
  })
  saveImgFun(BckUpFl)
  #loadFun(BckUpFl)
  source(parSrc, local = FALSE)
}

### Check that Cytoscape is installed and can run, then launch it.
Src <- paste0(libPath, "/extdata/R scripts/Sources/Cytoscape_init.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

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
  # Initialize ClueGO
  Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_init.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
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
            if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) {
              loadFun("GO_terms.RData")
            }
            #
            Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_enrich.R")
            #rstudioapi::documentOpen(Src)
            source(Src, local = FALSE)
            #
            clueGO_outDir <- dir
            clueGO_type <- "Enrichment (Right-sided hypergeometric test)"
            Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_enrich.R")
            #rstudioapi::documentOpen(Src)
            source(Src, local = FALSE)
            #
            # Cleanup - do it now, not within sources!
            for (i in allArgs) { try(rm(i), silent = TRUE) }
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
                temp$Offspring <- sapply(temp$Offspring, paste, collapse = ";")
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
              w <- which(sapply(colnames(temp), function(x) { class(temp[[x]]) }) == "character")
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
                  } else { w <- which(sapply(Reg_GO_terms[[tstbee]][,Kol2], function(x) {"+" %in% x})) }
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
                    temp[i, 2:(N+1)] <- sapply(Kol3, function(x) {
                      sum((GO_Plots[[tstbee]]$GO_terms[[x]] == "+")&(GO_Plots[[tstbee]]$GO_terms[[Kol3[i]]] == "+"),
                          na.rm = TRUE)
                    })
                  }
                  names(W) <- cleanNms(gsub(" [0-9]+%$", "", names(W)))
                  tst <- sapply(strsplit(names(W), " - "), length)
                  tst <- (min(tst) > 1)&(length(unique(tst)) == 1)
                  if (tst) {
                    tst <- as.data.frame(t(sapply(strsplit(names(W), " - "), unlist)))
                    l <- apply(tst, 2, function(x) { length(unique(x)) })
                    tst <- tst[,which(l > 1)]
                    names(W) <- apply(tst, 1, paste, collapse = " - ")
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
    Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_enrich.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    clueGO_outDir <- dir
    clueGO_type <- "Enrichment/Depletion (Two-sided hypergeometric test)"
    Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_enrich.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    # Cleanup - do it now, not within sources!
    for (i in allArgs) { try(rm(i), silent = TRUE) }
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
                           c("", "and ClueGO")[CytoScape+1], ".")
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
parLapply(parClust, 1:N.clust, function(x) {
  rm(list = ls())
  gc()
})
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)

#### Code chunk - Modified peptides analysis
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

#### Code chunk - Proteome ruler
if (protrul) {
  ProtRulRoot %<o% "log10(est. copies/cell) - "
  temp <- PG
  exprsRt <- paste0("Mean ", prtRfRoot)
  if (LocAnalysis) {
    temp <- temp[, grep(topattern(exprsRt), colnames(temp), invert = TRUE)]
    for (grp2 in SubCellFracAggr2$values) { #grp2 <- SubCellFracAggr2$values[1]
      em2 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == grp2),]
      for (grp in unique(em2[[SubCellFracAggr$column]])) { #grp <- unique(em2[[SubCellFracAggr$column]])[1]
        em <- em2[which(em2[[SubCellFracAggr$column]] == grp),]
        kol <- paste0(prtRfRoot, em$Ref.Sample.Aggregate)
        temp[[paste0(prtRfRoot, grp)]] <- apply(10^temp[, kol], 1, function(x) {
          log10(sum(is.all.good(x)))
        })
      }
      temp[[paste0(exprsRt, grp2)]] <- apply(temp[, paste0(prtRfRoot, unique(em2[[SubCellFracAggr$column]]))], 1, function(x) {
        mean(is.all.good(x))
      })
    }
  }
  temp <- try(Prot.Ruler(temp, db, exprsRt, NuclL = ProtRulNuclL), silent = TRUE)
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
  for (n in 1:length(plot_ly)) { #n <- 1
    i <- plot_ly[[n]]
    if (length(i)) {
      plot_ly_addresses[kount,] <- c(names(plot_ly)[n], "", "")
      for (j in 1:length(i)) { #j <- 1
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
                            raw.files = rawFiles, sc = sc, cl = parClust, MQtxt = indir)
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
if (length(GO_PG_col)) {
  if (!exists("GO_terms")) { loadFun("GO_terms.RData") }
  GO_PG_col <- GO_PG_col[which(GO_PG_col %in% GO_terms$ID)]
  if (length(GO_PG_col)) {
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
    w <- which(sapply(GO_PG_col, function(x) { sum(Offspring[[x]] %in% tmp$Term) }) == 0)
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
      tst <- setNames(sapply(GO_PG_col2, function(x) {
        sum(PG[[x]] == "+")
      }), GO_PG_col2)
      tst <- paste0(GO_PG_col2, " -> ", tst, collapse = "\n - ")
      tst <- paste0("Number of protein groups per GO term of interest\n - ", tst, "\n\n")
      cat(tst)
      write(tst, paste0(wd, "/Reg. analysis/GO enrich/Dataset/GO_terms_of_interest.txt"))
    }
  }
}

#### Code chunk - Create output tables
## PSMs
dir <- paste0(wd, "/Tables")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
w <- which(sapply(colnames(ev), function(x) { "list" %in% class(ev[[x]]) }))
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
  tst$L <- sapply(tst$x, length)
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
    w <- which(sapply(intColsTbl, nrow) > 0)
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
    w <- which(sapply(ratColsTbl, nrow) > 0)
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
    .obj <- unique(c(.obj, PepColList)) # Here easier than using a custom operator
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
    ColumnsTbl <- ColumnsTbl[which(sapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }) > 0)]
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
  KolNames <- gsub("Evidences?", c("PSMs", "Ev.")[(SearchSoft == "MAXQUANT")+1], KolNames)
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
  tst$L <- sapply(tst$x, length)
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
w <- which(sapply(intColsTbl, nrow) > 0)
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
w <- which(sapply(ratColsTbl, nrow) > 0)
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
  loadFun("GO_terms.RData") # Re-load object
  GOCC <- GO_terms$ID[which(GO_terms$Ontology == "CC")]
  tempData$"GO-ID (CC)" <- lapply(strsplit(tempData$`GO-ID`, ";"), function(x) { x[which(x %in% GOCC)] })
  w <- which(sapply(tempData$"GO-ID (CC)", length) > 0)
  tempData$"GO (CC)" <- ""
  tempData$"GO (CC)"[w] <- sapply(tempData$"GO-ID (CC)"[w], function(x) {
    paste(GO_terms$Term[match(x, GO_terms$ID)], collapse = ";")
  })
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
ColumnsTbl <- ColumnsTbl[which(sapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }) > 0)]
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
    tst$L <- sapply(tst$x, length)
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
  w <- which(sapply(ratColsTbl, nrow) > 0)
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
  ColumnsTbl <- ColumnsTbl[which(sapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }) > 0)]
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
  tmp$Razor <- sapply(strsplit(tmp$Razor, ";"), function(x) {
    as.logical(toupper(x))
  })
  tmp$IDs <- sapply(strsplit(tmp$IDs, ";"), as.numeric)
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
# cleanNms2 function specifically designed to clean names for Venn diagrams
cleanNms2 %<o% function(names, groups = VPAL, sep = "\n/vs/\n", simplify) {
  isList <- ("list" %in% class(names))
  if (missing(simplify)) { simplify <- !isList }
  if (!isList) {
    names <- sapply(strsplit(gsub("^\\(|\\)$", "", names), "\\) - \\("), unlist)
  }
  uNms <- unique(unlist(names))
  nuNms <- as.data.frame(t(as.data.frame(strsplit(uNms, "___"))))
  colnames(nuNms) <- groups$names
  nuNms$Full <- uNms
  w <- which(sapply(groups$names, function(x) { length(unique(nuNms[[x]])) }) > 1)
  nuNms$New <- do.call(paste, c(nuNms[, groups$names[w], drop = FALSE], sep = ""))
  nuNames <- lapply(names, function(x) { nuNms$New[match(x, nuNms$Full) ]})
  if (simplify) { nuNames <- sapply(nuNames, paste, collapse = sep) }
  return(nuNames)
}
#
msg <- "Venn diagrams"
ReportCalls <- AddMsg2Report(Space = FALSE)
HdrStlVenn <- createStyle(textDecoration = "bold", halign = "left", valign = "bottom", wrapText = TRUE,
                          numFmt = "TEXT", fontSize = 12, textRotation = 60)
II <- setNames(1, "Protein groups")
if ((exists("PTMs_pep"))&&(length(PTMs_pep))) {
  Mod2Venn <- names(PTMs_pep)
  II[paste0(Mod2Venn, "-mod. pept.")] <- 1+(seq_along(length(Mod2Venn)))
}
VennMx <- 7
for (ii in II) { #ii <- II[1] #ii <- II[2]
  if (ii == 1) {
    dir <- paste0(wd, "/Venn diagrams")
    ttest_Filt <- Reg_filters$"t-tests"$"By condition"
    vennRoot <- ""
    myData <- PG
    myRef <- paste0("Mean ", prtRfRoot)
    infoKol <- "Genes"
    idKol <- "Leading protein IDs"
    if (F.test) {
      Ftest_Filt <- Reg_filters$"F-tests"$"By condition"
      myFData <- F_test_data
    }
  } else {
    Ptm <- Mod2Venn[ii-1]
    dir <- paste0(wd, "/Reg. analysis/", Ptm, "/Venn diagrams")
    ttest_Filt <- PTMs_Reg_filters[[Ptm]]$"t-tests"$"By condition"
    vennRoot <- paste0(Ptm, " ")
    myData <- PTMs_pep[[Ptm]]
    myRef <- grep("^Mean ([^ ]+ )?log10\\(", colnames(myData), value = TRUE)
    myRef <- paste0(unique(gsub(" - .*", "", myRef)), " - ")
    myRef <- myRef[length(myRef)]
    infoKol <- "Proteins"
    idKol <- "Modified sequence"
    if (F.test) {
      Ftest_Filt <- PTMs_Reg_filters[[Ptm]]$"F-tests"$"By condition"
      myFData <- PTMs_F_test_data[[Ptm]]
    }
  }
  topTitle <- paste0("Venn diagram - ", names(II)[ii])
  infoKol <- infoKol[which(infoKol %in% colnames(myData))]
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  #
  wb <- createWorkbook()
  wbKount <- 0
  #
  dir2 <- paste0(dir, "/Observations")
  if (!dir.exists(dir2)) { dir.create(dir2, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir2))
  setwd(dir2)
  #
  comp_list <- setNames(lapply(VPAL$values, function(grp) {
    x <- myData[[paste0(myRef, grp)]]
    which(is.all.good(as.numeric(x), 2))
  }), cleanNms2(VPAL$values))
  w <- which(sapply(comp_list, length) > 0)
  VennExp <- names(comp_list)[w]
  OK <- length(w) > 1
  ttl <- paste0(vennRoot, "LFQ Venn diagram, all")
  if (length(w) > VennMx) {
    msg <- paste0("Too many groups, select at least 2 and up to ", VennMx,
                  " to include in Venn diagram ", ttl)
    opt <- sapply(VennExp, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") })
    VennExp <- VennExp[match(dlg_list(opt, opt[1:VennMx], TRUE, msg)$res, opt)]
    if (length(VennExp) == 1) {
      OK <- FALSE
      if (!is.na(VennExp)) {
        warning("Skipping per-sample-group observations Venn diagrams: you should have selected at least 2 groups!")
      }
    }
    if ((length(VennExp) > VennMx)&&(length(VennExp) < 1)) {
      msg <- paste0("Skipping per-sample-group observations Venn diagrams: you should have selected ", VennMx, " groups at most!")
      warning(msg)
      OK <- FALSE
    }
  }
  if (OK) {
    cat("Creating per-sample-group observations Venn diagrams...\n")
    comp_list <- comp_list[VennExp]
    plot <- venn(comp_list, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
    plot <- plot + ggtitle(topTitle, subtitle = "Global, LFQ") +
      theme(plot.title = element_text(size = 15),
            plot.subtitle = element_text(size = 10))
    poplot(plot)
    ggsave(paste0(dir2, "/", ttl, ".jpg"), plot, dpi = 150)
    ggsave(paste0(dir2, "/", ttl, ".pdf"), plot, dpi = 150)
    #system(paste0("open \"", dir2, "/", ttl, ".jpg", "\""))
    ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_img(Report, \"", dir2, "/", ttl, ".jpg\", height = 6, width = 6)"))
    wbKount <- wbKount+1
    SheetNm <- "Sample groups composition"
    addWorksheet(wb, SheetNm)
    writeData(wb, SheetNm, myData[, c("id", idKol, infoKol)], 1, 1)
    l <- length(comp_list)
    tmp <- sapply(names(comp_list), function(grp) {
      res <- rep("", nrow(myData))
      res[comp_list[[grp]]] <- "+"
      return(res)
    })
    writeData(wb, SheetNm, tmp, 4, 1)
    setRowHeights(wb, SheetNm, 1, 120)
    addStyle(wb, SheetNm, HdrStlVenn, 1, 1:(l+3))
  } else {
    msg <- paste0(ttl, ": not enough groups to compare!")
    ReportCalls <- AddMsg2Report(Space = FALSE)
  }
  #
  dir2 <- paste0(dir, "/Stat. tests")
  if (!dir.exists(dir2)) { dir.create(dir2, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir2))
  setwd(dir2)
  #
  vennlev <- c("up", "down")[1:(1+TwoSided)]
  if (("Venn.Groups" %in% colnames(Param))&&(Param$Venn.Groups != "")) { VennGrp <- Param$Venn.Groups } else { VennGrp <- "GLOBAL" }
  if (toupper(VennGrp) == "GLOBAL") { VennGrp <- "GLOBAL" }
  nmz <- names(ttest_Filt)
  if (length(nmz) > 1) {
    msg <- " - t-tests"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    for (r in vennlev) { #r <- "up"
      ReportCalls <- AddMsg2Report(Msg = r, Offset = TRUE, Space = FALSE)
      comp_list <- setNames(lapply(nmz, function(x) {
        ttest_Filt[[x]][[paste0("Filter_", r)]]
      }), cleanNms2(nmz))
      w <- which(sapply(comp_list, length) > 0)
      VennExp <- names(comp_list)[w]
      OK <- length(w) > 1
      ttl <- gsub("^g", "G", paste0(vennRoot, "t-tests Venn diagram, all, ", r))
      if (length(w) > VennMx) {
        msg <- paste0("Too many groups, select at least 2 and up to ", VennMx,
                      " to include in Venn diagram ", ttl)
        opt <- sapply(VennExp, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") })
        VennExp <- VennExp[match(dlg_list(opt, opt[1:VennMx], TRUE, msg)$res, opt)]
        if (length(VennExp) == 1) {
          OK <- FALSE
          if (!is.na(VennExp)) {
            warning(ttl, ": skipping per-sample-group observations. You should have selected at least 2 groups!")
          }
        }
        if ((length(VennExp) > VennMx)&&(length(VennExp) < 1)) {
          msg <- paste0(ttl, ": skipping per-sample-group observations. You should have selected ", VennMx, " groups at most!")
          warning(msg)
          OK <- FALSE
        }
      }
      if (OK) {
        comp_list <- comp_list[VennExp]
        kr <- paste0("Regulated - ", nmz)
        tmp <- myData[unique(unlist(comp_list)), c(idKol, infoKol, "id", kr)]
        if (r == "up") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
        if (r == "down") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
        tst <- apply(tmp[, kr], 1, function(x) { length(w[which(x %in% good)]) })
        tmp <- tmp[order(tst, decreasing = TRUE),]
        write.csv(tmp, paste0(dir2, "/", ttl, " - table.csv"), row.names = FALSE)
        plot <- venn(comp_list, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
        plot <- plot +
          ggtitle(topTitle, subtitle = paste0("t-tests, ", r)) +
          theme(plot.title = element_text(size = 15),
                plot.subtitle = element_text(size = 10))
        poplot(plot)
        ggsave(paste0(dir2, "/", ttl, ".jpg"), plot, dpi = 150)
        ggsave(paste0(dir2, "/", ttl, ".pdf"), plot, dpi = 150)
        #system(paste0("open \"", dir2, "/", ttl, ".jpg", "\""))
        ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_img(Report, \"", dir2, "/", ttl, ".jpg\", height = 6, width = 6)"))
        wbKount <- wbKount+1
        SheetNm <- paste0("t-test ", r)
        addWorksheet(wb, SheetNm)
        writeData(wb, SheetNm, myData[, c("id", idKol, infoKol)], 1, 1)
        l <- length(comp_list)
        tmp <- sapply(names(comp_list), function(grp) {
          res <- rep("", nrow(myData))
          res[comp_list[[grp]]] <- "+"
          return(res)
        })
        writeData(wb, SheetNm, tmp, 4, 1)
        setRowHeights(wb, SheetNm, 1, 120)
        addStyle(wb, SheetNm, HdrStlVenn, 1, 1:(l+3))
      } else {
        msg <- paste0(ttl, ": not enough groups with regulated proteins to compare!")
        ReportCalls <- AddMsg2Report(Space = FALSE)
      }
    }
    ReportCalls <- AddSpace2Report()
  }
  VennGrp2 <- parse.Param.aggreg(Param_filter(Param$Venn.Groups, "Rep"))
  if ((VennGrp != "GLOBAL")&&(length(VennGrp2$values) > 1)) {
    for (i in VennGrp2$values) { #i <- VennGrp2$values[1]
      nms <- nmz[which(nmz %in% Exp.map[which(Exp.map[[VennGrp2$column]] == i), VPAL$column])]
      if (length(nms) > 1) {
        j <- cleanNms(i)
        msg <- paste0(" - ", j)
        ReportCalls <- AddMsg2Report(Space = FALSE)
        for (r in vennlev) { #r <- "up"
          ReportCalls <- AddMsg2Report(Msg = r, Offset = TRUE, Space = FALSE)
          ttl <- paste0(vennRoot, "t-tests Venn diagram, ", j, ", ", r)
          comp_list <- setNames(lapply(nms, function(x) {
            ttest_Filt[[x]][[paste0("Filter_", r)]]
          }), cleanNms2(nms))
          w <- which(sapply(comp_list, length) > 0)
          VennExp <- names(comp_list)[w]
          OK <- length(w) > 1
          ttl <- gsub("^g", "G", paste0(vennRoot, "t-tests Venn diagram, all, ", r))
          if (length(w) > VennMx) {
            msg <- paste0("Too many groups, select at least 2 and up to ", VennMx,
                          " to include in Venn diagram ", ttl)
            opt <- sapply(VennExp, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") })
            VennExp <- VennExp[match(dlg_list(opt, opt[1:VennMx], TRUE, msg)$res, opt)]
            if (length(VennExp) == 1) {
              OK <- FALSE
              if (!is.na(VennExp)) {
                warning(ttl, ": skipping per-sample-group observations. You should have selected at least 2 groups!")
              }
            }
            if ((length(VennExp) > VennMx)&&(length(VennExp) < 1)) {
              msg <- paste0(ttl, ": skipping per-sample-group observations. You should have selected ", VennMx, " groups at most!")
              warning(msg)
              OK <- FALSE
            }
          }
          if (OK) {
            comp_list <- comp_list[VennExp]
            kr <- paste0("Regulated - ", nms)
            tmp <- myData[unique(unlist(comp_list)), c(idKol, infoKol, "id", kr)]
            if (r == "up") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
            if (r == "down") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
            tst <- apply(tmp[, kr], 1, function(x) { length(w[which(x %in% good)]) })
            tmp <- tmp[order(tst, decreasing = TRUE),]
            write.csv(tmp, paste0(dir2, "/", ttl, " - table.csv"), row.names = FALSE)
            names(comp_list) <- cleanNms(names(comp_list))
            plot <- venn(comp_list, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
            plot <- plot +
              ggtitle(topTitle, subtitle = paste0("t-tests, subgroup = ", j, ", ", r)) +
              theme(plot.title = element_text(size = 15),
                    plot.subtitle = element_text(size = 10))
            poplot(plot)
            ggsave(paste0(dir2, "/", ttl, ".jpg"), plot, dpi = 150)
            ggsave(paste0(dir2, "/", ttl, ".pdf"), plot, dpi = 150)
            #system(paste0("open \"", dir2, "/", ttl, ".jpg", "\""))
            ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_img(Report, \"", dir2, "/", ttl, ".jpg\", height = 6, width = 6)"))
            wbKount <- wbKount+1
            SheetNm <- paste0("t-test ", j, " ", r)
            addWorksheet(wb, SheetNm)
            writeData(wb, SheetNm, myData[, c("id", idKol, infoKol)], 1, 1)
            l <- length(comp_list)
            tmp <- sapply(names(comp_list), function(grp) {
              res <- rep("", nrow(myData))
              res[comp_list[[grp]]] <- "+"
              return(res)
            })
            writeData(wb, SheetNm, tmp, 4, 1)
            setRowHeights(wb, SheetNm, 1, 120)
            addStyle(wb, SheetNm, HdrStlVenn, 1, 1:(l+3))
          } else {
            msg <- paste0(ttl, ": not enough groups with regulated proteins to compare!")
            ReportCalls <- AddMsg2Report(Space = FALSE)
          }
        }
        ReportCalls <- AddSpace2Report()
      } else {
        msg <- paste0("t-tests Venn diagrams: ", i, " - could not draw: not enough groups to compare!")
        ReportCalls <- AddMsg2Report()
      }
    }
  }
  if (F.test) {
    dir2 <- paste0(dir, "/Stat. tests")
    if (!dir.exists(dir2)) { dir.create(dir2, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir2))
    setwd(dir2)
    #
    nmz <- names(Ftest_Filt)
    nmz <- nmz[which(!grepl(" VS ", nmz))]
    if (length(nmz) > 1) {
      msg <- " - F-tests"
      ReportCalls <- AddMsg2Report(Space = FALSE)
      for (r in vennlev) { #r <- "up"
        ReportCalls <- AddMsg2Report(Msg = r, Offset = TRUE, Space = FALSE)
        ttl <- paste0(vennRoot, "Venn diagram - F-test, global - ", r)
        nmz2 <- cleanNms2(nmz)
        comp_list <- setNames(lapply(nmz, function(x) {
          Ftest_Filt[[x]][[paste0("Filter_", r)]]
        }), nmz2)
        tst <- sapply(comp_list, length)
        w <- which(tst > 0)
        if (length(w) > 1) {
          if (length(w) <= VennMx) {
            comp_list <- comp_list[w]
            w <- unique(unlist(comp_list))
            wN <- which(!infoKol %in% colnames(myFData))
            if (length(wN)) {
              myFData[, infoKol[wN]] <- myData[match(myFData[[idKol]], myData[[idKol]]), infoKol]
            }
            tmp <- myFData[w, c(idKol, infoKol)]
            tmp$id <- myData$id[w]
            kr <- paste0("mod. F-test Regulated - ", nmz)
            tmp[, kr] <- myFData[w, kr]
            if (r == "up") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
            if (r == "down") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
            tst <- apply(tmp[, kr], 1, function(x) { length(w[which(x %in% good)]) })
            tmp <- tmp[order(tst, decreasing = TRUE),]
            write.csv(tmp, paste0(dir2, "/", ttl, " - table.csv"), row.names = FALSE)
            names(comp_list) <- cleanNms(names(comp_list))
            plot <- venn(comp_list, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
            plot <- plot +
              ggtitle(topTitle, subtitle = paste0(vennRoot, "F-test, ", r)) +
              theme(plot.title = element_text(size = 15),
                    plot.subtitle = element_text(size = 10))
            poplot(plot)
            ggsave(paste0(dir2, "/", ttl, ".jpg"), plot, dpi = 150)
            ggsave(paste0(dir2, "/", ttl, ".pdf"), plot, dpi = 150)
            #system(paste0("open \"", dir2, "/", ttl, ".jpg", "\""))
            ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_img(Report, \"", dir2, "/", ttl, ".jpg\", height = 6, width = 6)"))
            wbKount <- wbKount+1
            SheetNm <- paste0("F-test global ", r)
            addWorksheet(wb, SheetNm)
            writeData(wb, SheetNm, myData[, c("id", idKol, infoKol)], 1, 1)
            l <- length(comp_list)
            tmp <- sapply(names(comp_list), function(grp) {
              res <- rep("", nrow(myData))
              res[comp_list[[grp]]] <- "+"
              return(res)
            })
            writeData(wb, SheetNm, tmp, 4, 1)
            setRowHeights(wb, SheetNm, 1, 120)
            addStyle(wb, SheetNm, HdrStlVenn, 1, 1:(l+3))
          } else {
            msg <- paste0("Could not draw global F-tests Venn diagram: more than ", VennMx, " groups to compare!")
            ReportCalls <- AddMsg2Report(Space = FALSE)
          }
          ReportCalls <- AddSpace2Report()
        } else {
          msg <- "Could not draw global F-tests Venn diagram: not enough groups with regulated proteins to compare!"
          ReportCalls <- AddMsg2Report()
        }
      }
    } else {
      msg <- "Could not draw global F-tests Venn diagram: not enough groups to compare!"
      ReportCalls <- AddMsg2Report()
    }
  }
  if (wbKount) { saveWorkbook(wb, paste0(dir2, "/Venn diagrams.xlsx"), overwrite = TRUE) }
  setwd(wd)
  ReportCalls <- AddSpace2Report()
}

#### Code chunk - Coverage maps, XICs and heatmaps for proteins of interest
protlspep <- prot.list_pep
if (length(protlspep)) {
  test <- sapply(protlspep, function(i) { length(grsep2(i, PG$"Leading protein IDs")) })
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
  w <- which(sapply(tst, length) > 0)
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
          parLapply(parClust, u, function(sq) { #sq <- u[1] #sq <- u[2]
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
            ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".jpeg"), plot, dpi = 450, height = 10, width = 10)
            ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".pdf"), plot, height = 10, width = 10)
          })
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
      if (SearchSoft == "MAXQUANT") {
        if ((length(MQFold)==1)&&(dir.exists(MQFold))) {
          modFls <- paste0(MQFold, "/bin/conf/modifications", c("", ".local"), ".xml")
          modFls <- modFls[which(file.exists(modFls))]
        } else {
          if ((exists("mqFld"))&&(dir.exists(mqFld))) { dflt <- mqFld } else { dflt <- "C:" }
          dflt <- paste0(dflt, "/*.xml")
          modFls <- choose.files(dflt, "Select MaxQuant modifications file(s) as source of PTMs mass shifts:")
          if ((length(modFls) > 1)||(!is.na(modFls))) { mqFld <- unique(dirname(modFls)) }
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
            w <- which(sapply(colnames(Modifs), function(x) { class(Modifs[[x]]) }) == "list")
            for (i in w) { temp[[i]] <- sapply(temp[[i]],  paste, collapse = ", ") }
            write.csv(temp, "Workflow control/Modifications.csv", row.names = FALSE)
          }
          Modifs$"Mass shift" <- sapply(strsplit(Modifs$Composition, " "), function(x) {
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
          })
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
          tmp <- sapply(tmp, function(x) { #x <- tmp[1]
            x <- unlist(x)
            w <- grep("\\[.+\\]", x)
            x[w] <- paste0("[", round(Modifs$"Mass shift"[match(x[w], paste0("[", Modifs$Mark, "]"))], 0), "]")
            return(paste(x, collapse = ""))
          })
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
#
#
# Extend this to modified peptide filters!!! 
#
#
# NB: gets Ratios from the GO enrichment part of the script currently
# If you want to run this, you will thus first need to run the GO enrichment section first
setwd(wd)
# Packages for STRINGdb graphs
packs <- c("rbioapi", "png")
for (pack in packs) { require(pack, character.only = TRUE) }
# Create STRINGdb directory
GraphTypes %<o% c("Functional", "Physical")
dirs <- paste0(wd, "/STRINGdb/", GraphTypes)
for (dir in dirs) { if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) } }
dirlist <- unique(c(dirlist, dir))
# Process filters for STRINGdb
intNets <- list()
Tsts <- c("t-tests", "F-tests", "Localisation", "SAINTexpress")
WhTsts <- which(Tsts %in% names(Reg_filters))
msg <- "Starting STRINGdb network analysis"
ReportCalls <- AddMsg2Report(Space = FALSE)
tmpPG <- listMelt(strsplit(PG$"Leading protein IDs", ";"), PG$id)
tmpPG$L1 <- as.integer(tmpPG$L1)
IDs <- unique(unlist(lapply(names(Reg_filters), function(x) {
  i <- (x == "SAINTexpress")+1
  lapply(Reg_filters[[x]]$`By condition`, function(y) {
    w <- y$Filter
    dat <- get(c("PG", "allSAINTs")[i])
    w <- w[which(dat$"Potential contaminant"[w] != "+")]
    y <- dat[[c("Leading protein IDs", "Protein")[i]]][w]
    if (i == 1) { y <- unlist(strsplit(y, ";")) }
    return(y)
  })
})))
w <- which(db$`Protein ID` %in% IDs)
allTaxIDs <- unique(db$TaxID[w])
tmpPG <- tmpPG[which(tmpPG$value %in% IDs),]
source(parSrc, local = FALSE)
allProteins_mapped <- try(setNames(lapply(allTaxIDs, function(txid) { #txid <- allTaxIDs[1]
  kol <- c("Protein ID", "Common Name", "TAIR")
  kol <- kol[which(kol %in% colnames(db))]
  tmpDB <- db[which((db$`Protein ID` %in% IDs)&(db$TaxID == txid)), kol]
  x <- gsub("^cRAP[0-9]{3}", "", gsub("^CON__", "", tmpDB$`Protein ID`))
  y <- split(x, ceiling(seq_along(x)/20))
  n <- length(y)
  clusterExport(parClust, c("txid", "y", "tmpDB"), envir = environment())
  res <- parLapply(parClust, 1:n, function(i) {
    ids <- y[[i]]
    #a <- try(rbioapi::rba_string_map_ids(ids, txid), silent = TRUE)
    #if ("try-error" %in% class(a)) {
      try({
        rqst <- paste0("https://string-db.org/api/tsv/get_string_ids?identifiers=", paste(ids, collapse = "%0d"),
                       "&species=", txid)
        response <- httr::PUT(rqst, encode = "json")
        if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
        a <- httr::content(response, "parsed", show_col_types = FALSE)
        a$queryItem <- ids[a$queryIndex+1]
        m <- match(a$queryItem, tmpDB$`Protein ID`)
        a$Name <- tmpDB$`Common Name`[m]
        if ("TAIR" %in% colnames(tmpDB)) { # (For our dear plant-people)
          a$TAIR <- tmpDB$TAIR[m]
        }
      }, silent = TRUE)
    #}
    rs <- list(Outcome = (!"try-error" %in% class(a)))
    if (rs$Outcome) { rs$Result <- a }
    return(rs)
  })
  w <- which(sapply(res, function(rs) { rs$Outcome&&("data.frame" %in% class(rs$Result)) }))
  res <- res[w]
  if (length(res)) {
    res <- lapply(res, function(rs) { rs$Result })
    res <- plyr::rbind.fill(res)
  } else {
    res <- data.frame(Error = c())
  }
  return(res)
}), paste0("TaxID_", allTaxIDs)), silent = TRUE)
w <- which(!is.null(sapply(allProteins_mapped, nrow)))
allProteins_mapped <- allProteins_mapped[w]
if (4 %in% WhTsts) {
  tmpMap <- listMelt(strsplit(PG$`Leading protein IDs`, ";"), PG$id)
}
# See https://string-db.org/help/api/ for API help
if (length(WhTsts)&&length(allProteins_mapped)) {
  filtersDF <- lapply(WhTsts, function(tt) { #tt <- 1 #tt <- 2 #tt <- 3 #tt <- 4
    Filt <- Reg_filters[[Tsts[tt]]]$"By condition"
    nms <- sort(names(Filt))
    filtDF <- data.frame(Group = Tsts[tt],
                         Name = nms,
                         Type = c(1, 2)[(tt == 4)+1])
    filtDF$Filter <- lapply(filtDF$Name, function(nm) {
      Filt[[nm]]
    })
    filtDF$Reg <- filtDF$TaxID <- c()
    w <- which(sapply(filtDF$Filter, function(flt) {
      (sum(c("Ratios", "Filter") %in% names(flt)) == 2)
    }))
    filtDF <- filtDF[w,]
    nr <- nrow(filtDF)
    if (nr) {
      filtDF$W <- lapply(1:nr, function(x) { #x <- 1
        w <- filtDF$Filter[[x]]$Filter
        typ <- filtDF$Type[x]
        if (length(w)) {
          if (typ == 1) {
            w <- w[which(PG$"Potential contaminant"[w] != "+")]
          }
          if (typ == 2) {
            w <- w[which(allSAINTs$"Potential contaminant"[w] != "+")]
          }
        }
        return(w)
      })
      filtDF$Reg <- lapply(1:nr, function(x) { #x <- 1
        w <- filtDF$W[[x]]
        if (length(w)) {
          typ <- filtDF$Type[x]
          lfc <- filtDF$Filter[[x]]$Ratios[w]
          if (typ == 1) {
            reg <- data.frame("Leading protein IDs" = PG$"Leading protein IDs"[w],
                              "PG id" = PG$id[w],
                              "logFC" = lfc, # Here it would be neat to also have expression in the sample,
                              # but for this I need to rewrite filters!!!!!
                              check.names = FALSE)
            
            reg2 <- set_colnames(listMelt(strsplit(reg$"Leading protein IDs", ";"), 1:nrow(reg)), c("ID", "Row"))
            reg2$logFC <- reg$logFC[reg2$Row]
            reg2$"PG id" <- reg$"PG id"[reg2$Row]
            reg3 <- set_colnames(aggregate(reg2$logFC, list(reg2$ID), mean, na.rm = TRUE), c("ID", "logFC"))
            reg2 <- set_colnames(aggregate(reg2$"PG id", list(reg2$ID), unique), c("ID", "PG id")) # Should be only one value since we are using Leading protein IDs, which are unique
            reg2$logFC <- reg3$logFC[match(reg2$ID, reg3$ID)]
            reg2$Gene <- db$Gene[match(reg2$ID, db$`Protein ID`)]
          }
          if (typ == 2) {
            reg2 <- data.frame("ID" = allSAINTs$Protein[w],
                               "logFC" = lfc, # Here it would be neat to also have expression in the sample,
                               # but for this I need to rewrite filters!!!!!
                               "Gene" = allSAINTs$Gene[w],
                               check.names = FALSE)
            reg2$"PG id" <- tmpMap$L1[match(reg2$ID, tmpMap$value)]
          }
          reg2$TaxID <- db$TaxID[match(reg2$ID, db$"Protein ID")]
          reg2$ID <- gsub("^cRAP[0-9]{3}", "", gsub("^CON__", "", reg2$ID))
          reg2 <- aggregate(1:nrow(reg2), list(reg2$TaxID), function(x) {
            list(reg2[x,])
          })
          reg2 <- setNames(reg2$x, paste0("TaxID_", reg2$Group.1))
        } else {
          reg2 <- NA
        }
        return(reg2)
      })
      filtDF$W <- NULL
      filtDF$Filter <- NULL
      uTx <- unique(unlist(lapply(1:nr, function(x) { names(filtDF$Reg[[x]]) })))
      filtDF <- lapply(uTx, function(tx) { #tx <- uTx[1] #tx <- uTx[2]
        tmp <- filtDF
        nms <- lapply(tmp$Reg, names)
        tmp$Reg <- lapply(1:length(tmp$Reg), function(y) {
          if (tx %in% names(tmp$Reg[[y]])) {
            return(tmp$Reg[[y]][[tx]])
          } else { return() }
        })
        tmp$TaxID <- gsub("^TaxID_", "", tx)
        return(tmp)
      })
      filtDF <- plyr::rbind.fill(filtDF)
      if (!is.null(filtDF$Reg)) {
        filtDF <- filtDF[which(sapply(filtDF$Reg, function(x) { ("data.frame" %in% class(x))&&(nrow(x) > 0) })),]
      }
    }
    return(filtDF)
  })
  filtersDF <- plyr::rbind.fill(filtersDF)
  filtersDF <- filtersDF[which(paste0("TaxID_", filtersDF$TaxID) %in% names(allProteins_mapped)),]
  nr <- nrow(filtersDF)
  filtersDF <- rbind(filtersDF, filtersDF)
  filtersDF$GraphType <- GraphTypes[2]
  filtersDF$GraphType[1:nr] <- GraphTypes[1]
  txidsTst <- (length(unique(filtersDF$TaxID)) > 1)+1
  source(parSrc, local = FALSE)
  clusterExport(parClust, c("wd", "filtersDF", "allProteins_mapped", "tmpPG", "GraphTypes", "Exp", "txidsTst"), envir = environment())
  tstSTRINGs <- parLapply(parClust, 1:nrow(filtersDF), function(i) { #i <- 1 #i <- 4 #i <- 10 #i <- 11
    fltNm <- proteoCraft::cleanNms(filtersDF$Name[i], rep = "_")
    regTbl <- filtersDF$Reg[[i]]
    grphType <- filtersDF$GraphType[i]
    tstbee <- filtersDF$Group[i]
    txid <- filtersDF$TaxID[i]
    intNet_I <- NA
    STRINGplot_I <- NA
    img_I <- NA
    grphNm <- paste0(c("", paste0("taxID_", txid, "_"))[txidsTst],
                     gsub("[:\\*\\?<>\\|/| ,]+", "_", tstbee), "_", gsub("[:\\*\\?<>\\|/| ,]+", "_",
                                                                         gsub("\\(|\\)", "",
                                                                              fltNm)))
    imgpath <- paste0(wd, "/STRINGdb/", grphType, "/", grphNm, ".png")
    #cat(grphNm, "\n")
    #typ <- filtersDF$Type[i]
    proteins_mapped <- allProteins_mapped[[paste0("TaxID_", txid)]]
    proteins_mapped <- proteins_mapped[which(proteins_mapped$queryItem %in% regTbl$ID),]
    if (nrow(proteins_mapped) > 1) {
      m <- match(proteins_mapped$queryItem, regTbl$ID)
      proteins_mapped[, c("logFC", "PG id")] <- regTbl[m, c("logFC", "PG id")]
      #if (typ == 2) {
      #  proteins_mapped$PG <- tmpPG$L1[match(proteins_mapped$queryItem, tmpPG$value)]
      #}
      w <- which(nchar(proteins_mapped$Name) > 25)
      proteins_mapped$Name[w] <- paste0(substr(proteins_mapped$Name[w], 1, 22), "...")
      proteins_mapped$Name <- do.call(paste, c(proteins_mapped[, c("queryItem", "Name")],
                                               sep = "\n"))
      #if (typ == 1) {
      KOL <- "queryItem"
      #}
      #if (typ == 2) { KOL <- "PG" }
      KOL <- c(KOL, "Name")
      if ("TAIR" %in% colnames(proteins_mapped)) { # (For our dear plant-people)
        KOL <- c(KOL, "TAIR")
      }
      if (exists("intNet")) { rm(intNet) }
      # intNet <- try(rbioapi::rba_string_interactions_network(proteins_mapped$stringId,
      #                                                        txid,
      #                                                        500, network_type = tolower(grphType)),
      #               silent = TRUE)
      rqst <- paste0("https://string-db.org/api/tsv/network?identifiers=",
                     paste(proteins_mapped$stringId, collapse = "%0d"),
                     "&species=", txid,
                     "&network_type=", tolower(grphType),
                     "&required_score=500&network_flavor=confidence&show_query_node_labels=1&hide_disconnected_nodes=1")
      try({
        response <- httr::POST(rqst, encode = "json")
        if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
        intNet <- httr::content(response, "parsed", show_col_types = FALSE)
      })
      if (!"try-error" %in% class(intNet)) { # This may fail if we have too many nodes! The current limit is 2000.
        # They suggest to use their "Cytoscape stringApp", this may be worth looking into.
        if ((!is.null(intNet))&&("data.frame" %in% class(intNet))&&(nrow(intNet))) {
          for (ab in c("A", "B")) {
            m <- match(intNet[[paste0("stringId_", ab)]], proteins_mapped$stringId)
            intNet[, paste0(c("Original_", paste0(KOL, "_")), ab)] <- proteins_mapped[m, c("queryItem", KOL)]
          }
          stopifnot("character" %in% class(intNet$Original_A), "character" %in% class(intNet$Original_B))
          w <- which(proteins_mapped$stringId %in% c(intNet$stringId_A, intNet$stringId_A))
          if (length(w)) {
            intNet_I <- list("Network" = intNet,
                             "Mapping" = proteins_mapped)
            rqst <- paste0("https://string-db.org/api/highres_image/network?identifiers=",
                           paste(proteins_mapped$queryItem[w], collapse = "%0d"),
                           "&species=", txid,
                           "&network_type=", tolower(grphType),
                           "&required_score=500&network_flavor=confidence&show_query_node_labels=1&hide_disconnected_nodes=1")
            try({
              # graph <- rbioapi::rba_string_network_image(proteins_mapped$queryItem[w],
              #                                            "highres_image",
              #                                            imgpath,
              #                                            txid,
              #                                            required_score = 500,
              #                                            network_flavor = "confidence",
              #                                            use_query_labels = TRUE,
              #                                            flat_nodes = TRUE,
              #                                            hide_disconnected_nodes = TRUE # Not good enough on its own, 
              # )
              response <- httr::POST(rqst, encode = "json")
              if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
              img_I <- httr::content(response, "raw")
              writeBin(img_I, imgpath)
            }, silent = TRUE)
          }
        }
      }
    }
    return(list(Name = grphNm,
                Data = intNet_I,
                Bin = img_I,
                STRINGplot = imgpath))
  })
  STRINGplots %<o% setNames(lapply(GraphTypes, function(grphType) {
    w <- which(filtersDF$GraphType == grphType)
    sapply(w, function(x) { tstSTRINGs[[x]]$STRINGplot })
  }), GraphTypes)
  intNets %<o% setNames(lapply(GraphTypes, function(grphType) {
    w <- which(filtersDF$GraphType == grphType)
    w <- w[which(sapply(w, function(x) { ("list" %in% class(tstSTRINGs[[x]]$Data)) }))]
    setNames(lapply(w, function(x) { tstSTRINGs[[x]]$Data }),
             sapply(w, function(x) { tstSTRINGs[[x]]$Name }))
  }), GraphTypes)
  for (grphType in GraphTypes) { #grphType <- GraphTypes[1]
    if (length(STRINGplots[[grphType]])) {
      fls <- STRINGplots[[grphType]]
      fls <- fls[which(file.exists(fls))]
      if (length(fls) > 1) {
        PNGs <- parLapply(parClust, fls, readPNG)
        g <- parLapply(parClust, PNGs, function(png) { grid::rasterGrob(png, interpolate = TRUE) })
        dim <- ceiling(sqrt(length(g)))
        l <- length(gsub(".*/|\\.png$", "", fls))
        tst <- data.frame(tile = gsub(".*/|\\.png$", "", fls))
        tst$x <- rep(1:dim, dim)[1:l]
        tst$y <- as.numeric(sapply(dim:1, function(x) { rep(x, dim) }))[1:l]
        tst$Grob <- g
        tst$Lab <- gsub("\\)$", "", gsub("\\)_-_\\(", "\nvs\n", gsub("_by_condition_(\\()?", "\n", tst$tile)))
        ttl <- paste0(tolower(grphType), " STRING plots")
        plot <- ggplot(tst) +
          geom_tile(aes(x = x, y = y), width = 0.96, height = 0.96, fill = "white") + coord_equal() + theme_bw() +
          xlim(0.5, dim+0.5) + ylim(0.5, dim+1) +
          theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), plot.background = element_blank())
        for (i in 1:nrow(tst)) {
          plot <- plot + annotation_custom(tst$Grob[[i]], xmin = tst$x[i]-0.48, xmax = tst$x[i]+0.48, ymin = tst$y[i]-0.48, ymax = tst$y[i]+0.48)
        }
        plot <- plot +
          geom_text(aes(label = Lab, x = x, y = y + 0.48), cex = 2.5) +
          geom_text(label = paste0(grphType, " interactions"), x = dim/2+0.5, y = dim+0.7, cex = 5)
        poplot(plot)
        ggsave(paste0(wd, "/STRINGdb/", grphType, "/All ", ttl, ".jpeg"), plot,
               dpi = 120, width = 100, height = 100, units = "in", limitsize = FALSE)
        ggsave(paste0(wd, "/STRINGdb/", grphType, "/All ", ttl, ".pdf"), plot,
               dpi = 120, width = 100, height = 100, units = "in", limitsize = FALSE)
        ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg  = FALSE)
      }
    }
  }
  ReportCalls <- AddSpace2Report()
  #
  # Cytoscape
  if (CytoScape) {
    #if (clueGOahead) {
    #  rqst <- paste0(clueGO_URL, "/remove-all-cluego-analysis-results")
    #  try({
    #    response <- httr::DELETE(rqst, httr::timeout(10), encode = "json")
    #    if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
    #  }, silent = TRUE)
    #  clueGOahead <- FALSE
    #}
    ### Check that Cytoscape is installed and can run, then launch it.
    Src <- paste0(libPath, "/extdata/R scripts/Sources/Cytoscape_init.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    # Create directory for Cytoscape networks
    dirs <- paste0(wd, "/Cytoscape/", GraphTypes)
    for (dir in dirs) { if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) } }
    dirlist <- unique(c(dirlist, dir))
    #
    # Use Cytoscape to represent networks of gated proteins with overlayed logFC
    msg <- "Creating Cytoscape network .cx files..."
    ReportCalls <- AddMsg2Report()
    cat("Creating Cytoscape networks...\n")
    for (grphType in GraphTypes) { #grphType <- GraphTypes[1]
      cat(" ->", grphType, "\n")
      for (Nm in names(intNets[[grphType]])) { #Nm <- names(intNets[[grphType]])[1]
        cat("   -", Nm, "\n")
        intNet <- intNets[[grphType]][[Nm]]$Network
        mappings <- intNets[[grphType]][[Nm]]$Mapping
        # Legacy code for if doing by protein, not protein groups
        #mappings$name <- sapply(match(mappings$queryItem, db$"Protein ID"), function(m) {
        #  paste(c(paste0("Nm: ", db$"Common Name"[m]),
        #          paste0("Pr: ", db$"Protein ID"[m]),
        #          paste0("Gn: ", db$Gene[m])), collapse = "\n")
        #})
        #if (!grepl("^SAINTexpress_", Nm)) {
        mappings$name <- sapply(match(mappings$queryItem, db$`Protein ID`), function(m) {
          x <- c(paste0("Nm: ", db$`Common Name`[m]),
                 paste0("Pr: ", db$`Protein ID`[m]),
                 paste0("Gn: ", db$Gene[m]))
          if ("TAIR" %in% colnames(mappings)) { x <- c(paste0("TAIR: ", db$TAIR[m], " "), x) }
          return(paste(x, collapse = "\n"))
        })
        #  mappings$PG <- allSAINTs$PG_id[match(mappings$queryItem, allSAINTs$Protein)]
        #} else {
        #   mappings$name <- sapply(match(mappings$PG, PG$id), function(m) {
        #     pr <- unlist(strsplit(PG$"Leading protein IDs"[m], ";"))
        #     if (length(pr) > 1) { pr <- c(pr[1], "...") }
        #     gn <- unlist(strsplit(PG$Genes[m], ";"))
        #     if (length(gn) > 1) { gn <- c(gn[1], "...") }
        #     x <- c(paste0("Nm: ", PG$"Common Name (short)"[m]),
        #            paste0("Pr: ", paste(pr, collapse = ";")),
        #            paste0("Gn: ", paste(gn, collapse = ";")))
        #     if ("TAIR" %in% colnames(mappings)) { x <- c(paste0("TAIR: ", PG$TAIR[m], " "), x) }
        #     return(paste(x, collapse = "\n"))
        #   }) 
        #   mappings$PG <- allSAINTs$PG_id[match(mappings$queryItem, allSAINTs$Protein)]
        # }
        mappings$"Av. log10 expression" <- PG$"Av. log10 abundance"[match(mappings$PG, PG$id)]
        #if (!grepl("^SAINTexpress_", Nm)) {
        # intNet$Linkage <- do.call(paste, c(intNet[, paste0("PG_", c("A", "B"))], sep = "_"))
        # intNet <- Isapply(strsplit(unique(intNet$Linkage), "_"), unlist)
        # colnames(intNet) <- paste0("PG_", c("A", "B"))
        # for (ab in c("A", "B")) { #ab <- "A"
        #   intNet[[paste0("Name_", ab)]] <- mappings$name[match(intNet[[paste0("PG_", ab)]], mappings$PG)]
        # }
        #} else {
        intNet$Name_A <- mappings$name[match(intNet$Name_A, mappings$Name)]
        intNet$Name_B <- mappings$name[match(intNet$Name_B, mappings$Name)]
        #}
        kol <- c("Name_A", "Name_B")
        gD <- igraph::simplify(igraph::graph_from_data_frame(intNet[, kol], directed = FALSE))
        seq <- igraph::V(gD)
        logFC <- mappings$logFC[match(names(seq), mappings$name)]
        w <- which(is.na(logFC))
        logFC[w] <- 0
        gD <- igraph::set_vertex_attr(gD, "logFC", igraph::V(gD), logFC)
        Xprss <- mappings$"Av. log10 expression"[match(names(seq), mappings$name)]
        gD <- igraph::set_vertex_attr(gD, "Avg_expression", igraph::V(gD), Xprss)
        #igraph::vcount(gD)
        #igraph::ecount(gD)
        degAll <- igraph::degree(gD, v = igraph::V(gD), mode = "all")
        betAll <- igraph::betweenness(gD, v = igraph::V(gD), directed = FALSE) / (((igraph::vcount(gD) - 1) * (igraph::vcount(gD)-2)) / 2)
        betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
        rm(betAll)
        dsAll <- igraph::similarity(gD, vids = igraph::V(gD), mode = "all", method = "dice")
        gD <- igraph::set_vertex_attr(gD, "degree", index = igraph::V(gD), value = degAll)
        gD <- igraph::set_vertex_attr(gD, "betweenness", index = igraph::V(gD), value = betAll.norm)
        #summary(gD)
        F1 <- function(x) {
          data.frame(V4 = dsAll[which(igraph::V(gD)$name == as.character(x$V1)),
                                which(igraph::V(gD)$name == as.character(x$V2))])
        }
        dataSet.ext <- plyr::ddply(intNet[, kol], .variables = kol, function(x) data.frame(F1(x)))
        gD <- igraph::set_edge_attr(gD, "weight", index = igraph::E(gD), value = 0)
        gD <- igraph::set_edge_attr(gD, "similarity", index = igraph::E(gD), value = 0)
        for (i in 1:nrow(dataSet.ext)) {
          igraph::E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
          igraph::E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)
        }
        rm(dsAll, i, F1)
        #summary(gD)
        tst <- RCy3::createNetworkFromIgraph(gD, new.title = Nm)
        RCy3::layoutNetwork("force-directed defaultSpringLength=70 defaultSpringCoefficient=0.000003")
        rg <- max(abs(logFC))
        try(RCy3::setNodeColorMapping("logFC", c(-rg, 0, rg),
                                      c("#FF0000", "#999999", "#00FF00"), style.name = "default"), silent = TRUE)
        RCy3::setNodeShapeDefault("ellipse", style.name = "default")
        RCy3::setNodeFontSizeDefault(12, style.name = "default")
        RCy3::setNodeSizeMapping("Avg_expression", style.name = "default")
        RCy3::exportNetwork(paste0(wd, "/Cytoscape/", grphType, "/", Nm, "_", tolower(grphType), " network"), "CX")
        RCy3::deleteAllNetworks()
      }
    }
  }
  # Then close Cytoscape:
  RCy3::closeSession(save.before.closing = FALSE)
  cmd <- paste0("taskkill/im \"", gsub(".+/", "", CytoScVrs), "\" /f")
  #cat(cmd)
  shell(cmd)
}

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
#         sapply(unique(ev$"Raw file"), function(y) { length(which((temp$Type == x)&(temp$"Raw file" == y))) })
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
for (i in 1:length(MatMetCalls$Texts$DatAnalysis)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls, paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$DatAnalysis[", i,"], prop = WrdFrmt$",
                                                        c("Body", "Template_text")[(MatMetCalls$Texts$DatAnalysis[i] == "TEMPLATE")+1], "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#

# Finalize analysis
Src <- paste0(libPath, "/extdata/R scripts/Sources/Finalize_analysis.R")
#rstudioapi::documentOpen(Src)
#loadFun(BckUpFl)
source(Src, local = FALSE)

# End logging:
sink(NULL, type = "message")
#close(logcon)
rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
parLapply(parClust, 1:N.clust, function(x) {
  rm(list = ls())
  gc()
})
rm(ReportCalls) # Temporary fix until I figure out how to fix the grphtype bug - I thought I had
setwd(wd); saveImgFun(BckUpFl) # Leave an ultimate backup in the temporary folder
#loadFun(BckUpFl)

# Save final state of the environment
# This is done within the destination folder (outdir) because it will restart the session so has to be done last
# (this will interrupt the script flow so all commands queued after that are gone)
setwd(procdir)
pkgs <- gtools::loadedPackages()
dscrptFl <- paste0(procdir, "/DESCRIPTION")
tmp <- paste0("", do.call(paste, c(pkgs[, c("Name", "Version")], sep = " (")), ")")
tmp <- paste0("Depends: ", paste(tmp, collapse = ", "))
write(tmp, dscrptFl)
renv::snapshot(force = TRUE, prompt = FALSE, type = "explicit")
if ((exists("renv"))&&(renv)) { try(renv::deactivate(), silent = TRUE) }

### That's it, done!
#openwd(outdir)
#rm(list = ls())
