#### Code chunk - Initialization
if (!interactive()) { stop("This script should only be run within an interactive R session!") }
options(stringsAsFactors = FALSE)
options(install.packages.compile.from.source = "never")
#rm(list = ls()[which(!ls() %in% c("dtstNm", "wd", "indir", "outdir"))])

## The proteoCraft package can be re-installed at any time in the workflow (there is a specific script for this in the package's library folder),
## or just load it here:
if (exists(".obj")) { rm(".obj") }
require(proteoCraft)
dirlist %<o% c() # This should go!!!
ReUseAnsw %<o% FALSE
ReLoadPSMsBckp %<o% FALSE

RPath %<o% as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath %<o% paste0(RPath, "/proteoCraft")
homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
fls <- paste0(homePath, "/", c("Regulation analysis - master script.R",
                               "Regulation analysis - detailed script.R",
                               "Regulation analysis - detailed script_pepOnly.R",
                               "No replicates analysis - detailed script.R",
                               "Reload_renv_from_lock_file.R",
                               "Default_locations.xlsx",
                               "LC_columns.xlsx"))
tst <- sum(!file.exists(fls))
if (tst) { proteoCraft::Configure() }
scrptType %<o% "withReps"

# Parameters used by the master script:
###-|-### Workflows: setNames(c("Differential Protein Expression analysis", "Pull-Down (e.g. co-IP)", "Biotin-based Pull-Down (BioID, TurboID, APEX...)", "Time Course","SubCellular Localisation analysis"), c("REGULATION", "PULLDOWN", "BIOID", "TIMECOURSE", "LOCALISATION"))
###-|-### Replicates? TRUE
###-|-### External dependencies: Excel (loose); ScanHeadsman (loose); Cytoscape (loose); saintExpress (auto)

### Packages
## For convenience all (or most) of the packages used are loaded or installed here:
## CRAN packages:
if(!exists("cran_req")) { cran_req %<o% "pak" } else { cran_req %<o% cran_req }
if(!exists("bioc_req")) { bioc_req %<o% c() } else { bioc_req %<o% bioc_req }
cran_req <- unique(c(cran_req, "pak", "shiny", "renv", "R.utils", #"uchardet", # Should not be necessary anymore since Rcy3 replaced it with stringi in version 2.24.0
                     "qs2", "shinyWidgets", "DT", "shinyBS", "stringr", "gplots", "ggplot2",
                     "ggpubr", "gtools", "reshape", "reshape2", "compiler", "stats", "rgl", "ggrepel", "rstudioapi", "modeest", "minpack.lm",
                     "snow", "viridis", "pcaMethods", "impute", "imputeLCMD", "parallel", "coin", "openxlsx", "openxlsx2", "plotly", "Peptides",
                     "xml2", "pdftools", "statmod", "ggpolypath", "venn", "gridExtra", "svDialogs", "htmlwidgets", "magrittr", "tibble", "fs",
                     "officer", "hexbin", "igraph", "matlib", "umap", "plyr", "ggnewscale", "shinyjs", "shinyFiles", "TeachingDemos", "shinycssloaders",
                     "tidyr", "data.table", "ggplotify", "jpeg", "scattermore", "rpanel", "stringi", "lmtest", "ssh", "taxize", "arrow"))
bioc_req <- unique(c(bioc_req, "biomaRt", "GO.db", "UniProt.ws", "limma", "sva", "qvalue", "MSnbase", "DEP",
                     "Rgraphviz", "RCy3", "siggenes", "DEqMS", "rawrr"))
inst <- as.data.frame(installed.packages())
for (pack in cran_req) {
  if (!pack %in% inst$Package) {
    if (pack %in% c("pak", #"shiny",
                    "uchardet", "openxlsx2", "taxize")) {
      # Exceptions where for now we want a specific version to be installed,
      # or have to help the installer so it finds the right location
      if (pack == "pak") {
        install.packages("pak")
      }
      # if (pack == "shiny") { # Should be fixed now
      #   install.packages("https://cran.r-project.org/src/contrib/Archive/shiny/shiny_1.7.5.tar.gz", dependencies = TRUE)
      # }
      if (pack == "uchardet") {
        url <- "https://cran.r-project.org/src/contrib/Archive/uchardet/uchardet_1.1.1.tar.gz"
        destfile <- "uchardet_1.1.1.tar.gz"
        tst <- try(download.file(url, destfile, "curl"), silent = TRUE)
        if ("try-error" %in% class(tst)) { try(download.file(url, destfile, "wget"), silent = TRUE) }
        install.packages(destfile)
        unlink(destfile)
      }
      if (pack == "openxlsx2") {
        pak::pkg_install("JanMarvin/openxlsx2@v1.10", ask = FALSE) # ... until I can figure out what is happening...
      }
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
#devtools::install_github("cpanse/rawrr")
#rawrr::installRawFileReaderDLLs() # Deprecated
tst <- try(normalizePath(rawrr:::.rawrrAssembly(), winslash = "/"), silent = TRUE)
if (("try-error" %in% class(tst))||(!file.exists(tst))) {
  rawrr::installRawrrExe()
}

# Fast save and load functions
Src <- paste0(libPath, "/extdata/R scripts/Sources/Save_Load_fun.R")
#rstudioapi::documentOpen(Src)
source(Src)

# Set Shiny options, load functions for creating a Word report, create Excel styles
Src <- paste0(libPath, "/extdata/R scripts/Sources/ShinyOpt_Styles_and_Report.R")
#rstudioapi::documentOpen(Src)
source(Src)

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

# Update the proteoCraft package? (Usually a good idea)
msg <- "Should we update the proteoCraft package? (recommended)"
updt_proteoCraft %<o% c(TRUE, FALSE)[match(svDialogs::dlg_message(msg, "yesno")$res, c("yes", "no"))]

# Define input, output, project folder etc...
Src <- paste0(libPath, "/extdata/R scripts/Sources/Start_analysis.R")
#rstudioapi::documentOpen(Src)
source(Src)
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
source(parSrc)
setDTthreads(threads = N.clust)

# Load PSMs
Src <- paste0(libPath, "/extdata/R scripts/Sources/Load_PSMs.R")
#rstudioapi::documentOpen(Src)
source(Src)

# MS raw files map
tstFrMp <- FALSE
while(!tstFrMp) {
  Src <- paste0(libPath, "/extdata/R scripts/Sources/Fractions_Map_editor.R")
  #rstudioapi::documentOpen(Src)
  source(Src)
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
        if (is.character(vals)) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
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
        if (is.character(vals)) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
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
              if (is.character(vals)) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
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
              if (is.character(vals)) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
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
eval(parse(text = runApp))
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
  source(Src)
}
#

#### Code chunk - Load and process search database(s)
Src <- paste0(libPath, "/extdata/R scripts/Sources/Process_Fasta_DBs.R")
#rstudioapi::documentOpen(Src)
source(Src)

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
      source(parSrc)
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
    stop()
    ev <- ev[which(!is.na(ev2fr)),]
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
source(Src)
#
# Proteins of interest
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList.R")
#rstudioapi::documentOpen(Src)
source(Src)
#
# Targets
# Sometimes the user does not fill the Target factor with valid protein IDs... but this is what we would actually need.
# Here, if necessary, we will remap those to valid IDs:
Src <- paste0(libPath, "/extdata/R scripts/Sources/Targets.R")
#rstudioapi::documentOpen(Src)
source(Src)
#
# KnockOut, KnockIn of KnockDown
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
source(Src)
#
# Protein headers for shiny (update)
Src <- paste0(libPath, "/extdata/R scripts/Sources/protHeaders_for_shiny.R")
#rstudioapi::documentOpen(Src)
source(Src)

#### Code chunk - Define analysis parameters
# Create first PCA to check on sample relationships
if ((length(MQ.Exp) > 1)||(LabelType == "Isobaric")) { # Should be always TRUE
  source(parSrc)
  data <- ev
  colnames(data)[which(colnames(data) == "MQ.Exp")] <- "Parent sample"
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
  ReportCalls <- AddSpace2Report()
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              "body_add_fpar(Report, fpar(ftext(\"PSMs-level PCA plot:\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
  ReportCalls$Calls <- append(ReportCalls$Calls, list())
  dir <- paste0(wd, "/Dimensionality red. plots/PCA")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  LRepCalls <- length(ReportCalls$Calls)
  lsKl <- c("Modified sequence", Y)
  if (LabelType == "LFQ") { lsKl <- c(lsKl, X) }
  ls <- lapply(lsKl, function(kl) { data[[kl]] })
  tmp <- do.call(paste, c(data[, lsKl], sep = "---"))
  if (LabelType == "Isobaric") {
    kol2 <- gsub(topattern(ev.ref["Original"]), "", kol)
    data2 <- as.data.table(data[, kol])
    colnames(data2) <- kol
    data2$Group <- tmp
    data2 <- data2[, lapply(.SD, sum, na.rm = TRUE), keyby = Group]
    colnames(data2) <- c("Group", kol)
    data2 <- as.data.frame(data2)
    data2[, kol2] <- log10(data2[, kol])
    data2 <- data2[, which(!colnames(data2) %in% kol)]
    data2[, lsKl] <- data[match(data2$Group, tmp), lsKl]
    data2$Group <- NULL
  }
  if (LabelType == "LFQ") {
    if (X == "Parent sample") { kol2 <- get("MQ.Exp") } else { kol2 <- get(X) }
    data2 <- data.table(Intensity = data[[ev.col["Original"]]], Group = tmp)
    data2 <- data2[, list(`log10(Intensity)` = sum(Intensity, na.rm = TRUE)),
                   keyby = Group]
    data2$`log10(Intensity)` <- log10(data2$`log10(Intensity)`)
    data2 <- as.data.frame(data2)
    data2[, lsKl] <- data[match(data2$Group, tmp), lsKl]
    data2$Group <- NULL
    data <- spread(data2, X, "log10(Intensity)")
  }
  data <- data[, kol2]
  w <- which(is.na(data), arr.ind = TRUE)
  if (nrow(w)) {
    groups <- rep(1, length(MQ.Exp))
    temp <- Data_Impute2(data, groups, is.log = FALSE)
    data <- temp$Imputed_data
  }
  pcA <- prcomp(t(data[, kol2]), scale. = TRUE)
  if (length(pcA$rotation)) {
    scoresA <- as.data.frame(pcA$x)
    if ("PC2" %in% colnames(scoresA)) {
      scoresA$Sample <- rownames(scoresA)
      rownames(scoresA) <- NULL
      pvA <- round(100*(pcA$sdev)^2 / sum(pcA$sdev^2), 0)
      pvA <- pvA[which(pvA > 0)]
      pvA <- paste0("Original: ", paste(sapply(1:length(pvA), function(x) {
        paste0("PC", x, ": ", pvA[x], "%")
      }), collapse = ", "))
      scoresA$Colour <- scoresA$Label <- scoresA$Sample
      ttl <- "PCA plot - Samples (PSMs-level)"
      plot <- ggplot(scoresA) +
        geom_point(aes(x = PC1, y = PC2, colour = Colour)) +
        scale_color_viridis_d(begin = 0.25) +
        coord_fixed() + theme_bw() +
        geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
        ggtitle(ttl, subtitle = pvA) +
        geom_text_repel(aes(x = PC1, y = PC2, label = Label, colour = Colour),
                        size = 2.5, show.legend = FALSE)
      #poplot(plot)
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 20, height = 20, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 20, height = 20, units = "in")
      ReportCalls <- AddPlot2Report(Space = FALSE)
      Symb <- "circle"
      if ("PC3" %in% colnames(scoresA)) {
        plot_lyPSMsPCA <- plot_ly(scoresA, x = ~PC1, y = ~PC2, z = ~PC3,
                                  color = ~Colour, text = ~Label, type = "scatter3d", mode = "markers",
                                  symbol = I(Symb))
        #plot_lyPSMsPCA <- add_trace(plot_lyPSMsPCA, scoresA, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "text", showlegend = FALSE)
      } else {
        plot_lyPSMsPCA <- plot_ly(scoresA, x = ~PC1, y = ~PC2, color = ~Colour, text = ~Label, type = "scatter", mode = "markers",
                                  symbol = I(Symb))
        #plot_lyPSMsPCA <- add_trace(plot_lyPSMsPCA, scoresA, x = ~PC1, y = ~PC2, type = "scatter", mode = "text", showlegend = FALSE)
      }
      plot_lyPSMsPCA %<o% layout(plot_lyPSMsPCA, title = ttl)
      saveWidget(plot_lyPSMsPCA, paste0(dir, "/", ttl, ".html"), selfcontained = TRUE)
      #system(paste0("open \"", dir, "/", ttl, ".html"))
    } else {
      msg <- "Not enough valid data to draw a PSM-level PCA plot!"
      ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
    }
  }
  ReportCalls$Calls[[LRepCalls]] <- append(ReportCalls$Calls[[LRepCalls]],
                                           "body_add_par(Report, \"\", style = \"Normal\")")
} else {
  stop("Uh, I think you have the wrong analysis pipeline here...\nwhere are my sample groups and replicates!?!?!")
}
#
# Defaults
RefRat_Mode %<o% "2" # Values: RefRat_Mode = "2" or "1" # For now not a user-modifiable parameter, however this may change
WelchRoot %<o% "Welch's t-test -log10(Pvalue) - "
modRoot %<o% "Moderated t-test -log10(Pvalue) - "
permRoot %<o% "Permutations t-test -log10(Pvalue) - "
samRoot %<o% "SAM -log10(Pvalue) - "
pvalue.col %<o% c(WelchRoot, modRoot, permRoot, samRoot)
names(pvalue.col) <- sapply(pvalue.col, function(x) { unlist(strsplit(x, "\\.|\\'|\\ "))[1] })
ParamFls <- c(paste0(wd, "/Parameters.csv"),
              paste0(libPath, "/extData/Parameters_template.csv"))
ParamFl <- ParamFls[1]
if (!file.exists(ParamFl)) { ParamFl <- ParamFls[2] }
Param_Help <- read.csv(ParamFl, header = FALSE)
Param_Help <- Param_Help$V3
Param %<o% Param.load(ParamFl)
Param$vCPUs <- N.clust
Param$WD <- wd
if (ParamFl == ParamFls[2]) {
  Param$WD <- wd
  Param$Project <- dtstNm
  Param$Output <- SearchSoft
  Param$Fix.MQ.Isobaric.labels <- FALSE
  Param$Type <- WorkFlow
  Param$Label <- LabelType
  Param$MQ.Experiments <- paste(MQ.Exp, collapse = ";")
  Param$Search.DB <- paste(fastasTbl$Full, collapse = ";")
  Param$Search.DB.type <- paste(fastasTbl$Type, collapse = ";")
  Param$Search.DB.species <- paste(fastasTbl$Species, collapse = ";") 
  Param$Two.sided <- !(WorkFlow %in% c("PULLDOWN", "BIOID"))
  Param$Min.Pep.Size <- MinPepSz
  Param$PSMs <- paste(PSMsFl, collapse = ";")
  if (LabelType == "Isobaric") {
    Param$Label <- IsobarLab #
    Param$Label.Multiplicity <- length(get(IsobarLab)) #
    Param$Norma.Pep.Intens.IRS = "TF_TRUE"
    Param$Norma.Pep.Intens.IRS_Ref_channels = "SELECT"
    Param$Label.Purities.file <- "CHOOSEFL"
  }
  if ("Target" %in% Factors) {
    tmp <- FactorsLevels$Target
    tmp <- tmp[which((!is.na(tmp))&(tmp != "NA"))]
    Param$Prot.list <- paste(unique(c(tmp, unlist(strsplit(Param$Prot.list, ";")))), collapse = ";")
    Param$Prot.list_pep <- paste(unique(c(tmp, unlist(strsplit(Param$Prot.list_pep, ";")))), collapse = ";")
  }
  ptmDflt2 <- Param$PTM.analysis <- paste(grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE), collapse = ";")
  fTstDflt <- Param$F.test
  if (grepl("^TF_((TRUE)|(FALSE))$", fTstDflt)) { fTstDflt <- as.logical(gsub("^TF_", "", fTstDflt)) }
  if (!is.logical(fTstDflt)) { fTstDflt <- TRUE }
  Param$F.test <- fTstDflt
  Param$GO.enrichment <- TRUE
} else {
  ptmDflt2 <- ""
  if ("PTM.analysis" %in% colnames(Param)) {
    tmp <- unlist(strsplit(as.character(Param$PTM.analysis), ";"))
    tmp <- tmp[which(tmp %in% Modifs$`Full name`)]
    ptmDflt2 <- tmp
  }
}
if ((Param$Label == "LFQ")&&(isDIA)) { Param$Label <- "DIA" }
Param$WD <- wd # Super important!
ptmDflt1 <- grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE, invert = TRUE)
if (!"PTM.analysis" %in% colnames(Param)) { Param$PTM.analysis <- paste(ptmDflt2, collapse = ";") }
Mod4Quant %<o% Modifs$Mark[match(ptmDflt1, Modifs$`Full name`)]
Mod2Xclud %<o% set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                            c("Mark", "Where"))
goDflt <- suppressWarnings(as.logical(Param$GO.enrichment))
if ((!is.logical(goDflt))||(is.na(goDflt))) { goDflt <- TRUE }
fTstDflt <- suppressWarnings(as.logical(Param$F.test))
if ((!is.logical(fTstDflt))||(is.na(fTstDflt))) { fTstDflt <- TRUE }
# if ((!"Param_suppress_UI" %in% colnames(Param))||(!is.logical(Param$Param_suppress_UI))) {
#   # This is to allow bypassing UI-based parameters creation.
#   # Could be useful in cases you have created custom values for parameters which the script can handle but which,
#   # because I have made a mistake (those do happen), is getting overwritten by the UI.
#   Param$Param_suppress_UI <- FALSE
#   #
#   # Honestly, this is bloat and needs to go!
# }
if (!"Ratios.Groups_Nested" %in% colnames(Param)) { Param$Ratios.Groups_Nested <- WorkFlow != "Regulation" }
#Param$Ratios.Groups.Ref.Aggregate.Level <- "AUTOFACT"
#Param$Ratios.Groups.Ref.Aggregate.Level <- "Exp;Con;Rep"
tmp <- unlist(strsplit(Param$Ratios.Groups.Ref.Aggregate.Level, ";"))
tst <- (length(tmp)>1)||(!tmp %in% c("AUTOFACT", "MAP2FACTS"))
if (tst) {
  tst <- sum(!tmp %in% substr(Factors, 1, 3)) == 0
  if (!tst) { tmp <- "MAP2FACTS" }
}
if ((length(tmp) == 1)&&(tmp == "AUTOFACT")) {
  Param$Ratios.Groups.Ref.Aggregate.Level <- paste(substr(Factors, 1, 3), collapse = ";")
}
klustChoices %<o% c("K-means", "hierarchical")
KlustMeth %<o% 2
#
saintExprs %<o% (WorkFlow %in% c("PULLDOWN", "BIOID"))
if (Annotate) {
  allGO <- unique(unlist(strsplit(db$GO[which(!is.na(db$GO))], ";")))
  allGO2 <- paste0("GO:", gsub(".* \\[GO:|\\]$", "", allGO))
  dftlGO2 <- unique(unlist(strsplit(Param$GO.tabs, ";")))
  dftlGO2 <- dftlGO2[which(dftlGO2 %in% allGO2)]
  dftlGO <- allGO[match(dftlGO2, allGO2)]
  w <- c(which(allGO %in% dftlGO),
         which(!allGO %in% dftlGO))
  allGO <- allGO[w]; allGO2 <- allGO2[w]
}
if (!"GO.terms.for.proteins.of.interest" %in% colnames(Param)) { Param$GO.terms.for.proteins.of.interest <- FALSE }
if (!"Amica" %in% colnames(Param)) { Param$Amica <- TRUE }
if (!"Mirror.Ratios" %in% colnames(Param)) { Param$Mirror.Ratios <- FALSE }
if (!"Custom.PGs" %in% colnames(Param)) { Param$Custom.PGs <- "" }
if (!"TrueDisc_filter" %in% colnames(Param)) { Param$TrueDisc_filter <- "" }
DiscFiltModes %<o% c("Positive filter", "Negative filter", "Filter column")
DiscFiltModesHlp <- c(" - Positive filter: TRUE means values should not be set to NA for PGs matching this protein - values are set to NA for PGs without match or with match to FALSE)",
                      " - Negative filter: TRUE means values should be set to NA for PGs matching this protein - values are left unchanged for PGs without match or with match to FALSE)",
                      " - Filter column: Values are unaffected, but a new column is created marking with \"+\" proteins found in the provided filter.")
if (!"TrueDisc_filter_mode" %in% colnames(Param)) { Param$TrueDisc_filter_mode <- DiscFiltModes[1] }
if (!"CRAPome_file" %in% colnames(Param)) { Param$CRAPome_file <- "" }
shpDflt <- c("none", "vsn", "loess")
names(shpDflt) <- gsub("^none$", "FALSE", shpDflt)
shpDflt <- shpDflt[match(as.character(Param$Norma.Pep.Intens.Shape), names(shpDflt))]
if (is.na(shpDflt)) { shpDflt <- "none" }
pr <- c("Norma.Ev.Intens", "Norma.Pep.Intens", "Adv.Norma.Pep.Intens", "Norma.Pep.Intens.IRS", "Norma.Prot.Ratio", "Adv.Norma.Prot.Intens")
for (p in pr) { if ((!p %in% colnames(Param))||(!is.logical(Param[[p]]))||(is.na(Param[[p]]))) { Param[[p]] <- TRUE } }
p <- "Adv.Norma.Ev.Intens"
if ((!p %in% colnames(Param))||(!is.logical(Param[[p]]))||(is.na(Param[[p]]))) { Param[[p]] <- (length(unique(FracMap$Fraction)) > 1)|(length(unique(FracMap$`PTM-enriched`)) > 1) } 
if (!toupper(as.character(Param$Norma.Pep.Intens.Shape)) %in% c("FALSE", "VSN", "LOESS")) { Norma.Pep.Intens.Shape <- FALSE }
pr <- c("Norma.Pep.Ratio", "Adv.Norma.Pep.Ratio", "Norma.Prot.Ratio.to.Biot")
for (p in pr) { if ((!is.logical(Param[[p]]))||(is.na(Param[[p]]))) { Param[[p]] <- FALSE } }
QuantMethods %<o% setNames(c("Prot.Quant", "Prot.Quant + weights", "Prot.Quant.Unique", "Prot.Quant.Unique + weights",
                             "Prot.Quant2 + weights", "Prot.Quant2", "IQ_MaxLFQ", "Top3", "Top1"),
                           c(paste0("Profile_avg.", c("", ", weights = -log10(PEP)/CV", c(", unique peptides in priority", ", weights = -log10(PEP)/CV, unique peptides in priority"))),
                             paste0("Profile_avg.v2", c(", weights = -log10(PEP)/CV", "")), "MaxLFQ (iq)", "Top3", "Top1"))
QMdef <- "Prot.Quant.Unique"
if (("QuantMeth" %in% colnames(Param))&&(Param$QuantMeth %in% QuantMethods)) { QMdef <- Param$QuantMeth } else {
  Param$QuantMeth <- QMdef
}
if (!QMdef %in% QuantMethods[1:6]) { QMdef <- "Prot.Quant.Unique" }
QMdefnm <- names(QuantMethods)[match(QMdef, QuantMethods)]
dfltP4Q <- "Razor"
if (("Prot.Quant.Use" %in% colnames(Param))&&(!gsub(" |_|-|\\.", "", toupper(Param$Prot.Quant.Use)) %in% c("UNIQUE", "RAZOR", "ALL"))) {
  dfltP4Q <- Param$Prot.Quant.Use
} else {
  Param$Prot.Quant.Use <- dfltP4Q
}
if (("Update_Prot_matches" %in% colnames(Param))&&(is.logical(Param$Update_Prot_matches))&&(!is.na(Param$Update_Prot_matches))) {
  Update_Prot_matches %<o% Param$Update_Prot_matches
} else {
  Update_Prot_matches %<o% TRUE
  Param$Update_Prot_matches <- Update_Prot_matches
}
if (("Reuse_Prot_matches" %in% colnames(Param))&&(is.logical(Param$Reuse_Prot_matches))&&(!is.na(Param$Reuse_Prot_matches))) {
  Reuse_Prot_matches %<o% Param$Reuse_Prot_matches
} else {
  Reuse_Prot_matches %<o% ("evmatch.RData" %in% list.files(wd))
  Param$Reuse_Prot_matches <- Reuse_Prot_matches
}
if (("Pep.Impute" %in% colnames(Param))&&(is.logical(Param$Pep.Impute))&&(!is.na(Param$Pep.Impute))) {
  Impute %<o% as.logical(Param$Pep.Impute)
} else {
  Impute %<o% FALSE
  Param$Pep.Impute <- Impute
}
PepFoundInAtLeast %<o% 1
if ("PepFoundInAtLeast" %in% colnames(Param)) {
  PepFoundInAtLeast %<o% suppressWarnings(as.integer(Param$PepFoundInAtLeast))
  if ((is.na(PepFoundInAtLeast))||(PepFoundInAtLeast < 1)||(PepFoundInAtLeast > nrow(Exp.map))) {
    warning("Invalid \"PepFoundInAtLeast\" parameter, defaulting to 1")
    PepFoundInAtLeast <- 1
  }
}
Param$PepFoundInAtLeast <- PepFoundInAtLeast
mxRp <- max(Exp.map$Replicate)
mxN <- max(c(2, mxRp-1))
PepFoundInAtLeastGrp %<o% mxN
if ("PepFoundInAtLeastGrp" %in% colnames(Param)) {
  PepFoundInAtLeastGrp %<o% suppressWarnings(as.integer(Param$PepFoundInAtLeastGrp))
  if ((is.na(PepFoundInAtLeastGrp))||(PepFoundInAtLeastGrp < 1)||(PepFoundInAtLeastGrp > mxRp)) {
    warning(paste0("Invalid \"PepFoundInAtLeastGrp\" parameter, defaulting to ", mxN))
    PepFoundInAtLeastGrp <- mxN
  }
}
Param$PepFoundInAtLeastGrp <- PepFoundInAtLeastGrp
Param$PepFoundInAtLeastGrp <- PepFoundInAtLeastGrp
CytoScExe %<o% c()
tmp <- grep("cytoscape", list.dirs("C:/PROGRA~1", recursive = FALSE), value = TRUE, ignore.case = TRUE)
CytoScape %<o% (length(tmp) > 0)
if ("Cytoscape" %in% colnames(Param)) { CytoScape <- Param$Cytoscape }
if (length(tmp)) {
  CytoScExe <- sapply(tmp, function(x) { grep("Cytoscape\\.exe$", list.files(x, recursive = TRUE), value = TRUE) })
  CytoScExe <- paste0(tmp, "/", CytoScExe)
  if (length(CytoScExe) > 1) {
    tst <- sapply(CytoScExe, function(x) { file.info(x)$mtime })
    CytoScExe <- CytoScExe[order(tst, decreasing = TRUE)]
  }
} else {
  msg <- "Could not locate Cytoscape executable!"
  ReportCalls <- AddMsg2Report()
  CytoScape <- FALSE
}
if ("CytoScapePath" %in% colnames(Param)) {
  tmp <- normalizePath(Param$CytoScapePath, winslash = "/")
  if ((length(CytoScExe))&&(tmp %in% CytoScExe)) { CytoScExe <- tmp }
} else { Param$CytoScapePath <- CytoScExe[1] }
SpeciesTst %<o% "Unspecified"
if ("Taxonomy" %in% colnames(db)) {
  SpeciesTst <- unique(db$Taxonomy[which(gsub(" *(\\(|\\[).*", "", db[[dbOrgKol]]) == mainOrg)])
  SpeciesTst <- SpeciesTst[which(as.character(SpeciesTst) != "NA")][1]
}
KingdomTst %<o% aggregate(db$Kingdom, list(db$Kingdom), length)
KingdomTst <- KingdomTst[order(KingdomTst$x, decreasing = TRUE),]
KingdomTst <- KingdomTst$Group.1[1]
isEukaLike %<o% (KingdomTst %in% c("Eukaryota", "Archaea"))
if (("ProtRul" %in% colnames(Param))&&(is.logical(Param$ProtRul))&&(!is.na(Param$ProtRul))) {
  protrul %<o% Param$ProtRul
} else {
  protrul %<o% (!WorkFlow %in% c("PULLDOWN", "BIOID"))
  # Archaea and Eukaryotes have introns and histones, Bacteria do not
  protrul <- c(protrul, FALSE)[(!isEukaLike)+1]
  Param$ProtRul <- protrul
}
if (("ProtRulNuclL" %in% colnames(Param))&&(!is.na(as.integer(Param$ProtRulNuclL)))) {
  ProtRulNuclL <- as.integer(Param$ProtRulNuclL)
} else {
  ProtRulNuclL <- 196
  Param$ProtRulNuclL <- ProtRulNuclL
}
threshMsg %<o% paste0("% of ", c("control-to-control", "intra-sample group")[match(RefRat_Mode, c("1", "2"))], " ratios")
threshOpt %<o% c(threshMsg, "Absolute log2 FC threshold")
threshDflt <- threshOpt[1]
if (!"Ratios.Thresholds" %in% colnames(Param)) {
  Param$Ratios.Thresholds <- "% of intra-sample group ratios"
}
if (!Param$Ratios.Thresholds %in% threshOpt) {
  Param$Ratios.Thresholds <- "% of intra-sample group ratios"
}
if (Param$Ratios.Thresholds ==  "Absolute threshold") {
  Param$Ratios.Thresholds <- "Absolute log2 FC threshold"
}
threshDflt <- Param$Ratios.Thresholds
Mitch <- match(threshDflt, threshOpt)
if ("Ratios.Contamination.Rates" %in% colnames(Param)) {
  KontRt <- Param$Ratios.Contamination.Rates
  if ((!is.numeric(KontRt))||(KontRt < 0)) { KontRt <- c(0.05, 1)[Mitch] }
} else { KontRt <- c(0.05, 1)[Mitch] }
KontRt <- KontRt*c(100, 1)[Mitch]
KontGrps <- c("Ratio groups", "Experiments", "Whole dataset")
if ("Ratios.Contaminant.Groups" %in% colnames(Param)) {
  KontGrp <- Param$Ratios.Contaminant.Groups
} else { KontGrp <- KontGrps[1] }
if (!KontGrp %in% KontGrps) { KontGrp <- KontGrps[1] }
if (!exists("minInt")) { minInt <- 100 }
if ("Min.Intensity" %in% colnames(Param)) {
  tmp <- suppressWarnings(as.numeric(Param$Min.Intensity))
  if ((is.numeric(tmp))&&(is.finite(tmp))&&(tmp >= 0)) { minInt <- tmp }
}
minInt %<o% minInt
if ("BH.FDR.values" %in% colnames(Param)) {
  BH.FDR <- sort(as.numeric(unlist(unique(strsplit(as.character(Param$BH.FDR.values), ";")))))
} else { BH.FDR <- c(0.1, 0.2, 0.3) }
BH.FDR %<o% BH.FDR
# For PTMs normalisation to the parent PG, when no quantitation is available for the parent,
# how are we replacing NAs? 
NAsReplMeth %<o% 2
if ("PTM.analysis_NAsReplaceMethod" %in% colnames(Param)) {
  NAsReplMeth %<o% as.integer(Param$PTM.analysis_NAsReplaceMethod)
  if (!NAsReplMeth %in% 1:2) { NAsReplMeth <- 2 }
}
NAsReplMethods <- c("Impute", # Currently not recommended, until I can work out a much less random Imputation method
                    "Median" # Recommended method
                    )
#
TwoSidedDeflt <- c("Both directions", "Up-only")[(WorkFlow %in% c("PULLDOWN", "BIOID"))+1]
if (("Two.sided" %in% colnames(Param))&&(!is.na(Param$Two.sided))&&(is.logical(Param$Two.sided))) {
  TwoSidedDeflt <- Param$Two.sided
}
Mirror.Ratios %<o% FALSE
if ("Mirror.Ratios" %in% colnames(Param)) { Mirror.Ratios <- Param$Mirror.Ratios <- as.logical(Param$Mirror.Ratios) }
Param$Mirror.Ratios <- Mirror.Ratios
TwoSidedDeflt <- c(c("Up-only", "Down-only")[Mirror.Ratios+1], "Both directions")[TwoSidedDeflt+1]
#
mnFct <- c("Experiment", "Replicate")
if (WorkFlow == "TIMECOURSE") { mnFct <- c(mnFct, "Time.point") }
l <- length(mnFct)
coreNms <- setNames(c("Ratios.Groups.Ref.Aggregate.Level", "Norm.Groups", "Ratios.Groups", "Batch.correction", "GO.enrichment.Ref.Aggr"),
                    c(paste0("Individual Sample___Combination of Factors required to distinguish individual samples./nMust include ",
                             paste(mnFct[1:(l-1)], collapse = ", "), " and ", mnFct[l], "."),
                      "Normalisation groups___Factor(s) defining groups of samples to normalize to each-other. Note that which, if any, normalisations apply will be defined further down",
                      "Ratio groups___Factor(s) defining comparison groups, i.e. groups of samples, including at least some References (i.e. Controls), to compare to each others./nMust include Experiment, cannot include Replicate.",
                      "Batch___Optional: Factor(s) defining batches used for sva::ComBat-based correction.",
                      "GO enrichment___Optional: Only protein groups with at least one valid value in corresponding samples will be used as references for GO-terms enrichment tests."))
for (nm in coreNms) { if (!nm %in% colnames(Param)) { Param[[nm]] <- "Exp" } }
wMp <- c(which(colnames(Param) == coreNms[1]),
         which(colnames(Param) == coreNms[2]),
         which(colnames(Param) == coreNms[3]),
         which(colnames(Param) == coreNms[4]),
         which(colnames(Param) == coreNms[5]),
         which((Param[1,] == "MAP2FACTS")&(!colnames(Param) %in% coreNms)))
lstFct <- list()
dfltFct <- list()
for (w in wMp) {
  Opt <- Factors
  dflt <- dfdflt <- c("Experiment", "Replicate")[(colnames(Param)[w] == "Batch.correction")+1]
  if (!Param[1, w] %in% c("AUTOFACT", "MAP2FACTS")) {
    dflt <- Factors[match(unlist(strsplit(Param[1, w], ";")), names(Factors))]
    if ((length(dflt) == 1)&&(is.na(dflt))) { dflt <- dfdflt }
  }
  if (colnames(Param)[w] == "Ratios.Groups.Ref.Aggregate.Level") {
    dflt <- unique(c(dflt, mnFct))
  }
  if (colnames(Param)[w] == "Ratios.Groups") {
    Opt <- Factors[which(Factors != "Replicate")]
    dflt <- dflt[which(dflt != "Replicate")]
  }
  lbl <- gsub("\\.", " ", colnames(Param)[w])
  if (colnames(Param)[w] %in% coreNms) {
    lbl <- unlist(strsplit(names(coreNms)[match(colnames(Param)[w], coreNms)], "___"))
  } else { lbl <- c(gsub("\\.", " ", colnames(Param)[w]), Param_Help[w]) }
  names(Opt) <- NULL
  names(dflt) <- NULL
  dfltFct[[colnames(Param)[w]]] <- dflt
  blck <- list(list(br()),
               tags$table(
                 tags$tr(width = "80%",
                         tags$td(width = "25%",
                                 div(strong(lbl[1]))),
                         tags$td(width = "55%",
                                 selectInput(colnames(Param)[w],
                                             "",
                                             Opt,
                                             dflt,
                                             TRUE,
                                             TRUE)),
                         #addTooltip(session, colnames(Param)[w], Param_Help[w], "bottom", "hover", list(container = "body"))
                 )
               ))
  lbl2 <- unlist(strsplit(lbl[2], "/n"))
  for (lbl2a in lbl2) { blck <- append(blck, list(span(em(lbl2a)), br())) }
  if (colnames(Param)[w] == "Ratios.Groups") {
    dflt <- Param$Ratios.Groups_Nested
    if (!is.logical(dflt)) { dflt <- WorkFlow != "Regulation" }
    blck <- append(blck, list(radioButtons("IsNested", "Nested design? (i.e. are replicates paired?)",
                                           c(TRUE, FALSE), dflt, TRUE)))
  }
  blck <- append(blck, list(br()))
  lstFct <- append(lstFct, blck)
}
#
tstAdvOpt <- try(sum(file.exists(Param$Custom.PGs, Param$TrueDisc_filter, Param$CRAPome_file)) > 0)
if ("try-error" %in% class(tstAdvOpt)) { tstAdvOpt <- FALSE }
#if (!Param$Param_suppress_UI) {
  appNm <- paste0(dtstNm, " - Parameters")
  ui <- fluidPage(
    useShinyjs(),
    extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
    titlePanel(tag("u", "Parameters"),
               appNm),
    h2(dtstNm), 
    tags$hr(style = "border-color: black;"),
    br(),
    h4(strong("Factors")),
    span("Here we map diverse actions to experimental Factors:"),
    fluidRow(column(4, withSpinner(uiOutput("FactMappings"))),
             column(8, withSpinner(plotlyOutput("PSMsPCA", height = "600px")))),
    br(),
    tags$hr(style = "border-color: black;"),
    h4(strong("Proteins of interest")),
    pickerInput("IntProt", NULL, protHeads, protDeflt, TRUE, width = "600px",
                pickerOptions(title = "Search me",
                              `live-search` = TRUE,
                              actionsBox = TRUE,
                              deselectAllText = "Clear search",
                              showTick = TRUE)),
    br(),
    tags$hr(style = "border-color: black;"),
    h4(strong("Data processing")),
    fluidRow(
      column(3, numericInput("minInt",
                             "Exclude PSMs with intensity lower than...",
                             minInt,
                             0,
                             .Machine$double.xmax,
                             width = "100%")),
      column(3, checkboxInput("Impute", "Impute missing peptides-level values?", Impute, "100%")),
      column(3,
             checkboxInput("Update_Prot_matches", paste0("Update ", names(SearchSoft), "'s original protein-to-peptides assignments?"), Update_Prot_matches, "100%"),
             bsTooltip("Update_Prot_matches",
                       "Checking assignments may result in removal of some identifications. It is nonetheless recommended because we have observed occasional inconsistent peptides-to-protein assignments with some search software.",
                       placement = "right", trigger = "hover", options = list(container = "body")),
             withSpinner(uiOutput("ReloadMatches"))
      )),
    tags$hr(style="border-color: black;"),
    withSpinner(uiOutput("IsobarCorr")),
    withSpinner(uiOutput("Norm")),
    br(),
    # Quantitation
    ## Choice of algorithm + Proteomics ruler
    tags$hr(style="border-color: black;"),
    h4(strong("Protein Groups quantitation")),
    fluidRow(column(2, selectInput("QuantMeth",
                                   "Protein Groups-level quantitation algorithm:",
                                   names(QuantMethods)[1:6],
                                   QMdefnm,
                                   width = "100%")),
             column(2,
                    checkboxInput("ProtRul",
                                  "Apply Proteomic Ruler to estimate copy numbers per cell? (uses signal from all histones as reference; assumes inter-nucleosomal space = 196 bp, do not use if this assumption does not hold!)",
                                  protrul,
                                  "100%"),
                    numericInput("ProtRulNuclL",
                                 "Use inter-nucleosome length = ? (kb)",
                                 ProtRulNuclL,
                                 1,
                                 Inf,
                                 1,
                                 "100%")),
             column(2, selectInput("Prot.Quant.Use",
                                   "Peptides eligible for quantitation:",
                                   c("Unique", "Razor", "All"),
                                   dfltP4Q,
                                   width = "100%")),
             column(2,
                    h5("Use only peptidoforms found in at least how many samples in..."),
                    numericInput("PepFoundInAtLeast",
                                 " -> the whole dataset?", PepFoundInAtLeast, 1, nrow(Exp.map), 1, "100%"),
                    numericInput("PepFoundInAtLeastGrp",
                                 " -> one sample group at least? (overrides parameter above)",
                                 PepFoundInAtLeastGrp,
                                 1,
                                 mxRp,
                                 1,
                                 "100%"))
    ),
    # Note to self: I am for now excluding some methods, because I need to add code to calculate some columns for those, namely ratios.
    # This should be remedied asap, especially since there include such community favourites as IQ (= MaxLFQ) and Top3!!!
    br(),
    tags$hr(style = "border-color: black;"),
    fluidRow(column(2, radioButtons("Clustering", "Clustering method", klustChoices, klustChoices[1], TRUE, "100%"))),
    br(),
    tags$hr(style = "border-color: black;"),
    h4(strong("Statistical testing")),
    h5(strong("T-test(s)")),
    fluidRow(
      column(2,
             radioButtons("TwoSided", "Fold changes: test...", c("Both directions", "Up-only", "Down-only"),
                          TwoSidedDeflt, FALSE, "100%"),
             checkboxInput("Mirror", "Revert fold changes on plots? (default: log fold change = log2(Sample/Reference); revert: log2(Reference/Sample))",
                           Param$Mirror.Ratios, "100%")),
      column(4,
             h5(strong(" -> T-test: select default variant")),
             radioButtons("TtstPval", "", names(pvalue.col), names(pvalue.col)[2], TRUE, "100%"),
             h6(em("(you can change this later after comparing each test's power)")),
             br(),
             h5(strong(" -> F-test")),
             checkboxInput("Ftest", "run?", fTstDflt, "100%"),
             br(),
             uiOutput("sntXprs")),
      column(2,
             h5("Benjamini-Hochberg FDR thresholds"),
             withSpinner(uiOutput("FDR")),
             numericInput("FDR", "", 0.01, 0, 1),
             actionButton("AddFDR", "Add threshold"),
             actionButton("RemvFDR", "Remove threshold")),
      column(2,
             radioButtons("typeOfThresh", "Ratios thresholds: use...", threshOpt, threshDflt, FALSE, "100%"),
             numericInput("RatCont", "Threshold value = ", KontRt, 0, 100, 0.001, width = "100%"),
             selectInput("RatContGrp",
                         "Control ratios are grouped by:",
                         KontGrps,
                         KontGrp,
                         width = "100%"))
    ),
    h5(strong("F-test")),
    fluidRow(column(1, checkboxInput("Ftest", "", fTstDflt, "100%")),
    #          withSpinner(uiOutput("F_test_grps"))
    ),
    br(),
    tags$hr(style="border-color: black;"),
    withSpinner(uiOutput("GO")),
    tags$hr(style="border-color: black;"),
    withSpinner(uiOutput("CytoScape")),
    tags$hr(style="border-color: black;"),
    h4(strong("Post-translational modifications (PTMs)")),
    fluidRow(column(2, pickerInput("PTMsQuant", "Select PTM(s) eligible for use for Protein Groups quantitation:",
                                   Modifs$`Full name`[which(Modifs$Type == "Variable")], ptmDflt1, TRUE,
                                   pickerOptions(title = "Search me",
                                                 `live-search` = TRUE,
                                                 actionsBox = TRUE,
                                                 deselectAllText = "Clear search",
                                                 showTick = TRUE))),
             column(2, pickerInput("PTMsStats", "Select PTM(s) (if any) for which statistical tests will be performed and subtables written:",
                                   Modifs$`Full name`[which(Modifs$Type == "Variable")],
                                   unlist(strsplit(ptmDflt2, ";")), TRUE,
                                   pickerOptions(title = "Search me",
                                                 `live-search` = TRUE,
                                                 actionsBox = TRUE,
                                                 deselectAllText = "Clear search",
                                                 showTick = TRUE))),
             column(2, checkboxInput("PTMsReNorm",
                                     "Re-normalize modified peptides ratios to those of parent Protein Group(s)?",
                                     TRUE, "100%")),
             column(2, radioButtons("NAsReplMethod",
                                    "Some modified peptides do not have a quantified parent protein group to normalize to. Replaced missing values using:",
                                    NAsReplMethods,
                                    NAsReplMethods[NAsReplMeth],
                                    width = "100%"))
    ),
    h4(strong("Output tables")),
    fluidRow(column(2, checkboxInput("Amica", "Write tables compatible with https://bioapps.maxperutzlabs.ac.at/app/amica", Param$Amica, "100%"))),
    br(),
    tags$hr(style="border-color: black;"),
    checkboxInput("AdvOptOn", "Advanced options", tstAdvOpt),
    withSpinner(uiOutput("AdvOpt")),
    br(),
    actionButton("saveBtn", "Save"),
    br()
  )
  server <- function(input, output, session) {
    # Initialize variables to create in main environment
    PARAM <- reactiveVal(Param)
    m4Quant <- reactiveVal(Mod4Quant)
    m2Xclud <- reactiveVal(Mod2Xclud)
    ADVOPT <- reactiveVal(tstAdvOpt)
    #
    # Dynamic UI
    # Map Parameters to Factors
    output$PSMsPCA <- renderPlotly(plot_lyPSMsPCA)
    output$FactMappings <- renderUI({
      # lst <- list()
      # for (w in wMp) {
      #   Opt <- Factors
      #   dflt <- dfdflt <- c("Experiment", "Replicate")[(colnames(Param)[w] == "Batch.correction")+1]
      #   if (!Param[1, w] %in% c("AUTOFACT", "MAP2FACTS")) {
      #     dflt <- Factors[match(unlist(strsplit(Param[1, w], ";")), names(Factors))]
      #     if ((length(dflt) == 1)&&(is.na(dflt))) { dflt <- dfdflt }
      #   }
      #   if (colnames(Param)[w] == "Ratios.Groups.Ref.Aggregate.Level") {
      #     dflt <- unique(c(dflt, mnFct))
      #   }
      #   if (colnames(Param)[w] == "Ratios.Groups") {
      #     Opt <- Factors[which(Factors != "Replicate")]
      #     dflt <- dflt[which(dflt != "Replicate")]
      #   }
      #   lbl <- gsub("\\.", " ", colnames(Param)[w])
      #   if (colnames(Param)[w] %in% coreNms) {
      #     lbl <- unlist(strsplit(names(coreNms)[match(colnames(Param)[w], coreNms)], "___"))
      #   } else { lbl <- c(gsub("\\.", " ", colnames(Param)[w]), Param_Help[w]) }
      #   names(Opt) <- NULL
      #   names(dflt) <- NULL
      #   blck <- list(list(br()),
      #                tags$table(
      #                  tags$tr(width = "80%",
      #                          tags$td(width = "25%",
      #                                  div(strong(lbl[1]))),
      #                          tags$td(width = "55%",
      #                                  selectInput(colnames(Param)[w],
      #                                              "",
      #                                              Opt,
      #                                              dflt,
      #                                              TRUE,
      #                                              TRUE)),
      #                          #addTooltip(session, colnames(Param)[w], Param_Help[w], "bottom", "hover", list(container = "body"))
      #                  )
      #                ))
      #   lbl2 <- unlist(strsplit(lbl[2], "/n"))
      #   for (lbl2a in lbl2) { blck <- append(blck, list(span(em(lbl2a)), br())) }
      #   if (colnames(Param)[w] == "Ratios.Groups") {
      #     dflt <- Param$Ratios.Groups_Nested
      #     if (!is.logical(dflt)) { dflt <- WorkFlow != "Regulation" }
      #     blck <- append(blck, list(radioButtons("IsNested", "Nested design? (i.e. are replicates paired?)",
      #                                            c(TRUE, FALSE), dflt, TRUE)))
      #   }
      #   blck <- append(blck, list(br()))
      #   lst <- append(lst, blck)
      # }
      # return(lst)
      lstFct
    })
    output$RatioGroups <- renderUI({
      em(paste(unique(apply(ExpMap[, input[["Ratios.Groups"]], drop = FALSE], 1, paste, collapse = "___")), collapse = " / "))
    })
    # Factors
    sapply(wMp, function(w) {
      observeEvent(input[[colnames(Param)[w]]], {
        if (colnames(Param)[w] == "Ratios.Groups.Ref.Aggregate.Level") {
          l <- length(input[[colnames(Param)[w]]])
          if ((l < 3)||(sum(!c("Experiment", "Replicate") %in% input[[colnames(Param)[w]]]))) {
            shinyjs::disable("saveBtn")
          } else { shinyjs::enable("saveBtn") }
        }
        Par <- PARAM()
        Par[colnames(Par)[w]] <- paste0(input[[colnames(Param)[w]]], collapse = ";")
        PARAM(Par)
      }, ignoreNULL = FALSE)
    })
    # SAINTexpress
    output$sntXprs <- renderUI({
      lst <- list(list(br()))
      if (WorkFlow %in% c("PULLDOWN", "BIOID")) {
        lst <- list(list(h5(strong(" -> SAINTexpress")),
                         checkboxInput("saintExprs", "run?", saintExprs, "100%"))
        )
      }
      return(lst)
    })
    # CytoScape
    #suppress
    output$CytoScape <- renderUI({
      if (!length(CytoScExe)) {
        lst <- list(list(column(2, shinyFilesButton("CytoScVers2",
                                                    em("No CytoScape.exe detected, select one (optional)"),
                                                    "", FALSE))))
      }
      if (length(CytoScExe) == 1) {
        lst <- list(list(column(2, em(paste0("Version detected: ", CytoScExe)))))
      }
      if (length(CytoScExe) > 1) {
        lst <- list(list(fluidRow(column(2, pickerInput("CytoScVers1", "Choose CytoScape version:",
                                                        CytoScExe, CytoScExe[1], FALSE)))))
      }
      return(lst)
    })
    # Labels purity correction - only for isobarically labelled samples
    output$IsobarCorr <- renderUI({
      if (LabelType == "Isobaric") {
        lst <- list(shinyFilesButton("PurityFl", "Browse", em(paste0(IsobarLab, " purity table")), "", FALSE),
                    br())
      } else {
        lst <- list(HTML(""))
      }
      return(list(lst))
    })
    # Normalisations
    output$Norm <- renderUI({
      lst <- list(list(fluidRow(column(1, h4(strong("Normalisations"))),
                                column(2, checkboxInput("Norm", "Normalize data",
                                                        sum(Param$Norma.Ev.Intens,
                                                            Param$Norma.Pep.Intens,
                                                            Param$Norma.Prot.Ratio) > 0, "100%"))),
                       br(),
                       fluidRow(column(1, h5(strong(" -> PSMs-level: "))),
                                column(3, checkboxInput("evLM", "Levenberg-Marquardt",
                                                        Param$Adv.Norma.Ev.Intens, "100%"))),
                       br(),
                       fluidRow(column(1, h5(strong(" -> Peptidoforms-level: "))),
                                column(2, radioButtons("pepShape", "Corrections",
                                                       c("none", "vsn", "loess"),
                                                       shpDflt, TRUE, "100%")),
                                column(2, checkboxInput("pepLM", "Levenberg-Marquardt",
                                                        Param$Adv.Norma.Pep.Intens, "100%"))),
                       br()))
      if ((LabelType == "Isobaric")&&(length(Iso) > 1)) {
        lst <- append(lst, list(fluidRow(column(2, checkboxInput("IRS", "IRS", Param$Norma.Pep.Intens.IRS, "100%")))))
        dfltChan <- unlist(strsplit(Param$Norma.Pep.Intens.IRS_Ref_channels, ";"))
        dfltChan <- dfltChan[which(dfltChan %in% get(IsobarLab))]
        if (length(dfltChan) != length(Iso)) { dfltChan <- rep(get(IsobarLab)[1], length(Iso)) }
        dfltChan <- reactiveVal(dfltChan)
        for (i in length(Iso)) {
          m <- unique(Exp.map[which(Exp.map$Isobaric.set == Iso[i])],)
          lst <- append(lst, list(fluidRow(column(2, h5(strong(em(paste0(" -> ",
                                                                         paste(dfltChan(), collapse = " / "),
                                                                         ": "))))),
                                           column(2, pickerInput(paste0("intRef", Iso[i]),
                                                                 paste0("- Isobaric set ", Iso[i]),
                                                                 m$Isobaric.label, dfltChan()[i], TRUE,
                                                                 pickerOptions(title = "Search me",
                                                                               `live-search` = TRUE,
                                                                               actionsBox = TRUE,
                                                                               deselectAllText = "Clear search",
                                                                               showTick = TRUE)
                                           )))))
        }
      }
      # lst <- append(lst, list(fluidRow(column(1, h5(strong(" -> Protein Groups-level: "))),
      #                                  column(2, checkboxInput("prtLM", "Levenberg-Marquardt",
      #                                                          Param$Adv.Norma.Prot.Intens, "100%")),
      #                                  br()),
      #                         fluidRow(column(2,
      #                                         pickerInput("prt2PrtSlct", " -> To selected proteins ",
      #                                                     protHeads, prt2PrtSlctDflt, TRUE,
      #                                                     pickerOptions(title = "Search me",
      #                                                                   `live-search` = TRUE,
      #                                                                   actionsBox = TRUE,
      #                                                                   deselectAllText = "Clear search",
      #                                                                   showTick = TRUE)),
      #                                         br())),
      #                         if (WorkFlow == "BIOID") {
      #                           fluidRow(column(2,
      #                                           checkboxInput("prt2Biot", " -> To all biotinylated-proteins: ",
      #                                                         Param$Norma.Prot.Ratio.to.Biot, "100%")))
      #                         }
      # ))
      return(lst)
    })
    #
    # GO
    output$GO <- renderUI({
      lst <- list(list(br()))
      if (Annotate) {
        lst <- list(
          list(fluidRow(h4(strong("\tGO terms enrichment"))),
               fluidRow(column(1, checkboxInput("GOenrich", "GO enrichment", goDflt, "100%")),
                        column(2, pickerInput("GO.tabs", "GO terms of interest", allGO, dftlGO, TRUE,
                                              pickerOptions(title = "Search me",
                                                            `live-search` = TRUE,
                                                            actionsBox = TRUE,
                                                            deselectAllText = "Clear search",
                                                            showTick = TRUE))),
                        column(1, checkboxInput("GO2Int", "Use GO terms to define list of proteins of interest?",
                                                Param$GO.terms.for.proteins.of.interest , "100%"))))
        )
      }
      return(lst)
    })
    #
    # updtFTstUI <- function(reactive = TRUE) {
    #   # Update UI
    #   if (reactive) { FTst <- PARAM()$F.test } else { FTst <- Param$F.test }
    #   return(renderUI({
    #     if (FTst) {
    #       lst <- list()
    #       Opt <- Factors[which(Factors != "Replicate")]
    #       dflt <- dfdflt <- "Experiment"
    #       if (!Param$F.test_within %in% c("AUTOFACT", "MAP2FACTS")) {
    #         dflt <- Factors[match(unlist(strsplit(Param$F.test_within, ";")), names(Factors))]
    #         if ((!length(dflt))||((length(dflt) == 1)&&(is.na(dflt)))) { dflt <- dfdflt }
    #       }
    #       dflt <- dflt[which(dflt != "Replicate")]
    #       lbl <- "Groups within which to perform F-tests"
    #       names(Opt) <- NULL
    #       names(dflt) <- NULL
    #       w <- which(colnames(Param) == "F.test_within")
    #       blck <- list(list(br()),
    #                    tags$table(
    #                      tags$tr(width = "100%",
    #                              tags$td(width = "25%", div(strong("F-test groups"))),
    #                              addTooltip(session, colnames(Param)[w], Param_Help[w], "bottom", "hover", list(container = "body")),
    #                      ),
    #                      tags$tr(width = "100%", tags$td(width = "55%",
    #                                                      selectInput("F.test_within",
    #                                                                  "",
    #                                                                  Opt,
    #                                                                  dflt,
    #                                                                  TRUE,
    #                                                                  TRUE)))
    #                    ))
    #       blck <- append(blck, list(span(em("F-tests will be performed within groups defined by the levels of the factors you choose here.")), br()))
    #       blck <- append(blck, list(br()))
    #       lst <- append(lst, blck)
    #     } else { lst <- list(em("No analysis yet")) }
    #     return(lst)
    #   }))
    # }
    #output$F_test_grps <- updtFTstUI(reactive = FALSE)
    #
    output$ReloadMatches <- renderUI({
      if ("evmatch.RData" %in% list.files(wd)) {
        msg <- "Peptide-to-protein matches backup detected in folder: do you want to reload it?\n"
        lst <- list(list(list(br()),
                         tags$table(
                           tags$tr(width = "100%", tags$td(width = "55%", checkboxInput("Reuse_Prot_matches", msg, TRUE)))
                         )))
      } else {
        em(" ")
      }
    })
    #
    #
    # Event observers
    # Optional input files
    updtOptOn <- function(reactive = TRUE) {
      if (reactive) { tst <- ADVOPT() } else { tst <- tstAdvOpt }
      if (tst) {
        if (reactive) { rgKol <- input[["Ratios.Groups"]] } else { rgKol <- dfltFct[["Ratios.Groups"]] }
        RatGrps <- paste(unique(apply(ExpMap[, rgKol, drop = FALSE], 1, paste, collapse = " ")), collapse = " / ")
        lst <- vector("list", 1)
        lst[[1]] <- list(
          fluidRow(column(12, em("->"), 
                          shinyFilesButton("CustPG", em("Custom Protein Groups"), "", FALSE), br(),
                          em("Allows \"cheating\" with the naive Protein Groups assembly algorithm."), br(),
                          em("Useful when e.g. samples express from a custom construct to which matching peptides should be assigned in priority."), br(),
                          em("Should be a table with two columns: \"Leading protein IDs\" (\";\"-separated) and \"Priority\" (integer)."), br()
                          strong("Current selection = "), span(Param$Custom.PGs, style = "color:blue", .noWS = "outside"),
                          br(),
                          br()
          )),
          fluidRow(column(12, em("->"), 
                          shinyFilesButton("TrueDisc", em("TRUE/FALSE protein discovery filter"), "", FALSE), br(),
                          em("Select TRUE/FALSE protein discovery filter .csv file (useful when you know from another experiment which proteins are contaminants)."), br(),
                          em("Should contain a \"Protein ID\" column and one TRUE/FALSE column per \"Ratio group\" value."),
                          em("Current values: "), em(RatGrps), br(),
                          br(),
                          em(DiscFiltModesHlp[1]), br(),
                          em(DiscFiltModesHlp[2]), br(),
                          em(DiscFiltModesHlp[3]), br()
                          strong("Current selection = "), span(Param$TrueDisc_filter, style = "color:blue", .noWS = "outside"),
                          br(),
                          br()
          )),
          fluidRow(column(12, em("->"),
                          shinyFilesButton("CRAPome", em("CRAPome filter"), "", FALSE), br(),
                          em("CRAPome-like filter: 1 column table of protein accessions to mark as contaminants."), br(),
                          em("Column name = \"Protein ID\" or \"Protein IDs\", use \";\" if including more than one ID per row)."), br()
                          strong("Current selection = "), span(Param$CRAPome_file, style = "color:blue", .noWS = "outside"),
                          br(),
                          br()
          ))
        )
      } else { lst <- list(list(em(""))) }
      renderUI(lst)
    }
    output$AdvOpt <- updtOptOn(FALSE)
    observeEvent(input$AdvOptOn, {
      ADVOPT(input$AdvOptOn)
      output$AdvOpt <- updtOptOn()
    })
    observe({ shinyFileChoose(input, "CustPG", roots = getVolumes(), filetypes = "csv")
      {
        tmp <- input$CustPG
        if ((!is.null(tmp))&&(is.list(tmp))) {
          tmp <- parseFilePaths(getVolumes(), tmp)$datapath
          Par <- PARAM()
          Par$Custom.PGs <- normalizePath(tmp, winslash = "/")
          PARAM(Par)
        }
    }
    })
    observe({ shinyFileChoose(input, "TrueDisc", roots = getVolumes(), filetypes = "csv")
      {
        tmp <- input$TrueDisc
        if ((!is.null(tmp))&&(is.list(tmp))) {
          tmp <- parseFilePaths(getVolumes(), tmp)$datapath
          Par <- PARAM()
          Par$TrueDisc_filter <- normalizePath(tmp, winslash = "/")
          PARAM(Par)
          shinyjs::enable("TrueDiscMode")
        } else { shinyjs::disable("TrueDiscMode") }
    }
    })
    observeEvent(input$TrueDiscMode, {
      Par <- PARAM()
      Par$TrueDisc_filter_mode <- input$TrueDiscMode
      PARAM(Par)
    })
    observe({ shinyFileChoose(input, "CRAPome", roots = getVolumes(), filetypes = "csv" )
      {
        tmp <- input$CRAPome
        if ((!is.null(tmp))&&(is.list(tmp))) {
          tmp <- parseFilePaths(getVolumes(), tmp)$datapath
          Par <- PARAM()
          Par$CRAPome_file <- normalizePath(tmp, winslash = "/")
          PARAM(Par)
        }
    }
    })
    # Minimum intensity
    observeEvent(input$minInt, {
      assign("minInt", as.numeric(input$minInt), envir = .GlobalEnv)
      Par <- PARAM()
      Par$Min.Intensity <- as.numeric(input$minInt)
      PARAM(Par)
    })
    # Impute?
    observeEvent(input$Impute, {
      Par <- PARAM()
      Par$Pep.Impute <- input$Impute
      PARAM(Par)
    })
    # Mirror ratios
    observeEvent(input$Mirror, {
      Par <- PARAM()
      Par$Mirror.Ratios <- input$Mirror
      PARAM(Par)
    })
    # Update PSM-to-Protein matches?
    observeEvent(input$Update_Prot_matches, {
      Par <- PARAM()
      Par$Update_Prot_matches <- input$Update_Prot_matches
      PARAM(Par)
      if (input$Update_Prot_matches) { shinyjs::enable("Reuse_Prot_matches") }
      if (!input$Update_Prot_matches) { shinyjs::disable("Reuse_Prot_matches") }
    })
    observeEvent(input[["Reuse_Prot_matches"]], {
      Par <- PARAM()
      Par$Reuse_Prot_matches <- input$Reuse_Prot_matches
      PARAM(Par)
    })
    # Clustering method
    observeEvent(input$Clustering, {
      assign("KlustMeth", match(input$Clustering, klustChoices), envir = .GlobalEnv)
    })
    # Are analyses Two-sided?
    observeEvent(input$TwoSided, {
      Par <- PARAM()
      Par$Two.sided <- input$TwoSided == "Both directions"
      if (input$TwoSided == "Down-only") {
        shinyjs::disable("Mirror")
        Par$Mirror.Ratios <- TRUE
      }
      if (input$TwoSided != "Down-only") { shinyjs::enable("Mirror") }
      PARAM(Par)
    })
    # Normalisations
    observeEvent(input$Norm, {
      Par <- PARAM()
      Par$Norma.Ev.Intens <- input$Norm
      Par$Norma.Pep.Intens <- input$Norm
      Par$Norma.Prot.Ratio <- input$Norm
      PARAM(Par)
    })
    observeEvent(input$evLM, {
      Par <- PARAM()
      Par$Adv.Norma.Ev.Intens <- input$evLM
      PARAM(Par)
    })
    observeEvent(input$pepShape, {
      Par <- PARAM()
      Par$Norma.Pep.Intens.Shape <- input$pepShape
      PARAM(Par)
    })
    observeEvent(input$pepLM, {
      Par <- PARAM()
      Par$Adv.Norma.Pep.Intens <- input$pepLM
      PARAM(Par)
    })
    if ((LabelType == "Isobaric")&&(length(Iso) > 1)) {
      observeEvent(input$IRS, {
        Par <- PARAM()
        Par$Norma.Pep.Intens.IRS <- input$IRS
        PARAM(Par)
      })
      for (i in 1:length(Iso)) {
        observeEvent(input[[paste0("intRef", Iso[i])]], {
          Par <- PARAM()
          tmp <- dfltChan()
          tmp[i] <- input[[paste0("intRef", Iso[i])]]
          dfltChan(tmp)
          Par$Norma.Pep.Intens.IRS_Ref_channels <- paste(tmp, collapse = ";")
          PARAM(Par)
        }, ignoreNULL = FALSE)
      }
    }
    # observeEvent(input$prtLM, {
    #   Par <- PARAM()
    #   Par$Adv.Norma.Prot.Intens <- input$prtLM
    #   PARAM(Par)
    # })
    # observeEvent(input$prt2PrtSlct, {
    #   Par <- PARAM()
    #   Par$Norma.Prot.Ratio.to.proteins <- paste(db$`Protein ID`[dbOrd][match(input$prt2PrtSlct, protHeads)],
    #                                             collapse = ";")
    #   PARAM(Par)
    # }, ignoreNULL = FALSE)
    # if (WorkFlow == "BIOID") {
    #   observeEvent(input$prt2Biot, {
    #     Par <- PARAM()
    #     Par$Norma.Prot.Ratio.to.Biot <- input$prt2Biot
    #     PARAM(Par)
    #   })
    # }
    # Purity correction
    if (LabelType == "Isobaric") {
      observe({ shinyFileChoose(input, "PurityFl", roots = getVolumes(), filetypes = "csv")
        {
          tmp <- input$PurityFl
          if ((!is.null(tmp))&&(is.list(tmp))) {
            tmp <- parseFilePaths(getVolumes(), tmp)$datapath
            Par <- PARAM()
            Par$Label.Purities.file <- normalizePath(tmp, winslash = "/")
            PARAM(Par)
          }
      }
      })
    }
    # Quantitation
    observeEvent(input$QuantMeth, {
      Par <- PARAM()
      Par$QuantMeth <- QuantMethods[match(input$QuantMeth, names(QuantMethods))]
      PARAM(Par)
    })
    observeEvent(input$ProtRul, {
      Par <- PARAM()
      Par$ProtRul <- input$ProtRul
      PARAM(Par)
    })
    observeEvent(input$ProtRulNuclL, {
      Par <- PARAM()
      Par$ProtRulNuclL <- as.integer(input$ProtRulNuclL)
      PARAM(Par)
    })
    observeEvent(input$Prot.Quant.Use, {
      Par <- PARAM()
      Par$Prot.Quant.Use <- input$Prot.Quant.Use
      PARAM(Par)
    })
    observeEvent(input$PepFoundInAtLeast, {
      Par <- PARAM()
      Par$PepFoundInAtLeast <- as.integer(input$PepFoundInAtLeast)
      PARAM(Par)
    })
    observeEvent(input$PepFoundInAtLeastGrp, {
      Par <- PARAM()
      Par$PepFoundInAtLeastGrp <- as.integer(input$PepFoundInAtLeastGrp)
      PARAM(Par)
    })
    # Type of T-test
    observeEvent(input$TtstPval, {
      Par <- PARAM()
      Par$P.values.type <- input$TtstPval
      PARAM(Par)
    })
    # FDR values
    output$FDR <- renderUI(HTML(Param$BH.FDR.values))
    observeEvent(input$AddFDR, {
      Par <- PARAM()
      tmp <- sort(unique(c(unlist(strsplit(Par$BH.FDR.values, ";")), input$FDR)))
      tmp <- tmp[which((tmp <= 1)&(tmp > 0))]
      tmp <- paste(tmp, collapse = ";") 
      output$FDR <- renderUI(HTML(tmp))
      Par$BH.FDR.values <- tmp
      PARAM(Par)
    })
    observeEvent(input$RemvFDR, {
      Par <- PARAM()
      tmp <- sort(unique(unlist(strsplit(Par$BH.FDR.values, ";"))))
      tmp <- tmp[which(tmp != input$FDR)]
      tmp <- paste(tmp, collapse = ";") 
      output$FDR <- renderUI(HTML(tmp))
      Par$BH.FDR.values <- tmp
      PARAM(Par)
    })
    # Ratio-level thresholds
    observeEvent(input$typeOfThresh, {
      Par <- PARAM()
      Par$Ratios.Thresholds <- input$typeOfThresh
      PARAM(Par)
      assign("ratiosThresh", Par$Ratios.Thresholds, envir = .GlobalEnv)
    })
    observeEvent(input$RatCont, {
      Par <- PARAM()
      assign("KontRt", input$RatCont, envir = .GlobalEnv)
      PARAM(Par)
    })
    observeEvent(input$RatContGrp, {
      Par <- PARAM()
      Par$Ratios.Contaminant.Groups <- input$RatContGrp
      PARAM(Par)
    })
    # Is the design nested?
    observeEvent(input$IsNested, {
      Par <- PARAM()
      Par$Ratios.Groups_Nested <- as.logical(input$IsNested)
      PARAM(Par)
    })
    # Proteins of interest
    observeEvent(input$IntProt, {
      Par <- PARAM()
      Par$Prot.list_pep <- Par$Prot.list <- paste(db$`Protein ID`[dbOrd][match(input$IntProt, protHeads)], collapse = ";")
      PARAM(Par)
    }, ignoreNULL = FALSE)
    # F-test
    observeEvent(input$FTest, {
      Par <- PARAM()
      Par$F.test <- input$FTest
      PARAM(Par)
      #output$F_test_grps <- updtFTstUI()
    })
    # observeEvent(input$F.test_within, {
    #   Par <- PARAM()
    #   Par$F.test_within <- paste(substr(input$F.test_within, 1 , 3), collapse = ";")
    #   PARAM(Par)
    # }, ignoreNULL = FALSE)
    # GO enrichment
    observeEvent(input$GOenrich, {
      Par <- PARAM()
      Par$GO.enrichment <- paste(input$GOenrich, collapse = ";")
      PARAM(Par)
    })
    # GO terms of interest
    if (Annotate) {
      observeEvent(input$GO.tabs, {
        Par <- PARAM()
        tmpGO <- allGO2[match(input$GO.tabs, allGO)]
        Par$GO.tabs <- paste(tmpGO, collapse = ";")
        PARAM(Par)
      }, ignoreNULL = FALSE)
    }
    observeEvent(input$GO2Int, {
      Par <- PARAM()
      Par$GO.terms.for.proteins.of.interest <- input$GO2Int
      PARAM(Par)
    })
    #
    observeEvent(input$Amica, {
      Par <- PARAM()
      Par$Amica <- input$Amica
      PARAM(Par)
    })
    # SAINTexpress
    observeEvent(input$saintExprs, {
      Par <- PARAM()
      assign("saintExprs", Par$saintExpress, envir = .GlobalEnv)
      PARAM(Par)
    })
    # CytoScape
    observeEvent(input$CytoScVers1, {
      Par <- PARAM()
      Par$CytoScapePath <- input$CytoScVers1
      PARAM(Par)
    })
    observe({ shinyFileChoose(input, "CytoScVers2", roots = getVolumes(), filetypes = "exe")
      {
        tmp <- input$CytoScVers2
        if ((!is.null(tmp))&&(is.list(tmp))) {
          tmp <- parseFilePaths(getVolumes(), tmp)$datapath
          Par <- PARAM()
          Par$CytoScapePath <- normalizePath(tmp, winslash = "/")
          PARAM(Par)
        }
      }
    })
    # PTMs to use for PG Quant
    observeEvent(input$PTMsQuant, {
      m4Quant(Modifs$Mark[match(unlist(input$PTMsQuant), Modifs$`Full name`)])
      m2Xclud(set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                           c("Mark", "Where")))
    }, ignoreNULL = FALSE)
    # PTMs to test statistically
    observeEvent(input$PTMsStats, {
      Par <- PARAM()
      Par$PTM.analysis <- paste(input$PTMsStats, collapse = ";")
      PARAM(Par)
    }, ignoreNULL = FALSE)
    # Re-normalize PTM peptides
    observeEvent(input$PTMsReNorm, {
      Par <- PARAM()
      Par$PTM.analysis_Norm <- input$PTMsReNorm
      PARAM(Par)
    })
    observeEvent(input$NAsReplMethod, {
      Par <- PARAM()
      assign("NAsReplMeth", match(input$NAsReplMethod, NAsReplMethods), envir = .GlobalEnv)
      Par$PTM.analysis_NAsReplaceMethod <- NAsReplMeth
      PARAM(Par)
    })
    #
    # Save
    observeEvent(input$saveBtn, {
      Par <- PARAM()
      for (w in wMp) { # Extra verification
        if (nchar(Par[[w]])) { Par[[w]] <- paste(substr(unlist(strsplit(Par[[w]], ";")), 1, 3), collapse = ";") }
      }
      assign("Param", Par, envir = .GlobalEnv)
      assign("Mod4Quant", m4Quant(), envir = .GlobalEnv)
      assign("Mod2Xclud", m2Xclud(), envir = .GlobalEnv)
      stopApp()
    })
    #observeEvent(input$cancel, { stopApp() })
    session$onSessionEnded(function() { stopApp() })
  }
  eval(parse(text = runApp))
  #
  Param$Volcano.plots.Aggregate.Level <- Param_filter(Param$Ratios.Groups.Ref.Aggregate.Level, "Rep")
  Param$Ratios.Ref.Groups <- paste0(Param$Ratios.Groups, c("", ";Rep")[Param$Is])
  Param$Ratios.Plot.split <- "Exp"
  tmp <- unlist(strsplit(Param$Ratios.Groups.Ref.Aggregate.Level, ";"))
  tmp <- c(tmp[which(!tmp %in% c("Exp", "Rep"))], "Rep")
  Param$Ratios.Plot.wrap <- tmp[1]
  Param$Ratios.Plot.colour <- tmp[min(c(2, length(tmp)))]
  Param$Ratios.Contamination.Rates <- KontRt
  if (Param$Ratios.Thresholds == "% of intra-sample group ratios") {
    Param$Ratios.Contamination.Rates <- Param$Ratios.Contamination.Rates/100
  }
  # Apply defaults to non-UI parameters
  g <- grep("^TF_((TRUE)|(FALSE))$", Param[1,])
  for (w in g) { Param[[w]] <- as.logical(gsub("^TF_", "", Param[[w]])) }
  g <- grep("^((TRUE)|(FALSE))$", Param[1,])
  for (w in g) { Param[[w]] <- as.logical(Param[[w]]) }
#}
#if (!Param$Param_suppress_UI) {
  # Statistical tests
  dlg_message("Remember: replace this app with a \"Statistics: contrasts editor\" app!", "ok")
  #
  vpal <- unlist(strsplit(Param$Volcano.plots.Aggregate.Level, ";"))
  # We should select a reference level per ratio group, not ratio reference group!
  #rrg <- unlist(strsplit(Param$Ratios.Ref.Groups, ";"))
  rg <- unlist(strsplit(Param$Ratios.Groups, ";"))
  if (length(unique(Exp.map$Experiment)) == 1) {
    if (length(vpal) > 1)  { vpal <- vpal[which(vpal != "Exp")] }
    #if (length(rrg) > 1)  { rrg <- rrg[which(rrg != "Exp")] }
    if (length(rg) > 1)  { rg <- rg[which(rg != "Exp")] }
  }
  kolVPAL <- Factors[vpal]
  #kolRRG <- Factors[rrg]
  kolRG <- Factors[rg]
  factLevComb1 <- apply(Exp.map[, kolVPAL, drop = FALSE], 1 , paste, collapse = " ")
  #factLevComb2 <- apply(Exp.map[, kolRRG, drop = FALSE], 1 , paste, collapse = " ")
  factLevComb2 <- apply(Exp.map[, kolRG, drop = FALSE], 1 , paste, collapse = " ")
  if (("Reference" %in% colnames(Exp.map))&&("logical" %in% class(Exp.map$Reference))&&(sum(Exp.map$Reference))&&(sum(!Exp.map$Reference))) {
    w <- which(Exp.map$Reference)
    rfLev <- aggregate(factLevComb1[w], list(factLevComb2[w]), function(x) {
      x <- unique(x)
      stopifnot(length(x) == 1)
      return(x)
    })
    a <- aggregate(factLevComb1, list(factLevComb2), function(x) { list(unique(x)) })
    rfLev$y <- a$x[match(rfLev$Group.1, a$Group.1)]
    colnames(rfLev) <- c("Group", "Reference level", "All levels")
  } else {
    rfLev <- data.frame("Group" = unique(factLevComb2))
    rfLev$"Reference level" <- sapply(rfLev$Group, function(x) {
      w <- which(factLevComb2 == x)
      return(factLevComb1[w[1]])
    })
    rfLev$"All levels" <- lapply(rfLev$Group, function(x) {
      unique(factLevComb1[which(factLevComb2 == x)])
    })
  }
  makeSec <- FALSE
  # if (Param$F.test) {
  #   Factors3 <- Factors2[which(sapply(Factors2, function(Fact) {
  #     length(FactorsLevels[[Fact]]) > 1
  #   }))]
  #   if (WorkFlow == "TIMECOURSE") { Factors3 <- unique(c(Factors3, "Time.point")) }
  #   makeSec <- length(Factors3) > 1
  # }
  wdth <- paste0(30*max(c(nchar(unlist(rfLev$`All levels`)), 2)), "px")
  appNm <- paste0(dtstNm, " - Stat-tests")
  ui <- fluidPage(
    useShinyjs(),
    extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
    titlePanel(tag("u", "Statistical test - define comparisons"), 
               appNm),
    h2(dtstNm), 
    br(),
    actionButton("saveBtn", "Save"),
    br(),
    fluidPage(
      h3(" -> T-tests: define reference factors combination per comparisons group."),
      h5(em("Select a single reference level sample group per comparisons group.")),
      h5(em("Samples from this samples group will serve as the denominator of ratios and the controls for all t-tests performed within the comparisons group.")),
      br(),
      #withSpinner(
      DTOutput("refLevels")#)
      ,
      br(),
      br()
    ),
    # if (Param$F.test) {
    #   sidebarLayout(
    #     sidebarPanel(
    #       h3(" -> F-tests: design analyses"),
    #       h5(em("For F-tests, you can define multiple analyses, each with simple and/or double contrasts to different references.")),
    #       h5(em("For each analysis, select reference levels for each Factor relevant to the desired contrasts (comparisons) then add Analysis.")),
    #       h5(em("N.B.: Only combinations of levels corresponding to at least 1 valid sample group are allowed.")),
    #       br(),
    #       actionButton("addAnalysis", "Add new analysis"),
    #       withSpinner(uiOutput("NewAnalysis")),
    #       span(uiOutput("Msg"), style = "color:red"),
    #     ),
    #     mainPanel(
    #       h3("F-test analyses"),
    #       br(),
    #       uiOutput("Analyses")
    #     )
    #   )
    # },
    br()
  )
  #
  NuAn_tmplt <- data.frame(Analysis = "", Primary = "", Primary_Ref = "")
  if (makeSec) { NuAn_tmplt$Secondary <- NuAn_tmplt$Secondary_Ref <- "" }
  server <- function(input, output, session) {
    rfLev2 <- rfLev[, c("Group", "Reference level")]
    rfLev2$"Reference level" <- sapply(1:nrow(rfLev2), function(i) {
      as.character(selectInput(paste0("Ref_", as.character(i)),
                               "",
                               rfLev$`All levels`[[i]],
                               rfLev$`Reference level`[i],
                               selectize = FALSE,
                               width = wdth))
    })
    output$refLevels <- renderDT({ rfLev2 },
                                 FALSE,
                                 escape = FALSE,
                                 selection = "none",
                                 editable = TRUE,
                                 rownames = FALSE,
                                 options = list(dom = 't',
                                                paging = FALSE,
                                                ordering = FALSE),
                                 callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
    # sapply(1:nrow(rfLev2), function(x) {
    #   id <- paste0("Ref_", as.character(x))
    #   observeEvent(input[[id]],
    #                {
    #                  rfLev$`Reference level`[x] <- input[[id]]
    #                  assign("rfLev", rfLev, envir = .GlobalEnv)
    #                })
    # })
    # if (Param$F.test) { # Initialize variables to create in main environment
    #   #
    #   # Parse F-test Parameters
    #   fFact <- fRef <- FALSE
    #   #Fact <- Factors2[which(!Factors2 %in% c("Experiment", "Replicate"))]
    #   #w <- sapply(Fact, function(x) { length(FactorsLevels[[x]]) > 1 })
    #   #Fact[w]
    #   #tmp <- setNames(substr(Fact, 1, 3), Fact)
    #   #Param$F.test_factors <- paste0(paste(tmp[1:2], collapse = "___"), "_;_", tmp[3])
    #   #Param$F.test_factors_ref <- paste0(paste(sapply(names(tmp)[1:2], function(x) { FactorsLevels[[x]][1] }), collapse = "___"), "_;_", FactorsLevels[[names(tmp)[3]]][1])
    #   if ("F.test_factors" %in% colnames(Param)) {
    #     F_Fact <- unique(unlist(strsplit(Param$F.test_factors, "_\\|_")))
    #     F_Fact <- F_Fact[which(F_Fact != "")]
    #     if (length(F_Fact)) {
    #       fFact <- TRUE
    #     }
    #   }
    #   if ("F.test_factors_ref" %in% colnames(Param)) {
    #     F_Ref <- unique(unlist(strsplit(Param$F.test_factors_ref, "_\\|_")))
    #     F_Ref <- F_Ref[which(F_Ref != "")]
    #     if (length(F_Ref)) {
    #       fRef <- TRUE
    #     }
    #   }
    #   if (fFact && fRef) {
    #     F_Fact <- strsplit(F_Fact, "_;_")
    #     stopifnot(max(sapply(F_Fact, length)) <= 3) # We can have at most 3 factors aggregates: factors for primary and secondary contrasts, and blocking factors for nested designs.
    #     F_Fact <- lapply(F_Fact, function(x) { x[1:min(c(2, length(x)))] }) # We will add nesting later for simplicity using the UI. 
    #     F_Fact <- lapply(F_Fact, function(x) {
    #       x <- strsplit(x, "___")
    #       x <- x[1:2]
    #       return(x)
    #     })
    #     F_Fact <- setNames(F_Fact, paste0("Analysis_", 1:length(F_Fact)))
    #     F_Ref <- strsplit(F_Ref, "_;_")
    #     stopifnot(max(sapply(F_Ref, length)) <= 2) # We can have at most 2 here: no blocking factors.
    #     F_Ref <- lapply(F_Ref, function(x) { x[1:min(c(2, length(x)))] })
    #     F_Ref <- lapply(F_Ref, function(x) {
    #       x <- strsplit(x, "___")
    #       x <- x[1:2] 
    #       return(x)
    #     })
    #     F_Ref <- setNames(F_Ref, paste0("Analysis_", 1:length(F_Ref)))
    #     tst <- length(F_Fact) == length(F_Ref)
    #     if (tst) {
    #       tst <- sum(sapply(1:length(F_Fact), function(x) { length(F_Fact[[x]]) != length(F_Ref[[x]]) })) == 0
    #       if (tst) {
    #         tmp <- substr(Factors, 1, 3)
    #         tst <- sum(unlist(sapply(1:length(F_Fact), function(x) {
    #           sapply(1:length(F_Ref[[x]]), function(y) {
    #             if (length(F_Fact[[x]][[y]]) > 0) {
    #               sapply(1:length(F_Fact[[x]][[y]]), function(z) {
    #                 sum(!F_Ref[[x]][[y]][[z]] %in% FactorsLevels[[Factors[match(substr(F_Fact[[x]][[y]][[z]],1 , 3), tmp)]]])
    #               })  
    #             } else { 0 }
    #           })
    #         }))) == 0
    #         if (!tst) { fFact <- fRef <- FALSE }
    #       } else { fFact <- fRef <- FALSE }
    #     } else { fFact <- fRef <- FALSE }
    #   }
    #   if (fFact && fRef) {
    #     nms <- names(F_Fact)
    #     F_Analyses <- data.frame(Analysis = nms, row.names = nms)
    #     F_Analyses$Primary <- lapply(F_Fact, function(x) { x[[1]] })
    #     F_Analyses$Primary_Ref <- lapply(F_Ref, function(x) { x[[1]] })
    #     if (makeSec) {
    #       F_Analyses$Secondary <- lapply(F_Fact, function(x) { x[[2]] })
    #       F_Analyses$Secondary_Ref <- lapply(F_Ref, function(x) { x[[2]] })
    #     }
    #   } else {
    #     F_Analyses <- data.frame(Analysis = "", Primary = "", Primary_Ref = "")
    #     if (makeSec) { F_Analyses$Secondary <- F_Analyses$Secondary_Ref <- "" }
    #     F_Analyses <- F_Analyses[integer(0),]
    #   }
    #   #
    #   F_An <- reactiveVal(F_Analyses)
    #   nr <- reactiveVal(nrow(F_Analyses))
    #   mxNr <- reactiveVal(nrow(F_Analyses))
    #   #
    #   updtAnUI <- function(reactive = TRUE) {
    #     # Update UI
    #     if (reactive) { FA <- F_An() } else { FA <- F_Analyses }
    #     return(renderUI({
    #       if (nrow(FA)) {
    #         txts <- apply(FA, 1, function(r) {
    #           #r <- FA[1,]
    #           pr <- unlist(r[[2]])
    #           prRf <- unlist(r[[3]])
    #           txt <- paste0("-> Primary: ", paste(sapply(1:length(pr), function(x) {
    #             paste0(pr[x], " = ", prRf[x])
    #           }), collapse = ", "))
    #           if (makeSec) {
    #             sc <- unlist(r[[4]])
    #             if (length(sc)) {
    #               scRf <- unlist(r[[5]])
    #               txt <- paste0(txt, "_|_",
    #                             paste0("-> Secondary: ", paste(sapply(1:length(sc), function(x) {
    #                               paste0(sc[x], " = ", scRf[x])
    #                             }), collapse = ", ")))
    #             }
    #           }
    #           return(txt)
    #         })
    #         txts <- strsplit(txts, "_\\|_")
    #         lst <- vector("list", length(txts)*2)
    #         for (r in 1:length(txts)) {
    #           txt <- txts[[r]]
    #           if (length(txt) == 1) {
    #             lst[[r*2-1]] <- list(tags$table(
    #               tags$tr(width = "100%", tags$td(width = "100%", br())),
    #               tags$tr(width = "100%", tags$td(width = "100%", HTML(paste0("Analysis ", r, ":")))),
    #               tags$tr(width = "100%", tags$td(width = "100%", HTML(txt[[1]])))
    #             ))
    #           } 
    #           if (length(txt) == 2) {
    #             lst[[r*2-1]] <- list(tags$table(
    #               tags$tr(width = "100%", tags$td(width = "100%", br())),
    #               tags$tr(width = "100%", tags$td(width = "100%", HTML(paste0("Analysis ", r, ":")))),
    #               tags$tr(width = "100%", tags$td(width = "100%", HTML(txt[[1]]))),
    #               tags$tr(width = "100%", tags$td(width = "100%", HTML(txt[[2]])))
    #             ))
    #           }
    #           lst[[r*2]] <- list(actionButton(paste0("Rmv", r), "Remove"))
    #         }
    #       } else { lst <- list(em("No analysis yet")) }
    #       return(lst)
    #     }))
    #   }
    #   # Create original remove analysis observers
    #   if (nrow(F_Analyses)) {
    #     sapply(1:nrow(F_Analyses), function(r) {
    #       observeEvent(input[[paste0("Rmv", r)]], {
    #         FA <- F_An()
    #         rws <- 1:nrow(FA)
    #         rws <- rws[which(rws != r)]
    #         FA <- FA[rws,]
    #         nr(nrow(FA))
    #         if (nr()) { FA$Analysis <- paste0("Analysis_", 1:nr()) } # Rename analyses
    #         F_An(FA)
    #         output$Analyses <- updtAnUI()
    #       })
    #     })
    #     output$Analyses <- updtAnUI(reactive = FALSE)
    #   }
    #   #
    #   # Box to define a new analysis
    #   dfltsPrim <- reactive(setNames(rep("Not used", length(Factors3)), Factors3))
    #   if (makeSec) { dfltsSec <- reactive(setNames(rep("Not used", length(Factors3)), Factors3)) }
    #   output$NewAnalysis <- renderUI({
    #     tags$table(
    #       tags$tr(width = "100%",
    #               tags$td(width = "20%", div(em(""))),
    #               tags$td(width = "40%", div(em("Primary contrasts"))),
    #               if (makeSec) { tags$td(width = "40%", div(em("Secondary contrasts"))) }
    #       ),
    #       lapply(Factors3, function(Fact){
    #         tags$tr(width = "100%",
    #                 tags$td(width = "20%", div(strong(Fact))),
    #                 tags$td(width = "40%", div(selectInput(paste0("Prim", Fact), "",
    #                                                        c("Not used", FactorsLevels[[Fact]]), dfltsPrim()[[Fact]]))),
    #                 if (makeSec) { tags$td(width = "40%", div(selectInput(paste0("Sec", Fact), "",
    #                                                                       c("Not used", FactorsLevels[[Fact]]),
    #                                                                       dfltsSec()[[Fact]]))) }
    #         )
    #       })
    #     )
    #   })
    #   # Observers for factor reference selection
    #   # These will change the selected value for other factors in the column, if incompatible
    #   EM <- Exp.map
    #   for (Fact in Factors3) { EM[which(is.na(EM[[Fact]])), Fact] <- "NA" }
    #   output$Msg <- renderUI({ em(" ") })
    #   sapply(Factors3, function(Fact) {
    #     sapply(c("Prim", "Sec")[1:(makeSec+1)], function(i) {
    #       observeEvent(input[[paste0(i, Fact)]], {
    #         m <- EM[which(as.character(EM[[Fact]]) == input[[paste0(i, Fact)]]),]
    #         facTrs <- Factors3[which(Factors3 != Fact)]
    #         k <- 0
    #         for (fcT in facTrs) {
    #           if (!input[[paste0(i, fcT)]] %in% c("Not used", m[[fcT]])) {
    #             k <- k+1
    #             updateTextInput(inputId = paste0(i, fcT), value = "Not used")
    #           }
    #         }
    #         output$Msg <- renderUI({ em(c(" ", "Invalid levels combination!")[(k > 0)+1]) })
    #       })
    #     })
    #   })
    #   # Parse new analysis + update UI
    #   # Template new Analysis row
    #   kol <- c("Primary", "Primary_Ref", "Secondary", "Secondary_Ref")[1:(2^(1+makeSec))]
    #   observeEvent(input$addAnalysis, {
    #     Prim <- setNames(sapply(Factors3, function(x) { input[[paste0("Prim", x)]] }), Factors3)
    #     #Prim <- setNames(rep(NA, length(Factors3)), Factors3) ; Prim[Factors3[1]] <- FactorsLevels[[Factors3[1]]][1]
    #     #
    #     Prim <- Prim[which(Prim != "Not used")]
    #     if (length(Prim)) {
    #       NuAn <- NuAn_tmplt
    #       NuAn$Primary <- list(names(Prim))
    #       NuAn$Primary_Ref <- list(Prim)
    #       FA <- F_An()
    #       #FA <- F_Analyses
    #       if (makeSec) {
    #         Sec <- setNames(sapply(Factors3, function(x) { input[[paste0("Sec", x)]] }), Factors3)
    #         #Sec <- setNames(rep(NA, length(Factors3)), Factors3) ; Sec[Factors3[2]] <- FactorsLevels[[Factors3[2]]][1]
    #         Sec <- Sec[which(Sec != "Not used")]
    #         NuAn$Secondary <- list(names(Sec))
    #         NuAn$Secondary_Ref <- list(Sec)
    #       }
    #       tst1 <- apply(NuAn[, kol], 1, function(x) {
    #         x <- unlist(x)
    #         names(x) <- NULL
    #         return(paste(x, collapse = "@"))
    #       })
    #       tst2 <- apply(FA[, kol], 1, function(x) {
    #         x <- unlist(x)
    #         names(x) <- NULL
    #         return(paste(x, collapse = "@"))
    #       })
    #       if (!tst1 %in% tst2) {
    #         FA <- rbind(FA, NuAn)
    #         nr(nrow(FA))
    #         FA$Analysis <- paste0("Analysis_", 1:nr()) #FA$Analysis <- paste0("Analysis_", 1:nrow(FA))
    #         F_An(FA)
    #         # Create new analysis remove observers
    #         if (nr() > mxNr()) {
    #           sapply((mxNr()+1):nr(), function(r) {
    #             observeEvent(input[[paste0("Rmv", r)]], {
    #               FA <- F_An()
    #               rws <- 1:nrow(FA)
    #               rws <- rws[which(rws != r)]
    #               FA <- FA[rws,]
    #               nr(nrow(FA))
    #               if (nr()) { FA$Analysis <- paste0("Analysis_", 1:nr()) } # Rename analyses
    #               F_An(FA)
    #               output$Analyses <- updtAnUI()
    #             })
    #           })
    #           mxNr(max(c(mxNr()), nr()))
    #         }
    #         # Update UI
    #         output$Analyses <- updtAnUI()
    #       }
    #     }
    #   })
    #   #
    # }
    observeEvent(input$saveBtn, {
      # if (Param$F.test) {
      #   F_Analyses <<- F_An()
      #   #print(F_Analyses)
      #   # Update parameters here!
      #   if (nrow(F_Analyses)) {
      #     PAR <- Param
      #     tmpFct <- data.frame(Prim = sapply(F_Analyses$Primary, function(x) { paste(substr(x, 1, 3), collapse = "___") }))
      #     tmpRf <- data.frame(Prim = sapply(F_Analyses$Primary_Ref, paste, collapse = "___"))
      #     if (makeSec) {
      #       tmpFct$Sec <- sapply(F_Analyses$Secondary, function(x) { paste(substr(x, 1, 3), collapse = "___") })
      #       tmpRf$Sec <- sapply(F_Analyses$Secondary_Ref, paste, collapse = "___")
      #     } else {
      #       tmpFct$Sec <- ""
      #       tmpRf$Sec <- ""
      #     }
      #     if (Param$Ratios.Groups_Nested) { tmpFct$Block <- "Rep" }
      #     PAR$F.test_factors <- paste(gsub("_;_$", "", do.call(paste, c(tmpFct, sep = "_;_"))), collapse = "_|_")
      #     PAR$F.test_factors_ref <- paste(gsub("_;_$", "", do.call(paste, c(tmpRf, sep = "_;_"))), collapse = "_|_")
      #     Param <<- PAR
      #   } else { Param$F.test <- FALSE }
      # }
      rfLev$`Reference level` <- sapply(1:nrow(rfLev2), function(x) {
        input[[paste0("Ref_", as.character(x))]]
      })
      assign("rfLev", rfLev, envir = .GlobalEnv)
      stopApp()
    })
    #observeEvent(input$cancel, { stopApp() })
    session$onSessionEnded(function() { stopApp() })
  }
  eval(parse(text = runApp))
  #
  wRf <- lapply(1:nrow(rfLev), function(i) { #i <- 1
    w <- which(factLevComb2 == rfLev$Group[i])
    return(w[which(factLevComb1[w] == rfLev$"Reference level"[i])])
  })
  Exp.map$Reference <- 1:nrow(Exp.map) %in% unlist(wRf)
  write.csv(ExpMap, ExpMapPath, row.names = FALSE)
  tmp <- data.frame(x1 = colnames(Param), x2 = unlist(Param[1,]))
  colnames(tmp) <- NULL
  write.csv(tmp, paste0(wd, "/Parameters.csv"), row.names = FALSE)
#}
# if (Param$Param_suppress_UI) { # Candidate for deletion to simplify: the UI just works!
#   Param <- Param.load()
#   if (("Prot.Quant.Mod" %in% colnames(Param))&&(Param$Prot.Quant.Mod != "")&&(!is.na(Param$Prot.Quant.Mod))) {
#     warning("Parameter Prot.Quant.Mod is deprecated!")
#   }
#   if (("Prot.Quant.Mod.Excl" %in% colnames(Param))&&(Param$Prot.Quant.Mod.Excl != "")&&(!is.na(Param$Prot.Quant.Mod.Excl))) {
#     warning("Parameter Prot.Quant.Mod.Excl is deprecated!")
#   }
#   #Mod4Quant <- unlist(strsplit(Param$Prot.Quant.Mod, split = ";"))
#   ObjNm <- "Mod4Quant"
#   if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
#     msg <- "Which modifications are eligible for protein group quantitation?"
#     mods <- setNames(sapply(Modifs$`Full name`, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") }), Modifs$`Full name`)
#     pre <- mods[which(!grepl("phospho", Modifs$`Full name`, ignore.case = TRUE))]
#     tmp <- dlg_list(mods, pre, multiple = TRUE, title = msg)$res
#     if (length(tmp)) {
#       tmp <- Modifs$Mark[match(names(mods)[match(tmp, mods)], Modifs$`Full name`)]
#     }
#     assign(ObjNm, tmp)
#     AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
#     tmp <- AllAnsw[1,]
#     tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
#     tmp$Value <- list(get(ObjNm))
#     m <- match(ObjNm, AllAnsw$Parameter)
#     if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
#   }
#   # Note that modifications as indirectly inferred from data (when you use the Modifs object)
#   # do not currently for some cases reflect the exact context as defined in the search parameters.
#   # For instance, for "Gln->pyro-Glu" the position is "_" which means any N-terminus,
#   # but the modification as defined in the search is any N-terminal Q.
#   Mod2Xclud <- set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
#                             c("Mark", "Where"))
#   # Clustering method
#   ObjNm <- "KlustMeth"
#   if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
#     msg <- "Which clustering method do you want to use (used at samples and protein groups levels)?"
#     klustChoices <- c("hierarchical                                                                                                                                                                                ",
#                       "K-means                                                                                                                                                                           ")
#     tmp <- dlg_list(klustChoices, klustChoices[1], title = msg)$res
#     if (!length(tmp)) { tmp <- klustChoices[1] }
#     tmp <- c(2, 1)[match(tmp, klustChoices)] # 2,1: not a mistake
#     ObjNm %<c% tmp
#     AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
#     tmp <- AllAnsw[1,]
#     tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
#     tmp$Value <- list(get(ObjNm))
#     m <- match(ObjNm, AllAnsw$Parameter)
#     if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
#   }
# }
Param$FullQuant <- TRUE
FullQuant %<o% Param$FullQuant
protrul %<o% Param$ProtRul
Update_Prot_matches <- Param$Update_Prot_matches
Reuse_Prot_matches <- Param$Reuse_Prot_matches
Pep4QuantOpt %<o% c("Unique peptide IDs", "Razor peptide IDs", "Peptide IDs")
Pep4Quant %<o% Pep4QuantOpt[c(1, 2)[isEukaLike+1]]
if ("Prot.Quant.Use" %in% colnames(Param)) {
  Pep4Quant <- Pep4QuantOpt[which(c("UNIQUE", "RAZOR", "ALL") == toupper(Param$Prot.Quant.Use))]
}
ProtRulNuclL %<o% Param$ProtRulNuclL
CytoScVrs %<o% Param$CytoScapePath
CytoScape %<o% file.exists(CytoScVrs)
Nested %<o% as.logical(Param$Ratios.Groups_Nested)
if (!Nested) {
  Param$Ratios.Ref.Groups <- Param_filter(Param$Ratios.Ref.Groups, "Rep")
} else {
  Param$Ratios.Ref.Groups <- paste(unique(c(unlist(strsplit(Param$Ratios.Ref.Groups, ";")), "Rep")),
                                   collapse = ";")
}
#

# Defaults in case we missed a parameter
g <- grep("^TF_((TRUE)|(FALSE))$", Param[1,])
if (length(g)) {
  Param[1, g] <- gsub("^TF_", "",  Param[1, g])
  Param[, g] <- as.logical(Param[, g])
}
w <- grep("^ *((TRUE)|(FALSE)) *$", Param[1,], ignore.case = TRUE)
if (length(w)) {
  Param[1, w] <- gsub(" ", "", toupper(Param[1, w]))
  Param[, w] <- as.logical(Param[, w])
}

if (SearchSoft %in% c("DIANN", "FRAGPIPE")) {
  tmp <- unlist(strsplit(Param$PSMs, ";"))
  tmp <- tmp[which(file.exists(tmp))]
  tmp <- tmp[which(!tmp %in% PSMsFl)]
  if (length(tmp)) {
    msg <- paste0("Parameter \"PSMs\" is ignored in ", SearchSoft, " mode , instead the PSMs file is detected from the search parameters\n",
                  "The following file(s) will not be included:", paste0("\n - ", tmp), "\n")
    AddMsg2Report(Warning = TRUE, Space = TRUE)
  }
}
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> Parameters:\", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
for (i in 1:ncol(Param)) {
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              paste0("body_add_fpar(Report, fpar(ftext(\"   - ", colnames(Param)[i], ": ",
                                     Param[[i]], "\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))"))
}
ReportCalls <- AddSpace2Report()

# Create sub-directories vector:
dir <- c("Workflow control/MA plots", paste0("Workflow control/", evNm, "s", c("", "/Normalisation")),
         "Workflow control/Peptides", "Workflow control/Peptides/PCA plot", "Workflow control/Peptides/Intensities",
         "Workflow control/Peptides/Intensities/Normalisation - summary", "Workflow control/Peptides/Ratios",
         "Workflow control/Protein groups", "Workflow control/Protein groups/Expression",
         "Workflow control/Protein groups/Ratios", "Workflow control/Protein groups/P-values",
         "Reg. analysis", "Reg. analysis/t-tests",
         #"Reg. analysis/t-tests/pdf",
         #"Reg. analysis/t-tests/jpeg",
         #"Reg. analysis/t-tests/html",
         "PCA plots", "t-SNE plots", "Heatmaps", "Tables")
if (("Norma.Pep.Intens.Shape" %in% colnames(Param))&&(toupper(Param$Norma.Pep.Intens.Shape) %in% c("VSN", "LOESS"))) {
  dir <- c(dir, paste0("Workflow control/Peptides/Intensities/", toupper(Param$Norma.Pep.Intens.Shape), " normalisation"))
}
if (("Norma.Pep.Intens.IRS" %in% colnames(Param))&&(Param$Norma.Pep.Intens.IRS == TRUE)) {
  dir <- c(dir, "Workflow control/Peptides/Intensities/IRS normalisation")
}
if (("Batch.correction" %in% colnames(Param))&&(!as.character(Param$Batch.correction) %in% c("", "F", "FALSE"))) {
  dir <- c(dir, "Workflow control/Peptides/Intensities/Batch correction")
}
enrichGO %<o% (("GO.enrichment" %in% colnames(Param))&&(Param$GO.enrichment))
globalGO %<o% (("GO.enrichment_Whole_dataset" %in% colnames(Param))&&(Param$GO.enrichment_Whole_dataset))
if (enrichGO||globalGO) { dir <- c(dir, "Reg. analysis/GO enrich") }
dirlist <- unique(c(dirlist, paste0(wd, "/", dir)))

# In case we have enriched for a PTM, it helps to check how good the enrichment was:
# This should be before any PSMs are filtered out - we want to look at the data "straight out of the MS"
tstEnrich <- unique(FracMap$`PTM-enriched`)
tstEnrich <- tstEnrich[which((!is.na(tstEnrich))&(tstEnrich != "NA"))]
if (length(tstEnrich)) {
  dir <- paste0(wd, "/Summary plots")
  dirlist <- unique(c(dirlist, dir))
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  for (Mod in tstEnrich) { #Mod <- tstEnrich[1]
    pat <- paste0("\\(", Modifs$Mark[match(Mod, Modifs$`Full name`)], "\\)")
    ev[[Mod]] <- grepl(pat, ev$"Modified sequence")
    for (Type in c("PSMs", "Pep")) { #Type <- "PSMs"
      if (Type == "PSMs") {
        tst <- data.table(Mod = ev[[Mod]],
                          MS_file = ev$"Raw file path",
                          temp = 1)
        tst <- tst[, list(Count = sum(temp)), by = list(MS_file = MS_file, Mod = Mod)]
        tst <- as.data.frame(tst)
        Root <- "PSM"
      }
      if (Type == "Pep") {
        tst <- data.table(Mod = ev[[Mod]],
                          MS_file = ev$"Raw file path",
                          ModSeq = ev$`Modified sequence`)
        tst <- tst[, list(Count = length(unique(ModSeq))), by = list(MS_file = MS_file, Mod = Mod)]
        tst <- as.data.frame(tst)
        Root <- "peptidoform"
      }
      colnames(tst)[2] <- Mod
      m <- match(tst$MS_file, FracMap$"Raw file")
      tst$Sample <- FracMap$MQ.Exp[m]
      tst <- tst[which(!is.na(tst$Sample)),]
      tst <- tst[which(tst$Sample %in% unlist(Exp.map$MQ.Exp)),]
      tst <- tst[which(tst$Sample %in% unlist(Exp.map$MQ.Exp)),]
      #which(sapply(Exp.map$MQ.Exp, function(y) { x %in% unlist(y) }))
      tst2 <- reshape2::melt(tst)
      colnames(tst2) <- gsub("\\(|\\)|\\[|\\]", "", gsub(" ", "_", colnames(tst2)))
      frml <- as.formula(paste0("MS_file ~ `", gsub("\\(|\\)|\\[|\\]", "", gsub(" ", "_", Mod)), "`"))
      tst2 <- cast(tst2, frml, fun.aggregate = sum)
      kN <- paste0("Non-", Mod, "-modified")
      kY <- paste0(Mod, "-modified")
      colnames(tst2)[which(colnames(tst2) == "FALSE")] <- kN
      colnames(tst2)[which(colnames(tst2) == "TRUE")] <- kY
      tst2[[paste0(Mod, " [%]")]] <- signif(100*tst2[[kY]]/(tst2[[kY]]+tst2[[kN]]), 3)
      tst2$Sample <- tst$Sample[match(tst2$MS_file, tst$MS_file)]
      tst2 <- tst2[, c("Sample", "MS_file", kN, kY, paste0(Mod, " [%]"))]
      write.csv(tst2, paste0(dir, "/", Mod, "-", Root, "s per MS file.csv"), row.names = FALSE)
      rw <- unique(tst$MS_file)
      for (i in c(TRUE, FALSE)) {
        w <- which(tst[[Mod]] == i)
        w2 <- which(!rw %in% tst$MS_file[w])
        if (length(w2)) {
          tmp <- tst[which((tst$MS_file %in% rw[w2])&(tst[[Mod]] == !i)),]
          tmp[[Mod]] <- i
          tmp$Count <- 0
          tst <- rbind(tst, tmp)
        }
      }
      tst[[Mod]] <- factor(c("-", "+")[tst[[Mod]]+1], levels = c("-", "+"))
      ttl <- paste0(Mod, "-", Root, "s per MS file")
      plot <- ggplot(tst) + geom_bar(stat = "identity", position = "dodge",
                                     aes(x = MS_file, y = Count, fill = .data[[Mod]])) +
        theme_bw() + scale_fill_viridis(discrete = TRUE) + ggtitle(ttl) +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 5),
              plot.margin = unit(c(0, 0, 0, 3), "in"))
      #poplot(plot, 12, 20)
      ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300)
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)
      ReportCalls <- AddPlot2Report()
    }
  }
}

# Write PTMs table
temp <- Modifs
w <- which(sapply(colnames(Modifs), function(x) { class(Modifs[[x]]) }) == "list")
for (i in w) { temp[[i]] <- sapply(temp[[i]],  paste, collapse = ", ") }
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
write.csv(temp, paste0(dir, "/Modifications.csv"), row.names = FALSE)
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\"PTMs table:\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_table(Report, ReportCalls$Objects$AABiases)")
ReportCalls$Objects$AABiases <- temp
ReportCalls <- AddSpace2Report()

#### Code chunk - Create factor aggregates
a <- Aggregates
names(a) <- NULL
Aggregate.map %<o% data.frame(Aggregate.Name = names(Aggregates), Characteristics = a)
if (length(Aggregates) > 1) {
  temp1 <- c()
  temp2 <- list()
  kount <- 0
  for (i in 2:length(Aggregates)) {
    I <- combn(names(Aggregates), i)
    for (j in 1:ncol(I)) {
      kount <- kount + 1
      J <- I[,j]
      names(J) <- Aggregates[sapply(J, function(x) {which(names(Aggregates) == x)})]
      Exp.map[paste(J, collapse = "")] <- apply(Exp.map[, names(J)], 1, function(x) {
        paste(x, collapse = "___")
      })
      temp1 <- c(temp1, paste(J, collapse = ""))
      temp2[[kount]] <- names(J)
      paste(J, collapse = "") %<c% unique(Exp.map[[paste(J, collapse = "")]])
    }
  }
  temp1 <- data.frame(Aggregate.Name = temp1)
  temp1$Characteristics <- temp2
  Aggregate.map <- rbind(Aggregate.map, temp1)
}
Aggregate.list %<o% sapply(Aggregate.map$Aggregate.Name, function(x) { get(x) })
#This doesn't work for the master/detailed dual script approach (issue with environments!):

# Define reference sample aggregate, as well as ratio groups, ratio ref group and volcano plot aggregates:
if ((!"Norm.Groups" %in% colnames(Param))&&(Param$Norm.Groups == "")) { Param$Norm.Groups <- "Exp" }
Param.aggreg %<o% c()
parse.Param.aggreg.2("Ratios.Groups.Ref.Aggregate.Level", parsed.param.nm = "Ref.Sample.Aggregate")
for (i in c("Ratios.Groups", "Norm.Groups", "Ratios.Ref.Groups", "Volcano.plots.Aggregate.Level", "Ratios.Plot.split", "Ratios.Plot.wrap", "Ratios.Plot.colour")) {
  parse.Param.aggreg.2(i)
}
# Shorter synonyms
RSA %<o% Ref.Sample.Aggregate
VPAL %<o% Volcano.plots.Aggregate.Level
RRG %<o% Ratios.Ref.Groups
RG %<o% Ratios.Groups
#
if (("Adv.Norma.Pep.Ratio" %in% colnames(Param))&&(Param$Adv.Norma.Pep.Ratio != FALSE)) {
  parse.Param.aggreg.2("Adv.Norma.Pep.Ratio.Type.Group")
}
if (("Batch.correction" %in% colnames(Param))&&(!as.character(Param$Batch.correction) %in% c("", "F", "FALSE"))) {
  parse.Param.aggreg.2("Batch.correction")
}
a <- RSA$names
if (length(a) == 1) {
  Exp.map$Ref.Sample.Aggregate <- Exp.map[[a]]
} else { Exp.map$Ref.Sample.Aggregate <- apply(Exp.map[, a], 1, paste, collapse = "___") }
Exp.map <- Exp.map[order(Exp.map[[VPAL$column]], Exp.map$Replicate),]
tst <- setNames(lapply(MQ.Exp, function(x) { which(sapply(Exp.map$MQ.Exp, function(y) { sum(y %in% x) }) > 0) }), MQ.Exp)
tst2 <- sapply(tst, length)
if (LabelType == "LFQ") {
  stopifnot(max(tst2) == 1)
  tst <- tst[which(tst2 == 1)]
}
MQ.Exp <- MQ.Exp[which(MQ.Exp %in% names(tst))]
ev <- ev[which(ev$MQ.Exp %in% names(tst)),]

if (LabelType == "LFQ") {
  ev$RSA <- Exp.map$Ref.Sample.Aggregate[match(ev$MQ.Exp, names(tst))]
  tmp <- ev[, c("Proteins", "Intensity", "RSA")]
} else {
  # Isobaric case - subtle difference
  tmp <- ev[, c("Proteins", "Intensity")]
  tmp$RSA <- Exp.map$Ref.Sample.Aggregate[match(ev$MQ.Exp, names(tst))]
}
# Plot of contamination levels per sample
tmp2 <- listMelt(strsplit(tmp$Proteins, ";"), 1:nrow(tmp), c("Protein", "Row"))
m <- match(tmp2$Protein, db$`Protein ID`)
w <- which(!is.na(m))
tmp2 <- tmp2[w,]; m <- m[w]
# For our purpose here we must match contaminant proteins.
tmp2$Cont <- db$`Potential contaminant`[m]
if (tstorg) {
  tmp2$Organism <- db[m, dbOrgKol]
  tmp2 <- as.data.table(tmp2)
  f0 <- function(x) { c("Target", "Contaminant")[("Contaminant" %in% x)+1] }
  tmp2 <- tmp2[, list(x = f0(Organism)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
} else {
  tmp2 <- as.data.table(tmp2)
  f0 <- function(x) { c("Target", "Contaminant")[("+" %in% x)+1] }
  tmp2 <- tmp2[, list(x = f0(Cont)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
}
tmp$Organism <- tmp2$x[match(1:nrow(tmp), tmp2$Group.1)]
tmp <- tmp[which(!is.na(tmp$Organism)),]
tmp$Organism <- factor(tmp$Organism, levels = c("Contaminant", "Target"))
tmp$Intensity <- as.numeric(tmp$Intensity)
tmp <- aggregate(tmp$Intensity, list(tmp$RSA, tmp$Organism), sum, na.rm = TRUE)
colnames(tmp) <- c("Sample", "Organism", "Total intensity")
tmp$Sample <- factor(cleanNms(tmp$Sample), levels = cleanNms(Exp.map$Ref.Sample.Aggregate))
ttl <- "Contributions to TIC"
plot <- ggplot(tmp) +
  geom_bar(stat = "identity", aes(x = Sample, y = `Total intensity`, fill = Organism)) +
  theme_bw() + scale_fill_viridis(discrete = TRUE, begin = 0.8, end = 0.2) +
  ggtitle(ttl, subtitle = "Summed TIC for each class of identified peptides") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
ggsave(paste0(wd, "/Summary plots/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
ggsave(paste0(wd, "/Summary plots/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")

# Time points
if (exists("Tim")) {
  if (("Time.Points" %in% colnames(Param))&&(Param$Time.Points != "")) {
    tp %<o% as.character(sort(as.numeric(unlist(strsplit(Param$Time.Points, ";")))))
    if (("Time.Point.Names" %in% colnames(Param))&&(Param$Time.Point.Names != "")) {
      names(tp) <- gsub(" ", ".", unlist(strsplit(Param$Time.Point.Names, ";")))
    } else {
      names(tp) <- tp
    }
    if (sum(tp != Tim) > 0) {
      stop("Review your time points!")
    } else {
      Tim <- tp
    }
  } else {
    Tim <- as.character(sort(as.numeric(Tim)))
  }
}

# Refresh list of interesting proteins:
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList.R")
#rstudioapi::documentOpen(Src)
source(Src)
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList2.R")
#rstudioapi::documentOpen(Src)
source(Src)

if ("Prot.list_separate.plots" %in% colnames(Param)) {
  protsplit %<o% Param$Prot.list_separate.plots
  if (protsplit == "") { protsplit <- FALSE }
} else { protsplit %<o% FALSE }
#
if (Param$GO.terms.for.proteins.of.interest) {
  tmpGO <- unlist(strsplit(Param$GO.tabs, ";"))
  GO_prot.list <- list()
  GO_prot.list$Offspring <- lapply(tmpGO, function(x) {
    ont <- Ontology(x)
    x <- c(x, get(paste0("GO", ont, "OFFSPRING"))[[x]])
    x <- x[which(!is.na(x))]
    return(x)
  })
  tmpGO2 <- listMelt(strsplit(db$`GO-ID`, ";"), db$`Protein ID`)
  GO_prot.list$Proteins <- lapply(GO_prot.list$Offspring, function(x) {
    unique(tmpGO2$L1[which(tmpGO2$value %in% unlist(x))])
  })
  prot.list <- unique(c(prot.list, unlist(GO_prot.list$Proteins)))
}
prot.list.Cond %<o% (length(prot.list) > 0)
if (prot.list.Cond) {
  writeFasta(db[match(prot.list, db$`Protein ID`),], intPrtFst)
}

# Similar list as above: proteins for which a heatmap of peptides and a coverage map will be drawn:
if (!exists("prot.list_pep")) { prot.list_pep %<o% c() }
if ("Prot.list_pep" %in% colnames(Param)) {
  tmp <- as.character(Param$Prot.list_pep)
  if ((!as.character(tmp) %in% c("", "NA"))&&(rev(unlist(strsplit(tmp, "\\.")))[1] == "csv")) {
    if (tmp %in% list.files()) {
      tmp <- read.csv(tmp)
      if ("Protein.ID" %in% colnames(tmp)) { tmp <- tmp$Protein.ID } else {
        warning("There is no \"Protein.ID\" column in the file you provided!")
        tmp <- c()
      }
    } else {
      warning("I could not find the protein list for peptides coverage/heatmaps!")
      prot.list_pep <- c() }
  } else {
    if (is.na(tmp)||(tmp == "")) { tmp <- c() } else {
      tmp <- unlist(strsplit(tmp, ";"))
    }
  }
} else { tmp <- c() }
if (exists("TargetProteins")) { tmp <- c(tmp, TargetProteins) }
prot.list_pep %<o% unique(c(prot.list_pep, tmp))
# Filter lists to only keep existing ones
if (prot.list.Cond) {
  prot.list <- gsub("^CON_+", "", prot.list)
  prot.list <- prot.list[which(prot.list %in% db$"Protein ID")]
}
if (length(prot.list_pep)) {
  prot.list_pep <- gsub("^CON_+", "", prot.list_pep)
  prot.list_pep <- prot.list_pep[which(prot.list_pep %in% db$"Protein ID")]
}

# Custom protein groups
custPGs %<o% NA
if (("Custom.PGs" %in% colnames(Param))&&(!as.character(Param$Custom.PGs) %in% c("", "NA", "FALSE"))&&(file.exists(Param$Custom.PGs))) {
  custPGs_fl %<o% Param$Custom.PGs
  custPGs <- read.delim(custPGs_fl, check.names = FALSE)
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.delim(custPGs_fl, check.names = FALSE, sep = ",") }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_fl, check.names = FALSE, sep = ",") }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_fl, check.names = FALSE, sep = "\t") }
  if (colnames(custPGs)[1] != "Leading protein IDs") {
    warning("I could not make sense of the custom protein groups file provided, skipping!")
    custPGs <- NA
  } else {
    prot.list <- unique(c(prot.list, unlist(strsplit(custPGs$"Leading protein IDs", ";"))))
    prot.list_pep <- unique(c(prot.list_pep, unlist(strsplit(custPGs$"Leading protein IDs", ";"))))
  }
}

# Database of user-defined contaminants (not the default ones from CCP or MaxQuant)
## For cases where you searched with a second database of custom contaminants (e.g. E coli for C elegans samlpes)
if (("Cont.DB" %in% colnames(Param))&&(!toupper(as.character(Param$Cont.DB)) %in% c("", " ", "NA", "F", "FALSE"))) {
  temp <- unlist(strsplit(Param$Cont.DB, ";"))
  tst <- file.exists(temp)
  if (sum(!tst)) {
    msg <- paste0("The following contaminant fasta", c("", "s")[(sum(!tst) > 1)+1], "could not be found:", paste0(" - ", temp[which(!tst)], "\n"))
    stop(msg)
    #ReportCalls <- AddMsg2Report(Space = FALSE)
    #temp <- temp[which(tst)]
  }
  temp <- lapply(temp, Format.DB)
  temp <- plyr::rbind.fill(temp)
  w <- which(is.na(temp), arr.ind = TRUE)
  temp[w] <- ""
  temp$Organism_Full[which(temp$Organism_Full == "")] <- "Contaminant"
  temp$Organism[which(temp$Organism == "")] <- "Contaminant"
  temp$"Protein ID" <- paste0("CON__", gsub("^CON__", "",  temp$"Protein ID"))
  temp$"Potential contaminant" <- "+"
  # Remove all evidences which match one of these proteins:
  #test <- strsplit(ev$Proteins, ";")
  #test <- sapply(test, function(x) {sum(x %in% temp$"Protein ID")}) > 0
  #cont.ev %<o% ev[which(test),]
  #ev <- ev[which(!test),]
  if (exists("contDB")) { contDB <- plyr::rbind.fill(list(contDB, temp)) } else { contDB <- temp }
  w <- which(is.na(contDB), arr.ind = TRUE)
  contDB[w] <- ""
  db <- plyr::rbind.fill(list(db, temp))
  w <- which(is.na(db), arr.ind = TRUE)
  db[w] <- ""
}
w <- which(is.na(db$"Protein ID"))
stopifnot(length(w) == 0) # This would need immediate fixing => stop early!

# Write search database as backup
tmp <- rep("", nrow(db)*3)
tmp[3*(1:nrow(db))-2] <- db$Header
tmp[3*(1:nrow(db))-1] <- db$Sequence
write(tmp, paste0(wd, "/Concatenated search database.fasta"))

# Filters
# _______
#
# 1) Optional True/False-Discovery Filter
# Applied later at protein groups level.
# Any protein not in the filter sees its quantitative values set to 0.
#
DiscFilt %<o% ((!is.na(Param$TrueDisc_filter))&(Param$TrueDisc_filter != "")&(file.exists(Param$TrueDisc_filter)))
if (DiscFilt) {
  DiscFiltFl %<o% Param$TrueDisc_filter
  DiscFiltTbl %<o% read.csv(DiscFiltFl, check.names = FALSE)
  tst2 <- sum(!c("Protein ID", RG$values) %in% colnames(DiscFiltTbl)) == 0
  tst3 <- (!sum(!"logical" %in% sapply(RG$values, function(x) { class(DiscFiltTbl[[x]]) })))
  if (tst2+tst3 == 2) {
    DiscFiltMode %<o% Param$TrueDisc_filter_mode
    if (!DiscFiltMode %in% DiscFiltModes) { DiscFiltMode <- DiscFiltModes[1] }
    if (DiscFiltMode == "Filter column") {
      ObjNm <- "DiscFiltCol"
      if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
        msg <- "How should we name the filter column?"
        tmp <- dlg_input(msg, "Found in ...")$res
        ObjNm %<c% tmp
        AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
        tmp <- AllAnsw[1,]
        tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
        tmp$Value <- list(get(ObjNm))
        m <- match(ObjNm, AllAnsw$Parameter)
        if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
      }
    }
  } else {
    if (!tst2) { warning("Invalid filter (missing columns), skipping...") }
    if (!tst3) { warning("Invalid filter (except \"Protein ID\", all otherfilter columns should be logicals), skipping...") }
    DiscFilt <- FALSE
  }
}
#
# 2) Optional "CRAPome" filter
# For pull downs mostly.
# Provides a vector of accessions. Any protein group with one protein in this filter will be marked as a contaminant.
CRAPome %<o% ((!is.na(Param$CRAPome_file))&(Param$CRAPome_file != "")&(file.exists(Param$CRAPome_file)))
if (CRAPome) {
  CRAPomeFl %<o% Param$CRAPome_file
  CRAPomeProteins %<o% read.csv(CRAPomeFl, check.names = FALSE)
  stopifnot(sum(c("Protein ID", "Protein IDs") %in% colnames(CRAPomeProteins)) > 0)
  CRAPomeProteins <- c(CRAPomeProteins$`Protein ID`, CRAPomeProteins$`Protein IDs`)
  CRAPomeProteins <- unique(unlist(strsplit(CRAPomeProteins, ";")))
  if (length(CRAPomeProteins)) {
    w1 <- which(CRAPomeProteins %in% db$`Protein ID`)
    w2 <- which(!CRAPomeProteins %in% db$`Protein ID`)
    if (length(w2)) {
      warning(paste0(w2, " proteins (", round(100*length(w2)/nrow(db), 2), "%) in the provided CRAPome are not in the search database!"))
    }
  } else {
    warning("Empty CRAPome filter, skipping...")
    CRAPome <- FALSE
  }
}

# Pull-down specific parameters
IsPullDown %<o% (gsub(" |_|-|\\.", "", toupper(Param$Type)) %in% c("IP", "IMMUNOPRECIPITATION", "BIOID", "PULLDOWN"))
if (("Two.sided" %in% colnames(Param))&&(is.logical(Param$Two.sided))) {
  TwoSided %<o% Param$Two.sided
} else { TwoSided %<o% !IsPullDown }
# Impute
opt <- setNames(c(TRUE, FALSE), c("Yes", "No"))
dflt <- c("No", "Yes")[IsPullDown+1]
ObjNm <- "Impute"
if ("Pep.Impute" %in% colnames(Param)) { Impute %<o% as.logical(Param$Pep.Impute) } else {
  dflt <- "No"
  if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
    msg <- "Do you want to impute missing peptide values (recommended for pull-down experiments)?"
    tmp <- paste(rep(" ", 180), collapse = "")
    tmp <- opt[gsub(" +$", "", dlg_list(paste0(names(opt), tmp), paste0(dflt, tmp), title = msg)$res)]
    if (is.na(tmp)) { tmp <- FALSE }
    ObjNm %<c% tmp
    AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
    tmp <- AllAnsw[1,]
    tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
    tmp$Value <- list(get(ObjNm))
    m <- match(ObjNm, AllAnsw$Parameter)
    if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
  }
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

# Start writing Materials and Methods
cat("Writing Materials & Methods template...\n")
# Wet lab part
MatMetCalls %<o% list(Calls = list("read_docx()",
                                   paste0("body_add_fpar(MatMet, fpar(ftext(\"Sample", c("", "s")[(nrow(Exp.map) > 1)+1], " preparation\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")),
                      Texts = list())
if (ProcessedByUs) {
  if (LabelType == "LFQ") { WetLabMeth <- try(MatMet_WetLab(), silent = TRUE) }
  if (LabelType == "Isobaric") { WetLabMeth <- try(MatMet_WetLab(Label = IsobarLab), silent = TRUE) }
  if ("try-error" %in% class(WetLabMeth)) { WetLabMeth <- "TEMPLATE" }
} else { WetLabMeth <- "TEMPLATE" }
MatMetCalls$Texts$WetLab <- WetLabMeth
for (i in 1:length(WetLabMeth)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls, paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$WetLab[", i,"], prop = WrdFrmt$",
                                                        c("Body", "Template_text")[(WetLabMeth[i] == "TEMPLATE")+1], "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
# LCMS part
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_fpar(MatMet, fpar(ftext(\"LC-MS/MS analysis\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
LCMSMeth <- try(MatMet_LCMS(), silent = TRUE)
if ((("try-error" %in% class(LCMSMeth))||(is.null(LCMSMeth)))&&("mzML" %in% gsub(".*\\.", "", rawFiles))) {
  tmp <- gsub("\\.mzML", ".raw", rawFiles)
  LCMSMeth <- try(MatMet_LCMS(RawFiles = tmp), silent = TRUE)
}
if ((("try-error" %in% class(LCMSMeth)))||(is.null(LCMSMeth))) { LCMSMeth <- "TEMPLATE" }
L <- length(LCMSMeth)
if (L > 1) {
  LCMSMeth <- lapply(1:length(LCMSMeth), function(x) { c(paste0("Method ", x, ":"), LCMSMeth[[x]]) })
}
LCMSMeth <- unlist(LCMSMeth)
MatMetCalls$Texts$LCMS <- LCMSMeth
for (i in 1:length(LCMSMeth)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls, paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$LCMS[", i,"], prop = WrdFrmt$",
                                                        c("Body", "Template_text")[(LCMSMeth[i] == "TEMPLATE")+1], "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_fpar(MatMet, fpar(ftext(\"Data analysis\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
kol <- c("Organism_Full", "Organism")
kol <- kol[which(kol %in% colnames(db))]
tst <- sapply(kol, function(x) { length(unique(db[which(!as.character(db[[x]]) %in% c("", "NA")), x])) })
kol <- kol[order(tst, decreasing = TRUE)][1]
w <- which(db$`Potential contaminant` != "+")
Org %<o% aggregate(w, list(db[w, kol]), length)
colnames(Org) <- c("Organism", "Count")
Org <- Org[which(Org$Count == max(Org$Count)[1]),]
Org$Source <- aggregate(db$Source[which(db[[kol]] %in% Org$Organism)], list(db[which(db[[kol]] %in% Org$Organism), kol]), function(x) {
  unique(x[which(!is.na(x))])
})$x
Org$Source[which(is.na(Org$Source))] <- ""
tstOrg <- c("", "n")[(tolower(substr(Org$Organism, 1, 1)) %in% c("a", "e", "i", "o", "u"))+1]
if (nrow(Org) == 1) {
  dbTxt <- paste0("a", tstOrg, " ", Org$Organism, " proteome.")
  if (Org$Source != "") { dbTxt <- gsub("\\.$", paste0(" sourced from ", Org$Source, "."), dbTxt) }
} else {
  stop()
  l <- nrow(Org)
  tstSrc <- unique(Org$Source)
  if (length(tstSrc) == 1) {
    dbTxt <- paste0(paste(Org$Organism[1:(l-1)], collapse = ", "), " and ", Org$Organism[l], " proteomes",
                    c(".", paste0(" sourced from ", tstSrc, "."))[(tstSrc != "")+1])
  } else {
    dbTxt <- sapply(1:l, function(x) { paste0(Org$Organism[x], " (", Org$Source[x], ")") })
    dbTxt <- paste0(paste(dbTxt[1:(l-1)], collapse = ", "), " and ", dbTxt[l], " proteomes.")
  }
}
moult <- (length(rawFiles2) > 1)+1
DatAnalysisTxt %<o% "TEMPLATE"
if (SearchSoft == "MAXQUANT") {
  SearchSoftVers <- gsub(" *</?maxQuantVersion> *", "", grep("<maxQuantVersion>", mqpar, value = TRUE))
  gFx <- (grep("<fixedModifications>", mqpar)+1):(grep("</fixedModifications>", mqpar)-1)
  gVar <- (grep("<variableModifications>", mqpar)+1):(grep("</variableModifications>", mqpar)-1)
  FxMd <- gsub(" *</?string> *", "", mqpar[gFx])
  FxMdC <- grep("\\(C\\)", FxMd, value = TRUE)
  FxMd <- FxMd[which(!FxMd %in% FxMdC)]
  VarMd <- gsub(" *</?string> *", "", mqpar[gVar])
  DatAnalysisTxt <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in ",
                           names(SearchSoft))
  if ((length(SearchSoftVers))&&(nchar(SearchSoftVers))) {
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " version ", SearchSoftVers)
  }
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " against ", dbTxt)
  MBR <- gsub(" *</?matchBetweenRuns> *", "", grep("<matchBetweenRuns>", mqpar, value = TRUE))
  SecPep <- gsub(" *</?secondPeptide> *", "", grep("<secondPeptide>", mqpar, value = TRUE))
  DPep <- gsub(" *</?dependentPeptides> *", "", grep("<dependentPeptides>", mqpar, value = TRUE))
  tstFDRs <- grep("Fdr", mqpar, value = TRUE)
  FDRs <- gsub(" *</?[A-Z,a-z]*Fdr[A-Z,a-z]*> *", "", tstFDRs)
  w <- which(!FDRs %in% c("False", "True"))
  FDRs <- setNames(as.numeric(FDRs[w]), gsub(" *</?|>.+", "", tstFDRs[w]))
  FDRs <- FDRs[which(names(FDRs) != "psmFdrCrosslink")]
  if (DPep == "False") { FDRs <- FDRs[which(names(FDRs) != "dependentPeptideFdr")] }
  if (length(FxMdC)) {
    if (length(FxMdC) == 1) {
      txt <- paste0("Fixed cysteine modification was set to ", FxMdC, ".")
    } else {
      txt <- paste0(paste(FxMdC[1:(length(FxMdC)-1)], collapse = ", "), " and ", rev(FxMdC)[1], " were included as fixed cysteine modifications.")
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  }
  if (length(FxMd)) {
    if (length(FxMd) == 1) {
      txt <- paste0("An additional fixed modifications, ", FxMd, ", was also included.")
    } else {
      txt <- paste0("In addition, ", paste(FxMd[1:(length(FxMd)-1)], collapse = ", "), " and ", rev(FxMd)[1], " were included as fixed modifications.")
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  }
  if (length(VarMd)) {
    if (length(VarMd) == 1) {
      txt <- paste0(VarMd, " was set as variable modification.")
    } else {
      txt <- paste0("Variable modifications were set to ", paste(VarMd[1:(length(VarMd)-1)], collapse = ", "), " and ", rev(VarMd)[1], ".")
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  }
  tsts <- setNames(as.logical(toupper(c(MBR, SecPep, DPep))), c("Match-between-Runs", "Second Peptide Search", "Dependent Peptides"))
  y <- which(tsts); ly <- length(y)
  n <- which(!tsts); ln <- length(n)
  if (ly) {
    if (ly == 1) { txt <- paste0(names(tsts)[y], " was set to active") } else {
      txt <- paste0(paste(names(tsts)[y[1:(ly-1)]], collapse = ", "), " and ", names(tsts)[y[ly]], " were set to active")
    }
    if (ln) {
      if (ln == 1) { txt <- paste0(txt, ", while ", names(tsts)[n], " was turned off.") } else {
        txt <- paste0(txt, ", while ", paste(names(tsts)[n[1:(ln-1)]], collapse = ", "), " and ", names(tsts)[n[ln]], " were turned off.")
      }
    } else { txt <- paste0(txt, ".") }
  } else {
    txt <- paste0(paste(names(tsts)[n[1:(ln-1)]], collapse = ", "), " and ", names(tsts)[n[ln]], " were all turned off.")
  }
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  tst <- unique(FDRs)
  lf <- length(FDRs)
  nms <- gsub("Fdr$", "s", names(FDRs))
  if (length(tst) == 1) { txt <- paste0("All FDRs were set to ", tst*100, "%.") } else {
    txt <- paste0("FDRs were set to: ", paste(paste0(FDRs[1:(lf-1)]*100, "% (", nms[1:(lf-1)], ")"), collapse = ", "), " and ", paste0(FDRs[lf]*100, "% (", nms[lf], ")"), ".")
  }
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
}
if (SearchSoft == "DIANN") {
  SearchSoftVers <- gsub("^DIA-NN | .+$", "", grep("^DIA-NN [0-9]+\\.?[0-9]*", DIANNlog, value = TRUE))
  DIANNCall2 <- unlist(strsplit(DIANNCall, " +--"))
  FxMd <- gsub(" enabled as a fixed modification$", "", grep(" enabled as a fixed modification$", DIANNlog, value = TRUE))
  FxMdC <- grep("Cysteine", FxMd, value = TRUE)
  FxMd <- FxMd[which(!FxMd %in% FxMdC)]
  FxMdC <- gsub("cysteine[ ,-]?", "", FxMdC, ignore.case = TRUE)
  FxMdC <- paste0(toupper(substr(FxMdC, 1, 1)), substr(FxMdC, 2, nchar(FxMdC)))
  require(unimod)
  UniMod <- unimod::modifications
  VarMd <- gsub(" will be considered as variable$", "", grep(" will be considered as variable$", DIANNlog, value = TRUE))
  if (length(VarMd)) {
    VarMd <- set_colnames(Isapply(strsplit(gsub("^Modification ", "", VarMd), " with mass delta | at "), unlist), c("UniMod", "Delta mass", "AA"))
    VarMd <- set_colnames(aggregate(VarMd$AA, list(VarMd$UniMod, VarMd$`Delta mass`), list), c("UniMod", "Delta mass", "AA"))
    VarMd$UniMod <-  gsub("^UniMod:", "", VarMd$UniMod)
    VarMd$Type <- "Variable"
    VarMd$"Full name" <- apply(VarMd[, c("UniMod", "Delta mass")], 1, function(x) { #x <- VarMd[1, c("UniMod", "Delta mass")]
      unique(UniMod$Name[which((UniMod$UnimodId == x[[1]])&(round((UniMod$MonoMass == x[[2]])*2)/2 == 0))])
      
    })
    w <- which(lapply(VarMd$`Full name`, length) == 0)
    if (length(w)) {
      VarMd$`Full name`[w] <- apply(VarMd[w, c("Delta mass", "AA")], 1, function(x) {
        pos <- sort(unlist(x[[2]]))
        pos[which(pos == "*n")] <- "_" # Check, this may not be correct...
        # in fact, check whether this part of DIANN_to_MQ() should not be improved.
        Modifs$`Full name`[which((Modifs$`Mass shift` == x[[1]])&(Modifs$AA %in% pos))]
      })
    }
    VarMd$Text <- apply(VarMd[, c("Full name", "AA")], 1, function(x) {
      x[[2]][which(x[[2]] == "n")] <- "N-term"
      x[[2]][which(x[[2]] == "*n")] <- "protein N-term"
      x[[2]][which(x[[2]] == "c")] <- "C-term"
      x[[2]][which(x[[2]] == "*c")] <- "protein C-term"
      paste0(x[[1]], " (", paste(x[[2]], collapse = ""), ")")
    })
  }
  # Note that it is possible that a modification was searched but not found (so is present in VarMd but not Modifs)
  MBR <- grepl("--reanaly[s,z]e", DIANNCall)+1
  Libraries <- gsub("^lib +", "", grep("^lib", DIANNCall2, value = TRUE))
  nL <- length(Libraries)
  LibFree <- (nL == 0)+1
  DatAnalysisTxt <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in ",
                           names(SearchSoft))
  if ((length(SearchSoftVers))&&(nchar(SearchSoftVers))) {
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " version ", SearchSoftVers)
  }
  DatAnalysisTxt <- paste0(DatAnalysisTxt,
                           c(paste0(" against the following spectral librar", c("y", "ies")[(nL > 1)+1], ": ",
                                    paste(Libraries, collapse = "/"), " (see supplementary material), using "),
                             " in library-free mode against ")[LibFree], dbTxt)
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " Match-Between-Runs was turned ", c("on", "off")[MBR], ".")
  if (length(FxMdC)) {
    if (length(FxMdC) == 1) {
      txt <- paste0("Fixed cysteine modification was set to ", FxMdC, ".")
    } else {
      txt <- paste0(paste(FxMdC[1:(length(FxMdC)-1)], collapse = ", "), " and ", rev(FxMdC)[1], " were included as fixed cysteine modifications.")
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  }
  if (length(FxMd)) {
    if (length(FxMd) == 1) {
      txt <- paste0("An additional fixed modifications, ", FxMd, ", was also included.")
    } else {
      txt <- paste0("In addition, ", paste(FxMd[1:(length(FxMd)-1)], collapse = ", "), " and ", rev(FxMd)[1], " were included as fixed modifications.")
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  }
  if (("data.frame" %in% class(VarMd))&&(nrow(VarMd))) {
    if (nrow(VarMd) == 1) {
      txt <- paste0(VarMd$Text, " was set as variable modification.")
    } else {
      txt <- paste0("Variable modifications were set to ", paste(VarMd$Text[1:(nrow(VarMd)-1)], collapse = ", "), " and ", rev(VarMd$Text)[1], ".")
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  }
  tstFDR <- grep("^Output will be filtered at .+ FDR$", DIANNlog, value = TRUE)
  FDR <- gsub("^Output will be filtered at | FDR$", "", tstFDR)
  if (!grepl(" ?%$", FDR)) { FDR <- as.numeric(FDR)*100 }
  if ((FDR >= 0)&&(FDR <= 1)) { DatAnalysisTxt <- paste0(DatAnalysisTxt, " Data was filtered at ", FDR, "% FDR.") }
}
if (SearchSoft == "FRAGPIPE") {
  SearchSoftVers <- gsub("^# FragPipe \\(|\\).+$", "", grep("^# FragPipe \\(", FP_Workflow, value = TRUE))
  FxMdTbl <- gsub(topattern("msfragger.table.fix-mods="), "", grep(topattern("msfragger.table.fix-mods="), FP_Workflow, value = TRUE))
  FxMdTbl <- unlist(strsplit(FxMdTbl, ";"))
  FxMdTbl <- Isapply(strsplit(FxMdTbl, ","), unlist)
  colnames(FxMdTbl) <- c("Shift", "Site", "Enabled", "MaxOccurences")
  FxMdTbl <- FxMdTbl[which(FxMdTbl$Enabled == "true"),]
  FxMdTbl$AAName <- gsub(".+\\(|\\)", "", FxMdTbl$Site)
  FxMdTbl$AAName <- sapply(FxMdTbl$AAName, function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  })
  FxMdTbl$AA <- gsub(" \\(.+\\)", "", FxMdTbl$Site)
  FxMdTbl$Shift <- as.numeric(FxMdTbl$Shift)
  FxMdTbl <- FxMdTbl[which(FxMdTbl$Shift != 0),]
  FxMdTbl$Shift <- as.character(FxMdTbl$Shift)
  g <- grep("^-", FxMdTbl$Shift, invert = TRUE)
  FxMdTbl$Shift[g] <- paste0("+", FxMdTbl$Shift[g])
  FxMdTbl$Name <- apply(FxMdTbl[, c("Shift", "AAName")], 1, function(x) { paste0(x[[1]], " (", x[[2]], ")") })
  FxMdC <- FxMdTbl$Name[which(FxMdTbl$AA == "C")]
  FxMd <- FxMdTbl$Name[which(FxMdTbl$AA != "C")]
  VarMd <- gsub(topattern("msfragger.table.var-mods="), "", grep(topattern("msfragger.table.var-mods="), FP_Workflow, value = TRUE))
  VarMd <- unlist(strsplit(VarMd, ";"))
  VarMd <- Isapply(strsplit(VarMd, ","), unlist)
  colnames(VarMd) <- c("Shift", "Site", "Enabled", "MaxOccurences")
  VarMd <- VarMd[which(VarMd$Enabled == "true"),]
  for (i in 1:nrow(AA_table)) { VarMd$Site <- gsub(AA_table$AA[i], paste0(";", AA_table$AA[i], ";"), VarMd$Site) }
  VarMd$Site <- gsub("n", ";peptide N-term;", VarMd$Site)
  VarMd$Site <- gsub("c", ";peptide C-term;", VarMd$Site)
  VarMd$Site <- gsub("\\[\\^", ";protein N-term;", VarMd$Site)
  VarMd$Site <- gsub("\\^\\]", ";protein C-term;", VarMd$Site)
  VarMd$Site <- strsplit(gsub("^;|;$", "", VarMd$Site), ";;")
  VarMd$AAName <- sapply(VarMd$Site, function(x) { #
    x <- unlist(x)
    w <- which(nchar(x) == 1)
    x[w] <- as.character(AA_table$Name[match(x[w], as.character(AA_table$AA))])
    return(x)
  })
  VarMd$Shift <- gsub(" +", "", VarMd$Shift)
  g <- grep("^-", VarMd$Shift, invert = TRUE)
  VarMd$Shift[g] <- paste0("+", VarMd$Shift[g])
  VarMd$Text <- apply(VarMd[, c("Shift", "AAName", "Site")], 1, function(x) {
    x2 <- x[[2]]
    if (length(x2) > 1) { x2 <- paste(x[[3]], collapse = "") }
    return(paste0(x[[1]], " (", x2, ")"))
  })
  # Note that it is possible that a modification was searched but not found (so is present in VarMd but not Modifs)
  #
  # FragPipe is very versatile... which will not make this easy...
  # Unfortunately, the workflow file does not appear to store the type of Workflow which was preselected in the 2nd tab.
  # OpenSearch is indirectly identified as cases where unusually high values for tolerances are used
  PrTolUp <- gsub(topattern("msfragger.precursor_mass_upper="), "", grep(topattern("msfragger.precursor_mass_upper="), FP_Workflow, value = TRUE))
  PrTolDwn <- gsub(topattern("msfragger.precursor_mass_lower="), "", grep(topattern("msfragger.precursor_mass_lower="), FP_Workflow, value = TRUE))
  # Validations
  Chris <- as.logical(toupper(gsub(topattern("crystalc.run-crystalc="), "", grep(topattern("crystalc.run-crystalc="), FP_Workflow, value = TRUE))))
  Moses <- as.logical(toupper(gsub(topattern("peptide-prophet.run-peptide-prophet="), "", grep(topattern("peptide-prophet.run-peptide-prophet="), FP_Workflow, value = TRUE))))
  Illy <- as.logical(toupper(gsub(topattern("percolator.run-percolator="), "", grep(topattern("percolator.run-percolator="), FP_Workflow, value = TRUE))))
  Maud <- as.logical(toupper(gsub(topattern("ptmprophet.run-ptmprophet="), "", grep(topattern("ptmprophet.run-ptmprophet="), FP_Workflow, value = TRUE))))
  Camus <- as.logical(toupper(gsub(topattern("phi-report.run-report="), "", grep(topattern("phi-report.run-report="), FP_Workflow, value = TRUE))))
  David <- as.logical(toupper(gsub(topattern("ptmshepherd.run-shepherd="), "", grep(topattern("ptmshepherd.run-shepherd="), FP_Workflow, value = TRUE))))
  Kant <- as.logical(toupper(gsub(topattern("quantitation.run-label-free-quant="), "", grep(topattern("quantitation.run-label-free-quant="), FP_Workflow, value = TRUE))))
  Jon <- as.logical(toupper(gsub(topattern("ionquant.run-ionquant="), "", grep(topattern("ionquant.run-ionquant="), FP_Workflow, value = TRUE))))
  Matt <- as.logical(as.numeric(gsub(topattern("ionquant.mbr="), "", grep(topattern("ionquant.mbr="), FP_Workflow, value = TRUE))))
  Fritz <- as.logical(toupper(gsub(topattern("freequant.run-freequant="), "", grep(topattern("freequant.run-freequant="), FP_Workflow, value = TRUE))))
  Timothy <- as.logical(toupper(gsub(topattern("tmtintegrator.run-tmtintegrator="), "", grep(topattern("tmtintegrator.run-tmtintegrator="), FP_Workflow, value = TRUE))))
  #
  DatAnalysisTxt <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in ",
                           names(SearchSoft))
  if ((length(SearchSoftVers))&&(nchar(SearchSoftVers))) {
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " version ", SearchSoftVers)
  }
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " against ", dbTxt)
  if (Ump == 1) {
    DatAnalysisTxt <- gsub("\\.$",
                           " As a preliminary step, the DIAUmpire module was run to extract pseudo-MS/MS spectra from DIA spectra files.",
                           DatAnalysisTxt)
  }
  if (length(FxMdC)) {
    if (length(FxMdC) == 1) {
      txt <- paste0("Fixed cysteine modification was set to ", FxMdC, ".")
    } else {
      txt <- paste0(paste(FxMdC[1:(length(FxMdC)-1)], collapse = ", "), " and ", rev(FxMdC)[1], " were included as fixed cysteine modifications.")
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  }
  if (length(FxMd)) {
    if (length(FxMd) == 1) {
      txt <- paste0("An additional fixed modifications, ", FxMd, ", was also included.")
    } else {
      txt <- paste0("In addition, ", paste(FxMd[1:(length(FxMd)-1)], collapse = ", "), " and ", rev(FxMd)[1], " were included as fixed modifications.")
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  }
  if (nrow(VarMd)) {
    if (nrow(VarMd) == 1) {
      txt <- paste0(VarMd$Text, " was set as variable modification.")
    } else {
      txt <- paste0("Variable modifications were set to ", paste(VarMd$Text[1:(nrow(VarMd)-1)], collapse = ", "), " and ", rev(VarMd$Text)[1], ".")
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", txt)
  }
  w <- which(c(PrTolUp != "20", PrTolDwn != "-20"))
  l <- length(w)
  if (l) {
    PrTol <- c("lower", "upper")[w]
    PrTol[1] <- paste0(toupper(substr(PrTol[1], 1, 1)), substr(PrTol[1], 2 , nchar(PrTol[w])))
    PrTol <- paste(PrTol, collapse = " and ")
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", PrTol, " precursor mass tolerance", c("", "s")[l], " w", c("as", "ere")[l], " set to ",
                             paste(c(PrTolDwn, PrTolUp)[w], collapse = " and "), " ppm", c("", ", respectively")[l], ".")
  }
  # Validations
  w <- which(c(Chris, Moses, Illy)) # Normaly should only ever be length = 1
  l <- length(w)
  if (l) {
    txt <- c("Crystal-C", "PeptideProphet", "Percolator")[w]
    if (l > 1) { txt <- paste0(paste(txt[1:(l-1)], collapse = ", "), " and ", txt[l]) }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " Peptide identifications were validated using ", txt, ".")
    if (Maud) { DatAnalysisTxt <- gsub("\\.$", ", and", DatAnalysisTxt) }
  }
  if (Maud) { DatAnalysisTxt <- paste0(DatAnalysisTxt, " PTM sites", c(" were", "")[(l>0)+1], " localized using PTMProphet.") }
  # FDR levels
  if (Camus) {
    CamusPar <- toupper(gsub(topattern("phi-report.filter="), "", grep(topattern("phi-report.filter="), FP_Workflow, value = TRUE)))
    CamusPar <- unlist(strsplit(CamusPar, "--"))
    CamusPar <- CamusPar[2:length(CamusPar)]
    CamusPar <- data.frame(Par = CamusPar)
    CamusPar$Val <- gsub("[^ ]+ ", "", CamusPar$Par)
    CamusPar$Par <- gsub(" .+", "", CamusPar$Par)
    lvls <- setNames(c("ion", "psm", "peptide", "protein"), c("ION", "PSM", "PEP", "PROT"))
    w <- which(CamusPar$Par %in% names(lvls))
    l <- length(w)
    if (length(w)) {
      txt1 <- lvls[CamusPar$Par[w]]
      txt2 <- 100*as.numeric(CamusPar$Val[w])
      if (l > 1) {
        txt1 <- paste0(paste(txt1[1:(l-1)], collapse = ", "), " and ", txt1[l])
        txt2 <- paste0(paste(txt2[1:(l-1)], collapse = ", "), " and ", txt2[l])
      }
      DatAnalysisTxt <- paste0(DatAnalysisTxt, " Results were filtered in Philosopher at ", txt1, " level", c("", "s")[(l>1)+1],
                               " at FDR", c("", "s")[(l>1)+1], " ", txt2, "%", c("", ", respectively")[(l>1)+1], ".")
    }
  }
  # PTM Shepperd
  if (David) {
    DatAnalysisTxt <- gsub("\\.$", ", then PTM_Shepperd was run to localize Open search identified mass shifts." , DatAnalysisTxt)
    # Basic for now, of course eventually this should be expended.
  }
  # MS1 Quant
  if (Kant) {
    txt <- c()
    if (Jon) { txt <- c(txt, paste0("IonQuant", c("", " with match-between-runs turned on")[Matt+1])) }
    if (Fritz) { txt <- c(txt, "FreeQuant") }
    if (length(txt) > 1) { txt <- paste0("both ", paste(txt, collapse = paste0(" and ", c("", "with ")[Matt+1]))) }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " MS1-level peptide quantitation was performed using ", txt, ".")
  }
  # MS2 Quant
  if (Timothy) { DatAnalysisTxt <- paste0(DatAnalysisTxt, " MS2 reporter intensities were measured using TMT-integrator.") }
  # DiaNN
  if (Diane) {
    Libby <- gsub(topattern("diann.library="), "", grep(topattern("diann.library="), FP_Workflow, value = TRUE))
    Libby <- gsub("\\\\", "", gsub("\\\\\\\\", "/", Libby))
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " DIA quantitation was performed using DiaNN ",
                             c("in library-free mode",
                               paste0("using the following library: \"", Libby, "\" (see supplementary information)"))[(nchar(Libby) > 0)+1], ".")
  }
}
MatMetCalls$Texts$DatAnalysis <- DatAnalysisTxt
mSft <- match(SearchSoft, SearchSoftware)
DatAnalysisTxt <- paste0(names(SearchSoft), "'s output was re-processed using in-house R scripts, starting from the ",
                         c("evidence.txt", "main report", "psm.tsv")[mSft],
                         " table", c("", "s")[(length(PSMsFl)>1)+1], ".")

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
    source(parSrc)
    temp <- ProtMatch2(unique(ev$Sequence), db,
                       cl = parClust) # (ignore the warning for now until we remove contaminant evidences)
  }
  wh1 <- which(ev$Sequence %in% temp$Sequence)
  wh2 <- which(temp$Sequence %in% ev$Sequence)
  mtch1 <- match(ev$Sequence[wh1], temp$Sequence)
  if (!"Proteins" %in% colnames(ev)) { ev$Proteins <- "" } else {
    source(parSrc)
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

# Filter to keep only evidences with valid quantitative values:
if (LabelType == "LFQ") {
  source(parSrc)
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
if ((length(MQ.Exp) > 1)||(LabelType == "Isobaric")) { # Should be always
  source(parSrc)
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
      } else {
        if ((tst[1] >= tst[2]*3)||(tst[1] <= tst[2]/3)) {
          wrpKl <- paste0(paste0("`", X, "`"), "+", paste0("`", Y, "`"))
          tmp <- data2[, c(X, Y)]
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
          )# + scale_alpha_identity()
        if (length(unique(data$`PTM-enrich.`)) == 1) { plot <- plot + guides(colour = "none") }
        #if (1 %in% tst) {
        #  plot <- plot + facet_grid(as.formula(paste0("`PTM-enrich.`~", c(Y, X)[which(tst > 1)])))
        #} else {
        #  plot <- plot + facet_grid(as.formula(paste0("`PTM-enrich.`+", Y, "~", X)))
        #}
      } else {
        plot <- ggplot(data2) + geom_scattermore(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`,
                                                     colour = `Parent sample`), size = 1#, alpha = 0.1
        )
      }
      plot <- plot + scale_color_viridis(begin = 0.25, discrete = TRUE, option = "D") +
        geom_hline(yintercept = 0, colour = "grey") + xlab("A = mean log10(Intensity)") + ylab("M = sample log2(Ratio)") +
        geom_smooth(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`), color = "purple", formula = 'y~s(x, bs = "cs")',
                    linewidth = 0.1, linetype = "dashed", method = 'gam') +
        geom_text(data = annot2, aes(x = Amax, y = Y, label = Tag), hjust = 1, cex = 2) +
        coord_fixed(log10(2)/log2(10)) + theme_bw() + ggtitle(ttl) + ylim(-ylim, ylim) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"), axis.text = element_text(size = 3),
              strip.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 7),
              strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 7))
      if (length(wrpKl) == 1) {
        plot <- plot + facet_wrap(as.formula(paste0("~", "`", wrpKl, "`")), drop = TRUE)
      } else {
        plot <- plot + facet_grid(as.formula(paste0("`", Y, "`~`", X, "`")), drop = TRUE)
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
          source(parSrc)
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
source(Src)

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - Optional - Normalize evidence MS1 intensities, then, if applicable, MS2 reporter (Isobaric labelling) or fragment (DIA) intensities
# Step 0 for DIA measurements
# For DIA we have more accurate estimates of the ratio between samples from measurements of fragments (MS2-based quant)
# We will re-scale total precursor intensities across samples based on the relative amounts of each sample
#Param$Norma.Ev.Intens <- FALSE
if (Param$Norma.Ev.Intens) {
  msg <- paste0(evNm, "-level normalisations:\n------------------------------------\n")
  ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
  # Define groups - this will ensure that, if phospho (or other) -enrichment took place, these peptides will be normalized separately
  ev$"Normalisation group" <- "Standard"
  if ("PTM-enriched" %in% colnames(Frac.map)) {
    # (In column "PTM enriched", use NA to indicate no enrichment!!!)
    if (!"PTM-enriched" %in% colnames(Frac.map)) { Frac.map$"PTM-enriched" <- NA }
    ptmChck <- unique(Frac.map$"PTM-enriched")
    ptmChck <- ptmChck[which(!is.na(ptmChck))]
    if (sum(!ptmChck %in% Modifs$`Full name`)) { stop("Some of the modifications in column \"PTM-enriched\" of Fraction Map are invalid!") }
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
    dir <- paste0(wd, "/Workflow control/", evNm, "/Reporter intensities")
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
      ttl <- paste0(evNm, "s intensity normalisation")
      dir <- paste0(wd, "/Workflow control/", evNm, "/Normalisation")
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
    ttl <- paste0(evNm, "s intensity normalisation")
    dir <- paste0(wd, "/Workflow control/", evNm, "/Normalisation")
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

#### Code chunk - Create modified peptides table
tmp <- as.character(ev$id)
tmp2 <- data.table(id = tmp, mqxp = ev$MQ.Exp, mod = ev$"Modified sequence")
tmp2 <- tmp2[order(ev$id, decreasing = FALSE),]
source(parSrc)
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
      a <- unique(temp$"Modified sequence")
      clusterExport(parClust, list("temp", "j1"), envir = environment())
      b1 <- parLapply(parClust, a, function(x) { temp[which(temp$"Modified sequence" == x), j1] })
      b2 <- parLapply(parClust, a, function(x) { temp$PEP[which(temp$"Modified sequence" == x)] })
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
    tst <- aggregate(ev$"Normalisation group", list(ev$`Modified sequence`), unique)
    if ("character" %in% class(tst$x)) { # It could be that the same sequence is in different normalisation groups - if so we do not want to re-use PSM-level normalisation groups!
      pep$"Normalisation group" <- ev$"Normalisation group"[rvmtch2]
    } else { pep$"Normalisation group" <- "Standard" }
  } else { pep$"Normalisation group" <- "Standard" }
  pep$MQ.Exp <- ev$MQ.Exp[match(pep$`Modified sequence`, ev$`Modified sequence`)]
  nms <- Norm.Groups$names
  w <- which(!nms %in% colnames(pep))
  if (length(w)) { pep[, nms[w]] <- Exp.map[match(pep$MQ.Exp, Exp.map$MQ.Exp), nms[w]] }
  pep$"Normalisation group" <- do.call(paste, c(pep[, c(nms, "Normalisation group")], sep = "_"))
}

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
source(parSrc)
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
    if ("PC3" %in% colnames(scores1)) {
      plot_lyPCA <- plot_ly(scores1, x = ~PC1, y = ~PC2, z = ~PC3,
                            color = ~Samples_group, text = ~Label, type = "scatter3d", mode = "markers",
                            symbol = I(Symb))
      #plot_lyPCA <- add_trace(plot_lyPCA, scores1, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "text", showlegend = FALSE)
    } else {
      plot_lyPCA <- plot_ly(scores1, x = ~PC1, y = ~PC2, color = ~Samples_group, text = ~Label, type = "scatter", mode = "markers",
                            symbol = I(Symb))
      #plot_lyPCA <- add_trace(plot_lyPCA, scores1, x = ~PC1, y = ~PC2, type = "scatter", mode = "text", showlegend = FALSE)
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
#View(Exp.map[, c("Ref.Sample.Aggregate", "Use")])
appNm <- paste0(dtstNm, " - Check for outliers")
#Exp.map$Use <- "TRUE"
fctrs <- Factors[which(sapply(Factors, function(fct) { length(unique(Exp.map[[fct]])) }) > 1)]
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", "Check for outliers"),
             appNm),
  br(),
  h5(tags$div(
    "Check visually how samples are regrouped within expected sample groupings,", tags$br(),
    "then decide whether to remove any outliers (by unticking them in the \"Use\" column),", tags$br(),
    "finally click \"Save\".")),
  br(),
  sidebarLayout(
    sidebarPanel(
      actionButton("saveBtn", "Save"),
      span(uiOutput("Msg"), style = "color:red"),
      br(),
      DTOutput("Include"),
      br()
    ),
    mainPanel(
      plotlyOutput("PCA", height = "800px"),
      br(),
    )
  ))
server <- function(input, output, session) {
  Include <- Exp.map[, c(fctrs, "Use")]
  Include$Use <- shinyCheckInput(Exp.map$Use)
  tstGrps <- aggregate(Exp.map$Use, list(Exp.map[[VPAL$column]]), sum)
  if (min(tstGrps$x) < 2) { shinyjs::disable("saveBtn") }
  if (min(tstGrps$x) >= 2) { shinyjs::enable("saveBtn") }
  output$Include <- renderDT( { Include },
                              FALSE,
                              escape = FALSE,
                              selection = "none",
                              editable = FALSE,
                              rownames = FALSE,
                              options = list(dom = 't',
                                             paging = FALSE,
                                             ordering = FALSE),
                              callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
  output$Msg <- renderUI({ em(" ") })
  output$PCA <- renderPlotly(plot_lyPCA)
  # sapply(seq_len(nrow(Include)), function(x) {
  #   id <- paste0("check___", as.character(x))
  #   observeEvent(input[[id]], {
  #     Exp.map$Use[x] <- as.logical(input[[id]])
  #     assign("Exp.map", Exp.map, envir = .GlobalEnv)
  #     tstGrps <- aggregate(Exp.map$Use, list(Exp.map[[VPAL$column]]), sum)
  #     if (min(tstGrps$x) < 2) {
  #       shinyjs::disable("saveBtn")
  #       output$Msg <- renderUI({ em("Remember: this script requires at least 2 replicates per sample group!!!") })
  #     }
  #     if (min(tstGrps$x) >= 2) {
  #       shinyjs::enable("saveBtn")
  #       output$Msg <- renderUI({ em(" ") })
  #     }
  #   })
  # })
  observeEvent(input$saveBtn, {
    Exp.map$Use <- sapply(seq_len(nrow(Include)), function(x) { input[[paste0("check___", as.character(x))]] })
    tstGrps <- aggregate(Exp.map$Use, list(Exp.map[[VPAL$column]]), sum)
    if (min(tstGrps$x) < 2) {
      #shinyjs::disable("saveBtn")
      output$Msg <- renderUI({ em("Remember: this script requires at least 2 replicates per sample group!!!") })
    }
    if (min(tstGrps$x) >= 2) {
      #shinyjs::enable("saveBtn")
      output$Msg <- renderUI({ em(" ") })
      assign("Exp.map", Exp.map, envir = .GlobalEnv)
      write.csv(Exp.map, ExpMapPath, row.names = FALSE)
      stopApp()
    }
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
eval(parse(text = runApp))
#
#tmp <- read.csv(Param$Experiments.map)
#tmp$Use <- Exp.map$Use
#write.csv(tmp, file = paste0("backup_", Param$Experiments.map), row.names = FALSE)
Exp.map <- Exp.map[which(Exp.map$Use),]
if (LabelType == "Isobaric") { Iso <- sort(unique(Exp.map$Isobaric.set)) }
for (i in 1:nrow(Aggregate.map)) { #i <- 1
  n <- Aggregate.map$Aggregate.Name[i]
  if (nchar(n) == 3) { assign(n, unique(Exp.map[[unlist(Aggregate.map$Characteristics[i])]])) } else {
    assign(n, unique(Exp.map[[n]]))
  }
}
for (i in Param.aggreg) { #i <- Param.aggreg[1]
  a <- get(i)
  if (nchar(a$aggregate) == 3) { a$values <- sort(unique(Exp.map[[a$names]])) } else {
    a$values <- sort(unique(Exp.map[[a$aggregate]]))
  }
  assign(i, a)
}
# Update aliases
RSA <- Ref.Sample.Aggregate
RG <- Ratios.Groups
RRG <- Ratios.Ref.Groups
VPAL <- Volcano.plots.Aggregate.Level

#############################################
# Done!                                     #
#############################################

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
                   size = 0.01, alpha = 0.05) + ggtitle(ttl, subtitle = rfnm) + coord_fixed() + theme_bw() +
  facet_grid(`Strongest in...` ~ `Main charge`) + theme(strip.text.y.right = element_text(angle = 0)) +
  scale_colour_viridis_d(begin = 0.25)
#poplot(plot, 12, 20)
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
    if (shapenormtype == "LOESS") { pack <- "affy" }
    if (shapenormtype == "VSN") { pack <- "vsn" }
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
      geom_scattermore(aes(x = A, y = M, colour = Group), size = 0.5, alpha = 0.1) +
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
      geom_scattermore(aes(x = A, y = M, colour = Group), size = 0.5, alpha = 0.1) +
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
          if ("PC3" %in% colnames(scores1)) {
            plot_lyPCA1 <- plot_ly(scores1, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Iso, text = ~Label,
                                   type = "scatter3d", mode = "markers", symbol = I(Symb))
            plot_lyPCA1 <- add_trace(plot_lyPCA1, scores1, x = ~PC1, y = ~PC2, z = ~PC3,
                                     type = "scatter3d", mode = "text", showlegend = FALSE)
          } else {
            plot_lyPCA1 <- plot_ly(scores1, x = ~PC1, y = ~PC2, color = ~Iso, text = ~Label,
                                   type = "scatter", mode = "markers", symbol = I(Symb))
            plot_lyPCA1 <- add_trace(plot_lyPCA1, scores1, x = ~PC1, y = ~PC2,
                                     type = "scatter", mode = "text", showlegend = FALSE)
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
    source(parSrc)
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
          if ("Original_PC3" %in% colnames(scores2)) {
            plot_lyPCA2 <- plot_ly(scores2, x = ~Original_PC1, y = ~Original_PC2, z = ~Original_PC3, color = scores2$Batch, text = ~Label,
                                   type = "scatter3d", mode = "markers", symbol = I(Symb))
            plot_lyPCA2 <- add_trace(plot_lyPCA2, scores2, x = ~Original_PC1, y = ~Original_PC2, z = ~Original_PC3,
                                     type = "scatter3d", mode = "text", showlegend = FALSE)
          } else {
            plot_lyPCA2 <- plot_ly(scores2, x = ~Original_PC1, y = ~Original_PC2, color = scores2$Batch, text = ~Label,
                                   type = "scatter", mode = "markers", symbol = I(Symb))
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
          if ("PC3" %in% colnames(scores)) {
            plot_lyPCA3 <- plot_ly(scores3, x = ~Corrected_PC1, y = ~Corrected_PC2, z = ~Corrected_PC3, color = scores3$Batch, text = ~Label,
                                   type = "scatter3d", mode = "markers", symbol = I(Symb))
            plot_lyPCA3 <- add_trace(plot_lyPCA3, scores3, x = ~Corrected_PC1, y = ~Corrected_PC2, z = ~Corrected_PC3,
                                     type = "scatter3d", mode = "text", showlegend = FALSE)
          } else {
            plot_lyPCA3 <- plot_ly(scores3, x = ~Corrected_PC1, y = ~Corrected_PC2, color = scores3$Batch, text = ~Label,
                                   type = "scatter", mode = "markers", symbol = I(Symb))
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
                                 options = list(dom = 't',
                                                paging = FALSE,
                                                ordering = FALSE),
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
        eval(parse(text = runApp))
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
                   size = 0.01, alpha = 0.05) + ggtitle(ttl, subtitle = rfnm) + coord_fixed() + theme_bw() +
  scale_color_viridis_d(begin = 0.25) +
  facet_grid(`Strongest in...` ~ `Main charge`) + theme(strip.text.y.right = element_text(angle = 0))
#poplot(plot, 12, 20)
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
#poplot(plot, 12, 20)
print(plot)
ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
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
        tmp <- sapply(xr2, function(y) { log2(pep[[y]]/pep[[xr1]]) })
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
stopCluster(parClust)
source(parSrc)
saveImgFun(BckUpFl)
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
  source(parSrc)
  Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_prepare.R")
  #rstudioapi::documentOpen(Src)
  source(Src)
}

# Define design matrix and contrasts (limma) 
#  From Exp.map to design matrix
Coefficients %<o% Factors[which(!Factors %in% c("Experiment", "Replicate"))]
w <- which(sapply(Coefficients, function(x) { length(unique(FactorsLevels[[x]])) < nrow(Exp.map) }))
Coefficients <- Coefficients[w]
w <- which(sapply(Coefficients, function(x) { length(unique(Exp.map[[x]])) > 1 }))
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
eval(parse(text = contrCall))
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
    if (sum(sapply(RG$names, function(x) {! x %in% RSA$names })) > 0) {
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
}

### Check that CytoScape is installed and can run, then launch it.
Src <- paste0(libPath, "/extdata/R scripts/Sources/Cytoscape_init.R")
#rstudioapi::documentOpen(Src)
source(Src)

#### Code chunk - Modified peptides analysis
PepLabKol %<o% setNames(c("Mod. sequence", "Proteins", "Common protein names", "PEP"),
                        c("Mod. sequence", "Protein(s)", "Protein name(s)", "PEP"))
pep$"Mod. sequence" <- gsub("^_|_$", "", pep$"Modified sequence")
if (("Phospho.analysis" %in% colnames(Param))&&(Param$Phospho.analysis)) {
  warning("This argument is deprecated - use the \"PTM.analysis\" argument now.\nI understand you want to perform phosphopeptides analysis though so will be doing the editing for you... this time!\nCarry on please...")
  if ("PTM.analysis" %in% colnames(Param)) {
    a <- toupper(unlist(strsplit(Param$PTM.analysis, ";")))
    if (!"PHOSPHO" %in% a) { Param$PTM.analysis <- paste(c(a, "PHOSPHO"), collapse = ";") }
  } else { Param$PTM.analysis <- "PHOSPHO" }
}
if ("PTM.analysis" %in% colnames(Param)) {
  PTMs %<o% unlist(strsplit(Param$PTM.analysis, ";"))
  if (length(PTMs)) {
    msg <- "Modified peptides analysis"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    PTMs_ref.ratios %<o% list()
    PTMs_FDR.thresholds %<o% list()
    PTMs_pep %<o% list()
    PTMs_Reg_filters %<o% list()
    PTMs_int.ref %<o% list()
    PTMs_rat.ref %<o% list()
    PTM_normalize_NAsReplaceMethod %<o% list()
    if (F.test) {
      PTMs_F_test_data %<o% list()
      #PTMs_F_test_ref_ratios %<o% list() # Not needed
    }
    if (("PTM.analysis_Norm" %in% colnames(Param))&&(Param$PTM.analysis_Norm != "")) {
      tmp <- as.logical(unlist(strsplit(as.character(Param$PTM.analysis_Norm), ";")))
      if ((!length(tmp) %in% c(1, length(PTMs)))||(NA %in% tmp)) {
        warning("I cannot make sense of parameter PTM.analysis_Norm, defaulting to TRUE")
        tmp <- TRUE
      }
      if ((length(tmp) == 1)&&(length(PTMs) > 1)) {
        warning(paste0("Recycling PTM.analysis_Norm (value = ", as.character(tmp), ") over PTMs..."))
        tmp <- rep(tmp, length(PTMs))
      }
      PTM_normalize %<o% setNames(tmp, PTMs)
    } else { PTM_normalize %<o% setNames(rep(TRUE, length(PTMs)), PTMs) }
    if (enrichGO) {
      PTMs_GO_Plots %<o% list()
      PTMs_Reg_GO_terms %<o% list()
      PTMs_GO_enrich.dat %<o% list()
      PTMs_GO_enrich.FCRt %<o% list()
      PTMs_GO_enrich.tbl %<o% list()
      PTMs_GO_Plots %<o% list()
      PTMs_Reg_GO_terms %<o% list()
      #
      # Initialize ClueGO
      Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_init.R")
      #rstudioapi::documentOpen(Src)
      source(Src)
      #
    }
    tmp <- PTMs; l <- length(tmp)
    if (l > 1) { tmp <- paste0(paste(tmp[1:(l-1)], collapse = ", "), " and ", tmp[l], "") }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " Statistical analysis was performed for ", c("", "PTMs ")[(l > 1)+1], tmp, "-modified peptides.")
    for (ptm in PTMs) { #ptm <- PTMs[1]
      ptms.ref %<o% pep.ref[length(pep.ref)]
      ptms.ratios.ref %<o% setNames(pep.ratios.ref[length(pep.ratios.ref)],
                                    "Original")
      a <- unlist(strsplit(gsub("\\)$", "", ptm), "\\("))
      if (length(a) > 1) {
        Ptm <- paste0(toupper(substr(a[1], 1, 1)), substr(a[1], 2, nchar(ptm)), "(", a[2], ")")
      } else {
        Ptm <- paste0(toupper(substr(a[1], 1, 1)), substr(a[1], 2, nchar(ptm)))
      }
      w <- which(Modifs$"Full name" == Ptm)
      if (length(w) == 1) {
        p <- Modifs$Mark[w]
        ppat <- paste0("\\(", p, "\\)|\\(", p, ",|,", p, "\\)|,", p, ",") # Pattern to catch all instances of the mod
        pep[[Ptm]] <- grepl(ppat, pep$"Modified sequence")
        g <- which(pep[[Ptm]])
        if (length(g)) {
          source(parSrc)
          ReportCalls <- AddMsg2Report(Msg = paste0(" - ", ptm), Space = FALSE)
          dir <- c("", "/t-tests")
          if (F.test) { dir <- c(dir, "/F-tests") }
          dir <-  paste0(wd, "/Reg. analysis/", ptm, dir)
          for (d in dir) { if (!dir.exists(d)) { dir.create(d, recursive = TRUE) }}
          dirlist <- unique(c(dirlist, dir))
          pep[[paste0(Ptm, " ID")]] <- ""
          ptmpep <- pep[g,]
          pep[g, paste0(Ptm, " ID")] <- ptmpep$ModPep_ID <- 1:nrow(ptmpep)
          ptmsh <- substr(p, 1, 1)
          temp <- ptmpep[, c("Modified sequence", "Proteins")]
          temp$"Proteins" <- strsplit(temp$"Proteins", ";")
          temp$"Modified sequence" <- gsub(paste0("[^A-Z", ptmsh, "]"), "",
                                           gsub(ppat, ptmsh, temp$"Modified sequence"))
          #ptmpep[, c("Match(es)", paste0(Ptm, "-site(s)"))] <- ""
          dbsmall <- db[which(db$"Protein ID" %in% unique(unlist(temp$"Proteins"))), c("Protein ID", "Sequence")]
          dbsmall$"Seq*" <- gsub("I", "L", dbsmall$Sequence)
          temp$"ModSeq*" <- gsub("I", "L", temp$"Modified sequence")
          kol <- c("Proteins", "ModSeq*")
          temp$`ModSeq*` <- strsplit(temp$`ModSeq*`, "")
          temp2 <- temp[, kol]
          clusterExport(parClust, list("temp2", "dbsmall", "ptmsh"), envir = environment())
          temp3 <- parApply(parClust, temp2, 1, function(x) {
            #x <- temp[1, kol]
            #A <- sapply(1:nrow(temp), function(i) {
            #print(i)
            #x <- temp[i, kol]
            m <- unlist(x[[2]])
            m <- data.frame(Seq = m, Mod.seq = m, Test = FALSE)
            w1 <- which(m$Seq == ptmsh)
            w2 <- which(m$Seq != ptmsh)
            m$Mod.seq[w1-1] <- paste0(ptmsh, m$Mod.seq[w1-1])
            m$Test[w1-1] <- TRUE
            m <- m[w2,]
            l <- nrow(m)
            m$Offset <- 0:(l-1)
            q <- unlist(x[[1]])
            mtch <- match(q, dbsmall$"Protein ID")
            wN <- which(!is.na(mtch))
            mtch <- mtch[wN]
            q <- dbsmall$"Protein ID"[mtch[wN]]
            if (length(mtch)) {
              seq <- strsplit(dbsmall$"Seq*"[mtch], "")
              matches <- lapply(seq, function(S) { #S <- seq[1]
                S <- unlist(S)
                lS <- length(S)
                m1 <- m
                m1$Match <- apply(m1[, c("Seq", "Offset")], 1, function(y) { which(S == y[1]) - as.numeric(y[2]) })
                M <- unlist(m1$Match)
                M <- M[which(M > 0)]
                M <- aggregate(M, list(M), length)
                M <- M[order(-M$x),]
                M <- M$Group.1[which(M$x == l)]
                # Check that peptides are tryptic:
                #test <- sapply(M, function(y) {
                #  # r1: on the N-terminal end, is the peptide preceded by K, R or (if starting at position 2, M)?
                #  if (y > 1) {
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
              matches <- matches[which(sapply(matches, length) > 0)]
              if (length(matches)) {
                matches <- proteoCraft::listMelt(matches, ColNames = c("Match", "Protein"))
                matches <- aggregate(matches$Protein, list(matches$Match), paste, collapse = ";")
                colnames(matches) <- c("Match", "Proteins")
                w <- which(m$Test)
                matches$Sites <- sapply(matches$Match, function(y) {
                  y <- paste(sapply(w, function(z) {paste0(m$Mod.seq[z], y+m$Offset[z])}), collapse = "-")
                })
                matches$Match <- apply(matches[, c("Match", "Proteins")], 1, paste, collapse = " ")
                matches$Sites <- apply(matches[, c("Sites", "Proteins")], 1, paste, collapse = " ")
                matches <- apply(matches[, c("Match", "Sites")], 2, paste, collapse = "/")
              } else { matches <- c(NA, NA) }
            } else { matches <- c(NA, NA) }
            return(matches)
            #})
          })
          ptmpep[, c("Match(es)", paste0(Ptm, "-site(s)"))] <- as.data.frame(t(temp3))
          ptmpep[[paste0(Ptm, "-site")]] <- gsub(" .+", "", ptmpep[[paste0(Ptm, "-site(s)")]])
          ptmpep <- ptmpep[which(!is.na(ptmpep$`Match(es)`)),]
          ptmpep$tmp1 <- gsub("^_|_$", "", ptmpep$"Modified sequence")
          ptmpep$tmp2 <- ptmpep[[paste0(Ptm, "-site(s)")]]
          nc <- nchar(ptmpep$tmp2)
          w <- which(nc > 25)
          ptmpep$tmp2[w] <- paste0(substr(ptmpep$tmp2[w], 1, 22), "...")
          ptmpep$Code <- apply(ptmpep[, paste0("tmp", 1:2)], 1, paste, collapse = "\n")
          ptmpep$tmp1 <- NULL
          ptmpep$tmp2 <- NULL
          ptmpep$Name <- ""
          w <- which(ptmpep$"Proteins" != "")
          temp2 <- as.data.frame(t(sapply(strsplit(gsub("[/,;].+$", "", ptmpep[w, paste0(Ptm, "-site(s)")]), " "), unlist)))
          colnames(temp2) <- c("Site", "Protein")
          temp2$Protein <- db$"Common Name"[match(temp2$Protein, db$"Protein ID")]
          temp2$ModSeq <- gsub("^_|_$", "", ptmpep$`Modified sequence`[w])
          ptmpep$Name[w] <- do.call(paste, c(temp2[, c("Site", "ModSeq", "Protein")], sep = "\n"))
          ptmpep$Name[which(ptmpep$Name == "")] <- paste0("Unknown source ", ptm, "-modified peptide #", 1:length(which(ptmpep$Name == "")))
          #View(ptmpep[,c("Match(es)", "Modified sequence", "Code", paste0(Ptm, "-site(s)"))])
          if (grepl("^[P,p]hospho( \\([A-Z]+\\))?$", ptm)) {
            p_col <- paste0(gsub(" |\\(|\\)", ".", ptm), ".Probabilities")
            scd_col <- paste0(gsub(" |\\(|\\)", ".", ptm), ".Score.Diffs")
            if (sum(c(p_col, scd_col) %in% colnames(ptmpep)) == 2) {
              temp <- try(phos_QC(ptmpep, P_col = p_col, ScD_col = scd_col), silent = TRUE)
              if ("try-error" %in% class(temp)) {
                warning("No phospho QC performed, check colnames or the phos_QC function!")
              } else { ptmpep$High_Quality_ptmpep <- temp }
            }
          }
          # Convert intensities to log10
          # I've had it! I should've done this way earlier
          ref <- ptms.ref[1]
          ref2 <- paste0("log10(", gsub(" - $", ") - ", gsub("Evidence intensities", "PSMs. int.", ref)))
          kols <- grep(topattern(ref), colnames(ptmpep), value = TRUE)
          for (kol in kols) {
            smpl <- gsub(topattern(ref), "", kol)
            ptmpep[[paste0(ref2, smpl)]] <- log10(ptmpep[[kol]])
            ptmpep[[kol]] <- NULL
          }
          ptms.ref <- c("Original" = ref2)
          # Optional: normalize to parent protein group
          #
          pepRf <- ptms.ref[length(ptms.ref)]
          pepRatRf <- ptms.ratios.ref[length(ptms.ratios.ref)]
          #View(ptmpep[, grep(topattern(pepRf), colnames(ptmpep), value = TRUE)])
          #View(ptmpep[, grep(topattern(pepRatRf), colnames(ptmpep), value = TRUE)])
          #
          # Calculate average intensities and ratios, as well as Welch's t-test and moderated P-values;
          # For unpaired replicates a permutations t-test is also performed.
          samDir <- paste0(dir[1], "/", samSubDir)
          ebamDir <- paste0(dir[1], "/", ebamSubDir)
          dataType <- "modPeptides"
          #
          Src <- paste0(libPath, "/extdata/R scripts/Sources/Av_and_Stat_tests.R")
          #rstudioapi::documentOpen(Src)
          source(Src)
          #
          #kol <- grep(topattern(paste0("Mean ", pepRf)), colnames(ptmpep), value = TRUE)
          #tst <- apply(ptmpep[, kol], 2, function(x) { summary(is.all.good(log10(x))) })
          #kol2 <- grep(topattern(paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)])), colnames(ptmpep), value = TRUE)
          #tst2 <- apply(ptmpep[, kol2], 2, function(x) { summary(is.all.good(x)) })
          #
          # Also mean expression over whole dataset
          a <- grep(topattern(pepRf), colnames(ptmpep), value = TRUE)
          ptmpep$"Mean Expr." <- apply(ptmpep[, a], 1, function(x) { mean(is.all.good(unlist(x))) })
          # Create list of control ratio values for the purpose of identifying thresholds for plots:
          #
          if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
            ref.rat <- NULL
          }
          if (Param$Ratios.Thresholds == threshMsg) {
            PTMs_ref.ratios[[Ptm]] <- ref.rat <- setNames(lapply(VPAL$values, function(x) { #x <- VPAL$values[1]
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
              x <- grep(paste0(topattern(paste0(ptms.ratios.ref[length(ptms.ratios.ref)], x1, "_REF.to.REF_")), "[0-9]+"),
                        colnames(ptmpep), value = TRUE)
              if (length(x)) { x <- is.all.good(as.numeric(unlist(ptmpep[, x]))) } else { x <- NULL }
              return(x)
            }), VPAL$values)
          }
          #
          # Estimate P-value significance for a set of accepted FDRs:
          ## NB: For graphical reasons (volcano plots), there is only support for 4 different FDR values. This should suffice anyway.
          temp_thrsh <- c()
          A <- VPAL$values
          test <- sapply(A, function(x) { #x <- A[6]
            x <- paste0(pvalue.col[which(pvalue.use)], x)
            r <- x %in% colnames(ptmpep)
            if (r) { r <- length(is.all.good(as.numeric(ptmpep[[x]]))) > 0 }
            return(r)
          })
          A <- A[which(test)]
          for (a in A) { #a <- A[1]
            temp <- FDR(data = ptmpep,
                        aggregate = a,
                        pvalue_root = pvalue.col[which(pvalue.use)],
                        fdr = BH.FDR, returns = c(TRUE, TRUE), method = "BH")
            ptmpep[, colnames(temp$`Significance vector`)] <- temp$`Significance vector`
            temp_thrsh <- c(temp_thrsh, temp$Thresholds)
            PTMs_FDR.thresholds[[Ptm]] <- temp_thrsh
          }
          ptmpep$"1-PEP" <- 1 - ptmpep$PEP
          ptmpep$"log10(1-PEP)" <- log10(ptmpep$"1-PEP")
          a <- grep(topattern(pepRf), colnames(ptmpep), value = TRUE)
          a <- a[which(!grepl("\\.REF$", a))]
          ptmpep$"Av. log10 abundance" <- apply(ptmpep[, a], 1, function(x) { mean(is.all.good(unlist(x))) })
          ptmpep$"Rel. av. log10 abundance" <- ptmpep$"Av. log10 abundance"/max(is.all.good(ptmpep$"Av. log10 abundance"))
          # Arbitrary thresholds
          P <- Param
          P$Plot.labels <- "Name"
          P$Plot.metrics <- paste0("X:", paste0(ptms.ratios.ref[length(ptms.ratios.ref)], "Mean."),
                                   ";Y:", pvalue.col[which(pvalue.use)])
          ReportCalls <- AddMsg2Report(Msg = " - t-tests", Space = FALSE)
          #k1 <- grep(topattern(paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)])), colnames(ptmpep), value = TRUE)
          #df1 <- ptmpep[, k1]
          subDr <- gsub(topattern(wd), "", dir[2])
          tempVPptm <- Volcano.plot(ptmpep, "custom",
                                    paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)]),
                                    pvalue.col[which(pvalue.use)],
                                    experiments.map = Exp.map,
                                    aggregate.map = Aggregate.map,
                                    aggregate.list = Aggregate.list,
                                    aggregate.name = VPAL$aggregate,
                                    parameters = P,
                                    save = c("jpeg", "pdf"),
                                    labels = "FDR",
                                    Ref.Ratio.values = ref.rat,
                                    Ref.Ratio.method = paste0("obs", RefRat_Mode),
                                    ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                                    FDR.thresh = PTMs_FDR.thresholds[[Ptm]], arbitrary.lines = arbitrary.thr,
                                    proteins = prot.list, Proteins.col = "Code",
                                    return = TRUE, return.plot = TRUE,
                                    title = paste0(Ptm, " volcano plot "), subfolder = subDr,
                                    subfolderpertype = FALSE, Symmetrical = TwoSided,
                                    Size = "Rel. av. log10 abundance", Size.max = 2,
                                    plotly = create_plotly, plotly_local = create_plotly_local,
                                    plotly_labels = c(PepLabKol, paste0(Ptm, "-site")),
                                    cl = parClust)
          #k2 <- grep(topattern(paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)])), colnames(ptmpep), value = TRUE)
          #df2 <- ptmpep[, k2]
          #pepPlotFun(df1, df2, "Before VS after volc. plot", dir[1], FALSE)
          #
          #
          VP_list <- tempVPptm
          insrt <- ""
          Src <- paste0(libPath, "/extdata/R scripts/Sources/thresholds_Excel.R")
          #rstudioapi::documentOpen(Src)
          source(Src)
          #
          #
          ptmpep <- tempVPptm$Protein_groups_file
          volcano.plots[[Ptm]] <- tempVPptm$Plots
          n2 <- names(volcano.plots[[Ptm]]$Labelled)
          for (ttl in n2) {
            plot <- volcano.plots[[Ptm]]$Labelled[[ttl]]
            ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
          }
          if ((create_plotly)&&(!create_plotly_local)) {
            plot_ly[[paste0(Ptm, "_Volcano plots (t-tests)")]] <- tempVPptm$"Plotly plots"
          }
          # Specificity mark for untested proteins
          # We can assume that proteins with PSMs only in the specific pull-down samples are actually specifically enriched!
          # -> Label them as such!
          # NB: This used to be for pull-downs only but I am extending it to the whole workflow for now.
          #     This could also be parameter-controlled.
          # NB: For the F-test, this is done within the function
          #if (IsPullDown) {
          for (a in RG$values) { #a <- RG$values[1]
            e <- Exp.map[which(Exp.map[[RG$column]] == a),]
            grp1 <- unique(e[which(!e$Reference %in% c(TRUE, "TRUE")), VPAL$column])
            kolTR <- paste0("Regulated - ", grp1)
            w <- which(kolTR %in% colnames(ptmpep))
            if (length(w)) {
              grp1 <- grp1[w]
              kolTR <- kolTR[w]
              grp0 <- unique(e[which(e$Reference %in% c(TRUE, "TRUE")), VPAL$column])
              e0 <- unique(e$Ref.Sample.Aggregate[which(e[[VPAL$column]] == grp0)])
              kole0 <- paste0("Evidence IDs - ", e0)
              tst0 <- rowSums(ptmpep[, kole0] == "") == length(kole0)
              for (i in 1:length(grp1)) { #i <- 1
                e1i <- unique(e$Ref.Sample.Aggregate[which(e[[VPAL$column]] == grp1[i])])
                kole1i <- paste0("Evidence IDs - ", e1i)
                tst1i <- rowSums(ptmpep[, kole1i] != "") == length(kole1i)
                w1i <- which(tst0 & tst1i)
                if (length(w1i)) {
                  evcount <- apply(ptmpep[w1i, kole1i, drop = FALSE], 1, function(x) { length(unlist(strsplit(x, ";"))) })
                  txtup <- paste0("Specific: at least ", evcount, " PSMs/sample")
                  ptmpep[w1i, kolTR[i]] <- txtup
                }
              }
            } else { warning("I would expect \"Regulated ...\" columns in the modified peptides table by this stage!") }
          }
          #}
          # Create t-test filters:
          ## These can then be used for further steps down the line, such as volcano plots, etc...
          g <- grep("^Regulated - ", colnames(ptmpep), value = TRUE)
          g1 <- gsub("^Regulated - ", "", g)
          up <- grep("^up|^Specific", unique(unlist(ptmpep[, g])), value = TRUE)
          down <- grep("^down|^Anti-specific", unique(unlist(ptmpep[, g])), value = TRUE) # The "anti-specific" part will only become relevant if in future I disconnect symmetry and/or adding specific tags from pull-down experiments
          PTMs_Reg_filters[[Ptm]] <- list()
          PTMs_Reg_filters[[Ptm]]$"t-tests" <- list()
          if ("con" %in% filter_types) {
            PTMs_Reg_filters[[Ptm]]$"t-tests"$"By condition" <- list()
            for (i in 1:length(g)) {
              PTMs_Reg_filters[[Ptm]]$"t-tests"$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                                                Filter_up = sort(which(ptmpep[[g[i]]] %in% up)),
                                                                                Filter_down = sort(which(ptmpep[[g[i]]] %in% down)),
                                                                                Filter = sort(which(ptmpep[[g[i]]] %in% c(up, down))))
            }
          }
          if ("ref" %in% filter_types) {
            PTMs_Reg_filters[[Ptm]]$"t-tests"$"By reference" <- list()
            g2 <- as.data.frame(t(as.data.frame(strsplit(g1, "___"))))
            colnames(g2) <- VPAL$names
            tst <- apply(Exp.map[, RRG$names, drop = FALSE], 1, paste, collapse = "___")
            tmp <- apply(g2[,RRG$names, drop = FALSE], 1, paste, collapse = "___")
            g2$Ref <- sapply(tmp, function(x) { y <- unique(tst[which((Exp.map$Reference)&(tst == x))]) })
            for (i in unique(g2$Ref)) {
              w <- which(g2$Ref == i)
              PTMs_Reg_filters[[Ptm]]$"t-tests"$"By reference"[[i]] <- list(Columns = g[w],
                                                                            Filter_up = sort(which(apply(ptmpep[, g[w], drop = FALSE], 1, function(x) {
                                                                              length(which(x %in% up))
                                                                            }) > 0)),
                                                                            Filter_down = sort(which(apply(ptmpep[, g[w], drop = FALSE], 1, function(x) {
                                                                              length(which(x %in% down))
                                                                            }) > 0)),
                                                                            Filter = sort(which(apply(ptmpep[, g[w], drop = FALSE], 1, function(x) {
                                                                              length(which(x %in% c(up, down)))
                                                                            }) > 0)))
            }
          }
          if (sum(c("dat", "dat2") %in% filter_types)) {
            PTMs_Reg_filters[[Ptm]]$"t-tests"$"Whole dataset" <- list(Columns = g,
                                                                      Filter_up = sort(which(apply(ptmpep[, g, drop = FALSE], 1, function(x) {
                                                                        length(which(x %in% up))
                                                                      }) > 0)),
                                                                      Filter_down = sort(which(apply(ptmpep[, g, drop = FALSE], 1, function(x) {
                                                                        length(which(x %in% down))
                                                                      }) > 0)),
                                                                      Filter = sort(which(apply(ptmpep[, g, drop = FALSE], 1, function(x) {
                                                                        length(which(x %in% c(up, down)))
                                                                      }) > 0)))
          }
          if (F.test) {
            #kol <- grep(topattern(pepRf), colnames(ptmpep), value = TRUE)
            #View(ptmpep[, kol])
            ReportCalls <- AddMsg2Report(Msg = " - F-tests", Space = FALSE)
            # NB: id.col below should be unique!
            stopifnot(length(unique(ptmpep$Name)) == nrow(ptmpep))
            #
            dataType <- "modPeptides"
            #
            Src <- paste0(libPath, "/extdata/R scripts/Sources/run_F_test.R")
            #rstudioapi::documentOpen(Src)
            tstFtst <- try(source(Src), silent = TRUE)
            #
            if (!"try-error" %in% class(tstFtst)) {
              #F_test_ref_ratios %<o% F_volc$`Reference ratios` # Not needed
              volcano.plots[[Ptm]]$"F-tests_Unlabelled" <- F_volc$Plots$"Unlabelled"
              volcano.plots[[Ptm]]$"F-tests_Labelled" <- F_volc$Plots$"Labelled"
              n2 <- names(volcano.plots[[Ptm]]$"F-tests_Labelled")
              dir <- paste0(dir, "/Reg. analysis/F-tests")
              for (ttl in n2) {
                plot <- volcano.plots[[Ptm]]$"F-tests_Labelled"[[ttl]]
                ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
              }
              # Legacy code for web-hosted plotly plots:
              if ((create_plotly)&&(!create_plotly_local)) { plot_ly[[paste0(Ptm, "_Volcano plots (F-tests)")]] <- F_volc$"Plotly plots" }
              # Also a posteriori F-test P- or Q- values histogram:
              nbin <- 20
              bd <- (0:nbin)/nbin
              if (F_Root %in% colnames(ptmpep)) {
                temp <- data.frame(value = is.all.good(10^(-ptmpep[[F_Root]])))
                ttl <- paste0("Histogram: F-test moderated Pvalue (", ptm, ")")
                plot <- ggplot(temp, aes(x = value)) +
                  geom_histogram(bins = nbin, colour = "black", alpha = 0.25, fill = "green") +
                  guides(fill = "none") + theme_bw() + ggtitle(ttl)
                poplot(plot)
                dir <- paste0(wd, "/Workflow control/Protein groups/P-values")
                dirlist <- unique(c(dirlist, dir))
                if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
                ttla <- gsub(": *", " - ", ttl)
                ggsave(paste0(dir, "/", ttla, ".jpeg"), plot, dpi = 300)
                ggsave(paste0(dir, "/", ttla, ".pdf"), plot, dpi = 300)
                ReportCalls <- AddPlot2Report(Title = ttla)
              }
              # Create F-test filters:
              g <- grep("^mod\\. F-test Regulated - ", colnames(PTMs_F_test_data[[Ptm]]), value = TRUE)
              g1 <- gsub("^mod\\. F-test Regulated - ", "", g)
              PTMs_Reg_filters[[Ptm]]$"F-tests" <- list()
              if ("con" %in% filter_types) {
                PTMs_Reg_filters[[Ptm]]$"F-tests"$"By condition" <- list()
                for (i in 1:length(g)) {
                  PTMs_Reg_filters[[Ptm]]$"F-tests"$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                                                    Filter_up = sort(which(PTMs_F_test_data[[Ptm]][[g[i]]] %in% up)),
                                                                                    Filter_down = sort(which(PTMs_F_test_data[[Ptm]][[g[i]]] %in% down)),
                                                                                    Filter = sort(which(PTMs_F_test_data[[Ptm]][[g[i]]] %in% c(up, down))))
                }
              }
              if ("ref" %in% filter_types) {
                PTMs_Reg_filters[[Ptm]]$"F-tests"$"By reference" <- list()
                g2 <- as.data.frame(t(as.data.frame(strsplit(g1, " - "))))
                row.names(g2) <- NULL
                colnames(g2) <- c("Group", "Analysis", "Condition")
                g2$Ref <- ""
                w <- grep(" VS ", g2$Condition)
                if (length(w)) {
                  g2$tmp[w] <- sapply(strsplit(gsub("\\)$", "", g2$Condition[w]), split = " VS "), function(x) {x[[2]]})
                  g2$Ref[w] <- apply(g2[w, c("Group", "Analysis", "tmp")], 1, function(x) {
                    paste0(x[[1]], " - ", x[[2]], ", ref: ", x[[3]])
                  })
                }
                w <- which((grepl(" vs ", g2$Condition))&(!grepl(" VS ", g2$Condition)))
                if (length(w)) {
                  g2$tmp[w] <- sapply(strsplit(g2$Condition[w], split = " vs "), function(x) {x[[2]]})
                  g2$Ref[w] <- apply(g2[w, c("Group", "Analysis", "tmp")], 1, function(x) {
                    paste0(x[[1]], ", ref: ", x[[2]])
                  })
                }
                w <- which(!grepl(" vs ", g2$Condition))
                if (length(w)) { g2$Ref[w] <- g1[w] }
                for (i in unique(g2$Ref)) {
                  w <- which(g2$Ref == i)
                  PTMs_Reg_filters[[Ptm]]$"F-tests"$"By reference"[[i]] <- list(Columns = g[w],
                                                                                Filter_up = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                  length(which(x %in% up))
                                                                                }) > 0)),
                                                                                Filter_down = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                  length(which(x %in% down))
                                                                                }) > 0)),
                                                                                Filter = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                  length(which(x %in% c(up, down)))
                                                                                }) > 0)))
                }
              }
              if ("dat" %in% filter_types) {
                PTMs_Reg_filters[[Ptm]]$"F-tests"$"By analysis" <- list()
                g2 <- as.data.frame(t(as.data.frame(strsplit(g1, " - "))))
                colnames(g2) <- c("Group", "Analysis", "Condition")
                g2$Group_Analysis <- apply(g2[,c("Group", "Analysis")], 1, paste, collapse = "_")
                for (i in unique(g2$Group_Analysis)) {
                  w <- which(g2$Group_Analysis == i)
                  PTMs_Reg_filters[[Ptm]]$"F-tests"$"By analysis"[[i]] <- list(Columns = g[w],
                                                                               Filter_up = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                 length(which(x %in% up))
                                                                               }) > 0)),
                                                                               Filter_down = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                 length(which(x %in% down))
                                                                               }) > 0)),
                                                                               Filter = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                 length(which(x %in% c(up, down)))
                                                                               }) > 0)))
                }
              }
              if ("dat2" %in% filter_types) {
                PTMs_Reg_filters[[Ptm]]$"F-tests"$"Whole dataset" <- list(Columns = g,
                                                                          Filter_up = sort(which(apply(PTMs_F_test_data[[Ptm]][, g, drop = FALSE], 1, function(x) {
                                                                            length(which(x %in% up))
                                                                          }) > 0)),
                                                                          Filter_down = sort(which(apply(PTMs_F_test_data[[Ptm]][, g, drop = FALSE], 1, function(x) {
                                                                            length(which(x %in% down))
                                                                          }) > 0)),
                                                                          Filter = sort(which(apply(PTMs_F_test_data[[Ptm]][, g, drop = FALSE], 1, function(x) {
                                                                            length(which(x %in% c(up, down)))
                                                                          }) > 0)))
              }
            } else { warning("The F-test did not generate any plots! Investigate!") }
          }
          # Heatmap
          g <- paste0(pepRf, RSA$values)
          grps <- Exp.map[match(RSA$values, Exp.map$Ref.Sample.Aggregate), VPAL$column]
          w <- which(g %in% colnames(ptmpep))
          g <- g[w]; grps <- grps[w]
          temp <- ptmpep[, g]
          rownames(temp) <- gsub("\n", " ", ptmpep$Name)
          kol <- gsub(topattern(pepRf), "", g)
          colnames(temp) <- gsub(topattern(pepRf), "", colnames(temp))
          kol <- cleanNms(kol)
          colnames(temp) <- cleanNms(colnames(temp))
          w <- which(!is.finite(as.matrix(temp)), arr.ind = TRUE)
          temp[w] <- NA
          rwMns <- rowMeans(temp[, kol], na.rm = TRUE)
          w <- which(!is.na(rwMns))
          temp <- temp[w,]
          rwMns <- rwMns[w]
          temp[, kol] <- sweep(temp[, kol], 1, rwMns, "-")
          temp2 <- as.matrix(Data_Impute2(temp, grps)$Imputed_data)
          hcluster <- hclust(dist(t(temp2)))
          hdendro <- as.dendrogram(hcluster)
          hord <- order.dendrogram(hdendro)
          hord <- colnames(temp)[hord]
          vcluster <- hclust(dist(temp2))
          vdendro <- as.dendrogram(vcluster)
          vord <- order.dendrogram(vdendro)
          vord <- rownames(temp)[vord]
          require(ggdendro)
          hdendro.plot <- ggdendrogram(data = hdendro) + theme(axis.text.y = element_text(size = 0.1),
                                                               plot.margin = margin(0, 0, 0, 0, "cm"))
          vdendro.plot <- ggdendrogram(data = vdendro, rotate = TRUE, labels = FALSE, leaf_labels = FALSE) +
            theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                  panel.background = element_rect(fill = "transparent", colour = NA), 
                  plot.background = element_rect(fill = "transparent", colour = NA),
                  plot.margin = margin(0, 0, 0, 0, "cm"))
          # Data wrangling
          temp2 <- set_colnames(reshape::melt.data.frame(temp), c("Sample", "value"))
          temp2$Label <- rownames(temp)
          temp2$Sample <- as.character(temp2$Sample)
          temp2$Colour <- "grey"
          temp2$Size <- 1
          # Extract the order of the tips in the dendrograms
          # Order the levels according to their position in the clusters
          temp2$Xmin <- match(temp2$Sample, hord)-1
          temp2$Ymin <- match(temp2$Label, vord)-1
          Xscale <- length(unique(temp2$Sample))
          temp2$Label2 <- temp2$Label
          w <- which(nchar(temp2$Label2) > 25)
          temp2$Label2[w] <- paste0(substr(temp2$Label2[w], 1, 22), "...")
          # Create heatmap plot
          w1 <- which(temp2$Colour == "green")
          w2 <- which((temp2$Xmin == max(temp2$Xmin))&(temp2$Colour == "green"))
          nm <- paste0("Heatmap\n", Ptm, "-modified peptides")
          heatmap.plot <- ggplot(temp2) +
            geom_rect(aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
            geom_text(data = temp2, aes(y = Ymin+0.5, label = Label2),
                      x = length(unique(temp2$Sample)), colour = "grey", hjust = 0, vjust = 0.5, size = 2) +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                  axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                  panel.background = element_rect(fill = "transparent", colour = NA),
                  plot.margin = margin(0, 0, 0, 0, "cm")) +
            scale_fill_gradient2(low = "darkblue", mid = "lightgrey", high = "darkred") +
            xlab(NULL) + ylab(NULL)
          leg <- get_legend(heatmap.plot)
          heatmap.plot <- heatmap.plot + theme(legend.position = "none")
          htmp <- arrangeGrob(grobs = list(hdendro.plot, leg, heatmap.plot, vdendro.plot),
                              widths = c(1, 2*Xscale, 1, 3), heights = c(3, 5),
                              layout_matrix = rbind(c(NA, 1, NA, 2), c(3, 3, 3, 4)),
                              padding = 5)
          nm <- gsub("\n", " - ", nm)
          ggsave(paste0(dir[1], "/", nm, ".jpeg"), htmp, width = 20, height = 20, units = "in", dpi = 600)
          ggsave(paste0(dir[1], "/", nm, ".pdf"), htmp, width = 20, height = 20, units = "in", dpi = 600)
          ReportCalls <- AddPlot2Report(Title = nm, Dir = dir[1])
          #system(paste0("open \"", dir[1], "/", nm, ".jpeg", "\""))
          #system(paste0("open \"", dir[1], "/", nm, ".pdf", "\""))
          # Gene Ontology terms enrichment analysis
          if (enrichGO) {
            if ((!exists("GO_mappings"))&&(file.exists("GO_mappings.RData"))) {
              loadFun("GO_mappings.RData")
            }
            if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) {
              loadFun("GO_terms.RData")
            }
            p <- strsplit(ptmpep$Proteins, ";")
            test <- sapply(annot.col, function(x) { x %in% colnames(db) })
            annot.col2 <- annot.col[which(test)]
            tmpDB <- db[which(db$Observed), c("Protein ID", annot.col2)]
            clusterExport(parClust, c("tmpDB", "annot.col2"), envir = environment())
            temp <- parLapply(parClust, p, function(x) {
              m <- match(x, tmpDB$"Protein ID")
              y <- tmpDB[m, annot.col2]
              if (length(m) > 1) {
                y <- apply(y, 2, function(z) {
                  z <- z[which(!is.na(z))]
                  paste(sort(unique(unlist(strsplit(as.character(z), ";")))), collapse = ";")
                })
              }
              y <- setNames(as.character(unlist(y)), annot.col2)
              return(y)
            })
            temp <- as.data.frame((do.call(rbind, temp)))
            for (i in annot.col2) {
              w <- which(is.na(temp[[i]]))
              if (length(w)) { temp[w, i] <- "" }
            }
            ptmpep[, annot.col2] <- temp
            Kol <- c("Modified sequence", "Evidence IDs", "Sequence", "Proteins", "PEP", "id",
                     "Protein group IDs", "Razor protein group ID", "Proteins", "Leading razor proteins",
                     "Leading razor protein", "Gene names", "Protein names", "Gene names (all)",
                     "Protein names (all)", "Weights", "Match(es)", "Code", "Name", annot.col2,
                     "Potential contaminant")
            Kol <- Kol[which(Kol %in% colnames(ptmpep))]
            PTMs_GO_enrich.dat[[Ptm]] <- list()
            PTMs_GO_enrich.FCRt[[Ptm]] <- list()
            PTMs_GO_enrich.tbl[[Ptm]] <- list()
            PTMs_GO_Plots[[Ptm]] <- list()
            PTMs_Reg_GO_terms[[Ptm]] <- list()
            Tsts <- c("t-tests", "F-tests")
            WhTsts <- which(Tsts %in% names(PTMs_Reg_filters[[Ptm]]))
            for (tt in WhTsts) { #tt <- 1 #tt <- 2
              tstrt <- Tsts[tt]
              setwd(wd)
              dir <- paste0(wd, "/Reg. analysis/", ptm, "/GO enrich/", tstrt)
              if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
              dirlist <- unique(c(dirlist, dir))
              filt <- PTMs_Reg_filters[[Ptm]][[tstrt]]
              #By <- c("By condition", "By reference", "By analysis", "Whole dataset")
              By <- "By condition"
              By <- By[which(By %in% names(filt))]
              if (length(By)) {
                for (bee in By) { #bee <- By[1]
                  flt <- filt[[bee]]
                  if (bee == "Whole dataset") { flt <- list("Whole dataset" = flt) }
                  tstbee <- paste0(tstrt, "_", tolower(bee))
                  if (length(flt)) {
                    if (tt == 1) { temp <- ptmpep }
                    if (tt == 2) { temp <- PTMs_F_test_data[[Ptm]] }
                    for (kol in Kol) { temp[, kol] <- ptmpep[, kol] }
                    flt <- flt[order(names(flt))]
                    reg <- setNames(lapply(flt, function(x) { list(x$Columns) }), names(flt))
                    reg <- set_colnames(reshape2::melt(reg), c("Name", "Bleh", "For"))
                    reg$Bleh <- NULL
                    PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]] <- paste0("Mean ", c(ptms.ratios.ref[length(ptms.ratios.ref)], "log2(Ratio) - ")[tt])
                    reg$ParentFC <- gsub(".*Regulated - ", PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]], reg$Name)
                    reg$FCname <- paste0(PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]], reg$For)
                    UF <- unique(reg$For)
                    temPTM <- as.data.frame(sapply(UF, function(x) { #x <- UF[1]
                      x <- reg$ParentFC[which(reg$For == x)]
                      if (length(x) > 1) {
                        x <- apply(temp[, x], 1, log_ratio_av)
                      } else { x <- temp[[x]] }
                      return(x)
                    }))
                    colnames(temPTM) <- paste0(PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]], UF)
                    for (kol in Kol) { temPTM[, kol] <- ptmpep[, kol] }
                    #temPTM$"First protein" <- sapply(strsplit(temPTM$"Proteins", ";"), function(x) { unlist(x)[1] })
                    temp <- listMelt(strsplit(temPTM$Proteins, ";"), temPTM$id)
                    temp$Genes <- db$Gene[match(temp$value, db$`Protein ID`)]
                    temp <- aggregate(temp$Genes, list(temp$L1), function(x) {
                      paste(sort(unique(unlist(strsplit(x, ";")))), collapse = ";")
                    })
                    temPTM$Genes <- temp$x[match(temPTM$id, temp$Group.1)]
                    PTMs_GO_enrich.dat[[Ptm]][[tstbee]] <- temPTM
                    PTMs_GO_enrich.tbl[[Ptm]][[tstbee]] <- reg
                    flt <- setNames(lapply(UF, function(x) { flt[[x]]$Filter }), UF)
                    ttr <- btr <- ""
                    if (length(By) > 1) { ttr <- btr <- paste0(tolower(bee), "_") }
                    Pep.Ref.Filt <- tmpFilt <- setNames(lapply(names(flt), function(x) { 1:nrow(PTMs_GO_enrich.dat[[Ptm]][[tstbee]]) }), names(flt))
                    if (tt == 1) {
                      if (GO.enrich.MultiRefs) {
                        Pep.Ref.Filt <- try(setNames(lapply(names(flt), function(x) { #x <- names(flt)[1]
                          m <- Exp.map[which(Exp.map[[VPAL$column]] == x), , drop = FALSE]
                          y <- unique(m[[GO.enrichment.Ref.Aggr$column]])
                          stopifnot(length(y) == 1)
                          w <- which(Exp.map[[GO.enrichment.Ref.Aggr$column]] == y)
                          w1 <- which(apply(ptmpep[, paste0(ptms.ref, Exp.map$Ref.Sample.Aggregate[w])], 1, function(x) {
                            sum(is.all.good(x, 2))
                          }) > 0)
                          w2 <- flt[[x]] # Required for if we are imputing missing values
                          return(sort(unique(c(w1, w2))))
                        }), names(flt)), silent = TRUE)
                        if ("try-error" %in% class(Pep.Ref.Filt)) {
                          warning("Invalid \"GO.enrichment.Ref.Aggr\" argument: multiple references for GO enrichment are only feasible if each enrichment filter maps to a single reference! Skipping...")
                          GO.enrich.MultiRefs <- FALSE
                          Pep.Ref.Filt <- tmpFilt
                        }
                      }
                    }
                    if (tt == 2) {
                      cM <- as.data.frame(contrMatr_F)
                      Pep.Ref.Filt <- setNames(lapply(names(flt), function(x) { #x <- names(flt)[1]
                        x <- EM$Ref.Sample.Aggregate[which(EM$Group_ %in% unlist(expContrasts_F$All[match(x, expContrasts_F$name)]))]
                        kol <- paste0(ptms.ref[1], x)
                        which(parApply(parClust, ptmpep[, kol], 1, function(x) { length(proteoCraft::is.all.good(x)) }) > 0)
                      }), names(flt))
                    }
                    if ((length(Pep.Ref.Filt) > 1)||(!is.na(Pep.Ref.Filt))) {
                      flt <- setNames(lapply(names(flt), function(x) { flt[[x]][which(flt[[x]] %in% Pep.Ref.Filt[[x]])] }), names(flt))
                    }
                    ReportCalls <- AddMsg2Report(Msg = " - GO terms enrichment analysis", Space = FALSE)
                    #
                    Mode <- "regulated"
                    dataType <- "modPeptides"
                    #
                    Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_enrich.R")
                    #rstudioapi::documentOpen(Src)
                    source(Src)
                    #
                    clueGO_outDir <- dir
                    clueGO_type <- "Enrichment (Right-sided hypergeometric test)"
                    Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_enrich.R")
                    #rstudioapi::documentOpen(Src)
                    source(Src)
                    #
                    # Cleanup - do it now, not within sources!
                    for (i in allArgs) { try(rm(i), silent = TRUE) }
                    #
                    PTMs_GO_Plots[[Ptm]][[tstbee]] <- goRES
                    #
                    if (!is.null(names(PTMs_GO_Plots[[Ptm]][[tstbee]]))) {
                      if ("All_GO_terms" %in% names(PTMs_GO_Plots[[Ptm]][[tstbee]])) {
                        GO_terms <- PTMs_GO_Plots[[Ptm]][[tstbee]]$All_GO_terms
                        PTMs_GO_Plots[[Ptm]][[tstbee]]$All_GO_terms <- NULL  
                      }
                      n2 <- names(PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_plots)
                      for (ttl in n2) { #ttl <- n2[1]
                        plot <- PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_plots[[ttl]]
                        ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
                      }
                      if ((create_plotly)&&(!create_plotly_local)) {
                        plot_ly[[paste0(Ptm, "_GO plots - Regulated vs Observed - ", tstbee)]] <- PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_plot_ly
                      }
                      temp <- PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms
                      temp$Mapping <- NULL
                      if ("Offspring" %in% colnames(temp)) {
                        temp$Offspring <- parSapply(parClust, temp$Offspring, paste, collapse = ";")
                      }
                      temp$`Protein table row(s)` <- NULL
                      colnames(temp) <- cleanNms(colnames(temp), start = FALSE)
                      gn <- grep("^Genes", colnames(temp), value = TRUE)
                      pr <- grep("^Proteins", colnames(temp), value = TRUE)
                      pp <- grep("^Peptide IDs", colnames(temp), value = TRUE)
                      kn <- grep("^Count", colnames(temp), value = TRUE)
                      pv <- grep("^Pvalue", colnames(temp), value = TRUE)
                      zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
                      lf <- grep("^logFC", colnames(temp), value = TRUE)
                      si <- grep("^Significance", colnames(temp), value = TRUE)
                      #lp <- grep("^Leading protein IDs", colnames(temp), value = TRUE)
                      #kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pp, lp))]
                      kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pp, pr))]
                      #temp <- temp[, c(kl, si, gn, pp, lp, kn, pv, zs, lf)]
                      temp <- temp[, c(kl, si, gn, pp, pr, kn, pv, zs, lf)]
                      w <- apply(temp[, pv, drop = FALSE], 1, function(x) { sum(!is.na(x)) }) > 0
                      temp <- temp[w,]
                      tst <- apply(temp[, kn, drop = FALSE], 1, function(x) { sum(x[which(!is.na(x))]) })
                      temp <- temp[order(tst, decreasing = TRUE),]
                      temp <- temp[order(temp$Ontology, decreasing = FALSE),]
                      PTMs_Reg_GO_terms[[Ptm]][[tstbee]] <- temp
                      write.csv(temp, file = paste0(dir, "/", Ptm, " GO terms - ", tstbee, ".csv"), row.names = FALSE)
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
                      HdrStl <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                                            wrapText = TRUE, numFmt = "TEXT", fontSize = 12)
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
                          setColWidths(wb, ont, which(colnames(temp) %in% pp), 20)
                          addStyle(wb, ont, HdrStl, 1, 1:ncol(temp))
                          addStyle(wb, ont, createStyle(numFmt = "0"), 2:(length(w)+1),
                                   which(colnames(temp) %in% kn), gridExpand = TRUE)
                          addStyle(wb, ont, createStyle(numFmt = "0.000"), 2:(length(w)+1),
                                   which(colnames(temp) %in% c(zs, lf, pv)), gridExpand = TRUE)
                        }
                      }
                      if (kount) {
                        saveWorkbook(wb, paste0(dir, "/", Ptm, " GO terms - ", tstbee, ".xlsx"), overwrite = TRUE)
                        if (tt == 1) {
                          Kol2 <- paste0("Significance - ", cleanNms(VPAL$values), " ", max(BH.FDR)*100, "%")
                          Kol2 <- Kol2[which(Kol2 %in% colnames(PTMs_Reg_GO_terms[[Ptm]][[tstbee]]))]
                        }
                        if (tt == 2) {
                          Kol2 <- paste0("Significance - ", names(PTMs_Reg_filters[[Ptm]][[tstrt]][[bee]]), " ", max(BH.FDR)*100, "%")
                          Kol2 <- Kol2[which(Kol2 %in% colnames(PTMs_Reg_GO_terms[[Ptm]][[tstbee]]))]
                        }
                        if (length(Kol2)) {
                          if (length(Kol2) > 1) {
                            w <- which(apply(PTMs_Reg_GO_terms[[Ptm]][[tstbee]][, Kol2], 1, function(x) {"+" %in% x}))
                          } else { w <- which(sapply(PTMs_Reg_GO_terms[[Ptm]][[tstbee]][,Kol2], function(x) {"+" %in% x})) }
                          write.csv(PTMs_Reg_GO_terms[[Ptm]][[tstbee]][w,],
                                    file = paste0(dir, "/", Ptm, " regulated GO terms - ", tstbee, ".csv"),
                                    row.names = FALSE)
                        }
                        # Summary table and heatmap of number of co-regulated GO terms
                        Kol3 <- Kol2
                        N <- length(Kol3)
                        if (N > 1) {
                          temp <- as.data.frame(matrix(rep("", (N+1)^2), ncol = N+1))
                          W <- lapply(Kol3, function(x) { which(PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms[[x]] == "+") })
                          names(W) <- gsub(paste0("^Significance - | ", max(BH.FDR)*100, "%$"), "", Kol3)
                          temp[2:(N+1), 1] <- temp[1, 2:(N+1)] <- names(W)
                          for (i in 2:(N+1)) { #i <- 2
                            temp[i, 2:(N+1)] <- sapply(Kol3, function(x) {
                              sum((PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms[[x]] == "+")&(PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms[[Kol3[i]]] == "+"),
                                  na.rm = TRUE)
                            })
                          }
                          names(W) <- cleanNms(gsub(" [0-9]+%$", "", names(W)))
                          tst <- sapply(strsplit(names(W), " - "), length)
                          tst <- (min(tst) > 1)&(length(unique(tst)) == 1)
                          if (tst) {
                            tst <- as.data.frame(t(sapply(strsplit(names(W), " - "), unlist)))
                            l <- apply(tst, 2, function(x) { length(unique(x)) })
                            tst <- tst[, which(l > 1), drop = FALSE]
                            names(W) <- apply(tst, 1, paste, collapse = " - ")
                          }
                          temp[2:(N+1), 1] <- temp[1, 2:(N+1)] <- names(W)
                          nm <- paste0("N. of co-regulated GO terms\n", tstrt, "\n(", tolower(bee), ")")
                          write.csv(temp, file = paste0(dir, "/", gsub("\n", " - ", gsub("\n\\(", " (", nm)), ".csv"), row.names = FALSE)
                          temp2 <- temp[2:(N+1), 2:(N+1)]
                          colnames(temp2) <- temp[1, 2:(N+1)]
                          rownames(temp2) <-  temp[2:(N+1), 1]
                          for (i in 1:nrow(temp2)) { temp2[[i]] <- as.numeric(temp2[[i]]) }
                          if (max(is.all.good(unlist(temp2)))) {
                            temp2 <- as.matrix(temp2)
                            basic.heatmap(temp2, "N. of co-regulated GO terms", paste0(tstrt, "\n(", tolower(bee), ")"),
                                          save = c("pdf", "jpeg"), folder = dir)
                          } else { warning(paste0(tstrt, " co-regulated GO terms heatmap: nothing to plot")) }
                        } else { warning(paste0(tstrt, " performed for only one condition, skipping.")) }
                      }
                    }
                  } else { warning(paste0("Filter ", tstbee, " has length 0, skipping.")) }
                }
              } else { warning(paste0("No filters available for ", tstrt, ", skipping.")) }
            }
          }
          PTMs_pep[[Ptm]] <- ptmpep
          ReportCalls <- AddSpace2Report()
          PTMs_int.ref[[Ptm]] <- ptms.ref
          PTMs_rat.ref[[Ptm]] <- ptms.ratios.ref
          #write.csv(ptmpep, paste0("Tables/", Ptm, "-modified peptides.csv"), row.names = FALSE)
          #w <- grsep2(prot.list, ptmpep$Proteins)
          #if (length(w)) {
          #  write.csv(ptmpep[w,], paste0("Tables/", Ptm, "-modified peptide - proteins in list.csv"), row.names = FALSE)
          #}
        } else {
          warning(paste0("There are no ", ptm, "-modified peptides to report! Check your parameters!"))
        }
      } else {
        warning(paste0("PTM \"", ptm, "\" could not be recognized, please use names consistent with those specified in the search."))
      }
    }
  } else { stop("Really? There is no PTM-modified class of peptides to analyze? Why did you run this workflow then?") }
} else { stop("Really? There is no PTM-modified class of peptides to analyze? Why did you run this workflow then?") }

# Once we are done with using GO stuff, backup the final versions
if (Annotate) {
  # (Testing for existence in case we are rerunning part of the script)
  if (exists("GO_terms")) { saveFun(GO_terms, file = "GO_terms.RData") }
  if (exists("GO_mappings")) { saveFun(GO_mappings, file = "GO_mappings.RData") }
  # We can also remove our GO_mappings, we shan't be using them anymore
  rm(GO_terms)
  rm(GO_mappings)
  # otherwise they are also unused further down.
  .obj <- .obj[which(!.obj %in% c("GO_terms", "GO_mappings"))]
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
gc()
# It makes sense to close/re-create parallel clusters regularly to reduce memory usage
stopCluster(parClust)
source(parSrc)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - Create output tables
## PSMs
dir <- paste0(wd, "/Tables")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
w <- which(sapply(colnames(ev), function(x) { "list" %in% class(ev[[x]]) }))
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
    ColumnsTbl <- ColumnsTbl[which(sapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }) > 0)]
    ColumnsTbl <- set_colnames(reshape::melt.list(ColumnsTbl), c("Col", "Grp"))
    #tst <- aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), length); View(tst)
    #tst <- aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), unique); View(tst)
    #tst <- aggregate(1:nrow(ColumnsTbl), list(ColumnsTbl$Col), unique); w <- which(sapply(tst$x, length) > 1); setNames(tst$x[w], tst$Group.1[w])
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
    Src <- paste0(libPath, "/extdata/R scripts/Sources/parWrite_Excel_core_script.R")
    #rstudioapi::documentOpen(Src)
    source(Src)
  }
}
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/Write_Excel_end_script.R")
#rstudioapi::documentOpen(Src)
source(Src)
#WorkBook$get_active_sheet()
#xl_open(repFl)

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)



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
Mod2Venn <- names(PTMs_pep)
II <- setNames(1, "All peptidoforms")
#if ((exists("PTMs_pep"))&&(length(PTMs_pep))) {
  Mod2Write <- names(PTMs_pep)
  II[paste0(Mod2Write, "-mod. pept.")] <- 1+(seq_along(length(Mod2Write)))
#}
for (ii in II[2:length(II)]) { #ii <- II[2] #ii <- II[3]
  if (ii == 1) {
    # Placeholder
    # dir <- paste0(wd, "/Venn diagrams")
    # ttest_Filt <- Reg_filters$"t-tests"$"By condition"
    # vennRoot <- ""
    # myData <- pep
    # myRef <- grep("^Mean [^ ]+ log10\\(", colnames(myData), value = TRUE)
    # myRef <- paste0(unique(gsub(" - .*", "", myRef)), " - ")
    # myRef <- myRef[length(myRef)]
    # idKol <- "Modified sequence"
    # infoKol <- "Proteins"
    # if (F.test) {
    #   Ftest_Filt <- Reg_filters$"F-tests"$"By condition"
    #   myFData <- F_test_data
    # }
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
  #View(myData[, paste0(myRef, VPAL$values)])
  comp_list <- setNames(lapply(VPAL$values, function(grp) { #grp <- VPAL$values[1]
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
for (i in 1:length(MatMetCalls$Texts$DatAnalysis)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls, paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$DatAnalysis[", i,"], prop = WrdFrmt$",
                                                        c("Body", "Template_text")[(MatMetCalls$Texts$DatAnalysis[i] == "TEMPLATE")+1], "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#

# Finalize analysis
Src <- paste0(libPath, "/extdata/R scripts/Sources/Finalize_analysis.R")
#rstudioapi::documentOpen(Src)
source(Src)

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
