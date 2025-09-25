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
#ReUseAnsw %<o% FALSE
scrptType %<o% "noReps"
scrptTypeFull %<o% "noReps_PG_and_PTMs"

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

# Parameters used by the master script:
###-|-### Workflows: setNames(c("Discovery -> no comparisons, complex sample", "Band ID -> no comparisons, focus on coverage and the top proteins", "Regulation -> ratio analysis (up and down)", "Pull-down (incl. BioID) -> ratio analysis (choice between up only or up and down)"), c("Discovery", "Band ID", "Regulation", "Pull-down"))
###-|-### Replicates? FALSE
###-|-### External dependencies: Excel (loose); ScanHeadsman (loose)

### Packages
## For convenience all (or most) of the packages used are installed here:
## CRAN packages:
if(!exists("cran_req")) { cran_req %<o% "pak" } else { cran_req %<o% cran_req }
if(!exists("bioc_req")) { bioc_req %<o% c() } else { bioc_req %<o% bioc_req }
cran_req <- unique(c(cran_req, "pak", "shiny", "renv", "R.utils",
                     "qs2", "DT", "ggplot2", "ggpubr", "reshape2", "compiler", "stats",
                     "rgl", "ggrepel", "rstudioapi", "gtools", "minpack.lm", "parallel", "openxlsx", "openxlsx2", "openssl", "plotly",
                     "Peptides", "venn", "ggdendro", "ggpubr", "colorspace", "ggnewscale", "viridis", "factoextra", "NbClust", "gridExtra",
                     "svDialogs", "htmlwidgets", "magrittr", "tibble", "fs", "officer", "snow", "imputeLCMD", "ggplotify", "cowplot", "plyr",
                     "shinyjs", "shinyBS", "shinyFiles", "shinyWidgets", "TeachingDemos", "shinycssloaders", "jpeg", "data.table", "stringi",
                     "readr", "ssh", "taxize", "arrow"))
bioc_req <- unique(c(bioc_req, "UniProt.ws", "pcaMethods", "impute", "GO.db", "topGO", "pcaMethods"))
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
        if ("try-error" %in% class(tst)) { try(download.file(url, destfile, "wget"), silent = TRUE) }
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

# Load backup?
load_a_Bckp %<o% c(TRUE, FALSE)[match(svDialogs::dlg_message("Do you want to load a backup?", "yesno")$res, c("yes", "no"))]
if (load_a_Bckp) {
  tmp <- openxlsx2::read_xlsx(paste0(homePath, "/Default_locations.xlsx"))
  load_Bckp(startDir = tmp$Path[which(tmp$Folder == "Temporary folder")])
}

# Set Shiny options, load functions for creating a Word report, create Excel styles
Src <- paste0(libPath, "/extdata/R scripts/Sources/ShinyOpt_Styles_and_Report.R")
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
Src <- paste0(libPath, "/extdata/R scripts/Sources/Start_analysis.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
AnalysisParam %<o% list("Input folder" = inDirs,
                        "Temp. folder" = wd,
                        "N. of threads" = N.clust)

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

#### Code chunk - MS files map
if ("PTM-enriched" %in% colnames(FracMap)) {
  FracMap$"PTM-enriched"[which(!FracMap$"PTM-enriched" %in% Modifs$`Full name`)] <- NA
} else { FracMap$"PTM-enriched" <- NA }
# Try sorting automatically
tst <- grepl("_[0-9]+\\.d$", FracMap$`Raw file`)
if (sum(tst)) {
  FracMap$Bruker_run_ID <- NA
  FracMap$Bruker_run_ID[which(tst)] <- as.integer(gsub(".*_|\\.d$", "", FracMap$`Raw file`[which(tst)]))
  w1 <- which(!is.na(FracMap$Bruker_run_ID))
  w2 <- which(is.na(FracMap$Bruker_run_ID))
  w1 <- w1[order(FracMap$Bruker_run_ID[w1])]
  FracMap <- FracMap[c(w1, w2),]
}
FracMap$Use <- as.logical(FracMap$Use)
FracMap$Use[which(is.na(FracMap$Use))] <- TRUE
nr <- nrow(FracMap)
rws <- seq_len(nr)
# Original table column widths
wTest0 <- setNames(vapply(colnames(FracMap), function(k) { #k <- colnames(FracMap)[1]
  tmp <- FracMap[[k]]
  if ("logical" %in% class(tmp)) { tmp <- as.integer(tmp) }
  tmp <- as.character(tmp)
  x <- max(nchar(c(k, tmp)) + 3, na.rm = TRUE)
  x <- x*10
  if (is.na(x)) { x <- 15 } else { x <- max(c(ceiling(x/10)*10, 30)) }
  return(x)
}, 1), colnames(FracMap))
# Dummy table
frMap <- FracMap
frMap$"Raw files name" <- NULL
frMap$Use <- as.logical(toupper(frMap$Use))
frMap$Use[which(is.na(FracMap$Use))] <- TRUE
frMap$Use <- shinyCheckInput(frMap$Use, "Use")
frMap$"PTM-enriched" <- shinySelectInput(FracMap$"PTM-enriched",
                                         "PTMenriched",
                                         unique(c(Modifs$`Full name`, NA)),
                                         paste0(30*max(c(nchar(as.character(Modifs$`Full name`)), 2)), "px"))
frMap$"Parent sample" <- shinyTextInput(frMap$"Parent sample", "Sample", paste0(wTest0["Parent sample"], "px"))
k <- c("Raw file", "Parent sample", "Fraction", "PTM-enriched", "Bruker_run_ID", "Use")
k <- k[which(k %in% colnames(frMap))]
frMap <- frMap[, k]
# Estimate dummy table column widths
wTest1 <- vapply(colnames(frMap), function(k) { #k <- colnames(frMap)[1]
  if ((k == "Parent sample")&&(!k %in% names(wTest0))) { k <- "MQ.Exp" }
  if (k %in% names(wTest0)) { x <- wTest0[k] } else { x <- 30 }
  return(x)
}, 1)
wTest2 <- sum(wTest1) + 15 + ncol(frMap)*5
wTest1 <- paste0(as.character(wTest1), "px")
wTest1 <- aggregate((1:length(wTest1))-1, list(wTest1), c)
wTest1 <- apply(wTest1, 1, function(x) {
  x2 <- as.integer(x[[2]])
  list(width = x[[1]],
       targets = x2,
       names = colnames(frMap)[x2+1])
})
#
appNm <- paste0(dtstNm, " - ", FracMapNm)
if (exists("frMap2")) { rm(frMap2) }
ui <- shinyUI(fluidPage(titlePanel(tag("u", FracMapNm),
                                   appNm),
                        useShinyjs(),
                        setBackgroundColor( # Doesn't work
                          color = c(#"#F8F8FF",
                            "#EFE6F5"),
                          gradient = "linear",
                          direction = "bottom"
                        ),
                        extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
                        mainPanel(
                          withSpinner(DT::DTOutput("myFracMap", width = wTest2)),
                          br(),
                          actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                        )))
server <- function(input, output, session) {
  frMap2 <- FracMap
  output$myFracMap <- DT::renderDT({ frMap },
                                   FALSE,
                                   escape = FALSE,
                                   class = "compact",
                                   selection = "none",
                                   rownames = FALSE,
                                   editable = FALSE,
                                   options = list(dom = "t",
                                                  paging = FALSE,
                                                  ordering = FALSE,
                                                  autowidth = TRUE,
                                                  columnDefs = wTest1,
                                                  scrollX = FALSE),
                                   # the callback is essential to capture the inputs in each row
                                   callback = JS("table.rows().every(function(i, tab, row) {
  var $this = $(this.node());
  $this.attr('id', this.data()[0]);
  $this.addClass('shiny-input-container');
});
Shiny.unbindAll(table.table().node());
Shiny.bindAll(table.table().node());"))
  observeEvent(input$myFracMap_cell_edit, { # This is not goddamn working!!!
    kl <- colnames(frMap)[input$myFracMap_cell_edit$col+1]
    if (kl %in% colnames(frMap2)) {
      frMap2[input$myFracMap_cell_edit$row, kl] <<- input$myFracMap_cell_edit$value
    } else {
      warning(paste0("Could not find column ", kl, " from dummy table in the real table!"))
    }
  })
  observeEvent(input$saveBtn, {
    for (k in c("Use", "PTM-enriched", "Parent sample")) {
      root <- gsub("-", "", k)
      if (k == "Use") { tp <- TRUE }
      if (k %in% c("PTM-enriched", "Parent sample")) {
        tp <- ""
        if (k == "Parent sample") { root<- "Sample" }
      }
      frMap2[[k]] <- vapply(rws, function(x) { input[[paste0(root, "___", as.character(x))]] }, tp)
    }
    assign("frMap2", frMap2, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("frMap2"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
FracMap %<o% frMap2
#
tst <- try(write.csv(FracMap, file = FracMapPath, row.names = FALSE), silent = TRUE)
while ("try-error" %in% class(tst)) {
  dlg_message(paste0("File \"", FracMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(FracMap, file = FracMapPath, row.names = FALSE), silent = TRUE)
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

AnalysisParam$Type <- WorkFlow
MakeRatios %<o% FALSE
RatiosThresh_2sided %<o% TRUE
# No need to reload the local copy, values are updated in environment
if (WorkFlow %in% c("Regulation", "Pull-down")) {
  if (length(unique(FracMap$`Parent sample`)) == 1) {
    warning("Only one sample, skipping ratios analysis!")
    WorkFlow <- "Discovery"
  } else { MakeRatios <- TRUE }
}

#### Code chunk - Samples map
SamplesMapNm %<o% "Samples map"
SamplesMapPath %<o% paste0(wd, "/", SamplesMapNm, ".csv")
if ((!file.exists(SamplesMapPath))||(!nrow(SamplesMap))) {
  SamplesMap %<o% data.frame("MQ.Exp" = exp,
                             "Negative Filter" = FALSE,
                             check.names = FALSE)
}
if (MakeRatios) {
  if (!"Ratios group" %in% colnames(SamplesMap)) { SamplesMap$"Ratios group" <- 1 }
  if (!"Reference" %in% colnames(SamplesMap)) { SamplesMap$Reference <- FALSE }
}
nr <- nrow(SamplesMap)
rws <- 1:nr
# if ("Order" %in% colnames(SamplesMap)) {
#   u <- unique(SamplesMap$Order)
#   u <- u[which(!u %in% rws)]
#   if (length(u) < nr) { SamplesMap$Order <- rws }
# } else { SamplesMap$Order <- rws }
smplMapKol1 <- c("Reference", "Negative Filter", "Use")
for (kol in smplMapKol1) {
  if (kol %in% colnames(SamplesMap)) {
    SamplesMap[[kol]] <- as.logical(toupper(SamplesMap[[kol]]))
    SamplesMap[which(is.na(SamplesMap[[kol]])), kol] <- c(FALSE, TRUE)[(kol == "Use")+1]
  }
}
if (!"Use" %in% colnames(SamplesMap)) { SamplesMap$Use <- TRUE }
SamplesMap$Use <- suppressWarnings(as.logical(SamplesMap$Use))
SamplesMap$Use[which(is.na(SamplesMap$Use))] <- TRUE
smplMapKol <- smplMapKol1
if (MakeRatios) {
  smplMapKol <- c("Ratios group", smplMapKol1)
}
allIDs <- as.character(sapply(smplMapKol, function(x) {
  paste0(x, "___", rws)
}))
# Original table column widths
wTest0 <- setNames(vapply(colnames(SamplesMap), function(k) { #k <- colnames(SamplesMap)[1]
  l <- k
  if (k == "MQ.Exp") { l <- "Parent sample"}
  tmp <- SamplesMap[[k]]
  if ("logical" %in% class(tmp)) { tmp <- as.integer(tmp) }
  tmp <- as.character(tmp)
  x <- max(nchar(c(k, tmp)) + 3, na.rm = TRUE)
  x <- x*10
  if (is.na(x)) { x <- 30 } else { x <- max(c(ceiling(x/10)*10, 30)) }
  return(x)
}, 1), colnames(SamplesMap))
#
smplMap2 <- smplMap <- SamplesMap[, which(!colnames(SamplesMap) %in% "New name")] # This column is deprecated and ignored
colnames(smplMap2)[which(colnames(smplMap2) == "MQ.Exp")] <- "Parent sample"
# Estimate table column widths
wTest1 <- vapply(colnames(smplMap2), function(k) { #k <- colnames(smplMap2)[1]
  if (k == "Parent sample") { k <- "MQ.Exp" }
  if (k %in% names(wTest0)) { x <- wTest0[k] } else { x <- 30 }
  return(x)
}, 1)
wTest2 <- max(c(sum(wTest1) + 15 + ncol(smplMap2)*5, 600))
wTest1 <- paste0(as.character(wTest1), "px")
wTest1 <- aggregate((1:length(wTest1))-1, list(wTest1), c)
wTest1 <- apply(wTest1, 1, function(x) {
  x2 <- as.integer(x[[2]])
  list(width = x[[1]],
       targets = x2,
       names = colnames(smplMap2)[x2+1])
})
#
IDs <- c()
MSG <- reactiveVal("")
for (kol in c("Reference", "Negative Filter", "Use")) {
  if (kol %in% colnames(SamplesMap)) {
    smplMap2[[kol]] <- shinyCheckInput(SamplesMap[[kol]], kol)
    IDs <- c(IDs, paste0(kol, "___", seq_len(nrow(smplMap2))))
    stopifnot(length(IDs) == length(unique(IDs)))
  }
}
if (MakeRatios) {
  kol <- "Ratios group"
  smplMap2[[kol]] <- shinyNumInput(SamplesMap[[kol]], 1, Inf, 1, 1, root = kol)
  IDs <- c(IDs, paste0(kol, "___", seq_len(nrow(smplMap2))))
  stopifnot(length(IDs) == length(unique(IDs)))
}
#
if (exists("expOrder")) {
  tst <- sum(!exp %in% expOrder)
  if (sum(tst)) { rm(expOrder) }
}
if (!exists("expOrder")) {
  expOrder <- exp
}
#
msg <- ""
appNm <- paste0(dtstNm, " - Experiment map")
ui <- fluidPage(useShinyjs(),
                setBackgroundColor( # Doesn't work
                  color = c(#"#F8F8FF",
                    "#EBEFF7"),
                  gradient = "linear",
                  direction = "bottom"
                ),
                titlePanel(tag("u", "Experiment map"),
                           appNm),
                extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
                tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
                mainPanel(h4("Optional: define samples order"),
                          selectInput("expOrder", "",
                                      expOrder, expOrder, TRUE, TRUE, width = "1200px"),
                          br(),
                          br(),
                          span(uiOutput("Message"), style = "color:red"),
                          withSpinner(DT::DTOutput("mySampleMap", width = wTest2))),
                br(),
                actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"))
if (exists("smplMap3")) { rm(smplMap3) }
server <- function(input, output, session) {
  # Create copies table
  smplMap3 <- smplMap
  output$Message <- renderUI({ em(" ") })
  output$mySampleMap <- DT::renderDT({ smplMap2 },
                                     FALSE,
                                     escape = FALSE,
                                     class = "compact",
                                     selection = "none",
                                     rownames = FALSE,
                                     editable = FALSE,
                                     options = list(dom = "t",
                                                    paging = FALSE,
                                                    ordering = FALSE,
                                                    autowidth = TRUE,
                                                    columnDefs = wTest1,
                                                    scrollX = TRUE),
                                      callback = JS("table.rows().every(function(i, tab, row) {
  var $this = $(this.node());
  $this.attr('id', this.data()[0]);
  $this.addClass('shiny-input-container');
});
Shiny.unbindAll(table.table().node());
Shiny.bindAll(table.table().node());"))
  #
 # Observe events
  # if (length(IDs)) {
  #   sapply(IDs, function(id) {
  #     x <- gsub(".+_", "", id)
  #     kol <- substr(id, 1, nchar(id)-nchar(x)-1)
  #     x <- as.integer(x)
  #     observeEvent(input[[id]], {
  #       cat(input[[id]], "\n")
  #       val <- input[[id]]
  #       if (kol %in% c("Reference", "Negative Filter", "Use")) {
  #         val <<- as.logical(val)
  #       }
  #       SamplesMap[x, kol] <<- val
  #       if (kol == "Use") {
  #         tst <- sum(SamplesMap[kol])
  #         if (!tst) {
  #           msg <- "Include at least one sample!"
  #           shinyjs::disable("saveBtn")
  #         }
  #         if (tst) {
  #           msg <- " "
  #           shinyjs::enable("saveBtn")
  #         }
  #       }
  #       if ((MakeRatios)&&(kol == "Reference")) {
  #         l <- length(unique(SamplesMap[[kol]]))
  #         if (l != 2) {
  #           msg <- "This workflow compares different samples, select one sample as reference!"
  #           shinyjs::disable("saveBtn")
  #         }
  #         if (l == 2) {
  #           msg <- " "
  #           shinyjs::enable("saveBtn")
  #         }
  #       }
  #       output$Message <- renderUI({ em(msg) })
  #     })
  #   })
  # }
  observeEvent(input$mySampleMap_cell_edit, {
    smplMap3[input$mySampleMap_cell_edit$row,
             input$mySampleMap_cell_edit$col+1] <- input$mySampleMap_cell_edit$value
  })
  observeEvent(input$expOrder, {
    tmp <- input$expOrder
    tmp <- c(tmp, exp[which(!exp %in% tmp)])
    assign("expOrder", tmp, envir = .GlobalEnv)
  })
  # sapply(allIDs, function(id) {
  #   tmp <- unlist(strsplit(id, "___"))
  #   kol <- tmp[1]
  #   i <- tmp[2]
  #   observeEvent(input[[id]], {
  #     #cat("Event", id, "\n")
  #     #cat(input[[id]], "\n")
  #     smplMap3[i, kol] <<- as.logical(input[[id]])
  #   })
  # })
  observeEvent(input$saveBtn, {
    for (k in smplMapKol) {
      smplMap3[[k]] <- sapply(rws, function(x) { input[[paste0(k, "___", x)]] })
    }
    assign("smplMap3", smplMap3, envir = .GlobalEnv)
    tmp <- input$expOrder
    tmp <- c(tmp, exp[which(!exp %in% tmp)])
    assign("expOrder", tmp, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("smplMap3"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
SamplesMap %<o% smplMap3
exp <- expOrder
SamplesMap <- SamplesMap[match(SamplesMap$MQ.Exp, exp), ]
#
tst <- try(write.csv(SamplesMap, file = SamplesMapPath, row.names = FALSE), silent = TRUE)
while ("try-error" %in% class(tst)) {
  dlg_message(paste0("File \"", SamplesMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(SamplesMap, file = SamplesMapPath, row.names = FALSE), silent = TRUE)
}

if ("Negative Filter" %in% colnames(SamplesMap)) {
  SamplesMap$"Negative Filter" <- as.logical(toupper(SamplesMap$"Negative Filter"))
}
AnalysisParam$"Ratios analysis" <- MakeRatios
# Rename "MQ.Exp" to "Experiment"
colnames(SamplesMap)[which(colnames(SamplesMap) == "MQ.Exp")] <- "Experiment"

# Labelling
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
Src <- paste0(libPath, "/extdata/R scripts/Sources/Process_Fasta_DBs.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

evNm %<o% "PSM"
#evNm %<o% c("PSM", "Evidence")[(SearchSoft == "MAXQUANT")+1]

#### Code chunk - Load and process annotations
## This includes a QC step in case the database differs slightly from the one used by MQ, or if somehow some IDs have not been properly parsed.
Src <- paste0(libPath, "/extdata/R scripts/Sources/Load_Annotations.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
source(parSrc, local = FALSE)
Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_prepare.R") # Doing this earlier but also keep latter instance for now
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
AnalysisParam$Annotations <- Annotate

#### Code chunk - Define analysis parameters
Src <- paste0(libPath, "/extdata/R scripts/Sources/noRep_Parameters_editor_Main.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

ev$"Protein group IDs" <- NULL
ev$"Peptide ID" <- NULL

Script <- readLines(ScriptPath)
gc()
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

# Create Materials and Methods template
Src <- paste0(libPath, "/extdata/R scripts/Sources/autoMatMet.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Update peptide-to-protein mappings
Src <- paste0(libPath, "/extdata/R scripts/Sources/checkPep2Prot.R")
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
if (length(w) == 2) {
  for (i in genkol) { temp[[i]] <- strsplit(temp[[i]], ";") }
  temp <- apply(temp, 1, function(x) { paste(sort(unique(unlist(x))), collapse = ";") })
}
temp2 <- as.data.table(listMelt(strsplit(ev$Proteins, ";"), 1:nrow(ev)))
temp2$Gene <- temp[match(temp2$value, db$`Protein ID`)]
temp2 <- temp2[, list(Genes = paste(Gene, collapse = ";")), keyby = list(id = L1)]
temp2 <- as.data.frame(temp2)
ev$"Gene names" <- temp2$Genes[match(ev$id, temp2$id)]

# Deal with PTM-enriched data
FracMap$`PTM-enriched`[which(FracMap$`PTM-enriched` == "NA")] <- NA
EnrichedPTMs %<o% unique(FracMap$`PTM-enriched`[which(!is.na(FracMap$`PTM-enriched`))])
PTMriched %<o% (length(EnrichedPTMs) > 0)
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
    tmp[[ptm]] <- c("-", "+")[tmp[[ptm]] + 1]
    tmp[[ptm]] <- factor(tmp[[ptm]], levels = c("-", "+"))
    tmp$"Raw file" <- gsub(".*/|\\.[^\\.]+$", "", tmp$"Raw file")
    plot <- ggplot(tmp) + geom_bar(aes(x = `Raw file`, y = Count, fill = .data[[ptm]]), stat = "identity") +
      scale_fill_viridis(discrete = TRUE, option = "H", begin = 0.25, end = 0.8) + ggtitle(ttl) + 
      theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
    poplot(plot)
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150)
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150)
    #
    PTMev2Remov[[ptm]] <- lapply(Exp, function(exp) { #exp <- Exp[1]
      res <- c()
      fm <- FracMap[which(FracMap$`Parent sample` == exp),]
      tst <- unique(fm$`PTM-enriched`[which(fm$`Parent sample` == exp)])
      if ((ptm %in% tst)&&(length(tst) > 1)) {
        # We need to filter out peptides with the mark from all non-enriched samples
        smpls1 <- fm$`Raw file`[which((is.na(fm$`PTM-enriched`))|(fm$`PTM-enriched` != ptm))] # Non-enriched/flow through samples
        smpls2 <- fm$`Raw file`[which(fm$`PTM-enriched` == ptm)] # Enriched samples
        if (!smpls1 %in% unique(ev$`Raw file path`)) { # Pre-empt bugs
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
  msg <- paste0("Removing ", l, " PSMs (", signif(100*l/nrow(ev), 2), "%) with enriched PTMs from non-enriched samples.")
  ev <- ev[which(!ev$id %in% PTMev2Remov),]
}

# Negative filter
NegFilt %<o% (("Negative Filter" %in% colnames(SamplesMap))&&(sum(SamplesMap$`Negative Filter`)))

# Filter PSMs
# Optional (not used): remove "wrongly assigned" charge 1 evidences
RemovZ1 %<o% FALSE
w1 <- which(ev$Charge == 1)
wHt1 <- which(ev$Charge > 1)
if ((RemovZ1)&&(length(w1))) {
  cat("Removing the following presumably bogus identifications with Z=1:\n",
             paste(unique(ev$`Modified sequence`[w1]), collapse = "\n"))
  ev <- ev[wHt1,]
}
#
# Remove evidences with null intensity values
w1 <- which((is.all.good(ev$Intensity, 2))&(ev$Intensity > 0))
w2 <- which(!1:nrow(ev) %in% w1)
l2 <- length(w2)
if (l2) {
  warning(paste0("Removing ", l2, " (", round(100*l2/nrow(ev), 2), "%) peptide evidences with invalid or null intensity values!"))
  nullEv <- ev[w2,]
  ev <- ev[w1,]
}
w1 <- which((is.na(ev$Reverse))|(ev$Reverse != "+"))
w2 <- which(!1:nrow(ev) %in% w1)
l2 <- length(w2)
if (l2) {
  warning(paste0("Removing ", l2, " (", round(100*l2/nrow(ev), 2), "%) reverse peptide evidences!"))
  revEv <- ev[w2,]
  ev <- ev[w1, ]
}
w1 <- which((is.na(ev$"Potential contaminant"))|(ev$"Potential contaminant" != "+"))
w2 <- which(!1:nrow(ev) %in% w1)
l2 <- length(w2)
if (l2) {
  message(paste0(round(100*l2/nrow(ev), 2), "% of valid identifications are potential contaminants."))
  #  contEv <- ev[-w, ]
  #  ev <- ev[w, ]
}

Exp <- Exp[which(Exp %in% ev$Experiment)] # Update experiments

# DIA-only: MS2-based correction of MS1-based quantitative values
Src <- paste0(libPath, "/extdata/R scripts/Sources/MS2corr2MS1.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

#### Code chunk - Pepper correction
# (Maybe this should be done at peptides level?)
if (runPepper) {
  if (!require("qs", quietly = TRUE)) { install.packages("qs") }
  library(qs)
  PepScrptsDir %<o% paste0(RPath, "/proteoCraft/extdata/Pepper")
  stopifnot(dir.exists(PepScrptsDir))
  pyInit <- paste0(PepScrptsDir, "/Python_start.R")
  source(pyInit)
  #
  # Placeholder for when this gets ported to the main script
  # if (scrptType == "withReps") {
  #   if ("Pepper" %in% names(ev.col)) { ref <- ev.col[match("Pepper", names(ev.col))-1] } else {
  #     ref <- ev.col[length(ev.col)]
  #   }
  #   nuRef <- ev.col["Pepper"] <- paste0("Pepper ", ev.col["Original"])
  # }
  if (scrptType == "noReps") {
    if ("Pepper" %in% names(int.cols)) { ref <- int.cols[match("Pepper", names(int.cols))-1] } else {
      ref <- int.cols[length(int.cols)]
    }
    nuRef <- int.col <- int.cols["Pepper"] <- "Pepper-adj. intensity"
  }
  #
  # Detect recommended training parameters for all datasets
  trainParam <- suppressWarnings(readLines(paste0(PepScrptsDir, "/train_models.sh")))
  args <- c("peptide_file", "n_runs", "seq_length", "output_file", "filter_size", "n_filters", "n_layers", "n_nodes",
            "dropout", "learning_rate", "batch", "random_run")
  g <- grep("^python3 ", trainParam)
  l <- length(trainParam)
  trainParam <- sapply(args, function(arg) {
    sapply(1:length(g), function(x) {
      h <- grep(paste0("--", arg), trainParam)
      h <- gsub(paste0(".*'--", arg, "' +"), "", trainParam[h[which(h %in% g[x]:c(g, l)[x+1])]])
      pat1 <- paste0("^", substr(h, 1, 1))
      pat2 <- paste0(substr(h, 1, 1), ".*$")
      h <- gsub(pat2, "", gsub(pat1, "", h))
      return(h)
    })
  })
  trainParam <- as.data.frame(trainParam)
  rownames(trainParam) <- gsub(".+/|\\.tsv$", "", trainParam$peptide_file)
  #
  # Hyper-parameter tuning
  # a) Prepare dataset for one-hot-encoding
  smplKol <- "Experiment"
  cat("Samples (= Experiments):\n -", paste0(unique(SamplesMap[[smplKol]]), collapse = "\n - "), "\n")
  w <- which(ev[1,] %in% SamplesMap[[smplKol]])
  stopifnot(length(w) > 0)
  pepDir <- paste0(wd, "/Pepper")
  if (!dir.exists(pepDir)) { dir.create(pepDir, recursive = TRUE) }
  tempEv <- Pepper_TrainingData(ev, Modifs, SamplesMap, smplKol, smplKol, path = paste0(pepDir, "/ModPep4Pepper.tsv"), intCol = ref)
  Modifs <- tempEv$PTMs
  tempEv <- tempEv$Data
  #
  # Edit and run python script for One-Hot-Encoding
  #################################################
  Grps <- colnames(tempEv)[(max(grep("^Charge [0-9]+$", colnames(tempEv)))+1):length(colnames(tempEv))]
  lGrps <- length(Grps)
  cat(paste0("Dataset = ", dtstNm, ", number of samples = ", lGrps, "\n"))
  #
  OHfl <- paste0(pepDir, "/OneHotEncodedPepQuants.tsv")
  OHE <- readLines(paste0(PepScrptsDir, "/onehot_encode_peptide_sequences.py"))
  arg1 <- which(OHE == "input_filename = sys.argv[1]")
  arg2 <- which(OHE == "runs = int(sys.argv[2])")
  arg3 <- which(OHE == "output_filename = sys.argv[3]")
  OHE[arg1] <- paste0("input_filename = \"", paste0(pepDir, "/ModPep4Pepper.tsv"), "\"")
  OHE[arg2] <- paste0("runs = int(", length(Grps), ")")
  OHE[arg3] <- paste0("output_filename = \"", OHfl, "\"")
  write(OHE, paste0(pepDir, "/OneHotEncode.py"))
  #
  cmd <- paste0("python \"", pepDir, "/OneHotEncode.py\"")
  #cat(cmd)
  system(cmd)
  stopifnot(file.exists(OHfl))
  #
  # Hyperparameter tuning:
  tuneDir <- paste0(pepDir, "/Hyperparameter tuning")
  if (!dir.exists(tuneDir)) { dir.create(tuneDir, recursive = TRUE) }
  CleanUp <- FALSE
  GrdSrchFls <- list.files(tuneDir, "Gridsearch_Results_", full.names = TRUE)
  if (length(GrdSrchFls)) {
    opt <- setNames(c(FALSE, TRUE), c("No                                                                                                                                              ",
                                      "Yes                                                                                                                                             "))
    CleanUp <- opt[dlg_list(names(opt), names(opt)[1],
                            title = paste0(dtstNm, ": shall we start from scratch hyperparameter tuning?"))$res]
    names(CleanUp) <- NULL
  }
  # NB: This code is written in a way that it will always sample Step1Tries + Step2Tries parameter combinations.
  # 1 - randomly
  Step <- 1
  paramKol <- c("n_conv_layers", "n_filters", "filter_size", "n_layers", "n_nodes", "dropout", "batch", "learning_rate")
  if (CleanUp) { for (fl in GrdSrchFls) { unlink(fl) } } # Cleanup for fresh tuning
  HypTunFl <- paste0(PepScrptsDir, "/gridseach_parameters_neural_network.py")
  stopifnot(file.exists(HypTunFl))
  HypTunSrc <- paste0(PepScrptsDir, "/HyperTun.R")
  source(HypTunSrc)
  L <- length(GrdSrchFls)
  N <- Step1Tries - L
  kount <- 0
  while ((N > 0)&&(kount < N+1)) {
    clusterExport(parClust, "cmd", envir = environment())
    # Below, we are never running less tests than the number of cores, so we do not waste a round on a few tests.
    tst <- parSapply(parClust, 1:max(c(N.clust, N)), function(x) {
      try(system(cmd), silent = TRUE)
    })
    GrdSrchFls <- list.files(tuneDir, "Gridsearch_Results_", full.names = TRUE)
    L <- length(GrdSrchFls)
    N <- Step1Tries - L
    kount <- kount + 1
  }
  Step1Tries <- L # Update it if we got any bonus tests
  GrdSrch <- as.data.frame(t(sapply(strsplit(gsub("\\.tsv$", "", GrdSrchFls), "_"), function(x) {
    x <- unlist(x)
    x <- suppressWarnings(as.numeric(x))
    x <- x[which(!is.na(x))]
    return(x)
  })))
  colnames(GrdSrch) <- paramKol
  tst <- apply(GrdSrch, 2, function(x) { length(unique(x)) })
  GrdSrch <- GrdSrch[, which(tst > 1), drop = FALSE]
  paramKol2 <- paramKol[which(tst > 1)]
  rs <- plyr::rbind.fill(lapply(GrdSrchFls, read.delim))
  rs$X <- NULL
  GrdSrch[, gsub("\\.", "_", colnames(rs))] <- rs
  GrdSrch <- GrdSrch[order(GrdSrch$Val_percent_improvement, decreasing = TRUE),]
  bstImprov <- max(GrdSrch$Val_percent_improvement)
  # 2 - explore immediate neighborhood of best 10 results
  Step <- 2
  currParam <- GrdSrch[1:10, c(paramKol2, "No_of_epochs")]
  source(HypTunSrc)
  L <- length(GrdSrchFls)
  N <- Step1Tries + Step2Tries - L
  kount <- 0
  while ((N > 0)&&(kount < N+1)) {
    clusterExport(parClust, "cmd", envir = environment())
    # Below, we are never running less tests than the number of cores, so we do not waste a round on a few tests.
    tst <- parSapply(parClust, 1:max(c(N.clust, N)), function(x) {
      try(system(cmd), silent = TRUE)
    })
    GrdSrchFls <- list.files(tuneDir, "Gridsearch_Results_", full.names = TRUE)
    L <- length(GrdSrchFls)
    N <- Step1Tries + Step2Tries - L
    kount <- kount + 1
  }
  GrdSrch <- as.data.frame(t(sapply(strsplit(gsub("\\.tsv$", "", GrdSrchFls), "_"), function(x) {
    x <- unlist(x)
    x <- suppressWarnings(as.numeric(x))
    x <- x[which(!is.na(x))]
    return(x)
  })))
  colnames(GrdSrch) <- paramKol
  tst <- apply(GrdSrch, 2, function(x) { length(unique(x)) })
  GrdSrch <- GrdSrch[, which(tst > 1), drop = FALSE]
  paramKol2 <- paramKol[which(tst > 1)]
  rs <- plyr::rbind.fill(lapply(GrdSrchFls, read.delim))
  rs$X <- NULL
  GrdSrch[, gsub("\\.", "_", colnames(rs))] <- rs
  GrdSrch <- GrdSrch[order(GrdSrch$Val_percent_improvement, decreasing = TRUE),]
  bstImprov2 <- max(GrdSrch$Val_percent_improvement)
  print("")
  print(paste0(" - Step 1: best improvement = ", round(bstImprov, 2), "%"))
  print("")
  print(paste0(" - Step 2: best improvement = ", round(bstImprov2, 2), "%"))
  print("")
  #View(GrdSrch)
  #
  # Get final, best parameters
  bstParam <- GrdSrch[1, c(paramKol2, "No_of_epochs")]
  print(bstParam)
  bstParam$peptide_file <- OHfl
  bstParam$n_runs <- as.character(length(Grps))
  bstParam$seq_length <- "60"
  bstParam$random_run <- "0" # Obviously: we will do a proper training now!
  bstParam$output_file <- dtstNm
  #
  # Create model, using our own dataset and the above hypertuned parameters
  TrCoeff <- readLines(paste0(PepScrptsDir, "/peptide_coefficient_predictor.py"))
  
  # Parameters from train_models.sh"
  w <- which(TrCoeff %in% c("parser = argparse.ArgumentParser()", "args = parser.parse_args()", "print(args)"))
  TrCoeff[w] <- ""
  basePat <- paste0("^ *", topattern("parser.add_argument(", start = FALSE), " *")
  for (arg in colnames(trainParam)) { #arg <- colnames(trainParam)[1]
    pat <- paste0(basePat, "\'--", arg)
    g <- grep(pat, TrCoeff)
    stopifnot(length(g) == 1)
    val <- bstParam[[arg]]
    tst1 <- suppressWarnings(!is.na(as.numeric(val)))
    tst2 <- suppressWarnings(!is.na(as.integer(val)))
    if (tst1) {
      if (tst2) { TrCoeff[g] <- paste0(arg, " = float(", val, ")") } else {
        TrCoeff[g] <- paste0(arg, " = int(", val, ")")
      }
    } else { TrCoeff[g] <- paste0(arg, " = '", val, "'") }
  }
  # Default parameters
  G <- grep(basePat, TrCoeff)
  if (length(G)) {
    for (g in G) { #g <- G[1]
      tmp <- TrCoeff[g]
      arg <- gsub("\'.*", "", gsub(paste0(basePat, "\'--"), "", tmp))
      stopifnot(nchar(arg) > 0)
      val <- unlist(strsplit(tmp, "[,\\(] *default *= *"))
      if (!length(val) == 2) { stop(paste0("No default for argument ", arg, "! Give me a hand here, human!")) }
      val <- gsub(" *[,\\)].*", "", val[2])
      if (tst1) {
        if (tst2) { TrCoeff[g] <- paste0(arg, " = float(", val, ")") } else {
          TrCoeff[g] <- paste0(arg, " = int(", val, ")")
        }
      } else { TrCoeff[g] <- paste0(arg, " = '", val, "'") }
    }
  }
  TrCoeff <- gsub(" args\\.", " ", TrCoeff)
  TrCoeff <- gsub("\\(args\\.", "(", TrCoeff)
  arg1a <- which(TrCoeff == "data_df = pd.read_csv(peptide_file, sep = '\\t', index_col = 0)")
  TrCoeff[arg1a] <- "data_df = pd.read_csv(peptide_file, sep = '\t', index_col = 0)"
  #TrCoeff[arg4a] <- paste0("predicted_coefficients.to_csv(\"", output_file + \".tsv\", sep = '\\t')")
  g <- grep("trained_models/", TrCoeff) #TrCoeff[g]
  TrCoeff[g] <- gsub("trained_models/", paste0(pepDir, "/"), TrCoeff[g])
  # Small modification to adjust for a deprecation warning
  w <- which(TrCoeff == "init = tf.compat.v1.initialize_all_variables()")
  TrCoeff[w] <- "init = tf.compat.v1.global_variables_initializer()"
  # Bug fix
  w <- which(TrCoeff == "test_runs = np.array([[2*i, 2*i+1] for i in test_runs]).astype(int).ravel()")
  TrCoeff <- c(TrCoeff[1:w],
               "#################################################################################################",
               "### Insertion to remove cases where the above lines generate run indices outside of the range ###",
               "#################################################################################################",
               "train_runs = train_runs[train_runs+1 <= q_df.shape[1]]",
               "test_runs = test_runs[test_runs+1 <= q_df.shape[1]]",
               "#################################################################################################",
               "### Insertion ends                                                                            ###",
               "#################################################################################################",
               TrCoeff[(w+1):length(TrCoeff)])
  # Add code to save final layer dimensions
  w <- which(TrCoeff == "print(\"Saved model to disk\")")
  TrCoeff <- c(TrCoeff[1:w],
               "#################################################################################################",
               "### Insertion: this code saves the dimensions of the last weights layer to a small local file ###",
               "#################################################################################################",
               "weights = model.get_weights()",
               "dim = np.shape(weights[len(weights)-1])",
               "dim = np.matrix(dim)",
               #dim = np.reshape(dim, [2,1])",
               paste0("np.savetxt(\"",  pepDir, "/final_layer_dims.tsv\", dim, delimiter = '\t', fmt='%i')"),
               "#################################################################################################",
               "### Insertion ends                                                                            ###",
               "#################################################################################################",
               TrCoeff[(w+1):length(TrCoeff)])
  write(TrCoeff, paste0(pepDir, "/TrainCoeff.py"))
  
  py_clear_last_error()
  cat(" - Running coefficients training script...\n")
  cmd <- paste0("python \"", pepDir, "/TrainCoeff.py\"")
  #cat(cmd)
  tst <- try(system(cmd), silent = TRUE)
  stopifnot(tst == 0)
  #
  modlFl <- paste0(pepDir, "/", dtstNm, "_Coefficient_Predictor_Model.h5")
  if (tst == 0) { Outcome <- file.exists(modlFl) }
  if (Outcome) { cat("Success!!!\n") } else { cat("Failure!!!!!\n") }
  py_clear_last_error()
  #
  # Transfer to whole dataset
  # a) Prepare dataset for one-hot-encoding
  smplKol <- "Experiment"
  w <- which(ev[1,] %in% SamplesMap[[smplKol]])
  stopifnot(length(w) > 0)
  tempEv2 <- Pepper_ProcessData(ev, Modifs, SamplesMap, smplKol, smplKol, path = paste0(pepDir, "/ModPep4Pepper.tsv"), intCol = ref,
                                filter = FALSE)
  tempEv2 <- tempEv2$Data
  tmp <- data.frame(Peptides = nrow(tempEv2), Proteins = length(unique(tempEv2$Protein)), Samples = length(SamplesMap$Experiment))
  write.table(tmp, paste0(pepDir, "/Dataset dimensions.tsv"), sep = "\t", row.names = FALSE)
  Grps2 <- colnames(tempEv2)[(max(grep("^Charge [0-9]+$", colnames(tempEv2)))+1):length(colnames(tempEv2))]
  lGrps2 <- length(Grps2)
  #
  # Edit and run python script for One-Hot-Encoding
  #################################################
  #
  OHfl <- paste0(pepDir, "/OneHotEncodedPepQuants_Full.tsv")
  OHE <- readLines(paste0(PepScrptsDir, "/onehot_encode_peptide_sequences.py"))
  arg1 <- which(OHE == "input_filename = sys.argv[1]")
  arg2 <- which(OHE == "runs = int(sys.argv[2])")
  arg3 <- which(OHE == "output_filename = sys.argv[3]")
  OHE[arg1] <- paste0("input_filename = \"", paste0(pepDir, "/ModPep4Pepper.tsv"), "\"")
  OHE[arg2] <- paste0("runs = int(", length(Grps2), ")")
  OHE[arg3] <- paste0("output_filename = \"", OHfl, "\"")
  write(OHE, paste0(pepDir, "/OneHotEncode.py"))
  #
  cmd <- paste0("python \"", pepDir, "/OneHotEncode.py\"")
  #cat(cmd)
  system(cmd)
  #
  # Edit and run coefficients predictor transfer script
  #####################################################
  TrsfrCoeff <- readLines(paste0(PepScrptsDir, "/peptide_coefficient_predictor_transfer.py"))
  # Parameters from train_models.sh"
  w <- which(TrsfrCoeff %in% c("parser = argparse.ArgumentParser()", "args = parser.parse_args()", "print(args)"))
  TrsfrCoeff[w] <- ""
  basePat <- paste0("^ *", topattern("parser.add_argument(", start = FALSE), " *")
  bstParam2 <- bstParam
  tmp <- read.delim(paste0(pepDir, "/Dataset dimensions.tsv"))
  bstParam2$n_runs <- tmp$Samples
  bstParam2$output_file <- pepDir
  for (arg in colnames(trainParam)) { #arg <- colnames(trainParam)[1]
    pat <- paste0(basePat, "\'--", arg)
    g <- grep(pat, TrsfrCoeff)
    stopifnot(length(g) == 1)
    if (arg == "peptide_file") { val <- OHfl } else { val <- bstParam2[[arg]] }
    tst1 <- suppressWarnings(!is.na(as.numeric(val)))
    tst2 <- suppressWarnings(!is.na(as.integer(val)))
    if (tst1) {
      if (tst2) { TrsfrCoeff[g] <- paste0(arg, " = float(", val, ")") } else {
        TrsfrCoeff[g] <- paste0(arg, " = int(", val, ")")
      }
    } else { TrsfrCoeff[g] <- paste0(arg, " = '", val, "'") }
  }
  # Default parameters
  G <- grep(basePat, TrsfrCoeff)
  if (length(G)) {
    for (g in G) { #g <- G[1]
      tmp <- TrsfrCoeff[g]
      arg <- gsub("\'.*", "", gsub(paste0(basePat, "\'--"), "", tmp))
      stopifnot(nchar(arg) > 0)
      val <- unlist(strsplit(tmp, "[,\\(] *default *= *"))
      if (!length(val) == 2) { stop(paste0("No default for argument ", arg, "! Give me a hand here, human!")) }
      val <- gsub(" *[,\\)].*", "", val[2])
      if (tst1) {
        if (tst2) { TrsfrCoeff[g] <- paste0(arg, " = float(", val, ")") } else {
          TrsfrCoeff[g] <- paste0(arg, " = int(", val, ")")
        }
      } else { TrsfrCoeff[g] <- paste0(arg, " = '", val, "'") }
    }
  }
  TrsfrCoeff <- gsub(" args\\.", " ", TrsfrCoeff)
  TrsfrCoeff <- gsub("\\(args\\.", "(", TrsfrCoeff)
  arg1a <- which(TrsfrCoeff == "data_df = pd.read_csv(peptide_file, sep = '\\t', index_col = 0)")
  TrsfrCoeff[arg1a] <- "data_df = pd.read_csv(peptide_file, sep = '\t', index_col = 0)"
  TrsfrCoeff <- gsub("\"trained_models/\" +\\+ +output_file +\\+ +", "output_file + \"/\" + ", TrsfrCoeff)
  # Small modification to adjust for a deprecation warning
  w <- which(TrsfrCoeff == "init = tf.compat.v1.initialize_all_variables()")
  TrsfrCoeff[w] <- "init = tf.compat.v1.global_variables_initializer()"
  #
  # The most important:
  # The model path is hard-coded into the transfer script. Change the path to the one we want.
  w <- grep("^target_model\\.load_weights\\(", TrsfrCoeff)
  TrsfrCoeff[w] <- paste0("target_model.load_weights('", pepDir, "/", dtstNm, "_Coefficient_Predictor_Model.h5') ")
  #
  # Similarly, the dataset dimensions are hard coded too
  w1 <- which(TrsfrCoeff == "model = define_model(2438, 28)")
  tmp <- read.delim(paste0(pepDir, "/Dataset dimensions.tsv"))
  TrsfrCoeff[w1] <- paste0("model = define_model(", tmp$Proteins, ", len(train_runs))")
  w2 <- which(TrsfrCoeff == "target_model = define_model(5668, 612)")
  tmp <- unlist(read.table(paste0(pepDir, "/final_layer_dims.tsv"), col.names = FALSE))
  TrsfrCoeff[w2] <- paste0("target_model = define_model(", tmp[1], ", ", tmp[2], ")")
  w3 <- which(TrsfrCoeff == "            if run_count == 612:")
  TrsfrCoeff <- c(TrsfrCoeff[1:(w3-1)],
                  "            self.alphas = tf.Variable(np.random.rand(protein_count, run_count), trainable = True, dtype = 'float32')",
                  TrsfrCoeff[(w3+11):length(TrsfrCoeff)])
  #
  write(TrsfrCoeff, paste0(pepDir, "/TransferCoeff.py"))
  py_clear_last_error()
  cat(" - Running coefficients transfer script...\n")
  #saveImgFun(paste0(pepDir, "/Pepper_bckp.RData"))
  #loadFun(paste0(dtst, "/Pepper_bckp.RData"))
  #
  cmd <- paste0("python \"", pepDir, "/TransferCoeff.py\"")
  #cat(cmd)
  tst <- try(system(cmd), silent = TRUE)
  coeffFl <- paste0(pepDir, "/_transfer_inferred_coefficients_run0.tsv")
  if (tst == 0) { Outcome <- file.exists(modlFl) }
  if (Outcome) { cat("Success!!!\n") } else { cat("Failure!!!!!\n") }
  py_clear_last_error()
  #
  # Apply correction
  # Adj. int. = int. / coeff!!!
  tmp <- read.delim(coeffFl)
  colnames(tmp) <- c("ID", "Pepper coefficient")
  tmp$ID <- as.integer(gsub("^id_", "", tmp$ID))
  tmp$ModSeq <- ev$"Modified sequence"[match(tmp$ID, ev$id)]
  #
  ttl <- "Pepper - distribution of coefficients"
  plot <- ggplot(tmp) + geom_histogram(aes(x = log10(`Pepper coefficient`)), fill = "red", bins = 250) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() + ggtitle(ttl)
  poplot(plot)
  ggsave(paste0(pepDir, "/", ttl, ".jpeg"), plot, dpi = 150)
  #
  tmp2 <- ev[[ref]]/tmp$"Pepper coefficient"[match(ev$"Modified sequence", tmp$ModSeq)]
  tmp2 <- tmp2*median(ev[[ref]])/median(tmp2)
  ev[[nuRef]] <- tmp2
  #
  ttl <- "Pepper - corrected intensities distribution"
  tmp <- rbind(data.frame("Modified sequence"= ev$"Modified sequence",
                          "Intensity" = ev[[ref]],
                          "Type" = "Original",
                          check.names = FALSE),
               data.frame("Modified sequence"= ev$"Modified sequence",
                          "Intensity" = ev[[nuRef]],
                          "Type" = "Pepper-corrected",
                          check.names = FALSE))
  plot <- ggplot(tmp) + geom_density(stat = "density", aes(x = log10(`Intensity`), fill = Type), alpha = 0.3) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() + ggtitle(ttl)
  poplot(plot)
  ggsave(paste0(pepDir, "/", ttl, ".jpeg"), plot, dpi = 150)
  #
  # Test results ourselves
  g <- grep(";", ev$Proteins, invert = TRUE)
  tmp <- ev[g, c("Modified sequence", "Proteins", "Experiment", ref, nuRef)]
  colnames(tmp) <- c("ModSeq", "Proteins", "Experiment", "Original", "Pepper")
  tmp <- as.data.table(tmp)
  obj <- setNames(c("Original", "Pepper"), c("Before", "After"))
  call <- paste0(vapply(Exp, function(x) { paste0(x, " = sd(", x, ")") }, ""), collapse = ", ") # start building our call
  for (i in 1:2) {
    tmp2 <- copy(tmp)
    tmp2 <- dcast(tmp2, ModSeq + Proteins ~ Experiment, value.var = obj[i], fun.aggregate = sum, na.rm = TRUE)
    tmp2 <- as.data.frame(tmp2)
    w <- which(tmp2[, Exp] == 0, arr.ind = TRUE)
    tmp2[, Exp][w] <- NA
    tmp2 <- as.data.table(tmp2)
    call2 <- paste0("tmp2 <- tmp2[, list(", call, "), by = list(Proteins = Proteins)]")
    eval(parse(text = call2), envir = .GlobalEnv)
    tmp2 <- as.data.frame(tmp2)
    assign(names(obj)[i], tmp2)
    tmp2 <- colSums(tmp2[, Exp], na.rm = TRUE)
    assign(paste0(names(obj)[i], "_sums"), tmp2)
  }
  msg <- paste0("Pepper:\n\n#######\n   Final improvement per sample:\n\n\t -> Sum of intra-protein SDs, ratio before/after\n\t\t",
                paste(Exp, collapse = "\t"), "\n\t\t",
                paste(round(Before_sums/After_sums, 1), collapse = "\t"), "\n\n\t -> Overall: ",
                round(sum(is.all.good(unlist(Before[, Exp])))/sum(is.all.good(unlist(After[, Exp]))), 1), "\n\n")
  cat(msg)
  write(msg, paste0(pepDir, "/Final outcome.txt"))
}

Script <- readLines(ScriptPath)
gc()
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - Create modified peptides table
pep %<o% data.table(id = ev$id, seq = ev$"Modified sequence")
pep <- pep[, list(`Evidence IDs` = paste(id, collapse = ";")), by = list(`Modified sequence` = seq)]
pep <- as.data.frame(pep)
# Create matches between ev and pep
ev_to_pep <- match(pep$"Modified sequence", ev$"Modified sequence") # Only to be used when all matches would return the same value, e.g. sequence
pep_to_ev <- match(ev$"Modified sequence", pep$"Modified sequence") # In the reverse direction, works for all
# Proteins and sequence columns
if ("Modified sequence_verbose" %in% colnames(ev)) { pep$"Modified sequence_verbose" <- ev$"Modified sequence_verbose"[ev_to_pep] }
pep[, c("Sequence", "Proteins")] <- ev[ev_to_pep,  c("Sequence", "Proteins")]
# PEP: for a peptide, the lowest of the PEP of individual matching evidences
temp <- aggregate(ev$PEP, list(ev$"Modified sequence"), function(x) {
  x <- x[which(!is.na(x))]
  if (length(x)) { return(min(x)) } else { return(NA) }
})
pep$PEP <-  temp$x[match(pep$"Modified sequence", temp$Group.1)]
# IDs
pep$id <- c(1:nrow(pep)) # Note: these are new IDs, so will not match the peptide IDs in unprocessed MaxQuant output tables
ev$"Mod. peptide ID" <- NULL
ev$"Peptide ID" <- pep$id[pep_to_ev]
# Amino Acid counts
sq <- pep$Sequence
clusterExport(parClust, "sq", envir = environment())
tmp <- parSapply(parClust, proteoCraft::AA, function(aa) {
  nchar(sq) - nchar(gsub(aa, "", sq))
})
colnames(tmp) <- paste0(colnames(tmp), " Count")
pep[, colnames(tmp)] <- tmp
ev[, paste0(AA, " Count")] <- pep[pep_to_ev, paste0(AA, " Count")]
for (aa in c("O", "U")) { # Only keep the selenocysteine and pyrrolysine amino acid columns if they are non-empty (NB: pyrrolysine should really only be in some bacteria)
  if (sum(ev[[paste0(aa, " Count")]]) == 0) {
    ev[[paste0(aa, " Count")]] <- NULL
    pep[[paste0(aa, " Count")]] <- NULL
  }
}
# Length
pep$Length <- nchar(pep$Sequence)
ev$Length <- pep$Length[pep_to_ev]
# Intensity
temp <- aggregate(ev$Intensity, list(ev$"Modified sequence"), sum, na.rm = TRUE)
pep$Intensity <-  temp$x[match(pep$"Modified sequence", temp$Group.1)]
for (e in Exp) { #e <- Exp[1]
  temp <- ev[which(ev$Experiment == e),]
  temp <-  set_colnames(aggregate(temp[[int.col]], list(temp$"Modified sequence"), function(x) { sum(is.all.good(x)) }),
                        c("Modified sequence", int.col))
  pep[, paste0(int.col, " - ", e)] <- 0
  w <- which(pep$"Modified sequence" %in% temp$"Modified sequence")
  pep[w, paste0(int.col, " - ", e)] <- temp[match(pep$"Modified sequence"[w], temp$"Modified sequence"), int.col]
}
if (PTMriched) {
  for (ptm in EnrichedPTMs) { pep[[ptm]] <- ev[ev_to_pep, ptm] }
}

#### Code chunk - optionally impute missing expression values
if ((length(Exp) > 1)&&(ImputeMissData)) {
  kol <- grep(topattern(int.col), colnames(pep), value = TRUE)
  tst <- length(which(!is.all.good(log10(unlist(pep[, kol])), 2)))
  if (length(tst)) {
    cat("Incomplete peptides-level expression values matrix.\nImputing missing values with random draws from a gaussian distribution centered on the lowest observed value and with SD = 1/5 that of the data.\n")
    g <- grep(topattern(paste0(int.col, " - ")), colnames(pep), value = TRUE)
    temp <- log10(as.matrix(pep[, g]))
    w1 <- which((is.na(temp) | is.infinite(temp)), arr.ind = TRUE)
    w2 <- which(!(is.na(temp) | is.infinite(temp)), arr.ind = TRUE)
    Min <- min(unlist(temp[w2]))
    SD <- sd(unlist(temp[w2]))
    temp[w1] <- rnorm(nrow(w1), Min, SD/5) # Subsetting with arrays is sooo cool!
    kol2 <- gsub(topattern(int.col), paste0("Imput. ", int.cols["Original"]), g)
    colnames(temp) <- kol2
    temp <- 10^temp
    temp[w2] <- pep[, g][w2]
    temp <- as.data.frame(temp)
    temp[[paste0("Imput. ", int.cols["Original"])]] <- rowSums(temp[, kol2])
    pep[, colnames(temp)] <- temp
    int.cols["Imputed"] <- int.col <- paste0("Imput. ", int.cols["Original"])
  }
}

## Peptides level:
# Intensity distribution:
kol <- c("Modified sequence", paste0(int.col, " - ", Exp))
kol2 <- "Modified sequence"
form <- ".~variable"
if (PTMriched) {
  kol <- c(kol, EnrichedPTMs)
  kol2 <- c(kol2, EnrichedPTMs)
  form <- gsub("^\\.", "PTMs", form)
}
form <- as.formula(form)
temp <- pep[, kol]
temp <- reshape::melt.data.frame(temp, id.vars = kol2)
temp$variable <- gsub(topattern(paste0(int.col, " - ")), "", temp$variable)
temp$variable <- factor(temp$variable, levels = Exp)
temp$value <- log10(temp$value)
temp <- temp[which(is.all.good(temp$value, 2)),]
if (PTMriched) {
  for (ptm in EnrichedPTMs) { temp[[ptm]] <- c("", ptm)[temp[[ptm]]+1] }
  if (length(EnrichedPTMs) > 1) {
    temp$PTMs <- do.call(paste, c(temp[, EnrichedPTMs], sep = "-"))
  } else {
    temp$PTMs <- temp[[EnrichedPTMs]]
  }
}
ttl <- "Density plot - Peptides level"
plot <- ggplot(temp) + geom_histogram(aes(x = value, fill = variable), bins = 100) +
  facet_grid(form) +
  theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
  ggtitle(ttl) + xlab("log10(Peptides Intensity)")
poplot(plot)
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
# Correlation:
if (length(Exp) > 1) {
  temp <- pep[, c("Modified sequence", paste0(int.col, " - ", Exp))]
  for (e in Exp) {
    temp[[paste0("log10(Intensity) - ", e)]] <- log10(temp[[paste0(int.col, " - ", e)]])
  }
  comb <- as.data.frame(gtools::combinations(length(Exp), 2, Exp))
  temp2 <- temp[, grep(topattern("log10(Intensity) - "), colnames(temp), value = TRUE)]
  source(parSrc, local = FALSE)
  clusterExport(parClust, "temp2", envir = environment())
  temp2 <- parApply(parClust, comb, 1, function(x) {
    temp3 <- temp2[, paste0("log10(Intensity) - ", unlist(x))]
    temp3$X <- x[[1]]
    temp3$Y <- x[[2]]
    temp3$Comparison <- paste0(x[[1]], " (X) vs ", x[[2]], " (Y)")
    colnames(temp3)[1:2] <- c("log10(X intensity)", "log10(Y intensity)")
    return(temp3)
  })
  temp2 <- plyr::rbind.fill(temp2)
  temp2$"Modified sequence" <- temp$"Modified sequence"
  test <- parApply(parClust, temp2[,c("log10(X intensity)", "log10(Y intensity)")], 1, function(x) {
    length(proteoCraft::is.all.good(x))
  }) == 2
  temp2 <- temp2[which(test),]
  temp2$X <- factor(temp2$X, levels = Exp)
  temp2$Y <- factor(temp2$Y, levels = Exp)
  temp3 <- as.data.frame(t(sapply(unique(temp2$Comparison), function(x) { #x <- unique(temp2$Comparison)[1]
    x1 <- temp2[which(temp2$Comparison == unlist(x)), c("log10(X intensity)", "log10(Y intensity)")]
    x1 <- x1$"log10(Y intensity)"-x1$"log10(X intensity)"
    return(setNames(c(x, paste0("Median = ", round(median(x1), 3)), paste0("S.D. = ", round(sd(x1), 3))),
                    c("Comparison", "Median", "SD")))
  })))
  temp3$X <- gsub(" \\(X\\).+", "", temp3$Comparison)
  temp3$Y <- gsub(".+\\(X\\) vs ", "", gsub(" \\(Y\\)$", "", temp3$Comparison))
  temp3$X <- factor(temp3$X, levels = Exp)
  temp3$Y <- factor(temp3$Y, levels = Exp)
  temp3$R <- vapply(strsplit(gsub(" \\([XY]\\)", "", temp3$Comparison), " vs "), function(x) {
    return(paste0("R = ", round(cor(pep[[paste0(int.col, " - ", x[[1]])]],
                                    pep[[paste0(int.col, " - ", x[[2]])]]), 3)))
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
  clusterExport(parClust, list("get_density", "tmp2"), envir = environment())
  comps <- unique(temp2$Comparison)
  tmp2D <- setNames(parLapply(parClust, comps, function(cmp) { #cmp <- comps[1]
    w <- which(tmp2$Comparison == cmp)
    x <- try({
      y <- get_density(tmp2$`log10(X intensity)`[w], tmp2$`log10(Y intensity)`[w], n = 500)
      list(Success = TRUE, Density = data.frame(Density = y, Seq = tmp2$"Modified sequence"[w]))
    }, silent = TRUE)
    if ("try-error" %in% class(x)) { x <- list(Success = FALSE) }
    return(x)
  }), comps)
  tmp2D <- tmp2D[which(vapply(tmp2D, function(x) { x$Success }, TRUE))]
  tmp2D <- lapply(tmp2D, function(x) { x$Density })
  if (length(tmp2D)) {
    comps <- names(tmp2D)
    temp2 <- temp2[which(temp2$Comparison %in% comps),]
    temp2$Density <- 0
    for (cmp in comps) { #cmp <- comps[1]
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
      facet_grid(Y~X, drop = TRUE) + coord_fixed(1) +
      scale_color_viridis() +
      theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
                         strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) + ggtitle(ttl)
    poplot(plot)
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
  } else { warning("Do we really have enough data to continue? Investigate...") }
}
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))

# Calculate peptide ratios
rat.col %<o% "log2(Ratio)"
rat.cols %<o% c()
if (MakeRatios) {
  rat.cols["Original"] <- rat.col
  rat.grps %<o% unique(SamplesMap$`Ratios group`)
  rat.grps <- rat.grps[which(!is.na(rat.grps))]
  for (gr in rat.grps) { #gr <- rat.grps[1]
    m <- SamplesMap[which(SamplesMap$`Ratios group` == gr),]
    if (sum(c(TRUE, FALSE) %in% m$Reference) < 2) {
      warning(paste0("Ratios group ", gr, " - there should be reference (control) and non-reference samples!"))
    } else {
      ref <- apply(pep[, paste0(int.col, " - ", m$Experiment[which(m$Reference)]), drop = FALSE], 1, function(x) {
        x <- is.all.good(x)
        l <- length(x)
        if (!l) { x <- NA } else { x <- prod(x)^(1/l) }
        return(x)
      })
      w <- which(!m$Reference)
      for (x in w) { #x <- w[1]
        pep[[paste0(rat.col, " - ", m$Experiment[x])]] <- log2(pep[[paste0(int.col, " - ", m$Experiment[x])]]/ref)
      }
    }
  }
  # Peptide ratios distribution:
  dir <- paste0(wd, "/Workflow control")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  temp <- pep[, c("Modified sequence", grep(topattern(paste0(rat.col, " - ")), colnames(pep), value = TRUE))]
  temp <- reshape::melt.data.frame(temp, id.vars = "Modified sequence")
  temp$variable <- gsub(topattern(paste0(rat.col, " - ")), "", temp$variable)
  temp <- temp[which(is.all.good(temp$value, 2)),]
  ttl <- "Ratios density plot - Peptides level"
  plot <- ggplot(temp) + geom_histogram(aes(x = value, fill = variable), bins = 100) +
    geom_vline(xintercept = RatiosThresh, colour = "red") +
    ggtitle(ttl) + facet_grid(.~variable) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
    theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
    xlab("log2(Ratio)")
  if (RatiosThresh_2sided) { plot <- plot + geom_vline(xintercept = -RatiosThresh, colour = "red") }
  poplot(plot)
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
}

Script <- readLines(ScriptPath)
gc()
# It makes sense to close/re-create parallel clusters regularly to reduce memory usage + avoid corruption
stopCluster(parClust)
source(parSrc, local = FALSE)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - Assemble protein groups
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
# Basic fix, because I do not like the way I was doing Quality filters up to now
g <- grep("^Quality filter: ", colnames(PG), value = TRUE)
if (length(g)) {
  for (h in g) { #h <- g[1]
    PG[[h]] <- c("no -> dubious!", "")[match(PG[[h]], c("", "Keep"))]
  }
}
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
  test <- aggregate(test$Org, list(test$L1), function(x) { paste(unique(x), collapse = ";") })
  pgOrgKol %<o% c("Organism", "Organism(s)")[(sum(grepl(";", test))>0)+1]
  PG[[pgOrgKol]] <- test$x[match(PG$id, test$Group.1)]
}

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
poplot(plot)
dir <- paste0(wd, "/Summary plots")
#dirlist<- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300)
ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)

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
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
# Simplify Gene columns
genkol <- c("Genes", "Gene names")
w <- which(genkol %in% colnames(PG))
if (length(w) == 2) {
  temp <- PG[, genkol]
  for (i in genkol) { temp[[i]] <- strsplit(temp[[i]], ";") }
  PG$Genes <- apply(temp, 1, function(x) { paste(sort(unique(unlist(x))), collapse = ";") })
  PG$"Gene names" <- NULL
} else { if (length(w) == 1) { colnames(PG)[which(colnames(PG) %in% genkol)] <- "Genes" } }

# Here, if this is a BioID type experiment, we also want to mark protein groups which have Biotin peptides:
if (IsBioID) {
  wbiot <- grep("biot", Modifs$"Full name", ignore.case = TRUE)
  l <- length(wbiot)
  if (length(wbiot)) {
    if (l == 1) { tmp <- Modifs$"Full name"[wbiot] } else { tmp <- paste0(paste(Modifs$"Full name"[wbiot[1:(l-1)]], collapse = "\", \""), "\" and \"", Modifs$"Full name"[wbiot[l]]) }
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
    temp <- aggregate(temp$value, list(temp$L1), function(x) { paste(sort(x), collapse = ";") })
    PG[wpg, "Biot. peptide IDs"] <- temp$x[match(PG$id[wpg], temp$Group.1)]
    PG[["Biot. peptides count"]] <- vapply(strsplit(PG[["Biot. peptide IDs"]], ";"), length, 1)
    PG[["Biot. peptides [%]"]] <- round(100*PG[["Biot. peptides count"]]/PG$"Peptides count", 1)
    IsBioID2 <- TRUE
  } else {
    warning("I could not find any biotinylated peptides!")
    IsBioID <- FALSE
    AnalysisParam$"Type - advanced" <- "BioID... - except that no biotinylated peptides were detected!"
  }
}

# Number of spectra, evidences and peptides per sample:
source(parSrc, local = FALSE)
invisible(clusterCall(parClust, function() {
  library(proteoCraft)
  library(reshape)
  library(data.table)
  return()
}))
temp_PG <- data.frame(id = PG$id,
                      Accession1 = vapply(strsplit(PG$"Leading protein IDs", ";"), function(x) { unlist(x)[1] }, ""))
temp_PG$Pep <- parLapply(parClust, strsplit(PG$"Peptide IDs", ";"), as.integer)
tmp <- pep[, c("id", "Sequence")]
clusterExport(parClust, "tmp", envir = environment())
temp_PG$Pep <- parLapply(parClust, temp_PG$Pep, function(x) { tmp$Sequence[match(x, tmp$id)] })
temp_PG$Seq <- db$Sequence[match(temp_PG$Accession1, db$"Protein ID")]
if (!"Sequence coverage [%]" %in% colnames(PG)) {
  exports <- list("Coverage")
  clusterExport(parClust, exports, envir = environment())
  PG$"Sequence coverage [%]" <- round(100*parApply(parClust, temp_PG[, c("Seq", "Pep")], 1, function(x) {
    Coverage(x[[1]], x[[2]])
  }), 1)
}
CreateMSMSKol %<o% (("MS/MS IDs" %in% colnames(ev))&&(class(ev$"MS/MS IDs") %in% c("integer", "character")))
if (CreateMSMSKol) {
  # Somehow there are no MSMS IDs for DIA, but this may be a late tables writing bug.
  ev$temp <- parLapply(parClust, strsplit(as.character(ev$"MS/MS IDs"), ";"), as.integer)
  #PG[, paste0("Spectr", c("al count", "um IDs"))]
  temp <- reshape2::melt(setNames(lapply(strsplit(PG$`Evidence IDs`, ";"), as.integer), PG$id))
  temp$MSMSIDs <- ev$temp[match(temp$value, ev$id)]
  temp <- temp[which(vapply(temp$MSMSIDs, length, 1) > 0),] # Remove Match-Between-Runs evidences (no MS/MS)
  temp <- reshape2::melt(setNames(temp$MSMSIDs, temp$L1))
  temp <- do.call(data.frame, aggregate(temp$value, list(temp$L1), function(x) {
    x <- unique(x)
    return(c(Count = length(x), List = list(x)))
  }))
  temp$x.Count <- unlist(temp$x.Count)
  temp$Pasted <- vapply(temp$x.List, function(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") }, "")
  PG[, paste0("Spectr", c("al count", "um IDs"))] <- temp[match(PG$id, temp$Group.1), c("x.Count", "Pasted")]
  ev$temp <- NULL
}
temp_ev <- ev[, c("id", "Experiment", "Protein group IDs", "Peptide ID")]
if (CreateMSMSKol) { temp_ev$"MS/MS IDs" <- ev$"MS/MS IDs" }
temp_pep <- pep[, c("id", "Sequence")]
clusterExport(parClust, exports, envir = environment())
exports <- list("temp_ev", "temp_pep", "temp_PG", "IsBioID2", "Exp", "Modifs", "CreateMSMSKol")
if (IsBioID2) { exports <- append(exports, "wbiot") }
clusterExport(parClust, exports, envir = environment())
temp <- parLapply(parClust, Exp, function(exp) { #exp <- Exp[1]
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
  res[, kolk] <- 0
  res[, koli] <- ""
  res[[paste0("Sequence coverage [%] - ", exp)]] <- 0
  w <- which(temp_ev$Experiment == exp)
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
    res[w, paste0("Sequence coverage [%] - ", exp)] <- round(100*apply(temp_PG[w, c("Seq", "Pep")], 1, function(x) {
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
        temp1 <- reshape2::melt(setNames(lapply(strsplit(eB$"Protein group IDs", ";"), as.integer), eB$"Peptide ID"))
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

Script <- readLines(ScriptPath)
gc()
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - Calculate protein group-level quantitative values
PG.int.cols %<o% c()
PG.int.col %<o% "log10(Expr.) - "
PG.int.cols["Original"] <- PG.int.col
PG.rat.cols %<o% c()
PG.rat.col %<o% "log2(Ratio) - " #We need those defaults actually even if !MakeRatios
PG.rat.cols["Original"] <- PG.rat.col
if (ImputeMissData) {
  PG.int.cols["Imputed"] <- PG.int.col <- paste0("Imput. ", PG.int.cols["Original"])
  PG.rat.cols["Imputed"] <- PG.rat.col <- paste0("Imput. ", PG.rat.cols["Original"])
}

# Coverage columns
PG$"1st accession" <- vapply(strsplit(PG$`Leading protein IDs`, ";"), function(x) { unlist(x)[1] }, "")
PG$"Sequence (1st accession)" <- db$Sequence[match(PG$`1st accession`, db$`Protein ID`)]
source(parSrc, local = FALSE)
clusterExport(parClust, "Coverage", envir = environment())
for (exp in Exp) { #exp <- Exp[1]
  temp <- PG[, c("1st accession", "Sequence (1st accession)", paste0("Peptide IDs - ", exp))]
  w <- which(temp[[paste0("Peptide IDs - ", exp)]] != "")
  temp <- temp[w,]
  tmp <- listMelt(strsplit(temp[[paste0("Peptide IDs - ", exp)]], ";"), temp$`1st accession`)
  tmp$ModSeq <- pep$`Modified sequence`[match(as.integer(tmp$value), pep$id)]
  tmp <- aggregate(tmp$ModSeq, list(tmp$L1), c)
  temp$Peptides <- tmp$x[match(temp$`1st accession`, tmp$Group.1)]
  names(temp$"Sequence (1st accession)") <- temp$"1st accession"
  tmp <- temp[, c("Sequence (1st accession)", "Peptides")]
  clusterExport(parClust, "tmp", envir = environment())
  temp$Coverage <- parApply(parClust, tmp, 1, function(x) {
    round(100*Coverage(x[[1]], x[[2]]), 1)
  })
  PG[[paste0("Sequence coverage [%] - ", exp)]] <- 0
  PG[w, paste0("Sequence coverage [%] - ", exp)] <- temp$Coverage
}

# LFQ quant and ratios:
Exp.map <- SamplesMap
Exp.map$Ref.Sample.Aggregate <- SamplesMap$Experiment
Exp.map$Ratios_group <- paste0("Group", Exp.map$`Ratios group`)
RatGrp <- RefGrp <- list(aggregate = "Rat",
                         values = unique(Exp.map$Ratios_group),
                         names = "Ratios_group",
                         column = "Ratios_group")
SmplGrp <- list(aggregate = "Exp",
                values = Exp,
                names = "Experiment",
                column = "Experiment")
Aggregate.map <- list(Aggregate.name = c("Exp", "Rat"), Characteristics = c("Experiment", "Ratios_group"))
Aggregate.list <- list(Rat = unique(Exp.map$Ratios_group),
                       Exp = unique(Exp.map$Experiment))
Aggregates <- setNames(c("Ratios_group", "Experiment"), c("Rat", "Exp"))
ref.aggr.col <- RefGrp$names
quant.data %<o% Prot.Quant(PG, Pep4Quant, "PreferUnique", pep,
                           Summary.method = "mean", #Summary.weights = "Weights",
                           Intensity.weights = FALSE, Skip.ratios = !MakeRatios,
                           experiments.map = Exp.map, ref.groups = RefGrp, ratio.groups = RatGrp,
                           sample.groups = SmplGrp,
                           Pep.Intens.root = paste0(int.col, " - "), Pep.Ratios.root = paste0(rat.cols["Original"], " - "),
                           log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                           Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                           Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                           Min.N = NPep, N.clust = N.clust,
                           Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                           cl = parClust)
quant.data <- quant.data[, which(!grepl("\\.REF$", colnames(quant.data))), drop = FALSE]
colnames(quant.data) <- gsub("^log10\\(Expr\\.\\) - ", PG.int.col, colnames(quant.data))
colnames(quant.data) <- gsub("^log2\\(Ratio\\) - ", PG.rat.col, colnames(quant.data))
colnames(quant.data)[which(colnames(quant.data) == "Peptides IDs used for quantitation")] <- paste0("Peptide IDs used for quantitation - ",
                                                                                                    names(PG.int.cols)[match(PG.int.col, PG.int.cols)])
saveFun(quant.data, file = "quant.data.RData")
#loadFun("quant.data.RData")
#quant.data <- TopN(3, PG, Pep4Quant, pep, "id",
#             Pep.Intens.Nms = grep(topattern(paste0(int.col, " - ")), colnames(pep), value = TRUE),
#             log.Pep.Intens = FALSE,
#             Mods = Mod4Quant,
#             Min.Pep.Nb = NPep, corr = "global", Out.Norm = FALSE)
#PG.int.col <- paste0("log10(", int.col, ") - ")
#colnames(quant.data) <- gsub(paste0("^log10 - ", int.col, " - "), PG.int.col, colnames(temp))
#loadFun("quant.data.RData")

if (ImputeMissData) {
  # In that case we need to calculate expression values a second time, unfortunately... It doesn't take so long.
  refI <- int.cols[grep("^Imput\\.", int.cols)-1]
  refR <- rat.cols[grep("^Imput\\.", int.cols)-1] # This is correct!
  PGrefI <- PG.int.cols[grep("^Imput\\.", PG.int.cols)-1]
  PGrefR <- PG.rat.cols[grep("^Imput\\.", PG.int.cols)-1] # Again, intentional.
  quant.data2 %<o% Prot.Quant(PG, Pep4Quant, "PreferUnique", pep,
                              Summary.method = "mean", #Summary.weights = "Weights",
                              Intensity.weights = FALSE, Skip.ratios = !MakeRatios,
                              experiments.map = Exp.map, ref.groups = RefGrp, ratio.groups = RatGrp,
                              sample.groups = SmplGrp,
                              Pep.Intens.root = paste0(refI, " - "), Pep.Ratios.root = paste0(refR, " - "),
                              log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                              Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                              Mods = Mod4Quant,
                              Min.N = NPep, N.clust = N.clust,
                              Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              cl = parClust)
  #colnames(quant.data2)
  quant.data2 <- quant.data2[, which(!grepl("\\.REF$", colnames(quant.data2))), drop = FALSE]
  colnames(quant.data2) <- gsub("^log10\\(Expr\\.\\) - ", PGrefI, colnames(quant.data2))
  colnames(quant.data2) <- gsub("^log2\\(Ratio\\) - ", PGrefR, colnames(quant.data2))
  colnames(quant.data2)[which(colnames(quant.data2) == "Peptides IDs used for quantitation")] <- paste0("Peptide IDs used for quantitation - ",
                                                                                                        names(PG.int.cols)[match(PGrefI, PG.int.cols)])
  save(quant.data2, file = "quant.data2.RData")
  PG[, colnames(quant.data2)] <- quant.data2
}
PG[, colnames(quant.data)] <- quant.data

rm(Exp.map)

#### Code chunk - Re-normalize protein group expression values
# Normalize (Levenberg-Marquardt)
if ((length(Exp) > 1)&&(NormalizePG)) {
  g <- grep(topattern(PG.int.col), colnames(PG), value = TRUE)
  temp <- PG[, c("id", g)]
  m <- apply(temp[,g], 2, function(x) { median(is.all.good(x)) })
  M <- median(is.all.good(unlist(temp[,g])))
  temp[,g] <- sweep(temp[,g], 2, m, "-") + M
  temp <- AdvNorm.IL(temp[, c("id", g)], "id", exprs.col = g, exprs.log = TRUE)
  PG[, gsub(topattern(PG.int.col), paste0("Norm. ", PG.int.cols["Original"]), g)] <- temp[, paste0("AdvNorm.", g)]
  PG.int.cols["Normalized"] <- PG.int.col <- paste0("Norm. ", PG.int.cols["Original"])
  if (MakeRatios) {
    for (gr in unique(SamplesMap$`Ratios group`)) { #gr <- unique(SamplesMap$`Ratios group`)[1]
      m <- SamplesMap[which(SamplesMap$`Ratios group` == gr),]
      rf <- m$Experiment[which(m$Reference)]
      for (i in m$Experiment[which(!m$Reference)]) { #i <- m$Experiment[which(!m$Reference)][1]
        PG[[paste0("Norm. ", rat.cols["Original"], " - ", i)]] <- (PG[[paste0(PG.int.col, i)]] - PG[[paste0(PG.int.col, rf)]])/log10(2)
      }
    }
    PG.rat.cols["Normalized"] <- PG.rat.col <- paste0("Norm. ", PG.rat.cols["Original"])
  }
}

# Average expression columns:
if (length(Exp) > 1) {
  for (i in PG.int.cols) {
    PG[[paste0("Mean ", gsub(" - $", "", i))]] <- apply(PG[, grep(topattern(i), colnames(PG), value = TRUE)],
                                                        1, mean, na.rm = TRUE)
  }
}


#### Code chunk - Prepare Annotations and (if applicable) GO terms
# NB: I used to get functional annotations for all proteins in the protein group.
#     However we are now - and I think with reason - only using annotations from the leading protein(s)!
setwd(wd)
Src <- paste0(libPath, "/extdata/R scripts/Sources/Annotate_me.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# GO filters
if (GO_filt) {
  for (goID in GO_filter) { #goID <- GO_filter[1]
    # Get children terms
    gofilter <- unique(unlist(c(goID,
                                GOBPOFFSPRING[[goID]],
                                GOCCOFFSPRING[[goID]],
                                GOMFOFFSPRING[[goID]])))
    gofilter <- gofilter[which(!is.na(gofilter))]
    if (sum(gofilter %in% AllTerms)) {
      PG[[goID]] <- ""
      PG[grsep2(gofilter, PG$`GO-ID`), goID] <- "+"
      #which(vapply(strsplit(PG$`GO-ID`[-wtst],";"), function(x) { sum(x %in% gofilter) }, 1) != 0)
      #which(vapply(strsplit(PG$`GO-ID`[wtst],";"), function(x) { sum(x %in% gofilter) }, 1) == 0)
    }
  }
}

#### Code chunk - Correlation and distribution plots
# LFQ correlation plots:
# Idea for this plot: map color to protein type (list, organism or contaminant) as in profile plots
dir <- paste0(wd, "/Workflow control")
if (length(Exp) > 1) {
  comb <- as.data.frame(gtools::combinations(length(Exp), 2, Exp))
  kount <- 0
  temp <- list()
  for (klnm in names(PG.int.cols)) { #klnm <- names(PG.int.cols)[1]
    kol <- PG.int.cols[klnm]
    kolZ <- grep(topattern(kol), colnames(PG), value = TRUE)
    if (length(kolZ)) {
      kount <- kount + 1
      temp2 <- PG[, kolZ] 
      source(parSrc, local = FALSE)
      clusterExport(parClust, list("kol", "temp2", "klnm"), envir = environment())
      temp2 <- parApply(parClust, comb, 1, function(x) {
        temp3 <- temp2[, paste0(kol, unlist(x))]
        temp3$X <- x[[1]]
        temp3$Y <- x[[2]]
        temp3$Comparison <- paste0(x[[1]], " (X) vs ", x[[2]], " (Y)")
        colnames(temp3)[1:2] <- c("log10(X LFQ)", "log10(Y LFQ)")
        temp3$Type <- klnm
        return(temp3)
      })
      temp2 <- plyr::rbind.fill(temp2)
      temp2$"Common Name (short)" <- PG$"Common Name (short)"
      temp2$Type <- klnm
      temp[[kount]] <- temp2
    }
  }
  temp <- plyr::rbind.fill(temp)
  test <- parApply(parClust, temp[,c("log10(X LFQ)", "log10(Y LFQ)")], 1, function(x) {
    length(proteoCraft::is.all.good(x))
  }) == 2
  temp <- temp[which(test),]
  w <- aggregate(1:nrow(temp), list(temp$Type), list)
  temp2 <- data.frame(Comparison = rep(unique(temp$Comparison), length(PG.int.cols)))
  temp2$Type <- as.character(sapply(names(PG.int.cols), function(x) { rep(x, length(unique(temp$Comparison)))}))
  temp2[, c("Median", "SD", "R")] <- as.data.frame(t(apply(temp2[, c("Type", "Comparison")], 1, function(x) { #x <- temp2[1, c("Type", "Comparison")]
    x1 <- temp[which((temp$Type == x[[1]])&(temp$Comparison == x[[2]])), c("log10(X LFQ)", "log10(Y LFQ)")]
    x2 <- x1$"log10(Y LFQ)"-x1$"log10(X LFQ)"
    return(c(paste0("Median = ", round(median(x2), 3)),
             paste0("S.D. = ", round(sd(x2), 3)),
             paste0("R = ", round(cor(x1$"log10(X LFQ)", x1$"log10(Y LFQ)"), 3))))
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
  clusterExport(parClust, list("get_density", "tmp"), envir = environment())
  comps <- unique(temp$Comparison)
  tmpD <- setNames(parLapply(parClust, comps, function(cmp) { #cmp <- comps[1]
    w <- which(tmp$Comparison == cmp)
    x <- get_density(tmp$`log10(X LFQ)`[w], tmp$`log10(Y LFQ)`[w], n = 500)
    x <- data.frame(Density = x, Nm = tmp$"Common Name (short)"[w])
    return(x)
  }), comps)
  temp$Intensity <- 0
  for (cmp in comps) { #cmp <- comps[1]
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
    plot <- plot + geom_point(data = temp[which(temp$`In list` == "+"),], aes(x = `log10(X LFQ)`, y = `log10(Y LFQ)`), size = 1, shape = 1, color = "red") +
      geom_text(data = temp[which(temp$`In list` == "+"),], aes(label = `Common Name (short)`, x = `log10(X LFQ)`, y = `log10(Y LFQ)`),
                size = 2, hjust = 0, vjust = 0, color = "red")
  }
  poplot(plot, 12, 20)
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
}
g <- as.character(sapply(PG.int.cols, function(x) { grep(topattern(x), colnames(PG), value = TRUE) }))
test <- vapply(g, function(x) { length(is.all.good(PG[[x]])) }, 1)
print(test)
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))

# Intensities distribution:
long.dat <- list()
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
temp <- PG[, c("Common Name (short)", as.character(sapply(PG.int.cols, function(x) { paste0(x, Exp) })))]
temp <- reshape::melt.data.frame(temp, id.vars = "Common Name (short)")
colnames(temp)[which(colnames(temp) == "value")] <- "log10(Intensity)"
if (length(Exp) > 1) {
  temp2 <- PG[, c("Common Name (short)", paste0("Mean ", gsub(" - $", "", PG.int.cols)))]
  temp2 <- reshape::melt.data.frame(temp2, id.vars = "Common Name (short)")
  temp$"Mean log10(Intensity)" <- temp2$value
}
temp$Experiment <- gsub(paste0(".*", topattern(PG.int.cols["Original"], start = FALSE)), "", temp$variable)
temp$Experiment <- factor(temp$Experiment, levels = SamplesMap$Experiment)
temp$Type <- gsub(paste0(" ?", topattern(PG.int.cols["Original"], start = FALSE), ".*$"), "", temp$variable)
temp$Type[which(temp$Type == "")] <- "Orig."
temp$Type <- paste0("log10(", tolower(temp$Type), " LFQ)")
temp$Type <- factor(temp$Type, levels = paste0("log10(", c("orig", "imput", "norm"), ". LFQ)"))
long.dat$intens <- temp
temp <- temp[which(is.all.good(temp$"log10(Intensity)", 2)),]
ttl <- "LFQ density plot - PGs level"
plot <- ggplot(temp) + geom_histogram(aes(x = `log10(Intensity)`, fill = Type), bins = 100) +
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
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
  xlab("log10(Protein Groups LFQ)")
#poplot(plot)
M <- c()
if (GO_filt) {
  CompGOTerms %<o%setNames(c("GO:0005634", "GO:0005654", "GO:0000785", "GO:0005730", "GO:0005635", "GO:0005737", "GO:0005829", "GO:0005783", "GO:0005794", "GO:0031988", "GO:0005739", "GO:0009536", "GO:0005886", "GO:0031012", "GO:1903561"),
                           c("Nucleus", "Nucleoplasm", "Chromatin", "Nucleolus", "Nuclear envelope", "Cytoplasm", "Cytosol", "ER", "Golgi", "Vesicle", "Mitochondrion", "Plastid", "Plasma membrane", "Extracellular matrix", "Extracellular vesicle"))
  GO_filter1 %<o% unique(c(GO_filter, CompGOTerms))
  for (goID in GO_filter1) { #goID <- GO_filter1[1]
    # Get children terms
    gofilter <- unique(unlist(c(goID,
                                GOBPOFFSPRING[[goID]],
                                GOCCOFFSPRING[[goID]],
                                GOMFOFFSPRING[[goID]])))
    gofilter <- gofilter[which(!is.na(gofilter))]
    if (sum(gofilter %in% AllTerms)) {
      temp[[goID]] <- 0
      wtst <- grsep2(gofilter, PG$`GO-ID`)
      #which(vapply(strsplit(PG$`GO-ID`[-wtst],";"), function(x) { sum(x %in% gofilter) }, 1) != 0)
      #which(vapply(strsplit(PG$`GO-ID`[wtst],";"), function(x) { sum(x %in% gofilter) }, 1) == 0)
      wtst2 <- which(temp$"Common Name (short)" %in% PG$"Common Name (short)"[wtst])
      temp[wtst2, goID] <- 1
      w <- which(temp$Type == "log10(orig. LFQ)")
      wtst3 <- which(temp$"Common Name (short)"[w] %in% PG$"Common Name (short)"[wtst])
      if (length(wtst3) > 1) {
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
  temp2 <- reshape::melt.data.frame(temp2, id.vars = c("log10(Intensity)", "Experiment"))
  colnames(temp2) <- c("log10(Intensity)", "Experiment", "GO term", "+")
  temp2 <- temp2[which(temp2$"+" == 1),]
  if (globalGO) {
    temp2$`GO term` <- GO_terms$Term[match(temp2$`GO term`, GO_terms$ID)]
    temp2$`GO term` <- gsub("\\]$", "", gsub(" \\[", "\n",  temp2$`GO term`))
  }
  ttl2 <- "LFQ density plot - PGs level, GO terms"
  plot3 <- ggplot(temp2) + geom_density(stat = "density", aes(x = `log10(Intensity)`, fill = `GO term`),
                                        colour = "black", alpha = 0.5) +
    geom_density(data = temp, stat = "density", aes(x = `log10(Intensity)`), colour = "red",
                 linewidth = 1, linetype = "dotted", show.legend = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
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
  ggsave(paste0(dir, "/", ttl2, ".jpeg"), plot3a, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl2, ".pdf"), plot3a, dpi = 300, width = 10, height = 10, units = "in")
  leg3 <- get_legend(plot3)
  plot3 <- plot3 + theme(legend.position = "none")
  g2 <- ggplotGrob(plot2)
  g3 <- ggplotGrob(plot3)
  maxWidth <- grid::unit.pmax(g2$widths[2:5], g3$widths[2:5])
  g2$widths[2:5] <- as.list(maxWidth)
  g3$widths[2:5] <- as.list(maxWidth)
  plot4 <- arrangeGrob(grobs = list(g2, leg2, g3, leg3),
                       widths = c(4, 1),
                       heights = c(length(unique(temp$Type)), length(M)/5),
                       layout_matrix = rbind(c(1, 2), c(3, 4)),
                       padding = 5)
  plot4 <- as.ggplot(plot4)
  poplot(plot4)
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot4, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot4, dpi = 300, width = 10, height = 10, units = "in")
  #poplot(plot3)
} else {
  poplot(plot)
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
}

# Gene-Set Enrichment Analysis (GSEA)
if (runGSEA) {
  dataType <- "PG"
  Src <- paste0(libPath, "/extdata/R scripts/Sources/GSEA.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}
#

if (MakeRatios) {
  FC_filt %<o% c()
  FC_Smpls %<o% list()
  # Fold change filters:
  ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1]
  rat.grps <- unique(SamplesMap$`Ratios group`)
  rat.grps <- rat.grps[which(!is.na(rat.grps))]
  for (grp in rat.grps) { #grp <- rat.grps[1]
    SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
    smpl0 <- SmplMp$Experiment[which(SmplMp$Reference)]
    stopifnot(length(smpl0) == 1)
    smpl1 <- SmplMp$Experiment[which(!SmplMp$Reference)]
    e0 <- PG[[paste0(ref, smpl0)]]
    if (NegFilt) {
      nf <- SmplMp$Experiment[which(SmplMp$"Negative Filter")]
      if (length(nf)) {
        smpl1 <- smpl1[which(!smpl1 %in% nf)]
        nftst <- apply(PG[, paste0(ref, nf), drop = FALSE], 1, function(x) { length(is.all.good(x)) }) > 0
      }
    }
    FC_filt <- append(FC_filt, setNames(lapply(smpl1, function(x) {
      e1 <- PG[[paste0(ref, x)]]
      r1 <- PG[[paste0(PG.rat.col, x)]]
      if (RatiosThresh_2sided) { r1 <- abs(r1) }
      w <- which(((r1 >= RatiosThresh)&(is.all.good(e1, 2)))|(is.all.good(e1, 2)&(!is.all.good(e0, 2))))
      if ((NegFilt)&&(length(nf))) { w <- w[which(!nftst[w])] }
      return(w)
    }), smpl1))
    FC_Smpls[[grp]] <- list(Numerator = smpl1, Denominator = smpl0)
  }
  #g <- grep(topattern(PG.rat.col), colnames(PG), value = TRUE)
  #test <- setNames(vapply(g, function(x) { length(is.all.good(PG[[x]])) }, 1), gsub(topattern(PG.rat.col), "", g))
  #print(test)
  #
  # Ratios distribution:
  dir <- paste0(wd, "/Workflow control")
  temp <- PG[, c("Common Name (short)", as.character(sapply(PG.rat.cols, function(x) {
    grep(topattern(x), colnames(PG), value = TRUE)
  })))]
  temp <- reshape::melt.data.frame(temp, id.vars = "Common Name (short)")
  colnames(temp)[which(colnames(temp) == "value")] <- "log2(Ratio)"
  temp$Experiment <- gsub(paste0(".*", topattern(PG.rat.cols["Original"], start = FALSE)), "", temp$variable)
  temp$Experiment <- factor(temp$Experiment, levels = SamplesMap$Experiment)
  temp$Type <- gsub(paste0(" ?", topattern(PG.rat.cols["Original"], start = FALSE), ".*$"), "", temp$variable)
  temp$Type[which(temp$Type == "")] <- "Orig."
  temp$Type <- paste0("log2(", tolower(temp$Type), " ratio)")
  temp$Type <- factor(temp$Type, levels = paste0("log2(", c("orig", "norm", "imput"), ". ratio)"))
  long.dat$ratios <- temp
  temp <- temp[which(is.all.good(temp$"log2(Ratio)", 2)),]
  ttl <- "Ratios density plot - PGs level"
  plot <- ggplot(temp) + geom_histogram(aes(x = `log2(Ratio)`, fill = Type), bins = 100) +
    geom_vline(xintercept = RatiosThresh, colour = "red") +
    theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
    ggtitle(ttl) + facet_grid(Type~Experiment) +
    xlab("log2(Ratio)")
  if (RatiosThresh_2sided) { plot <- plot + geom_vline(xintercept = -RatiosThresh, colour = "red") }
  poplot(plot)
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
  # MA plots:
  temp <- long.dat$intens
  temp <- temp[which(temp$Experiment %in% long.dat$ratios$Experiment),]
  tst1 <- do.call(paste, c(temp[, c("Common Name (short)", "Experiment")], sep = "___"))
  tst2 <- do.call(paste, c(long.dat$ratios[, c("Common Name (short)", "Experiment")], sep = "___"))
  temp$"log2(Ratio)" <- long.dat$ratios$`log2(Ratio)`[match(tst1, tst2)]
  temp <- temp[which(is.all.good(temp$`log2(Ratio)`, 2)),]
  ttl <- "MA plots - PGs level"
  plot <- ggplot(temp) + geom_point(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`, colour = Type), size = 0.1) +
    geom_hline(yintercept = 0, linewidth = 0.8, linetype = "dashed") +
    ggtitle(ttl) + coord_fixed(log10(2)/log2(10)) + theme_bw() + facet_grid(Type~Experiment) +
    scale_color_viridis(option = "D", discrete = TRUE, begin = 0.25) +
    xlab("A = mean log10(Intensity)") + ylab("M = sample log2(Ratio)") +
    theme(strip.text.y = element_text(angle = 0))
  poplot(plot)
  ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
  # "Regulated/Enriched" columns
  ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1]
  for (grp in rat.grps) { #grp <- rat.grps[1]
    SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
    if (sum(c(TRUE, FALSE) %in% SmplMp$Reference) < 2) {
      stop("The reference column should include TRUE and FALSE values!")
    } else {
      w0 <- which(SmplMp$Reference)
      stopifnot(length(w0) == 1)
      w1 <- which(!SmplMp$Reference)
      xp0 <- SmplMp$Experiment[w0]
      for (x in w1) { #x <- w1[1]
        xp1 <- SmplMp$Experiment[x]
        wonly <- which(is.all.good(PG[[paste0(ref, xp1)]], 2)&(!is.all.good(PG[[paste0(ref, xp0)]], 2)))
        PG[[paste0("Regulated - ", xp1)]] <- ""
        if (RatiosThresh_2sided) {
          wup <- which(PG[[paste0(PG.rat.col, xp1)]] >= RatiosThresh)
          wdwn <- which(PG[[paste0(PG.rat.col, xp1)]] <= -RatiosThresh)
          PG[wup, paste0("Regulated - ", xp1)] <- "up"
          PG[wdwn, paste0("Regulated - ", xp1)] <- "down"
        } else {
          PG[[paste0("Enriched - ", xp1)]] <- c("", "+")[(PG[[paste0(PG.rat.col, xp1)]] >= RatiosThresh)+1]
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
    PTMs_intRf %<o% rev(int.cols[which(int.cols != paste0("Imput. ", int.cols["Original"]))])[1]
    PTMs_ratRf %<o% rev(rat.cols[which(rat.cols != paste0("Imput. ", rat.cols["Original"]))])[1]
    PTMs_intNm0 <- names(PTMs_intRf)
    PTMs_ratNm0 <- names(PTMs_ratRf)
    for (ptm in EnrichedPTMs) { #ptm <- EnrichedPTMs[1]
      p <- Modifs$Mark[match(ptm, Modifs$`Full name`)]
      ppat <- paste0("\\(", p, "\\)|\\(", p, ",|,", p, "\\)|,", p, ",") # Pattern to catch all instances of the mod
      PTMs_FC_filt[[ptm]] <- c()
      PTMs_FC_Smpls[[ptm]] <- list()
      PTM_normalize[[ptm]] <- TRUE
      ptmpep <- pep[which(pep[[ptm]]),]
      a <- unlist(strsplit(gsub("\\)$", "", ptm), "\\("))
      if (length(a) > 1) {
        Ptm <- paste0(toupper(substr(a[1], 1, 1)), substr(a[1], 2, nchar(ptm)), "(", a[2], ")")
      } else {
        Ptm <- paste0(toupper(substr(a[1], 1, 1)), substr(a[1], 2, nchar(ptm)))
      }
      p <- Modifs$Mark[match(ptm, Modifs$`Full name`)]
      ptmsh <- substr(p, 1, 1)
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
      ptmpep[, c("Match(es)", paste0(Ptm, "-site(s)"))] <- as.data.frame(t(apply(temp[, kol], 1, function(x) {
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
            matches <- set_colnames(reshape2::melt(matches), c("Match", "Protein"))
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
      })))
      ptmpep[[paste0(Ptm, "-site")]] <- gsub(" .+", "", ptmpep[[paste0(Ptm, "-site(s)")]])
      ptmpep <- ptmpep[which(!is.na(ptmpep$`Match(es)`)),]
      ptmpep$tmp1 <- gsub("^_|_$", "", ptmpep$`Modified sequence`)
      ptmpep$tmp2 <- ptmpep[[paste0(Ptm, "-site(s)")]]
      nc <- nchar(ptmpep$tmp2)
      w <- which(nc > 25)
      ptmpep$tmp2[w] <- paste0(substr(ptmpep$tmp2[w], 1, 23), "...")
      ptmpep$Code <- apply(ptmpep[, paste0("tmp", 1:2)], 1, paste, collapse = "\n")
      ptmpep$tmp1 <- NULL
      ptmpep$tmp2 <- NULL
      ptmpep$Name <- ""
      w <- which(ptmpep$"Leading proteins" != "")
      ptmpep$Name[w] <- vapply(strsplit(gsub("[/,;].+$", "", ptmpep[w, paste0(Ptm, "-site(s)")]), " "), function(x) {
        paste0(x[[1]], " ", db$"Common Name"[match(x[[2]], db$"Protein ID")])
      }, "")
      ptmpep$Name[which(ptmpep$Name == "")] <- paste0("Unknown ", ptm, "-modified peptide #", 1:length(which(ptmpep$Name == "")))
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
      #
      # Optional: normalize to parent protein group
      if (PTM_normalize[[ptm]]) {
        # Step 1: normalize ratios:
        a <- grep(topattern(PTMs_ratRf[PTMs_ratNm0]), colnames(ptmpep), value = TRUE)
        a1 <- gsub(topattern(paste0(PTMs_ratRf[PTMs_ratNm0], " - ")), PG.rat.col, a)
        stopifnot(sum(!a1 %in% colnames(PG)) == 0)
        temp <- Isapply(strsplit(ptmpep$"Protein group IDs", ";"), function(x) { #x <- strsplit(ptmpep$"Protein group IDs", ";")[1]
          x <- unlist(x)
          y <- PG[match(x, PG$id), a1, drop = FALSE]
          if (length(x) > 1) { y <- apply(y, 2, function(x) { mean(is.all.good(x)) }) }
          return(unlist(y))
        })
        ptmpep[, paste0("ReNorm. ", a)] <- ptmpep[, a] - temp # It's log data so "-", not "/"
        PTMs_ratRf["Re-normalized"] <- paste0("ReNorm. ", PTMs_ratRf[PTMs_ratNm0])
        # Step 2: adjust expression values to reflect the re-normalized ratios:
        # Intensities are not logged (for peptides) so we need to delog ratios and work in multiplicative, not additive, mode 
        PTMs_intRf["Re-normalized"] <- "ReNorm. int."
        temp <- set_colnames(data.frame(matrix(rep(0, nrow(ptmpep)*length(Exp)),
                                               ncol = length(Exp))),
                             paste0(PTMs_intRf["Re-normalized"], " - ", Exp))
        for (grp in rat.grps) { #grp <- rat.grps[1]
          e <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
          smpls0 <- e$Experiment[which(e$Reference %in% c("TRUE", TRUE))]
          stopifnot(length(smpls0) == 1) # This is a script without replicates! Only one reference is allowed per sample group!
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
          totB <- apply(ptmpep[, c(cole0, cole1)], 1, sum, na.rm = TRUE)
          if (length(colr0)) {
            av <- apply(ptmpep[, cole0, drop = FALSE], 1, mean, na.rm = TRUE)
            temp[, c(cols0, cols1)] <- sweep(2^ptmpep[, c(colr0, colr1)], 1, av, "*")
            #tst <- log10(temp[, cols1]/temp[, cols0])
            #tst2 <- apply(tst, 2, function(x) { summary(is.all.good(x)) })
          } else {
            temp[, cols0] <- ptmpep[, cole0] # This stays the same as before
            temp[, cols1] <- ptmpep[, cole0]*(2^ptmpep[, colr1])
          }
          #tst <- log10(temp[, cols1]/temp[, cols0])
          #tst2 <- apply(tst, 2, function(x) { summary(is.all.good(x)) })
          # Note: the price of normalisation is that often there is no valid parent protein so a lot of NAs are introduced
          #
          totA <- apply(temp[, c(cols0, cols1)], 1, sum, na.rm = TRUE)
          temp[, c(cols0, cols1)] <- temp[, c(cols0, cols1)]*totB/totA # Re-apply correct scale
        }
        ptmpep[, colnames(temp)] <- temp
        #kol <- grep("ReNorm. log2", colnames(ptmpep), value = TRUE)
        #tst <- apply(ptmpep[, kol], 2, function(x) { summary(is.all.good(x)) })
      }
      #
      # Fold change filters:
      Int <- PTMs_intRf[length(PTMs_intRf)]
      Rat <- PTMs_ratRf[length(PTMs_ratRf)]
      for (grp in rat.grps) { #grp <- rat.grps[1]
        SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
        smpl0 <- SmplMp$Experiment[which(SmplMp$Reference)]
        stopifnot(length(smpl0) == 1)
        smpl1 <- SmplMp$Experiment[which(!SmplMp$Reference)]
        e0 <- ptmpep[[paste0(Int, " - ", smpl0)]]
        if (NegFilt) {
          nf <- SmplMp$Experiment[which(SmplMp$"Negative Filter")]
          if (length(nf)) {
            smpl1 <- smpl1[which(!smpl1 %in% nf)]
            nftst <- apply(ptmpep[, paste0(Int, " - ", nf), drop = FALSE], 1, function(x) { length(is.all.good(x)) }) > 0
          }
        }
        PTMs_FC_filt[[ptm]] <- append(PTMs_FC_filt[[ptm]], setNames(lapply(smpl1, function(x) {
          e1 <- ptmpep[[paste0(Int, " - ", x)]]
          r1 <- ptmpep[[paste0(Rat, " - ", x)]]
          if (RatiosThresh_2sided) { r1 <- abs(r1) }
          w <- which(((r1 >= RatiosThresh)&(is.all.good(e1, 2)))|(is.all.good(e1, 2)&(!is.all.good(e0, 2))))
          if ((NegFilt)&&(length(nf))) { w <- w[which(!nftst[w])] }
          return(w)
        }), smpl1))
        PTMs_FC_Smpls[[ptm]][[grp]] <- list(Numerator = smpl1, Denominator = smpl0)
      }
      #
      # "Regulated/Enriched" columns
      Int <- PTMs_intRf[length(PTMs_intRf)]
      Rat <- PTMs_ratRf[length(PTMs_ratRf)]
      for (grp in rat.grps) { #grp <- rat.grps[1]
        SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
        if (sum(c(TRUE, FALSE) %in% SmplMp$Reference) < 2) {
          stop("The reference column should include TRUE and FALSE values!")
        } else {
          w0 <- which(SmplMp$Reference)
          stopifnot(length(w0) == 1)
          w1 <- which(!SmplMp$Reference)
          xp0 <- SmplMp$Experiment[w0]
          for (x in w1) { #x <- w1[1]
            xp1 <- SmplMp$Experiment[x]
            wonly <- which(is.all.good(ptmpep[[paste0(Int, " - ", xp1)]], 2)&(!is.all.good(ptmpep[[paste0(Int, " - ", xp0)]], 2)))
            PG[[paste0("Regulated - ", xp1)]] <- ""
            if (RatiosThresh_2sided) {
              wup <- which(ptmpep[[paste0(Rat, " - ", xp1)]] >= RatiosThresh)
              wdwn <- which(ptmpep[[paste0(Rat, " - ", xp1)]] <= -RatiosThresh)
              ptmpep[wup, paste0("Regulated - ", xp1)] <- "up"
              ptmpep[wdwn, paste0("Regulated - ", xp1)] <- "down"
            } else {
              ptmpep[[paste0("Enriched - ", xp1)]] <- c("", "+")[(ptmpep[[paste0(Rat, " - ", xp1)]] >= RatiosThresh)+1]
            }
            ptmpep[wonly, paste0("Regulated - ", xp1)] <- "specific"
          }
        }
      }
      # Ratios distribution:
      dir <- paste0(wd, "/Workflow control")
      temp <- ptmpep[, c("Code", as.character(sapply(paste0(PTMs_ratRf, " - "), function(x) {
        grep(topattern(x), colnames(ptmpep), value = TRUE)
      })))]
      temp <- reshape::melt.data.frame(temp, id.vars = "Code")
      colnames(temp)[which(colnames(temp) == "value")] <- "log2(Ratio)"
      temp$Experiment <- gsub(paste0(".*", topattern(paste0(PTMs_ratRf["Original"], " - "), start = FALSE)), "", temp$variable)
      temp$Experiment <- factor(temp$Experiment, levels = SamplesMap$Experiment)
      temp$Type <- gsub(paste0(" ?", topattern(PTMs_ratRf["Original"], start = FALSE), ".*$"), "", temp$variable)
      temp$Type[which(temp$Type == "")] <- "Orig."
      temp$Type <- paste0("log2(", tolower(temp$Type), " ratio)")
      temp$Type <- factor(temp$Type, levels = paste0("log2(", c("orig", "norm", "imput", "renorm"), ". ratio)"))
      long.dat$ratios <- temp
      temp <- temp[which(is.all.good(temp$"log2(Ratio)", 2)),]
      ttl <- paste0("Ratios density plot - ", ptm, "-modified peptides")
      plot <- ggplot(temp) + geom_histogram(aes(x = `log2(Ratio)`, fill = Type), bins = 100) +
        geom_vline(xintercept = RatiosThresh, colour = "red") +
        theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.25) +
        ggtitle(ttl) + facet_grid(Type~Experiment) +
        xlab("log2(Ratio)")
      if (RatiosThresh_2sided) { plot <- plot + geom_vline(xintercept = -RatiosThresh, colour = "red") }
      poplot(plot)
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
      #
      # Gene-Set Enrichment Analysis (GSEA)
      if (runGSEA) {
        dataType <- "modPeptides"
        Src <- paste0(libPath, "/extdata/R scripts/Sources/GSEA.R")
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
  for (exp in Exp) { #exp <- Exp[1]
    PG[[paste0("Biot. peptides count - ", exp)]] <- 0
    e <- ev[which(ev$Experiment == exp),]
    g <- grep(topattern(Modifs$Mark[wbiot], start = FALSE), e$"Modified sequence")
    if (length(g)) {
      e <- e[g,]
      temp2 <- listMelt(strsplit(e$"Protein group IDs", ";"), e$"Modified sequence")
      temp2 <- aggregate(temp2$L1, list(temp2$value), function(x) { length(unique(x)) })
      w2 <- which(PG$id %in% temp2$Group.1)
      PG[w2, paste0("Biot. peptides count - ", exp)] <- temp2$x[match(PG$id[w2], temp2$Group.1)]
    }
  }
}

#### Code chunk - Proteomic ruler
if (protrul) {
  ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1]
  if (length(Exp) > 1) { ref <- c(ref, paste0("Mean ", gsub(" - $", "", ref))) }
  temp <- try(Prot.Ruler(PG, db, ref, NuclL = ProtRulNuclL), silent = TRUE)
  if ((!"try-error" %in% class(temp))&&((!"logical" %in% class(temp)))) {
    PG <- temp$Protein.groups
    db <- temp$Database
  } else { protrul <- FALSE }
}

Script <- readLines(ScriptPath)
gc()
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - Summary table and QC plots
source(parSrc, local = FALSE)
Exp_summary %<o% MQ.summary(wd = wd, ev = ev, pg = PG, mods = setNames(Modifs$Mark, Modifs$"Full name"),
                            raw.files = rawFiles, sc = max(c(20, round(length(rawFiles2)/length(Exp)))),
                            save = c("jpeg", "pdf"), cl = parClust,
                            MQtxt = inDirs[which(SearchSoft == "MAXQUANT")])
write.csv(Exp_summary, paste0(wd, "/Workflow control/Summary.csv"), row.names = FALSE)
#Exp_summary <- read.csv(paste0(wd, "/Workflow control/Summary.csv"), check.names = FALSE)

# Plot of contamination levels per sample
tmp <- ev[, c("Proteins", "Intensity", "Experiment")]
tmp2 <- listMelt(strsplit(tmp$Proteins, ";"), 1:nrow(tmp), c("Protein", "Row"))
m <- match(gsub("^CON__", "", tmp2$Protein), gsub("^CON__", "", db$`Protein ID`))
w <- which(!is.na(m))
tmp2 <- tmp2[w,]; m <- m[w]
# For our purpose here we must match contaminant proteins.
tmp2$Cont <- db$`Potential contaminant`[m]
if (tstOrg) {
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
tmp <- aggregate(tmp$Intensity, list(tmp$Experiment, tmp$Organism), sum)
colnames(tmp) <- c("Experiment", "Organism", "Total intensity")
tmp$Experiment <- factor(tmp$Experiment, levels = SamplesMap$Experiment)
ttl <- "Contributions to TIC"
plot <- ggplot(tmp) + geom_bar(stat = "identity", aes(x = Experiment, y = `Total intensity`, fill = Organism)) +
  theme_bw() + scale_fill_viridis(discrete = TRUE, begin = 0.8, end = 0.2) +
  ggtitle(ttl, subtitle = "Summed TIC for each class of identified peptides")
poplot(plot)
ggsave(paste0(wd, "/Summary plots/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
ggsave(paste0(wd, "/Summary plots/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")

# Test for amino acid biases:
Src <- paste0(libPath, "/extdata/R scripts/Sources/AA_biases_test.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

Script <- readLines(ScriptPath)
gc()
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - Coverage maps for proteins of interest
Src <- paste0(libPath, "/extdata/R scripts/Sources/noRep_protPlots.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

Script <- readLines(ScriptPath)
gc()
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

### Check that CytoScape is installed and can run, then launch it.
#CytoScape <- TRUE #You may need to reset this
Src <- paste0(libPath, "/extdata/R scripts/Sources/Cytoscape_init.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Initialize ClueGO
if (enrichGO||globalGO) {
  Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_init.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
}
  
#### Code chunk - Gene Ontology terms enrichment analysis
source(parSrc, local = FALSE)
create_plotly %<o% TRUE
if (globalGO) {
  # Global dataset GO enrichment - expression per sample vs total proteome
  ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1]
  xprsFilt <- setNames(lapply(Exp, function(e) { which(PG[[paste0(ref, e)]] > 0) }), Exp)
  dir <- paste0(wd, "/GO enrichment analysis/Sample vs total proteome")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  setwd(dir)
  #temp <- PG
  #temp$`Leading protein IDs` <- gsub("CON__[^;]+", "", temp$`Leading protein IDs`) # Remove contaminants
  #
  Mode <- "dataset"
  dataType <- "PG"
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
  suppressWarnings(try(rm(list = allArgs), silent = TRUE))
  #
  temp <- goRES
  #
  setwd(wd)
  if (class(temp) != "logical") {
    temp <- temp$GO_terms
    temp$"Protein group IDs" <- vapply(temp$`Protein table row(s)`, function(x) {
      paste(PG$id[as.numeric(x)], collapse = ";")
    }, "")
    temp$Mapping <- NULL
    if ("Offspring" %in% colnames(temp)) {
      temp$Offspring <- vapply(temp$Offspring, paste, "", collapse = ";")
    }
    temp$`Protein table row(s)` <- NULL
    gn <- grep("^Genes", colnames(temp), value = TRUE)
    pr <- grep("^Proteins", colnames(temp), value = TRUE)
    pg <- grep("^Protein group IDs", colnames(temp), value = TRUE)
    kn <- grep("^Count", colnames(temp), value = TRUE)
    pv <- grep("^Pvalue", colnames(temp), value = TRUE)
    zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
    lf <- grep("^logFC", colnames(temp), value = TRUE)
    si <- grep("^Significance", colnames(temp), value = TRUE)
    kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, pr))]
    temp <- temp[, c(kl, si, gn, pg, pr, kn, pv, zs, lf)]
    w <- apply(temp[, pv, drop = FALSE], 1, function(x) { sum(!is.na(x)) }) > 0
    temp <- temp[w,]
    tst <- apply(temp[, kn, drop = FALSE], 1, function(x) { sum(x[which(!is.na(x))]) })
    temp <- temp[order(tst, decreasing = TRUE),]
    temp <- temp[order(temp$Ontology, decreasing = FALSE),]
    ExcelMax <- 32767
    w <- which(vapply(colnames(temp), function(x) { "character" %in% class(temp[[x]]) }, TRUE))
    if (length(w)) {
      for (i in w) { #i <- w[1]
        w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
        if (length(w1)) {
          temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1, ExcelMax-3), "...")
        }
      }
    }
    wb <- openxlsx2::wb_workbook()
    myFont <- openxlsx2::create_font(sz = "12", color = openxlsx2::wb_color("#FFFFFF"))
    wb$styles_mgr$add(myFont, "Header_font")
    HdrStl <- openxlsx2::create_cell_style(font_id = 1,
                                           num_fmt_id = "General",
                                           horizontal = "justify",
                                           vertical = "top",
                                           text_rotation = 60,
                                           wrap_text = TRUE)
    wb$styles_mgr$add(HdrStl, "Header_style")
    nTabs <- 0
    for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1]
      w <- which(temp$Ontology == Ontologies[ont])
      if (length(w)) {
        sheetNm <- ont
        if (sheetNm %in% wb_get_sheet_names(wb)) {
          wb <- openxlsx2::wb_remove_worksheet(wb, sheetNm)
          nTabs <- nTabs + 1
        }
        wb <- openxlsx2::wb_add_worksheet(wb, sheetNm)
        nTabs <- nTabs + 1
        wb <- openxlsx2::wb_add_data_table(wb, sheetNm, temp[w,], openxlsx2::wb_dims(2, 2))
        wb <- openxlsx2::wb_set_row_heights(wb, sheetNm, 2, 120)
        wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1, 1)
        wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1+1:ncol(temp), 12)
        wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1+which(colnames(temp) == "Term"), 45)
        wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1+which(colnames(temp) %in% c(gn, pr)), 20)
        dms <- openxlsx2::wb_dims(2, 1+1:ncol(temp))
        wb <- openxlsx2::wb_set_cell_style(wb, sheetNm, dms, wb$styles_mgr$get_xf_id("Header_style"))
        dms <- openxlsx2::wb_dims(1+1:length(w), 1+which(colnames(temp) %in% kn))
        wb <- openxlsx2::wb_add_numfmt(wb, sheetNm, dms, "0")
        dms <- openxlsx2::wb_dims(1+1:length(w), 1+which(colnames(temp) %in% c(zs, lf, pv)))
        wb <- openxlsx2::wb_add_numfmt(wb, sheetNm, dms, "0.000")
      }
    }
    if (nTabs) {
      goXLfl <- paste0(dir, "/GO terms - Samples vs Proteome.xlsx")
      openxlsx2::wb_save(wb, goXLfl)
      #openxlsx2::xl_open(goXLfl)
    }
  }
  if (enrichGO) {
    for (grp in rat.grps) { #grp <- rat.grps[1]
      wb <- openxlsx2::wb_workbook()
      myFont <- openxlsx2::create_font(sz = "12", color = openxlsx2::wb_color("#FFFFFF"))
      wb$styles_mgr$add(myFont, "Header_font")
      HdrStl <- openxlsx2::create_cell_style(font_id = 1,
                                             num_fmt_id = "General",
                                             horizontal = "justify",
                                             vertical = "top",
                                             text_rotation = 60,
                                             wrap_text = TRUE)
      wb$styles_mgr$add(HdrStl, "Header_style")
      nTabs <- 0
      SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
      tmp <- paste0("Comparisons to ", SmplMp$Experiment[which(SmplMp$Reference)])
      if (length(rat.grps) > 1) { tmp <- paste0("Group", grp, " (", tmp, ")") }
      dir <- paste0(wd, "/GO enrichment analysis/", tmp)
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      setwd(dir)
      db_obs %<o% db[which(db$"Protein ID" %in% unlist(strsplit(PG$"Leading protein IDs", ";"))),]
      fcFilt <- FC_filt[which(names(FC_filt) %in% SmplMp$Experiment)]
      nms <- names(fcFilt)
      filt2 <- setNames(lapply(nms, function(x) {
        unique(unlist(xprsFilt[SmplMp$Experiment]))
      }), nms)
      #
      Mode <- "regulated"
      dataType <- "PG"
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
      suppressWarnings(try(rm(list = allArgs), silent = TRUE))
      #
      temp <- goRES
      #
      setwd(wd)
      if (class(temp) != "logical") {
        temp <- temp$GO_terms
        temp$"Protein group IDs" <- vapply(temp$`Protein table row(s)`, function(x) {
          paste(PG$id[as.numeric(x)], collapse = ";")
        }, "")
        temp$Mapping <- NULL
        if ("Offspring" %in% colnames(temp)) {
          temp$Offspring <- vapply(temp$Offspring, paste, "", collapse = ";")
        }
        temp$`Protein table row(s)` <- NULL
        gn <- grep("^Genes", colnames(temp), value = TRUE)
        pr <- grep("^Proteins", colnames(temp), value = TRUE)
        pg <- grep("^Protein group IDs", colnames(temp), value = TRUE)
        kn <- grep("^Count", colnames(temp), value = TRUE)
        pv <- grep("^Pvalue", colnames(temp), value = TRUE)
        zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
        lf <- grep("^logFC", colnames(temp), value = TRUE)
        si <- grep("^Significance", colnames(temp), value = TRUE)
        kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, pr))]
        temp <- temp[, c(kl, si, gn, pg, pr, kn, pv, zs, lf)]
        w <- apply(temp[, pv, drop = FALSE], 1, function(x) { sum(!is.na(x)) }) > 0
        temp <- temp[w,]
        tst <- apply(temp[, kn, drop = FALSE], 1, function(x) { sum(x[which(!is.na(x))]) })
        temp <- temp[order(tst, decreasing = TRUE),]
        temp <- temp[order(temp$Ontology, decreasing = FALSE),]
        ExcelMax <- 32767
        w <- which(vapply(colnames(temp), function(x) { "character" %in% class(temp[[x]]) }, TRUE))
        if (length(w)) {
          for (i in w) { #i <- w[1]
            w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
            if (length(w1)) {
              temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1, ExcelMax-3), "...")
            }
          }
        }
        for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1]
          w <- which(temp$Ontology == Ontologies[ont])
          if (length(w)) {
            sheetNm <- ont
            if (sheetNm %in% wb_get_sheet_names(wb)) {
              wb <- openxlsx2::wb_remove_worksheet(wb, sheetNm)
              nTabs <- nTabs-1
            }
            wb <- openxlsx2::wb_add_worksheet(wb, sheetNm)
            nTabs <- nTabs + 1
            wb <- openxlsx2::wb_add_data_table(wb, sheetNm, temp[w,], openxlsx2::wb_dims(2, 2))
            wb <- openxlsx2::wb_set_row_heights(wb, sheetNm, 2, 120)
            wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1, 1)
            wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1+1:ncol(temp), 12)
            wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1+which(colnames(temp) == "Term"), 45)
            wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1+which(colnames(temp) %in% c(gn, pr)), 20)
            dms <- openxlsx2::wb_dims(2, 1+1:ncol(temp))
            wb <- openxlsx2::wb_set_cell_style(wb, sheetNm, dms, wb$styles_mgr$get_xf_id("Header_style"))
            dms <- openxlsx2::wb_dims(1+1:length(w), 1+which(colnames(temp) %in% kn))
            wb <- openxlsx2::wb_add_numfmt(wb, sheetNm, dms, "0")
            dms <- openxlsx2::wb_dims(1+1:length(w), 1+which(colnames(temp) %in% c(zs, lf, pv)))
            wb <- openxlsx2::wb_add_numfmt(wb, sheetNm, dms, "0.000")
          }
        }
        if (nTabs) {
          goXLfl <- paste0(dir, "/GO terms - Comp group ", grp, ".xlsx")
          openxlsx2::wb_save(wb, goXLfl)
          #openxlsx2::xl_open(goXLfl)
        }
      }
    }
    if (PTMriched) {
      for (ptm in EnrichedPTMs) {
        ptmpep <- PTMs_pep[[ptm]]
        for (grp in rat.grps) { #grp <- rat.grps[1]
          SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
          tmp <- paste0("Comparisons to ", SmplMp$Experiment[which(SmplMp$Reference)])
          if (length(rat.grps) > 1) { tmp <- paste0("Group", grp, " (", tmp, ")") }
          nms <- names(PTMs_FC_filt[[ptm]])
          filt3 <- setNames(lapply(nms, function(x) {
            # Filter: any expressed in the same group, i.e. at least with 1 non null value in the group
            which(rowSums(ptmpep[, paste0(PTMs_intRf[length(PTMs_intRf) # Not log!
                                                     ], " - ", SmplMp$Experiment)], na.rm = TRUE) > 0)
          }), nms)
          dir <- paste0(wd, "/GO enrichment analysis/", tmp)
          if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
          setwd(dir)
          #
          Mode <- "regulated"
          dataType <- "modPeptides"
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
          suppressWarnings(try(rm(list = allArgs), silent = TRUE))
          #
          PTMs_GO_Plots[[Ptm]][[tstbee]] <- goRES
          #
          setwd(wd)
          #
          setwd(wd)
          if (class(temp) != "logical") {
            temp <- temp$GO_terms
            temp$"Protein group IDs" <- vapply(temp$`Protein table row(s)`, function(x) {
              paste(PG$id[as.numeric(x)], collapse = ";")
            }, "")
            temp$Mapping <- NULL
            if ("Offspring" %in% colnames(temp)) {
              temp$Offspring <- vapply(temp$Offspring, paste, "", collapse = ";")
            }
            temp$`Protein table row(s)` <- NULL
            gn <- grep("^Genes", colnames(temp), value = TRUE)
            pr <- grep("^Proteins", colnames(temp), value = TRUE)
            pg <- grep("^Protein group IDs", colnames(temp), value = TRUE)
            kn <- grep("^Count", colnames(temp), value = TRUE)
            pv <- grep("^Pvalue", colnames(temp), value = TRUE)
            zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
            lf <- grep("^logFC", colnames(temp), value = TRUE)
            si <- grep("^Significance", colnames(temp), value = TRUE)
            kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, pr))]
            temp <- temp[, c(kl, si, gn, pg, pr, kn, pv, zs, lf)]
            w <- apply(temp[, pv, drop = FALSE], 1, function(x) { sum(!is.na(x)) }) > 0
            temp <- temp[w,]
            tst <- apply(temp[, kn, drop = FALSE], 1, function(x) { sum(x[which(!is.na(x))]) })
            temp <- temp[order(tst, decreasing = TRUE),]
            temp <- temp[order(temp$Ontology, decreasing = FALSE),]
            ExcelMax <- 32767
            w <- which(vapply(colnames(temp), function(x) { "character" %in% class(temp[[x]]) }, TRUE))
            if (length(w)) {
              for (i in w) { #i <- w[1]
                w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
                if (length(w1)) {
                  temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1, ExcelMax-3), "...")
                }
              }
            }
            for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1]
              w <- which(temp$Ontology == Ontologies[ont])
              if (length(w)) {
                sheetNm <- ont
                if (sheetNm %in% wb_get_sheet_names(wb)) {
                  wb <- openxlsx2::wb_remove_worksheet(wb, sheetNm)
                  nTabs <- nTabs-1
                }
                wb <- openxlsx2::wb_add_worksheet(wb, sheetNm)
                nTabs <- nTabs + 1
                wb <- openxlsx2::wb_add_data_table(wb, sheetNm, temp[w,], openxlsx2::wb_dims(2, 2))
                wb <- openxlsx2::wb_set_row_heights(wb, sheetNm, 2, 120)
                wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1, 1)
                wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1+1:ncol(temp), 12)
                wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1+which(colnames(temp) == "Term"), 45)
                wb <- openxlsx2::wb_set_col_widths(wb, sheetNm, 1+which(colnames(temp) %in% c(gn, pr)), 20)
                dms <- openxlsx2::wb_dims(2, 1+1:ncol(temp))
                wb <- openxlsx2::wb_set_cell_style(wb, sheetNm, dms, wb$styles_mgr$get_xf_id("Header_style"))
                dms <- openxlsx2::wb_dims(1+1:length(w), 1+which(colnames(temp) %in% kn))
                wb <- openxlsx2::wb_add_numfmt(wb, sheetNm, dms, "0")
                dms <- openxlsx2::wb_dims(1+1:length(w), 1+which(colnames(temp) %in% c(zs, lf, pv)))
                wb <- openxlsx2::wb_add_numfmt(wb, sheetNm, dms, "0.000")
              }
            }
            if (nTabs) {
              goXLfl <- paste0(dir, "/GO terms - ", ptm, " - Comp group ", grp, ".xlsx")
              openxlsx2::wb_save(wb, goXLfl)
              #openxlsx2::xl_open(goXLfl)
            }
          }
        }
      }
    }
  }
}

#### Code chunk - Venn diagrams
setwd(wd)
HdrStlVenn <- createStyle(textDecoration = "bold", halign = "left", valign = "bottom", wrapText = TRUE,
                          numFmt = "TEXT", fontSize = 12, textRotation = 60)
wb <- createWorkbook()
wbKount <- 0
VennMx <- 7
if (Venn_Obs) {
  dir <- paste0(wd, "/Venn diagrams")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  ref <- PG.int.cols["Original"]
  Grps <- unique(SamplesMap$`Ratios group`)
  for (grp in Grps) {
    sm <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
    Xp <- sm$Experiment
    test <- setNames(lapply(Xp, function(exp) {
      x <- PG[[paste0(ref, exp)]]
      which(is.all.good(x, 2))
    }), Xp)
    w <- which(vapply(test, length, 1) > 0)
    VennExp <- Xp[w]
    OK <- length(w) > 1
    if (OK) {
      if (length(w) > VennMx) {
        msg <- paste0("Too many samples, select at least 2 and up to ", VennMx,
                      " to include in the Venn diagram (comma-separated):")
        if (length(Grps) > 1) { msg <- paste0("Ratios group ", grp, ": ", msg) }
        tst <- !sm$Reference[match(VennExp, sm$Experiment)]
        tmp <- setNames(vapply(VennExp, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") }, ""), VennExp)
        VennExp <- names(tmp)[match(dlg_list(tmp, tmp[which(tst)[1:min(c(sum(tst, VennMx)))]], TRUE, msg)$res, tmp)]
        if (length(VennExp) < 2) {
          msg <- "Skipping per-sample observations Venn diagrams"
          if (length(Grps) > 1) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected at least 2 samples!")
          warning(msg)
          OK <- FALSE
        }
        if (length(VennExp) > VennMx) {
          msg <- "Skipping per-sample observations Venn diagrams"
          if (length(Grps) > 1) { msg <- paste0(msg, " for ratios group ", grp) }
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
      if (length(Grps) > 1) {
        ttl <- paste0(ttl, " - group ", grp)
        SheetNm <- paste0(SheetNm, "_grp", grp)
      }
      test <- test[VennExp]
      AnalysisParam$"Venn diagram (LFC) - samples" <- list(VennExp)
      plot <- venn(test, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
      plot <- plot + ggtitle("Venn diagram", subtitle = "Global, LFQ") +
        theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 10))
      poplot(plot)
      ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 150)
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150)
      #system(paste0("open \"", dir, "/", ttl, ".jpg", "\""))
      wbKount <- wbKount+1
      if (SheetNm %in% names(wb)) { removeWorksheet(wb, SheetNm) } 
      addWorksheet(wb, SheetNm)
      writeData(wb, SheetNm, PG[, c("id", "Leading protein IDs", "Genes")], 1, 1)
      l <- length(test)
      tmp <- sapply(names(test), function(smpl) {
        res <- rep("", nrow(PG))
        res[test[[smpl]]] <- "+"
        return(res)
      })
      writeData(wb, SheetNm, tmp, 4, 1)
      setRowHeights(wb, SheetNm, 1, 120)
      addStyle(wb, SheetNm, HdrStlVenn, 1, 1:(l+3))
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
    tst <- vapply(fc_filt, length, 1)
    w <- which(tst == 0)
    if (0 %in% tst) {
      l <- length(w)
      tmp <- names(fc_filt)[which(tst == 0)]
      if (l > 1) { tmp <- paste0(paste(tmp[1:(l-1)], collapse = ", "), " and ", tmp[l]) }
      warning(paste0("No filtered proteins for samples ", tmp, ", they will be skipped!"))
    }
    fc_filt <- fc_filt[which(tst > 0)]
    VennExp <- names(fc_filt)
    OK <- length(VennExp) > 1
    if (OK) {
      if (length(fc_filt) > VennMx) {
        msg <- paste0("Too many samples, select at least 2 and up to ", VennMx,
                      " to include in the Venn diagram (comma-separated):")
        if (length(Grps) > 1) { msg <- paste0("Ratios group ", grp, ": ", msg) }
        VennExp <- dlg_list(VennExp, VennExp[1:VennMx], TRUE, msg)$res
        if (length(VennExp) < 2) {
          msg <- "Skipping ratios Venn diagrams"
          if (length(Grps) > 1) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected at least 2 samples!")
          warning(msg)
          OK <- FALSE
        }
        if (length(VennExp) > VennMx) {
          msg <- "Skipping ratios Venn diagrams"
          if (length(Grps) > 1) { msg <- paste0(msg, " for ratios group ", grp) }
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
          updowntst <- setNames(lapply(VennExp, function(x) {
            #x <- VennExp[1]
            rs <- sign(PG[fc_filt[[x]], paste0(PG.rat.cols["Original"], x)])
            w <- which(is.na(rs))
            rs[w] <- c(1, -1)[is.na(PG[fc_filt[[x]][w], paste0(PG.int.cols["Original"], x)])+1]
            return(rs)
          }), VennExp)
        }
        setwd(wd); suppressWarnings(dir.create("Venn diagrams"))
        for (vt in VennTypes) { #vt <- ""
          ttl <- gsub("_\\(\\)$", "", paste0("Ratios_Venn_diagram_-_global", "_(", vt, ")"))
          SheetNm <- paste0(c("Up/down", "Up", "Down")[match(vt, VennTypes)], "-reg. PGs")
          if (length(Grps) > 1) {
            msg <- paste0(msg, " for ratios group ", grp)
            ttl <- paste0(ttl, " - group ", grp)
            SheetNm <- paste0(SheetNm, "_grp", grp)
          }
          sbttl <- rat.col
          flt <- fc_filt
          if (vt %in% c("up", "down")) {
            sbttl <- paste0(sbttl, ", ", vt)
            flt <- setNames(lapply(VennExp, function(x) {
              flt[[x]][which(updowntst[[x]] == c(1, -1)[match(vt, c("up", "down"))])]
            }), VennExp)
          }
          w <- which(vapply(flt, length, 1) > 0)
          if (length(w) > 1) {
            flt <- flt[w]
            plot <- venn(flt, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
            plot <- plot + ggtitle("Venn diagram", subtitle = sbttl) +
              theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 10))
            poplot(plot)
            ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 150)
            ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150)
            #system(paste0("open \"", dir, "/", ttl, ".jpg", "\""))
            wbKount <- wbKount+1
            if (SheetNm %in% names(wb)) { removeWorksheet(wb, SheetNm) } 
            addWorksheet(wb, SheetNm)
            writeData(wb, SheetNm, PG[, c("id", "Leading protein IDs", "Genes")], 1, 1)
            l <- length(flt)
            tmp <- sapply(names(flt), function(smpl) {
              res <- rep("", nrow(PG))
              res[flt[[smpl]]] <- "+"
              return(res)
            })
            writeData(wb, SheetNm, tmp, 4, 1)
            setRowHeights(wb, SheetNm, 1, 120)
            addStyle(wb, SheetNm, HdrStlVenn, 1, 1:(l+3))
          } else { message(gsub(" \\(\\)$", "", paste0("No overlaps to plot (", vt, ")"))) }
        }
      }
    } else {
      msg <- "Skipping ratios Venn diagrams"
      if (length(Grps) > 1) { msg <- paste0(msg, " for ratios group ", grp) }
      msg <- paste0(msg, ": not enough valid samples!")
      warning(msg)
    }
  }
}
if (wbKount) { saveWorkbook(wb, paste0(wd, "/Venn diagrams/Venn diagrams.xlsx"), overwrite = TRUE) }
setwd(wd)

#### Code chunk - Dimensionality reduction plots
## Process data
if (length(Exp) > 2) {
  setwd(wd); suppressWarnings(dir.create("PCA plots"))
  temp <- PG[, grep(topattern(PG.int.col), colnames(PG), value = TRUE)]
  colnames(temp) <- gsub(topattern(PG.int.col), "", colnames(temp))
  ## Impute data
  ## Here we have no way to decide between MAR/MCAR/MNAR,
  ## so we will instead replace every missing value with a random value drawn from a gaussian distribution of reduced m and sd
  m <- median(is.all.good(unlist(temp)))
  sd <- sd(is.all.good(unlist(temp)))
  for (i in colnames(temp)) {
    w <- which(!is.all.good(temp[[i]], 2))
    temp2 <- rnorm(length(w), m-3, sd/2)
    temp[w, i] <- temp2
  }
  nrm <- PG[[paste0("Mean ", gsub(" - $", "", PG.int.col))]]
  w <- which(is.all.good(nrm, 2))
  temp <- sweep(temp[w,], 1, nrm[w], "-")
  rownames(temp) <- PG$"Leading protein IDs"[w]
  temp <- temp + rnorm(nrow(temp)*ncol(temp), 0, 10^-24) # Add small random value in case we (very rarely) get non-unique values per row
  ## 1/ Samples level:
  pc <- prcomp(t(temp), scale. = TRUE)
  scores <- as.data.frame(pc$x)
  pv <- round(100*(pc$sdev)^2 / sum(pc$sdev^2), 0)
  pv <- pv[which(pv > 0)]
  pv <- paste0("Components: ", paste(vapply(1:length(pv), function(x) {
    paste0("PC", x, ": ", pv[x], "%")
  }, ""), collapse = ", "))
  scores$Sample <- rownames(scores)
  rownames(scores) <- NULL
  tst <- apply(scores[, c("PC1", "PC2")], 2, function(x) { length(unique(x)) })
  if (min(tst) > 1) {
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
      if ("try-error" %in% class(tst)) { tst <- try(saveWidget(plot_lyPCA, paste0(wd, "/PCA plots/", ttl, ".html")), silent = TRUE) }
      if (!"try-error" %in% class(tst)) { system(paste0("open \"", wd, "/PCA plots/", ttl, ".html")) }
    } else { poplot(plot, width = 18) }
    ggsave(paste0(wd, "/PCA plots/", ttl, ".jpeg"), plot, dpi = 150)
    ggsave(paste0(wd, "/PCA plots/", ttl, ".pdf"), plot, dpi = 150)
    ## 2/ Protein groups level:
    pc <- prcomp(temp, scale. = TRUE)
    scores <- as.data.frame(pc$x)
    pv <- round(100*(pc$sdev)^2 / sum(pc$sdev^2), 0)
    pv <- pv[which(pv > 0)]
    pv <- paste0("Components: ", paste(vapply(1:length(pv), function(x) {
      paste0("PC", x, ": ", pv[x], "%")
    }, ""), collapse = ", "))
    scores$"Leading protein IDs" <- rownames(scores)
    rownames(scores) <- NULL
    scores[, c("Protein IDs", "Protein group")] <- PG[match(scores$"Leading protein IDs", PG$"Leading protein IDs"),
                                                      c("Protein IDs", "Label")]
    scores$Alpha <- (scores$PC1^2 + scores$PC2^2)
    scores$Direction <- apply(temp[w,], 1, function(x) {
      wh <- which(is.all.good(10^x, 2))
      return(weighted.mean(c(1:length(colnames(temp)))[wh], 10^x[wh]))
    })
    scores$Class <- ""
    breaks <- 1:length(Exp)
    labels <- Exp
    if (prot.list.Cond) {
      g1 <- grsep2(prot.list, scores$"Protein IDs")
      if (length(g1)) {
        g2 <- grsep2(prot.list, scores$"Protein IDs", invert = TRUE)
        scores2 <- scores[g1,]
        scores <- scores[g2,]
        scores2$Alpha <- 1
        scores2$Class <- length(Exp)+1
      }
    }
    ttl <- "PCA plot - Protein groups (PG-level)"
    plot <- ggplot(scores) + geom_point(aes(x = PC1, y = PC2, alpha = Alpha, colour = Direction, text = `Protein group`), shape = 1) +
      coord_fixed() + ggtitle(ttl, subtitle = pv) +
      theme_bw() + guides(alpha = "none", shape = "none") +
      scale_color_gradientn(colors = hcl.colors(length(Exp)), breaks = breaks, labels = labels, guide = "legend")
    if ((prot.list.Cond)&&(length(g1)))
      plot <- plot + geom_point(data = scores2, colour = "red", shape = 2, aes(x = PC1, y = PC2, text = `Protein group`))
  }
  ggsave(paste0(wd, "/PCA plots/", ttl, ".jpeg"), plot, dpi = 150)
  ggsave(paste0(wd, "/PCA plots/", ttl, ".pdf"), plot, dpi = 150)
  if ("PC3" %in% colnames(scores)) {
    plot_lyPCAProt <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Direction, text = ~`Protein group`,
                              type = "scatter3d", mode = "markers", showlegend = FALSE, marker = list(size = 1))
    plot_lyPCAProt <- layout(plot_lyPCAProt, title = ttl)
    if ((prot.list.Cond)&&(length(g1))) {
      scores3 <- scores2
      scores3$"Protein group" <- gsub(" - | ?, ?", "<br>", scores3$"Protein group")
      plot_lyPCAProt <- add_trace(plot_lyPCAProt, data = scores3, x = ~PC1, y = ~PC2, z = ~PC3,
                                  type = "scatter3d", mode = "markers+text", color = I("red"), marker = list(size = 5),
                                  showlegend = FALSE, textposition = "bottom right")
    }
  } else { plot_lyPCAProt <- ggplotly(plot, tooltip = "text") }
  tst <- try(saveWidget(partial_bundle(plot_lyPCAProt), paste0(wd, "/PCA plots/", ttl, ".html")), silent = TRUE)
  if ("try-error" %in% class(tst)) { tst <- try(saveWidget(plot_lyPCAProt, paste0(wd, "/PCA plots/", ttl, ".html")), silent = TRUE) }
  if (!"try-error" %in% class(tst)) { system(paste0("open \"", wd, "/PCA plots/", ttl, ".html")) }
  #plot <- plot + geom_text_repel(data = scores, aes(x = PC1, y = PC2, label = `Protein group`, alpha = Alpha), size = 2.5)
  if ((prot.list.Cond)&&(length(g1))) {
    plot <- plot + geom_text_repel(data = scores2, colour = "red", size = 2.5,
                                   aes(x = PC1, y = PC2, label = `Protein group`))
  }
  #poplot(plot, width = 18)
  ggsave(paste0(wd, "/PCA plots/", ttl, " (labels).jpeg"), plot, dpi = 150)
  ggsave(paste0(wd, "/PCA plots/", ttl, " (labels).pdf"), plot)
} else { warning("No PCA plots drawn: samples are too similar!") }

#### Code chunk - Heatmaps with clustering at samples and protein groups level, highlighting proteins of interest
Src <- paste0(libPath, "/extdata/R scripts/Sources/cluster_Heatmap.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

Script <- readLines(ScriptPath)
gc()
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - Protein group profile plots and ranked abundance plots
if (runRankAbundPlots|runProfPlots) {
  source(parSrc, local = FALSE)
  require(RColorBrewer)
  require(colorspace)
  QuantTypes %<o% c("LFQ", "Coverage", "Spectra")
  tstOrg2 %<o% c()
  if (tstOrg) {
    PG$temp <- PG[[pgOrgKol]]
    PG$temp[which(PG$`Potential contaminant` == "+")] <- "Contaminant"
    tstOrg2 <- unique(PG$temp)
    tstOrg2 <- tstOrg2[which(tstOrg2 != "Contaminant")]
  }
  abbrFun %<o% function(x) { #x <- tst[[1]]
    g1 <- grep("[A-Z]", x)
    g2 <- grep("[a-z]", x)
    w1 <- g2[which((g2-1) %in% g1)]
    x[w1] <- "."
    gsub("\\.[^ ]+", ".", paste(x, collapse = ""))
  }
  tstOrg3 %<o% abbrFun(tstOrg2)
  for (QuantType in QuantTypes) { #QuantType <- "Coverage" #QuantType <- QuantTypes[1]
    myColors <- setNames("black", "-")
    myColors2 <- setNames(c("lightgrey", "brown"), c("-", "+"))
    if (length(tstOrg2)) {
      myColorsB <- myColors <- setNames(colorRampPalette(c("blue", "green"))(length(tstOrg2)), tstOrg2)
      w <- which(names(myColorsB) != "In list")
      names(myColorsB)[w] <- gsub("[a-z]+ ", ". ", names(myColorsB)[w])
      names(myColorsB) <- gsub(";", "\n", names(myColorsB))
      myColorsB[["potential contaminant"]] <- myColors[["potential contaminant"]] <- "grey"
    }
    myColorsB[[c("+", "In list")[(length(tstOrg2) > 0)+1]]] <- myColors[[c("+", "In list")[(length(tstOrg2) > 0)+1]]] <- "red"
    colScale <- scale_colour_manual(name = "Category", values = myColors)
    colScaleB <- scale_colour_manual(name = "Category", values = myColorsB)
    fillScale <- scale_fill_manual(name = "Category", values = myColors)
    colScale2 <- scale_colour_manual(name = "In list", values = myColors2)
    fillScale2 <- scale_fill_manual(name = "In list", values = myColors2)
    MainDir <- paste0(wd, "/Ranked abundance")
    SubDir <- paste0(MainDir, "/", QuantType)
    if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
    kolnm <- c("log10 LFQ", "Coverage [%]", "Spectral count")[match(QuantType, QuantTypes)]
    ref <- c(rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1],
             "Sequence coverage [%] - ",
             "Spectral count - ")[match(QuantType, QuantTypes)]
    g <- paste0(ref, Exp)
    Wh <- which(g %in% colnames(PG))
    g <- g[Wh]
    if (length(g)) {
      varkol <- c("Leading protein IDs", "Protein IDs", "Common Name (short)", "id", "Label")
      if (prot.list.Cond) { varkol <- c(varkol, "In list") }
      temp <- PG[, c(varkol, g)]
      if (prot.list.Cond) { temp$`In list`[which(temp$`In list` == "")] <- "-" }
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
      temp <- reshape::melt.data.frame(temp, id.vars = varkol)
      colnames(temp)[1:5] <- c(gsub(" ", "_", varkol[1:3]), "id", "Protein_Group")
      colnames(temp)[which(colnames(temp) == "variable")] <- "Sample"
      colnames(temp)[which(colnames(temp) == "value")] <- "Y"
      temp <- temp[which(is.all.good(temp$Y, 2)),]
      if (QuantType %in% c("Coverage", "Spectra")) { temp <- temp[which(temp$Y > 0),] }
      if (nrow(temp)) {
        if (length(tstOrg2)) { temp$Category <- PG$temp[match(temp$id, PG$id)] } else {
          temp$Category <- "-"
        }
        if (prot.list.Cond) {
          temp$Category[which(temp$"In list" == "+")] <- c("+", "In list")[(length(tstOrg2) > 0)+1]
        }
        lev <- c("In list", "+", tstOrg2, "-", "Contaminant")
        lev <- lev[which(lev %in% unique(temp$Category))]
        temp$Category <- factor(temp$Category, levels = lev)
        if ((QuantType == "LFQ")&&(AnalysisParam$Type == "Band ID")) {
          temp$Y <- 10^temp$Y
          kolnm <- "LFQ"
        }
        temp$Value <- paste0(kolnm, ": ", temp$Y)
        if (GO_filt) {
          for (goID in GO_filter) { #goID <- GO_filter[1]
            # Get children terms
            gofilter <- unique(unlist(c(goID,
                                        GOBPOFFSPRING[[goID]],
                                        GOCCOFFSPRING[[goID]],
                                        GOMFOFFSPRING[[goID]])))
            gofilter <- gofilter[which(!is.na(gofilter))]
            if (sum(gofilter %in% AllTerms)) {
              temp[[goID]] <- "-"
              wtst <- grsep2(gofilter, PG$`GO-ID`)
              #which(vapply(strsplit(PG$`GO-ID`[-wtst],";"), function(x) { sum(x %in% gofilter) }, 1) != 0)
              #which(vapply(strsplit(PG$`GO-ID`[wtst],";"), function(x) { sum(x %in% gofilter) }, 1) == 0)
              wtst2 <- which(temp$id %in% PG$id[wtst])
              temp[wtst2, goID] <- "+"
              #aggregate(temp$`GO:0043005`, list(temp$Sample), function(x) { setNames(aggregate(x, list(x), length)$x, c("-", "+")) })
              #View(PG[wtst, grep(topattern(PG.int.col), colnames(PG))])
            }
          }
        }
        invisible(clusterCall(parClust, function() {
          library(ggplot2)
          library(plotly)
          library(AnnotationDbi)
          library(htmlwidgets)
          return()
        }))
        clusterExport(parClust, c("QuantType", "QuantTypes", "Exp", "temp", "SubDir", "GO_filt", "GO_filter", "wd",
                                  "colScale", "fillScale", "colScale2", "fillScale2", "kolnm", "abbrFun"), envir = environment())
        if (runRankAbundPlots) {
          #for (exp in Exp[Wh]) { #exp <- Exp[Wh][1]
          tst <- parSapply(parClust, Exp[Wh], function(exp) {
            temp2 <- temp[which(temp$Sample == exp),]
            temp2 <- temp2[order(temp2$Y, decreasing = TRUE),]
            temp2$Protein_Group <- factor(temp2$Protein_Group, levels = temp2$Protein_Group)
            intmin <- floor(min(temp2$Y))
            intmax <- ceiling(max(temp2$Y))
            intscale <- intmax-intmin
            #xmax <- max(c(max(c(0, which(temp2$Category == "+")))+round(nrow(temp2)/15), nrow(temp2)))
            xmax <- nrow(temp2)*18/15
            PltTst <- setNames(c(TRUE, "+" %in% temp2$"In list"), c("All", "List"))
            if (GO_filt) { for (goID in GO_filter) {
              PltTst[goID] <- goID %in% colnames(temp2)
            } }
            for (i in 1:length(PltTst)) { #i <- 1
              if (PltTst[i]) {
                temp3 <- temp2
                ttl2 <- ttl <- paste0("Ranked abundance plots ",
                                      c("LFQ", "coverage", "spectral counts")[match(QuantType, QuantTypes)], " - ", exp)
                if (i == 1) {
                  catnm <- "Category"
                  filt <- 1:nrow(temp3)
                }
                if (i == 2) {
                  ttl2 <- ttl <- paste0(ttl, ", proteins of interest")
                  temp3$"In list" <- factor(temp3$"In list", levels = c("-", "+"))
                  catnm <- "In list"
                  filt <- which(temp3[[catnm]] == "+")
                }
                if (i > 2) {
                  goID <- names(PltTst)[i]
                  goID2 <- gsub("^GO:", "GO", goID)
                  nm <- AnnotationDbi::Term(goID)
                  temp3[[goID2]] <- factor(temp3[[goID]], levels = c("-", "+"))
                  ttl <- paste0(ttl, ", ", goID, " ", nm)
                  ttl2 <- gsub("/|:|\\*|\\?|<|>|\\|", "-", ttl)
                  myColors3 <- setNames(c("lightgrey", "purple"), c("-", "+"))
                  colScale3 <- scale_colour_manual(name = goID, values = myColors3)
                  fillScale3 <- scale_fill_manual(name = goID, values = myColors3)
                  catnm <- goID
                  filt <- which(temp3[[catnm]] == "+")
                }
                temp3$y <- temp3$Y + intmax*0.01
                plot <- ggplot(temp3) + geom_bar(stat = "identity",
                                                 aes(Protein_Group, Y, fill = .data[[catnm]], text = Value)) +
                  annotate("text", nrow(temp3)/2, intmax*1.2+intscale*0.04, label = paste0(length(unique(temp3$id)), " Protein Groups"),
                           hjust = 0.5, ) + theme_bw() + ggtitle(ttl) + ylab(kolnm) + 
                  coord_cartesian(xlim = c(1, xmax), ylim = c(intmin, intmax*1.2+intscale*0.05)) +
                  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks = element_blank(),
                        plot.margin = margin(r = 100))
                if (i == 1) { plot <- plot + colScale + fillScale }
                if (i == 2) { plot <- plot + colScale2 + fillScale2 }
                if (i > 2) { plot <- plot + colScale3 + fillScale3 }
                #poplot(plot, 12, 22)
                plot_ly <- ggplotly(plot, tooltip = c("Protein_Group", "text"))
                setwd(SubDir) # For some reason, unless I do this the default selfcontained = TRUE argument gets ignored and
                # a folder with external resources is created for each html plot!
                tst <- try(saveWidget(partial_bundle(plot_ly), paste0(SubDir, "/", ttl2, ".html")), silent = TRUE)
                if ("try-error" %in% class(tst)) { tst <- try(saveWidget(plot_ly, paste0(ttl2, ".html")), silent = TRUE) }
                #if ((i == 1)&&(!"try-error" %in% class(tst))) { system(paste0("open \"", SubDir, "/", ttl2, ".html")) }
                plot <- plot +
                  geom_text(data = temp3[filt,], angle = 45, hjust = 0, cex = 3.5,
                            aes(Protein_Group, y, colour = .data[[catnm]], label = Protein_Group))
                #poplot(plot, 12, 22)
                #setwd(SubDir)
                ggsave(filename = paste0(SubDir, "/", ttl2, ".jpeg"), plot, dpi = 600, width = 30, height = 10, units = "in")
                ggsave(filename = paste0(SubDir, "/", ttl2, ".pdf"), plot, dpi = 600, width = 30, height = 10, units = "in")
                setwd(wd)
              }
            }
          })
        }
        #}
        if ((runProfPlots)&&(length(g) > 1)) {
          PltTst <- setNames(c(TRUE, "+" %in% temp$`In list`), c("All", "List"))
          if (GO_filt) { for (goID in GO_filter) {
            PltTst[gsub(":", "", goID)] <- gsub(":", "", goID) %in% colnames(temp)
          } }
          for (i in 1:length(PltTst)) { #i <- 1
            MainDir <- paste0(wd, "/Profile plots")
            SubDir <- paste0(MainDir, "/", QuantType)
            if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
            if (PltTst[i]) {
              temp2 <- temp
              temp2$Category <- gsub(" *[\\(\\[].*", "", temp2$Category)
              temp2$Category <- as.factor(temp2$Category)
              ttl2 <- ttl <- paste0("Protein group ", c("LFQ", "coverage", "spectral count")[match(QuantType, QuantTypes)],
                                    " profiles")
              if (i == 1) {
                lev <- levels(temp2$Category)
                w <- which(!lev %in% c("In list", "Contaminant"))
                lev <- c(c("In list", "Contaminant"), lev[w])
                w <- which(!lev %in% c("In list", "Contaminant"))
                lev2 <- lev
                lev2[w] <- sapply(strsplit(as.character(lev[w]), ""), abbrFun)
                lev2[which(lev2 == "Contaminant")] <- "Cont."
                lev2 <- gsub(";", "\n", lev2)
                temp2$Category <- as.character(temp2$Category)
                temp2$Category <- lev2[match(temp2$Category, lev)]
                temp2$Category <- factor(temp2$Category, levels = lev2)
                catnm <- "Category"
                filt <- 1:nrow(temp2)
              }
              if (i == 2) {
                ttl2 <- ttl <- paste0(ttl, ", proteins of interest")
                temp2$"In list" <- factor(temp2$"In list", levels = c("-", "+"))
                catnm <- "In list"
                filt <- which(temp2[[catnm]] == "+")
              }
              if (i > 2) {
                goID <- names(PltTst)[i]
                goID2 <- gsub("^GO:", "GO", goID)
                nm <- AnnotationDbi::Term(goID)
                temp2[[goID]] <- factor(temp2[[goID]], levels = c("-", "+"))
                ttl <- paste0(ttl, ", ", goID, " ", nm)
                ttl2 <- gsub("/|:|\\*|\\?|<|>|\\|", "-", ttl)
                myColors3 <- setNames(c("lightgrey", "purple"), c("-", "+"))
                colScale3 <- scale_colour_manual(name = goID, values = myColors3)
                fillScale3 <- scale_fill_manual(name = goID, values = myColors3)
                catnm <- goID
                filt <- which(temp2[[catnm]] == "+")
              }
              m <- match(temp2$Category, c("-", "Contaminant", tstOrg3, "+", "In list"))
              temp2$LineType <- c(rep("dotted", 2), rep("dashed", length(tstOrg3)), rep("solid", 2))[m]
              #temp2$DotSize <- c(rep(0.2, 2), rep(0.3, length(tstOrg3)), rep(0.5, 2))[m]
              temp2$DotSize <- 1
              #temp2$Alpha <- c(rep(0.25, 2), rep(0.5, length(tstOrg3)), rep(1, 2))[m]
              temp2$Alpha <- 1
              Ngl <- c(0, 90)[(length(levels(temp2[[catnm]])) > 5) + 1]
              temp2 <- temp2[which(!is.na(temp2$Y)),]
              wTxt <- which(temp2$Sample == rev(levels(temp2$Sample))[1])
              frm <- as.formula(paste0("~`", catnm, "`"))
              plot <- ggplot(temp2, aes(text1 = Protein_Group, text2 = Value)) +
                geom_line(aes(x = Sample, y = Y, group = id, color = id), alpha = 0.1) +
                geom_point(aes(x = Sample, y = Y, size = DotSize, color = id)) + ggtitle(ttl) + ylab(kolnm) +
                theme_bw() + facet_grid(frm) +
                scale_size_identity(guide = "none") + scale_linetype_identity(guide = "none") +
                theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                      strip.text = element_text(face = "bold", size = 8, lineheight = 0.8, angle = Ngl),
                      strip.background = element_rect(fill = "lightblue", colour = "black", linewidth = 1)) +
                scale_alpha_identity(guide = "none") + scale_color_viridis(option = "D")
              plotxt <- plot + geom_text(data = temp2[wTxt,], aes(label = Protein_Group, x = Sample, y = Y, alpha = Alpha, color = id), hjust = 0, cex = 2)
              poplot(plotxt, 12, 20)
              setwd(SubDir)
              plotlyProfiles <- ggplotly(plot, tooltip = c("text1", "text2"))
              tst <- try(saveWidget(partial_bundle(plotlyProfiles), paste0(ttl2, ".html")), silent = TRUE)
              if ("try-error" %in% class(tst)) { tst <- try(saveWidget(plotlyProfiles, paste0(ttl2, ".html")), silent = TRUE) }
              if ((i == 1)&&(!"try-error" %in% class(tst))) { system(paste0("open \"", SubDir, "/", ttl2, ".html")) }
              ggsave(paste0(SubDir, "/", ttl2, ".jpeg"), plot, dpi = 150)
              ggsave(paste0(SubDir, "/", ttl2, ".pdf"), plot, dpi = 150)
              setwd(wd)
            }
          }
        }
      }
    }
  }
  PG$temp <- NULL
  invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
}

# Similar profiles but at peptides level
# This chunk has been vastly improved, and the others should be improved on the same model!!! 
source(parSrc, local = FALSE)
clusterExport(parClust, "abbrFun", envir = environment())
plotPepProf %<o% TRUE # For now, should come under control of a parameter eventually
if (plotPepProf) {
  #
  # This can be best parallelized thusly:
  # - Create data.frame with one row per combination of sample / quant_type / way to plot
  # - Parallelize over those
  #
  pepQuantTypes %<o% "intensities"
  PltTst <- c("All", "List")[1:(1+("+" %in% pep$`In list`))]
  plotTypes <- c("jpeg", "pdf", "html")
  l1 <- length(Exp)
  l2 <- length(pepQuantTypes)
  l3 <- length(PltTst)
  l4 <- length(plotTypes)
  plotsDF <- data.frame(Exp = rep(Exp, l2*l3*l4))
  plotsDF$Type <- as.character(sapply(pepQuantTypes, function(x) { rep(x, l1) }))
  plotsDF$Mode <- as.character(sapply(PltTst, function(x) {  rep(x, l1*l2) }))
  plotsDF$Ext <- as.character(sapply(plotTypes, function(x) {  rep(x, l1*l2*l3) }))
  #
  tmp <- listMelt(strsplit(pep$Proteins, ";"), 1:nrow(pep))
  m <- match(tmp$value, db$`Protein ID`)
  tmp$name <- db$`Common Name`[m]
  w <- which(!nchar(tmp$name) == 0)
  tmp$name[w] <- do.call(paste, c(tmp[w, c("name", "value")], sep = " "))
  w <- which(nchar(tmp$name) == 0)
  tmp$name[w] <- tmp$value[w]
  tmp3 <- aggregate(tmp$name, list(tmp$L1), paste, collapse = "\n") #This one is faster than with data.table!
  tmp3 <- data.frame(Seq = gsub("^_|_$", "", pep$`Modified sequence`), Name = tmp3$x[match(1:nrow(pep), tmp3$Group.1)])
  pep$Peptide_ID <- do.call(paste, c(tmp3, sep = "\n"))
  pep$"Peptide ID" <- gsub("\n", " ", pep$Peptide_ID)
  if (tstOrg) {
    tmp$Cont <- db$`Potential contaminant`[m]
    tmp$Org <- db[m, dbOrgKol]
    tmp1 <- aggregate(tmp$Cont, list(tmp$L1), function(x) {
      if ("+" %in% x) { x <- "+" } else { x <- "" }
      x
    })
    tmp2 <- aggregate(tmp$Org, list(tmp$L1), function(x) {
      if ("Contaminant" %in% x) { x <- "Contaminant" } else { x <- paste(unique(x), collapse = ";") }
      x
    })
    pep$`Potential contaminant` <- tmp1$x[match(1:nrow(pep), tmp1$Group.1)]
    pep$temp <- tmp2$x[match(1:nrow(pep), tmp2$Group.1)]
  }
  if (prot.list.Cond) {
    pep$"In list" <- ""
    pep$"In list"[grsep2(prot.list, pep$Proteins)] <- "+"
  }
  #
  plotsDF$kolnm <- plotsDF$root <- NA
  w <- which(plotsDF$Type == "intensities")
  plotsDF$kolnm[w] <- "Intensity"
  plotsDF$root[w] <- paste0(plotsDF$kolnm[w], " - ")
  w <- which(is.na(plotsDF$root))
  if (length(w)) { stop("This case is not addressed in the code yet!") }
  plotsDF$kol <- do.call(paste, c(plotsDF[, c("root", "Exp")], sep = ""))
  #
  # Prepare data
  temp_pep <- dfMelt(pep[, unique(plotsDF$kol)])
  colnames(temp_pep)[which(colnames(temp_pep) == "value")] <- "Y"
  temp_pep$variable <- as.character(temp_pep$variable)
  varkol <- c("Sequence", "Proteins", "Peptide ID", "Peptide_ID", "id")
  if (prot.list.Cond) { varkol <- c(varkol, "In list") }
  temp_pep[, varkol] <- pep[, varkol]
  if (prot.list.Cond) {
    temp_pep$`In list`[which(temp$`In list` == "")] <- "-"
    temp_pep$"In list" <- factor(temp_pep$"In list", levels = c("-", "+"))
  }
  if (length(tstOrg2)) { temp_pep$Category <- pep$temp[match(temp_pep$id, pep$id)] } else {
    temp_pep$Category <- "-"
  }
  temp_pep <- temp_pep[which(temp_pep$Y > 0),]
  temp_pep <- temp_pep[which(is.all.good(temp_pep$Y, 2)),]
  g <- grep("^Intensity - ", temp_pep$variable)
  temp_pep$Y[g] <- log10(temp_pep$Y[g])
  temp_pep$variable[g] <- gsub_Rep("^Intensity - ", "log10(Int.) - ", temp_pep$variable[g])
  g <- grep("^Intensity - ", plotsDF$kol)
  plotsDF$kol[g] <- gsub_Rep("^Intensity - ", "log10(Int.) - ", plotsDF$kol[g])
  plotsDF$root[g] <- "log10(Int.) - "
  temp_pep$Sample <- gsub_Rep(".* - ", "", temp_pep$variable)
  lev <- c("In list", "+", tstOrg2, "-", "Contaminant")
  lev <- lev[which(lev %in% unique(temp_pep$Category))]
  temp_pep$Category <- factor(temp_pep$Category, levels = lev)
  #
  # Serialize data
  tmp <- aggregate(1:nrow(temp_pep), list(temp_pep$Sample), list)
  sapply(Exp, function(exp) {
    write_rds(temp_pep[tmp$x[[which(tmp$Group.1 == exp)]],], paste0(wd, "/tmp_", exp, ".RDS"))
  })
  #
  MainDir <- paste0(wd, "/Ranked abundance")
  #
  # Scales
  myColorsB <- myColors <- setNames("black", "-")
  myColors2 <- setNames(c("lightgrey", "brown"), c("-", "+"))
  if (length(tstOrg2)) {
    myColorsB <- myColors <- setNames(grDevices::colorRampPalette(c("blue", "green"))(length(tstOrg2)), tstOrg2)
    w <- which(names(myColorsB) != "In list")
    names(myColorsB)[w] <- gsub("[a-z]+ ", ". ", names(myColorsB)[w])
    names(myColorsB) <- gsub(";", "\n", names(myColorsB))
    myColorsB[["potential contaminant"]] <- myColors[["potential contaminant"]] <- "grey"
  }
  myColorsB[[c("+", "In list")[(length(tstOrg2) > 0)+1]]] <- myColors[[c("+", "In list")[(length(tstOrg2) > 0)+1]]] <- "red"
  colScale <- ggplot2::scale_colour_manual(name = "Category", values = myColors)
  colScaleB <- ggplot2::scale_colour_manual(name = "Category", values = myColorsB)
  fillScale <- ggplot2::scale_fill_manual(name = "Category", values = myColors)
  colScale2 <- ggplot2::scale_colour_manual(name = "In list", values = myColors2)
  fillScale2 <- ggplot2::scale_fill_manual(name = "In list", values = myColors2)
  #
  clusterExport(parClust, list("plotsDF", "tstOrg2", "wd", "MainDir", "myColorsB", "myColors2", "colScale", "colScaleB",
                               "fillScale", "fillScale2"), envir = environment())
  cat("Plotting peptide ranked abundance plots\n")
  pepSortTst <- parLapply(parClust, 1:nrow(plotsDF), function(ii) { #ii <- 1
    exp <- plotsDF$Exp[ii]
    pepQuantType <- plotsDF$Type[ii]
    plotMode <- plotsDF$Mode[ii]
    kolnm <- plotsDF$kolnm[ii]
    root <- plotsDF$root[ii]
    fileType <- plotsDF$Ext[ii]
    #
    SubDir <- paste0(MainDir, "/Peptide ", pepQuantType)
    if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
    #
    # Load data
    temp <- readr::read_rds(paste0(wd, "/tmp_", exp, ".RDS"))
    nr <- nrow(temp)
    test <- aggregate(temp$Peptide_ID, list(temp$Peptide_ID), length)
    w <- which(test$x > 1)
    if (length(w)) {
      test <- test[w,]
      for (i in test$Group.1) {
        w <- which(temp$Peptide_ID == i)
        temp$Peptide_ID[w[2:length(w)]] <- paste0(temp$Peptide_ID[w[2:length(w)]], "_", 2:length(w))
        temp$"Peptide ID"[w[2:length(w)]] <- paste0(temp$"Peptide ID"[w[2:length(w)]], "_", 2:length(w))
      }
    }
    #
    temp$Value <- paste0(kolnm, ": ", temp$Y)
    temp <- temp[order(temp$Y, decreasing = TRUE),]
    temp$Peptide_ID <- factor(temp$Peptide_ID, levels = unique(temp$Peptide_ID))
    temp$"Peptide ID" <- factor(temp$"Peptide ID", levels = unique(temp$"Peptide ID"))
    intmin <- floor(min(temp$Y))
    intmax <- ceiling(max(temp$Y))
    intscale <- intmax-intmin
    xmax <- nr*18/15
    ttl <- paste0("Peptide ", gsub(" - $", "", root), " ranked abundance plots - ", exp)
    if (plotMode == "All") {
      catnm <- "Category"
      filt <- 1:nr
    }
    if (plotMode == "List") {
      ttl <- paste0(ttl, ", proteins of interest")
      catnm <- "In list"
      filt <- which(temp[[catnm]] == "+")
    }
    temp$y <- temp$Y + intmax*0.01
    flPath <- paste0(SubDir, "/", ttl, ".", fileType)
    plot <- ggplot2::ggplot(temp) +
      ggplot2::geom_bar(stat = "identity",
                        ggplot2::aes(x = Peptide_ID, y = Y, fill = .data[[catnm]], text = Value)) +
      ggplot2::annotate("text", nr/2, intmax*1.2+intscale*0.04, label = paste0(length(unique(temp$id)), " peptidoforms"),
                        hjust = 0.5) + ggplot2::theme_bw() + ggplot2::ggtitle(ttl) + ggplot2::ylab(kolnm) + 
      ggplot2::coord_cartesian(xlim = c(1, xmax),
                               ylim = c(intmin, intmax*1.2+intscale*0.05)) +
      ggplot2::theme(panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(r = 100))
    if (plotMode == "All") { plot <- plot + colScale + fillScale }
    if (plotMode == "List") { plot <- plot + colScale2 + fillScale2 }
    if (fileType %in% c("jpeg", "pdf")) {
      plot <- plot +
        ggplot2::geom_text(data = temp[filt,], angle = 45, hjust = 0, cex = 3.5, show.legend = FALSE,
                           ggplot2::aes(Peptide_ID, y, colour = .data[[catnm]], label = `Peptide ID`))
      #poplot(plot, 12, 22)
      ggplot2::ggsave(filename = flPath, plot, dpi = 600, width = 30, height = 10, units = "in")
    }
    if (fileType == "html") {
      plot_ly <- plotly::ggplotly(plot, tooltip = c("Peptide_ID", "text"))
      setwd(SubDir) # For some reason, unless I do this the default selfcontained = TRUE argument gets ignored and
      # a folder with external resources is created for each html plot!
      tst <- try(htmlwidgets::saveWidget(plotly::partial_bundle(plot_ly), flPath), silent = TRUE)
      if ("try-error" %in% class(tst)) { tst <- try(htmlwidgets::saveWidget(plot_ly, flPath), silent = TRUE) }
      setwd(wd)
    }
  })
  unlink(paste0(wd, "/tmp_", Exp, ".RDS"))
  if (length(Exp) > 1) {
    plotsDF2 <- plotsDF
    plotsDF2$Exp <- NULL
    plotsDF2$kol <- NULL
    tst <- do.call(paste, c(plotsDF2, sep = "_______"))
    tst <- aggregate(1:nrow(plotsDF2), list(tst), min)
    plotsDF2 <- plotsDF2[tst$x,]
    # temp_pep is all the data we need
    MainDir <- paste0(wd, "/Profile plots")
    clusterExport(parClust, list("plotsDF2", "MainDir", "tstOrg3"), envir = environment())
    write_rds(temp_pep, paste0(wd, "/tmp.RDS"))
    cat("Plotting peptide profile plots\n")
    pepProfTst <- parLapply(parClust, 1:nrow(plotsDF2), function(ii) { #ii <- 1
      pepQuantType <- plotsDF2$Type[ii]
      plotMode <- plotsDF2$Mode[ii]
      kolnm <- plotsDF2$kolnm[ii]
      root <- plotsDF2$root[ii]
      fileType <- plotsDF2$Ext[ii]
      #
      temp <- readr::read_rds(paste0(wd, "/tmp.RDS"))
      temp$Value <- temp$Y
      #
      SubDir <- paste0(MainDir, "/Peptide ", pepQuantType)
      if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
      #
      ttl <- paste0("Peptides ", gsub(" - $", "", root), " profiles")
      if (plotMode == "All") {
        lev <- levels(temp$Category)
        w <- which(!lev %in% c("In list", "Contaminant"))
        lev <- c(c("In list", "Contaminant"), lev[w])
        w <- which(!lev %in% c("In list", "Contaminant"))
        lev2 <- lev
        lev2[w] <- sapply(strsplit(as.character(lev[w]), ""), abbrFun)
        lev2[which(lev2 == "Contaminant")] <- "Cont."
        lev2 <- gsub(";", "\n", lev2)
        temp$Category <- as.character(temp$Category)
        temp$Category <- lev2[match(temp$Category, lev)]
        temp$Category <- factor(temp$Category, levels = lev2)
        catnm <- "Category"
        filt <- 1:nrow(temp)
      }
      if (plotMode == "List") {
        ttl <- paste0(ttl, ", proteins of interest")
        catnm <- "In list"
        filt <- which(temp[[catnm]] == "+")
      }
      m <- match(temp$Category, c("-", "Contaminant", tstOrg3, "+", "In list"))
      temp$LineType <- c(rep("dotted", 2), rep("dashed", length(tstOrg3)), rep("solid", 2))[m]
      #temp$DotSize <- c(rep(0.2, 2), rep(0.3, length(tstOrg3)), rep(0.5, 2))[m]
      temp$DotSize <- 1
      #temp$Alpha <- c(rep(0.25, 2), rep(0.5, length(tstOrg3)), rep(1, 2))[m]
      temp$Alpha <- 1
      Ngl <- c(0, 90)[(length(levels(temp[[catnm]])) > 5) + 1]
      wTxt <- which(temp$Sample == rev(levels(temp$Sample))[1])
      frm <- as.formula(paste0("~`", catnm, "`"))
      flPath <- paste0(SubDir, "/", ttl, ".", fileType)
      plot <- ggplot2::ggplot(temp, ggplot2::aes(text1 = Peptide_ID, text2 = Value)) +
        ggplot2::geom_line(ggplot2::aes(x = Sample, y = Y, group = id, color = id), alpha = 0.1, show.legend = FALSE) +
        ggplot2::geom_point(ggplot2::aes(x = Sample, y = Y, size = DotSize, color = id)) +
        ggplot2::ggtitle(ttl) + ggplot2::ylab(kolnm) +
        ggplot2::theme_bw() + ggplot2::facet_grid(frm) +
        ggplot2::scale_size_identity(guide = "none") +
        ggplot2::scale_linetype_identity(guide = "none") +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                       strip.text = ggplot2::element_text(face = "bold", size = 8, lineheight = 0.8, angle = Ngl),
                       strip.background = ggplot2::element_rect(fill = "lightblue", colour = "black", linewidth = 1)) +
        ggplot2::scale_alpha_identity(guide = "none") +
        viridis::scale_color_viridis(option = "D")
      #proteoCraft::poplot(plot, 12, 20)
      if (fileType %in% c("jpeg", "pdf")) {
        plot <- plot +
          ggplot2::geom_text(data = temp[wTxt,], ggplot2::aes(label = `Peptide ID`, x = Sample, y = Y, alpha = Alpha, color = id),
                             hjust = 0, cex = 2)
        #proteoCraft::poplot(plot, 12, 22)
        ggplot2::ggsave(filename = flPath, plot, dpi = 600, width = 30, height = 10, units = "in")
      }
      if (fileType == "html") {
        plot_ly <- plotly::ggplotly(plot, tooltip = c("Peptide_ID", "text")) # Super slowwwwww
        setwd(SubDir) # For some reason, unless I do this the default selfcontained = TRUE argument gets ignored and
        # a folder with external resources is created for each html plot!
        tst <- try(htmlwidgets::saveWidget(plotly::partial_bundle(plot_ly), flPath), silent = TRUE)
        if ("try-error" %in% class(tst)) { tst <- try(htmlwidgets::saveWidget(plot_ly, flPath), silent = TRUE) }
        setwd(wd)
      }
    })
  }
  unlink(paste0(wd, "/tmp.RDS"))
  pep$temp <- NULL
}
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))

# Negative filter
if (NegFilt) {
  e <- ev[grep("-MATCH$", ev$Type, invert = TRUE),]
  nf <- SamplesMap$Experiment[which(SamplesMap$"Negative Filter")]
  e <- e[which(e$Experiment %in% nf),]
  PG$"Direct identification in negative filter sample(s)" <- vapply(strsplit(PG$`Evidence IDs`, ";"), function(x) { sum(unlist(x) %in% e$id) }, 1) > 0
  PG$"Direct identification in negative filter sample(s)" <- c("", "+")[match(PG$"Direct identification in negative filter sample(s)", c(FALSE, TRUE))]
  if (MakeRatios) {
    exp <- SamplesMap$Experiment[which(!SamplesMap$Reference)]
    PG[which(PG$"Direct identification in negative filter sample(s)" == "+"),
       paste0(c("Enriched", "Regulated")[RatiosThresh_2sided+1], " - ", exp)] <- ""
  }
}

Script <- readLines(ScriptPath)
gc()
# It makes sense to close/re-create parallel clusters regularly to reduce memory usage + avoid corruption
stopCluster(parClust)
source(parSrc, local = FALSE)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

#### Code chunk - XML coverage columns
Src <- paste0(libPath, "/extdata/R scripts/Sources/xml_Coverage_columns.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
# Calculate maximum expected coverage per group
if (WorkFlow == "Band ID") {
  m <- match(frstProt, db$`Protein ID`)
  Dig <- data.frame(ID = frstProt, Seq = db$Sequence[m])
  Dig$Digest <- Digest(Dig$Seq, min = MinPepSz, missed = Missed, cl = parClust)
  Dig$MaxCov <- parApply(parClust, Dig[, c("Seq", "Digest")], 1, function(x) {
    Coverage(x[[1]], x[[2]])
  })
  PG$`Max. theoretical sequence coverage [%]` <- round(100*Dig$MaxCov, 1)
}

#### Code chunk - Create output tables
## PSMs
dir <- paste0(wd, "/Tables")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
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
for (nm in names(int.cols)) { #nm <- names(pep.ref[1])
  rpl <- intNms(nm, type = "pep")
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
for (nm in names(PG.int.cols)) { #nm <- names(pep.ref[1])
  rpl <- intNms(nm)
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
if (MakeRatios) {
  for (nm in unique(c(names(rat.cols), names(PG.rat.cols)))) { #nm <- names(rat.cols[1])
    rpl <- ratNms(nm)
    Styles[[paste0(rpl, ", avg.")]] <- "Summary Ratios"
    Styles[[paste0(rpl, ", indiv.")]] <- "Individual Ratios"
  }
}
fl <- system.file("extdata", "Report - column names - no replicates.xlsx", package = "proteoCraft")
styleNms <- openxlsx2::read_xlsx(fl, "tmp", colNames = FALSE)[,1]
# wb <- loadWorkbook(fl)
# addWorksheet(wb, "tmp")
# tmpFl <- temp_xlsx()
#w <- which(vapply(Styles, function(x) { "Style" %in% class(x) }, TRUE))
# styleNms %<o% names(Styles)[w]
# writeData(wb, "tmp", styleNms)
# for (i in 1:length(w)) { addStyle(wb, "tmp", Styles[[w[i]]], i, 1) }
# saveWorkbook(wb, tmpFl)
# openXL(tmpFl)
# WorkBook %<o% wb_load(tmpFl)
# Styles2 %<o% setNames(w, names(Styles)[w])
WorkBook %<o% wb_load(fl)
repFl <- paste0(wd, "/Tables/Report_", dtstNm, ".xlsx")
WorkBook <- wb_add_data(WorkBook, "Description", dtstNm, wb_dims(2, 5))
WorkBook <- wb_add_data(WorkBook, "Description", format(Sys.Date(), "%d/%m/%Y"), wb_dims(3, 5))
WorkBook <- wb_add_data(WorkBook, "Description", WhoAmI, wb_dims(4, 5))
tmp <- loadedPackages(TRUE)
WorkBook <- wb_add_data(WorkBook, "Description", tmp$Version[grep("proteoCraft", tmp$Name)], wb_dims(5, 5))
cat(" - Writing Excel report...\n")
# Function for editing the header
KolEdit <- function(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl) {
  #KolNames <- xlTabs[[sheetnm]]
  klnms <- KolNames
  KolNames[which(KolNames == "Evidence IDs")] <- paste0("All ", evNm, " IDs")
  for (nm in names(intTbl)) { #nm <- names(intTbl)[1] #nm <- names(intTbl)[2]
    m <- match(intTbl[[nm]]$Log, KolNames)
    w <- which(!is.na(m))
    if (length(w)) {
      rpl <- intNms(nm, type = "pep")
      KolNames[m[w]] <- paste0(rpl, " ", intTbl[[nm]]$Sample[w])
    }
  }
  if (MakeRatios) {
    for (nm in names(ratTbl)) {
      m <- match(ratTbl[[nm]]$Log, KolNames)
      w <- which(!is.na(m))
      if (length(w)) {
        rpl <- ratNms(nm)
        KolNames[m[w]] <- paste0(rpl, " ", ratTbl[[nm]]$Sample[w])
      }
    }
    KolNames <- gsub("^Enriched", "Enr.", KolNames)
    KolNames <- gsub("^Regulated", "Reg.", KolNames)
  }
  KntKol <- paste0(AA, " Count")
  KolNames[which(KolNames %in% KntKol)] <- gsub(" Count$", "", KolNames[which(KolNames %in% KntKol)])
  #
  KolNames <- gsub("( - )|(___)", " ", KolNames)
  #
  # Those names must be unique if the data is to be written as a table!
  tst <- aggregate(KolNames, list(KolNames), length)
  stopifnot(max(tst$x) == 1)
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
QualFilt %<o% c(#pgOrgKol,
  "Potential contaminant", "Only identified by site",
                grep("^Quality filter: ", colnames(PG), value = TRUE))
if (NegFilt) { QualFilt <- c(QualFilt, "Direct identification in negative filter sample(s)") }
II <- setNames(1, "All peptidoforms")
if ((length(Mod2Write))&&(PTMriched)) {
  II[paste0(Modifs$`Full name`[match(Mod2Write, Modifs$Mark)], "-mod. pept.")] <- 1+(seq_along(length(Mod2Write)))
}
for (ii in II) { #ii <- II[1] #ii <- II[2]
  tblMode <- tblMode2 <- "pep"
  TbNm <- names(II)[ii]
  tempData <- get(tblMode)
  intRf <- int.cols
  ratRf <- rat.cols
  if (ii > 1) {
    ptm <- Mod2Write[[ii-1]]
    Ptm <- Modifs$`Full name`[match(ptm, Modifs$Mark)]
    tempData <- PTMs_pep[[Ptm]]
    intRf <- PTMs_intRf
    ratRf <- PTMs_ratRf
    tblMode2 <- paste0(Ptm, "-modified pep")
  }
  # names(intRf) <- gsub("^Original ", "", paste0(names(intRf), " intensities"))
  # names(ratRf) <- gsub("^Original ", "", paste0(names(ratRf), " ratios"))
  # (In this workflow, as of this version all peptide intensities are non-log-transformed but ratios are!)
  if (nrow(tempData)) {
    CoreCol <- "Modified sequence"
    if ("Modified sequence_verbose" %in% colnames(tempData)) { CoreCol <- c(CoreCol, "Modified sequence_verbose") }
    CoreCol <- c(CoreCol, "Sequence", "id", "Proteins")
    CoreCol2 <- c("Leading proteins", "Leading razor proteins",
                  "Protein names", "Gene names", "Protein group IDs", "Razor protein group ID",
                  grep(" (Probabilities|Score Diffs)$", colnames(tempData), value = TRUE),
                  "Normalisation group")
    evcol <- "Evidence IDs"
    spcol <- "MS/MS count"
    intColsTbl <- setNames(lapply(names(intRf), function(nm) {
      res <- data.frame(nonLog = c(intRf[nm], paste0(intRf[nm], " - ", Exp)),
                        Log = c(paste0("log10(", intRf[nm], ")"), paste0("log10(", intRf[nm], ") - ", Exp)),
                        Type = c("Average", rep("Individual", length(Exp))),
                        Sample = c("Average", Exp))
      w <- which(res$nonLog %in% colnames(tempData))
      return(res[w,])
    }), names(intRf))
    quantCols <- intCols <- lapply(intColsTbl, function(x) { x$Log })
    gel0 <- unlist(lapply(intColsTbl, function(x) { x$nonLog }))
    gel <- unlist(intCols)
    tempData[, gel] <- log10(tempData[, gel0]) # Peptide expression values are not log-transformed
    tempData[, gel0] <- NULL
    for (gl in gel) {
      w <- which(is.infinite(tempData[[gl]]))
      tempData[w, gl] <- NA
    }
    quantcol <- gel
    for (nm in names(intRf)) { 
      if (!grepl("log[0-9]+", intRf[nm])) { # In case I am rerunning a bit
        intRf[nm] <- paste0("log10(", intRf[nm], ")")
      }
    }
    if (MakeRatios) {
      ratColsTbl <- setNames(lapply(names(ratRf), function(nm) {
        # Note: Peptide values are already log-transformed in this workflow!!!
        res <- data.frame(
          Log = c(paste0("Mean ", ratRf[nm]), paste0(ratRf[nm], " - ", Exp)),
          Type = c("Average", rep("Individual", length(Exp))),
          Sample = c("Average", Exp))
        w <- which(res$Log %in% colnames(tempData))
        return(res[w,])
      }), names(ratRf))
      ratCols <- lapply(ratColsTbl, function(x) { x$Log })
      grl <- unlist(ratCols)
      for (gr in grl) {
        w <- which(is.infinite(tempData[[gr]]))
        tempData[w, gr] <- NA
      }
      quantcol <- c(quantcol, grl)
      quantCols[names(ratRf)] <- ratCols
      regcol <- grep("^((Enriched)|(Regulated)) - ", colnames(tempData), value = TRUE)
    }
    aacol <- paste0(AA, " Count")
    qualFlt <- QualFilt[which(QualFilt %in% colnames(ev))]
    w <- which(!qualFlt %in% colnames(tempData))
    if (length(w)) {
      tempData[, qualFlt[w]] <- ev[match(tempData$"Modified sequence", ev$"Modified sequence"), qualFlt[w]]
    }
    dir <- paste0(wd, "/Tables")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    data.table::fwrite(tempData, paste0(dir, "/", TbNm, ".tsv"), sep = "\t", row.names = FALSE, na = "NA")
    w <- grsep2(prot.list, tempData$Proteins)
    if (length(w)) {
      data.table::fwrite(tempData[w,], paste0(wd, "/Tables/", TbNm, " - Proteins in list.tsv"),
                         sep = "\t", row.names = FALSE, na = "NA")
    }
    # Columns table
    # - IDs
    ColumnsTbl <- list(IDs = c(CoreCol, CoreCol2))
    # - Counts
    ColumnsTbl[["AA counts"]] <- aacol
    # - Evidence counts and IDs
    ColumnsTbl[["Global Ev. IDs"]] <- "Evidence IDs"
    ColumnsTbl[["Global Spec. counts"]] <- "MS/MS count"
    # - Expression values
    for (nm in names(intRf)) { #nm <- names(intRf[1])
      rpl <- intNms(nm, type = "pep")
      ColumnsTbl[[paste0(rpl, ", avg.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Average")]
      ColumnsTbl[[paste0(rpl, ", indiv.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Individual")]
    }
    if (MakeRatios) {
      for (nm in names(ratRf)) { #nm <- names(ratRf[1])
        rpl <- ratNms(nm)
        ColumnsTbl[[paste0(rpl, ", avg.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Average")]
        ColumnsTbl[[paste0(rpl, ", indiv.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Individual")]
      }
    }
    # - PEP
    ColumnsTbl[["PEP"]] <- "PEP"
    # - Filters
    ColumnsTbl[["Filters"]] <- qualFlt
    # Melt
    ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }, 1) > 0)]
    ColumnsTbl <- listMelt(ColumnsTbl, names(ColumnsTbl), c("Col", "Grp"))
    #aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), length)
    stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
    ColumnsTbl$Class <- ""
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Peptides information"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evcol))] <- "Evidence IDs"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(spcol))] <- "Spectral count"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% aacol)] <- "Amino Acid counts"
    ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters"))] <- "QC filters"
    for (nm in names(intRf)) { #nm <- names(intRf)[1]
      rpl <- intNms(nm, TRUE, type = "pep")
      ColumnsTbl$Class[grep(topattern(paste0(intRf[nm], " - ")), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[which(ColumnsTbl$Col == intRf[nm])] <- rpl
    }
    if (MakeRatios) {
      for (nm in names(ratRf)) { #nm <- names(ratRf)[1]
        rpl <- ratNms(nm, TRUE)
        ColumnsTbl$Class[grep(topattern(paste0(ratRf[nm], " - ")), ColumnsTbl$Col)] <- rpl
        ColumnsTbl$Class[which(ColumnsTbl$Col == ratRf[nm])] <- rpl
      }
      ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
    }
    ColumnsTbl$Class[grep("[Aa]nnotations", ColumnsTbl$Grp)] <- "Annotations"
    stopifnot(min(nchar(ColumnsTbl$Class)) > 0) #View(ColumnsTbl)
    w <- c(which(ColumnsTbl$Class == "General Peptides information"),
           unlist(lapply(names(intCols), function(nm) { which(ColumnsTbl$Class == intNms(nm, TRUE, type = "pep")) }))
    )
    if (MakeRatios) {
      w <- c(w,
             unlist(lapply(names(ratCols), function(nm) { which(ColumnsTbl$Class == ratNms(nm, TRUE)) })),
             which(ColumnsTbl$Class == "Regulated")
      )
    }
    w <- c(w, which(ColumnsTbl$Class == "QC filters"),
           which(ColumnsTbl$Class == "Evidence IDs"),
           which(ColumnsTbl$Class == "Spectral count"),
           which(ColumnsTbl$Class == "Amino Acid counts"))
    stopifnot(length(w) == nrow(ColumnsTbl))
    ColumnsTbl <- ColumnsTbl[w,]
    ColumnsTbl$Hide <- ColumnsTbl$Class %in% c("Spectral count", "Spectrum IDs",
                                               "Evidences count", "Evidence IDs",
                                               "Amino Acid counts", "Annotations")
    l <- length(intRf)
    if (l > 1) {
      for (nm in names(intRf)[1:(l-1)]) {
        ColumnsTbl$Hide[which(ColumnsTbl$Class == intNms(nm, TRUE, type = "pep"))] <- TRUE
      }
    }
    if (MakeRatios) {
      l <- length(ratRf)
      if (l > 1) {
        for (nm in names(ratRf)[1:(l-1)]) {
          ColumnsTbl$Hide[which(ColumnsTbl$Class == ratNms(nm, TRUE))] <- TRUE
        }
      }
    }
    #
    if (MakeRatios) { a <- KolEdit(ColumnsTbl$Col, intColsTbl, ratColsTbl) } else { a <- KolEdit(ColumnsTbl$Col, intColsTbl) }
    ColumnsTbl$edit_Col <- unlist(a)
    #
    Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
}
#saveFun(WorkBook, file = "WorkBook_bckp.RData")
#wb_save(WorkBook, repFl);xl_open(repFl)
#loadFun("WorkBook_bckp.RData")
TbNm <- "Protein groups"
tblMode <- tblMode2 <- "PG"
# Function for editing the header
KolEdit <- function(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl) {
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
    KolNames <- gsub(".*Regulated - ", "reg. ", KolNames)
  }
  #KolNames <- gsub(".*Pvalue\\)( - )?", "-log10 pval. ", KolNames)
  #KolNames <- gsub(".*Significant-", "signif. ", KolNames)
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
                  "Evidences count",
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
intRf <- PG.int.cols
names(intRf) <- paste0(names(intRf), " int.")
intColsTbl <- setNames(lapply(names(intRf), function(nm) { #nm <- names(intRf)[1]
  res <- data.frame(
    Log = c(paste0("Mean ", gsub(" - ", "", intRf[nm])), paste0(intRf[nm], Exp)),
    Type = c("Average", rep("Individual", length(Exp))),
    Sample = c("Average", Exp))
  w <- which(res$Log %in% colnames(tempData))
  return(res[w,])
}), names(intRf))
quantCols <- intCols <- lapply(intColsTbl, function(x) { x$Log })
gel <- unlist(intCols)
for (gl in gel) {
  w <- which(is.infinite(tempData[[gl]]))
  tempData[w, gl] <- NA
}
quantcol <- gel
if (MakeRatios) {
  ratRf <- PG.rat.cols
  names(ratRf) <- paste0(names(ratRf), " rat.")
  ratColsTbl <- setNames(lapply(names(ratRf), function(nm) {
    # Note: Peptide values are already log-transformed in this workflow!!!
    res <- data.frame(
      Log = c(paste0("Mean ", gsub(" - $", "", ratRf[nm])), paste0(ratRf[nm], Exp)),
      Type = c("Average", rep("Individual", length(Exp))),
      Sample = c("Average", Exp))
    w <- which(res$Log %in% colnames(tempData))
    return(res[w,])
  }), names(ratRf))
  ratCols <- lapply(ratColsTbl, function(x) { x$Log })
  grl <- unlist(ratCols)
  for (gr in grl) {
    w <- which(is.infinite(tempData[[gr]]))
    tempData[w, gr] <- NA
  }
  quantcol <- c(quantcol, grl)
  quantCols[names(ratRf)] <- ratCols
  regcol <- grep("^((Enriched)|(Regulated)) - ", colnames(tempData), value = TRUE)
}
if (protrul) { quantcol <- c(quantcol, grep(topattern("log10(est. copies/cell) - ", start = FALSE), colnames(tempData), value = TRUE)) }
covcol <- c(xmlCovCol,
            c("Sequence coverage [%]",
              "Unique + razor sequence coverage [%]",
              "Unique sequence coverage [%]")[1:c(1, 3)[isEukaLike+1]],
            grep(topattern("Sequence coverage [%] - "), colnames(tempData), value = TRUE)) # The complicated way, but ensures the order is correct
if (WorkFlow == "Band ID") {
  covcol <- c("Max. theoretical sequence coverage [%]", covcol)
}
kol <- c(kol, "Mol. weight [kDa]", covcol, "PEP", quantcol)
if ((exists("KlustKols"))&&(length(KlustKols))) { kol <- c(kol, KlustKols) }
qualFlt <- QualFilt
kol <- unique(c(kol, qualFlt))
if (Annotate) { kol <- c(kol, annot.col) }
kol <- unique(kol[which(kol %in% colnames(tempData))])
tempData <- tempData[, kol]
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
if (MakeRatios) {
  for (nm in names(ratRf)) { #nm <- names(ratRf[1])
    rpl <- ratNms(nm)
    ColumnsTbl[[paste0(rpl, ", avg.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Average")]
    ColumnsTbl[[paste0(rpl, ", indiv.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Individual")]
  }
  ColumnsTbl[["Regulated"]] <- regcol
}
# - Proteome Ruler
if (protrul) {
  ColumnsTbl[["Proteome Ruler"]] <- grep(topattern("log10(est. copies/cell) - ", start = FALSE), colnames(tempData), value = TRUE)
}
# - Annotations
if (Annotate) {
  annot.col2 <- gsub("_names$", " names", annot.col)
  AnnotTbl$Columns <- list(c("GO", "GO-ID"), c("Taxonomy", "TaxID"), NA, NA, NA, NA, "EMBL", NA)
  for (i in annot) { AnnotTbl$Columns[match(i, AnnotTbl$Name)] <- list(c(i, paste0(i, " names"))) }
  AnnotTbl$Columns[match("Other", AnnotTbl$Name)] <- list(annot.col2[which(!annot.col2 %in% unlist(AnnotTbl$Columns))])
  for (i in 1:nrow(AnnotTbl)) { ColumnsTbl[[paste0(AnnotTbl$Name[i], " annotations")]] <- AnnotTbl$Columns[[i]] }
}
# - PEP
ColumnsTbl[["PEP"]] <- "PEP"
# - Filters
ColumnsTbl[["Filters"]] <- qualFlt
ColumnsTbl[["In list"]] <- "In list"
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
ColumnsTbl$Class[which(ColumnsTbl$Grp == "In list")] <- "In list"
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
  kl <- c(paste0("Mean ", gsub(" - $", "", intRf[nm])), paste0(intRf[[nm]], Exp))
  kl <- kl[which(kl %in% ColumnsTbl$Col)]
  ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
}
ColumnsTbl$Class[grep("Ruler", ColumnsTbl$Grp)] <- "log10(est. copies/cell)"
if (MakeRatios) {
  for (nm in names(ratRf)) { #nm <- names(ratRf)[1]
    rpl <- ratNms(nm, TRUE)
    kl <- c(paste0("Mean ", gsub(" - $", "", ratRf[[nm]])), paste0(ratRf[[nm]], Exp))
    kl <- kl[which(kl %in% ColumnsTbl$Col)]
    ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
  }
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
}
ColumnsTbl$Class[grep("[Aa]nnotations", ColumnsTbl$Grp)] <- "Annotations"
ColumnsTbl$Class[grep("[Ss]equence coverage \\[%\\]", ColumnsTbl$Col)] <- "Sequence coverage [%]"
ColumnsTbl$Class[grep("^1st ID cov\\.", ColumnsTbl$Col)] <- "1st accession sequence coverage (peptides)"
if ((exists("KlustKols"))&&(length(KlustKols))) {
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "Cluster")] <- paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ")")
}
ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters", "Negative filter"))] <- "QC filters"
stopifnot(min(nchar(ColumnsTbl$Class)) > 0)
w <- c(which(ColumnsTbl$Class == "General Protein Group information"),
       unlist(lapply(names(intCols), function(nm) { which(ColumnsTbl$Class == intNms(nm, TRUE)) })),
       which(ColumnsTbl$Class == "log10(est. copies/cell)")
)
if (MakeRatios) {
  w <- c(w,
         unlist(lapply(names(ratCols), function(nm) { which(ColumnsTbl$Class == ratNms(nm, TRUE)) })),
         which(ColumnsTbl$Class == "Regulated")
  )
}
w <- c(w,
       which(ColumnsTbl$Class == "QC filters"),
       which(ColumnsTbl$Class == "In list"),
       which(ColumnsTbl$Class == "Sequence coverage [%]"),
       which(ColumnsTbl$Class == "1st accession sequence coverage (peptides)"),
       which(ColumnsTbl$Class == paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ")")),
       which(ColumnsTbl$Class == "Peptides count"),
       which(ColumnsTbl$Class == "Peptide IDs"),
       which(ColumnsTbl$Class == "Evidences count"),
       which(ColumnsTbl$Class == "Evidence IDs"),
       which(ColumnsTbl$Class == "Spectral count"),
       which(ColumnsTbl$Class == "Biotin peptides count"),
       which(ColumnsTbl$Class == "Biotin peptide IDs"),
       which(ColumnsTbl$Class == "Biotin evidences count"),
       which(ColumnsTbl$Class == "Biotin evidence IDs"),
       which(ColumnsTbl$Class == "Annotations"))
stopifnot(length(w) == nrow(ColumnsTbl))
ColumnsTbl <- ColumnsTbl[w,]
ColumnsTbl$Hide <- ColumnsTbl$Class %in% c("Peptide IDs", "Peptides count", "Evidence IDs", "Evidences count", "Spectral count", "Spectrum IDs",
                                           "Biotin peptides count", "Biotin peptide IDs", "Biotin evidences count", "Biotin evidence IDs",
                                           "Annotations")
if (length(intCols) > 1) {
  for (nm in names(intCols)[1:(length(intCols) - 1)]) {
    ColumnsTbl$Hide[which(ColumnsTbl$Class == intNms(nm, TRUE))] <- TRUE
  }
}
if (MakeRatios) {
  if (length(ratCols) > 1) {
    for (nm in names(ratCols)[1:(length(ratCols) - 1)]) {
      ColumnsTbl$Hide[which(ColumnsTbl$Class == ratNms(nm, TRUE))] <- TRUE
    }
  }
}
#
if (MakeRatios) { a <- KolEdit(ColumnsTbl$Col, intColsTbl, ratColsTbl) } else { a <- KolEdit(ColumnsTbl$Col, intColsTbl) }
ColumnsTbl$edit_Col <- unlist(a)
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/Write_Excel_end_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#xl_open(repFl)

# Save special quantitative table for proteins of interest
if ((length(Exp) > 1)&&(prot.list.Cond)) {
  dir <- paste0(wd, "/Proteins of interest")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  w <- grsep2(prot.list, PG$"Leading protein IDs")
  if (length(w)) {
    kols <- grep(topattern(PG.int.col), colnames(PG), value = TRUE)
    temp <- PG[w, c("Label", "Peptide IDs", kols)]
    for (kol in kols) {
      kol2 <- gsub(" log10\\(", " ", gsub("\\) - ", " - ", kol))
      temp[[kol2]] <- suppressWarnings(10^temp[[kol]])
      temp[which(!is.all.good(temp[[kol2]], 2)), kol] <- NA
      temp[[kol]] <- NULL
    }
    data.table::fwrite(temp, paste0(dir, "/Protein of interest profiles.tsv"),
                       sep = "\t", row.names = FALSE, na = "NA")
  }
}

#### Code chunk - Heatmap(s) of peptides mapping to proteins of interest
if ((length(Exp) > 1)&&(!is.null(prot.list))&&(length(prot.list))) {
  dir <- paste0(wd, "/Heatmaps")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  ref <- int.cols[which(names(int.cols) == "Imputed")-1]
  kol <- paste0(ref, " - ", Exp)
  kol <- kol[which(kol %in% colnames(pep))]
  StdWdth <- 6
  for (plp in prot.list) { #plp <- prot.list[1]
    Plp <- paste(db[which(db$"Protein ID" == plp), c("Common Name", "Protein ID")], collapse = " - ")
    grs <- grsep2(plp, pep$Proteins)
    if (length(grs)) {
      temp <- pep[grs, c("Sequence", "Modified sequence", kol)]
      Seq <- db$Sequence[match(plp, db$`Protein ID`)]
      temp$tst1 <- vapply(temp$Sequence, function(x) { nchar(unlist(strsplit(Seq, x))[1]) }, 1)
      temp$tst2 <- vapply(temp$Sequence, nchar, 1)
      temp$tst3 <- vapply(temp$"Modified sequence", nchar, 1)
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
        tst <- apply(temp[, Exp, drop = FALSE], 1, function(x) { mean(x[which(x > 0)]) })
        temp[, Exp] <- sweep(temp[, Exp, drop = FALSE], 1, tst, "/")
        temp2 <- set_colnames(reshape::melt.data.frame(temp, id.vars = "Modified sequence"),
                              c("Modified sequence", "Sample", "value"))
        temp2$"Modified sequence" <- gsub("^_|_$", "", temp2$"Modified sequence")
        temp2$Sample <- as.character(temp2$Sample)
        temp2$value <- suppressWarnings(log2(temp2$value))+StdWdth/2
        w <- which(!is.all.good(temp2$value, 2))
        temp2$value[w] <- NA
        temp2$value[which(temp2$value < -StdWdth/2)] <- -StdWdth/2
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
        temp2b <- data.frame(Xmin = c(Xscale/2 - 8*XScale2*((-Splits/2):(Splits/2)), 4),
                             Ymin = Yscale+10)
        temp2b$Xmax <- temp2b$Xmin + XScale2*8
        temp2b$Xmax[nrow(temp2b)] <- temp2b$Xmin[nrow(temp2b)]+0.5
        temp2b$value <- c((Splits:0)*StdWdth/Splits - StdWdth/2, NA)
        temp2b$Label <- ""
        temp2b$Label[c(1, Splits/2+1, Splits+1, Splits+2)] <- c(paste0("log10(max)", c("", -StdWdth/2, -StdWdth)), "0 or NA")
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
                plot.margin = margin(0, 0, 0, 0, "cm")) +
          coord_fixed(0.5) +
          scale_fill_gradient2(low = "#D55E00", mid = "black", high = "#009E73", na.value = "lightblue", breaks = c(-StdWdth/2, 0, StdWdth/2)) +
          xlab(NULL) + ylab(NULL) + theme(legend.position = "none") +
          xlim(Xlim[1], Xlim[2]) + ylim(Ylim[1], Ylim[2])
        poplot(heatmap.plot)
        ggsave(paste0(dir, "/", ttl, ".jpeg"), heatmap.plot, dpi = 600, width = 20, height = 12, units = "in")
        ggsave(paste0(dir, "/", ttl, ".pdf"), heatmap.plot, dpi = 600, width = 20, height = 12, units = "in")
        #system(paste0("open \"", dir, "/", ttl, ".jpeg", "\""))
        #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
      }
    }
  }
  setwd(wd)
}

#### Code chunk - peptide tables for visualizing the coverage of proteins of interest in 3D using SCV
kPBD <- kAlpha <- 0
if ((!is.null(prot.list))&&(length(prot.list))) {
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
        modFls <- paste0(MQFold, "/bin/conf/modifications", c("", ".local"), ".xml")
        modFls <- modFls[which(file.exists(modFls))]
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
        Modifs$Composition <- modFls$Composition[match(Modifs$`Full name`, modFls$Name)]
        Modifs$"Mass shift" <- sapply(strsplit(Modifs$Composition, " "), function(x) {
          #x <- strsplit(modifs$Formula, " ")[1]
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
    }
  }
  kount <- 0
  for (plp in prot.list) { #plp <- prot.list[1]
    grs <- grsep2(plp, pep$Proteins)
    if (length(grs)) {
      dir <- paste0(wd, "/Coverage/", plp)
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      # For each protein we want to download:
      # - all models for all fragments.
      # - latest version only!
      # PDB: we get PDB IDs from parsing the txt file
      if ("PDB" %in% colnames(db)) {
        tmp <- unlist(unlist(strsplit(db$PDB[match(plp, db$`Protein ID`)], ";")))
        tmp <- tmp[which(tmp != "")]
        kPBD <- length(tmp)
        if (kPBD) {
          for (i in tmp) {
            url <- paste0("https://files.rcsb.org/download/", i, ".pdb")
            try(download.file(url, paste0(dir, "/PDB-", i, ".pdb")), silent = TRUE)
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
        kPBD <- kAlpha <- 0 # Re-initialize
        # We have found at least one model which can be used to visualize coverage
        # Let's write peptidoforms
        tmp <- pep$`Modified sequence`[grs]
        if (SCV_PTMs) {
          tmp <- gsub("_", "", tmp)
          tmp <- gsub("\\)", "]_",gsub("\\(", "_[", tmp))
          tmp <- strsplit(tmp, "_")
          tmp <- vapply(tmp, function(x) { #x <- tmp[1]
            x <- unlist(x)
            w <- grep("\\[.+\\]", x)
            x[w] <- paste0("[", round(Modifs$`Mass shift`[match(x[w], paste0("[", Modifs$Mark, "]"))], 0), "]")
            return(paste(x, collapse = ""))
          }, "")
        } else { tmp <- gsub("[^A-Z]", "", tmp) }
        write(tmp, paste0(dir, "/SCV - observed peptides.txt"))
        kount <- kount+1
      }
    }
  }
  # If this worked, we write a guide in the Coverage folder
  if (kount) {
    Guide <- c("Visualing protein coverage in 3D using SCV",
               "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
               "",
               "If, for protein accession of interest \"PROTEIN\", a 3D model is available, then the following files are created in its subfolder:",
               " - \"peptidoforms.txt\": sequences of identified peptides",
               " - For each fragment number i for which a structure is available (ranging from 1 to n), \"AF-PROTEIN-F1-model_v#.pdb\", where \"#\" is the latest valid version.",
               "",
               "To visualise the peptides onto the folded protein structure, navigate to https://scv.lab.gy and:",
               " - paste the peptides into the \"PSM/peptide list\" field",
               " - load the pdb ",
               "")
    write(Guide, paste0(wd, "/Coverage/SCV - how to visualise protein coverage in 3D.txt"))
  }
}

# Save parameters
temp <- data.frame(Parameter = names(AnalysisParam),
                   Value = vapply(AnalysisParam, function(x) { paste(unlist(x), collapse = ";") }, ""),
                   row.names = NULL)
write.csv(temp, "Analysis parameters.csv", row.names = FALSE)
saveFun(AnalysisParam, "AnalysisParam.RData")
#
# Save copy of this script to local work directory
fs::file_copy(ScriptPath, wd, overwrite = TRUE)

Script <- readLines(ScriptPath)
gc()
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
saveImgFun(BckUpFl)
#loadFun(BckUpFl)

# Write PTMs table
temp <- Modifs
w <- which(vapply(colnames(Modifs), function(x) { "list" %in% class(Modifs[[x]]) }, TRUE))
for (i in w) { temp[[i]] <- vapply(temp[[i]],  paste, "", collapse = ", ") }
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
write.csv(temp, paste0(dir, "/Modifications.csv"), row.names = FALSE)

# Write SDRF file in case you want to submit to PRIDE
Src <- paste0(libPath, "/extdata/R scripts/Sources/SDRF_4_PRIDE.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

# Finalize analysis
Src <- paste0(libPath, "/extdata/R scripts/Sources/Finalize_analysis.R")
#rstudioapi::documentOpen(Src)
#loadFun(BckUpFl)
source(Src, local = FALSE)

### That's it, done!
#openwd(outDir)
#rm(list = ls())
