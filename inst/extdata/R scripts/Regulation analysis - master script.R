# Note:
# This master script hasn't been used in a long time and probably isn't too reliable.
# In particular, it seems heavily affected by the issues affecting the order of calls ("leaks", as I call them) which I have observed when shiny and/or svDialogs are being used.
# Proceed cautiously...

#### Code chunk - Initialization
if (!interactive()) { stop("This script should only be run within an interactive R session!") }
options(stringsAsFactors = FALSE)
options(install.packages.compile.from.source = "never")
#rm(list = ls()[which(!ls() %in% c("dtstNm", "wd", "indir", "outdir"))])

## The proteoCraft package can be re-installed at any time in the workflow (there is a specific script for this in the package's library folder),
## or just load it here:
require(proteoCraft)

homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
fls <- paste0(homePath, "/", c("Regulation analysis - master script.R",
                               "Regulation analysis - detailed script.R",
                               "No replicates analysis - detailed script.R",
                               "Default_locations.xlsx",
                               "LC_columns.xlsx"))
tst <- sum(!file.exists(fls))
if (tst) { proteoCraft::Configure() }

dirlist %<o% c()

# Create behind-the-scenes environment
BehindTheScenes %<o% new.env()

### Packages
## CRAN packages:
cran_req %<o% "rstudioapi"
bioc_req %<o% c()
cran_req <- unique(c(cran_req, "gtools", "svDialogs", "shiny", "shinyFiles", "rstudioapi", "TeachingDemos", "openxlsx", "tibble"))
if ((as.numeric(R.Version()$major) >= 3)||((as.numeric(R.Version()$major) == 3)&&(as.numeric(unlist(strsplit(R.Version()$minor, "\\."))[1]) >= 6))) {
  R.older.than.3.6 %<o% TRUE
} else { R.older.than.3.6 %<o% FALSE }
for (i in cran_req) { if (!require(i, character.only = TRUE)) { install.packages(i) } }
## Bioconductor packages:
biocInstall %<o% function(pack, load = TRUE, newR = R.older.than.3.6, recursive = TRUE) {
  if (newR) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
      #BiocManager::install()
    }
    if (!require(pack, character.only = TRUE)) { BiocManager::install(pack, update = FALSE) }
  } else {
    source("https://bioconductor.org/biocLite.R")
    biocLite(suppressUpdates = TRUE)
    if (!require(pack, character.only = TRUE)) { biocLite(pack, suppressUpdates = TRUE) }
  }
  if (load) { library(pack, character.only = TRUE) }
}
for (pack in bioc_req) { biocInstall(pack) }
# Now load packages
for (pack in c(cran_req, bioc_req)) { library(pack, character.only = TRUE) }


Rpath %<o% as.data.frame(library()$results)
Rpath <- normalizePath(Rpath$LibPath[match("proteoCraft", Rpath$Package)], winslash = "/")
fl <- paste0(homePath, "/Default_locations.xlsx")
flTst <- file.exists(fl)
if (!flTst) { file.copy(paste0(Rpath, "/proteoCraft/extdata/Default_locations.xlsx"), fl) }
inRoot <- read.xlsx(fl)
if (is.na(inRoot$Path[match("Temporary folder", inRoot$Folder)])) {
  inRoot$Path[match("Temporary folder", inRoot$Folder)] <- homePath
  wb <- createWorkbook()
  addWorksheet(wb, "Default_folders")
  writeData(wb, "Default_folders", inRoot, 1, 1)
  setColWidths(wb, "Default_folders", match(c("Folder", "Path", "Help"), colnames(inRoot)), c(25, 60, 150))
  HdStl <- createStyle(textDecoration = c("bold", "underline"))
  PthStl <- createStyle(fontName = "Consolas", textDecoration = "italic")
  HlpStl <- createStyle(textDecoration = "italic")
  addStyle(wb, "Default_folders", HdStl, 1, 1:ncol(inRoot), stack = TRUE)
  addStyle(wb, "Default_folders", PthStl, 1+(1:nrow(inRoot)), match("Path", colnames(inRoot)), stack = TRUE)
  addStyle(wb, "Default_folders", HlpStl, 1+(1:nrow(inRoot)), match("Help", colnames(inRoot)), stack = TRUE)
  saveWorkbook(wb, fl, overwrite = TRUE)
  #system(paste0("open \"", fl, "\""))
}
m <- match(c("Search folder", "Results delivery folder", "Archive folder"), inRoot$Folder)
inRoot2 <- setNames(inRoot$Path[m], paste0(inRoot$Folder[m], " (", inRoot$Path[m], ")"))
inRoot2 <- c(inRoot2, getVolumes()()) # this makes the directory at the base of your computer.
if (exists("indir")) {
  tst <- sapply(inRoot2, function(rt) { grepl(topattern(rt), indir) })
  inRoot2 <- c(inRoot2[which(tst)][order(nchar(inRoot2[which(tst)]), decreasing = TRUE)],
               inRoot2[which(!tst)])
}
if (!exists("WhoAmI")) { WhoAmI <- Sys.getenv("USERNAME") }
usesRep %<o% TRUE # Should be used to check if the custom script is compatible, see below!!!

outRoot <- inRoot2[c(grep("^Results delivery folder ", names(inRoot2)),
                     grep("^Archive folder ", names(inRoot2)),
                     grep("^Results delivery folder |^Archive folder ", names(inRoot2), invert = TRUE))]
scrpts <- setNames(c("Regulation analysis - detailed script.R", "No replicates analysis - detailed script.R", NA), c("Replicates", "No replicates", "Custom"))
scrptPaths <- setNames(paste0(homePath, "/", scrpts[which(!is.na(scrpts))]), names(scrpts)[which(!is.na(scrpts))])
WorkFlows %<o% setNames(lapply(scrptPaths, function(pth) { #pth <- scrptPaths[1]
  x <- readLines(pth)
  tmp <- gsub("^###-\\|-### *Workflows: *", "", grep("^###-\\|-### *Workflows: *", x, value = TRUE))
  eval(parse(text = tmp))
}), names(scrptPaths))
ScriptPath %<o% scrptPaths["Replicates"]
if (!exists("WorkFlow")) { WorkFlow <- WorkFlows["Replicates"][1] }

if ((!exists("indir"))||(!dir.exists(indir))) { indir <- inRoot2[grep("^Search folder ", names(inRoot2))[1]] }
if ((!exists("outdir"))||(!dir.exists(outdir))) { outdir <- inRoot2[grep("^Results delivery folder ", names(inRoot2))[1]] }
if (!exists("dtstNm")) { dtstNm <- "My_latest_proteomics_dataset" }
indir %<o% indir
outdir %<o% outdir
dtstNm %<o% dtstNm
WhoAmI %<o% WhoAmI
WorkFlow %<o% WorkFlow
append <- FALSE
defRt <- inRoot2[match(gsub("/.*", "/", inRoot2[grep("^Search folder ", names(inRoot2))]), inRoot2)]
tmp <- c("MaxQuant", "DiaNN", "FragPipe", "Proteome Discoverer")
SearchSoftware %<o% setNames(gsub(" ", "", toupper(tmp)), tmp)
ui <- shinyUI(fluidPage(tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
                        shinyjs::useShinyjs(),
                        titlePanel(tag("u", "Start analysis"),
                                   "Start analysis"),
                        mainPanel(
                          sidebarLayout(
                            sidebarPanel(
                              textInput("Who", "Enter your name?", WhoAmI, "100%"),
                              h5(em("For traceability only: your name will be written into a tab of the output Excel table created by this workflow, so as to keep a trace of who ran their analysis.")),
                              br(),
                              textInput("dtstNm", "Enter a name for this project or analysis:", dtstNm, "100%"),
                              br(),
                              fluidRow(column(12,
                                              selectInput("Script", "Select analysis script...", names(scrpts)[1:2], names(scrpts)[1], width = "100%"),
                                              em("or..."),
                                              shinyFilesButton("CustScriptFl", em("optionally select a compatible custom data analysis script"), "", FALSE),
                                              strong("Selected script:"), br(),
                                              uiOutput("ScrptMsg"),
                                              br(),
                                              uiOutput("Workflows"),
                                              br()
                              )),
                              br(),
                              br()
                            ),
                            mainPanel(
                              br(),
                              #actionButton("indir", "Select input directory"),
                              HTML("Select an input folder."), br(),
                              shinyDirButton("indir", "Input MS search directory", "Select input directory"),
                              HTML("This should contain the output from a DiaNN, FragPipe, MaxQuant or Proteome Discoverer search."), br(),
                              strong("Input directory = "), span(uiOutput("Indir"), style = "color:blue", .noWS = "outside"),
                              br(),
                              uiOutput("inputType"),
                              br(), br(),
                              HTML("Select an output folder, where the results will be copied from the temporary directory once processing finishes."),
                              br(),
                              #actionButton("outdir", "Select input directory"),
                              fluidRow(column(12,
                                              shinyDirButton("outdir", "Output directory", ""), br(),
                                              strong("Output directory = "), span(uiOutput("Outdir"), style = "color:blue", .noWS = "outside")),
                                       column(4, checkboxInput("append", "Append dataset name?", FALSE))),
                              br()
                            )),
                          br(),
                          em("Once you are finished, click "), actionButton("saveBtn", "Save"), em(" to continue."), br(),
                          em("Saving is only possible once:"), br(),
                          em(" - a valid dataset name is entered,"), br(),
                          em(" - valid in- and output directories are selected."),
                          br()
                        )))
server <- shinyServer(function(input, output, session) {
  outdirRoot <- outdir
  while ((!dir.exists(outdirRoot))&&(grepl("/", outdirRoot))) { outdirRoot <- gsub("/[^/]+$", "", outdirRoot) }
  #
  # Reactive values
  WHO <- reactiveVal(WhoAmI)
  SCRIPT <- reactiveVal("Replicates")
  SCRIPTfl <- reactiveVal(scrptPaths["Replicates"])
  SCRIPTpar <- reactiveVal("")
  USECUST <- reactiveVal(FALSE)
  WORKFLOWS <- reactiveVal(setNames(WorkFlows[["Replicates"]], NULL))
  WORKFLOW <- reactiveVal(setNames(WorkFlow, NULL))
  INDIR <- reactiveVal(indir)
  INDIRtype <- reactiveVal("diaNN")
  OUTDIRRoot <- reactiveVal(outdirRoot)
  OUTDIR <- reactiveVal(outdir)
  DATASETNAME <- reactiveVal(dtstNm)
  APPEND <- reactiveVal(append)
  #
  # Reactive functions to update UI
  prtInDr <- function(reactive = TRUE) {
    if (reactive) { msg = INDIR() } else { msg <- indir }
    renderUI({ em(msg) })
  }
  prtOutDr <- function(reactive = TRUE) {
    if (reactive) {
      if (APPEND()) {
        dtstNm2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", DATASETNAME())
        if (nchar(dtstNm2)) { OUTDIR(paste0(OUTDIRRoot(), "/", dtstNm2)) } else { OUTDIR(OUTDIRRoot()) }
      } else { OUTDIR(OUTDIRRoot()) }
      renderUI({ em(OUTDIR()) })
    } else {
      if (append) {
        dtstNm2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", dtstNm)
        if (nchar(dtstNm2)) { outdir <- (paste0(outdirRoot, "/", dtstNm2)) } else { outdir <- outdirRoot }
      } else { outdir <- outdirRoot }
      renderUI({ em(outdir) })
    }
  }
  updtOutDr <- function(reactive = TRUE) {
    if (reactive) {
      if (APPEND()) {
        dtstNm2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", DATASETNAME())
        if (nchar(dtstNm2)) { OUTDIR(paste0(OUTDIRRoot(), "/", dtstNm2)) } else { OUTDIR(OUTDIRRoot()) }
      } else { OUTDIR(OUTDIRRoot()) }
    } else {
      if (append) {
        dtstNm2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", dtstNm)
        if (nchar(dtstNm2)) { outdir <- (paste0(outdirRoot, "/", dtstNm2)) } else { outdir <- outdirRoot }
      } else { outdir <- outdirRoot }
    }
  }
  updtWrkfls <- function() {
    tmp <- grep("^Workflows: *", SCRIPTpar(), value = TRUE)
    if (length(tmp)) {
      tmp <- gsub("^Workflows: *", "", tmp)
      tmp <- setNames(eval(parse(text = tmp)), NULL)
      WORKFLOWS(tmp)
      renderUI({
        lst <- vector("list", 1)
        lst[[1]] <- list(selectInput("Workflow", "Select data analysis workflow", WORKFLOWS(), WORKFLOWS()[1], width = "100%"))
        return(lst)
      })
    } else { renderUI({ em("") }) }
  }
  updtScrptMsg <- function() { # This function as written can only be used in a reactive context!
    pth <- SCRIPTfl()
    if (file.exists(pth)) {
      tmp <- SCRIPTpar()
      tmp <- as.logical(gsub("^Replicates\\? *", "", grep("^Replicates\\? *", tmp, value = TRUE)))
      msg <- list(span(HTML(paste0("   ", pth)), style = "color:blue", .noWS = "outside"),
                  br(),
                  span(em(c("This script does not support multiple-replicates and compares individual samples without replicates-dependent statistical tests.",
                            "This script supports statistical tests on sample groups as long as each pair of groups includes 2 replicates or more.")[tmp+1]),
                       style = "color:orange", .noWS = "outside"))
    } else { msg <- HTML("") }
    renderUI(msg)
  }
  updtType <- function(reactive = TRUE) {
    if ((exists("SearchSoft"))&&(!is.na(SearchSoft))) { dflt <- names(SearchSoft) } else {
      if (reactive) { fls <- list.files(INDIR()) } else { fls <- list.files(indir) }
      tsv2txt <- length(grep("\\.tsv$", fls))/length(grep("\\.txt$", fls))
      dflt <- c(FALSE,
                length(grep("\\.log\\.txt$", fls)) > 0,
                length(grep("\\.fp-manifest$", fls)) > 0,
                length(grep("\\.xlsx$", fls)) > 0)
      if (!sum(dflt)) { dflt[1] <- TRUE }
      if (sum(dflt) > 1) {
        w <- which(dflt)
        dflt[w[2:length(w)]] <- FALSE
      }
      dflt <- names(SearchSoftware)[which(dflt)]
    }
    lst <- vector("list", 1)
    lst[[1]] <- list(selectInput("SearchSoft", "Select search software used", names(SearchSoftware), dflt, width = "200px"))
    return(renderUI(lst))
  }
  #
  # Initial UI renders
  output$Indir <- prtInDr(reactive = FALSE)
  output$Outdir <- prtOutDr(reactive = FALSE)
  output$ScrptMsg <- renderUI(em(""))
  output$Workflows <- renderUI({
    lst <- vector("list", 1)
    lst[[1]] <- list(selectInput("Workflow", "Select data analysis workflow", setNames(WorkFlows[["Replicate"]], NULL), setNames(WorkFlow, NULL), width = "100%"))
    return(lst)
  })
  output$inputType <- updtType(reactive = FALSE)
  #
  # Observers
  observeEvent(input$Who, { WHO(input$Who) })
  observeEvent(input$Script, {
    SCRIPT(input$Script)
    USECUST(FALSE)
    SCRIPTfl(scrptPaths[SCRIPT()])
    tmp <- readLines(SCRIPTfl())
    tmp <- gsub("^###-\\|-### *", "", grep("^###-\\|-### *", tmp, value = TRUE))
    SCRIPTpar(tmp)
    output$ScrptMsg <- updtScrptMsg()
    output$Workflows <- updtWrkfls()
    #}
  })
  observe({
    shinyFileChoose(input, "CustScriptFl", roots = getVolumes(), filetypes = "R"
                    # Note: do not use default root and paths until they start behaving properly!!!
    )
    {
      tmp <<- input$CustScriptFl
      print(class(tmp))
      if ((!is.null(tmp))&&(is.list(tmp))) {
        USECUST(TRUE)
        tmp <- parseFilePaths(getVolumes(), tmp)$datapath
        SCRIPTfl(tmp)
        tmp2 <- readLines(SCRIPTfl())
        tmp2 <- gsub("^###-\\|-### *", "", grep("^###-\\|-### *", tmp2, value = TRUE))
        SCRIPTpar(tmp2)
        output$ScrptMsg <- updtScrptMsg()
        output$Workflows <- updtWrkfls()
      } else {
        USECUST(FALSE)
      }
    }
  })
  observeEvent(input$Workflow, { WORKFLOW(input$Workflow) })
  observe({
    shinyDirChoose(input, "indir", roots = inRoot2,
                   # Note: do not use default root and paths until they start behaving properly!!!
    )
    {
      tmp <- input$indir
      if ("list" %in% class(tmp)) {
        INDIR(paste0(inRoot2[tmp$root], paste(tmp$path, collapse = "/")))
        output$Indir <- prtInDr()
        if ((dir.exists(OUTDIRRoot()))&&(nchar(DATASETNAME()))) { shinyjs::enable("saveBtn") }
        output$inputType <- updtType()
      } else { shinyjs::disable("saveBtn") }
    }
  })
  observeEvent(input$SearchSoft, { INDIRtype(input$SearchSoft) })
  observe({
    shinyDirChoose(input, "outdir", roots = outRoot)
    {
      tmp <- input$outdir
      if ("list" %in% class(tmp)) {
        OUTDIRRoot(paste0(outRoot[tmp$root], paste(tmp$path, collapse = "/")))
        OUTDIR(OUTDIRRoot())
        updtOutDr()
        output$Outdir <- prtOutDr()
        if ((dir.exists(INDIR()))&&(nchar(DATASETNAME()))) { shinyjs::enable("saveBtn") }
      } else { shinyjs::disable("saveBtn") }
    }
  })
  observeEvent(input$dtstNm, {
    if (nchar(input$dtstNm)) {
      DATASETNAME(input$dtstNm)
      updtOutDr()
      output$Outdir <- prtOutDr()
      if ((dir.exists(INDIR()))&&(dir.exists(OUTDIRRoot()))) { shinyjs::enable("saveBtn") }
    } else { shinyjs::disable("saveBtn") }
  })
  observeEvent(input$append, {
    APPEND(input$append)
    updtOutDr()
    output$Outdir <- prtOutDr()
  }, ignoreInit = TRUE)
  observeEvent(input$saveBtn, {
    WhoAmI <- WHO()
    if (WhoAmI == "Your name here") { WhoAmI <- NA }
    WhoAmI <<- WhoAmI
    indir <<- INDIR()
    outdir <<- OUTDIR()
    dtstNm <<- DATASETNAME()
    append <<- APPEND()
    WorkFlows <<- WORKFLOWS()
    WorkFlow <<- WORKFLOW()
    ScriptPath <<- SCRIPTfl()
    SearchSoft <<- SearchSoftware[INDIRtype()]
    stopApp()
  })
  session$onSessionEnded(function() { stopApp() })
})
print(shinyApp(ui, server, options = list(launch.browser = TRUE)))
SearchSoft %<o% SearchSoft
#
tst <- gsub("[^A-Z,a-z,0-9]", "", unlist(sapply(c(dtstNm, indir, outdir), strsplit, "/")))
tst <- sum(sapply(tst, char2seed, set = FALSE))
tst2 <- lapply(unlist(strsplit(WhoAmI, " +")), function(x) {
  x <- unlist(strsplit(x, ""))
  return(x[1:(min(c(2, length(x))))])
})
tst2 <- c(sapply(tst2, function(x) { x[[1]] }),
          sapply(tst2, function(x) { if (length(x) >= 2) { x[[2]] } else { NA } }))
tst2 <- tst2[which(!is.na(tst2))]
tst2 <- paste(tst2[1:(min(c(2, length(tst2))))], collapse = "")
tst3 <- gsub("[a-z, ]", "", dtstNm)
wd %<o% paste0(inRoot$Path[match("Temporary folder", inRoot$Folder)], "/", tst2, "_", tst3, "_", tst)
if (!dir.exists(wd)) { dir.create(wd, recursive = TRUE) }
setwd(wd)
#
cat(c(paste0("Script path:\n -> ", ScriptPath, "\n\n"), 
      paste0("Input directory:\n -> ", indir, "\n\n"),
      paste0("Final output directory:\n -> ", outdir, "\n(folder will be created at the end of the workflow if it doesn't exist yet)\n\n")),
    paste0("Temporary work directory:\n -> ", wd, "\n\n"),
    "The data is processed within the temporary work directory, which is meant to have a short path to avoid issues with saving too long paths.\n",
    "The script then attempts to copy the files to the final output directory.\n",
    "If this fails (e.g. because of too long paths), then you can always manually copy the analysis from the temporary folder.\n",
    sep = "")
#

#load(paste0(wd, "/Backup.RData"))
#openwd()

# Should we prompt the user for going forward after each chunk or just go through all without prompting? (unless specified in chunk itself)
BehindTheScenes$ProceedWithCaution <- FALSE
#BehindTheScenes$ProceedWithCaution <- c(FALSE, TRUE)[match(dlg_message("Do you want to proceed step-by-step?", "yesno")$res, c("no", "yes"))]

BehindTheScenes$Workflow <- WorkFlow
BehindTheScenes$ScriptFile <- ScriptPath

# (When updating the detailed script, reload it here:)
BehindTheScenes$Script <- suppressWarnings(readLines(BehindTheScenes$ScriptFile))
BehindTheScenes$Chunks <- data.frame(Start = 1,
                                     End = length(BehindTheScenes$Script))
# Identify chunks which follow an environment backup:
#g0 <- grep("^ *save\\.image\\(\"Backup\\.RData\"\\)$", BehindTheScenes$Script)
g1 <- grep("^#{4} Code chunk - ", BehindTheScenes$Script)
g1 <- g1[which(g1 > 1)]
if (length(g1)) {
  #g1 <- sapply(c(1, g0), function(x) { g1[which(g1 > x)][1] })
  #g0 <- sapply(g1, function(x) { rev(g0[which(g0 < x)])[1] })
  g1 <- unique(g1[which(!is.na(g1))])
  #g0 <- unique(g0[which(!is.na(g0))])
  #if (length(g1) > 1) {
  BehindTheScenes$Chunks <- data.frame(Start = c(1, g1),
                                       End = c(g1-1, length(BehindTheScenes$Script)))
  #}
}
g1 <- c(1, g1)
tst <- apply(BehindTheScenes$Chunks[, c("Start", "End")], 1, function(x) { x[[2]] - x[[1]] })
BehindTheScenes$Chunks$Text <- gsub("^#{4} Code chunk - ", "", BehindTheScenes$Script[g1])
BehindTheScenes$Chunks$Action <- c("run", "check")[BehindTheScenes$ProceedWithCaution+1]
g2 <- grep(" //", BehindTheScenes$Chunks$Text)
BehindTheScenes$Chunks$Action[g2] <- gsub(".* //", "", BehindTheScenes$Chunks$Text[g2])
BehindTheScenes$Chunks$Text[g2] <- gsub(" //.*", "", BehindTheScenes$Chunks$Text[g2])

#load(paste0(wd, "/Backup.RData"))
BehindTheScenes$Steps <- paste0(" - ", 1:length(BehindTheScenes$Chunks$Text), ": ", BehindTheScenes$Chunks$Text, ", action = ", BehindTheScenes$Chunks$Action)
BehindTheScenes$Chunks$Action[which(BehindTheScenes$Chunks$Action %in% c("run", "check"))] <- c("run", "check")[BehindTheScenes$ProceedWithCaution+1]
if ((!"Pos" %in% names(BehindTheScenes))||(is.na(BehindTheScenes$Pos))||(BehindTheScenes$Pos > nrow(BehindTheScenes$Chunks))) { BehindTheScenes$Pos <- 1 }
if (file.exists(paste0(wd, "/Backup.RData"))) {
  if (file.exists(paste0(wd, "/Last completed step.txt"))) {
    msg <- "An analysis has been performed in this directory before. Do you want to..."
    opt <- c("Start from scratch?                                                                                                                                                                                                     ",
             "Perform partial processing?                                                                                                                                                                                             ")
  } else {
    msg <- "An analysis has been performed in this directory before, but no record of the last completed step was found (you probably ran the detailed script manually). Should we..."
    opt <- c("Start from scratch?                                                                                                                                                                                                                                                      ",
             "Perform partial processing anyway? (you know where to start from)                                                                                                                                                                                                        ")
  }
  BehindTheScenes$Partial_processing <- c(FALSE, TRUE)[match(dlg_list(opt, opt[1], title = msg, "yesno")$res, opt)]
  if (BehindTheScenes$Partial_processing) {
    #dflt <- 1 # For when executing line-by-line outside of loop
    if (!file.exists(paste0(wd, "/Last completed step.txt"))) { dflt <- 1 } else {
      dflt <- as.integer(gsub("^ -> Last completed step: ", "", readLines(paste0(wd, "/Last completed step.txt"))))+1
      dflt <- min(c(dflt, length(BehindTheScenes$Steps)))
    }
    msg <- "Which step do you want to process from?"
    opt <- sapply(BehindTheScenes$Steps, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") })
    BehindTheScenes$Pos <- match(dlg_list(opt, opt[dflt], title = msg, "yesno")$res, opt)
    dir <- paste0(wd, "/Workflow control")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    save(BehindTheScenes, file = paste0(dir, "/BehindTheScenes_env.RData"))
    load(paste0(wd, "/Backup.RData"))
    load(paste0(dir, "/BehindTheScenes_env.RData"))
    # Check that packages are loaded
    for (pack in cran_req) { if (!require(pack, character.only = TRUE)) { install.packages(pack) } }
    for (pack in bioc_req) { biocInstall(pack) }
    for (pack in c(cran_req, bioc_req)) { require(pack, character.only = TRUE) }
    # It is good to reload parameters by default here, as I often update them when stopping the analysis
    suppressWarnings(try(Param <- Param.load(), silent = TRUE))
  }
}
cat(paste0("Starting from: ", gsub("^ - ", "", BehindTheScenes$Steps[BehindTheScenes$Pos]), "\n"))

#load("All_decisions.RData")
#m <- which(AllAnsw$Parameter == "Impute")
#AllAnsw$Value[[m]] <- FALSE
#
#save("AllAnsw", file = "All_decisions.RData")
#load("All_decisions.RData")
#AllAnsw$Value[[m]]

#BehindTheScenes$Pos <- BehindTheScenes$Pos + 1
#BehindTheScenes$Pos <- BehindTheScenes$Pos - 1
#BehindTheScenes$Chunks$Action[BehindTheScenes$Pos] <- "run" # To restart, if checking each chunk and having stopped
#BehindTheScenes$Chunks$Action[which(BehindTheScenes$Chunks$Text == "Assemble protein groups")] <- "check"
#Param <- Param.load()
BehindTheScenes$Code <- list()
BehindTheScenes$OneChunkToRuleThemAll <- TRUE
with(BehindTheScenes, {
  Debug <- FALSE
  #Debug <- TRUE
  if (Debug) { for (nm in names(BehindTheScenes)) { assign(nm, BehindTheScenes[[nm]]) } }
  if (OneChunkToRuleThemAll) {
    if (Chunks$Action[Pos] == "ask") {
      Pos_Mess <- paste0("Code chunk #", Pos, ": \"", Chunks$Text[Pos], "\" -> action?")
      opt <- c("run                                                                                                                                                                                                                                                                  ",
               "stop workflow                                                                                                                                                                                                                                                        ")
      tmp <- dlg_list(opt, opt[1], title = Pos_Mess)$res
      if (!length(tmp)) { tmp <- opt[3] }
      Chunks$Action[Pos] <- gsub(" .+", "", tmp)
    }
    if (Chunks$Action[Pos] == "check") {
      Pos_Mess <- paste0("Run code chunk #", Pos, ": \"", Chunks$Text[Pos], "\"?")
      Chunks$Action[Pos] <- c("run", "stop")[match(dlg_message(Pos_Mess, "yesno")$res, c("yes", "no"))]
    }
    if (Chunks$Action[Pos] == "run") {
      if (Pos != 1) { cat(paste0("Executing script, starting at step ", Chunks$Text[[Pos]], "\n")) }
      Code <- Script[Chunks$Start[Pos]:length(Script)]
      write(Code, "temp_chunk.R")
      Sys.sleep(1)
      source("temp_chunk.R")
      Sys.sleep(1)
      unlink("temp_chunk.R")
      write(paste0(" -> Last completed step: ", nrow(Chunks)), paste0(wd, "/Last completed step.txt"))
    }
    if (Chunks$Action[Pos] %in% c("run", "skip")) { Pos <- Pos+1 }
  } else {
    #BehindTheScenes$Steps[BehindTheScenes$Pos]
    while ((Pos <= nrow(Chunks))&&(Chunks$Action[Pos] %in% c("check", "run", "skip", "ask"))) {
      if (Chunks$Action[Pos] == "ask") {
        Pos_Mess <- paste0("Code chunk #", Pos, ": \"", Chunks$Text[Pos], "\" -> action?")
        opt <- c("run                                                                                                                                                                                                                                                                  ",
                 "skip                                                                                                                                                                                                                                                                 ",
                 "stop workflow                                                                                                                                                                                                                                                        ")
        tmp <- dlg_list(opt, opt[1], title = Pos_Mess)$res
        if (!length(tmp)) { tmp <- opt[3] }
        Chunks$Action[Pos] <- gsub(" .+", "", tmp)
      }
      if (Chunks$Action[Pos] == "check") {
        Pos_Mess <- paste0("Run code chunk #", Pos, ": \"", Chunks$Text[Pos], "\"?")
        Chunks$Action[Pos] <- c("run", "stop")[match(dlg_message(Pos_Mess, "yesno")$res, c("yes", "no"))]
      }
      if (Chunks$Action[Pos] == "run") {
        cat(paste0(Chunks$Text[[Pos]], "\n"))
        Code[[Pos]] <- Script[Chunks$Start[Pos]:Chunks$End[Pos]]
        #Code <- BehindTheScenes$Script[BehindTheScenes$Chunks$Start[BehindTheScenes$Pos]:BehindTheScenes$Chunks$End[BehindTheScenes$Pos]]
        #if (grepl("^ +save\\.image\\(\"Backup\\.RData\"\\)", Code[1])) {
        #  Code <- Code[(grep("^\\}", Code)[1]+1):length(Code)]
        #}
        #while ((Code[1] == "")||(grepl("^save\\.image\\(\"Backup\\.RData\"\\)", Code[1]))) { Code <- Code[2:length(Code)] }
        #eval(parse(text = Code), .GlobalEnv)
        write(Code[[Pos]], "temp_chunk.R")
        Sys.sleep(1)
        tst <- try(source("temp_chunk.R"), silent = TRUE)
        if (!"try-error" %in% class(tst)) {
          Sys.sleep(1)
          unlink("temp_chunk.R")
          write(paste0(" -> Last completed step: ", Steps[Pos]), paste0(wd, "/Last completed step.txt"))          
          Pos <- Pos+1
        } else { stop() }
        rm(tst)
      }
    }
    if (Chunks$Action[Pos] == "skip") {
      Pos_Mess <- paste0("Skipping code chunk #", Pos, ": \"", Chunks$Text[Pos], "\"")
      cat(paste0(Pos_Mess, "\n"))
      Pos <- Pos+1
    }
  }
})
if (file.exists("temp_chunk.R")) { unlink("temp_chunk.R") }
# Useful command lines
#View(BehindTheScenes$Chunks) #-> View chunks table
#system(paste0("open \"", ScriptFile, "\"")) #-> open detailed script
#closeAllConnections() #-> close all currently connections, especially the log file
#Param <- Param.load(filter.deprecated = TRUE) #-> reload parameters file
