######################
### Start analysis ###
######################

if (!exists("renv")) { renv %<o% FALSE }

# dflt <- "My_latest_proteomics_dataset"
# if (exists("dtstNm")) { dflt <- dtstNm }
# dtstNm %<o% svDialogs::dlg_input("Enter a name for this project or analysis:", dflt)$res
if (!exists("dtstNm")) { dtstNm <- "My_latest_proteomics_dataset" }
dtstNm %<o% dtstNm

# From now on we will set seed to avoid results reproducibility issues
if (!exists("mySeed")) { mySeed <- 1234567 }
# Ultimately this should be made a user-defined input in the App below and be written into dataset details
mySeed %<o% mySeed

# Define input, output, project folder etc...
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1, 1)) }
if (!exists("writeRaws")) { writeRaws <- TRUE }
if (!exists("writeSearch")) { writeSearch <- TRUE }

ObjNm <- "ProcessedByUs"
if (!exists(ObjNm)) {
  if ((scrptType == "withReps")&&(ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) {
    ProcessedByUs <- AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]]
  } else {
    ProcessedByUs <- TRUE
  }
}
ProcessedByUs %<o% as.logical(ProcessedByUs)
if (is.na(ProcessedByUs)) { ProcessedByUs <- TRUE }
if (scrptType == "withReps") {
  AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
  tmp <- AllAnsw[1,]
  tmp[, c("Parameter", "Message")] <- c(ObjNm, "No questions asked!")
  tmp$Value <- list(get(ObjNm))
  m <- match(ObjNm, AllAnsw$Parameter)
  if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
}

N.clust %<o% N.clust
writeRaws %<o% writeRaws
writeSearch %<o% writeSearch
ProcessedByUs %<o% ProcessedByUs
require(openxlsx)
require(shiny)
require(shinyjs)
require(shinyFiles)
require(TeachingDemos)

if (!RunByMaster) {
  fl <- paste0(homePath, "/Default_locations.xlsx")
  flTst <- file.exists(fl)
  if (!flTst) { file.copy(paste0(RPath, "/proteoCraft/extdata/Default_locations.xlsx"), fl) }
  inRoot <- openxlsx::read.xlsx(fl)
  if (is.na(inRoot$Path[match("Temporary folder", inRoot$Folder)])) {
    inRoot$Path[match("Temporary folder", inRoot$Folder)] <- homePath
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Default_folders")
    openxlsx::writeData(wb, "Default_folders", inRoot, 1, 1)
    openxlsx::setColWidths(wb, "Default_folders", match(c("Folder", "Path", "Help"), colnames(inRoot)), c(25, 60, 150))
    HdStl <- openxlsx::createStyle(textDecoration = c("bold", "underline"))
    PthStl <- openxlsx::createStyle(fontName = "Consolas", textDecoration = "italic")
    HlpStl <- openxlsx::createStyle(textDecoration = "italic")
    openxlsx::addStyle(wb, "Default_folders", HdStl, 1, 1:ncol(inRoot), stack = TRUE)
    openxlsx::addStyle(wb, "Default_folders", PthStl, 1+(1:nrow(inRoot)), match("Path", colnames(inRoot)), stack = TRUE)
    openxlsx::addStyle(wb, "Default_folders", HlpStl, 1+(1:nrow(inRoot)), match("Help", colnames(inRoot)), stack = TRUE)
    openxlsx::saveWorkbook(wb, fl, overwrite = TRUE)
    #system(paste0("open \"", fl, "\""))
  }
  m <- match(c("Search folder", "Results delivery folder", "Archive folder"), inRoot$Folder)
  inRoot2 <- setNames(inRoot$Path[m], paste0(inRoot$Folder[m], " (", inRoot$Path[m], ")"))
  inRoot2 <- c(inRoot2, shinyFiles::getVolumes()()) # this makes the directory at the base of your computer.
  if (exists("indir")) {
    tst <- sapply(inRoot2, function(rt) { grepl(proteoCraft::topattern(rt), indir) })
    inRoot2 <- c(inRoot2[which(tst)][order(nchar(inRoot2[which(tst)]), decreasing = TRUE)],
                 inRoot2[which(!tst)])
  }
  if (!exists("WhoAmI")) { WhoAmI <- Sys.getenv("USERNAME") }
  #
  outRoot <- inRoot2[c(grep("^Results delivery folder ", names(inRoot2)),
                       grep("^Archive folder ", names(inRoot2)),
                       grep("^Results delivery folder |^Archive folder ", names(inRoot2), invert = TRUE))]
  scrptPaths <- setNames(ScriptPath,
                         c("Replicates",
                           "No replicates")[match(scrptType, c("withReps", "noReps"))])
  WorkFlows %<o% eval(parse(text = gsub("^###-\\|-### *Workflows: *", "", grep("^###-\\|-### *Workflows: *", readLines(scrptPaths), value = TRUE))),
                      envir = .GlobalEnv)
  if ((!exists("WorkFlow"))||(is.null(WorkFlow))) {
    WorkFlow <- WorkFlows[1]
  } else {
    #if (scrptType == "withReps") { # Commented because this should apply to both, if incorrect change to if (scrptType == "noReps") {
    if ((!WorkFlow %in% WorkFlows)&&(WorkFlow %in% names(WorkFlows))) {
      WorkFlow <- WorkFlows[WorkFlow]
    } # Because I was stupid when I created Workflow and Workflows
    # The former matches a name, not a value, of the latter
    if (!WorkFlow %in% WorkFlows) { WorkFlow <- WorkFlows[1] }
    #}
  }
  #
  if ((!exists("indir"))||(!dir.exists(indir))) { indir <- inRoot2[grep("^Search folder ", names(inRoot2))[1]] }
  if ((!exists("outdir"))||(!dir.exists(outdir))) { outdir <- inRoot2[grep("^Results delivery folder ", names(inRoot2))[1]] }
  append <- FALSE
  defRt <- inRoot2[match(gsub("/.*", "/", inRoot2[grep("^Search folder ", names(inRoot2))]), inRoot2)]
  tmp <- c("MaxQuant", "DiaNN", "FragPipe", "Proteome Discoverer")
  SearchSoftware %<o% setNames(gsub(" ", "", toupper(tmp)), tmp)
  scriptPar <- readLines(ScriptPath)
  scriptPar <- gsub("^###-\\|-### *", "", grep("^###-\\|-### *", scriptPar, value = TRUE))
  appNm <- "Start analysis"
  dtstNm2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", dtstNm)
  ui <- shiny::shinyUI(shiny::fluidPage(shiny::tags$head(shiny::tags$style(shiny::HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
                                        shinyjs::useShinyjs(),
                                        shinyjs::extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
                                        shiny::titlePanel(shiny::tag("u", "Start analysis"),
                                                          appNm),
                                        shiny::mainPanel(
                                          shiny::sidebarLayout(
                                            shiny::sidebarPanel(
                                              shiny::textInput("Who", "Enter your name?", WhoAmI, "100%"),
                                              shiny::h5(shiny::em("For traceability only: your name will be written into a tab of the output Excel table created by this workflow, so as to keep a trace of who ran their analysis.")),
                                              shiny::br(),
                                              shiny::textInput("dtstNm", "Enter a name for this project or analysis:", dtstNm, "100%"),
                                              shiny::br(),
                                              shiny::uiOutput("Workflows"),
                                              shiny::numericInput("vCPUs", "Number of available vCPUs (threads):", N.clust, 1,  N.clust, 1),
                                              shiny::numericInput("Seed", "Set a seed for reproducible random processes:", mySeed, -Inf, Inf, 1),
                                              shiny::br(),
                                              shiny::br()
                                            ),
                                            shiny::mainPanel(
                                              shiny::br(),
                                              #shiny::actionButton("indir", "Select input directory"),
                                              shiny::HTML("Select an input folder."),
                                              shiny::br(),
                                              shinyFiles::shinyDirButton("indir", "Input MS search directory", "Select input directory"),
                                              shiny::HTML("This should contain the output from a DiaNN, FragPipe, MaxQuant or Proteome Discoverer search."),
                                              shiny::br(),
                                              shiny::strong("Input directory = "),
                                              shiny::span(shiny::uiOutput("Indir"), style = "color:blue", .noWS = "outside"),
                                              shiny::br(),
                                              shiny::uiOutput("inputType"),
                                              shiny::br(),
                                              shiny::br(),
                                              shiny::HTML("Select an output folder, where the results will be copied from the temporary directory once processing finishes."),
                                              shiny::br(),
                                              #shiny::actionButton("outdir", "Select input directory"),
                                              shiny::fluidRow(shiny::column(12,
                                                                            shinyFiles::shinyDirButton("outdir", "Output directory", ""),
                                                                            shiny::br(),
                                                                            shiny::strong("Output directory = "),
                                                                            shiny::span(shiny::uiOutput("Outdir"), style = "color:blue", .noWS = "outside")),
                                                              shiny::column(4,
                                                                            shiny::checkboxInput("append", "Append dataset name?", FALSE))),
                                              shiny::br(),
                                              shiny::column(4,
                                                            shiny::checkboxInput("writeRaws", "Copy raw files to delivery folder?", writeRaws)),
                                              shiny::column(4,
                                                            shiny::checkboxInput("writeSearch", "Copy search results to delivery folder?", writeSearch)),
                                              shiny::column(4,
                                                            shiny::checkboxInput("ProcessedByUs", "Attempt to write interactively a template sample prep text?", ProcessedByUs)),
                                              shiny::br())),
                                          shiny::br(),
                                          shiny::em("Once you are finished, click "),
                                          shinyWidgets::actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                                          shiny::em(" to continue."),
                                          shiny::br(),
                                          shiny::em("Saving is only possible once:"),
                                          shiny::br(),
                                          shiny::em(" - a valid dataset name is entered,"),
                                          shiny::br(),
                                          shiny::em(" - valid in- and output directories are selected."),
                                          shiny::br()
                                        )))
  server <- shiny::shinyServer(function(input, output, session) {
    outdirRoot <- outdir
    while ((!dir.exists(outdirRoot))&&(grepl("/", outdirRoot))) { outdirRoot <- gsub("/[^/]+$", "", outdirRoot) }
    #
    # Reactive values
    WHO <- shiny::reactiveVal(WhoAmI)
    USECUST <- shiny::reactiveVal(FALSE)
    INDIR <- shiny::reactiveVal(indir)
    INDIRtype <- shiny::reactiveVal("DiaNN")
    OUTDIRRoot <- shiny::reactiveVal(outdirRoot)
    OUTDIR <- shiny::reactiveVal(outdir)
    DATASETNAME <- shiny::reactiveVal(dtstNm)
    APPEND <- shiny::reactiveVal(append)
    # Non-reactive values
    wrkFlws <- setNames(WorkFlows, NULL)
    wrkFlw <- setNames(WorkFlow, NULL)
    #
    # Reactive functions to update UI
    prtInDr <- function(reactive = TRUE) {
      if (reactive) { msg = INDIR() } else { msg <- indir }
      shiny::renderUI({ shiny::em(msg) })
    }
    prtOutDr <- function(reactive = TRUE) {
      if (reactive) {
        if (APPEND()) {
          dtstNm2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", DATASETNAME())
          if (nchar(dtstNm2)) { OUTDIR(paste0(OUTDIRRoot(), "/", dtstNm2)) } else { OUTDIR(OUTDIRRoot()) }
        } else { OUTDIR(OUTDIRRoot()) }
        shiny::renderUI({ shiny::em(OUTDIR()) })
      } else {
        if (append) {
          if (nchar(dtstNm2)) { outdir <- (paste0(outdirRoot, "/", dtstNm2)) } else { outdir <- outdirRoot }
        } else { outdir <- outdirRoot }
        shiny::renderUI({ shiny::em(outdir) })
      }
    }
    updtOutDr <- function(reactive = TRUE) {
      dtstNm2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", dtstNm)
      if (reactive) {
        if (APPEND()) {
          dtstNm2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", DATASETNAME())
          if (nchar(dtstNm2)) { OUTDIR(paste0(OUTDIRRoot(), "/", dtstNm2)) } else { OUTDIR(OUTDIRRoot()) }
        } else { OUTDIR(OUTDIRRoot()) }
      } else {
        if (append) {
          if (nchar(dtstNm2)) { outdir <- (paste0(outdirRoot, "/", dtstNm2)) } else { outdir <- outdirRoot }
        } else { outdir <- outdirRoot }
      }
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
      lst[[1]] <- list(selectInput("SearchSoft",
                                   "Select search software used",
                                   names(SearchSoftware),
                                   dflt,
                                   width = "200px"))
      return(shiny::renderUI(lst))
    }
    #
    # Initial UI renders
    output$Indir <- prtInDr(reactive = FALSE)
    output$Outdir <- prtOutDr(reactive = FALSE)
    output$ScrptMsg <- shiny::renderUI(list(shiny::span(shiny::em("This script supports statistical tests on sample groups as long as each pair of groups includes 2 replicates or more."),
                                                        style = "color:red", .noWS = "outside")))
    output$Workflows <- shiny::renderUI({
      lst <- vector("list", 1)
      lst[[1]] <- list(selectInput("Workflow",
                                   "Select data analysis workflow",
                                   wrkFlws,
                                   wrkFlw,
                                   width = "100%"))
      return(lst)
    })
    output$inputType <- updtType(reactive = FALSE)
    #
    # Observers
    shiny::observeEvent(input$Who, { WHO(input$Who) })
    shiny::observeEvent(input$Workflow, { WorkFlow <<- input$Workflow })
    shiny::observe({
      shinyFiles::shinyDirChoose(input, "indir", roots = inRoot2)
      # Note: do not use default root and paths until they start behaving properly!!!
      {
        tmp <- input$indir
        if ("list" %in% class(tmp)) {
          INDIR(paste0(inRoot2[tmp$root], paste(tmp$path, collapse = "/")))
          output$Indir <- prtInDr()
          if ((dir.exists(OUTDIRRoot()))&&(nchar(DATASETNAME()))) { shinyjs::enable("saveBtn") }
          if ((dir.exists(OUTDIRRoot()))&&(nchar(dtstNm))) { shinyjs::enable("saveBtn") }
          output$inputType <- updtType()
        } else { shinyjs::disable("saveBtn") }
      }
    })
    shiny::observeEvent(input$SearchSoft, { INDIRtype(input$SearchSoft) })
    shiny::observe({
      shinyFiles::shinyDirChoose(input, "outdir", roots = outRoot)
      {
        tmp <- input$outdir
        if ("list" %in% class(tmp)) {
          OUTDIRRoot(paste0(outRoot[tmp$root], paste(tmp$path, collapse = "/")))
          OUTDIR(OUTDIRRoot())
          updtOutDr()
          output$Outdir <- prtOutDr()
          if ((dir.exists(INDIR()))&&(nchar(DATASETNAME()))) { shinyjs::enable("saveBtn") }
          if ((dir.exists(INDIR()))&&(nchar(dtstNm))) { shinyjs::enable("saveBtn") }
        } else { shinyjs::disable("saveBtn") }
      }
    })
    shiny::observeEvent(input$dtstNm, {
      if (nchar(input$dtstNm)) {
        DATASETNAME(input$dtstNm)
        updtOutDr()
        output$Outdir <- prtOutDr()
        if ((dir.exists(INDIR()))&&(dir.exists(OUTDIRRoot()))) { shinyjs::enable("saveBtn") }
      } else { shinyjs::disable("saveBtn") }
    })
    shiny::observeEvent(input$append, {
      APPEND(input$append)
      updtOutDr()
      output$Outdir <- prtOutDr()
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$vCPUs, {
      N.clust <<- input$vCPUs
    })
    shiny::observeEvent(input$writeRaws, {
      updtOutDr()
      writeRaws <<- input$writeRaws
    })
    shiny::observeEvent(input$Seed, {
      updtOutDr()
      mySeed <<- input$Seed
    })
    shiny::observeEvent(input$writeSearch, {
      writeSearch <<- input$writeSearch
    })
    shiny::observeEvent(input$ProcessedByUs, {
      ProcessedByUs <<- input$ProcessedByUs
      ObjNm <- "ProcessedByUs"
      if ((scrptType == "withReps")&&(ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) {
        AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
        tmp <- AllAnsw[1,]
        tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
        tmp$Value <- list(get(ObjNm))
        m <- match(ObjNm, AllAnsw$Parameter)
        if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
        AllAnsw <<- AllAnsw
      }
    })
    shiny::observeEvent(input$saveBtn, {
      WhoAmI <- WHO()
      if (WhoAmI == "Your name here") { WhoAmI <- NA }
      WhoAmI <<- WhoAmI
      indir <<- INDIR()
      outdir <<- OUTDIR()
      dtstNm <<- DATASETNAME()
      append <<- APPEND()
      SearchSoft <<- SearchSoftware[INDIRtype()]
      WorkFlow <<- input$Workflow
      shiny::stopApp()
    })
    shiny::observeEvent(input$cancel, { shiny::stopApp() })
    session$onSessionEnded(function() { shiny::stopApp() })
  })
  eval(parse(text = runApp), envir = .GlobalEnv)
  #
  SearchSoft %<o% SearchSoft
  #
  indir <- gsub("/+", "/", indir)
  outdir <- gsub("/+", "/", outdir)
  #
  tst <- gsub("[^A-Z,a-z,0-9]", "", unlist(sapply(c(dtstNm, indir, outdir), strsplit, "/")))
  tst <- sum(sapply(tst, TeachingDemos::char2seed, set = FALSE))
  tst2 <- lapply(unlist(strsplit(WhoAmI, " +")), function(x) {
    x <- unlist(strsplit(x, ""))
    return(x[1:(min(c(2, length(x))))])
  })
  tst2 <- c(sapply(tst2, function(x) { x[[1]] }),
            sapply(tst2, function(x) { if (length(x) >= 2) { x[[2]] } else { NA } }))
  tst2 <- tst2[which(!is.na(tst2))]
  tst2 <- paste(tst2[1:(min(c(2, length(tst2))))], collapse = "")
  tst3 <- gsub("[a-z, ]", "", dtstNm)
  projDir %<o% paste0(inRoot$Path[match("Temporary folder", inRoot$Folder)], "/", tst2, "_", tst3, "_", tst)
  wd %<o% paste0(projDir, "/Proc")
  if (!dir.exists(wd)) { dir.create(wd, recursive = TRUE) }
  setwd(projDir)
  #
  msg <- c(paste0("Dataset/project name:\n -> ", dtstNm), "",
           paste0("Script path:\n -> ", ScriptPath), "",
           paste0("Input directory:\n -> ", indir), "",
           paste0("Final output directory:\n -> ", outdir, "\n(folder will be created at the end of the workflow if it doesn't exist yet)"), "",
           paste0("Temporary work directory:\n -> ", wd), "",
           paste0("Seed:\n -> ", mySeed), "",
           "The data is processed within the temporary work directory, which is meant to have a short path to avoid issues with saving too long paths.",
           "The script then attempts to copy the files to the final output directory.",
           "If this fails (e.g. because of too long paths), then you can always manually copy the analysis from the temporary folder.")
  #cat(paste0(msg, "\n"))
  svDialogs::dlg_message(gsub("\n -> ", "\n>", msg), "ok", rstudio = FALSE)
  write(c(msg, ""), paste0(wd, "/Dataset details.txt"))
  #
  WorkFlow <- names(WorkFlows)[match(WorkFlow, WorkFlows)]
}
indir %<o% indir
outdir %<o% outdir
WhoAmI %<o% WhoAmI
WorkFlow %<o% WorkFlow
#create_project(path = wd, open = FALSE, rstudio = TRUE)
#usethis::proj_set(wd)
BckUpFl %<o% paste0(wd, "/Backup.RData")
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
#
# Below: this would apply if running a project with renv active, after reloading a lock
# if (renv) {
#   cran_req <- unique(c(cran_req, "renv"))
#   if (!is.na(RPath)) {
#     libPath <- paste0(RPath, "/proteoCraft")
#     Src <- paste0(libPath, "/extdata/R scripts/Sources/Save_Load_fun.R")
#     #system(paste0("open \"", Src, "\""))

#     source(Src, local = FALSE)
#   } else {
#     saveImgFun <- function(file) { # This one adapted from https://github.com/qsbase/qs2/issues/new?template=Blank+issue
#       obj <- base::ls(envir = .GlobalEnv)
#       if (exists(".obj", envir = .GlobalEnv)) {
#         obj <- unique(c(".obj", obj))
#         obj <- obj[which(sapply(obj, exists, envir = .GlobalEnv))]
#       }
#       obj <- grep("^[A-Za-z\\.][A-Za-z\\.0-9_]*$", obj, value = TRUE)
#       do.call(qs2::qs_savem,
#               c(lapply(obj, as.symbol),
#                 file = file,
#                 nthreads = parallel::detectCores()-1)
#       )
#     }
#   }
#   saveImgFun(BckUpFl)
#   inst <- as.data.frame(installed.packages())
#   if (!"renv" %in% inst$Package) { pak::pkg_install("renv", ask = FALSE, upgrade = TRUE, dependencies = TRUE) }
#   require("renv")
#   nuEnv <- (!file.exists("renv.lock"))
#   if (nuEnv) {
#     # This should never be the case: normally you would've run Reload_renv_from_lock_file.R already so there would be a project!
#     renv::init(force = TRUE)
#     # Session restarts here!!!
#   } else {
#     renv::load()
#     renv::restore()
#   }
# }
options(stringsAsFactors = FALSE)
options(install.packages.compile.from.source = "never")
#
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
# if (renv) {
#   if (!is.na(RPath)) {
#     libPath <- paste0(RPath, "/proteoCraft")
#     Src <- paste0(libPath, "/extdata/R scripts/Sources/Save_Load_fun.R")
#     #system(paste0("open \"", Src, "\""))

#     source(Src, local = FALSE)
#   } else {
#     loadFun %<o% function(file) {
#       tst <- try(qs2::qs_readm(file, env = globalenv(), nthreads = max(c(parallel::detectCores()-1, 1))), silent = TRUE)
#       if ("try-error" %in% class(tst)) { load(file, envir = globalenv()) }
#       require(proteoCraft) # This is because we have to remove the 
#     }
#   }
#   loadFun(paste0(getwd(), "/Proc/Backup.RData"))
# }
#loadFun(paste0(getwd(), "/Backup.RData"))
setwd(wd)
# if (renv) {
#   stopifnot(grepl("/Proc$", wd))
#   stopifnot("renv" %in% list.dirs(projDir, full.names = FALSE, recursive = FALSE))
#   #if (nuEnv) {
#   inst <- as.data.frame(installed.packages())
#   w <- which(!c(cran_req, bioc_req) %in% inst$Package)
#   while (length(w)) {
#     pack <- c(cran_req, bioc_req)[w][1]
#     if (pack %in% c("pak", "shiny", "uchardet", "openxlsx2", "taxize", "unimod")) {
#       if (pack == "pak") {
#         install.packages("pak", dependencies = TRUE)
#       }
#       if (pack == "shiny") {
#         # This is the recommended version!!! Others cause some issues which I have not managed to fix yet.
#         install.packages("https://cran.r-project.org/src/contrib/Archive/shiny/shiny_1.7.5.tar.gz", dependencies = TRUE)
#       }
#       if (pack == "uchardet") {
#         url <- "https://cran.r-project.org/src/contrib/Archive/uchardet/uchardet_1.1.1.tar.gz"
#         destfile <- "uchardet_1.1.1.tar.gz"
#         tst <- try(download.file(url, destfile, "curl"), silent = TRUE)
#         if ("try-error" %in% class(tst)) { try(download.file(url, destfile, "wget"), silent = TRUE) }
#         install.packages(destfile, dependencies = TRUE)
#         unlink(destfile)
#       }
#       if (pack == "openxlsx2") {
#         pak::pkg_install("JanMarvin/openxlsx2@v1.10", ask = FALSE, upgrade = TRUE, dependencies = TRUE) # ... until I can figure out what is happening...
#       }
#       # if (pack == "myTAI") {
#       #   pak::pkg_install("drostlab/myTAI@v0.9.3", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
#       # }
#       if (pack == "taxize") {
#         pak::pkg_install("ropensci/bold", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
#         pak::pkg_install("ropensci/taxize", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
#       }
#       if (pack == "unimod") {
#         pak::pkg_install("rformassspectrometry/unimod", ask = FALSE, upgrade = TRUE, dependencies = TRUE)
#       }
#     } else {
#       tst <- try(pak::pkg_install(pack, ask = FALSE, upgrade = TRUE, dependencies = TRUE), silent = TRUE)
#       if ("try-error" %in% class(tst)) {
#         tst <- try(install.packages(pack, dependencies = TRUE), silent = TRUE)
#       }
#       if ("try-error" %in% class(tst)) {
#         warning(paste0("Package ", pack, " wasn't installed properly, skipping..."))
#         cran_req <- cran_req[which(cran_req != pack)]
#         bioc_req <- bioc_req[which(bioc_req != pack)]
#       }
#     }
#     inst <- as.data.frame(installed.packages())
#     w <- which(!c(cran_req, bioc_req) %in% inst$Package)
#   }
#   #}
# }
for (pack in c(cran_req, bioc_req, "proteoCraft")) {
  try(library(pack, character.only = TRUE), silent = TRUE)
  if ("try-error" %in% class(tst)) {
    warning(paste0("Package ", pack, " could not be loaded!"))
  }
  # 
  # add something here to catch issues with packages which cannot be unloaded...
}
tst <- try(normalizePath(rawrr:::.rawrrAssembly(), winslash = "/"), silent = TRUE)
if (("try-error" %in% class(tst))||(!file.exists(tst))) {
  rawrr::installRawrrExe()
}
data.table::setDTthreads(threads = detectCores()-1)
#
inst <- as.data.frame(installed.packages())
if ((!"proteoCraft" %in% inst$Package)||((exists("updt_proteoCraft"))&&(updt_proteoCraft))) {
  locFl <- paste0(homePath, "/Default_locations.xlsx")
  locs <- openxlsx2::read_xlsx(locFl)
  pckgloc <- paste0(locs$Path[match("Server share", locs$Folder)], "/proteoCraft_package")
  goOn <- dir.exists(pckgloc)
  if (goOn) {
    pckgs <- list.files(pckgloc, "\\.tar\\.gz$", full.names = TRUE)
    goOn <- length(pckgs) > 0
    if (goOn) {
      pckgs <- data.frame(Name = pckgs, check.names = FALSE)
      temp <- strsplit(gsub(".*/proteoCraft_|\\.tar\\.gz$", "", pckgs$Name), "\\.")
      pckgs[, paste0("V", 1:4)] <- as.data.frame(t(sapply(temp, function(x) {
        x <- as.numeric(unlist(x))
        if (length(x) < 4) { x <- c(x, rep(0, 4-length(x))) }
        return(x)
      })))
      k <- 1
      pckgs <- pckgs[which(pckgs[[paste0("V", k)]] == max(pckgs[[paste0("V", k)]])), , drop = FALSE]
      while (nrow(pckgs) > 1) {
        k <- k+1  
        pckgs <- pckgs[which(pckgs[[paste0("V", k)]] == max(pckgs[[paste0("V", k)]])), , drop = FALSE]
      }
      nuPack <- pckgs$Name
    } else {
      nuPack <- rstudioapi::selectFile("Select tarball (\\.tar.gz) for updating the proteoCraft package", filter = "tarball (*.tar.gz)")
    }
  }
  if ((!is.na(nuPack))&&(file.exists(nuPack))) {
    cat("Updating proteoCraft package...\n")
    unloadNamespace("proteoCraft")
    remove.packages("proteoCraft")
    install.packages(nuPack, dependencies = TRUE)
  }
}
require(proteoCraft)
set.seed(mySeed)
if (renv) {
  if (nuEnv) {
    renv::snapshot(force = TRUE, prompt = FALSE, exclude = "fastSave")
  }
}
LocAnalysis %<o% (WorkFlow %in% c("LOCALISATION", "LOCALIZATION"))
IsBioID %<o% (gsub(" |_|-|\\.", "", toupper(WorkFlow)) == "BIOID")
#
# Update values
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
parSrc %<o% paste0(libPath, "/extdata/R scripts/Sources/make_check_Cluster.R")
#
# Detect backups and decide whether to keep them
#labelMode <- match(LabelType, c("LFQ", "Isobaric"))
FracMapNm %<o% "Fractions map"
FracMapPath %<o% paste0(wd, "/", FracMapNm, ".csv")
intPrtFst %<o% paste0(wd, "/Proteins of interest.fasta")
allBckps %<o% data.frame(Full = FracMapPath,
                         File = basename(FracMapPath),
                         Role = "Map of MS files to biological samples")
allBckps$ObjNm <- list("FracMap_reloaded")
if (scrptType == "withReps") {
  ExpMapNm %<o% "Experiment map"
  ExpMapPath %<o% paste0(wd, "/", ExpMapNm, ".csv")
  tmpDF <- data.frame(File = c(basename(ExpMapPath),
                               "Factors.RData",
                               "Parameters.csv"),
                      Role = c("Experimental structure map",
                               "Experimental factors",
                               "Analysis parameters"))
  tmpDF$ObjNm <- list("ExpMap",
                      c("Factors", "FactorsLevels"),
                      "Param")
  tmpDF$Full <- paste0(wd, "/", tmpDF$File)
  allBckps <- rbind(allBckps, tmpDF)
}
if (scrptType == "noReps") {
  SamplesMapNm %<o% "Samples map"
  SamplesMapPath %<o% paste0(wd, "/", SamplesMapNm, ".csv")
  allBckps <- rbind(allBckps,
                    data.frame(Full = SamplesMapPath,
                               File = basename(SamplesMapPath),
                               Role = "Experimental structure map",
                               ObjNm = "SamplesMap"))
}
tmpDF <- data.frame(File = c(basename(intPrtFst),
                             "Parsed_annotations.RData",
                             "evmatch.RData"),
                    Role = c("FASTA of proteins of special interest",
                             "Parsed functional annotations",
                             "Matches of peptide sequences to parent proteins"),
                    ObjNm = c("prot.list",
                              "Parsed_annotations",
                              "evmatch"))
tmpDF$Full <- paste0(wd, "/", tmpDF$File)
allBckps <- rbind(allBckps, tmpDF)
if (SearchSoft %in% c("DIANN", "FRAGPIPE")) {
  m <- match(SearchSoft, c("DIANN", "FRAGPIPE"))
  PSMsBckp %<o% paste0(c("diaNN", "FragPipe")[m], " PSMs converted to MQ-like format.RData")
  tmpDF <- data.frame(File = PSMsBckp,
                      Role = "Processed PSMs",
                      ObjNm = paste0("ev_", c("DIANN", "FP")[m], "2MQ"))
  tmpDF$Full <- paste0(wd, "/", tmpDF$File)
  allBckps <- rbind(allBckps, tmpDF)
}
#View(allBckps[which(!file.exists(allBckps$Full)),])
allBckps <- allBckps[which(file.exists(allBckps$Full)),]
reloadedBckps %<o% allBckps[NULL,]
if (nrow(allBckps)) {
  allBckps$Dir <- dirname(allBckps$Full)
  allBckps$Value <- allBckps$File
  w <- which(allBckps$Dir != wd)
  if (length(w)) { allBckps$Value[w] <- allBckps$Full[w] }
  allBckps$Value <- paste0(do.call(paste, c(allBckps[, c("Role", "Value")], sep = " (file = ")), ")")
  bckps2Reload %<o% dlg_list(allBckps$Value, allBckps$Value, TRUE, title = "Backups detected: which should we reload?")$res
  if (length(bckps2Reload)) {
    reloadedBckps <- allBckps[match(bckps2Reload, allBckps$Value),]
    for (i in 1:nrow(reloadedBckps)) {
      ext <- tolower(gsub(".*\\.", "", reloadedBckps$File[i]))
      if (ext == "fasta") { fastas_reloaded <- reloadedBckps$Full[i] }
      if (ext == "rdata") { loadFun(reloadedBckps$Full[i]) }
      if (ext == "csv") {
        tmp <- read.csv(reloadedBckps$Full[i], check.names = FALSE)
        areUok <- TRUE
        if (reloadedBckps$Role[i] == "Map of MS files to biological samples") {
          colnames(tmp)[which(colnames(tmp) == "Raw.file")] <- "Raw file" # Backwards compatibility
          if (sum(!c("Parent sample", "MQ.Exp")#[labelMode]
                  %in% colnames(tmp)) == 2) {
            warning("Invalid Fractions map reloaded, ignoring...")
            areUok <- FALSE
          }
        }
        if ((reloadedBckps$Role[i] == "Experimental structure map")&&(scrptType == "withReps")) {
          # Backwards compatibility
          colnames(tmp)[which(colnames(tmp) == "Sample.name")] <- "Sample name" 
          #colnames(tmp)[which(colnames(tmp) == "Isobaric.set")] <- "Isobaric set"
          colnames(tmp)[which(colnames(tmp) == "Isobaric.label")] <- "Isobaric label"
          colnames(tmp)[which(colnames(tmp) == "Isobaric.label.details")] <- "Isobaric label details"
          if (exists("FracMap")) {
            expKl <- c("MQ.Exp", "Parent sample")
            expKl <- expKl[which(expKl %in% colnames(FracMap))[1]]
            # tst <- unique(FracMap[[expKl]][which(FracMap$Use)])
            # if (sum(!tst %in% tmp[[c("Sample name", "MQ.Exp")[labelMode]]])) {
            #   warning("Invalid Experiment map reloaded, ignoring...")
            #   areUok <- FALSE
            # }
          }
        }
        rm(list = c(reloadedBckps$ObjNm[[i]])) # We want to make sure that no invalid version of the object lingers
        if (areUok) { assign(reloadedBckps$ObjNm[[i]], tmp) }
      }
    }
    .obj <- unique(c(.obj, unlist(reloadedBckps$ObjNm)))
  }
}
if ((!nrow(reloadedBckps))||(!"FASTA of proteins of special interest" %in% reloadedBckps$Role)) {
  loadInt %<o% c(TRUE, FALSE)[match(dlg_message("Load a fasta of proteins of interest?", "yesno")$res, c("yes", "no"))]
  if (loadInt) {
    intFast <- selectFile(paste0("Select proteins of interest fasta", intPrtFst), path = wd)
    if (!is.null(intFast)) {
      intFast <- gsub("^~", normalizePath(Sys.getenv("HOME"), winslash = "/"), intFast)
      if (intFast != intPrtFst) {
        file.copy(intFast, intPrtFst, TRUE)
        cat(paste0("   FYI: a copy of your input fasta has been saved at \"", intPrtFst, "\"..."))
      }
      fastas <- unique(c(fastas, intPrtFst))
    }
  }
}
Reuse_Prot_matches %<o% ("Matches of peptide sequences to parent proteins" %in% reloadedBckps$Role)
ReLoadPSMsBckp %<o% ("Processed PSMs" %in% reloadedBckps$Role)
#
saveImgFun(BckUpFl)
