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
  if (exists("inDirs")) {
    tst <- vapply(inRoot2, function(rt) {
      sum(vapply(inDirs, function(indir) {
        grepl(proteoCraft::topattern(rt), indir)
      }, TRUE)) > 0
    }, TRUE)
    inRoot2 <- c(inRoot2[which(tst)][order(nchar(inRoot2[which(tst)]), decreasing = TRUE)],
                 inRoot2[which(!tst)])
  }
  if (!exists("WhoAmI")) { WhoAmI <- Sys.getenv("USERNAME") }
  #
  # outRoot <- inRoot2[c(grep("^Results delivery folder ", names(inRoot2)),
  #                      grep("^Archive folder ", names(inRoot2)),
  #                      grep("^Results delivery folder |^Archive folder ", names(inRoot2), invert = TRUE))]
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
  if (exists("inDirs")) {
    inDirs <- inDirs[which(dir.exists(inDirs))]
    if (!length(inDirs)) { rm(inDirs) }
  }
  if (!exists("inDirs")) { inDirs <- inRoot2[grep("^Search folder ", names(inRoot2))[1]] }
  defRt <- inRoot2[match(gsub("/.*", "/", inRoot2[grep("^Search folder ", names(inRoot2))]), inRoot2)]
  tmp <- c("MaxQuant", "DiaNN", "FragPipe", "Proteome Discoverer")
  SearchSoftware %<o% setNames(gsub(" ", "", toupper(tmp)), tmp)
  srchSoftOpt <- names(SearchSoftware)
  #
  scriptPar <- readLines(ScriptPath)
  scriptPar <- gsub("^###-\\|-### *", "", grep("^###-\\|-### *", scriptPar, value = TRUE))
  appNm <- "Start analysis"
  dtstNm2 <- gsub(":|\\*|\\?|<|>|\\||/", "-", dtstNm)
  updt_Type0 <- function(searchDir, default = NULL) {
    if ((is.null(default))||(!default %in% srchSoftOpt)) {
      fls <- list.files(searchDir)
      tsv2txt <- length(grep("\\.tsv$", fls))/length(grep("\\.txt$", fls))
      default <- c(FALSE,
                   length(grep("\\.log\\.txt$", fls)) > 0,
                   length(grep("\\.fp-manifest$", fls)) > 0,
                   length(grep("\\.xlsx$", fls)) > 0)
      if (!sum(default)) { default[1] <- TRUE }
      if (sum(default) > 1) {
        w <- which(default)
        default[w[2:length(w)]] <- FALSE
      }
      default <- srchSoftOpt[which(default)]
    }
    return(default)
  }
  if ((!exists("SearchSoft"))||(sum(!SearchSoft %in% names(SearchSoftware)))||(length(SearchSoft) != length(inDirs))) {
    SearchSoft <- vapply(inDirs, updt_Type0, "")
  }
  nr0 <- length(inDirs)
  inputTblDflt <- data.frame("Input search folder" = "",
                             "Value" = gsub("/[^/]+$", "", inDirs[1]),
                             "Search engine" = SearchSoftware[1],
                             "Remove" = "",
                             check.names = FALSE)
  rownames(inputTblDflt) <- NULL
  inputTbl <- data.frame("Input search folder" = rep("", nr0),
                         "Value" = inDirs,
                         "Search engine" = SearchSoft,
                         "Remove" = "",
                         check.names = FALSE)
  rownames(inputTbl) <- NULL
  nr0 <- nrow(inputTbl)
  inputTbl2 <- inputTbl
  inputTbl2$`Input search folder` <- vapply(paste0("selectDir___", 1:nr0), function(id) {
    as.character(shiny::actionButton(id, "Select input"))
  }, "")
  fSlct0 <- function(i, data = inputTbl, opt = srchSoftOpt) {
    srchSoftOpt2 <- paste0("<option value=\"", opt, "\"",
                           c("", " selected")[(opt == data$`Search engine`[i])+1],
                           ">", opt, "</option>", collapse = "")
    return(paste0("<select id=\"searchSoft___", i, "\" style=\"width:200px;\">", srchSoftOpt2, "</select>"))
  }
  inputTbl2$`Search engine` <- vapply(1:nr0, fSlct0, "")
  inputTbl2$Remove <- vapply(paste0("removeMe___", 1:nr0), function(id) {
    as.character(shiny::actionButton(id, "Remove dataset"))
  }, "")
  colnames(inputTbl2)[4] <- ""
  #
  ui <- shiny::shinyUI(
    shiny::fluidPage(
      shinyjs::useShinyjs(),
      shinyWidgets::setBackgroundColor( # Doesn't work
        color = c(#"#F8F8FF",
          "#EBEFF7"),
        gradient = "linear",
        direction = "bottom"
      ),
      shinyjs::extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
      shiny::tags$head(shiny::tags$style(shiny::HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
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
            shiny::h4(shiny::strong("INPUTS")),
            shiny::br(),
            shiny::h5("Select output folder from a supported search engine:"),
            DT::DTOutput("inDirs"),
            shiny::actionButton("addBtn", "+ add input folder"),
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
  #eval(parse(text = runApp), envir = .GlobalEnv)
  slctDirXprs <- expression({
    dat <- INPUTTBL()
    nr <- nrow(dat)
    drs <- dat$Value
    dflt <- dat$Value[i]
    if (!dir.exists(dflt)) { dflt <- rev(drs[which(dir.exists(drs))]) }
    if (!length(dflt)) { dflt <- "C:/" }
    if (length(dflt) > 1) { dflt <- dflt[1] }
    dr <- rstudioapi::selectDirectory(paste0("Select ", c("", "additional ")[(i>1)+1], "input folder"),
                                      path = dflt)
    if ((length(dr) == 1)&&(!is.na(dr))&&(dir.exists(dr))) {
      if (dr %in% drs) {
        warning("You already selected this directory! Ignoring...")
      } else {
        dat2 <- INPUTTBL2()
        dat$Value[i] <- dat2$Value[i] <- dr
        tp <- updt_Type0(dr)
        dat$`Search engine`[i] <- tp
        INPUTTBL(dat)
        inputTbl <<- dat
        inputTbl2$`Search engine`[i] <- fSlct0(i, dat)
        INPUTTBL2(dat2)
        inputTbl2 <<- dat2
        output$inDirs <- updt_inDirs()
      }
    }
  })
  #eval(parse(text = runApp), envir = .GlobalEnv)
  srchSoftXprs <- expression({
    dat <- INPUTTBL()
    dat$`Search engine`[i] <- ev$value
    INPUTTBL(dat)
    if ("" %in% colnames(dat)) { stop() }
    inputTbl <<- dat
    dat2 <- INPUTTBL2()
    tmp <<- fSlct0(i, dat)
    dat2$`Search engine`[i] <- tmp
    INPUTTBL2(dat2)
    inputTbl2 <<- dat2
    output$inDirs <- updt_inDirs()
  })
  #eval(parse(text = runApp), envir = .GlobalEnv)
  rmvDirXprs <- expression({
    dat <- INPUTTBL()
    nr <- nrow(dat)
    w <- which(1:nr != i)
    if (length(w)) { # Can't remove all!!!
      dat <- dat[w,]
      INPUTTBL(dat)
      if ("" %in% colnames(dat)) { stop() }
      inputTbl <<- dat
      dat2 <- INPUTTBL2()[w,]
      nr2 <- nrow(dat)
      # Re-generate table with IDs from updated row position
      dat2$`Input search folder` <- vapply(paste0("selectDir___", 1:nr2), function(id) {
        as.character(shiny::actionButton(id, "Select input"))
      }, "")
      dat2$`Search engine` <- vapply(1:nr2, fSlct0, "", dat)
      dat2[[4]] <- vapply(paste0("removeMe___", 1:nr2), function(id) {
        as.character(shiny::actionButton(id, "Remove dataset"))
      }, "")
      INPUTTBL2(dat2)
      inputTbl2 <<- dat2
      output$inDirs <- updt_inDirs()
    }
  })
  #eval(parse(text = runApp), envir = .GlobalEnv)
  server <- shiny::shinyServer(function(input, output, session) {
    WHO <- shiny::reactiveVal(WhoAmI)
    USECUST <- shiny::reactiveVal(FALSE)
    DATASETNAME <- shiny::reactiveVal(dtstNm)
    INPUTTBL <- shiny::reactiveVal(inputTbl)
    INPUTTBL2 <- shiny::reactiveVal(inputTbl2)
    # Non-reactive values
    wrkFlws <- setNames(WorkFlows, NULL)
    wrkFlw <- setNames(WorkFlow, NULL)
    #
    # Reactive functions to update UI
    updt_inDirs <- function(reactive = TRUE) {
      if (reactive) { dat2 <- INPUTTBL2() } else { dat2 <- inputTbl2 }
      return(DT::renderDT({ dat2 },
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
                                         columnDefs = list(list(width = "100px", targets = 0),
                                                           list(width = "400px", targets = 1),
                                                           list(width = "200px", targets = 2),
                                                           list(width = "100px", targets = 3)),
                                         scrollX = FALSE),
                          # the callback is essential to capture the inputs in each row
                          callback = DT::JS("
// Buttons
table.on('click', 'button', function() {
  var id = this.id;
  Shiny.setInputValue('dt_event', {type: 'button', id: id}, {priority: 'event'});
});
// Dropdowns
table.on('change', 'select', function() {
  var id = this.id;
  var val = this.value;
  Shiny.setInputValue('dt_event', {type: 'select', id: id, value: val}, {priority: 'event'});
});
")))
    }
    observeEvent(input$dt_event, {
      ev <- input$dt_event
      id <- ev$id
      i <- as.integer(gsub(".*___", "", id))
      root <- gsub("___.*", "", id)
      if (ev$type == "button") {
        # Directory selection
        if(root == "selectDir") { eval(slctDirXprs) }
        # Row removal?
        if(root == "removeMe") { eval(rmvDirXprs) }
      }
      if (ev$type == "select") {
        # Software selection
        if(root == "searchSoft") { eval(srchSoftXprs) }
      }
    })
    observeEvent(input$addBtn, {
      dat <- INPUTTBL()
      dat2 <- INPUTTBL2()
      drs <- dat$Value
      nr <- length(drs)
      i <- nr+1
      datDflt <- data.frame("Input search folder" = "",
                            "Value" = gsub("/[^/]+$", "", dat$Value[nr]),
                            "Search engine" = dat$`Search engine`[nr],
                            "Remove" = "",
                            check.names = FALSE)
      rownames(datDflt) <- NULL
      dat <- rbind(dat, datDflt)
      INPUTTBL(dat)
      if ("" %in% colnames(dat)) { stop() }
      inputTbl <<- dat
      datDflt2 <- datDflt
      datDflt2$"Input search folder" <- as.character(shiny::actionButton(paste0("selectDir___", i), "Select input"))
      datDflt2$`Search engine` <- fSlct0(i, datDflt)
      datDflt2$Remove <- as.character(shiny::actionButton(paste0("removeMe___", i), "Remove dataset"))
      colnames(datDflt2)[4] <- ""
      dat2 <- rbind(dat2, datDflt2)
      INPUTTBL2(dat2)
      inputTbl2 <<- dat2
      output$inDirs <- updt_inDirs()
    })
    #
    # Initial UI renders
    output$inDirs <- updt_inDirs(reactive = FALSE)
    output$Workflows <- shiny::renderUI({
      lst <- vector("list", 1)
      lst[[1]] <- list(selectInput("Workflow",
                                   "Select data analysis workflow",
                                   wrkFlws,
                                   wrkFlw,
                                   width = "100%"))
      return(lst)
    })
    #
    # Observers
    shiny::observeEvent(input$Who, { WHO(input$Who) })
    shiny::observeEvent(input$Workflow, { WorkFlow <<- input$Workflow })
    #
    shiny::observeEvent(input$dtstNm, {
      if (nchar(input$dtstNm)) {
        DATASETNAME(input$dtstNm)
        dat <- INPUTTBL()
        if (length(sum(dir.exists(dat$Value)))) { shinyjs::enable("saveBtn") } else { shinyjs::disable("saveBtn") }
      } else { shinyjs::disable("saveBtn") }
    })
    #
    shiny::observeEvent(input$vCPUs, {
      N.clust <<- input$vCPUs
    })
    shiny::observeEvent(input$writeRaws, {
      #updt_OutDr()
      writeRaws <<- input$writeRaws
    })
    shiny::observeEvent(input$Seed, {
      #updt_OutDr()
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
      dtstNm <<- DATASETNAME()
      WorkFlow <<- input$Workflow
      inputTbl <<- INPUTTBL()
      inDirs <<- inputTbl$Value
      SearchSoft <<- inputTbl$`Search engine`
      shiny::stopApp()
    })
    shiny::observeEvent(input$cancel, { shiny::stopApp() })
    session$onSessionEnded(function() { shiny::stopApp() })
  })
  eval(parse(text = runApp), envir = .GlobalEnv)
  #
  #
  inDirs <- gsub("/+", "/", inDirs)
  #
  tst <- gsub("[^A-Z,a-z,0-9]", "", unlist(sapply(c(dtstNm, inDirs), strsplit, "/")))
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
           paste0("Input director", c("y", "ies")[(length(inDirs) > 1) + 1], ":\n -> ", paste(inDirs, collapse = "\n -> ")), "",
           paste0("Temporary work directory:\n -> ", wd), "",
           paste0("Seed:\n -> ", mySeed), "",
           "The data is processed within the temporary work directory, which is meant to have a short path to avoid issues with saving too long paths.",
           "The script then attempts to copy the files to the final output directory.",
           "If this fails (e.g. because of too long paths), then you can always manually copy the analysis from the temporary folder.")
  #cat(paste0(msg, "\n"))
  svDialogs::dlg_message(msg, "ok", rstudio = TRUE)
  write(c(msg, ""), paste0(wd, "/Dataset details.txt"))
  #
  WorkFlow <- names(WorkFlows)[match(WorkFlow, WorkFlows)]
}
SearchSoft %<o% SearchSoftware[SearchSoft]
if (sum(SearchSoft %in% names(SearchSoftware))) { SearchSoft %<o% SearchSoftware[SearchSoft] } # Because I am stupid!
inDirs %<o% inDirs
WhoAmI %<o% WhoAmI
WorkFlow %<o% WorkFlow
BckUpFl %<o% paste0(wd, "/Backup.RData")
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
#
options(stringsAsFactors = FALSE)
options(install.packages.compile.from.source = "never")
#
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")

setwd(wd)


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
    renv::snapshot(force = TRUE, prompt = FALSE#, exclude = "fastSave" # We are not including fastSave anymore
    )
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
tmp <- lapply(1:length(inDirs), function(dir_i) {
  # No need for MQ, it's faster and easier to just reload and do the minimal processing we do
  if (SearchSoft[dir_i] %in% c("DIANN", "FRAGPIPE")) {
    m <- match(SearchSoft[dir_i], c("DIANN", "FRAGPIPE"))
    psmsBckpFl_i <- paste0(c("diaNN", "FragPipe")[m], " PSMs converted to MQ-like format_", dir_i, ".RData")
    tmpDF <- data.frame(File = psmsBckpFl_i,
                        Role = "Processed PSMs",
                        ObjNm = paste0("ev_", c("DIANN", "FP")[m], "2MQ_", dir_i))
    tmpDF$Full <- paste0(wd, "/", tmpDF$File)
    return(tmpDF)
  }
})
tmp <- plyr::rbind.fill(tmp)
allBckps <- rbind(allBckps, tmp)
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
  #
  unusedBckps <- allBckps[which(!allBckps$Value %in% bckps2Reload),]
  if (nrow(unusedBckps)) { # Remove the object if it already exists and we do not want to reload it!
    suppressWarnings(rm(list = unlist(unusedBckps$ObjNm)))
  }
  if (length(bckps2Reload)) {
    reloadedBckps <- allBckps[match(bckps2Reload, allBckps$Value),]
    for (i in 1:nrow(reloadedBckps)) { #i <- 1
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
          colnames(tmp)[which(colnames(tmp) == "Isobaric.label")] <- "Isobaric label"
          colnames(tmp)[which(colnames(tmp) == "Isobaric.label.details")] <- "Isobaric label details"
          if (exists("FracMap")) {
            expKl <- c("MQ.Exp", "Parent sample")
            expKl <- expKl[which(expKl %in% colnames(FracMap))[1]]
          }
        }
        if (exists(reloadedBckps$ObjNm[[i]])) {
          # We want to make sure that no invalid version of the object lingers
          rm(list = c(reloadedBckps$ObjNm[[i]]))
        }
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
#
saveImgFun(BckUpFl)
