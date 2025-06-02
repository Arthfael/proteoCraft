#
options(stringsAsFactors = FALSE)
options(buildtools.check = function(action) TRUE )
rm(list = ls()[!ls() %in% c("fls")])

if (!require(pak)) { install.packages("pak") }

getInt <- FALSE
cran_req <- c("pak",
              "officer",
              "devtools",
              "data.table",
              "XML",
              "plyr",
              "svDialogs",
              "parallel",
              "snow",
              "ggplot2",
              "rpanel",
              "shiny",
              "shinyjs",
              "DT",
              "shinycssloaders",
              "viridis",
              "ggnewscale",
              #"pracma",
              "ggformula",
              "splines",
              "curl",
              "gtools")
tst <- sapply(cran_req, function(pck) { require(pck, character.only = TRUE, quietly = TRUE) })
w <- which(!tst)
if (length(w)) { pak::pkg_install(cran_req[w]) }
for (pck in cran_req) { library(pck, character.only = TRUE) }
if (!require(proteoCraft)) {
  tst <- try(pak::pkg_install("Arthfael/proteoCraft"))
  if ("try-error" %in% clas(tst)) {
    tst <- try(devtools::install_github("Arthfael/proteoCraft"))
  }
  library(proteoCraft)
}

# Work-in-Progress only behaviour
WIP <- FALSE
#WIP <- TRUE
# From https://stackoverflow.com/questions/2643719/find-contiguous-stretches-of-equal-data-in-a-vector
contigEq <- function(data, value, position) {
  id <- cumsum(c(1, as.numeric(diff(data) != 0)))
  id == id[position]
}

RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/aRmel")
# Set Shiny options, load functions for creating a Word report, create Excel styles
Src <- paste0(libPath, "/extdata/R scripts/Sources/ShinyOpt_Styles_and_Report.R")
#rstudioapi::documentOpen(Src)
#rstudioapi::documentOpen(paste0("H:/aRmel_package/aRmel/inst/extdata/R scripts/Sources/", basename(Src)))
source(Src)

# Create cluster
parSrc <- paste0(libPath, "/extdata/R scripts/Sources/make_check_Cluster.R")
#rstudioapi::documentOpen(parSrc)
source(parSrc)

# Some useful functions
cleanRawNm <- function(rawFileNm) { gsub("\\.raw$", "", rawFileNm, ignore.case = TRUE) }
environment(cleanRawNm) <- .GlobalEnv
clusterExport(parClust, "cleanRawNm", envir = environment())

# Install/load packages
rawrrVers <- c("github",
               "bioc",
               "bioc_1.11.14")
tst <- try({
  if (!require(rawrr, quietly = TRUE)) {
    rawrrVers <- dlg_list(rawrrVers, rawrrVers[3], title = "Select which rawrr version should be installed")$res
    if (rawrrVers == "github") {
      if(!require(devtools)) { install.packages("devtools") }
      devtools::install_github("cpanse/rawrr")
      #install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawrr_0.2.1.tar.gz', repo = NULL) # Old address, now on github
    }
    if (rawrrVers == "bioc") {
      if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
      BiocManager::install("rawrr", version = "1.11")
      BiocManager::install("rawrr")
    }
    if (rawrrVers == "bioc_1.11.14") {
      url <- "https://bioconductor.org/packages/3.19/bioc/src/contrib/Archive/rawrr/rawrr_1.11.14.tar.gz"
      require(curl)
      dstfl <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/rawrr_1.11.14.tar.gz")
      curl_download(url, dstfl)
      install.packages(dstfl)
    }
    # 
    require(rawrr)
    yesRawFileReaderLicenseIsAccepted <- function () {
      licenseFile <- file.path(system.file(package = "rawrr"), 
                               "rawrrassembly", "RawFileReaderLicense.txt")
      stopifnot(file.exists(licenseFile))
      eulaFile <- file.path(rawrrAssemblyPath(), "eula.txt")
      msg <- c("# By changing the setting below to TRUE you are accepting ", 
               "the Thermo License agreement.")
      if (!file.exists(eulaFile)) {
        file.show(licenseFile)
        response <- "y"
        if (tolower(response) == "y") {
          if (isFALSE(dir.exists(dirname(eulaFile)))) {
            dir.create(dirname(eulaFile), recursive = TRUE)
          }
          fileConn <- file(eulaFile)
          writeLines(paste(msg, paste0("# ", date()), "eula=true", 
                           sep = "\n"), fileConn)
          close(fileConn)
          return(TRUE %in% grepl("eula=true", tolower(readLines(eulaFile))))
        }
      } else { return(TRUE %in% grepl("eula=true", tolower(readLines(eulaFile)))) }
      msg <- "Yes, we accept Thermo's License agreement, get on with it!"
      cat(msg, "\n")
    }
    installRawFileReaderDLLsNoAcpt <- function (sourceUrl = .thermofisherlsmsUrl(), ...) {
      rawfileReaderDLLsPath <- rawrrAssemblyPath()
      if (isTRUE(dir.exists(rawfileReaderDLLsPath))) {
        msg <- sprintf("removing files in directory '%s'", rawfileReaderDLLsPath)
        message(msg)
        file.remove(file.path(rawrrAssemblyPath(), list.files(rawrrAssemblyPath())))
      }
      if (isFALSE(dir.exists(rawfileReaderDLLsPath))) {
        dir.create(rawfileReaderDLLsPath, recursive = TRUE)
      }
      #
      stopifnot(yesRawFileReaderLicenseIsAccepted())
      .rawfileReaderDLLs <- getAnywhere(.rawfileReaderDLLs)
      .rawfileReaderDLLs <- .rawfileReaderDLLs$objs[[1]]
      rv <- vapply(.rawfileReaderDLLs(), function(dll) {
        destfile <- file.path(rawfileReaderDLLsPath, dll)
        download.file(file.path(sourceUrl, dll), destfile = destfile, 
                      mode = "wb", ...)
      }, 0)
      rv
    }
    installRawFileReaderDLLsNoAcpt(sourceUrl = .thermofisherlsmsUrl())
    #rawrr::installRawFileReaderDLLs(sourceUrl = .thermofisherlsmsUrl())
    #require(rawrr); getAnywhere(".isRawFileReaderLicenseAccepted")
    rawrr::installRawrrExe()
  }
  require(rawrr)
}, silent = TRUE)
getInt <- !("try-error" %in% class(tst))
#
if (getInt) {
  ticFun <- function(fl) { #fl <- fls[1]
    x <- try(readChromatogram(fl, type = "tic"), silent = TRUE)
    if (!"try-error" %in% class(x)) {
      res <- list(Outcome = TRUE,
                  Output = data.frame("Raw file" = fl,
                                      "Raw file name" = cleanRawNm(basename(fl)),
                                      "Retention time" = as.numeric(x$times),
                                      "Intensity" = as.numeric(x$intensities),
                                      check.names = FALSE))
    } else { res <- list(Outcome = FALSE,
                         Error = x) }
    return(res)
  }
  bpcFun <- function(fl) { #fl <- fls[1]
    x <- try(readChromatogram(fl, type = "bpc"), silent = TRUE)
    if (!"try-error" %in% class(x)) {
      res <- list(Outcome = TRUE,
                  Output = data.frame("Raw file" = fl,
                                      "Raw file name" = cleanRawNm(basename(fl)),
                                      "Retention time" = as.numeric(x$times),
                                      "Intensity" = as.numeric(x$intensities),
                                      check.names = FALSE))
    } else { res <- list(Outcome = FALSE,
                         Error = x) }
    return(res)
  }
  xicFun <- function(fl, mass, tol, filter) { #fl <- fls[1]
    x <- try(readChromatogram(fl, type = "xic", mass, tol, filter), silent = TRUE)
    if (!"try-error" %in% class(x)) {
      res <- list(Outcome = TRUE,
                  Output = data.frame("Raw file" = fl,
                                      "Raw file name" = cleanRawNm(basename(fl)),
                                      "M/Z" = mass,
                                      "Tolerance" = tol,
                                      "Retention time" = as.numeric(x[[1]]$times),
                                      "Intensity" = as.numeric(x[[1]]$intensities),
                                      check.names = FALSE))
    } else { res <- list(Outcome = FALSE,
                         Error = x) }
    return(res)
  }
  environment(ticFun) <- .GlobalEnv
  environment(bpcFun) <- .GlobalEnv
  environment(xicFun) <- .GlobalEnv
}
#require(RaMS)
dflt <- "B:/group/lsfgrp/Mass_Spec/Acquired_data_v2"
if (!dir.exists(dflt)) { dflt <- "B:/archive/lsfgrp/MS/Acquired_data" }
if (!dir.exists(dflt)) { dflt <- "D:/Data" }
if (!dir.exists(dflt)) { dflt <- "D:/groups_temp" }

# Raw files
if (exists("fls")) {
  fls <- fls[which(file.exists(fls))]
} else { fls <- c() }
l <- length(fls)
if (l) {
  msg <- paste0("There are ", l, " input Raw files already present in environment from a previous run:")
  opt <- c("Remove them                                                                                                             ",
           "Reprocess them                                                                                                          ")
  startFresh <- c(TRUE, FALSE)[match(dlg_list(opt, opt[1], title = msg)$res, opt)]
  if (startFresh) { fls <- c() }
}
l <- length(fls)
if (!l) {
  filt <- matrix(data = c("Thermo raw file", "*.raw;*.RAW"), ncol = 2,
                 dimnames = list("raw file"))
  #fls <- "B:/archive/lsfgrp/MS/Acquired data/frimlgrp/LCMS_JiFrATeplova1_m"
  fls <- normalizePath(choose.files(paste0(dflt, "/*.raw"), filters = filt), winslash = "/")
  fls <- fls[order(file.info(fls)$mtime, decreasing = FALSE)]
  if (getInt) {
    tst <- parSapply(parClust, fls, function(fl) {
      rawrr::readFileHeader(fl)$`Creation date`
    }) # Get created time directly from raw file!
    tst <- as.POSIXct.default(tst, tz = Sys.timezone(), format = "%d/%m/%Y %H:%M:%S")
    fls <- fls[order(tst, decreasing = FALSE)]
  }
}
#a <- rawrr::readFileHeader(rawfile = fl)

# Work directory
wd <- unique(dirname(fls))[1]
dtstNm <- gsub(".*/", "", wd)
tst <- try(suppressWarnings(write("Test", paste0(wd, "/test.txt"))), silent = TRUE)
while ("try-error" %in% class(tst)) {
  wd <- rstudioapi::selectDirectory("Choose a work directory where we have write permission!", path = "D:/")
  tst <- try(suppressWarnings(write("Test", paste0(wd, "/test.txt"))), silent = TRUE)
}
unlink(paste0(wd, "/test.txt"))
setwd(wd)
clusterExport(parClust, "wd", envir = environment())

# Convert to mzML
deer <- list()
ParsDirs <- grep("/ThermoRawFileParser",
                 c(list.dirs("C:/ThermoRawFileParser", full.names = TRUE, recursive = FALSE),
                   list.dirs("C:/Program Files", full.names = TRUE, recursive = FALSE),
                   list.dirs(paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads"), full.names = TRUE, recursive = FALSE)),
                 value = TRUE)
# tst <- sapply(ParsDirs, function(x) { sum(c("ThermoRawFileParser.exe", "ThermoFisher.CommonCore.RawFileReader.dll") %in% list.files(x)) == 2 })
# ParsDirs <- ParsDirs[which(tst)]
tst <- sapply(ParsDirs, function(x) { "ThermoRawFileParser" %in% list.files(x) })
ParsDirs <- ParsDirs[which(tst)]
if (!length(ParsDirs)) {
  url <- "https://github.com/compomics/ThermoRawFileParser/archive/refs/heads/master.zip"
  require(curl)
  dstfl <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/ThermoRawFileParser-master.zip")
  curl_download(url, dstfl)
  ParsDirs <- "C:/ThermoRawFileParser/Master"
  if (!dir.exists(ParsDirs)) { dir.create(ParsDirs, recursive = TRUE) }
  unzip(dstfl, exdir = ParsDirs)
}
if (length(ParsDirs) > 1) {
  tst <- sapply(ParsDirs, function(x) { file.info(x)$ctime })
  ParsDirs <- ParsDirs[which(tst == max(tst))[1]]
}
deer$ParsDir <- ParsDirs
MSConvertInst <- ("C:/Program Files/ProteoWizard"%in% list.dirs("C:/Program Files", recursive = FALSE))
if (!MSConvertInst) {
  tst <- grep("ProteoWizard ", list.dirs("C:/Users/Thermo/AppData/Local/Apps", full.names = FALSE), value = TRUE)
  if (length(tst)) {
    MSConvertInst <- TRUE
    MSConvertDir <- paste0("C:/Users/Thermo/AppData/Local/Apps/", tst[1])
    deer$MSConvertDir <- MSConvertDir
  }
} else {
  MSConvertDirs <- grep("/ProteoWizard/ProteoWizard [^/]+",
                        list.dirs("C:/Program Files/ProteoWizard", recursive = FALSE), value = TRUE)
  if (!length(MSConvertDirs)) { MSConvertInst <- FALSE } else {
    #MSConvertDirs <- c("C:/Program Files/ProteoWizard/ProteoWizard 3.0.19172.57d620127",
    #                   "C:/Program Files/ProteoWizard/ProteoWizard 3.0.22099.89b871a",
    #                  "C:/Program Files/ProteoWizard/ProteoWizard 3.0.22317.1e024d4")
    MSConvertVers <- as.data.frame(t(sapply(strsplit(gsub(".*/ProteoWizard ", "", MSConvertDirs), "\\."), unlist)))
    MSConvertVers <- MSConvertVers[order(MSConvertVers$V1, MSConvertVers$V2, MSConvertVers$V3, MSConvertVers$V4, decreasing = TRUE),]
    MSConvertVers <- MSConvertVers[1,]
    MSConvertDir <- paste0("C:/Program Files/ProteoWizard/ProteoWizard ", paste(MSConvertVers, collapse = "."))
    deer$MSConvertDir <- MSConvertDir
  }
}
Convert_mode <- "thermorawfileparser"
zlib <- TRUE
PeakPicking <- TRUE

# Files and Factors
fls0 <- fls
l0 <- l1 <- length(fls)
tst <- sum(grepl("/", fls0))
while ((l0 == l1)&&(tst == l0)) {
  fls0 <- gsub("^[^/]+/", "", fls0)
  l0 <- length(fls0)
  tst <- sum(grepl("/", fls0))
}
fls0 <- gsub("\\.raw$", "", fls0)

# Type of analysis
analysisTypes <- data.frame(Type = c("AC", "GC", "PDE_A", "PDE_G", "Nuc_A", "Nuc_G"),
                            Quantified = c("3'5'cAMP", "3'5'cGMP", "AMP", "GMP", "2'3'cAMP", "2'3'cGMP"))
analysisTypes$Nucleotides <- list(c("3'5'cAMP", "2'3'cAMP", "ATP"),
                                  c("3'5'cGMP", "2'3'cGMP", "GTP"),
                                  c("3'5'cAMP", "AMP"),
                                  c("3'5'cGMP", "GMP"),
                                  c("2'3'cAMP", "RNA"),
                                  c("2'3'cGMP", "RNA"))
msg <- "What type(s) of analysis are we running?"
opt <- sapply(analysisTypes$Type, function(x) { paste(c(x, rep(" ", max(c(100, 250-nchar(x))))), collapse = "") })
analysisType <- analysisTypes$Type[match(dlg_list(opt, opt[1], title = msg, multiple = TRUE)$res, opt)]

# Get Experimental Factors
minFact <- c("Replicate", "Analysis_group", "Samples_group", "Role", "Protein", "Nucleotide")
minFactDesc <- setNames(c("maximum number of replicates",
                          "group of samples to analyze together - which includes controls (buffer blank, standards, no ATP control...) and samples",
                          "group of samples which are replicates of the same condition",
                          "role in the experiment",
                          "name of the protein or protein variant which you are testing",
                          "for standards only: exact nucleotide present in the standard; choose \"none\" for samples"),
                        minFact)
if (file.exists("Factors.RData")) {
  load("Factors.RData")
  Factors <- setNames(tmp$Factors, substr(tmp$Factors, 1, 3))
  Factors <- Factors[which(!is.na(Factors))]
  if (length(Factors)) {
    FactorsLevels <- tmp$Levels[Factors]
    FactorsLevels <- lapply(FactorsLevels, function(x) { grep(" ", x, invert = TRUE, value = TRUE) }) # No spaces allowed in factor levels!!!
    # This is strict because we use spaces in the shiny app to separate multiple levels when entered together
  } else { rm(Factors) }
}
if (!exists("Factors")) {
  Factors <- minFact
  FactorsLevels <- setNames(lapply(Factors, function(x) { c("") }), Factors)
} else { Factors <- unique(c(Factors, minFact)) }
if (!exists("FactorsLevels")) {
  FactorsLevels <- setNames(lapply(Factors, function(x) { c("") }), Factors)
}
w <- which(!Factors %in% names(FactorsLevels))
if (length(w)) {
  FactorsLevels[Factors[w]] <- c()
}
Factors <- Factors[which(!is.na(Factors))]
FactorsLevels <- FactorsLevels[Factors]
FactorsLevels$"Analysis_group" <- unique(c(analysisType, FactorsLevels$"Analysis_group"))
tmp <- FactorsLevels$Replicate
tmp <- suppressWarnings(as.integer(tmp))
tmp <- tmp[which(!is.na(tmp))]
if (!length(tmp)) { tmp <- 1 }
FactorsLevels$Replicate <- 1:max(tmp)
FactorsLevels$Role <- unique(c("Blank", "Buffer_control", "Tag_control", "Standard", "Sample", FactorsLevels$Role))
FactorsLevels$Nucleotide <- unique(c(unlist(analysisTypes$Nucleotides[which(analysisTypes$Type %in% analysisType)]),
                                     "none", FactorsLevels$Nucleotide))
tmp <- FactorsLevels$Nucleotide[which(!FactorsLevels$Nucleotide %in% "none")]
FactorsLevels$Samples_group <- unique(c(paste0("Standard_", tmp), FactorsLevels$Samples_group))

#rm(Factors, FactorsLevels)
intFact <- c("Replicate", "Isobaric.set")
dfltInt <- c(1, 1)
if (!exists("blnksPat")) { blnksPat <<- "^blank(_[0-9]+)?$" }
appNm <- paste0(dtstNm, " - Experimental Factors")
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  titlePanel(tag("u", "Experimental structure editor"),
             appNm),
  h2(dtstNm), 
  br(),
  tags$hr(style = "border-color: black;"),
  h3("Select MS files to analyze and plot, in the desired order!"),
  textInput("negStrng",
            "Blanks identification: negative filter (regular expression) to apply",
            blnksPat),
  h5(em("(case-sensitive; use ^ and $ for a string's beginning and end, respectively, [0-9] for any number...")),
  h5(em("see "),
     #a("this link", href = "https://www.tutorialspoint.com/tcl-tk/tcl_regular_expressions.htm"), # Disappointingly, the link worked but clicking it killed the app, which, well, is kinda stupid...
     a("https://www.tutorialspoint.com/tcl-tk/tcl_regular_expressions.htm"),
     em(" for full regex syntax)")),
  span(uiOutput("Order")),
  span(uiOutput("ordrMsg"), style = "color:grey"),
  br(),
  tags$hr(style = "border-color: black;"),
  h3("Enter the names and levels of Experiment Factors."),
  h4("These are, so to speak, the columns in the experiment map. Each may have different levels, for instance, Replicate: 1-2, Protein: MBP;TIR1;SLY1..."),
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
      h5(em(" - use underscores, no spaces or hyphens within level names")),
      h5(em(" - use spaces as separators to add multiple levels at a time")),
      br(),
      span(uiOutput("fctrMsg"), style = "color:red"),
      h3("Factor levels:"),
      uiOutput("Factors"),
      br(),
      br(),
      actionButton("saveBtn", "Save")
    )
  )
)
# Pre-processing for server
# - Default order
if (!exists("flsOrd0_dflt")) { flsOrd0_dflt <- fls0 }
flsOrd0 <- c(flsOrd0_dflt,
             fls0[which(!fls0 %in% flsOrd0_dflt)])
# - Default blanks pattern
grpOK <- FALSE
if (nchar(blnksPat)) {
  tst <- try(grep(blnksPat, fls0, value = TRUE), silent = TRUE)
  grpOK <- (!"try-error" %in% class(tst))
}
if (grpOK) {
  blnkFls0 <- tst
  flsOrd0 <- flsOrd0[which(!flsOrd0 %in% blnkFls0)]
  flsOrd0 <- c(flsOrd0,
               fls0[which(!fls0 %in% c(flsOrd0, blnkFls0))])
} else {
  blnkFls0 <- c()
}
server <- function(input, output, session) {
  # Initialize reactive variables
  blanksPAT <- reactiveVal(blnksPat) # Regex pattern
  smplFILES <- reactiveVal(flsOrd0) # Files defined as sample files based on the regex
  flsORD <- reactiveVal(flsOrd0) # Order of sample files; may exclude some
  FACT <- reactiveVal(Factors)
  FACTLevels <- reactiveVal(FactorsLevels)
  #
  output$ordrMsg <- renderUI({ em(" ") })
  output$Order <- renderUI(selectInput("Order",
                                       "",
                                       flsOrd0,
                                       flsOrd0,
                                       TRUE,
                                       TRUE,
                                       width = "1200px"))
  updtOpt <- function(reactive = TRUE) {
    if (reactive) {
      blnksPat <<- blanksPAT()
      fls0_a <- flsORD()
    } else {
      fls0_a <- flsOrd0
    }
    fls0_a <- c(fls0_a,
                fls0[which(!flsOrd0 %in% fls0_a)])
    grpOK <- FALSE
    if (nchar(blnksPat)) {
      tst <- try(grep(blnksPat, fls0_a, value = TRUE), silent = TRUE)
      grpOK <- (!"try-error" %in% class(tst))
    }
    if (grpOK) {
      blnk0_b <- tst
      myFls0_b <- fls0_a[which(!fls0_a %in% blnk0_b)]
    } else {
      myFls0_b <- fls0_a
      blnk0_b <- c()
    }
    if (reactive) {
      smplFILES(myFls0_b)
      slct <- input$Order[which(input$Order %in% myFls0_b)]
    } else{
      slct <- flsOrd0[which(flsOrd0 %in% myFls0_b)]
    }
    updateSelectInput(inputId = "Order",
                      choices = myFls0_b,
                      selected = slct)
  }
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
  output$minFact <- renderUI({ HTML(paste0(" - ", minFact, " => ", minFactDesc, collapse = "<br>")) })
  output$fctrMsg <- renderUI({ em(" ") })
  # Initialize
  output$Factors <- updtFactUI(FALSE)
  #
  # Observers for files selection
  observeEvent(input$Order, {
    flsORD(input$Order)
    flsOrd0_dflt <<- flsOrd0 <<- input$Order
    n <- length(fls0) - length(flsOrd0) - length(blnkFls0)
    if (n) {
      output$ordrMsg <- renderUI({ em(paste0("Files currectly not selected: ", n)) })
    } else {
      output$ordrMsg <- renderUI({ em(" ") })
    }
  })
  observeEvent(input$negStrng, {
    blanksPAT(input$negStrng)
    updtOpt()
  })
  #
  # Observers for already extent factors
  #  - Add new Factor level
  sapply(Factors, function(Fact) {
    observeEvent({ input[[paste0(Fact, "_levAdd")]] }, {
      vals <- input[[paste0(Fact, "_lev")]]
      vals <- vals[which(!is.na(vals))]
      if (length(vals)) {
        if (is.character(vals)) { vals <- unlist(strsplit(vals, " ")) }
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
    observeEvent({ input[[paste0(Fact, "_levRmv")]] }, {
      vals <- input[[paste0(Fact, "_lev")]]
      vals <- vals[which(!is.na(vals))]
      if (length(vals)) {
        if (is.character(vals)) { vals <- unlist(strsplit(vals, " ")) }
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
      #if ((nchar(Fact) < 3)||(grepl("^[0-9]", Fact))||(substr(Fact, 1, 3) %in% substr(FACT(), 1, 3))) {
      #  msg <- "Invalid Factor name! Must be at least 3 characters long and start with a capital letter! The first 3 characters must be unique to this Factor!"
      #} else {
        #Fact <- paste0(toupper(substr(Fact, 1, 1)), tolower(substr(Fact, 2, nchar(Fact))))
        if (!Fact %in% FACT()) {
          FACT(c(FACT(), Fact))
          tmp <- FACTLevels()
          tmp[Fact] <- list(NULL)
          FACTLevels(tmp)
          # Also, ABSOLUTELY crucial: create new observers for level addition/removal!!!
          observeEvent({ input[[paste0(Fact, "_levAdd")]] }, {
            vals <- input[[paste0(Fact, "_lev")]]
            vals <- vals[which(!is.na(vals))]
            if (length(vals)) {
              if (is.character(vals)) { vals <- unlist(strsplit(vals, " ")) }
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
          observeEvent({ input[[paste0(Fact, "_levRmv")]] }, {
            vals <- input[[paste0(Fact, "_lev")]]
            vals <- vals[which(!is.na(vals))]
            if (length(vals)) {
              if (is.character(vals)) { vals <- unlist(strsplit(vals, " ")) }
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
      #}
    })
    output$fctrMsg <- renderUI({ em(msg) })
    output$Factors <- updtFactUI()
  })
  # Remove extant Factor
  observeEvent(input$rmvFactor, {
    msg <- " "
    tmp <- gsub("[^A-Z,a-z,0-9]", "\\.", input$nuFact)
    if ((nchar(tmp) < 3)||(grepl("^[0-9]", tmp))) {
      msg <- "Invalid Factor name! Must be at least 3 characters long and start with a capital letter!"
    } else {
      #tmp <- paste0(toupper(substr(tmp, 1, 1)), tolower(substr(tmp, 2, nchar(tmp))))
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
    output$fctrMsg <- renderUI({ em(msg) })
    output$Factors <- updtFactUI()
  })
  #
  #
  observeEvent(input$saveBtn, {
    Factors <<- FACT()
    FactorsLevels <<- FACTLevels()
    flsOrd0_dflt <<- flsOrd0 <<- flsORD()
    grpOK <- FALSE
    if (nchar(blnksPat)) {
      tst <- try(grep(blnksPat, fls0, value = TRUE), silent = TRUE)
      grpOK <- (!"try-error" %in% class(tst))
    }
    if (grpOK) { blnkFls0 <<- tst } else {
      stop("No blank files selected!!!")
    }
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
eval(parse(text = runApp))
FactorsLevels <- setNames(lapply(Factors, function(fct) {
  x <- FactorsLevels[[fct]]
  x[which(!is.na(x))]
}), Factors)
Factors <- Factors[which(sapply(FactorsLevels[Factors], length) > 0)]
stopifnot(sum(!minFact %in% Factors) == 0)
Factors <- c(Factors[which(Factors != "Replicate")], "Replicate")
FactorsLevels <- FactorsLevels[Factors]
FactorsLevels$Protein <- unique(c(NA, FactorsLevels$Protein))
tmp <- list(Factors = Factors, Levels = FactorsLevels)
save(tmp, file = "Factors.RData")
for (fct in Factors) {
  assign(fct, FactorsLevels[[fct]])
}
slctFls <- fls[match(flsOrd0, fls0)]
slctFls0 <- flsOrd0

# Files map
ExpMapNm <- "Experiment map"
ExpMapPath <- paste0(wd, "/", ExpMapNm, ".csv")
Factors2 <- Factors[which(!Factors %in% minFact)]
if (file.exists(ExpMapPath)) {
  ExpMap2 <- read.csv(ExpMapPath, check.names = FALSE)
  colnames(ExpMap2)[which(colnames(ExpMap2) == "Sample.name")] <- "Sample name" # Backwards compatibility
  reLoad <- sum(!slctFls %in% ExpMap2$`MS raw file`) == 0
  if ((reLoad)&&(exists("ExpMap"))) {
    reLoad <- c(TRUE, FALSE)[match(dlg_message("Reload Experiment map from disk?\n(you will still be able to edit it)", "yesno")$res, c("yes", "no"))]
  }
  if (reLoad) {
    ExpMap <- ExpMap2[match(slctFls, ExpMap2$`MS raw file`), ]
  }
}
if (!exists("ExpMap")) {
  ExpMap <- data.frame("MS raw file" = slctFls,
                       "MS raw file name" = slctFls0,
                       "Role" = "Sample",
                       "Replicate" = "?",
                       check.names = FALSE)
}
for (Fact in Factors2) { ExpMap[[Fact]] <- "?" }
# Edit map
ExpData <- ExpMap#read.csv(ExpMapPath, check.names = FALSE)
for (Fact in Factors[which(!Factors %in% colnames(ExpData))]) { ExpData[[Fact]] <- "?" }
tst <- sapply(FactorsLevels, length)
Fact1 <- Factors[which(tst == 1)]
Fact2 <- Factors[which(tst > 1)]
nr <- nrow(ExpData)
rws <- seq_len(nr)
AllIDs <- setNames(lapply(Fact2, function(x) { paste0(x, "___", rws)} ), Fact2)
ALLIDS <- setNames(unlist(AllIDs), NULL)
Fact2IDs <- setNames(lapply(Fact2, function(x) { paste0(x, "___", rws)} ), Fact2)
facLevels2 <- lapply(FactorsLevels, function(x) {
  unique(c(x, NA))
})
facLevels2$Use <- c(TRUE, FALSE) # Not used, but just in case
appNm <- paste0(dtstNm, " - Experiment map")
if (length(Fact1)) { for (fct in Fact1) { #fct <- Fact1[1]
  ExpData[[fct]] <- FactorsLevels[[fct]]
} }
# Pre-processing for server:
xpDat <- ExpData[, c("MS raw file", Factors)]
# - Simplify file names
nc <- nchar(xpDat$"MS raw file")
nr <- nrow(xpDat)
m <- min(nc) # Should always be larger or equal to 2 (much larger)
tst <- sapply(2:m, function(x) {
  length(unique(substr(xpDat$"MS raw file", 1, x)))
}) == 1
if (sum(!tst)) {
  w <- which(!tst)[1]+1
  xpDat$"MS raw file" <- substr(xpDat$"MS raw file", w, nc)
}
# - Original table column widths
wTest0 <- setNames(sapply(colnames(ExpData), function(k) { #k <- colnames(ExpData)[1]
  tst <- k %in% Fact2
  if (k == "MS raw file") {
    x <- max(c(nchar(k),
               nchar(as.character(xpDat[[k]])) + 5 + 3*tst), na.rm = TRUE)
  } else {
    x <- max(c(nchar(k),
               nchar(as.character(ExpData[[k]])) + 5 + 3*tst), na.rm = TRUE)
  }
  if (tst) {
    x <- max(c(x, nchar(FactorsLevels[[k]]) + 6), na.rm = TRUE)
  }
  x <- x*7
  if (is.na(x)) { x <- 15 } else { x <- max(c(ceiling(x/10)*10, 30)) }
  return(x)
}), colnames(ExpData))
#xpDat$Reference <- shinyCheckInput(ExpData$Reference, "Reference")
# - Edit table
kol <- c()
for (fct in Fact2) { #fct <- Fact2[1]
  IDs <- Fact2IDs[[fct]]
  lvls <- FactorsLevels[[fct]]
  lvls2 <- lvls[which(!is.na(lvls))]
  xpDat[[fct]] <- shinySelectInput(ExpData[[fct]],
                                   fct,
                                   lvls,
                                   paste0(30*max(c(nchar(c(as.character(lvls2), fct)), 2)), "px")#,
                                   #c(TRUE, FALSE)[(fct %in% c("Analysis_group", "Role", "Replicate", "Protein"))+1]
  )
  fdNm <- paste0(fct, "___FD")
  wTest0[fdNm] <- 15
  xpDat[[fdNm]] <- shinyFDInput(fct, nr, TRUE, paste0(wTest0[fdNm], "px"))
  kol <- c(kol, fct, fdNm)
  if (fct == "Replicate") {
    incrNm <- paste0(fct, "___INCR")
    wTest0[incrNm] <- 15
    xpDat[[incrNm]] <- shinyPlusFD(fct, nr, width = paste0(wTest0[incrNm], "px"))
    kol <- c(kol, incrNm)
  }
}
kol2 <- colnames(xpDat)[which(!colnames(xpDat) %in% c(kol))]
xpDat <- xpDat[, c(kol2, kol)]
# - Estimate table column widths
wTest1 <- sapply(colnames(xpDat), function(k) { #k <- colnames(xpDat)[1]
  if (k == "Parent sample") { k <- "MQ.Exp" }
  if (k %in% names(wTest0)) { x <- wTest0[k] } else { x <- 30 }
  return(x)
})
wTest2 <- sum(wTest1)
wTest1 <- paste0(as.character(wTest1), "px")
wTest1 <- aggregate((1:length(wTest1))-1, list(wTest1), c)
wTest1 <- apply(wTest1, 1, function(x) {
  list(width = x[[1]],
       targets = x[[2]],
       names = colnames(xpDat)[x[[2]]+1])
})
oTHerz <- which(colnames(xpDat) != "Parent sample")
g <- grep("___((FD)|(INCR))$", colnames(xpDat))
colnames(xpDat)[g] <- "" 
edith <- list(target = "column",
              disable = list(columns = c(0, match(Fact1, colnames(xpDat))-1)))
tmp <- c(0:(ncol(xpDat)-1))
tmp <- tmp[which(!tmp %in% edith$disable$columns)]
edith$enable <- list(columns = tmp)
L <- length(FactorsLevels[["Replicate"]])
dflt_Rpl <- 1:(nr+L) %% L
dflt_Rpl[which(dflt_Rpl == 0)] <- L
dflt_Rpl <- FactorsLevels[["Replicate"]][dflt_Rpl] # Easy template
if (exists("ExpData3")) { rm(ExpData3) }
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  # Dummy, hidden div to load arrow-down icon
  tags$div(
    class = "container",
    style = "display: none;",
    tags$div(
      style = "margin-top: 50xp;",
      actionButton(
        "add_thing",
        label = "do it",
        class = "btn-success",
        icon = icon("arrow-down")
      )
    )
  ),
  #
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", ExpMapNm),
             appNm),
  h2(dtstNm), 
  h3("Define the experiment's structure:"),
  br(),
  actionButton("saveBtn", "Save"),
  withSpinner(DT::DTOutput("ExpTbl", width = wTest2))
)
server <- function(input, output, session) {
  ExpData3 <- ExpData # Output table
  output$ExpTbl <- renderDT({ xpDat },
                            FALSE,
                            escape = FALSE,
                            selection = "none",
                            editable = edith,
                            rownames = FALSE,
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
  #
  # Fill down
  sapply(1:length(ALLIDS), function(x) {
    id1 <- ALLIDS[x]
    id2 <- paste0(ALLIDS[x], "___FD")
    tmp <- unlist(strsplit(id1, "___"))
    fct <- tmp[[1]]
    i <- as.integer(tmp[[2]])
    # observeEvent(input[[id1]],
    observeEvent(input[[id2]],
                 {
                   if (i < nr) {
                     x <- input[[id1]]
                     tp <- typeof(facLevels2[[fct]])
                     if (typeof(x) == tp) { x <- get(paste0("as.", tp))(x) }
                     for (k in (i+1):nr) {
                       idK <- paste0(fct, "___", as.character(k))
                      # if (fct == "Use") {
                      #   updateCheckboxInput(session, idK, NULL, x)
                      # } else {
                         updateSelectInput(session, idK, NULL, facLevels2[[fct]], x)
                      # }
                     }
                   }
                 }
    )
  })
  # Incremental fill-down for replicates
  sapply(rws, function(i) {
    if (i < nr) {
      id1 <- paste0("Replicate___", i)
      id2 <- paste0("Replicate___", i, "___INCR")
      observeEvent(input[[id2]],
                   {
                     if (i < nr) {
                       x <- input[[id1]]
                       tp <- typeof(facLevels2[[fct]])
                       if (typeof(x) == tp) { x <- get(paste0("as.", tp))(x) }
                       rplRg <- (i+1):nr
                       l <- length(rplRg)
                       m <- match(x, dflt_Rpl)+1
                       rplVal <- dflt_Rpl[m:(m+l-1)]
                       for (k in rplRg) {
                         y <- rplVal[k-i]
                         idK <- paste0("Replicate___", as.character(k))
                         updateSelectInput(session, idK, NULL, facLevels2[["Replicate"]], y)
                       }
                     }
                   }
      )
    }
  })
  # Manual cell edit (sample names)
  observeEvent(input$ExpTbl_cell_edit, {
    ExpData3[input$ExpTbl_cell_edit$row,
             input$ExpTbl_cell_edit$col+1] <- input$ExpTbl_cell_edit$value
  })
  # Save
  observeEvent(input$saveBtn, {
    for (fct in #c(
         Fact2#, "Use")
         ) {
      ExpData3[[fct]] <- sapply(rws, function(i) {
        input[[paste0(fct, "___", i)]]
      })
      #if (fct == "Use") {
      #  ExpData3[[fct]] <- as.logical(ExpData3[[fct]])
      #} else {
        typ <- typeof(FactorsLevels[[fct]])
        ExpData3[[fct]] <- get(paste0("as.", typ))(ExpData3[[fct]])
        # Consider here detecting if a factor needs conversion to numeric/integer...
      #}
    }
    assign("ExpData3", ExpData3, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("ExpData3"))) {
  eval(parse(text = runApp))
  runKount <- runKount+1
}
runKount <- 0
tmpTbl <- ExpMap <- ExpData3
tst <- lapply(colnames(tmpTbl), function(x) { typeof(tmpTbl[[x]]) })
w <- which(tst == "list")
if (length(w)) { for (i in w) { tmpTbl[[i]] <- sapply(tmpTbl[[i]], paste, collapse = ";") }}
tst <- try(write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE), silent = TRUE)
if ("try-error" %in% class(tst)) {
  dlg_message(paste0("File \"", ExpMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE)
}
ExpMap$"MS raw file name" <- slctFls0[match(ExpMap$`MS raw file`, slctFls)]

#
Quant <- TRUE
blnkFls <- fls[which(fls0 %in% blnkFls0)]
refFls <- ExpMap$`MS raw file`[which(ExpMap$Role == "Buffer_control")]
if (!length(refFls)) {
  warning("No Buffer_control references defined: we will not be able to quantify!")
  Quant <- FALSE
}

# Files ranges
rgL <- 24
ctrlRanges <- flRanges <- list()
for (grp in Analysis_group) { #grp <- Analysis_group[1]
  w0 <- which((ExpMap$"Analysis_group" == grp)
              &(!ExpMap$Role %in% c("Blank", "Sample")))
  l0 <- length(w0)
  if (l0) { ctrlRanges[[grp]] <- w0 }
  for (prot in Protein) { #prot <- Protein[1]
    if (is.na(prot)) {
      w1 <- which((ExpMap$"Analysis_group" == grp)
                  &(!ExpMap$Role %in% c("Buffer_control", "Tag_control"))
                  &(is.na(ExpMap$Protein)))
    } else {
      w1 <- which((ExpMap$"Analysis_group" == grp)
                  &(!ExpMap$Role %in% c("Buffer_control", "Tag_control"))
                  &(ExpMap$Protein == prot))
    }
    l1 <- length(w1)
    if (l1) {
      N <- ceiling(l1/rgL)
      for (i in 1:N) {
        y <- (1:rgL)+(i-1)*rgL
        y <- y[which(y <= l1)]
        y <- y[which(!w1[y] %in% ctrlRanges[[grp]])]
        flRanges[[paste0(grp, " ", prot, " ", i)]] <- w1[y]
      }
    }
  }
}
stopifnot(identical(setNames(sort(c(unlist(flRanges), unlist(ctrlRanges))), NULL),
                    which(ExpMap$Role != "Blank")))

# Convert to mzMLs
mzMLs <- gsub("\\.raw$", ".mzML", slctFls, ignore.case = TRUE)
pressFls <- gsub("\\.mzML$", ".csv", mzMLs)
clusterExport(parClust, list("slctFls", "mzMLs", "pressFls"), envir = environment())
#w <- 1:length(mzMLs)
w <- which((!file.exists(mzMLs))|(file.size(mzMLs) <= 2000))
if (length(w)) {
  if (tolower(Convert_mode) == "thermorawfileparser") { # Mode 1: using ThermoRawFileParser
    clusterExport(parClust, list("wd", "deer", "zlib", "PeakPicking"), envir = environment())
    tst <- parSapply(parClust, w, function(i) { #i <- w[1]
      cmd <- paste0("ThermoRawFileParser -i=\"",
                    slctFls[i], "\" -b=\"", gsub(".*/", paste0(wd, "/"), mzMLs[i]), "\" -f=2 -a",
                    c(" -z", "")[zlib+1], c(" -p", "")[PeakPicking+1])
      cmd <- c("C:", paste0("cd \"", deer$ParsDir, "\""), cmd)
      write(cmd, paste0(wd, "/tmp", i, ".bat"))
      #cat(cmd)
      #system(cmd)
      cmd2 <- paste0("\"", wd, "/tmp", i, ".bat\"") # I have to go through an intermediate batch file to run cmd,
      # because it is one of those which works in Windows command line but not when passed to system() or shell() in R.
      #cat(cmd2)
      #writeClipboard(cmd)
      shell(cmd2)
      unlink(paste0(wd, "/tmp", i, ".bat"))
    })
  }
  if (tolower(Convert_mode) == "msconvert") { # Mode 2: using msconvert
    write(slctFls[w], file = paste0(wd, "/tmp_MS_files.txt"))
    precRecal <- FALSE
    cmd <- paste0("\"", deer$MSConvertDir, "/msconvert.exe\" -f \"", wd, "/tmp_MS_files.txt\" -o \"",
                  wd, "\" --mzML --64 -z --filter \"peakPicking true 1-\"")
    #cat(cmd)
    system(cmd)
    unlink(paste0(wd, "/tmp_MS_files.txt"))
  }
}

# Getting the pressure profile is unfortunately not doable using rawrr
# This must be done from mzMLs
#tst <- grabMzmlData(mzMLs[1], "everything")
#View(tst$metadata)
pressScript <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/mzML_pressure_to_csv.py")
if (!file.exists(pressScript)) {
  download.file("https://gist.githubusercontent.com/caetera/0921b33f0c6201a538436906cc965fff/raw/d1af134fe228ce6a23be3e5bc3c49d20a8447ab2/mzML_pressure_to_csv.py",
                pressScript)
}
stopifnot(file.exists(pressScript))
# Create temporary python pressure script edited to create a more predictably named output
tmp <- suppressWarnings(readLines(pressScript))
w <- which(tmp == "        csv_path = file.replace('.mzML', f'_{convert_to_safe_filename(chrom_name)}.csv')")
tmp[w] <- "        csv_path = file.replace('.mzML', f'.csv')"
pressScript2 <- paste0(wd, "/pressScript.py")
write(tmp, pressScript2)
w <- which(file.exists(mzMLs))
clusterExport(parClust, list("pressScript2", "mzMLs"), envir = environment())
tst <- parSapply(parClust, w, function(i) { #i <- w[1]
  cmd <- paste0("python \"", pressScript2, "\" ", paste0("\"", mzMLs[i], "\"", collapse = " "))
  #cat(cmd)
  system(cmd)
})
#
cmd <- paste0("open \"", wd, "\"")
#cat(cmd)
#system(cmd)
tstPress <- file.exists(pressFls)
#stopifnot(length(w) == 0)
if (sum(!tstPress)) {
  warning(paste0("Couldn't get pressure profile from the following file(s):",
             paste0("\n - ", slctFls0[which(!tstPress)], collapse = ""), "\n"))
}
unlink(pressScript2) # Remove temporary python script 
#
clusterExport(parClust, "slctFls", envir = environment())
press <- setNames(parLapply(parClust, which(tstPress), function(i) {
  x <- data.table::fread(pressFls[i], integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  fl <- slctFls[i]
  colnames(x) <- c("Retention time", "Pressure")
  x$"Raw file" <- factor(fl, levels = slctFls)
  x$"Raw file name" <- factor(cleanRawNm(basename(fl)), levels = cleanRawNm(basename(slctFls)))
  return(x)
}), slctFls[which(tstPress)])
allChroms <- list(Pressure = press)
#allChroms$Pressure <- press
#
#chromtypes <- setNames(c("tic", "bpc", "xic"), c("TIC", "Base peak", "XIC"))
allFls <- c(blnkFls, slctFls)
if (getInt) {
  tol <- as.numeric(dlg_input("Enter mass tolerance (ppm)", 20)$res)
  clusterCall(parClust, function() library(rawrr))
  MS2s <- parLapply(parClust, allFls, function(fl) { #fl <- slctFls[1]
    x <- rawrr::readIndex(fl)
    x <- x[which(x$MSOrder == "Ms2"),]
    return(unique(gsub("@.*", "", unique(x$scanType))))
  })
  ms2s <- unique(unlist(MS2s))
  if (length(ms2s) == 1) {
    cat(paste0("Single MS2 filter detected: \"", ms2s, "\"\n"))
  } else {
    tmp <- sapply(ms2s, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") })
    tmp <- dlg_list(tmp, tmp[1], TRUE, "Select filter(s) for which you want to get an XIC")$res
    ms2s <- gsub(" *$", "", tmp)
  }
  dflts <- c("134.0472", "150.0421", NA)
  masses <- setNames(lapply(ms2s, function(ms2) { #ms2 <- ms2s[1]
    msg <- paste0(ms2, " XIC: enter M/Z to extract (if multiple M/Zs, separate them with \" \")")
    m <- match(ms2, c("FTMS - p ESI Full ms2 328.0452",
                      "FTMS - p ESI Full ms2 344.0402"))
    if (is.na(m)) { m <- 3 }
    dflt <- dflts[m]
    mass <- as.numeric(unlist(strsplit(dlg_input(msg, dflt)$res, " ")))
    return(mass[which(!is.na(mass))])
  }), ms2s)
  # TIC
  tic <- setNames(parLapply(parClust, allFls, ticFun), slctFls)
  w <- sapply(tic, function(x) { x$Outcome })
  tic <- lapply(tic[w], function(x) { x$Output })
  if (length(tic)) { allChroms$TIC <- tic }
  # BPC
  bpc <- setNames(parLapply(parClust, allFls, bpcFun), allFls)
  w <- sapply(bpc, function(x) { x$Outcome })
  bpc <- lapply(bpc[w], function(x) { x$Output })
  if (length(bpc)) { allChroms$BPC <- bpc }
  # XICs
  clusterExport(parClust, list("ms2s", "tol", "masses", "xicFun"), envir = environment())
  fltms2 <- listMelt(masses, names(masses))
  clusterExport(parClust, list("tol", "xicFun"), envir = environment())
  xicNms <- do.call(paste, c(fltms2, sep = " for "))
  xic <- setNames(apply(fltms2, 1, function(x) {
    mass <- x[[1]]
    ms2 <- x[[2]]
    clusterExport(parClust, list("ms2", "mass"), envir = environment())
    setNames(parLapply(parClust, allFls, function(fl) { xicFun(fl, mass, tol, ms2) }), allFls)
  }), xicNms)
  w <- setNames(lapply(xicNms, function(x) { #x <- names(xic)[1]
    which(sapply(xic[[x]], function(y) { y$Outcome }))
  }), xicNms)
  xic <- setNames(lapply(xicNms, function(x) { #x <- names(xic)[1]
    setNames(lapply(allFls[w[[x]]], function(fl) { xic[[x]][[fl]]$Output }), allFls[w[[x]]])
  }), xicNms)
  w <- which(setNames(lapply(xicNms, length), xicNms) > 0)
  xic <- xic[w]
  for (nm in names(xic)) { allChroms[[nm]] <- xic[[nm]] }
}
#
for (nm in names(allChroms)) { #nm <- names(allChroms)[1]
  allChroms[[nm]] <- plyr::rbind.fill(allChroms[[nm]])
  if (nm == "Pressure") { #nm <- "Pressure"
    tmp <- as.data.table(allChroms[[nm]])
    tmp$Bin <- round(as.numeric(tmp$`Retention time`, 1))
    tmp <- tmp[, list(Pressure = mean(Pressure)),
               by = list(`Raw file` = `Raw file`,
                         `Raw file name` = `Raw file name`,
                         `Retention time` = Bin)]
    allChroms[[nm]] <- as.data.frame(tmp)
  }
}
#

# Define retention time ranges
fullRTRange <- unlist(lapply(allChroms, function(x) { x$`Retention time` }))
fullRTRange <- c(min(fullRTRange), max(fullRTRange))
getRTRange <- TRUE; kount <- 1
dfltRanges <- data.frame(Start = c(fullRTRange[1], 0, 2, 4.5, 0, 3, 5),
                         End = c(fullRTRange[2], 1, 3.5, 6, 1, 4.5, 6),
                         Type = c("Full", rep("Sub-range", 6)),
                         Nucleotide = c(NA,
                                        "ATP", "2'3'cAMP", "3'5'cAMP",
                                        "GTP", "2'3'cGMP", "3'5'cGMP"))
dfltRanges <- dfltRanges[c(1, which(dfltRanges$Nucleotide %in% Nucleotide)),]
dfltRanges$"Analysis_group" <- lapply(dfltRanges$Nucleotide, function(x) {
  rs <- Analysis_group[which(Analysis_group %in% analysisTypes$Type[which(sapply(analysisTypes$Nucleotides, function(y) { x %in% y }))])]
  if (!length(rs)) { rs <- Analysis_group }
  return(rs)
})
kount <- 0
msg <- paste0("Extract ", c("", "another ")[(kount > 0)+1], "custom RT sub-range?")
opt <- c("No                                                                           ",
         "Yes                                                                          ")
getRTRange <- c(FALSE, TRUE)[match(dlg_list(opt, opt[1], title = msg)$res, opt)]
while (getRTRange) {
  if (length(Analysis_group) > 1) {
    dflt <- c(Analysis_group, "All")
    dflt <- setNames(sapply(dflt, function(x) {
      paste0(c(x, rep(" ", 250-nchar(x))), collapse = "")
    }), dflt)
    wh <- dlg_list(dflt, dflt[1], title = "To which group(s) will this range apply?")$res
    wh <- names(dflt)[match(wh, dflt)]
    if (wh == "All") { wh <- Analysis_group }
  } else { wh <- Analysis_group }
  dflt <- fullRTRange
  if ((exists("RTRng"))&&(!is.null(RTRng))&&(is.numeric(RTRng))&&(length(RTRng) == 2)) { dflt <- RTRng } else { dflt <- fullRTRange }
  RTRng <- as.numeric(unlist(strsplit(dlg_input("Enter min. and max. RT separated by \"-\"",
                                                paste(dflt, collapse = "-"))$res, " *- *")))
  RTRng <- RTRng[which(!is.na(RTRng))]
  if (length(RTRng) != 2) {
    warning("Invalid RT range! Cancelling...")
  } else {
    RTRng <- sort(RTRng)
    if (RTRng[1] < fullRTRange[1]) {
      warning("Too low RT range!")
      RTRng[1] <- fullRTRange[1]
    }
    if (RTRng[2] > fullRTRange[2]) {
      warning("Too high RT range!")
      RTRng[2] <- fullRTRange[2]
    }
    RTRng <- data.frame(Start = RTRng[1],
                        End = RTRng[2],
                        Type = "Sub-range")
    if ((RTRng$Start == fullRTRange[1])&&(RTRng$End == fullRTRange[2])) { RTRng$Type <- "Full" }
    RTRng$"Analysis_group" <- list(wh)
    msg <- paste0("Which nucleotide is relevant for this range? (RT = ", RTRng$Start, "-", RTRng$End, ")")
    opt <- sapply(Nucleotide, function(x) { paste(c(x, rep(" ", max(c(100, 250-nchar(x))))), collapse = "") })
    RTRng$Nucleotide <- Nucleotide[match(dlg_list(opt, opt[1], title = msg)$res, opt)]
    dfltRanges <- rbind(dfltRanges, RTRng)
    kount <- kount+1
    msg <- paste0("Extract ", c("", "another ")[(kount > 0)+1], "custom RT sub-range?")
    getRTRange <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  }
}
rtRanges <- dfltRanges
dfltRanges$"Analysis_group" <- sapply(dfltRanges$"Analysis_group", paste, collapse = ",")
tst <- do.call(paste, c(dfltRanges, sep = "_"))
tst <- aggregate(1:length(tst), list(tst), min)
tst <- tst$x[order(tst$x)]
rtRanges <- rtRanges[tst,]
#
lnTypes <- setNames(rep("solid", length(Role)), Role)
lnTypes["Buffer_control"] <- "dotted"
lnTypes["Tag_control"] <- "dotdash"
lnTypes["Standard"] <- "dashed"

# Slightly edited map (NA is not nice in the plot facet labels)
flsOrd <- ExpMap$`MS raw file`[match(flsOrd0, ExpMap$`MS raw file name`)]
ExpMap2 <- ExpMap
for (k in colnames(ExpMap2)) {
  ExpMap2[[k]] <- as.character(ExpMap2[[k]])
  w <- which((is.na(ExpMap2[[k]]))|(ExpMap2[[k]] == "NA"))
  ExpMap2[w, k] <- "_"
  if (k == "MS raw file name") {
    ExpMap2[[k]] <- factor(ExpMap2[[k]], levels = flsOrd0)
  }
  if (k == "MS raw file") {
    ExpMap2[[k]] <- factor(ExpMap2[[k]], levels = flsOrd)
  }
  if (k %in% Factors) {
    ExpMap2[[k]] <- factor(ExpMap2[[k]], levels = c(FactorsLevels[[k]], "_"))
  }
}
fact <- Factors[which(!Factors %in% c("Replicate", "Samples_group", "Samples_group2"))]
fact <- fact[which(sapply(fact, function(x) {
  length(unique(ExpMap2[[x]]))
}) > 1)]
if (length(fact) > 1) {
  # Check for synonymous factors
  tmp <- ExpMap2[, fact]
  for (fct in fact) {
    tmp[[fct]] <- as.numeric(as.factor(tmp[[fct]]))
  }
  comb <- gtools::combinations(length(fact), 2, fact)
  comb <- as.data.frame(comb)
  tst <- apply(comb, 1, function(x) {
    identical(tmp[[x[1]]], tmp[[x[2]]])
  })
  w <- which(tst)
  if (length(w)) {
    fct2Remov <- unique(comb[w, 2])
    fact <- fact[which(!fact %in% fct2Remov)]
  }
}
if (!length(fact)) {
  fact <- "Protein"
}
fact <- unique(c(fact, "Role"))
ExpMap2$Samples_group2 <- do.call(paste, c(ExpMap2[, fact, drop = FALSE], sep = " "))
tst <- aggregate(ExpMap2$Samples_group, list(ExpMap2$Samples_group2), function(x) { length(unique(x)) })
if (max(tst$x) > 1) {
  fact <- unique(c(fact, "Samples_group"))
  ExpMap2$Samples_group2 <- do.call(paste, c(ExpMap2[, fact, drop = FALSE], sep = " "))
}
Samples_group2 <- unique(ExpMap2$Samples_group2)
FactorsLevels[["Samples_group2"]] <- Samples_group2
Factors <- unique(c("Samples_group2", Factors))
ExpMap2$Samples_group2 <- factor(ExpMap2$Samples_group2, levels = FactorsLevels[["Samples_group2"]])

# Peak extract
propFlt <- function(x,
                    p = 0.1, # Proportion of values to average
                    upper = TRUE, # If TRUE, the proportion applies to the upper tail of the data, if not to the lower tail.
                    average = TRUE # If TRUE, returns mean of filtered data, if FALSE returns data 
                    ) {
  x <- is.all.good(x) 
  l <- length(x)
  x <- sort(x, decreasing = upper)
  x <- x[1:round(l*p)]
  if (average) { x <- mean(x) }
  return(x)
}
peakXtract <- function(fileName, # File name
                       data = chr,
                       WIP = FALSE, # Work-in-Progress? If so we will be plotting stuff to check behaviour
                       binSz = 10, # Default RT bin size (s)
                       targRg,
                       bgMode = "loess") { #fileName <- f[1] #fileName <- f[2] #fileName <- f[3] #fileName <- f[4] #fileName <- f[5] #fileName <- f[6]
  rtPrec <- 0.01/60 # 10 ms
  #fileName <- "B:/group/lsfgrp/Mass_Spec/Acquired_data_v2/frimlgrp/LCMS_JiFrATeplova13_m/blank.raw"
  #WIP <- FALSE
  #WIP <- TRUE
  #
  flNm2 <- ExpMap$`MS raw file name`[match(fileName, ExpMap$`MS raw file`)]
  # Default outputs
  qnt <- 0
  peakXtrFine <- setNames(rep(NA, 3),
                          c("Start", "Apex", "End"))
  #
  w <- which(RTchrm01$`Raw file` == fileName)
  if (length(w)) {
    chr <- RTchrm01[w,]
    chr <- chr[which((chr$`Retention time` > RTlim$Start)&(chr$`Retention time` < RTlim$End)),]
    #
    # Quantitative column to use:
    kol <- "Intensity" # We should ALWAYS use raw intensity to identify the peak!
    # Similarly, we should quantify un-corrected intensities, then potentially subtract at a LATER stage
    #
    if (missing("targRg")) {
      mid <- mean(c(RTlim$Start, RTlim$End))
      chr$RTdist <- abs(chr$`Retention time` - mid)
      mx <- max(chr[[kol]][which((chr$`Retention time` >= RTlim$Start)&
                                   (chr$`Retention time` <= RTlim$End))])
    } else { # Here we bias for an expected retention time
      mid <- targRg["Mean"]
      chr$RTdist <- abs(chr$`Retention time` - mid)
      wdth <- max(c(2*targRg["SD"], 0.01), na.rm = TRUE)
      mx <- max(chr[[kol]][which((chr$`Retention time` >= targRg["lWidth"]-wdth)&
                                   (chr$`Retention time` <= targRg["rWidth"]+wdth))])
    }
    wMx <- which(chr[[kol]] == mx)
    if (length(wMx) > 1) { # Tie breaker, for unexpected non-unicity: these are real measurements
      wMx <- wMx[which(chr$RTdist[wMx] == min(chr$RTdist[wMx]))[1]]
    }
    peakXtrFine["Apex"] <- chr$`Retention time`[wMx]
    if (mx > 0) {
      try({
        # First, find background
        # Below: idea which did not work: it assumes a flat baseline
        #chr2 <- chr
        #chr2 <- chr2[order(chr2$Intensity, decreasing = FALSE),]
        #chr2 <- chr2[1:round(nrow(chr2)/10),]
        chrRg <- c(min(chr$`Retention time`),
                   max(chr$`Retention time`))
        chrWdth <- chrRg[2]-chrRg[1]
        nBgBinz <- 5
        bgBins <- data.frame(Bin = 1:nBgBinz)
        bgBins$Start <- chrRg[1] + (bgBins$Bin-1)*chrWdth/nBgBinz
        bgBins$End <- chrRg[1] + bgBins$Bin*chrWdth/nBgBinz
        bgBins$"lowest 10th percentile" <- lapply(1:nBgBinz, function(y) {
          w <- which((chr$`Retention time` > bgBins$Start[y])&
                       (chr$`Retention time` <= bgBins$End[y]))
          y <- setNames(chr$Intensity[w], chr$`Retention time`[w])
          propFlt(y, upper = FALSE, average = FALSE)
        })
        myBg <- plyr::rbind.fill(lapply(bgBins$"lowest 10th percentile", function(y) {
          data.frame(Int = y,
                     RT = as.character(names(y)))
        }))
        myBg <- rbind(data.frame(RT = chr$`Retention time`[1], Int = chr$Intensity[1]),
                      myBg,
                      data.frame(RT = chr$`Retention time`[nrow(chr)], Int = chr$Intensity[nrow(chr)]))
        # Two background identification methods are implemented
        if (bgMode == "loess") {
          bgMod <- loess(Int~RT, data = myBg)
          chr$bgPred <- predict(bgMod, chr$`Retention time`)
        }
        if (bgMode == "spline") {
          bgMod <- smooth.spline(myBg$RT, myBg$Int, nknots = 4) # Fit a spline
          chr$bgPred <- predict(bgMod, chr$`Retention time`)$y
        }
        if (WIP) {
          plot <- ggplot(chr, aes(x = `Retention time`, y = Intensity)) +
            geom_line(linewidth = 0.3, color = "red", aes(y = bgPred)) +
            geom_line() +
            theme_bw() + ggtitle(flNm2)
          #poplot(plot)
        }
        # Identify peak limits
        chr$"Retention time bin" <- round(chr$"Retention time"*60/binSz)
        chrSmooth <- data.table(chr[, c("Retention time bin", "Intensity")])
        chrSmooth <- chrSmooth[,
                               list(Intensity = mean(Intensity)),
                               by = list(`Retention time` = `Retention time bin`)]
        chrSmooth <- as.data.frame(chrSmooth)
        # - Find rough extent of smoothed peak using bins
        #   Assumption: the peak will be the single local maximum in the given range (we are only doing this in narrow ranges)
        #    - Either we have a pure standard and this should be correct,
        #    - or we just transfer the learned RT/width from standards (if targRg isn't missing),
        #    - or we just do the best we can without a standard.
        tMxBin <- round(peakXtrFine["Apex"]*60/binSz)
        wMxBin <- which(chrSmooth$`Retention time` == tMxBin)
        #
        if (WIP) {
          plot2 <- ggplot(chrSmooth, aes(x = `Retention time`)) +
            geom_line(aes(y = Intensity)) +
            theme_bw() + ggtitle(flNm2)
          #poplot(plot2)
        }
        #
        n <- nrow(chrSmooth)
        chrSmooth$Larger_than_prev <- chrSmooth[[kol]] - c(0, chrSmooth[[kol]][1:(n-1)])*1.01
        chrSmooth$Larger_than_next <- chrSmooth[[kol]] - c(chrSmooth[[kol]][2:n], 0)*1.01
        lftExt <- c(wMxBin:1)[which(chrSmooth$Larger_than_prev[wMxBin:1] > 0)]
        rgtExt <- c(wMxBin:n)[which(chrSmooth$Larger_than_next[wMxBin:n] > 0)]
        lftExtTst <- lftExt-c(lftExt[2:length(lftExt)], 0)
        rgtExtTst <- rgtExt-c(rgtExt[2:length(rgtExt)], 0)
        lftExtTst <- contigEq(lftExtTst, 1, 1)
        rgtExtTst <- contigEq(rgtExtTst, -1, 1)
        lftExtTst <- max(which(lftExtTst))
        rgtExtTst <- max(which(rgtExtTst))
        peakXtrI <- c(wMxBin-lftExtTst,
                      wMxBin+rgtExtTst)
        peakXtr <- chrSmooth$`Retention time`[peakXtrI]
        if (length(peakXtr) == 2) {
          if (WIP) {
            plot2 <- plot2 + geom_vline(xintercept = peakXtr, linetype = "dashed")
            poplot(plot2)
          }
          # Then take the local minimum in that bin as refined peak extremity
          tmp <- sapply(peakXtr, function(y) {
            w <- which(chr$"Retention time bin" == y)
            if (y == peakXtr[1]) { w <- rev(w) }
            int <- chr$Intensity[w]
            rt <- chr$`Retention time`[w]
            return(mean(rt[which(int == min(int))]))
          })
          peakXtrFine["Start"] <- tmp[1]
          peakXtrFine["End"] <- tmp[2]
          if (WIP) {
            plot <- plot + geom_vline(xintercept = peakXtrFine, linetype = "dotted")
            #poplot(plot)
          }
          # - Now we want to refine the peaks using a spline model
          wdth <- peakXtrFine["End"]-peakXtrFine["Start"]
          mlt <- 1
          doGoOn <- TRUE
          kount <- 0
          while ((doGoOn)&&(kount < 5)) {
            kount <- kount + 1
            tstFindPeak <- try({
              chr2 <- chr[which((chr$`Retention time` >= peakXtrFine["Start"]-wdth*mlt)&
                                  (chr$`Retention time` <= peakXtrFine["End"]+wdth*mlt)),
                          c("Retention time", "Intensity")]
              colnames(chr2) <- c("RT", "Int")
              pkMod <- smooth.spline(chr2$RT, chr2$Int) # Fit a spline
              pkModPred <- data.frame("Retention time" = chr2$RT,
                                      Intensity = predict(pkMod, chr2$RT)$y,
                                      check.names = FALSE)
              # Try to use derivatives (NB: after smoothing is better)
              pkModPred$logx <- log(pkModPred$"Retention time")
              pkModPred$dlogx <- c(0, diff(pkModPred$logx))
              pkModPred$dy <- c(0, diff(pkModPred$Intensity))
              pkModPred$dy_dlogx <- pkModPred$dy / pkModPred$dlogx
              pkModPred <- pkModPred[which(!is.na(pkModPred$dy_dlogx)),]
              w <- which((pkModPred$`Retention time` < min(pkModPred$`Retention time`) + wdth*0.1*mlt)|
                           (pkModPred$`Retention time` > max(pkModPred$`Retention time`) - wdth*0.1*mlt))
              pkDerivLin <- nls(dy_dlogx ~ a*`Retention time` + b,
                                data = pkModPred[w,],
                                start = list(a = 0,
                                             b = mean(c(pkModPred$dy_dlogx[w][1],
                                                        rev(pkModPred$dy_dlogx[w])[1])))) # Fit a line
              if (WIP) {
                plot3 <- ggplot(pkModPred) +
                  geom_line(aes(x = `Retention time`, y = dy_dlogx)) +
                  geom_abline(slope = pkDerivLin$m$getPars()["a"], intercept = pkDerivLin$m$getPars()["b"]) +
                  theme_bw() + ggtitle("Derivative", subtitle = flNm2)
                #poplot(plot3)
              }
              xtrm2 <- c(min(pkModPred$`Retention time`), max(pkModPred$`Retention time`))
              wdth2 <- xtrm2[2]-xtrm2[1]
              n <- round(wdth2/rtPrec)
              tstDF <- data.frame(step = 1:n,
                                  RT = xtrm2[1] + ((1:n)-1)*rtPrec)
              tstDF$line <- pkDerivLin$m$getPars()["a"]*tstDF$RT + pkDerivLin$m$getPars()["b"]
              tstDF$spline <- vapply(tstDF$RT, function(y) {
                z <- which(pkModPred$`Retention time` == y)
                if (length(z)) {
                  v <- pkModPred$dy_dlogx[z]
                } else {
                  z1 <- max(which(pkModPred$`Retention time` < y))
                  z2 <- min(which(pkModPred$`Retention time` > y))
                  t1 <- pkModPred$`Retention time`[z1]
                  t2 <- pkModPred$`Retention time`[z2]
                  v1 <- pkModPred$dy_dlogx[z1]
                  v2 <- pkModPred$dy_dlogx[z2]
                  v <- v1 + (v2-v1)*(y-t1)/(t2-t1)
                }
                return(v)
              }, 1)
              tstDF$Test <- tstDF$spline-tstDF$line
              tstDF$cross <- NA
              tstDF$cross[2:n] <- (tstDF$Test[1:(n-1)] >= 0)&(tstDF$Test[2:n] <= 0)
              wMxx <- which(tstDF$cross[2:n])
              if (length(wMxx) > 1) {
                tstDF$dist <- abs(tstDF$RT - peakXtrFine["Apex"])
                wMxx <- round(mean(wMxx[which(tstDF$dist[wMxx] == min(tstDF$dist[wMxx]))]))
              }
              peakXtrFine["Apex"] <- tstDF$RT[wMxx] # Refined peak apex (may not be worth doing)
              xtrm3a <- c(min(tstDF$RT), max(tstDF$RT))
              xtrm3b <- c(min(tstDF$Test), max(tstDF$Test))
              #tstDF$Test2 <- tstDF$Test - (xtrm3b[2] - xtrm3b[1])*(tstDF$RT - xtrm3a[1])/(xtrm3a[2] - xtrm3a[1])
              #tstDF$Test2 <- tstDF$Test2 - tstDF$Test2[wMxx]
              if (WIP) {
                plot3 <- plot3 +
                  geom_line(data = tstDF, aes(x = RT, y = Test), color = "cyan")
                #geom_line(data = tstDF, aes(x = RT, y = Test2), color = "blue")
                #poplot(plot3)
              }
              # Detect peak boundaries
              tstDF$high <- NA
              #tstDF$Test3[2:n] <- tstDF$Test2[2:n] > tstDF$Test2[1:(n-1)]
              tstDF$Test3[2:n] <- tstDF$Test[2:n] > tstDF$Test[1:(n-1)]
              wOK <- which((((tstDF$Test > 0)&(tstDF$RT < peakXtrFine["Apex"]))|#Only consider parts going up...
                              ((tstDF$Test < 0)&(tstDF$RT > peakXtrFine["Apex"])))& # ... and down on either sides, and
                             (abs(tstDF$Test/max(xtrm3b)) < 0.01)) # smallest local slope allowed for cutoff (gotta be within range of the plain)
              flt1 <- c(min(which(!c(1:n) %in% wOK)),
                        max(which(!c(1:n) %in% wOK)))
              wOK <- wOK[which(!wOK %in% 1:flt1[1])]
              wOK <- wOK[which(!wOK %in% flt1[2]:n)]
              minOK <- min(wOK)
              maxOK <- max(wOK)
              tst <- vapply(c(-1, 1), function(y) {
                suppressWarnings({
                  #tstDF$RT[w1]
                  #tstDF$RT[w2]
                  if (y == -1) { # move along
                    w1 <- max(c(minOK:wMxx)[which(tstDF$Test3[minOK:wMxx])], na.rm = TRUE)
                    w2 <- max(c(minOK:w1)[which(!tstDF$Test3[minOK:w1])], na.rm = TRUE)
                    goOn <- tstDF$Test[w2] > 0
                    while ((goOn)&&(!w2 %in% wOK)) {
                      r <- w1 <- max(c(minOK:w2)[which(tstDF$Test3[minOK:w2])], na.rm = TRUE)
                      if (is.finite(w1)) {
                        r <- max(c(minOK:w1)[which(!tstDF$Test3[minOK:w1])], na.rm = TRUE)
                      }
                      if (is.finite(r)) {
                        w2 <- r
                        goOn <- tstDF$Test[w2] > 0
                      } else { goOn <- FALSE }
                    }
                  }
                  if (y == 1) {
                    w1 <- suppressWarnings(min(c(wMxx:maxOK)[which(tstDF$Test3[wMxx:maxOK])], na.rm = TRUE))
                    w2 <- min(c(w1:maxOK)[which(!tstDF$Test3[w1:maxOK])], na.rm = TRUE)
                    goOn <- tstDF$Test[w2] < 0
                    while ((goOn)&&(!w2 %in% wOK)) {
                      r <- w1 <- min(c(w2:maxOK)[which(tstDF$Test3[w2:maxOK])], na.rm = TRUE)
                      if (is.finite(w1)) {
                        r <- min(c(w1:maxOK)[which(!tstDF$Test3[w1:maxOK])], na.rm = TRUE)
                      }
                      
                      if (is.finite(r)) {
                        w2 <- r
                        goOn <- tstDF$Test[w2] < 0
                      } else { goOn <- FALSE }
                    }
                  }
                })
                return(w2)
              }, 1)
              #tstDF$RT[tst]
              peakXtrFine["Start"] <- tstDF$RT[tst[1]]
              peakXtrFine["End"] <- tstDF$RT[tst[2]]
              if (WIP) {
                plot3 <- plot3 +
                  geom_vline(xintercept = peakXtrFine, color = "magenta", linetype = "dashed")
                poplot(plot3)
                plot <- plot +
                  geom_vline(xintercept = peakXtrFine, color = "magenta", linetype = "dashed")
                poplot(plot)
              }
            }, silent = TRUE)
            doGoOn <- "try-error" %in% class(tstFindPeak)
            if (doGoOn) {
              mlt <- mlt*1.2
            }
          }
          peakXtrFine <- vapply(peakXtrFine, function(x) {
            tmp <- abs(chr$`Retention time` - x)
            chr$`Retention time`[which(tmp == min(tmp))[1]]
          }, 1)
          #
          m <- match(peakXtrFine[c("Start", "End")], chr$`Retention time`)
          pkRg <- sort(m[1]:(m[2]-1)) # Shouldn't ever need sorting, but jic
          qnt <- sum(vapply(pkRg, function(i) {
            mean(chr$Intensity[i+0:1]-chr$bgPred[i+0:1])*(chr$`Retention time`[i+1]-chr$`Retention time`[i])
          }, 1))
        }
      }, silent = TRUE)
    }
  }
  if ((length(qnt) != 1)||(is.null(qnt))) {
    warning(paste0(x, ": something unexpected happened, investigate!"))
    qnt <- NA
  } # Should never happen
  return(list(Quant = qnt,
              PeakRT = peakXtrFine))
}
for (nm in names(allChroms)) { #nm <- names(allChroms)[1] #nm <- names(allChroms)[2] #nm <- names(allChroms)[4]
  chrm <- allChroms[[nm]]
  chrm <- chrm[which((chrm$`Retention time` >= fullRTRange[1])&(chrm$`Retention time` <= fullRTRange[2])),]
  Ykol <- c("Intensity", "Pressure")[(nm == "Pressure")+1]
  if ("data.frame" %in% class(chrm)) {
    if (nm == "Pressure") { # We only plot the first minute to check for bubbles
      chrm <- chrm[which(chrm$`Retention time` <= 1),]
    }
    if (!nm %in% c("TIC", "BPC", "Pressure")) {
      tmp <- unlist(strsplit(nm, " for "))
      mass <- as.numeric(tmp[1])
      ms2 <- tmp[2]
      ms2Trg <- as.numeric(gsub(".* ", "", ms2))
    }
    chrm <- chrm[which(chrm$`Raw file name` %in% ExpMap$`MS raw file name`),]
    chrm$`Raw file` <- factor(chrm$`Raw file`, levels = c(slctFls, blnkFls))
    chrm$`Raw file name` <- factor(chrm$`Raw file name`, levels = c(slctFls0, blnkFls0))
    chrm$Sort <- match(chrm$`Raw file`, ExpMap2$`MS raw file`)
    chrm <- chrm[order(chrm$Sort),]
    chrm[, Factors] <- ExpMap2[chrm$Sort, Factors]
    chrm$Sample <- factor(chrm$`Raw file`, levels = slctFls)
    chrm$Role[which(chrm$`Raw file` %in% blnkFls)] <- "Blank"
    chrm$Linetype <- "solid"
    chrm$Linetype[which(chrm$Role == "Buffer_control")] <- "dotted"
    chrm$Linetype[which(chrm$Role == "Tag_control")] <- "dotdash"
    chrm$Linetype[which(chrm$Role == "Standard")] <- "dashed"
    #offset <- max(chrm$Int)/3
    #nc <- nchar(floor(offset))
    #offset <- ceiling(offset/10^(nc-1))*10^(nc-1)
    #if (!is.finite(offset)) { offset <- 0 }
    #chrm$Int <- chrm[[Ykol]] + offset*(match(chrm$`Raw file`, slctFls)-1)
    #unique(chrm$`Raw file`)
    #
    wRef <- which(chrm$Role == "Buffer_control")
    quantThisOne <- (Ykol != "Pressure")&(Quant)&(length(wRef) > 0)#&(!nm %in% c("Pressure", "TIC", "BPC"))
    #
    availGrps <- unique(chrm$"Analysis_group")
    availGrps <- availGrps[which(!is.na(availGrps))]
    for (grp in availGrps) { #grp <- availGrps[1]
      availRTRgs <- 1
      if (!nm %in% names(allChroms)[1:3]) {
        availRTRgs <- 1:nrow(rtRanges)
      }
      availRTRgs <- availRTRgs[which(sapply(rtRanges$"Analysis_group"[availRTRgs], function(x) {
        grp %in% x
      }))] # Available RT ranges to plot for this analysis group
      #
      if (length(availRTRgs)) {
        # Look at the different RT ranges
        for (I in availRTRgs) { #I <- availRTRgs[1] #I <- availRTRgs[3] #I <- availRTRgs[4]
          rg0 <- ctrlRanges[[grp]] # Controls file range
          RTlim <- rtRanges[I, c("Start", "End")] # Retention range range
          RTchrm <- chrm[which((chrm$`Retention time` >= RTlim$Start)&(chrm$`Retention time` <= RTlim$End)),]
          Targ <- rtRanges$Nucleotide[I]
          w0 <- which(RTchrm$`Raw file` %in% ExpMap2$`MS raw file`[rg0])
          for (prot in Protein) { #prot <- Protein[1] #prot <- Protein[2] #prot <- Protein[3] # Let's plot each protein separately
            pat <- paste0("^", grp, " ", prot, " ")
            flRgs <- grep(paste0(pat, "[0-9]+$"), names(flRanges), value = TRUE) # Here these are file ranges, not RT ranges
            if (length(flRgs)) {
              flRgs <- gsub(pat, "", flRgs)
              #
              # Create formula for facet_grid:
              if (is.na(prot)) {
                em01 <- ExpMap2[which((is.na(ExpMap2$Protein))|(ExpMap2$Role %in% c("Standard", "Buffer_control", "Tag_control"))),]
              } else {
                em01 <- ExpMap2[which((ExpMap2$Protein == prot)|(ExpMap2$Role %in% c("Standard", "Buffer_control", "Tag_control"))),]
              }
              em1 <- em01[which(!em01$Role %in% c("Standard", "Buffer_control")),]
              Factors3 <- Factors[which(!Factors %in% c("Replicate", "Samples_group"))] # Not using Samples_group because we replace it with Samples_group2!!!
              tst <- sapply(Factors3, function(fct) {
                (fct %in% colnames(em1))&&(length(unique(em1[[fct]])) > 1)
              })
              Factors3 <- Factors3[which(tst)]
              if (length(Factors3) > 1) {
                # Check for synonymous factors
                tmp <- em01[, Factors3]
                for (fct in Factors3) {
                  tmp[[fct]] <- as.numeric(as.factor(tmp[[fct]]))
                }
                comb <- gtools::combinations(length(Factors3), 2, Factors3)
                comb <- as.data.frame(comb)
                tst <- apply(comb, 1, function(x) {
                  identical(tmp[[x[1]]], tmp[[x[2]]])
                })
                w <- which(tst)
                if (length(w)) {
                  fct2Remov <- unique(comb[w, 2])
                  Factors3 <- Factors3[which(!Factors3 %in% fct2Remov)]
                }
              }
              if (!length(Factors3)) {
                Factors3 <- "Protein"
              }
              Factors3 <- unique(c(Factors3, "Role")) # Role must ALWAYS be present!
              #if (length(Factors3)) {
              Form0 <- as.formula(paste0("`MS raw file name` ~ ", paste(c("Replicate", Factors3), collapse = " + "))) # This formula is used to check that the factors used to create the fact grid formula (Form, below) allow for unique samples per facet
              tst <- aggregate(Form0, data = em01, function(x) { length(unique(x)) })
              tst <- max(tst$`MS raw file name`)
              if (tst > 1) {
                Factors4 <- Factors[which(!Factors %in% c("Replicate", Factors3))] # Should always be at least length 1
                l <- length(Factors4)
                stopifnot(l > 1) # This means something unexpected happened: several raw files have been given the same factors and will be indistinguishable,
                # which should never happen
                # If this becomes a problem we will have to build in a check into the experiment map app
                goOn <- TRUE
                i <- 1
                while ((goOn)&&(i <= l)) {
                  comb <- gtools::combinations(length(Factors4), i, Factors4)
                  tst <- lapply(1:nrow(comb), function(x) { c(Factors3, comb[x, ]) })
                  Form0 <- sapply(tst, function(x) { as.formula(paste0("`MS raw file name` ~ ", paste(c("Replicate", x), collapse = " + "))) })
                  tst <- sapply(Form0, function(x) { max(aggregate(x, data = em01, function(y) { length(unique(y)) })$`MS raw file name`) })
                  if (1 %in% tst) {
                    goOn <- FALSE
                    w <- which(tst == 1)[1]
                    Factors3 <- c(Factors3, comb[w,])
                  } else {
                    i <- i+1
                  }
                }
              }
              Form <- as.formula(paste0(paste(Factors3, collapse = " + "), " ~ Replicate")) # Formula for how to layout our facets on the ggplot
              em01$Color_class <- do.call(paste, c(em01[, Factors3, drop = FALSE], sep = " "))
              # } else {
              #   Form <- as.formula(paste0(". ~ Replicate")) # Formula for how to layout our facets on the ggplot
              #   em01$Color_class <- NA
              # }
              em01$Color_class[which((em01$Color_class == "NA")&(em01$Role == "Standard"))] <- "Standard"
              em01$Color_class[which((em01$Color_class == "NA")&(em01$Role == "Buffer_control"))] <- "Buffer_control"
              #
              for (rgI in flRgs) { #rgI <- flRgs[1]
                rg1 <- flRanges[[paste0(grp, " ", prot, " ", rgI)]] # Samples file range
                w1 <- which(as.character(RTchrm$`Raw file`) %in% ExpMap2$`MS raw file`[rg1]) # not em01!
                RTchrm01 <- RTchrm[c(w0, w1),]
                if (nm %in% c("TIC", "BPC", "Pressure")) {
                  ttl <- paste0(dtstNm, " - ", nm)
                } else {
                  ttl <- paste0(dtstNm, " - XIC, transition ", ms2Trg, " to ", mass)
                }
                if (length(Analysis_group) > 1) { ttl <- paste0(ttl, " ", grp) }
                if (length(Protein) > 1) { ttl <- paste0(ttl, " ", prot) }
                if (I > 1) {
                  ttl <- paste0(ttl, ", RT ", rtRanges$Start[I], "-", rtRanges$End[I], " min")
                }
                if (length(flRgs) > 1) {
                  ttl <- paste0(ttl, "_", rgI)
                }
                quantThisOne2 <- quantThisOne&(I > 1)&(!is.na(Targ))
                if (quantThisOne2) {
                  f <- as.character(unique(RTchrm01$`Raw file`))
                  m <- match(f, em01$`MS raw file`)
                  stdF <- f[which((em01$Role[m] == "Standard")
                                  &(em01$Nucleotide[m] == Targ))]
                  quantThisOne2 <- length(stdF)
                }
                if (quantThisOne2) {
                  # Case 1 - We have standards
                  f2 <- f[which(!f %in% stdF)]
                  # Find peak of standards
                  quantLst <- setNames(lapply(stdF, peakXtract), stdF)
                  RTchar <- as.data.frame(t(sapply(stdF, function(x) { quantLst[[x]]$PeakRT })))
                  RTchar <- c(Mean = mean(RTchar$Apex, na.rm = TRUE),
                              SD = sd(RTchar$Apex, na.rm = TRUE),
                              lWidth = min(RTchar$Start, na.rm = TRUE),
                              rWidth = max(RTchar$End, na.rm = TRUE))
                  # Find peak of other samples using information from standards
                  quantLst <- append(quantLst, lapply(f2, peakXtract, targRg = RTchar))
                  names(quantLst) <- c(stdF, f2)
                  #
                  quant <- sapply(quantLst, function(x) { x$Quant })
                  quant <- data.frame("Raw file" = names(quantLst),
                                      "Peak intensity" = quant,
                                      check.names = FALSE)
                  quant[, c("Peak start", "Peak apex", "Peak end")] <- as.data.frame(t(sapply(quantLst, function(x) { x$PeakRT })))
                  quant <- quant[which((!is.na(quant$`Peak intensity`))&(quant$`Peak intensity` > 0)),]
                  quant[, c(Factors, "Raw file name")] <- ExpMap2[match(quant$"Raw file", ExpMap2$`MS raw file`), c(Factors, "MS raw file name")]
                  quant$"Peak intensity (rounded)" <- round(quant$"Peak intensity")
                  quant$Color_class <- em01$Color_class[match(quant$"Raw file", em01$`MS raw file`)]
                  fwrite(quant, paste0(wd, "/", ttl, "_quant.csv"), quote = FALSE, sep = ",", col.names = TRUE, na = "NA")
                }
                #
                for (fct in Factors) {
                  RTchrm01[[fct]] <- factor(RTchrm01[[fct]], unique(em01[[fct]]))
                }
                RTchrm01$Color_class <- em01$Color_class[match(RTchrm01$Samples_group2, em01$Samples_group2)]
                #RTchrm1 <- RTchrm01[which(!RTchrm01$Role %in% c("Standard", "Buffer_control")),]
                #RTchrm0 <- RTchrm01[which(RTchrm01$Role %in% c("Standard", "Buffer_control")),]
                #smpls <- unique(RTchrm1$Sample)
                # We will plot raw intensities, but note that quantitation includes background correction to buffer control for samples
                yMax <- max(RTchrm01[[Ykol]])
                xLim <- c(min(RTchrm01$`Retention time`),
                          max(RTchrm01$`Retention time`))
                xRange <- xLim[2]-xLim[1]
                plot <- ggplot(RTchrm01,
                               aes(x = `Retention time`, y = .data[[Ykol]], #linetype = Role,
                                   group = `Raw file name`, color = Color_class)) +
                  geom_line() +
                  ggtitle(ttl, subtitle = c("", Targ)[(!is.na(Targ))+1]) +
                  coord_fixed(xRange/(yMax*3)) + facet_grid(Form) +
                  #scale_linetype_manual(values = lnTypes, guide = "legend") +
                  theme_bw() +
                  theme(strip.text.y = element_text(angle = 0),
                        legend.position = "bottom") +
                  ylim(0, yMax)
                #
                if (quantThisOne2) {
                  plot <- plot +
                    geom_vline(data = quant, aes(xintercept = `Peak start`), linetype = "dotted") +
                    geom_vline(data = quant, aes(xintercept = `Peak apex`), linetype = "dotted") +
                    geom_vline(data = quant, aes(xintercept = `Peak end`), linetype = "dotted") +
                    geom_text(data = quant, aes(label = `Peak intensity (rounded)`,
                                                x = min(RTchrm01$`Retention time`) + xRange*0.05),
                              y = yMax*0.8, hjust = 0, size = 2.5, show.legend = FALSE)
                  poplot(plot, 12, 22)
                }
                wdth <- (max(Replicate)+2)*2
                ggsave(paste0(wd, "/", ttl, ".jpeg"), plot, dpi = 300, width = wdth, units = "in")
                ggsave(paste0(wd, "/", ttl, ".pdf"), plot, dpi = 300, width = wdth, units = "in")
                system(paste0("open \"", wd, "/", ttl, ".jpeg\""))
              }
            }
          }
        }
      }
    }
  } else {
    if (nm %in% c("TIC", "BPC")) {
      warning(paste0("Failure to extract ", nm))
    } else {
      warning(paste0("Failure to extract XICs @ mass = ", mass, " Da, tol. = ", tol, ", filter = '", ms2, "'"))
    }
  }
}
#nms <- names(allChroms)[which(names(allChroms) != "Pressure")] # No need to write pressure, we have it already
nms <- names(allChroms) # Actually let's write it, it's spread over different files
sapply(nms, function(nm) { #nm <- "BPC"
  nm2 <- paste0(nm, ".csv")
  if (!nm %in% c("TIC", "BPC", "Pressure")) { nm2 <- paste0("XIC ", nm2) }
  x <- allChroms[[nm]]
  x$`Raw file name` <- NULL
  x <- x[, c("Raw file", colnames(x)[which(colnames(x) != "Raw file")])]
  data.table::fwrite(x, paste0(wd, "/", nm2), row.names = FALSE, na = "NA")
})
openwd()

# To do:
# - Parallelize plots printing bit
# - Shiny-based app to display single graph for file of interest
# - Propagate lessons learned in this script to "Convert Raw files..." script
