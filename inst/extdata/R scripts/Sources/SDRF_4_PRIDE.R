# This is an embryo script to write SDRF files for PRIDE submissions.
# For now it is incomplete, and only works for label-free datasets.
# Assumes Trypsin
# Assumes...
# Just...
# Check...
# The output!!!
#
# Check the specs of the SDRF format here:
# https://github.com/bigbio/proteomics-sample-metadata/blob/master/sdrf-proteomics/README.adoc

require(shiny)
require(shinyWidgets)
require(shinyjs)
require(shinycssloaders)
require(DT)
require(proteoCraft)
if (!exists("homePath")) {
  homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
}
# For validations, make sure that python is installed and during Configure() make sure the parse_sdrf.exe script is added!!!
locDirsFl <- paste0(homePath, "/Default_locations.xlsx")
if (file.exists(locDirsFl)) {
  locDirs <- openxlsx2::read_xlsx(locDirsFl)
  if ("Python" %in% locDirs$Folder) {
    pyTest <- !as.logical(try({
      #pyPath0 <- locDirs$Path[match("Python", locDirs$Folder)]
      #cat(pyPath <- gsub("/", "\\\\", pyPath0))
      # Update pip
      cmd <- "pip install --upgrade pip"
      system(cmd, show.output.on.console = FALSE)
      # Install sdrf-pipelines
      cmd <- "pip install sdrf-pipelines"
      system(cmd, show.output.on.console = FALSE)
    }, silent = TRUE))
  }
}
#load_Bckp()

setwd(wd)
sdrfPth <- paste0(wd, "/SDRF.tsv")
myOntDF <- data.frame(Name = c("CellType",
                               "Tissue",
                               "Disease",
                               "DevStages",
                               "MSInstr",
                               "AcqMeth",
                               "QuantMeth",
                               "LabelMeth",
                               "PTMs"),
                      Column = c("characteristics[cell type]",
                                 "characteristics[organism part]",
                                 "characteristics[disease]",
                                 "characteristics[developmental stage]",
                                 "comment[instrument]",
                                 "comment[proteomics data acquisition method]",
                                 NA,
                                 "comment[label]",
                                 "characteristics[enrichment process]"))
myOntDF$myOntology <- paste0("my", myOntDF$Name)

reload_SDRF <- file.exists(sdrfPth)
if (reload_SDRF) {
  auld_SDRF <- read.delim(sdrfPth, sep = "\t", check.names = FALSE)
  filtOnt <- function(ontNm, colNm) {
    if ((!is.na(colNm))&&(colNm %in% colnames(auld_SDRF))) {
      tmp <- unique(auld_SDRF[[colNm]])
      tmp <- tmp[which(!is.na(tmp))]
      tmp <- tmp[which(!tolower(tmp) %in% c("unknown", "not available"))]
      tmp <- tmp[which(nchar(tmp) > 0)]
      assign(ontNm, tmp, envir = .GlobalEnv)
    }
  }
  for (i in 1:nrow(myOntDF)) {
    filtOnt(myOntDF$myOntology[i], myOntDF$Column[i])
  }
}
#
w <- which(c("Organism_Full", "Organism") %in% colnames(db))
if (!exists("tstOrg")) {
  tstOrg <- (length(w) > 0)
}
if (tstOrg) {
  dbOrgKol <- c("Organism_Full", "Organism")[w[1]]
  tst <- gsub(" *(\\(|\\[).*", "", db[[dbOrgKol]])
  tst <- aggregate(tst, list(tst), length)
  tst <- tst[order(tst$x, decreasing = TRUE),]
  mainOrg %<o% tst$Group.1[1]
} else {
  dbOrgKol <- c()
  mainOrg <- "[UNKNOWN_ORGANISM]"
}
dbOrgKol %<o% dbOrgKol
tstOrg %<o% tstOrg
mainOrg %<o% mainOrg
#
ontoFls <- data.frame(File = paste0(homePath, "/", c("Tissues",
                                                     "Modifications",
                                                     "Diseases",
                                                     "Cell_types",
                                                     "PRIDE",
                                                     "MS_acq_meth",
                                                     "MS_quant_meth",
                                                     "MS_label_meth",
                                                     "MS_models",
                                                     paste0("DevStages_", c("Human",
                                                                            "Mouse",
                                                                            "Zebrafish",
                                                                            "Fly",
                                                                            "Worms"))),
                                    ".csv"))
ontoFls$Name <- gsub("\\.csv$", "", basename(ontoFls$File))
ontoFls$Exists <- file.exists(ontoFls$File)
#sum(!ontoFls$Exists)
#View(ontoFls)
isFnd <- setNames(ontoFls$Exists, ontoFls$Name)
if (isFnd["Cell_types"]) {
  tmp <- data.table::fread(ontoFls$File[match("Cell_types", ontoFls$Name)])
  availCellTypes <- sort(unique(tmp$label))
  cellTypeTxt <- " the cell type analyzed"
  if ((!exists("myCellType"))||(!is.character(myCellType))||(sum(!myCellType %in% availCellTypes))) {
    myCellType <- c()
    myCellType2 <- "not available"
  } else {
    myCellType2 <- myCellType
  }
  availCellTypes <- c(myCellType, availCellTypes[which(!availCellTypes %in% myCellType)])
}
if (isFnd["Tissues"]) {
  tmp <- data.table::fread(ontoFls$File[match("Tissues", ontoFls$Name)])
  availTissues <- sort(unique(tmp$label))
  tissueTxt <- " the tissue or organism part analyzed"
  if ((!exists("myTissue"))||(!is.character(myTissue))||(sum(!myTissue %in% availTissues))) {
    myTissue <- c()
    myTissue2 <- "not available"
  } else {
    myTissue2 <- myTissue
  }
  availTissues <- c(myTissue, availTissues[which(!availTissues %in% myTissue)])
}
if (isFnd["Modifications"]) {
  tmp <- data.table::fread(ontoFls$File[match("Modifications", ontoFls$Name)])
  availMods <- sort(unique(tmp$label))
  modsTxt <- " any Post-Translational Modifications which were included in the analysis"
  if ((!exists("myPTMs"))||(!is.character(myPTMs))||(sum(!myPTMs %in% availMods))) {
    myPTMs <- c()
    myPTMs2 <- "not available"
  } else {
    myPTMs2 <- myPTMs
  }
  availMods <- c(myPTMs, availMods[which(!availMods %in% myPTMs)])
}
if (isFnd["Diseases"]) {
  tmp <- data.table::fread(ontoFls$File[match("Diseases", ontoFls$Name)])
  availDiseases <- sort(unique(tmp$label))
  diseasesTxt <- " any diseases relevant to the dataset"
  if ((!exists("myDisease"))||(!is.character(myDisease))||(sum(!myDisease %in% availDiseases))) {
    myDisease <- c()
    myDisease2 <- "not available"
  } else {
    myDisease2 <- myDisease
  }
  availDiseases <- c(myDisease, availDiseases[which(!availDiseases %in% myDisease)])
}
if (isFnd["MS_models"]) {
  availInstrDF <- data.table::fread(ontoFls$File[match("MS_models", ontoFls$Name)])
  availInstr <- availInstrDF$Instrument
  Vendors2Instr <- aggregate(availInstrDF$Instrument, list(availInstrDF$Vendor), function(x) { sort(unique(x)) })
  colnames(Vendors2Instr) <- c("Vendor", "Instruments")
  availMSVend <- c("All vendors", sort(Vendors2Instr$Vendor))
  tmp2 <- data.frame(Vendor = "All vendors")
  tmp2$Instruments <- list(sort(unique(unlist(Vendors2Instr$Instruments))))
  Vendors2Instr <- rbind(Vendors2Instr, tmp2)
  MSInstrTxt <- " the Mass Spectrometer model(s) on which the data was acquired"
  if ((!exists("myVendor"))||(!is.character(myVendor))||(sum(!myVendor %in% Vendors2Instr$Vendor))) { myVendor <- "All vendors" }
  myMSInstr <- c()
  if ((!exists("myMSInstr"))||(!is.character(myMSInstr))||(sum(!myMSInstr %in% availInstr))) {
    if ((exists("LCMS_instr"))&&("list" %in% class(LCMS_instr))&&("MS" %in% names(LCMS_instr))) {
      myMSInstr <- availInstr[which(availInstr %in% LCMS_instr$MS)]
    }
  }
  myMSInstr2 <- myMSInstr
  if (!length(myMSInstr2)) { myMSInstr2 <- "not available" }
  availInstr <- c(myMSInstr, availInstr[which(!availInstr %in% myMSInstr)])
}
if (isFnd["MS_acq_meth"]) {
  tmp <- data.table::fread(ontoFls$File[match("MS_acq_meth", ontoFls$Name)])
  availAcqMeth <- sort(unique(tmp$Name))
  acqMethTxt <- " relevant MS acquisition method(s)"
  if ((!exists("myAcqMeth"))||(!is.character(myAcqMeth))||(!sum(myAcqMeth %in% availAcqMeth))) {
    myAcqMeth <- c()
    myAcqMeth2 <- "not available"
  } else {
    myAcqMeth2 <- myAcqMeth
  }
  availAcqMeth <- c(myAcqMeth, availAcqMeth[which(!availAcqMeth %in% myAcqMeth)])
}
if (isFnd["MS_quant_meth"]) {
  tmp <- data.table::fread(ontoFls$File[match("MS_quant_meth", ontoFls$Name)])
  availQuantMeth <- sort(unique(tmp$Name))
  quantMethTxt <- " relevant MS quantification method(s)"
  if ((!exists("myQuantMeth"))||(!is.character(myQuantMeth))||(!sum(myQuantMeth %in% availQuantMeth))) {
    myQuantMeth <- c()
    myQuantMeth2 <- "not available"
  } else {
    myQuantMeth2 <- myQuantMeth
  }
  availQuantMeth <- c(myQuantMeth, availQuantMeth[which(!availQuantMeth %in% myQuantMeth)])
}
if (isFnd["MS_label_meth"]) {
  tmp <- data.table::fread(ontoFls$File[match("MS_label_meth", ontoFls$Name)])
  availLabelMeth <- sort(unique(tmp$Name))
  labelMethTxt <- " relevant MS labelling method(s)"
  if ((!exists("myLabelMeth"))||(!is.character(myLabelMeth))||(!sum(myLabelMeth %in% availLabelMeth))) {
    myLabelMeth <- c()
    myLabelMeth2 <- "not available"
  } else {
    myLabelMeth2 <- myLabelMeth
  }
  availLabelMeth <- c(myLabelMeth, availLabelMeth[which(!availLabelMeth %in% myLabelMeth)])
}
myDS1 <- NULL
myDS2 <- ""
appNm <- paste0(dtstNm, " - SDRF editor")
uiA <- shiny::fluidPage(
  shinyjs::useShinyjs(),
  shinyWidgets::setBackgroundColor(
    color = "#E6F7F4",
    gradient = "linear",
    direction = "bottom"
  ),
  shinyjs::extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  shiny::titlePanel(shiny::tag("u", "SDRF editor"),
                    appNm),
  shiny::h2(dtstNm),
  h5("An SDRF file is a controlled vocabulary map of MS runs to biological samples, and can be requested when uploading a dataset to public repositories such as https://www.ebi.ac.uk/pride/"),
  if (!tstOrg) {
    shiny::fluidRow(shiny::column(4,
                                  shiny::tags$hr(style = "border-color: black;"),
                                  shiny::textInput("Organism",
                                                   "Enter parent organism (ontology = https://www.ebi.ac.uk/ols4/ontologies/NCBITaxon)",
                                                   "Homo sapiens"),
                                  shiny::br()))
  },
  shiny::tags$hr(style = "border-color: black;"),
  shiny::br(),
  shiny::fluidRow(
    if (isFnd["Cell_types"]) {
      shiny::column(4,
                    shinyWidgets::pickerInput("Cell_type",
                                              paste0("Select", cellTypeTxt),
                                              availCellTypes,
                                              myCellType,
                                              TRUE,
                                              shinyWidgets::pickerOptions(title = "Search me",
                                                                          `live-search` = TRUE,
                                                                          actionsBox = TRUE,
                                                                          deselectAllText = "Clear search",
                                                                          showTick = TRUE)),
                    shiny::br())
    } else {
      shiny::column(4,
                    shiny::textInput("Cell_type",
                                     paste0("Input", cellTypeTxt, " (ontology = https://www.ebi.ac.uk/ols4/ontologies/CL)"),
                                     myCellType2),
                    shiny::br())
    },
    if (isFnd["Tissues"]) {
      shiny::column(4,
                    shinyWidgets::pickerInput("Tissue_type",
                                              paste0("Select", tissueTxt),
                                              availTissues,
                                              myTissue,
                                              TRUE,
                                              shinyWidgets::pickerOptions(title = "Search me",
                                                                          `live-search` = TRUE,
                                                                          actionsBox = TRUE,
                                                                          deselectAllText = "Clear search",
                                                                          showTick = TRUE)),
                    shiny::br())
    } else {
      shiny::column(4,
                    shiny::textInput("Tissue_type",
                                     paste0("Input", tissueTxt, " (ontology = https://www.ebi.ac.uk/ols4/ontologies/BTO)"),
                                     myTissue2),
                    shiny::br())
    },
    if (isFnd["Diseases"]) {
      shiny::column(4,
                    shinyWidgets::pickerInput("Diseases",
                                              paste0("Select", diseasesTxt),
                                              availDiseases,
                                              myDisease,
                                              TRUE,
                                              shinyWidgets::pickerOptions(title = "Search me",
                                                                          `live-search` = TRUE,
                                                                          actionsBox = TRUE,
                                                                          deselectAllText = "Clear search",
                                                                          showTick = TRUE)),
                    shiny::br())
    } else {
      shiny::column(4,
                    shiny::textInput("Diseases",
                                     paste0("Input", diseasesTxt, " (ontology = https://www.ebi.ac.uk/ols4/ontologies/DOID)"),
                                     myDisease2),
                    shiny::br())
    },
  ),
  shiny::tags$hr(style = "border-color: black;"),
  strong("Developmental stage(s)"),
  shiny::uiOutput("myDevStages"),
  shiny::br(),
  shiny::tags$hr(style = "border-color: black;"),
  shiny::br(),
  shiny::fluidRow(
    shiny::column(2,
                  if (isFnd["MS_models"]) {
                    shinyWidgets::pickerInput("Vendor",
                                              "Select instrument vendor(s)",
                                              availMSVend,
                                              myVendor,
                                              TRUE,
                                              shinyWidgets::pickerOptions(title = "Search me",
                                                                          `live-search` = TRUE,
                                                                          actionsBox = TRUE,
                                                                          deselectAllText = "Clear search",
                                                                          showTick = TRUE))
                  },
                  shiny::uiOutput("myMS_instr"),
                  shiny::br()),
    if (isFnd["MS_acq_meth"]) {
      shiny::column(2,
                    shinyWidgets::pickerInput("MS_acq_meth",
                                              paste0("Select", acqMethTxt),
                                              availAcqMeth,
                                              myAcqMeth,
                                              TRUE,
                                              shinyWidgets::pickerOptions(title = "Search me",
                                                                          `live-search` = TRUE,
                                                                          actionsBox = TRUE,
                                                                          deselectAllText = "Clear search",
                                                                          showTick = TRUE)),
                    shiny::br())
    } else {
      shiny::column(2,
                    shiny::textInput("MS_acq_meth",
                                     paste0("Input", acqMethTxt, " (ontology = https://www.ebi.ac.uk/ols4/ontologies/pride)"),
                                     myAcqMeth2),
                    shiny::br())
    },
    if (isFnd["MS_quant_meth"]) {
      shiny::column(2,
                    shinyWidgets::pickerInput("MS_quant_meth",
                                              paste0("Select", quantMethTxt),
                                              availQuantMeth,
                                              myQuantMeth,
                                              TRUE,
                                              shinyWidgets::pickerOptions(title = "Search me",
                                                                          `live-search` = TRUE,
                                                                          actionsBox = TRUE,
                                                                          deselectAllText = "Clear search",
                                                                          showTick = TRUE)),
                    shiny::br())
    } else {
      shiny::column(2,
                    shiny::textInput("MS_quant_meth",
                                     paste0("Input", quantMethTxt, " (ontology = https://www.ebi.ac.uk/ols4/ontologies/pride)"),
                                     myQuantMeth2),
                    shiny::br())
    },
    if (isFnd["MS_label_meth"]) {
      shiny::column(2,
                    shinyWidgets::pickerInput("MS_label_meth",
                                              paste0("Select", labelMethTxt),
                                              availLabelMeth,
                                              myLabelMeth,
                                              TRUE,
                                              shinyWidgets::pickerOptions(title = "Search me",
                                                                          `live-search` = TRUE,
                                                                          actionsBox = TRUE,
                                                                          deselectAllText = "Clear search",
                                                                          showTick = TRUE)),
                    shiny::br())
    } else {
      shiny::column(2,
                    shiny::textInput("MS_label_meth",
                                     paste0("Input", labelMethTxt, " (ontology = https://www.ebi.ac.uk/ols4/ontologies/pride)"),
                                     myLabelMeth2),
                    shiny::br())
    },
    if (isFnd["Modifications"]) {
      shiny::column(2,
                    shinyWidgets::pickerInput("Modifications",
                                              paste0("Select", modsTxt),
                                              availMods,
                                              myPTMs,
                                              TRUE,
                                              shinyWidgets::pickerOptions(title = "Search me",
                                                                          `live-search` = TRUE,
                                                                          actionsBox = TRUE,
                                                                          deselectAllText = "Clear search",
                                                                          showTick = TRUE)),
                    shiny::br())
    } else {
      shiny::column(2,
                    shiny::textInput("Modifications",
                                     paste0("Input", modsTxt, " (ontology = https://www.ebi.ac.uk/ols4/ontologies/MOD)"),
                                     myPTMs2),
                    shiny::br())
    },
  ),
  shinyWidgets::actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
  shiny::br(),
  shiny::br()
)
serverA <- function(input, output, session) {
  VENDOR <- shiny::reactiveVal(myVendor)
  MSINSTR <- shiny::reactiveVal(myMSInstr)
  if (!tstOrg) {
    ORGA <- shiny::reactiveVal(mainOrg)
  }
  # MS instrument update function
  updtMSInstr <- function(reactive = TRUE) {
    if (reactive) {
      myVend <- VENDOR()
      myInstr <- MSINSTR()
    } else {
      myVend <- myVendor
      myInstr <- myMSInstr
    }
    availMS <- unlist(Vendors2Instr$Instruments[match(myVend, Vendors2Instr$Vendor)])
    myInstr <- myInstr[which(myInstr %in% availMS)]
    shiny::renderUI(
      if (isFnd["MS_models"]) {
        shinyWidgets::pickerInput("MS_instr",
                                  paste0("Select", MSInstrTxt),
                                  availMS,
                                  myInstr,
                                  TRUE,
                                  shinyWidgets::pickerOptions(title = "Search me",
                                                              `live-search` = TRUE,
                                                              actionsBox = TRUE,
                                                              deselectAllText = "Clear search",
                                                              showTick = TRUE))
      } else {
        shiny::textInput("MS_instr",
                         paste0("Input", MSInstrTxt, " (ontology = https://www.ebi.ac.uk/ols4/ontologies/ms/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FMS_1000031)"),
                         myMSInstr2)
      }
    )
  }
  updtDevStgs <- function(reactive = TRUE) {
    if (reactive) { myOrg <- ORGA() } else { myOrg <- mainOrg }
    devOnt <- devOntVal <- NA
    if (grepl("^D((\\.)|(anio))[ _\\-\\.]?rerio", myOrg)) { devOnt <- "Zebrafish" }
    if (grepl("^D((\\.)|(rosophila))[ _\\-\\.]?melanogaster", myOrg)) { devOnt <- "Fly" }
    if (grepl("^M((\\.)|(us))[ _\\-\\.]?musculus", myOrg)) { devOnt <- "Mouse" }
    if (grepl("^H((\\.)|(omo))[ _\\-\\.]?sapiens", myOrg)) { devOnt <- "Human" }
    if (grepl("^C((\\.)|(aenorhabditis))[ _\\-\\.]?elegans", myOrg)) { devOnt <- "Worms" }
    if (!is.na(devOnt)) {
      dvStgDF <- read.csv(paste0(homePath, "/DevStages_", devOnt, ".csv"))
      opt <- dvStgDF$Name
      if (exists("myDevStages")) {
        myDevStages <- unique(myDevStages)
        myDS1 <- myDevStages[which(myDevStages %in% opt)]
        myDS2 <- paste(myDevStages[which(!myDevStages %in% opt)], collapse = "|")
      }
    }
    shiny::renderUI(list(list(shiny::fluidRow(
      shiny::column(1),
      if (!is.na(devOnt)) {
        shiny::column(2,
                      shinyWidgets::pickerInput("devStages1",
                                                "Select stage(s)...",
                                                opt,
                                                myDS1,
                                                TRUE,
                                                shinyWidgets::pickerOptions(title = "Search me",
                                                                            `live-search` = TRUE,
                                                                            actionsBox = TRUE,
                                                                            deselectAllText = "Clear search",
                                                                            showTick = TRUE)),
                      shiny::br())
      },
      shiny::column(6,
                    shiny::textInput("devStages2",
                                     paste0("... or manually enter stage(s), separating multiple entries with '|'"),
                                     myDS2,
                                     "100%"))
    ))))
  }
  output$myMS_instr <- updtMSInstr(FALSE)
  output$myDevStages <- updtDevStgs(FALSE)
  if (isFnd["Cell_types"]) {
    shiny::observeEvent(input$Cell_type, {
      cltp <- unique(input$Cell_type)
      tmp <- c(cltp, availCellTypes[which(!availCellTypes %in% cltp)])
      shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                      "Cell_type",
                                      paste0("Select", cellTypeTxt),
                                      cltp,
                                      tmp)
    })
  }
  if (isFnd["Tissues"]) {
    shiny::observeEvent(input$Tissue_type, {
      tstp <- unique(input$Tissue_type)
      tmp <- c(tstp, availTissues[which(!availTissues %in% tstp)])
      shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                      "Tissue_type",
                                      paste0("Select", tissueTxt),
                                      tstp,
                                      tmp)
    })
  }
  if (isFnd["Modifications"]) {
    shiny::observeEvent(input$Modifications, {
      mds <- unique(input$Modifications)
      tmp <- c(mds, availMods[which(!availMods %in% mds)])
      shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                      "Modifications",
                                      paste0("Select", modsTxt),
                                      mds,
                                      tmp)
    })
  }
  if (isFnd["Diseases"]) {
    shiny::observeEvent(input$Diseases, {
      dis <- unique(input$Diseases)
      tmp <- c(dis, availDiseases[which(!availDiseases %in% dis)])
      shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                      "Diseases",
                                      paste0("Select", diseasesTxt),
                                      dis,
                                      tmp)
    })
  }
  if (isFnd["MS_models"]) {
    shiny::observeEvent(input$Vendor, {
      vnd <- unique(input$Vendor)
      if ("All vendors" %in% vnd) {
        tmp <- vnd <- availMSVend
      } else {
        tmp <- c(vnd, availMSVend[which(!availMSVend %in% vnd)])
      }
      VENDOR(vnd)
      if (isFnd["MS_models"]) {
        shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                        "Vendor",
                                        "Select instrument vendor(s)",
                                        vnd,
                                        tmp)
      }
      #}
      output$myMS_instr <- updtMSInstr()
    })
  }
  shiny::observeEvent(input$MS_instr, {
    instr <- unique(input$MS_instr)
    MSINSTR(instr)
    if (isFnd["MS_models"]) {
      tmp <- c(instr, availInstr[which(!availInstr %in% instr)])
      shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                      "MS_instr",
                                      paste0("Select", MSInstrTxt),
                                      instr,
                                      tmp)
    }
  })
  if (isFnd["MS_acq_meth"]) {
    shiny::observeEvent(input$MS_acq_meth, {
      acqmtd <- unique(input$MS_acq_meth)
      tmp <- c(acqmtd, availAcqMeth[which(!availAcqMeth %in% acqmtd)])
      shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                      "MS_acq_meth",
                                      paste0("Select", acqMethTxt),
                                      acqmtd,
                                      tmp)
    })
  }
  if (isFnd["MS_quant_meth"]) {
    shiny::observeEvent(input$MS_quant_meth, {
      qtmtd <- unique(input$MS_quant_meth)
      tmp <- c(qtmtd, availQuantMeth[which(!availQuantMeth %in% qtmtd)])
      shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                      "MS_quant_meth",
                                      paste0("Select", quantMethTxt),
                                      qtmtd,
                                      tmp)
    })
  }
  if (isFnd["MS_label_meth"]) {
    shiny::observeEvent(input$MS_label_meth, {
      lblmtd <- unique(input$MS_label_meth)
      tmp <- c(lblmtd, availLabelMeth[which(!availLabelMeth %in% lblmtd)])
      shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                      "MS_label_meth",
                                      paste0("Select", labelMethTxt),
                                      lblmtd,
                                      tmp)
    })
  }
  if (!tstOrg) {
    shiny::observeEvent(input$Organism, {
      ORGA(input$Organism)
      output$myDevStages <- updtDevStgs()
    })
  }
  #
  # Save
  shiny::observeEvent(input$saveBtn, {
    if (!tstOrg) {
      assign("mainOrg", unique(input$Organism), envir = .GlobalEnv)
    }
    assign("myCellType", unique(input$Cell_type), envir = .GlobalEnv)
    assign("myTissue", unique(input$Tissue_type), envir = .GlobalEnv)
    assign("myDisease", unique(input$Diseases), envir = .GlobalEnv)
    tmp <- unique(c(input$devStages1,
                    unlist(strsplit(input$devStages2, "\\|"))))
    tmp <- tmp[which(!is.na(tmp))]
    tmp <- tmp[which(nchar(tmp) > 0)]
    assign("myDevStages", tmp, envir = .GlobalEnv)
    assign("myVendor", unique(input$Vendor), envir = .GlobalEnv)
    assign("myMSInstr", unique(input$MS_instr), envir = .GlobalEnv)
    assign("myAcqMeth", unique(input$MS_acq_meth), envir = .GlobalEnv)
    assign("myQuantMeth", unique(input$MS_quant_meth), envir = .GlobalEnv)
    assign("myLabelMeth", unique(input$MS_label_meth), envir = .GlobalEnv)
    assign("myPTMs", unique(input$Modifications), envir = .GlobalEnv)
    assign("appRunTest", TRUE, envir = .GlobalEnv)
    shiny::stopApp()
  })
}
appTxtA <- gsub("myApp", "myAppA", gsub("\\(ui", "(uiA", gsub(", server", ", serverA", runApp)))
runKount <- 0
while ((!runKount)||(!exists("appRunTest"))) {
  eval(parse(text = appTxtA), envir = .GlobalEnv)
  runKount <- runKount+1
}
# Update myVendor
tmpVnd <- availMSVend[which(availMSVend != "All vendors")]
myVendor <- tmpVnd[which(vapply(tmpVnd, function(x) {
  w <- which(availInstrDF$Vendor == x)
  sum(myMSInstr %in% availInstrDF$Instrument[w])
}, 1) > 0)]
if ((length(myVendor))&&(sum(!myVendor %in% c("Bruker Daltonics", "Thermo Fisher Scientific")))) {
  stop("Vendor not yet supported!")
}
# Create duplicates which we will use now to map to our data
myOntologies <- c("mainOrg", myOntDF$myOntology)
for (i in myOntologies) {
  tmp <- unique(get(i))
  tmp <- tmp[which(tolower(tmp) != "unknown")]
  tmp <- tmp[which(!is.na(tmp))]
  tmp <- c(tmp, "not available")
  assign(paste0(i, 2), tmp)
}
#
if (scrptType == "withReps") {
  Frac.map2 <- Frac.map
}
if (scrptType == "noReps") {
  Frac.map2 <- FracMap
}
reload_SDRF <- exists("auld_SDRF")
if (reload_SDRF) {
  SDRF <- auld_SDRF
  reload_SDRF <- (nrow(SDRF) >= nrow(Frac.map2))&&("comment[data file]" %in% colnames(SDRF))
}
if (reload_SDRF) {
  tmp <- gsub(".*/", "", Frac.map2$`Raw file`)
  g <- grep("\\.d$", tmp)
  if (length(g)) {
    tmp[g] <- paste0(tmp[g], ".zip")
  }
  m <- match(tmp, SDRF$"comment[data file]")
  reload_SDRF <- sum(is.na(m)) == 0
  if (reload_SDRF) {
    SDRF <- SDRF[m,]
  } else {
    warning("Invalid SDRF backup reloaded, ignoring...")
    reload_SDRF <- FALSE
  }
}
nr <- nrow(Frac.map2)
if (scrptType == "withReps") {
  if (LabelType == "LFQ") {
    myKol <- colnames(Exp.map)[which(!colnames(Exp.map) %in% names(Aggregate.list))]
    myKol <- myKol[which(!myKol %in% c("Fractions", "Reference", "MQ.Exp", "Sample.name", "Use", "Ref.Sample.Aggregate"))]
    if (length(unique(Exp.map$Experiment)) == 1) { myKol <- myKol[which(myKol != "Experiment")] }
    tmp <- listMelt(Exp.map$MQ.Exp, 1:nrow(Exp.map))
    tmp[, myKol] <- Exp.map[tmp$L1, myKol]
    Frac.map2[, myKol] <- tmp[match(Frac.map2$MQ.Exp, tmp$value), myKol]
  }
  if (LabelType == "Isobaric") {
    tmp <- listMelt(Exp.map$MQ.Exp, 1:nrow(Exp.map), c("MQ.Exp", "ExpRow"))
    tmp <- aggregate(tmp$ExpRow, list(tmp$MQ.Exp), list)
    colnames(tmp) <- c("MQ.Exp", "ExpRow")
    m <- match(Frac.map$MQ.Exp, tmp$MQ.Exp)
    tmp <- tmp[m,]
    tmp$FracRow <- 1:nrow(tmp)
    tmp <- listMelt(tmp$ExpRow, tmp$FracRow, c("ExpRow", "FracRow"))
    k <- colnames(Frac.map2)
    k <- k[which(!k %in% c("Unique.Frac.ID", "Use", "MQ.Exp"))]
    tmp[, k] <- Frac.map2[tmp$FracRow, k]
    k <- colnames(Exp.map)
    k <- grep("^([A-Z][a-z]{2}){2,}$", k, value = TRUE, invert = TRUE)
    k <- k[which(!k %in% c("Fractions", "Ref.Sample.Aggregate", "Use", "MQ.Exp"))]
    tmp[, k] <- Exp.map[tmp$ExpRow, k]
    tmp$FracRow <- NULL
    tmp$ExpRow <- NULL
    tmp$`Isobaric label details` <- paste0(IsobarLab, gsub("^[A-Za-z]+", "", tmp$`Isobaric label details`))
    # Hard fix
    w <- which(tmp$`Isobaric label details` == "TMT126C")
    if (length(w)) { tmp$`Isobaric label details`[w] <- "TMT126" }
    #
    Frac.map2 <- tmp
  }
}
if (scrptType == "noReps") {
  myKol <- colnames(SamplesMap)
  myKol <- myKol[which(!myKol %in% c("Fractions", "Reference", "MQ.Exp", "Sample.name", "Use"))]
  myKol <- grep("__$", myKol, value = TRUE, invert = TRUE)
  tmp <- listMelt(SamplesMap$Experiment, 1:nrow(SamplesMap))
  tmp[, myKol] <- SamplesMap[tmp$L1, myKol]
  Frac.map2[, myKol] <- tmp[match(Frac.map2$`Parent sample`, tmp$value), myKol]
}
if ((reload_SDRF)&&(nrow(SDRF) != nrow(Frac.map2))) {
  rm(SDRF)
  reload_SDRF <- FALSE
}
nr <- nrow(Frac.map2) # Update just in case
if (!reload_SDRF) {
  SDRF <- data.frame("source name" = paste0("sample ", 1:nr),
                     "characteristics[organism]" = tolower(mainOrg),
                     check.names = FALSE)
}
if (!"assay name" %in% colnames(SDRF)) { SDRF$"assay name" <- paste0("Run ",  1:nrow(SDRF)) }
if (!"comment[data file]" %in% colnames(SDRF)) {
  SDRF$"comment[data file]" <- gsub(".*/", "", Frac.map2$`Raw file`)
  g <- grep("\\.d$", SDRF$"comment[data file]")
  if (length(g)) {
    SDRF$"comment[data file]"[g] <- paste0(SDRF$"comment[data file]"[g], ".zip")
  }
}
if (!"comment[technical replicate]" %in% colnames(SDRF)) {
  SDRF$"comment[technical replicate]" <- 1 # Please note so far none of our workflows support technical replicates!!!
}
if (!"comment[fraction identifier]" %in% colnames(SDRF)) { SDRF$"comment[fraction identifier]" <- Frac.map2$Fraction }
if (!"characteristics[biological replicate]" %in% colnames(SDRF)) {
  if ("Replicate" %in% colnames(Frac.map2)) {
    SDRF$"characteristics[biological replicate]" <- Frac.map2$Replicate
  } else {
    SDRF$"characteristics[biological replicate]" <- 1
  }
}
# if (!"comment[label]" %in% colnames(SDRF)) {
#   SDRF$"comment[label]" <- "label free" # Change this if doing TMT/SILAC!!!
# }
if (!"comment[cleavage agent details]" %in% colnames(SDRF)) {
  SDRF$"comment[cleavage agent details]" <- "NT=Trypsin; AC=MS:1001251; CS=(?⇐[KR])(?!P)" # Only trypsin supported so far!
}
#
SDRF2 <- SDRF[, c("comment[data file]", "comment[fraction identifier]"), drop = FALSE]
colnames(SDRF2) <- c("MS raw file", "Fraction")
if (LabelType == "Isobaric") {
  SDRF2 <- cbind(Frac.map2[, "Sample name", drop = FALSE],
                 SDRF2)
}
rws <- seq_len(nr)
wTest0 <- setNames(vapply(1:nrow(myOntDF), function(i) { #i <- 1
  x <- max(nchar(c(myOntDF$myOntology[i],
                   get(paste0(myOntDF$myOntology[i], 2))))) + 6
  x <- x*10
  if (is.na(x)) { x <- 15 } else { x <- max(c(ceiling(x/10)*10, 30)) }
  return(x)
}, 1), myOntDF$Name)
#
myOntDF$Candidate_col_root <- list(c("Cell Type", "Cell type"),
                                   c("Tissue", "Organ"),
                                   c("Disease", "Illness"),
                                   c("Developmental Stage", "Dev. Stage", "Developmental stage", "Dev. stage"),
                                   NA,
                                   NA,
                                   NA,
                                   c(NA, "Isobaric label details")[(LabelType == "Isobaric") + 1],
                                   "PTM-enriched")
myOntDF$Candidate_col <- lapply(1:nrow(myOntDF), function(i) { #i <- 8
  candKols <- myOntDF$Candidate_col_root[[i]]
  candKols <- candKols[which(!is.na(candKols))]
  if (length(candKols)) {
    candKols <- unique(c(candKols,
                         gsub(" ", "", candKols),
                         gsub(" ", ".", candKols),
                         gsub(" ", "_", candKols),
                         gsub(" ", "-", candKols)))
    candKols <- c(candKols, toupper(candKols), tolower(candKols))
    candKols <- candKols[which(candKols %in% colnames(Frac.map2))]
  }
  if (length(candKols)) {
    candKols <- candKols[which(vapply(candKols, function(k) { #k <- candKols[1]
      sum(!unique(Frac.map2[[k]]) %in% c(get(myOntDF$myOntology[i]), get(paste0(myOntDF$myOntology[i], "2"))))
    }, 1) == 0)]
  }
  return(candKols)
})
ALLIDS <- editOnts <- c()
wNotNA <- which(!is.na(myOntDF$Column))
#paste0("my", myOntDF$Name[wNotNA])
for (i in wNotNA) { #i <- wNotNA[1]
  nm <- myOntDF$Name[i]
  myOntVal <- get(paste0(myOntDF$myOntology[i], "2"))
  if (length(myOntVal) == 1) {
    SDRF[[myOntDF$Column[i]]] <- myOntVal
  } else {
    editOnts <- c(editOnts, nm)
    strtVal <- c()
    if ((reload_SDRF)&&(myOntDF$Column[i] %in% colnames(SDRF))) {
      strtVal <- SDRF[[myOntDF$Column[i]]]
    } else {
      if (length(myOntDF$Candidate_col[[i]]) == 1) {
        strtVal <- Frac.map2[[myOntDF$Candidate_col[[i]]]]
      }
    }
    strtVal <- strtVal[which(strtVal %in% myOntVal)]
    if (length(strtVal) != nr) { strtVal <- rep(myOntVal[1], nr) }
    SDRF2[[nm]] <- shinySelectInput(strtVal,
                                    nm,
                                    myOntVal,
                                    paste0(wTest0[nm], "px"),
                                    FALSE)
    SDRF2[[paste0(nm, "___FD")]] <- shinyFDInput(nm, nr)
    ALLIDS <- c(ALLIDS,
                paste0(nm, "___", rws))
  }
}
idsL <- length(ALLIDS)
#
if (idsL) {
  wTest1 <- vapply(colnames(SDRF2), function(k) { #k <- colnames(SDRF2)[1]
    if (nchar(k)) {
      if (k %in% names(wTest0)) { x <- wTest0[k] } else {
        k1 <- k
        if (k == "Sample name") { k1 <- "factor value[sample name]" }
        if (k == "MS raw file") { k1 <- "comment[data file]" }
        if (k == "Fraction") { k1 <- "comment[fraction identifier]" }
        if (k1 %in% colnames(SDRF)) {
          if (k == "MS raw file") {
            x <- max(nchar(c(k, gsub(".*/", "", SDRF[[k1]]))))*8 + 6
          } else {
            x <- max(nchar(c(k, SDRF[[k1]])))*8 + 6
          }
        } else {
          x <- 30
        }
      }
    } else { x <- 30 }
    return(x)
  }, 1)
  wTest2 <- sum(wTest1) + 15 + ncol(SDRF2)*5
  wTest1 <- paste0(as.character(wTest1), "px")
  wTest1 <- aggregate((1:length(wTest1))-1, list(wTest1), c)
  wTest1 <- apply(wTest1, 1, function(x) {
    x2 <- as.integer(x[[2]])
    list(width = x[[1]],
         targets = x2,
         names = colnames(SDRF2)[x2+1])
  })
  g <- grep("___((FD)|(INCR))$", colnames(SDRF2))
  colnames(SDRF2)[g] <- ""
  appNm <- paste0(dtstNm, " - SDRF editor")
  uiB <- fluidPage(
    useShinyjs(),
    setBackgroundColor( # Doesn't work
      color = c(#"#F8F8FF",
        "#EBEFF7"),
      gradient = "linear",
      direction = "bottom"
    ),
    extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
    #Dummy, hidden div to load arrow-down icon
    # tags$div(
    #   class = "container",
    #   style = "display: none;",
    #   tags$div(
    #     style = "margin-top: 50xp;",
    #     actionButton(
    #       "add_thing",
    #       label = "do it",
    #       class = "btn-success",
    #       icon = icon("arrow-down")
    #     )
    #   )
    # ),
    tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
    titlePanel(tag("u", "SDRF editor"),
               appNm),
    h2(dtstNm), 
    h3("Assign controlled vocabulary descriptors to MS runs."),
    br(),
    actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
    withSpinner(DT::DTOutput("SDRF_tbl", width = wTest2))
  )
  serverB <- function(input, output, session) {
    # Initialize output table
    SDRF3 <- SDRF # Output table
    #
    # Render dummy table
    output$SDRF_tbl <- DT::renderDT({ SDRF2 },
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
    #
    # Fill down
    sapply(1:idsL, function(x) { #x <- 1
      id1 <- ALLIDS[x]
      id2 <- paste0(ALLIDS[x], "___FD")
      tmp <- unlist(strsplit(id1, "___"))
      fct <- tmp[[1]]
      i <- as.integer(tmp[[2]])
      if (i < nr) {
        observeEvent(input[[id2]],
                     {
                       x <- input[[id1]]
                       myOntVal <- get(paste0("my", fct, "2"))
                       tp <- typeof(myOntVal)
                       if (typeof(x) != tp) { x <- get(paste0("as.", tp))(x) }
                       for (k in (i+1):nr) {
                         idK <- paste0(fct, "___", as.character(k))
                         updateSelectInput(shiny::getDefaultReactiveDomain(), idK, NULL, myOntVal, x)
                       }
                     })
        
      }
    })
    # Save
    observeEvent(input$saveBtn, {
      for (ont in editOnts) {
        typ <- typeof(get(paste0("my", ont, "2")))
        as.Fun <- get(paste0("as.", typ))
        SDRF3[[ont]] <- sapply(rws, function(i) {
          as.Fun(input[[paste0(ont, "___", i)]])
        })
      }
      assign("SDRF3", SDRF3, envir = .GlobalEnv)
      stopApp()
    })
    #observeEvent(input$cancel, { stopApp() })
    session$onSessionEnded(function() { stopApp() })
  }
  if (exists("SDRF3")) { rm(SDRF3) }
  appTxtB <- gsub("myApp", "myAppB", gsub("\\(ui", "(uiB", gsub(", server", ", serverB", runApp)))
  runKount <- 0
  while ((!runKount)||(!exists("appRunTest"))) {
    eval(parse(text = appTxtB), envir = .GlobalEnv)
    runKount <- runKount+1
  }
  SDRF[, myOntDF$Column[match(editOnts, myOntDF$Name)]] <- SDRF3[, editOnts]
}
#
for (k in myKol) {
  k1 <- tolower(k)
  if (k == "Replicate") { k1 <- "biological replicate" }
  SDRF[[paste0("characteristics[", k1, "]")]] <- Frac.map2[[k]]
  SDRF[[paste0("factor value[", k1, "]")]] <- Frac.map2[[k]]
}
#
# TO DO:
#   For DIA, it is recommended to add a "comment[MS1 scan range]" (example value = "400m/z - 1200m/z")
#   -> get this from the method!!!
#
# Column order matters...
kol1 <- colnames(SDRF)
kol2 <- c("source name",
          grep("^characteristics\\[", kol1, value = TRUE),
          "assay name",
          grep("^comment\\[", kol1, value = TRUE),
          grep("^factor value\\[", kol1, value = TRUE))
stopifnot(length(kol1) == length(kol2))
SDRF <- SDRF[, kol2]
if (LabelType == "Isobaric") {
  SDRF$`factor value[isobaric label]` <- NULL
  SDRF$`factor value[isobaric label details]` <- NULL
  SDRF$`characteristics[isobaric label]` <- NULL
  SDRF$`characteristics[isobaric label details]` <- NULL
}
#
#View(SDRF)
#writeClipboard(sdrfPth)
tst <- try(write.table(SDRF, sdrfPth, sep = "\t", quote = FALSE, row.names = FALSE), silent = TRUE)
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) { # We only want this to happen if the file is locked for editing
  dlg_message(paste0("File \"", sdrfPth, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.table(SDRF, sdrfPth, sep = "\t", quote = FALSE, row.names = FALSE), silent = TRUE)
}
#openxlsx::openXL(sdrfPth)
#
if (pyTest) {
  cat("Validating SDRF file...\n")
  cmd <- paste0("parse_sdrf validate-sdrf --sdrf_file \"", sdrfPth, "\"")
  #cat(cmd)
  system(cmd)
  #openwd(dirname(sdrfPth))
  cat("\n... but remember: it is always a good idea to also check it manually!\n")
}
