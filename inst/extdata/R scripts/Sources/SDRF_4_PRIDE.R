# This is an embryo script to write SDRF files for PRIDE submissions.
# For now it is probably incomplete, and only works for label-free datasets
# There will probably have to be additions for phospho-enriched samples too
# This could get included in the two main workflows, though I think it's better to keep it separate for now...
# Also the path to the validation tool is hard coded for now
# Assumes Trypsin
# Assumes single instrument type for all
# Assumes...
# Just... check the output!!!

require(svDialogs)
require(proteoCraft)

#load_Bckp()

setwd(wd)
if (LabelType %in% c("LFQ", "DIA")) { # Actually for DIA experiments the value should be "LFQ", not "DIA"!
  w <- which(c("Organism_Full", "Organism") %in% colnames(db))
  tstorg %<o% (length(w) > 0)
  if (tstorg) {
    OrgKol %<o% c("Organism_Full", "Organism")[w[1]]
    tst <- gsub(" *(\\(|\\[).*", "", db[[OrgKol]])
    tst <- aggregate(tst, list(tst), length)
    tst <- tst[order(tst$x, decreasing = TRUE),]
    mainOrg %<o% tst$Group.1[1]
  }
  fls <- data.frame(File = paste0(homePath, "/", c("Tissues",
                                                   "Modifications",
                                                   "Diseases",
                                                   "Cell_types",
                                                   "PRIDE",
                                                   "MS_Quant_meth",
                                                   "MS_models"), ".csv"))
  fls$Name <- gsub("\\.csv$", "", basename(fls$File))
  fls$Exists <- file.exists(fls$File)
  isFnd <- setNames(fls$Exists, fls$Name)
  if (isFnd["Cell_types"]) {
    tmp <- data.table::fread(fls$File[match("Cell_types", fls$Name)])
    availCellTypes <- sort(unique(tmp$label))
    cellTypeTxt <- " the cell type analyzed"
    if ((!exists("myTissue"))||(!is.character(myTissue))||(sum(!myTissue %in% availCellTypes))) {
      myCellType <- c()
      myCellType2 <- "not available"
    } else {
      myCellType2 <- myCellType
    }
    availCellTypes <- c(myCellType, availCellTypes[which(!availCellTypes %in% myCellType)])
  }
  if (isFnd["Tissues"]) {
    tmp <- data.table::fread(fls$File[match("Tissues", fls$Name)])
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
    tmp <- data.table::fread(fls$File[match("Modifications", fls$Name)])
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
    tmp <- data.table::fread(fls$File[match("Diseases", fls$Name)])
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
    tmp <- data.table::fread(fls$File[match("MS_models", fls$Name)])
    availInstr <- tmp$Instrument
    Vendors2Instr <- aggregate(tmp$Instrument, list(tmp$Vendor), function(x) { sort(unique(x)) })
    colnames(Vendors2Instr) <- c("Vendor", "Instruments")
    availMSVend <- c("All vendors", sort(Vendors2Instr$Vendor))
    tmp2 <- data.frame(Vendor = "All vendors")
    tmp2$Instruments <- list(sort(unique(unlist(Vendors2Instr$Instruments))))
    Vendors2Instr <- rbind(Vendors2Instr, tmp2)
    MSInstrTxt <- " the Mass Spectrometer model(s) on which the data was acquired"
    if ((!exists("myVendor"))||(!is.character(myVendor))||(sum(!myVendor %in% Vendors2Instr$Vendor))) { myVendor <- "All vendors" }
    if ((!exists("myMSInstr"))||(!is.character(myMSInstr))||(sum(!myMSInstr %in% availInstr))) {
      myMSInstr <- c()
      myMSInstr2 <- "not available"
    } else {
      myMSInstr2 <- myMSInstr
    }
    availInstr <- c(myMSInstr, availInstr[which(!availInstr %in% myMSInstr)])
  }
  if (isFnd["MS_Quant_meth"]) {
    tmp <- data.table::fread(fls$File[match("MS_Quant_meth", fls$Name)])
    availQuantMeth <- sort(unique(tmp$label))
    quantMethTxt <- " relevant MS quantification method(s)"
    if ((!exists("myQuantMeth"))||(!is.character(myQuantMeth))||(!sum(myQuantMeth %in% availQuantMeth))) {
      myQuantMeth <- c()
      myQuantMeth2 <- "not available"
    } else {
      myQuantMeth2 <- myQuantMeth
    }
    availQuantMeth <- c(myQuantMeth, availQuantMeth[which(!availQuantMeth %in% myQuantMeth)])
  }
  appNm <- paste0(dtstNm, " - SDRF editor")
  ui2 <- shiny::fluidPage(
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
    if (!tstorg) {
      shiny::column(4,
                    shiny::tags$hr(style = "border-color: black;"),
                    shiny::textInput("Organism",
                                     "Enter parent organism (ontology = https://www.ebi.ac.uk/ols4/ontologies/NCBITaxon)",
                                     "Homo sapiens"),
                    shiny::br())
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
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
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
      if (isFnd["MS_Quant_meth"]) {
        shiny::column(4,
                      shinyWidgets::pickerInput("MS_Quant_meth",
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
        shiny::column(4,
                      shiny::textInput("MS_Quant_meth",
                                       paste0("Input", quantMethTxt, " (ontology = https://www.ebi.ac.uk/ols4/ontologies/pride)"),
                                       myQuantMeth2),
                      shiny::br())
      },
      if (isFnd["Modifications"]) {
        shiny::column(4,
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
        shiny::column(4,
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
  server2 <- function(input, output, session) {
    VENDOR <- shiny::reactiveVal(myVendor)
    MSINSTR <- shiny::reactiveVal(myMSInstr)
    # MS instrument update function
    updtMSInstr <- function(reactive = TRUE) {
      if (reactive) { myVend <- VENDOR() } else { myVend <- myVendor }
      if (reactive) { myInstr <- MSINSTR() } else { myInstr <- myMSInstr }
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
    output$myMS_instr <- updtMSInstr(FALSE)
    if (isFnd["Cell_types"]) {
      shiny::observeEvent(input$Cell_type, {
        tmp <- c(input$Cell_type, availCellTypes[which(!availCellTypes %in% input$Cell_type)])
        shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                        "Cell_type",
                                        paste0("Select", cellTypeTxt),
                                        input$Cell_type,
                                        tmp)
      })
    }
    if (isFnd["Tissues"]) {
      shiny::observeEvent(input$Tissue_type, {
        tmp <- c(input$Tissue_type, availTissues[which(!availTissues %in% input$Tissue_type)])
        shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                        "Tissue_type",
                                        paste0("Select", tissueTxt),
                                        input$Tissue_type,
                                        tmp)
      })
    }
    if (isFnd["Modifications"]) {
      shiny::observeEvent(input$Modifications, {
        tmp <- c(input$Modifications, availMods[which(!availMods %in% input$Modifications)])
        shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                        "Modifications",
                                        paste0("Select", modsTxt),
                                        input$Modifications,
                                        tmp)
      })
    }
    if (isFnd["Diseases"]) {
      shiny::observeEvent(input$Diseases, {
        tmp <- c(input$Diseases, availDiseases[which(!availDiseases %in% input$Diseases)])
        shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                        "Diseases",
                                        paste0("Select", diseasesTxt),
                                        input$Diseases,
                                        tmp)
      })
    }
    if (isFnd["MS_models"]) {
      shiny::observeEvent(input$Vendor, {
        vnd <- input$Vendor
        #vnd <- vnd[which(vnd != "All vendors")]
        #avVnd <- availMSVend[which(availMSVend != "All vendors")]
        #w <- which(!avVnd %in% vnd)
        #if ((!length(w))||(!length(vnd))) { vnd <- "All vendors" }
        #print(vnd)
        VENDOR(vnd)
        #if (!length(vnd)) { vnd <- NULL }
        #if (!identical(vnd, input$Vendor)) {
        if (isFnd["MS_models"]) {
          shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                          "Vendor",
                                          "Select instrument vendor(s)",
                                          vnd,
                                          availMSVend)
        }
        #}
        output$myMS_instr <- updtMSInstr()
      })
    }
    shiny::observeEvent(input$MS_instr, {
      if (isFnd["MS_models"]) {
        tmp <- c(input$MS_instr, availInstr[which(!availInstr %in% input$MS_instr)])
        shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                        "MS_instr",
                                        paste0("Select", MSInstrTxt),
                                        input$MS_instr,
                                        tmp)
      }
      MSINSTR(input$MS_instr)
    })
    if (isFnd["MS_Quant_meth"]) {
      shiny::observeEvent(input$MS_Quant_meth, {
        tmp <- c(input$MS_Quant_meth, availQuantMeth[which(!availQuantMeth %in% input$MS_Quant_meth)])
        shinyWidgets::updatePickerInput(shiny::getDefaultReactiveDomain(),
                                        "MS_Quant_meth",
                                        paste0("Select", quantMethTxt),
                                        input$MS_Quant_meth,
                                        tmp)
      })
    }
    #
    # Save
    shiny::observeEvent(input$saveBtn, {
      if (!tstorg) {
        assign("mainOrg", input$Organism, envir = .GlobalEnv)
      }
      assign("myCellType", input$Cell_type, envir = .GlobalEnv)
      assign("myTissue", input$Tissue_type, envir = .GlobalEnv)
      assign("myDisease", input$Diseases, envir = .GlobalEnv)
      assign("myVendor", input$Vendor, envir = .GlobalEnv)
      assign("myMSInstr", input$MS_instr, envir = .GlobalEnv)
      assign("myQuantMeth", input$MS_Quant_meth, envir = .GlobalEnv)
      assign("myPTMs", input$Modifications, envir = .GlobalEnv)
      assign("appRunTest", TRUE, envir = .GlobalEnv)
      shiny::stopApp()
    })
  }
  appTxt2 <- gsub("myApp", "myApp2", gsub("\\(ui", "(ui2", gsub(", server", ", server2", runApp)))
  runKount <- 0
  while ((!runKount)||(!exists("appRunTest"))) {
    eval(parse(text = appTxt2), envir = .GlobalEnv)
    runKount <- runKount+1
  }
  res <- c("mainOrg", "myCellType", "myTissue", "myDisease", "myMSInstr", "myQuantMeth", "myPTMs")
  for (i in res) {
    tmp <- get(i)
    if (!length(tmp)) { tmp <- "not available" }
    assign(paste0(i, 2), tmp)
  }
  #
  if (scrptType == "withReps") {
    L <- nrow(Frac.map)
    Frac.map2 <- Frac.map
    kol <- colnames(Exp.map)[which(!colnames(Exp.map) %in% names(Aggregate.list))]
    kol <- kol[which(!kol %in% c("Fractions", "Reference", "MQ.Exp", "Sample.name", "Use", "Ref.Sample.Aggregate"))]
    if (length(unique(Exp.map$Experiment)) == 1) { kol <- kol[which(kol != "Experiment")] }
    tmp <- listMelt(Exp.map$MQ.Exp, 1:nrow(Exp.map))
    tmp[, kol] <- Exp.map[tmp$L1, kol]
    Frac.map2[, kol] <- tmp[match(Frac.map2$MQ.Exp, tmp$value), kol]
    SDRF <- data.frame("source name" = paste0("sample ", 1:L),
                       "characteristics[organism]" = tolower(mainOrg),
                       check.names = FALSE)
  }
  if (scrptType == "noReps") {
    L <- nrow(FracMap)
    Frac.map2 <- FracMap
    kol <- colnames(SamplesMap)
    kol <- kol[which(!kol %in% c("Fractions", "Reference", "MQ.Exp", "Sample.name", "Use"))]
    kol <- grep("__$", kol, value = TRUE, invert = TRUE)
    tmp <- listMelt(SamplesMap$Experiment, 1:nrow(SamplesMap))
    tmp[, kol] <- SamplesMap[tmp$L1, kol]
    Frac.map2[, kol] <- tmp[match(Frac.map2$`Parent sample`, tmp$value), kol]
    SDRF <- data.frame("source name" = paste0("sample ", 1:L),
                       "characteristics[organism]" = tolower(mainOrg),
                       check.names = FALSE)
  }
  #
  # TO DO:
  # For now we can deal with only one value, however we should create a shiny app (based on Exp Map editor) allowing,
  # wherever there are multiple possibilities, to assign one to each row dynamically!!!
  #
  # PTMs and quant meth are not used for now: do they need to be in the SDRF?
  #
  if ("Tissue" %in% kol) { SDRF$"characteristics[organism part]" <- Frac.map2$Tissue } else {
    SDRF$"characteristics[organism part]" <- myTissue2
  }
  kol <- kol[which(!kol %in% c("Tissue"))]
  if ("Replicate" %in% kol) { SDRF$"characteristics[biological replicate]" <- Frac.map2$Replicate } else {
    SDRF$"characteristics[biological replicate]" <- 1
  }
  kol <- kol[which(!kol %in% c("Replicate"))]
  if ("Disease" %in% kol) { SDRF$"characteristics[disease]" <- Frac.map2$Disease } else {
    SDRF$"characteristics[disease]" <- myDisease2
  }
  kol <- kol[which(!kol %in% c("Disease"))]
  if ("Cell type" %in% kol) { SDRF$"characteristics[cell type]" <- Frac.map2$"Cell type" } else {
    SDRF$"characteristics[cell type]" <- myCellType2
  }
  kol <- kol[which(!kol %in% c("Cell type"))]
  if ("Developmental stage" %in% kol) { SDRF$"characteristics[developmental stage]" <- Frac.map2$"Developmental stage" } else {
    #
    # TO DO: add to the detected ontologies + 1st shiny app
    #
    SDRF$"characteristics[developmental stage]" <- dlg_input("At which development stage are the samples?", "not available")$res
  }
  kol <- kol[which(!kol %in% c("Developmental stage"))]
  SDRF$"assay name" <- paste0("Run ", 1:L)
  SDRF$"comment[technical replicate]" <- 1
  SDRF$"comment[fraction identifier]" <- Frac.map2$Fraction
  #
  #
  #
  #
  #
  SDRF$"comment[label]" <- "label free" # Change this if doing TMT/SILAC!!!
  #
  #
  #
  #
  #
  SDRF$"comment[data file]" <- paste0(Frac.map2$`Raw files name`, ".d.zip")
  if (!exists("homePath")) {
    homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
  }
  SDRF$"comment[instrument]" <- myMSInstr2
  #
  SDRF$"comment[cleavage agent details]" <- "NT=Trypsin; AC=MS:1001251; CS=(?⇐[KR])(?!P)"
  if (length(kol)) { for (k in kol) { SDRF[[paste0("factor value[", tolower(k), "]")]] <- Frac.map2[[k]] } }
  #View(SDRF)
  sdrfPth <- paste0(outdir, "/SDRF.tsv")
  #writeClipboard(sdrfPth)
  write.table(SDRF, sdrfPth, sep = "\t", quote = FALSE, row.names = FALSE)
  #openxlsx::openXL(sdrfPth)
  #
  #
  # TO DO: for validations, make sure that python is installed and during Configure() make sure the parse_sdrf.exe script is added!!!
  fl <- paste0("", gsub("/Documents$", "", gsub("\\\\", "/", Sys.getenv("HOME"))),
               "/AppData/Roaming/Python/Python310/Scripts/parse_sdrf.exe")
  if (file.exists(fl)) {
    cat("Validating SDRF file...\n")
    cmd <- paste0("\"", fl, "\" validate-sdrf --sdrf_file \"", sdrfPth, "\"")
    #cat(cmd)
    system(cmd)
    #openwd(dirname(sdrfPth))
  } else {
    cat("To allow this workflow to automatically validate this file, make sure that a) python is installed and b) follow installation instructions for https://github.com/bigbio/sdrf-pipelines - the current script will then be able to run the validation tool directly...")
  }
  cat("\n... but remember: it's always a good idea to also check it manually!\n")
} else {
  # That part would be easy to write but isn't done yet
  warning("Writing SDRF files for is only available for LFQ experiments for the time being...")  
}
