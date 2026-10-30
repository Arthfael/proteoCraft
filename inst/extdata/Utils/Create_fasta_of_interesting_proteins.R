require(proteoCraft)

msg <- "Select fasta file"
#fastaFl <- normalizePath(choose.files("D:/Fasta_databases/*.fa", msg, FALSE), winslash = "/")
fastaFl <- rstudioapi::selectFile(msg,
                                  path = "D:/Fasta_databases/*",
                                  filter = "fasta file (*.fasta|*.fas|*.fa|*.faa|*.fasta.fas|*.txt)")
wd <- dirname(fastaFl)
setwd(wd)
Sources <- c("UniProtKB", "Ensembl", "RefSeq-RNA", "RefSeq-Protein", "RefSeq-CDS", "NCBI", "TAIR", "custom")
dbSource <- svDialogs::dlg_list(Sources, Sources[1L], title = "Where is the Fasta file from?")$res
db <- Format.DB(fastaFl, dbSource)
db <- db[grep("^>rev_", db$Header, invert = TRUE),]
contRmv <- c(TRUE, FALSE)[match(svDialogs::dlg_message("Remove contaminant entries (\"CON__...\")?", "yesno")$res,
                                c("yes", "no"))]
if (contRmv) {
  db <- db[grep("^CON__", db$"Protein ID", invert = TRUE),]
}
frgRmv <- c(TRUE, FALSE)[match(svDialogs::dlg_message("Remove incomplete (\"Fragment\") entries?", "yesno")$res,
                               c("yes", "no"))]
if (frgRmv) {
  db <- db[grep("\\(Fragment\\)", db$Header, invert = TRUE),]
}

dbOrd <- 1L:nrow(db)
protDeflt <- NULL
protHeads <- sub("^>", "", db$Header[dbOrd])
ui <- shiny::fluidPage(
  shinyjs::useShinyjs(),
  shiny::titlePanel("Select proteins to extract and write into a smaller fasta", "Proteins of interest"),
  shinyWidgets::pickerInput("IntProt", NULL, protHeads, protDeflt, TRUE,
                            shinyWidgets::pickerOptions(title = "Search me",
                                                        `live-search` = TRUE,
                                                        actionsBox = TRUE,
                                                        deselectAllText = "Clear search")),
  shiny::br(),
  shiny::textInput("Destination", "Destination file name", "Proteins_of_interest.fasta", "600px"),
  shiny::checkboxInput("Overwrite", "Overwrite conflicts?", FALSE),
  shiny::br(),
  shiny::actionButton("saveBtn", "Save"),
  shiny::span(shiny::uiOutput("Msg"), style = "color:red"),
  shiny::br()
)
server <- \(input, output, session) {
  # Proteins of interest
  output$Msg <- shiny::renderUI({ shiny::em(" ") })
  shiny::observeEvent(input$IntProt, {
    assign("Prot.list", db$`Protein ID`[dbOrd][match(input$IntProt, protHeads)], envir = .GlogbalEnv)
  })
  shiny::observeEvent(input$Destination, {
    if (nchar(input$Destination)) {
      if (input$Overwrite || (!file.exists(paste0(wd, "/", input$Destination)))) {
        shinyjs::enable("saveBtn")
        output$Msg <- shiny::renderUI({ shiny::em(" ") })
      }
      if ((!input$Overwrite) && file.exists(paste0(wd, "/", input$Destination))) {
        shinyjs::disable("saveBtn")
        output$Msg <- shiny::renderUI({ shiny::em("File already exist! Change file name or tick overwrite!") })
      }
    }
    if (!nchar(input$Destination)) {
      shinyjs::disable("saveBtn")
      output$Msg <- shiny::renderUI({ shiny::em("Enter a valid file name!") })
    }
  })
  shiny::observeEvent(input$Overwrite, {
    if (input$Overwrite && nchar(input$Destination)) {
      shinyjs::enable("saveBtn")
      output$Msg <- shiny::renderUI({ shiny::em(" ") })
    }
    if ((!input$Overwrite) && file.exists(paste0(wd, "/", input$Destination))) {
      shinyjs::disable("saveBtn")
      output$Msg <- shiny::renderUI({ shiny::em("File already exist! Change file name or tick overwrite!") })
    }
  })
  shiny::observeEvent(input$saveBtn, {
    assign("flNm", input$Destination, envir = .GlogbalEnv)
    dbFilt <- db[match(Prot.list, db$`Protein ID`),]
    msg <- apply(dbFilt[, c("Protein ID", "Common Name")], 1L, paste, collapse = " = ")
    msg <- paste0("\nProteins selected:\n", paste0(" - ", msg, collapse = "\n"), "\n")
    cat(msg)
    writeFasta(dbFilt, paste0(wd, "/", flNm))
    shiny::stopApp()
  })
  session$onSessionEnded(\() { shiny::stopApp() })
}
print(shiny::shinyApp(ui, server, options = list(launch.browser = TRUE)))
openwd()
