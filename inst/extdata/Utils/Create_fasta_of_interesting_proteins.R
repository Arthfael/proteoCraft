require(proteoCraft)
require(svDialogs)
require(shiny)
require(shinyWidgets)
require(shinyjs)

fastaFl <- normalizePath(choose.files("D:/Fasta_databases/*.fa", "Select fasta file", FALSE), winslash = "/")
wd <- dirname(fastaFl)
setwd(wd)
Sources <- c("UniProtKB", "Ensembl", "RefSeq-RNA", "RefSeq-Protein", "RefSeq-CDS", "NCBI", "TAIR", "custom")
dbSource <- dlg_list(Sources, Sources[1], title = "Where is the Fasta file from?")$res
db <- Format.DB(fastaFl, dbSource)
db <- db[grep("^>rev_", db$Header, invert = TRUE),]
contRmv <- c(TRUE, FALSE)[match(dlg_message("Remove contaminant entries (\"CON__...\")?", "yesno")$res, c("yes", "no"))]
if (contRmv) {
  db <- db[grep("^CON__", db$"Protein ID", invert = TRUE),]
}
frgRmv <- c(TRUE, FALSE)[match(dlg_message("Remove incomplete (\"Fragment\") entries?", "yesno")$res, c("yes", "no"))]
if (frgRmv) {
  db <- db[grep("\\(Fragment\\)", db$Header, invert = TRUE),]
}

dbOrd <- 1:nrow(db)
protDeflt <- NULL
protHeads <- gsub("^>", "", db$Header[dbOrd])
ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel("Select proteins to extract and write into a smaller fasta", "Proteins of interest"),
  pickerInput("IntProt", NULL, protHeads, protDeflt, TRUE,
              pickerOptions(title = "Search me",
                            `live-search` = TRUE,
                            actionsBox = TRUE,
                            deselectAllText = "Clear search")),
  br(),
  textInput("Destination", "Destination file name", "Proteins_of_interest.fasta", "600px"),
  checkboxInput("Overwrite", "Overwrite conflicts?", FALSE),
  br(),
  actionButton("saveBtn", "Save"),
  span(uiOutput("Msg"), style = "color:red"),
  br()
)
server <- function(input, output, session) {
  # Proteins of interest
  output$Msg <- renderUI({ em(" ") })
  observeEvent(input$IntProt, {
    Prot.list <<- db$`Protein ID`[dbOrd][match(input$IntProt, protHeads)]
  })
  observeEvent(input$Destination, {
    if (nchar(input$Destination) > 0) {
      if ((input$Overwrite)||(!file.exists(paste0(wd, "/", input$Destination)))) {
        shinyjs::enable("saveBtn")
        output$Msg <- renderUI({ em(" ") })
      }
      if ((!input$Overwrite)&&(file.exists(paste0(wd, "/", input$Destination)))) {
        shinyjs::disable("saveBtn")
        output$Msg <- renderUI({ em("File already exist! Change file name or tick overwrite!") })
      }
    }
    if (nchar(input$Destination) == 0) {
      shinyjs::disable("saveBtn")
      output$Msg <- renderUI({ em("Enter a valid file name!") })
    }
  })
  observeEvent(input$Overwrite, {
    if ((input$Overwrite)&&(nchar(input$Destination) > 0)) {
      shinyjs::enable("saveBtn")
      output$Msg <- renderUI({ em(" ") })
    }
    if ((!input$Overwrite)&&(file.exists(paste0(wd, "/", input$Destination)))) {
      shinyjs::disable("saveBtn")
      output$Msg <- renderUI({ em("File already exist! Change file name or tick overwrite!") })
    }
  })
  observeEvent(input$saveBtn, {
    flNm <<- input$Destination
    dbFilt <- db[match(Prot.list, db$`Protein ID`),]
    msg <- apply(dbFilt[, c("Protein ID", "Common Name")], 1, paste, collapse = " = ")
    msg <- paste0("\nProteins selected:\n", paste0(" - ", msg, collapse = "\n"), "\n")
    cat(msg)
    writeFasta(dbFilt, paste0(wd, "/", flNm))
    stopApp()
  })
  session$onSessionEnded(function() { stopApp() })
}
print(shinyApp(ui, server, options = list(launch.browser = TRUE)))
openwd()
