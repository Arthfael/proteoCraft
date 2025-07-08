
HEIGHT <- "500px"
allProt <- PG$Label
dfltPrt <- PG$Label[grsep2(prot.list[1], PG$"Protein IDs")]
appNm <- paste0(dtstNm, " - Data exploration")
if (!exists("plotLeatMaps")) { loadFun(paste0(wd, "/Clustering/HeatMaps.RData")) }
if (!exists("dimRedPlotLy")) { loadFun(paste0(wd, "/Dimensionality red. plots/DimRedPlots.RData")) }
if (!exists("QuantLy")) { loadFun(paste0(wd, "/Sorting plots/QuantPlots.RData")) }
if (!exists("ProfLy")) { loadFun(paste0(wd, "/Profile plots/ProfilePlots.RData")) }
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#F1F5DF"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", "Data exploration"),
             appNm),
  br(),
  em("Here you can explore you data."), br(),
  em("Click "), actionButton("exitBtn", "exit"), em(" to continue along the workflow."), br(),
  br(),
  #
  h4("Clustering heatmap"),
  selectInput("HeatMap", "Select heatmap", names(plotLeatMaps), "Global"),
  radioButtons("NormType", "Normalisation method:", c("Norm. by row", "Z-scored"), "Norm. by row"),
  fluidRow(withSpinner(plotlyOutput("MYplotLeatMap", height = HEIGHT))),
  #
  h4("Dimensionality reduction plots"),
  selectInput("DimRedPlot", "Select method",
              names(dimRedPlotLy)[which(names(dimRedPlotLy) != "Samples PCA")],
              "PCA"),
  fluidRow(
    column(8, withSpinner(plotlyOutput("PGDimRed", height = HEIGHT))),
    column(4, withSpinner(plotlyOutput("SmplsDimRed", height = HEIGHT)))
  ),
  #
  h4("Sorted Quantification plots"),
  fluidRow(column(2, selectInput("QuantType", "Select Quant method", QuantTypes, "Expression")),
           column(4, pickerInput("ProfPlotProt", "Select protein(s) to display", allProt, dfltPrt, TRUE,
                                 pickerOptions(title = "Search me",
                                               `live-search` = TRUE,
                                               actionsBox = TRUE,
                                               deselectAllText = "Clear search",
                                               showTick = TRUE)))),
  withSpinner(plotlyOutput("ProfPlot", height = HEIGHT)),
  # fluidRow(withSpinner(uiOutput("SortedPGPlot", height = HEIGHT))),
  br(),
  br()
)
server <- function(input, output, session) {
  #
  HEATMAP <- reactiveVal("Global")
  NORMTYPE <- reactiveVal("Norm. by row")
  DIMRED <- reactiveVal("PCA")
  QUANTTYPE <- reactiveVal("Expression")
  PROFPROT <- reactiveVal(prot.list[1])
  # SORTDSMPL <- reactiveVal(names(QuantLy[["Expression"]])[1])
  #
  updtHtMp <- function(reactive = TRUE) {
    if (reactive) {
      renderPlotly(plotLeatMaps[[HEATMAP()]][[NORMTYPE()]])
    } else {
      renderPlotly(plotLeatMaps[["Global"]][["Norm. by row"]])
    }
  }
  updtDimRed <- function(reactive = TRUE) {
    if (reactive) {
      renderPlotly(dimRedPlotLy[[DIMRED()]])
    } else {
      renderPlotly(dimRedPlotLy[["PCA"]])
    }
  }
  updtProfPlot <- function(reactive = TRUE) {
    if (reactive) {
      qunt <- QUANTTYPE()
      prt <- PROFPROT()
    } else {
      qunt <- "Expression"
      prt <- dfltPrt
    }
    pltly <- ProfLy[[qunt]]$grey
    if (length(prt)) {
      dat <- ProfLy[[qunt]]$data
      dat <- dat[which(dat$`Protein Group` %in% prt),]
      if (nrow(dat)) {
        #if (reactive) { print(PROFPROT()) }
        pltly <- add_trace(pltly, data = dat, x = ~Sample, y = ~Y,
                           split = ~`Protein Group`,
                           color = ~`Protein Group`, type = "scatter",
                           mode = "lines+markers", text = ~`Protein Group`, connectgaps = FALSE,
                           name = "")
      }
    }
    renderPlotly(pltly)
  }
  # updtSortPlot <- function(reactive = TRUE) {
  #   if (reactive) {
  #     ls <- list(list(
  #       selectInput("SortedSample", "Select sample", names(QuantLy[[QUANTTYPE()]]), SORTDSMPL()),
  #       fluidRow(withSpinner(plotlyOutput(renderPlotly(QuantLy[[QUANTTYPE()]][[SORTDSMPL()]]), height = HEIGHT)))
  #     ))
  #   } else {
  #     s <- list(list(
  #       selectInput("QuantPlot", "Select sample", names(QuantLy[["Expression"]]), names(QuantLy[["Expression"]])[1]),
  #       fluidRow(withSpinner(plotlyOutput(renderPlotly(QuantLy[["Expression"]][[1]]), height = HEIGHT)))
  #     ))
  #   }
  #   return(renderUI(ls))
  # }
  output$MYplotLeatMap <- updtHtMp(FALSE)
  output$PGDimRed <- updtDimRed(FALSE)
  output$SmplsDimRed <- renderPlotly(dimRedPlotLy[["Samples PCA"]])
  output$ProfPlot <- updtProfPlot(FALSE)
  # output$SortedPGPlot <- updtSortPlot(FALSE)
  observeEvent(input$HeatMap, {
    HEATMAP(input$HeatMap)
    output$MYplotLeatMap <- updtHtMp()
  })
  observeEvent(input$NormType, {
    NORMTYPE(input$NormType)
    output$MYplotLeatMap <- updtHtMp()
  })
  observeEvent(input$DimRedPlot, {
    DIMRED(input$DimRedPlot)
    output$PGDimRed <- updtDimRed()
  })
  observeEvent(input$QuantType, {
    QUANTTYPE(input$QuantType)
    output$ProfPlot <- updtProfPlot()
  })
  observeEvent(input$ProfPlotProt, {
    PROFPROT(input$ProfPlotProt)
    output$ProfPlot <- updtProfPlot()
  })
  # observeEvent(input$SortedSample, {
  #   SORTDSMPL(input$SortedSample)
  #   output$SortedPGPlot <- updtSortPlot()
  # })
  # Save
  observeEvent(input$exitBtn, { stopApp() })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
eval(parse(text = runApp), envir = .GlobalEnv)
#
