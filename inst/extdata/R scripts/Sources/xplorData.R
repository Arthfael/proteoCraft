####################
# Explore the data #
####################
#
# Opens a shiny app where clustering heatmaps, dimensionality reduction plots, protein profile plots or sorting plots may be visualized.
#
HEIGHT <- "500px"
# Defaults
allProt <- PG$Label
dfltPrt <- PG$Label[grsep2(prot.list, PG$"Protein IDs")]
if (!length(dfltPrt)) { dfltPrt <- PG$Label[1] }
appNm <- paste0(dtstNm, " - Data exploration")
if (!exists("plotLeatMaps")) { loadFun(paste0(wd, "/Clustering/HeatMaps.RData")) }
if (!exists("dimRedPlotLy")) { loadFun(paste0(wd, "/Dimensionality red. plots/DimRedPlots.RData")) }
if (!exists("QuantLy")) { loadFun(paste0(wd, "/Sorting plots/QuantPlots.RData")) }
if (!exists("ProfLy")) { loadFun(paste0(wd, "/Profile plots/ProfilePlots.RData")) }
nmsDmRds <- names(dimRedPlotLy)[which(names(dimRedPlotLy) != "Samples PCA")]
dfltHtMp <- c("Global", names(plotLeatMaps))
dfltHtMp <- dfltHtMp[which(dfltHtMp %in% names(plotLeatMaps))[1]]
dfltNrmTp <- c("Norm. by row", names(plotLeatMaps[[dfltHtMp]]))
dfltNrmTp <- dfltNrmTp[which(dfltNrmTp %in% names(plotLeatMaps[[dfltHtMp]]))[1]]
dfltDmRd <- c("PCA", nmsDmRds)
dfltDmRd <- dfltDmRd[which(dfltDmRd %in% nmsDmRds)[1]]
dfltQntTp <- c("Expression", names(QuantLy))
dfltQntTp <- dfltQntTp[which(dfltQntTp %in% names(QuantLy))[1]]
#
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
  selectInput("HeatMap", "Select heatmap", names(plotLeatMaps), dfltHtMp),
  radioButtons("NormType", "Normalisation method:", c("Norm. by row", "None"), "None"),
  fluidRow(withSpinner(plotlyOutput("MYplotLeatMap", height = HEIGHT))),
  #
  h4("Dimensionality reduction plots"),
  selectInput("DimRedPlot", "Select method", nmsDmRds, dfltDmRd),
  fluidRow(
    column(8, withSpinner(plotlyOutput("PGDimRed", height = HEIGHT))),
    column(4, withSpinner(plotlyOutput("SmplsDimRed", height = HEIGHT)))
  ),
  #
  h4("Sorted Quantification plots"),
  fluidRow(column(2, selectInput("QuantType", "Select Quant method", QuantTypes, dfltQntTp)),
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
  HEATMAP <- reactiveVal(dfltHtMp)
  NORMTYPE <- reactiveVal(dfltNrmTp)
  DIMRED <- reactiveVal(dfltDmRd)
  QUANTTYPE <- reactiveVal(dfltQntTp)
  PROFPROT <- reactiveVal(prot.list[1])
  # SORTDSMPL <- reactiveVal(names(QuantLy[[dfltQntTp]])[1])
  #
  updtHtMp <- function(reactive = TRUE) {
    if (reactive) {
      renderPlotly(plotLeatMaps[[HEATMAP()]][[NORMTYPE()]]$Plot)
    } else {
      renderPlotly(plotLeatMaps[[dfltHtMp]][[dfltNrmTp]]$Plot)
    }
  }
  updtDimRed <- function(reactive = TRUE) {
    if (reactive) {
      renderPlotly(suppressMessages(dimRedPlotLy[[DIMRED()]]))
    } else {
      renderPlotly(suppressMessages(dimRedPlotLy[[dfltDmRd]]))
    }
  }
  updtProfPlot <- function(reactive = TRUE) {
    if (reactive) {
      qunt <- QUANTTYPE()
      prt <- PROFPROT()
    } else {
      qunt <- dfltQntTp
      prt <- dfltPrt
    }
    if (length(prt)) {
      ggCall <- gsub(", +alpha += +[^\\)]+\\)", ")", ProfLy[[qunt]]$ggCall)
      plCall <- ProfLy[[qunt]]$plCall
      profData <- ProfLy[[qunt]]$data
      profData <- profData[which(profData$`Protein Group` %in% prt),]
      if (nrow(profData)) {
        ttl <- ProfLy[[qunt]]$Ttl
        kolnm <- ProfLy[[qunt]]$data2$colName
        myFacets <- ProfLy[[qunt]]$data2$facets
        Ngl <- ProfLy[[qunt]]$data2$angle
        # Now that the original objects have their proper values again, we can evaluate
        eval(parse(text = ggCall))
        eval(parse(text = plCall))
        renderPlotly(plotlyProfiles)
      }
    }
  }
  # updtSortPlot <- function(reactive = TRUE) {
  #   if (reactive) {
  #     ls <- list(list(
  #       selectInput("SortedSample", "Select sample", names(QuantLy[[QUANTTYPE()]]), SORTDSMPL()),
  #       fluidRow(withSpinner(plotlyOutput(renderPlotly(QuantLy[[QUANTTYPE()]][[SORTDSMPL()]]), height = HEIGHT)))
  #     ))
  #   } else {
  #     s <- list(list(
  #       selectInput("QuantPlot", "Select sample", names(QuantLy[[dfltQntTp]]), names(QuantLy[[dfltQntTp]])[1]),
  #       fluidRow(withSpinner(plotlyOutput(renderPlotly(QuantLy[[dfltQntTp]][[1]]), height = HEIGHT)))
  #     ))
  #   }
  #   return(renderUI(ls))
  # }
  output$MYplotLeatMap <- updtHtMp(FALSE)
  output$PGDimRed <- updtDimRed(FALSE)
  output$SmplsDimRed <- renderPlotly(suppressMessages(dimRedPlotLy[["Samples PCA"]]))
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
