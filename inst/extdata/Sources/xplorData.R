####################
# Explore the data #
####################
#
# Opens a shiny app where clustering heatmaps, dimensionality reduction plots, protein profile plots or sorting plots may be visualized.
#
Height <- 500L
HEIGHT <- paste0(as.character(Height, "px"))
# Defaults
loadFun(paste0(wd, "/Clustering/HeatMaps.RData"))
loadFun(paste0(wd, "/Dimensionality red. plots/DimRedPlots.RData"))
loadFun(paste0(wd, "/Sorting plots/quantPlots.RData"))
loadFun(paste0(wd, "/Profile plots/profilePlots.RData"))
heatMaps_ON <- exists("plotLeatMaps")&&(length(plotLeatMaps) > 0L)
dimRed_ON <- exists("dimRedPlotLy")&&(length(dimRedPlotLy) > 0L)
quant_ON <- exists("ggQuantLy")&&(length(ggQuantLy) > 0L)
profile_ON <- exists("ggProfLy")&&(length(ggProfLy) > 0L)
if (sum(c(heatMaps_ON, dimRed_ON, quant_ON, profile_ON))) {
  HEIGHT2 <- paste0(as.character(Height*2L/(quant_ON+profile_ON)), "px")
  if (dimRed_ON) {
    nmsDmRds <- names(dimRedPlotLy)[which(names(dimRedPlotLy) != "Samples PCA")]
    dfltDmRd <- c("PCA", nmsDmRds)
    dfltDmRd <- dfltDmRd[which(dfltDmRd %in% nmsDmRds)[1L]]
  }
  if (heatMaps_ON) {
    dfltHtMp <- c("Global", names(plotLeatMaps))
    dfltHtMp <- dfltHtMp[which(dfltHtMp %in% names(plotLeatMaps))[1L]]
    dfltNrmTp <- c("Norm. by row", names(plotLeatMaps[[dfltHtMp]]))
    dfltNrmTp <- dfltNrmTp[which(dfltNrmTp %in% names(plotLeatMaps[[dfltHtMp]]))[1L]]
  }
  if (quant_ON||profile_ON) {
    dfltQuant <- c("Expression", "LFQ")
    if (profile_ON) {
      allProt <- PG$Label
      if (prot.list.Cond) {
        g <- grsep(prot.list, x = PG$"Protein IDs")
        allProt <- allProt[c(g,
                             which(!1L:nrow(PG) %in% g))]
      }
      dfltProt <- allProt[1L]
      dfltQuant <- union(dfltQuant, names(ggProfLy))
      dfltQuant <- dfltQuant[which(dfltQuant %in% names(ggProfLy))[1L]]
    }
    if (quant_ON) {
      dfltQuant <- union(dfltQuant, names(ggQuantLy))
      dfltQuant <- dfltQuant[which(dfltQuant %in% names(ggQuantLy))[1L]]
      allSamples <- unique(unlist(lapply(names(ggQuantLy), \(x) { names(ggQuantLy[[x]]) })))
    }
  }
  #
  ui <- fluidPage(
    useShinyjs(),
    setBackgroundColor(color = c("#F8F8FF", "#F1F5DF"),
                       gradient = "linear", direction = "bottom"),
    titlePanel(tagList(tags$u("Data exploration"), dtstNm)),
    br(),
    em("Here you can explore the structure of the dataset."), br(),
    em("Click "), actionButton("exitBtn", "exit"), em(" to continue."), br(), br(),
    #
    if (heatMaps_ON) {
      fluidRow(column(12L,
                      h4("Clustering heatmap"),
                      selectInput("HeatMap", "Select heatmap", names(plotLeatMaps), dfltHtMp),
                      radioButtons("NormType", "Normalisation method:", c("Norm. by row", "None"), "None"),
                      withSpinner(plotlyOutput("myHeatmap_plot"))))
    },
    #
    if (dimRed_ON) {
      fluidRow(column(8L,
                      h4("Dimensionality reduction plots"),
                      selectInput("DimRedPlot", "Select method", nmsDmRds, dfltDmRd),
                      withSpinner(plotlyOutput("myPGDimRed_plot", height = HEIGHT))),
               column(4L,
                      br(),
                      br(),
                      withSpinner(plotlyOutput("mySmplsDimRed_plot", height = HEIGHT))))
    },
    #
    if (quant_ON||profile_ON) {
      h4("Protein plots")
    },
    if (quant_ON||profile_ON) {
      fluidRow(
        column(3L,
               selectInput("QuantType", "Select Quant method", QuantTypes, dfltQuant)),
        if (profile_ON) {
          column(3L,
                 pickerInput("ProfPlotProt", "Select protein(s) to display", allProt, dfltProt, TRUE,
                             pickerOptions(title = "Search me",
                                           `live-search` = TRUE,
                                           actionsBox = TRUE,
                                           deselectAllText = "Clear search",
                                           showTick = TRUE)))
        },
        if (quant_ON) {
          column(3L,
                 selectInput("Sample", "Select Sample to plot", allSamples, allSamples[1L]))
        },
      )
    },
    if (quant_ON||profile_ON) {
      fluidRow(
        if (profile_ON) {
          column(12L/(quant_ON+profile_ON),
                 withSpinner(plotlyOutput("myProfilePlot", height = HEIGHT2)))
        },
        if (quant_ON) {
          column(12L/(quant_ON+profile_ON),
                 withSpinner(plotlyOutput("SortedPGPlot", height = HEIGHT2)))
        },
      )
    },
    #
    br(), br()
  )
  server <- \(input, output, session) {
    HEATMAP <- reactiveVal()
    NORMTYPE <- reactiveVal()
    DIMRED <- reactiveVal()
    QUANTTYPE <- reactiveVal()
    PROFPROT <- reactiveVal()
    SAMPLE <- reactiveVal()
    #
    if (heatMaps_ON) {
      HEATMAP <- reactiveVal(dfltHtMp)
      NORMTYPE <- reactiveVal(dfltNrmTp)
    }
    if (dimRed_ON) { DIMRED <- reactiveVal(dfltDmRd) }
    if (quant_ON||profile_ON) { QUANTTYPE <- reactiveVal(dfltQuant) }
    if (profile_ON) { PROFPROT <- reactiveVal(dfltProt) }
    if (quant_ON) { SAMPLE <- reactiveVal(allSamples[1L]) }
    #
    observeEvent(input$HeatMap, HEATMAP(input$HeatMap))
    observeEvent(input$NormType, NORMTYPE(input$NormType))
    observeEvent(input$DimRedPlot, DIMRED(input$DimRedPlot))
    observeEvent(input$QuantType, QUANTTYPE(input$QuantType))
    observeEvent(input$ProfPlotProt, PROFPROT(input$ProfPlotProt))
    observeEvent(input$Sample, SAMPLE(input$Sample))
    #
    output$myHeatmap_plot <- renderPlotly(suppressMessages({
      req(heatMaps_ON)
      req(plotLeatMaps[[HEATMAP()]][[NORMTYPE()]]$Plot)
      return(plotLeatMaps[[HEATMAP()]][[NORMTYPE()]]$Plot)
    }))
    #
    output$mySmplsDimRed_plot <- renderPlotly(suppressMessages({
      req(dimRed_ON)
      req(dimRedPlotLy[["Samples PCA"]])
      dimRedPlotLy[["Samples PCA"]]
    }))
    output$myPGDimRed_plot <- renderPlotly(suppressMessages({
      req(dimRed_ON)
      req(dimRedPlotLy[[DIMRED()]])
      return(dimRedPlotLy[[DIMRED()]])
    }))
    #
    output$myProfilePlot <- renderPlotly(suppressMessages({
      req(profile_ON)
      req(ggProfLy[[QUANTTYPE()]])
      req(PROFPROT())
      myProfPlot <- ggProfLy[[QUANTTYPE()]]
      ggCall <- gsub(", +alpha += +[^\\)]+\\)", ")", myProfPlot$ggCall)
      plCall <- myProfPlot$plCall
      profData <- myProfPlot$data
      srchCol <- c("Protein IDs", "Proteins")[grepl("^peptides ", QUANTTYPE())+1L]
      #if (grepl("^peptides ", QUANTTYPE())) { } else { }
      toolTip <- paste0("text", as.character(1L:4L))
      yKol <- myProfPlot$data2$yCol
      colKol <- myProfPlot$data2$color
      txt4Kol <- myProfPlot$data2$labCol
      w <- grsep(gsub(" - .*", "", PROFPROT()), x = profData[[srchCol]])
      req(w)
      profData <- profData[w,]
      ttl <- myProfPlot$Ttl
      kolnm <- myProfPlot$data2$colName
      myFacets <- myProfPlot$data2$facets
      Ngl <- myProfPlot$data2$angle
      # Now that the original objects have their proper values again, we can evaluate
      eval(parse(text = ggCall))
      eval(parse(text = plCall))
      return(plotlyProfiles)
    }))
    output$SortedPGPlot <- renderPlotly(suppressMessages({
      req(quant_ON)
      req(ggQuantLy[[QUANTTYPE()]][[SAMPLE()]]$plotly)
      return(ggQuantLy[[QUANTTYPE()]][[SAMPLE()]]$plotly)
    }))
    #
    observeEvent(input$exitBtn, stopApp())
    session$onSessionEnded(\() stopApp())
  }
  eval(parse(text = runApp), envir = .GlobalEnv)
  shinyCleanup()
}
