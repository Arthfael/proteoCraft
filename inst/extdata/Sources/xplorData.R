####################
# Explore the data #
####################
#
# Opens a shiny app where clustering heatmaps, dimensionality reduction plots, protein profile plots or sorting plots may be visualized.
#
Height <- 500L
HEIGHT <- paste0(as.character(Height, "px"))
# Defaults
try({
  loadFun(heatMaps_fl)
  loadFun(dimRed_fl)
  loadFun(paste0(wd, "/Ranked abundance/quantPlots.RDS"))
  loadFun(paste0(wd, "/Profile plots/profilePlots.RDS"))
}, silent = TRUE)
heatMaps_ON <- exists("plotLeatMaps") && (length(plotLeatMaps) > 0L)
dimRed_ON <- exists("dimRedPlotLy") && (length(dimRedPlotLy) > 0L)
quant_ON <- exists("ggQuantLy") && (length(ggQuantLy) > 0L)
profile_ON <- exists("ggProfLy") && (length(ggProfLy) > 0L)
if (sum(c(heatMaps_ON, dimRed_ON, quant_ON, profile_ON))) {
  HEIGHT2 <- paste0(as.character(Height*2L/(quant_ON+profile_ON)), "px")
  if (dimRed_ON) {
    nmsDmRds <- names(dimRedPlotLy$PG)[names(dimRedPlotLy$PG) != "Samples PCA"]
    dfltDmRd <- c("PCA", nmsDmRds)
    dfltDmRd <- base::intersect(dfltDmRd, nmsDmRds)[1L]
  }
  if (heatMaps_ON) {
    dfltHtMp <- c("Global", names(plotLeatMaps))
    dfltHtMp <- base::intersect(dfltHtMp, names(plotLeatMaps))[1L]
    dfltNrmTp <- c("Norm. by row", names(plotLeatMaps[[dfltHtMp]]))
    dfltNrmTp <- base::intersect(dfltNrmTp, names(plotLeatMaps[[dfltHtMp]]))[1L]
  }
  if (quant_ON || profile_ON) {
    dfltQuant <- c("Expression", "LFQ")
    if (profile_ON) {
      allProt <- PG$Label
      if (prot.list.Cond) {
        g <- grsep(prot.list, x = PG$"Protein IDs")
        allProt <- allProt[c(g,
                             which(!1L:nrow(PG) %in% g))]
      }
      dfltProt <- allProt[1L]
      dfltQuant <- base::union(dfltQuant, names(ggProfLy))
      dfltQuant <- base::intersect(dfltQuant, names(ggProfLy))[1L]
    }
    if (quant_ON) {
      dfltQuant <- base::union(dfltQuant, names(ggQuantLy))
      dfltQuant <- base::intersect(dfltQuant, names(ggQuantLy))[1L]
      allSamples <- unique(unlist(lapply(names(ggQuantLy), \(x) { names(ggQuantLy[[x]]) })))
    }
  }
  #
  ui <- shiny::fluidPage(
    shinyjs::useShinyjs(),
    shinyWidgets::setBackgroundColor(color = c("#F8F8FF", "#F1F5DF"),
                                     gradient = "linear", direction = "bottom"),
    shiny::titlePanel(shiny::tagList(shiny::tags$u("Data exploration"), dtstNm)),
    shiny::br(),
    shiny::em("Here you can explore the structure of the dataset."), shiny::br(),
    shiny::em("Click "), shiny::actionButton("exitBtn", "exit"), shiny::em(" to continue."), shiny::br(), shiny::br(),
    #
    if (heatMaps_ON) {
      shiny::div(
        shiny::fluidRow(shiny::column(12L,
                                      shiny::h4("Clustering heatmap"))),
        shiny::fluidRow(shiny::column(2L,
                                      shiny::selectInput("HeatMap", "Select heatmap", names(plotLeatMaps), dfltHtMp)),
                        shiny::column(2L,
                                      shiny::radioButtons("NormType", "Normalisation method:", c("Norm. by row", "None"), "None"))),
        shiny::fluidRow(shiny::column(12L,
                                      shinycssloaders::withSpinner(shiny::plotlyOutput("myHeatmap_plot"))))
      )
    },
    #
    if (dimRed_ON) {
      shiny::fluidRow(shiny::column(8L,
                                    shiny::h4("Dimensionality reduction plots"),
                                    shiny::selectInput("DimRedPlot", "Select method", nmsDmRds, dfltDmRd),
                                    shinycssloaders::withSpinner(shiny::plotlyOutput("myPGDimRed_plot", height = HEIGHT))),
                      shiny::column(4L,
                                    shiny::br(),
                                    shiny::br(),
                                    shinycssloaders::withSpinner(shiny::plotlyOutput("mySmplsDimRed_plot", height = HEIGHT))))
    },
    #
    if (quant_ON || profile_ON) {
      shiny::h4("Protein plots")
    },
    if (quant_ON || profile_ON) {
      shiny::fluidRow(shiny::column(3L,
                                    shiny::selectInput("QuantType", "Select Quant method", QuantTypes, dfltQuant)),
                                    if (profile_ON) {
                                      shiny::column(3L,
                                                    shinyWidgets::pickerInput("ProfPlotProt", "Select protein(s) to display", allProt, dfltProt, TRUE,
                                                                              shinyWidgets::pickerOptions(title = "Search me",
                                                                              `live-search` = TRUE,
                                                                              actionsBox = TRUE,
                                                                              deselectAllText = "Clear search",
                                                                              showTick = TRUE)))
                                    },
                                    if (quant_ON) {
                                      shiny::column(3L,
                                                    shiny::selectInput("Sample",
                                                                       paste0("Select Sample ", c("", "Group ")[(scrptType == "withReps")+1L], "to plot"),
                                                                       allSamples, allSamples[1L]))
                                    },
      )
    },
    if (quant_ON || profile_ON) {
      shiny::fluidRow(
        if (profile_ON) {
          shiny::column(12L/(quant_ON + profile_ON),
                        shinycssloaders::withSpinner(shiny::plotlyOutput("myProfilePlot", height = HEIGHT2)))
        },
        if (quant_ON) {
          shiny::column(12L/(quant_ON + profile_ON),
                        shinycssloaders::withSpinner(shiny::plotlyOutput("SortedPGPlot", height = HEIGHT2)))
        },
      )
    },
    #
    shiny::br(), shiny::br()
  )
  server <- \(input, output, session) {
    HEATMAP <- shiny::reactiveVal()
    NORMTYPE <- shiny::reactiveVal()
    DIMRED <- shiny::reactiveVal()
    QUANTTYPE <- shiny::reactiveVal()
    PROFPROT <- shiny::reactiveVal()
    SAMPLE <- shiny::reactiveVal()
    #
    if (heatMaps_ON) {
      HEATMAP <- shiny::reactiveVal(dfltHtMp)
      NORMTYPE <- shiny::reactiveVal(dfltNrmTp)
    }
    if (dimRed_ON) { DIMRED <- shiny::reactiveVal(dfltDmRd) }
    if (quant_ON || profile_ON) { QUANTTYPE <- shiny::reactiveVal(dfltQuant) }
    if (profile_ON) { PROFPROT <- shiny::reactiveVal(dfltProt) }
    if (quant_ON) { SAMPLE <- shiny::reactiveVal(allSamples[1L]) }
    #
    shiny::observeEvent(input$HeatMap, HEATMAP(input$HeatMap))
    shiny::observeEvent(input$NormType, NORMTYPE(input$NormType))
    shiny::observeEvent(input$DimRedPlot, DIMRED(input$DimRedPlot))
    shiny::observeEvent(input$QuantType, QUANTTYPE(input$QuantType))
    shiny::observeEvent(input$ProfPlotProt, PROFPROT(input$ProfPlotProt))
    shiny::observeEvent(input$Sample, SAMPLE(input$Sample))
    #
    output$myHeatmap_plot <- shiny::renderPlotly(suppressMessages({
      shiny::req(heatMaps_ON)
      shiny::req(plotLeatMaps[[HEATMAP()]][[NORMTYPE()]]$Plot)
      return(plotLeatMaps[[HEATMAP()]][[NORMTYPE()]]$Plot)
    }))
    #
    output$mySmplsDimRed_plot <- shiny::renderPlotly(suppressMessages({
      req(dimRed_ON)
      req(dimRedPlotLy$PG[["Samples PCA"]])
      dimRedPlotLy$PG[["Samples PCA"]]
    }))
    output$myPGDimRed_plot <- shiny::renderPlotly(suppressMessages({
      shiny::req(dimRed_ON)
      shiny::req(dimRedPlotLy$PG[[DIMRED()]])
      return(dimRedPlotLy$PG[[DIMRED()]])
    }))
    #
    output$myProfilePlot <- shiny::renderPlotly(suppressMessages({
      shiny::req(profile_ON)
      shiny::req(ggProfLy[[QUANTTYPE()]])
      shiny::req(PROFPROT())
      myProfPlot <- ggProfLy[[QUANTTYPE()]]
      ggCall <- gsub(", +alpha += +[^)]+\\)", ")", myProfPlot$ggCall)
      plCall <- myProfPlot$plCall
      profData <- myProfPlot$data
      srchCol <- c("Protein IDs", "Proteins")[grepl("^peptides ", QUANTTYPE())+1L]
      #if (grepl("^peptides ", QUANTTYPE())) { } else { }
      toolTip <- paste0("text", as.character(1L:4L))
      yKol <- myProfPlot$data2$yCol
      colKol <- myProfPlot$data2$color
      txt4Kol <- myProfPlot$data2$labCol
      w <- grsep(sub(" - .*", "", PROFPROT()), x = profData[[srchCol]])
      shiny::req(w)
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
    output$SortedPGPlot <- shiny::renderPlotly(suppressMessages({
      shiny::req(quant_ON)
      shiny::req(ggQuantLy[[QUANTTYPE()]][[SAMPLE()]]$plotly)
      return(ggQuantLy[[QUANTTYPE()]][[SAMPLE()]]$plotly)
    }))
    #
    shiny::observeEvent(input$exitBtn, shiny::stopApp())
    session$onSessionEnded(\() shiny::stopApp())
  }
  g <- shiny:::.globals
  g$appState <- NULL
  eval(parse(text = run_App), envir = .GlobalEnv)
  shinyCleanup()
}
try({
  rm(plotLeatMaps,
     dimRedPlotLy,
     ggQuantLy,
     ggProfLy)
}, silent = TRUE)
