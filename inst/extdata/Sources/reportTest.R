library(shiny)
library(shinyjs)
library(bslib)
library(htmltools)
library(DT)
library(plotly)

htmlRprtFl <- paste0(wd, "/Report_", dtstNm, ".html")
htmlSCRprtFl <- paste0(wd, "/Report_", dtstNm, "_SC.html")

# Reload processed data from the report
allNms <- openxlsx2::wb_get_sheet_names(openxlsx2::wb_load(repFl))
nms <- setdiff(allNms, c("Description", "Quality control"))
xlDat <- setNames(lapply(nms, \(nm) {
  dat <- openxlsx2::read_xlsx(repFl, match(nm, allNms), 2L)
  if ("Potential contaminant" %in% colnames(dat)) {
    w <- which(is.na(dat$"Potential contaminant"))
    if (length(w)) { dat$"Potential contaminant"[w] <- "" }
  }
  return(dat)
}), nms)

#
plotHght <- paste0(round(screenRes$height*0.75), "px")

# UI functions
make_prot_tab <- \(prot,
                   prots = allProt,
                   dflt = dfltProt,
                   shiny = TRUE) {
  # - show:
  # Proteins tab
  #################################
  #     dropdown for protein      #
  #################################
  # ->
  #################################
  #      comment for protein      #
  #################################
  ################## ##############
  #samples dropdown# #            #
  ################## #            #
  ################## #   Ratios   #
  #                # #    plot    #
  #    Coverage    # #            #
  #                # #            #
  ################## ##############
  # Peptides table
  myTags <- if (shiny) {
    tags$div(
      if (prot.list.Cond && (length(prots) > 1L)) {
        tags$div(selectInput("myProtein", "Select protein", allProt, dflt),
                 tags$br())
      },
      uiOutput("protComment"),
      tags$br(),
      fluidRow(column(6L,
                      plotlyOutput("ratioPlot", height = plotHght)),
               column(6L,
                      if (length(Exp) > 1L) {
                        selectInput("mySample", "", Exp, Exp[1L]) 
                      },
                      plotlyOutput("coverPlot", height = plotHght))),
      br(),
      tags$hr(style = "border-color: black;"),
      uiOutput("protPep"),
      br())
  } else {
    ## Coverage plots ###################################################
    cov_plots <- list()
    dfltSmpl <- Exp[1L]
    for (prot in prots) {
      for (exp in Exp) {
        visible <- (prot == dflt) && (exp == dfltSmpl)
        cov_plots <- append(cov_plots,
                            tags$div(id = paste0("cov_", prot, "_", exp),
                                     style = paste("width:100%;",
                                                   if (visible) { "display:block;" } else { "display:none;" } ),
                                     covPlots[[prot]]$logInt[[exp]]))
      }
    }
    ## Ratio plots ######################################################
    ratioPlotsUI <- list()
    if (exists("ratioPlots")) {
      for (prot in intersect(prots, names(ratioPlots))) {
        visible <- (prot == dflt)
        ratioPlotsUI <- append(ratioPlotsUI,
                               tags$div(id = paste0("rat_", prot),
                                        style = if (visible) { "display:block;" } else { "display:none;" },
                                        ratioPlots[[prot]]))
      }
    }
    ## UI ###############################################################
    tagList(
      if (length(prots) > 1L) {
        fluidRow(column(12L,
                        tags$label("Protein"),
                        tags$select(id = "myProtein",
                                    lapply(prots, \(pr) {
                                      tags$option(value = pr,
                                                  selected = pr == dflt[1L],
                                                  pr)
                                    })),
                        tags$br()))
      },
      make_comment_ui(dflt, shiny),
      tags$br(),
      fluidRow(column(6L,
                      tags$label("Sample"),
                      if (length(Exp) > 1L) {
                        tags$select(id = "mySample",
                                    lapply(Exp, \(x) {
                                      tags$option(value = x,
                                                  selected = x == Exp[1L],
                                                  x)
                                    }))
                      },
                      tags$br(),
                      cov_plots),
               if (length(ratioPlotsUI)) {
                 column(6L, ratioPlotsUI)
               },
      ),
      tags$script(HTML("function updateProteinPlots() {
  var prot = document.getElementById('myProtein').value;
  var sample = document.getElementById('mySample').value;
  document.querySelectorAll('[id^=\"cov_\"]')
    .forEach(function(el) {
      el.style.display = 'none';
    });
  var cov = document.getElementById('cov_' + prot + '_' + sample);
  if (cov)
    cov.style.display = 'block';
  document.querySelectorAll('[id^=\"rat_\"]')
    .forEach(function(el) {
      el.style.display = 'none';
    });
  var rat = document.getElementById('rat_' + prot);
  if (rat)
    rat.style.display = 'block';
  window.dispatchEvent(
    new Event('resize')
  );
}
document.getElementById('mySample')
  .addEventListener('change',
                    updateProteinPlots);
var prot = document.getElementById('myProtein');
if (prot)
  prot.addEventListener('change',
                        updateProteinPlots);")),
      tags$br(),
      tags$hr(style = "border-color: black;"))
  }
return(myTags)
}
make_comment_ui <- \(id,
                     shiny = TRUE,
                     values = allComments,
                     toggle = FALSE) {
  if (shiny) {
    textAreaInput(inputId = paste0("comment_", id),
                  label = NULL,
                  value = values[id],
                  width = "100%",
                  height = "150px")
  } else {
    style <- "white-space: pre-wrap; padding: 10px;"
    if (toggle) {
      m <- match(id, names(values))
      if (m > 1L) { style <- "display:none;"}
    }
    tags$div(class = "comment-box",
             `data-id` = paste0("comment_", id),
             style = style,
             values[id])
  }
}
make_bar <- \(x) {
  sprintf(
    '<div style="position:relative; width:100%%; background:#eee; height:16px; border-radius:4px;">
        <div style="width:%s%%; background:#4CAF50; height:100%%; border-radius:4px;"></div>
        <div style="position:absolute; top:0; left:50%%; transform:translateX(-50%%);
                    font-size:11px; line-height:16px; color:black;">
            %.1f%%
        </div>
     </div>',
    x, x)
}
make_tbl_ui <- \(exp = Exp, #exp <- Exp[1L] #exp <- Exp[2L]
                 tab = "Protein groups", # can also be "All peptidoforms"; we will eventually add "`PTM`-modified", where `PTM` can be any PTM of interest
                 filt = NULL, #filt = allProt[1L]
                 dat = xlDat) {
  df <- dat[[tab]]
  smplCols_lst <- setNames(lapply(exp, \(xp) {
    grep(topattern(paste0(" ", xp), FALSE, TRUE), colnames(df), value = TRUE)
  }), exp)
  smplCols <- setNames(unlist(smplCols_lst), NULL)
  coreCols <- c("PEP", "Potential contaminant")
  if (tab %in% c("Protein groups", "All peptidoforms")) {
    if (tab == "Protein groups") {
      filtCol <- "Protein IDs"
      coreCols <- union(c("Leading protein IDs", filtCol, "Common Names", "Genes", "Mol. weight [kDa]"), coreCols)
      intRoot <- "expr"
    }
    if (tab == "All peptidoforms") {
      filtCol <- "Proteins"
      coreCols <- union(c("Modified sequence_verbose", #"Sequence",
                          filtCol), coreCols)
      intRoot <- "int"
    }
  } else {
    stop("TO DO!")
    filtCol <- "Proteins"
    coreCols <- union(c("Modified sequence_verbose", #"Sequence",
                        filtCol), coreCols) # Check before use...
    intRoot <- "int" # Presumably...
  }
  #
  xprCols <- grep(paste0("log10\\(([^\\)]+ )?", intRoot, "\\.\\) "), smplCols, value = TRUE)
  fullIntRoot <- rev(paste0(vapply(strsplit(xprCols, "\\)"), `[[`, "", 1L), ") "))[1L]
  xprCols <- paste0(fullIntRoot, exp)
  repXprCols <- if (length(exp) == 1L) { fullIntRoot } else { xprCols }
  repXprCols <- sub("\\([Nn]or\\. ", "(norm. ", repXprCols)
  #
  ratCols <- grep(paste0("log2\\(([^\\)]+ )?rat\\.\\) "), smplCols, value = TRUE)
  useRat <- length(ratCols) > 0L
  if (useRat) {
    fullRatRoot <- rev(paste0(vapply(strsplit(ratCols, "\\)"), `[[`, "", 1L), ") "))[1L]
    ratCols <- paste0(fullRatRoot, exp)
    ratCols <- intersect(ratCols, colnames(df))
    repRatCols <- ratCols
    repRatCols <- sub("\\([Nn]or\\. ", "(norm. ", repRatCols)
  }
  #
  colNms <- c(coreCols, xprCols)
  repColNms <- c(sub("_verbose$", "", coreCols), repXprCols)
  if (useRat) {
    colNms <- c(colNms, ratCols)
    repColNms <- c(repColNms, repRatCols)
  }
  if (tab == "Protein groups") {
    k <- unlist(lapply(c("Pep. count ", "PSMs count "), \(k) { paste0(k, exp) }))
    colNms <- c(colNms, k)
    repColNms <- if (length(exp) == 1L) { c(repColNms, "Pep. count", "PSMs count") } else { c(repColNms, k) }
  }
  df <- df[, colNms]
  flt <- if (is.null(filt)) {
    1L:nrow(df)
  } else {
    grsep(db$`Protein ID`[match(filt, db$`Common Name`)], x = df[[filtCol]])
  }
  df <- df[flt,]
  colnames(df) <- colNms <- repColNms
  if ("Modified sequence" %in% colNms) {
    df$"Modified sequence" <- gsub("^_|_$", "",  df$"Modified sequence")
  }
  xprCols <- repXprCols
  if (useRat) { ratCols <- repRatCols }
  if (tab == "Protein groups") {
    covCols <- paste0("Cov. ", exp)
    repCovCols <- if (length(exp) == 1L) { "Coverage" } else { covCols }
    if (covCols %in% colnames(dat[[tab]])) {
      df[, repCovCols] <- dat[[tab]][flt, covCols]
    } else {
      if (("Coverage" %in% names(dat)) && (covCols %in% colnames(dat$Coverage))) {
        m <- match(df$"Protein IDs"[flt], dat$Coverage$`Protein IDs`)
        df[, repCovCols] <- dat$Coverage[m, covCols]
      }
    }
    covCols <- repCovCols
  }
  xprRng <- range(df[, xprCols], na.rm = TRUE)
  #
  # Make sure this re-ordering is done after any other data is added from dat to df!
  orderVect <- df[, xprCols]
  if (length(exp) > 1L) { orderVect <- apply(orderVect, 1L, \(x) { mean(x[which(is.finite(x))]) }) }
  df <- df[order(orderVect, decreasing = TRUE),]
  #
  quantCols <- xprCols
  if (useRat) {
    ratRng <- range(df[[ratCols]], na.rm = TRUE)
    quantCols <- c(quantCols, ratCols)
  }
  col2 <- setdiff(colnames(df), c("PEP", quantCols))
  df[, col2] <- sapply(col2, \(k) {
    if ((tab == "Protein groups") && (k %in% covCols)) {
      make_bar(df[[k]])
    } else {
      sprintf("<div class=\"cell-wrap\">%s</div>", df[[k]])
    }
  })
  wTest1 <- setNames(vapply(colnames(df), \(k) { #k <- colnames(df)[1L]
    tmp <- as.character(df[[k]])
    x <- min(c(250L, max(nchar(c(k, tmp))+3L, na.rm = TRUE)*10L))
    if (is.na(x)) { x <- 50L }
    return(as.integer(x))
  }, 1L), colnames(df))
  wTest2 <- sum(wTest1) + 15L + ncol(df)*5L
  wTest1 <- paste0(as.character(wTest1), "px")
  wTest1 <- aggregate((1L:length(wTest1))-1L, list(wTest1), c)
  wTest1 <- apply(wTest1, 1L, \(x) {
    x2 <- as.integer(x[[2L]])
    list(width = x[[1L]],
         targets = x2,
         names = colnames(df)[x2+1L])
  })
  df <- DT::datatable(df,
                      rownames = FALSE,
                      class = "compact",
                      escape = FALSE,
                      options = list(#scrollX = TRUE,
                        scrollY = "500px",
                        autoWidth = FALSE,
                        columnDefs = wTest1))
  df <- DT::formatRound(df, c("PEP", quantCols), digits = 5L)
  df <- DT::formatStyle(df, "PEP",
                        backgroundColor = DT::styleInterval(10^-seq(10, 0, length.out = 99L),
                                                            colorRampPalette(rev(ColScaleList$PEP))(100L)))
  df <- DT::formatStyle(df, xprCols,
                        backgroundColor = DT::styleInterval(seq(xprRng[1L], xprRng[2L], length.out = 99L),
                                                            colorRampPalette(ColScaleList$`Individual Expr`)(100L)))
  if (useRat) {
    df <- DT::formatStyle(df, ratCols,
                          backgroundColor = DT::styleInterval(seq(ratRng[1L], ratRng[2L], length.out = 99L),
                                                              colorRampPalette(ColScaleList$`Individual Ratios`)(100L)))
  }
  return(df)
}
# make_plot_div <- function(widget, id, visible = TRUE) {
#   tags$div(id = id,
#            style = sprintf(
#              "display:%s; width:100%%;",
#              if (visible) "block" else "none"
#            ),
#            widget)
# }
make_smpl_tab <- \(exp,
                   shiny = TRUE,
                   dflt = dfltQuant) {
  if (shiny) {
    tagList(make_comment_ui(exp, shiny),
            selectInput(paste0("quant_", exp), "", c("LFQ", "Coverage"), dflt[exp]),
            plotlyOutput(paste0("quantLy_", exp), height = plotHght),
            br(),
            br(),
            tags$hr(style = "border-color: black;"),
            make_tbl_ui(exp))
  } else {
    js <- sprintf(
      "document.getElementById('%s').addEventListener(
  'change',
  function() {
    var selected = document.getElementById('%s').value;
    document.querySelectorAll('[id^=\"%s\"]').forEach(function(div) {
      div.style.display = 'none';
    });
    document.getElementById('%s' + selected).style.display = 'block';
    window.dispatchEvent(new Event('resize'));
  }
);",
      paste0("quant_", exp),
      paste0("quant_", exp),
      paste0("Quant_", exp, "_"),
      paste0("Quant_", exp, "_")
    )
    tagList(make_comment_ui(exp, shiny),
            tags$div(id = paste0("Quant_", exp, "_", 1L),
                     style = if (dflt[exp] == "LFQ") { "display:block;" } else { "display:none;" },
                     ggQuantLy$LFQ[[exp]]$plotly),
            tags$div(id = paste0("Quant_", exp, "_", 2L),
                     style = if (dflt[exp] == "Coverage") { "display:block;" } else { "display:none;" },
                     ggQuantLy$Coverage[[exp]]$plotly),
            br(),
            tags$select(id = paste0("quant_", exp),
                        lapply(1L:2L, \(i) {
                          tags$option(value = i,
                                      c("LFQ", "Coverage")[i])
                        })),
            br(),
            tags$hr(style = "border-color: black;"),
            tags$script(HTML(js)),
            make_tbl_ui(exp))
  }
}
make_strt_tab <- \(type = "Z-scored",
                   shiny = TRUE) {
  if (shiny) {
    tagList(make_comment_ui("Dataset overview", shiny),
            fluidRow(column(6L,
                            plotlyOutput("heatMap", height = plotHght)),
                     column(6L,
                            plotlyOutput("PCA", height = plotHght))),
            br())
  } else {
    tagList(make_comment_ui("Dataset overview", shiny),
            tags$div(style = "
        display: grid;
        grid-template-columns: repeat(2, 1fr);
        gap: 12px;
      ",
                     list(plotLeatMaps$Global[[type]]$Render,
                          height = plotHght)))
  }
  
}
make_QC_ui <- \(shiny = TRUE,
                plotsList = QC_plotLys) {
  if (shiny) {
    tabPanel("QC",
             tagList(selectInput("QC1", "", names(plotsList), names(plotsList)[1L]),
                     fluidRow(column(8L,
                                     plotlyOutput("QCplotLy", height = plotHght)),
                              column(4L,
                                     uiOutput("QCtxt"))),
                     br()))
  } else {
    tabPanel("QC",
             tagList(tags$select(id = "myQC",
                                 lapply(seq_along(plotsList), \(i) {
                                   tags$option(value = i,
                                               names(plotsList)[i])
                                 })),
                     lapply(seq_along(plotsList), \(i) {
                       nm <- names(plotsList)[i]
                       fluidRow(column(8L,
                                       tags$div(id = paste0("QC_", i),
                                                style = if (i == 1L) { "display:block;" } else { "display:none;" },
                                                plotsList[[i]])),
                                column(4L,
                                       make_comment_ui(nm, FALSE, toggle = TRUE)))
                     }),
                     tags$script(HTML("document.getElementById('myQC').addEventListener(
    'change',
    function() {
        var selected = document.getElementById('myQC').value;
        document.querySelectorAll('[id^=\"QC_\"]').forEach(function(div) { div.style.display = 'none'; });
        document.getElementById('QC_' + selected).style.display = 'block';
        document.querySelectorAll('[id^=\"comment_\"]').forEach(function(div) { div.style.display = 'none'; });
        var div = document.getElementById('comment_' + selected);
        div.style.display = 'block';
        div.style.whiteSpace = 'pre-wrap';
        div.style.padding = '10px';
        window.dispatchEvent(new Event('resize'));
    }
);")),
                     br()))
  }
}
make_ui <- \(shiny = TRUE) {
  tabs <- lapply(myTabs, function(x) {
    if (x == "Dataset overview") {
      return(tabPanel(x,
                      make_strt_tab(shiny = shiny)))
    }
    if (x %in% Exp) {
      return(tabPanel(paste("sample = ", x),
                      make_smpl_tab(x,
                                    shiny = shiny)))
    }
    if (x == "Proteins of interest") {
      return(tabPanel(x,
                      tagList(
                        if (ratiosTest) {
                          make_prot_tab(dfltProt,
                                        shiny = shiny)
                        },
                        tags$br()
                      )))
    }
    if (x == "QC") {
      return(make_QC_ui(shiny = shiny))
    }
  })
  return(bslib::navset_tab(!!!tabs))
}

# If necessary reload plots data
loadFun(paste0(wd, "/Clustering/HeatMaps.RData"))
loadFun(paste0(wd, "/Sorting plots/quantPlots.RData"))
flPCA <- paste0(wd, "/Dimensionality red. plots/DimRedPlots.RData")
tstPCA <- file.exists(flPCA)
if (tstPCA) {
  loadFun(flPCA)
  tstPCA <- exists("dimRedPlotLy") && ("PCA" %in% names(dimRedPlotLy))
}

# Fix to plotly autoscaling
if (exists("plotLeatMaps")) {
  for (x in names(plotLeatMaps)) { #x <- names(plotLeatMaps)[1L]
    for (y in names(plotLeatMaps[[x]])) { #y <- names(plotLeatMaps[[x]])[1L]
      plotLeatMaps[[x]][[y]]$Render$x$layout$xaxis$autorange <- TRUE
      plotLeatMaps[[x]][[y]]$Render$x$layout$yaxis$autorange <- TRUE
    }
  }
}
if (exists("ggQuantLy")) {
  for (x in names(ggQuantLy)) { #x <- names(ggQuantLy)[1L]
    for (y in names(ggQuantLy[[x]])) { #y <- names(ggQuantLy[[x]])[1L]
      ggQuantLy[[x]][[y]]$plotly$x$layout$xaxis$autorange <- TRUE
      ggQuantLy[[x]][[y]]$plotly$x$layout$yaxis$autorange <- TRUE
    }
  }
}
if (tstPCA) {
  for (x in names(dimRedPlotLy)) { #x <- names(dimRedPlotLy)[1L]
    dimRedPlotLy[[x]]$x$layout$xaxis$autorange <- TRUE
    dimRedPlotLy[[x]]$x$layout$yaxis$autorange <- TRUE
  }
}
if (exists("covPlots")) {
  for (x in names(covPlots)) { #x <- names(covPlots)[1L]
    for (y in names(covPlots[[x]])) { #y <- names(covPlots[[x]])[1L]
      for (z in names(covPlots[[x]][[y]])) { #z <- names(covPlots[[x]][[y]])[1L]
        covPlots[[x]][[y]][[z]]$x$layout$xaxis$autorange <- TRUE
        covPlots[[x]][[y]][[z]]$x$layout$yaxis$autorange <- TRUE
      }
    }
  }
}

# Plot HTML paths
#myPlots <- list.files(paste0(wd, "/Ranked abundance/LFQ"), "\\.html$", full.names = TRUE)
#names(myPlots) <- gsub(".* - |\\.html$", "", myPlots)
#nPl <- length(myPlots)
myTabs <- nms <- c("Dataset overview", Exp)
if (prot.list.Cond) {
  allProt <- names(covPlots)
  dfltProt <- allProt[1L]
  myTabs <- union(myTabs, "Proteins of interest")
  nms <- union(nms, allProt)
}
myTabs <- union(myTabs, "QC")
nms <- union(nms, c("QC", names(QC_plotLys)))
dfltComment <- paste0(nrow(PG), " protein groups were identified from ", nrow(ev), " PSMs", " corresponding to ", nrow(pep),
                      " distinct peptidoforms. ...")
if ((!exists("allComments")) || (!is.character(allComments))) {
  allComments <- setNames(vapply(nms, \(nm) {
    if (nm == "Dataset overview") {
      dfltComment
    } else { "" }
  }, ""), nms)
}
if (sum(!nms %in% names(allComments))) {
  nms_ <- setdiff(nms, names(allComments))
  allComments[nms_] <- ""
  if ("Dataset overview" %in% nms_) { allComments$"Dataset overview" <- dfltComment }
}
allComments <- allComments[nms]
ratiosTest <- exists("ratioPlots") && (length(ratioPlots) > 0L)

#
quantTst <- setNames(lapply(names(ggQuantLy), \(tp) { names(ggQuantLy[[tp]]) }), names(ggQuantLy))
quantTst <- listMelt(quantTst, ColNames = c("Sample", "Type"))
quantTst <- aggregate(quantTst$Type, list(quantTst$Sample), list)
colnames(quantTst) <- c("Sample", "Types")
dfltQuant <- quantTst[, "Sample", drop = FALSE]
dfltQuant$Type <- vapply(quantTst$Types, \(x) { x[[1L]] }, "")
quantTst <- setNames(quantTst$Types, quantTst$Sample)
dfltQuant <- setNames(dfltQuant$Type, dfltQuant$Sample)
#
appPage <- 1L
appNm <- "Edit report"
ui <- fluidPage(useShinyjs(),
                extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
                tags$head(tags$style(HTML("table.dataTable td {
    white-space: normal !important;
    vertical-align: top !important;
}")),
                          tags$style(HTML(".cell-wrap {
    max-height: 80px !important;
    overflow: hidden !important;
    white-space: normal;
    vertical-align: top !important;
    overflow-wrap: break-word !important;
    word-break: break-word !important;
}"))),
                titlePanel(tag("u", appNm),
                           appNm),
                br(),
                fluidRow(column(4L,
                                h2(dtstNm)),
                         column(8L,
                                actionBttn("xprtBtn", " export final html report", icon = icon("file-export"), color = "success", style = "pill"),
                                uiOutput("xprtMsg"))),
                br(),
                uiOutput("myUI"),
                br(),
                br())
server <- \(input, output, session) {
  # if (prot.list.Cond) {
  #   PROT <- reactiveVal(dfltProt)
  # }
  QUANT <- reactiveVal(dfltQuant)
  # Render UI
  output$xprtMsg <- renderUI(HTML(""))
  output$myUI <- renderUI({ make_ui() })
  output$heatMap <- renderPlotly(plotLeatMaps$Global$`Norm. by row`$Render)
  if (tstPCA) {
    output$PCA <- renderPlotly(dimRedPlotLy$PCA)
  }
  #
  for (exp in Exp) {
    idQ <- paste0("quant_", exp)
    idQLy <- paste0("quantLy_", exp)
    output[[idQLy]] <- renderPlotly(ggQuantLy[[input[[idQ]]]][[exp]]$plotly)
  }
  #
  # Event observers
  #  - Comments
  sapply(names(allComments), \(nm) {
    observeEvent(input[[paste0("comment_", nm)]], {
      allComments[[nm]] <- input[[paste0("comment_", nm)]]
      allComments <<- allComments
    })
  })
  #  - Quant method
  sapply(Exp, \(exp) {
    idQ <- paste0("quant_", exp)
    idQLy <- paste0("quantLy_", exp)
    observeEvent(input[[idQ]], {
      dfltQuant <- QUANT()
      dfltQuant[exp] <- input[[idQ]]
      QUANT(dfltQuant)
      assign("dfltQuant", dfltQuant, envir = .GlobalEnv)
      # Update plot
      output[[idQLy]] <- renderPlotly(ggQuantLy[[input[[idQ]]]][[exp]]$plotly)
    })
  })
  #  - Proteins tab
  if (prot.list.Cond) {
    if (length(allProt) > 1L) {
      observeEvent(input$myProtein, {
        output$ratioPlot <- renderPlotly(ratioPlots[[input$myProtein]])
        output$coverPlot <- if (length(Exp) > 1L) {
          renderPlotly(covPlots[[input$myProtein]]$logInt[[input$mySample]])
        } else {
          renderPlotly(covPlots[[input$myProtein]]$logInt[[Exp]])
        }
        output$protComment <- renderUI(make_comment_ui(input$myProtein))
        output$protPep <- renderUI(make_tbl_ui(tab = "All peptidoforms",
                                               filt = input$myProtein))
      })
    } else {
      output$ratioPlot <- renderPlotly(ratioPlots[[allProt]])
      output$coverPlot <- if (length(Exp) > 1L) {
        renderPlotly(covPlots[[allProt]]$logInt[[input$mySample]])
      } else {
        renderPlotly(covPlots[[allProt]]$logInt[[Exp]])
      }
      output$protComment <- renderUI(make_comment_ui(allProt))
      output$protPep <- renderUI(make_tbl_ui(tab = "All peptidoforms",
                                             filt = allProt))
    }
  }
  if (length(Exp) > 1L) {
    observeEvent(input$mySample, {
      output$coverPlot <- if (prot.list.Cond) {
        renderPlotly(covPlots[[input$myProtein]]$logInt[[input$mySample]])
      } else {
        renderPlotly(covPlots[[allProt]]$logInt[[input$mySample]])
      }
    })
  }
  #  - QC tab
  observeEvent(input$QC1, {
    output$QCplotLy <- renderPlotly(QC_plotLys[[input$QC1]])
    output$QCtxt <- renderUI(make_comment_ui(input$QC1))
  })
  #  - Render final report
  observeEvent(input$xprtBtn, {
    output$xprtMsg <- renderUI(em("Exporting .html report, this will take a few seconds...",
                                  style = "color:red",
                                  .noWS = "outside"))
    # 1. Rebuild the SAME UI we use in the app
    page <- # tagList(div(h3(paste0(dtstNm, " - report"), br(),
      #                "Analysis run by: ", em(WhoAmI), br(),
      #                "Date: ", em(Sys.Date()), br(),
      #                "Package: ", em("proteoCraft v", package.version("proteoCraft")), br())),
      make_ui(FALSE)#)
    # 2. Wrap as browsable HTML
    page <- htmltools::browsable(page)
    # 3. Save to disk
    htmltools::save_html(page, htmlRprtFl)
    assign("appRunTst", TRUE, envir = .GlobalEnv)
    stopApp()
  })
  session$onSessionEnded(\() { stopApp() })
}
runKount <- 0L
if (exists("appRunTst")) { rm(appRunTst) }
while ((!runKount) || (!exists("appRunTst")) || (!file.exists(htmlRprtFl))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  shinyCleanup()
  runKount <- runKount + 1L
  appRunTst <- TRUE
}

# We now have our html... but it depends on local libraries... we want those embedded in it!
h2 <- h1 <- readLines(htmlRprtFl)
rg1 <- grep("</?head>", h1) + c(1L, -1L)
rg1 <- rg1[1L]:rg1[2L]
hd1 <- h1[rg1]
hd1 <- data.frame(original = hd1)
hd1$new <- hd1$original
g <- grep("^ *<((style)|(script)|(link))( *[^>]+)?>", hd1$original)
hd1$original[g]
require(base64enc)
read_file <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}
file_to_data_uri <- function(path) {
  ext <- tools::file_ext(path)
  mime <- switch(tolower(ext),
                 "woff2" = "font/woff2",
                 "woff"  = "font/woff",
                 "ttf"   = "font/ttf",
                 "png"   = "image/png",
                 "jpg"   = "image/jpeg",
                 "jpeg"  = "image/jpeg",
                 "svg"   = "image/svg+xml",
                 "gif"   = "image/gif",
                 "application/octet-stream")
  paste0("data:",
         mime,
         ";base64,",
         base64enc::base64encode(path))
}
# - embed scripts
read_asset <- function(path) {
  readChar(path, # Do not use readLines, which isn't binary-safe!
           nchars = file.info(path)$size,
           useBytes = TRUE)
}
inline_script <- function(path) {
  txt <- paste(read_asset(path), collapse = "")
  txt <- gsub("</script",
              "<\\/script",
              txt,
              ignore.case = TRUE)
  paste0("<script>\n",
         txt,
         "\n</script>")
}
gs <- grep("^ *<script src=\"", hd1$original)
hd1$new[gs] <- vapply(sub("\".*", "", sub("^ *<script src=\"", paste0(wd, "/"), hd1$original[gs])), inline_script, "")
# - embed css
inline_css <- function(path) {
  paste0(  "<style>\n",
           read_asset(path),
           "\n</style>")
}
gc <- grepl("^ *<link href=\"", hd1$original)
hd1$new[gc] <- vapply(sub("\".*", "", sub("^ *<link href=\"", paste0(wd, "/"), hd1$original[gc])), inline_css, "")
h2[rg1] <- hd1$new
write(h2, htmlSCRprtFl)



# TO DO
# - remove local /lib folder
# - overwrite report instead of keeping two versions
