####################################################
# Optional: choose whether to remove any outliers  #
####################################################
#View(Exp.map[, c("Ref.Sample.Aggregate", "Use")])
appNm <- paste0(dtstNm, " - Check for outliers")
#Exp.map$Use <- "TRUE"
fctrs <- Factors
Include <- Exp.map[, c(fctrs, "Use")]
#colnames(Include) <- gsub("\\.", "_", colnames(Include))
#fctrs <- gsub("\\.", "_", fctrs)
Include$Use <- shinyCheckInput(Exp.map$Use)
# Table width
wTest0 <- setNames(sapply(colnames(Exp.map), function(k) { #k <- colnames(Exp.map)[1]
  tmp <- Exp.map[[k]]
  if ("logical" %in% class(tmp)) { tmp <- as.integer(tmp) }
  tmp <- as.character(tmp)
  x <- max(nchar(c(k, tmp)) + 3, na.rm = TRUE)
  x <- x*10
  return(x)
}), #gsub("\\.", "_",
colnames(Exp.map))#)
wTest1 <- vapply(colnames(Include), function(k) { #k <- colnames(Include)[1]
  if (k %in% names(wTest0)) { x <- wTest0[k] } else { x <- 30 }
  return(x)
}, 1)
wTest1 <- round(wTest1*screenRes$width*0.45/sum(wTest1))
#sum(wTest1)
wTest1 <- paste0(as.character(wTest1), "px")
wTest1 <- aggregate((1:length(wTest1))-1, list(wTest1), c)
wTest1 <- apply(wTest1, 1, function(x) {
  x1 <- x[[1]]
  x2 <- as.integer(x[[2]])
  list(width = x1,
       targets = x2,
       names = colnames(Include)[x2+1])
})
#
if (exists("appRunTest")) { rm(appRunTest) }
dotLab <- RSA$names
colLab <- VPAL$names
if (length(Exp) == 1) {
  dotLab <- dotLab[which(dotLab != "Experiment")]
  colLab <- colLab[which(colLab != "Experiment")]
}
dotLab <- paste(dotLab, collapse = " ")
colLab <- paste(colLab, collapse = " ")
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#E5EDE1"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", "Check for outliers"),
             appNm),
  br(),
  h5(tags$div(
    "Check visually how samples are regrouped within expected sample groupings,", tags$br(),
    "then decide whether to remove any outliers (by unticking them in the \"Use\" column),", tags$br(),
    "finally click \"Save\".")),
  br(),
  sidebarLayout(
    sidebarPanel(
      actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
      span(uiOutput("Msg"), style = "color:red"),
      br(),
      DTOutput("Include"),
      br(),
      width = 6
    ),
    mainPanel(
      h4(em("Plot mappings:")),
      h5(em(paste0("Dot labels = ", dotLab))),
      h5(em(paste0("Colour = ", outlierAnnot_color))),
      h5(em(paste0("Colour labels = ", colLab))),
      h5(em(paste0("Shape = ", outlierAnnot_shape))),
      plotlyOutput("PCA", height = paste0(screenRes$width*0.4, "px")),
      br(),
      width = 6
    )
  ))
server <- function(input, output, session) {
  tstGrps <- aggregate(Exp.map$Use, list(Exp.map[[VPAL$column]]), sum)
  if (min(tstGrps$x) < 2) { shinyjs::disable("saveBtn") }
  if (min(tstGrps$x) >= 2) { shinyjs::enable("saveBtn") }
  output$Include <- renderDT({ Include },
                             FALSE,
                             escape = FALSE,
                             class = "compact",
                             selection = "none",
                             editable = FALSE,
                             rownames = FALSE,
                             options = list(
                               dom = "t",
                               paging = FALSE,
                               ordering = FALSE,
                               autowidth = TRUE,
                               columnDefs = wTest1,
                               scrollX = FALSE
                             ),
                             callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
  output$Msg <- renderUI({ em(" ") })
  output$PCA <- renderPlotly(plot_lyPCA)
  # sapply(seq_len(nrow(Include)), function(x) {
  #   id <- paste0("check___", as.character(x))
  #   observeEvent(input[[id]], {
  #     Exp.map$Use[x] <- as.logical(input[[id]])
  #     assign("Exp.map", Exp.map, envir = .GlobalEnv)
  #     tstGrps <- aggregate(Exp.map$Use, list(Exp.map[[VPAL$column]]), sum)
  #     if (min(tstGrps$x) < 2) {
  #       shinyjs::disable("saveBtn")
  #       output$Msg <- renderUI({ em("Remember: this script requires at least 2 replicates per sample group!!!") })
  #     }
  #     if (min(tstGrps$x) >= 2) {
  #       shinyjs::enable("saveBtn")
  #       output$Msg <- renderUI({ em(" ") })
  #     }
  #   })
  # })
  observeEvent(input$saveBtn, {
    Exp.map$Use <- sapply(seq_len(nrow(Include)), function(x) { input[[paste0("check___", as.character(x))]] })
    tstGrps <- aggregate(Exp.map$Use, list(Exp.map[[VPAL$column]]), sum)
    if (min(tstGrps$x) < 2) {
      #shinyjs::disable("saveBtn")
      output$Msg <- renderUI({ em("Remember: this script requires at least 2 replicates per sample group!!!") })
    }
    if (min(tstGrps$x) >= 2) {
      #shinyjs::enable("saveBtn")
      output$Msg <- renderUI({ em(" ") })
      #
      assign("Exp.map", Exp.map, envir = .GlobalEnv)
      #
      tmpTbl <- Exp.map
      tst <- lapply(colnames(tmpTbl), function(x) { typeof(tmpTbl[[x]]) })
      w <- which(tst == "list")
      if (length(w)) { for (i in w) { tmpTbl[[i]] <- sapply(tmpTbl[[i]], paste, collapse = ";") }}
      tst <- try(write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE), silent = TRUE)
      while ("try-error" %in% class(tst)) {
        dlg_message(paste0("File \"", ExpMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
        tst <- try(write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE), silent = TRUE)
      }
      #
      assign("appRunTest", TRUE, envir = .GlobalEnv)
      stopApp()
    }
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("appRunTest"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
#
#tmp <- read.csv(Param$Experiments.map)
#tmp$Use <- Exp.map$Use
#write.csv(tmp, file = paste0("backup_", Param$Experiments.map), row.names = FALSE)
Exp.map <- Exp.map[which(Exp.map$Use),]
if (LabelType == "Isobaric") { Iso <- sort(unique(Exp.map$Isobaric.set)) }
for (i in 1:nrow(Aggregate.map)) { #i <- 1
  n <- Aggregate.map$Aggregate.Name[i]
  if (nchar(n) == 3) { assign(n, unique(Exp.map[[unlist(Aggregate.map$Characteristics[i])]])) } else {
    assign(n, unique(Exp.map[[n]]))
  }
}
for (i in Param.aggreg) { #i <- Param.aggreg[1]
  a <- get(i)
  if (nchar(a$aggregate) == 3) { a$values <- sort(unique(Exp.map[[a$names]])) } else {
    a$values <- sort(unique(Exp.map[[a$aggregate]]))
  }
  assign(i, a)
}
# Update aliases
RSA <- Ref.Sample.Aggregate
RG <- Ratios.Groups
RRG <- Ratios.Ref.Groups
VPAL <- Volcano.plots.Aggregate.Level

#############################################
# Done!                                     #
#############################################
