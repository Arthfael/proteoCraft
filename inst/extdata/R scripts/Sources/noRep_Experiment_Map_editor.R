###############
# Samples map #
###############
#
# Map samples to comparison groups, define reference or negative filter samples...

SamplesMapNm %<o% "Samples map"
SamplesMapPath %<o% paste0(wd, "/", SamplesMapNm, ".csv")
if (file.exists(SamplesMapPath)) {
  SamplesMap <- data.table::fread(SamplesMapPath, check.names = FALSE, data.table = FALSE)
  if (("MQ.Exp" %in% colnames(SamplesMap))&&(!"Parent sample" %in% colnames(SamplesMap))) {
    SamplesMap$"Parent sample" <- SamplesMap$"MQ.Exp"
  }
}
exp <- unique(FracMap$`Parent sample`) # Update
if ((!file.exists(SamplesMapPath))||((exists("SamplesMap"))&&(nrow(SamplesMap) < length(unique(FracMap$`Parent sample`))))) {
  SamplesMap <- data.frame("Parent sample" = exp,
                           "Negative Filter" = FALSE,
                           check.names = FALSE)
}
if (MakeRatios) {
  if (!"Ratios group" %in% colnames(SamplesMap)) { SamplesMap$"Ratios group" <- 1 }
  if (!"Reference" %in% colnames(SamplesMap)) { SamplesMap$Reference <- FALSE }
} else {
  SamplesMap$"Ratios group" <- NULL
  SamplesMap$Reference <- NULL
}
SamplesMap %<o% SamplesMap
nr <- nrow(SamplesMap)
rws <- 1:nr
# if ("Order" %in% colnames(SamplesMap)) {
#   u <- unique(SamplesMap$Order)
#   u <- u[which(!u %in% rws)]
#   if (length(u) < nr) { SamplesMap$Order <- rws }
# } else { SamplesMap$Order <- rws }
smplMapKol1 <- c("Negative Filter", "Use")
if (MakeRatios) {
  smplMapKol1 <- c("Reference", smplMapKol1)
}
for (kol in smplMapKol1) {
  if (kol %in% colnames(SamplesMap)) {
    SamplesMap[[kol]] <- as.logical(toupper(SamplesMap[[kol]]))
    SamplesMap[which(is.na(SamplesMap[[kol]])), kol] <- c(FALSE, TRUE)[(kol == "Use")+1]
  }
}
if (!"Use" %in% colnames(SamplesMap)) { SamplesMap$Use <- TRUE }
SamplesMap$Use <- suppressWarnings(as.logical(SamplesMap$Use))
SamplesMap$Use[which(is.na(SamplesMap$Use))] <- TRUE
smplMapKol <- smplMapKol1
if (MakeRatios) {
  smplMapKol <- c("Ratios group", smplMapKol)
}
allIDs <- as.character(sapply(smplMapKol, function(x) {
  paste0(x, "___", rws)
}))
# Original table column widths
wTest0 <- setNames(vapply(colnames(SamplesMap), function(k) { #k <- colnames(SamplesMap)[1]
  l <- k
  if (k == "MQ.Exp") { l <- "Parent sample"}
  tmp <- SamplesMap[[k]]
  if ("logical" %in% class(tmp)) { tmp <- as.integer(tmp) }
  tmp <- as.character(tmp)
  x <- max(nchar(c(k, tmp)) + 3, na.rm = TRUE)
  x <- x*10
  if (is.na(x)) { x <- 30 } else { x <- max(c(ceiling(x/10)*10, 30)) }
  return(x)
}, 1), colnames(SamplesMap))
#
if (exists("expOrder")) {
  tst <- sum(!exp %in% expOrder)
  if (sum(tst)) { rm(expOrder) }
}
if (!exists("expOrder")) {
  expOrder <- exp
}
SamplesMap <- SamplesMap[match(expOrder, SamplesMap$"Parent sample"),]
#
smplMap2 <- smplMap <- SamplesMap[, which(!colnames(SamplesMap) %in% "New name")] # This column is deprecated and ignored
#colnames(smplMap2)[which(colnames(smplMap2) == "MQ.Exp")] <- "Parent sample"
# Estimate table column widths
wTest1 <- vapply(colnames(smplMap2), function(k) { #k <- colnames(smplMap2)[1]
  if (k == "Parent sample") { k <- "MQ.Exp" }
  if (k %in% names(wTest0)) { x <- wTest0[k] } else { x <- 30 }
  return(x)
}, 1)
wTest2 <- max(c(sum(wTest1) + 15 + ncol(smplMap2)*5, 600))
wTest1 <- paste0(as.character(wTest1), "px")
wTest1 <- aggregate((seq_along(wTest1))-1, list(wTest1), c)
wTest1 <- apply(wTest1, 1, function(x) {
  x2 <- as.integer(x[[2]])
  list(width = x[[1]],
       targets = x2,
       names = colnames(smplMap2)[x2+1])
})
#
IDs <- c()
MSG <- reactiveVal("")
for (kol in c("Reference", "Negative Filter", "Use")) {
  if (kol %in% colnames(SamplesMap)) {
    smplMap2[[kol]] <- shinyCheckInput(SamplesMap[[kol]], kol)
    IDs <- c(IDs, paste0(kol, "___", seq_len(nrow(smplMap2))))
    stopifnot(length(IDs) == length(unique(IDs)))
  }
}
if (MakeRatios) {
  kol <- "Ratios group"
  smplMap2[[kol]] <- shinyNumInput(SamplesMap[[kol]], 1, Inf, 1, 1, root = kol)
  IDs <- c(IDs, paste0(kol, "___", seq_len(nrow(smplMap2))))
  stopifnot(length(IDs) == length(unique(IDs)))
}

#
msg <- ""
appNm <- paste0(dtstNm, " - Experiment map")
ui <- fluidPage(useShinyjs(),
                setBackgroundColor( # Doesn't work
                  color = c(#"#F8F8FF",
                    "#EBEFF7"),
                  gradient = "linear",
                  direction = "bottom"
                ),
                titlePanel(tag("u", "Experiment map"),
                           appNm),
                extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
                tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
                mainPanel(h4("Optional: define samples order"),
                          selectInput("expOrder", "",
                                      expOrder, expOrder, TRUE, TRUE, width = "1200px"),
                          br(),
                          br(),
                          span(uiOutput("Message"), style = "color:red"),
                          withSpinner(DT::DTOutput("mySampleMap", width = wTest2))),
                br(),
                actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"))
if (exists("smplMap3")) { rm(smplMap3) }
server <- function(input, output, session) {
  # Create copies table
  smplMap3 <- smplMap
  output$Message <- renderUI({ em(" ") })
  output$mySampleMap <- DT::renderDT({ smplMap2 },
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
                                                    scrollX = TRUE),
                                     callback = JS("table.rows().every(function(i, tab, row) {
  var $this = $(this.node());
  $this.attr('id', this.data()[0]);
  $this.addClass('shiny-input-container');
});
Shiny.unbindAll(table.table().node());
Shiny.bindAll(table.table().node());"))
  #
  # Observe events
  # if (length(IDs)) {
  #   sapply(IDs, function(id) {
  #     x <- gsub(".+_", "", id)
  #     kol <- substr(id, 1, nchar(id)-nchar(x)-1)
  #     x <- as.integer(x)
  #     observeEvent(input[[id]], {
  #       cat(input[[id]], "\n")
  #       val <- input[[id]]
  #       if (kol %in% c("Reference", "Negative Filter", "Use")) {
  #         val <<- as.logical(val)
  #       }
  #       SamplesMap[x, kol] <<- val
  #       if (kol == "Use") {
  #         tst <- sum(SamplesMap[kol])
  #         if (!tst) {
  #           msg <- "Include at least one sample!"
  #           shinyjs::disable("saveBtn")
  #         }
  #         if (tst) {
  #           msg <- " "
  #           shinyjs::enable("saveBtn")
  #         }
  #       }
  #       if ((MakeRatios)&&(kol == "Reference")) {
  #         l <- length(unique(SamplesMap[[kol]]))
  #         if (l != 2) {
  #           msg <- "This workflow compares different samples, select one sample as reference!"
  #           shinyjs::disable("saveBtn")
  #         }
  #         if (l == 2) {
  #           msg <- " "
  #           shinyjs::enable("saveBtn")
  #         }
  #       }
  #       output$Message <- renderUI({ em(msg) })
  #     })
  #   })
  # }
  observeEvent(input$mySampleMap_cell_edit, {
    smplMap3[input$mySampleMap_cell_edit$row,
             input$mySampleMap_cell_edit$col+1] <- input$mySampleMap_cell_edit$value
  })
  observeEvent(input$expOrder, {
    tmp <- input$expOrder
    tmp <- c(tmp, exp[which(!exp %in% tmp)])
    assign("expOrder", tmp, envir = .GlobalEnv)
  })
  # sapply(allIDs, function(id) {
  #   tmp <- unlist(strsplit(id, "___"))
  #   kol <- tmp[1]
  #   i <- tmp[2]
  #   observeEvent(input[[id]], {
  #     #cat("Event", id, "\n")
  #     #cat(input[[id]], "\n")
  #     smplMap3[i, kol] <<- as.logical(input[[id]])
  #   })
  # })
  observeEvent(input$saveBtn, {
    for (k in smplMapKol) {
      smplMap3[[k]] <- sapply(rws, function(x) { input[[paste0(k, "___", x)]] })
    }
    assign("smplMap3", smplMap3, envir = .GlobalEnv)
    tmp <- input$expOrder
    tmp <- c(tmp, exp[which(!exp %in% tmp)])
    assign("expOrder", tmp, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("smplMap3"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
SamplesMap %<o% smplMap3
exp <- expOrder
SamplesMap <- SamplesMap[match(expOrder, SamplesMap$"Parent sample"),]
#
tmp <- SamplesMap
w <- which(vapply(colnames(tmp), function(x) { "list" %in% class(tmp[[x]]) }, TRUE))
if (length(w)) { for (i in w) { tmp[[i]] <- vapply(tmp[[i]], paste, "", collapse = ";") } }
tst <- try(write.csv(tmp, file = SamplesMapPath, row.names = FALSE), silent = TRUE)
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) { # We only want this to happen if the file is locked for editing
  dlg_message(paste0("File \"", SamplesMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(tmp, file = SamplesMapPath, row.names = FALSE), silent = TRUE)
}

if ("Negative Filter" %in% colnames(SamplesMap)) {
  SamplesMap$"Negative Filter" <- as.logical(toupper(SamplesMap$"Negative Filter"))
}
AnalysisParam$"Ratios analysis" <- MakeRatios
# Experiment column, synonym of "Parent sample" (but keep both)
SamplesMap$Experiment <- SamplesMap$"Parent sample"
#
msg <- "Are you happy with the edits? (click no to try again)"
tstXpMp <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno", rstudio = TRUE)$res, c("yes", "no"))]
