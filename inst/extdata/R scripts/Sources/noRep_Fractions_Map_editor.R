################
# MS files map #
################
#
# Edit a map defining the relationship between MS runs and biological samples

if ("PTM-enriched" %in% colnames(FracMap)) {
  FracMap$"PTM-enriched"[which(!FracMap$"PTM-enriched" %in% Modifs$`Full name`)] <- NA
} else { FracMap$"PTM-enriched" <- NA }
# Try sorting automatically
tst <- grepl("_[0-9]+\\.d$", FracMap$`Raw file`)
if (sum(tst)) {
  FracMap$Bruker_run_ID <- NA
  FracMap$Bruker_run_ID[which(tst)] <- as.integer(gsub(".*_|\\.d$", "", FracMap$`Raw file`[which(tst)]))
  w1 <- which(!is.na(FracMap$Bruker_run_ID))
  w2 <- which(is.na(FracMap$Bruker_run_ID))
  w1 <- w1[order(FracMap$Bruker_run_ID[w1])]
  FracMap <- FracMap[c(w1, w2),]
}
FracMap$Use <- as.logical(FracMap$Use)
FracMap$Use[which(is.na(FracMap$Use))] <- TRUE
nr <- nrow(FracMap)
rws <- seq_len(nr)
# Original table column widths
wTest0 <- setNames(vapply(colnames(FracMap), function(k) { #k <- colnames(FracMap)[1]
  tmp <- FracMap[[k]]
  if ("logical" %in% class(tmp)) { tmp <- as.integer(tmp) }
  tmp <- as.character(tmp)
  x <- max(nchar(c(k, tmp)) + 3, na.rm = TRUE)
  x <- x*10
  if (is.na(x)) { x <- 15 } else { x <- max(c(ceiling(x/10)*10, 30)) }
  return(x)
}, 1), colnames(FracMap))
# Dummy table
frMap <- FracMap
frMap$"Raw files name" <- NULL
frMap$Use <- as.logical(toupper(frMap$Use))
frMap$Use[which(is.na(FracMap$Use))] <- TRUE
frMap$Use <- shinyCheckInput(frMap$Use, "Use")
frMap$"PTM-enriched" <- shinySelectInput(FracMap$"PTM-enriched",
                                         "PTMenriched",
                                         unique(c(Modifs$`Full name`, NA)),
                                         paste0(30*max(c(nchar(as.character(Modifs$`Full name`)), 2)), "px"))
frMap$"Parent sample" <- shinyTextInput(frMap$"Parent sample", "Sample", paste0(wTest0["Parent sample"], "px"))
k <- c("Raw file", "Parent sample", "Fraction", "PTM-enriched", "Bruker_run_ID", "Use")
k <- k[which(k %in% colnames(frMap))]
frMap <- frMap[, k]
# Estimate dummy table column widths
wTest1 <- vapply(colnames(frMap), function(k) { #k <- colnames(frMap)[1]
  if ((k == "Parent sample")&&(!k %in% names(wTest0))) { k <- "MQ.Exp" }
  if (k %in% names(wTest0)) { x <- wTest0[k] } else { x <- 30 }
  return(x)
}, 1)
wTest2 <- sum(wTest1) + 15 + ncol(frMap)*5
wTest1 <- paste0(as.character(wTest1), "px")
wTest1 <- aggregate((seq_along(wTest1))-1, list(wTest1), c)
wTest1 <- apply(wTest1, 1, function(x) {
  x2 <- as.integer(x[[2]])
  list(width = x[[1]],
       targets = x2,
       names = colnames(frMap)[x2+1])
})
#
appNm <- paste0(dtstNm, " - ", FracMapNm)
if (exists("frMap2")) { rm(frMap2) }
ui <- shinyUI(fluidPage(titlePanel(tag("u", FracMapNm),
                                   appNm),
                        useShinyjs(),
                        setBackgroundColor( # Doesn't work
                          color = c(#"#F8F8FF",
                            "#EFE6F5"),
                          gradient = "linear",
                          direction = "bottom"
                        ),
                        extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
                        mainPanel(
                          withSpinner(DT::DTOutput("myFracMap", width = wTest2)),
                          br(),
                          actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                        )))
server <- function(input, output, session) {
  frMap2 <- FracMap
  output$myFracMap <- DT::renderDT({ frMap },
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
                                                  scrollX = FALSE),
                                   # the callback is essential to capture the inputs in each row
                                   callback = JS("table.rows().every(function(i, tab, row) {
  var $this = $(this.node());
  $this.attr('id', this.data()[0]);
  $this.addClass('shiny-input-container');
});
Shiny.unbindAll(table.table().node());
Shiny.bindAll(table.table().node());"))
  observeEvent(input$myFracMap_cell_edit, { # This is not goddamn working!!!
    kl <- colnames(frMap)[input$myFracMap_cell_edit$col+1]
    if (kl %in% colnames(frMap2)) {
      frMap2[input$myFracMap_cell_edit$row, kl] <<- input$myFracMap_cell_edit$value
    } else {
      warning(paste0("Could not find column ", kl, " from dummy table in the real table!"))
    }
  })
  observeEvent(input$saveBtn, {
    for (k in c("Use", "PTM-enriched", "Parent sample")) {
      root <- gsub("-", "", k)
      if (k == "Use") { tp <- TRUE }
      if (k %in% c("PTM-enriched", "Parent sample")) {
        tp <- ""
        if (k == "Parent sample") { root<- "Sample" }
      }
      frMap2[[k]] <- vapply(rws, function(x) { input[[paste0(root, "___", as.character(x))]] }, tp)
    }
    assign("frMap2", frMap2, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("frMap2"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
FracMap %<o% frMap2
#
tst <- try(write.csv(FracMap, file = FracMapPath, row.names = FALSE), silent = TRUE)
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) { # We only want this to happen if the file is locked for editing
  dlg_message(paste0("File \"", FracMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(FracMap, file = FracMapPath, row.names = FALSE), silent = TRUE)
}
#
msg <- "Are you happy with the edits? (click no to try again)"
tstFrMp <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno", rstudio = TRUE)$res, c("yes", "no"))]
