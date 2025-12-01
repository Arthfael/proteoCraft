################################################
# Edit comparison groups map - Histones script #
################################################

appNm <- paste0(dtstNm, " - Groups map")
grpsMapFl <- paste0(dstDir, "/Groups_map.csv")
reLoad <- file.exists(grpsMapFl)
if (reLoad) {
  grpsMap <- read.csv(grpsMapFl,
                      check.names = FALSE)
  #
  reLoad <- (sum(c("Group", "Reference", "Comparison_group") %in% colnames(grpsMap)) == 3)&&
    (sum(!Groups %in% grpsMap$Group) == 0)
}
if (!reLoad) {
  grpsMap <- data.frame(Group = Groups,
                        Reference = FALSE,
                        Comparison_group = 1)
}
grpsMap$Reference <- as.logical(grpsMap$Reference)
w <- which(is.na(grpsMap$Reference))
grpsMap$Reference[w] <- TRUE
#
nr <- nrow(grpsMap)
rws <- seq_len(nr)
#
grpsMap2 <- grpsMap
grpsMap2$Comparison_group <- shinyTextInput(grpsMap2$Comparison_group, "Comparison_group", width = "100%")
grpsMap2$Comparison_group___FD <- shinyFDInput("Comparison_group", nr)
grpsMap2$Reference <- shinyCheckInput(grpsMap2$Reference, "Reference")
grpsMap2$Reference___FD <- shinyFDInput("Reference", nr)
grpsMap2 <- grpsMap2[, c("Group", "Comparison_group", "Comparison_group___FD",
                         "Reference", "Reference___FD")]
ALLIDS <- as.character(sapply(c("Comparison_group", "Reference"), function(x) {
  paste0(x, "___", as.character(rws))
}))
idsL <- length(ALLIDS)
# Table width
wTest0 <- setNames(vapply(colnames(grpsMap), function(k) { #k <- colnames(grpsMap2)[1]
  tmp <- grpsMap[[k]]
  if ("logical" %in% class(tmp)) { tmp <- as.integer(tmp) }
  tmp <- as.character(tmp)
  x <- max(nchar(c(k, tmp)) + 3, na.rm = TRUE)
  x <- x*10
  return(x)
}, 1), colnames(grpsMap))
wTest1 <- vapply(colnames(grpsMap2), function(k) { #k <- colnames(grpsMap2)[1]
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
       names = colnames(grpsMap2)[x2+1])
})
#
g <- grep("___((FD)|(INCR))$", colnames(grpsMap2))
colnames(grpsMap2)[g] <- ""
#
if (exists("appRunTest")) { rm(appRunTest) }
ui <- fluidPage(useShinyjs(),
                setBackgroundColor( # Doesn't work
                  color = c(#"#F8F8FF",
                    "#E5EDE1"),
                  gradient = "linear",
                  direction = "bottom"
                ),
                extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
                tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
                titlePanel(tag("u", "Samples map"),
                           appNm),
                br(),
                h5("Edit the relationship between sample groups, defining comparison groups which must include exactly one reference group and any number of non-reference groups, then click \"Save\"."),
                br(),
                actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                br(),
                DTOutput("grpsMap2"),
                br(),
                br())
server <- function(input, output, session) {
  grpsMap3 <- grpsMap
  output$grpsMap2 <- renderDT({ grpsMap2 },
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
  # Fill down
  sapply(1:idsL, function(x) {
    id1 <- ALLIDS[x]
    id2 <- paste0(ALLIDS[x], "___FD")
    tmp <- unlist(strsplit(id1, "___"))
    fct <- tmp[[1]]
    i <- as.integer(tmp[[2]])
    if (i < nr) {
      observeEvent(input[[id2]],
                   {
                     x <- input[[id1]]
                     for (k in (i+1):nr) {
                       idK <- paste0(fct, "___", as.character(k))
                       if (fct == "Comparison_group") {
                         updateTextInput(session, idK, NULL, x)
                       }
                       if (fct == "Reference") {
                         updateCheckboxInput(session, idK, NULL, as.logical(x))
                       }
                     }
                   })
    }
  })
  observeEvent(input$saveBtn, {
    kls <- c("Comparison_group", "Reference")
    for (kl in kls) {
      grpsMap3[[kl]] <- sapply(rws, function(i) {
        input[[paste0(kl, "___", i)]]
      })
      if (kl == "Reference") {
        grpsMap3[[kl]] <- as.logical(grpsMap3[[kl]])
      }
    }
    assign("grpsMap3", grpsMap3, envir = .GlobalEnv)
    #
    assign("appRunTest", TRUE, envir = .GlobalEnv)
    stopApp()
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
#
tmpTbl <- grpsMap3
tst <- lapply(colnames(grpsMap3), function(x) { typeof(grpsMap3[[x]]) })
w <- which(tst == "list")
if (length(w)) { for (i in w) { grpsMap3[[i]] <- vapply(grpsMap3[[i]], paste, "", collapse = ";") }}
tst <- try(write.csv(tmpTbl, file = grpsMapFl, row.names = FALSE, quote = TRUE), silent = TRUE)
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
  dlg_message(paste0("File \"", grpsMapFl, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(tmpTbl, file = grpsMapFl, row.names = FALSE, quote = TRUE), silent = TRUE)
}
#
groupsMap <- grpsMap
groupsMap[, colnames(grpsMap3)] <- grpsMap3
#groupsMap <- groupsMap[which(groupsMap$Use),]
groupsMap <- groupsMap[which(groupsMap$Group %in% samplesMap$Group),]
#
#############################################
# Done!                                     #
#############################################
