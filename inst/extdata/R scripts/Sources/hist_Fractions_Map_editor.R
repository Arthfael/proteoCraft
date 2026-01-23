######################################
# Edit samples map - Histones script #
######################################

appNm <- paste0(dtstNm, " - Samples map")
smplsMapFl <- paste0(dstDir, "/Samples_map.csv")
reLoad <- file.exists(smplsMapFl)
dfltKol <- c("Sample", "Group", "Replicate")
if (reLoad) {
  smplsMap <- read.csv(smplsMapFl,
                         check.names = FALSE)
  # Backwards compatibility
  if (("Old" %in% colnames(smplsMap))&&(!"Sample" %in% colnames(smplsMap))) {
    colnames(smplsMap)[which(colnames(smplsMap) == "Old")] <- "Sample"
  }
  if (("New" %in% colnames(smplsMap))&&(!"Sample name" %in% colnames(smplsMap))) {
    colnames(smplsMap)[which(colnames(smplsMap) == "New")] <- "Sample name"
  }
  #
  reLoad <- sum(dfltKol %in% colnames(smplsMap)) == 3
}
if (!reLoad) {
  smplsMap <- data.frame(Sample = unique(ev$`Raw file`),
                         Group = "?",
                         Replicate = "?",
                         Use = TRUE)
}
dfltFact <- c()
if ((exists("myFact"))&&(length(myFact))&&(is.character(myFact))) {
  dfltFact <- unique(c(dfltFact, myFact))
}
kol <- colnames(smplsMap)
kol <- kol[which(!kol %in% c(dfltKol, "Use", "Sample name"))]
if (length(kol)) {
  dfltFact <- unique(c(dfltFact, kol))
}
if (!length(dfltFact)) {
  dfltFact <- c("Fact1", "Fact2")
}
dfltFact <- paste(dfltFact, collapse = " / ")
myFact <- dlg_input("Enter name of columns you would like to add to the table (in addition to Sample, Replicate and Group; use \" / \" as separator)",
                    dfltFact)$res
myFact <- unlist(strsplit(myFact, " +/ +"))
myFact <- unique(myFact)
myFact <- myFact[which(myFact != "")]
myFact <- myFact[which(!myFact %in% c(dfltKol, "Use", "Sample name"))]

if (!"Sample name" %in% colnames(smplsMap)) {
  smplsMap$"Sample name" <- smplsMap$Sample
}
if (!"Use" %in% colnames(smplsMap)) {
  smplsMap$Use <- TRUE
}
smplsMap$Use <- as.logical(smplsMap$Use)
w <- which(is.na(smplsMap$Use))
smplsMap$Use[w] <- TRUE

nr <- nrow(smplsMap)
rws <- seq_len(nr)
chRws <- as.character(rws)
if (!exists("nRep")) { nRep <- nr }
dflt_Rpl <- 1:(nr+nRep) %% nRep
dflt_Rpl[which(dflt_Rpl == 0)] <- nRep
#
if (length(myFact)) {
  w <- which(!myFact %in% colnames(smplsMap))
  if (length(w)) {
    smplsMap[, myFact[w]] <- ""
  }
  # for (fct in myFact) {
  #   wNA <- which(is.na(smplsMap[[fct]]))
  #   smplsMap[wNA, fct] <- ""
  # }
}
#
smplsMap2 <- smplsMap
smplsMap2$`Sample name` <- shinyTextInput(smplsMap2$`Sample name`, "Sample name", width = "100%")
smplsMap2$Group <- shinyTextInput(smplsMap2$Group, "Group", width = "100%")
smplsMap2$Group___FD <- shinyFDInput("Group", nr)
smplsMap2$Replicate <- shinyNumInput(smplsMap2$Replicate, 1, Inf, 1, root = "Replicate", width = "100%")
smplsMap2$Replicate___FD <- shinyFDInput("Replicate", nr)
smplsMap2$Replicate___INCR <- shinyPlusFD("Replicate", nr)
smplsMap2$Use <- shinyCheckInput(smplsMap2$Use, "Use")
smplsMap2$Use___FD <- shinyFDInput("Use", nr)
kol <- c("Sample", "Sample name", "Group", "Group___FD")
if (length(myFact)) {
  for (fct in myFact) {
    smplsMap2[[fct]] <- shinyTextInput(smplsMap2[[fct]], fct, width = "100%")
    fdKl <- paste0(fct, "___FD")
    smplsMap2[[fdKl]] <- shinyFDInput(fct, nr)
    kol <- c(kol, fct, fdKl)
  }
}
kol <- c(kol, "Replicate", "Replicate___FD", "Replicate___INCR", "Use", "Use___FD")
smplsMap2 <- smplsMap2[, kol]
ALLIDS <- as.character(sapply(c("Group", myFact, "Replicate", "Use"), function(x) {
  paste0(x, "___", chRws))
}))
idsL <- length(ALLIDS)
# Table width
wTest0 <- setNames(vapply(colnames(smplsMap), function(k) { #k <- colnames(smplsMap2)[1]
  tmp <- smplsMap[[k]]
  if ("logical" %in% class(tmp)) { tmp <- as.integer(tmp) }
  tmp <- as.character(tmp)
  x <- max(nchar(c(k, tmp)) + 3, na.rm = TRUE)
  x <- x*10
  return(x)
}, 1), colnames(smplsMap))
wTest1 <- vapply(colnames(smplsMap2), function(k) { #k <- colnames(smplsMap2)[1]
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
       names = colnames(smplsMap2)[x2+1])
})
#
g <- grep("___((FD)|(INCR))$", colnames(smplsMap2))
colnames(smplsMap2)[g] <- ""
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
                h5("Edit the relationship between different samples, defining sample groups with replicate numbers, then click \"Save\"."),
                br(),
                actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                br(),
                DTOutput("smplsMap2"),
                br(),
                #plotlyOutput("PCA", height = paste0(screenRes$width*0.4, "px")),
                br())
server <- function(input, output, session) {
  smplsMap3 <- smplsMap
  output$smplsMap2 <- renderDT({ smplsMap2 },
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
                       if (fct %in% c("Group", myFact)) {
                         updateTextInput(session, idK, NULL, x)
                       }
                       if (fct == "Replicate") {
                         updateNumericInput(session, idK, NULL, x, 1, Inf, 1)
                       }
                       if (fct == "Use") {
                         updateCheckboxInput(session, idK, NULL, as.logical(x))
                       }
                     }
                   })
    }
  })
  # Incremental fill-down for replicates
  sapply(rws, function(i) {
    if (i < nr) {
      iChr <- as.character(i)
      id1 <- paste0("Replicate___", iChr)
      id2 <- paste0("Replicate___", iChr, "___INCR")
      observeEvent(input[[id2]],
                   {
                     x <- input[[id1]]
                     rplRg <- (i+1):nr
                     l <- length(rplRg)
                     m <- match(x, dflt_Rpl)+1
                     rplVal <- dflt_Rpl[m:(m+l-1)]
                     for (k in rplRg) {
                       y <- rplVal[k-i]
                       idK <- paste0("Replicate___", as.character(k))
                       updateSelectInput(session, idK, NULL, 1:nRep, y)
                     }
                   }
      )
    }
  })
  #output$PCA <- renderPlotly(plot_lyPCA)
  observeEvent(input$saveBtn, {
    kls <- c("Sample name", "Group", "Replicate", "Use")
    for (kl in kls) {
      smplsMap3[[kl]] <- sapply(rws, function(i) {
        input[[paste0(kl, "___", as.character(i))]]
      })
      if (kl == "Replicate") {
        smplsMap3[[kl]] <- as.integer(smplsMap3[[kl]])
      }
      if (kl == "Use") {
        smplsMap3[[kl]] <- as.logical(smplsMap3[[kl]])
      }
    }
    assign("smplsMap3", smplsMap3, envir = .GlobalEnv)
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
samplesMap <- smplsMap
samplesMap[, colnames(smplsMap3)] <- smplsMap3
samplesMap <- samplesMap[which(samplesMap$Use),]
#
tmpTbl <- smplsMap3
tst <- lapply(colnames(smplsMap3), function(x) { typeof(smplsMap3[[x]]) })
w <- which(tst == "list")
if (length(w)) { for (i in w) { smplsMap3[[i]] <- vapply(smplsMap3[[i]], paste, "", collapse = ";") }}
tst <- try(write.csv(tmpTbl, file = smplsMapFl, row.names = FALSE, quote = TRUE), silent = TRUE)
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
  dlg_message(paste0("File \"", smplsMapFl, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(tmpTbl, file = smplsMapFl, row.names = FALSE, quote = TRUE), silent = TRUE)
}
#############################################
# Done!                                     #
#############################################
