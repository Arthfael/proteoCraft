# A sourced script to edit the Experiment map (table defining experimental structure)
#
library(shiny)
library(shinyjs)
library(DT)
#
frKl <- c("MQ.Exp", "Parent sample")
frKl <- frKl[which(frKl %in% colnames(FracMap))[1]]
Factors2 %<o% Factors[which(!Factors %in% c("Experiment",
                                            "Replicate",
                                            "Isobaric.set",
                                            "Time.point"))]
#
ExpMapNm %<o% "Experiment map"
ExpMapPath %<o% paste0(wd, "/", ExpMapNm, ".csv")
#
if (file.exists(ExpMapPath)) {
  reLoad <- TRUE
  if (exists("ExpMap")) {
    reLoad <- c(TRUE, FALSE)[match(dlg_message("Reload Experiment map from disk?\n(you will still be able to edit it)", "yesno")$res, c("yes", "no"))]
  }
  if (reLoad) {
    ExpMap <- read.csv(ExpMapPath, check.names = FALSE)
    colnames(ExpMap)[which(colnames(ExpMap) == "Sample.name")] <- "Sample name" # Backwards compatibility
    tst <- unique(FracMap[[frKl]][which(FracMap$Use)])
    if (sum(!tst %in% ExpMap$"Sample name")) {
      warning("Invalid Experiment map reloaded, ignoring...")
      rm(ExpMap)
    }
  }
}
if (!exists("ExpMap")) {
  if (LabelType == "LFQ") {
    ExpMap <- data.frame("Experiment" = "?",
                         "Replicate" = "?",
                         "MQ.Exp" = MQ.Exp,
                         "Reference" = FALSE,
                         "Sample name" = MQ.Exp,
                         "Use" = TRUE,
                         check.names = FALSE)
  }
  if (LabelType == "Isobaric") {
    ExpMap <- data.frame("Experiment" = "?",
                         "Replicate" = "?",
                         "MQ.Exp" = rep(MQ.Exp, length(get(IsobarLab))),
                         "Reference" = FALSE,
                         "Sample name" = rep(MQ.Exp, length(get(IsobarLab))),
                         "Isobaric.label" = rep(get(IsobarLab), length(MQ.Exp)),
                         "Isobaric.label.details" = rep(IsobarLabDet, length(MQ.Exp)),
                         "Isobaric.set" = "?",
                         "Use" = TRUE,
                         check.names = FALSE)
  }
}
if ((LocAnalysis)&&(!"Proportion" %in% colnames(ExpMap))) { ExpMap$Proportion <- 1 }
for (Fact in Factors2) { if (!Fact %in% colnames(ExpMap)) { ExpMap[[Fact]] <- "?" } }
tmp <- aggregate(FracMap$Fraction, list(FracMap[[frKl]]), function(x) { paste(sort(unique(x)), collapse = ";") })
ExpMap$Fractions <- tmp$x[match(ExpMap$MQ.Exp, tmp$Group.1)]
ExpMap$Use <- as.logical(sapply(strsplit(ExpMap$MQ.Exp, ";"), function(x) { max(FracMap$Use[match(x, FracMap[[frKl]])]) }))
tmpTbl <- ExpMap
tst <- lapply(colnames(tmpTbl), function(x) { typeof(tmpTbl[[x]]) })
w <- which(tst == "list")
if (length(w)) { for (i in w) { tmpTbl[[i]] <- sapply(tmpTbl[[i]], paste, collapse = ";") }}
tst <- try(write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE), silent = TRUE)
if ("try-error" %in% class(tst)) {
  dlg_message(paste0("File \"", ExpMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE)
}
ExpMap <- ExpMap[which(ExpMap$MQ.Exp %in% FracMap[[frKl]]),]
#
# Edit map
ExpData <- read.csv(ExpMapPath, check.names = FALSE)
ExpData <- ExpData[which(ExpData$MQ.Exp %in% FracMap[[frKl]]),]
for (Fact in Factors[which(!Factors %in% colnames(ExpData))]) { ExpData[[Fact]] <- "?" }
ExpData$Use <- as.logical(sapply(strsplit(ExpData$MQ.Exp, ";"), function(x) { max(FracMap$Use[match(x, FracMap[[frKl]])]) }))
if (LocAnalysis) {
  if (!"Proportion" %in% colnames(ExpData)) {
    ExpData$Proportion <- 1
  } else {
    ExpData$Proportion <- suppressWarnings(as.numeric(ExpData$Proportion))
    w <- which(is.na(ExpData$Proportion))
    if (length(w)) {
      warning("Invalid values in \"Proportion\" column!")
      ExpData$Proportion[w] <- 1
    }
  }
}
tst <- sapply(FactorsLevels, length)
Fact1 <- Factors[which(tst == 1)]
Fact2 <- Factors[which(tst > 1)]
Others <- c("MQ.Exp", "Sample name")
Others <- Others[which(Others %in% colnames(ExpData))]
nr <- nrow(ExpData)
rws <- seq_len(nr)
OtherIDs <- setNames(lapply(Others, function(x) { paste0(x, "___", rws)} ), Others)
Fact2IDs <- setNames(lapply(Fact2, function(x) { paste0(x, "___", rws)} ), Fact2)
AllIDs <- append(OtherIDs, Fact2IDs)
AllIDs$Use <- paste0("Use___", rws)
ALLIDS <- setNames(unlist(AllIDs), NULL)
for (fct in Fact1) { #fct <- Fact1[1]
  # (Since the table can be manually edited too, we want to make sure to correct any typos in our single level factor columns)
  ExpData[[fct]] <- rep(FactorsLevels[[fct]], nr)
}
L <- length(FactorsLevels[["Replicate"]])
dflt_Rpl <- 1:(nr+L) %% L
dflt_Rpl[which(dflt_Rpl == 0)] <- L
dflt_Rpl <- FactorsLevels[["Replicate"]][dflt_Rpl] # Easy template
k1 <- c("MQ.Exp", "Experiment", "Replicate")
k2 <- colnames(ExpData)
k2 <- k2[which(!k2 %in% c(k1, "Sample name", "Use"))]
ExpData <- ExpData[which(ExpData$MQ.Exp %in% MQ.Exp), ]
L <- length(ALLIDS)
#
# Original table width
wTest0 <- setNames(sapply(colnames(ExpData), function(k) { #k <- colnames(ExpData)[1]
  tst <- k %in% Fact2
  x <- max(c(nchar(k),
             nchar(as.character(ExpData[[k]])) + 3 + 3*tst), na.rm = TRUE)
  if (tst) {
    x <- max(c(x, nchar(FactorsLevels[[k]]) + 6))
  }
  x <- x*10
  if (is.na(x)) { x <- 15 } else { x <- max(c(ceiling(x/10)*10, 30)) }
  return(x)
}), colnames(ExpData))
#
# Dummy table for app
tst <- sapply(FactorsLevels, length)
Fact1 <- Factors[which(tst == 1)]
Fact2 <- Factors[which(tst > 1)]
#Fact2 <- "Developmental.stage"
myKol <- unique(c(Others, Factors, "Sample name"))
ExpData2 <- ExpData[, myKol]
wMQExp <- which(colnames(ExpData2) == "MQ.Exp")
colnames(ExpData2)[wMQExp] <- "Parent sample"
ExpData2 <- ExpData2[, c("Parent sample", colnames(ExpData2)[which(colnames(ExpData2) != "Parent sample")])]
kol <- c()
for (fct in Fact2) { #fct <- Fact2[1]
  IDs <- Fact2IDs[[fct]]
  lvls <- FactorsLevels[[fct]]
  lvls2 <- lvls[which(!is.na(lvls))]
  ExpData2[[fct]] <- shinySelectInput(ExpData[[fct]],
                                      fct,
                                      lvls,
                                      paste0(wTest0[fct], "px"),
                                      fct != "Replicate")
  fdNm <- paste0(fct, "___FD")
  wTest0[fdNm] <- 15
  ExpData2[[fdNm]] <- shinyFDInput(fct, nr, TRUE, paste0(wTest0[fdNm], "px"))
  kol <- c(kol, fct, fdNm)
  if (fct == "Replicate") {
    incrNm <- paste0(fct, "___INCR")
    wTest0[incrNm] <- 15
    ExpData2[[incrNm]] <- shinyPlusFD(fct, nr, width = paste0(wTest0[incrNm], "px"))
    kol <- c(kol, incrNm)
  }
}
if (LocAnalysis) {
  ExpData2$Proportion <- shinyPropInput(ExpData$Proportion)
  fdNm <- "Proportion___FD"
  wTest0[fdNm] <- 15
  ExpData2[[fdNm]] <- shinyFDInput("Proportion", nr, TRUE, paste0(wTest0[fdNm], "px"))
  kol <- c(kol, "Proportion", fdNm)
}
wMQExpWdth <- paste0(as.character(max(nchar(ExpData2$"Parent sample"))*10), "px")
kol2 <- colnames(ExpData2)[which(!colnames(ExpData2) %in% c(kol, "Use", "Sample name"))]
ExpData2$Use <- shinyCheckInput(ExpData$Use,
                                "Use")
ExpData2$Use___FD <- shinyFDInput("Use", nr, TRUE)
ExpData2 <- ExpData2[, c(kol2, kol, "Use", "Use___FD", "Sample name")]
# Update table width
wTest1 <- sapply(colnames(ExpData2), function(k) { #k <- colnames(ExpData2)[1]
  if (k == "Parent sample") { k <- "MQ.Exp" }
  if (k %in% names(wTest0)) { x <- wTest0[k] } else { x <- 30 }
  return(x)
})
wTest2 <- sum(wTest1)
wTest1 <- paste0(as.character(wTest1), "px")
wTest1 <- aggregate((1:length(wTest1))-1, list(wTest1), c)
wTest1 <- apply(wTest1, 1, function(x) {
  list(width = x[[1]],
       targets = x[[2]],
       names = colnames(ExpData2)[x[[2]]+1])
})
#
g <- grep("___((FD)|(INCR))$", colnames(ExpData2))
colnames(ExpData2)[g] <- "" 
#
facLevels2 <- lapply(FactorsLevels, function(x) {
  unique(c(x, NA))
})
facLevels2$Use <- c(TRUE, FALSE)
#
edith <- list(target = "column",
              disable = list(columns = c(0, match(Fact1, colnames(ExpData2))-1)))
tmp <- c(0:(ncol(ExpData2)-1))
tmp <- tmp[which(!tmp %in% edith$disable$columns)]
edith$enable <- list(columns = tmp)
#
appNm <- paste0(dtstNm, " - Exp. map")
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  # Dummy, hidden div to load arrow-down icon
  tags$div(
    class = "container",
    style = "display: none;",
    tags$div(
      style = "margin-top: 50xp;",
      actionButton(
        "add_thing",
        label = "do it",
        class = "btn-success",
        icon = icon("arrow-down")
      )
    )
  ),
  #
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", ExpMapNm),
             appNm),
  h2(dtstNm), 
  h3("Define the relationship between biological samples and experimental factors."),
  br(),
  h4(em(tags$div("\"Parent sample\" is equivalent with \"Experiments\" in the sense of MaxQuant's \"Raw data\" or FragPipe's \"Workflow\" tabs.", tags$br(),
                 "Note that this is different from the meaning of \"Experiment\" we use here (that of a group of samples to process together and compare to each other)."))),
  br(),
  actionButton("saveBtn", "Save"),
  withSpinner(DT::DTOutput("ExpTbl", width = wTest2))
)
if (exists("ExpData3")) { rm(ExpData3) }
server <- function(input, output, session) {
  # Initialize output table
  ExpData3 <- ExpData # Output table
  #
  # Render dummy table
  output$ExpTbl <- DT::renderDT({ ExpData2 },
    FALSE,
    escape = FALSE,
    selection = "none",
    rownames = FALSE,
    editable = edith,
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
  #
  # Fill down
  sapply(1:L, function(x) {
    id1 <- ALLIDS[x]
    id2 <- paste0(ALLIDS[x], "___FD")
    tmp <- unlist(strsplit(id1, "___"))
    fct <- tmp[[1]]
    i <- as.integer(tmp[[2]])
    observeEvent(input[[id2]],
                 {
                   if (i < nr) {
                     x <- input[[id1]]
                     tp <- typeof(facLevels2[[fct]])
                     if (typeof(x) == tp) { x <- get(paste0("as.", tp))(x) }
                     for (k in (i+1):nr) {
                       idK <- paste0(fct, "___", as.character(k))
                       if (fct == "Use") {
                         updateCheckboxInput(session, idK, NULL, x)
                       } else {
                         updateSelectInput(session, idK, NULL, facLevels2[[fct]], x)
                       }
                     }
                   }
                 })
  })
  # Incremental fill-down for replicates
  sapply(rws, function(i) {
    if (i < nr) {
      id1 <- paste0("Replicate___", i)
      id2 <- paste0("Replicate___", i, "___INCR")
      observeEvent(input[[id2]],
                   {
                     if (i < nr) {
                       x <- input[[id1]]
                       tp <- typeof(facLevels2[[fct]])
                       if (typeof(x) == tp) { x <- get(paste0("as.", tp))(x) }
                       rplRg <- (i+1):nr
                       l <- length(rplRg)
                       m <- match(x, dflt_Rpl)+1
                       rplVal <- dflt_Rpl[m:(m+l-1)]
                       for (k in rplRg) {
                         y <- rplVal[k-i]
                         idK <- paste0("Replicate___", as.character(k))
                         updateSelectInput(session, idK, NULL, facLevels2[["Replicate"]], y)
                       }
                     }
                   }
      )
    }
  })
  # Manual cell edit (sample names)
  observeEvent(input$ExpTbl_cell_edit, {
    ExpData3[input$ExpTbl_cell_edit$row,
            input$ExpTbl_cell_edit$col+1] <- input$ExpTbl_cell_edit$value
  })
  # Save
  observeEvent(input$saveBtn, {
    for (fct in c(Fact2, "Use")) {
      ExpData3[[fct]] <- sapply(rws, function(i) {
        input[[paste0(fct, "___", i)]]
      })
      if (fct == "Use") {
        ExpData3[[fct]] <- as.logical(ExpData3[[fct]])
      } else {
        typ <- typeof(FactorsLevels[[fct]])
        ExpData3[[fct]] <- get(paste0("as.", typ))(ExpData3[[fct]])
        # Consider here detecting if a factor needs conversion to numeric/integer...
      }
    }
    ExpData3$"Sample name" <- ExpData2$"Sample name"
    assign("ExpData3", ExpData3, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("ExpData3"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
#
Exp.map %<o% ExpData3
k0 <- colnames(ExpMap)
k0 <- k0[which(!k0 %in% colnames(Exp.map))]
if (length(k0)) {
  Exp.map[, k0] <- Exp.map[match(Exp.map$MQ.Exp, Exp.map$MQ.Exp), k0]
}
tmpTbl <- Exp.map
tst <- lapply(colnames(tmpTbl), function(x) { typeof(tmpTbl[[x]]) })
w <- which(tst == "list")
if (length(w)) { for (i in w) { tmpTbl[[i]] <- sapply(tmpTbl[[i]], paste, collapse = ";") }}
tst <- try(write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE), silent = TRUE)
if ("try-error" %in% class(tst)) {
  dlg_message(paste0("File \"", ExpMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE)
}
#
#system(paste0("open \"", wd, "/", ExpMapNm, ".csv\""))
#Exp.map2 %<o% read.csv(ExpMapPath, check.names = FALSE)
Exp.map$Use <- as.logical(Exp.map$Use)
Exp.map$Use[which(is.na(Exp.map$Use))] <- FALSE
Exp.map <- Exp.map[which(Exp.map$Use),]
kol <- colnames(Exp.map)
kol1 <- c("MQ.Exp", "Sample name", "Use")
kol2 <- c("Reference", "Fractions", "Proportion", Factors)
kol2 <- kol2[which(kol2 %in% colnames(Exp.map))]
Exp.map <- aggregate(Exp.map[, kol1],
                     lapply(kol2, function(k) {
                       x <- Exp.map[[k]]
                       w <- which(is.na(x))
                       if (length(w)) { x[w] <- "NA" }
                       return(x)
                     }), function(x) { paste(unique(x), collapse = ";") })
colnames(Exp.map) <- c(kol2, kol1)
for (k in kol2) {
  x <- Exp.map[[k]]
  w <- which(x == "NA")
  if (length(w)) { Exp.map[w, k] <- NA }
}
if (!"list" %in% class(Exp.map$MQ.Exp)) { Exp.map$MQ.Exp <- strsplit(Exp.map$MQ.Exp, ";") }
tst <- setNames(lapply(Factors, function(Fact) {
  magrittr::set_colnames(aggregate(Exp.map[[Fact]], list(Exp.map[[Fact]]), length), c("Level", "Count"))
}), Factors)
msg2 <- msg <- paste(c("Check the number of samples per factor level below. Is everything ok? If not, click \"no\" to get back to editing the table.\n\n   -----\n", unlist(lapply(names(tst), function(nm) {
  c(paste0("-> ", nm),
    paste0("\t", c(paste(colnames(tst[[nm]]), collapse = "\t"),
                   apply(tst[[nm]], 1, paste, collapse = "\t"))),
    c("\n   -----\n"))
})), "\n"), collapse = "\n")
if (nchar(msg) > 1000) {
  cat(msg)
  msg2 <- paste0(substr(msg, 1, 996), "...")
}
tstXpMp <- c(TRUE, FALSE)[match(dlg_message(msg2, "yesno", rstudio = FALSE)$res, c("yes", "no"))]
