# A sourced script to edit the MS files map
#
library(shiny)
library(shinyjs)
library(DT)
#
labelMode <- match(LabelType, c("LFQ", "Isobaric"))
# Important distinction here:
# - If LabelType == "LFQ", a priori "Parent sample" == MQ.Exp, but this can be changed!
# - If LabelType == "Isobaric", at this stage each raw file maps to a single MQ.Exp but normally to several "Parent sample" values!
#   Thus, in that case we use MQ.Exp as the basis here and do not have a "Parent sample" column!!!
#
if ((LabelType == "LFQ")&&(!"Parent sample" %in% colnames(FracMap))&&("MQ.Exp" %in% colnames(FracMap))) {
  FracMap$"Parent sample" <- FracMap$MQ.Exp
}
if (LabelType == "Isobaric") {
  if ("Isobaric set" %in% colnames(FracMap)) {
    if (!"Isobaric.set" %in% colnames(FracMap)) {
      FracMap$Isobaric.set <- FracMap$"Isobaric set"
    }
    FracMap$"Isobaric set" <- NULL
  }
  if (!"Isobaric.set" %in% colnames(FracMap)) {
    FracMap$"Isobaric.set" <- 1
  }
}
# if (("MQ.Exp" %in% colnames(FracMap))&&("Experiment" %in% colnames(ev))) {
#   FracMap$Exp_Backup <- FracMap$MQ.Exp
# }
if ("Use" %in% colnames(FracMap)) {
  FracMap$Use <- as.logical(FracMap$Use)
  FracMap$Use[which(is.na(FracMap$Use))] <- FALSE
} else { FracMap$Use <- TRUE }
if ("PTM-enriched" %in% colnames(FracMap)) {
  FracMap$"PTM-enriched"[which(!FracMap$"PTM-enriched" %in% Modifs$`Full name`)] <- NA
} else { FracMap$"PTM-enriched" <- NA }
FracMap$Use <- as.logical(toupper(FracMap$Use))
FracMap$Use[which(is.na(FracMap$Use))] <- TRUE
nr <- nrow(FracMap)
rws <- seq_len(nr)
#
parKol <- c("Parent sample", "MQ.Exp")[labelMode]
parKol2 <- c("Sample", "Experiment")[labelMode]
kol1 <- c(parKol2, "Fraction", "Use", "PTMenriched")
kol0 <- c(parKol, "Fraction", "Use", "PTM-enriched")
if (LabelType == "Isobaric") {
  kol1 <- c(kol1, "IsobaricSet")
  kol0 <- c(kol0, "Isobaric.set")
}
ALLIDS <- unlist(lapply(kol1, function(x) { paste0(x, "___", rws) }))
allPTMs <- unique(c(Modifs$`Full name`, NA))
#
# Dummy for shiny app
FracMap2a <- FracMap
if ((LabelType == "LFQ")&&("MQ.Exp" %in% colnames(FracMap2a))) { FracMap2a$MQ.Exp <- NULL }
colnames(FracMap2a) <- gsub("-", "", gsub("\\.", " ", colnames(FracMap2a)))
if (LabelType == "Isobaric") {
  # Exceptions
  colnames(FracMap2a)[which(colnames(FracMap2a) == "Isobaric set")] <- "Isobaric.set"
  colnames(FracMap2a)[which(colnames(FracMap2a) == "MQ Exp")] <- "MQ.Exp"
}
lu <- length(unique(FracMap2a$"Raw files name"))
if (lu == nr) {
  FracMap2a$"Raw file" <- NULL
} else {
  tst <- gsub("/[^/]+$", "", FracMap2a$"Raw file")
  nc <- nchar(tst)
  mnc <- min(nc)
  if (mnc >= 5) {
    w <- max(which(sapply(1:mnc, function(x) { length(unique(substr(tst, 1, x))) }) == 1))
    FracMap2a$"Raw file" <- paste0("...", substr(FracMap2a$"Raw file", w+1, nchar(FracMap2a$"Raw file")))
  }
  FracMap2a$"Raw files name" <- NULL
}
# Original table width
wTest0 <- setNames(sapply(colnames(FracMap2a), function(k) { #k <- colnames(FracMap2a)[1]
  tmp <- FracMap2a[[k]]
  if ("logical" %in% class(tmp)) { tmp <- as.integer(tmp) }
  tmp <- as.character(tmp)
  x <- max(nchar(c(k, tmp)) + 3, na.rm = TRUE)
  if (k == "PTMenriched") { x <- max(c(x, nchar(allPTMs) + 3), na.rm = TRUE) }
  x <- x*10
  if (is.na(x)) { x <- 15 } else { x <- max(c(ceiling(x/10)*10, 30)) }
  return(x)
}), colnames(FracMap2a))
FracMap2 <- FracMap2a[, which(colnames(FracMap2a) != "Exp_Backup")]
FracMap2$Use <- shinyCheckInput(FracMap2a$Use,
                                "Use")
FracMap2[[parKol2]] <- shinyTextInput(FracMap[[parKol]],
                                      parKol2,
                                      paste0(wTest0[parKol], "px"))
FracMap2$Fraction <- shinyNumInput(FracMap2a$Fraction,
                                   1,
                                   Inf,
                                   1,
                                   1,
                                   paste0(wTest0["Fraction"], "px"),
                                   "Fraction")
FracMap2$"PTMenriched" <- shinySelectInput(FracMap$"PTM-enriched",
                                         "PTMenriched",
                                         allPTMs,
                                         paste0(wTest0["PTMenriched"], "px"))
if (LabelType == "Isobaric") {
  FracMap2$"IsobaricSet" <- shinyNumInput(FracMap2$"Isobaric.set",
                                          1,
                                          Inf,
                                          1,
                                          1,
                                          paste0(wTest0["Isobaric.set"], "px"),
                                          "IsobaricSet")
  FracMap2$"Isobaric.set" <- NULL
}
kol <- colnames(FracMap2)[which(!colnames(FracMap2) %in% c(kol0, kol1))]
#kol %in% names(wTest0)
#wTest0[kol]
ALLFDIDS <- c()
for (k0 in kol0) { #k0 <- parKol #k0 <-  "Isobaric.set"
  k1 <- kol1[match(k0, kol0)]
  fdNm <- paste0(k1, "___FD")
  FracMap2[[fdNm]] <- shinyFDInput(k1, nr, TRUE)
  ALLFDIDS <- c(ALLFDIDS, paste0(k1, "___", rws, "___FD"))
  kol <- c(kol, k1, fdNm)
}
kol <- grep("^MQ[ \\.]Exp", kol, value = TRUE, invert = TRUE)
FracMap2 <- FracMap2[, kol]
wTest1 <- sapply(colnames(FracMap2), function(k1) { #k1 <- colnames(FracMap2)[1] #k1 <- "IsobaricSet"
  if (k1 %in% names(wTest0)) {
    x <- wTest0[k1]
  } else {
    k0 <- kol0[match(k1, kol1)]
    if ((!is.na(k0))&&(k0 %in% names(wTest0))) { # Should not happen, currently we use kol1 names for wTest1
      x <- wTest0[k0]
    } else {
      x <- 30
    }
  }
  if (is.na(x)) { x <- 30 }
  return(x)
})
wTest2 <- sum(wTest1)
wTest1 <- paste0(as.character(wTest1), "px")
wTest1 <- aggregate((1:length(wTest1))-1, list(wTest1), c)
wTest1 <- apply(wTest1, 1, function(x) {
  list(width = x[[1]],
       targets = x[[2]],
       names = colnames(FracMap2)[x[[2]]+1])
})
#
g <- grep("___((FD)|(INCR))$", colnames(FracMap2))
colnames(FracMap2)[g] <- ""
#
kN <- c("Raw files name", "Parameter group", "PTMs")
wN <- which(colnames(FracMap2) %in% kN) - 1
wY <- which(!colnames(FracMap2) %in% kN) - 1
edith <- list(target = "column",
              disable = list(columns = wN),
              enable = list(columns = wY))
#
if (exists("FracMap3")) { rm(FracMap3) }
appNm <- paste0(dtstNm, " - MS files map")
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#EFE6F5"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  # Dummy, hidden div to load arrow-down icon
  # tags$div(
  #   class = "container",
  #   style = "display: none;",
  #   tags$div(
  #     style = "margin-top: 50xp;",
  #     actionButton(
  #       "add_thing",
  #       label="do it",
  #       class="btn-success",
  #       icon=icon("arrow-down")
  #     )
  #   )
  # ),
  # #
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", FracMapNm),
             appNm),
  h2(dtstNm), 
  br(),
  h3(paste0("Enter details of MS raw files to \"", parKol2, "\" and \"Fraction\" mappings, then click \"Save\".")),
  br(),
  h4(em(c(paste0("NB: Here we use \"", parKol2,
                 "\" in the same sense as what MaxQuant's \"Raw data\" and FragPipe's \"Workflow\" tabs call \"Experiment\"\n(i.e. group of MS files, such as fractions, derived from a same parent biological sample)."),
          paste0("NB: Here \"", parKol2,
                 "\" does not mean parent biological sample but combined isobarically-labelled sample - we will deal with individual samples later"))[labelMode])),
  br(),
  mainPanel(#tabsetPanel(tabPanel(
    actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
    DT::DTOutput("FracTbl", width = paste0(wTest2, "px")),
    br(),
    width = 12
  ))
server <- function(input, output, session) {
  # Initialize output table
  FracMap3 <- FracMap
  #
  # Render dummy table
  output$FracTbl <- DT::renderDT({ FracMap2 },
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
      Shiny.bindAll(table.table().node());")
  )
  # Fill down
  sapply(ALLFDIDS, function(id2) {
    id1 <- gsub("___FD$", "", id2)
    #cat(id2, "\n")
    tmp <- unlist(strsplit(id1, "___"))
    k1 <- tmp[[1]]
    #k0 <- kol0[match(k1, kol1)]
    i <- as.integer(tmp[[2]])
    observeEvent(input[[id2]],
                 {
                   if (i < nr) {
                     x <- input[[id1]]
                     if (k1 == "Sample") { x <- as.character(x) }
                     if (k1 %in% c("Fraction", "IsobaricSet")) { x <- as.integer(x) }
                     if (k1 == "Use") { x <- as.logical(x) }
                     if (k1 == "ISo") { x <- as.integer(x) }
                     for (k in (i+1):nr) {
                       idK <- paste0(k1, "___", as.character(k))
                       if (k1 == "PTMenriched") { updateSelectInput(session, idK, NULL, allPTMs, x) }
                       if (k1 %in% c("Fraction", "IsobaricSet")) { updateNumericInput(session, idK, NULL, x, 1, Inf, 1) }
                       if (k1 == "Use") { updateCheckboxInput(session, idK, NULL, x) }
                       if (k1 == "Sample") { updateTextInput(session, idK, NULL, x) }
                     }
                   }
                 }, ignoreInit = TRUE)
  })
  # Manual cell edit (sample names)
  observeEvent(input$FracTbl_cell_edit, {
    k2 <- colnames(FracTbl_cell_edit$col+1)
    m3 <- match(k2, colnames(FracMap3))
    FracMap3[input$FracTbl_cell_edit$row, m3] <- input$FracTbl_cell_edit$value
  })
  # Save
  observeEvent(input$saveBtn, {
    for (k1 in kol1) {
      k0 <- kol0[match(k1, kol1)]
      FracMap3[[k0]] <- sapply(rws, function(i) { input[[paste0(k1, "___", as.character(i))]] })
    }
    assign("FracMap3", FracMap3, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("FracMap3"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
#
FracMap <- FracMap3
# Below: inactivated for now, only useful for reruns and may cause issues because MQ.Exp/"Parent sample" are ambiguous at this stage
# if ("Exp_Backup" %in% colnames(FracMap)) {
#   ev$Experiment <- FracMap$MQ.Exp[match(ev$Experiment, FracMap$Exp_Backup)]
# }
#
if (LabelType == "Isobaric") {
  FracMap$Isobaric.set <- as.integer(FracMap$Isobaric.set)
  Iso %<o% sort(unique(FracMap$Isobaric.set))
}
for (k in colnames(FracMap)) { FracMap[[k]] <- gsub("\t|\n|^ +| +$", "", FracMap[[k]]) } # Avoid issues with incorrect entries (hidden spaces, newlines...)
FracMap$Use <- as.logical(toupper(FracMap$Use))
tst <- try(write.csv(FracMap, file = FracMapPath, row.names = FALSE), silent = TRUE)
if ("try-error" %in% class(tst)) {
  dlg_message(paste0("File \"", FracMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  write.csv(FracMap, file = FracMapPath, row.names = FALSE)
}
if (LabelType == "LFQ") {
  colnames(FracMap)[which(colnames(FracMap) == "Parent sample")] <- "MQ.Exp"
}
FracMap <- FracMap[which(FracMap$Use),]
MQ.Exp %<o% sort(unique(FracMap$MQ.Exp))
#
kol <- setNames(c("Raw file", "Fraction", "PTM-enriched"), c("Raw files", "Fractions", "Type of PTM-enrichment sample"))
if (LabelType == "Isobaric") {
  kol["Isobaric.set"] <- "Isobaric.set"
}
kol <- kol[which(kol %in% colnames(FracMap))]
stopifnot(length(kol) > 0)
tst <- aggregate(FracMap[, kol], list(FracMap$MQ.Exp), function(x) { length(unique(x)) })
colnames(tst)[1] <- "Biological sample"
kount <- 0
strt <- unique(substr(tst$"Biological sample", 1, 1))
while (length(strt) == 1) {
  kount <- kount + 1
  tst$"Biological sample" <- substr(tst$"Biological sample", 2, nchar(tst$"Biological sample"))
  strt <- unique(substr(tst$"Biological sample", 1, 1))
}
nc <- nchar(c("Biological sample", tst[[1]]))
mx <- max(nc)
w <- which(nc < mx)-1
if (length(w)) {
  if (0 %in% w) { colnames(tst)[1] <- paste0(c(colnames(tst)[1], rep(" ", mx-nc[1])), collapse = "") }
  w <- w[which(w > 0)]
  tst[w, 1] <- sapply(w, function(x) { paste0(c(tst[x, 1], rep(" ", mx-nc[x+1])), collapse = "") })
}
msg <- c(paste(colnames(tst), collapse = "\t"), do.call(paste, c(tst, sep = "\t")))
msg2 <- msg <- paste(c(paste0("Check the number of:\n - MS files\n - fractions\n - enriched sample types",
                              c("", "\n - isobaric sets")[labelMode],
                              "\nper parent biological sample below. Is everything ok? If not, click \"no\" to get back to editing the table.\n\n   -----\n"),
                       msg, "\n   -----\n"), collapse = "\n")
if (nchar(msg) > 1000) {
  cat(msg)
  msg2 <- paste0(substr(msg, 1, 996), "...")
}
tstFrMp <- c(TRUE, FALSE)[match(dlg_message(msg2, "yesno", rstudio = FALSE)$res, c("yes", "no"))]
