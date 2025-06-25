# Initial values for outputs
tmpDat2 <- NA
wAG2 <- wAG1
Outcome <- TRUE
txt2 <- ""
#
if (length(Iso) <= 1) {
  if (length(Iso) == 1) {
    warning("Skipping IRS normalisation: there is only 1 Isobarically labelled sample (aka \"Isobaric Set\")")
  } else { stop("\"LabelType\" is \"Isobaric\" but \"Iso\" has invalid length, investigate!") }
} else {
  irsDr <- paste0(nrmDr, "/Step ", nrmStp, " - IRS normalisation")
  if (!dir.exists(irsDr)) { dir.create(irsDr, recursive = TRUE) }
  dirlist <- unique(c(dirlist, irsDr))
  #
  currSamples <- allSamples[which(allSamples %in% colnames(tmpDat1))]
  m <- match(currSamples, Exp.map$Ref.Sample.Aggregate)
  design <- data.frame(Sample = currSamples,
                       Iso = Exp.map$Isobaric.set[m],
                       Samples_group = Exp.map[m, VPAL$column])
  design$Original.order <- 1:nrow(design)
  design <- design[order(design$Samples_group),]
  design <- design[order(design$Iso),]
  design$New.order <- 1:nrow(design)
  #
  # Define reference channels
  #Exp.map$Isobaric.set <- c(1, 1, 2, 1, 1, 2)
  IsoMapNm %<o% "Isobaric sets map"
  IsoMapPath %<o% paste0(wd, "/", IsoMapNm, ".csv")
  if (file.exists("IsoMap")) {
    reLoad <- TRUE
    if (exists("IsoMap")) {
      reLoad <- c(TRUE, FALSE)[match(dlg_message("Reload Isobaric sets map from disk?\n(you will still be able to edit it)", "yesno")$res, c("yes", "no"))]
    }
    if (reLoad) {
      IsoMap <- read.csv(IsoMapPath, check.names = FALSE)
      IsoMap$`All channels` <- strsplit(IsoMap$`All channels`, ";")
      IsoMap$`Reference channel(s)` <- strsplit(IsoMap$`Reference channel(s)`, ";")
      # Column "Set" MUST be absent AND valid
      tst1 <- TRUE
      if ("Set" %in% colnames(IsoMap)) {
        tst1 <- sum(!unique(Exp.map$Isobaric.set[which(Exp.map$Use)]) %in% IsoMap$Set)
      }
      # Column "All channels" MUST be present and valid
      tst2 <- sum((!"All channels" %in% colnames(IsoMap))&&
                    (sum(sapply(IsoMap$"All channels", function(x) {
                      !x %in% unique(unlist(Exp.map$"Isobaric label details"))
                    }))))
      # IF column "Reference channel(s)" is present, it MUST be valid
      tst3 <- FALSE
      if ("Reference channel(s)" %in% colnames(IsoMap)) {
        tst3 <- aggregate(Exp.map$"Isobaric label details", list(Exp.map$Isobaric.set), list)
        tst3 <- sum(vapply(IsoMap$Set, function(i) {
          x <- IsoMap$"Reference channel(s)"[[i]]
          x <- x[which(x != "")]
          sum(!x %in% tst3$x[[i]])
        }, 1) > 0)
      }
      if (sum(tst1, tst2, tst3)) {
        warning("Invalid Isobaric sets map reloaded, ignoring...")
        rm(IsoMap)
      }
    }
  }
  if (!exists("IsoMap")) {
    IsoMap <- aggregate(Exp.map$"Isobaric label details", list(Exp.map$Isobaric.set), list)
    colnames(IsoMap) <- c("Set", "All channels")
    IsoMap$"Reference channel(s)" <- lapply(1:nrow(IsoMap), function(x) { })
  }
  rws <- 1:nrow(IsoMap)
  IsoMap2 <- IsoMap[, c("Set", "Reference channel(s)")]
  IsoMap2$"Reference channel(s)" <- vapply(rws, function(i) {
    as.character(shiny::selectInput(paste0("ID___", as.character(i)),
                                    "",
                                    unlist(IsoMap$`All channels`[i]),
                                    unlist(IsoMap$`Reference channel(s)`[i])[1],
                                    TRUE,
                                    FALSE,
                                    width = "300px"))
  }, "a")
  wTest1 <- setNames(vapply(colnames(IsoMap2), function(k) { #k <- colnames(IsoMap2)[1]
    x <- max(c(nchar(k),
               nchar(as.character(IsoMap[[k]])) + 6), na.rm = TRUE)*10
    if (is.na(x)) { x <- 15 } else { x <- max(c(ceiling(x/10)*10, 30)) }
    return(x)
  }, 1), colnames(IsoMap2))
  wTest2 <- sum(wTest1)
  wTest1 <- paste0(as.character(wTest1), "px")
  wTest1 <- aggregate((1:length(wTest1))-1, list(wTest1), c)
  wTest1 <- apply(wTest1, 1, function(x) {
    x2 <- as.integer(x[[2]])
    list(width = x[[1]],
         targets = x2,
         names = colnames(IsoMap2)[x2+1])
  })
  appNm <- "IRS reference channels"
  if (exists("IHAVERUN")) { rm(IHAVERUN) }
  ui <- fluidPage(
    useShinyjs(),
    setBackgroundColor( # Doesn't work
      color = c(#"#F8F8FF",
        "#F2E4F5"),
      gradient = "linear",
      direction = "bottom"
    ),
    extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
    tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
    titlePanel(tag("u", appNm),
               appNm),
    br(),
    fluidRow(column(5,
                    actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                    h5(HTML("Select 1 reference channel per Isobaric set for Internal Reference Scaling normalisation")),
                    h5(HTML("We define an \"Isobaric set\" as a combined sample made up of pooled samples individually labelled isobarically.")),
                    h5(HTML("If no reference channel was included, select all (or any you want to average as your reference for that set) to simulate a reference!")),
                    withSpinner(DT::DTOutput("IsoMap", width = wTest2)))),
    br(),
    br()
  )
  server <- function(input, output, session) {
    # Initialize output table
    IsoMap3 <- IsoMap # Output table
    #
    # Render dummy table
    output$IsoMap <- DT::renderDT({ IsoMap2 },
                                  FALSE,
                                  escape = FALSE,
                                  selection = "none",
                                  editable = FALSE,
                                  rownames = FALSE,
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
    # Save
    observeEvent(input$saveBtn, {
      IsoMap3$"Reference channel(s)" <- sapply(rws, function(i) {
        id <- paste0("ID___", as.character(i))
        return(input[[id]])
      })
      assign("IsoMap3", IsoMap3, envir = .GlobalEnv)
      assign("IHAVERUN", TRUE, .GlobalEnv)
      stopApp()
    })
    #observeEvent(input$cancel, { stopApp() })
    session$onSessionEnded(function() { stopApp() })
  }
  runKount <- 0
  while ((!runKount)||(!exists("IHAVERUN"))) {
    eval(parse(text = runApp), envir = .GlobalEnv)
    runKount <- runKount+1
  }
  IsoMap %<o% IsoMap3
  tmpTbl <- IsoMap
  tst <- lapply(colnames(tmpTbl), function(x) { typeof(tmpTbl[[x]]) })
  w <- which(tst == "list")
  if (length(w)) { for (i in w) { tmpTbl[[i]] <- sapply(tmpTbl[[i]], paste, collapse = ";") }}
  tst <- try(write.csv(tmpTbl, file = IsoMapPath, row.names = FALSE), silent = TRUE)
  if ("try-error" %in% class(tst)) {
    dlg_message(paste0("File \"", IsoMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
    write.csv(tmpTbl, file = IsoMapPath, row.names = FALSE)
  }
  #
  # Temporarily impute (unavoidable here if we want to show PCAs)
  ImpGrps <- Exp.map[match(currSamples, Exp.map$Ref.Sample.Aggregate),
                     VPAL$column]
  tmp <- proteoCraft::Data_Impute2(tmpDat1[, currSamples], ImpGrps)
  tmpDat1Imp <- tmp$Imputed_data
  Pos <- tmp$Positions_Imputed
  #
  # A more complex description of IRS can be found here
  # https://pwilmart.github.io/IRS_normalization/understanding_IRS.html
  # However there is a mistake in the code there, as the geometric row means (instead of row sums) should be used, as for the average
  # Here the code is also adapted for log scale
  setNms <- paste0("Set", as.character(IsoMap$Set))
  # All samples per set
  mixed <- unique(unlist(lapply(Factors, function(x) { which(Exp.map[[x]] == "Mixed_IRS") })))
  mixedSets <- paste0("Set", as.character(unique(Exp.map$Isobaric.set[mixed])))
  irsSamples <- setNames(lapply(1:nrow(IsoMap), function(i) {
    Exp.map$Ref.Sample.Aggregate[which((! 1:nrow(Exp.map) %in% mixed)&(Exp.map$Isobaric.set == IsoMap$Set[i]))]
  }), setNms)
  # Either intensity from IRS channel, or failing that log10 row means for the group
  irsRowMeans <- setNames(lapply(setNms, function(nm) {
    if (nm %in% mixedSets) {
      x <- tmpDat1Imp[, Exp.map$Ref.Sample.Aggregate[mixed[which(mixedSets == nm)]], drop = FALSE]
    } else {
      x <- tmpDat1Imp[, irsSamples[[nm]]]
    }
    #log10(parApply(parClust, tmpDat1Imp[, x], 1, function(y) { sum(10^proteoCraft::is.all.good(y)) }))
    #parApply(parClust, tmpDat1Imp[, x], 1, function(y) { sum(proteoCraft::is.all.good(y)) })
    x <- parApply(parClust, x, 1, function(y) { mean(proteoCraft::is.all.good(y)) })
    return(x)
  }), setNms)
  #View(do.call(cbind, irsRowMeans))
  # log10 geometric mean
  allMeans <- parApply(parClust, tmpDat1Imp[, currSamples], 1, function(x) {
    #log10(mean(10^proteoCraft::is.all.good(x)))
    mean(proteoCraft::is.all.good(x))
  })
  #View(data.frame(meanOfAll = allMeans))
  # log10 scaling factors:
  irsFact <- setNames(lapply(setNms, function(x) {
    allMeans - irsRowMeans[[x]]
  }), setNms)
  #View(do.call(cbind, irsFact))
  # log10 new data
  tmpDat2Imp <- lapply(setNms, function(x) { #x <- setNms[1]
    sweep(tmpDat1Imp[, irsSamples[[x]]], 1, irsFact[[x]], "+")
  })
  tmpDat2Imp <- do.call(cbind, c(tmpDat2Imp))
  #View(tmpDat2Imp)
  # Remove imputed values
  tmpDat2 <- tmpDat2Imp
  w <- which(!currSamples %in% colnames(tmpDat2))
  tmpDat2[, currSamples[w]] <- NA
  wImp <- which(Pos, arr.ind = TRUE)
  if (nrow(wImp)) { tmpDat2[, currSamples][wImp] <- tmpDat1[, currSamples][wImp] }
  tmpDat2 <- tmpDat2[, colnames(tmpDat2Imp)]
  #
  wAG2 <- wAG1
  #wAGstrict <- which(parApply(parClust, tmpDat2, 1, function(x) { length(proteoCraft::is.all.good(x)) }) == length(currSamples))
  # Plot
  map <- Exp.map[, c("Ref.Sample.Aggregate", "Isobaric.set")]
  w1 <- which(colnames(tmpDat1Imp) %in% map$Ref.Sample.Aggregate) # Should be all
  PCAlyLst <- scoresLst <- PCsLst <- list()
  prev <- "before"
  tmp <- pcaBatchPlots(tmpDat1Imp[, currSamples],
                       prev,
                       "Isobaric.set",
                       map = Exp.map,
                       SamplesCol = "Ref.Sample.Aggregate",
                       intRoot = "",
                       dir = irsDr,
                       ttl = "PCA plot - IRS batch corr.",
                       isRef = !currSamples %in% unlist(irsSamples))
  PCAlyLst[[prev]] <- tmp$PlotLy
  scoresLst[[prev]] <- tmp$Scores
  PCsLst[[prev]] <- tmp$PCs
  curr <- "IRS"
  w <- which(currSamples %in% colnames(tmpDat2Imp))
  tmp <- pcaBatchPlots(tmpDat2Imp[, currSamples[w]],
                       curr,
                       "Isobaric.set",
                       map = Exp.map,
                       SamplesCol = "Ref.Sample.Aggregate",
                       intRoot = "",
                       dir = irsDr,
                       ttl = "PCA plot - IRS batch corr.",
                       isRef = !currSamples[w] %in% unlist(irsSamples))
  PCAlyLst[[curr]] <- tmp$PlotLy
  scoresLst[[curr]] <- tmp$Scores
  PCsLst[[curr]] <- tmp$PCs
  #
  appNm <- "IRS batch correction"
  msg <- "Keep results from IRS batch correction? (untick to cancel correction)"
  if ((!exists("KeepIRSRes"))||(length(KeepIRSRes) != 1)||(!is.logical(KeepIRSRes))||(is.na(KeepIRSRes))) {
    KeepIRSRes <- TRUE
  }
  if (exists("IHAVERUN")) { rm(IHAVERUN) }
  ui <- fluidPage(
    useShinyjs(),
    setBackgroundColor( # Doesn't work
      color = c(#"#F8F8FF",
        "#EDE1F0"),
      gradient = "linear",
      direction = "bottom"
    ),
    extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
    tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
    titlePanel(tag("u", appNm),
               appNm),
    br(),
    fluidRow(column(12,
                    checkboxInput("KeepIRSRes", msg, KeepIRSRes),
                    actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                    h4("Recommended criteria:"),
                    h5(HTML(paste0("&nbsp;Does the original grouping follow known ", IsobarLab, " sets (it usually does)?"))),
                    h5(HTML(paste0("&nbsp;&nbsp;-> If yes: accept the correction if it removes the original grouping by ", IsobarLab, " sets."))),
                    h5(HTML(paste0("&nbsp;&nbsp;-> If no: you can reject the correction, congratulations: you were lucky, this time there was no ",
                                   IsobarLab, " set-related batch effect.")))),
             br(),
             fluidRow(column(5, withSpinner(plotlyOutput("Before", height = "600px"))),
                      column(5, withSpinner(plotlyOutput("After", height = "600px"))))),
    br(),
    br()
  )
  server <- function(input, output, session) {
    output$Before <- renderPlotly(PCAlyLst[[prev]]$Isobaric.set)
    output$After <- renderPlotly(PCAlyLst[[curr]]$Isobaric.set)
    observeEvent(input[["KeepIRSRes"]], {
      assign("KeepIRSRes", as.logical(input[["KeepIRSRes"]]), envir = .GlobalEnv)
    })
    observeEvent(input$saveBtn, {
      assign("KeepIRSRes", as.logical(input[["KeepIRSRes"]]), envir = .GlobalEnv)
      assign("IHAVERUN", TRUE, .GlobalEnv)
      stopApp()
    })
    #observeEvent(input$cancel, { stopApp() })
    session$onSessionEnded(function() { stopApp() })
  }
  runKount <- 0
  while ((!runKount)||(!exists("IHAVERUN"))) {
    eval(parse(text = runApp), envir = .GlobalEnv)
    runKount <- runKount+1
  }
  msg <- paste0(" -> Internal Reference Scaling-based correction of ", IsobarLab, "-associated batch effect ", c("rejec", "accep")[KeepIRSRes+1], "ted.\n")
  if (KeepIRSRes) {
    txt2 <- paste0("corrected against the ", IsobarLab, "-associated batch effect using the Internal Reference Scaling method (P. Wilmarth)")
  }
  Outcome <- KeepComBatRes
  cat(msg)
}
  