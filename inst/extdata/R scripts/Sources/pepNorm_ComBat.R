# These normalisations are more like corrections: they will change the shape of the vector
#
# Initial values for outputs
tmpDat2 <- NA
wAG2 <- wAG1
Outcome <- TRUE
txt2 <- ""
#
pack <- "sva"
bioc_req <- unique(c(bioc_req, pack))
if (!require(pack, character.only = TRUE, quietly = TRUE)) { pak::pkg_install(p) }
require(pack, character.only = TRUE)
#
btchDir <- paste0(nrmDr, "/Step ", nrmStp, " - Batch corr.")
if (!dir.exists(btchDir)) { dir.create(btchDir, recursive = TRUE) }
dirlist <- unique(c(dirlist, btchDir))
#
# Let's do some imputation:
# (remember to remove those afterwards!)
currSamples <- allSamples[which(allSamples %in% colnames(tmpDat1))]
ImpGrps <- Exp.map[match(currSamples, Exp.map$Ref.Sample.Aggregate),
                   VPAL$column]
tmpDat2 <- tmpDat1[, currSamples]*NA
#
myBatch3 <- myBatch2 <- myBatch <- normSequence[[nrmStp]]$Batch #myBatch <- "Lit"
if (length(myBatch) > 1) {
  myBatch2 <- paste(substr(myBatch, 1, 3), collapse = "")
  myBatch3 <- paste(myBatch, collapse = "/")
}
mod0a <- model.matrix(~1, data = Exp.map[match(currSamples, Exp.map$Ref.Sample.Aggregate),])
#
# Impute
tmp <- proteoCraft::Data_Impute2(tmpDat1[, currSamples], ImpGrps)
tmpDat2Imp <- tmpDat1Imp <- tmp$Imputed_data
Pos <- tmp$Positions_Imputed
#
for (lGrp in NormGrps$Group) { #lGrp <- NormGrps$Group[1] # Longitudinal group (peptide class)
  grpMtch <- match(NormGrps$IDs[[match(lGrp, NormGrps$Group)]],
                   tmpDat1$id[wAG1])
  grpMtch <- grpMtch[which(!is.na(grpMtch))]
  #
  # For ComBat we only use the longitudinal groups (peptide normalisation group),
  # not the transversal groups (comparison/ratio groups).
  # Indeed, batches will often intersect with the latter
  btchs <- Exp.map[match(currSamples, Exp.map$Ref.Sample.Aggregate), myBatch2]
  tmpDat2Imp[grpMtch, currSamples] <- ComBat(dat = tmpDat1Imp[grpMtch,],
                                            batch = btchs,
                                            mod = mod0a,
                                            par.prior = TRUE)
}
tmpDat2 <- tmpDat2Imp[, currSamples]
#
wAG2 <- wAG1
# Remove imputed data:
w <- which(Pos, arr.ind = TRUE)
if (nrow(w)) { tmpDat2[, currSamples][w] <- tmpDat1[, currSamples][w] }
#
# Plot
PCAlyLst <- scoresLst <- PCsLst <- list()
tmp <- pcaBatchPlots(tmpDat1Imp[wAG1, currSamples],
                     "original",
                     myBatch2,
                     Exp.map,
                     intRoot = "")
PCAlyLst[["original"]] <- tmp$PlotLy
scoresLst[["original"]] <- tmp$Scores
PCsLst[["original"]] <- tmp$PCs
tmp <- pcaBatchPlots(tmpDat2Imp[wAG1, currSamples],
                     "original",
                     myBatch2,
                     Exp.map,
                     intRoot = "")
PCAlyLst[[myBatch2]] <- tmp$PlotLy
scoresLst[[myBatch2]] <- tmp$Scores
PCsLst[[myBatch2]] <- tmp$PCs
#
appNm <- paste0("Batch corr.: original -> ", myBatch3)
msg <- "Keep results from ComBat batch correction? (untick to cancel correction)"
if ((!exists("KeepComBatRes"))||(length(KeepComBatRes) != 1)||(!is.logical(KeepComBatRes))||(is.na(KeepComBatRes))) {
  KeepComBatRes <- TRUE
}
PCs <- data.frame("Component" = paste0("PC", as.character(1:length(PCsLst[["original"]]$sdev))),
                  "Before (%)" = round(100*(PCsLst[["original"]]$sdev)^2 / sum(PCsLst[["original"]]$sdev^2), 0),
                  "After (%)" = round(100*(PCsLst[[myBatch2]]$sdev)^2 / sum(PCsLst[[myBatch2]]$sdev^2), 0))
if (exists("IHAVERUN")) { rm(IHAVERUN) }
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#EEFAE6"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", appNm),
             appNm),
  br(),
  fluidRow(column(5,
                  checkboxInput("KeepResults", msg, KeepComBatRes),
                  actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                  h4("Recommended criteria:"),
                  h5(HTML("&nbsp;Does the original grouping follow known batches?")),
                  h5(HTML("&nbsp;&nbsp;-> If no: only accept the correction if it improves the apparent grouping of expected sample groups.")),
                  h5(HTML("&nbsp;&nbsp;-> If yes: accept the correction if...")),
                  h5(HTML("&nbsp;&nbsp;&nbsp;- ... it removes the original grouping by batches...")),
                  h5(HTML("&nbsp;&nbsp;&nbsp;- ... or it improves the apparent grouping of expected sample groups.")),
                  withSpinner(DTOutput("PCs")))),
  br(),
  fluidRow(column(6, withSpinner(plotlyOutput("Before", height = "600px"))),
           column(6, withSpinner(plotlyOutput("After", height = "600px")))),
  br(),
  br()
)
server <- function(input, output, session) {
  output$Before <- renderPlotly(PCAlyLst[["original"]][[myBatch2]])
  output$After <- renderPlotly(PCAlyLst[[myBatch2]][[myBatch2]])
  output$PCs <- renderDT({ PCs },
                         FALSE,
                         escape = FALSE,
                         selection = "none",
                         editable = FALSE,
                         rownames = FALSE,
                         options = list(
                           dom = 't',
                           paging = FALSE,
                           ordering = FALSE
                         ),
                         callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
  observeEvent(input[["KeepResults"]], {
    assign("KeepComBatRes", as.logical(input[["KeepResults"]]), envir = .GlobalEnv)
  })
  observeEvent(input$saveBtn, {
    assign("KeepComBatRes", as.logical(input[["KeepResults"]]), envir = .GlobalEnv)
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
msg <- paste0(" -> correction of ", myBatch3, "-associated batch effect ", c("rejec", "accep")[KeepComBatRes+1], "ted.\n")
if (KeepComBatRes) {
  txt2 <- paste0("corrected against the ", myBatch3, "-associated batch effect")
}
Outcome <- KeepComBatRes
cat(msg)
