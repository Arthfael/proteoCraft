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
if (!require(pack, character.only = TRUE, quietly = TRUE)) { pak::pak(p) }
require(pack, character.only = TRUE)
#
btchDir <- paste0(nrmDr, "/Step ", nrmStp, " - Batch corr.")
if (!dir.exists(btchDir)) { dir.create(btchDir, recursive = TRUE) }
dirlist <- unique(c(dirlist, btchDir))
#
# Let's do some imputation:
# (remember to remove those afterwards!)
currSamples <- allSamples[which(allSamples %in% colnames(tmpDat1))]
smplsMtch <- match(currSamples, Exp.map$Ref.Sample.Aggregate)
ImpGrps <- Exp.map[smplsMtch,
                   VPAL$column]
tmpDat2 <- tmpDat1[, currSamples]*NA
#
myBatch3 <- myBatch2 <- myBatch <- normSequence[[nrmStp]]$Batch #myBatch <- "Lit"
if (length(myBatch) > 1L) {
  myBatch2 <- paste(substr(myBatch, 1L, 3L), collapse = "")
  myBatch3 <- paste(myBatch, collapse = "/")
}
m <- match(gsub("___", "_", as.character(expMap[smplsMtch, RSA$limmaCol])), row.names(designMatr))
mdlMtr <- designMatr_noBatch[m,] # Note: it is important for good batch correction to also make ComBat aware of other covariates so it doesn't remove them too!
# Hence why the "~1" model should be avoided. We now pre-make the ComBat model matrix in advance. Important: it cannot be based on a formula with ~ 0  intercept, otherwise ComBat will fail!
#mdlMtr <- model.matrix(~1, data = Exp.map[match(currSamples, Exp.map$Ref.Sample.Aggregate),]) # Old code for reference, do NOT use!
#
# Impute
tmp <- Data_Impute2(tmpDat1[, currSamples], ImpGrps)
tmpDat2Imp <- tmpDat1Imp <- tmp$Imputed_data
Pos <- tmp$Positions_Imputed
#
for (lGrp in NormGrps$Group) { #lGrp <- NormGrps$Group[1L] # Longitudinal group (peptide class)
  w <- which(tmpDat1$id[wAG1] %in% NormGrps$IDs[[match(lGrp, NormGrps$Group)]])
  if (length(w)) {
    # For ComBat we only use longitudinal groups (peptide normalisation group),
    # not transversal groups (comparison/ratio groups).
    # Indeed, batches will often intersect with the latter
    btchs <- Exp.map[smplsMtch, myBatch2]
    tmpDat2Imp[w, currSamples] <- ComBat(tmpDat1Imp[w,],
                                         btchs,
                                         mdlMtr,
                                         par.prior = TRUE)
  }
  
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
                     myBatch2,
                     myBatch2,
                     Exp.map,
                     intRoot = "")
PCAlyLst[[myBatch2]] <- tmp$PlotLy
scoresLst[[myBatch2]] <- tmp$Scores
PCsLst[[myBatch2]] <- tmp$PCs
#
appNm <- paste0("Batch corr.: original -> ", myBatch3)
msg <- "Keep results from ComBat batch correction? (untick to cancel correction)"
corrTst <- t.test(unlist(tmpDat2Imp[wAG1, currSamples]),
                  unlist(tmpDat1Imp[wAG1, currSamples]),
                  paired = TRUE)
if ("Decision" %in% (normSequence[[nrmStp]])) {
  KeepComBatRes <- normSequence[[nrmStp]]$Decision
}
if (!validLogicPar("KeepComBatRes")) {
  KeepComBatRes <- corrTst$p.value < 0.01 # Default: keep results only if data is very significantly different
}
l <- min(c(length(PCsLst$original$sdev), length(PCsLst[[myBatch2]]$sdev)))
l2 <- min(c(5L, l))
PCs <- data.frame("Component" = paste0("PC", as.character(1L:l2)),
                  "Before (%)" = round(100*(PCsLst[["original"]]$sdev[1L:l2])^2L / sum(PCsLst[["original"]]$sdev^2L), 0L),
                  "After (%)" = round(100*(PCsLst[[myBatch2]]$sdev[1L:l2])^2L / sum(PCsLst[[myBatch2]]$sdev^2L), 0L),
                  check.names = FALSE)
if (l2 < l) {
  PCs <- rbind(PCs, data.frame(Component = "...",
                               "Before (%)" = "...",
                               "After (%)" = "...",
                               check.names = FALSE))
}
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
                         class = "compact",
                         selection = "none",
                         editable = FALSE,
                         rownames = FALSE,
                         options = list(
                           dom = 't',
                           paging = FALSE,
                           ordering = FALSE
                         ),
      #                    callback = JS("table.rows().every(function(i, tab, row) {
      #   var $this = $(this.node());
      #   $this.attr('id', this.data()[0]);
      #   $this.addClass('shiny-input-container');
      # });
      # Shiny.unbindAll(table.table().node());
      # Shiny.bindAll(table.table().node());")
      )
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
runKount <- 0L
while ((!runKount)||(!exists("IHAVERUN"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  shinyCleanup()
  runKount <- runKount+1L
}
msg <- paste0(" -> correction of ", myBatch3, "-associated batch effect ", c("rejec", "accep")[KeepComBatRes+1], "ted.\n")
if (KeepComBatRes) {
  txt2 <- paste0("corrected against the ", myBatch3, "-associated batch effect")
}
normSequence[[nrmStp]]$Decision <- Outcome <- KeepComBatRes
cat(msg)
