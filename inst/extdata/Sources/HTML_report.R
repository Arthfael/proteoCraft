# Create html report to:
# - organize the most important tables and plots into a coherent html file
# - capture user comments on each section

library(shiny)
library(shinyjs)
library(bslib)
library(htmltools)
library(DT)
library(plotly)
library(jsonlite)

htmlRprtFl <- paste0(wd, "/Report_", dtstNm, ".html")
if (scrptType == "withReps") {
  smplGrps <- setNames(cleanNms(VPAL$values), VPAL$values)
  if (!"clean_Group_name" %in% colnames(Exp.map)) {
    Exp.map$clean_Group_name <- smplGrps[Exp.map[[VPAL$column]]]
  }
}

# Reload processed data from the report
if (!exists("xlDat")) { loadFun(paste0(wd, "/Tables/xlDat.RDS")) }
peptidoTst <- "All peptidoforms" %in% names(xlDat)

# Reload materials and methods
if (!exists("matmethTxt")) {
  matmethTxt <- list("Samples preparation" = MatMetCalls$Texts$WetLab,
                     "LC-MS/MS analysis" = MatMetCalls$Texts$LCMS,
                     "Data analysis" = MatMetCalls$Texts$DatAnalysis)
}
matmethSections <- names(matmethTxt)

# Reload plots data
tstRat <- (scrptType == "noReps") && MakeRatios && exists("ratioPlots_fl") && file.exists(ratioPlots_fl)
if (tstRat) {
  loadFun(ratioPlots_fl)
  tstRat <- length(ratioPlots) > 0L
}
heatMaps_ON <- file.exists(heatMaps_fl)
if (heatMaps_ON) {
  loadFun(heatMaps_fl)
  heatMaps_ON <- exists("plotLeatMaps") && length(plotLeatMaps)
}
dimRed_fl <- paste0(wd, "/Dimensionality red. plots/DimRedPlots.RDS")
PCA_ON <- file.exists(dimRed_fl)
if (PCA_ON) {
  loadFun(dimRed_fl)
  PCA_ON <- exists("dimRedPlotLy") && ("PCA" %in% names(dimRedPlotLy$PG))
}
Venn_ON <- exists("Venn_fl") && file.exists(Venn_fl)
if (Venn_ON) {
  loadFun(Venn_fl)
  Venn_ON <- exists("plotly_Venn") && ("Global, LFQ" %in% names(plotly_Venn))
}
strtColWdth <- 12L/max(c(1L, Venn_ON + PCA_ON))
tstCov <- exists("covPlots_fl") && file.exists(covPlots_fl)
if (tstCov) {
  loadFun(covPlots_fl)
  tstCov <- length(covPlots) > 0L
}

# Fix to plotly autoscaling + remove some Modebar tools (redundant: they should already be gone, but in case we reload old data) + fix warnings
global_autorange <- "function(el, x) {
  var gd = el;
  function globalRange(axisPrefix) {
    var axes = Object.keys(gd._fullLayout).filter(function(k) {
      return k.match(new RegExp('^' + axisPrefix + 'axis[0-9]*$'));
    });
    if (axes.length <= 1)
      return;
    var minVal = Infinity;
    var maxVal = -Infinity;
    axes.forEach(function(name) {
      var axis = gd._fullLayout[name];
      if (axis && axis.range) {
        minVal = Math.min(minVal, axis.range[0], axis.range[1]);
        maxVal = Math.max(maxVal, axis.range[0], axis.range[1]);
      }
    });
    if (isFinite(minVal) && isFinite(maxVal)) {
      axes.forEach(function(name) {
        Plotly.relayout(gd, name + '.range', [minVal, maxVal]);
      });
    }
  }
  globalRange('x');
  globalRange('y');
}"
if (heatMaps_ON) {
  for (x in names(plotLeatMaps)) { #x <- names(plotLeatMaps)[1L]
    for (y in names(plotLeatMaps[[x]])) { #y <- names(plotLeatMaps[[x]])[1L]
      p <- plotLeatMaps[[x]][[y]]$Plot
      p$x$layout$xaxis$autorange <- TRUE
      p$x$layout$yaxis$autorange <- TRUE
      p <- htmlwidgets::onRender(p, global_autorange)
      plotLeatMaps[[x]][[y]]$Plot <- plotly::config(p,
                                                    modeBarButtonsToRemove = c("select2d", "lasso2d"))
    }
  }
}
loadFun(paste0(wd, "/Ranked abundance/quantPlots.RDS"))
for (x in names(ggQuantLy)) { #x <- names(ggQuantLy)[1L]
  for (y in names(ggQuantLy[[x]])) { #y <- names(ggQuantLy[[x]])[1L]
    p <- ggQuantLy[[x]][[y]]$plotly
    p$plotly$x$layout$xaxis$autorange <- TRUE
    p$plotly$x$layout$yaxis$autorange <- TRUE
    p <- htmlwidgets::onRender(p, global_autorange)
    ggQuantLy[[x]][[y]]$plotly <- plotly::config(p,
                                                 modeBarButtonsToRemove = c("select2d", "lasso2d"))
  }
}
if (PCA_ON) {
  for (x in names(dimRedPlotLy)) { #x <- names(dimRedPlotLy)[1L]
    for (y in names(dimRedPlotLy[[x]])) { #x <- names(dimRedPlotLy[[x]])[1L]
      p <- dimRedPlotLy$PG[[x]][[y]]
      p$x$layout$xaxis$autorange <- TRUE
      p$x$layout$yaxis$autorange <- TRUE
      p <- htmlwidgets::onRender(p, global_autorange)
      dimRedPlotLy[[x]][[y]] <- plotly::config(p,
                                               modeBarButtonsToRemove = c("select2d", "lasso2d"))
    }
  }
}
if (tstCov) {
  for (x in names(covPlots)) { #x <- names(covPlots)[1L]
    for (y in names(covPlots[[x]])) { #y <- names(covPlots[[x]])[1L]
      for (z in names(covPlots[[x]][[y]])) { #z <- names(covPlots[[x]][[y]])[1L]
        p <- covPlots[[x]][[y]][[z]]
        p$x$layout$xaxis$autorange <- TRUE
        p$x$layout$yaxis$autorange <- TRUE
        p <- htmlwidgets::onRender(p, global_autorange)
        covPlots[[x]][[y]][[z]] <- plotly::config(p,
                                                  modeBarButtonsToRemove = c("select2d", "lasso2d"))
      }
    }
  }
}
if (tstRat) {
  for (x in names(ratioPlots)) { #x <- names(ratioPlots)[1L]
    p <- ratioPlots[[x]]
    p$x$layout$xaxis$autorange <- TRUE
    p$x$layout$yaxis$autorange <- TRUE
    p <- htmlwidgets::onRender(p, global_autorange)
    ratioPlots[[x]] <- plotly::config(p,
                                      modeBarButtonsToRemove = c("select2d", "lasso2d"))
  }
}
loadFun(qcBckUpFl)
for (x in names(QC_plotLys)) { #x <- names(QC_plotLys)[1L]
  p <- QC_plotLys[[x]]
  p$x$layout$xaxis$autorange <- TRUE
  p$x$layout$yaxis$autorange <- TRUE
  p <- htmlwidgets::onRender(p, global_autorange)
  QC_plotLys[[x]] <- plotly::config(p,
                                    modeBarButtonsToRemove = c("select2d", "lasso2d"))
}
if (Venn_ON) {
  for (x in names(plotly_Venn)) { #x <- names(plotly_Venn)[1L]
    p <- plotly_Venn[[x]]
    p$x$layout$xaxis$autorange <- TRUE
    p$x$layout$yaxis$autorange <- TRUE
    plotly_Venn[[x]] <- plotly::config(p,
                                       modeBarButtonsToRemove = c("select2d", "lasso2d"))
  }
}
loadFun(GO_plot_ly_fl)
for (tt in names(GO_plot_ly$PG)) {
  for (x in names(GO_plot_ly$PG[[tt]])) { #x <- names(GO_plot_ly$PG[[tt]])[1L]
    for (nm in c("Bar", "Bubble")) {
      p <- GO_plot_ly$PG[[tt]][[x]][[nm]]
      p$x$layout$xaxis$autorange <- TRUE
      p$x$layout$yaxis$autorange <- TRUE
      GO_plot_ly$PG[[tt]][[x]][[nm]] <- plotly::config(p,
                                                                modeBarButtonsToRemove = c("select2d", "lasso2d"))
    }
  }
}
if (!is.null(GO_plot_ly$Prot$SAINTexpress)) {
  for (x in names(GO_plot_ly$Prot$SAINTexpress)) { #x <- names(GO_plot_ly$Prot$SAINTexpress)[1L]
    for (nm in c("Bar", "Bubble")) {
      p <- GO_plot_ly$Prot$SAINTexpress[[x]][[nm]]
      p$x$layout$xaxis$autorange <- TRUE
      p$x$layout$yaxis$autorange <- TRUE
      GO_plot_ly$Prot$SAINTexpress[[x]][[nm]] <- plotly::config(p,
                                                                modeBarButtonsToRemove = c("select2d", "lasso2d"))
    }
  }
}

#
#plotHght <- "400px"
plotHght <- paste0(round(screenRes$height*0.75), "px")
nmsHtMp <- names(plotLeatMaps$Global)
if ((scrptType == "noReps") && (length(Exp) == 2L)) {
  nmsHtMp <- setdiff(nmsHtMp, "Z-scored")
}
nmsHtMp <- intersect(union("None", nmsHtMp), nmsHtMp)
plotHtMpHght <- paste0(round(min(c(400L, vapply(nmsHtMp, \(nm) { plotLeatMaps$Global[[nm]]$Plot$sizingPolicy$defaultHeight }, 1)))), "px")
plotPCAHght <- "700px"

# UI functions
tbl_css <- tags$style(HTML("table.dataTable td {
  white-space: normal !important;
  vertical-align: top !important;
}
.cell-wrap {
  max-height: 80px !important;
  overflow: hidden !important;
  white-space: normal;
  vertical-align: top !important;
  overflow-wrap: break-word !important;
  word-break: break-word !important;
}
table.dataTable thead th {
    background-color: #4472c4 !important;
    color: white !important;
}"))
logoFl <- list.files(homePath,
                     "^logo\\.((gif)|(tiff?)|(jpe?g)|(png))$", full.names = TRUE)[1L]
report_header <- tags$header(
  fluidRow(if (length(logoFl)) {
    column(2L,
           tags$img(src = knitr::image_uri(logoFl),
                    style = "max-width: 100%; height: auto;"))
  },
  column(10L,
         br(),
         paste0(dtstNm, " - report"),
         br(),
         "Analysis run by: ", em(WhoAmI),
         br(),
         "Date: ", em(Sys.Date()),
         br(),
         "Package: ", em(paste0("proteoCraft v", package.version("proteoCraft"))),
         br(),
         br())),
  style = "background: linear-gradient(to right, #e8f1ff, #ffffff); padding: 1rem; margin-bottom: 1rem;")


# Functions
make_comment_ui <- \(id,
                     shiny = TRUE,
                     values = allComments,
                     ON = TRUE,
                     root = "comment_") {
  if (shiny) {
    textAreaInput(inputId = paste0(root, id),
                  label = NULL,
                  value = values[id],
                  width = "100%",
                  height = "150px")
  } else {
    style <- if (ON) { "display: block; white-space: pre-wrap; padding: 10px;" } else { "display: none;" }
    tags$div(class = "comment-box",
             id = paste0(root, match(id, names(values))),
             style = style,
             values[id])
  }
}
make_bar <- \(x) {
  sprintf("<div style=\"position: relative; width: 100%%; background: #eee; height: 16px; border-radius: 4px;\">
  <div style=\"width: %s%%; background: #4CAF50; height: 100%%; border-radius: 4px;\"></div>
  <div style=\"position: absolute; top: 0; left: 50%%; transform: translateX(-50%%); font-size: 11px; line-height: 16px; color: black;\">
    %.1f%%
  </div>
</div>", x, x)
}
make_ctrst_tbl_ui <- \(contr, #contr <- myContrasts$Contrast[1L] #contr <- myContrasts$Contrast[2L] #contr <- myContrasts$Contrast[3L]
                       tab = "Protein groups", # can also be "All peptidoforms"; we will eventually add "`PTM`-modified", where `PTM` can be any PTM of interest
                       filt = NULL, #filt = allProt[1L] # Filter by "Common Name"
                       dat = xlDat,
                       minN = 1L) {
  m <- match(contr, myContrasts$Contrast)
  exp <- setNames(myContrasts[m, c("A_samples", "B_samples")],
                  myContrasts[m, c("A", "B")])
  exp <- lapply(exp, \(x) { Exp.map$Clean_name[match(unlist(x), Exp.map[[RSA$column]])] })
  grps <- names(exp)
  pgTest <- (tab == "Protein groups")
  df <- dat[[tab]]
  tmp <- grep("\n", colnames(df), value = TRUE)
  tmp2 <- gsub(" /\n.*", "", tmp)
  tmp2 <- gsub(".*\n", "", tmp2)
  smplCols_lst <- setNames(lapply(grps, \(xp) {
    tmp[tmp2 %in% c(xp, exp[[xp]])]
  }), grps)
  smplCols <- setNames(unlist(smplCols_lst), NULL)
  coreCols <- "PEP"
  if (tab %in% c("Protein groups", "All peptidoforms")) {
    if (pgTest) {
      filtCol <- "Protein IDs"
      coreCols <- union(c("Leading protein IDs", filtCol, "Common Names", "Genes", "Mol. weight [kDa]", "Potential contaminant"), coreCols)
      intRoot <- "expr"
    }
    if (tab == "All peptidoforms") {
      filtCol <- "Proteins"
      coreCols <- union(c("Modified sequence_verbose", #"Sequence",
                          filtCol), coreCols)
      intRoot <- "int"
    }
  } else {
    stop("TO DO!")
    filtCol <- "Proteins"
    coreCols <- union(c("Modified sequence_verbose", #"Sequence",
                        filtCol), coreCols) # Check before use...
    intRoot <- "int" # Presumably...
  }
  #
  xprCols <- grep(paste0("log10\\(([^\\)]+ )?", intRoot, "\\.\\) "), smplCols, value = TRUE)
  fullIntRoot <- rev(paste0(vapply(strsplit(xprCols, "\n"), `[[`, "", 1L), "\n"))[1L]
  smpls_and_grps <- unlist(lapply(grps, \(grp) { c(grp, exp[[grp]]) }))
  xprCols <- paste0(fullIntRoot, smpls_and_grps)
  xprCols <- intersect(xprCols, colnames(df))
  repXprCols <- if (length(exp) == 1L) { sub(" *\n$", "", fullIntRoot) } else { xprCols }
  #
  ratCols <- grep("log2\\(.*rat\\.\\) \n", smplCols, value = TRUE)
  stopifnot(length(ratCols) > 0L) # For contrasts we always MUST have a logFC column!
  fullRatRoot <- rev(paste0(vapply(strsplit(ratCols, "\n"), `[[`, "", 1L), "\n"))[1L]
  ratCol <- paste0(fullRatRoot, sub(" - ", " /\n", contr))
  ratCol <- repRatCol <- intersect(ratCol, colnames(df)) 
  stopifnot(length(ratCol) == 1L) # Again
  #
  colNms <- c(coreCols, xprCols, ratCol)
  repColNms <- c(sub("_verbose$", "",
                     sub("^Potential contaminant$", "Cont.",
                         sub("^Mol\\. weight \\[kDa\\]$", "MW (kDa)",
                             sub("^Common Names$", "Common names", coreCols)))),
                 repXprCols, repRatCol)
  if (pgTest) {
    pepCountCols <- intersect(paste0("Pep. count \n", grps), colnames(df))
    psmCountCols <- intersect(paste0("PSMs count \n", grps), colnames(df))
    k <- union(pepCountCols, psmCountCols)
    if (length(k)) {
      colNms <- union(colNms, k)
      repColNms <- if (length(exp) == 1L) { union(repColNms, c("Pep. count", "PSMs count")) } else { union(repColNms, k) }
    }
  }
  df <- df[, colNms]
  flt <- if (is.null(filt)) {
    1L:nrow(df)
  } else {
    grsep(db$`Protein ID`[match(filt, db$`Common Name`)], x = df[[filtCol]])
  }
  if (pgTest && is.integer(minN) && (minN > 0L) && length(pepCountCols)) {
    flt <- flt[apply(df[flt, pepCountCols, drop = FALSE], 1L, max, na.rm = TRUE) >= minN]
  }
  if (!length(flt)) { return() }
  df <- df[flt,]
  colnames(df) <- colNms <- repColNms
  if ("Modified sequence" %in% colNms) {
    df$"Modified sequence" <- gsub("^_|_$", "",  df$"Modified sequence")
  }
  xprCols <- repXprCols
  ratCol <- repRatCol
  covCols <- NULL
  if (pgTest) {
    covCols <- paste0("Cov. \n", smpls_and_grps)
    repCovCols <- if (length(exp) == 1L) { "Coverage" } else { covCols }
    if (sum(covCols %in% colnames(dat[[tab]]))) {
      df[, repCovCols] <- dat[[tab]][flt, covCols]
    } else {
      if (("Coverage" %in% names(dat)) && (sum(covCols %in% colnames(dat$Coverage)))) {
        m <- match(df$"Protein IDs", dat$Coverage$`Protein IDs`)
        df[, repCovCols] <- dat$Coverage[m, covCols]
      }
    }
    covCols <- repCovCols
  }
  xprRng <- range(df[, xprCols], na.rm = TRUE)
  #
  # Make sure this re-ordering is done after any other data is added from dat to df!
  orderVect <- df[, xprCols]
  if (length(exp) > 1L) { orderVect <- apply(orderVect, 1L, \(x) { mean(x[is.finite(x)]) }) }
  df <- df[order(orderVect, decreasing = TRUE),]
  #
  quantCols <- xprCols
  ratRng <- range(df[, ratCol], na.rm = TRUE)
  quantCols <- c(quantCols, ratCol)
  covSortCols <- character(0L)
  covSortVals <- NULL
  if (pgTest && (!is.null(covCols)) && length(covCols)) {
    covSortCols <- paste0(covCols, "__sort")
    covSortVals <- setNames(lapply(covCols, \(k) {
      suppressWarnings(as.numeric(df[[k]]))
    }), covSortCols)
  }
  col2 <- setdiff(colnames(df), c("PEP", quantCols))
  df[, col2] <- sapply(col2, \(k) {
    if (pgTest && (k %in% covCols)) {
      make_bar(df[[k]])
    } else {
      sprintf("<div class=\"cell-wrap\">%s</div>", df[[k]])
    }
  })
  if (length(covSortCols)) {
    df[covSortCols] <- covSortVals
  }
  wTest1 <- setNames(vapply(colnames(df), \(k) { #k <- colnames(df)[1L]
    # tmp <- as.character(df[[k]])
    # x <- min(c(250L, max(nchar(c(k, tmp)) + 3L, na.rm = TRUE)*8L))
    # if (is.na(x)) { x <- 50L }
    # return(as.integer(x))
    max(c(min(c(nchar(k)*8L + 24L,
                250L)),
          50L))
  }, 1L), colnames(df))
  wTest2 <- sum(wTest1) + 15L + ncol(df)*5L
  wTest1 <- paste0(as.character(wTest1), "px")
  wTest1 <- aggregate((1L:length(wTest1)) - 1L, list(wTest1), c)
  wTest1 <- apply(wTest1, 1L, \(x) {
    x2 <- as.integer(x[[2L]])
    list(width = x[[1L]],
         targets = x2,
         names = colnames(df)[x2 + 1L])
  })
  covOrderDefs <- list()
  if (length(covSortCols)) {
    covVisibleTargets <- match(covCols, colnames(df)) - 1L
    covSortTargets <- match(covSortCols, colnames(df)) - 1L
    covOrderDefs <- c(Map(\(visible_col, sort_col) {
      list(targets = visible_col,
           orderData = sort_col)
    },
    covVisibleTargets,
    covSortTargets),
    list(list(targets = covSortTargets,
              visible = FALSE,
              searchable = FALSE)))
  }
  columnDefs_all <- c(unname(wTest1), covOrderDefs)
  header_help <- c("PEP" = c("Posterior Error Probability = estimate of the probability that the peptide is a false discovery (local FDR) based on the local density of decoy hits",
                             "Posterior Error Probability = estimate of the probability that the protein group is a false discovery (i.e. all its assigned peptides are false discoveries), calculated as the product of the PEPs of individual peptides)")[pgTest + 1L],
                   "Leading protein IDs" = "Accessions of the minimum number of proteins required to explain the observed peptides assigned to the protein group",
                   "Protein IDs" = "Accessions of all proteins the peptides in the group could originate from",
                   "Proteins" = "Accessions of all proteins the peptide could originate from",
                   "Modified sequence" = "Peptide sequence with any post-translational modification(s) detected",
                   "Cont." = paste0(c("Peptides matching", "Proteins from")[pgTest + 1L],
                                    " a list of common environmental and laboratory contaminants, inc. e.g. Trypsin, Keratins, BSA,... are marked with a \"+\""),
                   "Pep. count" = "Number of peptidoforms",
                   "PSMs count" = "Number of individual identifications (Peptide-to-Spectrum matches)",
                   "MW (kDa)" = "Molecular weight of the first protein in column \"Leading protein IDs\"",
                   "Coverage" = "Sequence coverage of the first protein in column \"Leading protein IDs\"",
                   "log10(" = c("log10 intensity-based estimated abundance",
                                "log10 intensity-based estimated abundance")[pgTest + 1L], # For now the same, but should become more taylor-made based on which specific intensity/expression column we choose to show (e.g. normalised, imputed, corrected etc...)
                   "log2(" = "Estimated log2 Fold Change",
                   "pval." = "P-value (0 is good, 1 is bad) of the statistical test for this contrast",
                   "-log10 pval." = "-log10-transformed P-value (infinite is good, 0 is bad) of the statistical test for this contrast",
                   "reg." = "Final decision on differential expression")
  df <- DT::datatable(df,
                      rownames = FALSE,
                      class = "compact",
                      escape = FALSE,
                      options = list(scrollX = TRUE,
                                     scrollY = "500px",
                                     pageLength = 100L,
                                     lengthMenu = list(c(10L, 25L, 50L, 100L, -1L),
                                                       c("10", "25", "50", "100", "All")),
                                     autoWidth = FALSE,
                                     columnDefs = columnDefs_all,
                                     initComplete = JS(sprintf("function(settings, json) {
  const tips = %s;
  const api = this.api();
  // Sort prefixes longest-first, so more specific matches win
  const entries = Object.entries(tips).sort(([a], [b]) => b.length - a.length);
  api.columns().every(function(i) {
    const header = api.column(i).header();
    const colName = header.textContent.trim();
    for (const [prefix, helpText] of entries) {
      if (colName.startsWith(prefix)) {
        header.setAttribute('title', helpText);
        break;
      }
    }
  });
}",
                                                               jsonlite::toJSON(as.list(header_help), auto_unbox = TRUE)))))
  df <- DT::formatRound(df, c("PEP", quantCols), digits = 5L)
  df <- DT::formatStyle(df, "PEP",
                        backgroundColor = DT::styleInterval(10L^-seq(10, 0, length.out = 99L),
                                                            colorRampPalette(rev(ColScaleList$PEP))(100L)))
  df <- DT::formatStyle(df, xprCols,
                        backgroundColor = DT::styleInterval(seq(xprRng[1L], xprRng[2L], length.out = 99L),
                                                            colorRampPalette(ColScaleList$`Individual Expr`)(100L)))
  df <- DT::formatStyle(df, ratCol,
                        backgroundColor = DT::styleInterval(seq(ratRng[1L], ratRng[2L], length.out = 99L),
                                                            colorRampPalette(ColScaleList$`Summary Ratios`)(100L)))
  return(tags$div(df,
                  style = paste0("background: #ffffff;")))
}
make_smpl_tbl_ui <- \(exp, #exp <- Exp[1L] #exp <- Exp[2L] #exp <- smplGrps[1L]
                      tab = "Protein groups", # can also be "All peptidoforms"; we will eventually add "`PTM`-modified", where `PTM` can be any PTM of interest
                      filt = NULL, #filt = allProt[1L] # Filter by "Common Name"
                      dat = xlDat,
                      minN = 1L) {
  pgTest <- (tab == "Protein groups")
  df <- dat[[tab]]
  if (missing(exp)) {
    exp <- if (scrptType == "noReps") { Exp } else { smplGrps }
  }
  if ((scrptType == "withReps") && (tab == "All peptidoforms")) {
    exp <- unique(unlist(lapply(exp, \(x) {
      Exp.map$Clean_name[Exp.map$clean_Group_name == x]
    })))
  }
  smplCols_lst <- setNames(lapply(exp, \(xp) {
    grep(topattern(paste0("\n", xp), FALSE, TRUE), colnames(df), value = TRUE)
  }), exp)
  smplCols <- setNames(unlist(smplCols_lst), NULL)
  coreCols <- "PEP"
  if (tab %in% c("Protein groups", "All peptidoforms")) {
    if (pgTest) {
      filtCol <- "Protein IDs"
      coreCols <- union(c("Leading protein IDs", filtCol, "Common Names", "Genes", "Mol. weight [kDa]", "Potential contaminant"), coreCols)
      intRoot <- "expr"
    }
    if (tab == "All peptidoforms") {
      filtCol <- "Proteins"
      coreCols <- union(c("Modified sequence_verbose", #"Sequence",
                          filtCol), coreCols)
      intRoot <- "int"
    }
  } else {
    stop("TO DO!")
    filtCol <- "Proteins"
    coreCols <- union(c("Modified sequence_verbose", #"Sequence",
                        filtCol), coreCols) # Check before use...
    intRoot <- "int" # Presumably...
  }
  #
  xprCols <- grep(paste0("log10\\(([^\\)]+ )?", intRoot, "\\.\\) "), smplCols, value = TRUE)
  fullIntRoot <- rev(paste0(vapply(strsplit(xprCols, "\n"), `[[`, "", 1L), "\n"))[1L]
  xprCols <- paste0(fullIntRoot, exp)
  repXprCols <- if (length(exp) == 1L) { sub(" *\n$", "", fullIntRoot) } else { xprCols }
  #
  ratCols <- grep(paste0("log2\\(.*rat\\.\\) \n"), smplCols, value = TRUE)
  useRat <- length(ratCols) > 0L
  if (useRat) {
    fullRatRoot <- rev(paste0(vapply(strsplit(ratCols, "\n"), `[[`, "", 1L), "\n"))[1L]
    if (scrptType == "noReps") {
      ratCols <- paste0(fullRatRoot, exp)
    } else {
      ratCols <- unlist(lapply(exp, \(xp) { grep(topattern(paste0(fullRatRoot, xp, " /\n")), smplCols, value = TRUE) }))
    }
    repRatCols <- ratCols <- intersect(ratCols, colnames(df))
    useRat <- length(ratCols) > 0L
  }
  #
  colNms <- c(coreCols, xprCols)
  repColNms <- c(sub("_verbose$", "",
                     sub("^Potential contaminant$", "Cont.",
                         sub("^Mol\\. weight \\[kDa\\]$", "MW (kDa)",
                             sub("^Common Names$", "Common names", coreCols)))),
                 repXprCols)
  if (useRat) {
    colNms <- c(colNms, ratCols)
    repColNms <- c(repColNms, repRatCols)
  }
  if (pgTest) {
    pepCountCols <- intersect(paste0("Pep. count \n", exp), colnames(df))
    psmCountCols <- intersect(paste0("PSMs count \n", exp), colnames(df))
    k <- union(pepCountCols, psmCountCols)
    if (length(k)) {
      colNms <- union(colNms, k)
      repColNms <- if (length(exp) == 1L) { union(repColNms, c("Pep. count", "PSMs count")) } else { union(repColNms, k) }
    }
  }
  df <- df[, colNms]
  flt <- if (is.null(filt)) {
    1L:nrow(df)
  } else {
    grsep(db$`Protein ID`[match(filt, db$`Common Name`)], x = df[[filtCol]])
  }
  if (pgTest && is.integer(minN) && (minN > 0L) && length(pepCountCols)) {
    flt <- flt[apply(df[flt, pepCountCols, drop = FALSE], 1L, max, na.rm = TRUE) >= minN]
  }
  if (!length(flt)) { return() }
  df <- df[flt,]
  colnames(df) <- colNms <- repColNms
  if ("Modified sequence" %in% colNms) {
    df$"Modified sequence" <- gsub("^_|_$", "",  df$"Modified sequence")
  }
  xprCols <- repXprCols
  if (useRat) { ratCols <- repRatCols }
  covCols <- NULL
  if (pgTest) {
    covCols <- paste0("Cov. \n", exp)
    repCovCols <- if (length(exp) == 1L) { "Coverage" } else { covCols }
    if (sum(covCols %in% colnames(dat[[tab]]))) {
      df[, repCovCols] <- dat[[tab]][flt, covCols]
    } else {
      if (("Coverage" %in% names(dat)) && (sum(covCols %in% colnames(dat$Coverage)))) {
        m <- match(df$"Protein IDs", dat$Coverage$`Protein IDs`)
        df[, repCovCols] <- dat$Coverage[m, covCols]
      }
    }
    covCols <- repCovCols
  }
  xprRng <- range(df[, xprCols], na.rm = TRUE)
  #
  # Make sure this re-ordering is done after any other data is added from dat to df!
  orderVect <- df[, xprCols]
  if (length(exp) > 1L) { orderVect <- apply(orderVect, 1L, \(x) { mean(x[is.finite(x)]) }) }
  df <- df[order(orderVect, decreasing = TRUE),]
  #
  quantCols <- xprCols
  if (useRat) {
    ratRng <- range(df[, ratCols], na.rm = TRUE)
    quantCols <- c(quantCols, ratCols)
  }
  covSortCols <- character(0L)
  covSortVals <- NULL
  if (pgTest && (!is.null(covCols)) && length(covCols)) {
    covSortCols <- paste0(covCols, "__sort")
    covSortVals <- setNames(lapply(covCols, \(k) {
      suppressWarnings(as.numeric(df[[k]]))
    }), covSortCols)
  }
  col2 <- setdiff(colnames(df), c("PEP", quantCols))
  df[, col2] <- sapply(col2, \(k) {
    if (pgTest && (k %in% covCols)) {
      make_bar(df[[k]])
    } else {
      sprintf("<div class=\"cell-wrap\">%s</div>", df[[k]])
    }
  })
  if (length(covSortCols)) {
    df[covSortCols] <- covSortVals
  }
  wTest1 <- setNames(vapply(colnames(df), \(k) { #k <- colnames(df)[1L]
    # tmp <- as.character(df[[k]])
    # x <- min(c(250L, max(nchar(c(k, tmp)) + 3L, na.rm = TRUE)*8L))
    # if (is.na(x)) { x <- 50L }
    # return(as.integer(x))
    max(c(min(c(nchar(k)*8L + 24L,
                250L)),
          50L))
  }, 1L), colnames(df))
  wTest2 <- sum(wTest1) + 15L + ncol(df)*5L
  wTest1 <- paste0(as.character(wTest1), "px")
  wTest1 <- aggregate((1L:length(wTest1)) - 1L, list(wTest1), c)
  wTest1 <- apply(wTest1, 1L, \(x) {
    x2 <- as.integer(x[[2L]])
    list(width = x[[1L]],
         targets = x2,
         names = colnames(df)[x2 + 1L])
  })
  covOrderDefs <- list()
  if (length(covSortCols)) {
    covVisibleTargets <- match(covCols, colnames(df)) - 1L
    covSortTargets <- match(covSortCols, colnames(df)) - 1L
    covOrderDefs <- c(Map(\(visible_col, sort_col) {
      list(targets = visible_col,
           orderData = sort_col)
    },
    covVisibleTargets,
    covSortTargets),
    list(list(targets = covSortTargets,
              visible = FALSE,
              searchable = FALSE)))
  }
  columnDefs_all <- c(unname(wTest1), covOrderDefs)
  header_help <- c("PEP" = c("Posterior Error Probability = estimate of the probability that the peptide is a false discovery (local FDR) based on the local density of decoy hits",
                             "Posterior Error Probability = estimate of the probability that the protein group is a false discovery (i.e. all its assigned peptides are false discoveries), calculated as the product of the PEPs of individual peptides)")[pgTest + 1L],
                   "Leading protein IDs" = "Accessions of the minimum number of proteins required to explain the observed peptides assigned to the protein group",
                   "Protein IDs" = "Accessions of all proteins the peptides in the group could originate from",
                   "Proteins" = "Accessions of all proteins the peptide could originate from",
                   "Modified sequence" = "Peptide sequence with any post-translational modification(s) detected",
                   "Cont." = paste0(c("Peptides matching", "Proteins from")[pgTest + 1L],
                                    " a list of common environmental and laboratory contaminants, inc. e.g. Trypsin, Keratins, BSA,... are marked with a \"+\""),
                   "Pep. count" = "Number of peptidoforms",
                   "PSMs count" = "Number of individual identifications (Peptide-to-Spectrum matches)",
                   "MW (kDa)" = "Molecular weight of the first protein in column \"Leading protein IDs\"",
                   "Coverage" = "Sequence coverage of the first protein in column \"Leading protein IDs\"",
                   "log10(" = c("log10 intensity-based estimated abundance",
                                "log10 intensity-based estimated abundance")[pgTest + 1L], # For now the same, but should become more taylor-made based on which specific intensity/expression column we choose to show (e.g. normalised, imputed, corrected etc...)
                   "log2(" = "estimated log2 Fold Change")
  df <- DT::datatable(df,
                      rownames = FALSE,
                      class = "compact",
                      escape = FALSE,
                      options = list(scrollX = TRUE,
                                     scrollY = "500px",
                                     pageLength = 100L,
                                     lengthMenu = list(c(10L, 25L, 50L, 100L, -1L),
                                                       c("10", "25", "50", "100", "All")),
                                     autoWidth = FALSE,
                                     columnDefs = columnDefs_all,
                                     initComplete = JS(sprintf("function(settings, json) {
  const tips = %s;
  const api = this.api();
  // Sort prefixes longest-first, so more specific matches win
  const entries = Object.entries(tips).sort(([a], [b]) => b.length - a.length);
  api.columns().every(function(i) {
    const header = api.column(i).header();
    const colName = header.textContent.trim();
    for (const [prefix, helpText] of entries) {
      if (colName.startsWith(prefix)) {
        header.setAttribute('title', helpText);
        break;
      }
    }
  });
}",
                                                               jsonlite::toJSON(as.list(header_help), auto_unbox = TRUE)))))
  df <- DT::formatRound(df, c("PEP", quantCols), digits = 5L)
  df <- DT::formatStyle(df, "PEP",
                        backgroundColor = DT::styleInterval(10L^-seq(10, 0, length.out = 99L),
                                                            colorRampPalette(rev(ColScaleList$PEP))(100L)))
  df <- DT::formatStyle(df, xprCols,
                        backgroundColor = DT::styleInterval(seq(xprRng[1L], xprRng[2L], length.out = 99L),
                                                            colorRampPalette(ColScaleList$`Individual Expr`)(100L)))
  if (useRat) {
    df <- DT::formatStyle(df, ratCols,
                          backgroundColor = DT::styleInterval(seq(ratRng[1L], ratRng[2L], length.out = 99L),
                                                              colorRampPalette(ColScaleList$`Individual Ratios`)(100L)))
  }
  return(tags$div(df,
                  style = paste0("background: #ffffff;")))
}
make_prot_tab <- \(dflt = dfltProt,
                   prots = allProt,
                   shiny = TRUE) {
  myCol <- tolower(viridis::viridis(6L, alpha = 0.2)[4L])
  myExp <- if (scrptType == "noReps") { Exp } else { setNames(smplGrps, NULL) }
  # - show:
  # Proteins tab
  #################################
  #     dropdown for protein      #
  #################################
  # ->
  #################################
  #      comment for protein      #
  #################################
  ################## ##############
  #samples dropdown# #            #
  ################## #            #
  ################## #   Ratios   #
  #                # #    plot    #
  #    Coverage    # #            #
  #                # #            #
  ################## ##############
  # Peptides table
  if (shiny) {
    tagList(tags$div(
      if (prot.list.Cond && (length(prots) > 1L)) {
        tags$div(selectInput("myProtein", "Select protein", prots, dflt),
                 br())
      },
      uiOutput("protComment"),
      br(),
      fluidRow(column(6L,
                      if (length(myExp) > 1L) {
                        selectInput("mySample", "", myExp, myExp[1L]) 
                      },
                      plotlyOutput("coverPlot", height = plotHght)),
               if (scrptType == "noReps") {
                 column(6L,
                        plotlyOutput("ratioPlot", height = plotHght))
               },
      ),
      br(),
      br(),
      tags$hr(style = "border-color: black;"),
      if (peptidoTst) { uiOutput("protPep") },
      br(),
      style = paste0("background: ", myCol, ";")))
  } else {
    ## Coverage plots ###################################################
    if (tstCov) {
      dfltExp <- myExp[1L]
      exp2smpl <- listMelt(lapply(prots, \(pr) { myExp }), prots, ColNames = c("Sample", "Protein"))
      cov_plots <- lapply(1L:nrow(exp2smpl), \(i) {
        exp <- exp2smpl$Sample[i]
        pr <- exp2smpl$Protein[i]
        tags$div(id = paste0("cov_", pr, "_", exp),
                 style = paste("width: 100%; display: ",
                               if ((pr == dflt) && (exp == dfltExp)) { "block" } else { "none" },
                               ";"),
                 covPlots[[pr]]$logInt[[exp]])
      })
    }
    ## Ratio plots ######################################################
    ratio_plots_ui <- NULL
    if (tstRat) {
      prots2 <- intersect(prots, names(ratioPlots))
      if (length(prots2)) {
        ratio_plots_ui <- lapply(prots2, \(pr) {
          tags$div(id = paste0("rat_", pr),
                   style = paste("width: 100%; display: ",
                                 if (pr == dflt) { "block" } else { "display" },
                                 ";"),
                   ratioPlots[[pr]])
        })
      }
    }
    #
    ## Comments #########################################################
    prot_comments <- allComments[prots]
    prComments <- lapply(prots, \(pr) {
      make_comment_ui(pr,
                      FALSE,
                      prot_comments,
                      pr == dflt,
                      "prComment_")
    })
    #
    ## Peptide tables  ##################################################
    if (peptidoTst) {
      pepTables <- lapply(prots, \(pr) {
        m <- match(pr, prots)
        tags$div(id = paste0("pepTable_", m),
                 style = if (pr == dflt) { "width: 100%; display: block;" } else { "display: none;" },
                 make_smpl_tbl_ui(tab = "All peptidoforms",
                                  filt = pr))
      })
    }
    ## UI ###############################################################
    tagList(tags$div(
      if (length(prots) > 1L) {
        fluidRow(column(12L,
                        make_select_tag("myProtein",
                                        "",
                                        "myProtein",
                                        prots,
                                        dflt),
                        br()))
      },
      prComments,
      br(),
      fluidRow(column(6L,
                      if (length(myExp) > 1L) {
                        make_select_tag("mySample",
                                        "",
                                        "mySample",
                                        myExp,
                                        myExp[1L])
                      },
                      br(),
                      cov_plots),
               if (!is.null(ratio_plots_ui)) {
                 column(6L, ratio_plots_ui)
               },
      ),
      br(),
      br(),
      tags$hr(style = "border-color: black;"),
      if (peptidoTst) { pepTables },
      tags$script(HTML(paste0("function updateProteinTab() {
  const protEl = document.getElementById('myProtein');
  const sampleEl = document.getElementById('mySample');
  const singleSample = ",
                              jsonlite::toJSON(if (length(Exp) == 1L) { myExp[1L] } else { NULL },
                                               auto_unbox = TRUE),
                              ";
  // Nothing useful to update if there isn't even a protein
  if (!protEl) {
    return;
  }
  const prot = protEl.value;
  const ind = protEl.selectedIndex + 1;
  // If there is no sample selector, use the single sample
  // encoded by the first/only option, if available.
  const sample = sampleEl ? sampleEl.value : singleSample;
  const comm = 'prComment_' + ind;
  document.querySelectorAll('[id^=\"cov_\"]').forEach(function(el) {
    el.style.display = 'none';
  });
  if (sample !== null) {
    const cov = document.getElementById('cov_' + prot + '_' + sample);
    if (cov) {
      cov.style.display = 'block';
    }
  }
  document.querySelectorAll('[id^=\"pepTable_\"]').forEach(function(el) {
    el.style.display = 'none';
  });
  const pepTblID = document.getElementById('pepTable_' + ind);
  if (pepTblID) {
    pepTblID.style.width = '100%';
    pepTblID.style.display = 'block';
  }
  document.querySelectorAll('[id^=\"rat_\"]').forEach(function(el) {
    el.style.display = 'none';
  });
  const rat = document.getElementById('rat_' + prot);
  if (rat) {
    rat.style.display = 'block';
  }
  document.querySelectorAll('[id^=\"prComment_\"]').forEach(function(el) {
    el.style.display = 'none';
  });
  const comment = document.getElementById(comm);
  if (comment) {
    comment.style.display = 'block';
    comment.style.whiteSpace = 'pre-wrap';
    comment.style.padding = '10px';
  }
  window.dispatchEvent(new Event('resize'));
}
const mySample = document.getElementById('mySample');
if (mySample) {
  mySample.addEventListener('change', updateProteinTab);
}
const myProt = document.getElementById('myProtein');
if (myProt) {
  myProt.addEventListener('change', updateProteinTab);
}
"))),
      br(),
      style = paste0("background: ", myCol, ";")))
  }
}
make_summTbl_ui <- \() {
  df <- t(Exp_summary[, grep(" - % ", colnames(Exp_summary), invert = TRUE, value = TRUE)])
  colnames(df) <- df[1L,]
  df <- df[2L:nrow(df),]
  if (scrptType == "withReps") {
    tmp <- listMelt(Exp.map$MQ.Exp, Exp.map$Clean_name)
    m <- unlist(lapply(unlist(Exp.map$MQ.Exp), \(x) { which(Frac.map$MQ.Exp == x) }))
    df <- df[, c("Whole dataset", Frac.map$`Raw files name`[m])]
    df[1L, 2L:ncol(df)] <- tmp$L1[match(Frac.map$MQ.Exp[m], tmp$value)]
  } 
  #
  # Drop fixed PTMs (alkylation)
  fxdMods <- Modifs$`Full name`[Modifs$Type == "Fixed"]
  l <- length(fxdMods)
  if (l) {
    if (l > 1L) {
      fxdMods <- paste0("(", paste0("(", fxdMods, ")", collapse = "|"), ")")
    }
    pat <- paste0("^", fxdMods, " - ")
    df <- df[grep(fxdMods, rownames(df), invert = TRUE),]
  }
  #
  wdth <- paste0(160L*(ncol(df) + 1L) + 40*ncol(df), "px")
  #hght <- paste0(100L*(nrow(df)+1L), "px")
  rownames(df) <- sub("eptides$", "eptidoforms", rownames(df))
  rowHelp <- setNames(rep("", nrow(df)), rownames(df))
  rowHelp["PSMs"] <- "\"Peptide-Spectrum-Matches\" = individual identifications by the search engine"
  rowHelp["Peptidoforms"] <- "Peptides in a specific post-translationally-modified state"
  rowHelp["Protein groups"] <- "Groups of sequence-related proteins whose presence in the dataset is inferred from a collected of observed peptide sequences.\nA protein group includes:\n - one or more \"leading\" protein(s), which explain all peptide sequences assigned to the group and, if multiple, are indistinguishable based on observations,\n - ... as well as any other proteins which have no unique (= proteotypic) peptide but can produce some of the peptides in the group."
  df <- DT::datatable(df,
                      rownames = TRUE,
                      class = "compact",
                      escape = FALSE,
                      width = wdth,
                      #height = hght,
                      height = "auto",
                      fillContainer = FALSE,
                      options = list(rowCallback = JS(sprintf("function(row, data, displayNum, displayIndex, dataIndex) {
  const tips = %s;
  // With rownames = TRUE, the row name is usually in the first cell
  const rowNameCell = row.cells[0];
  const rowName = rowNameCell.textContent.trim();
  if (Object.prototype.hasOwnProperty.call(tips, rowName)) {
    rowNameCell.setAttribute('title', tips[rowName]);
  }
}",
                                                              jsonlite::toJSON(as.list(rowHelp), auto_unbox = TRUE))),
                                     scrollX = TRUE,
                                     paging = FALSE,
                                     lengthChange = FALSE,
                                     autoWidth = FALSE,
                                     dom = "ft",
                                     columnDefs = list(list(width = "160px",
                                                            targets = 1L:ncol(df) - 1L))))
  return(tags$div(df,
                  style = "background: #ffffff;"))
}
make_select_tag <- \(id,
                     label,
                     name,
                     values,
                     selected) {
  tags$div(if (nchar(label)) {
    tags$label(`for` = id,
               label)
  },
  tags$select(id = id,
              name = name,
              `data-default` = selected,
              lapply(values, \(x) {
                tags$option(value = x,
                            selected = if (x == selected) { "selected" } else { NULL },
                            x)
              })))
}
make_smpl_tab <- \(exp,
                   shiny = TRUE,
                   quant = quantMeth,
                   dflt = dfltQuant) {
  myCol <- tolower(viridis::viridis(6L, alpha = 0.2)[2L])
  exp2 <- if (scrptType == "noReps") { exp } else { names(smplGrps)[match(exp, smplGrps)] }
  lQ <- length(quant)
  if (shiny) {
    tagList(tags$div(
      make_comment_ui(exp, shiny),
      selectInput(paste0("quant_", exp), "", quant, dflt[exp]),
      plotlyOutput(paste0("quantLy_", exp), height = plotHght),
      br(),
      br(),
      tags$hr(style = "border-color: black;"),
      make_smpl_tbl_ui(exp),
      style = paste0("background: ", myCol, ";")))
  } else {
    id1 <- paste0("quant_", exp)
    id2 <- paste0("quant_", exp, "_")
    js <- sprintf("document.getElementById('%s').addEventListener('change', function() {
  const selected = document.getElementById('%s').selectedIndex + 1;
  document.querySelectorAll('[id^=\"%s\"]').forEach(function(div) {
    div.style.display = 'none';
  });
  document.getElementById('%s' + selected).style.display = 'block';
  window.dispatchEvent(new Event('resize'));
});",
                  id1,
                  id1,
                  id2,
                  id2)
    tagList(tags$div(
      make_comment_ui(exp, shiny),
      lapply(1L:lQ, \(i) {
        tags$div(id = paste0("quant_", exp, "_", as.character(i)),
                 style = if (quant[i] == dflt[exp]) { "display: block;" } else { "display: none;" },
                 ggQuantLy[[quant[i]]][[exp2]]$plotly)
      }),
      br(),
      make_select_tag(id1,
                      "",
                      id1,
                      quant,
                      dflt[exp]),
      br(),
      br(),
      tags$hr(style = "border-color: black;"),
      tags$script(HTML(js)),
      make_smpl_tbl_ui(exp),
      style = paste0("background: ", myCol, ";")))
  }
}
make_ctrst_tab <- \(contr,
                    shiny = TRUE) {
  myCol <- tolower(viridis::viridis(6L, alpha = 0.2)[3L])
  styleOn <- paste0("display: block; height: ", plotHtMpHght)
  contr2 <- gsub(" ", "_", contr)
  saintIDs <- c(paste0("SAINTexpress volcano plot ", contr),
                paste0(contr2, c("_SAINT_volcPlot", "_SAINT_GObars")))
  if (runGSEA) {
    GSEA_IDs <- paste0(contr2, "_GSEA", as.character(1L:4L))
  }
  saintXPRS <- saintExprs && (saintIDs[1L] %in% names(volcPlotly$SAINTexpress))
  if (shiny) {
    tagList(tags$div(
      make_comment_ui(contr, shiny),
      if (saintXPRS) {
        div(h3("SAINTexpress"),
            fluidRow(column(6L,
                            plotlyOutput(paste0(contr2, "_SAINT_volcPlot"), height = "600px")),
                     if (enrichGO) {
                       column(6L,
                              br(),
                              br(),
                              plotlyOutput(paste0(contr2, "_SAINT_GObars"), height = "600px"))
                     },
            ),
            style = "background: #ffffff;")
      },
      if (!saintXPRS) {
        div(h3("t-test"),
            fluidRow(column(6L,
                            plotlyOutput(paste0(contr2, "_volcPlot"), height = "600px")),
                     if (enrichGO) {
                       column(6L,
                              br(),
                              br(),
                              plotlyOutput(paste0(contr2, "_GObars"), height = "600px"))
                     },
            ),
            style = "background: #ffffff;")
      },
      tags$hr(style = "border-color: black;"),
      if (F.test) {
        # Add F-test part here... or maybe dropdown to choose f-/F-test... or drop F-test altogether?
      },
      if (runGSEA) {
        div(
          div(h3("GSEA"),
              fluidRow(column(6L,
                              plotlyOutput(GSEA_IDs[1L]),
                              plotlyOutput(GSEA_IDs[2L])),
                       column(6L,
                              plotlyOutput(GSEA_IDs[3L]),
                              plotlyOutput(GSEA_IDs[4L]))),
              style = "background: #ffffff;"),
          tags$hr(style = "border-color: black;"))
      },
      make_ctrst_tbl_ui(contr),
      style = paste0("background: ", myCol, ";")))
  } else {
    styleOn6 <- "display: block; height: 600px"
    styleOn4 <- "display: block; height: 400px"
    tagList(tags$div(
      make_comment_ui(contr, shiny),
      if (saintXPRS) {
        div(h3("SAINTexpress"),
            fluidRow(column(6L,
                            tags$div(id = saintIDs[2L],
                                     style = styleOn6,
                                     volcPlotly$SAINTexpress[[saintIDs[1L]]]$Plot)),
                     if (enrichGO) {
                       column(6L,
                              br(),
                              br(),
                              tags$div(id = saintIDs[3L],
                                       style = styleOn6,
                                       GO_plot_ly$Prot$SAINTexpress[[contr]]$Bar))
                     },
            ),
            style = "background: #ffffff;")
      },
      if (!saintXPRS) {
        div(h3("t-test"),
            fluidRow(column(6L,
                            tags$div(id = paste0(contr2, "_volcPlot"),
                                     style = styleOn6,
                                     volcPlotly$"t-test"[[paste0("Volcano plot ", contr)]]$Plot)),
                     if (enrichGO) {
                       column(6L,
                              br(),
                              br(),
                              tags$div(id = paste0(contr2, "_GObars"),
                                       style = styleOn6,
                                       GO_plot_ly$PG$"t-test"[[contr]]$Bar))
                     },
            ),
            style = "background: #ffffff;")
      },
      tags$hr(style = "border-color: black;"),
      if (F.test) {
        # Add F-test part here... or maybe dropdown to choose f-/F-test... or drop F-test altogether?
      },
      if (runGSEA) {
        div(
          div(h3("GSEA"),
              # NB: I also tried the plotly::subplot() approach to displaying the plots together in one,
              # but this fails (subplots look corrupted, possibly because they are slightly hacky)
              fluidRow(column(6L,
                              tags$div(id = GSEA_IDs[1L],
                                       style = styleOn4,
                                       GSEA_plots$standard$PG$`GSEA dotplot`[[contr]]),
                              tags$div(id = GSEA_IDs[2L],
                                       style = styleOn4,
                                       GSEA_plots$standard$PG$`GSEA enrichment map`[[contr]])),
                       column(6L,
                              tags$div(id = GSEA_IDs[3L],
                                       style = styleOn4,
                                       GSEA_plots$standard$PG$`GSEA ridge plot`[[contr]]),
                              tags$div(id = GSEA_IDs[4L],
                                       style = styleOn4,
                                       GSEA_plots$standard$PG$`GSEA category net plot`[[contr]]))),
              style = "background: #ffffff;"),
          tags$hr(style = "border-color: black;"))
      },
      make_ctrst_tbl_ui(contr),
      style = paste0("background: ", myCol, ";")))
  }
}
make_strt_tab <- \(shiny = TRUE) {
  myCol <- tolower(viridis::viridis(6L, alpha = 0.2)[1L])
  # TO DO
  # Add Dataset GO terms enrichment
  if (shiny) {
    tagList(tags$div(
      make_comment_ui("Dataset overview", shiny),
      br(),
      h4(strong(tags$ul(em("Summary table")))),
      make_summTbl_ui(),
      br(),
      br(),
      if (heatMaps_ON) {
        column(12L,
               selectInput("myHeatMap",
                           "",
                           nmsHtMp,
                           nmsHtMp[1L]),
               plotlyOutput("heatMap", height = plotHtMpHght))
      },
      if (globalGO && ("Observed dataset" %in% names(GO_plot_ly))) {
        fluidRow(column(12L,
                        plotlyOutput("GO_enrich_Dataset", height = plotHtMpHght)))
      },
      fluidRow(
        if (PCA_ON) {
          column(4L*(3L-Venn_ON),
                 plotlyOutput("PCA", height = plotPCAHght))
        },
        if (Venn_ON) {
          column(4L,
                 plotlyOutput("Venn", height = plotHtMpHght))
        },
      ),
      br(),
      style = paste0("background: ", myCol, ";")))
  } else {
    styleOn <- paste0("display: block; height: ", plotHtMpHght)
    tagList(tags$div(
      make_comment_ui("Dataset overview", shiny),
      br(),
      h4(strong(tags$ul(em("Summary table")))),
      make_summTbl_ui(),
      br(),
      br(),
      if (heatMaps_ON) {
        fluidRow(column(strtColWdth,
                        make_select_tag("myHeatMap",
                                        "",
                                        "myHeatMap",
                                        nmsHtMp,
                                        nmsHtMp[1L]),
                        lapply(nmsHtMp, \(nm) {
                          i <- match(nm, nmsHtMp)
                          tags$div(id = paste0("HeatMap_", i),
                                   style = if (i == 1L) { styleOn } else { "display: none;" },
                                   plotLeatMaps$Global[[nm]]$Plot)
                        })))
      },
      if (globalGO && ("Observed dataset" %in% names(GO_plot_ly))) {
        fluidRow(column(12L,
                        tags$div(id = "GO_enrich_Dataset",
                                 style = if (i == 1L) { styleOn } else { "display: none;" },
                                 GO_plot_ly$`Observed dataset`$Bar)))
      },
      fluidRow(
        if (PCA_ON) {
          column(4L*(3L-Venn_ON),
                 tags$div(id = "PCA",
                          style = styleOn,
                          dimRedPlotLy$PG$PCA))
        },
        if (Venn_ON) {
          column(4L,
                 tags$div(id = "Venn",
                          style = styleOn,
                          plotly_Venn$`Global, LFQ`))
        },
      ),
      br(),
      tags$script(HTML(paste0("document.getElementById('myHeatMap').addEventListener('change', function() {
  const HtMpID = document.getElementById('myHeatMap').selectedIndex + 1;
  const HtMp = document.getElementById('HeatMap_' + HtMpID);
  document.querySelectorAll('[id^=\"HeatMap_\"]').forEach(function(el) {
    el.style.display = 'none';
  });
  HtMp.style.display = 'block';
  HtMp.style.height = '", plotHtMpHght, "';
  window.dispatchEvent(
    new Event('resize')
  );
});"))),
      style = paste0("background: ", myCol, ";")))
  }
}
make_QC_tab <- \(shiny = TRUE,
                 plotsList = QC_plotLys) {
  myCol <- tolower(viridis::viridis(6L, alpha = 0.2)[5L])
  if (shiny) {
    tagList(tags$div(
      selectInput("QC1", "", names(plotsList), names(plotsList)[1L]),
      fluidRow(column(8L,
                      plotlyOutput("QCplotLy", height = plotHght)),
               column(4L,
                      uiOutput("QCtxt"))),
      br(),
      style = paste0("background: ", myCol, ";")))
  } else {
    QC_comments <- allComments[names(plotsList)]
    tagList(tags$div(
      make_select_tag("myQC",
                      "",
                      "myQC",
                      names(plotsList),
                      names(plotsList)[1L]),
      lapply(seq_along(plotsList), \(i) { #i <-1L #i <- i+1L
        nm <- names(plotsList)[i]
        fluidRow(column(8L,
                        tags$div(id = paste0("QC_", i),
                                 style = if (i == 1L) { "display: block;" } else { "display: none;" },
                                 plotsList[[nm]])),
                 column(4L,
                        make_comment_ui(nm,
                                        FALSE,
                                        QC_comments,
                                        nm == names(plotsList)[1L],
                                        "QCcomment_")))
      }),
      br(),
      tags$script(HTML("document.getElementById('myQC').addEventListener('change', function() {
  const selected = document.getElementById('myQC').selectedIndex + 1;
  document.querySelectorAll('[id^=\"QC_\"]').forEach(function(div) { div.style.display = 'none'; });
  document.getElementById('QC_' + selected).style.display = 'block';
  document.querySelectorAll('[id^=\"QCcomment_\"]').forEach(function(div) { div.style.display = 'none'; });
  const div = document.getElementById('QCcomment_' + selected);
  div.style.display = 'block';
  div.style.whiteSpace = 'pre-wrap';
  div.style.padding = '10px';
  window.dispatchEvent(new Event('resize'));
});")),
      style = paste0("background: ", myCol, ";")))
  }
}
make_matmet_tab <- \(matmeth = matmethTxt,
                     shiny = TRUE) {
  myCol <- tolower(viridis::viridis(6L, alpha = 0.2)[6L])
  # We want to load the processed materials and methods (potentially edited by the user)
  # ============> This should be ideally run as part of the finalization script, after the materials and method edition stage
  #
  hght <- vapply(strsplit(matmeth, "\n"), \(x) {
    paste0(as.character(20L*(sum(ceiling(nchar(unlist(x))/ceiling(screenRes$width/5.75)))+2L)), "px")
  }, "")
  if (shiny) {
    tagList(tags$div(
      textAreaInput("MatMet_SamplePrep", matmethSections[1L], matmeth[1L], "100%", hght[1L]),
      br(),
      textAreaInput("MatMet_LCMS", matmethSections[2L], matmeth[2L], "100%", hght[2L]),
      br(),
      textAreaInput("MatMet_DataAnalysis", matmethSections[3L], matmeth[3L], "100%", hght[3L]),
      br(),
      style = paste0("background: ", myCol, ";")))
  } else {
    tagList(tags$div(
      h5(matmethSections[1L]),
      tags$p(matmeth[1L]),
      br(),
      h5(matmethSections[2L]),
      tags$p(matmeth[2L]),
      br(),
      h5(matmethSections[3L]),
      tags$p(matmeth[3L]),
      br(),
      style = paste0("background: ", myCol, ";")))
  }
}
make_ui_noReps <- \(tabNames = myTabs,
                    shiny = TRUE) {
  tabs <- lapply(tabNames, \(x) {
    if (x == "Dataset overview") {
      return(tabPanel(x,
                      make_strt_tab(shiny = shiny)))
    }
    if (x %in% Exp) {
      return(tabPanel(paste0("sample = ", x),
                      make_smpl_tab(x,
                                    shiny = shiny)))
    }
    if (x == "Proteins of interest") {
      return(tabPanel(x,
                      make_prot_tab(dfltProt,
                                    shiny = shiny)))
    }
    if (x == "QC") {
      return(tabPanel(x,
                      make_QC_tab(shiny = shiny)))
    }
    if (x == "Materials and methods") {
      return(tabPanel(x,
                      make_matmet_tab(shiny = shiny)))
    }
  })
  return(bslib::navset_tab(!!!tabs))
}
make_ui_Reps <- \(tabNames = myTabs,
                  shiny = TRUE) {
  tabs <- lapply(tabNames, \(x) {
    if (x == "Dataset overview") {
      return(tabPanel(x,
                      make_strt_tab(shiny = shiny)))
    }
    if (x %in% smplGrps) {
      return(tabPanel(paste0("sample group = ", x),
                      make_smpl_tab(x,
                                    shiny = shiny)))
    }
    if (x %in% myContrasts$Contrast) {
      return(tabPanel(paste0("contrast = ", x),
                      make_ctrst_tab(contr = x,
                                     shiny = shiny)))
    }
    if (x == "Proteins of interest") {
      return(tabPanel(x,
                      make_prot_tab(dfltProt,
                                    shiny = shiny)))
    }
    if (x == "QC") {
      return(tabPanel(x,
                      make_QC_tab(shiny = shiny)))
    }
    if (x == "Materials and methods") {
      return(tabPanel(x,
                      make_matmet_tab(shiny = shiny)))
    }
  })
  return(bslib::navset_tab(!!!tabs))
}
make_ui <- if (scrptType == "noReps") { make_ui_noReps } else { make_ui_Reps }

# Plot HTML paths
#myPlots <- list.files(paste0(wd, "/Ranked abundance/LFQ"), "\\.html$", full.names = TRUE)
#names(myPlots) <- gsub(".* - |\\.html$", "", myPlots)
#nPl <- length(myPlots)
myTabs <- nms <- if (scrptType == "noReps") {
  c("Dataset overview", Exp)
} else {
  c("Dataset overview", smplGrps, myContrasts$Contrast)
}
if (prot.list.Cond) {
  if (tstCov) {
    allProt <- names(covPlots)[vapply(names(covPlots), \(x) { length(covPlots[[x]]$logInt) > 0L }, TRUE)]
    dfltProt <- allProt[1L]
    nms <- union(nms, allProt)
  } else {
    allProt <- do.call(paste, c(db[match(prot.list, db$`Protein ID`), c("Protein ID", "Common Name")], sep = "_"))
    dfltProt <- allProt[1L]
  }
  myTabs <- union(myTabs, "Proteins of interest")
} else {
  dfltProt <- c()
}
myTabs <- union(myTabs, c("QC", "Materials and methods"))
nms <- union(nms, c("QC", names(QC_plotLys)))
dfltComment <- paste0(nrow(PG), " protein groups were identified from ", nrow(ev), " PSMs", " corresponding to ", nrow(pep),
                      " distinct peptidoforms. ...")
if ((!exists("allComments")) || (!is.character(allComments))) {
  allComments <- setNames(vapply(nms, \(nm) {
    if (nm == "Dataset overview") {
      dfltComment
    } else { "" }
  }, ""), nms)
}
nms_ <- setdiff(nms, names(allComments))
if (length(nms_)) { # Generate defaults
  allComments[nms_] <- ""
  if ("Dataset overview" %in% nms_) { allComments$"Dataset overview" <- dfltComment }
}
allComments %<o% allComments[nms]

#
quantLst <- setNames(lapply(names(ggQuantLy), \(tp) { names(ggQuantLy[[tp]]) }), names(ggQuantLy))
quantLst <- listMelt(quantLst, ColNames = c("Sample", "Type"))
quantLst <- aggregate(quantLst$Type, list(quantLst$Sample), list)
colnames(quantLst) <- c("Sample", "Types")
if (scrptType == "withReps") {
  quantLst$Sample <- cleanNms(quantLst$Sample)
}
dfltQuant <- quantLst[, "Sample", drop = FALSE]
dfltQuant$Type <- vapply(quantLst$Types, \(x) { x[[1L]] }, "")
quantLst <- setNames(quantLst$Types, quantLst$Sample)
dfltQuant <- setNames(dfltQuant$Type, dfltQuant$Sample)
quantMeth <- unique(unlist(quantLst))
#
appPage <- 1L
appNm <- "Edit report"
ui <- fluidPage(useShinyjs(),
                extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
                tags$head(tbl_css),
                titlePanel(tag("u", appNm),
                           appNm),
                br(),
                fluidRow(column(4L,
                                h2(dtstNm),
                                br()),
                         column(8L,
                                actionBttn("xprtBtn", " export final html report", icon = icon("file-export"), color = "success", style = "pill"),
                                br(),
                                br(),
                                uiOutput("xprtMsg"))),
                br(),
                uiOutput("myUI"),
                br(),
                br())
server <- \(input, output, session) {
  # if (prot.list.Cond) {
  #   PROT <- reactiveVal(dfltProt)
  # }
  myExp <- if (scrptType == "noReps") { Exp } else { setNames(smplGrps, NULL) }
  QUANT <- reactiveVal(dfltQuant)
  XPRTMSG <- reactiveVal(NULL)
  MYPROT <- reactiveVal(dfltProt)
  SAMPLE <- reactiveVal(myExp[1L])
  NORMMETH <- reactiveVal("None")
  # Render UI
  output$xprtMsg <- renderUI(XPRTMSG())
  output$myUI <- renderUI(make_ui())
  if (heatMaps_ON) {
    output$heatMap <- renderPlotly(plotLeatMaps$Global[[NORMMETH()]]$Plot)
  }
  if (PCA_ON) {
    output$PCA <- renderPlotly(dimRedPlotLy$PG$PCA)
  }
  if (Venn_ON) {
    output$Venn <- renderPlotly(plotly_Venn$`Global, LFQ`)
  }
  #
  lapply(myExp, \(exp) {
    exp2 <- if (scrptType == "noReps") { exp } else { names(smplGrps)[match(exp, smplGrps)] }
    idQ <- paste0("quant_", exp)
    idQLy <- paste0("quantLy_", exp)
    output[[idQLy]] <- renderPlotly(ggQuantLy[[input[[idQ]]]][[exp2]]$plotly)
  })
  #
  if (scrptType == "withReps") {
    lapply(myContrasts$Contrast, \(contr) {
      # Don't use a for loop here! In absence of a reactive component to how we access the plotly plot in the list,
      # this would display the plots from the last contrast in all contrast tabs!
      contr2 <- gsub(" ", "_", contr)
      output[[paste0(contr2, "_volcPlot")]] <- renderPlotly(volcPlotly$"t-test"[[paste0("Volcano plot ", contr)]]$Plot)
      if (enrichGO) {
        output[[paste0(contr2, "_GObars")]] <- renderPlotly(GO_plot_ly$PG$"t-test"[[contr]]$Bar)
      }
      if (F.test) {
        # Add F-test part here... or maybe dropdown to choose f-/F-test... or drop F-test altogether?
      }
      if (runGSEA) {
        GSEA_IDs <- paste0(contr2, "_GSEA", as.character(1L:4L))
        output[[GSEA_IDs[1L]]] <- renderPlotly(GSEA_plots$standard$PG$`GSEA dotplot`[[contr]])
        output[[GSEA_IDs[2L]]] <- renderPlotly(GSEA_plots$standard$PG$`GSEA enrichment map`[[contr]])
        output[[GSEA_IDs[3L]]] <- renderPlotly(GSEA_plots$standard$PG$`GSEA ridge plot`[[contr]])
        output[[GSEA_IDs[4L]]] <- renderPlotly(GSEA_plots$standard$PG$`GSEA category net plot`[[contr]])
      }
      saintIDs <- c(paste0("SAINTexpress volcano plot ", contr),
                    paste0(contr2, c("_SAINT_volcPlot", "_SAINT_GObars")))
      if (saintExprs && (saintIDs[1L] %in% names(volcPlotly$SAINTexpress))) {
        output[[saintIDs[2L]]] <- renderPlotly(volcPlotly$SAINTexpress[[saintIDs[1L]]]$Plot)
        if (enrichGO) {
          output[[saintIDs[3L]]] <- renderPlotly(GO_plot_ly$Prot$SAINTexpress[[contr]]$Bar)
        }
      }
    })
  }
  if (globalGO && (!is.null(GO_plot_ly$PG$Dataset$`Observed dataset`$Bar))) {
    output$GO_enrich_Dataset <- renderPlotly(GO_plot_ly$PG$Dataset$`Observed dataset`$Bar)
  }
  #
  # Event observers
  observeEvent(input$myHeatMap, { NORMMETH(input$myHeatMap) })
  #  - Comments
  sapply(names(allComments), \(nm) {
    observeEvent(input[[paste0("comment_", nm)]], {
      allComments[[nm]] <- input[[paste0("comment_", nm)]]
      allComments <<- allComments
    })
  })
  #  - Quant method
  sapply(myExp, \(exp) {
    exp2 <- if (scrptType == "noReps") { exp } else { names(smplGrps)[match(exp, smplGrps)] }
    idQ <- paste0("quant_", exp)
    idQLy <- paste0("quantLy_", exp)
    observeEvent(input[[idQ]], {
      dfltQuant <- QUANT()
      dfltQuant[exp] <- input[[idQ]]
      QUANT(dfltQuant)
      assign("dfltQuant", dfltQuant, envir = .GlobalEnv)
      # Update plot
      output[[idQLy]] <- renderPlotly(ggQuantLy[[input[[idQ]]]][[exp2]]$plotly)
    })
  })
  #  - Proteins tab
  if (prot.list.Cond) {
    if (scrptType == "noReps") {
      output$ratioPlot <- renderPlotly(ratioPlots[[MYPROT()]])
    }
    output$coverPlot <- renderPlotly({
      p <- covPlots[[MYPROT()]]$logInt[[SAMPLE()]]
      if (is.null(p)) {
        return(plot_ly(type = "scatter",
                       mode = "markers") |>
                 layout(xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        annotations = list(list(text = "No identifications for this protein in this sample!",
                                                x = 0.5,
                                                y = 0.5,
                                                xref = "paper",
                                                yref = "paper",
                                                showarrow = FALSE))))
      }
      return(p)
    })
    output$protComment <- renderUI(make_comment_ui(MYPROT()))
    if (peptidoTst) {
      output$protPep <- renderUI({
        make_smpl_tbl_ui(tab = "All peptidoforms",
                         filt = MYPROT())
      })  
    }
    if (tstCov && (length(allProt) > 1L)) {
      observeEvent(input$myProtein, { MYPROT(input$myProtein) })
    }
    if (length(myExp) > 1L) {
      observeEvent(input$mySample, {
        SAMPLE(input$mySample)
      })
    }
  }
  #  - QC tab
  observeEvent(input$QC1, {
    output$QCplotLy <- renderPlotly(QC_plotLys[[input$QC1]])
    output$QCtxt <- renderUI(make_comment_ui(input$QC1))
  })
  #  - Materials and methods
  observeEvent(input$MatMet_SamplePrep, {
    txt <- matmethTxt
    txt[matmethSections[1L]] <- input$MatMet_SamplePrep
    assign("matmethTxt", txt, envir = .GlobalEnv)
  })
  observeEvent(input$MatMet_LCMS, {
    txt <- matmethTxt
    txt[matmethSections[2L]] <- input$MatMet_LCMS
    assign("matmethTxt", txt, envir = .GlobalEnv)
  })
  observeEvent(input$MatMet_DataAnalysis, {
    txt <- matmethTxt
    txt[matmethSections[3L]] <- input$MatMet_DataAnalysis
    assign("matmethTxt", txt, envir = .GlobalEnv)
  })
  #  - Render final report
  observeEvent(input$xprtBtn, {
    XPRTMSG(em("Exporting .html report, this will take a few seconds...",
               style = "color:green",
               .noWS = "outside"))
    later::later(\() {
      # Wrapping in this allows displaying the message before export completes
      # 1. Rebuild the SAME UI we use in the app
      page <- bslib::page_fluid(tags$head(tbl_css),
                                tags$script(HTML("document.addEventListener('DOMContentLoaded', function() {
  document.querySelectorAll('select[data-default]').forEach(function(sel) {
    sel.value = sel.dataset.default;
    sel.dispatchEvent(new Event('change'));
  });
});")),
                                report_header,
                                make_ui(shiny = FALSE))
      # 2. Wrap as browsable HTML
      page <- htmltools::browsable(page)
      # 3. Save to disk
      htmltools::save_html(page, htmlRprtFl)
      assign("appRunTst", TRUE, envir = .GlobalEnv)
      stopApp()
    }, 0.1)
  })
  session$onSessionEnded(\() { stopApp() })
}
#eval(parse(text = run_App), envir = .GlobalEnv)
runKount <- 0L
if (exists("appRunTst")) { rm(appRunTst) }
while ((!runKount) || (!exists("appRunTst")) || (!file.exists(htmlRprtFl))) {
  eval(parse(text = run_App), envir = .GlobalEnv)
  shinyCleanup()
  runKount <- runKount + 1L
}

# We now have our html... but it depends on local libraries...
# ---> We want those embedded in it so it is fully portable!
h2 <- h1 <- readr::read_lines(htmlRprtFl)
rg1 <- grep("</?head>", h1) + c(1L, -1L)
rg1 <- rg1[1L]:rg1[2L]
hd1 <- h1[rg1]
hd1 <- data.frame(original = hd1)
hd1$new <- hd1$original
g <- grep("^ *<((style)|(script)|(link))( *[^>]+)?>", hd1$original)
hd1$original[g]
require(base64enc)
read_file <- \(path) {
  paste(readr::read_lines(path, warn = FALSE), collapse = "\n")
}
file_to_data_uri <- \(path) {
  ext <- tools::file_ext(path)
  mime <- switch(tolower(ext),
                 "woff2" = "font/woff2",
                 "woff"  = "font/woff",
                 "ttf"   = "font/ttf",
                 "png"   = "image/png",
                 "jpg"   = "image/jpeg",
                 "jpeg"  = "image/jpeg",
                 "svg"   = "image/svg+xml",
                 "gif"   = "image/gif",
                 "application/octet-stream")
  paste0("data:",
         mime,
         ";base64,",
         base64enc::base64encode(path))
}
# - embed scripts
read_asset <- \(path) {
  readChar(path, # Do not use readLines, which isn't binary-safe!
           nchars = file.info(path)$size,
           useBytes = TRUE)
}
inline_script <- \(path) {
  txt <- paste(read_asset(path), collapse = "")
  txt <- gsub("</script",
              "<\\/script",
              txt,
              ignore.case = TRUE)
  paste0("<script>\n",
         txt,
         "\n</script>")
}
gs <- grep("^ *<script src=\"", hd1$original)
hd1$new[gs] <- vapply(sub("\".*", "", sub("^ *<script src=\"", paste0(wd, "/"), hd1$original[gs])), inline_script, "")
# - embed css
inline_css <- \(path) {
  paste0(  "<style>\n",
           read_asset(path),
           "\n</style>")
}
gc <- grepl("^ *<link href=\"", hd1$original)
hd1$new[gc] <- vapply(sub("\".*", "", sub("^ *<link href=\"", paste0(wd, "/"), hd1$original[gc])), inline_css, "")
h2[rg1] <- hd1$new
write(h2, htmlRprtFl)
removeDirectory(paste0(wd, "/lib"), TRUE, FALSE)

# Write Mat Meth template as separate file
MatMetCalls$Texts$WetLab <- matmethTxt["Samples preparation"]
MatMetCalls$Texts$LCMS <- matmethTxt["LC-MS/MS analysis"]
MatMetCalls$Texts$DatAnalysis <- matmethTxt["Data analysis"]
setwd(wd)
tmp <- paste0("MatMet <- ", unlist(MatMetCalls$Calls))
tmpSrc <- paste0(wd, "/tmp.R")
write(tmp, tmpSrc)
MatMetFl <- paste0(wd, "/Materials and methods_WIP.docx")
tst <- try({
  source(tmpSrc)
  #rstudioapi::documentOpen(tmpSrc)
  MatMet %<o% MatMet
  print(MatMet, target = MatMetFl)
}, silent = TRUE)
if (inherits(tst, "try-error")) {
  warning("Couldn't write materials and methods template, investigate...")
}
unlink(tmpSrc)

try({
  rm(ggQuantLy,
     plotLeatMaps,
     dimRedPlotLy,
     plotly_Venn,
     QC_plotLys)
}, silent = TRUE)
