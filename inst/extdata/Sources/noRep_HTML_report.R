library(shiny)
library(shinyjs)
library(bslib)
library(htmltools)
library(DT)
library(plotly)
library(jsonlite)

htmlRprtFl <- paste0(wd, "/Report_", dtstNm, ".html")

# Reload processed data from the report
allNms <- openxlsx2::wb_get_sheet_names(openxlsx2::wb_load(repFl))
nms <- setdiff(allNms, c("Description", "Quality control"))
xlDat <- setNames(lapply(nms, \(nm) {
  dat <- openxlsx2::read_xlsx(repFl, match(nm, allNms), 2L)
  if ("Potential contaminant" %in% colnames(dat)) {
    w <- which(is.na(dat$"Potential contaminant"))
    if (length(w)) { dat$"Potential contaminant"[w] <- "" }
  }
  return(dat)
}), nms)
#
# Reload materials and methods
matmethSections <- c("Samples preparation", "LC-MS/MS analysis", "Data analysis")
fl <- paste0(wd, "/Materials and methods_WIP.docx")
if (!exists("matmethTxt")) {
  matmethTxt <- setNames(rep("", 3L),
                         matmethSections)
}
if (file.exists(fl)) {
  matmethTxt_fromFl <- officer::read_docx(fl)
  matmethTxt_fromFl <- officer::docx_summary(matmethTxt_fromFl)
  w <- match(matmethSections, matmethTxt_fromFl$text)
  matmethTxt_fromFl <- setNames(vapply(1L:3L, \(i) {
    rg <- (w[i]+1L):(c(w, nrow(matmethTxt_fromFl)+1L)[i+1L]-1L)
    paste(setdiff(matmethTxt_fromFl$text[rg], ""), collapse = "\n")
  }, ""), matmethSections)
  w <- which(matmethTxt == "")
  if (length(w)) {
    matmethTxt[w] <- matmethTxt_fromFl[w]
  }
}
#
tstRat <- MakeRatios && exists("ratioPlots") && (length(ratioPlots) > 0L)
# If necessary reload plots data
flHtMp <- paste0(wd, "/Clustering/HeatMaps.RData")
tstHtMp <- file.exists(flHtMp)
if (tstHtMp) {
  loadFun(flHtMp)
  tstHtMp <- exists("plotLeatMaps") && length(plotLeatMaps)
}
loadFun(paste0(wd, "/Sorting plots/quantPlots.RData"))
flPCA <- paste0(wd, "/Dimensionality red. plots/DimRedPlots.RData")
tstPCA <- file.exists(flPCA)
if (tstPCA) {
  loadFun(flPCA)
  tstPCA <- exists("dimRedPlotLy") && ("PCA" %in% names(dimRedPlotLy))
}
flVenn <- paste0(wd, "/Venn diagrams/Venn_plotly.RDS")
tstVenn <- file.exists(flVenn)
if (tstVenn) {
  loadFun(flVenn)
  tstVenn <- exists("plotly_Venn") && ("Global, LFQ" %in% names(plotly_Venn))
}
strtColWdth <- 12L/(tstVenn + tstPCA)

# Fix to plotly autoscaling + remove some Modebar tools (redundant: they should already be gone, but in case we reload old data)
if (tstHtMp) {
  for (x in names(plotLeatMaps)) { #x <- names(plotLeatMaps)[1L]
    for (y in names(plotLeatMaps[[x]])) { #y <- names(plotLeatMaps[[x]])[1L]
      plotLeatMaps[[x]][[y]]$Render$x$layout$xaxis$autorange <- TRUE
      plotLeatMaps[[x]][[y]]$Render$x$layout$yaxis$autorange <- TRUE
      plotLeatMaps[[x]][[y]]$Render <- plotly::config(plotLeatMaps[[x]][[y]]$Render,
                                                      modeBarButtonsToRemove = c("select2d", "lasso2d"))
    }
  }
}
if (!exists("ggQuantLy")) { stop("Ugh... really? No ranked abundance plots? Shouldn't we always have those by now?") }
for (x in names(ggQuantLy)) { #x <- names(ggQuantLy)[1L]
  for (y in names(ggQuantLy[[x]])) { #y <- names(ggQuantLy[[x]])[1L]
    ggQuantLy[[x]][[y]]$plotly$x$layout$xaxis$autorange <- TRUE
    ggQuantLy[[x]][[y]]$plotly$x$layout$yaxis$autorange <- TRUE
    ggQuantLy[[x]][[y]]$plotly <- plotly::config(ggQuantLy[[x]][[y]]$plotly,
                                                 modeBarButtonsToRemove = c("select2d", "lasso2d"))
  }
}
if (tstPCA) {
  for (x in names(dimRedPlotLy)) { #x <- names(dimRedPlotLy)[1L]
    dimRedPlotLy[[x]]$x$layout$xaxis$autorange <- TRUE
    dimRedPlotLy[[x]]$x$layout$yaxis$autorange <- TRUE
    dimRedPlotLy[[x]] <- plotly::config(dimRedPlotLy[[x]],
                                        modeBarButtonsToRemove = c("select2d", "lasso2d"))
  }
}
if (exists("covPlots")) {
  for (x in names(covPlots)) { #x <- names(covPlots)[1L]
    for (y in names(covPlots[[x]])) { #y <- names(covPlots[[x]])[1L]
      for (z in names(covPlots[[x]][[y]])) { #z <- names(covPlots[[x]][[y]])[1L]
        covPlots[[x]][[y]][[z]]$x$layout$xaxis$autorange <- TRUE
        covPlots[[x]][[y]][[z]]$x$layout$yaxis$autorange <- TRUE
        covPlots[[x]][[y]][[z]] <- plotly::config(covPlots[[x]][[y]][[z]],
                                                  modeBarButtonsToRemove = c("select2d", "lasso2d"))
      }
    }
  }
}
if (tstRat) {
  for (x in names(ratioPlots)) { #x <- names(ratioPlots)[1L]
    ratioPlots[[x]]$x$layout$xaxis$autorange <- TRUE
    ratioPlots[[x]]$x$layout$yaxis$autorange <- TRUE
    ratioPlots[[x]] <- plotly::config(ratioPlots[[x]] ,
                                      modeBarButtonsToRemove = c("select2d", "lasso2d"))
  }
}
for (x in names(QC_plotLys)) { #x <- names(QC_plotLys)[1L]
  QC_plotLys[[x]]$x$layout$xaxis$autorange <- TRUE
  QC_plotLys[[x]]$x$layout$yaxis$autorange <- TRUE
  QC_plotLys[[x]] <- plotly::config(QC_plotLys[[x]] ,
                                    modeBarButtonsToRemove = c("select2d", "lasso2d"))
}
if (tstVenn) {
  for (x in names(plotly_Venn)) { #x <- names(plotly_Venn)[1L]
    plotly_Venn[[x]]$x$layout$xaxis$autorange <- TRUE
    plotly_Venn[[x]]$x$layout$yaxis$autorange <- TRUE
    plotly_Venn[[x]] <- plotly::config(plotly_Venn[[x]] ,
                                       modeBarButtonsToRemove = c("select2d", "lasso2d"))
  }
}

#
#plotHght <- "400px"
plotHght <- paste0(round(screenRes$height*0.75), "px")
nmsHtMp <- names(plotLeatMaps$Global)
nmsHtMp <- intersect(union("None", nmsHtMp), nmsHtMp)
plotHtMpHght <- paste0(round(min(c(400L, vapply(nmsHtMp, \(nm) { plotLeatMaps$Global[[nm]]$Render$sizingPolicy$defaultHeight }, 1)))), "px")

# UI functions
tbl_css <- tags$style(HTML("table.dataTable th,
table.dataTable td {
  white-space: normal !important;
  vertical-align: top;
}
table.dataTable td .cell-wrap {
  display: block;
  white-space: normal !important;
  overflow-wrap: anywhere;
  word-break: break-word;
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
make_prot_tab <- \(dflt = dfltProt,
                   prots = allProt,
                   shiny = TRUE) {
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
  myTags <- if (shiny) {
    tags$div(
      if (prot.list.Cond && (length(prots) > 1L)) {
        tags$div(selectInput("myProtein", "Select protein", prots, dflt),
                 br())
      },
      uiOutput("protComment"),
      br(),
      fluidRow(column(6L,
                      if (length(Exp) > 1L) {
                        selectInput("mySample", "", Exp, Exp[1L]) 
                      },
                      plotlyOutput("coverPlot", height = plotHght)),
               column(6L,
                      plotlyOutput("ratioPlot", height = plotHght))),
      br(),
      br(),
      tags$hr(style = "border-color: black;"),
      uiOutput("protPep"),
      br())
  } else {
    ## Coverage plots ###################################################
    dfltSmpl <- Exp[1L]
    exp2smpl <- listMelt(lapply(prots, \(pr) { Exp }), prots, ColNames = c("Sample", "Protein"))
    cov_plots <- lapply(1L:nrow(exp2smpl), \(i) {
      exp <- exp2smpl$Sample[i]
      pr <- exp2smpl$Protein[i]
      tags$div(id = paste0("cov_", pr, "_", exp),
               style = paste("width: 100%; display: ",
                             if ((pr == dflt) && (exp == dfltSmpl)) { "block" } else { "none" },
                             ";"),
               covPlots[[pr]]$logInt[[exp]])
    })
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
    pepTables <- lapply(prots, \(pr) {
      m <- match(pr, prots)
      tags$div(id = paste0("pepTable_", m),
               style = if (pr == dflt) { "width: 100%; display: block;" } else { "display: none;" },
               make_tbl_ui(tab = "All peptidoforms",
                           filt = pr))
    })
    ## UI ###############################################################
    tagList(
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
                      if (length(Exp) > 1L) {
                        make_select_tag("mySample",
                                        "",
                                        "mySample",
                                        Exp,
                                        Exp[1L])
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
      pepTables,
      tags$script(HTML("function updateProteinTab() {
  const prot = document.getElementById('myProtein').value;
  const ind = document.getElementById('myProtein').selectedIndex + 1;
  const sample = document.getElementById('mySample').value;
  const comm = 'prComment_' + ind
  document.querySelectorAll('[id^=\"cov_\"]').forEach(function(el) {
      el.style.display = 'none';
  });
  const cov = document.getElementById('cov_' + prot + '_' + sample);
  if (cov)
    cov.style.display = 'block';
  document.querySelectorAll('[id^=\"pepTable_\"]').forEach(function(el) {
      el.style.display = 'none';
  });
  const pepTblID = document.getElementById('pepTable_' + ind);
  if (pepTblID)
    pepTblID.style.width = '100%';
    pepTblID.style.display = 'block';
  document.querySelectorAll('[id^=\"rat_\"]').forEach(function(el) {
    el.style.display = 'none';
  });
  const rat = document.getElementById('rat_' + prot);
  if (rat)
    rat.style.display = 'block';
  document.querySelectorAll('[id^=\"prComment_\"]').forEach(function(el) {
    el.style.display = 'none';
  });
  document.getElementById(comm).style.display = 'block';
  document.getElementById(comm).style.whiteSpace = 'pre-wrap';
  document.getElementById(comm).style.padding = '10px';
  window.dispatchEvent(
    new Event('resize')
  );
}
document.getElementById('mySample').addEventListener('change', updateProteinTab);
const myProt = document.getElementById('myProtein');
if (myProt)
  myProt.addEventListener('change', updateProteinTab);")),
      br(),
      tags$hr(style = "border-color: black;"))
  }
return(myTags)
}
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
make_tbl_ui <- \(exp = Exp, #exp <- Exp[1L] #exp <- Exp[2L]
                 tab = "Protein groups", # can also be "All peptidoforms"; we will eventually add "`PTM`-modified", where `PTM` can be any PTM of interest
                 filt = NULL, #filt = allProt[1L] # Filter by "Common Name"
                 dat = xlDat,
                 minN = 1L) {
  pgTest <- (tab == "Protein groups")
  df <- dat[[tab]]
  smplCols_lst <- setNames(lapply(exp, \(xp) {
    grep(topattern(paste0(" ", xp), FALSE, TRUE), colnames(df), value = TRUE)
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
  fullIntRoot <- rev(paste0(vapply(strsplit(xprCols, "\\)"), `[[`, "", 1L), ") "))[1L]
  xprCols <- paste0(fullIntRoot, exp)
  repXprCols <- if (length(exp) == 1L) { fullIntRoot } else { xprCols }
  repXprCols <- sub("\\([Nn]or\\. ", "(norm. ", repXprCols)
  #
  ratCols <- grep(paste0("log2\\(([^\\)]+ )?rat\\.\\) "), smplCols, value = TRUE)
  useRat <- length(ratCols) > 0L
  if (useRat) {
    fullRatRoot <- rev(paste0(vapply(strsplit(ratCols, "\\)"), `[[`, "", 1L), ") "))[1L]
    ratCols <- paste0(fullRatRoot, exp)
    ratCols <- intersect(ratCols, colnames(df))
    repRatCols <- ratCols
    repRatCols <- sub("\\([Nn]or\\. ", "(norm. ", repRatCols)
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
    pepCountCols <- intersect(paste0("Pep. count ", exp), colnames(df))
    psmCountCols <- intersect(paste0("PSMs count ", exp), colnames(df))
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
    flt <- flt[which(apply(df[flt, pepCountCols, drop = FALSE], 1L, max, na.rm = TRUE) >= minN)]
  }
  df <- df[flt,]
  colnames(df) <- colNms <- repColNms
  if ("Modified sequence" %in% colNms) {
    df$"Modified sequence" <- gsub("^_|_$", "",  df$"Modified sequence")
  }
  xprCols <- repXprCols
  if (useRat) { ratCols <- repRatCols }
  covCols <- NULL
  if (pgTest) {
    covCols <- paste0("Cov. ", exp)
    repCovCols <- if (length(exp) == 1L) { "Coverage" } else { covCols }
    if (covCols %in% colnames(dat[[tab]])) {
      df[, repCovCols] <- dat[[tab]][flt, covCols]
    } else {
      if (("Coverage" %in% names(dat)) && (covCols %in% colnames(dat$Coverage))) {
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
  if (length(exp) > 1L) { orderVect <- apply(orderVect, 1L, \(x) { mean(x[which(is.finite(x))]) }) }
  df <- df[order(orderVect, decreasing = TRUE),]
  #
  quantCols <- xprCols
  if (useRat) {
    ratRng <- range(df[[ratCols]], na.rm = TRUE)
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
  return(df)
}
make_summTbl_ui <- \() {
  df <- t(Exp_summary[, grep(" - % ", colnames(Exp_summary), invert = TRUE, value = TRUE)])
  colnames(df) <- df[1L,]
  df <- df[2L:nrow(df),]
  #
  # Drop fixed PTMs (alkylation)
  fxdMods <- Modifs$`Full name`[which(Modifs$Type == "Fixed")]
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
  lQ <- length(quant)
  if (shiny) {
    tagList(make_comment_ui(exp, shiny),
            selectInput(paste0("quant_", exp), "", quant, dflt[exp]),
            plotlyOutput(paste0("quantLy_", exp), height = plotHght),
            br(),
            br(),
            tags$hr(style = "border-color: black;"),
            make_tbl_ui(exp))
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
    tagList(make_comment_ui(exp, shiny),
            lapply(1L:lQ, \(i) {
              tags$div(id = paste0("quant_", exp, "_", as.character(i)),
                       style = if (quant[i] == dflt[exp]) { "display: block;" } else { "display: none;" },
                       ggQuantLy[[quant[i]]][[exp]]$plotly)
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
            make_tbl_ui(exp))
  }
}
make_strt_tab <- \(shiny = TRUE) {
  if (shiny) {
    tagList(make_comment_ui("Dataset overview", shiny),
            br(),
            h4(strong(tags$ul(em("Summary table")))),
            make_summTbl_ui(),
            br(),
            br(),
            if (tstHtMp) {
              column(12L,
                     selectInput("myHeatMap",
                                 "",
                                 nmsHtMp,
                                 nmsHtMp[1L]),
                     plotlyOutput("heatMap", height = plotHtMpHght))
            },
            fluidRow(
              if (tstPCA) {
                column(4L*(3L-tstVenn),
                       plotlyOutput("PCA", height = plotHtMpHght))
              },
              if (tstVenn) {
                column(4L,
                       plotlyOutput("Venn", height = plotHtMpHght))
              },
            ),
            br())
  } else {
    styleOn <- paste0("display: block; height: ", plotHtMpHght)
    tagList(make_comment_ui("Dataset overview", shiny),
            br(),
            h4(strong(tags$ul(em("Summary table")))),
            make_summTbl_ui(),
            br(),
            br(),
            if (tstHtMp) {
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
                                         plotLeatMaps$Global[[nm]]$Render)
                              })))
            },
            fluidRow(
              if (tstPCA) {
                column(4L*(3L-tstVenn),
                       tags$div(id = "PCA",
                                style = styleOn,
                                dimRedPlotLy$PCA))
              },
              if (tstVenn) {
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
});"))))
  }
}
make_QC_ui <- \(shiny = TRUE,
                plotsList = QC_plotLys) {
  if (shiny) {
    tabPanel("QC",
             tagList(selectInput("QC1", "", names(plotsList), names(plotsList)[1L]),
                     fluidRow(column(8L,
                                     plotlyOutput("QCplotLy", height = plotHght)),
                              column(4L,
                                     uiOutput("QCtxt"))),
                     br()))
  } else {
    QC_comments <- allComments[names(plotsList)]
    tabPanel("QC",
             tagList(make_select_tag("myQC",
                                     "",
                                     "myQC",
                                     names(plotsList),
                                     names(plotsList)[1L]),
                     lapply(seq_along(plotsList), \(i) {
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
});"))))
  }
}
make_matmet_ui <- \(matmeth = matmethTxt,
                    shiny = TRUE) {
  # We want to load the processed materials and methods (potentially edited by the user)
  # ============> This should be ideally run as part of the finalization script, after the materials and method edition stage
  #
  hght <- vapply(strsplit(matmeth, "\n"), \(x) {
    paste0(as.character(20L*(sum(ceiling(nchar(unlist(x))/ceiling(screenRes$width/5.75)))+2L)), "px")
  }, "")
  if (shiny) {
    tabPanel("Materials and methods",
             tagList(textAreaInput("MatMet_SamplePrep", matmethSections[1L], matmeth[1L], "100%", hght[1L]),
                     br(),
                     textAreaInput("MatMet_LCMS", matmethSections[2L], matmeth[2L], "100%", hght[2L]),
                     br(),
                     textAreaInput("MatMet_DataAnalysis", matmethSections[3L], matmeth[3L], "100%", hght[3L]),
                     br()))
  } else {
    tabPanel("Materials and methods",
             tagList(h5(matmethSections[1L]),
                     tags$p(matmeth[1L]),
                     br(),
                     h5(matmethSections[2L]),
                     tags$p(matmeth[2L]),
                     br(),
                     h5(matmethSections[3L]),
                     tags$p(matmeth[3L]),
                     br()))
  }
}
make_ui <- \(tabNames = myTabs,
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
                      tagList(
                        if (tstRat) {
                          make_prot_tab(dfltProt,
                                        shiny = shiny)
                        },
                        br()
                      )))
    }
    if (x == "QC") {
      return(make_QC_ui(shiny = shiny))
    }
    if (x == "Materials and methods") {
      return(make_matmet_ui(shiny = shiny))
    }
  })
  return(bslib::navset_tab(!!!tabs))
}

# Plot HTML paths
#myPlots <- list.files(paste0(wd, "/Ranked abundance/LFQ"), "\\.html$", full.names = TRUE)
#names(myPlots) <- gsub(".* - |\\.html$", "", myPlots)
#nPl <- length(myPlots)
myTabs <- nms <- c("Dataset overview", Exp)
if (prot.list.Cond) {
  allProt <- names(covPlots)
  dfltProt <- allProt[1L]
  myTabs <- union(myTabs, "Proteins of interest")
  nms <- union(nms, allProt)
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
                tags$head(tags$style(HTML("table.dataTable td {
    white-space: normal !important;
    vertical-align: top !important;
}")),
                          tags$style(HTML(".cell-wrap {
    max-height: 80px !important;
    overflow: hidden !important;
    white-space: normal;
    vertical-align: top !important;
    overflow-wrap: break-word !important;
    word-break: break-word !important;
}"))),
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
  QUANT <- reactiveVal(dfltQuant)
  XPRTMSG <- reactiveVal(NULL)
  MYPROT <- reactiveVal(dfltProt)
  SAMPLE <- reactiveVal(Exp[1L])
  NORMMETH <- reactiveVal("None")
  # Render UI
  output$xprtMsg <- renderUI(XPRTMSG())
  output$myUI <- renderUI({ make_ui() })
  if (tstHtMp) {
    output$heatMap <- renderPlotly(plotLeatMaps$Global[[NORMMETH()]]$Render)
  }
  if (tstPCA) {
    output$PCA <- renderPlotly(dimRedPlotLy$PCA)
  }
  if (tstVenn) {
    output$Venn <- renderPlotly(plotly_Venn$`Global, LFQ`)
  }
  #
  for (exp in Exp) {
    idQ <- paste0("quant_", exp)
    idQLy <- paste0("quantLy_", exp)
    output[[idQLy]] <- renderPlotly(ggQuantLy[[input[[idQ]]]][[exp]]$plotly)
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
  sapply(Exp, \(exp) {
    idQ <- paste0("quant_", exp)
    idQLy <- paste0("quantLy_", exp)
    observeEvent(input[[idQ]], {
      dfltQuant <- QUANT()
      dfltQuant[exp] <- input[[idQ]]
      QUANT(dfltQuant)
      assign("dfltQuant", dfltQuant, envir = .GlobalEnv)
      # Update plot
      output[[idQLy]] <- renderPlotly(ggQuantLy[[input[[idQ]]]][[exp]]$plotly)
    })
  })
  #  - Proteins tab
  if (prot.list.Cond) {
    output$ratioPlot <- renderPlotly(ratioPlots[[MYPROT()]])
    output$coverPlot <- renderPlotly(covPlots[[MYPROT()]]$logInt[[SAMPLE()]])
    output$protComment <- renderUI(make_comment_ui(MYPROT()))
    output$protPep <- renderUI(make_tbl_ui(tab = "All peptidoforms",
                                           filt = MYPROT()))
    if (length(allProt) > 1L) {
      observeEvent(input$myProtein, { MYPROT(input$myProtein) })
    }
    if (length(Exp) > 1L) {
      observeEvent(input$mySample, { SAMPLE(input$mySample) })
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
runKount <- 0L
if (exists("appRunTst")) { rm(appRunTst) }
while ((!runKount) || (!exists("appRunTst")) || (!file.exists(htmlRprtFl))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  shinyCleanup()
  runKount <- runKount + 1L
}

# We now have our html... but it depends on local libraries...
# ---> We want those embedded in it so it is fully portable!
h2 <- h1 <- readLines(htmlRprtFl)
rg1 <- grep("</?head>", h1) + c(1L, -1L)
rg1 <- rg1[1L]:rg1[2L]
hd1 <- h1[rg1]
hd1 <- data.frame(original = hd1)
hd1$new <- hd1$original
g <- grep("^ *<((style)|(script)|(link))( *[^>]+)?>", hd1$original)
hd1$original[g]
require(base64enc)
read_file <- \(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
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

# To do: write Mat Meth locally
