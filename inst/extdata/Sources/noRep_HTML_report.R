library(shiny)
library(shinyjs)
library(bslib)
library(htmltools)
library(DT)
library(plotly)

htmlRprtFl <- paste0(wd, "/Report_", dtstNm, ".html")
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
report_header <- tags$header(paste0(dtstNm, " - report"),
    br(),
    "Analysis run by: ", em(WhoAmI),
    br(),
    "Date: ", em(Sys.Date()),
    br(),
    "Package: ", em(paste0("proteoCraft v", package.version("proteoCraft"))),
    br(),
    br())

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
#plotHght <- "400px"
plotHght <- paste0(round(screenRes$height*0.75), "px")
nmsHtMp <- names(plotLeatMaps$Global)
nmsHtMp <- intersect(union("None", nmsHtMp), nmsHtMp)
plotHtMpHght <- paste0(round(min(c(400, vapply(nmsHtMp, \(nm) { plotLeatMaps$Global[[nm]]$Render$sizingPolicy$defaultHeight }, 1)))), "px")

# UI functions
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
    if (exists("ratioPlots")) {
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
      pepTables,
      tags$script(HTML("function updateProteinTab() {
  var prot = document.getElementById('myProtein').value;
  var ind = document.getElementById('myProtein').selectedIndex + 1;
  var sample = document.getElementById('mySample').value;
  var comm = 'prComment_' + ind
  document.querySelectorAll('[id^=\"cov_\"]').forEach(function(el) {
      el.style.display = 'none';
  });
  var cov = document.getElementById('cov_' + prot + '_' + sample);
  if (cov)
    cov.style.display = 'block';
  document.querySelectorAll('[id^=\"pepTable_\"]').forEach(function(el) {
      el.style.display = 'none';
  });
  var pepTblID = document.getElementById('pepTable_' + ind);
  if (pepTblID)
    pepTblID.style.width = '100%';
    pepTblID.style.display = 'block';
  document.querySelectorAll('[id^=\"rat_\"]').forEach(function(el) {
    el.style.display = 'none';
  });
  var rat = document.getElementById('rat_' + prot);
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
var prot = document.getElementById('myProtein');
if (prot)
  prot.addEventListener('change', updateProteinTab);")),
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
  df <- dat[[tab]]
  smplCols_lst <- setNames(lapply(exp, \(xp) {
    grep(topattern(paste0(" ", xp), FALSE, TRUE), colnames(df), value = TRUE)
  }), exp)
  smplCols <- setNames(unlist(smplCols_lst), NULL)
  coreCols <- c("PEP", "Potential contaminant")
  if (tab %in% c("Protein groups", "All peptidoforms")) {
    if (tab == "Protein groups") {
      filtCol <- "Protein IDs"
      coreCols <- union(c("Leading protein IDs", filtCol, "Common Names", "Genes", "Mol. weight [kDa]"), coreCols)
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
  repColNms <- c(sub("_verbose$", "", coreCols), repXprCols)
  if (useRat) {
    colNms <- c(colNms, ratCols)
    repColNms <- c(repColNms, repRatCols)
  }
  if (tab == "Protein groups") {
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
  if ((tab == "Protein groups") && is.integer(minN) && (minN > 0L) && length(pepCountCols)) {
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
  if (tab == "Protein groups") {
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
  covSortCols <- character(0)
  covSortVals <- NULL
  if (tab == "Protein groups" && (!is.null(covCols)) && length(covCols)) {
    covSortCols <- paste0(covCols, "__sort")
    covSortVals <- setNames(lapply(covCols, \(k) {
      suppressWarnings(as.numeric(df[[k]]))
    }), covSortCols)
  }
  col2 <- setdiff(colnames(df), c("PEP", quantCols))
  df[, col2] <- sapply(col2, \(k) {
    if ((tab == "Protein groups") && (k %in% covCols)) {
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
  df <- DT::datatable(df,
                      rownames = FALSE,
                      class = "compact",
                      escape = FALSE,
                      options = list(scrollX = TRUE,
                                     scrollY = "500px",
                                     autoWidth = FALSE,
                                     columnDefs = columnDefs_all))
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
  wdth <- paste0(250L*(ncol(df) + 1L), "px")
  #hght <- paste0(100L*(nrow(df)+1L), "px")
  df <- DT::datatable(df,
                      rownames = TRUE,
                      class = "compact",
                      escape = FALSE,
                      #autoHideNavigation = TRUE,
                      width = wdth,
                      #height = hght,
                      options = list(scrollX = TRUE,
                                     scrollY = "500px",
                                     autoWidth = FALSE,
                                     columnDefs = list(list(width = "200px",
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
  var selected = document.getElementById('%s').selectedIndex + 1;
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
            fluidRow(column(strtColWdth,
                            selectInput("myHeatMap",
                                        "",
                                        nmsHtMp,
                                        nmsHtMp[1L]),
                            plotlyOutput("heatMap", height = plotHtMpHght)),
                     column(strtColWdth,
                            plotlyOutput("PCA", height = plotHtMpHght))),
            br())
  } else {
    styleOn <- paste0("display: block; height: ", plotHtMpHght)
    tagList(make_comment_ui("Dataset overview", shiny),
            br(),
            h4(strong(tags$ul(em("Summary table")))),
            make_summTbl_ui(),
            br(),
            br(),
            fluidRow(
              if (tstHtMp) {
                column(strtColWdth,
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
                       }))
              },
              if (tstPCA) {
                column(strtColWdth,
                       tags$div(id = "PCA",
                                style = styleOn,
                                dimRedPlotLy$PCA))
              },
            ),
            br(),
            tags$script(HTML(paste0("document.getElementById('myHeatMap').addEventListener('change', function() {
  var HtMpID = document.getElementById('myHeatMap').selectedIndex + 1;
  var HtMp = document.getElementById('HeatMap_' + HtMpID);
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
  var selected = document.getElementById('myQC').selectedIndex + 1;
  document.querySelectorAll('[id^=\"QC_\"]').forEach(function(div) { div.style.display = 'none'; });
  document.getElementById('QC_' + selected).style.display = 'block';
  document.querySelectorAll('[id^=\"QCcomment_\"]').forEach(function(div) { div.style.display = 'none'; });
  var div = document.getElementById('QCcomment_' + selected);
  div.style.display = 'block';
  div.style.whiteSpace = 'pre-wrap';
  div.style.padding = '10px';
  window.dispatchEvent(new Event('resize'));
});"))))
  }
}
make_ui <- \(shiny = TRUE) {
  tabs <- lapply(myTabs, function(x) {
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
                        if (ratiosTest) {
                          make_prot_tab(dfltProt,
                                        shiny = shiny)
                        },
                        br()
                      )))
    }
    if (x == "QC") {
      return(make_QC_ui(shiny = shiny))
    }
  })
  return(bslib::navset_tab(!!!tabs))
}

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
strtColWdth <- 12L/(tstHtMp + tstPCA)

# Fix to plotly autoscaling
if (tstHtMp) {
  for (x in names(plotLeatMaps)) { #x <- names(plotLeatMaps)[1L]
    for (y in names(plotLeatMaps[[x]])) { #y <- names(plotLeatMaps[[x]])[1L]
      plotLeatMaps[[x]][[y]]$Render$x$layout$xaxis$autorange <- TRUE
      plotLeatMaps[[x]][[y]]$Render$x$layout$yaxis$autorange <- TRUE
    }
  }
}
if (exists("ggQuantLy")) {
  for (x in names(ggQuantLy)) { #x <- names(ggQuantLy)[1L]
    for (y in names(ggQuantLy[[x]])) { #y <- names(ggQuantLy[[x]])[1L]
      ggQuantLy[[x]][[y]]$plotly$x$layout$xaxis$autorange <- TRUE
      ggQuantLy[[x]][[y]]$plotly$x$layout$yaxis$autorange <- TRUE
    }
  }
}
if (tstPCA) {
  for (x in names(dimRedPlotLy)) { #x <- names(dimRedPlotLy)[1L]
    dimRedPlotLy[[x]]$x$layout$xaxis$autorange <- TRUE
    dimRedPlotLy[[x]]$x$layout$yaxis$autorange <- TRUE
  }
}
if (exists("covPlots")) {
  for (x in names(covPlots)) { #x <- names(covPlots)[1L]
    for (y in names(covPlots[[x]])) { #y <- names(covPlots[[x]])[1L]
      for (z in names(covPlots[[x]][[y]])) { #z <- names(covPlots[[x]][[y]])[1L]
        covPlots[[x]][[y]][[z]]$x$layout$xaxis$autorange <- TRUE
        covPlots[[x]][[y]][[z]]$x$layout$yaxis$autorange <- TRUE
      }
    }
  }
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
myTabs <- union(myTabs, "QC")
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
if (sum(!nms %in% names(allComments))) {
  nms_ <- setdiff(nms, names(allComments))
  allComments[nms_] <- ""
  if ("Dataset overview" %in% nms_) { allComments$"Dataset overview" <- dfltComment }
}
allComments %<o% allComments[nms]
ratiosTest <- exists("ratioPlots") && (length(ratioPlots) > 0L)

#
quantTst <- setNames(lapply(names(ggQuantLy), \(tp) { names(ggQuantLy[[tp]]) }), names(ggQuantLy))
quantTst <- listMelt(quantTst, ColNames = c("Sample", "Type"))
quantTst <- aggregate(quantTst$Type, list(quantTst$Sample), list)
colnames(quantTst) <- c("Sample", "Types")
dfltQuant <- quantTst[, "Sample", drop = FALSE]
dfltQuant$Type <- vapply(quantTst$Types, \(x) { x[[1L]] }, "")
quantTst <- setNames(quantTst$Types, quantTst$Sample)
dfltQuant <- setNames(dfltQuant$Type, dfltQuant$Sample)
quantMeth <- unique(unlist(quantTst))
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
  output$heatMap <- renderPlotly(plotLeatMaps$Global[[NORMMETH()]]$Render)
  if (tstPCA) {
    output$PCA <- renderPlotly(dimRedPlotLy$PCA)
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
        make_ui(FALSE))
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

# We now have our html... but it depends on local libraries... we want those embedded in it!
h2 <- h1 <- readLines(htmlRprtFl)
rg1 <- grep("</?head>", h1) + c(1L, -1L)
rg1 <- rg1[1L]:rg1[2L]
hd1 <- h1[rg1]
hd1 <- data.frame(original = hd1)
hd1$new <- hd1$original
g <- grep("^ *<((style)|(script)|(link))( *[^>]+)?>", hd1$original)
hd1$original[g]
require(base64enc)
read_file <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}
file_to_data_uri <- function(path) {
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


# TO DO
# - remove local /lib folder
# - overwrite report instead of keeping two versions
