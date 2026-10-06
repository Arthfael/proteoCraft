# UI functions used by the HTML report
safe_id <- \(x) { gsub("[: ]", "_", x) }
myCSS <- tags$style(HTML("table.dataTable td {
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
}
.plot-container {
     border-radius: 10px;
     overflow: hidden;
}
"))

fix_plotly_sizing <- \(x, px_height = 700L, depth = 0L) {
  if (depth > 40L) { return(x) }
  if (inherits(x, "plotly")) {
    x$x$layout$xaxis$autorange <- TRUE
    x$x$layout$yaxis$autorange <- TRUE
    x <- htmlwidgets::onRender(x, global_autorange)
    x$sizingPolicy$browser$fill <- FALSE   # stop depending on a flex ancestor for height
    x$height <- px_height                  # fixed pixel height, independent of CSS context
    return(plotly::config(x, modeBarButtonsToRemove = c("select2d", "lasso2d")))
  }
  if (is.list(x) && identical(class(x), "list")) { return(lapply(x, fix_plotly_sizing, px_height = px_height, depth = depth + 1L)) }
  x
}

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
                     ON = TRUE) {
  if (!id %in% names(values)) { stop("Invalid comment name!") }
  if (shiny) {
    textAreaInput(inputId = paste0("comment_", safe_id(id)),
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
                       minN = 1L,
                       regMsg,
                       filtMsg) {
  stopifnot(tab %in% names(dat))
  m <- match(contr, myContrasts$Contrast)
  contr2 <- sub(" - ", " /\n", contr)
  exp <- setNames(myContrasts[m, c("A_samples", "B_samples")],
                  myContrasts[m, c("A", "B")])
  exp <- lapply(exp, \(x) { Exp.map$Clean_name[match(unlist(x), Exp.map[[RSA$column]])] })
  grps <- names(exp)
  pgTest <- (tab == "Protein groups")
  df <- dat[[tab]]
  tmp <- grep("\n", colnames(df), value = TRUE)
  tmp2 <- gsub(".*\n", "", sub(" /\n.*", "", tmp))
  smplCols_lst <- setNames(lapply(grps, \(xp) {
    tmp[tmp2 %in% #c(
          xp#, exp[[xp]])
    ]
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
    filtCol <- "Proteins"
    coreCols <- union(c("Modified sequence_verbose", #"Sequence",
                        filtCol), coreCols) # Check before use...
    intRoot <- "int" # Presumably...
  }
  #
  xprCols <- grep(paste0("log10\\(([^\\)]+ )?", intRoot, "\\.\\) "), colnames(df), value = TRUE)
  fullIntRoot <- rev(paste0(vapply(strsplit(xprCols, "\n"), `[[`, "", 1L), "\n"))[1L]
  smpls <- unlist(exp[grps])
  smpls_and_grps <- unlist(lapply(grps, \(grp) { c(grp, exp[[grp]]) }))
  xprCols <- intersect(paste0(fullIntRoot, smpls),
                       colnames(df))
  if (!length(xprCols)) { # In this case we only have an average column
    xprCols <- intersect(paste0(smpls_and_grps, smpls),
                         colnames(df))
  }
  repXprCols <- if (length(grps) == 1L) { sub(" *\n$", "", fullIntRoot) } else { xprCols }
  #
  ratCols <- grep("log2\\(.*rat\\.\\) \n", smplCols, value = TRUE)
  stopifnot(length(ratCols) > 0L) # For contrasts we always MUST have a logFC column!
  fullRatRoot <- rev(paste0(vapply(strsplit(ratCols, "\n"), `[[`, "", 1L), "\n"))[1L]
  ratCol <- paste0(fullRatRoot, contr2)
  ratCol <- repRatCol <- intersect(ratCol, colnames(df)) 
  stopifnot(length(ratCol) == 1L) # Again
  #
  PValCol <- paste0(sub(" -log10\\(Pvalue\\) - ", "\n-log10 pval. \n", pvalue.col[pvalue.use]), contr2)
  repPValCol <- sub("\n-log10 pval\\. \n", "\npval. \n", PValCol)
  decCol <- repDecCol <- paste0("reg. \n", contr2)
  #
  colNms <- c(coreCols, xprCols, ratCol, PValCol, decCol)
  repColNms <- c(sub("_verbose$", "",
                     sub("^Potential contaminant$", "Cont.",
                         sub("^Mol\\. weight \\[kDa\\]$", "MW (kDa)",
                             sub("^Common Names$", "Common names", coreCols)))),
                 repXprCols, repRatCol, repPValCol, repDecCol)
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
  df[[PValCol]] <- 10L^-df[[PValCol]]
  flt <- if (is.null(filt)) {
    1L:nrow(df)
  } else {
    grsep(db$`Protein ID`[match(filt, db$`Common Name`)], x = df[[filtCol]])
  }
  if (pgTest && is.integer(minN) && (minN > 0L) && length(pepCountCols)) {
    flt <- flt[apply(df[flt, pepCountCols, drop = FALSE], 1L, max, na.rm = TRUE) >= minN]
  }
  if (!length(flt)) {
    if (missing(filtMsg)) { filtMsg <- paste0("No matching ", sub("^All", "", tab), " to show...") }
    return(div(em(HTML(filtMsg)),
               br(),
               br()))
  }
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
  # Filter now to only return regulated proteins
  rg <- grep("^((up)|(down)),|specific", df[[repDecCol]])
  if (!length(rg)) {
    if (missing(regMsg)) { regMsg <- paste0("No significant ", tolower(gsub("^All|s$", "", tab)), " to show...") }
    return(div(em(HTML(regMsg)),
               br(),
               br()))
  }
  df <- df[rg, ]
  df <- df[, setdiff(colnames(df), repDecCol)]
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
                                     scrollY = "1000px",
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
                   shiny = TRUE,
                   tables) {
  if (missing(tables)) { tables <- !shiny }
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
                      withSpinner(plotlyOutput("coverPlot", height = plotHght))),
               if (scrptType == "noReps") {
                 column(6L,
                        withSpinner(plotlyOutput("ratioPlot", height = plotHght)))
               },
      ),
      br(),
      br(),
      tags$hr(style = "border-color: black;"),
      if (tables && peptidoTst) { uiOutput("protPep") },
      br(),
      style = paste0("background: ", myCol, ";")))
  } else {
    ## Coverage plots ###################################################
    if (cov_ON) {
      dfltExp <- myExp[1L]
      exp2smpl <- listMelt(lapply(prots, \(pr) { myExp }), prots, ColNames = c("Sample", "Protein"))
      cov_plots <- lapply(1L:nrow(exp2smpl), \(i) {
        exp <- exp2smpl$Sample[i]
        pr <- exp2smpl$Protein[i]
        tags$div(id = paste0("cov_", pr, "_", exp),
                 style = paste("width: 100%; display: ",
                               if ((pr == dflt) && (exp == dfltExp)) { "block" } else { "none" },
                               ";"),
                 class = "plot-container",
                 fix_plotly_sizing(covPlots[[pr]]$logInt[[exp]]))
      })
    }
    ## Ratio plots ######################################################
    ratio_plots_ui <- NULL
    if (rat_ON) {
      prots2 <- intersect(prots, names(ratioPlots))
      if (length(prots2)) {
        ratio_plots_ui <- lapply(prots2, \(pr) {
          tags$div(id = paste0("rat_", pr),
                   style = paste("width: 100%; display: ",
                                 if (pr == dflt) { "block" } else { "display" },
                                 ";"),
                   class = "plot-container",
                   fix_plotly_sizing(ratioPlots[[pr]]))
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
                      pr == dflt)
    })
    #
    ## Peptide tables  ##################################################
    if (tables && peptidoTst) {
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
      if (tables && peptidoTst) { pepTables },
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
                      width = "100%",
                      #width = wdth,
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
  return(tags$div(h4(strong(em("Summary table"))),
                  df,
                  style = "background: #ffffff; overflow: auto;"))
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
                   quant,
                   dflt = dfltQuant,
                   tables) {
  if (runRankAbundPlots) {
    if (missing(quant)) { quant <- quantMeth }
    if (missing(dflt)) { dflt <- dfltQuant }
  }
  if (missing(tables)) { tables <- !shiny }
  myCol <- tolower(viridis::viridis(6L, alpha = 0.2)[2L])
  exp2 <- if (scrptType == "noReps") { exp } else { names(smplGrps)[match(exp, smplGrps)] }
  exp_ <- safe_id(exp)
  lQ <- length(quant)
  id1 <- paste0("quant_", exp)
  if (shiny) {
    tagList(tags$div(
      uiOutput(paste0("cmmnt_", exp_)),
      if (runRankAbundPlots) {
        div(selectInput(id1, "", quant, dflt[exp]),
            withSpinner(plotlyOutput(paste0("quantLy_", exp), height = plotHght)),
            br())
      },
      br(),
      tags$hr(style = "border-color: black;"),
      if (tables) {
        uiOutput(paste0("PG_table_", exp))
      },
      style = paste0("background: ", myCol, ";")))
  } else {
    id2 <- paste0(id1, "_")
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
      if (runRankAbundPlots) {
        div(lapply(1L:lQ, \(i) {
          tags$div(id = paste0("quant_", exp, "_", as.character(i)),
                   style = if (quant[i] == dflt[exp]) { "display: block;" } else { "display: none;" },
                   class = "plot-container",
                   fix_plotly_sizing(ggQuantLy[[quant[i]]][[exp2]]$plotly))
        }),
        br(),
        make_select_tag(id1,
                        "",
                        id1,
                        quant,
                        dflt[exp]),
        br())
      },
      br(),
      tags$hr(style = "border-color: black;"),
      tags$script(HTML(js)),
      if (tables) {
        make_smpl_tbl_ui(exp)
      },
      style = paste0("background: ", myCol, ";")))
  }
}
make_ctrst_tab <- \(contr,
                    shiny = TRUE,
                    tables,
                    ptm) {
  if (missing(tables)) { tables <- !shiny }
  myCol <- tolower(viridis::viridis(6L, alpha = 0.2)[3L])
  styleOn <- paste0("display: block; height: ", plotHtMpHght)
  contr_ <- safe_id(contr)
  volcID <- paste0(contr_, "_volcPlot")
  commentID <- paste0("cmmnt_", contr_)
  goID <- paste0(contr_, "_GObars")
  saintIDs <- c(paste0("SAINTexpress volcano plot ", contr),
                paste0(contr_, c("_SAINT_volcPlot", "_SAINT_GObars")))
  saintXPRS <- saintExprs & (saintIDs[1L] %in% names(volcPlotly$SAINTexpress))
  GSEA_tst <- runGSEA
  if (runGSEA) {
    GSEA_IDs <- paste0(contr_, "_GSEA", as.character(1L:4L))
    GSEA_tst <- { if (missing(ptm)) { "PG" } else { ptm } } %in% wh_GSEA[[contr]]
  }
  tblID <- paste0(contr_, "_tbl")
  if (!missing(ptm)) {
    volcID <- paste0(ptm, "_", volcID)
    commentID <- paste0("cmmnt_", ptm, "_", contr_)
    goID <- paste0(ptm, "_", goID)
    if (runGSEA) {
      GSEA_IDs <- paste0(ptm, "_", GSEA_IDs)
    }
    tblID <- paste0(ptm, "_", tblID)
    # SAINTexpress and PTMs are not compatible right now - do they even make sense to combine?
  }
  volcNm1 <- if (missing(ptm)) { "t-test" } else { paste0(ptm, " t-test") }
  volcNm2 <- if (missing(ptm)) { paste0("Volcano plot ", contr) } else { paste0(ptm, " volcano plot ", contr) }
  slotNm <- if (missing(ptm)) { "PG" } else { ptm }
  cmmntNm <- if (missing(ptm)) { contr} else { paste0(Ptm, ": ", contr) }
  if (shiny) {
    tagList(tags$div(
      uiOutput(commentID),
      if (missing(ptm) && saintXPRS) {
        div(h3("SAINTexpress"),
            fluidRow(column(6L,
                            withSpinner(plotlyOutput(saintIDs[2L], height = "600px"))),
                     if (enrichGO) {
                       column(6L,
                              br(),
                              br(),
                              withSpinner(plotlyOutput(saintIDs[3L], height = "600px")))
                     },
            ),
            style = "background: #ffffff;")
      } else {
        div(h3("t-test"),
            fluidRow(column(6L,
                            withSpinner(plotlyOutput(volcID, height = "600px"))),
                     if (enrichGO) {
                       column(6L,
                              br(),
                              br(),
                              withSpinner(plotlyOutput(goID, height = "600px")))
                     },
            ),
            style = "background: #ffffff;")
      },
      tags$hr(style = "border-color: black;"),
      if (F.test) {
        # Add F-test part here... or maybe dropdown to choose f-/F-test... or drop F-test altogether?
      },
      if (GSEA_tst) {
        div(
          div(h3("GSEA"),
              fluidRow(column(6L,
                              h5("dot plot"),
                              withSpinner(plotlyOutput(GSEA_IDs[1L], height = "600px")),
                              h5("enrichment map"),
                              withSpinner(plotlyOutput(GSEA_IDs[2L], height = "600px"))),
                       column(6L,
                              h5("ridge plot"),
                              withSpinner(plotlyOutput(GSEA_IDs[3L], height = "600px")),
                              h5("net plot"),
                              withSpinner(plotlyOutput(GSEA_IDs[4L], height = "600px")))),
              style = "background: #ffffff;"),
          tags$hr(style = "border-color: black;"))
      },
      if (tables) {
        uiOutput(tblID)
      },
      style = paste0("background: ", myCol, ";")))
  } else {
    styleOn6 <- "display: block; height: 600px"
    styleOn5 <- "display: block; height: 500px"
    tagList(tags$div(
      make_comment_ui(cmmntNm, shiny),
      if (missing(ptm) && saintXPRS) {
        div(h3("SAINTexpress"),
            fluidRow(column(6L,
                            tags$div(id = saintIDs[2L],
                                     style = styleOn6,
                                     class = "plot-container",
                                     fix_plotly_sizing(volcPlotly$SAINTexpress[[saintIDs[1L]]]$Plot, 600L))),
                     if (enrichGO && (!is.null(GO_plot_ly$Prot$SAINTexpress[[contr]]$Bar))) {
                       column(6L,
                              br(),
                              br(),
                              tags$div(id = saintIDs[3L],
                                       style = styleOn6,
                                       class = "plot-container",
                                       fix_plotly_sizing(GO_plot_ly$Prot$SAINTexpress[[contr]]$Bar, 600L)))
                     },
            ),
            style = "background: #ffffff;")
      } else {
        div(h3("t-test"),
            fluidRow(column(6L,
                            tags$div(id = volcID,
                                     style = styleOn6,
                                     class = "plot-container",
                                     fix_plotly_sizing(volcPlotly[[volcNm1]][[volcNm2]]$Plot, 600L))),
                     if (enrichGO && (!is.null(GO_plot_ly[[slotNm]]$"t-test"[[contr]]$Bar))) {
                       column(6L,
                              br(),
                              br(),
                              tags$div(id = goID,
                                       style = styleOn6,
                                       class = "plot-container",
                                       fix_plotly_sizing(GO_plot_ly[[slotNm]]$"t-test"[[contr]]$Bar, 600L)))
                     },
            ),
            style = "background: #ffffff;")
      },
      tags$hr(style = "border-color: black;"),
      if (F.test) {
        # Add F-test part here... or maybe dropdown to choose f-/F-test... or drop F-test altogether?
      },
      if (GSEA_tst) {
        div(
          div(h3("GSEA"),
              # NB: I also tried the plotly::subplot() approach to displaying the plots together in one,
              # but this fails (subplots look corrupted, possibly because they are slightly hacky)
              fluidRow(column(6L,
                              h5("dot plot"),
                              tags$div(id = GSEA_IDs[1L],
                                       style = styleOn5,
                                       class = "plot-container",
                                       fix_plotly_sizing(GSEA_plotly$standard[[slotNm]][[GSEA_plotNms[1L]]][[contr]], 500L)),
                              h5("enrichment map"),
                              tags$div(id = GSEA_IDs[2L],
                                       style = styleOn5,
                                       class = "plot-container",
                                       fix_plotly_sizing(GSEA_plotly$standard[[slotNm]][[GSEA_plotNms[2L]]][[contr]], 500L))),
                       column(6L,
                              h5("ridge plot"),
                              tags$div(id = GSEA_IDs[3L],
                                       style = styleOn5,
                                       class = "plot-container",
                                       fix_plotly_sizing(GSEA_plotly$standard[[slotNm]][[GSEA_plotNms[3L]]][[contr]], 500L)),
                              h5("net plot"),
                              tags$div(id = GSEA_IDs[4L],
                                       style = styleOn5,
                                       class = "plot-container",
                                       fix_plotly_sizing(GSEA_plotly$standard[[slotNm]][[GSEA_plotNms[4L]]][[contr]], 500L)))),
              style = "background: #ffffff;"),
          tags$hr(style = "border-color: black;"))
      },
      if (tables) {
        if (missing(ptm)) {
          make_ctrst_tbl_ui(contr)
        } else {
          make_ctrst_tbl_ui(contr,
                            paste0(ptm, "-mod. pept."))
        }
      },
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
      make_summTbl_ui(),
      br(),
      if (heatMaps_ON) {
        div(fluidRow(column(1L,
                            em("Heatmap type:")),
                     column(11L,
                            selectInput("myHeatMap",
                                        "",
                                        nmsHtMp,
                                        nmsHtMp[1L]))),
            fluidRow(column(12L,
                            withSpinner(plotlyOutput("heatMap", height = "700px")))))
      },
      br(),
      fluidRow(
        if (globalGO && (!is.null(GO_plot_ly$PG$Dataset$`Observed dataset`$Bar))) {
          column(5L,
                 withSpinner(plotlyOutput("GO_enrich_Dataset", height = "700px")))
        },
        if (PCA_ON) {
          column(4L,
                 withSpinner(plotlyOutput("PCA", height = "700px")))
        },
        if (Venn_ON) {
          column(3L,
                 withSpinner(plotlyOutput("Venn", height = "700px")))
        },
      ),
      br(),
      style = paste0("background: ", myCol, ";")))
  } else {
    styleOn7 <- "display: block; height: 700px"
    tagList(tags$div(
      make_comment_ui("Dataset overview", shiny),
      br(),
      make_summTbl_ui(),
      br(),
      if (heatMaps_ON) {
        div(fluidRow(column(1L,
                            em("Heatmap type:")),
                     column(11L,
                            make_select_tag("myHeatMap",
                                            "",
                                            "myHeatMap",
                                            nmsHtMp,
                                            nmsHtMp[1L]))),
            fluidRow(column(12L,
                            lapply(nmsHtMp, \(nm) {
                              i <- match(nm, nmsHtMp)
                              tags$div(id = paste0("HeatMap_", i),
                                       style = if (i == 1L) { styleOn7 } else { "display: none;" },
                                       class = "plot-container",
                                       fix_plotly_sizing(plotLeatMaps$Global[[nm]]$Plot))
                            }))))
      },
      br(),
      fluidRow(
        if (globalGO && (!is.null(GO_plot_ly$PG$Dataset$`Observed dataset`$Bar))) {
          column(5L,
                 tags$div(id = "GO_enrich_Dataset",
                          style = styleOn7,
                          class = "plot-container",
                          fix_plotly_sizing(GO_plot_ly$PG$Dataset$`Observed dataset`$Bar)))
        },
        if (PCA_ON) {
          column(4L,
                 tags$div(id = "PCA",
                          style = styleOn7,
                          class = "plot-container",
                          fix_plotly_sizing(dimRedPlotLy$PG$PCA)))
        },
        if (Venn_ON) {
          column(3L,
                 tags$div(id = "Venn",
                          style = styleOn7,
                          class = "plot-container",
                          fix_plotly_sizing(plotly_Venn$`Global, LFQ`)))
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
  QC_comments <- allComments[names(plotsList)]
  if (shiny) {
    tagList(tags$div(
      selectInput("QC", "", names(plotsList), names(plotsList)[1L]),
      fluidRow(column(8L,
                      withSpinner(plotlyOutput("QCplotLy", height = plotHght))),
               column(4L,
                      uiOutput("QCtxt"))),
      br(),
      style = paste0("background: ", myCol, ";")))
  } else {
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
                                 class = "plot-container",
                                 fix_plotly_sizing(plotsList[[nm]], plotHght))),
                 column(4L,
                        make_comment_ui(nm,
                                        shiny,
                                        QC_comments,
                                        nm == names(plotsList)[1L])))
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
  hght <- vapply(lapply(matmeth, strsplit, split = "\n"), \(x) {
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
      em(HTML("This is an automatically-generated materials and methods template and may contain inaccuracies. Please remember to check with us the details before including this in a publication!")),
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
  #tabNames <- tabNames[c(1L, 8L, 13L)]
  tabs <- lapply(tabNames, \(x) {
    if (x == "Dataset overview") {
      return(tabPanel(x,
                      make_strt_tab(shiny = shiny)))
    }
    if (showSmplGrpTabs && (x %in% smplGrps)) {
      return(tabPanel(paste0("sample group = ", x),
                      make_smpl_tab(x,
                                    shiny = shiny)))
    }
    if (x %in% myContrasts$Contrast) {
      return(tabPanel(paste0("contrast = ", x),
                      make_ctrst_tab(contr = x,
                                     shiny = shiny)))
    }
    if (length(PTMs)) {
      tmp <- unlist(lapply(PTMs, \(Ptm) { paste0(Ptm, ": ", myContrasts$Contrast)}))
      if (x %in% tmp) { #x <- tmp[1L]
        Ptm <- sub(": .*", "", x)
        contr <- sub(topattern(paste0(Ptm, ": ")), "", x)
        return(tabPanel(paste0(Ptm, " contrast = ", contr),
                        make_ctrst_tab(contr = contr,
                                       shiny = shiny,
                                       ptm = Ptm)))
      }
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
