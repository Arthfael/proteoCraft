# Create a html report to:
# - organize the most important tables and plots into a coherent set of html files
# - capture user comments on each section
#
library(shiny)
library(shinyjs)
library(bslib)
library(htmltools)
library(DT)
library(plotly)
library(jsonlite)
library(svDialogs)
require(base64enc)

showSmplGrpTabs <- TRUE
if (scrptType == "withReps") {
  smplGrps <- setNames(cleanNms(VPAL$values), VPAL$values)
  if (!"clean_Group_name" %in% colnames(Exp.map)) {
    Exp.map$clean_Group_name <- smplGrps[Exp.map[[VPAL$column]]]
  }
  # For now, as long as we stick to individual HTMLs per tabs as oppose to a single ginormous HTML,
  # also saving sample group tabs in all circumstances is fine.
  # if ((!exists("showSmplGrpTabs")) || (!is.logical(showSmplGrpTabs))) {
  #   showSmplGrpTabs <- ((WorkFlow %in% c("PULLDOWN", "BIOID", "LOCALISATION"))
  #                       | sum(c("TISSUE", "CELLLINE", "COMPARTMENT") %in% toupper(gsub(" _-\\.", "", Factors))))
  # }
  # showSmplGrpTabs %<o% showSmplGrpTabs[1L]
  # if (!exists("showSmplGrpTabs_askedOnce")) { showSmplGrpTabs_askedOnce %<o% FALSE }
  # if (!showSmplGrpTabs_askedOnce) {
  #   m <- match(showSmplGrpTabs, c(TRUE, FALSE))
  #   opt <- c("Yes                                                                                                ",
  #            "No                                                                                                 ")
  #   showSmplGrpTabs <- c(TRUE, FALSE)[match(sub(" +$", "", dlg_list(opt, opt[m],
  #                                                                   title = "Include sample group tabs in the HTML report?")$res),
  #                                           c("Yes", "No"))]
  #   showSmplGrpTabs_askedOnce <- TRUE
  # }
}

funFile <- paste0(libPath, "/extdata/Sources/HTML_report_Fun.R")
source(funFile)
# Below: leftover from an aborted attempt at parallelisation
# obj2Xport <- readLines(funFile)
# obj2Xport <- grep("^[^ ]+ *((<-)|(%<[co]%)) *", obj2Xport, value = TRUE)
# obj2Xport <- gsub(" *((<-)|(%<[co]%)) *.*", "", obj2Xport)
# obj2Xport <- intersect(obj2Xport, ls())

if (!exists("matmethTxt")) {
  matmethTxt <- c("Samples preparation" = paste(MatMetCalls$Texts$WetLab, collapse = "\n"),
                  "LC-MS/MS analysis" = paste(MatMetCalls$Texts$LCMS, collapse = "\n"),
                  "Data analysis" = paste(MatMetCalls$Texts$DatAnalysis, collapse = "\n"))
}
matmethSections <- names(matmethTxt)

if (!exists("xlDat")) { loadFun(paste0(wd, "/Tables/xlDat.RDS")) }
peptidoTst <- "All peptidoforms" %in% names(xlDat)
#
if (!exists("volcPlotly_fl")) { volcPlotly_fl <- paste0(wd, "/Reg. analysis/volcPlotly.RDS") }
if (!exists("volcPlotly")) { loadFun(volcPlotly_fl) }
#
rat_ON <- exists("ratioPlots")
if (!rat_ON) {
  if (!exists("ratioPlots_fl")) { ratioPlots_fl <- paste0(wd, "/Protein plots/ratioPlots.RDS") }
  rat_ON <- (scrptType == "noReps") && MakeRatios && file.exists(ratioPlots_fl)
  if (rat_ON) { loadFun(ratioPlots_fl) }
}
if (rat_ON) { rat_ON <- length(ratioPlots) > 0L }
#
if (!exists("heatMaps_fl")) { heatMaps_fl <- paste0(wd, "/Reg. analysis/volcPlotly.RDS") }
heatMaps_ON <- exists("plotLeatMaps")
if (!heatMaps_ON) {
  heatMaps_ON <- file.exists(heatMaps_fl)
  if (heatMaps_ON) { loadFun(heatMaps_fl) }
}
if (heatMaps_ON) {
  heatMaps_ON <- exists("plotLeatMaps") && length(plotLeatMaps)
}
#
if (!exists("dimRed_fl")) { dimRed_fl <- paste0(wd, "/Dimensionality red. plots/DimRedPlots.RDS") }
PCA_ON <- exists("dimRedPlotLy")
if (!PCA_ON) {
  PCA_ON <- file.exists(dimRed_fl)
  if (PCA_ON) { loadFun(dimRed_fl) }
}
if (PCA_ON) {
  PCA_ON <- exists("dimRedPlotLy") && (!is.null(dimRedPlotLy$PG$PCA))
}
#
if (!exists("Venn_fl")) { Venn_fl <- paste0(wd, "/Venn diagrams/Venn_plotly.RDS") }
Venn_ON <- exists("plotly_Venn")
if (!Venn_ON) {
  Venn_ON <- file.exists(Venn_fl)
  if (Venn_ON) { loadFun(Venn_fl) }
}
if (Venn_ON) {
  Venn_ON <- exists("plotly_Venn") && ("Global, LFQ" %in% names(plotly_Venn))
}
#
if (!exists("covPlots_fl")) { covPlots_fl <- paste0(wd, "/Protein plots/covPlots.RDS") }
cov_ON <- exists("covPlots")
if (!cov_ON) {
  cov_ON <- file.exists(covPlots_fl)
  if (cov_ON) { loadFun(covPlots_fl) }
}
if (cov_ON) {
  cov_ON <- length(covPlots) > 0L
}
#
if (runGSEA) {
  if (!exists("GSEA_plotly_fl")) { GSEA_plotly_fl <- paste0(wd, "/Reg. analysis/GSEA/GSEA_plotly.RDS") }
  if (!exists("GSEA_plotly")) { loadFun(GSEA_plotly_fl) }
}
#
if (runRankAbundPlots) {
  if (!exists("ggQuantLy_fl")) { ggQuantLy_fl <- paste0(wd, "/Ranked abundance/quantPlots.RDS") }
  if (!exists("ggQuantLy")) { loadFun(ggQuantLy_fl) }
}
#
if (!exists("qcBckUpFl")) { qcBckUpFl <- paste0(wd, "/QC.RDS") }
if (!exists("QC_plotLys")) { loadFun(qcBckUpFl) }
#
GO_ON <- Annotate && (enrichGO || globalGO)
if (GO_ON) {
  if (!exists("GO_plot_ly_fl")) { GO_plot_ly_fl <- paste0(wd, "/Reg. analysis/GO enrich/GO_plot_ly.RDS") }
  if (!exists("GO_plot_ly")) { loadFun(GO_plot_ly_fl) }
}

# -----------------------------------------------------------------------------------------------------
# Re-apply plotly::partial_bundle() in THIS session.
# Every object above was loaded back in via loadFun() from an RDS file written in an earlier R session.
# partial_bundle() caches its trimmed JS bundle under that earlier session's tempdir(), which no longer
# so save_html()'s copyDependencyToDir() fails trying to copy from a path like
# "AppData/Local/Temp/2/RtmpXXXXXX/plotly-cartesian-....min.js".
# Re-bundling now writes fresh files into the CURRENT tempdir(), which save_html() can find.
# NOTE: this block is currently commented out in your copy -- if that's intentional (debugging
# something else?), fine, but if it was left off by accident, save_html() will hit the same stale-
# tempdir error as before once any of these objects were loaded via loadFun() in an earlier session.
# -----------------------------------------------------------------------------------------------------
# partial_bundle_list <- \(x) {
#   if (inherits(x, "plotly")) { return(plotly::partial_bundle(x)) }
#   if (is.list(x) && identical(class(x), "list")) { return(lapply(x, partial_bundle_list)) }
#   return(x)
# }
# volcPlotly <- partial_bundle_list(volcPlotly)
# if (rat_ON) { ratioPlots <- partial_bundle_list(ratioPlots) }
# if (heatMaps_ON) { plotLeatMaps <- partial_bundle_list(plotLeatMaps) }
# if (PCA_ON) { dimRedPlotLy <- partial_bundle_list(dimRedPlotLy) }
# if (Venn_ON) { plotly_Venn <- partial_bundle_list(plotly_Venn) }
# if (cov_ON) { covPlots <- partial_bundle_list(covPlots) }
# if (GO_ON) { GO_plot_ly <- partial_bundle_list(GO_plot_ly) }
# if (runGSEA) { GSEA_plotly <- partial_bundle_list(GSEA_plotly) }
# ggQuantLy <- partial_bundle_list(ggQuantLy)
# QC_plotLys <- partial_bundle_list(QC_plotLys)

plotHght <- paste0(round(screenRes$height*0.75), "px")
nmsHtMp <- names(plotLeatMaps$Global)
if ((scrptType == "noReps") && (length(Exp) == 2L)) {
  nmsHtMp <- setdiff(nmsHtMp, "Z-scored")
}
nmsHtMp <- intersect(union("None", nmsHtMp), nmsHtMp)
plotHtMpHght <- paste0(max(c(500L, as.integer(round(vapply(nmsHtMp, \(nm) {
  plotLeatMaps$Global[[nm]]$Plot$sizingPolicy$defaultHeight
}, 1))))), "px")
plotPCAHght <- "700px"

myTabs <- if (scrptType == "noReps") {
  c("Dataset overview", Exp)
} else {
  c("Dataset overview", smplGrps, myContrasts$Contrast)
}
PTMs <- c()
if (exists("PTMs_pep")) { PTMs <- names(PTMs_pep) }
if (length(PTMs)) {
  for (Ptm in PTMs) {
    myTabs <- union(myTabs, paste0(Ptm, ": ", myContrasts$Contrast))
  }
}
nms <- myTabs
if (prot.list.Cond) {
  if (cov_ON) {
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
    if (nm == "Dataset overview") { dfltComment } else { "" }
  }, ""), nms)
}
nms_ <- setdiff(nms, names(allComments))
if (length(nms_)) {
  allComments[nms_] <- ""
  if ("Dataset overview" %in% nms_) { allComments$"Dataset overview" <- dfltComment }
}
allComments %<o% allComments[nms]

if (runRankAbundPlots) {
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
}

myExp <- if (scrptType == "noReps") { Exp } else { setNames(smplGrps, NULL) }


# Define split-file plan
# ======================
sectionFileName <- \(label) { paste0(wd, "/", dtstNm, " - ", gsub(":", "", label), ".html") }

sections <- list()
sections[["Dataset Overview"]] <- list(type = "summary", label = "Dataset Overview", file = sectionFileName("Dataset Overview"))
if ((scrptType == "noReps") || showSmplGrpTabs) {
  sections[myExp] <- lapply(myExp, \(exp) {
    list(type = "sample", label = exp, file = sectionFileName(exp), key = exp)
  })
  # for (exp in myExp) {
  #   sections[[exp]] <- list(type = "sample", label = exp, file = sectionFileName(exp), key = exp)
  # }
}
if (scrptType == "withReps") {
  sections[myContrasts$Contrast] <- lapply(myContrasts$Contrast, \(contr) {
    list(type = "contrast", label = sub(" - ", " /\n", contr), file = sectionFileName(contr), key = contr)
  })
  # for (contr in myContrasts$Contrast) {
  #   sections[[contr]] <- list(type = "contrast", label = sub(" - ", " /\n", contr), file = sectionFileName(contr), key = contr)
  # }
  for (Ptm in PTMs) {
    # NB: key is the plain contrast name (contr), not the combined id --
    # make_ctrst_tab(contr, ...) does match(contr, myContrasts$Contrast) internally,
    # which fails silently/wrongly if fed "Ptm: contr" instead of just "contr".
    sections[paste0(Ptm, ": ", myContrasts$Contrast)] <-  lapply( myContrasts$Contrast, \(contr) {
      id <- paste0(Ptm, ": ", contr)
      list(type = Ptm, label = sub(" - ", " /\n", id), file = sectionFileName(id), key = contr)
    })
    # for (contr in myContrasts$Contrast) {
    #   id <- paste0(Ptm, ": ", contr)
    #   sections[[id]] <- list(type = Ptm, label = sub(" - ", " /\n", id), file = sectionFileName(id), key = contr)
    # }
  }
}
if (prot.list.Cond) {
  sections[["Proteins of interest"]] <- list(type = "proteins", label = "Proteins of interest",
                                             file = sectionFileName("Proteins of interest"))
}
sections[["QC"]] <- list(type = "qc", label = "QC", file = sectionFileName("QC"))
sections[["Materials & Methods"]] <- list(type = "matmeth", label = "Materials & Methods", file = sectionFileName("Materials & Methods"))
sectionFiles <- vapply(sections, `[[`, "", "file")   # named vector: label -> full path

# =============================================================================
#    nav bar + asset-embedding helpers, and the per-file page builder
#    (this calls your real make_*_tab / make_*_ui functions with shiny = FALSE,
#    tables = FALSE)
# =============================================================================
make_nav_bar <- \(current_label) { #current_label <- sec$label
  # One row per group; a row with no items renders nothing (e.g. no "Proteins
  # of interest" row-1 entry when prot.list.Cond is FALSE).
  nav_row <- \(label, items, itemTexts) {
    if (!length(label)) { return(NULL) }
    if (missing(itemTexts)) { itemTexts <- items }
    tags$div(style = "margin-bottom:6px;",
             if (nzchar(label)) { tags$strong(label) },
             lapply(1L:length(items), \(i) { #i <- 1L
               itm <- items[i]
               txt <- itemTexts[i]
               fn <- basename(sections[[itm]]$file)
               if (identical(itm, current_label)) {
                 tags$span(txt, style = "font-weight:bold; margin-right:16px;")
               } else {
                 tags$a(txt, href = fn, style = "margin-right:16px;")
               }
             }))
  }
  allLbls <- names(sections)
  types <- vapply(sections, `[[`, "", "type")
  # Row 1: the fixed, un-labelled set -- order matters more than lookup speed here.
  topOrder <- c("Dataset Overview", "QC", "Proteins of interest", "Materials & Methods")
  topLbls <- intersect(topOrder, allLbls[types %in% c("summary", "qc", "proteins", "matmeth")])
  sampleLbls <- allLbls[types == "sample"]
  ptmTypes <- setdiff(unique(types), c("summary", "sample", "contrast", "proteins", "qc", "matmeth"))
  contrastLbls <- allLbls[types == "contrast"]
  contrastLbls2 <- lapply(sections, \(sec){
    if ("key" %in% names(sec)) { return(sec$key) }
    return()
  })
  # Any type that isn't one of the six fixed kinds is a PTM name (see sections
  # construction above); unique() preserves the PTMs' original order since the
  # sections loop nests PTM outside contrast.
  tags$nav(style = "padding:10px 0 16px 0; border-bottom:1px solid #ddd; margin-bottom:20px;",
           nav_row("", topLbls),
           nav_row(paste0("Sample", c(" group", "")[(scrptType == "noReps")+1L], "s: "),
                   sampleLbls),
           nav_row(paste0(c("C", "Protein group c"), "ontrasts: ")[(length(ptmTypes) > 0L)+1L],
                   contrastLbls),
           lapply(ptmTypes, \(ptm) { nav_row(paste0(ptm, " contrasts: "), allLbls[types == ptm],
                                             contrastLbls2[types == ptm]) }))
}
# Re-run of the same "reset <select data-default> on load" defensive script the original single-file export used.
select_default_script <- tags$script(HTML(
  "document.addEventListener('DOMContentLoaded', function() {
  document.querySelectorAll('select[data-default]').forEach(function(sel) {
    sel.value = sel.dataset.default;
    sel.dispatchEvent(new Event('change'));
  });
});"))

# Self-contained-ness:
# inline every local <script src=...> / <link href=...> referenced from a saved html file's <head>.
# Refactored from the original bottom-of-script embedding logic so it can be called once per output file.
embed_assets <- \(html_path) {
  h1 <- readr::read_lines(html_path)
  h2 <- h1
  rg1 <- grep("</?head>", h1) + c(1L, -1L)
  rg1 <- rg1[1L]:rg1[2L]
  hd1 <- data.frame(original = h1[rg1])
  hd1$new <- hd1$original
  read_asset <- \(path) { readChar(path, nchars = file.info(path)$size, useBytes = TRUE) }
  inline_script <- \(path) {
    txt <- paste(read_asset(path), collapse = "")
    txt <- gsub("</script", "<\\/script", txt, ignore.case = TRUE)
    paste0("<script>\n", txt, "\n</script>")
  }
  inline_css <- \(path) { paste0("<style>\n", read_asset(path), "\n</style>") }
  gs <- grep("^ *<script src=\"", hd1$original)
  if (length(gs)) {
    hd1$new[gs] <- vapply(sub("\".*", "", sub("^ *<script src=\"", paste0(wd, "/"), hd1$original[gs])),
                          inline_script, "")
  }
  gc <- grepl("^ *<link href=\"", hd1$original)
  if (any(gc)) {
    hd1$new[gc] <- vapply(sub("\".*", "", sub("^ *<link href=\"", paste0(wd, "/"), hd1$original[gc])),
                          inline_css, "")
  }
  h2[rg1] <- hd1$new
  write(h2, html_path)
  return()
}

# NB fix (this was the do.call(switch, ...) bug): do.call(switch, list(...)) is
# NOT equivalent to switch() -- list(...) evaluates every element eagerly before
# do.call ever runs, so every branch's make_*_tab() call fired regardless of
# sec$type (that's what produced "attempt to select less than one element" --
# make_ctrst_tab() running with sec$key == NULL whenever sec$type wasn't
# "contrast"). A plain switch() only evaluates its one matched branch, so it's
# lazy on its own -- no dispatch table needed. Since sections already sets
# `type = Ptm` for PTM-contrast entries (see the sections loop above), a single
# if/else fallback for "sec$type is a PTM name" handles any number of PTMs
# without a dynamic branch list.
build_section <- \(sec) { #sec <- sections[[1L]]
  content <- if (sec$type %in% c("summary", "sample", "contrast", "proteins", "qc", "matmeth")) {
    switch(sec$type,
           summary  = make_strt_tab(shiny = FALSE),
           sample   = make_smpl_tab(sec$key, shiny = FALSE),
           contrast = make_ctrst_tab(sec$key, shiny = FALSE),
           proteins = make_prot_tab(dfltProt, allProt, shiny = FALSE),
           qc       = make_QC_tab(shiny = FALSE),
           matmeth  = make_matmet_tab(matmethTxt, shiny = FALSE))
  } else {
    # sec$type is a PTM name here (e.g. "Phospho"); sec$key is the plain
    # contrast name after the sections-loop fix above.
    make_ctrst_tab(sec$key, shiny = FALSE, ptm = sec$type)
  }
  bslib::page_fluid(tags$head(myCSS),
                    select_default_script,
                    report_header,
                    make_nav_bar(sec$label),
                    content)
}
if (runGSEA) {
  GSEA_plotNms <- c("GSEA dotplot", "GSEA enrichment map", "GSEA ridge plot", "GSEA category net plot")
  wh_GSEA <- setNames(lapply(myContrasts$Contrast, \(contr) {
    c("PG", PTMs)[which(vapply(c("PG", PTMs), \(GSEA_nm) { #GSEA_nm <- c("PG", PTMs)[1L] #GSEA_nm <- c("PG", PTMs)[2L]
      sum(vapply(GSEA_plotNms, \(plotNm) { #plotNm <- GSEA_plotNms[1L]
        p <- GSEA_plotly$standard[[GSEA_nm]][[plotNm]][[contr]]
        (!is.null(p)) && (inherits(p, "plotly"))
      }, TRUE))
    }, 1L) > 0L)]
  }), myContrasts$Contrast)
}

appPage <- 1L
appNm <- "Edit report"
ui <- fluidPage(useShinyjs(),
                extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
                tags$head(myCSS),
                titlePanel(tag("u", appNm), appNm),
                br(),
                fluidRow(column(4L, h2(dtstNm), br()),
                         column(8L,
                                actionBttn("xprtBtn", " export final html report", icon = icon("file-export"), color = "success", style = "pill"),
                                br(), br(),
                                uiOutput("xprtMsg"))),
                br(),
                uiOutput("myUI"),
                br(), br())
server <- \(input, output, session) {
  QUANT <- reactiveVal(dfltQuant)
  XPRTMSG <- reactiveVal(NULL)
  MYPROT <- reactiveVal(dfltProt)
  SAMPLE <- reactiveVal(myExp[1L])
  NORMMETH <- reactiveVal("None")
  ALLCOMMENTS <- reactiveVal(allComments)
  #
  output$xprtMsg <- renderUI(XPRTMSG())
  output$myUI <- renderUI(make_ui())
  if (heatMaps_ON) { output$heatMap <- renderPlotly(plotLeatMaps$Global[[NORMMETH()]]$Plot) }
  if (PCA_ON) { output$PCA <- renderPlotly(dimRedPlotLy$PG$PCA) }
  if (Venn_ON) { output$Venn <- renderPlotly(plotly_Venn$`Global, LFQ`) }
  #
  lapply(myExp, \(exp) {
    exp_ <- safe_id(exp)
    exp2 <- if (scrptType == "noReps") { exp } else { names(smplGrps)[match(exp, smplGrps)] }
    output[[paste0("cmmnt_", exp_)]] <- renderUI(make_comment_ui(exp))
    idQ <- paste0("quant_", exp)
    idQLy <- paste0("quantLy_", exp)
    idSmplPGTbl <- paste0("PG_table_", exp)
    output[[idQLy]] <- renderPlotly(ggQuantLy[[input[[idQ]]]][[exp2]]$plotly)
    output[[idSmplPGTbl]] <- renderUI(make_smpl_tbl_ui(exp))
  })
  #
  if (scrptType == "withReps") {
    lapply(myContrasts$Contrast, \(contr) { #contr <- myContrasts$Contrast[1L]
      contr_ <- safe_id(contr)
      volcIDs <- list(PG = paste0(contr_, "_volcPlot"))
      commentIDs <- list(PG = paste0("cmmnt_", contr_))
      goIDs <- list(PG = paste0(contr_, "_GObars"))
      #tblIDs <- list(PG = ... # For now we do not show tables in the Shiny app
      if (runGSEA) {
        GSEA_IDs <- list(PG = paste0(contr_, "_GSEA", as.character(1L:4L)))
      }
      if (length(PTMs)) {
        volcIDs[PTMs] <- lapply(PTMs, \(Ptm) { paste0(Ptm, "_", volcIDs$PG) })
        commentIDs[PTMs] <- lapply(PTMs, \(Ptm) { paste0("cmmnt_", Ptm, "_", contr_) })
        goIDs[PTMs] <- lapply(PTMs, \(Ptm) { paste0(Ptm, "_", goIDs$PG) })
        #tblIDs[PTMs] <- lapply(PTMs, \(Ptm) { ... # For now we do not show tables in the Shiny app
        if (runGSEA) {
          GSEA_IDs[PTMs] <- lapply(PTMs, \(Ptm) { paste0(Ptm, "_", GSEA_IDs$PG) })
        }
      }
      output[[volcIDs$PG]] <- renderPlotly(volcPlotly$"t-test"[[paste0("Volcano plot ", contr)]]$Plot)
      output[[commentIDs$PG]] <- renderUI(make_comment_ui(contr))
      if (enrichGO) { output[[goIDs$PG]] <- renderPlotly(GO_plot_ly$PG$"t-test"[[contr]]$Bar) }
      saintIDs <- c(paste0("SAINTexpress volcano plot ", contr), paste0(contr_, c("_SAINT_volcPlot", "_SAINT_GObars")))
      if (saintExprs && (saintIDs[1L] %in% names(volcPlotly$SAINTexpress))) {
        output[[saintIDs[2L]]] <- renderPlotly(volcPlotly$SAINTexpress[[saintIDs[1L]]]$Plot)
        if (enrichGO) { output[[saintIDs[3L]]] <- renderPlotly(GO_plot_ly$Prot$SAINTexpress[[contr]]$Bar) }
      }
      lapply(PTMs, \(Ptm) { #Ptm <- PTMs[1L]
        output[[volcIDs[[Ptm]]]] <- renderPlotly(volcPlotly[[paste0(Ptm, " t-test")]][[paste0(Ptm, " volcano plot ", contr)]]$Plot)
        output[[commentIDs[[Ptm]]]] <- renderUI(make_comment_ui(paste0(Ptm, ": ", contr)))
        if (enrichGO) { output[[goIDs[[Ptm]]]] <- renderPlotly(GO_plot_ly[[Ptm]]$"t-test"[[contr]]$Bar) }
        # No SAINTexpress analysis with PTMs
      })
      if (runGSEA && length(wh_GSEA[[contr]])) {
        # We will deal with PG and PTMs in one loop
        # We only want to draw space for the plots if they are available
        lapply(wh_GSEA[[contr]], \(GSEA_nm) {
          output[[GSEA_IDs[[GSEA_nm]][[1L]]]] <- renderPlotly(GSEA_plotly$standard[[GSEA_nm]][[GSEA_plotNms[1L]]][[contr]])
          output[[GSEA_IDs[[GSEA_nm]][[2L]]]] <- renderPlotly(GSEA_plotly$standard[[GSEA_nm]][[GSEA_plotNms[2L]]][[contr]])
          output[[GSEA_IDs[[GSEA_nm]][[3L]]]] <- renderPlotly(GSEA_plotly$standard[[GSEA_nm]][[GSEA_plotNms[3L]]][[contr]])
          output[[GSEA_IDs[[GSEA_nm]][[4L]]]] <- renderPlotly(GSEA_plotly$standard[[GSEA_nm]][[GSEA_plotNms[4L]]][[contr]])
        })
      }
    })
  }
  if (globalGO && (!is.null(GO_plot_ly$PG$Dataset$`Observed dataset`$Bar))) {
    output$GO_enrich_Dataset <- renderPlotly(GO_plot_ly$PG$Dataset$`Observed dataset`$Bar)
  }
  #
  observeEvent(input$myHeatMap, { NORMMETH(input$myHeatMap) })
  sapply(names(allComments), \(id) {
    id_ <- paste0("comment_", safe_id(id))
    observeEvent(input[[id_]], {
      #cat(input[[id_]], "\n", tmp[[id]], "\n\n")
      tmp <- ALLCOMMENTS()
      tmp[[id]] <- input[[id_]]
      ALLCOMMENTS(tmp)
    })
  })
  sapply(myExp, \(exp) {
    exp2 <- if (scrptType == "noReps") { exp } else { names(smplGrps)[match(exp, smplGrps)] }
    idQ <- paste0("quant_", exp)
    idQLy <- paste0("quantLy_", exp)
    observeEvent(input[[idQ]], {
      dq <- QUANT()
      dq[exp] <- input[[idQ]]
      QUANT(dq)
      assign("dfltQuant", dq, envir = .GlobalEnv)
      output[[idQLy]] <- renderPlotly(ggQuantLy[[input[[idQ]]]][[exp2]]$plotly)
    })
  })
  if (prot.list.Cond) {
    if (scrptType == "noReps") { output$ratioPlot <- renderPlotly(ratioPlots[[MYPROT()]]) }
    output$coverPlot <- renderPlotly({
      p <- covPlots[[MYPROT()]]$logInt[[SAMPLE()]]
      if (is.null(p)) {
        return(plot_ly(type = "scatter", mode = "markers") |>
                 layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
                        annotations = list(list(text = "No identifications for this protein in this sample!",
                                                x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE))))
      }
      return(p)
    })
    output$protComment <- renderUI(make_comment_ui(MYPROT()))
    if (peptidoTst) {
      output$protPep <- renderUI({ make_smpl_tbl_ui(tab = "All peptidoforms", filt = MYPROT()) })
    }
    if (cov_ON && (length(allProt) > 1L)) { observeEvent(input$myProtein, { MYPROT(input$myProtein) }) }
    if (length(myExp) > 1L) { observeEvent(input$mySample, { SAMPLE(input$mySample) }) }
  }
  observeEvent(input$QC, {
    output$QCplotLy <- renderPlotly(QC_plotLys[[input$QC]])
    output$QCtxt <- renderUI(make_comment_ui(input$QC))
  })
  observeEvent(input$MatMet_SamplePrep, {
    txt <- matmethTxt; txt[matmethSections[1L]] <- input$MatMet_SamplePrep
    assign("matmethTxt", txt, envir = .GlobalEnv)
  })
  observeEvent(input$MatMet_LCMS, {
    txt <- matmethTxt; txt[matmethSections[2L]] <- input$MatMet_LCMS
    assign("matmethTxt", txt, envir = .GlobalEnv)
  })
  observeEvent(input$MatMet_DataAnalysis, {
    txt <- matmethTxt; txt[matmethSections[3L]] <- input$MatMet_DataAnalysis
    assign("matmethTxt", txt, envir = .GlobalEnv)
  })
  # ------------------------------------------------------------------------------
  # Export writes one static, table-free file per section (via build_section(),
  # which calls the make_*_tab functions with shiny = FALSE, tables = FALSE),
  # instead of one combined file previously
  # ------------------------------------------------------------------------------
  observeEvent(input$xprtBtn, {
    assign("allComments", ALLCOMMENTS(), envir = .GlobalEnv)
    # lapply(names(allComments), \(id) {
    #   id_ <- paste0("comment_", safe_id(id))
    #   print(paste0(id, " = ", input[[id_]], " | ", allComments[id]))
    # })
    XPRTMSG(em("Exporting .html report(s), this will take a few minutes...", style = "color:green", .noWS = "outside"))
    later::later(\() {
      # Commented attempt at parallelization: this won't work!
      # We are already extending time in this session to reload the plotly widget lists stored as .RDS, and these are massive.
      # We risk not saving much time AND exploding our RAM usage beyong possible.
      # source(parSrc)
      # parallel::clusterExport(parClust, obj2Xport, envir = environment()) # This would have been insufficient, we need to export more objects
      # invisible(parallel::parLapply(parClust, sections, \(sec) {
      #   page <- htmltools::browsable(build_section(sec))
      #   htmltools::save_html(page, sec$file, libdir = "lib")
      # }))
      cat("Writing html report pages:\n")
      for (sec in sections) { #sec <- sections[1L]
        cat(" -", sec$file, "\n")
        page <- htmltools::browsable(build_section(sec))
        htmltools::save_html(page, sec$file, libdir = "lib")
      }
      assign("appRunTst", TRUE, envir = .GlobalEnv)
      stopApp()
    }, 0.1)
  })
  session$onSessionEnded(\() { stopApp() })
}
runKount <- 0L
if (exists("appRunTst")) { rm(appRunTst) }
while ((!runKount) || (!exists("appRunTst")) || (!all(file.exists(sectionFiles)))) {
  g <- shiny:::.globals
  g$appState <- NULL
  eval(parse(text = run_App), envir = .GlobalEnv)
  shinyCleanup()
  runKount <- runKount + 1L
}

# Make each file fully portable (embed its local lib/ assets), then clean up
# ==========================================================================
cat("Embedding dependencies...\n")
source(parSrc)
parallel::clusterExport(parClust, "wd", envir = environment())
invisible(parallel::parLapply(parClust, sectionFiles, embed_assets))
removeDirectory(paste0(wd, "/lib"), TRUE, FALSE)

# Save the final materials and methods to a local file
# ====================================================
cat("Writing materials and methods template...\n")
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
  MatMet %<o% MatMet
  print(MatMet, target = MatMetFl)
}, silent = TRUE)
if (inherits(tst, "try-error")) {
  warning("Couldn't write materials and methods template, investigate...")
}
unlink(tmpSrc)

source(bckpSrc)

cat("Done!\n\n")

# try({
#   rm(ggQuantLy, plotLeatMaps, dimRedPlotLy, plotly_Venn, QC_plotLys)
# }, silent = TRUE)
