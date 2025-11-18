# Use browser or RStudio for Shiny apps?
screenRes %<o% rpanel::rp.screenresolution()
shinyOpts %<o% c("RStudio", "System default browser")
#shinyOpt %<o% svDialogs::dlg_list(shinyOpts, title = "Open Shiny apps in...")$res
shinyOpt %<o% "RStudio"
if (shinyOpt == shinyOpts[1]) {
  #runApp %<o% "print(shiny::shinyApp(ui, server, options = list(height = screenRes$height, width = screenRes$width)))"
  runApp %<o% c("myApp <- shiny::shinyApp(ui, server, options = list(height = screenRes$height, width = \"100%\"))",
                "shiny::runApp(myApp)")
  #myViewer %<o% shiny::dialogViewer("Viewer", width = screenRes$width, height = screenRes$height)
  # runApp %<o% c("myApp <- shiny::shinyApp(ui, server)",
  #               "shiny::runGadget(myApp, viewer = myViewer, stopOnCancel = FALSE)")
}
if (shinyOpt == shinyOpts[2]) {
  runApp %<o% c("myApp <- shiny::shinyApp(ui, server, options = list(height = \"100%\", width = \"100%\", launch.browser = TRUE))",
                "print(myApp)")
}
jsToggleFS %<o% "shinyjs.toggleFullScreen = function() {
        var element = document.documentElement,
          enterFS = element.requestFullscreen || element.msRequestFullscreen || element.mozRequestFullScreen || element.webkitRequestFullscreen,
          exitFS = document.exitFullscreen || document.msExitFullscreen || document.mozCancelFullScreen || document.webkitExitFullscreen;
        if (!document.fullscreenElement && !document.msFullscreenElement && !document.mozFullScreenElement && !document.webkitFullscreenElement) {
          enterFS.call(document);
        } else {
          exitFS.call(document);
        }
      }"

# Create Word formatting styles (package used = officer):
WrdFrmt %<o% list()
WrdFrmt$Main_title <- officer::fp_text(bold = TRUE, font.family = "Calibri", font.size = 12, underlined = TRUE) # For Main title
WrdFrmt$Section_title <- officer::fp_text(bold = TRUE, font.family = "Calibri") # For section title
WrdFrmt$Section_title_ital <- officer::fp_text(bold = TRUE, italic = TRUE, font.family = "Calibri") # For italic title
WrdFrmt$Body_text <- officer::fp_text(font.family = "Calibri") # For normal body text
WrdFrmt$Body_text_ital <- officer::fp_text(italic = TRUE, font.family = "Calibri") # For italic body text
WrdFrmt$Body_text_under <- officer::fp_text(underlined = TRUE, font.family = "Calibri") # For underlined body text
WrdFrmt$Template_text <- officer::fp_text(font.family = "Calibri", color = "red") # For empty template body text
WrdFrmt$just <- officer::fp_par(text.align = "justify")
WrdFrmt$left <- officer::fp_par(text.align = "left")
ReportCalls %<o% list(Calls = list("read_docx()",
                                   "body_add_par(Report, \"\", style = \"Normal\")",
                                   "body_add_fpar(Report, fpar(ftext(paste0(\"Data analysis - \", start_date), prop = WrdFrmt$Main_title), fp_p = WrdFrmt$just))",
                                   "body_add_par(Report, \"\", style = \"Normal\")",
                                   "body_add_par(Report, \"\", style = \"Normal\")",
                                   "body_add_fpar(Report, fpar(ftext(\"Inputs\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$left))"),
                      Objects = list(),
                      Plots = list())
AddPlot2Report %<o% function(RCName = "ReportCalls", Plot = plot, Title = ttl, Space = TRUE, Jpeg = TRUE, Dir = dir) {
  try({
    reportCalls <- get(RCName, envir = .GlobalEnv)
    if (Jpeg) {
      Path <- paste0(gsub("/+$", "", normalizePath(Dir, winslash = "/")), "/", Title, ".jp", c("", "e"), "g")
      w <- which(file.exists(Path))
      stopifnot(length(w) > 0)
      Path <- Path[w[1]]
    } else {
      Title2 <- Title
      kount <- 0
      while (Title2 %in% names(reportCalls$Plots)) {
        kount <- kount + 1
        Title2 <- paste0(Title, "_", kount)
      }
    }
    reportCalls$Calls <- append(reportCalls$Calls,
                                paste0("body_add_fpar(Report, fpar(ftext(\"", gsub(":", "_", Title),
                                       "\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))"))
    if (Jpeg) {
      reportCalls$Calls <- append(reportCalls$Calls, paste0("body_add_img(Report, \"", Path, "\", height = 6, width = 6)"))
    } else {
      reportCalls$Calls <- append(reportCalls$Calls,
                                  paste0("body_add_gg(Report, ", RCName, "$Plots[[\"", Title2, "\"]])"))
      reportCalls$Plots[[Title2]] <- Plot
    }
    if (Space) { reportCalls$Calls <- append(reportCalls$Calls, "body_add_par(Report, \"\", style = \"Normal\")") }
    return(reportCalls)
  }, silent = TRUE)
}
AddMsg2Report %<o% function(RCName = "ReportCalls", Msg = msg, Print = TRUE, Warning = FALSE, Space = TRUE, Offset = FALSE) {
  try({
    reportCalls <- get(RCName, envir = .GlobalEnv)
    if (Print+Warning) { if (Warning) { warning(Msg) } else { cat(paste0(Msg, "\n")) } }
    ReportCalls$Calls <- append(ReportCalls$Calls,
                                paste0("body_add_fpar(Report, fpar(ftext(\"", c("", "     ")[Offset+1], Msg,
                                       "\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))"))
    if (Space) { reportCalls$Calls <- append(reportCalls$Calls, "body_add_par(Report, \"\", style = \"Normal\")") }
    return(reportCalls)
  }, silent = TRUE)
}
AddSpace2Report %<o% function(RCName = "ReportCalls") {
  try({
    reportCalls <- get(RCName, envir = .GlobalEnv)
    reportCalls <- get(RCName, envir = .GlobalEnv)
    reportCalls$Calls <- append(reportCalls$Calls, "body_add_par(Report, \"\", style = \"Normal\")")
    return(reportCalls)
  }, silent = TRUE) 
}

# Create Excel formatting styles (package used = openxlsx):
# These styles are used later down when writing Excel tables
Styles %<o% list(IDs = openxlsx::createStyle(textDecoration = c("bold", "italic"), numFmt = "TEXT",
                                             fontSize = 11, valign = "center"))
# - Peptide and evidence counts and IDs
pepcolours %<o% colorRampPalette(c("thistle2", "thistle1"))(8)
evcolours %<o% colorRampPalette(c("slategray2", "slategray1"))(8)
spcolours %<o% colorRampPalette(c("wheat2", "wheat1"))(8)
Styles[["Global Pep. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = pepcolours[1], fgFill = pepcolours[1],
                                                     numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Global Pep. counts"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = pepcolours[2], fgFill = pepcolours[2],
                                                        numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Pep. IDs"]] <- openxlsx::createStyle(bgFill = pepcolours[3], fgFill = pepcolours[3],
                                              numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Pep. counts"]] <- openxlsx::createStyle(bgFill = pepcolours[4], fgFill = pepcolours[4],
                                                 numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Global Ev. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = evcolours[1], fgFill = evcolours[1],
                                                    numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Ev. IDs"]] <- openxlsx::createStyle(bgFill = evcolours[3], fgFill = evcolours[3],
                                             numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Ev. counts"]] <- openxlsx::createStyle(bgFill = evcolours[4], fgFill = evcolours[4],
                                                numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Global Spec. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = spcolours[1], fgFill = spcolours[1],
                                                      numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Global Spec. counts"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = spcolours[2], fgFill = spcolours[2],
                                                         numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Spec. IDs"]] <- openxlsx::createStyle(bgFill = spcolours[3], fgFill = spcolours[3],
                                               numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Spec. counts"]] <- openxlsx::createStyle(bgFill = spcolours[4], fgFill = spcolours[4],
                                                  numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Global Biot. Pep. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = pepcolours[5], fgFill = pepcolours[1],
                                                           numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Global Biot. Pep. counts"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = pepcolours[6], fgFill = pepcolours[2],
                                                              numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Biot. Pep. IDs"]] <- openxlsx::createStyle(bgFill = pepcolours[5], fgFill = pepcolours[7],
                                                    numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Biot. Pep. counts"]] <- openxlsx::createStyle(bgFill = pepcolours[6], fgFill = pepcolours[8],
                                                       numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Global Biot.Ev. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = evcolours[5], fgFill = evcolours[1],
                                                         numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Biot. Ev. IDs"]] <- openxlsx::createStyle(bgFill = evcolours[5], fgFill = evcolours[7],
                                                   numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Biot. Ev. counts"]] <- openxlsx::createStyle(bgFill = evcolours[6], fgFill = evcolours[8],
                                                      numFmt = "COMMA", fontSize = 11, valign = "center")
# (Placeholders, not used currently)
Styles[["Global Biot. Spec. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = spcolours[5], fgFill = spcolours[1],
                                                            numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Global Biot. Spec. counts"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = spcolours[6], fgFill = spcolours[2],
                                                               numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Biot. Spec. IDs"]] <- openxlsx::createStyle(bgFill = spcolours[5], fgFill = spcolours[7],
                                                     numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Biot. Spec. counts"]] <- openxlsx::createStyle(bgFill = spcolours[6], fgFill = spcolours[8],
                                                        numFmt = "COMMA", fontSize = 11, valign = "center")
#
Styles[["Biot. Pep. %"]] <- "Biot. Pep. %"
## - Individual Expr
Styles[["Individual Expr"]] <- "Individual Expr"
## - Summary Expr
Styles[["Summary Expr"]] <- "Summary Expr"
## - Summary Expr: standard errors
Styles[["Summary Expr: standard errors"]] <- "Summary Expr: standard errors"
## - Proteome Ruler
Styles[["Proteome Ruler"]] <- "Proteome Ruler"
## - Individual Ratios
Styles[["Individual Ratios"]] <- "Individual Ratios"
## - Summary Ratios
Styles[["Summary Ratios"]] <- "Summary Ratios"
## - Summary Ratios: P-values, standard errors and significance
Styles[["Summary Ratios: standard errors"]] <- "Summary Ratios: standard errors"
Styles[["P-values"]] <- "P-values"
# - Significant
Styles[["Significant"]] <- openxlsx::createStyle(textDecoration = "italic", bgFill = "bisque", fgFill = "bisque",
                                                 numFmt = "TEXT", fontSize = 11, valign = "center")
# - Regulated
Styles[["Regulated"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = "lightgoldenrod1", fgFill = "lightgoldenrod1",
                                               numFmt = "TEXT", fontSize = 11, valign = "center")
# F-test
Styles[["F-test summary Ratios"]] <- Styles[["Summary Ratios"]]
Styles[["F-test P-values"]] <- Styles[["P-values"]]
Styles[["F-test significant"]] <- Styles[["Significant"]]
Styles[["F-test regulated"]] <- Styles[["Regulated"]]
# Localisation
Styles[["Localisation"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = "azure", fgFill = "azure",
                                                  numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["SSDs"]] <- "SSDs"
Styles[["Mean SSDs"]] <- "Mean SSDs"
Styles[["SSDs P-values"]] <- Styles[["P-values"]]
Styles[["SSDs significant"]] <- Styles[["Significant"]]
Styles[["Re-localized"]] <- Styles[["Regulated"]]
Styles[["Marker"]] <- openxlsx::createStyle(textDecoration = "bold", numFmt = "TEXT", fontSize = 11, valign = "center")
# - Annotations
annot %<o% c("InterPro", "Pfam", "PIRSF", "PROSITE")
AnnotTbl %<o% data.frame(Name = c("GO", "Taxonomy", annot, "EMBL", "Other"),
                         Colour = colorRampPalette(c("cyan", "cyan3"))(8), check.names = FALSE)
for (i in 1:nrow(AnnotTbl)) {
  Styles[[paste0(AnnotTbl$Name[i], " annotations")]] <- openxlsx::createStyle(bgFill = AnnotTbl$Colour[i], fgFill = AnnotTbl$Colour[i],
                                                                              numFmt = "TEXT", fontSize = 11, valign = "center")
}
# - PEP
Styles[["PEP"]] <- "PEP"
# - Filters (contaminants, only identified by site)
Styles[["Filters"]] <- openxlsx::createStyle(textDecoration = "italic", halign = "center", valign = "center", wrapText = TRUE,
                                             numFmt = "TEXT", fontSize = 12)
# - Clusters
Styles[["Cluster"]] <- "Cluster"
# - Header
Styles[["Header 1"]] <- openxlsx::createStyle(fontColour = "white", bgFill = "#5b9bd5", fgFill = "#5b9bd5",
                                              border = "left", borderStyle = "thin", borderColour = "dimgrey",
                                              textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE,
                                              numFmt = "TEXT", fontSize = 12)
Styles[["Header 2"]] <- openxlsx::createStyle(fontColour = "white", bgFill = "#5b9bd5", fgFill = "#5b9bd5",
                                              border = "left", borderStyle = "thin", borderColour = "dimgrey",
                                              textDecoration = c("bold", "underline"), halign = "center", valign = "center", wrapText = TRUE,
                                              numFmt = "TEXT", fontSize = 12)
Styles[["LeftBorder"]] <- openxlsx::createStyle(border = "left", borderStyle = "thick")
Styles[["RightBorderThin"]] <- openxlsx::createStyle(border = "right", borderStyle = "thin", borderColour = "dimgrey")
Styles[["UpperBorder"]] <- openxlsx::createStyle(border = "top", borderStyle = "thick")
Styles[["NoUpperBorder"]] <- openxlsx::createStyle(border = "top", borderStyle = "none")
# - Proteins in list
Styles[["Protein list"]] <- openxlsx::createStyle(textDecoration = "bold", fontColour = "brown")
# - Negative filter
Styles[["Negative filter"]] <- "Negative filter"
# - Summary tab
Styles[["Interlines"]] <- openxlsx::createStyle(fontSize = 12, halign = "left", valign = "center", textDecoration = "underline2", wrapText = FALSE)
Styles[["Small table - TEXT"]] <- openxlsx::createStyle(fontSize = 10, halign = "center", valign = "center", wrapText = FALSE, numFmt = "TEXT")
Styles[["Small table - NUMBER"]] <- openxlsx::createStyle(fontSize = 10, halign = "center", valign = "center", wrapText = FALSE, numFmt = "NUMBER")
Styles[["Header 3"]] <- openxlsx::createStyle(bgFill = "darkgrey", fgFill = "blue",
                                              fontSize = 11, valign = "center", halign = "center", wrapText = TRUE,
                                              border = "TopBottomLeftRight", borderColour = "black", borderStyle = "thick", numFmt = "TEXT")
# Also colour scales and text decorations for conditional formatting
ColScaleList %<o% list(`Individual Expr` = c("red", "green"),
                       `Summary Expr` = c("red", "green"),
                       `Summary Expr: standard errors` = c("firebrick1", "white"),
                       `Proteome Ruler` = c("red", "green"),
                       `Individual Ratios` = c("blue", "grey95", "red"),
                       `Summary Ratios` = c("blue", "grey95", "red"),
                       `Summary Ratios: standard errors` = c("firebrick1", "white"),
                       `P-values` = c("cyan", "brown1"),
                       `F-test summary Ratios` = c("blue", "grey95", "red"),
                       `F-test P-values` = c("cyan", "brown1"),
                       `PEP` = c("brown1", "cyan"),
                       `Cluster` = c("darkolivegreen1", "darkorange"),
                       `SSDs` = c("grey95", "green"),
                       `Mean SSDs` = c("grey95", "green"),
                       `FDR` = c("cyan", "brown1"),
                       `AvgP` = c("grey95", "green"),
                       `MaxP` = c("grey95", "green"),
                       `topoAvgP` = c("grey95", "green"),
                       `topoMaxP` = c("grey95", "green"),
                       `SaintScore` = c("grey85", "aquamarine1"),
                       `OddsScore` = c("grey90", "aquamarine3"))
ContainsList %<o% list(`Negative filter` = list(rule = "+",
                                                style = openxlsx::createStyle(bgFill = "red",
                                                                              fgFill = "red",
                                                                              fontColour = "black",
                                                                              textDecoration = "bold",
                                                                              halign = "center")),
                       `GO-columns` = list(rule = "+",
                                           style = openxlsx::createStyle(bgFill = "lightblue",
                                                                         fgFill = "lightblue",
                                                                         fontColour = "black",
                                                                         textDecoration = "bold",
                                                                         halign = "center")) # Currently not used... but actually that's ok
)
DecoList %<o% list(`Summary Expr` = "bold",
                   `Summary Expr: standard errors` = "italic",
                   `Summary Ratios` = "bold",
                   `Summary Ratios: standard errors` = "italic",
                   `P-values` = "italic",
                   `F-test summary Ratios` = "bold",
                   `F-test P-values` = "italic",
                   `PEP` = "italic",
                   `Cluster` = "bold",
                   `Mean SSDs` = "bold",
                   `FDR` = "italic",
                   `SaintScore` = "italic",
                   `OddsScore` = "italic")
SignifList %<o% list(`Individual Expr` = 5,
                     `Summary Expr` = 5,
                     `Summary Expr: standard errors` = 5,
                     `Proteome Ruler` = 5,
                     `Individual Ratios` = 5,
                     `Summary Ratios` = 5,
                     `Summary Ratios: standard errors` = 5,
                     `P-values` = 5,
                     `F-test summary Ratios` = 5,
                     `F-test P-values` = 5,
                     `PEP` = 5,
                     `SSDs` = 5,
                     `Mean SSDs` = 5,
                     `FDR` = 5,
                     `AvgP` = 5,
                     `MaxP` = 5,
                     `topoAvgP` = 5,
                     `topoMaxP` = 5,
                     `SaintScore` = 5,
                     `OddsScore` = 5)
ColList %<o% list(`Cluster` = "black")
HAlignList %<o% list(`Cluster` = "center")
ExcelMax %<o% 32767
