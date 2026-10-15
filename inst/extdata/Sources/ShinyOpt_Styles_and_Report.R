# Use browser or RStudio for Shiny apps?
screenRes %<o% rpanel::rp.screenresolution()
shinyOpts %<o% c("RStudio", "System default browser")
#shinyOpt %<o% svDialogs::dlg_list(shinyOpts, title = "Open Shiny apps in...")$res
shinyOpt %<o% "RStudio"
if (shinyOpt == shinyOpts[1L]) {
  #run_App %<o% "print(shiny::shinyApp(ui, server, options = list(height = screenRes$height, width = screenRes$width)))"
  run_App %<o% c("myApp <- shiny::shinyApp(ui, server, options = list(height = screenRes$height, width = \"100%\"))",
                 "shiny::runApp(myApp)")
  # myViewer %<o% shiny::dialogViewer("Viewer", width = screenRes$width, height = screenRes$height)
  # run_App %<o% c("myApp <- shiny::shinyApp(ui, server)",
  #                "shiny::runGadget(myApp, viewer = myViewer, stopOnCancel = FALSE)")
}
if (shinyOpt == shinyOpts[2L]) {
  run_App %<o% c("myApp <- shiny::shinyApp(ui, server, options = list(height = \"100%\", width = \"100%\", launch.browser = TRUE))",
                 "print(myApp)")
}

# Create Excel formatting styles (package used = openxlsx):
# These styles are used later down when writing Excel tables
Styles %<o% list(IDs = openxlsx::createStyle(textDecoration = c("bold", "italic"), numFmt = "TEXT",
                                             fontSize = 11L, valign = "center"))
# - Peptide and evidence counts and IDs
pepcolours %<o% colorRampPalette(c("thistle2", "thistle1"))(8L)
evcolours %<o% colorRampPalette(c("slategray2", "slategray1"))(8L)
spcolours %<o% colorRampPalette(c("wheat2", "wheat1"))(8L)
Styles[["Global Pep. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = pepcolours[1L], fgFill = pepcolours[1L],
                                                     numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Global Pep. counts"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = pepcolours[2L], fgFill = pepcolours[2L],
                                                        numFmt = "COMMA", fontSize = 11L, valign = "center")
Styles[["Pep. IDs"]] <- openxlsx::createStyle(bgFill = pepcolours[3L], fgFill = pepcolours[3L],
                                              numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Pep. counts"]] <- openxlsx::createStyle(bgFill = pepcolours[4L], fgFill = pepcolours[4L],
                                                 numFmt = "COMMA", fontSize = 11L, valign = "center")
Styles[["Global Ev. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = evcolours[1L], fgFill = evcolours[1L],
                                                    numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Ev. IDs"]] <- openxlsx::createStyle(bgFill = evcolours[3L], fgFill = evcolours[3L],
                                             numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Ev. counts"]] <- openxlsx::createStyle(bgFill = evcolours[4L], fgFill = evcolours[4L],
                                                numFmt = "COMMA", fontSize = 11L, valign = "center")
Styles[["Global Spec. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = spcolours[1L], fgFill = spcolours[1L],
                                                      numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Global Spec. counts"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = spcolours[2L], fgFill = spcolours[2L],
                                                         numFmt = "COMMA", fontSize = 11L, valign = "center")
Styles[["Spec. IDs"]] <- openxlsx::createStyle(bgFill = spcolours[3L], fgFill = spcolours[3L],
                                               numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Spec. counts"]] <- openxlsx::createStyle(bgFill = spcolours[4L], fgFill = spcolours[4L],
                                                  numFmt = "COMMA", fontSize = 11L, valign = "center")
Styles[["Global Biot. Pep. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = pepcolours[5L], fgFill = pepcolours[1L],
                                                           numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Global Biot. Pep. counts"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = pepcolours[6L], fgFill = pepcolours[2L],
                                                              numFmt = "COMMA", fontSize = 11L, valign = "center")
Styles[["Biot. Pep. IDs"]] <- openxlsx::createStyle(bgFill = pepcolours[5L], fgFill = pepcolours[7L],
                                                    numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Biot. Pep. counts"]] <- openxlsx::createStyle(bgFill = pepcolours[6L], fgFill = pepcolours[8L],
                                                       numFmt = "COMMA", fontSize = 11L, valign = "center")
Styles[["Global Biot.Ev. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = evcolours[5L], fgFill = evcolours[1L],
                                                         numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Biot. Ev. IDs"]] <- openxlsx::createStyle(bgFill = evcolours[5L], fgFill = evcolours[7L],
                                                   numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Biot. Ev. counts"]] <- openxlsx::createStyle(bgFill = evcolours[6L], fgFill = evcolours[8L],
                                                      numFmt = "COMMA", fontSize = 11L, valign = "center")
# (Placeholders, not used currently)
Styles[["Global Biot. Spec. IDs"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = spcolours[5L], fgFill = spcolours[1L],
                                                            numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Global Biot. Spec. counts"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = spcolours[6L], fgFill = spcolours[2L],
                                                               numFmt = "COMMA", fontSize = 11L, valign = "center")
Styles[["Biot. Spec. IDs"]] <- openxlsx::createStyle(bgFill = spcolours[5L], fgFill = spcolours[7L],
                                                     numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["Biot. Spec. counts"]] <- openxlsx::createStyle(bgFill = spcolours[6L], fgFill = spcolours[8L],
                                                        numFmt = "COMMA", fontSize = 11L, valign = "center")
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
                                                 numFmt = "TEXT", fontSize = 11L, valign = "center")
# - Regulated
Styles[["Regulated"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = "lightgoldenrod1", fgFill = "lightgoldenrod1",
                                               numFmt = "TEXT", fontSize = 11L, valign = "center")
# F-test
Styles[["F-test summary Ratios"]] <- Styles[["Summary Ratios"]]
Styles[["F-test P-values"]] <- Styles[["P-values"]]
Styles[["F-test significant"]] <- Styles[["Significant"]]
Styles[["F-test regulated"]] <- Styles[["Regulated"]]
# Localisation
Styles[["Localisation"]] <- openxlsx::createStyle(textDecoration = "bold", bgFill = "azure", fgFill = "azure",
                                                  numFmt = "TEXT", fontSize = 11L, valign = "center")
Styles[["SSDs"]] <- "SSDs"
Styles[["Mean SSDs"]] <- "Mean SSDs"
Styles[["SSDs P-values"]] <- Styles[["P-values"]]
Styles[["SSDs significant"]] <- Styles[["Significant"]]
Styles[["Re-localized"]] <- Styles[["Regulated"]]
Styles[["Marker"]] <- openxlsx::createStyle(textDecoration = "bold", numFmt = "TEXT", fontSize = 11L, valign = "center")
# - Annotations
annot %<o% c("InterPro", "Pfam", "PIRSF", "PROSITE")
AnnotTbl %<o% data.frame(Name = c("GO", "Taxonomy", annot, "EMBL", "Other"),
                         Colour = colorRampPalette(c("cyan", "cyan3"))(8L), check.names = FALSE)
for (i in 1L:nrow(AnnotTbl)) {
  Styles[[paste0(AnnotTbl$Name[i], " annotations")]] <- openxlsx::createStyle(bgFill = AnnotTbl$Colour[i], fgFill = AnnotTbl$Colour[i],
                                                                              numFmt = "TEXT", fontSize = 11L, valign = "center")
}
# - PEP
Styles[["PEP"]] <- "PEP"
# - Filters (contaminants, only identified by site)
Styles[["Filters"]] <- openxlsx::createStyle(textDecoration = "italic", halign = "center", valign = "center", wrapText = TRUE,
                                             numFmt = "TEXT", fontSize = 12L)
# - Clusters
Styles[["Cluster"]] <- "Cluster"
# - Header
Styles[["Header 1"]] <- openxlsx::createStyle(fontColour = "white", bgFill = "#5b9bd5", fgFill = "#5b9bd5",
                                              border = "left", borderStyle = "thin", borderColour = "dimgrey",
                                              textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE,
                                              numFmt = "TEXT", fontSize = 12L)
Styles[["Header 2"]] <- openxlsx::createStyle(fontColour = "white", bgFill = "#5b9bd5", fgFill = "#5b9bd5",
                                              border = "left", borderStyle = "thin", borderColour = "dimgrey",
                                              textDecoration = c("bold", "underline"), halign = "center", valign = "center", wrapText = TRUE,
                                              numFmt = "TEXT", fontSize = 12L)
Styles[["LeftBorder"]] <- openxlsx::createStyle(border = "left", borderStyle = "thick")
Styles[["RightBorderThin"]] <- openxlsx::createStyle(border = "right", borderStyle = "thin", borderColour = "dimgrey")
Styles[["UpperBorder"]] <- openxlsx::createStyle(border = "top", borderStyle = "thick")
Styles[["NoUpperBorder"]] <- openxlsx::createStyle(border = "top", borderStyle = "none")
# - Proteins in list
Styles[["Protein list"]] <- openxlsx::createStyle(textDecoration = "bold", fontColour = "brown")
# - Negative filter
Styles[["Negative filter"]] <- "Negative filter"
# - Summary tab
Styles[["Interlines"]] <- openxlsx::createStyle(fontSize = 12L, halign = "left", valign = "center", textDecoration = "underline2", wrapText = FALSE)
Styles[["Small table - TEXT"]] <- openxlsx::createStyle(fontSize = 10L, halign = "center", valign = "center", wrapText = FALSE, numFmt = "TEXT")
Styles[["Small table - NUMBER"]] <- openxlsx::createStyle(fontSize = 10L, halign = "center", valign = "center", wrapText = FALSE, numFmt = "NUMBER")
Styles[["Header 3"]] <- openxlsx::createStyle(bgFill = "darkgrey", fgFill = "blue",
                                              fontSize = 11L, valign = "center", halign = "center", wrapText = TRUE,
                                              border = "TopBottomLeftRight", borderColour = "black", borderStyle = "thick", numFmt = "TEXT")
