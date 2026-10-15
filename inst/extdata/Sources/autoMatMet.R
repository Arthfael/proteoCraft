# Start writing Materials and Methods
if (ProcessedByUs) { cat("Materials & Methods: (the user will be prompted for information...)\n") }
#
# 1) Wet lab
if (scrptType == "noReps") { nr <- nrow(SamplesMap) }
if (scrptType == "withReps") { nr <- nrow(Exp.map) }
MatMetCalls %<o% list(Calls = list("read_docx()",
                                   paste0("body_add_fpar(MatMet, fpar(ftext(\"Sample", c("", "s")[(nr>1L)+1L],
                                          " preparation\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")),
                      Texts = list())
if (ProcessedByUs) {
  if (scrptType == "noReps") { WetLabMeth <- try(MatMet_WetLab(exp.map = SamplesMap), silent = TRUE) }
  if (scrptType == "withReps") {
    if (LabelType == "LFQ") { WetLabMeth <- try(MatMet_WetLab(), silent = TRUE) }
    if (LabelType == "Isobaric") { WetLabMeth <- try(MatMet_WetLab(Label = IsobarLab), silent = TRUE) }
  }
  if (inherits(WetLabMeth, "try-error")) { WetLabMeth <- "TEMPLATE" }
} else { WetLabMeth <- "TEMPLATE" }
MatMetCalls$Texts$WetLab <- WetLabMeth
for (i in 1L:length(WetLabMeth)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls,
                              paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$WetLab[", i,"], prop = WrdFrmt$",
                                     c("Body", "Template_text")[(WetLabMeth[i] == "TEMPLATE")+1L], "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#
# 2) LCMS
if ((!exists("LCMS_instr")) || (!inherits(LCMS_instr, "list")) || (sum(c("LC", "MS") %in% names(LCMS_instr)) != 2L)) {
  LCMS_instr <- list(LC = c(),
                     MS = c())
}
MatMetCalls$Calls <- append(MatMetCalls$Calls,
                            "body_add_fpar(MatMet, fpar(ftext(\"LC-MS/MS analysis\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
LCMS_meth_lst <- try(MatMet_LCMS(cl = parClust), silent = TRUE)
mzMLtst <- ("mzML" %in% gsub(".*\\.", "", rawFiles))
if (mzMLtst) { # Fix for when we searched mzML-converted files
  if (inherits(LCMS_meth_lst, "try-error") || is.null(LCMS_meth_lst)) {
    tmp <- gsub("\\.mzML", ".raw", rawFiles)
    LCMS_meth_lst <- try(MatMet_LCMS(RawFiles = tmp, cl = parClust), silent = TRUE)
    if (inherits(LCMS_meth_lst, "try-error") || is.null(LCMS_meth_lst)) {
      tmp <- gsub("\\.mzML", ".d", rawFiles)
      LCMS_meth_lst <- try(MatMet_LCMS(RawFiles = tmp, cl = parClust), silent = TRUE)
    }
  }
}
if ((inherits(LCMS_meth_lst, "try-error")) || is.null(LCMS_meth_lst)) {
  LCMS_meth <- "TEMPLATE"
} else {
  LCMS_meth <- LCMS_meth_lst$Text
  LCMS_instr <- LCMS_meth_lst$Instruments
}
LCMS_meth %<o% LCMS_meth
LCMS_instr %<o% LCMS_instr
MatMetCalls$Texts$LCMS <- LCMS_meth
for (i in 1L:length(LCMS_meth)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls,
                              paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$LCMS[", i,
                                     "], prop = WrdFrmt$",
                                     c("Body", "Template_text")[(LCMS_meth[i] == "TEMPLATE")+1L],
                                     "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#
# 3) Search
MatMetCalls$Texts$DatAnalysis <- MatMet_Search
#
# 4) Post-processing
if ((!exists("DatAnalysisTxt")) || (!inherits(DatAnalysisTxt, "list"))) { DatAnalysisTxt %<o% list() }
DatAnalysisTxt$Start <- if (length(SearchSoft) == 1L) {
    paste0(names(SearchSoft), "'s output was re-processed using in-house R scripts, starting from the ",
           c("evidence.txt", "main report", "psm.tsv")[match(SearchSoft, SearchSoftware)],
           " table", c("", "s")[(length(PSMsFls)>1L)+1L], ".")
  } else {
    "The output PSMs tables from individual search engines were converted to similar formats and combined, then re-processed using in-house R scripts."
  }
MatMetCalls$Calls <- append(MatMetCalls$Calls,
                            "body_add_fpar(MatMet, fpar(ftext(\"Data analysis\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
