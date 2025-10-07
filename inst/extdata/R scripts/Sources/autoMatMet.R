# Start writing Materials and Methods
cat("Writing Materials & Methods template...\n")
#
# 1) Wet lab
if (scrptType == "noReps") { nr <- nrow(SamplesMap) }
if (scrptType == "withReps") { nr <- nrow(Exp.map) }
MatMetCalls %<o% list(Calls = list("read_docx()",
                                   paste0("body_add_fpar(MatMet, fpar(ftext(\"Sample", c("", "s")[(nr > 1)+1],
                                          " preparation\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")),
                      Texts = list())
if (ProcessedByUs) {
  if (scrptType == "noReps") { WetLabMeth <- try(MatMet_WetLab(exp.map = SamplesMap), silent = TRUE) }
  if (scrptType == "withReps") {
    if (LabelType == "LFQ") { WetLabMeth <- try(MatMet_WetLab(), silent = TRUE) }
    if (LabelType == "Isobaric") { WetLabMeth <- try(MatMet_WetLab(Label = IsobarLab), silent = TRUE) }
  }
  if ("try-error" %in% class(WetLabMeth)) { WetLabMeth <- "TEMPLATE" }
} else { WetLabMeth <- "TEMPLATE" }
MatMetCalls$Texts$WetLab <- WetLabMeth
for (i in 1:length(WetLabMeth)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls,
                              paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$WetLab[", i,"], prop = WrdFrmt$",
                                     c("Body", "Template_text")[(WetLabMeth[i] == "TEMPLATE")+1], "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#
# 2) LCMS
if ((!exists("LCMS_instr"))||(!"list" %in% class(LCMS_instr))||(sum(c("LC", "MS") %in% names(LCMS_instr)) != 2)) {
  LCMS_instr <- list(LC = c(),
                     MS = c())
}
MatMetCalls$Calls <- append(MatMetCalls$Calls,
                            "body_add_fpar(MatMet, fpar(ftext(\"LC-MS/MS analysis\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
LCMS_meth_lst <- try(MatMet_LCMS(cl = parClust), silent = TRUE)
mzMLtst <- ("mzML" %in% gsub(".*\\.", "", rawFiles))
if (mzMLtst) { # Fix for when we searched mzML-converted files
  if (("try-error" %in% class(LCMS_meth_lst))||(is.null(LCMS_meth_lst))) {
    tmp <- gsub("\\.mzML", ".raw", rawFiles)
    LCMS_meth_lst <- try(MatMet_LCMS(RawFiles = tmp, cl = parClust), silent = TRUE)
    if (("try-error" %in% class(LCMS_meth_lst))||(is.null(LCMS_meth_lst))) {
      tmp <- gsub("\\.mzML", ".d", rawFiles)
      LCMS_meth_lst <- try(MatMet_LCMS(RawFiles = tmp, cl = parClust), silent = TRUE)
    }
  }
}
if ((("try-error" %in% class(LCMS_meth_lst)))||(is.null(LCMS_meth_lst))) {
  LCMS_meth <- "TEMPLATE"
} else {
  LCMS_meth <- LCMS_meth_lst$Text
  LCMS_instr <- LCMS_meth_lst$Instruments
}
LCMS_meth %<o% LCMS_meth
LCMS_instr %<o% LCMS_instr
MatMetCalls$Texts$LCMS <- LCMS_meth
for (i in 1:length(LCMS_meth)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls,
                              paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$LCMS[", i,
                                     "], prop = WrdFrmt$",
                                     c("Body", "Template_text")[(LCMS_meth[i] == "TEMPLATE")+1],
                                     "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#
# 3) Search
MatMetCalls$Texts$DatAnalysis <- MatMet_Search
#
# 4) Post-processing
DatAnalysisTxt %<o% "TEMPLATE"
if (length(SearchSoft) == 1) {
  DatAnalysisTxt <- paste0(names(SearchSoft), "'s output was re-processed using in-house R scripts, starting from the ",
                           c("evidence.txt", "main report", "psm.tsv")[match(SearchSoft, SearchSoftware)],
                           " table", c("", "s")[(length(PSMsFls)>1)+1], ".")
} else {
  DatAnalysisTxt <- "The output PSMs tables from the individual engines were converted to similar formats and combined, then re-processed using in-house R scripts."
}
# For the scripts with reps we will expand this later, during subsequent parts of the script...
if (scrptType == "noReps") {
  # ... but for the no-reps script this is done in one block here:
  DatAnalysisTxt <- paste0(DatAnalysisTxt,
                           c("", " Peptide-to-protein assignments were checked, then")[Update_Prot_matches+1],
                           paste0(" Protein Groups were assembled and quantified using an algorithm which: ",
                                  "i) computes a mean protein group-level profile across samples using individual, normalized peptidoform profiles (\"relative quantitation\" step), then ",
                                  "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step).",
                                  c(" Only unique peptidoforms were used.",
                                    " Only unique and razor peptidoforms were used.",
                                    "")[match(Pep4Quant, Pep4QuantOpt)], collapse = " "))
  if (nrow(Mod2Xclud)) {
    tmp <- Modifs$`Full name`[match(Mod2Xclud$Mark, Modifs$Mark)]
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " ", tmp, "-peptidoforms and their unmodified counterparts were excluded from the calculations.")
  }
  DatAnalysisTxt <- paste0(DatAnalysisTxt,
                           c("", " Protein group-level quantitative values were normalized using the Levenberg-Marquardt procedure.")[NormalizePG+1])
  MatMetCalls$Texts$DatAnalysis <- c(MatMetCalls$Texts$DatAnalysis, DatAnalysisTxt)
  L <- length(MatMetCalls$Texts$DatAnalysis)
  for (i in 1:L) {
    MatMetCalls$Calls <- append(MatMetCalls$Calls,
                                paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$DatAnalysis[", i,
                                       "], prop = WrdFrmt$",
                                       c("Body", "Template_text")[(MatMetCalls$Texts$DatAnalysis[i] == "TEMPLATE")+1],
                                       "_text), fp_p = WrdFrmt$just))"))
  }
  MatMetCalls$Calls <- append(MatMetCalls$Calls,
                              "body_add_par(MatMet, \"\", style = \"Normal\")")
} # For reps, this is done later
