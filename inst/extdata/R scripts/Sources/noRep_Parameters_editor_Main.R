#### Define analysis parameters
require(shiny)
require(shinyFiles)
require(shinycssloaders)
require(shinyjs)
require(shinyBS)
require(htmlwidgets)
#
Exp %<o% unique(SamplesMap$Experiment)
nSmpls <- length(Exp)
moreThan1Sample %<o% (length(Exp) > 1)
#
# Boolean functions to check parameter values
Src <- paste0(libPath, "/extdata/R scripts/Sources/parBooleans.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# Species
SpeciesTst %<o% "Unspecified"
if ("Taxonomy" %in% colnames(db)) {
  SpeciesTst <- unique(db$Taxonomy[which(gsub(" *(\\(|\\[).*", "", db[[dbOrgKol]]) == mainOrg)])
  SpeciesTst <- SpeciesTst[which(as.character(SpeciesTst) != "NA")][1]
}
if ("Kingdom" %in% colnames(db)) {
  KingdomTst <- aggregate(db$Kingdom, list(db$Kingdom), length)
  KingdomTst <- KingdomTst[order(KingdomTst$x, decreasing = TRUE),]
  KingdomTst <- KingdomTst$Group.1[1]
} else {
  KingdomTst <- "unknown" # Could be user prompted
}
KingdomTst %<o% KingdomTst
isEukaLike %<o% (KingdomTst %in% c("Eukaryota", "Archaea"))
#
if (file.exists("AnalysisParam.RData")) {
  tmp <- AnalysisParam
  loadFun("AnalysisParam.RData")
  for (nm in names(AnalysisParam)) {
    if (!nm %in% names(tmp)) { tmp[[nm]] <- AnalysisParam[[nm]] }
  }
  AnalysisParam <- tmp
}
#
ptmDflt1 <- grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE, invert = TRUE)
Mod4Quant %<o% Modifs$Mark[match(ptmDflt1, Modifs$`Full name`)]
if ("PTMs eligible for quantitation" %in% names(AnalysisParam)) {
  Mod4Quant <- AnalysisParam$"PTMs eligible for quantitation"
}
ptmDflt1 <- Modifs$`Full name`[match(Mod4Quant, Modifs$Mark)]
Mod2Xclud %<o% set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                            c("Mark", "Where"))
if ("PTMs to exclude from quantitation" %in% names(AnalysisParam)) {
  Mod2Xclud <- AnalysisParam$"PTMs to exclude from quantitation"
}
ptmDflt2 <- paste(grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE), collapse = ";")
#
klustChoices %<o% c("K-means", "hierarchical")
KlustMeth %<o% 1 # Changed from 2
#
if (Annotate) {
  allGO <- unique(unlist(strsplit(db$GO[which(!is.na(db$GO))], ";")))
  allGO2 <- paste0("GO:", gsub(".* \\[GO:|\\]$", "", allGO))
}
if (!"GO.terms.for.proteins.of.interest" %in% names(AnalysisParam)) { AnalysisParam$GO.terms.for.proteins.of.interest <- FALSE }
if (!"Custom.PGs" %in% names(AnalysisParam)) { AnalysisParam$Custom.PGs <- "" }
if (!"CRAPome_file" %in% names(AnalysisParam)) { AnalysisParam$CRAPome_file <- "" }
#
pr <- c("Norma.Ev.Intens", "Norma.Pep.Intens", "Adv.Norma.Pep.Intens", "Norma.Prot.Ratio",
        "Adv.Norma.Prot.Intens", "Adv.Norma.Ev.Intens",
        "Norma.Pep.Ratio", "Adv.Norma.Pep.Ratio", "Norma.Prot.Ratio.to.Biot")
prDflt <- setNames(rep(TRUE, length(pr)), pr)
prDflt["Adv.Norma.Ev.Intens"] <- (length(unique(FracMap$Fraction)) > 1)|
  (length(unique(FracMap$`PTM-enriched`)) > 1)
prDflt[c("Norma.Pep.Ratio", "Adv.Norma.Pep.Ratio", "Norma.Prot.Ratio.to.Biot")] <- FALSE
for (p in pr) { #p <- pr[1]
  dflt <- prDflt[p]
  if (!p %in% names(AnalysisParam)) { AnalysisParam[[p]] <- dflt } else {
    tmp <- AnalysisParam[[p]]
    if (!validLogicPar("tmp")) { AnalysisParam[[p]] <- dflt }
  }
}
#
# Peptide classes to use for quantitation
Pep4QuantOpt %<o% setNames(c("Unique peptide IDs", "Razor peptide IDs", "Peptide IDs"),
                           c("Unique", "Razor", "All"))
if ((!validCharPar("Pep4Quant", Pep4QuantOpt))&&("Peptide classes eligible for quantitation" %in% names(AnalysisParam))) {
  tmp1 <- toupper(gsub(" |_|-|\\.", "", as.character(AnalysisParam[["Peptide classes eligible for quantitation"]])))
  tmp2 <- Pep4QuantOpt
  names(tmp2) <- toupper(names(tmp2))
  m12 <- match(tmp1, names(tmp2))
  if (!is.na(m12)) {
    Pep4Quant <- Pep4QuantOpt[m12]
  }
}
if (!validCharPar("Pep4Quant", Pep4QuantOpt)) {
  Pep4Quant <- Pep4QuantOpt[c("Unique", "Razor")[isEukaLike+1]]
}
Pep4Quant %<o% Pep4Quant
AnalysisParam$Prot.Quant.Use <- names(Pep4Quant)
#
# How many peptide do we need to quantify a protein?
if ((!validIntegPar("N_Pep"))&&("N. of peptidoforms for quantitation" %in% names(AnalysisParam))) {
  tmp1 <- AnalysisParam$"N. of peptidoforms for quantitation"
  if (validIntegPar("tmp1")) { N_Pep <- as.integer(tmp1) }
}
if (!validIntegPar("N_Pep")) {
  # The default is to use ALL of the information, even if we only have a single peptide
  # We will flag proteins flagged with 1 peptide as dodgy
  N_Pep <- 1
}
N_Pep %<o% N_Pep
AnalysisParam$"N. of peptidoforms for quantitation" <- N_Pep
#
# How many unique peptides (if available) to use to the exclusivity of any others
# If set to 0, we will use peptides regardless of whether unique or not
if ((!validIntegPar("N_unique_Pep", 0))&&("Use.N.unique" %in% names(AnalysisParam))) {
  tmp1 <- AnalysisParam$Use.N.unique
  if (validIntegPar("tmp1", 0)) { N_unique_Pep <- as.integer(tmp1) }
}
if (!validIntegPar("N_unique_Pep", 0)) {
  N_unique_Pep <- 3 # ... which means if we have 3 unique peptides we will not use any non-unique ones
}
N_unique_Pep %<o% N_unique_Pep
AnalysisParam$Use.N.unique <- N_unique_Pep
#
ptmDflt1 <- grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE, invert = TRUE)
Mod4Quant %<o% Modifs$Mark[match(ptmDflt1, Modifs$`Full name`)]
if ("Prot.Quant.Mod.Excl" %in% names(AnalysisParam)) {
  Mod4Quant <- Mod4Quant[which(!Mod4Quant %in% unlist(strsplit(AnalysisParam$Prot.Quant.Mod.Excl, ";")))]
}
ptmDflt1 <- Modifs$`Full name`[match(Mod4Quant, Modifs$Mark)]
Mod2Xclud %<o% set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                            c("Mark", "Where"))
#
allQuantAlgos %<o% data.frame(Algorithm = c("limpa",
                                            "LM",
                                            "MaxLFQ (iq)",
                                            "QFeatures"),
                              Details = c("limpa: https://www.biorxiv.org/content/10.1101/2025.04.28.651125v1",
                                          "in-house MaxLFQ-like algorithm",
                                          "MaxLFQ as implemented in iq: https://pubmed.ncbi.nlm.nih.gov/31909781/",
                                          "QFeatures: Gatto L, Vanderaa C (2025). QFeatures: Quantitative features for mass spectrometry data. doi:10.18129/B9.bioc.QFeatures"),
                              Help = c("limpa uses a probabilistic approach to model protein group abundance, is meant to be used upstream of limma and requires at least two samples. When this method is chosen, Imputation and Protein Group-level re-normalization are turned off by default, and the type of t-test is set to \"Moderated\".",
                                       "An in-house MaxLFQ-like algorithm using the Levenberg-Marquardt procedure to align peptide profiles prior to averaging.",
                                       "The fast version of MaxLFQ as implemented in package iq.",
                                       "QFeatures aggregateFeatures() uses a robust summarization procedure and is meant to be used upstream of msqrob2 statistics."))
quantAlgoOpt %<o% allQuantAlgos$Algorithm
if (scrptType == "noReps") { # limpa needs at least 2 samples
  quantAlgoOpt <- quantAlgoOpt[which(!quantAlgoOpt == "limpa")]
}
if (("QuantMeth" %in% names(AnalysisParam))&&(!"Quant_algorithm" %in% names(AnalysisParam))) { # Old AnalysisParameter name
  AnalysisParam$Quant_algorithm <- AnalysisParam$QuantMeth
}
if ((!validCharPar("quantAlgo", quantAlgoOpt))&&("Quant_algorithm" %in% names(AnalysisParam))) {
  tmp1 <- AnalysisParam$Quant_algorithm
  if (validCharPar("tmp1", quantAlgoOpt)) { quantAlgo <- tmp1 }
}
if (!validCharPar("quantAlgo", quantAlgoOpt)) {
  quantAlgo <- c("LM", "limpa")[match(scrptType, c("noReps", "withReps"))]
}
quantAlgo %<o% quantAlgo
AnalysisParam$Quant_algorithm <- quantAlgo
quantMsg <- function(Quant) { allQuantAlgos$Details[match(Quant, allQuantAlgos$Algorithm)] }
# Re-scaling
allReScAlgoOpt %<o% data.frame(Algorithm = c("limpa",
                                             "QFeatures",
                                             "median",
                                             "topN",
                                             "weighted.mean",
                                             "max",
                                             "sum",
                                             "MaxLFQ (iq)"),
                               Details = c("Use limpa's probabilistic model-based scale",
                                           "Use QFeatures' scale, based on robustSummary(peptide intensities)",
                                           "Re-scale to the median of quantitative peptide intensities",
                                           "Re-scale to the mean of up to N quantitative peptide intensities",
                                           "Re-scale to a weighted mean of peptide intensities",
                                           "Re-scale to the most intense peptide",
                                           "Re-scale to the sum of all quantitative peptide intensities (makes no sense, but you can do it anyway)",
                                           "Use MaxLFQ's scale"))
reScAlgoOpt %<o% allReScAlgoOpt$Algorithm
if (scrptType == "noReps") { # limpa needs at least 2 samples
  reScAlgoOpt <- reScAlgoOpt[which(!reScAlgoOpt == "limpa")]
}
if ((!validCharPar("reScAlgo", reScAlgoOpt))&&("ReScaling_algorithm" %in% names(AnalysisParam))) {
  tmp1 <- AnalysisParam$ReScaling_algorithm
  if (validCharPar("tmp1", reScAlgoOpt)) { reScAlgo <- tmp1 }
}
if (!validCharPar("reScAlgo", reScAlgoOpt)) {
  if (quantAlgo %in% reScAlgoOpt) { reScAlgo <- quantAlgo }
  if (quantAlgo == "LM") { reScAlgo <- "max" }
}
reScAlgo %<o% reScAlgo
AnalysisParam$ReScaling_algorithm <- reScAlgo
reScMsg <- function(Rescale, Quant) {
  if (Rescale == Quant) {
    return("Use scale from quantitation algorithm")
  }
  if (Rescale %in% allReScAlgoOpt$Algorithm) {
    return(allReScAlgoOpt$Details[match(Rescale, allReScAlgoOpt$Algorithm)])
  }
  return("???... WHAT?!")
}
if (!validLogicPar("reScale")) {
  reScale <- (quantAlgo != reScAlgo)|(quantAlgo == "LM")
}
# Top-N
if ((!validIntegPar("topN"))&&("topN" %in% names(AnalysisParam))) {
  tmp1 <- AnalysisParam$topN
  if (validIntegPar("tmp1")) { topN <- as.integer(tmp1) }
}
if (!validIntegPar("topN")) { topN <- 3 } # Usually top3 
topN %<o% topN
AnalysisParam$topN <- topN
#
# Top-N correction: when averaging first/second/third peptides, should we 
# Let's say we have 2 proteins to re-scale:
# - for the first, we have 3 peptides... swell!
# - for the second, we have 1... shucks!
# Now, we want a value in all cases. So we'll also work with 1 peptide.
# But rank 3 peptides are systematically lower intensity than rank 2 peptides which are lower intensity than rank 1 peptide.
# So if we want to average 3-or-failing-that-2-or-failing-that-1-peptide(s),
# and still get a similar treatment for proteins with 1, 2 or 3+ peptides,
# we need to correct for that difference.
if ((!validLogicPar("topN_correct"))&&("topN_correct" %in% names(AnalysisParam))) {
  tmp1 <- as.logical(AnalysisParam$topN_correct)
  if (validLogicPar("tmp1")) { topN_correct <- tmp1 }
}
if (!validLogicPar("topN_correct")) { topN_correct <- TRUE }
topN_correct %<o% topN_correct
AnalysisParam$topN_correct <- topN_correct
#
# Check peptide-to-protein matches
if ((!validLogicPar("Update_Prot_matches"))&&("Update_Prot_matches" %in% names(AnalysisParam))) {
  tmp1 <- as.logical(AnalysisParam$Update_Prot_matches)
  if (validLogicPar("tmp1")) { Update_Prot_matches <- tmp1 }
}
if (!validLogicPar("Update_Prot_matches")) { Update_Prot_matches <- TRUE }
Update_Prot_matches %<o% Update_Prot_matches # See https://github.com/vdemichev/DiaNN/discussions/1631
AnalysisParam$Update_Prot_matches <- Update_Prot_matches
#
if ((!validLogicPar("Impute"))&&("ImputeMissData" %in% names(AnalysisParam))) {
  tmp1 <- as.logical(AnalysisParam$ImputeMissData)
  if (validLogicPar("tmp1")) { Impute <- tmp1 }
}
if (!validLogicPar("Impute")) { Impute <- FALSE }
Impute %<o% Impute
AnalysisParam$ImputeMissData <- Impute
#
if ((!validIntegPar("PepFoundInAtLeast", 1))&&("PepFoundInAtLeast" %in% names(AnalysisParam))) {
  tmp1 <- AnalysisParam$PepFoundInAtLeast
  if (validIntegPar("tmp1", 0)) { PepFoundInAtLeast <- as.integer(tmp1) }
}
if (!validIntegPar("PepFoundInAtLeast", 0)) {
  PepFoundInAtLeast <- 1 # ... which means if we have 3 unique peptides we will not use any non-unique ones
}
PepFoundInAtLeast %<o% PepFoundInAtLeast
AnalysisParam$PepFoundInAtLeast <- PepFoundInAtLeast
#
# Cytoscape
CytoScExe %<o% c()
cytoTst <- grep("cytoscape", list.dirs("C:/PROGRA~1", recursive = FALSE), value = TRUE, ignore.case = TRUE)
CytoScape %<o% (length(cytoTst) > 0)
if (CytoScape) {
  CytoScExe <- unlist(lapply(cytoTst, function(x) {
    grep("/Cytoscape\\.exe$", list.files(x, recursive = TRUE, full.names = TRUE), value = TRUE)
  }))
}
CytoScape %<o% (length(CytoScExe) > 0)
if (length(CytoScExe) > 1) {
  tst <- vapply(CytoScExe, function(x) { file.info(x)$mtime }, 1)
  CytoScExe <- CytoScExe[order(tst, decreasing = TRUE)]
}
if ("CytoScapePath" %in% names(AnalysisParam)) {
  AnalysisParam$CytoscapePath <- AnalysisParam$CytoScapePath
  AnalysisParam$CytoScapePath <- NULL
}
if ("CytoscapePath" %in% names(AnalysisParam)) {
  tmp <- normalizePath(AnalysisParam$CytoscapePath, winslash = "/")
  if ((length(CytoScExe))&&(tmp %in% CytoScExe)) { CytoScExe <- tmp }
}
if (CytoScape&&("Cytoscape" %in% names(AnalysisParam))) {
  tmp1 <- as.logical(AnalysisParam$Cytoscape)
  if (validLogicPar("tmp1")&&(!CytoScape)) {
    msg <- "Sorry, can't run Cytoscape: couldn't find a valid executable!"
    ReportCalls <- AddMsg2Report(Warning = TRUE)
  } else { CytoScape <- tmp1 }
}
AnalysisParam$Cytoscape <- CytoScape
#
# Proteome ruler
if ((!validLogicPar("protrul"))&&("Proteomic ruler calculated" %in% names(AnalysisParam))) {
  tmp1 <- as.logical(AnalysisParam$"Proteomic ruler calculated")
  if (validLogicPar("tmp1")) { protrul <- tmp1 }
}
if (!validLogicPar("protrul")) {
  # Proteome ruler re-scaling makes no sense if we are dealing with a pull down
  protrul %<o% (!WorkFlow %in% c("PULLDOWN", "BIOID"))
  # Archaea and Eukaryotes have introns and histones, Bacteria do not
  protrul <- c(protrul, FALSE)[(!isEukaLike)+1]
}
protrul %<o% protrul
AnalysisParam$ProtRul <- protrul
if ((!validIntegPar("ProtRulNuclL", 1))&&("ProtRulNuclL" %in% names(AnalysisParam))) {
  tmp1 <- AnalysisParam$ProtRulNuclL
  if (validIntegPar("tmp1")) { ProtRulNuclL <- as.integer(tmp1) }
}
if (!validIntegPar("ProtRulNuclL")) { ProtRulNuclL <- 196 } # Usually top3 
ProtRulNuclL %<o% ProtRulNuclL
AnalysisParam$ProtRulNuclL <- ProtRulNuclL
#
Update_Prot_matches %<o% TRUE # See https://github.com/vdemichev/DiaNN/discussions/1631
if ("Update_Prot_matches" %in% names(AnalysisParam)) {
  Update_Prot_matches <- as.logical(AnalysisParam$Update_Prot_matches)
  if ((is.na(Update_Prot_matches))||(is.null(Update_Prot_matches))) { Update_Prot_matches <- TRUE }
}
if ((!is.logical(Update_Prot_matches))||(length(Update_Prot_matches) != 1)||(is.na(Update_Prot_matches))) {
  Update_Prot_matches <- TRUE
}
AnalysisParam$Update_Prot_matches <- Update_Prot_matches
#
PepFoundInAtLeast %<o% 1
if ("PepFoundInAtLeast" %in% names(AnalysisParam)) {
  PepFoundInAtLeast <- suppressWarnings(as.integer(AnalysisParam$PepFoundInAtLeast))
  if ((is.na(PepFoundInAtLeast))||(PepFoundInAtLeast < 1)||(!moreThan1Sample)||(PepFoundInAtLeast > length(Exp))) {
    PepFoundInAtLeast <- 1
  }
}
AnalysisParam$PepFoundInAtLeast <- PepFoundInAtLeast
#
dbOrd <- 1:nrow(db)
protDflt <- NULL
tmp <- c()
if ("Protein of interest" %in% names(db)) {
  tmp <- unique(c(tmp, db$`Protein ID`[which(db$"Protein of interest")]))
}
if ("Prot.list" %in% names(AnalysisParam)) {
  tmp <- unique(c(tmp, unlist(strsplit(AnalysisParam$Prot.list, ";"))))
}
if ("Prot.list_pep" %in% names(AnalysisParam)) {
  tmp <- unique(c(tmp, unlist(strsplit(AnalysisParam$Prot.list_pep, ";"))))
}
if (length(tmp)) {
  m <- match(tmp, db$`Protein ID`)
  m <- m[which(!is.na(m))]
  dbOrd <- c(m, which(!db$`Protein ID` %in% tmp))
}
protHeads <- gsub("^>", "", db$Header[dbOrd])
if (length(tmp)) { protDflt <- protHeads[1:length(m)] }
#
tstAdvOpt <- try(sum(file.exists(AnalysisParam$Custom.PGs, AnalysisParam$CRAPome_file)) > 0)
if ("try-error" %in% class(tstAdvOpt)) { tstAdvOpt <- FALSE }
if (MakeRatios) {
  RatiosThresh %<o% 1
  if ("Ratios analysis - threshold" %in% names(AnalysisParam)) {
    RatiosThresh <- AnalysisParam$"Ratios analysis - threshold"
  } else {
    AnalysisParam$"Ratios analysis - threshold" <-  RatiosThresh
  }
  if ((is.na(RatiosThresh))||(!is.numeric(RatiosThresh))||(!RatiosThresh < 0)) {
    RatiosThresh <- 1
  }
  RatiosThresh_2sided %<o% c(TRUE, FALSE)[(WorkFlow == "Pull-down")+1]
  if ("Ratios analysis - threshold is two-sided" %in% names(AnalysisParam)) {
    RatiosThresh_2sided <- AnalysisParam$"Ratios analysis - threshold is two-sided"
  } else {
    AnalysisParam$"Ratios analysis - threshold is two-sided" <- RatiosThresh_2sided
  }
  RatiosThresh_2sided <- as.logical(RatiosThresh_2sided)
  if ((is.na(RatiosThresh_2sided))||(is.null(RatiosThresh_2sided))) {
    RatiosThresh_2sided <- TRUE
  }
}
#
removeMBR %<o% FALSE
if ("Proteins list: remove match-between-runs" %in% names(AnalysisParam)) {
  removeMBR <- as.logical(AnalysisParam$"Proteins list: remove match-between-runs")
  if ((is.na(removeMBR))||(is.null(removeMBR))||(!moreThan1Sample)) { removeMBR <- FALSE }
}
AnalysisParam$"Proteins list: remove match-between-runs" <- removeMBR
#
NormalizePG %<o% moreThan1Sample
if ("NormalizePG" %in% names(AnalysisParam)) {
  NormalizePG <- as.logical(AnalysisParam$NormalizePG)
  if ((is.na(NormalizePG))||(is.null(NormalizePG))||(!moreThan1Sample)) { NormalizePG <- moreThan1Sample }
}
AnalysisParam$NormalizePG <- NormalizePG
#
# GO enrichment analysis
# - vs total proteome
if ("globalGO" %in% names(AnalysisParam)) {
  tmp <- as.logical(AnalysisParam$globalGO)
  if (validLogicPar("tmp")) {
    globalGO <- tmp
  }
}
if (!validLogicPar("globalGO")) {
  globalGO <- Annotate & (WorkFlow %in% c("Discovery", "Regulation", "Pull-down"))
}
if (is.na(globalGO)) { globalGO <- FALSE }
globalGO %<o% globalGO
# - vs reference sample
enrichGO %<o% (globalGO&MakeRatios)
if (Annotate) {
  if (!exists("GO_filter")) { GO_filter %<o% c() }
  if ("GO terms of interest" %in% names(AnalysisParam)) {
    GO_filter <- AnalysisParam$"GO terms of interest"
  } else { AnalysisParam$"GO terms of interest" <- GO_filter }
  GO_filter1 <- allGO[match(GO_filter, allGO2)]
}
#
Mod2Write %<o% c()
if ("Mod2Write" %in% names(AnalysisParam)) {
  Mod2Write <- AnalysisParam$Mod2Write
} else { AnalysisParam$Mod2Write <- Mod2Write }
if ((length(Mod2Write))&&(!Mod2Write %in% Modifs$Mark)) { Mod2Write <- c() }
#
# Venn diagrams
Venn_Obs %<o% moreThan1Sample
if ("Venn diagrams: observed" %in% names(AnalysisParam)) {
  Venn_Obs <- AnalysisParam$"Venn diagrams: observed"
}
Venn_Obs <- moreThan1Sample # No matter the parameters, if there is only one sample, we DO NOT draw Venn diagrams, period
Venn_Obs <- as.logical(Venn_Obs)
if ((is.na(Venn_Obs))||(is.null(Venn_Obs))) { Venn_Obs <- FALSE }
AnalysisParam$"Venn diagrams: observed" <- Venn_Obs
Venn_Obs %<o% Venn_Obs
#
# Some default parameters which nicely follow the same structure
myPar <- c("Pepper", "GSEA", "ProfPlots", "RankAbundPlots", "ClueGO")
for (parI in myPar) {
  parNm <- paste0("run", parI)
  # Lowest level default: defined by context
  par_dflt <- par_dflt2 <- c(FALSE,
                             WorkFlow %in% c("Discovery", "Regulation", "Pull-down"),
                             moreThan1Exp,
                             TRUE,
                             FALSE)[match(parI, myPar)]
  # Level 1 default: defined by AnalysisParam
  if (parNm %in% names(AnalysisParam)) {
    par_dflt <- as.logical(AnalysisParam[[parNm]])
  }
  # Level 2 default: defined by existing value
  parOK <- (exists(parNm))
  if (parOK) {
    tmpPar <- get(parNm)
    parOK <- (length(tmpPar) == 1)&&(is.logical(tmpPar))&&(!is.na(tmpPar))
  }
  if (parOK) { par_dflt <- tmpPar }
  # Backup value
  if (is.na(par_dflt)) { par_dflt <- par_dflt2 }
  #
  if (!parOK) { assign(parNm, par_dflt) }
  assign(parNm, par_dflt)
  .obj <- unique(c(parNm, .obj))
}
#
appNm <- paste0(dtstNm, " - Parameters")
make_ui <- function() {
  fluidPage(
    useShinyjs(),
    setBackgroundColor( # Doesn't work
      color = c(#"#F8F8FF",
        "#E6F7F4"),
      gradient = "linear",
      direction = "bottom"
    ),
    tags$head(
      tags$style(
        '.inner.open {
            overflow-y: hidden !important;
        }'
      )
    ),
    extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
    titlePanel(tag("u", "Parameters"),
               appNm),
    h2(dtstNm), 
    br(),
    tags$hr(style = "border-color: black;"),
    h4("Proteins of interest"),
    pickerInput("IntProt", NULL, protHeads, protDflt, TRUE,
                pickerOptions(title = "Search me",
                              `live-search` = TRUE,
                              actionsBox = TRUE,
                              deselectAllText = "Clear search")),
    br(),
    tags$hr(style = "border-color: black;"),
    h4("Data processing"),
    fluidRow(
      if (moreThan1Sample) {
        column(2,
               checkboxInput("Impute",
                             "Impute missing peptides-level values?",
                             Impute, "100%"))
      },
      column(2,
             checkboxInput("Update_Prot_matches",
                           paste0("Update ", names(SearchSoft),
                                  "'s original protein-to-peptides assignments?"),
                           Update_Prot_matches, "100%"),
             bsTooltip("Update_Prot_matches",
                       "Checking assignments may result in removal of some identifications. It is nonetheless recommended because we have observed occasional inconsistent peptides-to-protein assignments with some search software.",
                       placement = "right", trigger = "hover",
                       options = list(container = "body"))#,
             #withSpinner(uiOutput("ReloadMatches"))
      ),
      if (moreThan1Sample) {
        column(2, checkboxInput("prtNorm", "Normalize data",
                                NormalizePG, "100%")) 
      },
      column(2, checkboxInput("runPepper", "Run Pepper ML-based PSMs intensity correction?",
                              runPepper, "100%"))
    ),
    tags$hr(style = "border-color: black;"),
    br(),
    # Quantitation
    ## Choice of algorithm + Proteomics ruler
    tags$hr(style = "border-color: black;"),
    h4(strong("Protein Groups quantitation")),
    fluidRow(column(2,
                    h4("Peptides eligible for quantitation:"),
                    selectInput("Prot.Quant.Use",
                                "Peptide classes",
                                names(Pep4QuantOpt),
                                names(Pep4QuantOpt)[match(setNames(Pep4Quant, NULL),
                                                          setNames(Pep4QuantOpt, NULL))],
                                width = "100%"),
                    pickerInput("PTMsQuant",
                                "PTM(s)",
                                Modifs$"Full name", # Keep all mods here, otherwise we may mistakenly exclude carbamidomethyl peptides!!!
                                ptmDflt1,
                                TRUE,
                                pickerOptions(title = "Search me",
                                              `live-search` = TRUE,
                                              actionsBox = TRUE,
                                              deselectAllText = "Clear search",
                                              showTick = TRUE)),
                    numericInput("N_Pep",
                                 "Minimum number of peptides for quantitation", N_Pep, 1, Inf, 1, "100%"),
                    numericInput("N_unique_Pep",
                                 "Number of unique peptides which - if available - will be used to the exclusivity of non-unique ones",
                                 N_unique_Pep, 1, Inf, 1, "100%"),
                    strong(em("Use only peptidoforms found in at least how many samples in...")),
                    numericInput("PepFoundInAtLeast",
                                 " -> the whole dataset?", PepFoundInAtLeast, 1, nSmpls, 1, "100%")),
             column(3,
                    h4("Main LFQ algorithm"),
                    selectInput("Quant_algorithm",
                                "Quantitation algorithm",
                                quantAlgoOpt,
                                quantAlgo,
                                width = "100%"),
                    uiOutput("QuantMsg"),
                    br(),
                    h4("Re-Scaling"),
                    selectInput("ReScaling",
                                "Re-scaling algorithm",
                                reScAlgoOpt,
                                reScAlgo,
                                width = "100%"),
                    em("Re-scaling has no effect on individual intra-protein (row-wise) fold changes but defines inter-protein (column-wise) relative abundance estimates."),
                    br(),
                    br(),
                    strong(em("TopN parameters")),
                    fillRow(em(HTML("&nbsp;&nbsp;&nbsp;-> N = ")),
                            numericInput("topN", NULL, topN, 1, Inf, 1, "100%"),
                            height = "40px"),
                    fillRow(em(HTML("&nbsp;&nbsp;&nbsp;-> Correct systematic intensity shift between ranks?")),
                            checkboxInput("topN_correct", NULL, topN_correct, "100%"),
                            height = "40px")),
             column(2,
                    h4("Proteomic Ruler"),
                    checkboxInput("ProtRul",
                                  "Estimate copy numbers per cell using signal from all histones as reference?",
                                  protrul,
                                  "100%"),
                    numericInput("ProtRulNuclL",
                                 "Use inter-nucleosome length = ? (kb)",
                                 ProtRulNuclL,
                                 1,
                                 Inf,
                                 1,
                                 "100%"))
    ),
    br(),
    fluidRow(
      if (moreThan1Sample) {
        column(2,
               checkboxInput("removeMBR",
                             "Exclude match-between-runs peptides from coverage analysis?",
                             removeMBR, "100%")) 
      }
    ),
    tags$hr(style = "border-color: black;"),
    br(),
    withSpinner(uiOutput("Ratios")),
    # Note to self: I am for now excluding some methods, because I need to add code to calculate some columns for those, namely ratios.
    # This should be remedied asap, especially since there include such community favourites as IQ (= MaxLFQ) and Top3!!!
    br(),
    tags$hr(style = "border-color: black;"),
    fluidRow(
      column(2,
             if (moreThan1Sample) {
               checkboxInput("runProfPlots", "Draw protein profile plots?",
                             runProfPlots, "100%")
             },
             checkboxInput("runRankAbundPlots", "Draw protein ranked abundance plots?",
                           runRankAbundPlots, "100%")),
      if (moreThan1Sample) {
        column(2, radioButtons("Clustering",
                               "Heatmaps: clustering method",
                               klustChoices,
                               klustChoices[1], TRUE, "100%"))
      },
      if (Annotate) {
        column(2, checkboxInput("runGSEA", "Run Gene Set Enrichment Analysis (GSEA)?",
                                runGSEA, "100%"))
      }
    ),
    if (moreThan1Sample) {
      checkboxInput("Venn_Obs", "Draw Venn diagrams?", Venn_Obs, "100%")
    },
    br(),
    tags$hr(style = "border-color: black;"),
    withSpinner(uiOutput("GO")),
    tags$hr(style = "border-color: black;"),
    withSpinner(uiOutput("CytoScape")),
    h4("Post-translational modifications (PTMs)"),
    fluidRow(column(2, pickerInput("Mod2Write",
                                   "Select PTM(s) for which to write a specific tab in the report",
                                   Modifs$`Full name`, ptmDflt2, TRUE,
                                   pickerOptions(title = "Search me",
                                                 `live-search` = TRUE,
                                                 actionsBox = TRUE,
                                                 deselectAllText = "Clear search"))),
             #column(2, checkboxInput("PTMsReNorm", "Re-normalize modified peptides ratios to those of parent Protein Group(s)?", TRUE, "100%"))
    ),
    br(),
    tags$hr(style = "border-color: black;"),
    checkboxInput("AdvOptOn", "Advanced options", tstAdvOpt),
    withSpinner(uiOutput("AdvOpt")),
    br(),
    actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
    br(),
    br()
  )
}
if (exists("appRunTest")) { rm(appRunTest) }
server1 <- function(input, output, session) {
  # Initialize variables to create in main environment
  m4Quant <- reactiveVal(Mod4Quant)
  m2Xclud <- reactiveVal(Mod2Xclud)
  ADVOPT <- reactiveVal(tstAdvOpt)
  PARAM <- reactiveVal(AnalysisParam)
  #
  updtQuantMsg <- function(reactive = TRUE) {
    myQuant <- { if (reactive) { input$Quant_algorithm } else { quantAlgo } }
    quantMsg <- allQuantAlgos$Help[match(myQuant, allQuantAlgos$Algorithm)]
    renderUI(em(quantMsg))
  }
  #
  # Ratios
  output$Ratios <- renderUI({
    if (MakeRatios) {
      lst <- list(list(h4("Ratios analysis")),
                  list(fluidRow(column(2,
                                       numericInput("RatiosThresh", "Fold change threshold (log2)",
                                                    RatiosThresh, 0, width = "100%")),
                                column(2,
                                       checkboxInput("RatiosThresh_2sided",
                                                     "Are you interested in down-regulated protein groups as well?",
                                                     RatiosThresh_2sided, "100%")))))
    } else { lst <- list() }
    return(lst)
  })
  #
  # GO
  output$GO <- renderUI({
    lst <- list(list(br()))
    if (Annotate) {
      lst <- list(
        list(h4("GO terms enrichment"),
             fluidRow(column(1,
                             checkboxInput("GOenrich", "GO enrichment", globalGO, "100%"),
                             checkboxInput("runClueGO", "run ClueGO enrichment (NB: this is rather slow)", runClueGO, "100%")),
                      column(2,
                             pickerInput("GO.tabs", "GO terms of interest",
                                         allGO, GO_filter1, TRUE,
                                         pickerOptions(title = "Search me",
                                                       `live-search` = TRUE,
                                                       actionsBox = TRUE,
                                                       deselectAllText = "Clear search"))),
                      column(1,
                             checkboxInput("GO2Int", "Use GO terms to define list of proteins of interest?",
                                           AnalysisParam$GO.terms.for.proteins.of.interest , "100%"))))
      )
    }
    return(lst)
  })
  #
  # CytoScape
  output$CytoScape <- renderUI({
    if (!length(CytoScExe)) {
      lst <- list(list(column(2, shinyFilesButton("CytoScVers2",
                                                  em("No Cytoscape.exe detected, select one (optional)"),
                                                  "", FALSE))))
    }
    if (length(CytoScExe) == 1) {
      lst <- list(list(column(2, em(paste0("Version detected: ", CytoScExe)))))
    }
    if (length(CytoScExe) > 1) {
      lst <- list(list(fluidRow(column(2, pickerInput("CytoScVers1",
                                                      "Choose Cytoscape version:",
                                                      CytoScExe,
                                                      CytoScExe[1],
                                                      FALSE)))))
    }
    return(lst)
  })
  #
  # Optional input files
  updtOptOn <- function(reactive = TRUE) {
    if (reactive) { tst <- ADVOPT() } else { tst <- tstAdvOpt }
    if (tst) {
      lst <- vector("list", 1)
      lst[[1]] <- list(fluidRow(column(2,
                                       em("->"), 
                                       shinyFilesButton("CustPG",
                                                        em("Custom Protein Groups"), "", FALSE),
                                       br(),
                                       em("Allows \"cheating\" with the naive Protein Groups assembly algorithm."), br(),
                                       em("Useful when e.g. samples express from a custom construct to which matching peptides should be assigned in priority."), br(),
                                       em("Should be a table with two columns: \"Leading protein IDs\" (\";\"-separated) and \"Priority\" (integer)."),
                                       br(), br(),
                                       em("Current selection = "),
                                       span(AnalysisParam$Custom.PGs, style = "color:blue", .noWS = "outside"),
                                       br()
      ),
      column(2,
             em("->"),
             shinyFilesButton("CRAPome", em("CRAPome filter"), "", FALSE), br(),
             em("CRAPome-like filter: 1 column table of protein accessions to mark as contaminants."), br(),
             em("Column name = \"Protein ID\" or \"Protein IDs\", use \";\" if including more than one ID per row)."),
             br(), br(),
             em("Current selection = "),
             span(AnalysisParam$CRAPome_file, style = "color:blue", .noWS = "outside"),
             br()
      )
      ))
    } else { lst <- list(list(em(""))) }
    renderUI(lst)
  }
  output$AdvOpt <- updtOptOn(FALSE)
  # Event observers
  observeEvent(input$AdvOptOn, {
    ADVOPT(input$AdvOptOn)
    output$AdvOpt <- updtOptOn()
  })
  observe({ shinyFileChoose(input, "CustPG", roots = getVolumes(), filetypes = "csv")
    {
      tmp <- input$CustPG
      if ((!is.null(tmp))&&(is.list(tmp))) {
        tmp <- parseFilePaths(getVolumes(), tmp)$datapath
        Par <- PARAM()
        Par$Custom.PGs <- normalizePath(tmp, winslash = "/")
        PARAM(Par)
      }
  }
  })
  observe({ shinyFileChoose(input, "CRAPome", roots = getVolumes(), filetypes = "csv" )
    {
      tmp <- input$CRAPome
      if ((!is.null(tmp))&&(is.list(tmp))) {
        tmp <- parseFilePaths(getVolumes(), tmp)$datapath
        Par <- PARAM()
        Par$CRAPome_file <- normalizePath(tmp, winslash = "/")
        PARAM(Par)
      }
  }
  })
  # Pepper
  observeEvent(input$runPepper, {
    runPepper <- as.logical(input$runPepper)
    assign("runPepper", runPepper, envir = .GlobalEnv)
    Par <- PARAM()
    Par$runPepper <- runPepper
    PARAM(Par)
  })
  # Impute?
  observeEvent(input$Impute, {
    Impute <- as.logical(input$Impute)
    assign("Impute", input$Impute, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Pep.Impute <- Impute
    PARAM(Par)
  })
  # Update PSM-to-Protein matches?
  observeEvent(input$Update_Prot_matches, {
    Update_Prot_matches <- input$Update_Prot_matches
    assign("Update_Prot_matches", Update_Prot_matches, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Update_Prot_matches <- Update_Prot_matches
    PARAM(Par)
  })
  # Clustering method
  observeEvent(input$Clustering, {
    assign("KlustMeth", match(input$Clustering, klustChoices), envir = .GlobalEnv)
  })
  observeEvent(input$removeMBR, {
    removeMBR <- as.logical(input$removeMBR)
    assign("removeMBR", removeMBR, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"Proteins list: remove match-between-runs" <- removeMBR
    PARAM(Par)
  })
  # Quantitation
  observeEvent(input$Quant_algorithm, {
    assign("quantAlgo", input$Quant_algorithm, envir = .GlobalEnv)
    if (quantAlgo == "limpa") {
      updateCheckboxInput(inputId = "Impute", value = FALSE)
      # updatePickerInput(inputId = "Norma.Prot.Ratio.to.proteins", selected = NULL)
      # updatePickerInput(inputId = "Norma.Prot.Ratio.to.GO", selected = NULL)
      # updateCheckboxInput(inputId = "Norma.Prot.Ratio.to.Biot", value = FALSE)
      updateRadioButtons(inputId = "TtstPval", selected = "Moderated")
      shinyjs::disable("Impute")
      # shinyjs::disable("Norma.Prot.Ratio.to.proteins")
      # shinyjs::disable("Norma.Prot.Ratio.to.GO")
      # shinyjs::disable("Norma.Prot.Ratio.to.Biot")
      #
      # Where to deal with batch-correction?
      # Not here: it should be encoded in the design matrix in general rather than corrected for, whether using limpa or not for quant.
    } else {
      shinyjs::enable("Impute")
      # shinyjs::enable("Norma.Prot.Ratio.to.proteins")
      # shinyjs::enable("Norma.Prot.Ratio.to.GO")
      # shinyjs::enable("Norma.Prot.Ratio.to.Biot")
    }
    Par <- PARAM()
    Par$Quant_algorithm <- quantAlgo
    PARAM(Par)
    output$QuantMsg <- updtQuantMsg()
  })
  observeEvent(input$ReScaling, {
    assign("reScAlgo", input$ReScaling, envir = .GlobalEnv)
    if (input$ReScaling == "topN") {
      shinyjs::enable("topN")
      shinyjs::enable("topN_correct")
    } else {
      shinyjs::disable("topN")
      shinyjs::disable("topN_correct")
    }
    Par <- PARAM()
    Par$ReScaling_algorithm <- reScAlgo
    PARAM(Par)
  })
  # - Peptides eligible for quant
  #   - Classes
  observeEvent(input$Prot.Quant.Use, {
    assign("Pep4Quant", Pep4QuantOpt[input$Prot.Quant.Use])
    Par <- PARAM()
    Par$Prot.Quant.Use <- input$Prot.Quant.Use
    PARAM(Par)
  })
  #   - PTMs to use for PG Quant
  observeEvent(input$PTMsQuant, {
    assign("ptmDflt1", input$PTMsQuant, envir = .GlobalEnv)
    m4Quant(Modifs$Mark[match(unlist(input$PTMsQuant), Modifs$`Full name`)])
    m2Xclud(set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                         c("Mark", "Where")))
  }, ignoreNULL = FALSE)
  #   - Number of samples in which observed
  observeEvent(input$PepFoundInAtLeast, {
    PepFoundInAtLeast <- as.integer(input$PepFoundInAtLeast)
    assign("PepFoundInAtLeast", PepFoundInAtLeast, envir = .GlobalEnv)
    Par <- PARAM()
    Par$PepFoundInAtLeast <- PepFoundInAtLeast
    PARAM(Par)
  })
  #   - Number of samples in at least one sample group in which observed
  observeEvent(input$PepFoundInAtLeastGrp, {
    PepFoundInAtLeastGrp <- as.integer(input$PepFoundInAtLeastGrp)
    assign("PepFoundInAtLeastGrp", PepFoundInAtLeastGrp, envir = .GlobalEnv)
    Par <- PARAM()
    Par$PepFoundInAtLeastGrp <- PepFoundInAtLeastGrp
    PARAM(Par)
  })
  #   - Number of peptides to use for quantitation
  observeEvent(input$N_Pep, {
    N_Pep <- as.integer(input$N_Pep)
    assign("N_Pep", N_Pep, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"N. of peptidoforms for quantitation" <- N_Pep
    PARAM(Par)
  })
  #   - Number of unique peptides to use, if available, to the exclusion of non-unique ones
  observeEvent(input$N_unique_Pep, {
    N_unique_Pep <- as.integer(input$N_unique_Pep)
    assign("N_unique_Pep", N_unique_Pep, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Use.N.unique <- N_unique_Pep
    PARAM(Par)
  })
  # - Top N parameters
  #   - N = ?
  observeEvent(input$topN, {
    topN <- as.integer(input$topN)
    assign("topN", topN, envir = .GlobalEnv)
    Par <- PARAM()
    Par$topN <- topN
    PARAM(Par)
  })
  #   - topN correct?
  observeEvent(input$topN_correct, {
    topN_correct <- as.logical(input$topN_correct)
    assign("topN_correct", topN_correct, envir = .GlobalEnv)
    Par <- PARAM()
    Par$topN_correct <- topN_correct
    PARAM(Par)
  })
  # - Proteome ruler
  #   - Yes/no?
  observeEvent(input$ProtRul, {
    protrul <- as.logical(input$ProtRul)
    if (protrul) { shinyjs::enable("ProtRulNuclL") }
    if (!protrul) { shinyjs::disable("ProtRulNuclL") }
    assign("protrul", protrul, envir = .GlobalEnv)
    Par <- PARAM()
    Par$ProtRul <- protrul
    PARAM(Par)
  })
  #   - Nucleosome length
  observeEvent(input$ProtRulNuclL, {
    ProtRulNuclL <- as.integer(input$ProtRulNuclL)
    assign("ProtRulNuclL", ProtRulNuclL, envir = .GlobalEnv)
    Par <- PARAM()
    Par$ProtRulNuclL <- as.integer(input$ProtRulNuclL)
    PARAM(Par)
  })
  # # WGCNA
  # observeEvent(input$runWGCNA, {
  #   runWGCNA <- as.logical(input$runWGCNA)
  #   assign("runWGCNA", runWGCNA, envir = .GlobalEnv)
  #   Par <- PARAM()
  #   Par$runWGCNA <- runWGCNA
  #   PARAM(Par)
  # })
  # GSEA
  observeEvent(input$runGSEA, {
    runGSEA <- as.logical(input$runGSEA)
    assign("runGSEA", runGSEA, envir = .GlobalEnv)
    Par <- PARAM()
    Par$runGSEA <- runGSEA
    PARAM(Par)
  })
  # Profile plots
  observeEvent(input$runProfPlots, {
    runProfPlots <- as.logical(input$runProfPlots)
    assign("runProfPlots", runProfPlots, envir = .GlobalEnv)
    Par <- PARAM()
    Par$runProfPlots <- runProfPlots
    PARAM(Par)
  })
  # Ranked abundance plots
  observeEvent(input$runRankAbundPlots, {
    runRankAbundPlots <- as.logical(input$runRankAbundPlots)
    assign("runRankAbundPlots", runRankAbundPlots, envir = .GlobalEnv)
    Par <- PARAM()
    Par$runRankAbundPlots <- runRankAbundPlots
    PARAM(Par)
  })
  #
  observeEvent(input$RatiosThresh, {
    RatiosThresh <- as.numeric(input$RatiosThresh)
    assign("RatiosThresh", RatiosThresh, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"Ratios analysis - threshold" <- RatiosThresh
    PARAM(Par)
  })
  observeEvent(input$RatiosThresh_2sided, {
    RatiosThresh_2sided <- as.logical(input$RatiosThresh_2sided)
    assign("RatiosThresh_2sided", RatiosThresh_2sided, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"Ratios analysis - threshold is two-sided" <- RatiosThresh_2sided
    PARAM(Par)
  })
  # Normalizations
  observeEvent(input$prtNorm, {
    NormalizePG <- as.logical(input$prtNorm)
    assign("NormalizePG", NormalizePG, envir = .GlobalEnv)
    Par <- PARAM()
    Par$NormalizePG <- NormalizePG
    PARAM(Par)
  })
  # Proteins of interest
  observeEvent(input$IntProt, {
    prot.list <- db$`Protein ID`[dbOrd][match(input$IntProt, protHeads)]
    assign("prot.list", prot.list, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Prot.list <- prot.list
    PARAM(Par)
  }, ignoreNULL = FALSE)
  # GO enrichment
  observeEvent(input$GOenrich, {
    if (input$GOenrich) { shinyjs::enable("runClueGO") }
    if (!input$GOenrich) { shinyjs::disable("runClueGO") }
    globalGO <- as.logical(input$GOenrich)
    enrichGO <- globalGO&MakeRatios
    assign("globalGO", globalGO, envir = .GlobalEnv)
    assign("enrichGO", enrichGO, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"GO terms enrichment analysis" <- globalGO
    PARAM(Par)
  })
  observeEvent(input$runClueGO, {
    runClueGO <- input$runClueGO&input$GOenrich
    assign("runClueGO", runClueGO, envir = .GlobalEnv)
    Par <- PARAM()
    Par$runClueGO <- runClueGO
    PARAM(Par) 
  })
  # GO terms of interest
  if (Annotate) {
    observeEvent(input$GO.tabs, {
      GO_filter <- allGO2[match(input$GO.tabs, allGO)]
      assign("GO_filter", GO_filter, envir = .GlobalEnv)
      if (length(GO_filter1)) {
        mY <- match(GO_filter1, allGO)
        mN <- which(!allGO %in% GO_filter1)
        assign("allGO", allGO[c(mY, mN)], envir = .GlobalEnv)
        assign("allGO2", allGO2[c(mY, mN)], envir = .GlobalEnv)
      }
      Par <- PARAM()
      Par$"GO terms of interest" <- GO_filter
      PARAM(Par)
    }, ignoreNULL = FALSE)
    observeEvent(input$GO2Int, {
      Par <- PARAM()
      Par$GO.terms.for.proteins.of.interest <- input$GO2Int
      PARAM(Par)
    })
  }
  #
  observeEvent(input$Venn_Obs, {
    Venn_Obs <- as.logical(input$Venn_Obs)
    assign("Venn_Obs", Venn_Obs, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"Venn diagrams: observed" <- Venn_Obs
    PARAM(Par)
  })
  # PTMs to use for PG Quant
  observeEvent(input$PTMsQuant, {
    m4Quant(Modifs$Mark[match(unlist(input$PTMsQuant), Modifs$`Full name`)])
    m2Xclud(set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                         c("Mark", "Where")))
  }, ignoreNULL = FALSE)
  # PTMs to write a tab for
  observeEvent(input$Mod2Write, {
    Mod2Write <- Modifs$Mark[match(input$Mod2Write, Modifs$`Full name`)]
    assign("Mod2Write", Mod2Write, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Mod2Write <- Mod2Write
    PARAM(Par)
  }, ignoreNULL = FALSE)
  # CytoScape
  if (!length(CytoScExe)) {
    observeEvent(input$CytoScVers1, {
      CytoScExe <- input$CytoScVers1
      assign("CytoScExe", CytoScExe, envir = .GlobalEnv)
      Par <- PARAM()
      Par$CytoscapePath <- CytoScExe
      PARAM(Par)
    })
  }
  if (length(CytoScExe) > 1) {
    observe({ shinyFileChoose(input, "CytoScVers2", roots = getVolumes(), filetypes = "exe")
      {
        tmp <- input$CytoScVers2
        if ((!is.null(tmp))&&(is.list(tmp))) {
          assign("CytoScExe", input$CytoScVers2, envir = .GlobalEnv)
          tmp <- parseFilePaths(getVolumes(), tmp)$datapath
          Par <- PARAM()
          Par$CytoscapePath <- normalizePath(tmp, winslash = "/")
          PARAM(Par)
        }
    }
    })
  }
  # Save
  observeEvent(input$saveBtn, {
    Par <- PARAM()
    assign("AnalysisParam", Par, envir = .GlobalEnv)
    assign("Mod4Quant", m4Quant(), envir = .GlobalEnv)
    assign("Mod2Xclud", m2Xclud(), envir = .GlobalEnv)
    assign("appRunTest", TRUE, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
appTxt1 <- gsub("myApp", "myApp1", gsub("\\(ui", "(ui1", gsub(", server", ", server1", runApp)))
runKount <- 0
while ((!runKount)||(!exists("appRunTest"))) {
  ui1 <- make_ui() # Update ui with new values
  eval(parse(text = appTxt1), envir = .GlobalEnv)
  runKount <- runKount+1
}
#
# Post-processing
prot.list %<o% AnalysisParam$Prot.list
if (AnalysisParam$GO.terms.for.proteins.of.interest) {
  tmpGO <- AnalysisParam$"GO terms of interest"
  GO_prot.list <- list()
  GO_prot.list$Offspring <- lapply(tmpGO, function(x) {
    ont <- Ontology(x)
    x <- c(x, get(paste0("GO", ont, "OFFSPRING"))[[x]])
    x <- x[which(!is.na(x))]
    return(x)
  })
  tmpGO2 <- listMelt(strsplit(db$`GO-ID`, ";"), db$`Protein ID`)
  GO_prot.list$Proteins <- lapply(GO_prot.list$Offspring, function(x) {
    unique(tmpGO2$L1[which(tmpGO2$value %in% unlist(x))])
  })
  prot.list <- unique(c(prot.list, unlist(GO_prot.list$Proteins)))
}
prot.list.Cond %<o% (length(prot.list) > 0)
if (prot.list.Cond) {
  m <- match(prot.list, db$`Protein ID`)
  x <- aggregate(db$`Protein ID`[m], list(db$`Common Name`[m]), c)
  IDs.list %<o% setNames(as.list(x$x), x$Group.1)
  prot.names %<o% names(IDs.list)
  db$"Potential contaminant"[which(db$`Protein ID` %in% prot.list)] <- ""
}
if (runClueGO&&!enrichGO) { runClueGO <- FALSE}
#
custPGs_file %<o% AnalysisParam$Custom.PGs
custPGsTst <- (!is.na(custPGs_file))&(file.exists(custPGs_file))
if (custPGsTst) {
  if (dirname(custPGs_file) != wd) {
    file.copy(custPGs_file, wd)
    custPGs_file <- paste0(wd, "/", basename(custPGs_file))
  }
  custPGs %<o% read.delim(custPGs_file, check.names = FALSE)
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.delim(custPGs_file, check.names = FALSE, sep = ",") }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_file, check.names = FALSE) }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_file, check.names = FALSE, sep = "\t") }
  if (colnames(custPGs)[1] != "Leading protein IDs") {
    warning("I could not make sense of that file! Skipping.")
    custPGs <- NA
  }
} else { custPGs %<o% NA }
#
if (prot.list.Cond) {
  temp <- db[which(db$`Protein ID` %in% prot.list),]
  writeFasta(temp, paste0(wd, "/Proteins of interest.fasta"))
  AnalysisParam$"Proteins list: proteins" <- IDs.list
  AnalysisParam$"Proteins list: names" <- prot.names
}

# GO terms of interest (for profile plots, LFQ plots and similar, GO-specific tabs...)
GO_filt %<o% FALSE
if (exists("GO_filter")) { GO_filt <- length(GO_filter) > 0 }
GO_filter %<o% GO_filter
if (GO_filt) {
  library(GO.db)
  AllTerms %<o% unique(unlist(strsplit(db$`GO-ID`, ";")))
  AllTermNames %<o% unique(unlist(strsplit(db$GO, ";")))
}

# Venn diagrams
AnalysisParam$"Venn diagrams: observed" <- Venn_Obs
Venn_Ratios %<o% (Venn_Obs & MakeRatios)
AnalysisParam$"Venn diagrams: regulated (up)" <- Venn_Ratios
AnalysisParam$"Venn diagrams: regulated (down)" <- Venn_Ratios & RatiosThresh_2sided
#
AnalysisParam$"PTMs eligible for quantitation" <- Mod4Quant
AnalysisParam$"PTMs to exclude from quantitation" <- Mod2Xclud

# Experiments
SamplesMapPath %<o% paste0(wd, "/SamplesMap.csv")
tst <- (("Experiment" %in% colnames(ev))&&(sum(is.na(ev$Experiment)) == 0)) 
if (tst) {
  if ("Reference" %in% colnames(SamplesMap)) {
    SamplesMap <- SamplesMap[which(!is.na(SamplesMap$Reference)),]
  }
  Exp %<o% SamplesMap$Experiment
  w1 <- which(ev$Experiment %in% Exp)
  w2 <- which(!ev$Experiment %in% Exp)
  if (length(w2)) {
    warning(paste0("Removing ", length(w2), " evidences from undefined experiments (check \"Reference\" column of Samples map)..."))
    ev <- ev[w1,]
  }
} else {
  Exp %<o% "Exp1"
  ev$Experiment <- Exp
  SamplesMap <- data.frame("Experiment" = Exp,
                           "Ratios group" = 1,
                           "Reference" = TRUE,
                           "Negative Filter" = FALSE,
                           "Use" = TRUE,
                           check.names = FALSE)
  tst <- try(write.csv(SamplesMap, file = SamplesMapPath, row.names = FALSE), silent = TRUE)
  while ("try-error" %in% class(tst)) {
    dlg_message(paste0("File \"", SamplesMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
    tst <- try(write.csv(SamplesMap, file = SamplesMapPath, row.names = FALSE), silent = TRUE)
  }
}
if (identical(c("MQ.Exp", "Experiment") %in% colnames(SamplesMap), c(TRUE, FALSE))) {
  # Here we use Experiment as synonym for MQ.Exp
  colnames(SamplesMap)[which(colnames(SamplesMap) == "MQ.Exp")] <- "Experiment"
}
