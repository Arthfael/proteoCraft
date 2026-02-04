#### Define analysis parameters
require(shiny)
require(shinyFiles)
require(shinycssloaders)
require(shinyjs)
require(shinyBS)
require(htmlwidgets)
#
moreThan1Sample %<o% TRUE
#
# Boolean functions to check parameter values
Src <- paste0(libPath, "/extdata/R scripts/Sources/parBooleans.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# PCA prior to shiny app
Src <- paste0(libPath, "/extdata/R scripts/Sources/rep_Parameters_editor_PCA.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# Protein headers for shiny
Src <- paste0(libPath, "/extdata/R scripts/Sources/protHeaders_for_shiny.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# Proteins of interest
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# Targets
# Sometimes the user does not fill the Target factor with valid protein IDs... but this is what we would actually need.
# Here, if necessary, we will remap those to valid IDs:
Src <- paste0(libPath, "/extdata/R scripts/Sources/Targets.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# KnockOut, KnockIn or KnockDown
tst <- tolower(gsub("[- _]", "", Factors))
if (sum(c("knockout", "knockin", "knockdown") %in% tst)) {
  w <- which(c("knockout", "knockin", "knockdown") %in% tst)
  # There should be only one for now, because all three share the same 3-characters root = "Kno"
  # This should evolve but will be difficult, knowing how complex this script is now.
  prot.list %<o% unique(c(unique(Exp.map[[Factors["Kno"]]]), prot.list))
  prot.list_pep %<o% unique(c(unique(Exp.map[[Factors["Kno"]]]), prot.list_pep))
}
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList2.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# Protein headers for shiny (update)
Src <- paste0(libPath, "/extdata/R scripts/Sources/protHeaders_for_shiny.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
# Defaults
nSmpls <- nrow(Exp.map)
RefRat_Mode %<o% "2" # Values: RefRat_Mode = "2" or "1" # For now not a user-modifiable parameter, however this may change
StudentRoot %<o% "Student's t-test -log10(Pvalue) - "
WelchRoot %<o% "Welch's t-test -log10(Pvalue) - "
modRoot %<o% "Moderated t-test -log10(Pvalue) - "
deqmsRoot %<o% "DEqMS mod. t-test -log10(Pvalue) - "
permRoot %<o% "Permutations t-test -log10(Pvalue) - "
samRoot %<o% "SAM -log10(Pvalue) - "
odpRoot %<o% "ODP -log10(Pvalue) - "
lrtRoot %<o% "LRT -log10(Pvalue) - "
#
pvalue.col %<o% c(StudentRoot, WelchRoot, modRoot, permRoot, samRoot, odpRoot, lrtRoot)
names(pvalue.col) <- vapply(pvalue.col, function(x) { unlist(strsplit(x, "\\.|\\'|\\ "))[1] }, "")
ParamFls <- c(paste0(wd, "/Parameters.csv"),
              paste0(libPath, "/extData/Parameters_template.csv"))
ParamPath %<o% ParamFls[1] # The parameters file we will save to!
ParamFl %<o% ParamFls[1] # The parameters file we (re-)load!
if (!file.exists(ParamFl)) { ParamFl <- ParamFls[2] }
Param %<o% Param.load(ParamFl)
tmp <- read.csv(ParamFl, header = FALSE)
if (ncol(tmp) == 3) {
  Param_Help <- tmp$V3
  names(Param_Help) <- tmp$V1
} else {
  tmp <- read.csv(ParamFls[2], header = FALSE)
  Param_Help <- vapply(colnames(Param), function(x) {
    if (x %in% tmp$V1) { return(tmp$V3[match(x, tmp$V1)]) }
    return("")
  }, "")
}
Param$vCPUs <- N.clust
Param$WD <- wd
if (ParamFl == ParamFls[2]) {
  Param$WD <- wd
  Param$Project <- dtstNm
  Param$Fix.MQ.Isobaric.labels <- FALSE
  Param$Type <- WorkFlow
  Param$Label <- LabelType
  Param$MQ.Experiments <- paste(MQ.Exp, collapse = ";")
  Param$Search.DB <- paste(fastasTbl$Full, collapse = ";")
  Param$Search.DB.type <- paste(fastasTbl$Type, collapse = ";")
  Param$Search.DB.species <- paste(fastasTbl$Species, collapse = ";") 
  Param$Two.sided <- !(WorkFlow %in% c("PULLDOWN", "BIOID"))
  Param$Min.Pep.Size <- MinPepSz
  Param$PSMs <- paste(PSMsFls, collapse = ";")
  if (LabelType == "Isobaric") {
    Param$Label <- IsobarLab #
    Param$Label.Multiplicity <- length(get(IsobarLab)) #
    Param$Label.Purities.file <- ""
  }
  if ("Target" %in% Factors) {
    tmp <- FactorsLevels$Target
    tmp <- tmp[which((!is.na(tmp))&(tmp != "NA"))]
    Param$Prot.list <- paste(unique(c(tmp, unlist(strsplit(Param$Prot.list, ";")))), collapse = ";")
    Param$Prot.list_pep <- paste(unique(c(tmp, unlist(strsplit(Param$Prot.list_pep, ";")))), collapse = ";")
  }
  ptmDflt2 <- Param$PTM.analysis <- paste(grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE), collapse = ";")
  #
  Param$GO.enrichment <- TRUE
} else {
  ptmDflt2 <- ""
  if ("PTM.analysis" %in% colnames(Param)) {
    tmp <- unlist(strsplit(as.character(Param$PTM.analysis), ";"))
    tmp <- tmp[which(tmp %in% Modifs$`Full name`)]
    ptmDflt2 <- tmp
  }
}
if (!"PTM.analysis_Norm" %in% colnames(Param)) { Param$PTM.analysis_Norm <- TRUE } 
tmp <- Param$PTM.analysis_Norm
if (!validLogicPar("tmp")) { Param$PTM.analysis_Norm <- TRUE }
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
if (!"PTM.analysis" %in% colnames(Param)) { Param$PTM.analysis <- paste(ptmDflt2, collapse = ";") }
if ("Output" %in% colnames(Param)) { Param$Output <- NULL } # Deprecated
tmp <- PSMsFls
Param$Search.dir <- paste(inDirs, collapse = ";")
#if ((Param$Label == "LFQ")&&(isDIA)) { Param$Label <- "DIA" } # Nope! isDIA can be length > 1 now! 
Param$WD <- wd
#
goDflt <- suppressWarnings(as.logical(Param$GO.enrichment))
if ((!is.logical(goDflt))||(is.na(goDflt))) { goDflt <- TRUE }
#
fTstDflt <- { if (exists("F.test")) { F.test } else { Param$F.test } }
if (!length(fTstDflt)) { fTstDflt <- TRUE }
fTstDflt <- suppressWarnings(as.logical(fTstDflt[1]))
if (is.na(fTstDflt)) { fTstDflt <- TRUE }
#
if ((!validLogicPar("Nested"))&&("Ratios.Groups_Nested" %in% colnames(Param))) {
  tmp1 <- as.logical(Param$Ratios.Groups_Nested)
  if (validLogicPar("tmp1")) { Nested <- tmp1 }
}
if (!validLogicPar("Nested")) { Nested <- WorkFlow != "Regulation" }
Nested %<o% Nested
Param$Ratios.Groups_Nested <- Nested
#
#Param$Ratios.Groups.Ref.Aggregate.Level <- "AUTOFACT"
#Param$Ratios.Groups.Ref.Aggregate.Level <- "Exp;Con;Rep"
tmp <- unlist(strsplit(Param$Ratios.Groups.Ref.Aggregate.Level, ";"))
tst <- (length(tmp) > 1)||(!tmp %in% c("AUTOFACT", "MAP2FACTS"))
if (tst) {
  tst <- sum(!tmp %in% substr(Factors, 1, 3)) == 0
  if (!tst) { tmp <- "MAP2FACTS" }
}
if ((length(tmp) == 1)&&(tmp == "AUTOFACT")) {
  Param$Ratios.Groups.Ref.Aggregate.Level <- paste(substr(Factors, 1, 3), collapse = ";")
}
klustChoices %<o% c("K-means", "hierarchical")
KlustMeth %<o% 1 # Changed from 2
#
if ((!validLogicPar("saintExprs"))&&("saintExprs" %in% colnames(Param))) {
  tmp1 <- as.logical(Param$saintExprs)
  if (validLogicPar("tmp1")) { saintExprs <- (WorkFlow %in% c("PULLDOWN", "BIOID")) }
}
if (!validLogicPar("saintExprs")) { saintExprs <- WorkFlow != "Regulation" }
saintExprs %<o% saintExprs
Param$saintExprs <- saintExprs
#
if (Annotate) {
  allGO <- unique(unlist(strsplit(db$GO[which(!is.na(db$GO))], ";")))
  allGO2 <- paste0("GO:", gsub(".* \\[GO:|\\]$", "", allGO))
  dftlGO2 <- unique(unlist(strsplit(Param$GO.tabs, ";")))
  dftlGO2 <- dftlGO2[which(dftlGO2 %in% allGO2)]
  dftlGO <- allGO[match(dftlGO2, allGO2)]
  w <- c(which(allGO %in% dftlGO),
         which(!allGO %in% dftlGO))
  allGO <- allGO[w]; allGO2 <- allGO2[w]
  # allGO = Names
  # allGO2 = IDs
  if (("Norma.Prot.Ratio.to.GO" %in% colnames(Param))&&
      ("character" %in% class(Param$Norma.Prot.Ratio.to.GO))&&
      (nchar(Param$Norma.Prot.Ratio.to.GO))) {
    tmp <- unlist(strsplit(Param$Norma.Prot.Ratio.to.GO, ";"))
    tmp <- tmp[which(tmp %in% allGO2)]
    if (length(tmp)) {
      nrm2GO <- allGO[match(tmp, allGO2)]
      nrm2GOall <- c(nrm2GO,
                     allGO[which(!allGO %in% nrm2GO)])
    }
  } else {
    nrm2GOall <- allGO
    nrm2GO <- c()
  }
}
if (!"GO.terms.for.proteins.of.interest" %in% colnames(Param)) { Param$GO.terms.for.proteins.of.interest <- FALSE }
if (!"Amica" %in% colnames(Param)) { Param$Amica <- TRUE }
#if (!"Mirror.Ratios" %in% colnames(Param)) { Param$Mirror.Ratios <- FALSE }
if (!"Custom.PGs" %in% colnames(Param)) { Param$Custom.PGs <- "" }
if (!"TrueDisc_filter" %in% colnames(Param)) { Param$TrueDisc_filter <- "" }
DiscFiltModes %<o% c("Positive filter", "Negative filter", "Filter column")
DiscFiltModesHlp <- c(" - Positive filter: TRUE means values should not be set to NA for PGs matching this protein - values are set to NA for PGs without match or with match to FALSE)",
                      " - Negative filter: TRUE means values should be set to NA for PGs matching this protein - values are left unchanged for PGs without match or with match to FALSE)",
                      " - Filter column: Values are unaffected, but a new column is created marking with \"+\" proteins found in the provided filter.")
if (!"TrueDisc_filter_mode" %in% colnames(Param)) { Param$TrueDisc_filter_mode <- DiscFiltModes[1] }
if (!"CRAPome_file" %in% colnames(Param)) { Param$CRAPome_file <- "" }
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
  if (!p %in% colnames(Param)) { Param[[p]] <- dflt } else {
    tmp <- Param[[p]]
    if (!validLogicPar("tmp")) { Param[[p]] <- dflt }
  }
}
#
# Peptide classes to use for quantitation
Pep4QuantOpt %<o% setNames(c("Unique peptide IDs", "Razor peptide IDs", "Peptide IDs"),
                           c("Unique", "Razor", "All"))
if ((!validCharPar("Pep4Quant", Pep4QuantOpt))&&("Prot.Quant.Use" %in% colnames(Param))) {
  tmp1 <- toupper(gsub(" |_|-|\\.", "", as.character(Param$Prot.Quant.Use)))
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
Param$Prot.Quant.Use <- names(Pep4Quant)
#
# How many peptide do we need to quantify a protein?
if ((!validIntegPar("N_Pep"))&&("N. of peptidoforms for quantitation" %in% colnames(Param))) {
  tmp1 <- Param$"N. of peptidoforms for quantitation"
  if (validIntegPar("tmp1")) { N_Pep <- as.integer(tmp1) }
}
if (!validIntegPar("N_Pep")) {
  # The default is to use ALL of the information, even if we only have a single peptide
  # We will flag proteins flagged with 1 peptide as dodgy
  N_Pep <- 1
}
N_Pep %<o% N_Pep
Param$"N. of peptidoforms for quantitation" <- N_Pep
#
# How many unique peptides (if available) to use to the exclusivity of any others
# If set to 0, we will use peptides regardless of whether unique or not
if ((!validIntegPar("N_unique_Pep", 0))&&("Use.N.unique" %in% colnames(Param))) {
  tmp1 <- Param$Use.N.unique
  if (validIntegPar("tmp1", 0)) { N_unique_Pep <- as.integer(tmp1) }
}
if (!validIntegPar("N_unique_Pep", 0)) {
  N_unique_Pep <- 3 # ... which means if we have 3 unique peptides we will not use any non-unique ones
}
N_unique_Pep %<o% N_unique_Pep
Param$Use.N.unique <- N_unique_Pep
# Minimum number of peptides for quantitation?
if ((!validIntegPar("N_Pep"))&&("Min.N.pep" %in% colnames(Param))) {
  tmp1 <- Param$Min.N.pep
  if (validIntegPar("tmp1")) { N_Pep <- as.integer(tmp1) }
}
if (!validIntegPar("N_Pep")) {
  N_Pep <- 1 # The default is to use ALL of the information, even if we only have a single peptide
}
N_Pep %<o% N_Pep
Param$Min.N.pep <- N_Pep
#
ptmDflt1 <- grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE, invert = TRUE)
Mod4Quant %<o% Modifs$Mark[match(ptmDflt1, Modifs$`Full name`)]
if ("Prot.Quant.Mod.Excl" %in% colnames(Param)) {
  Mod4Quant <- Mod4Quant[which(!Mod4Quant %in% unlist(strsplit(Param$Prot.Quant.Mod.Excl, ";")))]
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
if (("QuantMeth" %in% colnames(Param))&&(!"Quant_algorithm" %in% colnames(Param))) { # Old parameter name
  Param$Quant_algorithm <- Param$QuantMeth
}
if ((!validCharPar("quantAlgo", quantAlgoOpt))&&("Quant_algorithm" %in% colnames(Param))) {
  tmp1 <- Param$Quant_algorithm
  if (validCharPar("tmp1", quantAlgoOpt)) { quantAlgo <- tmp1 }
}
if (!validCharPar("quantAlgo", quantAlgoOpt)) {
  quantAlgo <- c("LM", "limpa")[match(scrptType, c("noReps", "withReps"))]
}
quantAlgo %<o% quantAlgo
Param$Quant_algorithm <- quantAlgo
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
if ((!validCharPar("reScAlgo", reScAlgoOpt))&&("ReScaling_algorithm" %in% colnames(Param))) {
  tmp1 <- Param$ReScaling_algorithm
  if (validCharPar("tmp1", reScAlgoOpt)) { reScAlgo <- tmp1 }
}
if (!validCharPar("reScAlgo", reScAlgoOpt)) {
  if (quantAlgo %in% reScAlgoOpt) { reScAlgo <- quantAlgo }
  if (quantAlgo == "LM") { reScAlgo <- "max" }
}
reScAlgo %<o% reScAlgo
Param$ReScaling_algorithm <- reScAlgo
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
if ((!validIntegPar("topN"))&&("topN" %in% colnames(Param))) {
  tmp1 <- Param$topN
  if (validIntegPar("tmp1")) { topN <- as.integer(tmp1) }
}
if (!validIntegPar("topN")) { topN <- 3 } # Usually top3 
topN %<o% topN
Param$topN <- topN
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
if ((!validLogicPar("topN_correct"))&&("topN_correct" %in% colnames(Param))) {
  tmp1 <- as.logical(Param$topN_correct)
  if (validLogicPar("tmp1")) { topN_correct <- tmp1 }
}
if (!validLogicPar("topN_correct")) { topN_correct <- TRUE }
topN_correct %<o% topN_correct
Param$topN_correct <- topN_correct
#
# Check peptide-to-protein matches
if ((!validLogicPar("Update_Prot_matches"))&&("Update_Prot_matches" %in% colnames(Param))) {
  tmp1 <- as.logical(Param$Update_Prot_matches)
  if (validLogicPar("tmp1")) { Update_Prot_matches <- tmp1 }
}
if (!validLogicPar("Update_Prot_matches")) { Update_Prot_matches <- TRUE }
Update_Prot_matches %<o% Update_Prot_matches # See https://github.com/vdemichev/DiaNN/discussions/1631
Param$Update_Prot_matches <- Update_Prot_matches
#
if ((!validLogicPar("Impute"))&&("Pep.Impute" %in% colnames(Param))) {
  tmp1 <- as.logical(Param$Pep.Impute)
  if (validLogicPar("tmp1")) { Impute <- tmp1 }
}
if (!validLogicPar("Impute")) { Impute <- FALSE }
Impute %<o% Impute
Param$Pep.Impute <- Impute
#
if ((!validIntegPar("PepFoundInAtLeast", 1))&&("PepFoundInAtLeast" %in% colnames(Param))) {
  tmp1 <- Param$PepFoundInAtLeast
  if (validIntegPar("tmp1", 0)) { PepFoundInAtLeast <- as.integer(tmp1) }
}
if (!validIntegPar("PepFoundInAtLeast", 0)) {
  PepFoundInAtLeast <- 1 # ... which means if we have 3 unique peptides we will not use any non-unique ones
}
PepFoundInAtLeast %<o% PepFoundInAtLeast
Param$PepFoundInAtLeast <- PepFoundInAtLeast
#
Exp.map$Replicate <- as.integer(Exp.map$Replicate)
mxRp <- max(as.integer(Exp.map$Replicate), na.rm = TRUE)
mxN <- max(c(2, mxRp-1))
if ((!validIntegPar("PepFoundInAtLeastGrp", 1))&&("PepFoundInAtLeastGrp" %in% colnames(Param))) {
  tmp1 <- Param$PepFoundInAtLeastGrp
  if (validIntegPar("tmp1", 0)) { PepFoundInAtLeastGrp <- as.integer(tmp1) }
}
if (!validIntegPar("PepFoundInAtLeastGrp", 0)) {
  PepFoundInAtLeastGrp <- mxN # ... which means if we have 3 unique peptides we will not use any non-unique ones
}
PepFoundInAtLeastGrp %<o% PepFoundInAtLeastGrp
Param$PepFoundInAtLeastGrp <- PepFoundInAtLeastGrp
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
if ("CytoScapePath" %in% colnames(Param)) {
  Param$CytoscapePath <- Param$CytoScapePath
  Param$CytoScapePath <- NULL
}
if ("CytoscapePath" %in% colnames(Param)) {
  tmp <- normalizePath(Param$CytoscapePath, winslash = "/")
  if ((length(CytoScExe))&&(tmp %in% CytoScExe)) { CytoScExe <- tmp }
}
if (CytoScape&&("Cytoscape" %in% colnames(Param))) {
  tmp1 <- as.logical(Param$Cytoscape)
  if (validLogicPar("tmp1")&&(!CytoScape)) {
    msg <- "Sorry, can't run Cytoscape: couldn't find a valid executable!"
    ReportCalls <- AddMsg2Report(Warning = TRUE)
  } else { CytoScape <- tmp1 }
}
Param$Cytoscape <- CytoScape
#
# Proteome ruler
if ((!validLogicPar("protrul"))&&("ProtRul" %in% colnames(Param))) {
  tmp1 <- as.logical(Param$ProtRul)
  if (validLogicPar("tmp1")) { protrul <- tmp1 }
}
if (!validLogicPar("protrul")) {
  # Proteome ruler re-scaling makes no sense if we are dealing with a pull down
  protrul %<o% (!WorkFlow %in% c("PULLDOWN", "BIOID"))
  # Archaea and Eukaryotes have introns and histones, Bacteria do not
  protrul <- c(protrul, FALSE)[(!isEukaLike)+1]
}
protrul %<o% protrul
Param$ProtRul <- protrul
if ((!validIntegPar("ProtRulNuclL", 1))&&("ProtRulNuclL" %in% colnames(Param))) {
  tmp1 <- Param$ProtRulNuclL
  if (validIntegPar("tmp1")) { ProtRulNuclL <- as.integer(tmp1) }
}
if (!validIntegPar("ProtRulNuclL")) { ProtRulNuclL <- 196 } # Usually top3 
ProtRulNuclL %<o% ProtRulNuclL
Param$ProtRulNuclL <- ProtRulNuclL
#
threshMsg %<o% paste0("% of ", c("control-to-control", "intra-sample group")[match(RefRat_Mode, c("1", "2"))], " ratios")
threshOpt %<o% c(threshMsg, "Absolute log2 FC threshold")
threshDflt <- threshOpt[1]
if (!"Ratios.Thresholds" %in% colnames(Param)) {
  Param$Ratios.Thresholds <- "% of intra-sample group ratios"
}
if (!Param$Ratios.Thresholds %in% threshOpt) {
  Param$Ratios.Thresholds <- "% of intra-sample group ratios"
}
if (Param$Ratios.Thresholds ==  "Absolute threshold") {
  Param$Ratios.Thresholds <- "Absolute log2 FC threshold"
}
threshDflt <- Param$Ratios.Thresholds
Mitch <- match(threshDflt, threshOpt)
if ("Ratios.Contamination.Rates" %in% colnames(Param)) {
  KontRt <- Param$Ratios.Contamination.Rates
  if ((!is.numeric(KontRt))||(KontRt < 0)) { KontRt <- c(0.05, 1)[Mitch] }
} else { KontRt <- c(0.05, 1)[Mitch] }
KontRt <- KontRt*c(100, 1)[Mitch]
KontGrps <- c("Ratio groups", "Experiments", "Whole dataset")
if ((!"Ratios.Contaminant.Groups" %in% colnames(Param))||(!Param$Ratios.Contaminant.Groups %in% KontGrps)) {
  Param$Ratios.Contaminant.Groups <- KontGrps[1]
}
#
if (!exists("minInt")) { minInt <- 0 }
if ("Min.Intensity" %in% colnames(Param)) {
  tmp <- suppressWarnings(as.numeric(Param$Min.Intensity))
  if ((is.numeric(tmp))&&(is.finite(tmp))&&(tmp >= 0)) { minInt <- tmp }
}
minInt %<o% minInt
#
BH.FDR %<o% {
  if ("BH.FDR.values" %in% colnames(Param)) {
    sort(as.numeric(unlist(unique(strsplit(as.character(Param$BH.FDR.values), ";")))))
  } else { c(0.1, 0.2, 0.3) }
}
# For PTMs normalisation to the parent PG, when no quantitation is available for the parent,
# how are we replacing NAs? 
NAsReplMeth %<o% 2
if ("PTM.analysis_NAsReplaceMethod" %in% colnames(Param)) {
  NAsReplMeth %<o% as.integer(Param$PTM.analysis_NAsReplaceMethod)
  if (!NAsReplMeth %in% 1:2) { NAsReplMeth <- 2 }
}
NAsReplMethods <- c("Impute", # Currently not recommended, until I can work out a much less random Imputation method
                    "Median" # Recommended method
)
#
TwoSidedDeflt <- c("Both directions", "Up-only")[(WorkFlow %in% c("PULLDOWN", "BIOID"))+1]
if (("Two.sided" %in% colnames(Param))&&(!is.na(Param$Two.sided))&&(is.logical(Param$Two.sided))) {
  TwoSidedDeflt <- Param$Two.sided
}
Mirror.Ratios %<o% FALSE
#if ("Mirror.Ratios" %in% colnames(Param)) { Mirror.Ratios <- Param$Mirror.Ratios <- as.logical(Param$Mirror.Ratios) }
Param$Mirror.Ratios <- Mirror.Ratios
TwoSidedDeflt <- c(c("Up-only", "Down-only")[Mirror.Ratios+1], "Both directions")[TwoSidedDeflt+1]
# ROC
ROC_GOterms %<o% c()
ROCfilt_GOterms_Pos %<o% c()
ROCfilt_GOterms_Neg %<o% c()
annotRep <- ((Annotate)&&(scrptType == "withReps"))
if (annotRep) {
  if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) { loadFun("GO_terms.RData") }
  if ("ROC.GO.terms" %in% colnames(Param)) {
    Param$ROC_GOterms <- Param$ROC.GO.terms
    Param$ROC.GO.terms <- NULL
  }
  if ("ROC_GOterms" %in% colnames(Param)) { 
    tmp <- unlist(strsplit(Param$ROC_GOterms, ";"))
    tmp <- tmp[which(tmp %in% GO_terms$ID)]
    if (length(tmp)) { ROC_GOterms <- tmp }
  }
  if ("ROCfilt_GOterms_Pos" %in% colnames(Param)) {
    tmp <- unlist(strsplit(Param$ROCfilt_GOterms_Pos, ";"))
    tmp <- tmp[which(tmp %in% GO_terms$ID)]
    if (length(tmp)) { ROCfilt_GOterms_Pos <- tmp }
  }
  if ("ROCfilt_GOterms_Neg" %in% colnames(Param)) {
    tmp <- unlist(strsplit(Param$ROCfilt_GOterms_Neg, ";"))
    tmp <- tmp[which(tmp %in% GO_terms$ID)]
    if (length(tmp)) { ROCfilt_GOterms_Neg <- tmp }
  }
}
if (length(ROC_GOterms)) {
  ROC2_GOterms <- allGO[match(ROC_GOterms, allGO2)]
  w <- c(which(allGO %in% ROC2_GOterms),
         which(!allGO %in% ROC2_GOterms))
  ROC2_allGO1 <- allGO[w]
  ROC2_allGO2 <- allGO2[w]
} else {
  ROC2_allGO1 <- ROC2_allGO2 <- ROC2_GOterms <- c()
}
if (length(ROCfilt_GOterms_Pos)) {
  ROC1_GOterms_Pos <- allGO[match(ROCfilt_GOterms_Pos, allGO2)]
  w <- c(which(allGO %in% ROC1_GOterms_Pos),
         which(!allGO %in% ROC1_GOterms_Pos))
  ROC1_allGOPos1 <- allGO[w]
  ROC1_allGOPos2 <- allGO2[w]
} else {
  ROC1_allGOPos1 <- ROC1_allGOPos2 <- ROC1_GOterms_Pos <- c()
}
if (length(ROCfilt_GOterms_Neg)) {
  ROC1_GOterms_Neg <- allGO[match(ROCfilt_GOterms_Neg, allGO2)]
  w <- c(which(allGO %in% ROC1_GOterms_Neg),
         which(!allGO %in% ROC1_GOterms_Neg))
  ROC1_allGONeg1 <- allGO[w]
  ROC1_allGONeg2 <- allGO2[w]
} else {
  ROC1_allGONeg1 <- ROC1_allGONeg2 <- ROC1_GOterms_Neg <- c()
}
#
if (!exists("normDat")) {
  normDat <- sum(Param$Norma.Ev.Intens,
                 Param$Norma.Pep.Intens,
                 Param$Norma.Prot.Ratio) > 0
}
normDat %<o% normDat
# - Peptide normalisation methods
#   Currently not implemented:
#      - Quantile: what's the point of being quantitative if we are going to replace measurements with quantiles? I do not see for now a rationale to implement this.
#      - RUV normalization (could be interesting, see https://www.bioconductor.org/packages/release/bioc/vignettes/RUVnormalize/inst/doc/RUVnormalize.pdf)
pepNormMethods %<o% list(list(Method = "median",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { median(x, na.rm = TRUE) }"),
                         list(Method = "mean",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { mean(x, na.rm = TRUE) }"),
                         list(Method = "Levenberg-Marquardt",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { median(x, na.rm = TRUE) }" # used for the surrounding steps
                         ),
                         list(Method = "sum",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { log10(sum(10^x, na.rm = TRUE)) }"),
                         list(Method = "logSum",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { median(x, na.rm = TRUE) }"),
                         list(Method = "max",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { max(x, na.rm = TRUE) }"),
                         list(Method = "mode",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { modeest::mlv(is.all.good(x), method = \"Parzen\") }"),
                         list(Method = "proteins",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { mean(x, na.rm = TRUE) }", # used for the surrounding steps
                              Proteins = NA),
                         list(Method = "biotinylated proteins",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { mean(x, na.rm = TRUE) }" # used for the surrounding steps
                         ),
                         list(Method = "GO terms",
                              Source = "pepNorm_General.R",
                              funCall = "normFun <- function(x) { mean(x, na.rm = TRUE) }", # used for the surrounding steps
                              Terms = NA),
                         list(Method = "LOESS",
                              Source = "pepNorm_Shape.R"),
                         list(Method = "VSN",
                              Source = "pepNorm_Shape.R"),
                         list(Method = "IRS",
                              Source = "pepNorm_IRS.R"),
                         list(Method = "ComBat",
                              Source = "pepNorm_ComBat.R",
                              Batch = NA
                         ))
L <- length(pepNormMethods)
pepNormMethodsDF <- data.frame(Method = vapply(1:L, function(i) { pepNormMethods[[i]]$Method }, ""),
                               Source = vapply(1:L, function(i) { pepNormMethods[[i]]$Source }, ""))
# General normalisation methods
genPepNormMeth <- pepNormMethodsDF$Method[which(pepNormMethodsDF$Source == "pepNorm_General.R")]
# genPepNormMeth <- genPepNormMeth[which(!genPepNormMeth %in% c("proteins", "GO terms"))] # These have their own box // NO THEY DON'T, NOT AS STEPS!
if (!IsBioID) {
  genPepNormMeth <- genPepNormMeth[which(genPepNormMeth != "biotinylated proteins")]
}
# "Shape" (variance stabilisation) normalisation methods
shapePepNormMeth <- pepNormMethodsDF$Method[which(pepNormMethodsDF$Source == "pepNorm_Shape.R")]
# Default normalisation sequence
dfltNormSeq <- list(list(Method = "median"),
                    list(Method = "IRS"),
                    list(Method = "ComBat",
                         Batch = "Replicate"),
                    list(Method = "ComBat",
                         Batch = "Isobaric.set"),
                    list(Method = "Levenberg-Marquardt"))
if (LabelType == "Isobaric") {
  if (length(Iso) == 1) {
    dfltNormSeq <- dfltNormSeq[which(vapply(dfltNormSeq, function(x) { !((x$Method == "ComBat")&(x$Batch = "Isobaric.set")) }, TRUE))]
  }
} else {
  dfltNormSeq <- dfltNormSeq[which(vapply(dfltNormSeq, function(x) { x$Method }, "") != "IRS")]
  dfltNormSeq <- dfltNormSeq[which(vapply(dfltNormSeq, function(x) {
    if ("Batch" %in% names(x)) { return(x$Batch) }
    return("")
  }, "") != "Isobaric.set")]
}
if (!Nested) {
  dfltNormSeq <- dfltNormSeq[which(vapply(dfltNormSeq, function(x) { x$Method }, "") != "ComBat")]
}
if (exists("pepNormSeq")) { dfltNormSeq <- pepNormSeq }
#
# 2 functions to convert between different normalisation sequence formats:
normSeqProc12 <- function(seq) { #seq <- dfltNormSeq
  l <- length(seq)
  stopifnot(l > 0)
  vapply(1:l, function(i) {
    x <- seq[[i]]$Method
    if ("Batch" %in% names(seq[[i]])) {
      x <- paste0(x, ": ", paste(seq[[i]]$Batch, collapse = ";"))
    }
    return(paste0(i, " - ", x))
  }, "")
}
normSeqProc21 <- function(seq2) { #seq2 <- dfltNormSeq2
  l <- length(seq2)
  stopifnot(l > 0)
  seq2 <- gsub("^[0-9]+ - ", "", seq2)
  dict <- vapply(pepNormMethods, function(i) { i$Method }, "")
  lapply(1:l, function(i) {
    x <- unlist(strsplit(seq2[[i]], ": "))
    rs <- list(Method = x[[1]])
    m <- match(x[[1]], dict)
    nms <- names(pepNormMethods[[m]])
    nms <- nms[which(nms != "Method")]
    if (length(nms)) {
      rs[nms] <- pepNormMethods[[m]][nms]
    }
    if (length(x) == 2) {
      rs$Batch <- unlist(strsplit(x[[2]], ";"))
    }
    return(rs)
  })
}
if (!exists("normSequence")) { normSequence <- dfltNormSeq }
normSequence %<o% normSequence
dfltNormSeq2 <- normSeqProc12(normSequence)
#
# Some default parameters which nicely follow the same structure
# (NB: Pepper isn't yet integrated in the workflows with replicates)
myPar <- c("Pepper", "GSEA", "ProfPlots", "RankAbundPlots", "ClueGO", "WGCNA")
for (parI in myPar) {
  parNm <- paste0("run", parI)
  # Lowest level default: defined by context
  par_dflt <- par_dflt2 <- c(FALSE,
                             TRUE, # Should be (Annotate)|(Org %in% ... ), cf. GSEA source: there are 20-ish organisms for which a specific annotations package is usable
                             TRUE,
                             TRUE,
                             FALSE,
                             nrow(Exp.map >= 15))[match(parI, myPar)]
  
  # Level 1 default: defined by Param
  if (parNm %in% colnames(Param)) {
    par_dflt <- as.logical(Param[[parNm]])
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
mnFct <- c("Experiment", "Replicate")
if (WorkFlow == "TIMECOURSE") { mnFct <- c(mnFct, "Time.point") }
l <- length(mnFct)
coreNms <- setNames(c("Ratios.Groups.Ref.Aggregate.Level",
                      "Norm.Groups",
                      "Ratios.Groups",
                      "GO.enrichment.Ref.Aggr"
                      #, "Batch.correction"
),
c(paste0("Individual Sample___Combination of Factors required to distinguish individual samples./nMust include ",
         paste(mnFct[1:(l-1)], collapse = ", "), " and ", mnFct[l], "."),
  "Normalisation groups___Factor(s) defining groups of samples to normalize to each-other. Note that which, if any, normalisations apply will be defined further down",
  "Ratio groups___Factor(s) defining comparison groups, i.e. groups of samples, including at least some References (i.e. Controls), to compare to each others./nMust include Experiment, cannot include Replicate.",
  "GO enrichment___Optional: Only protein groups with at least one valid value in corresponding samples will be used as references for GO-terms enrichment tests."
  #, "Batch___Optional: Factor(s) defining batches used for sva::ComBat-based correction."
))
for (nm in coreNms) { if (!nm %in% colnames(Param)) { Param[[nm]] <- "Exp" } }
wMp <- c(which(colnames(Param) == coreNms[1]),
         which(colnames(Param) == coreNms[2]),
         which(colnames(Param) == coreNms[3]),
         which(colnames(Param) == coreNms[4]),
         #which(colnames(Param) == coreNms[5]),
         which((Param[1,] == "MAP2FACTS")&(!colnames(Param) %in% c(coreNms, "Batch.correction"))))
lstFct <- list()
dfltFct <- list()
for (w in wMp) { #w <- wMp[3]
  myFct <- colnames(Param)[w]
  Opt <- Factors
  dflt <- dfdflt <- c("Experiment", "Replicate")[(myFct == "Batch.correction")+1]
  if (!Param[1, w] %in% c("AUTOFACT", "MAP2FACTS")) {
    dflt <- Factors[match(unlist(strsplit(Param[1, w], ";")), names(Factors))]
    if ((length(dflt) == 1)&&(is.na(dflt))) { dflt <- dfdflt }
  }
  if (myFct == "Ratios.Groups.Ref.Aggregate.Level") {
    Opt <- Factors
    dflt <- unique(c(dflt, mnFct))
    dflt <- dflt[which(!dflt %in% c("Batch"))]
  }
  if (myFct %in% c("Ratios.Groups", "GO.enrichment.Ref.Aggr")) {
    Opt <- Factors[which(Factors != "Replicate")]
    dflt <- dflt[which(!dflt %in% c("Batch", "Replicate"))]
  }
  lbl <- gsub("\\.", " ", myFct)
  if (myFct %in% coreNms) {
    lbl <- unlist(strsplit(names(coreNms)[match(myFct, coreNms)], "___"))
  } else {
    lbl <- gsub("\\.", " ", myFct)
    if (myFct %in% names(Param_Help)) { lbl <- c(lbl, Param_Help[myFct]) }
  }
  names(Opt) <- NULL
  names(dflt) <- NULL
  dfltFct[[myFct]] <- dflt
  blck <- list(list(br()),
               tags$table(
                 tags$tr(width = "80%",
                         tags$td(width = "25%",
                                 div(strong(lbl[1]))),
                         tags$td(width = "55%",
                                 selectInput(myFct,
                                             "",
                                             Opt,
                                             dflt,
                                             TRUE,
                                             TRUE)),
                         #addTooltip(session, myFct, Param_Help[w], "bottom", "hover", list(container = "body"))
                 )
               ))
  lbl2 <- unlist(strsplit(lbl[2], "/n"))
  for (lbl2a in lbl2) { blck <- append(blck, list(span(em(lbl2a)), br())) }
  if (myFct == "Ratios.Groups") {
    nstDflt <- Param$Ratios.Groups_Nested
    if (!is.logical(nstDflt)) { nstDflt <- WorkFlow != "Regulation" }
    blck <- append(blck, list(radioButtons("IsNested", "Nested design? (i.e. are replicates paired?)",
                                           c(TRUE, FALSE), nstDflt, TRUE)))
  }
  blck <- append(blck, list(br()))
  lstFct <- append(lstFct, blck)
}
#
useSAM_thresh %<o% FALSE
tstAdvOpt <- try(sum(file.exists(Param$Custom.PGs, Param$TrueDisc_filter, Param$CRAPome_file)) > 0)
if ("try-error" %in% class(tstAdvOpt)) { tstAdvOpt <- FALSE }
#
mtchCheckMsg1 <- "Not all search software will map peptides to protein IDs in the search database the same way. Using this function ensures consistent results regardless of search engine."
mtchCheckMsg2 <- "!Checking assignments may result in removal of some identifications!"
F_test_override <- FALSE
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
    extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
    titlePanel(tag("u", "Parameters"),
               #em(appNm)), # Doesn't work =(
               appNm),
    h2(dtstNm),
    tags$hr(style = "border-color: black;"),
    br(),
    h4(strong("Factors")),
    span("Here we map specific actions to experimental Factors:"),
    fluidRow(column(4,
                    uiOutput("RSA_msg"),
                    withSpinner(uiOutput("FactMappings"))),
             column(8, withSpinner(plotlyOutput("PSMsPCA", height = "600px")))),
    br(),
    tags$hr(style = "border-color: black;"),
    h4(strong("Proteins of interest")),
    pickerInput("IntProt",
                NULL,
                protHeads,
                protDeflt,
                TRUE,
                pickerOptions(title = "Search me",
                              `live-search` = TRUE,
                              actionsBox = TRUE,
                              deselectAllText = "Clear search",
                              showTick = TRUE),
                width = "600px"),
    br(),
    tags$hr(style = "border-color: black;"),
    h4(strong("Data processing")),
    fluidRow(
      column(3,
             numericInput("minInt",
                          "Exclude PSMs with intensity lower than...",
                          minInt,
                          0,
                          .Machine$double.xmax,
                          width = "100%")),
      # Imputation: add a new parameter here to define the counts limit between MCAR and MNAR 
      column(3,
             checkboxInput("Impute", "Impute missing peptides-level values?", Impute, "100%"),
             em("Imputation is done at peptides (not PSMs or protein groups) level, after removing any outliers, before re-normalization."),
             em("This is independent of any temporary imputations done by default to run specific procedures without data loss."),
             em(HTML("&nsbp;- MAR/MCAR data is imputed with a KNN method.")),
             em(HTML("&nsbp;- MNAR data is imputed with the QRILC method."))),
      column(3,
             checkboxInput("Update_Prot_matches", paste0("Update peptide-to-protein assignments?"), Update_Prot_matches, "100%"),
             bsTooltip("Update_Prot_matches",
                       paste0(mtchCheckMsg1, "\n", mtchCheckMsg2),
                       placement = "right", trigger = "hover", options = list(container = "body")),
             h5(em(mtchCheckMsg1)),
             h5(em(mtchCheckMsg2))#, withSpinner(uiOutput("ReloadMatches"))
      )
    ),
    br(),
    if (annotRep) {
      fluidRow(column(6,
                      checkboxInput("ROC1on", "PSMs-level ROC filter", length(ROC_GOterms) > 0, "100%"),
                      withSpinner(uiOutput("ROC1"))))
    },
    tags$hr(style = "border-color: black;"),
    withSpinner(uiOutput("IsobarCorr")),
    h4(strong("Normalisation")),
    withSpinner(uiOutput("Norm")),
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
                                 " -> the whole dataset?", PepFoundInAtLeast, 1, nSmpls, 1, "100%"),
                    numericInput("PepFoundInAtLeastGrp",
                                 " -> one sample group at least? (overrides parameter above)",
                                 PepFoundInAtLeastGrp,
                                 1,
                                 mxRp,
                                 1,
                                 "100%")),
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
    tags$hr(style = "border-color: black;"),
    fluidRow(column(2, radioButtons("Clustering",
                                    "Heatmaps: clustering method",
                                    klustChoices,
                                    klustChoices[1], TRUE, "100%"))),
    br(),
    tags$hr(style = "border-color: black;"),
    h4(strong("Statistical tests")),
    h5(strong(" -> t-test(s)")),
    fluidRow(
      column(2,
             radioButtons("TwoSided", "Fold changes: test...", c("Both directions", "Up-only", "Down-only"),
                          TwoSidedDeflt, FALSE, "100%"),
             # checkboxInput("Mirror", "Revert fold changes on plots? (default: log fold change = log2(Sample/Reference); revert: log2(Reference/Sample))",
             #               Param$Mirror.Ratios, "100%")
      ),
      column(4,
             strong("Volcano plot: select default variant"),
             radioButtons("TtstPval", "", names(pvalue.col), "Moderated", TRUE, "100%"),
             h6(em(" - Welch's t-test is a modified form of Student's original version which is more robust to variance inequality.")),
             h6(em(" - Moderated t-test (limma): re-samples individual row variances using global dataset variance to provide a more robust estimate.")),
             h6(em(" - Permutation t-test (coin): based on permutations of the samples from each group.")),
             h6(em(" - SAM's modified t-test (siggenes): corrects for poor variance estimates for low replicate number data by optimizing a small constant s0 added to the denominator of the test statistic.")),
             h6(em(" - LRT (Likelihood Ratio Test //edge): tests for the ratio of the likelihoods of observed data under two models (explanatory variable has an effect /vs/ no effect)")),
             h6(em(" - ODP (Optimal Discovery Procedure, Storey et al., 2007 //edge): uses all relevant information from all genes in order to test each one for differential expression; has been demonstrated to have optimal power.")),
             if (scrptTypeFull == "withReps_PG_and_PTMs") {
               # Not available for now for the peptides-only workflow,
               # but this could be done by turning the code for this into an app which is called individually for each PTM class.
               # This would mean allowing for a different selection for each PTM class.
               # We will anyway want to allow for multiple tests selected in the future, so go out of the "choose one t-test" variant approach...
               h6(em("(you can change which one will be used for Volcano plots later, after we compare each test's power)"))
             },
             checkboxInput("useSAM_thresh", "For Student's t-test, plot SAM-based curved significance thresholds?", useSAM_thresh, "100%")),
      column(2,
             h5("Benjamini-Hochberg FDR thresholds"),
             withSpinner(uiOutput("FDR")),
             numericInput("FDR", "", 0.01, 0, 1),
             actionBttn("AddFDR", "Add threshold", color = "primary", size = "xs", style = "pill"),
             actionBttn("RemvFDR", "Remove threshold", color = "primary", size = "xs", style = "pill")),
      column(2,
             radioButtons("typeOfThresh", "Ratios thresholds: use...", threshOpt, threshDflt, FALSE, "100%"),
             numericInput("RatCont", "Threshold value = ", KontRt, 0, 100, 0.001, width = "100%"),
             selectInput("RatContGrp",
                         "Control ratios are grouped by:",
                         KontGrps,
                         Param$Ratios.Contaminant.Groups,
                         width = "100%"))),
    br(),
    fluidRow(column(2,
                    h5(strong(" -> ANOVA (moderated F-test //limma)")),
                    checkboxInput("run_F_test", "Run?", fTstDflt, "100%"),
                    em("(only makes sense if N(sample groups) > 2)")),
             column(2,
                    uiOutput("sntXprs")),
             if (annotRep) {
               column(3,
                      checkboxInput("ROC2on", "ROC analysis of P-values (experimental)", length(ROC_GOterms) > 0, "100%"),
                      withSpinner(uiOutput("ROC2")))
             }
    ),
    tags$hr(style = "border-color: black;"),
    withSpinner(uiOutput("GO")),
    if (scrptTypeFull == "withReps_PG_and_PTMs") {
      list(
        tags$hr(style = "border-color: black;"),
        h4(strong("Optional analyses")),
        fluidRow(column(2,
                        checkboxInput("runProfPlots", "Draw protein profile plots?",
                                      runProfPlots, "100%"),
                        checkboxInput("runRankAbundPlots", "Draw protein ranked abundance plots?",
                                      runRankAbundPlots, "100%")),
                 column(2,
                        checkboxInput("runWGCNA", "Run Weighted Gene Correlation Network Analysis (WGCNA)?",
                                      runWGCNA, "100%"),
                        checkboxInput("runGSEA", "Run Gene Set Enrichment Analysis (GSEA)?",
                                      runGSEA, "100%")))
      )
    },
    tags$hr(style = "border-color: black;"),
    withSpinner(uiOutput("CytoScape")),
    tags$hr(style = "border-color: black;"),
    h4(strong("Post-translational modifications (PTMs)")),
    fluidRow(column(2, pickerInput("PTMsStats",
                                   "Select PTM(s) (if any) for which statistical tests will be performed and subtables written:",
                                   Modifs$`Full name`[which(Modifs$Type == "Variable")],
                                   unlist(strsplit(ptmDflt2, ";")),
                                   TRUE,
                                   pickerOptions(title = "Search me",
                                                 `live-search` = TRUE,
                                                 actionsBox = TRUE,
                                                 deselectAllText = "Clear search",
                                                 showTick = TRUE))),
             column(2, checkboxInput("PTMsReNorm",
                                     "Re-normalize modified peptides ratios to those of parent Protein Group(s)?",
                                     Param$PTM.analysis_Norm, "100%")),
             column(2, radioButtons("NAsReplMethod",
                                    "Some modified peptides do not have a quantified parent protein group to re-normalize to. Replace missing values using:",
                                    NAsReplMethods,
                                    NAsReplMethods[NAsReplMeth],
                                    width = "100%"))
    ),
    h4(strong("Output tables")),
    fluidRow(column(3, checkboxInput("Amica", "Write tables compatible with https://bioapps.maxperutzlabs.ac.at/app/amica", Param$Amica, "100%"))),
    br(),
    tags$hr(style = "border-color: black;"),
    checkboxInput("AdvOptOn", "Advanced options", tstAdvOpt),
    withSpinner(uiOutput("AdvOpt")),
    br(),
    actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
    br(),
    br()
    # ,
    # setBackgroundColor( # Doesn't work
    #    color = c("#F8F8FF", "#91E0E3"),
    #    gradient = "linear",
    #    direction = "bottom"
    # )
  )
} 
if (exists("appRunTest")) { rm(appRunTest) }
server1 <- function(input, output, session) {
  NORMALIZE <- reactiveVal(normDat)
  NORMSEQ <- reactiveVal(dfltNormSeq2)
  PURFL <- reactiveVal(Param$Label.Purities.file)
  RSA_Msg <- reactiveVal("")
  F_TEST <- reactiveVal(fTstDflt)
  #
  # Server sub-functions
  # (unfortunately, it does not look like these can be pre-created outside the server - the reactive objects do not work if they are)
  updtNorm <- function(reactive = TRUE) {
    if (reactive) {
      nrmSq2 <- NORMSEQ()
      nrm <- NORMALIZE()
    } else {
      nrmSq2 <- dfltNormSeq2
      nrm <- normDat
    }
    nrmSq <- normSeqProc21(nrmSq2)
    if (nrm) {
      lst <- list(list(
        #
        # Parameters shared by multiple normalisations
        h5("Re-normalize to ..."),
        fluidRow(column(3,
                        pickerInput("Norma.Prot.Ratio.to.proteins",
                                    " -> select proteins",
                                    protHeads,
                                    nrm2PrtSlctDflt,
                                    TRUE,
                                    pickerOptions(title = "Search me",
                                                  `live-search` = TRUE,
                                                  actionsBox = TRUE,
                                                  deselectAllText = "Clear search",
                                                  showTick = TRUE))),
                 if (Annotate) {
                   column(3,
                          pickerInput("Norma.Prot.Ratio.to.GO",
                                      " -> select GO terms",
                                      nrm2GOall,
                                      nrm2GO,
                                      TRUE,
                                      pickerOptions(title = "Search me",
                                                    `live-search` = TRUE,
                                                    actionsBox = TRUE,
                                                    deselectAllText = "Clear search",
                                                    showTick = TRUE)))
                 },
                 if (WorkFlow == "BIOID") {
                   column(3,
                          checkboxInput("Norma.Prot.Ratio.to.Biot", " -> all biotinylated proteins: ",
                                        Param$Norma.Prot.Ratio.to.Biot, "100%"))
                 },
        ),
        br(),
        #
        # PSMs normalisations
        fluidRow(column(1,
                        strong(" -> PSMs-level:")),
                 column(2,
                        checkboxInput("evLM", "Levenberg-Marquardt",
                                      Param$Adv.Norma.Ev.Intens, "100%"))),
        br(),
        #
        # Peptidoforms
        fluidRow(column(6,
                        selectInput("PepNormSeq",
                                    " -> Peptidoforms-level:",
                                    nrmSq2,
                                    nrmSq2,
                                    TRUE,
                                    TRUE,
                                    width = "90%"))),
        #h5(strong(em(paste(nrmSq2, collapse = " -> ")))),
        fluidRow(column(2,
                        actionBttn("addClassicNorm", "Add classic normalisation", color = "primary", size = "xs", style = "pill"),
                        pickerInput("ClassicNormMeth",
                                    "methods:",
                                    genPepNormMeth,
                                    genPepNormMeth[1],
                                    width = "80%")),
                 column(2,
                        actionBttn("addVarCorr", "Add variance intensity-bias correction", color = "primary", size = "xs", style = "pill"),
                        pickerInput("VarCorrMeth",
                                    "methods:",
                                    shapePepNormMeth,
                                    shapePepNormMeth[1],
                                    width = "80%")),
                 column(2,
                        actionBttn("addComBat", "Add ComBat batch correction", color = "primary", size = "xs", style = "pill"),
                        pickerInput("ComBatBatch",
                                    "batches (select one or more):",
                                    Factors,
                                    c("Replicate", "Batch")[("Batch" %in% Factors)+1],
                                    TRUE,
                                    width = "80%")),
                 if (LabelType == "Isobaric") {
                   column(2,
                          actionBttn("addIRS", "Add IRS batch correction", color = "primary", size = "xs", style = "pill"),
                          h5(em("IRS = Internal Reference Scaling")))
                 },
        ),
        h5(em("(It goes without saying that you should aim for a REASONABLE normalisation scheme!!!)")),
        br()
      ))
      #
      # Protein groups
      if (scrptTypeFull == "withReps_PG_and_PTMs") {
        lst <- append(lst, list(fluidRow(column(2,
                                                h5(strong(" -> Protein Groups-level:"))),
                                         column(2,
                                                checkboxInput("prtLM", "Levenberg-Marquardt",
                                                              Param$Adv.Norma.Prot.Intens, "100%")))
        ))
      }
      lst <- append(lst, list(br()))
    } else {
      lst <- list(list(br()))
    }
    renderUI(lst)
  }
  updtOptOn <- function(reactive = TRUE) {
    tst <- { if (reactive) { ADVOPT() } else { tstAdvOpt } }
    if (tst) {
      rgKol <- { if (reactive) { input[["Ratios.Groups"]] } else { dfltFct[["Ratios.Groups"]] } }
      RatGrps <- paste(unique(apply(ExpMap[, rgKol, drop = FALSE], 1, paste, collapse = " ")), collapse = " / ")
      lst <- vector("list", 1)
      lst[[1]] <- list(
        fluidRow(column(12, em("->"), 
                        shinyFilesButton("CustPG", em("Custom Protein Groups"), "", FALSE), br(),
                        em("Allows \"cheating\" with the naive Protein Groups assembly algorithm."), br(),
                        em("Useful when for instance genetically-motified/transduced samples synthesize a custom protein (e.g. fusion construct) to which matching peptides should be assigned in priority over the canonical form."), br(),
                        em("Should be a table with two columns: \"Leading protein IDs\" (\";\"-separated) and \"Priority\" (integer)."), br(),
                        strong("Current selection = "), span(Param$Custom.PGs, style = "color:blue", .noWS = "outside"),
                        br(),
                        br()
        )),
        fluidRow(column(12, em("->"), 
                        shinyFilesButton("TrueDisc", em("TRUE/FALSE protein discovery filter"), "", FALSE), br(),
                        em("Select TRUE/FALSE protein discovery filter .csv file (useful when you know from another experiment which proteins are contaminants)."), br(),
                        em("Should contain a \"Protein ID\" column and one TRUE/FALSE column per \"Ratio group\" value."),
                        em("Current values: "), em(RatGrps), br(),
                        br(),
                        em(DiscFiltModesHlp[1]), br(),
                        em(DiscFiltModesHlp[2]), br(),
                        em(DiscFiltModesHlp[3]), br(),
                        strong("Current selection = "), span(Param$TrueDisc_filter, style = "color:blue", .noWS = "outside"),
                        br(),
                        br()
        )),
        fluidRow(column(12, em("->"),
                        shinyFilesButton("CRAPome", em("CRAPome filter"), "", FALSE), br(),
                        em("CRAPome-like filter: 1 column table of protein accessions to mark as contaminants."), br(),
                        em("Column name = \"Protein ID\" or \"Protein IDs\", use \";\" if including more than one ID per row)."), br(),
                        strong("Current selection = "), span(Param$CRAPome_file, style = "color:blue", .noWS = "outside"),
                        br(),
                        br()
        ))
      )
    } else { lst <- list(list(em(""))) }
    renderUI(lst)
  }
  if (annotRep) { # Not implemented yet for script without Reps
    updtROC1 <- function(reactive = TRUE) {
      tst <- { if (reactive) { ROC1ON() } else { length(c(ROCfilt_GOterms_Pos, ROCfilt_GOterms_Neg)) > 0 } }
      if (tst) {
        lst <- list(
          fluidRow(column(4,
                          pickerInput("ROC1_pos",
                                      "Positive GO term(s)",
                                      ROC1_allGOPos1,
                                      ROC1_GOterms_Pos,
                                      TRUE,
                                      pickerOptions(title = "Search me",
                                                    `live-search` = TRUE,
                                                    actionsBox = TRUE,
                                                    deselectAllText = "Clear search",
                                                    showTick = TRUE))),
                   column(4,
                          pickerInput("ROC1_neg",
                                      "Negative GO term(s)",
                                      ROC1_allGONeg1,
                                      ROC1_GOterms_Neg,
                                      TRUE,
                                      pickerOptions(title = "Search me",
                                                    `live-search` = TRUE,
                                                    actionsBox = TRUE,
                                                    deselectAllText = "Clear search",
                                                    showTick = TRUE)),
                          em("Optional: these can be used for an inverted ROC analysis."))))
      } else { lst <- list(list(em(""))) }
      renderUI(lst)
    }
    updtROC2 <- function(reactive = TRUE) {
      tst <- { if (reactive) { ROC2ON() } else { length(ROC_GOterms) > 0 } }
      if (tst) {
        lst <- list(
          fluidRow(column(4,
                          em("GO terms to evaluate acceptable P-value significance thresholds, e.g. if you know that all proteins with a given term should be significant. Currently this is only used for personal interpretation of the plots and not fed back into the automated significance analysis."),
                          pickerInput("ROC2_terms",
                                      "",
                                      ROC2_allGO1,
                                      ROC2_GOterms,
                                      TRUE,
                                      pickerOptions(title = "Search me",
                                                    `live-search` = TRUE,
                                                    actionsBox = TRUE,
                                                    deselectAllText = "Clear search",
                                                    showTick = TRUE)))))
      } else { lst <- list(list(em(""))) }
      renderUI(lst)
    }
  }
  updtRSAmsg <- function(reactive = TRUE) {
    myMsg <- { if (reactive) { RSA_Msg() } else { "" } }
    renderUI(h4(strong(em(myMsg,
                          style = paste0("color:", c("black", "red")[(nchar(myMsg) > 0)+1]),
                          .noWS = "outside"))))
  }
  updtQuantMsg <- function(reactive = TRUE) {
    myQuant <- { if (reactive) { input$Quant_algorithm } else { quantAlgo } }
    quantMsg <- allQuantAlgos$Help[match(myQuant, allQuantAlgos$Algorithm)]
    renderUI(em(quantMsg))
  }
  #
  # Initialize variables to create in main environment
  PARAM <- reactiveVal(Param)
  m4Quant <- reactiveVal(Mod4Quant)
  m2Xclud <- reactiveVal(Mod2Xclud)
  ADVOPT <- reactiveVal(tstAdvOpt)
  if (annotRep) {
    ROC1ON <- reactiveVal(length(c(ROCfilt_GOterms_Pos,
                                   ROCfilt_GOterms_Neg)) > 0)
    ROC2ON <- reactiveVal(length(ROC_GOterms) > 0)
  }
  #
  # Dynamic UI
  # Map Parameters to Factors
  output$QuantMsg <- updtQuantMsg(FALSE)
  output$PSMsPCA <- renderPlotly(plot_lyPSMsPCA)
  output$RSA_msg <- updtRSAmsg(FALSE)
  output$FactMappings <- renderUI({ lstFct })
  output$RatioGroups <- renderUI({
    em(paste(unique(do.call(paste, c(ExpMap[, input[["Ratios.Groups"]], drop = FALSE], sep = "___"))),
             collapse = " / "))
  })
  # Factors
  sapply(wMp, function(w) {
    myFct <- colnames(Param)[w]
    observeEvent(input[[myFct]], {
      tmpVal <- input[[myFct]] # Current aggregate
      if (myFct == "Ratios.Groups.Ref.Aggregate.Level") {
        # For ref samples aggregate - the aggregate which defines individual samples /// what a terrible name I chose back then... ///
        # we have to do some extra homework:
        l <- length(tmpVal)
        tst2 <- tst3 <- tst4 <- FALSE
        l2 <- 2
        if (length(tmpVal)) {
          tst2 <- sum(!c("Experiment", "Replicate") %in% tmpVal)
          tmpRSA <- do.call(paste, c(Exp.map[, tmpVal, drop = FALSE], sep = "___"))
          UtmpRSA <- unique(tmpRSA)
          tst3 <- length(UtmpRSA) < nSmpls
          tmpVal2 <- tmpVal[which(tmpVal != "Replicate")] # Defines sample groups...
          if (length(tmpVal2)) {
            xpRws <- 1:nSmpls
            tmpVPAL <- do.call(paste, c(Exp.map[, tmpVal2, drop = FALSE], sep = "___"))
            tst4 <- aggregate(xpRws, list(tmpVPAL, Exp.map$Replicate), length)
            tst4 <- max(tst4$x) > 1
            l2 <- length(unique(tmpVPAL))
          }
        }
        msg <- ""
        if ((l < 3)||(tst2)||(tst3)||(tst4)) {
          msg <- c()
          if ((l < 3)||(tst2)) { msg <- c(msg, "You MUST always include Experiment and Replicate here, as well as at least one other contrasted factor!") }
          if (tst3) { msg <- c(msg, "The chosen experimental Factors do not discriminate fully between some samples (rows in the experiment map), add more!") }
          if (tst4) { msg <- c(msg, "It is not allowed for several samples to have the same replicate number within a group!") }
          msg <- paste(msg, collapse = " / ")
          shinyjs::disable("saveBtn")
        } else {
          shinyjs::enable("saveBtn")
        }
        RSA_Msg(msg)
        output$RSA_msg <- updtRSAmsg()
        #
        dflt <- F_TEST()
        tmpVal <- input$Ratios.Groups.Ref.Aggregate.Level
        if (!length(dflt)) { dflt <- TRUE }
        dflt <- suppressWarnings(as.logical(dflt[1]))
        if (is.na(dflt)) { dflt <- TRUE }
        updateCheckboxInput(session = getDefaultReactiveDomain(),
                            inputId = "run_F_test",
                            label = "Run?",
                            value = fTstDflt)
        if (l2 < 3) {
          shinyjs::disable("run_F_test")
          assign("F_test_override", TRUE, envir = .GlobalEnv)
        } else {
          shinyjs::enable("run_F_test")
          assign("F_test_override", FALSE, envir = .GlobalEnv)
        }
      }
      Par <- PARAM()
      Par[colnames(Par)[w]] <- paste0(tmpVal, collapse = ";")
      PARAM(Par)
    }, ignoreNULL = FALSE)
  })
  # SAINTexpress
  output$sntXprs <- renderUI({
    lst <- list(list(br()))
    if (WorkFlow %in% c("PULLDOWN", "BIOID")) {
      lst <- list(list(h5(strong(" -> SAINTexpress")),
                       checkboxInput("saintExprs", "run?", saintExprs, "100%"))
      )
    }
    return(lst)
  })
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
  # Labels purity correction - only for isobarically labelled samples
  updtIsoPur <- function(reactive = TRUE, lblType = LabelType) {
    if (lblType == "Isobaric") {
      fl <- { if (reactive) { PURFL() } else { Param$Label.Purities.file } }
      lst <- list(list(fluidRow(column(1,
                                       shinyFilesButton("PurityFl", "Browse", paste0(IsobarLab, " purity table"), "", FALSE)),
                                column(3, em(fl))),
                       br()))
      
    } else {
      lst <- list(list(HTML("")))
    }
    renderUI(lst)
  }
  output$IsobarCorr <- updtIsoPur(FALSE)
  #
  # Normalisations
  output$Norm <- updtNorm(FALSE)
  #
  # GO
  output$GO <- renderUI({
    lst <- list(list(br()))
    if (Annotate) {
      lst <- list(
        list(h4(strong("GO terms enrichment")),
             fluidRow(column(1,
                             checkboxInput("GOenrich", "GO enrichment", goDflt, "100%"),
                             checkboxInput("runClueGO", "run ClueGO enrichment (NB: this is rather slow)", runClueGO, "100%")),
                      column(2,
                             pickerInput("GO.tabs",
                                         "GO terms of interest",
                                         allGO,
                                         dftlGO,
                                         TRUE,
                                         pickerOptions(title = "Search me",
                                                       `live-search` = TRUE,
                                                       actionsBox = TRUE,
                                                       deselectAllText = "Clear search",
                                                       showTick = TRUE))),
                      column(1,
                             checkboxInput("GO2Int", "Use GO terms to define list of proteins of interest?",
                                           as.character(toupper(Param$GO.terms.for.proteins.of.interest)) == "TRUE", # Clumsy but possibly more backwards proof
                                           "100%"))))
      )
    }
    return(lst)
  })
  # # Pepper
  # observeEvent(input$runPepper, {
  #   runPepper <- as.logical(input$runPepper)
  #   assign("runPepper", runPepper, envir = .GlobalEnv)
  #   Par <- PARAM()
  #   Par$runPepper <- runPepper
  #   PARAM(Par)
  # })
  # Impute?
  observeEvent(input$Impute, {
    Impute <- as.logical(input$Impute)
    assign("Impute", Impute, envir = .GlobalEnv)
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
    KlustMeth <- match(input$Clustering, klustChoices)
    assign("KlustMeth", KlustMeth, envir = .GlobalEnv)
  })
  # Quantitation
  observeEvent(input$Quant_algorithm, {
    quantAlgo <- input$Quant_algorithm
    assign("quantAlgo", quantAlgo, envir = .GlobalEnv)
    if (quantAlgo == "limpa") {
      updateCheckboxInput(inputId = "Impute", value = FALSE)
      updatePickerInput(inputId = "Norma.Prot.Ratio.to.proteins", selected = NULL)
      updatePickerInput(inputId = "Norma.Prot.Ratio.to.GO", selected = NULL)
      updateCheckboxInput(inputId = "Norma.Prot.Ratio.to.Biot", value = FALSE)
      updateRadioButtons(inputId = "TtstPval", selected = "Moderated")
      shinyjs::disable("Impute")
      shinyjs::disable("Norma.Prot.Ratio.to.proteins")
      shinyjs::disable("Norma.Prot.Ratio.to.GO")
      shinyjs::disable("Norma.Prot.Ratio.to.Biot")
      #
      # Where to deal with batch-correction?
      # Not here: it should be encoded in the design matrix in general rather than corrected for, whether using limpa or not for quant.
    } else {
      shinyjs::enable("Impute")
      shinyjs::enable("Norma.Prot.Ratio.to.proteins")
      shinyjs::enable("Norma.Prot.Ratio.to.GO")
      shinyjs::enable("Norma.Prot.Ratio.to.Biot")
    }
    Par <- PARAM()
    Par$Quant_algorithm <- quantAlgo
    PARAM(Par)
    output$QuantMsg <- updtQuantMsg()
  })
  observeEvent(input$ReScaling, {
    reScAlgo <- input$ReScaling
    assign("reScAlgo", reScAlgo, envir = .GlobalEnv)
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
    Pep4QuantOpt[input$Prot.Quant.Use]
    assign("Pep4Quant", Pep4Quant, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Prot.Quant.Use <- input$Prot.Quant.Use
    PARAM(Par)
  })
  #   - PTMs to use for PG Quant
  observeEvent(input$PTMsQuant, {
    ptmDflt1 <- input$PTMsQuant
    assign("ptmDflt1", ptmDflt1, envir = .GlobalEnv)
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
    assign("protrul", protrul, envir = .GlobalEnv)
    if (protrul) { shinyjs::enable("ProtRulNuclL") }
    if (!protrul) { shinyjs::disable("ProtRulNuclL") }
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
  # WGCNA
  observeEvent(input$runWGCNA, {
    runWGCNA <- as.logical(input$runWGCNA)
    assign("runWGCNA", runWGCNA, envir = .GlobalEnv)
    Par <- PARAM()
    Par$runWGCNA <- runWGCNA
    PARAM(Par)
  })
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
  # Are analyses Two-sided?
  observeEvent(input$TwoSided, {
    TwoSided <- input$TwoSided == "Both directions"
    assign("TwoSided", TwoSided, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Two.sided <- TwoSided
    # if (input$TwoSided == "Down-only") {
    #   shinyjs::disable("Mirror")
    #   Par$Mirror.Ratios <- TRUE
    # }
    # if (input$TwoSided != "Down-only") { shinyjs::enable("Mirror") }
    PARAM(Par)
  })
  #
  # Optional input files
  output$AdvOpt <- updtOptOn(FALSE)
  # Event observers
  # Optional input files
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
  observe({ shinyFileChoose(input, "TrueDisc", roots = getVolumes(), filetypes = "csv")
    {
      tmp <- input$TrueDisc
      if ((!is.null(tmp))&&(is.list(tmp))) {
        tmp <- parseFilePaths(getVolumes(), tmp)$datapath
        Par <- PARAM()
        Par$TrueDisc_filter <- normalizePath(tmp, winslash = "/")
        PARAM(Par)
        shinyjs::enable("TrueDiscMode")
      } else { shinyjs::disable("TrueDiscMode") }
  }
  })
  observeEvent(input$TrueDiscMode, {
    Par <- PARAM()
    Par$TrueDisc_filter_mode <- input$TrueDiscMode
    PARAM(Par)
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
  # ROC
  if (annotRep) { # Not implemented yet for script without Reps
    output$ROC1 <- updtROC1(FALSE)
    output$ROC2 <- updtROC1(FALSE)
    observeEvent(input$ROC1on, {
      ROC1ON(input$ROC1on)
      output$ROC1 <- updtROC1()
    })
    observeEvent(input$ROC2on, {
      ROC2ON(input$ROC2on)
      output$ROC2 <- updtROC2()
    })
    observeEvent(input$ROC1_pos, {
      ROCfilt_GOterms_Pos <- ROC1_allGOPos2[match(input$ROC1_pos, ROC1_allGOPos1)]
      assign("ROCfilt_GOterms_Pos", ROCfilt_GOterms_Pos, envir = .GlobalEnv)
      Par <- PARAM()
      Par$ROCfilt_GOterms_Pos <- paste(ROCfilt_GOterms_Pos, collapse = ";")
      PARAM(Par)
    })
    observeEvent(input$ROC1_neg, {
      ROCfilt_GOterms_Neg <- ROC1_allGONeg2[match(input$ROC1_neg, ROC1_allGONeg1)]
      assign("ROCfilt_GOterms_Neg", ROCfilt_GOterms_Neg, envir = .GlobalEnv)
      Par <- PARAM()
      Par$ROCfilt_GOterms_Neg <- paste(ROCfilt_GOterms_Neg, collapse = ";")
      PARAM(Par)
    })
    observeEvent(input$ROC2_terms, {
      ROC_GOterms <- ROC2_allGO2[match(input$ROC2_terms, ROC2_allGO1)]
      assign("ROC_GOterms", ROC_GOterms, envir = .GlobalEnv)
      Par <- PARAM()
      Par$ROC_GOterms <- paste(ROC_GOterms, collapse = ";")
      PARAM(Par)
    })
  }
  #
  # Minimum intensity
  observeEvent(input$minInt, {
    minInt <- as.numeric(input$minInt)
    assign("minInt", minInt, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Min.Intensity <- minInt
    PARAM(Par)
  })
  # Mirror ratios
  # observeEvent(input$Mirror, {
  #   Mirror.Ratios <- as.logical(input$Mirror)
  #   assign("Mirror.Ratios", Mirror.Ratios, envir = .GlobalEnv)
  #   Par <- PARAM()
  #   Par$Mirror.Ratios <- Mirror.Ratios
  #   PARAM(Par)
  # })
  # Normalisations
  observeEvent(input$Norm, {
    output$Norm <- updtNorm()
    NORMALIZE(input$Norm)
    Par <- PARAM()
    Par$Norma.Ev.Intens <- input$Norm
    Par$Norma.Pep.Intens <- input$Norm
    Par$Norma.Prot.Ratio <- input$Norm & (input$Quant_algorithm != "limpa")
    PARAM(Par)
  })
  observeEvent(input$evLM, {
    Par <- PARAM()
    Par$Adv.Norma.Ev.Intens <- input$evLM
    PARAM(Par)
  })
  observeEvent(input$Norma.Prot.Ratio.to.proteins, {
    Par <- PARAM()
    Par$Norma.Prot.Ratio.to.proteins <- paste(db$`Protein ID`[dbOrd][match(input$Norma.Prot.Ratio.to.proteins, protHeads)],
                                              collapse = ";")
    PARAM(Par)
  }, ignoreNULL = FALSE)
  if (Annotate) {
    observeEvent(input$Norma.Prot.Ratio.to.GO, {
      Par <- PARAM()
      tmpGO <- allGO2[match(input$Norma.Prot.Ratio.to.GO, allGO)]
      Par$Norma.Prot.Ratio.to.GO <- paste(tmpGO, collapse = ";")
      PARAM(Par)
    }, ignoreNULL = FALSE)
  }
  observeEvent(input$addClassicNorm, {
    tmp <- NORMSEQ()
    tmp <- c(tmp, paste0(length(tmp)+1, " - ", input$ClassicNormMeth))
    NORMSEQ(tmp)
    updateSelectInput(inputId = "PepNormSeq",
                      choices = tmp,
                      selected = tmp)
    normSequence <- normSeqProc21(tmp)
    assign("normSequence", normSequence, envir = .GlobalEnv)
  }, ignoreInit = TRUE)
  observeEvent(input$addVarCorr, {
    tmp <- NORMSEQ()
    tmp <- c(tmp, paste0(length(tmp)+1, " - ", input$VarCorrMeth))
    NORMSEQ(tmp)
    updateSelectInput(inputId = "PepNormSeq",
                      choices = tmp,
                      selected = tmp)
    normSequence <- normSeqProc21(tmp)
    assign("normSequence", normSequence, envir = .GlobalEnv)
  }, ignoreInit = TRUE)
  observeEvent(input$addComBat, {
    tmp <- NORMSEQ()
    tmp <- c(tmp, paste0(length(tmp)+1, " - ", paste0("ComBat: ", paste(input$ComBatBatch, collapse = ";"))))
    NORMSEQ(tmp)
    updateSelectInput(inputId = "PepNormSeq",
                      choices = tmp,
                      selected = tmp)
    normSequence <- normSeqProc21(tmp)
    assign("normSequence", normSequence, envir = .GlobalEnv)
  }, ignoreInit = TRUE)
  if (LabelType == "Isobaric") {
    observeEvent(input$addIRS, {
      tmp <- NORMSEQ()
      tmp <- c(tmp, paste0(length(tmp)+1, " - IRS"))
      NORMSEQ(tmp)
      updateSelectInput(inputId = "PepNormSeq",
                        choices = tmp,
                        selected = tmp)
      normSequence <- normSeqProc21(tmp)
      assign("normSequence", normSequence, envir = .GlobalEnv)
    }, ignoreInit = TRUE)
  }
  observeEvent(input$PepNormSeq, {
    tmp <- input$PepNormSeq
    tmp <- paste0(1:length(tmp), gsub("^[0-9]+ - ", " - ", input$PepNormSeq))
    NORMSEQ(tmp)
    updateSelectInput(inputId = "PepNormSeq",
                      choices = tmp,
                      selected = tmp)
    normSequence <- normSeqProc21(tmp)
    assign("normSequence", normSequence, envir = .GlobalEnv)
  }, ignoreInit = TRUE)
  #
  if (scrptTypeFull == "withReps_PG_and_PTMs") {
    observeEvent(input$prtLM, {
      Par <- PARAM()
      Par$Adv.Norma.Prot.Intens <- input$prtLM
      PARAM(Par)
    })
    if (WorkFlow == "BIOID") {
      observeEvent(input$Norma.Prot.Ratio.to.Biot, {
        Par <- PARAM()
        Par$Norma.Prot.Ratio.to.Biot <- input$Norma.Prot.Ratio.to.Biot
        PARAM(Par)
      })
    }
  }
  # Purity correction
  if (LabelType == "Isobaric") {
    observe({ shinyFileChoose(input, "PurityFl", roots = getVolumes(), filetypes = "csv")
      {
        tmp <- input$PurityFl
        if ((!is.null(tmp))&&(is.list(tmp))) {
          tmp <- parseFilePaths(getVolumes(), tmp)$datapath
          tmp <- normalizePath(tmp, winslash = "/")
          PURFL(tmp)
          Par <- PARAM()
          Par$Label.Purities.file <- tmp
          PARAM(Par)
          output$IsobarCorr <- updtIsoPur()
        }
    }
    })
  }
  # Default type of T-test
  observeEvent(input$TtstPval, {
    Par <- PARAM()
    Par$P.values.type <- input$TtstPval
    PARAM(Par)
  })
  # FDR values
  output$FDR <- renderUI(HTML(Param$BH.FDR.values))
  observeEvent(input$AddFDR, {
    Par <- PARAM()
    tmp <- sort(unique(c(unlist(strsplit(Par$BH.FDR.values, ";")), input$FDR)))
    tmp <- tmp[which((tmp <= 1)&(tmp > 0))]
    tmp <- paste(tmp, collapse = ";") 
    output$FDR <- renderUI(HTML(tmp))
    Par$BH.FDR.values <- tmp
    PARAM(Par)
  })
  observeEvent(input$RemvFDR, {
    Par <- PARAM()
    tmp <- sort(unique(unlist(strsplit(Par$BH.FDR.values, ";"))))
    tmp <- tmp[which(tmp != input$FDR)]
    tmp <- paste(tmp, collapse = ";") 
    output$FDR <- renderUI(HTML(tmp))
    Par$BH.FDR.values <- tmp
    PARAM(Par)
  })
  # Ratio-level thresholds
  observeEvent(input$typeOfThresh, {
    ratiosThresh <- input$typeOfThresh
    assign("ratiosThresh", ratiosThresh, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Ratios.Thresholds <- ratiosThresh
    PARAM(Par)
  })
  observeEvent(input$RatCont, {
    KontRt <- as.numeric(input$RatCont)
    assign("KontRt", KontRt, envir = .GlobalEnv)
    Par <- PARAM()
    Par$Ratios.Contamination.Rates <- KontRt
    PARAM(Par)
  })
  observeEvent(input$RatContGrp, {
    Par <- PARAM()
    Par$Ratios.Contaminant.Groups <- input$RatContGrp
    PARAM(Par)
  })
  # Is the design nested?
  observeEvent(input$IsNested, {
    Par <- PARAM()
    Par$Ratios.Groups_Nested <- as.logical(input$IsNested)
    PARAM(Par)
  })
  # Proteins of interest
  observeEvent(input$IntProt, {
    prot.list <- db$`Protein ID`[dbOrd][match(input$IntProt, protHeads)]
    assign("prot.list", prot.list, envir = .GlobalEnv)
    Par <- PARAM()
    ##Par$Prot.list_pep <-
    Par$Prot.list <- paste(prot.list, collapse = ";")
    PARAM(Par)
  }, ignoreNULL = FALSE)
  # Use curved SAM thresholds for Student's t-test
  observeEvent(input$useSAM_thresh, {
    useSAM_thresh <- as.logical(input$useSAM_thresh)
    assign("useSAM_thresh", useSAM_thresh, envir = .GlobalEnv)
    #output$F_test_grps <- updtFTstUI()
  })
  # F-test
  observeEvent(input$run_F_test, {
    Par <- PARAM()
    Par$F.test <- tmp <- as.logical(input$run_F_test)
    F_TEST(tmp)
    PARAM(Par)
    #output$F_test_grps <- updtFTstUI()
  })
  # observeEvent(input$F.test_within, {
  #   Par <- PARAM()
  #   Par$F.test_within <- paste(substr(input$F.test_within, 1 , 3), collapse = ";")
  #   PARAM(Par)
  # }, ignoreNULL = FALSE)
  # GO enrichment
  observeEvent(input$GOenrich, {
    if (input$GOenrich) { shinyjs::enable("runClueGO") }
    if (!input$GOenrich) { shinyjs::disable("runClueGO") }
    Par <- PARAM()
    Par$GO.enrichment <- paste(input$GOenrich, collapse = ";")
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
      Par <- PARAM()
      tmpGO <- allGO2[match(input$GO.tabs, allGO)]
      Par$GO.tabs <- paste(tmpGO, collapse = ";")
      PARAM(Par)
    }, ignoreNULL = FALSE)
  }
  observeEvent(input$GO2Int, {
    Par <- PARAM()
    Par$GO.terms.for.proteins.of.interest <- paste(input$GO2Int, collapse = ";")
    PARAM(Par)
  })
  #
  observeEvent(input$Amica, {
    Par <- PARAM()
    Par$Amica <- input$Amica
    PARAM(Par)
  })
  # SAINTexpress
  observeEvent(input$saintExprs, {
    saintExprs <- as.logical(input$saintExprs)
    assign("saintExprs", saintExprs, envir = .GlobalEnv)
    Par <- PARAM()
    Par$saintExpress <- saintExprs
    PARAM(Par)
  })
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
          CytoScExe <- input$CytoScVers2
          assign("CytoScExe", CytoScExe, envir = .GlobalEnv)
          tmp <- parseFilePaths(getVolumes(), tmp)$datapath
          Par <- PARAM()
          Par$CytoscapePath <- normalizePath(tmp, winslash = "/")
          Par$CytoScape
          PARAM(Par)
        }
    }
    })
  }
  # PTMs to test statistically
  observeEvent(input$PTMsStats, {
    Par <- PARAM()
    Par$PTM.analysis <- paste(input$PTMsStats, collapse = ";")
    PARAM(Par)
  }, ignoreNULL = FALSE)
  # Re-normalize PTM peptides
  observeEvent(input$PTMsReNorm, {
    Par <- PARAM()
    Par$PTM.analysis_Norm <- input$PTMsReNorm
    PARAM(Par)
  })
  observeEvent(input$NAsReplMethod, {
    NAsReplMeth <- match(input$NAsReplMethod, NAsReplMethods)
    assign("NAsReplMeth", NAsReplMeth, envir = .GlobalEnv)
    Par <- PARAM()
    Par$PTM.analysis_NAsReplaceMethod <- NAsReplMeth
    PARAM(Par)
  })
  #
  # Save
  observeEvent(input$saveBtn, {
    Par <- PARAM()
    for (w in wMp) { # Extra verification
      if (nchar(Par[[w]])) { Par[[w]] <- paste(substr(unlist(strsplit(Par[[w]], ";")), 1, 3), collapse = ";") }
    }
    assign("Param", Par, envir = .GlobalEnv)
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
if (F_test_override) {
  Param$F.test <- FALSE
}
Param$Volcano.plots.Aggregate.Level <- Param_filter(Param$Ratios.Groups.Ref.Aggregate.Level, "Rep")
Param$Ratios.Ref.Groups <- paste0(Param$Ratios.Groups, c("", ";Rep")[Param$Is])
Param$Ratios.Plot.split <- "Exp"
tmp <- unlist(strsplit(Param$Ratios.Groups.Ref.Aggregate.Level, ";"))
tmp <- c(tmp[which(!tmp %in% c("Exp", "Rep"))], "Rep")
Param$Ratios.Plot.wrap <- tmp[1]
Param$Ratios.Plot.colour <- tmp[min(c(2, length(tmp)))]
if (Param$Ratios.Thresholds == "% of intra-sample group ratios") {
  Param$Ratios.Contamination.Rates <- Param$Ratios.Contamination.Rates/100
}
# Apply defaults to non-UI parameters
g <- grep("^TF_((TRUE)|(FALSE))$", Param[1,])
for (w in g) { Param[[w]] <- as.logical(gsub("^TF_", "", Param[[w]])) }
g <- grep("^((TRUE)|(FALSE))$", Param[1,])
for (w in g) { Param[[w]] <- as.logical(Param[[w]]) }
#}
if (runClueGO&&!enrichGO) { runClueGO <- FALSE}
#
Param$FullQuant <- TRUE
FullQuant %<o% Param$FullQuant
protrul %<o% Param$ProtRul
Update_Prot_matches <- Param$Update_Prot_matches
ProtRulNuclL %<o% Param$ProtRulNuclL
CytoScVrs %<o% Param$CytoscapePath
CytoScape %<o% file.exists(CytoScVrs)
Nested %<o% as.logical(Param$Ratios.Groups_Nested)
Param$Ratios.Ref.Groups <- { if (!Nested) {
  Param_filter(Param$Ratios.Ref.Groups, "Rep")
} else {
  paste(unique(c(unlist(strsplit(Param$Ratios.Ref.Groups, ";")), "Rep")), collapse = ";")
} }
#
# Defaults in case we missed a parameter
g <- grep("^TF_((TRUE)|(FALSE))$", Param[1,])
if (length(g)) {
  Param[1, g] <- gsub("^TF_", "",  Param[1, g])
  Param[, g] <- as.logical(Param[, g])
}
w <- grep("^ *((TRUE)|(FALSE)) *$", Param[1,], ignore.case = TRUE)
if (length(w)) {
  Param[1, w] <- gsub(" ", "", toupper(Param[1, w]))
  Param[, w] <- as.logical(Param[, w])
}

ReportCalls$Calls <- AddTxt2Report(" -> Parameters:")
for (i in 1:ncol(Param)) {
  ReportCalls$Calls <- AddTxt2Report(paste0("   - ", colnames(Param)[i], ": ", Param[[i]]))
}
ReportCalls <- AddSpace2Report()
#
# Create sub-directories vector:
# (this should go...)
dir <- c("Workflow control/MA plots", paste0("Workflow control/", evNm, "s", c("", "/Normalisation")),
         "Workflow control/Peptides", "Workflow control/Peptides/PCA plot", "Workflow control/Peptides/Intensities",
         "Workflow control/Peptides/Intensities/Normalisation - summary", "Workflow control/Peptides/Ratios",
         "Workflow control/Protein groups", "Workflow control/Protein groups/Expression",
         "Workflow control/Protein groups/Ratios", "Workflow control/Protein groups/P-values",
         "Reg. analysis", "Reg. analysis/t-tests",
         "PCA plots", "t-SNE plots", "Heatmaps", "Tables")
enrichGO %<o% (("GO.enrichment" %in% colnames(Param))&&(Param$GO.enrichment))
globalGO %<o% (("GO.enrichment_Whole_dataset" %in% colnames(Param))&&(Param$GO.enrichment_Whole_dataset))
if (enrichGO||globalGO) { dir <- c(dir, "Reg. analysis/GO enrich") }
dirlist <- unique(c(dirlist, paste0(wd, "/", dir)))

# Refresh list of interesting proteins:
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList2.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)

if ("Prot.list_separate.plots" %in% colnames(Param)) {
  protsplit %<o% Param$Prot.list_separate.plots
  if (protsplit == "") { protsplit <- FALSE }
} else { protsplit %<o% FALSE }
#
if (Param$GO.terms.for.proteins.of.interest) {
  tmpGO <- unlist(strsplit(Param$GO.tabs, ";"))
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
  writeFasta(db[match(prot.list, db$`Protein ID`),], intPrtFst)
}

# Similar list as above: proteins for which a heatmap of peptides and a coverage map will be drawn:
if (!exists("prot.list_pep")) { prot.list_pep %<o% c() }
if ("Prot.list_pep" %in% colnames(Param)) {
  tmp <- as.character(Param$Prot.list_pep)
  if ((!as.character(tmp) %in% c("", "NA"))&&(rev(unlist(strsplit(tmp, "\\.")))[1] == "csv")) {
    if (tmp %in% list.files()) {
      tmp <- read.csv(tmp)
      if ("Protein.ID" %in% colnames(tmp)) { tmp <- tmp$Protein.ID } else {
        warning("There is no \"Protein.ID\" column in the file you provided!")
        tmp <- c()
      }
    } else {
      warning("I could not find the protein list for peptides coverage/heatmaps!")
      prot.list_pep <- c() }
  } else {
    tmp <- { if (is.na(tmp)||(tmp == "")) { c() } else { unlist(strsplit(tmp, ";")) } }
  }
} else { tmp <- c() }
if (exists("TargetProteins")) { tmp <- c(tmp, TargetProteins) }
prot.list_pep %<o% unique(c(prot.list_pep, tmp))
# Filter lists to only keep existing ones
if (prot.list.Cond) {
  prot.list <- gsub("^CON_+", "", prot.list)
  prot.list <- prot.list[which(prot.list %in% db$"Protein ID")]
}
if (length(prot.list_pep)) {
  prot.list_pep <- gsub("^CON_+", "", prot.list_pep)
  prot.list_pep <- prot.list_pep[which(prot.list_pep %in% db$"Protein ID")]
}

# Custom protein groups
custPGs %<o% NA
if (("Custom.PGs" %in% colnames(Param))&&(!as.character(Param$Custom.PGs) %in% c("", "NA", "FALSE"))&&(file.exists(Param$Custom.PGs))) {
  custPGs_fl %<o% Param$Custom.PGs
  custPGs <- read.delim(custPGs_fl, check.names = FALSE)
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.delim(custPGs_fl, check.names = FALSE, sep = ",") }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_fl, check.names = FALSE, sep = ",") }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_fl, check.names = FALSE, sep = "\t") }
  if (colnames(custPGs)[1] != "Leading protein IDs") {
    warning("I could not make sense of the custom protein groups file provided, skipping!")
    custPGs <- NA
  } else {
    prot.list <- unique(c(prot.list, unlist(strsplit(custPGs$"Leading protein IDs", ";"))))
    prot.list_pep <- unique(c(prot.list_pep, unlist(strsplit(custPGs$"Leading protein IDs", ";"))))
  }
}

# Database of user-defined contaminants (not the default ones from CCP or MaxQuant)
## For cases where you searched with a second database of custom contaminants (e.g. E. coli for C. elegans samples)
if (("Cont.DB" %in% colnames(Param))&&(!toupper(as.character(Param$Cont.DB)) %in% c("", " ", "NA", "F", "FALSE"))) {
  temp <- unlist(strsplit(Param$Cont.DB, ";"))
  tst <- file.exists(temp)
  if (sum(!tst)) {
    msg <- paste0("The following contaminant fasta", c("", "s")[(sum(!tst) > 1)+1], "could not be found:", paste0(" - ", temp[which(!tst)], "\n"))
    stop(msg)
    #ReportCalls <- AddMsg2Report(Space = FALSE)
    #temp <- temp[which(tst)]
  }
  temp <- lapply(temp, Format.DB)
  temp <- plyr::rbind.fill(temp)
  w <- which(is.na(temp), arr.ind = TRUE)
  temp[w] <- ""
  temp$Organism_Full[which(temp$Organism_Full == "")] <- "Contaminant"
  temp$Organism[which(temp$Organism == "")] <- "Contaminant"
  temp$"Protein ID" <- paste0("CON__", gsub("^CON__", "",  temp$"Protein ID"))
  temp$"Potential contaminant" <- "+"
  # Remove all evidences which match one of these proteins:
  #test <- strsplit(ev$Proteins, ";")
  #test <- vapply(test, function(x) { sum(x %in% temp$"Protein ID") }, 1) > 0
  #cont.ev %<o% ev[which(test),]
  #ev <- ev[which(!test),]
  contDB <- { if (exists("contDB")) { plyr::rbind.fill(list(contDB, temp)) } else { temp } }
  w <- which(is.na(contDB), arr.ind = TRUE)
  contDB[w] <- ""
  db <- plyr::rbind.fill(list(db, temp))
  w <- which(is.na(db), arr.ind = TRUE)
  db[w] <- ""
}
w <- which(is.na(db$"Protein ID"))
stopifnot(length(w) == 0) # This would need immediate fixing => stop early!

# Write search database as backup
tmp <- rep("", nrow(db)*3)
tmp[3*(1:nrow(db))-2] <- db$Header
tmp[3*(1:nrow(db))-1] <- db$Sequence
write(tmp, paste0(wd, "/Concatenated search database.fasta"))


# Filters
# _______
#
# 1) Optional True/False-Discovery Filter
# Applied later at protein groups level.
# Any protein not in the filter sees its quantitative values set to 0.
#
DiscFilt %<o% ((!is.na(Param$TrueDisc_filter))&(Param$TrueDisc_filter != "")&(file.exists(Param$TrueDisc_filter)))
if (DiscFilt) {
  DiscFiltFl %<o% Param$TrueDisc_filter
  DiscFiltTbl %<o% read.csv(DiscFiltFl, check.names = FALSE)
  tst2 <- sum(!c("Protein ID", RG$values) %in% colnames(DiscFiltTbl)) == 0
  tst3 <- sum(!vapply(RG$values, function(x) { !"logical" %in% class(DiscFiltTbl[[x]]) }, TRUE))
  if (tst2+tst3 == 2) {
    DiscFiltMode %<o% Param$TrueDisc_filter_mode
    if (!DiscFiltMode %in% DiscFiltModes) { DiscFiltMode <- DiscFiltModes[1] }
    if (DiscFiltMode == "Filter column") {
      ObjNm <- "DiscFiltCol"
      if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
        msg <- "How should we name the filter column?"
        tmp <- dlg_input(msg, "Found in ...")$res
        ObjNm %<c% tmp
        AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
        tmp <- AllAnsw[1,]
        tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
        tmp$Value <- list(get(ObjNm))
        m <- match(ObjNm, AllAnsw$Parameter)
        if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
      }
    }
  } else {
    if (!tst2) { warning("Invalid filter (missing columns), skipping...") }
    if (!tst3) { warning("Invalid filter (except \"Protein ID\", all otherfilter columns should be logicals), skipping...") }
    DiscFilt <- FALSE
  }
}
#
# 2) Optional "CRAPome" filter
# For pull downs mostly.
# Provides a vector of accessions. Any protein group with one protein in this filter will be marked as a contaminant.
CRAPome %<o% ((!is.na(Param$CRAPome_file))&(Param$CRAPome_file != "")&(file.exists(Param$CRAPome_file)))
if (CRAPome) {
  CRAPomeFl %<o% Param$CRAPome_file
  CRAPomeProteins %<o% read.csv(CRAPomeFl, check.names = FALSE)
  stopifnot(sum(c("Protein ID", "Protein IDs") %in% colnames(CRAPomeProteins)) > 0)
  CRAPomeProteins <- c(CRAPomeProteins$`Protein ID`, CRAPomeProteins$`Protein IDs`)
  CRAPomeProteins <- unique(unlist(strsplit(CRAPomeProteins, ";")))
  if (length(CRAPomeProteins)) {
    w1 <- which(CRAPomeProteins %in% db$`Protein ID`)
    w2 <- which(!CRAPomeProteins %in% db$`Protein ID`)
    if (length(w2)) {
      warning(paste0(w2, " proteins (", round(100*length(w2)/nrow(db), 2), "%) in the provided CRAPome are not in the search database!"))
    }
  } else {
    warning("Empty CRAPome filter, skipping...")
    CRAPome <- FALSE
  }
}

# Pull-down specific parameters
IsPullDown %<o% (gsub(" |_|-|\\.", "", toupper(Param$Type)) %in% c("IP", "IMMUNOPRECIPITATION", "BIOID", "PULLDOWN"))
TwoSided %<o% { if (("Two.sided" %in% colnames(Param))&&(is.logical(Param$Two.sided))) {
  Param$Two.sided
} else { !IsPullDown } }
# Impute
opt <- setNames(c(TRUE, FALSE), c("Yes", "No"))
dflt <- c("No", "Yes")[IsPullDown+1]
ObjNm <- "Impute"
if ("Pep.Impute" %in% colnames(Param)) { Impute %<o% as.logical(Param$Pep.Impute) } else {
  dflt <- "No"
  if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
    msg <- "Do you want to impute missing peptide values (recommended for pull-down experiments)?"
    tmp <- paste(rep(" ", 180), collapse = "")
    tmp <- opt[gsub(" +$", "", dlg_list(paste0(names(opt), tmp), paste0(dflt, tmp), title = msg)$res)]
    if (is.na(tmp)) { tmp <- FALSE }
    ObjNm %<c% tmp
    AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
    tmp <- AllAnsw[1,]
    tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
    tmp$Value <- list(get(ObjNm))
    m <- match(ObjNm, AllAnsw$Parameter)
    if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
  }
}
# Save parameters
tmp <- data.frame(Param = colnames(Param),
                  Value = vapply(colnames(Param), function(x) {
                    x <- as.character(Param[[x]])
                    l <- length(x)
                    if (l > 1) {
                      #g <- grep("^[A-Z][a-z]{2}$", x)
                      x <- paste(x, collapse = #c(
                                   ";"#, "")[(length(g) == l) + 1]
                                 )
                    }
                    return(x)
                  }, ""))
tmp$Help <- vapply(tmp$Param, function(x) {
  if (x %in% names(Param_Help)) { return(Param_Help[x]) }
  return("")
}, "")
colnames(tmp) <- NULL
tst <- try(write.csv(tmp, file = ParamPath, row.names = FALSE), silent = TRUE)
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
  dlg_message(paste0("File \"", ParamPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(tmp, file = ParamPath, row.names = FALSE), silent = TRUE)
}

# Important!
# We often copy and paste factor levels -> they must be stored as characters
# Otherwise, they may get pasted improperly, e.g. if there are numerics with different number of characters,
# such as 1 and 10, then this gets pasted as "... 1" and "...10" which is SUPER annoying!!!
w <- which(vapply(colnames(Exp.map), function(x) { sum(c("numeric", "integer") %in% class(Exp.map[[x]])) > 0 }, TRUE))
for (i in w) {
  Exp.map[[i]] <- gsub(" ", "", # Overkill, just in case!
                       as.character(Exp.map[[i]]))
}

#### Code chunk - Create factor aggregates
a <- Aggregates
names(a) <- NULL
Aggregate.map %<o% data.frame(Aggregate.Name = names(Aggregates), Characteristics = a)
if (length(Aggregates) > 1) {
  temp1 <- c()
  temp2 <- list()
  kount <- 0
  for (i in 2:length(Aggregates)) {
    I <- combn(names(Aggregates), i)
    for (j in 1:ncol(I)) {
      kount <- kount + 1
      J <- I[,j]
      names(J) <- Aggregates[match(J, names(Aggregates))]
      Exp.map[paste(J, collapse = "")] <- apply(Exp.map[, names(J)], 1, function(x) {
        gsub(" ", "", paste(x, collapse = "___"))
      })
      temp1 <- c(temp1, paste(J, collapse = ""))
      temp2[[kount]] <- names(J)
      paste(J, collapse = "") %<c% unique(Exp.map[[paste(J, collapse = "")]])
    }
  }
  temp1 <- data.frame(Aggregate.Name = temp1)
  temp1$Characteristics <- temp2
  Aggregate.map <- rbind(Aggregate.map, temp1)
}
Aggregate.list %<o% setNames(lapply(Aggregate.map$Aggregate.Name, function(x) { get(x) }), Aggregate.map$Aggregate.Name)
#This doesn't work for the master/detailed dual script approach (issue with environments!)
#
# Define reference sample aggregate, as well as ratio groups, ratio ref group and volcano plot aggregates:
if ((!"Norm.Groups" %in% colnames(Param))&&(Param$Norm.Groups == "")) { Param$Norm.Groups <- "Exp" }
Param.aggreg %<o% c()
parse.Param.aggreg.2("Ratios.Groups.Ref.Aggregate.Level", parsed.param.nm = "Ref.Sample.Aggregate")
for (i in c("Ratios.Groups", "Norm.Groups", "Ratios.Ref.Groups", "Volcano.plots.Aggregate.Level", "Ratios.Plot.split", "Ratios.Plot.wrap", "Ratios.Plot.colour")) {
  parse.Param.aggreg.2(i)
}
# Shorter synonyms
RSA %<o% Ref.Sample.Aggregate
VPAL %<o% Volcano.plots.Aggregate.Level
RRG %<o% Ratios.Ref.Groups
RG %<o% Ratios.Groups
#
if (("Adv.Norma.Pep.Ratio" %in% colnames(Param))&&(Param$Adv.Norma.Pep.Ratio != FALSE)) {
  parse.Param.aggreg.2("Adv.Norma.Pep.Ratio.Type.Group")
}
# if (("Batch.correction" %in% colnames(Param))&&(!as.character(Param$Batch.correction) %in% c("", "F", "FALSE"))) { # Deprecated!!!
#   parse.Param.aggreg.2("Batch.correction")
# }
a <- RSA$names
Exp.map$Ref.Sample.Aggregate <- { 
  if (length(a) == 1) { Exp.map[[a]] } else { do.call(paste, c(Exp.map[, a], sep = "___")) }
}
Exp.map <- Exp.map[order(Exp.map[[VPAL$column]], Exp.map$Replicate),]
if (!"list" %in% class(Exp.map$MQ.Exp)) { Exp.map$MQ.Exp <- strsplit(Exp.map$MQ.Exp, ";") }
tstMQXp <- listMelt(Exp.map$MQ.Exp, 1:nrow(Exp.map), c("MQ.Exp", "Row"))
tstMQXp <- aggregate(tstMQXp$Row, list(tstMQXp$MQ.Exp), list)
tstMQXp <- setNames(tstMQXp$x, tstMQXp$Group.1) 
MQ.Exp <- MQ.Exp[which(MQ.Exp %in% names(tstMQXp))]
ev <- ev[which(ev$MQ.Exp %in% names(tstMQXp)),]
