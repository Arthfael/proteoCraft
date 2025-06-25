### Analysis parameters
# Defaults
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
if (!exists("Nested")) { Nested <- FALSE }
pvalue.col %<o% c(StudentRoot, WelchRoot, modRoot, permRoot, samRoot, odpRoot, lrtRoot)
names(pvalue.col) <- vapply(pvalue.col, function(x) { unlist(strsplit(x, "\\.|\\'|\\ "))[1] }, "")
ParamFls <- c(paste0(wd, "/Parameters.csv"),
              paste0(libPath, "/extData/Parameters_template.csv"))
ParamPath %<o% ParamFls[1]
ParamFl %<o% ParamFls[1]
if (!file.exists(ParamFl)) { ParamFl <- ParamFls[2] }
Param_Help <- read.csv(ParamFl, header = FALSE)
Param_Help <- Param_Help$V3
Param %<o% Param.load(ParamFl)
Param$vCPUs <- N.clust
Param$WD <- wd
if (ParamFl == ParamFls[2]) {
  Param$WD <- wd
  Param$Project <- dtstNm
  Param$Output <- SearchSoft
  Param$Fix.MQ.Isobaric.labels <- FALSE
  Param$Type <- WorkFlow
  Param$Label <- LabelType
  Param$MQ.Experiments <- paste(MQ.Exp, collapse = ";")
  Param$Search.DB <- paste(fastasTbl$Full, collapse = ";")
  Param$Search.DB.type <- paste(fastasTbl$Type, collapse = ";")
  Param$Search.DB.species <- paste(fastasTbl$Species, collapse = ";") 
  Param$Two.sided <- !(WorkFlow %in% c("PULLDOWN", "BIOID"))
  Param$Min.Pep.Size <- MinPepSz
  Param$PSMs <- paste(PSMsFl, collapse = ";")
  if (LabelType == "Isobaric") {
    Param$Label <- IsobarLab #
    Param$Label.Multiplicity <- length(get(IsobarLab)) #
    Param$Norma.Pep.Intens.IRS = "TF_TRUE"
    Param$Norma.Pep.Intens.IRS_Ref_channels = "SELECT"
    Param$Label.Purities.file <- ""
  }
  if ("Target" %in% Factors) {
    tmp <- FactorsLevels$Target
    tmp <- tmp[which((!is.na(tmp))&(tmp != "NA"))]
    Param$Prot.list <- paste(unique(c(tmp, unlist(strsplit(Param$Prot.list, ";")))), collapse = ";")
    Param$Prot.list_pep <- paste(unique(c(tmp, unlist(strsplit(Param$Prot.list_pep, ";")))), collapse = ";")
  }
  ptmDflt2 <- Param$PTM.analysis <- paste(grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE), collapse = ";")
  fTstDflt <- Param$F.test
  if (grepl("^TF_((TRUE)|(FALSE))$", fTstDflt)) { fTstDflt <- as.logical(gsub("^TF_", "", fTstDflt)) }
  if (!is.logical(fTstDflt)) { fTstDflt <- TRUE }
  Param$F.test <- fTstDflt
  Param$GO.enrichment <- TRUE
} else {
  ptmDflt2 <- ""
  if ("PTM.analysis" %in% colnames(Param)) {
    tmp <- unlist(strsplit(as.character(Param$PTM.analysis), ";"))
    tmp <- tmp[which(tmp %in% Modifs$`Full name`)]
    ptmDflt2 <- tmp
  }
}
if ((Param$Label == "LFQ")&&(isDIA)) { Param$Label <- "DIA" }
Param$WD <- wd # Super important!
ptmDflt1 <- grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE, invert = TRUE)
if (!"PTM.analysis" %in% colnames(Param)) { Param$PTM.analysis <- paste(ptmDflt2, collapse = ";") }
Mod4Quant %<o% Modifs$Mark[match(ptmDflt1, Modifs$`Full name`)]
Mod2Xclud %<o% set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                            c("Mark", "Where"))
goDflt <- suppressWarnings(as.logical(Param$GO.enrichment))
if ((!is.logical(goDflt))||(is.na(goDflt))) { goDflt <- TRUE }
fTstDflt <- suppressWarnings(as.logical(Param$F.test))
if ((!is.logical(fTstDflt))||(is.na(fTstDflt))) { fTstDflt <- TRUE }
Param$F.test <- fTstDflt
# if ((!"Param_suppress_UI" %in% colnames(Param))||(!is.logical(Param$Param_suppress_UI))) {
#   # This is to allow bypassing UI-based parameters creation.
#   # Could be useful in cases you have created custom values for parameters which the script can handle but which,
#   # because I have made a mistake (those do happen), is getting overwritten by the UI.
#   Param$Param_suppress_UI <- FALSE
#   #
#   # Honestly, this is bloat and needs to go!
# }
if (!"Ratios.Groups_Nested" %in% colnames(Param)) { Param$Ratios.Groups_Nested <- WorkFlow != "Regulation" }
#Param$Ratios.Groups.Ref.Aggregate.Level <- "AUTOFACT"
#Param$Ratios.Groups.Ref.Aggregate.Level <- "Exp;Con;Rep"
tmp <- unlist(strsplit(Param$Ratios.Groups.Ref.Aggregate.Level, ";"))
tst <- (length(tmp)>1)||(!tmp %in% c("AUTOFACT", "MAP2FACTS"))
if (tst) {
  tst <- sum(!tmp %in% substr(Factors, 1, 3)) == 0
  if (!tst) { tmp <- "MAP2FACTS" }
}
if ((length(tmp) == 1)&&(tmp == "AUTOFACT")) {
  Param$Ratios.Groups.Ref.Aggregate.Level <- paste(substr(Factors, 1, 3), collapse = ";")
}
klustChoices %<o% c("K-means", "hierarchical")
KlustMeth %<o% 2
#
saintExprs %<o% (WorkFlow %in% c("PULLDOWN", "BIOID"))
if ("saintExprs" %in% colnames(Param)) {
  tmp <- Param$saintExprs <- as.logical(Param$saintExprs)
  if (!is.na(tmp)) { saintExprs <- tmp }
} else { Param$saintExprs <- saintExprs }
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
if (!"Mirror.Ratios" %in% colnames(Param)) { Param$Mirror.Ratios <- FALSE }
if (!"Custom.PGs" %in% colnames(Param)) { Param$Custom.PGs <- "" }
if (!"TrueDisc_filter" %in% colnames(Param)) { Param$TrueDisc_filter <- "" }
DiscFiltModes %<o% c("Positive filter", "Negative filter", "Filter column")
DiscFiltModesHlp <- c(" - Positive filter: TRUE means values should not be set to NA for PGs matching this protein - values are set to NA for PGs without match or with match to FALSE)",
                      " - Negative filter: TRUE means values should be set to NA for PGs matching this protein - values are left unchanged for PGs without match or with match to FALSE)",
                      " - Filter column: Values are unaffected, but a new column is created marking with \"+\" proteins found in the provided filter.")
if (!"TrueDisc_filter_mode" %in% colnames(Param)) { Param$TrueDisc_filter_mode <- DiscFiltModes[1] }
if (!"CRAPome_file" %in% colnames(Param)) { Param$CRAPome_file <- "" }
shpDflt <- c("none", "vsn", "loess")
names(shpDflt) <- gsub("^none$", "FALSE", shpDflt)
shpDflt <- shpDflt[match(as.character(Param$Norma.Pep.Intens.Shape), names(shpDflt))]
if (is.na(shpDflt)) { shpDflt <- "none" }
pr <- c("Norma.Ev.Intens", "Norma.Pep.Intens", "Adv.Norma.Pep.Intens", "Norma.Pep.Intens.IRS", "Norma.Prot.Ratio", "Adv.Norma.Prot.Intens")
for (p in pr) { if ((!p %in% colnames(Param))||(!is.logical(Param[[p]]))||(is.na(Param[[p]]))) { Param[[p]] <- TRUE } }
p <- "Adv.Norma.Ev.Intens"
if ((!p %in% colnames(Param))||(!is.logical(Param[[p]]))||(is.na(Param[[p]]))) { Param[[p]] <- (length(unique(FracMap$Fraction)) > 1)|(length(unique(FracMap$`PTM-enriched`)) > 1) } 
if (!toupper(as.character(Param$Norma.Pep.Intens.Shape)) %in% c("FALSE", "VSN", "LOESS")) { Norma.Pep.Intens.Shape <- FALSE }
pr <- c("Norma.Pep.Ratio", "Adv.Norma.Pep.Ratio", "Norma.Prot.Ratio.to.Biot")
for (p in pr) { if ((!is.logical(Param[[p]]))||(is.na(Param[[p]]))) { Param[[p]] <- FALSE } }
QuantMethods %<o% setNames(c("Prot.Quant", "Prot.Quant + weights", "Prot.Quant.Unique", "Prot.Quant.Unique + weights",
                             "Prot.Quant2 + weights", "Prot.Quant2", "IQ_MaxLFQ", "Top3", "Top1"),
                           c(paste0("Profile_avg.", c("", ", weights = -log10(PEP)/CV", c(", unique peptides in priority", ", weights = -log10(PEP)/CV, unique peptides in priority"))),
                             paste0("Profile_avg.v2", c(", weights = -log10(PEP)/CV", "")), "MaxLFQ (iq)", "Top3", "Top1"))
QMdef <- "Prot.Quant.Unique"
if (("QuantMeth" %in% colnames(Param))&&(Param$QuantMeth %in% QuantMethods)) { QMdef <- Param$QuantMeth } else {
  Param$QuantMeth <- QMdef
}
if (!QMdef %in% QuantMethods[1:6]) { QMdef <- "Prot.Quant.Unique" }
QMdefnm <- names(QuantMethods)[match(QMdef, QuantMethods)]
dfltP4Q <- "Razor"
if (("Prot.Quant.Use" %in% colnames(Param))&&(!gsub(" |_|-|\\.", "", toupper(Param$Prot.Quant.Use)) %in% c("UNIQUE", "RAZOR", "ALL"))) {
  dfltP4Q <- Param$Prot.Quant.Use
} else {
  Param$Prot.Quant.Use <- dfltP4Q
}
if (("Update_Prot_matches" %in% colnames(Param))&&(is.logical(Param$Update_Prot_matches))&&(!is.na(Param$Update_Prot_matches))) {
  Update_Prot_matches %<o% Param$Update_Prot_matches
} else {
  Update_Prot_matches %<o% TRUE
  Param$Update_Prot_matches <- Update_Prot_matches
}
if (("Reuse_Prot_matches" %in% colnames(Param))&&(is.logical(Param$Reuse_Prot_matches))&&(!is.na(Param$Reuse_Prot_matches))) {
  Reuse_Prot_matches %<o% Param$Reuse_Prot_matches
} else {
  Reuse_Prot_matches %<o% ("evmatch.RData" %in% list.files(wd))
  Param$Reuse_Prot_matches <- Reuse_Prot_matches
}
if (("Pep.Impute" %in% colnames(Param))&&(is.logical(Param$Pep.Impute))&&(!is.na(Param$Pep.Impute))) {
  Impute %<o% as.logical(Param$Pep.Impute)
} else {
  Impute %<o% FALSE
  Param$Pep.Impute <- Impute
}
PepFoundInAtLeast %<o% 1
if ("PepFoundInAtLeast" %in% colnames(Param)) {
  PepFoundInAtLeast %<o% suppressWarnings(as.integer(Param$PepFoundInAtLeast))
  if ((is.na(PepFoundInAtLeast))||(PepFoundInAtLeast < 1)||(PepFoundInAtLeast > nrow(Exp.map))) {
    warning("Invalid \"PepFoundInAtLeast\" parameter, defaulting to 1")
    PepFoundInAtLeast <- 1
  }
}
Param$PepFoundInAtLeast <- PepFoundInAtLeast
Exp.map$Replicate <- as.integer(Exp.map$Replicate)
mxRp <- max(as.integer(Exp.map$Replicate), na.rm = TRUE)
mxN <- max(c(2, mxRp-1))
PepFoundInAtLeastGrp %<o% mxN
if ("PepFoundInAtLeastGrp" %in% colnames(Param)) {
  PepFoundInAtLeastGrp %<o% suppressWarnings(as.integer(Param$PepFoundInAtLeastGrp))
  if ((is.na(PepFoundInAtLeastGrp))||(PepFoundInAtLeastGrp < 1)||(PepFoundInAtLeastGrp > mxRp)) {
    warning(paste0("Invalid \"PepFoundInAtLeastGrp\" parameter, defaulting to ", mxN))
    PepFoundInAtLeastGrp <- mxN
  }
}
Param$PepFoundInAtLeastGrp <- PepFoundInAtLeastGrp
Param$PepFoundInAtLeastGrp <- PepFoundInAtLeastGrp
CytoScExe %<o% c()
tmp <- grep("cytoscape", list.dirs("C:/PROGRA~1", recursive = FALSE), value = TRUE, ignore.case = TRUE)
CytoScape %<o% (length(tmp) > 0)
if ("Cytoscape" %in% colnames(Param)) { CytoScape <- Param$Cytoscape }
if (length(tmp)) {
  CytoScExe <- unlist(lapply(tmp, function(x) { grep("/Cytoscape\\.exe$", list.files(x, recursive = TRUE, full.names = TRUE), value = TRUE) }))
  if (length(CytoScExe) > 1) {
    tst <- sapply(CytoScExe, function(x) { file.info(x)$mtime })
    CytoScExe <- CytoScExe[order(tst, decreasing = TRUE)]
  }
} else {
  msg <- "Could not locate Cytoscape executable!"
  ReportCalls <- AddMsg2Report()
  CytoScape <- FALSE
}
if ("CytoScapePath" %in% colnames(Param)) {
  tmp <- normalizePath(Param$CytoScapePath, winslash = "/")
  if ((length(CytoScExe))&&(tmp %in% CytoScExe)) { CytoScExe <- tmp }
} else { Param$CytoScapePath <- CytoScExe[1] }
SpeciesTst %<o% "Unspecified"
if ("Taxonomy" %in% colnames(db)) {
  SpeciesTst <- unique(db$Taxonomy[which(gsub(" *(\\(|\\[).*", "", db[[dbOrgKol]]) == mainOrg)])
  SpeciesTst <- SpeciesTst[which(as.character(SpeciesTst) != "NA")][1]
}
KingdomTst %<o% aggregate(db$Kingdom, list(db$Kingdom), length)
KingdomTst <- KingdomTst[order(KingdomTst$x, decreasing = TRUE),]
KingdomTst <- KingdomTst$Group.1[1]
isEukaLike %<o% (KingdomTst %in% c("Eukaryota", "Archaea"))
if (("ProtRul" %in% colnames(Param))&&(is.logical(Param$ProtRul))&&(!is.na(Param$ProtRul))) {
  protrul %<o% Param$ProtRul
} else {
  protrul %<o% (!WorkFlow %in% c("PULLDOWN", "BIOID"))
  # Archaea and Eukaryotes have introns and histones, Bacteria do not
  protrul <- c(protrul, FALSE)[(!isEukaLike)+1]
  Param$ProtRul <- protrul
}
if (("ProtRulNuclL" %in% colnames(Param))&&(!is.na(as.integer(Param$ProtRulNuclL)))) {
  ProtRulNuclL <- as.integer(Param$ProtRulNuclL)
} else {
  ProtRulNuclL <- 196
  Param$ProtRulNuclL <- ProtRulNuclL
}
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
if ("Ratios.Contaminant.Groups" %in% colnames(Param)) {
  KontGrp <- Param$Ratios.Contaminant.Groups
} else { KontGrp <- KontGrps[1] }
if (!KontGrp %in% KontGrps) { KontGrp <- KontGrps[1] }
if (!exists("minInt")) { minInt <- 0 }
if ("Min.Intensity" %in% colnames(Param)) {
  tmp <- suppressWarnings(as.numeric(Param$Min.Intensity))
  if ((is.numeric(tmp))&&(is.finite(tmp))&&(tmp >= 0)) { minInt <- tmp }
}
minInt %<o% minInt
if ("BH.FDR.values" %in% colnames(Param)) {
  BH.FDR <- sort(as.numeric(unlist(unique(strsplit(as.character(Param$BH.FDR.values), ";")))))
} else { BH.FDR <- c(0.1, 0.2, 0.3) }
BH.FDR %<o% BH.FDR
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
if ("Mirror.Ratios" %in% colnames(Param)) { Mirror.Ratios <- Param$Mirror.Ratios <- as.logical(Param$Mirror.Ratios) }
Param$Mirror.Ratios <- Mirror.Ratios
TwoSidedDeflt <- c(c("Up-only", "Down-only")[Mirror.Ratios+1], "Both directions")[TwoSidedDeflt+1]
# ROC
ROC_GOterms %<o% c()
ROCfilt_GOterms_Pos %<o% c()
ROCfilt_GOterms_Neg %<o% c()
annotRep <- ((Annotate)&&(scrptType == "withReps"))
if (annotRep) {
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
  ROC2_GOterms_dflt <- allGO[match(ROC_GOterms, allGO2)]
} else {
  ROC2_GOterms_dflt <- c()
}
w <- c(which(allGO %in% ROC2_GOterms_dflt),
       which(!allGO %in% ROC2_GOterms_dflt))
ROC2_allGO1 <- allGO[w]
ROC2_allGO2 <- allGO2[w]
if (length(ROCfilt_GOterms_Pos)) {
  ROC1_GOterms_Pos_dflt <- allGO[match(ROCfilt_GOterms_Pos, allGO2)]
} else {
  ROC1_GOterms_Pos_dflt <- c()
}
w <- c(which(allGO %in% ROC1_GOterms_Pos_dflt),
       which(!allGO %in% ROC1_GOterms_Pos_dflt))
ROC1_allGOPos1 <- allGO[w]
ROC1_allGOPos2 <- allGO2[w]
if (length(ROCfilt_GOterms_Neg)) {
  ROC1_GOterms_Neg_dflt <- allGO[match(ROCfilt_GOterms_Neg, allGO2)]
} else {
  ROC1_GOterms_Neg_dflt <- c()
}
w <- c(which(allGO %in% ROC1_GOterms_Neg_dflt),
       which(!allGO %in% ROC1_GOterms_Neg_dflt))
ROC1_allGONeg1 <- allGO[w]
ROC1_allGONeg2 <- allGO2[w]
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
#      - RUV normalisation (see https://www.bioconductor.org/packages/release/bioc/vignettes/RUVnormalize/inst/doc/RUVnormalize.pdf)
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
pepNormMethodsDF <- data.frame(Method = vapply(1:L, function(i) { pepNormMethods[[i]]$Method }, "a"),
                               Source = vapply(1:L, function(i) { pepNormMethods[[i]]$Source }, "a"))
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
  dfltNormSeq <- dfltNormSeq[which(vapply(dfltNormSeq, function(x) { x$Method }, "a") != "IRS")]
}
if (!Nested) {
  dfltNormSeq <- dfltNormSeq[which(vapply(dfltNormSeq, function(x) { x$Method }, "a") != "ComBat")]
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
  }, "a")
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
         which((Param[1,] == "MAP2FACTS")&(!colnames(Param) %in% coreNms)))
lstFct <- list()
dfltFct <- list()
for (w in wMp) {
  Opt <- Factors
  dflt <- dfdflt <- c("Experiment", "Replicate")[(colnames(Param)[w] == "Batch.correction")+1]
  if (!Param[1, w] %in% c("AUTOFACT", "MAP2FACTS")) {
    dflt <- Factors[match(unlist(strsplit(Param[1, w], ";")), names(Factors))]
    if ((length(dflt) == 1)&&(is.na(dflt))) { dflt <- dfdflt }
  }
  if (colnames(Param)[w] == "Ratios.Groups.Ref.Aggregate.Level") {
    dflt <- unique(c(dflt, mnFct))
  }
  if (colnames(Param)[w] == "Ratios.Groups") {
    Opt <- Factors[which(Factors != "Replicate")]
    dflt <- dflt[which(dflt != "Replicate")]
  }
  lbl <- gsub("\\.", " ", colnames(Param)[w])
  if (colnames(Param)[w] %in% coreNms) {
    lbl <- unlist(strsplit(names(coreNms)[match(colnames(Param)[w], coreNms)], "___"))
  } else { lbl <- c(gsub("\\.", " ", colnames(Param)[w]), Param_Help[w]) }
  names(Opt) <- NULL
  names(dflt) <- NULL
  dfltFct[[colnames(Param)[w]]] <- dflt
  blck <- list(list(br()),
               tags$table(
                 tags$tr(width = "80%",
                         tags$td(width = "25%",
                                 div(strong(lbl[1]))),
                         tags$td(width = "55%",
                                 selectInput(colnames(Param)[w],
                                             "",
                                             Opt,
                                             dflt,
                                             TRUE,
                                             TRUE)),
                         #addTooltip(session, colnames(Param)[w], Param_Help[w], "bottom", "hover", list(container = "body"))
                 )
               ))
  lbl2 <- unlist(strsplit(lbl[2], "/n"))
  for (lbl2a in lbl2) { blck <- append(blck, list(span(em(lbl2a)), br())) }
  if (colnames(Param)[w] == "Ratios.Groups") {
    dflt <- Param$Ratios.Groups_Nested
    if (!is.logical(dflt)) { dflt <- WorkFlow != "Regulation" }
    blck <- append(blck, list(radioButtons("IsNested", "Nested design? (i.e. are replicates paired?)",
                                           c(TRUE, FALSE), dflt, TRUE)))
  }
  blck <- append(blck, list(br()))
  lstFct <- append(lstFct, blck)
}
#
useSAM_thresh %<o% TRUE
tstAdvOpt <- try(sum(file.exists(Param$Custom.PGs, Param$TrueDisc_filter, Param$CRAPome_file)) > 0)
if ("try-error" %in% class(tstAdvOpt)) { tstAdvOpt <- FALSE }
#if (!Param$Param_suppress_UI) {
#

#
appNm <- paste0(dtstNm, " - Parameters")
ui1 <- fluidPage(
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
                  withSpinner(uiOutput("FactMappings")),
                  uiOutput("RSA_msg")),
           column(8, withSpinner(plotlyOutput("PSMsPCA", height = "600px")))),
  br(),
  tags$hr(style = "border-color: black;"),
  h4(strong("Proteins of interest")),
  pickerInput("IntProt", NULL, protHeads, protDeflt, TRUE, width = "600px",
              pickerOptions(title = "Search me",
                            `live-search` = TRUE,
                            actionsBox = TRUE,
                            deselectAllText = "Clear search",
                            showTick = TRUE)),
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
    column(3,
           checkboxInput("Impute", "Impute missing peptides-level values?", Impute, "100%")),
    column(3,
           checkboxInput("Update_Prot_matches", paste0("Update ", names(SearchSoft), "'s original protein-to-peptides assignments?"), Update_Prot_matches, "100%"),
           bsTooltip("Update_Prot_matches",
                     "Checking assignments may result in removal of some identifications. It is nonetheless recommended because we have observed occasional inconsistent peptides-to-protein assignments with some search software.",
                     placement = "right", trigger = "hover", options = list(container = "body")),
           withSpinner(uiOutput("ReloadMatches"))
    )),
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
  fluidRow(column(2, selectInput("QuantMeth",
                                 "Protein Groups-level quantitation algorithm:",
                                 names(QuantMethods)[1:6],
                                 QMdefnm,
                                 width = "100%")),
           column(2,
                  checkboxInput("ProtRul",
                                "Apply Proteomic Ruler to estimate copy numbers per cell? (uses signal from all histones as reference; assumes inter-nucleosomal space = 196 bp, do not use if this assumption does not hold!)",
                                protrul,
                                "100%"),
                  numericInput("ProtRulNuclL",
                               "Use inter-nucleosome length = ? (kb)",
                               ProtRulNuclL,
                               1,
                               Inf,
                               1,
                               "100%")),
           column(2, selectInput("Prot.Quant.Use",
                                 "Peptides eligible for quantitation:",
                                 c("Unique", "Razor", "All"),
                                 dfltP4Q,
                                 width = "100%")),
           column(2,
                  h5("Use only peptidoforms found in at least how many samples in..."),
                  numericInput("PepFoundInAtLeast",
                               " -> the whole dataset?", PepFoundInAtLeast, 1, nrow(Exp.map), 1, "100%"),
                  numericInput("PepFoundInAtLeastGrp",
                               " -> one sample group at least? (overrides parameter above)",
                               PepFoundInAtLeastGrp,
                               1,
                               mxRp,
                               1,
                               "100%"))
  ),
  # Note to self: I am for now excluding some methods, because I need to add code to calculate some columns for those, namely ratios.
  # This should be remedied asap, especially since there include such community favourites as IQ (= MaxLFQ) and Top3!!!
  br(),
  tags$hr(style = "border-color: black;"),
  fluidRow(column(2, radioButtons("Clustering", "Clustering method", klustChoices, klustChoices[1], TRUE, "100%"))),
  br(),
  tags$hr(style = "border-color: black;"),
  h4(strong("Statistical testing")),
  h5(strong(" -> t-test(s)")),
  fluidRow(
    column(2,
           radioButtons("TwoSided", "Fold changes: test...", c("Both directions", "Up-only", "Down-only"),
                        TwoSidedDeflt, FALSE, "100%"),
           checkboxInput("Mirror", "Revert fold changes on plots? (default: log fold change = log2(Sample/Reference); revert: log2(Reference/Sample))",
                         Param$Mirror.Ratios, "100%")),
    column(4,
           h5(strong("Volcano plot: select default variant")),
           radioButtons("TtstPval", "", names(pvalue.col), "Moderated", TRUE, "100%"),
           h6(em(" - Welch's t-test is a modified form of Student's original version which is more robust to variance inequality.")),
           h6(em(" - Moderated t-test (//limma): re-samples individual row variances using global dataset variance to provide a more robust estimate.")),
           h6(em(" - Permutation t-test (//coin): based on permutations of the samples from each group.")),
           h6(em(" - SAM's modified t-test (//siggenes): corrects for poor variance estimates for very reproducible data by optimizing a small constant s0 added to the denominator of the test statistic.")),
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
                       KontGrp,
                       width = "100%"))),
  br(),
  fluidRow(column(2,
                  h5(strong(" -> ANOVA (moderated F-test //limma)")),
                  checkboxInput("Run?", "", fTstDflt, "100%"),
                  em("(only makes sense if Nb. sample groups > 2)")),
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
  tags$hr(style = "border-color: black;"),
  withSpinner(uiOutput("CytoScape")),
  tags$hr(style = "border-color: black;"),
  h4(strong("Post-translational modifications (PTMs)")),
  fluidRow(column(2, pickerInput("PTMsQuant", "Select PTM(s) eligible for use for Protein Groups quantitation:",
                                 Modifs$`Full name`[which(Modifs$Type == "Variable")], ptmDflt1, TRUE,
                                 pickerOptions(title = "Search me",
                                               `live-search` = TRUE,
                                               actionsBox = TRUE,
                                               deselectAllText = "Clear search",
                                               showTick = TRUE))),
           column(2, pickerInput("PTMsStats", "Select PTM(s) (if any) for which statistical tests will be performed and subtables written:",
                                 Modifs$`Full name`[which(Modifs$Type == "Variable")],
                                 unlist(strsplit(ptmDflt2, ";")), TRUE,
                                 pickerOptions(title = "Search me",
                                               `live-search` = TRUE,
                                               actionsBox = TRUE,
                                               deselectAllText = "Clear search",
                                               showTick = TRUE))),
           column(2, checkboxInput("PTMsReNorm",
                                   "Re-normalize modified peptides ratios to those of parent Protein Group(s)?",
                                   TRUE, "100%")),
           column(2, radioButtons("NAsReplMethod",
                                  "Some modified peptides do not have a quantified parent protein group to normalize to. Replaced missing values using:",
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
if (exists("appRunTest")) { rm(appRunTest) }
server1 <- function(input, output, session) {
  NORMALIZE <- reactiveVal(normDat)
  NORMSEQ <- reactiveVal(dfltNormSeq2)
  PURFL <- reactiveVal(Param$Label.Purities.file)
  RSA_Msg <- reactiveVal("")
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
                        pickerInput("nrm2PrtSlct", " -> select proteins",
                                    protHeads, nrm2PrtSlctDflt, TRUE,
                                    pickerOptions(title = "Search me",
                                                  `live-search` = TRUE,
                                                  actionsBox = TRUE,
                                                  deselectAllText = "Clear search",
                                                  showTick = TRUE))),
                 if (Annotate) {
                   column(3,
                          pickerInput("nrm2GOSlct", " -> select GO terms",
                                      nrm2GOall, nrm2GO, TRUE,
                                      pickerOptions(title = "Search me",
                                                    `live-search` = TRUE,
                                                    actionsBox = TRUE,
                                                    deselectAllText = "Clear search",
                                                    showTick = TRUE)))
                 },
                 if (WorkFlow == "BIOID") {
                   column(3,
                          checkboxInput("prt2Biot", " -> all biotinylated proteins: ",
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
                        pickerInput("ClassicNormMeth", "methods:", genPepNormMeth, genPepNormMeth[1], width = "80%")),
                 column(2,
                        actionBttn("addVarCorr", "Add variance intensity-bias correction", color = "primary", size = "xs", style = "pill"),
                        pickerInput("VarCorrMeth", "methods:", shapePepNormMeth, shapePepNormMeth[1], width = "80%")),
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
    if (reactive) { tst <- ADVOPT() } else { tst <- tstAdvOpt }
    if (tst) {
      if (reactive) { rgKol <- input[["Ratios.Groups"]] } else { rgKol <- dfltFct[["Ratios.Groups"]] }
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
      if (reactive) { tst <- ROC1ON() } else { tst <- length(c(ROCfilt_GOterms_Pos, ROCfilt_GOterms_Neg)) > 0 }
      if (tst) {
        lst <- list(
          fluidRow(column(4,
                          pickerInput("ROC1_pos", "Positive GO term(s)", ROC1_allGOPos1, ROC1_GOterms_Pos_dflt, TRUE,
                                      pickerOptions(title = "Search me",
                                                    `live-search` = TRUE,
                                                    actionsBox = TRUE,
                                                    deselectAllText = "Clear search",
                                                    showTick = TRUE))),
                   column(4,
                          pickerInput("ROC1_neg", "Negative GO term(s)", ROC1_allGONeg1, ROC1_GOterms_Neg_dflt, TRUE,
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
      if (reactive) { tst <- ROC2ON() } else { tst <- length(ROC_GOterms) > 0 }
      if (tst) {
        lst <- list(
          fluidRow(column(4,
                          em("GO terms to evaluate acceptable P-value significance thresholds, e.g. if you know that all proteins with a given term should be significant. Currently this is only used for personal interpretation of the plots and not fed back into the automated significance analysis."),
                          pickerInput("ROC2_terms", "", ROC2_allGO1, ROC2_GOterms_dflt, TRUE,
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
    if (reactive) { myMsg <- RSA_Msg() } else { myMsg <- "" }
    renderUI(h4(strong(em(myMsg, style = "color:red", .noWS = "outside"))))
  }
  output$RSA_msg <- updtRSAmsg(FALSE)
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
  output$PSMsPCA <- renderPlotly(plot_lyPSMsPCA)
  output$FactMappings <- renderUI({
    lstFct
  })
  output$RatioGroups <- renderUI({
    em(paste(unique(apply(ExpMap[, input[["Ratios.Groups"]], drop = FALSE], 1, paste, collapse = "___")), collapse = " / "))
  })
  # Factors
  sapply(wMp, function(w) {
    observeEvent(input[[colnames(Param)[w]]], {
      tmpVal <- input[[colnames(Param)[w]]]
      if (colnames(Param)[w] == "Ratios.Groups.Ref.Aggregate.Level") {
        l <- length(tmpVal)
        tst2 <- tst3 <- tst4 <- FALSE
        if (length(tmpVal)) {
          tst2 <- sum(!c("Experiment", "Replicate") %in% tmpVal)
          tmpRSA <- do.call(paste, c(Exp.map[, tmpVal], sep = "___"))
          UtmpRSA <- unique(tmpRSA)
          tst3 <- length(UtmpRSA) < nrow(Exp.map)
          tmpVal2 <- tmpVal[which(tmpVal != "Replicate")]
          if (length(tmpVal2)) {
            tmpVPAL <- do.call(paste, c(Exp.map[, tmpVal2], sep = "___"))
            tst4 <- aggregate(1:nrow(Exp.map), list(tmpVPAL, Exp.map$Replicate), length)
            tst4 <- max(tst4$x) > 1
          }
        }
        if ((l < 3)||(tst2)||(tst3)||(tst4)) {
          msg <- c()
          if ((l < 3)||(tst2)) { msg <- c(msg, "You MUST always include Experiment and Replicate here, as well as at least one other contrasted factor!") }
          if (tst3) { msg <- c(msg, "The chosen experimental Factors do not discriminate fully between some samples (rows in the experiment map), add more!") }
          if (tst4) { msg <- c(msg, "It is not allowed for several samples to have the same replicate number within a group!") }
          RSA_Msg(paste(msg, collapse = " / "))
          shinyjs::disable("saveBtn")
        } else {
          RSA_Msg("")
          shinyjs::enable("saveBtn")
        }
        output$RSA_msg <- updtRSAmsg()
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
  #suppress
  output$CytoScape <- renderUI({
    if (!length(CytoScExe)) {
      lst <- list(list(column(2, shinyFilesButton("CytoScVers2",
                                                  em("No CytoScape.exe detected, select one (optional)"),
                                                  "", FALSE))))
    }
    if (length(CytoScExe) == 1) {
      lst <- list(list(column(2, em(paste0("Version detected: ", CytoScExe)))))
    }
    if (length(CytoScExe) > 1) {
      lst <- list(list(fluidRow(column(2, pickerInput("CytoScVers1", "Choose CytoScape version:",
                                                      CytoScExe, CytoScExe[1], FALSE)))))
    }
    return(lst)
  })
  #
  # Labels purity correction - only for isobarically labelled samples
  updtIsoPur <- function(reactive = TRUE, lblType = LabelType) {
    if (lblType == "Isobaric") {
      if (reactive) {
        fl <- PURFL()
      } else {
        fl <- Param$Label.Purities.file
      }
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
             fluidRow(column(1, checkboxInput("GOenrich", "GO enrichment", goDflt, "100%")),
                      column(2, pickerInput("GO.tabs", "GO terms of interest", allGO, dftlGO, TRUE,
                                            pickerOptions(title = "Search me",
                                                          `live-search` = TRUE,
                                                          actionsBox = TRUE,
                                                          deselectAllText = "Clear search",
                                                          showTick = TRUE))),
                      column(1, checkboxInput("GO2Int", "Use GO terms to define list of proteins of interest?",
                                              as.character(toupper(Param$GO.terms.for.proteins.of.interest)) == "TRUE", # Clumsy but possibly more backwards proof
                                              "100%"))))
      )
    }
    return(lst)
  })
  #
  # updtFTstUI <- function(reactive = TRUE) {
  #   # Update UI
  #   if (reactive) { FTst <- PARAM()$F.test } else { FTst <- Param$F.test }
  #   return(renderUI({
  #     if (FTst) {
  #       lst <- list()
  #       Opt <- Factors[which(Factors != "Replicate")]
  #       dflt <- dfdflt <- "Experiment"
  #       if (!Param$F.test_within %in% c("AUTOFACT", "MAP2FACTS")) {
  #         dflt <- Factors[match(unlist(strsplit(Param$F.test_within, ";")), names(Factors))]
  #         if ((!length(dflt))||((length(dflt) == 1)&&(is.na(dflt)))) { dflt <- dfdflt }
  #       }
  #       dflt <- dflt[which(dflt != "Replicate")]
  #       lbl <- "Groups within which to perform F-tests"
  #       names(Opt) <- NULL
  #       names(dflt) <- NULL
  #       w <- which(colnames(Param) == "F.test_within")
  #       blck <- list(list(br()),
  #                    tags$table(
  #                      tags$tr(width = "100%",
  #                              tags$td(width = "25%", div(strong("F-test groups"))),
  #                              addTooltip(session, colnames(Param)[w], Param_Help[w], "bottom", "hover", list(container = "body")),
  #                      ),
  #                      tags$tr(width = "100%", tags$td(width = "55%",
  #                                                      selectInput("F.test_within",
  #                                                                  "",
  #                                                                  Opt,
  #                                                                  dflt,
  #                                                                  TRUE,
  #                                                                  TRUE)))
  #                    ))
  #       blck <- append(blck, list(span(em("F-tests will be performed within groups defined by the levels of the factors you choose here.")), br()))
  #       blck <- append(blck, list(br()))
  #       lst <- append(lst, blck)
  #     } else { lst <- list(em("No analysis yet")) }
  #     return(lst)
  #   }))
  # }
  #output$F_test_grps <- updtFTstUI(reactive = FALSE)
  #
  output$ReloadMatches <- renderUI({
    if ("evmatch.RData" %in% list.files(wd)) {
      msg <- "Peptide-to-protein matches backup detected in folder: do you want to reload it?\n"
      lst <- list(list(list(br()),
                       tags$table(
                         tags$tr(width = "100%", tags$td(width = "55%", checkboxInput("Reuse_Prot_matches", msg, TRUE)))
                       )))
    } else {
      em(" ")
    }
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
  # Optional input files
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
      Par <- PARAM()
      ROCfilt_GOterms_Pos <<- ROC1_allGOPos2[match(input$ROC1_pos, ROC1_allGOPos1)]
      Par$ROCfilt_GOterms_Pos <- paste(ROCfilt_GOterms_Pos, collapse = ";")
      PARAM(Par)
    })
    observeEvent(input$ROC1_neg, {
      Par <- PARAM()
      ROCfilt_GOterms_Neg <<- ROC1_allGONeg2[match(input$ROC1_neg, ROC1_allGONeg1)]
      Par$ROCfilt_GOterms_Neg <- paste(ROCfilt_GOterms_Neg, collapse = ";")
      PARAM(Par)
    })
    observeEvent(input$ROC2_terms, {
      Par <- PARAM()
      ROC_GOterms <<- ROC2_allGO2[match(input$ROC2_terms, ROC2_allGO1)]
      Par$ROC_GOterms <- paste(ROC_GOterms, collapse = ";")
      PARAM(Par)
    })
  }
  #
  # Minimum intensity
  observeEvent(input$minInt, {
    assign("minInt", as.numeric(input$minInt), envir = .GlobalEnv)
    Par <- PARAM()
    Par$Min.Intensity <- as.numeric(input$minInt)
    PARAM(Par)
  })
  # Impute?
  observeEvent(input$Impute, {
    Par <- PARAM()
    Par$Pep.Impute <- input$Impute
    PARAM(Par)
  })
  # Mirror ratios
  observeEvent(input$Mirror, {
    Par <- PARAM()
    Par$Mirror.Ratios <- input$Mirror
    PARAM(Par)
  })
  # Update PSM-to-Protein matches?
  observeEvent(input$Update_Prot_matches, {
    Par <- PARAM()
    Par$Update_Prot_matches <- input$Update_Prot_matches
    PARAM(Par)
    if (input$Update_Prot_matches) { shinyjs::enable("Reuse_Prot_matches") }
    if (!input$Update_Prot_matches) { shinyjs::disable("Reuse_Prot_matches") }
  })
  observeEvent(input[["Reuse_Prot_matches"]], {
    Par <- PARAM()
    Par$Reuse_Prot_matches <- input$Reuse_Prot_matches
    PARAM(Par)
  })
  # Clustering method
  observeEvent(input$Clustering, {
    assign("KlustMeth", match(input$Clustering, klustChoices), envir = .GlobalEnv)
  })
  # Are analyses Two-sided?
  observeEvent(input$TwoSided, {
    Par <- PARAM()
    Par$Two.sided <- input$TwoSided == "Both directions"
    if (input$TwoSided == "Down-only") {
      shinyjs::disable("Mirror")
      Par$Mirror.Ratios <- TRUE
    }
    if (input$TwoSided != "Down-only") { shinyjs::enable("Mirror") }
    PARAM(Par)
  })
  # Normalisations
  observeEvent(input$Norm, {
    output$Norm <- updtNorm()
    NORMALIZE(input$Norm)
    Par <- PARAM()
    Par$Norma.Ev.Intens <- input$Norm
    Par$Norma.Pep.Intens <- input$Norm
    Par$Norma.Prot.Ratio <- input$Norm
    PARAM(Par)
  })
  observeEvent(input$evLM, {
    Par <- PARAM()
    Par$Adv.Norma.Ev.Intens <- input$evLM
    PARAM(Par)
  })
  # observeEvent(input$pepShape, {
  #   Par <- PARAM()
  #   Par$Norma.Pep.Intens.Shape <- input$pepShape
  #   PARAM(Par)
  # })
  observeEvent(input$nrm2PrtSlct, {
    Par <- PARAM()
    Par$Norma.Prot.Ratio.to.proteins <- paste(db$`Protein ID`[dbOrd][match(input$nrm2PrtSlct, protHeads)],
                                              collapse = ";")
    PARAM(Par)
  }, ignoreNULL = FALSE)
  if (Annotate) {
    observeEvent(input$nrm2GOSlct, {
      Par <- PARAM()
      tmpGO <- allGO2[match(input$nrm2GOSlct, allGO)]
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
    assign("normSequence", normSeqProc21(tmp), envir = .GlobalEnv)
  }, ignoreInit = TRUE)
  observeEvent(input$addVarCorr, {
    tmp <- NORMSEQ()
    tmp <- c(tmp, paste0(length(tmp)+1, " - ", input$VarCorrMeth))
    NORMSEQ(tmp)
    updateSelectInput(inputId = "PepNormSeq",
                      choices = tmp,
                      selected = tmp)
    assign("normSequence", normSeqProc21(tmp), envir = .GlobalEnv)
  }, ignoreInit = TRUE)
  observeEvent(input$addComBat, {
    tmp <- NORMSEQ()
    tmp <- c(tmp, paste0(length(tmp)+1, " - ", paste0("ComBat: ", paste(input$ComBatBatch, collapse = ";"))))
    NORMSEQ(tmp)
    updateSelectInput(inputId = "PepNormSeq",
                      choices = tmp,
                      selected = tmp)
    assign("normSequence", normSeqProc21(tmp), envir = .GlobalEnv)
  }, ignoreInit = TRUE)
  if (LabelType == "Isobaric") {
    observeEvent(input$addIRS, {
      tmp <- NORMSEQ()
      tmp <- c(tmp, paste0(length(tmp)+1, " - IRS"))
      NORMSEQ(tmp)
      updateSelectInput(inputId = "PepNormSeq",
                        choices = tmp,
                        selected = tmp)
      assign("normSequence", normSeqProc21(tmp), envir = .GlobalEnv)
    }, ignoreInit = TRUE)
  }
  observeEvent(input$PepNormSeq, {
    tmp <- input$PepNormSeq
    tmp <- paste0(1:length(tmp), gsub("^[0-9]+ - ", " - ", input$PepNormSeq))
    NORMSEQ(tmp)
    updateSelectInput(inputId = "PepNormSeq",
                      choices = tmp,
                      selected = tmp)
    assign("normSequence", normSeqProc21(tmp), envir = .GlobalEnv)
  }, ignoreInit = TRUE)
  #
  # observeEvent(input$pepLM, {
  #   Par <- PARAM()
  #   Par$Adv.Norma.Pep.Intens <- input$pepLM
  #   PARAM(Par)
  # })
  # if ((LabelType == "Isobaric")&&(length(Iso) > 1)) {
  #   observeEvent(input$IRS, {
  #     Par <- PARAM()
  #     Par$Norma.Pep.Intens.IRS <- input$IRS
  #     PARAM(Par)
  #   })
  #   for (i in 1:length(Iso)) {
  #     observeEvent(input[[paste0("intRef", Iso[i])]], {
  #       Par <- PARAM()
  #       tmp <- dfltChan()
  #       tmp[i] <- input[[paste0("intRef", Iso[i])]]
  #       dfltChan(tmp)
  #       Par$Norma.Pep.Intens.IRS_Ref_channels <- paste(tmp, collapse = ";")
  #       PARAM(Par)
  #     }, ignoreNULL = FALSE)
  #   }
  # }
  if (scrptTypeFull == "withReps_PG_and_PTMs") {
    observeEvent(input$prtLM, {
      Par <- PARAM()
      Par$Adv.Norma.Prot.Intens <- input$prtLM
      PARAM(Par)
    })
    if (WorkFlow == "BIOID") {
      observeEvent(input$prt2Biot, {
        Par <- PARAM()
        Par$Norma.Prot.Ratio.to.Biot <- input$prt2Biot
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
  # Quantitation
  observeEvent(input$QuantMeth, {
    Par <- PARAM()
    Par$QuantMeth <- QuantMethods[match(input$QuantMeth, names(QuantMethods))]
    PARAM(Par)
  })
  observeEvent(input$ProtRul, {
    Par <- PARAM()
    Par$ProtRul <- input$ProtRul
    PARAM(Par)
  })
  observeEvent(input$ProtRulNuclL, {
    Par <- PARAM()
    Par$ProtRulNuclL <- as.integer(input$ProtRulNuclL)
    PARAM(Par)
  })
  observeEvent(input$Prot.Quant.Use, {
    Par <- PARAM()
    Par$Prot.Quant.Use <- input$Prot.Quant.Use
    PARAM(Par)
  })
  observeEvent(input$PepFoundInAtLeast, {
    Par <- PARAM()
    Par$PepFoundInAtLeast <- as.integer(input$PepFoundInAtLeast)
    PARAM(Par)
  })
  observeEvent(input$PepFoundInAtLeastGrp, {
    Par <- PARAM()
    Par$PepFoundInAtLeastGrp <- as.integer(input$PepFoundInAtLeastGrp)
    PARAM(Par)
  })
  # Type of T-test
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
    Par <- PARAM()
    Par$Ratios.Thresholds <- input$typeOfThresh
    PARAM(Par)
    assign("ratiosThresh", Par$Ratios.Thresholds, envir = .GlobalEnv)
  })
  observeEvent(input$RatCont, {
    Par <- PARAM()
    assign("KontRt", input$RatCont, envir = .GlobalEnv)
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
    Par <- PARAM()
    ##Par$Prot.list_pep <-
    Par$Prot.list <- paste(db$`Protein ID`[dbOrd][match(input$IntProt, protHeads)], collapse = ";")
    PARAM(Par)
  }, ignoreNULL = FALSE)
  # Use curved SAM thresholds for Student's t-test
  observeEvent(input$useSAM_thresh, {
    useSAM_thresh <<- input$useSAM_thresh
    #output$F_test_grps <- updtFTstUI()
  })
  # F-test
  observeEvent(input$FTest, {
    Par <- PARAM()
    Par$F.test <- input$FTest
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
    Par <- PARAM()
    Par$GO.enrichment <- paste(input$GOenrich, collapse = ";")
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
    Par <- PARAM()
    Par$saintExpress <- input$saintExprs
    assign("saintExprs", Par$saintExpress, envir = .GlobalEnv)
    PARAM(Par)
  })
  # CytoScape
  observeEvent(input$CytoScVers1, {
    Par <- PARAM()
    Par$CytoScapePath <- input$CytoScVers1
    PARAM(Par)
  })
  observe({ shinyFileChoose(input, "CytoScVers2", roots = getVolumes(), filetypes = "exe")
    {
      tmp <- input$CytoScVers2
      if ((!is.null(tmp))&&(is.list(tmp))) {
        tmp <- parseFilePaths(getVolumes(), tmp)$datapath
        Par <- PARAM()
        Par$CytoScapePath <- normalizePath(tmp, winslash = "/")
        PARAM(Par)
      }
  }
  })
  # PTMs to use for PG Quant
  observeEvent(input$PTMsQuant, {
    m4Quant(Modifs$Mark[match(unlist(input$PTMsQuant), Modifs$`Full name`)])
    m2Xclud(set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                         c("Mark", "Where")))
  }, ignoreNULL = FALSE)
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
    Par <- PARAM()
    assign("NAsReplMeth", match(input$NAsReplMethod, NAsReplMethods), envir = .GlobalEnv)
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
  session$onSessionEnded(function() { 
    stopApp()
  })
}
appTxt1 <- gsub("myApp", "myApp1", gsub("\\(ui", "(ui1", gsub(", server", ", server1", runApp)))
runKount <- 0
while ((!runKount)||(!exists("appRunTest"))) {
  eval(parse(text = appTxt1), envir = .GlobalEnv)
  runKount <- runKount+1
}
#
Param$Volcano.plots.Aggregate.Level <- Param_filter(Param$Ratios.Groups.Ref.Aggregate.Level, "Rep")
Param$Ratios.Ref.Groups <- paste0(Param$Ratios.Groups, c("", ";Rep")[Param$Is])
Param$Ratios.Plot.split <- "Exp"
tmp <- unlist(strsplit(Param$Ratios.Groups.Ref.Aggregate.Level, ";"))
tmp <- c(tmp[which(!tmp %in% c("Exp", "Rep"))], "Rep")
Param$Ratios.Plot.wrap <- tmp[1]
Param$Ratios.Plot.colour <- tmp[min(c(2, length(tmp)))]
Param$Ratios.Contamination.Rates <- KontRt
if (Param$Ratios.Thresholds == "% of intra-sample group ratios") {
  Param$Ratios.Contamination.Rates <- Param$Ratios.Contamination.Rates/100
}
# Apply defaults to non-UI parameters
g <- grep("^TF_((TRUE)|(FALSE))$", Param[1,])
for (w in g) { Param[[w]] <- as.logical(gsub("^TF_", "", Param[[w]])) }
g <- grep("^((TRUE)|(FALSE))$", Param[1,])
for (w in g) { Param[[w]] <- as.logical(Param[[w]]) }
#}
