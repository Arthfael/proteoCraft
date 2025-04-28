#### Code chunk - Define analysis parameters
# Create first PCA to check on sample relationships
if ((length(MQ.Exp) > 1)||(LabelType == "Isobaric")) { # Should be always TRUE
  source(parSrc)
  data <- ev
  colnames(data)[which(colnames(data) == "MQ.Exp")] <- "Parent sample"
  data <- data[which(data$Reverse != "+"),]
  data <- data[which((is.na(data$"Potential contaminant"))|(data$"Potential contaminant" != "+")),]
  if (LabelType == "Isobaric") {
    kol <- grep(paste0(topattern(ev.ref["Original"]), "[0-9]+$"), colnames(data), value = TRUE)
  } else { kol <- ev.col["Original"] }
  w <- which(rowSums(data[, kol, drop = FALSE], na.rm = TRUE) > 0)
  data <- data[w,]
  if (!"Fraction" %in% colnames(data)) { data$Fraction <- 1 }
  Fraction <- sort(unique(data$Fraction), decreasing = FALSE)
  Experiment <- Exp
  kols <- c("Parent sample", "Fraction", "Experiment")
  if (LabelType == "Isobaric") {
    X <- "Label"
    kols <- c("Fraction", "Parent sample", "Experiment") # The order matters!
    tst <- sapply(kols, function(x) { length(unique(data[[x]])) })
    w1 <- which(tst > 1)
    w2 <- which(tst >= 1) 
    if (length(w1)) { Y <- kols[w1[1]] } else { Y <- kols[w2[1]] }
  }
  if (LabelType == "LFQ") {
    kols <- c("Parent sample", "Fraction", "Experiment") # The order matters!
    tst <- sapply(kols, function(x) { length(unique(data[[x]])) })
    w1 <- which(tst > 1)
    w2 <- which(tst >= 1) 
    X <- kols[w1[1]]
    if (length(w1) > 1) { Y <- kols[w1[2]] } else { Y <- kols[w2[2]] }
  }
  kols <- kols[which(!kols %in% c(X, Y))]
  ReportCalls <- AddSpace2Report()
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              "body_add_fpar(Report, fpar(ftext(\"PSMs-level PCA plot:\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
  ReportCalls$Calls <- append(ReportCalls$Calls, list())
  dir <- paste0(wd, "/Dimensionality red. plots/PCA")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  LRepCalls <- length(ReportCalls$Calls)
  lsKl <- c("Modified sequence", Y)
  if (LabelType == "LFQ") { lsKl <- c(lsKl, X) }
  ls <- lapply(lsKl, function(kl) { data[[kl]] })
  tmp <- do.call(paste, c(data[, lsKl], sep = "---"))
  if (LabelType == "Isobaric") {
    kol2 <- gsub(topattern(ev.ref["Original"]), "", kol)
    data2 <- as.data.table(data[, kol])
    colnames(data2) <- kol
    data2$Group <- tmp
    data2 <- data2[, lapply(.SD, sum, na.rm = TRUE), keyby = Group]
    colnames(data2) <- c("Group", kol)
    data2 <- as.data.frame(data2)
    data2[, kol2] <- log10(data2[, kol])
    data2 <- data2[, which(!colnames(data2) %in% kol)]
    data2[, lsKl] <- data[match(data2$Group, tmp), lsKl]
    data2$Group <- NULL
  }
  if (LabelType == "LFQ") {
    if (X == "Parent sample") { kol2 <- get("MQ.Exp") } else { kol2 <- get(X) }
    data2 <- data.table(Intensity = data[[ev.col["Original"]]], Group = tmp)
    data2 <- data2[, list(`log10(Intensity)` = sum(Intensity, na.rm = TRUE)),
                   keyby = Group]
    data2$`log10(Intensity)` <- log10(data2$`log10(Intensity)`)
    data2 <- as.data.frame(data2)
    data2[, lsKl] <- data[match(data2$Group, tmp), lsKl]
    data2$Group <- NULL
    data <- spread(data2, X, "log10(Intensity)")
  }
  data <- data[, kol2]
  w <- which(is.na(data), arr.ind = TRUE)
  if (nrow(w)) {
    groups <- rep(1, length(MQ.Exp))
    temp <- Data_Impute2(data, groups, is.log = FALSE)
    data <- temp$Imputed_data
  }
  pcA <- prcomp(t(data[, kol2]), scale. = TRUE)
  if (length(pcA$rotation)) {
    scoresA <- as.data.frame(pcA$x)
    if ("PC2" %in% colnames(scoresA)) {
      scoresA$Sample <- rownames(scoresA)
      rownames(scoresA) <- NULL
      pvA <- round(100*(pcA$sdev)^2 / sum(pcA$sdev^2), 0)
      pvA <- pvA[which(pvA > 0)]
      pvA <- paste0("Original: ", paste(sapply(1:length(pvA), function(x) {
        paste0("PC", x, ": ", pvA[x], "%")
      }), collapse = ", "))
      scoresA$Label <- scoresA$Sample
      m <- match(scoresA$Sample, Exp.map$MQ.Exp)
      if (sum(is.na(m))) {
        warning("Mapping samples through MQ.Exp to sample groups failed, check code!\nMapping colors to samples instead of sample groups...")
        scoresA$Colour <- scoresA$Sample
      } else {
        tmp <- Exp.map[m, Factors[which(Factors != "Replicate")]]
        tmp <- tmp[, which(sapply(colnames(tmp), function(x) { length(unique(tmp[[x]])) > 1 })), drop = FALSE]
        scoresA$Colour <- do.call(paste, c(tmp, sep = " "))
      }
      ttl <- "PCA plot - Samples (PSMs-level)"
      plot <- ggplot(scoresA) +
        geom_point(aes(x = PC1, y = PC2, colour = Colour)) +
        scale_color_viridis_d(begin = 0.25) +
        coord_fixed() + theme_bw() +
        geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
        ggtitle(ttl, subtitle = pvA) +
        geom_text_repel(aes(x = PC1, y = PC2, label = Label, colour = Colour),
                        size = 2.5, show.legend = FALSE)
      #poplot(plot)
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 20, height = 20, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 20, height = 20, units = "in")
      ReportCalls <- AddPlot2Report(Space = FALSE)
      Symb <- "circle"
      # Custom color scale
      scoresA$`Samples group` <- factor(scoresA$Colour)
      if ("PC3" %in% colnames(scoresA)) {
        plot_lyPSMsPCA <- plot_ly(scoresA, x = ~PC1, y = ~PC2, z = ~PC3,
                                  text = ~Label, type = "scatter3d", mode = "markers",
                                  color = ~`Samples group`, colors = "viridis",
                                  symbol = I(Symb))
      } else {
        plot_lyPSMsPCA <- plot_ly(scoresA, x = ~PC1, y = ~PC2,
                                  text = ~Label, type = "scatter", mode = "markers",
                                  color = ~`Samples group`, colors = "viridis",
                                  symbol = I(Symb))
      }
      plot_lyPSMsPCA %<o% layout(plot_lyPSMsPCA, title = ttl)
      saveWidget(plot_lyPSMsPCA, paste0(dir, "/", ttl, ".html"), selfcontained = TRUE)
      #system(paste0("open \"", dir, "/", ttl, ".html"))
    } else {
      msg <- "Not enough valid data to draw a PSM-level PCA plot!"
      ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
    }
  }
  ReportCalls$Calls[[LRepCalls]] <- append(ReportCalls$Calls[[LRepCalls]],
                                           "body_add_par(Report, \"\", style = \"Normal\")")
} else {
  stop("Uh, I think you have the wrong analysis pipeline here...\nwhere are my sample groups and replicates!?!?!")
}
#
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
pvalue.col %<o% c(StudentRoot, WelchRoot, modRoot, permRoot, samRoot, odpRoot, lrtRoot)
names(pvalue.col) <- sapply(pvalue.col, function(x) { unlist(strsplit(x, "\\.|\\'|\\ "))[1] })
ParamFls <- c(paste0(wd, "/Parameters.csv"),
              paste0(libPath, "/extData/Parameters_template.csv"))
ParamFl <- ParamFls[1]
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
    Param$Label.Purities.file <- "CHOOSEFL"
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
mxRp <- max(Exp.map$Replicate)
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
  CytoScExe <- sapply(tmp, function(x) { grep("Cytoscape\\.exe$", list.files(x, recursive = TRUE), value = TRUE) })
  CytoScExe <- paste0(tmp, "/", CytoScExe)
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
if (!exists("minInt")) { minInt <- 100 }
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
#
mnFct <- c("Experiment", "Replicate")
if (WorkFlow == "TIMECOURSE") { mnFct <- c(mnFct, "Time.point") }
l <- length(mnFct)
coreNms <- setNames(c("Ratios.Groups.Ref.Aggregate.Level", "Norm.Groups", "Ratios.Groups", "Batch.correction", "GO.enrichment.Ref.Aggr"),
                    c(paste0("Individual Sample___Combination of Factors required to distinguish individual samples./nMust include ",
                             paste(mnFct[1:(l-1)], collapse = ", "), " and ", mnFct[l], "."),
                      "Normalisation groups___Factor(s) defining groups of samples to normalize to each-other. Note that which, if any, normalisations apply will be defined further down",
                      "Ratio groups___Factor(s) defining comparison groups, i.e. groups of samples, including at least some References (i.e. Controls), to compare to each others./nMust include Experiment, cannot include Replicate.",
                      "Batch___Optional: Factor(s) defining batches used for sva::ComBat-based correction.",
                      "GO enrichment___Optional: Only protein groups with at least one valid value in corresponding samples will be used as references for GO-terms enrichment tests."))
for (nm in coreNms) { if (!nm %in% colnames(Param)) { Param[[nm]] <- "Exp" } }
wMp <- c(which(colnames(Param) == coreNms[1]),
         which(colnames(Param) == coreNms[2]),
         which(colnames(Param) == coreNms[3]),
         which(colnames(Param) == coreNms[4]),
         which(colnames(Param) == coreNms[5]),
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
appNm <- paste0(dtstNm, " - Parameters")
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  titlePanel(tag("u", "Parameters"),
             appNm),
  h2(dtstNm), 
  tags$hr(style = "border-color: black;"),
  br(),
  h4(strong("Factors")),
  span("Here we map diverse actions to experimental Factors:"),
  fluidRow(column(4, withSpinner(uiOutput("FactMappings"))),
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
    column(3, numericInput("minInt",
                           "Exclude PSMs with intensity lower than...",
                           minInt,
                           0,
                           .Machine$double.xmax,
                           width = "100%")),
    column(3, checkboxInput("Impute", "Impute missing peptides-level values?", Impute, "100%")),
    column(3,
           checkboxInput("Update_Prot_matches", paste0("Update ", names(SearchSoft), "'s original protein-to-peptides assignments?"), Update_Prot_matches, "100%"),
           bsTooltip("Update_Prot_matches",
                     "Checking assignments may result in removal of some identifications. It is nonetheless recommended because we have observed occasional inconsistent peptides-to-protein assignments with some search software.",
                     placement = "right", trigger = "hover", options = list(container = "body")),
           withSpinner(uiOutput("ReloadMatches"))
    )),
  tags$hr(style="border-color: black;"),
  withSpinner(uiOutput("IsobarCorr")),
  withSpinner(uiOutput("Norm")),
  br(),
  # Quantitation
  ## Choice of algorithm + Proteomics ruler
  tags$hr(style="border-color: black;"),
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
  h5(strong("T-test(s)")),
  fluidRow(
    column(2,
           radioButtons("TwoSided", "Fold changes: test...", c("Both directions", "Up-only", "Down-only"),
                        TwoSidedDeflt, FALSE, "100%"),
           checkboxInput("Mirror", "Revert fold changes on plots? (default: log fold change = log2(Sample/Reference); revert: log2(Reference/Sample))",
                         Param$Mirror.Ratios, "100%")),
    column(4,
           h5(strong(" -> Volcano plot: select default variant")),
           radioButtons("TtstPval", "", names(pvalue.col), "Moderated", TRUE, "100%"),
           h6(em(" - Welch's t-test is a modified form of Student's original version which is more robust to variance inequality.")),
           h6(em(" - Moderated t-test (limma package): re-samples individual row variances using global dataset variance to provide a more robust estimate.")),
           h6(em(" - Permutation t-test (coin package): based on permutations of the samples from each group.")),
           h6(em(" - SAM's modified t-test (siggenes package): corrects for poor variance estimates for very reproducible data by optimizing and adding a small constant s0 to the denominator of the test.")),
           h6(em(" - LRT (Likelihood Ratio Test, edge package): tests for the ratio of the likelihoods of observed data under two models (explanatory variable has an effect /vs/ no effect)")),
           h6(em(" - ODP (Optimal Discovery Procedure, Storey et al., 2007 - edge package): uses all relevant information from all genes in order to test each one for differential expression; has been demonstrated to have optimal power.")),
           if (scrptTypeFull == "withReps_PG_and_PTMs") {
             # Not available for now for the peptides-only workflow,
             # but this could be done by turning the code for this into an app which is called individually for each PTM class.
             # This would mean allowing for a different selection for each PTM class.
             # We will anyway want to allow for multiple tests selected in the future, so go out of the "choose one t-test" variant approach...
             h6(em("(you can change which one will be used for Volcano plots later, after we compare each test's power)")),
             checkboxInput("useSAM_thresh", "For Student's t-test, plot SAM-based curved significance thresholds?", useSAM_thresh, "100%"),
             br(),
           }
           checkboxInput("useSAM_thresh", "For Student's t-test, plot SAM-based curved significance thresholds?", useSAM_thresh, "100%"),
           br(),
           h5(strong(" -> ANOVA (moderated F-test - limma)")),
           checkboxInput("Ftest", "run?", fTstDflt, "100%"),
           br(),
           uiOutput("sntXprs")),
    column(2,
           h5("Benjamini-Hochberg FDR thresholds"),
           withSpinner(uiOutput("FDR")),
           numericInput("FDR", "", 0.01, 0, 1),
           actionButton("AddFDR", "Add threshold"),
           actionButton("RemvFDR", "Remove threshold")),
    column(2,
           radioButtons("typeOfThresh", "Ratios thresholds: use...", threshOpt, threshDflt, FALSE, "100%"),
           numericInput("RatCont", "Threshold value = ", KontRt, 0, 100, 0.001, width = "100%"),
           selectInput("RatContGrp",
                       "Control ratios are grouped by:",
                       KontGrps,
                       KontGrp,
                       width = "100%"))
  ),
  tags$hr(style="border-color: black;"),
  withSpinner(uiOutput("GO")),
  tags$hr(style="border-color: black;"),
  withSpinner(uiOutput("CytoScape")),
  tags$hr(style="border-color: black;"),
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
  fluidRow(column(2, checkboxInput("Amica", "Write tables compatible with https://bioapps.maxperutzlabs.ac.at/app/amica", Param$Amica, "100%"))),
  br(),
  tags$hr(style="border-color: black;"),
  checkboxInput("AdvOptOn", "Advanced options", tstAdvOpt),
  withSpinner(uiOutput("AdvOpt")),
  br(),
  actionButton("saveBtn", "Save"),
  br()
)
server <- function(input, output, session) {
  # Initialize variables to create in main environment
  PARAM <- reactiveVal(Param)
  m4Quant <- reactiveVal(Mod4Quant)
  m2Xclud <- reactiveVal(Mod2Xclud)
  ADVOPT <- reactiveVal(tstAdvOpt)
  #
  # Dynamic UI
  # Map Parameters to Factors
  output$PSMsPCA <- renderPlotly(plot_lyPSMsPCA)
  output$FactMappings <- renderUI({
    # lst <- list()
    # for (w in wMp) {
    #   Opt <- Factors
    #   dflt <- dfdflt <- c("Experiment", "Replicate")[(colnames(Param)[w] == "Batch.correction")+1]
    #   if (!Param[1, w] %in% c("AUTOFACT", "MAP2FACTS")) {
    #     dflt <- Factors[match(unlist(strsplit(Param[1, w], ";")), names(Factors))]
    #     if ((length(dflt) == 1)&&(is.na(dflt))) { dflt <- dfdflt }
    #   }
    #   if (colnames(Param)[w] == "Ratios.Groups.Ref.Aggregate.Level") {
    #     dflt <- unique(c(dflt, mnFct))
    #   }
    #   if (colnames(Param)[w] == "Ratios.Groups") {
    #     Opt <- Factors[which(Factors != "Replicate")]
    #     dflt <- dflt[which(dflt != "Replicate")]
    #   }
    #   lbl <- gsub("\\.", " ", colnames(Param)[w])
    #   if (colnames(Param)[w] %in% coreNms) {
    #     lbl <- unlist(strsplit(names(coreNms)[match(colnames(Param)[w], coreNms)], "___"))
    #   } else { lbl <- c(gsub("\\.", " ", colnames(Param)[w]), Param_Help[w]) }
    #   names(Opt) <- NULL
    #   names(dflt) <- NULL
    #   blck <- list(list(br()),
    #                tags$table(
    #                  tags$tr(width = "80%",
    #                          tags$td(width = "25%",
    #                                  div(strong(lbl[1]))),
    #                          tags$td(width = "55%",
    #                                  selectInput(colnames(Param)[w],
    #                                              "",
    #                                              Opt,
    #                                              dflt,
    #                                              TRUE,
    #                                              TRUE)),
    #                          #addTooltip(session, colnames(Param)[w], Param_Help[w], "bottom", "hover", list(container = "body"))
    #                  )
    #                ))
    #   lbl2 <- unlist(strsplit(lbl[2], "/n"))
    #   for (lbl2a in lbl2) { blck <- append(blck, list(span(em(lbl2a)), br())) }
    #   if (colnames(Param)[w] == "Ratios.Groups") {
    #     dflt <- Param$Ratios.Groups_Nested
    #     if (!is.logical(dflt)) { dflt <- WorkFlow != "Regulation" }
    #     blck <- append(blck, list(radioButtons("IsNested", "Nested design? (i.e. are replicates paired?)",
    #                                            c(TRUE, FALSE), dflt, TRUE)))
    #   }
    #   blck <- append(blck, list(br()))
    #   lst <- append(lst, blck)
    # }
    # return(lst)
    lstFct
  })
  output$RatioGroups <- renderUI({
    em(paste(unique(apply(ExpMap[, input[["Ratios.Groups"]], drop = FALSE], 1, paste, collapse = "___")), collapse = " / "))
  })
  # Factors
  sapply(wMp, function(w) {
    observeEvent(input[[colnames(Param)[w]]], {
      if (colnames(Param)[w] == "Ratios.Groups.Ref.Aggregate.Level") {
        l <- length(input[[colnames(Param)[w]]])
        if ((l < 3)||(sum(!c("Experiment", "Replicate") %in% input[[colnames(Param)[w]]]))) {
          shinyjs::disable("saveBtn")
        } else { shinyjs::enable("saveBtn") }
      }
      Par <- PARAM()
      Par[colnames(Par)[w]] <- paste0(input[[colnames(Param)[w]]], collapse = ";")
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
  # Labels purity correction - only for isobarically labelled samples
  output$IsobarCorr <- renderUI({
    if (LabelType == "Isobaric") {
      lst <- list(shinyFilesButton("PurityFl", "Browse", em(paste0(IsobarLab, " purity table")), "", FALSE),
                  br())
    } else {
      lst <- list(HTML(""))
    }
    return(list(lst))
  })
  # Normalisations
  output$Norm <- renderUI({
    lst <- list(list(fluidRow(column(1, h4(strong("Normalisations"))),
                              column(2, checkboxInput("Norm", "Normalize data",
                                                      sum(Param$Norma.Ev.Intens,
                                                          Param$Norma.Pep.Intens,
                                                          Param$Norma.Prot.Ratio) > 0, "100%"))),
                     br(),
                     fluidRow(column(1, h5(strong(" -> PSMs-level: "))),
                              column(3, checkboxInput("evLM", "Levenberg-Marquardt",
                                                      Param$Adv.Norma.Ev.Intens, "100%"))),
                     br(),
                     fluidRow(column(1, h5(strong(" -> Peptidoforms-level: "))),
                              column(2, radioButtons("pepShape", "Corrections",
                                                     c("none", "vsn", "loess"),
                                                     shpDflt, TRUE, "100%")),
                              column(2, checkboxInput("pepLM", "Levenberg-Marquardt",
                                                      Param$Adv.Norma.Pep.Intens, "100%"))),
                     br()))
    if ((LabelType == "Isobaric")&&(length(Iso) > 1)) {
      lst <- append(lst, list(fluidRow(column(2, checkboxInput("IRS", "IRS", Param$Norma.Pep.Intens.IRS, "100%")))))
      dfltChan <- unlist(strsplit(Param$Norma.Pep.Intens.IRS_Ref_channels, ";"))
      dfltChan <- dfltChan[which(dfltChan %in% get(IsobarLab))]
      if (length(dfltChan) != length(Iso)) { dfltChan <- rep(get(IsobarLab)[1], length(Iso)) }
      dfltChan <- reactiveVal(dfltChan)
      for (i in length(Iso)) {
        m <- unique(Exp.map[which(Exp.map$Isobaric.set == Iso[i])],)
        lst <- append(lst, list(fluidRow(column(2, h5(strong(em(paste0(" -> ",
                                                                       paste(dfltChan(), collapse = " / "),
                                                                       ": "))))),
                                         column(2, pickerInput(paste0("intRef", Iso[i]),
                                                               paste0("- Isobaric set ", Iso[i]),
                                                               m$Isobaric.label, dfltChan()[i], TRUE,
                                                               pickerOptions(title = "Search me",
                                                                             `live-search` = TRUE,
                                                                             actionsBox = TRUE,
                                                                             deselectAllText = "Clear search",
                                                                             showTick = TRUE)
                                         ))),
                                br()))
      }
    }
    if (scrptTypeFull == "withReps_PG_and_PTMs") {
      lst <- append(lst, list(fluidRow(column(1, h5(strong(" -> Protein Groups-level: "))),
                                       column(2, checkboxInput("prtLM", "Levenberg-Marquardt",
                                                               Param$Adv.Norma.Prot.Intens, "100%")),
                                       br()),
                              fluidRow(column(2,
                                              pickerInput("prt2PrtSlct", " -> To selected proteins ",
                                                          protHeads, prt2PrtSlctDflt, TRUE,
                                                          pickerOptions(title = "Search me",
                                                                        `live-search` = TRUE,
                                                                        actionsBox = TRUE,
                                                                        deselectAllText = "Clear search",
                                                                        showTick = TRUE)),
                                              br())),
                              if (WorkFlow == "BIOID") {
                                fluidRow(column(2,
                                                checkboxInput("prt2Biot", " -> To all biotinylated-proteins: ",
                                                              Param$Norma.Prot.Ratio.to.Biot, "100%")))
                              }
      ))
    }
    return(lst)
  })
  #
  # GO
  output$GO <- renderUI({
    lst <- list(list(br()))
    if (Annotate) {
      lst <- list(
        list(fluidRow(h4(strong("\tGO terms enrichment"))),
             fluidRow(column(1, checkboxInput("GOenrich", "GO enrichment", goDflt, "100%")),
                      column(2, pickerInput("GO.tabs", "GO terms of interest", allGO, dftlGO, TRUE,
                                            pickerOptions(title = "Search me",
                                                          `live-search` = TRUE,
                                                          actionsBox = TRUE,
                                                          deselectAllText = "Clear search",
                                                          showTick = TRUE))),
                      column(1, checkboxInput("GO2Int", "Use GO terms to define list of proteins of interest?",
                                              Param$GO.terms.for.proteins.of.interest , "100%"))))
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
  #
  # Event observers
  # Optional input files
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
                        em("Useful when e.g. samples express from a custom construct to which matching peptides should be assigned in priority."), br(),
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
  output$AdvOpt <- updtOptOn(FALSE)
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
  observeEvent(input$pepShape, {
    Par <- PARAM()
    Par$Norma.Pep.Intens.Shape <- input$pepShape
    PARAM(Par)
  })
  observeEvent(input$pepLM, {
    Par <- PARAM()
    Par$Adv.Norma.Pep.Intens <- input$pepLM
    PARAM(Par)
  })
  if ((LabelType == "Isobaric")&&(length(Iso) > 1)) {
    observeEvent(input$IRS, {
      Par <- PARAM()
      Par$Norma.Pep.Intens.IRS <- input$IRS
      PARAM(Par)
    })
    for (i in 1:length(Iso)) {
      observeEvent(input[[paste0("intRef", Iso[i])]], {
        Par <- PARAM()
        tmp <- dfltChan()
        tmp[i] <- input[[paste0("intRef", Iso[i])]]
        dfltChan(tmp)
        Par$Norma.Pep.Intens.IRS_Ref_channels <- paste(tmp, collapse = ";")
        PARAM(Par)
      }, ignoreNULL = FALSE)
    }
  }
  if (scrptTypeFull == "withReps_PG_and_PTMs") {
    observeEvent(input$prtLM, {
      Par <- PARAM()
      Par$Adv.Norma.Prot.Intens <- input$prtLM
      PARAM(Par)
    })
    observeEvent(input$prt2PrtSlct, {
      Par <- PARAM()
      Par$Norma.Prot.Ratio.to.proteins <- paste(db$`Protein ID`[dbOrd][match(input$prt2PrtSlct, protHeads)],
                                                collapse = ";")
      PARAM(Par)
    }, ignoreNULL = FALSE)
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
          Par <- PARAM()
          Par$Label.Purities.file <- normalizePath(tmp, winslash = "/")
          PARAM(Par)
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
    Par$GO.terms.for.proteins.of.interest <- input$GO2Int
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
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
eval(parse(text = runApp))
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
#if (!Param$Param_suppress_UI) {
# Statistical tests
#dlg_message("Remember: replace this app with a \"Statistics: contrasts editor\" app!", "ok")
#
vpal <- unlist(strsplit(Param$Volcano.plots.Aggregate.Level, ";"))
# We should select a reference level per ratio group, not ratio reference group!
#rrg <- unlist(strsplit(Param$Ratios.Ref.Groups, ";"))
rg <- unlist(strsplit(Param$Ratios.Groups, ";"))
if (length(unique(Exp.map$Experiment)) == 1) {
  if (length(vpal) > 1)  { vpal <- vpal[which(vpal != "Exp")] }
  #if (length(rrg) > 1)  { rrg <- rrg[which(rrg != "Exp")] }
  if (length(rg) > 1)  { rg <- rg[which(rg != "Exp")] }
}
kolVPAL <- Factors[vpal]
#kolRRG <- Factors[rrg]
kolRG <- Factors[rg]
factLevComb1 <- apply(Exp.map[, kolVPAL, drop = FALSE], 1 , paste, collapse = " ")
#factLevComb2 <- apply(Exp.map[, kolRRG, drop = FALSE], 1 , paste, collapse = " ")
factLevComb2 <- apply(Exp.map[, kolRG, drop = FALSE], 1 , paste, collapse = " ")
refTst <- ("Reference" %in% colnames(Exp.map))&&("logical" %in% class(Exp.map$Reference))&&(sum(Exp.map$Reference))&&(sum(!Exp.map$Reference))
rfLev <- data.frame("Group" = unique(factLevComb2))
rfLev$"All levels" <- lapply(rfLev$Group, function(x) {
  unique(factLevComb1[which(factLevComb2 == x)])
})
rfLev$"Reference level" <- sapply(rfLev$"All levels", function(x) { #x <- unlist(rfLev$"All levels"[2])
  if (refTst) {
    w <- which(Exp.map$Reference[match(x, factLevComb1)])
    if (length(w)) { w <- w[1] } else { w <- 1 }
    rs <- x[w]
  } else {
    rs <- factLevComb1[w[1]]
  }
  return(rs)
})
makeSec <- FALSE
# if (Param$F.test) {
#   Factors3 <- Factors2[which(sapply(Factors2, function(Fact) {
#     length(FactorsLevels[[Fact]]) > 1
#   }))]
#   if (WorkFlow == "TIMECOURSE") { Factors3 <- unique(c(Factors3, "Time.point")) }
#   makeSec <- length(Factors3) > 1
# }
wdth <- paste0(30*max(c(nchar(unlist(rfLev$`All levels`)), 2)), "px")
appNm <- paste0(dtstNm, " - Stat-tests")
ui <- fluidPage(
  useShinyjs(),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  titlePanel(tag("u", "Statistical test - define comparisons"), 
             appNm),
  h2(dtstNm), 
  br(),
  actionButton("saveBtn", "Save"),
  br(),
  fluidPage(
    h3(" -> T-tests: define reference factors combination per comparisons group."),
    h5(em("Select a single reference level sample group per comparisons group.")),
    h5(em("Samples from this samples group will serve as the denominator of ratios and the controls for all t-tests performed within the comparisons group.")),
    br(),
    #withSpinner(
    DTOutput("refLevels")#)
    ,
    br(),
    br()
  ),
  # if (Param$F.test) {
  #   sidebarLayout(
  #     sidebarPanel(
  #       h3(" -> F-tests: design analyses"),
  #       h5(em("For F-tests, you can define multiple analyses, each with simple and/or double contrasts to different references.")),
  #       h5(em("For each analysis, select reference levels for each Factor relevant to the desired contrasts (comparisons) then add Analysis.")),
  #       h5(em("N.B.: Only combinations of levels corresponding to at least 1 valid sample group are allowed.")),
  #       br(),
  #       actionButton("addAnalysis", "Add new analysis"),
  #       withSpinner(uiOutput("NewAnalysis")),
  #       span(uiOutput("Msg"), style = "color:red"),
  #     ),
  #     mainPanel(
  #       h3("F-test analyses"),
  #       br(),
  #       uiOutput("Analyses")
  #     )
  #   )
  # },
  br()
)
#
NuAn_tmplt <- data.frame(Analysis = "", Primary = "", Primary_Ref = "")
if (makeSec) { NuAn_tmplt$Secondary <- NuAn_tmplt$Secondary_Ref <- "" }
rfLev2 <- rfLev[, c("Group", "Reference level")]
rfLev2$"Reference level" <- sapply(1:nrow(rfLev2), function(i) {
  as.character(selectInput(paste0("Ref_", as.character(i)),
                           "",
                           rfLev$`All levels`[[i]],
                           rfLev$`Reference level`[i],
                           selectize = FALSE,
                           width = wdth))
})
server <- function(input, output, session) {
  output$refLevels <- renderDT({ rfLev2 },
                               FALSE,
                               escape = FALSE,
                               selection = "none",
                               editable = TRUE,
                               rownames = FALSE,
                               options = list(dom = 't',
                                              paging = FALSE,
                                              ordering = FALSE
                               ),
                               callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
  # sapply(1:nrow(rfLev2), function(x) {
  #   id <- paste0("Ref_", as.character(x))
  #   observeEvent(input[[id]],
  #                {
  #                  rfLev$`Reference level`[x] <- input[[id]]
  #                  assign("rfLev", rfLev, envir = .GlobalEnv)
  #                })
  # })
  # if (Param$F.test) { # Initialize variables to create in main environment
  #   #
  #   # Parse F-test Parameters
  #   fFact <- fRef <- FALSE
  #   #Fact <- Factors2[which(!Factors2 %in% c("Experiment", "Replicate"))]
  #   #w <- sapply(Fact, function(x) { length(FactorsLevels[[x]]) > 1 })
  #   #Fact[w]
  #   #tmp <- setNames(substr(Fact, 1, 3), Fact)
  #   #Param$F.test_factors <- paste0(paste(tmp[1:2], collapse = "___"), "_;_", tmp[3])
  #   #Param$F.test_factors_ref <- paste0(paste(sapply(names(tmp)[1:2], function(x) { FactorsLevels[[x]][1] }), collapse = "___"), "_;_", FactorsLevels[[names(tmp)[3]]][1])
  #   if ("F.test_factors" %in% colnames(Param)) {
  #     F_Fact <- unique(unlist(strsplit(Param$F.test_factors, "_\\|_")))
  #     F_Fact <- F_Fact[which(F_Fact != "")]
  #     if (length(F_Fact)) {
  #       fFact <- TRUE
  #     }
  #   }
  #   if ("F.test_factors_ref" %in% colnames(Param)) {
  #     F_Ref <- unique(unlist(strsplit(Param$F.test_factors_ref, "_\\|_")))
  #     F_Ref <- F_Ref[which(F_Ref != "")]
  #     if (length(F_Ref)) {
  #       fRef <- TRUE
  #     }
  #   }
  #   if (fFact && fRef) {
  #     F_Fact <- strsplit(F_Fact, "_;_")
  #     stopifnot(max(sapply(F_Fact, length)) <= 3) # We can have at most 3 factors aggregates: factors for primary and secondary contrasts, and blocking factors for nested designs.
  #     F_Fact <- lapply(F_Fact, function(x) { x[1:min(c(2, length(x)))] }) # We will add nesting later for simplicity using the UI. 
  #     F_Fact <- lapply(F_Fact, function(x) {
  #       x <- strsplit(x, "___")
  #       x <- x[1:2]
  #       return(x)
  #     })
  #     F_Fact <- setNames(F_Fact, paste0("Analysis_", 1:length(F_Fact)))
  #     F_Ref <- strsplit(F_Ref, "_;_")
  #     stopifnot(max(sapply(F_Ref, length)) <= 2) # We can have at most 2 here: no blocking factors.
  #     F_Ref <- lapply(F_Ref, function(x) { x[1:min(c(2, length(x)))] })
  #     F_Ref <- lapply(F_Ref, function(x) {
  #       x <- strsplit(x, "___")
  #       x <- x[1:2] 
  #       return(x)
  #     })
  #     F_Ref <- setNames(F_Ref, paste0("Analysis_", 1:length(F_Ref)))
  #     tst <- length(F_Fact) == length(F_Ref)
  #     if (tst) {
  #       tst <- sum(sapply(1:length(F_Fact), function(x) { length(F_Fact[[x]]) != length(F_Ref[[x]]) })) == 0
  #       if (tst) {
  #         tmp <- substr(Factors, 1, 3)
  #         tst <- sum(unlist(sapply(1:length(F_Fact), function(x) {
  #           sapply(1:length(F_Ref[[x]]), function(y) {
  #             if (length(F_Fact[[x]][[y]]) > 0) {
  #               sapply(1:length(F_Fact[[x]][[y]]), function(z) {
  #                 sum(!F_Ref[[x]][[y]][[z]] %in% FactorsLevels[[Factors[match(substr(F_Fact[[x]][[y]][[z]],1 , 3), tmp)]]])
  #               })  
  #             } else { 0 }
  #           })
  #         }))) == 0
  #         if (!tst) { fFact <- fRef <- FALSE }
  #       } else { fFact <- fRef <- FALSE }
  #     } else { fFact <- fRef <- FALSE }
  #   }
  #   if (fFact && fRef) {
  #     nms <- names(F_Fact)
  #     F_Analyses <- data.frame(Analysis = nms, row.names = nms)
  #     F_Analyses$Primary <- lapply(F_Fact, function(x) { x[[1]] })
  #     F_Analyses$Primary_Ref <- lapply(F_Ref, function(x) { x[[1]] })
  #     if (makeSec) {
  #       F_Analyses$Secondary <- lapply(F_Fact, function(x) { x[[2]] })
  #       F_Analyses$Secondary_Ref <- lapply(F_Ref, function(x) { x[[2]] })
  #     }
  #   } else {
  #     F_Analyses <- data.frame(Analysis = "", Primary = "", Primary_Ref = "")
  #     if (makeSec) { F_Analyses$Secondary <- F_Analyses$Secondary_Ref <- "" }
  #     F_Analyses <- F_Analyses[integer(0),]
  #   }
  #   #
  #   F_An <- reactiveVal(F_Analyses)
  #   nr <- reactiveVal(nrow(F_Analyses))
  #   mxNr <- reactiveVal(nrow(F_Analyses))
  #   #
  #   updtAnUI <- function(reactive = TRUE) {
  #     # Update UI
  #     if (reactive) { FA <- F_An() } else { FA <- F_Analyses }
  #     return(renderUI({
  #       if (nrow(FA)) {
  #         txts <- apply(FA, 1, function(r) {
  #           #r <- FA[1,]
  #           pr <- unlist(r[[2]])
  #           prRf <- unlist(r[[3]])
  #           txt <- paste0("-> Primary: ", paste(sapply(1:length(pr), function(x) {
  #             paste0(pr[x], " = ", prRf[x])
  #           }), collapse = ", "))
  #           if (makeSec) {
  #             sc <- unlist(r[[4]])
  #             if (length(sc)) {
  #               scRf <- unlist(r[[5]])
  #               txt <- paste0(txt, "_|_",
  #                             paste0("-> Secondary: ", paste(sapply(1:length(sc), function(x) {
  #                               paste0(sc[x], " = ", scRf[x])
  #                             }), collapse = ", ")))
  #             }
  #           }
  #           return(txt)
  #         })
  #         txts <- strsplit(txts, "_\\|_")
  #         lst <- vector("list", length(txts)*2)
  #         for (r in 1:length(txts)) {
  #           txt <- txts[[r]]
  #           if (length(txt) == 1) {
  #             lst[[r*2-1]] <- list(tags$table(
  #               tags$tr(width = "100%", tags$td(width = "100%", br())),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(paste0("Analysis ", r, ":")))),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(txt[[1]])))
  #             ))
  #           } 
  #           if (length(txt) == 2) {
  #             lst[[r*2-1]] <- list(tags$table(
  #               tags$tr(width = "100%", tags$td(width = "100%", br())),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(paste0("Analysis ", r, ":")))),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(txt[[1]]))),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(txt[[2]])))
  #             ))
  #           }
  #           lst[[r*2]] <- list(actionButton(paste0("Rmv", r), "Remove"))
  #         }
  #       } else { lst <- list(em("No analysis yet")) }
  #       return(lst)
  #     }))
  #   }
  #   # Create original remove analysis observers
  #   if (nrow(F_Analyses)) {
  #     sapply(1:nrow(F_Analyses), function(r) {
  #       observeEvent(input[[paste0("Rmv", r)]], {
  #         FA <- F_An()
  #         rws <- 1:nrow(FA)
  #         rws <- rws[which(rws != r)]
  #         FA <- FA[rws,]
  #         nr(nrow(FA))
  #         if (nr()) { FA$Analysis <- paste0("Analysis_", 1:nr()) } # Rename analyses
  #         F_An(FA)
  #         output$Analyses <- updtAnUI()
  #       })
  #     })
  #     output$Analyses <- updtAnUI(reactive = FALSE)
  #   }
  #   #
  #   # Box to define a new analysis
  #   dfltsPrim <- reactive(setNames(rep("Not used", length(Factors3)), Factors3))
  #   if (makeSec) { dfltsSec <- reactive(setNames(rep("Not used", length(Factors3)), Factors3)) }
  #   output$NewAnalysis <- renderUI({
  #     tags$table(
  #       tags$tr(width = "100%",
  #               tags$td(width = "20%", div(em(""))),
  #               tags$td(width = "40%", div(em("Primary contrasts"))),
  #               if (makeSec) { tags$td(width = "40%", div(em("Secondary contrasts"))) }
  #       ),
  #       lapply(Factors3, function(Fact){
  #         tags$tr(width = "100%",
  #                 tags$td(width = "20%", div(strong(Fact))),
  #                 tags$td(width = "40%", div(selectInput(paste0("Prim", Fact), "",
  #                                                        c("Not used", FactorsLevels[[Fact]]), dfltsPrim()[[Fact]]))),
  #                 if (makeSec) { tags$td(width = "40%", div(selectInput(paste0("Sec", Fact), "",
  #                                                                       c("Not used", FactorsLevels[[Fact]]),
  #                                                                       dfltsSec()[[Fact]]))) }
  #         )
  #       })
  #     )
  #   })
  #   # Observers for factor reference selection
  #   # These will change the selected value for other factors in the column, if incompatible
  #   EM <- Exp.map
  #   for (Fact in Factors3) { EM[which(is.na(EM[[Fact]])), Fact] <- "NA" }
  #   output$Msg <- renderUI({ em(" ") })
  #   sapply(Factors3, function(Fact) {
  #     sapply(c("Prim", "Sec")[1:(makeSec+1)], function(i) {
  #       observeEvent(input[[paste0(i, Fact)]], {
  #         m <- EM[which(as.character(EM[[Fact]]) == input[[paste0(i, Fact)]]),]
  #         facTrs <- Factors3[which(Factors3 != Fact)]
  #         k <- 0
  #         for (fcT in facTrs) {
  #           if (!input[[paste0(i, fcT)]] %in% c("Not used", m[[fcT]])) {
  #             k <- k+1
  #             updateTextInput(inputId = paste0(i, fcT), value = "Not used")
  #           }
  #         }
  #         output$Msg <- renderUI({ em(c(" ", "Invalid levels combination!")[(k > 0)+1]) })
  #       })
  #     })
  #   })
  #   # Parse new analysis + update UI
  #   # Template new Analysis row
  #   kol <- c("Primary", "Primary_Ref", "Secondary", "Secondary_Ref")[1:(2^(1+makeSec))]
  #   observeEvent(input$addAnalysis, {
  #     Prim <- setNames(sapply(Factors3, function(x) { input[[paste0("Prim", x)]] }), Factors3)
  #     #Prim <- setNames(rep(NA, length(Factors3)), Factors3) ; Prim[Factors3[1]] <- FactorsLevels[[Factors3[1]]][1]
  #     #
  #     Prim <- Prim[which(Prim != "Not used")]
  #     if (length(Prim)) {
  #       NuAn <- NuAn_tmplt
  #       NuAn$Primary <- list(names(Prim))
  #       NuAn$Primary_Ref <- list(Prim)
  #       FA <- F_An()
  #       #FA <- F_Analyses
  #       if (makeSec) {
  #         Sec <- setNames(sapply(Factors3, function(x) { input[[paste0("Sec", x)]] }), Factors3)
  #         #Sec <- setNames(rep(NA, length(Factors3)), Factors3) ; Sec[Factors3[2]] <- FactorsLevels[[Factors3[2]]][1]
  #         Sec <- Sec[which(Sec != "Not used")]
  #         NuAn$Secondary <- list(names(Sec))
  #         NuAn$Secondary_Ref <- list(Sec)
  #       }
  #       tst1 <- apply(NuAn[, kol], 1, function(x) {
  #         x <- unlist(x)
  #         names(x) <- NULL
  #         return(paste(x, collapse = "@"))
  #       })
  #       tst2 <- apply(FA[, kol], 1, function(x) {
  #         x <- unlist(x)
  #         names(x) <- NULL
  #         return(paste(x, collapse = "@"))
  #       })
  #       if (!tst1 %in% tst2) {
  #         FA <- rbind(FA, NuAn)
  #         nr(nrow(FA))
  #         FA$Analysis <- paste0("Analysis_", 1:nr()) #FA$Analysis <- paste0("Analysis_", 1:nrow(FA))
  #         F_An(FA)
  #         # Create new analysis remove observers
  #         if (nr() > mxNr()) {
  #           sapply((mxNr()+1):nr(), function(r) {
  #             observeEvent(input[[paste0("Rmv", r)]], {
  #               FA <- F_An()
  #               rws <- 1:nrow(FA)
  #               rws <- rws[which(rws != r)]
  #               FA <- FA[rws,]
  #               nr(nrow(FA))
  #               if (nr()) { FA$Analysis <- paste0("Analysis_", 1:nr()) } # Rename analyses
  #               F_An(FA)
  #               output$Analyses <- updtAnUI()
  #             })
  #           })
  #           mxNr(max(c(mxNr()), nr()))
  #         }
  #         # Update UI
  #         output$Analyses <- updtAnUI()
  #       }
  #     }
  #   })
  #   #
  # }
  observeEvent(input$saveBtn, {
    # if (Param$F.test) {
    #   F_Analyses <<- F_An()
    #   #print(F_Analyses)
    #   # Update parameters here!
    #   if (nrow(F_Analyses)) {
    #     PAR <- Param
    #     tmpFct <- data.frame(Prim = sapply(F_Analyses$Primary, function(x) { paste(substr(x, 1, 3), collapse = "___") }))
    #     tmpRf <- data.frame(Prim = sapply(F_Analyses$Primary_Ref, paste, collapse = "___"))
    #     if (makeSec) {
    #       tmpFct$Sec <- sapply(F_Analyses$Secondary, function(x) { paste(substr(x, 1, 3), collapse = "___") })
    #       tmpRf$Sec <- sapply(F_Analyses$Secondary_Ref, paste, collapse = "___")
    #     } else {
    #       tmpFct$Sec <- ""
    #       tmpRf$Sec <- ""
    #     }
    #     if (Param$Ratios.Groups_Nested) { tmpFct$Block <- "Rep" }
    #     PAR$F.test_factors <- paste(gsub("_;_$", "", do.call(paste, c(tmpFct, sep = "_;_"))), collapse = "_|_")
    #     PAR$F.test_factors_ref <- paste(gsub("_;_$", "", do.call(paste, c(tmpRf, sep = "_;_"))), collapse = "_|_")
    #     Param <<- PAR
    #   } else { Param$F.test <- FALSE }
    # }
    rfLev$`Reference level` <- sapply(1:nrow(rfLev2), function(x) {
      input[[paste0("Ref_", as.character(x))]]
    })
    assign("rfLev", rfLev, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
eval(parse(text = runApp))
#
wRf <- lapply(1:nrow(rfLev), function(i) { #i <- 1
  w <- which(factLevComb2 == rfLev$Group[i])
  return(w[which(factLevComb1[w] == rfLev$"Reference level"[i])])
})
Exp.map$Reference <- 1:nrow(Exp.map) %in% unlist(wRf)
tmp <- Exp.map
tmp$MQ.Exp <- sapply(tmp$MQ.Exp, paste, collapse = ";")
tst <- try(write.csv(tmp, file = ExpMapPath, row.names = FALSE), silent = TRUE)
if ("try-error" %in% class(tst)) {
  dlg_message(paste0("File \"", ExpMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  write.csv(tmp, file = ExpMapPath, row.names = FALSE)
}
tmp <- data.frame(x1 = colnames(Param), x2 = unlist(Param[1,]))
colnames(tmp) <- NULL
write.csv(tmp, paste0(wd, "/Parameters.csv"), row.names = FALSE)
#}
# if (Param$Param_suppress_UI) { # Candidate for deletion to simplify: the UI just works!
#   Param <- Param.load()
#   if (("Prot.Quant.Mod" %in% colnames(Param))&&(Param$Prot.Quant.Mod != "")&&(!is.na(Param$Prot.Quant.Mod))) {
#     warning("Parameter Prot.Quant.Mod is deprecated!")
#   }
#   if (("Prot.Quant.Mod.Excl" %in% colnames(Param))&&(Param$Prot.Quant.Mod.Excl != "")&&(!is.na(Param$Prot.Quant.Mod.Excl))) {
#     warning("Parameter Prot.Quant.Mod.Excl is deprecated!")
#   }
#   #Mod4Quant <- unlist(strsplit(Param$Prot.Quant.Mod, split = ";"))
#   ObjNm <- "Mod4Quant"
#   if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
#     msg <- "Which modifications are eligible for protein group quantitation?"
#     mods <- setNames(sapply(Modifs$`Full name`, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") }), Modifs$`Full name`)
#     pre <- mods[which(!grepl("phospho", Modifs$`Full name`, ignore.case = TRUE))]
#     tmp <- dlg_list(mods, pre, multiple = TRUE, title = msg)$res
#     if (length(tmp)) {
#       tmp <- Modifs$Mark[match(names(mods)[match(tmp, mods)], Modifs$`Full name`)]
#     }
#     assign(ObjNm, tmp)
#     AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
#     tmp <- AllAnsw[1,]
#     tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
#     tmp$Value <- list(get(ObjNm))
#     m <- match(ObjNm, AllAnsw$Parameter)
#     if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
#   }
#   # Note that modifications as indirectly inferred from data (when you use the Modifs object)
#   # do not currently for some cases reflect the exact context as defined in the search parameters.
#   # For instance, for "Gln->pyro-Glu" the position is "_" which means any N-terminus,
#   # but the modification as defined in the search is any N-terminal Q.
#   Mod2Xclud <- set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
#                             c("Mark", "Where"))
#   # Clustering method
#   ObjNm <- "KlustMeth"
#   if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
#     msg <- "Which clustering method do you want to use (used at samples and protein groups levels)?"
#     klustChoices <- c("hierarchical                                                                                                                                                                                ",
#                       "K-means                                                                                                                                                                           ")
#     tmp <- dlg_list(klustChoices, klustChoices[1], title = msg)$res
#     if (!length(tmp)) { tmp <- klustChoices[1] }
#     tmp <- c(2, 1)[match(tmp, klustChoices)] # 2,1: not a mistake
#     ObjNm %<c% tmp
#     AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
#     tmp <- AllAnsw[1,]
#     tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
#     tmp$Value <- list(get(ObjNm))
#     m <- match(ObjNm, AllAnsw$Parameter)
#     if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
#   }
# }
Param$FullQuant <- TRUE
FullQuant %<o% Param$FullQuant
protrul %<o% Param$ProtRul
Update_Prot_matches <- Param$Update_Prot_matches
Reuse_Prot_matches <- Param$Reuse_Prot_matches
Pep4QuantOpt %<o% c("Unique peptide IDs", "Razor peptide IDs", "Peptide IDs")
Pep4Quant %<o% Pep4QuantOpt[c(1, 2)[isEukaLike+1]]
if ("Prot.Quant.Use" %in% colnames(Param)) {
  Pep4Quant <- Pep4QuantOpt[which(c("UNIQUE", "RAZOR", "ALL") == toupper(Param$Prot.Quant.Use))]
}
ProtRulNuclL %<o% Param$ProtRulNuclL
CytoScVrs %<o% Param$CytoScapePath
CytoScape %<o% file.exists(CytoScVrs)
Nested %<o% as.logical(Param$Ratios.Groups_Nested)
if (!Nested) {
  Param$Ratios.Ref.Groups <- Param_filter(Param$Ratios.Ref.Groups, "Rep")
} else {
  Param$Ratios.Ref.Groups <- paste(unique(c(unlist(strsplit(Param$Ratios.Ref.Groups, ";")), "Rep")),
                                   collapse = ";")
}
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

if (SearchSoft %in% c("DIANN", "FRAGPIPE")) {
  tmp <- unlist(strsplit(Param$PSMs, ";"))
  tmp <- tmp[which(file.exists(tmp))]
  tmp <- tmp[which(!tmp %in% PSMsFl)]
  if (length(tmp)) {
    msg <- paste0("Parameter \"PSMs\" is ignored in ", SearchSoft, " mode , instead the PSMs file is detected from the search parameters\n",
                  "The following file(s) will not be included:", paste0("\n - ", tmp), "\n")
    ReportCalls <- AddMsg2Report(Warning = TRUE, Space = TRUE)
  }
}
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\" -> Parameters:\", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$left))")
for (i in 1:ncol(Param)) {
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              paste0("body_add_fpar(Report, fpar(ftext(\"   - ", colnames(Param)[i], ": ",
                                     Param[[i]], "\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$left))"))
}
ReportCalls <- AddSpace2Report()

# Create sub-directories vector:
dir <- c("Workflow control/MA plots", paste0("Workflow control/", evNm, "s", c("", "/Normalisation")),
         "Workflow control/Peptides", "Workflow control/Peptides/PCA plot", "Workflow control/Peptides/Intensities",
         "Workflow control/Peptides/Intensities/Normalisation - summary", "Workflow control/Peptides/Ratios",
         "Workflow control/Protein groups", "Workflow control/Protein groups/Expression",
         "Workflow control/Protein groups/Ratios", "Workflow control/Protein groups/P-values",
         "Reg. analysis", "Reg. analysis/t-tests",
         #"Reg. analysis/t-tests/pdf",
         #"Reg. analysis/t-tests/jpeg",
         #"Reg. analysis/t-tests/html",
         "PCA plots", "t-SNE plots", "Heatmaps", "Tables")
if (("Norma.Pep.Intens.Shape" %in% colnames(Param))&&(toupper(Param$Norma.Pep.Intens.Shape) %in% c("VSN", "LOESS"))) {
  dir <- c(dir, paste0("Workflow control/Peptides/Intensities/", toupper(Param$Norma.Pep.Intens.Shape), " normalisation"))
}
if (("Norma.Pep.Intens.IRS" %in% colnames(Param))&&(Param$Norma.Pep.Intens.IRS == TRUE)) {
  dir <- c(dir, "Workflow control/Peptides/Intensities/IRS normalisation")
}
if (("Batch.correction" %in% colnames(Param))&&(!as.character(Param$Batch.correction) %in% c("", "F", "FALSE"))) {
  dir <- c(dir, "Workflow control/Peptides/Intensities/Batch correction")
}
enrichGO %<o% (("GO.enrichment" %in% colnames(Param))&&(Param$GO.enrichment))
globalGO %<o% (("GO.enrichment_Whole_dataset" %in% colnames(Param))&&(Param$GO.enrichment_Whole_dataset))
if (enrichGO||globalGO) { dir <- c(dir, "Reg. analysis/GO enrich") }
dirlist <- unique(c(dirlist, paste0(wd, "/", dir)))

# In case we have enriched for a PTM, it helps to check how good the enrichment was:
# This should be before any PSMs are filtered out - we want to look at the data "straight out of the MS"
tstEnrich <- unique(FracMap$`PTM-enriched`)
tstEnrich <- tstEnrich[which((!is.na(tstEnrich))&(tstEnrich != "NA"))]
if (length(tstEnrich)) {
  dir <- paste0(wd, "/Summary plots")
  dirlist <- unique(c(dirlist, dir))
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  for (Mod in tstEnrich) { #Mod <- tstEnrich[1]
    pat <- paste0("\\(", Modifs$Mark[match(Mod, Modifs$`Full name`)], "\\)")
    ev[[Mod]] <- grepl(pat, ev$"Modified sequence")
    for (Type in c("PSMs", "Pep")) { #Type <- "PSMs"
      if (Type == "PSMs") {
        tst <- data.table(Mod = ev[[Mod]],
                          MS_file = ev$"Raw file path",
                          temp = 1)
        tst <- tst[, list(Count = sum(temp)), by = list(MS_file = MS_file, Mod = Mod)]
        tst <- as.data.frame(tst)
        Root <- "PSM"
      }
      if (Type == "Pep") {
        tst <- data.table(Mod = ev[[Mod]],
                          MS_file = ev$"Raw file path",
                          ModSeq = ev$`Modified sequence`)
        tst <- tst[, list(Count = length(unique(ModSeq))), by = list(MS_file = MS_file, Mod = Mod)]
        tst <- as.data.frame(tst)
        Root <- "peptidoform"
      }
      colnames(tst)[2] <- Mod
      m <- match(tst$MS_file, FracMap$"Raw file")
      tst$Sample <- FracMap$MQ.Exp[m]
      tst <- tst[which(!is.na(tst$Sample)),]
      tst <- tst[which(tst$Sample %in% unlist(Exp.map$MQ.Exp)),]
      tst <- tst[which(tst$Sample %in% unlist(Exp.map$MQ.Exp)),]
      #which(sapply(Exp.map$MQ.Exp, function(y) { x %in% unlist(y) }))
      tst2 <- reshape2::melt(tst)
      colnames(tst2) <- gsub("\\(|\\)|\\[|\\]", "", gsub(" ", "_", colnames(tst2)))
      frml <- as.formula(paste0("MS_file ~ `", gsub("\\(|\\)|\\[|\\]", "", gsub(" ", "_", Mod)), "`"))
      tst2 <- cast(tst2, frml, fun.aggregate = sum)
      kN <- paste0("Non-", Mod, "-modified")
      kY <- paste0(Mod, "-modified")
      colnames(tst2)[which(colnames(tst2) == "FALSE")] <- kN
      colnames(tst2)[which(colnames(tst2) == "TRUE")] <- kY
      tst2[[paste0(Mod, " [%]")]] <- signif(100*tst2[[kY]]/(tst2[[kY]]+tst2[[kN]]), 3)
      tst2$Sample <- tst$Sample[match(tst2$MS_file, tst$MS_file)]
      tst2 <- tst2[, c("Sample", "MS_file", kN, kY, paste0(Mod, " [%]"))]
      write.csv(tst2, paste0(dir, "/", Mod, "-", Root, "s per MS file.csv"), row.names = FALSE)
      rw <- unique(tst$MS_file)
      for (i in c(TRUE, FALSE)) {
        w <- which(tst[[Mod]] == i)
        w2 <- which(!rw %in% tst$MS_file[w])
        if (length(w2)) {
          tmp <- tst[which((tst$MS_file %in% rw[w2])&(tst[[Mod]] == !i)),]
          tmp[[Mod]] <- i
          tmp$Count <- 0
          tst <- rbind(tst, tmp)
        }
      }
      tst[[Mod]] <- factor(c("-", "+")[tst[[Mod]]+1], levels = c("-", "+"))
      ttl <- paste0(Mod, "-", Root, "s per MS file")
      plot <- ggplot(tst) + geom_bar(stat = "identity", position = "dodge",
                                     aes(x = MS_file, y = Count, fill = .data[[Mod]])) +
        theme_bw() + scale_fill_viridis(discrete = TRUE) + ggtitle(ttl) +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 5),
              plot.margin = unit(c(0, 0, 0, 3), "in"))
      #poplot(plot, 12, 20)
      ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300)
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)
      ReportCalls <- AddPlot2Report()
    }
  }
}

# Write PTMs table
temp <- Modifs
w <- which(sapply(colnames(Modifs), function(x) { class(Modifs[[x]]) }) == "list")
for (i in w) { temp[[i]] <- sapply(temp[[i]],  paste, collapse = ", ") }
dir <- paste0(wd, "/Workflow control")
dirlist <- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
write.csv(temp, paste0(dir, "/Modifications.csv"), row.names = FALSE)
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\"PTMs table:\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_table(Report, ReportCalls$Objects$AABiases)")
ReportCalls$Objects$AABiases <- temp
ReportCalls <- AddSpace2Report()

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
      names(J) <- Aggregates[sapply(J, function(x) {which(names(Aggregates) == x)})]
      Exp.map[paste(J, collapse = "")] <- apply(Exp.map[, names(J)], 1, function(x) {
        paste(x, collapse = "___")
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
Aggregate.list %<o% sapply(Aggregate.map$Aggregate.Name, function(x) { get(x) })
#This doesn't work for the master/detailed dual script approach (issue with environments!):

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
if (("Batch.correction" %in% colnames(Param))&&(!as.character(Param$Batch.correction) %in% c("", "F", "FALSE"))) {
  parse.Param.aggreg.2("Batch.correction")
}
a <- RSA$names
if (length(a) == 1) {
  Exp.map$Ref.Sample.Aggregate <- Exp.map[[a]]
} else { Exp.map$Ref.Sample.Aggregate <- apply(Exp.map[, a], 1, paste, collapse = "___") }
Exp.map <- Exp.map[order(Exp.map[[VPAL$column]], Exp.map$Replicate),]
tst <- setNames(lapply(MQ.Exp, function(x) { which(sapply(Exp.map$MQ.Exp, function(y) { sum(y %in% x) }) > 0) }), MQ.Exp)
tst2 <- sapply(tst, length)
if (LabelType == "LFQ") {
  stopifnot(max(tst2) == 1)
  tst <- tst[which(tst2 == 1)]
}
MQ.Exp <- MQ.Exp[which(MQ.Exp %in% names(tst))]
ev <- ev[which(ev$MQ.Exp %in% names(tst)),]
if (LabelType == "LFQ") {
  ev$RSA <- Exp.map$Ref.Sample.Aggregate[match(ev$MQ.Exp, names(tst))]
  tmp <- ev[, c("Proteins", "Intensity", "RSA")]
} else {
  # Isobaric case - subtle difference
  tmp <- ev[, c("Proteins", "Intensity")]
  tmp$RSA <- Exp.map$Ref.Sample.Aggregate[match(ev$MQ.Exp, names(tst))]
}
# Plot of contamination levels per sample
dir <- paste0(wd, "/Summary plots")
dirlist <- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
tmp2 <- listMelt(strsplit(tmp$Proteins, ";"), 1:nrow(tmp), c("Protein", "Row"))
m <- match(tmp2$Protein, db$`Protein ID`)
w <- which(!is.na(m))
tmp2 <- tmp2[w,]; m <- m[w]
# For our purpose here we must match contaminant proteins.
tmp2$Cont <- db$`Potential contaminant`[m]
if (tstorg) {
  tmp2$Organism <- db[m, dbOrgKol]
  tmp2 <- as.data.table(tmp2)
  f0 <- function(x) { c("Target", "Contaminant")[("Contaminant" %in% x)+1] }
  tmp2 <- tmp2[, list(x = f0(Organism)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
} else {
  tmp2 <- as.data.table(tmp2)
  f0 <- function(x) { c("Target", "Contaminant")[("+" %in% x)+1] }
  tmp2 <- tmp2[, list(x = f0(Cont)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
}
tmp$Organism <- tmp2$x[match(1:nrow(tmp), tmp2$Group.1)]
tmp <- tmp[which(!is.na(tmp$Organism)),]
tmp$Organism <- factor(tmp$Organism, levels = c("Contaminant", "Target"))
tmp$Intensity <- as.numeric(tmp$Intensity)
tmp <- aggregate(tmp$Intensity, list(tmp$RSA, tmp$Organism), sum, na.rm = TRUE)
colnames(tmp) <- c("Sample", "Organism", "Total intensity")
tmp$Sample <- factor(cleanNms(tmp$Sample), levels = cleanNms(Exp.map$Ref.Sample.Aggregate))
ttl <- "Contributions to TIC"
plot <- ggplot(tmp) +
  geom_bar(stat = "identity", aes(x = Sample, y = `Total intensity`, fill = Organism)) +
  theme_bw() + scale_fill_viridis(discrete = TRUE, begin = 0.8, end = 0.2) +
  ggtitle(ttl, subtitle = "Summed TIC for each class of identified peptides") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
ggsave(paste0(wd, "/Summary plots/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
ggsave(paste0(wd, "/Summary plots/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")

# Time points
if (exists("Tim")) {
  if (("Time.Points" %in% colnames(Param))&&(Param$Time.Points != "")) {
    tp %<o% as.character(sort(as.numeric(unlist(strsplit(Param$Time.Points, ";")))))
    if (("Time.Point.Names" %in% colnames(Param))&&(Param$Time.Point.Names != "")) {
      names(tp) <- gsub(" ", ".", unlist(strsplit(Param$Time.Point.Names, ";")))
    } else {
      names(tp) <- tp
    }
    if (sum(tp != Tim) > 0) {
      stop("Review your time points!")
    } else {
      Tim <- tp
    }
  } else {
    Tim <- as.character(sort(as.numeric(Tim)))
  }
}

# Refresh list of interesting proteins:
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList.R")
#rstudioapi::documentOpen(Src)
source(Src)
Src <- paste0(libPath, "/extdata/R scripts/Sources/protList2.R")
#rstudioapi::documentOpen(Src)
source(Src)

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
    if (is.na(tmp)||(tmp == "")) { tmp <- c() } else {
      tmp <- unlist(strsplit(tmp, ";"))
    }
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
  #test <- sapply(test, function(x) {sum(x %in% temp$"Protein ID")}) > 0
  #cont.ev %<o% ev[which(test),]
  #ev <- ev[which(!test),]
  if (exists("contDB")) { contDB <- plyr::rbind.fill(list(contDB, temp)) } else { contDB <- temp }
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
  tst3 <- (!sum(!"logical" %in% sapply(RG$values, function(x) { class(DiscFiltTbl[[x]]) })))
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
if (("Two.sided" %in% colnames(Param))&&(is.logical(Param$Two.sided))) {
  TwoSided %<o% Param$Two.sided
} else { TwoSided %<o% !IsPullDown }
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

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc)
