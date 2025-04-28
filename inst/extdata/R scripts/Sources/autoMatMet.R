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
  MatMetCalls$Calls <- append(MatMetCalls$Calls, paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$WetLab[", i,"], prop = WrdFrmt$",
                                                        c("Body", "Template_text")[(WetLabMeth[i] == "TEMPLATE")+1], "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#
# 2) LCMS
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_fpar(MatMet, fpar(ftext(\"LC-MS/MS analysis\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
LCMSMeth <- try(MatMet_LCMS(), silent = TRUE)
if ((("try-error" %in% class(LCMSMeth))||(is.null(LCMSMeth)))&&("mzML" %in% gsub(".*\\.", "", rawFiles))) {
  tmp <- gsub("\\.mzML", ".raw", rawFiles)
  LCMSMeth <- try(MatMet_LCMS(RawFiles = tmp), silent = TRUE)
  if ((("try-error" %in% class(LCMSMeth))||(is.null(LCMSMeth)))&&("mzML" %in% gsub(".*\\.", "", rawFiles))) {
    tmp <- gsub("\\.mzML", ".d", rawFiles)
    LCMSMeth <- try(MatMet_LCMS(RawFiles = tmp), silent = TRUE)
  }
}
if ((("try-error" %in% class(LCMSMeth)))||(is.null(LCMSMeth))) { LCMSMeth <- "TEMPLATE" }
L <- length(LCMSMeth)
if (L > 1) {
  LCMSMeth <- lapply(1:length(LCMSMeth), function(x) { c(paste0("Method ", x, ":"), LCMSMeth[[x]]) })
}
LCMSMeth <- unlist(LCMSMeth)
MatMetCalls$Texts$LCMS <- LCMSMeth
for (i in 1:length(LCMSMeth)) {
  MatMetCalls$Calls <- append(MatMetCalls$Calls, paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts$LCMS[", i,"], prop = WrdFrmt$",
                                                        c("Body", "Template_text")[(LCMSMeth[i] == "TEMPLATE")+1], "_text), fp_p = WrdFrmt$just))"))
}
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
#
# 3) Search
MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_fpar(MatMet, fpar(ftext(\"Data analysis\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
kol <- c("Organism_Full", "Organism")
kol <- kol[which(kol %in% colnames(db))]
tst <- sapply(kol, function(x) { length(unique(db[which(!as.character(db[[x]]) %in% c("", "NA")), x])) })
kol <- kol[order(tst, decreasing = TRUE)][1]
w <- which(db$`Potential contaminant` != "+")
Org %<o% aggregate(w, list(db[w, kol]), length)
colnames(Org) <- c("Organism", "Count")
Org <- Org[which(Org$Count == max(Org$Count)[1]),]
Org$Source <- aggregate(db$Source[which(db[[kol]] %in% Org$Organism)], list(db[which(db[[kol]] %in% Org$Organism), kol]), function(x) {
  unique(x[which(!is.na(x))])
})$x
Org$Source[which(is.na(Org$Source))] <- ""
tstOrg <- c("", "n")[(tolower(substr(Org$Organism, 1, 1)) %in% c("a", "e", "i", "o", "u"))+1]
if (nrow(Org) == 1) {
  dbTxt <- paste0("a", tstOrg, " ", Org$Organism, " proteome.")
  if (Org$Source != "") { dbTxt <- gsub("\\.$", paste0(" sourced from ", Org$Source, "."), dbTxt) }
} else {
  stop()
  l <- nrow(Org)
  tstSrc <- unique(Org$Source)
  if (length(tstSrc) == 1) {
    dbTxt <- paste0(paste(Org$Organism[1:(l-1)], collapse = ", "), " and ", Org$Organism[l], " proteomes",
                    c(".", paste0(" sourced from ", tstSrc, "."))[(tstSrc != "")+1])
  } else {
    dbTxt <- sapply(1:l, function(x) { paste0(Org$Organism[x], " (", Org$Source[x], ")") })
    dbTxt <- paste0(paste(dbTxt[1:(l-1)], collapse = ", "), " and ", dbTxt[l], " proteomes.")
  }
}
moult <- (length(rawFiles2) > 1)+1
DatAnalysisTxt %<o% "TEMPLATE"
if (SearchSoft == "MAXQUANT") {
  SearchSoftVers <- gsub(" *</?maxQuantVersion> *", "", grep("<maxQuantVersion>", mqpar, value = TRUE))
  gFx <- (grep("<fixedModifications>", mqpar)+1):(grep("</fixedModifications>", mqpar)-1)
  gVar <- (grep("<variableModifications>", mqpar)+1):(grep("</variableModifications>", mqpar)-1)
  FxMd <- gsub(" *</?string> *", "", mqpar[gFx])
  FxMdC <- grep("\\(C\\)", FxMd, value = TRUE)
  FxMd <- FxMd[which(!FxMd %in% FxMdC)]
  VarMd <- gsub(" *</?string> *", "", mqpar[gVar])
  SearchTxt <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in ",
                      names(SearchSoft))
  if ((length(SearchSoftVers))&&(nchar(SearchSoftVers))) {
    SearchTxt <- paste0(SearchTxt, " version ", SearchSoftVers)
  }
  SearchTxt <- paste0(SearchTxt, " against ", dbTxt)
  MBR <- gsub(" *</?matchBetweenRuns> *", "", grep("<matchBetweenRuns>", mqpar, value = TRUE))
  SecPep <- gsub(" *</?secondPeptide> *", "", grep("<secondPeptide>", mqpar, value = TRUE))
  DPep <- gsub(" *</?dependentPeptides> *", "", grep("<dependentPeptides>", mqpar, value = TRUE))
  tstFDRs <- grep("Fdr", mqpar, value = TRUE)
  FDRs <- gsub(" *</?[A-Z,a-z]*Fdr[A-Z,a-z]*> *", "", tstFDRs)
  w <- which(!FDRs %in% c("False", "True"))
  FDRs <- setNames(as.numeric(FDRs[w]), gsub(" *</?|>.+", "", tstFDRs[w]))
  FDRs <- FDRs[which(names(FDRs) != "psmFdrCrosslink")]
  if (DPep == "False") { FDRs <- FDRs[which(names(FDRs) != "dependentPeptideFdr")] }
  if (length(FxMdC)) {
    if (length(FxMdC) == 1) {
      txt <- paste0("Fixed cysteine modification was set to ", FxMdC, ".")
    } else {
      txt <- paste0(paste(FxMdC[1:(length(FxMdC)-1)], collapse = ", "), " and ", rev(FxMdC)[1], " were included as fixed cysteine modifications.")
    }
    SearchTxt <- paste0(SearchTxt, " ", txt)
  }
  if (length(FxMd)) {
    if (length(FxMd) == 1) {
      txt <- paste0("An additional fixed modifications, ", FxMd, ", was also included.")
    } else {
      txt <- paste0("In addition, ", paste(FxMd[1:(length(FxMd)-1)], collapse = ", "), " and ", rev(FxMd)[1], " were included as fixed modifications.")
    }
    SearchTxt <- paste0(SearchTxt, " ", txt)
  }
  if (length(VarMd)) {
    if (length(VarMd) == 1) {
      txt <- paste0(VarMd, " was set as variable modification.")
    } else {
      txt <- paste0("Variable modifications were set to ", paste(VarMd[1:(length(VarMd)-1)], collapse = ", "), " and ", rev(VarMd)[1], ".")
    }
    SearchTxt <- paste0(SearchTxt, " ", txt)
  }
  tsts <- setNames(as.logical(toupper(c(MBR, SecPep, DPep))), c("Match-between-Runs", "Second Peptide Search", "Dependent Peptides"))
  y <- which(tsts); ly <- length(y)
  n <- which(!tsts); ln <- length(n)
  if (ly) {
    if (ly == 1) { txt <- paste0(names(tsts)[y], " was set to active") } else {
      txt <- paste0(paste(names(tsts)[y[1:(ly-1)]], collapse = ", "), " and ", names(tsts)[y[ly]], " were set to active")
    }
    if (ln) {
      if (ln == 1) { txt <- paste0(txt, ", while ", names(tsts)[n], " was turned off.") } else {
        txt <- paste0(txt, ", while ", paste(names(tsts)[n[1:(ln-1)]], collapse = ", "), " and ", names(tsts)[n[ln]], " were turned off.")
      }
    } else { txt <- paste0(txt, ".") }
  } else {
    txt <- paste0(paste(names(tsts)[n[1:(ln-1)]], collapse = ", "), " and ", names(tsts)[n[ln]], " were all turned off.")
  }
  SearchTxt <- paste0(SearchTxt, " ", txt)
  tst <- unique(FDRs)
  lf <- length(FDRs)
  nms <- gsub("Fdr$", "s", names(FDRs))
  if (length(tst) == 1) { txt <- paste0("All FDRs were set to ", tst*100, "%.") } else {
    txt <- paste0("FDRs were set to: ", paste(paste0(FDRs[1:(lf-1)]*100, "% (", nms[1:(lf-1)], ")"), collapse = ", "), " and ", paste0(FDRs[lf]*100, "% (", nms[lf], ")"), ".")
  }
  SearchTxt <- paste0(SearchTxt, " ", txt)
}
if (SearchSoft == "DIANN") {
  SearchSoftVers <- gsub("^DIA-NN | .+$", "", grep("^DIA-NN [0-9]+\\.?[0-9]*", DIANNlog, value = TRUE))
  DIANNCall2 <- unlist(strsplit(DIANNCall, " +--"))
  FxMd <- gsub(" enabled as a fixed modification$", "", grep(" enabled as a fixed modification$", DIANNlog, value = TRUE))
  FxMdC <- grep("Cysteine", FxMd, value = TRUE)
  FxMd <- FxMd[which(!FxMd %in% FxMdC)]
  FxMdC <- gsub("cysteine[ ,-]?", "", FxMdC, ignore.case = TRUE)
  FxMdC <- paste0(toupper(substr(FxMdC, 1, 1)), substr(FxMdC, 2, nchar(FxMdC)))
  require(unimod)
  UniMod <- unimod::modifications
  VarMd <- gsub(" will be considered as variable$", "", grep(" will be considered as variable$", DIANNlog, value = TRUE))
  if (length(VarMd)) {
    VarMd <- set_colnames(Isapply(strsplit(gsub("^Modification ", "", VarMd), " with mass delta | at "), unlist), c("UniMod", "Delta mass", "AA"))
    VarMd <- set_colnames(aggregate(VarMd$AA, list(VarMd$UniMod, VarMd$`Delta mass`), list), c("UniMod", "Delta mass", "AA"))
    VarMd$UniMod <-  gsub("^UniMod:", "", VarMd$UniMod)
    VarMd$Type <- "Variable"
    VarMd$"Full name" <- apply(VarMd[, c("UniMod", "Delta mass")], 1, function(x) { #x <- VarMd[1, c("UniMod", "Delta mass")]
      unique(UniMod$Name[which((UniMod$UnimodId == x[[1]])&(round((UniMod$MonoMass == x[[2]])*2)/2 == 0))])
    })
    w <- which(lapply(VarMd$`Full name`, length) == 0)
    if (length(w)) {
      VarMd$`Full name`[w] <- apply(VarMd[w, c("Delta mass", "AA")], 1, function(x) {
        pos <- sort(unlist(x[[2]]))
        pos[which(pos == "*n")] <- "_" # Check, this may not be correct...
        # in fact, check whether this part of DIANN_to_MQ() should not be improved.
        Modifs$`Full name`[which((Modifs$`Mass shift` == x[[1]])&(Modifs$AA %in% pos))]
      })
    }
    VarMd$Text <- apply(VarMd[, c("Full name", "AA")], 1, function(x) {
      x[[2]][which(x[[2]] == "n")] <- "N-term"
      x[[2]][which(x[[2]] == "*n")] <- "protein N-term"
      x[[2]][which(x[[2]] == "c")] <- "C-term"
      x[[2]][which(x[[2]] == "*c")] <- "protein C-term"
      paste0(x[[1]], " (", paste(x[[2]], collapse = ""), ")")
    })
  }
  # Note that it is possible that a modification was searched but not found (so is present in VarMd but not Modifs)
  MBR <- grepl("--reanaly[s,z]e", DIANNCall)+1
  Libraries <- gsub("^lib +", "", grep("^lib", DIANNCall2, value = TRUE))
  nL <- length(Libraries)
  LibFree <- (nL == 0)+1
  SearchTxt <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in ",
                      names(SearchSoft))
  if ((length(SearchSoftVers))&&(nchar(SearchSoftVers))) {
    SearchTxt <- paste0(SearchTxt, " version ", SearchSoftVers)
  }
  SearchTxt <- paste0(SearchTxt,
                      c(paste0(" against the following spectral librar", c("y", "ies")[(nL > 1)+1], ": ",
                               paste(Libraries, collapse = "/"), " (see supplementary material), using "),
                        " in library-free mode against ")[LibFree], dbTxt)
  SearchTxt <- paste0(SearchTxt, " Match-Between-Runs was turned ", c("on", "off")[MBR], ".")
  if (length(FxMdC)) {
    if (length(FxMdC) == 1) {
      txt <- paste0("Fixed cysteine modification was set to ", FxMdC, ".")
    } else {
      txt <- paste0(paste(FxMdC[1:(length(FxMdC)-1)], collapse = ", "), " and ", rev(FxMdC)[1], " were included as fixed cysteine modifications.")
    }
    SearchTxt <- paste0(SearchTxt, " ", txt)
  }
  if (length(FxMd)) {
    if (length(FxMd) == 1) {
      txt <- paste0("An additional fixed modifications, ", FxMd, ", was also included.")
    } else {
      txt <- paste0("In addition, ", paste(FxMd[1:(length(FxMd)-1)], collapse = ", "), " and ", rev(FxMd)[1], " were included as fixed modifications.")
    }
    SearchTxt <- paste0(SearchTxt, " ", txt)
  }
  if (("data.frame" %in% class(VarMd))&&(nrow(VarMd))) {
    if (nrow(VarMd) == 1) {
      txt <- paste0(VarMd$Text, " was set as variable modification.")
    } else {
      txt <- paste0("Variable modifications were set to ", paste(VarMd$Text[1:(nrow(VarMd)-1)], collapse = ", "), " and ", rev(VarMd$Text)[1], ".")
    }
    SearchTxt <- paste0(SearchTxt, " ", txt)
  }
  tstFDR <- grep("^Output will be filtered at .+ FDR$", DIANNlog, value = TRUE)
  FDR <- gsub("^Output will be filtered at | FDR$", "", tstFDR)
  if (!grepl(" ?%$", FDR)) { FDR <- as.numeric(FDR)*100 }
  if ((FDR >= 0)&&(FDR <= 1)) { SearchTxt <- paste0(SearchTxt, " Data was filtered at ", FDR, "% FDR.") }
}
if (SearchSoft == "FRAGPIPE") {
  SearchSoftVers <- gsub("^# FragPipe \\(|\\).+$", "", grep("^# FragPipe \\(", FP_Workflow, value = TRUE))
  FxMdTbl <- gsub(topattern("msfragger.table.fix-mods="), "", grep(topattern("msfragger.table.fix-mods="), FP_Workflow, value = TRUE))
  FxMdTbl <- unlist(strsplit(FxMdTbl, ";"))
  FxMdTbl <- Isapply(strsplit(FxMdTbl, ","), unlist)
  colnames(FxMdTbl) <- c("Shift", "Site", "Enabled", "MaxOccurences")
  FxMdTbl <- FxMdTbl[which(FxMdTbl$Enabled == "true"),]
  FxMdTbl$AAName <- gsub(".+\\(|\\)", "", FxMdTbl$Site)
  FxMdTbl$AAName <- sapply(FxMdTbl$AAName, function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  })
  FxMdTbl$AA <- gsub(" \\(.+\\)", "", FxMdTbl$Site)
  FxMdTbl$Shift <- as.numeric(FxMdTbl$Shift)
  FxMdTbl <- FxMdTbl[which(FxMdTbl$Shift != 0),]
  FxMdTbl$Shift <- as.character(FxMdTbl$Shift)
  g <- grep("^-", FxMdTbl$Shift, invert = TRUE)
  FxMdTbl$Shift[g] <- paste0("+", FxMdTbl$Shift[g])
  FxMdTbl$Name <- apply(FxMdTbl[, c("Shift", "AAName")], 1, function(x) { paste0(x[[1]], " (", x[[2]], ")") })
  FxMdC <- FxMdTbl$Name[which(FxMdTbl$AA == "C")]
  FxMd <- FxMdTbl$Name[which(FxMdTbl$AA != "C")]
  VarMd <- gsub(topattern("msfragger.table.var-mods="), "", grep(topattern("msfragger.table.var-mods="), FP_Workflow, value = TRUE))
  VarMd <- unlist(strsplit(VarMd, ";"))
  VarMd <- Isapply(strsplit(VarMd, ","), unlist)
  colnames(VarMd) <- c("Shift", "Site", "Enabled", "MaxOccurences")
  VarMd <- VarMd[which(VarMd$Enabled == "true"),]
  for (i in 1:nrow(AA_table)) { VarMd$Site <- gsub(AA_table$AA[i], paste0(";", AA_table$AA[i], ";"), VarMd$Site) }
  VarMd$Site <- gsub("n", ";peptide N-term;", VarMd$Site)
  VarMd$Site <- gsub("c", ";peptide C-term;", VarMd$Site)
  VarMd$Site <- gsub("\\[\\^", ";protein N-term;", VarMd$Site)
  VarMd$Site <- gsub("\\^\\]", ";protein C-term;", VarMd$Site)
  VarMd$Site <- strsplit(gsub("^;|;$", "", VarMd$Site), ";;")
  VarMd$AAName <- sapply(VarMd$Site, function(x) { #
    x <- unlist(x)
    w <- which(nchar(x) == 1)
    x[w] <- as.character(AA_table$Name[match(x[w], as.character(AA_table$AA))])
    return(x)
  })
  VarMd$Shift <- gsub(" +", "", VarMd$Shift)
  g <- grep("^-", VarMd$Shift, invert = TRUE)
  VarMd$Shift[g] <- paste0("+", VarMd$Shift[g])
  VarMd$Text <- apply(VarMd[, c("Shift", "AAName", "Site")], 1, function(x) {
    x2 <- x[[2]]
    if (length(x2) > 1) { x2 <- paste(x[[3]], collapse = "") }
    return(paste0(x[[1]], " (", x2, ")"))
  })
  # Note that it is possible that a modification was searched but not found (so is present in VarMd but not Modifs)
  #
  # FragPipe is very versatile... which will not make this easy...
  # Unfortunately, the workflow file does not appear to store the type of Workflow which was preselected in the 2nd tab.
  Ump <- as.logical(toupper(gsub(topattern("diaumpire.run-diaumpire="), "", grep(topattern("diaumpire.run-diaumpire="), FP_Workflow, value = TRUE))))
  Diane <- as.logical(toupper(gsub(topattern("diann.run-dia-nn="), "", grep(topattern("diann.run-dia-nn="), FP_Workflow, value = TRUE))))
  isDIA <- Diane|Ump # I think theoretically you could use either without the other... although probably for any DIA dataset you will use DiaNN
  #stopifnot(((Param$Label == "DIA")&isDIA)|((Param$Label != "DIA")&!isDIA)) # Sanity check
  # OpenSearch is indirectly identified as cases where unusually high values for tolerances are used
  PrTolUp <- gsub(topattern("msfragger.precursor_mass_upper="), "", grep(topattern("msfragger.precursor_mass_upper="), FP_Workflow, value = TRUE))
  PrTolDwn <- gsub(topattern("msfragger.precursor_mass_lower="), "", grep(topattern("msfragger.precursor_mass_lower="), FP_Workflow, value = TRUE))
  # Validations
  Chris <- as.logical(toupper(gsub(topattern("crystalc.run-crystalc="), "", grep(topattern("crystalc.run-crystalc="), FP_Workflow, value = TRUE))))
  Moses <- as.logical(toupper(gsub(topattern("peptide-prophet.run-peptide-prophet="), "", grep(topattern("peptide-prophet.run-peptide-prophet="), FP_Workflow, value = TRUE))))
  Illy <- as.logical(toupper(gsub(topattern("percolator.run-percolator="), "", grep(topattern("percolator.run-percolator="), FP_Workflow, value = TRUE))))
  Maud <- as.logical(toupper(gsub(topattern("ptmprophet.run-ptmprophet="), "", grep(topattern("ptmprophet.run-ptmprophet="), FP_Workflow, value = TRUE))))
  Camus <- as.logical(toupper(gsub(topattern("phi-report.run-report="), "", grep(topattern("phi-report.run-report="), FP_Workflow, value = TRUE))))
  David <- as.logical(toupper(gsub(topattern("ptmshepherd.run-shepherd="), "", grep(topattern("ptmshepherd.run-shepherd="), FP_Workflow, value = TRUE))))
  Kant <- as.logical(toupper(gsub(topattern("quantitation.run-label-free-quant="), "", grep(topattern("quantitation.run-label-free-quant="), FP_Workflow, value = TRUE))))
  Jon <- as.logical(toupper(gsub(topattern("ionquant.run-ionquant="), "", grep(topattern("ionquant.run-ionquant="), FP_Workflow, value = TRUE))))
  Matt <- as.logical(as.numeric(gsub(topattern("ionquant.mbr="), "", grep(topattern("ionquant.mbr="), FP_Workflow, value = TRUE))))
  Fritz <- as.logical(toupper(gsub(topattern("freequant.run-freequant="), "", grep(topattern("freequant.run-freequant="), FP_Workflow, value = TRUE))))
  Timothy <- as.logical(toupper(gsub(topattern("tmtintegrator.run-tmtintegrator="), "", grep(topattern("tmtintegrator.run-tmtintegrator="), FP_Workflow, value = TRUE))))
  #
  SearchTxt <- paste0(c("The r", "R")[moult], "aw file", c(" was", "s were")[moult], " searched in ",
                      names(SearchSoft))
  if ((length(SearchSoftVers))&&(nchar(SearchSoftVers))) {
    SearchTxt <- paste0(SearchTxt, " version ", SearchSoftVers)
  }
  SearchTxt <- paste0(SearchTxt, " against ", dbTxt)
  if (Ump == 1) {
    SearchTxt <- gsub("\\.$",
                      " As a preliminary step, the DIAUmpire module was run to extract pseudo-MS/MS spectra from DIA spectra files.",
                      SearchTxt)
  }
  if (length(FxMdC)) {
    if (length(FxMdC) == 1) {
      txt <- paste0("Fixed cysteine modification was set to ", FxMdC, ".")
    } else {
      txt <- paste0(paste(FxMdC[1:(length(FxMdC)-1)], collapse = ", "), " and ", rev(FxMdC)[1], " were included as fixed cysteine modifications.")
    }
    SearchTxt <- paste0(SearchTxt, " ", txt)
  }
  if (length(FxMd)) {
    if (length(FxMd) == 1) {
      txt <- paste0("An additional fixed modifications, ", FxMd, ", was also included.")
    } else {
      txt <- paste0("In addition, ", paste(FxMd[1:(length(FxMd)-1)], collapse = ", "), " and ", rev(FxMd)[1], " were included as fixed modifications.")
    }
    SearchTxt <- paste0(SearchTxt, " ", txt)
  }
  if (nrow(VarMd)) {
    if (nrow(VarMd) == 1) {
      txt <- paste0(VarMd$Text, " was set as variable modification.")
    } else {
      txt <- paste0("Variable modifications were set to ", paste(VarMd$Text[1:(nrow(VarMd)-1)], collapse = ", "), " and ", rev(VarMd$Text)[1], ".")
    }
    SearchTxt <- paste0(SearchTxt, " ", txt)
  }
  w <- which(c(PrTolUp != "20", PrTolDwn != "-20"))
  l <- length(w)
  if (l) {
    PrTol <- c("lower", "upper")[w]
    PrTol[1] <- paste0(toupper(substr(PrTol[1], 1, 1)), substr(PrTol[1], 2 , nchar(PrTol[w])))
    PrTol <- paste(PrTol, collapse = " and ")
    SearchTxt <- paste0(SearchTxt, " ", PrTol, " precursor mass tolerance", c("", "s")[l], " w", c("as", "ere")[l], " set to ",
                        paste(c(PrTolDwn, PrTolUp)[w], collapse = " and "), " ppm", c("", ", respectively")[l], ".")
  }
  # Validations
  w <- which(c(Chris, Moses, Illy)) # Normally should only ever be length = 1
  l <- length(w)
  if (l) {
    txt <- c("Crystal-C", "PeptideProphet", "Percolator")[w]
    if (l > 1) { txt <- paste0(paste(txt[1:(l-1)], collapse = ", "), " and ", txt[l]) }
    SearchTxt <- paste0(SearchTxt, " Peptide identifications were validated using ", txt, ".")
    if (Maud) { SearchTxt <- gsub("\\.$", ", and", SearchTxt) }
  }
  if (Maud) { SearchTxt <- paste0(SearchTxt, " PTM sites", c(" were", "")[(l>0)+1], " localized using PTMProphet.") }
  # FDR levels
  if (Camus) {
    CamusPar <- toupper(gsub(topattern("phi-report.filter="), "", grep(topattern("phi-report.filter="), FP_Workflow, value = TRUE)))
    CamusPar <- unlist(strsplit(CamusPar, "--"))
    CamusPar <- CamusPar[2:length(CamusPar)]
    CamusPar <- data.frame(Par = CamusPar)
    CamusPar$Val <- gsub("[^ ]+ ", "", CamusPar$Par)
    CamusPar$Par <- gsub(" .+", "", CamusPar$Par)
    lvls <- setNames(c("ion", "psm", "peptide", "protein"), c("ION", "PSM", "PEP", "PROT"))
    w <- which(CamusPar$Par %in% names(lvls))
    l <- length(w)
    if (length(w)) {
      txt1 <- lvls[CamusPar$Par[w]]
      txt2 <- 100*as.numeric(CamusPar$Val[w])
      if (l > 1) {
        txt1 <- paste0(paste(txt1[1:(l-1)], collapse = ", "), " and ", txt1[l])
        txt2 <- paste0(paste(txt2[1:(l-1)], collapse = ", "), " and ", txt2[l])
      }
      SearchTxt <- paste0(SearchTxt, " Results were filtered in Philosopher at ", txt1, " level", c("", "s")[(l>1)+1],
                          " at FDR", c("", "s")[(l>1)+1], " ", txt2, "%", c("", ", respectively")[(l>1)+1], ".")
    }
  }
  # PTM Shepherd
  if (David) {
    SearchTxt <- gsub("\\.$", ", then PTM_Shepperd was run to localize Open search identified mass shifts." , SearchTxt)
    # Basic for now, of course eventually this should be expended.
  }
  # MS1 Quant
  if (Kant) {
    txt <- c()
    if (Jon) { txt <- c(txt, paste0("IonQuant", c("", " with match-between-runs turned on")[Matt+1])) }
    if (Fritz) { txt <- c(txt, "FreeQuant") }
    if (length(txt) > 1) { txt <- paste0("both ", paste(txt, collapse = paste0(" and ", c("", "with ")[Matt+1]))) }
    SearchTxt <- paste0(SearchTxt, " MS1-level peptide quantitation was performed using ", txt, ".")
  }
  # MS2 Quant
  if (Timothy) { SearchTxt <- paste0(SearchTxt, " MS2 reporter intensities were measured using TMT-integrator.") }
  # DiaNN
  if (Diane) {
    Libby <- gsub(topattern("diann.library="), "", grep(topattern("diann.library="), FP_Workflow, value = TRUE))
    Libby <- gsub("\\\\", "", gsub("\\\\\\\\", "/", Libby))
    SearchTxt <- paste0(SearchTxt, " DIA quantitation was performed using DiaNN ",
                        c("in library-free mode",
                          paste0("using the following library: \"", Libby, "\" (see supplementary information)"))[(nchar(Libby) > 0)+1], ".")
  }
}
MatMetCalls$Texts$DatAnalysis <- SearchTxt
mSft <- match(SearchSoft, SearchSoftware)
#
# 4) Post-processing
# ...
if (scrptType == "noReps") {
  # ... for the no-reps script this is done in one block...
  DatAnalysisTxt <- paste0(names(SearchSoft), "'s output was re-processed using in-house R scripts, starting from the ",
                           c("evidence.txt", "main report", "psm.tsv")[match(SearchSoft, SearchSoftware)],
                           " table", c("", "s")[(length(PSMsFl)>1)+1], ".",
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
  MatMetCalls$Texts <- c(MatMetCalls$Texts, DatAnalysisTxt)
  L <- length(MatMetCalls$Texts)
  MatMetCalls$Calls <- append(MatMetCalls$Calls, paste0("body_add_fpar(MatMet, fpar(ftext(MatMetCalls$Texts[", L, "], prop = WrdFrmt$Body_text), fp_p = WrdFrmt$just))"))
  MatMetCalls$Calls <- append(MatMetCalls$Calls, "body_add_par(MatMet, \"\", style = \"Normal\")")
}
if (scrptType == "withReps") {
  # ... whereas for the reps script this is started here and expanded during subsequent parts of the script...
  DatAnalysisTxt <- paste0(names(SearchSoft), "'s output was re-processed using in-house R scripts, starting from the ",
                           c("evidence.txt", "main report", "psm.tsv")[mSft],
                           " table", c("", "s")[(length(PSMsFl)>1)+1], ".")
}
