##############################
# Create output Excel tables #
##############################
#
# Write the main peptidoforms- and protein groups-level, multi-tabs report.
#
# Create openxlsx2 styles
#   It may make sense from the way Excel works, but I HATE how openxlsx2 deals with styles!
#   Anyway... 
#   So. We will. CHEAT!
#   I have saved a dummy tab with my old openxlsx styles,
#   which I will load in openxlsx2 to get and copy the styles from.
MakeRatios <- TRUE # (Used by the sourced, core sub-script)
intNms <- function(nms, topLvl = FALSE, type = "PG") {
  m <- match(type, c("pep", "PG"))
  root <- c("Intensity", "Expression")[m]
  mode <- topLvl+1
  sapply(nms, function(nm) {
    if (nm %in% c("Original", "Intensity", "Expression", "Original int.", "Intensity int.", "Expression int.")) {
      nm <- c(c("int.", "expr.")[m], root)[mode]
    } else {
      if (nm %in% c("ReNorm.", "ReNorm. int.")) { nm <- "re-norm" } else { nm <- substr(nm, 1, min(c(3, nchar(nm)))) }
      nm <- paste0(nm, ". ", c(c("int.", "expr.")[m], root)[mode])
    }
    paste0("log10(", nm, ")")
  })
}
ratNms <- function(nms, topLvl = FALSE) {
  mode <- topLvl+1
  sapply(nms, function(nm) {
    if (nm %in% c("Original", "Ratios", "Original rat.", "Ratios rat.")) {
      nm <- c("rat.", "Ratio")[mode]
    } else {
      if (nm %in% c("ReNorm.", "ReNorm. rat.")) { nm <- "Re-norm" } else { nm <- substr(nm, 1, min(c(3, nchar(nm)))) }
      nm <- paste0(nm, ". ", c("rat.", "ratios")[mode])
    }
    paste0("log2(", nm, ")")
  })
}
for (nm in names(pep.ref)) { #nm <- names(pep.ref[1])
  rpl <- intNms(nm, type = "pep")
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
for (nm in names(Prot.Expr.Root)) { #nm <- names(pep.ref[1])
  rpl <- intNms(nm)
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
for (nm in unique(c(names(pep.ratios.ref), names(Prot.Rat.Root)))) { #nm <- names(pep.ratios.ref[1])
  rpl <- ratNms(nm)
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Ratios"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Ratios"
}
fl <- system.file("extdata", "Report - column names - with replicates.xlsx", package = "proteoCraft")
styleNms <- openxlsx2::read_xlsx(fl, "tmp", colNames = FALSE)[,1]
WorkBook %<o% wb_load(fl)
repFl <- paste0(wd, "/Tables/Report_", dtstNm, ".xlsx")
WorkBook <- wb_add_data(WorkBook, "Description", dtstNm, wb_dims(2, 5))
WorkBook <- wb_add_data(WorkBook, "Description", format(Sys.Date(), "%d/%m/%Y"), wb_dims(3, 5))
WorkBook <- wb_add_data(WorkBook, "Description", WhoAmI, wb_dims(4, 5))
tmp <- loadedPackages(TRUE)
WorkBook <- wb_add_data(WorkBook, "Description", tmp$Version[grep("proteoCraft", tmp$Name)], wb_dims(5, 5))
WorkBook <- wb_set_base_font(WorkBook, 11, font_name = "Calibri")
cat(" - Writing Excel report...\n")
#
# Function for editing our header
KolEdit <- function(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl) {
  #KolNames <- xlTabs[[sheetnm]]
  klnms <- KolNames
  KolNames <- gsub("Peptides?", "Pep.", KolNames)
  KolNames <- gsub("Evidences?", "PSMs", KolNames)
  KolNames <- gsub("Spectr((al)|(um))", "Spec.", KolNames)
  KolNames <- gsub("Razor", "Raz.", KolNames)
  KolNames <- gsub("Unique", "Uniq.", KolNames)
  KolNames <- gsub("MS/MS", "MS2", KolNames)
  # This would be the place to edit PER sample evidence and MS/MS columns,
  # currently not needed because those do not exist yet
  for (nm in names(intTbl)) { #nm <- names(intTbl)[1] #nm <- names(intTbl)[2]
    m <- match(intTbl[[nm]]$Log, KolNames)
    w <- which(!is.na(m))
    if (length(w)) {
      rpl <- intNms(nm, type = "pep")
      KolNames[m[w]] <- paste0(rpl, " ", intTbl[[nm]]$Sample[w])
    }
  }
  if (!missing("ratTbl")) { # Actually, should never be missing in this workflow!
    for (nm in names(ratTbl)) {
      m <- match(ratTbl[[nm]]$Log, KolNames)
      w <- which(!is.na(m))
      if (length(w)) {
        rpl <- ratNms(nm)
        KolNames[m[w]] <- paste0(rpl, " ", ratTbl[[nm]]$Sample[w])
      }
    }
  }
  wNF <- grep("^mod\\. F-test ", KolNames, invert = TRUE)
  if (F.test) {
    wF <- grep("^mod\\. F-test ", KolNames)
    KolNames[wF] <- gsub("^mod\\. F-test +", "F-test ", KolNames[wF]) # Shorter F-test tag
    KolNames[wF] <- gsub(" +-log10\\(Pvalue\\)( - )?", " -log10 pval. ", KolNames[wF])
    KolNames[wF] <- gsub(" +Regulated - ", " reg. ", KolNames[wF])
    KolNames[wF] <- gsub(" +Significant-", " signif. ", KolNames[wF])
  }
  KolNames[wNF] <- gsub(".*-log10\\(Pvalue\\)( - )?", "-log10 pval. ", KolNames[wNF])
  KolNames[wNF] <- gsub(".*Significant-", "signif. ", KolNames[wNF])
  KolNames[wNF] <- gsub(".*Regulated - ", "reg. ", KolNames[wNF])
  KolNames[wNF] <- gsub(".*Significant-", "signif. ", KolNames[wNF])
  KntKol <- paste0(AA, " Count")
  KolNames[which(KolNames %in% KntKol)] <- gsub(" Count$", "", KolNames[which(KolNames %in% KntKol)])
  #
  KolNames <- gsub("( - )|(___)", " ", KolNames)
  #
  # F-test
  g <- grep("F-test: ", klnms)
  if (length(g)) {
    KolNames[g] <- paste0("F-test ", KolNames[g])
    g <- grep("F-test .*F(-| |\\.|_)?test", KolNames)
    stopifnot(length(g) == 0)
  }
  #
  # Those names must be unique if the data is to be written as a table!
  # Which is annoying, because this limits how much fat we can cut
  tst <- aggregate(KolNames, list(KolNames), c)
  tst$L <- vapply(tst$x, length, 1)
  tst <- tst[which(tst$Group.1 != ""),]
  stopifnot(max(tst$L) == 1)
  #tst$x[which(tst$L > 1)]
  #
  KolNames <- as.data.frame(t(KolNames))
  colnames(KolNames) <- klnms
  return(KolNames)
}
#KolEdit(xlTabs[[sheetnm]], intColsTbl, ratColsTbl)
if ((prot.list.Cond)&&(!"In list" %in% colnames(ev))) {
  g <- grsep2(prot.list, ev$Proteins)
  w <- rep(FALSE, nrow(ev))
  w[g] <- TRUE
  ev$"In list" <- c("", "+")[w+1]
}
QualFilt %<o% c("In list", "Potential contaminant", "Only identified by site",
                grep("^Quality filter: ", colnames(PG), value = TRUE),
                grep("^Quantity Quality$", colnames(pep), value = TRUE))
if ((DiscFilt)&&(DiscFiltMode == DiscFiltModes[3])) { QualFilt <- c(QualFilt, DiscFiltCols) }
II <- setNames(1, "All peptidoforms")
if ((exists("PTMs_pep"))&&(length(PTMs_pep))) {
  Mod2Write <- names(PTMs_pep)
  II[paste0(Mod2Write, "-mod. pept.")] <- 1+(seq_along(length(Mod2Write)))
}
for (ii in II) { #ii <- II[1] #ii <- II[2]
  tblMode <- tblMode2 <- "pep"
  TbNm <- names(II)[ii]
  tempData <- get(tblMode)
  intRf <- pep.ref
  ratRf <- pep.ratios.ref
  if (ii > 1) {
    Ptm <- Mod2Write[[ii-1]]
    ptm <- Modifs$Mark[match(Ptm, Modifs$`Full name`)]
    tempData <- PTMs_pep[[Ptm]]
    tempData$Name <- gsub("\n", " ", tempData$Name)
    intRf <- PTMs_int.ref[[Ptm]]
    ratRf <- PTMs_rat.ref[[Ptm]]
    tblMode2 <- paste0(Ptm, "-modified pep")
  }
  names(intRf) <- paste0(names(intRf), " int.")
  names(ratRf) <- paste0(names(ratRf), " rat.")
  if (nrow(tempData)) {
    CoreCol %<o% c("id", "Modified sequence")
    if ("Modified sequence_verbose" %in% colnames(tempData)) { CoreCol <- c(CoreCol, "Modified sequence_verbose") }
    CoreCol <- c(CoreCol, "Sequence", "Proteins")
    if (ii > 1) { CoreCol <- c(CoreCol, "Name", paste0(Ptm, "-site(s)")) }
    CoreCol2 %<o% c("Leading proteins", "Leading razor proteins",
                    "Protein names", "Gene names", "Protein group IDs", "Razor protein group ID",
                    grep(" (Probabilities|Score Diffs)$", colnames(tempData), value = TRUE),
                    "Normalisation group")
    evcol <- "Evidence IDs"
    spcol <- "MS/MS count"
    gel <- setNames(lapply(intRf, function(rf) {
      x <- c(paste0(rf, RSA$values),
             paste0("Mean ", rf, VPAL$values))
      return(x[which(x %in% colnames(tempData))])
    }), intRf)
    if (ii == 1) {
      # Log transform for normal tables - not necessary for PTM-modified tables as we already transformed
      gel2 <- setNames(lapply(intRf, function(rf) {
        x <- c(paste0(rf, RSA$values),
               paste0("Mean ", rf, VPAL$values))
        w <- which(x %in% colnames(tempData))
        rf2 <- gsub("Evidence intensities - ", "log10(Int.) - ", rf)
        x <- c(paste0(rf2, RSA$values),
               paste0("Mean ", rf, VPAL$values))
        return(x[w])
      }), intRf)
      for (rf in intRf) { #rf <- intRf[1]
        kol1 <- gel[[rf]]
        kol2 <- gel2[[rf]]
        if (length(kol1)) {
          tempData[, kol2] <- suppressWarnings(log10(tempData[, kol1]))
          for (kl in kol2) {
            w <- which(is.infinite(tempData[[kl]]))
            tempData[w, kl] <- NA
          }
        }
      }
      intRf <- sapply(intRf, function(rf) {
        gsub("Evidence intensities - ", "log10(Int.) - ", rf)
      })
      gel <- unlist(gel2)
    }
    #
    smpls <- c(VPAL$values, RSA$values)
    if (length(Exp) == 1) { smpls <- gsub(topattern(paste0(Exp, "___")), "", smpls) }
    smpls <- gsub("___", " ", smpls)
    intColsTbl <- setNames(lapply(names(intRf), function(nm) { #nm <- names(intRf)[1]
      res <- data.frame(Log = c(paste0("Mean ", intRf[nm], VPAL$values),
                                paste0(intRf[nm], RSA$values)),
                        Type = c(rep("Average", length(VPAL$values)),
                                 rep("Individual", length(RSA$values))),
                        Sample = smpls)
      w <- which(res$Log %in% colnames(tempData))
      return(res[w,])
    }), names(intRf))
    w <- which(vapply(intColsTbl, nrow, 1) > 0)
    intColsTbl <- intColsTbl[w]; intRf <- intRf[w]
    quantCols <- intCols <- lapply(intColsTbl, function(x) { x$Log })
    ratColsTbl <- setNames(lapply(names(ratRf), function(nm) {
      res <- data.frame(Log = c(paste0("Mean ", ratRf[nm], VPAL$values),
                                paste0(ratRf[nm], RSA$values)),
                        Type = c(rep("Average", length(VPAL$values)),
                                 rep("Individual", length(RSA$values))),
                        Sample = smpls)
      w <- which(res$Log %in% colnames(tempData))
      return(res[w,])
    }), names(ratRf))
    w <- which(vapply(ratColsTbl, nrow, 1) > 0)
    ratColsTbl <- ratColsTbl[w]; ratRf <- ratRf[w]
    ratCols <- lapply(ratColsTbl, function(x) { x$Log })
    grl <- unlist(ratCols)
    for (gr in grl) {
      w <- which(is.infinite(tempData[[gr]]))
      tempData[w, gr] <- NA
    }
    quantCols[names(ratRf)] <- ratCols
    regcol <- grep("^((Enriched)|(Regulated)) - ", colnames(tempData), value = TRUE)
    signcol <- grep("^Significant-FDR=[1-9][0-9]*\\.*[0-9]*% - ", colnames(tempData), value = TRUE)
    signcol <- grep(" - Analysis_[0-9]+", signcol, invert = TRUE, value = TRUE)
    quantcol <- unlist(quantCols)
    PepColList %<o% c("gel", "grl", "quantcol", "signcol", "regcol") # These are any column for which we want to gsub "___" to " "
    .obj <- unique(c(PepColList, .obj)) # Here easier than using a custom operator
    if (ii > 1) {
      gpl <- grep(topattern(pvalue.col[which(pvalue.use)]), colnames(tempData), value = TRUE)
      quantcol <- c(quantcol, gpl)
      PepColList <- c(PepColList, "gpl")
    }
    aacol <- paste0(AA, " Count")
    qualFlt <- QualFilt[which(QualFilt %in% colnames(ev))]
    w <- which(!qualFlt %in% colnames(tempData))
    if (length(w)) {
      tempData[, qualFlt[w]] <- ev[match(tempData$"Modified sequence", ev$"Modified sequence"), qualFlt[w]]
    }
    kol <- c(CoreCol, "In list", CoreCol2, evcol, spcol, "PEP", quantcol, signcol, regcol, qualFlt[which(qualFlt != "In list")], aacol)
    if (ii > 1) { kol <- c(kol, "Code") }
    if (Annotate) {
      PepAnnotCol %<o% annot.col[which(annot.col %in% colnames(tempData))]
      kol <- c(kol, PepAnnotCol)
    }
    #tst <- data.frame(Names = names(kol), Column = setNames(kol, NULL), Found = kol %in% colnames(tempData));View(tst)
    kol <- kol[which(kol %in% colnames(tempData))]
    #kol[which(!kol %in% colnames(tempData))]
    #colnames(tempData)[which(!colnames(tempData) %in% kol)]
    tempData <- tempData[, kol]
    # If there is only one experiment, remove it from the names here...
    colnames(tempData) <- cleanNms(colnames(tempData), start = FALSE)
    for (i in PepColList) { assign(i, cleanNms(get(i), start = FALSE)) }
    intColsTbl <- lapply(intColsTbl, function(x) {
      x$Log <- cleanNms(x$Log, start = FALSE)
      x
    })
    ratColsTbl <- lapply(ratColsTbl, function(x) {
      x$Log <- cleanNms(x$Log, start = FALSE)
      x
    })
    intCols <- lapply(intCols, cleanNms, start = FALSE)
    ratCols <- lapply(ratCols, cleanNms, start = FALSE)
    quantCols <- lapply(quantCols, cleanNms, start = FALSE)
    colnames(tempData) <- gsub("_names$", " names", colnames(tempData))
    if (Annotate) { PepAnnotCol <- gsub("_names$", " names", PepAnnotCol) }
    for (k in regcol) {
      tempData[which(tempData[[k]] == "non significant"), k] <- "n.s."
      tempData[which(tempData[[k]] == ""), k] <- "n.t."
    }
    if ((ii > 1)&&(F.test)) {
      tempPepF <- PTMs_F_test_data[[Ptm]]
      tempPepF <- tempPepF[, which(!colnames(tempPepF) %in% c(Param$Plot.labels, "Rel. log10(Peptides count)", "Av. log10 abundance"))]
      colnames(tempPepF) <- cleanNms(colnames(tempPepF), start = FALSE)
      #mnratcolF %<o% grep("Mean log2\\(Ratio\\) - ", colnames(tempPepF), value = TRUE)
      #m <- match(mnratcolF, colnames(tempPepF))
      #colnames(tempPepF)[m] <- paste0("mod. F-test ", mnratcolF)
      #mnratcolF <- paste0("mod. F-test ", mnratcolF)
      pvalcolF %<o% F_Root
      signcolF %<o% grep("^mod\\. F-test Significant", colnames(tempPepF), value = TRUE)
      regcolF %<o% grep("^mod\\. F-test Regulated", colnames(tempPepF), value = TRUE)
      Fkol %<o% c(regcolF, #mnratcolF,
                  pvalcolF, signcolF)
      for (k in regcolF) {
        tempPepF[which(tempPepF[[k]] == "non significant"), k] <- "n.s."
        tempPepF[which(tempPepF[[k]] == ""), k] <- "n.t."
      }
      tempData[, Fkol] <- tempPepF[match(tempData$"Modified sequence", tempPepF$"Modified sequence"), Fkol]
    }
    dir <- paste0(wd, "/Tables")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    data.table::fwrite(tempData, paste0(dir, "/", TbNm, ".tsv"), sep = "\t", row.names = FALSE, na = "NA")
    w <- grsep2(prot.list, tempData$Proteins)
    if (length(w)) {
      data.table::fwrite(tempData[w,], paste0(wd, "/Tables/", TbNm, " - Proteins in list.tsv"),
                         sep = "\t", row.names = FALSE, na = "NA")
    }
    # Which columns are affected by each style
    # - IDs
    ColumnsTbl %<o% list(IDs = c(CoreCol, CoreCol2))
    # - Counts
    ColumnsTbl[["AA counts"]] <- aacol
    # - Evidence counts and IDs
    ColumnsTbl[["Global Ev. IDs"]] <- "Evidence IDs"
    ColumnsTbl[["Global Spec. counts"]] <- "MS/MS count"
    # - Individual Expr
    ColumnsTbl[["Individual Expr"]] <- grep("^Mean ", unlist(intCols), value = TRUE, invert = TRUE)
    # - Summary Expr
    ColumnsTbl[["Summary Expr"]] <- grep("^Mean ", unlist(intCols), value = TRUE)
    # - Individual Ratios
    ColumnsTbl[["Individual Ratios"]] <- grep("^Mean ", unlist(ratCols), value = TRUE, invert = TRUE)
    # - Summary Ratios
    ColumnsTbl[["Summary Ratios"]] <- grep("^Mean ", unlist(ratCols), value = TRUE)
    if (ii > 1) {
      # - Summary Ratios: P-values and significance
      ColumnsTbl[["P-values"]] <- gpl
      # - Significant
      ColumnsTbl[["Significant"]] <- signcol
      # - Regulated
      ColumnsTbl[["Regulated"]] <- regcol
      # F-test
      if ((ii > 1)&&(F.test)) {
        #ColumnsTbl[["F-test summary Ratios"]] <- mnratcolF
        ColumnsTbl[["F-test P-values"]] <- pvalcolF
        ColumnsTbl[["F-test significant"]] <- signcolF
        ColumnsTbl[["F-test regulated"]] <- regcolF
      }
    }
    # - Annotations
    if (Annotate) {
      AnnotTbl$Columns <- list(c("GO", "GO-ID"), c("Taxonomy", "TaxID"), NA, NA, NA, NA, "EMBL", NA)
      for (i in annot) { AnnotTbl$Columns[match(i, AnnotTbl$Name)] <- list(c(i, paste0(i, " names"))) }
      AnnotTbl$Columns[match("Other", AnnotTbl$Name)] <- list(annot.col2[which(!annot.col2 %in% unlist(AnnotTbl$Columns))])
      for (i in 1:nrow(AnnotTbl)) { ColumnsTbl[[paste0(AnnotTbl$Name[i], " annotations")]] <- AnnotTbl$Columns[[i]] }
    }
    # - PEP
    ColumnsTbl[["PEP"]] <- "PEP"
    # - Filters
    ColumnsTbl[["Filters"]] <- qualFlt
    # Melt
    ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }, 1) > 0)]
    ColumnsTbl <- set_colnames(reshape::melt.list(ColumnsTbl), c("Col", "Grp"))
    #tst <- aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), length); View(tst)
    #aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), unique)
    stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
    ColumnsTbl$Class <- ""
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Peptides information"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evcol))] <- "Evidence IDs"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(spcol))] <- "Spectral count"
    for (nm in names(intRf)) { #nm <- names(intRf)[1]
      rpl <- intNms(nm, TRUE, type = "pep")
      ColumnsTbl$Class[grep(topattern(intRf[nm]), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[grep(topattern(paste0("Mean ", intRf[nm])), ColumnsTbl$Col)] <- rpl
    }
    for (nm in names(ratRf)) { #nm <- names(ratRf)[1]
      rpl <- ratNms(nm, TRUE)
      ColumnsTbl$Class[grep(topattern(ratRf[nm]), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[grep(topattern(paste0("Mean ", ratRf[nm])), ColumnsTbl$Col)] <- rpl
    }
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "P-values")] <- gsub(" - $", "", pvalue.col[which(pvalue.use)])
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% signcol)] <- "Significant"
    if ((ii > 1)&&(F.test)) {
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test summary Ratios")] <- "F-test"
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test P-values")] <- "F-test"
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test significant")] <- "F-test"
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test regulated")] <- "F-test"
    }
    ColumnsTbl$Class[grep("[Aa]nnotations", ColumnsTbl$Grp)] <- "Annotations"
    ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters"))] <- "QC filters"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% aacol)] <- "Amino Acid counts"
    ColumnsTbl$Hide <- ColumnsTbl$Class %in% c("Spectral count", "Spectrum IDs", "Amino Acid counts", "Annotations", "Cluster (hierarch.)")
    #
    if (MakeRatios) { a <- KolEdit(ColumnsTbl$Col, intColsTbl, ratColsTbl) } else { a <- KolEdit(ColumnsTbl$Col, intColsTbl) }
    ColumnsTbl$edit_Col <- unlist(a)
    #
    Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #saveFun(WorkBook, file = "WorkBook_bckp.RData")
    #wb_save(WorkBook, paste0(wd, "/tst.xlsx")); xl_open(paste0(wd, "/tst.xlsx"))
    #loadFun("WorkBook_bckp.RData")
  }
}
TbNm <- "Protein groups"
tblMode <- tblMode2 <- "PG"
# Function for editing the header
KolEdit <- function(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl, locTbl = SubCellLocTbl) {
  #KolNames <- xlTabs[[sheetnm]]
  klnms <- KolNames
  KolNames <- gsub("Peptides?", "Pep.", KolNames)
  KolNames <- gsub("Evidences?", "PSMs", KolNames)
  KolNames <- gsub("Spectr((al)|(um))", "Spec.", KolNames)
  KolNames <- gsub("Razor", "Raz.", KolNames)
  KolNames <- gsub("Unique", "Uniq.", KolNames)
  KolNames <- gsub("MS/MS", "MS2", KolNames)
  for (nm in names(intTbl)) { #nm <- names(intTbl)[1] #nm <- names(intTbl)[2]
    m <- match(intTbl[[nm]]$Log, KolNames)
    w <- which(!is.na(m))
    if (length(w)) {
      rpl <- intNms(nm)
      KolNames[m[w]] <- paste0(rpl, " ", intTbl[[nm]]$Sample[w])
    }
  }
  if (!missing("ratTbl")) { # Actually, should never be missing in this workflow!
    for (nm in names(ratTbl)) {
      m <- match(ratTbl[[nm]]$Log, KolNames)
      w <- which(!is.na(m))
      if (length(w)) {
        rpl <- ratNms(nm)
        KolNames[m[w]] <- paste0(rpl, " ", ratTbl[[nm]]$Sample[w])
      }
    }
  }
  if (!missing("locTbl")) {
    for (nm in names(locTbl)) {
      m <- match(locTbl[[nm]]$Columns, KolNames)
      w <- which(!is.na(m))
      if (length(w)) {
        KolNames[m[w]] <- paste0(nm, locTbl[[nm]]$Sample[w])
      }
    }
  }
  wNF <- grep("^mod\\. F-test ", KolNames, invert = TRUE)
  if (F.test) {
    wF <- grep("^mod\\. F-test ", KolNames)
    KolNames[wF] <- gsub("^mod\\. F-test +", "F-test ", KolNames[wF]) # Shorter F-test tag
    KolNames[wF] <- gsub(" +-log10\\(Pvalue\\)( - )?", " -log10 pval. ", KolNames[wF])
    KolNames[wF] <- gsub(" +Regulated - ", " reg. ", KolNames[wF])
    KolNames[wF] <- gsub(" +Significant-", " signif. ", KolNames[wF])
  }
  KolNames <- gsub(".*Pvalue\\)( - )?", "-log10 pval. ", KolNames)
  KolNames <- gsub(".*Significant-", "signif. ", KolNames)
  KolNames <- gsub(".*Regulated - ", "reg. ", KolNames)
  KolNames <- gsub("log10\\(est\\. copies/cell\\) (- )?", "ProtRul. ", KolNames)
  KolNames <- gsub(paste0(topattern("Sequence coverage [%] ", start = FALSE), "(- )?"), "Cov. ", KolNames)
  KolNames[which(KolNames == "Max. theoretical sequence coverage [%]")] <- "Theoretical max."
  KolNames[which(KolNames == "Sequence coverage [%]")] <- "All peptides"
  KolNames[which(KolNames == "Uniq. + razor sequence coverage [%]")] <- "Unique + razor"
  KolNames[which(KolNames == "Uniq. sequence coverage [%]")] <- "Unique"
  KolNames <- gsub("^Cluster \\([^\\)]+\\) - ", "Clust. ", KolNames)
  #
  KolNames <- gsub("( - )|(___)", " ", KolNames)
  #
  # F-test
  g <- grep("F-test: ", klnms)
  if (length(g)) {
    KolNames[g] <- paste0("F-test ", KolNames[g])
    g <- grep("F-test .*F(-| |\\.|_)?test", KolNames)
    stopifnot(length(g) == 0)
  }
  #
  # Those names must be unique if the data is to be written as a table!
  # Which is annoying, because this limits how much fat we can cut
  tst <- aggregate(KolNames, list(KolNames), c)
  tst$L <- vapply(tst$x, length, 1)
  tst <- tst[which(tst$Group.1 != ""),]
  stopifnot(max(tst$L) == 1)
  #tst$x[which(tst$L > 1)]
  #
  KolNames <- as.data.frame(t(KolNames))
  colnames(KolNames) <- klnms
  return(KolNames) #View(KolNames)
}
#
tempData <- get(tblMode)
CoreCol <- c("Leading protein IDs", "Common Names", "Genes")
CoreCol2 <- c("Protein IDs", "Names", "id", "Mol. weight [kDa]")
pepevspeccol <- c("Peptides count",
                  grep("^Peptides count - ", colnames(tempData), value = TRUE),
                  "Peptide IDs", "Razor peptide IDs", "Unique peptide IDs",
                  grep("^Peptide IDs - ", colnames(tempData), value = TRUE),
                  "psmS count",
                  grep("^Evidences count - ", colnames(tempData), value = TRUE),
                  "Evidence IDs",
                  grep("^Evidence IDs - ", colnames(tempData), value = TRUE),
                  "Spectral count",
                  grep("^Spectral count - ", colnames(tempData), value = TRUE),
                  "Spectrum IDs",
                  grep("^Spectrum IDs - ", colnames(tempData), value = TRUE),
                  "Biot. peptides count",
                  grep("^Biot\\. peptides count - ", colnames(tempData), value = TRUE),
                  "Biot. peptide IDs",
                  grep("^Biot\\. peptides IDs - ", colnames(tempData), value = TRUE),
                  "Biot. peptides [%]",
                  "Biot. evidences count",
                  grep("^Biot\\. evidences count - ", colnames(tempData), value = TRUE),
                  "Biot. evidence IDs",
                  grep("^Biot\\. evidence IDs - ", colnames(tempData), value = TRUE),
                  "Biot. spectral count",
                  grep("^Biot\\. spectral count - ", colnames(tempData), value = TRUE),
                  "Biot. spectrum IDs",
                  grep("^Biot\\. spectrum IDs - ", colnames(tempData), value = TRUE))
pepevspeccol <- pepevspeccol[which(pepevspeccol %in% colnames(tempData))]
pepcountcol1 <- grep("[Pp]eptides count$", pepevspeccol, value = TRUE)
pepcountcol2 <- grep("[Pp]eptides count - ", pepevspeccol, value = TRUE)
pepidcol1 <- grep("[Pp]eptide IDs$", pepevspeccol, value = TRUE)
pepidcol2 <- grep("[Pp]eptide IDs - ", pepevspeccol, value = TRUE)
evcountcol1 <- grep("[Ev]vidences count$", pepevspeccol, value = TRUE)
evcountcol2 <- grep("[Ev]vidences count - ", pepevspeccol, value = TRUE)
evidcol1 <- grep("[Ev]vidence IDs$", pepevspeccol, value = TRUE)
evidcol2 <- grep("[Ev]vidence IDs - ", pepevspeccol, value = TRUE)
speccountcol1 <- grep("[Ss]pectral count$", pepevspeccol, value = TRUE)
speccountcol2 <- grep("[Ss]pectral count - ", pepevspeccol, value = TRUE)
specidcol1 <- grep("[Ss]pectrum IDs$", pepevspeccol, value = TRUE)
specidcol2 <- grep("[Ss]pectrum IDs - ", pepevspeccol, value = TRUE)
pepevspeccol <- c(pepcountcol1, pepcountcol2, pepidcol1, pepidcol2,
                  evcountcol1, evcountcol2, evidcol1, evidcol2,
                  speccountcol1, speccountcol2, specidcol1, specidcol2)
pepevspeccola <- grep("((([Ss]pectral|[Pp]eptides|[Ee]vidences) count)|(([Ss]pectrum|[Pp]eptide|[Ee]vidence) IDs))$", pepevspeccol, value = TRUE)
pepevspeccolb <- grep("((([Ss]pectral|[Pp]eptides|[Ee]vidences) count)|(([Ss]pectrum|[Pp]eptide|[Ee]vidence) IDs)) - ", pepevspeccol, value = TRUE)
kol <- c(CoreCol, "In list", CoreCol2, pepevspeccol)
intRf <- Prot.Expr.Root
names(intRf) <- paste0(names(intRf), " int.")
intColsTbl <- setNames(lapply(names(intRf), function(nm) { #nm <- names(intRf)[1]
  res <- data.frame(Log = c(paste0("Mean ", intRf[nm], VPAL$values),
                            paste0(intRf[nm], RSA$values)),
                    Type = c(rep("Average", length(VPAL$values)),
                             rep("Individual", length(RSA$values))),
                    Sample = smpls)
  w <- which(res$Log %in% colnames(tempData))
  return(res[w,])
}), names(intRf))
w <- which(vapply(intColsTbl, nrow, 1) > 0)
intColsTbl <- intColsTbl[w]; intRf <- intRf[w]
quantCols <- intCols <- lapply(intColsTbl, function(x) { x$Log })
gel <- unlist(intCols)
for (gl in gel) {
  w <- which(is.infinite(tempData[[gl]]))
  tempData[w, gl] <- NA
}
quantcol <- gel
ratRf <- Prot.Rat.Root
names(ratRf) <- paste0(names(ratRf), " rat.")
ratColsTbl <- setNames(lapply(names(ratRf), function(nm) {
  res <- data.frame(Log = c(paste0("Mean ", ratRf[nm], VPAL$values),
                            paste0(ratRf[nm], RSA$values)),
                    Type = c(rep("Average", length(VPAL$values)),
                             rep("Individual", length(RSA$values))),
                    Sample = smpls)
  w <- which(res$Log %in% colnames(tempData))
  return(res[w,])
}), names(ratRf))
w <- which(vapply(ratColsTbl, nrow, 1) > 0)
ratColsTbl <- ratColsTbl[w]; ratRf <- ratRf[w]
ratCols <- lapply(ratColsTbl, function(x) { x$Log })
grl <- unlist(ratCols)
for (gr in grl) {
  w <- which(is.infinite(tempData[[gr]]))
  tempData[w, gr] <- NA
}
quantcol <- c(quantcol, grl)
quantCols[names(ratRf)] <- ratCols
if (protrul) {
  tmp <- grep(topattern("log10(est. copies/cell) - ", start = FALSE), colnames(tempData), value = TRUE)
  quantcol <- c(quantcol, tmp)
  quantCols[["Proteomic ruler"]] <- tmp
}
pvalcol <- grep(topattern(pvalue.col[pvalue.use]), colnames(tempData), value = TRUE)
regcol <- grep("^((Enriched)|(Regulated)) - ", colnames(tempData), value = TRUE)
signcol <- grep("^Significant-FDR=[1-9][0-9]*\\.*[0-9]*% - ", colnames(tempData), value = TRUE)
signcol <- grep(" - Analysis_[0-9]+", signcol, invert = TRUE, value = TRUE)
covcol <- c(xmlCovCol,
            c("Sequence coverage [%]",
              "Unique + razor sequence coverage [%]",
              "Unique sequence coverage [%]")[1:c(1, 3)[isEukaLike+1]],
            grep(topattern("Sequence coverage [%] - "), colnames(tempData), value = TRUE)) # The complicated way, but ensures the order is correct
kol <- c(kol, "Mol. weight [kDa]", covcol, "PEP", covcol, quantcol, pvalcol, regcol, signcol)
if ((exists("KlustKols"))&&(length(KlustKols))) { kol <- c(kol, KlustKols) }
qualFlt <- QualFilt
if (length(GO_PG_col)) { qualFlt <- c(qualFlt, GO_PG_col2) }
kol <- c(kol, qualFlt[which(qualFlt != "In list")])
if (Annotate) { kol <- c(kol, annot.col) }
if (Annotate&&LocAnalysis) {
  PG$Marker <- ""
  w <- which(PG$Label %in% names(SubCellMark2))
  PG$Marker[w] <- SubCellMark2[match(PG$Label[w], names(SubCellMark2))]
  SubCellLocTbl <- list()
  lokol1 <- paste0("Localisation - ", SubCellFracAggr$values)
  w <- which(lokol1 %in% colnames(PG))
  lokol1 <- lokol1[w]
  SubCellLocTbl$Localisation <- list(Columns = lokol1,
                                     Sample = SubCellFracAggr$values[w])
  kol <- c(kol, "Marker", lokol1)
  if (LocAnalysis2) {
    rt2 <- SSD.Root
    lokol2 <- grep(topattern(rt2), colnames(PG), value = TRUE)
    rt3 <- paste0("Mean ", SSD.Root)
    lokol3 <- grep(topattern(rt3), colnames(PG), value = TRUE)
    SubCellLocTbl$SSDs <- list(Columns = c(lokol3, lokol2),
                               Sample = c(gsub(topattern(rt3), "", lokol3),
                                          gsub(topattern(rt2), "", lokol2)))
    rt <- SSD.Pval.Root
    lokol4 <- grep(topattern(rt), colnames(PG), value = TRUE)
    SubCellLocTbl$"SSD P-vals." <- list(Columns = lokol4,
                                        Sample = gsub(topattern(rt), "", lokol4))
    rt <- "Signif. SSDs-FDR="
    lokol5 <- grep(topattern(rt), colnames(PG), value = TRUE)
    SubCellLocTbl$"Signif. SSDs" <- list(Columns = lokol5,
                                         Sample = gsub(topattern(rt), "", lokol5))
    rt <- "Re-localized"
    lokol6 <- grep(topattern(rt), colnames(PG), value = TRUE)
    SubCellLocTbl$"Re-Loc." <- list(Columns = lokol6,
                                    Sample = gsub(topattern(rt), "", lokol6))
    kol <- c(kol, lokol6, lokol2, lokol3, lokol4, lokol5)
  }
} else { SubCellLocTbl <- NULL}
kol <- unique(kol[which(kol %in% colnames(tempData))])
tempData <- tempData[, kol]
#
if (F.test) {
  tmpPGf <- F_test_data
  tmpPGf <- tmpPGf[, which(!colnames(tmpPGf) %in% c(Param$Plot.labels, "Rel. log10(Peptides count)", "Av. log10 abundance"))]
  colnames(tmpPGf) <- cleanNms(colnames(tmpPGf), start = FALSE)
  #mnratcolF %<o% grep("Mean log2\\(Ratio\\) - ", colnames(tmpPGf), value = TRUE)
  #m <- match(mnratcolF, colnames(tmpPGf))
  #colnames(tmpPGf)[m] <- paste0("mod. F-test ", mnratcolF)
  #mnratcolF <- paste0("mod. F-test ", mnratcolF)
  pvalcolF %<o% F_Root
  signcolF %<o% grep("^mod\\. F-test Significant", colnames(tmpPGf), value = TRUE)
  regcolF %<o% grep("^mod\\. F-test Regulated", colnames(tmpPGf), value = TRUE)
  for (k in regcolF) {
    tmpPGf[which(tmpPGf[[k]] == "non significant"), k] <- "n.s."
    tmpPGf[which(tmpPGf[[k]] == ""), k] <- "n.t."
  }
  Fkol %<o% c(regcolF, #mnratcolF,
              pvalcolF, signcolF)
  tempData[, Fkol] <- tmpPGf[match(tempData$`Protein IDs`, tmpPGf$`Protein IDs`), Fkol]
}
if (Annotate&&LocAnalysis) {
  if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) { loadFun("GO_terms.RData") }
  GOCC <- GO_terms$ID[which(GO_terms$Ontology == "CC")]
  tempData$"GO-ID (CC)" <- lapply(strsplit(tempData$`GO-ID`, ";"), function(x) { x[which(x %in% GOCC)] })
  w <- which(vapply(tempData$"GO-ID (CC)", length, 1) > 0)
  tempData$"GO (CC)" <- ""
  tempData$"GO (CC)"[w] <- vapply(tempData$"GO-ID (CC)"[w], function(x) {
    paste(GO_terms$Term[match(x, GO_terms$ID)], collapse = ";")
  }, "")
  if (LocAnalysis2) {
    for (k in lokol6) {
      tempData[which(tempData[[k]] == "non significant"), k] <- "n.s."
      tempData[which(tempData[[k]] == ""), k] <- "n.t."
      tempData[[k]] <- gsub("^up, FDR = ", "re-loc., FDR = ", tempData[[k]])
    }
  }
}
#
# Which columns are affected by each style
# - IDs
ColumnsTbl <- list(IDs = c(CoreCol, CoreCol2))
# - Peptide and evidence counts and IDs
ColumnsTbl[["Global Pep. IDs"]] <- pepidcol1
ColumnsTbl[["Global Pep. counts"]] <- pepcountcol1
ColumnsTbl[["Pep. IDs"]] <- pepidcol2
ColumnsTbl[["Pep. counts"]] <- pepcountcol2
ColumnsTbl[["Global Ev. IDs"]] <- evidcol1
ColumnsTbl[["Global Ev. counts"]] <- evcountcol1
ColumnsTbl[["Ev. IDs"]] <- evidcol2
ColumnsTbl[["Ev. counts"]] <- evcountcol2
ColumnsTbl[["Global Spec. IDs"]] <- specidcol1
ColumnsTbl[["Global Spec. counts"]] <- speccountcol1
ColumnsTbl[["Spec. IDs"]] <- specidcol2
ColumnsTbl[["Spec. counts"]] <- speccountcol2
if (IsBioID2) {
  # Some of these may currently be empty - we don't want to overload the files with columns
  ColumnsTbl[["Global Biot. Pep. IDs"]] <- biotpepidcol1
  ColumnsTbl[["Global Biot. Pep. counts"]] <- biotpepcountcol1
  ColumnsTbl[["Biot. Pep. IDs"]] <- biotpepidcol2
  ColumnsTbl[["Biot. Pep. counts"]] <- biotpepcountcol2
  ColumnsTbl[["Global Biot. Ev. IDs"]] <- biotevidcol1
  ColumnsTbl[["Global Biot. Ev. counts"]] <- biotevcountcol1
  ColumnsTbl[["Biot. Ev. IDs"]] <- biotevidcol2
  ColumnsTbl[["Biot. Ev. counts"]] <- biotevcountcol2
  ColumnsTbl[["Biot. Pep. %"]] <- "Biot. peptides [%]"
}
# Quantitation
# - Expression values
for (nm in names(intRf)) { #nm <- names(intRf[1])
  rpl <- intNms(nm)
  ColumnsTbl[[paste0(rpl, ", avg.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Average")]
  ColumnsTbl[[paste0(rpl, ", indiv.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Individual")]
}
# - Ratios
for (nm in names(ratRf)) { #nm <- names(ratRf[1])
  rpl <- ratNms(nm)
  ColumnsTbl[[paste0(rpl, ", avg.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Average")]
  ColumnsTbl[[paste0(rpl, ", indiv.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Individual")]
}
# - P-values
ColumnsTbl[["P-values"]] <- pvalcol
# - Significant
ColumnsTbl[["Significant"]] <- signcol
# - Regulated
ColumnsTbl[["Regulated"]] <- regcol
# - Proteome Ruler
if (protrul) {
  ColumnsTbl[["Proteome Ruler"]] <- grep(topattern("log10(est. copies/cell) - ", start = FALSE), colnames(tempData), value = TRUE)
}
#
# F-test
if (F.test) {
  #ColumnsTbl[["F-test summary Ratios"]] <- mnratcolF
  ColumnsTbl[["F-test P-values"]] <- pvalcolF
  ColumnsTbl[["F-test significant"]] <- signcolF
  ColumnsTbl[["F-test regulated"]] <- regcolF
}
# - Annotations
if (Annotate) {
  annot.col2 <- gsub("_names$", " names", annot.col)
  AnnotTbl$Columns <- list(c("GO", "GO-ID"), c("Taxonomy", "TaxID"), NA, NA, NA, NA, "EMBL", NA)
  annot <- c("InterPro", "Pfam", "PIRSF", "PROSITE")
  for (i in annot) { AnnotTbl$Columns[match(i, AnnotTbl$Name)] <- list(c(i, paste0(i, " names"))) }
  AnnotTbl$Columns[match("Other", AnnotTbl$Name)] <- list(annot.col2[which(!annot.col2 %in% unlist(AnnotTbl$Columns))])
  for (i in 1:nrow(AnnotTbl)) { ColumnsTbl[[paste0(AnnotTbl$Name[i], " annotations")]] <- AnnotTbl$Columns[[i]] }
  if (LocAnalysis) {
    ColumnsTbl[["GO annotations"]] <- c(ColumnsTbl[["GO"]], "GO (CC)")
    ColumnsTbl[["GO-ID annotations"]] <- c(ColumnsTbl[["GO-ID"]], "GO-ID (CC)")
    ColumnsTbl[["Localisation"]] <- lokol1
    ColumnsTbl[["Marker"]] <- "Marker"
    if (LocAnalysis2) {
      ColumnsTbl[["SSDs"]] <- lokol2
      ColumnsTbl[["Mean SSDs"]] <- lokol3
      ColumnsTbl[["SSDs P-values"]] <- lokol4
      ColumnsTbl[["SSDs significant"]] <- lokol5
      ColumnsTbl[["Re-localized"]] <- lokol6
    }
  }
}
# - PEP
ColumnsTbl[["PEP"]] <- "PEP"
# - Filters
ColumnsTbl[["Filters"]] <- qualFlt
# - Clusters
if ((exists("KlustKols"))&&(length(KlustKols))) { ColumnsTbl[["Cluster"]] <- KlustKols }
# - Coverage
ColumnsTbl[["Coverage"]] <- covcol
# Melt
ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }, 1) > 0)]
ColumnsTbl <- listMelt(ColumnsTbl, names(ColumnsTbl), c("Col", "Grp"))
stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
#tst <- aggregate(ColumnsTbl$Col, list(ColumnsTbl$Col), length)
#tst[which(tst$x > 1),]
ColumnsTbl$Class <- ""
ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Protein Group information"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(pepcountcol1, pepcountcol2))] <- "Peptides count"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(pepidcol1, pepidcol2))] <- "Peptide IDs"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evcountcol1, evcountcol2))] <- "Evidences count"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evidcol1, evidcol2))] <- "Evidence IDs"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(speccountcol1, speccountcol2))] <- "Spectral count"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(specidcol1, specidcol2))] <- "Spectrum IDs"
if (IsBioID2) {
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(biotpepcountcol1, biotpepcountcol2))] <- "Biotin peptides count"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(biotpepidcol1, biotpepidcol2))] <- "Biotin peptide IDs"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(biotevcountcol1, biotevcountcol2))] <- "Biotin evidences count"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(biotevidcol1, biotevidcol2))] <- "Biotin evidence IDs"
}
for (nm in names(intRf)) { #nm <- names(intRf)[1]
  rpl <- intNms(nm, TRUE)
  kl <- c(paste0("Mean ", intRf[[nm]], VPAL$values), paste0(intRf[[nm]], RSA$values))
  kl <- kl[which(kl %in% ColumnsTbl$Col)]
  ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
}
ColumnsTbl$Class[grep("Ruler", ColumnsTbl$Grp)] <- "log10(est. copies/cell)"
for (nm in names(ratRf)) { #nm <- names(ratRf)[1]
  rpl <- ratNms(nm, TRUE)
  kl <- c(paste0("Mean ", ratRf[[nm]], VPAL$values), paste0(ratRf[[nm]], RSA$values))
  kl <- kl[which(kl %in% ColumnsTbl$Col)]
  ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
}
ColumnsTbl$Class[which(ColumnsTbl$Grp == "Proteome Ruler")] <- "log10(est. copies/cell)"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% pvalcol)] <- "P-value"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% signcol)] <- "Significant"
ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
if (F.test) {
  #ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test summary Ratios")] <- "F-test"
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test P-values")] <- "F-test"
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test significant")] <- "F-test"
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "F-test regulated")] <- "F-test"
}
if (Annotate) {
  ColumnsTbl$Class[grep("[Aa]nnotations", ColumnsTbl$Grp)] <- "Annotations"
  if (LocAnalysis) {
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "Localisation")] <- "Localisation"
    if (LocAnalysis2) {
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "SSDs P-values")] <- gsub(" - $", "", SSD.Pval.Root)
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "SSDs significant")] <- paste0(gsub(" - $", "", SSD.Pval.Root))
      ColumnsTbl$Class[which(ColumnsTbl$Grp == "Re-localized")] <- "Re-localized"
    }
  }
}
ColumnsTbl$Class[grep("[Ss]equence coverage \\[%\\]", ColumnsTbl$Col)] <- "Sequence coverage [%]"
ColumnsTbl$Class[grep("^1st ID cov\\.", ColumnsTbl$Col)] <- "1st accession sequence coverage (peptides)"
if ((exists("KlustKols"))&&(length(KlustKols))) {
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "Cluster")] <- paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ")")
}
ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters", "Negative filter"))] <- "QC filters"
stopifnot(min(nchar(ColumnsTbl$Class)) > 0)
#View(ColumnsTbl[which(nchar(ColumnsTbl$Class) == 0),])
w <- c(which(ColumnsTbl$Class == "General Protein Group information"),
       as.integer(unlist(lapply(names(intCols), function(nm) {
         rpl <- intNms(nm, TRUE)
         which(ColumnsTbl$Class == rpl)
       }))),
       which(ColumnsTbl$Class == "log10(est. copies/cell)"),
       as.integer(unlist(lapply(names(ratCols), function(nm) {
         rpl <- ratNms(nm, TRUE)
         which(ColumnsTbl$Class == rpl)
       }))),
       which(ColumnsTbl$Class == "P-value"),
       which(ColumnsTbl$Class == "Significant"),
       which(ColumnsTbl$Class == "Regulated"),
       #which(ColumnsTbl$Class == "F-test: Mean log2(Ratio)"),
       #which(ColumnsTbl$Class == "F-test: -log10(P-values)"),
       #which(ColumnsTbl$Class == "F-test: Significant"),
       #which(ColumnsTbl$Class == "F-test: Regulated"),
       which(ColumnsTbl$Class == "F-test"))
if (Annotate&&LocAnalysis) {
  w <- c(w,
         which(ColumnsTbl$Class == gsub(" - $", "", SSD.Pval.Root)),
         which(ColumnsTbl$Class == paste0(gsub(" - $", "", SSD.Pval.Root))),
         which(ColumnsTbl$Class == "Re-localized"))
}
w <- c(w,
       which(ColumnsTbl$Class == "QC filters"),
       which(ColumnsTbl$Class == "Sequence coverage [%]"),
       which(ColumnsTbl$Class == "1st accession sequence coverage (peptides)"),
       which(ColumnsTbl$Class == paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ")")),
       which(ColumnsTbl$Class == "Peptides count"),
       which(ColumnsTbl$Class == "Peptide IDs"),
       which(ColumnsTbl$Class == "Evidences count"),
       which(ColumnsTbl$Class == "Evidence IDs"),
       which(ColumnsTbl$Class == "Spectral count"),
       which(ColumnsTbl$Class == "Spectrum IDs"),
       which(ColumnsTbl$Class == "Biotin peptides count"),
       which(ColumnsTbl$Class == "Biotin peptide IDs"),
       which(ColumnsTbl$Class == "Biotin evidences count"),
       which(ColumnsTbl$Class == "Biotin evidence IDs"),
       which(ColumnsTbl$Class == "Annotations"))
stopifnot(length(w) == nrow(ColumnsTbl))
#View(ColumnsTbl[which(!1:nrow(ColumnsTbl) %in% w),])
ColumnsTbl <- ColumnsTbl[w,]
ColumnsTbl$Hide <- ColumnsTbl$Class %in% c("Peptide IDs", "Peptides count", "Evidence IDs", "Evidences count", "Spectral count", "Spectrum IDs",
                                           "Biotin peptides count", "Biotin peptide IDs", "Biotin evidences count", "Biotin evidence IDs",
                                           "Annotations",
                                           unique(grep("^F-test: ", ColumnsTbl$Class, value = TRUE)) # Let's not show too much stuff
)
if (length(intCols) > 1) {
  for (nm in names(intCols)[1:(length(intCols) - 1)]) {
    rpl <- intNms(nm)
    ColumnsTbl$Hide[which(ColumnsTbl$Class == rpl)] <- TRUE
  }
}
if (length(ratCols) > 1) {
  for (nm in names(ratCols)[1:(length(ratCols) - 1)]) {
    rpl <- ratNms(nm)
    ColumnsTbl$Hide[which(ColumnsTbl$Class == rpl)] <- TRUE
  }
}
#
if (MakeRatios) { a <- KolEdit(ColumnsTbl$Col, intColsTbl, ratColsTbl) } else { a <- KolEdit(ColumnsTbl$Col, intColsTbl) }
ColumnsTbl$edit_Col <- unlist(a)
#wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#saveFun(WorkBook, file = "WorkBook_bckp.RData")
#wb_save(WorkBook, paste0(wd, "/tst.xlsx")); xl_open(paste0(wd, "/tst.xlsx"))
#loadFun("WorkBook_bckp.RData")
if (saintExprs) {
  tblMode <- tblMode2 <- TbNm <- "SAINTexpress"
  # Function for editing the header
  KolEdit <- function(KolNames, #intTbl = intColsTbl,
                      ratTbl = ratColsTbl) {
    #KolNames <- xlTabs[[sheetnm]]
    klnms <- KolNames
    KolNames <- gsub("^log2\\(FC\\) - ", "log2(FC) ", KolNames)
    KolNames <- gsub("^AvgP - ", "AvgP ", KolNames)
    KolNames <- gsub("^MaxP - ", "MaxP ", KolNames)
    KolNames <- gsub("^TopoAvgP - ", "TopoAvgP ", KolNames)
    KolNames <- gsub("^TopoMaxP - ", "TopoMaxP ", KolNames)
    KolNames <- gsub("^SaintScore - ", "SaintScore ", KolNames)
    KolNames <- gsub("^OddsScore - ", "OddsScore ", KolNames)
    KolNames <- gsub("^boosted_by - ", "Boosted by... ", KolNames)
    KolNames <- gsub("^BFDR - ", "BFDR ", KolNames)
    KolNames <- cleanNms(KolNames, start = FALSE)
    #
    # Those names must be unique if the data is to be written as a table!
    # Which is annoying, because this limits how much fat we can cut
    tst <- aggregate(KolNames, list(KolNames), c)
    tst$L <- vapply(tst$x, length, 1)
    tst <- tst[which(tst$Group.1 != ""),]
    stopifnot(max(tst$L) == 1)
    #tst$x[which(tst$L > 1)]
    #
    KolNames <- as.data.frame(t(KolNames))
    colnames(KolNames) <- klnms
    return(KolNames) #View(KolNames)
  }
  #
  tempData <- allSAINTs
  CoreCol <- c("Protein", "Gene", "Common Name")
  quantCols <- quantcol <- "Av. log10 abundance"
  ratRf <- setNames("log2(FC) - ", "log2(rat.), avg.")
  ratColsTbl <- setNames(lapply(names(ratRf), function(nm) {
    res <- data.frame(Log = paste0(ratRf[nm], VPAL$values),
                      Type = "Average",
                      Sample = VPAL$values)
    w <- which(res$Log %in% colnames(tempData))
    return(res[w,])
  }), names(ratRf))
  w <- which(vapply(ratColsTbl, nrow, 1) > 0)
  ratColsTbl <- ratColsTbl[w]; ratRf <- ratRf[w]
  ratCols <- lapply(ratColsTbl, function(x) { x$Log })
  grl <- unlist(ratCols)
  for (gr in grl) {
    w <- which(is.infinite(tempData[[gr]]))
    tempData[w, gr] <- NA
  }
  quantcol <- c(quantcol, grl)
  quantCols[names(ratRf)] <- ratCols
  fdrcol <- grep("^BFDR - ", colnames(tempData), value = TRUE)
  avgPcol <- grep("^AvgP - ", colnames(tempData), value = TRUE)
  maxPcol <- grep("^MaxP - ", colnames(tempData), value = TRUE)
  topoAvgPcol <- grep("^TopoAvgP - ", colnames(tempData), value = TRUE)
  topoMaxPcol <- grep("^TopoMaxP - ", colnames(tempData), value = TRUE)
  saintcol <- grep("^SaintScore - ", colnames(tempData), value = TRUE)
  oddscol <- grep("^OddsScore - ", colnames(tempData), value = TRUE)
  boostcol <- grep("^boosted_by - ", colnames(tempData), value = TRUE)
  qualFlt <- "Potential contaminant"
  kol <- c(CoreCol, "In list", quantcol, fdrcol, avgPcol, maxPcol, topoAvgPcol, topoMaxPcol,
           #saintcol, oddscol, # Those 2 columns look a bit useless... 
           boostcol, qualFlt[which(qualFlt != "In list")])
  kol <- unique(kol)
  kol <- kol[which(kol %in% colnames(tempData))]
  tempData <- tempData[, kol]
  #
  # Re-order
  for (k in rev(fdrcol)) {
    tempData <- tempData[order(tempData[[k]], decreasing = FALSE),]
  }
  m <- match(prot.list, tempData$Protein)
  m <- m[which(!is.na(m))]
  w <- c(m, which(!tempData$Protein %in% prot.list))
  tempData <- tempData[w,]
  #
  # Which columns are affected by each style
  # - IDs
  ColumnsTbl <- list(IDs = CoreCol)
  # Quantitation
  # - Ratios
  for (nm in names(ratRf)) { #nm <- names(ratRf[1])
    ColumnsTbl[[names(ratRf)]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Average")]
  }
  # SAINTexpress stats
  ColumnsTbl$"P-values" <- fdrcol
  ColumnsTbl$AvgP <- avgPcol
  ColumnsTbl$MaxP <- maxPcol
  ColumnsTbl$topoAvgP <- topoAvgPcol
  ColumnsTbl$topoMaxP <- topoMaxPcol
  ColumnsTbl$SaintScore <- saintcol
  ColumnsTbl$OddsScore <- oddscol
  ColumnsTbl$boosted_by <- boostcol
  # - Filters
  ColumnsTbl$Filters <- qualFlt
  # Melt
  ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, function(x) { length(x[which(!is.na(x))]) }, 1) > 0)]
  ColumnsTbl <- listMelt(ColumnsTbl, names(ColumnsTbl), c("Col", "Grp"))
  stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
  #tst <- aggregate(ColumnsTbl$Col, list(ColumnsTbl$Col), length)
  #tst[which(tst$x > 1),]
  ColumnsTbl$Class <- ""
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Protein Group information"
  for (nm in names(ratRf)) { #nm <- names(ratRf)[1]
    rpl <- "log2(rat.), avg."
    kl <- paste0(ratRf[[nm]], VPAL$values)
    kl <- kl[which(kl %in% ColumnsTbl$Col)]
    ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
  }
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% fdrcol)] <- "P-values"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(avgPcol, maxPcol, topoAvgPcol, topoMaxPcol))] <- "SAINTexpress P"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% saintcol)] <- "SAINTexpress score"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% oddscol)] <- "SAINTexpress odds score"
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% boostcol)] <- "SAINTexpress boost"
  ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters", "Negative filter"))] <- "QC filters"
  stopifnot(min(nchar(ColumnsTbl$Class)) > 0)
  #View(ColumnsTbl[which(nchar(ColumnsTbl$Class) == 0),])
  w <- c(which(ColumnsTbl$Class == "General Protein Group information"),
         as.integer(unlist(lapply(names(ratCols), function(nm) {
           which(ColumnsTbl$Class == "log2(rat.), avg.")
         }))),
         which(ColumnsTbl$Grp == "P-values"),
         which(ColumnsTbl$Grp == "AvgP"),
         which(ColumnsTbl$Grp == "MaxP"),
         which(ColumnsTbl$Grp == "topoAvgP"),
         which(ColumnsTbl$Grp == "topoMaxP"),
         which(ColumnsTbl$Grp == "SaintScore"),
         which(ColumnsTbl$Grp == "OddsScore"),
         which(ColumnsTbl$Grp == "boosted_by"),
         which(ColumnsTbl$Class == "QC filters"))
  stopifnot(length(w) == nrow(ColumnsTbl))
  #View(ColumnsTbl[which(!1:nrow(ColumnsTbl) %in% w),])
  ColumnsTbl <- ColumnsTbl[w,]
  ColumnsTbl$Hide <- ColumnsTbl$Class %in% "SAINTexpress P" # Let's not show too much stuff
  #
  a <- KolEdit(ColumnsTbl$Col, ratColsTbl)
  ColumnsTbl$edit_Col <- unlist(a)
  Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
  #saveFun(WorkBook, file = "WorkBook_bckp.RData")
  #wb_save(WorkBook, paste0(wd, "/tst.xlsx")); xl_open(paste0(wd, "/tst.xlsx"))
  #loadFun("WorkBook_bckp.RData")
}
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/Write_Excel_end_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#WorkBook$get_active_sheet()
#xl_open(repFl)
