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
intNms <- \(nms, topLvl = FALSE, type = "PG") {
  m <- match(type, c("pep", "PG"))
  root <- c("Intensity", "Expression")[m]
  mode <- topLvl+1L
  sapply(nms, \(nm) {
    if (nm %in% c("Original", "Intensity", "Expression", "Original int.", "Intensity int.", "Expression int.")) {
      nm <- c(c("int.", "expr.")[m], root)[mode]
    } else {
      nm <- if (tolower(gsub("-", "", nm)) %in% c("renorm.", "renorm. int.", "renormalized", "renormalized int.")) {
        "re-norm"
      } else { substr(nm, 1L, min(c(3L, nchar(nm)))) }
      nm <- paste0(nm, ". ", c(c("int.", "expr.")[m], root)[mode])
    }
    paste0("log10(", nm, ")")
  })
}
ratNms <- \(nms, topLvl = FALSE) {
  mode <- topLvl+1L
  sapply(nms, \(nm) {
    if (nm %in% c("Original", "Ratios", "Original rat.", "Ratios rat.")) {
      nm <- c("rat.", "Ratio")[mode]
    } else {
      nm <- if (tolower(gsub("-", "", nm)) %in% c("renorm.", "renorm. rat.", "renormalized", "renormalized rat.")) {
        "re-norm"
      } else { substr(nm, 1L, min(c(3L, nchar(nm)))) }
      nm <- paste0(nm, ". ", c("rat.", "ratios")[mode])
    }
    paste0("log2(", nm, ")")
  })
}
for (nm in names(int.cols)) { #nm <- names(int.cols)[1L]
  rpl <- intNms(nm, type = "pep")
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
for (nm in names(PG.int.cols)) { #nm <- names(PG.int.cols)[1L]
  rpl <- intNms(nm)
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
if (MakeRatios) {
  for (nm in unique(c(names(rat.cols), names(PG.rat.cols)))) { #nm <- unique(c(names(rat.cols), names(PG.rat.cols)))[1L]
    rpl <- ratNms(nm)
    Styles[[paste0(rpl, ", avg.")]] <- "Summary Ratios"
    Styles[[paste0(rpl, ", indiv.")]] <- "Individual Ratios"
  }
  if (exists("PTMs_intRf")) {
    for (nm in names(PTMs_intRf)) { #nm <- names(PTMs_intRf)[1L]
      rpl <- intNms(nm, type = "pep")
      Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
      Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
    }
  }
  if (exists("PTMs_ratRf")) {
    for (nm in names(PTMs_ratRf)) { #nm <- names(PTMs_ratRf)[1L]
      rpl <- ratNms(nm)
      Styles[[paste0(rpl, ", avg.")]] <- "Summary Ratios"
      Styles[[paste0(rpl, ", indiv.")]] <- "Individual Ratios"
    }
  }
}
fl <- system.file("extdata", "Report - column names - no replicates.xlsx", package = "proteoCraft")
styleNms <- openxlsx2::read_xlsx(fl, "tmp", colNames = FALSE)[, 1L]
# wb <- loadWorkbook(fl)
# addWorksheet(wb, "tmp")
# tmpFl <- temp_xlsx()
#w <- which(vapply(Styles, \(x) { "Style" %in% class(x) }, TRUE))
# styleNms %<o% names(Styles)[w]
# writeData(wb, "tmp", styleNms)
# for (i in seq_along(w)) { addStyle(wb, "tmp", Styles[[w[i]]], i, 1L) }
# saveWorkbook(wb, tmpFl)
# openXL(tmpFl)
# WorkBook %<o% wb_load(tmpFl)
# Styles2 %<o% setNames(w, names(Styles)[w])
WorkBook %<o% wb_load(fl)
repFl <- paste0(wd, "/Tables/Report_", dtstNm, ".xlsx")
WorkBook <- wb_add_data(WorkBook, "Description", dtstNm, wb_dims(2L, 5L))
WorkBook <- wb_add_data(WorkBook, "Description", format(Sys.Date(), "%d/%m/%Y"), wb_dims(3L, 5L))
WorkBook <- wb_add_data(WorkBook, "Description", WhoAmI, wb_dims(4L, 5L))
tmp <- loadedPackages(TRUE)
WorkBook <- wb_add_data(WorkBook, "Description", tmp$Version[grep("proteoCraft", tmp$Name)], wb_dims(5L, 5L))
cat(" - Writing Excel report...\n")
# Function for editing the header
KolEdit <- \(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl) {
  #KolNames <- xlTabs[[sheetnm]]
  klnms <- KolNames
  KolNames[which(KolNames == "Evidence IDs")] <- paste0("All ", evNm, " IDs")
  for (nm in names(intTbl)) { #nm <- names(intTbl)[1L] #nm <- names(intTbl)[2L]
    m <- match(intTbl[[nm]]$Log, KolNames)
    w <- which(!is.na(m))
    if (length(w)) {
      rpl <- intNms(nm, type = "pep")
      KolNames[m[w]] <- paste0(rpl, " ", intTbl[[nm]]$Sample[w])
    }
  }
  if (MakeRatios) {
    for (nm in names(ratTbl)) {
      m <- match(ratTbl[[nm]]$Log, KolNames)
      w <- which(!is.na(m))
      if (length(w)) {
        rpl <- ratNms(nm)
        KolNames[m[w]] <- paste0(rpl, " ", ratTbl[[nm]]$Sample[w])
      }
    }
    KolNames <- gsub("^Enriched", "Enr.", KolNames)
    KolNames <- gsub("^Regulated", "Reg.", KolNames)
  }
  KntKol <- paste0(AA, " Count")
  KolNames[which(KolNames %in% KntKol)] <- gsub(" Count$", "", KolNames[which(KolNames %in% KntKol)])
  #
  KolNames <- gsub("( - )|(___)", " ", KolNames)
  #
  # Those names must be unique if the data is to be written as a table!
  tst <- aggregate(KolNames, list(KolNames), length)
  stopifnot(max(tst$x) == 1L)
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
  ev$"In list" <- c("", "+")[w+1L]
}
QualFilt %<o% c(#pgOrgKol,
  "Potential contaminant", "Only identified by site",
  grep("^Quality filter: ", colnames(PG), value = TRUE))
if (NegFilt) { QualFilt <- c(QualFilt, "Direct identification in negative filter sample(s)") }
II <- setNames(1L, "All peptidoforms")
if ((length(Mod2Write))&&(PTMriched)) {
  II[paste0(Modifs$`Full name`[match(Mod2Write, Modifs$Mark)], "-mod. pept.")] <- 1L+(seq_along(length(Mod2Write)))
}
for (ii in II) { #ii <- II[1L] #ii <- II[2L]
  tblMode <- tblMode2 <- "pep"
  TbNm <- names(II)[ii]
  tempData <- get(tblMode)
  intRf <- int.cols
  ratRf <- rat.cols
  if (ii > 1L) {
    ptm <- Mod2Write[[ii-1L]]
    Ptm <- Modifs$`Full name`[match(ptm, Modifs$Mark)]
    tempData <- PTMs_pep[[Ptm]]
    intRf <- PTMs_intRf
    ratRf <- PTMs_ratRf
    tblMode2 <- paste0(Ptm, "-modified pep")
  }
  # names(intRf) <- gsub("^Original ", "", paste0(names(intRf), " intensities"))
  # names(ratRf) <- gsub("^Original ", "", paste0(names(ratRf), " ratios"))
  # (In this workflow, as of this version all peptide intensities are non-log-transformed but ratios are!)
  if (nrow(tempData)) {
    CoreCol <- "Modified sequence"
    if ("Modified sequence_verbose" %in% colnames(tempData)) { CoreCol <- c(CoreCol, "Modified sequence_verbose") }
    CoreCol <- c(CoreCol, "Sequence", "id", "Proteins")
    CoreCol2 <- c("Leading proteins", "Leading razor proteins",
                  "Protein names", "Gene names", "Protein group IDs", "Razor protein group ID",
                  grep(" (Probabilities|Score Diffs)$", colnames(tempData), value = TRUE),
                  "Normalisation group")
    evcol <- "Evidence IDs"
    spcol <- "MS/MS count"
    intColsTbl <- setNames(lapply(names(intRf), \(nm) {
      res <- data.frame(nonLog = c(intRf[nm], paste0(intRf[nm], " - ", Exp)),
                        Log = c(paste0("log10(", intRf[nm], ")"), paste0("log10(", intRf[nm], ") - ", Exp)),
                        Type = c("Average", rep("Individual", length(Exp))),
                        Sample = c("Average", Exp))
      w <- which(res$nonLog %in% colnames(tempData))
      return(res[w,])
    }), names(intRf))
    quantCols <- intCols <- lapply(intColsTbl, \(x) { x$Log })
    gel0 <- unlist(lapply(intColsTbl, \(x) { x$nonLog }))
    gel <- unlist(intCols)
    tempData[, gel] <- log10(tempData[, gel0]) # Peptide expression values are not log-transformed
    tempData[, gel0] <- NULL
    for (gl in gel) {
      w <- which(is.infinite(tempData[[gl]]))
      tempData[w, gl] <- NA
    }
    quantcol <- gel
    for (nm in names(intRf)) { 
      if (!grepl("log[0-9]+", intRf[nm])) { # In case I am rerunning a bit
        intRf[nm] <- paste0("log10(", intRf[nm], ")")
      }
    }
    if (MakeRatios) {
      ratColsTbl <- setNames(lapply(names(ratRf), \(nm) {
        # Note: Peptide values are already log-transformed in this workflow!!!
        res <- data.frame(
          Log = c(paste0("Mean ", ratRf[nm]), paste0(ratRf[nm], " - ", Exp)),
          Type = c("Average", rep("Individual", length(Exp))),
          Sample = c("Average", Exp))
        w <- which(res$Log %in% colnames(tempData))
        return(res[w,])
      }), names(ratRf))
      ratCols <- lapply(ratColsTbl, \(x) { x$Log })
      grl <- unlist(ratCols)
      for (gr in grl) {
        w <- which(is.infinite(tempData[[gr]]))
        tempData[w, gr] <- NA
      }
      quantcol <- c(quantcol, grl)
      quantCols[names(ratRf)] <- ratCols
      regcol <- grep("^((Enriched)|(Regulated)) - ", colnames(tempData), value = TRUE)
    }
    aacol <- paste0(AA, " Count")
    qualFlt <- QualFilt[which(QualFilt %in% colnames(ev))]
    w <- which(!qualFlt %in% colnames(tempData))
    if (length(w)) {
      tempData[, qualFlt[w]] <- ev[match(tempData$"Modified sequence", ev$"Modified sequence"), qualFlt[w]]
    }
    dir <- paste0(wd, "/Tables")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    data.table::fwrite(tempData, paste0(dir, "/", TbNm, ".tsv"), sep = "\t", row.names = FALSE, na = "NA")
    w <- grsep2(prot.list, tempData$Proteins)
    if (length(w)) {
      data.table::fwrite(tempData[w,], paste0(wd, "/Tables/", TbNm, " - Proteins in list.tsv"),
                         sep = "\t", row.names = FALSE, na = "NA")
    }
    # Columns table
    # - IDs
    ColumnsTbl <- list(IDs = c(CoreCol, CoreCol2))
    # - Counts
    ColumnsTbl[["AA counts"]] <- aacol
    # - Evidence counts and IDs
    ColumnsTbl[["Global Ev. IDs"]] <- "Evidence IDs"
    ColumnsTbl[["Global Spec. counts"]] <- "MS/MS count"
    # - Expression values
    for (nm in names(intRf)) { #nm <- names(intRf[1L])
      rpl <- intNms(nm, type = "pep")
      ColumnsTbl[[paste0(rpl, ", avg.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Average")]
      ColumnsTbl[[paste0(rpl, ", indiv.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Individual")]
    }
    if (MakeRatios) {
      for (nm in names(ratRf)) { #nm <- names(ratRf[1L])
        rpl <- ratNms(nm)
        ColumnsTbl[[paste0(rpl, ", avg.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Average")]
        ColumnsTbl[[paste0(rpl, ", indiv.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Individual")]
      }
    }
    # - PEP
    ColumnsTbl[["PEP"]] <- "PEP"
    # - Filters
    ColumnsTbl[["Filters"]] <- qualFlt
    # Melt
    ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, \(x) { sum(!is.na(x)) }, 1L) > 0L)]
    ColumnsTbl <- listMelt(ColumnsTbl, names(ColumnsTbl), c("Col", "Grp"))
    #aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), length)
    stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
    ColumnsTbl$Class <- ""
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Peptides information"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evcol))] <- "Evidence IDs"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(spcol))] <- "Spectral count"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% aacol)] <- "Amino Acid counts"
    ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters"))] <- "QC filters"
    for (nm in names(intRf)) { #nm <- names(intRf)[1L]
      rpl <- intNms(nm, TRUE, type = "pep")
      ColumnsTbl$Class[grep(topattern(paste0(intRf[nm], " - ")), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[which(ColumnsTbl$Col == intRf[nm])] <- rpl
    }
    if (MakeRatios) {
      for (nm in names(ratRf)) { #nm <- names(ratRf)[1L]
        rpl <- ratNms(nm, TRUE)
        ColumnsTbl$Class[grep(topattern(paste0(ratRf[nm], " - ")), ColumnsTbl$Col)] <- rpl
        ColumnsTbl$Class[which(ColumnsTbl$Col == ratRf[nm])] <- rpl
      }
      ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
    }
    ColumnsTbl$Class[grep("[Aa]nnotations", ColumnsTbl$Grp)] <- "Annotations"
    stopifnot(min(nchar(ColumnsTbl$Class)) > 0L) #View(ColumnsTbl)
    w <- c(which(ColumnsTbl$Class == "General Peptides information"),
           unlist(lapply(names(intCols), \(nm) { which(ColumnsTbl$Class == intNms(nm, TRUE, type = "pep")) }))
    )
    if (MakeRatios) {
      w <- c(w,
             unlist(lapply(names(ratCols), \(nm) { which(ColumnsTbl$Class == ratNms(nm, TRUE)) })),
             which(ColumnsTbl$Class == "Regulated")
      )
    }
    w <- c(w, which(ColumnsTbl$Class == "QC filters"),
           which(ColumnsTbl$Class == "Evidence IDs"),
           which(ColumnsTbl$Class == "Spectral count"),
           which(ColumnsTbl$Class == "Amino Acid counts"))
    stopifnot(length(w) == nrow(ColumnsTbl))
    ColumnsTbl <- ColumnsTbl[w,]
    ColumnsTbl$Hide <- ColumnsTbl$Class %in% c("Spectral count", "Spectrum IDs",
                                               "Evidences count", "Evidence IDs",
                                               "Amino Acid counts", "Annotations")
    l <- length(intRf)
    if (l > 1L) {
      for (nm in names(intRf)[1L:(l-1L)]) {
        ColumnsTbl$Hide[which(ColumnsTbl$Class == intNms(nm, TRUE, type = "pep"))] <- TRUE
      }
    }
    if (MakeRatios) {
      l <- length(ratRf)
      if (l > 1L) {
        for (nm in names(ratRf)[1L:(l-1L)]) {
          ColumnsTbl$Hide[which(ColumnsTbl$Class == ratNms(nm, TRUE))] <- TRUE
        }
      }
    }
    #
    a <- if (MakeRatios) { KolEdit(ColumnsTbl$Col, intColsTbl, ratColsTbl) } else { KolEdit(ColumnsTbl$Col, intColsTbl) }
    ColumnsTbl$edit_Col <- unlist(a)
    #
    Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
}
#saveFun(WorkBook, file = "WorkBook_bckp.RData")
#wb_save(WorkBook, repFl);xl_open(repFl)
#loadFun("WorkBook_bckp.RData")
TbNm <- "Protein groups"
tblMode <- tblMode2 <- "PG"
# Function for editing the header
KolEdit <- \(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl) {
  #KolNames <- xlTabs[[sheetnm]]
  klnms <- KolNames
  KolNames <- gsub("Peptides?", "Pep.", KolNames)
  KolNames <- gsub("Evidences?", "PSMs", KolNames)
  KolNames <- gsub("Spectr((al)|(um))", "Spec.", KolNames)
  KolNames <- gsub("Razor", "Raz.", KolNames)
  KolNames <- gsub("Unique", "Uniq.", KolNames)
  KolNames <- gsub("MS/MS", "MS2", KolNames)
  for (nm in names(intTbl)) { #nm <- names(intTbl)[1L] #nm <- names(intTbl)[2L]
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
    KolNames <- gsub(".*Regulated - ", "reg. ", KolNames)
  }
  #KolNames <- gsub(".*Pvalue\\)( - )?", "-log10 pval. ", KolNames)
  #KolNames <- gsub(".*Significant-", "signif. ", KolNames)
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
  # Those names must be unique if the data is to be written as a table!
  # Which is annoying, because this limits how much fat we can cut
  tst <- aggregate(KolNames, list(KolNames), c)
  tst$L <- lengths(tst$x)
  tst <- tst[which(tst$Group.1 != ""),]
  stopifnot(max(tst$L) == 1L)
  #tst$x[which(tst$L > 1L)]
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
                  "Evidences count",
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
intRf <- PG.int.cols
names(intRf) <- paste0(names(intRf), " int.")
intColsTbl <- setNames(lapply(names(intRf), \(nm) { #nm <- names(intRf)[1L]
  res <- data.frame(
    Log = c(paste0("Mean ", gsub(" - ", "", intRf[nm])), paste0(intRf[nm], Exp)),
    Type = c("Average", rep("Individual", length(Exp))),
    Sample = c("Average", Exp))
  w <- which(res$Log %in% colnames(tempData))
  return(res[w,])
}), names(intRf))
quantCols <- intCols <- lapply(intColsTbl, \(x) { x$Log })
gel <- unlist(intCols)
for (gl in gel) {
  w <- which(is.infinite(tempData[[gl]]))
  tempData[w, gl] <- NA
}
quantcol <- gel
if (MakeRatios) {
  ratRf <- PG.rat.cols
  names(ratRf) <- paste0(names(ratRf), " rat.")
  ratColsTbl <- setNames(lapply(names(ratRf), \(nm) {
    # Note: Peptide values are already log-transformed in this workflow!!!
    res <- data.frame(
      Log = c(paste0("Mean ", gsub(" - $", "", ratRf[nm])), paste0(ratRf[nm], Exp)),
      Type = c("Average", rep("Individual", length(Exp))),
      Sample = c("Average", Exp))
    w <- which(res$Log %in% colnames(tempData))
    return(res[w,])
  }), names(ratRf))
  ratCols <- lapply(ratColsTbl, \(x) { x$Log })
  grl <- unlist(ratCols)
  for (gr in grl) {
    w <- which(is.infinite(tempData[[gr]]))
    tempData[w, gr] <- NA
  }
  quantcol <- c(quantcol, grl)
  quantCols[names(ratRf)] <- ratCols
  regcol <- grep("^((Enriched)|(Regulated)) - ", colnames(tempData), value = TRUE)
}
if (protrul) { quantcol <- c(quantcol, grep(topattern("log10(est. copies/cell) - ", start = FALSE), colnames(tempData), value = TRUE)) }
covcol <- c(#xmlCovCol,
            c("Sequence coverage [%]",
              "Unique + razor sequence coverage [%]",
              "Unique sequence coverage [%]")[1L:c(1L, 3L)[isEukaLike+1L]],
            grep(topattern("Sequence coverage [%] - "), colnames(tempData), value = TRUE))
if (WorkFlow == "Band ID") {
  covcol <- c("Max. theoretical sequence coverage [%]", covcol)
}
kol <- c(kol, "Mol. weight [kDa]", covcol, "PEP", quantcol)
if ((exists("KlustKols"))&&(length(KlustKols))) { kol <- c(kol, KlustKols) }
qualFlt <- QualFilt
kol <- unique(c(kol, qualFlt))
if (Annotate) { kol <- c(kol, annot.col) }
kol <- unique(kol[which(kol %in% colnames(tempData))])
tempData <- tempData[, kol]
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
for (nm in names(intRf)) { #nm <- names(intRf[1L])
  rpl <- intNms(nm)
  ColumnsTbl[[paste0(rpl, ", avg.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Average")]
  ColumnsTbl[[paste0(rpl, ", indiv.")]] <- intColsTbl[[nm]]$Log[which(intColsTbl[[nm]]$Type == "Individual")]
}
# - Ratios
if (MakeRatios) {
  for (nm in names(ratRf)) { #nm <- names(ratRf[1L])
    rpl <- ratNms(nm)
    ColumnsTbl[[paste0(rpl, ", avg.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Average")]
    ColumnsTbl[[paste0(rpl, ", indiv.")]] <- ratColsTbl[[nm]]$Log[which(ratColsTbl[[nm]]$Type == "Individual")]
  }
  ColumnsTbl[["Regulated"]] <- regcol
}
# - Proteome Ruler
if (protrul) {
  ColumnsTbl[["Proteome Ruler"]] <- grep(topattern("log10(est. copies/cell) - ", start = FALSE),
                                         colnames(tempData), value = TRUE)
}
# - Annotations
if (Annotate) {
  annot.col2 <- gsub("_names$", " names", annot.col)
  AnnotTbl$Columns <- list(c("GO", "GO-ID"), c("Taxonomy", "TaxID"), NA, NA, NA, NA, "EMBL", NA)
  for (i in annot) { AnnotTbl$Columns[match(i, AnnotTbl$Name)] <- list(c(i, paste0(i, " names"))) }
  AnnotTbl$Columns[match("Other", AnnotTbl$Name)] <- list(annot.col2[which(!annot.col2 %in% unlist(AnnotTbl$Columns))])
  for (i in 1L:nrow(AnnotTbl)) { ColumnsTbl[[paste0(AnnotTbl$Name[i], " annotations")]] <- AnnotTbl$Columns[[i]] }
}
# - PEP
ColumnsTbl[["PEP"]] <- "PEP"
# - Filters
ColumnsTbl[["Filters"]] <- qualFlt
ColumnsTbl[["In list"]] <- "In list"
# - Clusters
if ((exists("KlustKols"))&&(length(KlustKols))) { ColumnsTbl[["Cluster"]] <- KlustKols }
# - Coverage
ColumnsTbl[["Coverage"]] <- covcol
# Melt
ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, \(x) { sum(!is.na(x)) }, 1L) > 0L)]
ColumnsTbl <- listMelt(ColumnsTbl, names(ColumnsTbl), c("Col", "Grp"))
stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
#tst <- aggregate(ColumnsTbl$Col, list(ColumnsTbl$Col), length)
#tst[which(tst$x > 1L),]
ColumnsTbl$Class <- ""
ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Protein Group information"
ColumnsTbl$Class[which(ColumnsTbl$Grp == "In list")] <- "In list"
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
for (nm in names(intRf)) { #nm <- names(intRf)[1L] #nm <- names(intRf)[2L]
  rpl <- intNms(nm, TRUE)
  kl <- c(paste0("Mean ", gsub(" - $", "", intRf[nm])), paste0(intRf[[nm]], Exp))
  kl <- kl[which(kl %in% ColumnsTbl$Col)]
  ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
}
ColumnsTbl$Class[grep("Ruler", ColumnsTbl$Grp)] <- "log10(est. copies/cell)"
if (MakeRatios) {
  for (nm in names(ratRf)) { #nm <- names(ratRf)[1L]
    rpl <- ratNms(nm, TRUE)
    kl <- c(paste0("Mean ", gsub(" - $", "", ratRf[[nm]])), paste0(ratRf[[nm]], Exp))
    kl <- kl[which(kl %in% ColumnsTbl$Col)]
    ColumnsTbl$Class[match(kl, ColumnsTbl$Col)] <- rpl
  }
  ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
}
ColumnsTbl$Class[grep("[Aa]nnotations", ColumnsTbl$Grp)] <- "Annotations"
ColumnsTbl$Class[grep("[Ss]equence coverage \\[%\\]", ColumnsTbl$Col)] <- "Sequence coverage [%]"
ColumnsTbl$Class[grep("^1st ID cov\\.", ColumnsTbl$Col)] <- "1st accession sequence coverage (peptides)"
if ((exists("KlustKols"))&&(length(KlustKols))) {
  ColumnsTbl$Class[which(ColumnsTbl$Grp == "Cluster")] <- paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ")")
}
ColumnsTbl$Class[which(ColumnsTbl$Grp %in% c("PEP", "Filters", "Negative filter"))] <- "QC filters"
stopifnot(min(nchar(ColumnsTbl$Class)) > 0L)
w <- c(which(ColumnsTbl$Class == "General Protein Group information"),
       unlist(lapply(names(intCols), \(nm) { which(ColumnsTbl$Class == intNms(nm, TRUE)) })),
       which(ColumnsTbl$Class == "log10(est. copies/cell)")
)
if (MakeRatios) {
  w <- c(w,
         unlist(lapply(names(ratCols), \(nm) { which(ColumnsTbl$Class == ratNms(nm, TRUE)) })),
         which(ColumnsTbl$Class == "Regulated")
  )
}
w <- c(w,
       which(ColumnsTbl$Class == "QC filters"),
       which(ColumnsTbl$Class == "In list"),
       which(ColumnsTbl$Class == "Sequence coverage [%]"),
       which(ColumnsTbl$Class == "1st accession sequence coverage (peptides)"),
       which(ColumnsTbl$Class == paste0("Cluster (", c("K-means", "hierarch.")[KlustMeth], ")")),
       which(ColumnsTbl$Class == "Peptides count"),
       which(ColumnsTbl$Class == "Peptide IDs"),
       which(ColumnsTbl$Class == "Evidences count"),
       which(ColumnsTbl$Class == "Evidence IDs"),
       which(ColumnsTbl$Class == "Spectral count"),
       which(ColumnsTbl$Class == "Biotin peptides count"),
       which(ColumnsTbl$Class == "Biotin peptide IDs"),
       which(ColumnsTbl$Class == "Biotin evidences count"),
       which(ColumnsTbl$Class == "Biotin evidence IDs"),
       which(ColumnsTbl$Class == "Annotations"))
stopifnot(length(w) == nrow(ColumnsTbl))
ColumnsTbl <- ColumnsTbl[w,]
ColumnsTbl$Hide <- ColumnsTbl$Class %in% c("Peptide IDs", "Peptides count", "Evidence IDs", "Evidences count", "Spectral count", "Spectrum IDs",
                                           "Biotin peptides count", "Biotin peptide IDs", "Biotin evidences count", "Biotin evidence IDs",
                                           "Annotations")
if (length(intCols) > 1L) {
  for (nm in names(intCols)[1L:(length(intCols) - 1L)]) {
    ColumnsTbl$Hide[which(ColumnsTbl$Class == intNms(nm, TRUE))] <- TRUE
  }
}
if (MakeRatios) {
  if (length(ratCols) > 1L) {
    for (nm in names(ratCols)[1L:(length(ratCols) - 1L)]) {
      ColumnsTbl$Hide[which(ColumnsTbl$Class == ratNms(nm, TRUE))] <- TRUE
    }
  }
}
#
a <- if (MakeRatios) { KolEdit(ColumnsTbl$Col, intColsTbl, ratColsTbl) } else { KolEdit(ColumnsTbl$Col, intColsTbl) }
ColumnsTbl$edit_Col <- unlist(a)
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/fstWrite_Excel_core_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
Src <- paste0(libPath, "/extdata/R scripts/Sources/Write_Excel_end_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#xl_open(repFl)
