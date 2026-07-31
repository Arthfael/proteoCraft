##############################
# Create output Excel tables #
##############################
#
# Version for peptides-only script with replicates
#
# Write the main peptidoforms- and protein groups-level, multi-tabs report.
#
## Main peptidoforms-level, multi-tabs report
# Create openxlsx2 styles
#   It may make more sense from the way Excel works, but I HATE how openxlsx2 deals with styles!
#   This is not a statement on openxlsx2 or Jan Marvin, which are both awesome!
#   But the one thing I preferred in openxlsx over openxlsx2 was how it dealt with style.
#   Anyway... 
#   So. We will. CHEAT!
#   I have saved a dummy tab with my old openxlsx styles,
#   which I will load in openxlsx2 to get and copy the styles from.
MakeRatios <- TRUE # (Used by the sourced, core sub-script)
intNms <- \(nms, topLvl = FALSE, type = "PG", newLine = FALSE) {
  m <- match(type, c("pep", "PG"))
  root <- c("Intensity", "Expression")[m]
  mode <- topLvl+1L
  vapply(nms, \(nm) {
    if (nm %in% c("Original", "Intensity", "Expression", "Original int.", "Intensity int.", "Expression int.")) {
      nm <- c(c("int.", "expr.")[m], root)[mode]
    } else {
      nmTst <- tolower(gsub("-|\\.", "", nm))
      nm <- if (grepl("((re)|(back))norm", nmTst)) {
        "re-norm" 
      } else { tolower(substr(nm, 1L, min(c(3L, nchar(nm))))) }
      nm <- paste0(nm, ". ", c(c("int.", "expr.")[m], root)[mode])
    }
    nm <- paste0("log10(", nm, ")")
    if (newLine) { nm <- paste0(nm, " ///NL///") }
    return(nm)
  }, "")
}
ratNms <- \(nms, topLvl = FALSE, newLine = FALSE) {
  mode <- topLvl+1L
  vapply(nms, \(nm) {
    if (nm %in% c("Original", "Ratios", "Original rat.", "Ratios rat.")) {
      nm <- c("rat.", "Ratio")[mode]
    } else {
      nmTst <- tolower(gsub("-|\\.", "", nm))
      nm <- if (grepl("((re)|(back))norm", nmTst)) {
        "re-norm"
      } else { tolower(substr(nm, 1L, min(c(3L, nchar(nm))))) }
      nm <- paste0(nm, ". ", c("rat.", "ratios")[mode])
    }
    paste0("log2(", nm, ")")
    if (newLine) { nm <- paste0(nm, " ///NL///") }
    return(nm)
  }, "")
}
for (nm in names(pep.ref)) { #nm <- names(pep.ref[1L])
  rpl <- intNms(nm, type = "pep")
  Styles[[paste0(rpl, ", avg.")]] <- "Summary Expr"
  Styles[[paste0(rpl, ", indiv.")]] <- "Individual Expr"
}
fl <- system.file("extdata", "Report - column names - with replicates.xlsx", package = "proteoCraft")
styleNms <- openxlsx2::read_xlsx(fl, "tmp", colNames = FALSE)[, 1L]
WorkBook %<o% wb_load(fl)
repFl <- paste0(wd, "/Tables/Report_", dtstNm, ".xlsx")
WorkBook <- wb_add_data(WorkBook, "Description", dtstNm, wb_dims(2L, 5L))
WorkBook <- wb_add_data(WorkBook, "Description", format(Sys.Date(), "%d/%m/%Y"), wb_dims(3L, 5L))
WorkBook <- wb_add_data(WorkBook, "Description", WhoAmI, wb_dims(4L, 5L))
tmp <- loadedPackages(TRUE)
WorkBook <- wb_add_data(WorkBook, "Description", tmp$Version[grep("proteoCraft", tmp$Name)], wb_dims(5L, 5L))
WorkBook <- wb_set_base_font(WorkBook, 11L, font_name = "Calibri")
cat(" - Writing Excel report...\n")
#
# Function for editing our header
replacements_Map <- data.frame(in_PG = c(myContrasts$Contrast, RSA$values, VPAL$values),
                               Name = c(sub(" - ", " ///VS/// ", myContrasts$Contrast), cleanNms(c(RSA$values, VPAL$values))))
replacements_Map <- cbind(replacements_Map,
                          data.frame(in_PG = replacements_Map$Name,
                                     Namw = replacements_Map$Name))
replacements_Map$nchar <- nchar(replacements_Map$in_PG)
replacements_Map <- replacements_Map[order(replacements_Map$nchar, decreasing = TRUE),]
replFun <- \(colNames,
             replMap = replacements_Map) {
  replColNames <- colNames
  nc <- nchar(colNames)
  replMap$matches <- lapply(1L:nrow(replMap), \(x) {
    which(substr(colNames, nc-replMap$nchar[[x]]+1L, nc) == replMap$in_PG[[x]])
  })
  replMap <- replMap[which(lengths(replMap$matches) > 0L),]
  if (!nrow(replMap)) { return(replColNames) }
  col2val <- listMelt(replMap$matches, replMap$in_PG, c("match", "val"))
  col2val$nchar <- nchar(col2val$match)
  col2val <- aggregate(1L:nrow(col2val), list(col2val$match), \(x) {
    v <- col2val$val[[x]]
    l <- col2val$nchar[[x]]
    return(v[which(l == max(l))])
  })
  colnames(col2val) <- c("column", "match")
  m <- match(col2val$match, replMap$in_PG)
  col2val$repl <- replMap$Name[m]
  col2val$nchar <- replMap$nchar[m]
  #
  w <- col2val$column
  replColNames[w] <- paste0(substr(replColNames[w], 1L, nc[w] - col2val$nchar),
                            col2val$repl)
  return(replColNames)
}
KolEdit <- \(KolNames, intTbl = intColsTbl, ratTbl = ratColsTbl) {
  #KolNames <- xlTabs[[sheetnm]]
  klnms <- KolNames
  KolNames <- sub("Peptides?", "Pep.", KolNames)
  KolNames <- sub("Evidences?", "PSMs", KolNames)
  KolNames <- sub("Spectr((al)|(um))", "Spec.", KolNames)
  #KolNames <- sub("Razor", "Raz.", KolNames)
  KolNames <- sub("Unique", "Uniq.", KolNames)
  KolNames <- sub("MS/MS", "MS2", KolNames)
  # This would be the place to edit PER sample evidence and MS/MS columns,
  # currently not needed because those do not exist yet
  for (nm in names(intTbl)) { #nm <- names(intTbl)[1L] #nm <- names(intTbl)[2L]
    m <- match(intTbl[[nm]]$Log, KolNames)
    w <- which(!is.na(m))
    if (length(w)) {
      rpl <- intNms(nm, type = "pep", newLine = TRE)
      KolNames[m[w]] <- paste0(rpl, intTbl[[nm]]$Sample[w])
    }
  }
  if ((!missing("ratTbl")) && (!is.null(ratTbl))) {
    for (nm in names(ratTbl)) {
      m <- match(ratTbl[[nm]]$Log, KolNames)
      w <- which(!is.na(m))
      if (length(w)) {
        rpl <- ratNms(nm, newLine = TRUE)
        KolNames[m[w]] <- paste0(rpl, ratTbl[[nm]]$Sample[w])
      }
    }
  }
  #
  if (F.test) {
    wF <- grep("mod\\. F[-\\.]?test:? +", KolNames)
    KolNames[wF] <- sub("mod\\. F[-\\.]?test:? +", "F-test ", KolNames[wF]) # Shorter F-test tag
  }
  KolNames <- sub(" +-log10.*Pvalue\\)", " -log10 pval.", KolNames)
  KolNames <- sub(".*Pvalue", "pval.", KolNames)
  KolNames <- sub(".*Regulated", " reg.", KolNames)
  gSg <- grep("Significant-", KolNames)
  KolNames[gSg] <- sub(" +Significant-", " ", KolNames[gSg])
  #
  KolNames <- sub("^Cluster \\([^\\)]+\\)", "Clust.", KolNames)
  #
  KntKol <- paste0(AA, " Count")
  wI <- which(KolNames %in% KntKol)
  KolNames[wI] <- sub(" Count$", "", KolNames[wI])
  #
  KolNames <- replFun(KolNames)
  KolNames <- gsub(" - ", " ///NL///", KolNames)
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
  return(KolNames)
}
#KolEdit(xlTabs[[sheetnm]], intColsTbl, ratColsTbl)
if ((prot.list.Cond)&&(!"In list" %in% colnames(ev))) {
  g <- grsep2(prot.list, ev$Proteins)
  w <- rep(FALSE, nrow(ev))
  w[g] <- TRUE
  ev$"In list" <- c("", "+")[w+1L]
}
QualFilt %<o% c("In list", "Potential contaminant", "Only identified by site")
if ((DiscFilt)&&(DiscFiltMode == DiscFiltModes[3])) { QualFilt <- c(QualFilt, DiscFiltCols) }
II <- setNames(1, "All peptidoforms")
if ((exists("PTMs_pep"))&&(length(PTMs_pep))) {
  stop("Really? There is no PTM-modified class of peptides to analyze? Why did you run this workflow then?")
} else {
  Mod2Write <- names(PTMs_pep)
  II[paste0(Mod2Write, "-mod. pept.")] <- 1L+(seq_along(length(Mod2Write)))
}
for (ii in II) { #ii <- II[1L] #ii <- II[2L]
  tblMode <- tblMode2 <- "pep"
  TbNm <- names(II)[ii]
  tempData <- get(tblMode)
  intRf <- pep.ref
  ratRf <- pep.ratios.ref
  if (ii > 1L) {
    Ptm <- Mod2Write[[ii-1L]]
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
    if (ii > 1L) { CoreCol <- c(CoreCol, "Name", paste0(Ptm, "-site(s)")) }
    CoreCol2 %<o% c("Leading razor proteins",
                    "Protein names", "Gene names", "Protein group IDs", "Razor protein group ID",
                    grep(" (Probabilities|Score Diffs)$", colnames(tempData), value = TRUE),
                    "Normalisation group")
    evcol <- "Evidence IDs"
    spcol <- "MS/MS count"
    gel <- setNames(lapply(intRf, \(rf) {
      x <- c(paste0(rf, RSA$values),
             paste0("Mean ", rf, VPAL$values))
      return(x[which(x %in% colnames(tempData))])
    }), intRf)
    if (ii == 1L) {
      # Log transform for normal tables - not necessary for PTM-modified tables as we already transformed
      gel2 <- setNames(lapply(intRf, \(rf) {
        x <- c(paste0(rf, RSA$values),
               paste0("Mean ", rf, VPAL$values))
        w <- which(x %in% colnames(tempData))
        rf2 <- sub("Int\\. - ", "log10(Int.) - ", rf)
        x <- c(paste0(rf2, RSA$values),
               paste0("Mean ", rf, VPAL$values))
        return(x[w])
      }), intRf)
      for (rf in intRf) { #rf <- intRf[1L]
        kol1 <- gel[[rf]]
        kol2 <- gel2[[rf]]
        if (length(kol1)) {
          tempData[, kol2] <- suppressWarnings(log10(tempData[, kol1]))
          for (kl in kol2) {
            w <- which(is.infinite(tempData[[kl]]))
            tempData[w, kl] <- NA_real_
          }
        }
      }
      intRf <- sapply(intRf, \(rf) {
        sub("Int\\. - ", "log10(Int.) - ", rf)
      })
      gel <- unlist(gel2)
    }
    #
    smpls <- c(VPAL$values, RSA$values)
    if (length(Exp) == 1L) { smpls <- gsub(topattern(paste0(Exp, "___")), "", smpls) }
    smpls <- gsub("___", " ", smpls)
    intColsTbl <- setNames(lapply(names(intRf), \(nm) { #nm <- names(intRf)[1L]
      res <- data.frame(Log = c(paste0("Mean ", intRf[nm], VPAL$values),
                                paste0(intRf[nm], RSA$values)),
                        Type = c(rep("Average", length(VPAL$values)),
                                 rep("Individual", length(RSA$values))),
                        Sample = smpls)
      w <- which(res$Log %in% colnames(tempData))
      return(res[w,])
    }), names(intRf))
    w <- which(vapply(intColsTbl, nrow, 1L) > 0L)
    intColsTbl <- intColsTbl[w]; intRf <- intRf[w]
    quantCols <- intCols <- lapply(intColsTbl, \(x) { x$Log })
    ratColsTbl <- setNames(lapply(names(ratRf), \(nm) {
      res <- data.frame(Log = c(paste0("Mean ", ratRf[nm], VPAL$values),
                                paste0(ratRf[nm], RSA$values)),
                        Type = c(rep("Average", length(VPAL$values)),
                                 rep("Individual", length(RSA$values))),
                        Sample = smpls)
      w <- which(res$Log %in% colnames(tempData))
      return(res[w,])
    }), names(ratRf))
    w <- which(vapply(ratColsTbl, nrow, 1L) > 0L)
    ratColsTbl <- ratColsTbl[w]; ratRf <- ratRf[w]
    ratCols <- lapply(ratColsTbl, \(x) { x$Log })
    grl <- unlist(ratCols)
    for (gr in grl) {
      w <- which(is.infinite(tempData[[gr]]))
      tempData[w, gr] <- NA_real_
    }
    quantCols[names(ratRf)] <- ratCols
    regcol <- grep("^((Enriched)|(Regulated)) - ", colnames(tempData), value = TRUE)
    signcol <- grep("^Significant-FDR=[1-9][0-9]*\\.*[0-9]*% - ", colnames(tempData), value = TRUE)
    signcol <- grep(" - Analysis_[0-9]+", signcol, invert = TRUE, value = TRUE)
    quantcol <- unlist(quantCols)
    PepColList %<o% c("gel", "grl", "quantcol", "signcol", "regcol") # These are any column for which we want to gsub "___" to " "
    .obj <- unique(c(PepColList, .obj)) # Here easier than using a custom operator
    if (ii > 1L) {
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
    kol <- c(CoreCol, CoreCol2, evcol, spcol, "PEP", quantcol, signcol, regcol, qualFlt, aacol)
    if (ii > 1L) { kol <- c(kol, "Code") }
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
    intColsTbl <- lapply(intColsTbl, \(x) {
      x$Log <- cleanNms(x$Log, start = FALSE)
      x
    })
    ratColsTbl <- lapply(ratColsTbl, \(x) {
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
    if ((ii > 1L)&&(F.test)) {
      tempPepF <- PTMs_F_test_data[[Ptm]]
      tempPepF <- tempPepF[, which(!colnames(tempPepF) %in% c(Param$Plot.labels, "Rel. log10(Peptides count)", "Av. log10 abundance"))]
      colnames(tempPepF) <- cleanNms(colnames(tempPepF), start = FALSE)
      mnratcolF %<o% grep("Mean log2\\(Ratio\\) - ", colnames(tempPepF), value = TRUE)
      m <- match(mnratcolF, colnames(tempPepF))
      colnames(tempPepF)[m] <- paste0("mod. F-test ", mnratcolF)
      mnratcolF <- paste0("mod. F-test ", mnratcolF)
      pvalcolF %<o% F_Root
      signcolF %<o% grep("^mod\\. F-test Significant", colnames(tempPepF), value = TRUE)
      regcolF %<o% grep("^mod\\. F-test Regulated", colnames(tempPepF), value = TRUE)
      Fkol %<o% c(regcolF, mnratcolF, pvalcolF, signcolF)
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
    if (ii > 1L) {
      # - Summary Ratios: P-values and significance
      ColumnsTbl[["P-values"]] <- gpl
      # - Significant
      ColumnsTbl[["Significant"]] <- signcol
      # - Regulated
      ColumnsTbl[["Regulated"]] <- regcol
      # F-test
      if ((ii > 1L)&&(F.test)) {
        ColumnsTbl[["F-test summary Ratios"]] <- mnratcolF
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
      for (i in 1L:nrow(AnnotTbl)) { ColumnsTbl[[paste0(AnnotTbl$Name[i], " annotations")]] <- AnnotTbl$Columns[[i]] }
    }
    # - PEP
    ColumnsTbl[["PEP"]] <- "PEP"
    # - Filters
    ColumnsTbl[["Filters"]] <- qualFlt
    # Melt
    ColumnsTbl <- ColumnsTbl[which(vapply(ColumnsTbl, \(x) { length(x[which(!is.na(x))]) }, 1L) > 0L)]
    ColumnsTbl <- set_colnames(reshape::melt.list(ColumnsTbl), c("Col", "Grp"))
    #tst <- aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), length); View(tst)
    #tst <- aggregate(ColumnsTbl$Grp, list(ColumnsTbl$Col), unique); View(tst)
    #tst <- aggregate(1L:nrow(ColumnsTbl), list(ColumnsTbl$Col), unique); w <- which(lengths(tst$x) > 1L); setNames(tst$x[w], tst$Group.1[w])
    stopifnot(nrow(ColumnsTbl) == length(unique(ColumnsTbl$Col)))
    ColumnsTbl$Class <- ""
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "IDs")] <- "General Peptides information"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(evcol))] <- "Evidence IDs"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% c(spcol))] <- "Spectral count"
    for (nm in names(intRf)) { #nm <- names(intRf)[1L]
      rpl <- intNms(nm, TRUE, type = "pep")
      ColumnsTbl$Class[grep(topattern(intRf[nm]), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[grep(topattern(paste0("Mean ", intRf[nm])), ColumnsTbl$Col)] <- rpl
    }
    for (nm in names(ratRf)) { #nm <- names(ratRf)[1L]
      rpl <- ratNms(nm, TRUE)
      ColumnsTbl$Class[grep(topattern(ratRf[nm]), ColumnsTbl$Col)] <- rpl
      ColumnsTbl$Class[grep(topattern(paste0("Mean ", ratRf[nm])), ColumnsTbl$Col)] <- rpl
    }
    ColumnsTbl$Class[which(ColumnsTbl$Grp == "P-values")] <- gsub(" - $", "", pvalue.col[which(pvalue.use)])
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% regcol)] <- "Regulated"
    ColumnsTbl$Class[which(ColumnsTbl$Col %in% signcol)] <- "Significant"
    if ((ii > 1L)&&(F.test)) {
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
    Src <- paste0(libPath, "/extdata/Sources/fstWrite_Excel_core_script.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
}
#
Src <- paste0(libPath, "/extdata/Sources/Write_Excel_end_script.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#WorkBook$get_active_sheet()
#xl_open(repFl)
