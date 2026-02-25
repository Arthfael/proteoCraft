# Filter for cases where number of characters is more than the Excel limit:
## Should only ever be character columns (lists of spectra IDs, sequences...)
w <- which(sapply(colnames(tempData), function(x) { class(tempData[[x]]) }) == "character")
if (length(w)) {
  for (i in w) { #i <- w[1]
    w1 <- which(nchar(tempData[[colnames(tempData)[i]]]) > ExcelMax)
    if (length(w1)) {
      tempData[[colnames(tempData)[i]]][w1] <- paste0(substr(tempData[[colnames(tempData)[i]]][w1], 1, ExcelMax-3), "...")
    }
  }
}
#
tempData$id <- as.numeric(tempData$id) # Safety to avoid issues with opening Excel file later
tempData <- tempData[order(tempData$id, decreasing = FALSE),]
#
LS <- c("sheetnmsA", "sheetnmsB", "xlTabs", "whLevel", "data_filt", "data_order")
suppressWarnings(rm(list = LS))
# Peptidoform file with 1 tab for each class of information
## Full Protein Groups file (1 tab for all columns)
xlTabs <- list()
xlTabs[[TbNm]] <- colnames(tempData)
if (tblMode == "PG") {
  # We now want coverage columns in another tab - they are otherwise slowing the main tab too much
  xlTabs$Coverage <- c(CoreCol, CoreCol2, covcol)
  tmp <- xlTabs[[TbNm]]
  tmp <- tmp[which(!tmp %in% covcol)]
  xlTabs[[TbNm]] <- tmp
}
# Write tab, apply formatting etc...
for (l in LS) { if (!exists(l)) { assign(l, list()) } }
sheetnmsA <- names(xlTabs)
sheetnmsB <- sheetnmsA[which(!sheetnmsA %in% c("Description", "Quality control"))]
cat(paste0("   Building ", tblMode2, " tabs...\n"))
if (tblMode == "pep") {
  datCol <- ColumnsTbl$Col[unique(which(ColumnsTbl$Class %in% c(intNms(names(intRf), TRUE, "pep"),
                                                                ratNms(names(ratRf), TRUE))))]
}
if (tblMode == "PG") {
  datCol <- ColumnsTbl$Col[unique(c(which(ColumnsTbl$Class %in% c(intNms(names(intColsTbl), TRUE),
                                                                  ratNms(names(ratColsTbl), TRUE))),
                                    grep("est\\. copies/cell", ColumnsTbl$Class)))]
}
stopifnot(length(datCol) > 0)
for (sheetnm in sheetnmsB) { #sheetnm <- sheetnmsB[1] #sheetnm <- sheetnmsB[2]
  # Very important early filter, to avoid discrepancies between header and columns!
  # Should be done at least before the header is written
  xlTabs[[sheetnm]] <- xlTabs[[sheetnm]][which(xlTabs[[sheetnm]] %in% colnames(tempData))]
  xlTabs[[sheetnm]] <- ColumnsTbl$Col[which(ColumnsTbl$Col %in% xlTabs[[sheetnm]])] # The order in ColumnsTbl should make sense!
  #
  if (!sheetnm %in% names(data_filt)) { data_filt[[sheetnm]] <- 1:nrow(tempData) }
  if (!sheetnm %in% names(data_order)) { data_order[[sheetnm]] <- 1:length(data_filt[[sheetnm]]) }
  if (sheetnm %in% wb_get_sheet_names(WorkBook)) { WorkBook <- wb_remove_worksheet(WorkBook, sheetnm) }
  #WorkBook <- wb_add_worksheet(WorkBook, sheetnm)
  #saveFun(WorkBook, file = "WorkBook_bckp.RData")
  #loadFun("WorkBook_bckp.RData")
  # Vertically center values in cells
  wtst <- data_filt[[sheetnm]]
  sheetrows <- length(wtst)+2
  ord <- data_order[[sheetnm]]
  rowRg <- 1:sheetrows
  colRg <- 1:length(xlTabs[[sheetnm]])
  #
  # Write the data
  # Styles list which applies to this tab
  ColumnsTbl2 <- data.frame(Col = xlTabs[[sheetnm]], Grp = "", Class = "")
  w <- which(ColumnsTbl2$Col %in% ColumnsTbl$Col)
  ColumnsTbl2[w, c("Grp", "Class")] <- ColumnsTbl[match(xlTabs[[sheetnm]][w], ColumnsTbl$Col), c("Grp", "Class")] 
  # Test for series of consecutive columns from the same class
  tstKol <- lapply(1:nrow(ColumnsTbl2), function(x) { #x <- 1
    res <- x
    w <- which(ColumnsTbl2$Class == ColumnsTbl2$Class[x])
    y <- x-1
    while (y %in% w) {
      res <- c(y, res)
      y <- y-1
    }
    y <- x+1
    while (y %in% w) {
      res <- c(res, y)
      y <- y+1
    }
    return(paste(c(min(res), max(res)), collapse = ";"))
  })
  tstKol <- as.data.frame(t(sapply(strsplit(unique(unlist(tstKol)), ";"), as.numeric)))
  colnames(tstKol) <- c("First", "Last")
  tstKol$Group <- ColumnsTbl2$Class[tstKol$First]
  tstKol <- tstKol[which(tstKol$Group != ""), , drop = FALSE]
  whLevel[[sheetnm]] <- rep(1, length(xlTabs[[sheetnm]])) # For now all columns are to write in the 1st header row; this may change
  # Prepare 1st level header
  m <- match(xlTabs[[sheetnm]], ColumnsTbl$Col) #m <- match(ColumnsTbl$Col, ColumnsTbl$Col)
  tmpKol2 <- ColumnsTbl$edit_Col[m]
  w <- which(tmpKol2 == "")
  if (length(w)) {
    clss <- ColumnsTbl$Class[match(xlTabs[[sheetnm]][w], ColumnsTbl$Col)]
    w2 <- match(clss, tstKol$Group)
    if (tstKol$First[w2] == tstKol$Last[w2]) {
      tmpKol2[w] <- clss
      tstKol <- tstKol[-w2]
    }
  }
  tmpKol2 <- setNames(tmpKol2, NULL)
  temp <- tempData[wtst[ord], xlTabs[[sheetnm]]]
  colnames(temp) <- tmpKol2 #Edit column names
  #
  if (sheetnm %in% wb_get_sheet_names(WorkBook)) { WorkBook <- wb_remove_worksheet(WorkBook, sheetnm) }
  WorkBook <- wb_add_worksheet(WorkBook, sheetnm)
  #
  tblNm <- tolower(gsub(" |\\$", ".", sheetnm))
  WorkBook <- wb_add_data_table(WorkBook, sheetnm, temp,
                                dims = wb_dims(rows = 2, cols = 1),
                                table_name = tblNm,
                                col_names = TRUE, table_style = "TableStyleMedium2",
                                banded_rows = TRUE, banded_cols = FALSE)
  #wb_save(WorkBook, repFl);xl_open(repFl)
  #
  # Attempt to override the banding in the data columns,
  # because #N/A values are not affected by conditional formatting
  # and seeing them in alternate colours looks... weird.
  tst <- data.frame(IsQuant = xlTabs[[sheetnm]] %in% datCol)
  tst$IsChange <- c(FALSE, sapply(2:nrow(tst), function(x) {
    as.numeric(tst$IsQuant[x]!= tst$IsQuant[x-1])
  }))
  tst$Block <- cumsum(tst$IsChange)
  blocks <- aggregate(1:nrow(tst), list(tst$Block), unique)
  colnames(blocks) <- c("Block", "Columns")
  blocks$IsQuant <- tst$IsQuant[match(blocks$Block, tst$Block)]
  w <- which(blocks$IsQuant)
  for (i in w) {
    WorkBook <- wb_add_fill(WorkBook, sheetnm,
                            wb_dims(rows = (1:nrow(temp))+2, cols = blocks$Columns[[i]]),
                            color = wb_color("yellow"))
  }
  #
  # Write 1st level header
  if (nrow(tstKol)) {
    for (i in 1:nrow(tstKol)) {
      WorkBook <- wb_add_data(WorkBook, sheetnm, tstKol$Group[i], wb_dims(1, tstKol$First[i])) # Write 1st header level
      if (tstKol$Last[i] > tstKol$First[i]) { # Merge if 1st header level is same for several columns
        WorkBook <- wb_merge_cells(WorkBook, sheetnm, wb_dims(1, tstKol$First[i]:tstKol$Last[i]))
      }
    }
  } else { stop() }
  #wb_save(WorkBook, repFl);xl_open(repFl)
  #
  # Style the data
  # - vertical align to center
  WorkBook <- wb_add_cell_style(WorkBook, sheetnm, wb_dims(rowRg, colRg), vertical = "center")
  # - Header borders
  WorkBook <- wb_add_border(WorkBook, sheetnm, wb_dims(1:2, colRg))
  # 
  stl <- wb_get_cell_style(WorkBook, "tmp", wb_dims(match("Header 1", styleNms), 1))
  if (stl != "") { WorkBook <- wb_set_cell_style(WorkBook, sheetnm, wb_dims(2, colRg), stl) }
  stl <- wb_get_cell_style(WorkBook, "tmp", wb_dims(match("Header 2", styleNms), 1))
  if (stl != "") { WorkBook <- wb_set_cell_style(WorkBook, sheetnm, wb_dims(1, colRg), stl) }
  kl <- "id"
  if (tblMode == "pep") { kl <- c(kl, aacol) }
  w <- which(xlTabs[[sheetnm]] %in% kl) # Peptide ID column, Amino Acid counts
  if (length(w)) { WorkBook <- wb_set_col_widths(WorkBook, sheetnm, w, 5) }
  # Expression and annotations columns, PEP, some others...
  kl <- c(quantcol, "PEP", CoreCol[which(CoreCol != "id")])
  if (Annotate) { kl <- c(kl, annot.col2) }
  w <- unique(c(which(xlTabs[[sheetnm]] %in% kl),
                grep("Regulated - ", xlTabs[[sheetnm]]))) # Specific columns
  if (length(w)) { WorkBook <- wb_set_col_widths(WorkBook, sheetnm, w, 15) }
  #
  if (tblMode == "PG") {
    w <- which(xlTabs[[sheetnm]] %in% xmlCovCol)
    if (length(w)) {
      WorkBook <- wb_set_col_widths(WorkBook, sheetnm, w, 100)
      WorkBook <- wb_add_cell_style(WorkBook, sheetnm, wb_dims(3:sheetrows, w), vertical = "top")
    }
    w <- which(xlTabs[[sheetnm]] %in% covcol[which(!covcol %in% xmlCovCol)])
    if (length(w)) { WorkBook <- wb_set_col_widths(WorkBook, sheetnm, w, 11) }
  }
  #
  w <- grep("([Ee]vidence|[Pp]eptide) IDs", xlTabs[[sheetnm]]) # Peptide and protein IDs columns
  if (length(w)) { WorkBook <- wb_set_col_widths(WorkBook, sheetnm, w, 10) }
  #
  w <- which(xlTabs[[sheetnm]] %in% qualFlt) # Quality filter, contaminants
  if (length(w)) { WorkBook <- wb_set_col_widths(WorkBook, sheetnm, w, 14) }
  WorkBook <- wb_set_row_heights(WorkBook,
                                 sheetnm,
                                 rowRg[1:2],
                                 ceiling(c(max(nchar(tstKol$Group)),
                                           max(nchar(tstKol)))/20)*20)
  WorkBook <- wb_add_cell_style(WorkBook, sheetnm, "A1", vertical = "top")
  WorkBook <- wb_freeze_pane(WorkBook, sheetnm, first_active_row  = 3, firstActiveCol = length(CoreCol)+1)
  #wb_save(WorkBook, repFl);xl_open(repFl)
  w1 <- which(!is.na(sapply(xlTabs[[sheetnm]], function(x) { #x <- xlTabs[[sheetnm]][1]
    x <- match(ColumnsTbl$Grp[match(x, ColumnsTbl$Col)], styleNms)
  })))
  w2 <- colRg
  w2 <- w2[which(!w2 %in% w1)]
  CS1 <- setNames(lapply(xlTabs[[sheetnm]][w1], function(x) { #x <- xlTabs[[sheetnm]][w1][1]
    m <- match(x, ColumnsTbl$Col)
    res <- list(Group = ColumnsTbl$Grp[m],
                Class = ColumnsTbl$Class[m])
    res$tmpCell <- wb_dims(match(res$Group, styleNms), 1)
    res$Style <- wb_get_cell_style(WorkBook, "tmp", res$tmpCell)
    return(res)
  }), xlTabs[[sheetnm]][w1])
  CS2 <- setNames(lapply(xlTabs[[sheetnm]][w2], function(x) { #x <- xlTabs[[sheetnm]][w2][1]
    m <- match(x, ColumnsTbl$Col)
    res <- list(Group = ColumnsTbl$Grp[m],
                Class = ColumnsTbl$Class[m])
    res$Style <- Styles[[res$Group]]
    return(res)
  }), xlTabs[[sheetnm]][w2])
  CS1tst <- aggregate(1:length(CS1), list(sapply(CS1, paste, collapse = "-")), function(x) {
    c(min(x), max(x))
  })
  CS1tst[, c("Min", "Max")] <- as.data.frame(CS1tst$x)
  CS1tst$x <- NULL
  CS1tst <- CS1tst[order(CS1tst$Min),]
  CS1tst$Range <- apply(CS1tst[, c("Min", "Max")], 1, function(x) { list(w1[x[[1]]:x[[2]]]) })
  CS2tst <- aggregate(1:length(CS2), list(sapply(CS2, paste, collapse = "-")), function(x) {
    c(min(x), max(x))
  })
  CS2tst[, c("Min", "Max")] <- as.data.frame(CS2tst$x)
  CS2tst$x <- NULL
  CS2tst <- CS2tst[order(CS2tst$Min),]
  CS2tst$Range <- apply(CS2tst[, c("Min", "Max")], 1, function(x) { list(w2[x[[1]]:x[[2]]]) })
  for (rw in 1:nrow(CS1tst)) {
    mn <- CS1tst$Min[rw]
    mx <- CS1tst$Max[rw]
    rg <- unlist(CS1tst$Range[[rw]])
    xlRg <- wb_dims(3:sheetrows, rg)
    WorkBook <- wb_set_cell_style(WorkBook, sheetnm, xlRg, CS1[[mn]]$Style)
  }
  #unlist(lapply(names(CS2), function(nm) { as.character(CS2[[nm]]$Style) }))
  #saveFun(WorkBook, file = "WorkBook_bckp.RData")
  #wb_save(WorkBook, repFl);xl_open(repFl)
  #loadFun("WorkBook_bckp.RData")
  for (rw in 1:nrow(CS2tst)) {
    mn <- CS2tst$Min[rw]
    mx <- CS2tst$Max[rw]
    rg <- unlist(CS2tst$Range[[rw]])
    xlRg <- wb_dims(3:sheetrows, rg)
    if ("Style" %in% names(CS2[[mn]])) {
      if (CS2[[mn]]$Style %in% names(DecoList)) {
        deco <- DecoList[[CS2[[mn]]$Style]]
        WorkBook <- wb_add_font(WorkBook, sheetnm, xlRg,
                                italic = c("", "true")[("italic" %in% deco)+1],
                                bold = c("", "true")[("bold" %in% deco)+1])
      }
      if (CS2[[mn]]$Style %in% names(SignifList)) {
        signi <- SignifList[[CS2[[mn]]$Style]]
        WorkBook <- wb_add_numfmt(WorkBook, sheetnm, xlRg,
                                  paste(c("0.", rep(0, signi)), collapse = ""))
      }
      if (CS2[[mn]]$Style %in% names(ColScaleList)) {
        condFrmt <- ColScaleList[[CS2[[mn]]$Style]]
        WorkBook <- wb_add_conditional_formatting(WorkBook, sheetnm, xlRg,
                                                  type = "colorScale", style = condFrmt)
        # WorkBook <- wb_add_font(WorkBook, sheetnm, xlRg, color = wb_color(hex = "FF000000"))
      }
      if (CS2[[mn]]$Style %in% names(ColList)) {
        clr <- ColList[[CS2[[mn]]$Style]]
        WorkBook <- wb_add_font(WorkBook, sheetnm, xlRg, color = wb_color(hex = "white"))
      }
      if (CS2[[mn]]$Style %in% names(HAlignList)) {
        hal <- HAlignList[[CS2[[mn]]$Style]]
        WorkBook <- wb_add_cell_style(WorkBook, sheetnm, xlRg
                                      , horizontal = hal)
      }
      if (CS2[[mn]]$Style %in% names(ContainsList)) {
        cntn <- ContainsList[[CS2[[mn]]$Style]]
        stl <- capture.output(cntn$style)
        deco <- tolower(gsub(".* ", "", gsub(" *$", "", grep("Font decoration", stl, value = TRUE))))
        cntnFrt <- gsub(".* ", "", gsub(" *$", "", grep("Font colour", stl, value = TRUE)))
        cntnBgrd <- gsub(".* ", "", gsub(" *$", "", grep("Cell fill background", stl, value = TRUE)))
        cntnHrz <- gsub(".* ", "", gsub(" *$", "", grep("Cell horz. align", stl, value = TRUE)))
        WorkBook <- wb_add_dxfs_style(WorkBook, CS2[[mn]]$Style, font_color = wb_color(cntnFrt),
                                      bg_fill = wb_color(cntnBgrd),
                                      text_bold = c("", "single")[("bold" %in% deco)+1],
                                      text_italic = c("", "italic")[("italic" %in% deco)+1])
        WorkBook <- wb_add_conditional_formatting(WorkBook, sheetnm, xlRg,
                                                  type = "containsText", style = CS2[[mn]]$Style, rule = cntn$rule)
      }
    }
  }
  #wb_save(WorkBook, repFl);xl_open(repFl)
  #
  if (prot.list.Cond) {
    kol <- paste0("Leading protein", c("", " ID")[match(tblMode, c("pep", "PG"))], "s")
    w <- grsep2(prot.list, tempData[wtst[ord], kol])+2
    if (length(w)) {
      stl <- wb_get_cell_style(WorkBook, "tmp", wb_dims(match("Protein list", styleNms), 1))
      if (stl != "") {
        for (i in w) {
          dims <- wb_dims(i, colRg)
          WorkBook <- wb_set_cell_style(WorkBook, sheetnm, dims, stl)
          if (tblMode == "PG") {
            WorkBook <- wb_set_row_heights(WorkBook, sheetnm, dims, 50)
          }
        }
      }
    }
  }
  #
  # - Groupings
  tst <- aggregate(ColumnsTbl$Hide, list(ColumnsTbl$Class), unique)
  stopifnot("logical" %in% class(tst$x))
  colnames(tst) <- c("Group", "Hide")
  tst$Col <- lapply(tst$Group, function(x) {
    m <- which(ColumnsTbl$Class == x)
    return(which(xlTabs[[sheetnm]] %in% ColumnsTbl$Col[m]))
  })
  tst <- tst[which(sapply(tst$Col, length) > 0),]
  filt <- c("General Peptides information", unique(ColumnsTbl$Class[match(intCols[[rev(names(intCols))[1]]], ColumnsTbl$Col)]), "QC filters")
  if (MakeRatios) {
    filt <- c(filt, unique(ColumnsTbl$Class[match(ratCols[[rev(names(ratCols))[1]]], ColumnsTbl$Col)]), "Regulated")
  }
  tst <- tst[which(!tst$Group %in% filt),]
  if (nrow(tst)) {
    tst$Min <- sapply(tst$Col, min)
    tst$Max <- sapply(tst$Col, max)
    tst$Max <- tst$Max+1 # For openxlsx2
    tst <- tst[order(tst$Min, decreasing = FALSE),]
    for (i in 1:nrow(tst)) { WorkBook <- wb_group_cols(WorkBook, sheetnm, tst$Min[i]:tst$Max[i], tst$Hide[i]) }
  }
  #
  if (tblMode == "PG") {
    # - Coverage columns:
    #   - Xml encoded coverage
    kols <- xmlCovCol
    kols <- kols[which(kols %in% xlTabs[[sheetnm]])]
    if (length(kols)) {
      for (kol in kols) {
        dms <- wb_dims(3, match(kol, xlTabs[[sheetnm]]))
        WorkBook <- wb_add_data(WorkBook, sheetnm, tempData[wtst[ord], kol], dms)
      }
    }
    #   - Data bar summaries
    kols <- grep("\\[%\\]", covcol, value = TRUE)
    kols <- kols[which(kols %in% xlTabs[[sheetnm]])]
    if (length(kols)) {
      for (kol in kols) {
        dms <- wb_dims(3:sheetrows, match(kol, xlTabs[[sheetnm]]))
        WorkBook <- wb_add_conditional_formatting(WorkBook, sheetnm, dms, type = "dataBar",
                                                  style = c("FF63C384", "0063C384"), rule = c(0, 100))
      }
    }
  }
}
