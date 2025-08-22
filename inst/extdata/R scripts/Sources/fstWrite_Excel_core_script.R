
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
if (tblMode == "SAINTexpress") { idKl <- "Protein" } else {
  idKl <- "id"
  tempData[[idKl]] <- as.numeric(tempData$id) # Safety to avoid issues with opening Excel file later
  tempData <- tempData[order(tempData[[idKl]], decreasing = FALSE),]
}
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
cat(paste0("   Building ", tblMode2, " tab...\n"))
if (tblMode == "pep") {
  nms <- intNms(names(intRf), TRUE, "pep")
  if (MakeRatios) {
    nms <- c(nms, ratNms(names(ratRf), TRUE))
  }
  datCol <- ColumnsTbl$Col[unique(which(ColumnsTbl$Class %in% nms))]
}
if (tblMode == "PG") {
  nms <- intNms(names(intColsTbl), TRUE)
  if (MakeRatios) {
    nms <- c(nms, ratNms(names(ratColsTbl), TRUE))
  }
  datCol <- ColumnsTbl$Col[unique(c(which(ColumnsTbl$Class %in% nms),
                                    grep("est\\. copies/cell", ColumnsTbl$Class)))]
}
if (tblMode == "SAINTexpress") {
  datCol <- ColumnsTbl$Col[which(ColumnsTbl$Class == "log2(rat.), avg.")]
}
stopifnot(length(datCol) > 0)
wrtHeader <- TRUE
for (sheetnm in sheetnmsB) { #sheetnm <- sheetnmsB[1] #sheetnm <- sheetnmsB[2]
  # Very important early filter, to avoid discrepancies between header and columns!
  # Should be done at least before the header is written
  xlTabs[[sheetnm]] <- xlTabs[[sheetnm]][which(xlTabs[[sheetnm]] %in% colnames(tempData))]
  xlTabs[[sheetnm]] <- ColumnsTbl$Col[which(ColumnsTbl$Col %in% xlTabs[[sheetnm]])] # The order in ColumnsTbl should make sense!
  #xlTabs[[sheetnm]] <- colnames(tempData)[which(colnames(tempData) %in% xlTabs[[sheetnm]])] # We now define the order not in ColumnsTbl but when we build tempData
  # This may be reverted, but is the current trend!!!
  #
  if (!sheetnm %in% names(data_filt)) { data_filt[[sheetnm]] <- 1:nrow(tempData) }
  if (!sheetnm %in% names(data_order)) { data_order[[sheetnm]] <- 1:length(data_filt[[sheetnm]]) }
  #WorkBook <- wb_add_worksheet(WorkBook, sheetnm)
  #saveFun(WorkBook, file = "WorkBook_bckp.RData")
  #loadFun("WorkBook_bckp.RData")
  # Vertically center values in cells
  wtst <- data_filt[[sheetnm]]
  sheetrows <- length(wtst)+2
  ord <- data_order[[sheetnm]]
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
  #xlTabs[[sheetnm]][which(is.na(m))]
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
  if (!nrow(tstKol)) { stop() }
  tmpKol2 <- setNames(tmpKol2, NULL)
  myData <- tempData[wtst[ord], xlTabs[[sheetnm]]]
  colnames(myData) <- tmpKol2 #Edit column names
  nRws <- nrow(myData)
  nCol <- ncol(myData)
  tblRws <- c(1, nRws)
  hdRg <- c(1, 2)
  tblRws <- tblRws + 2
  # Round PEPs, we don't need all of those zeros!!!
  if ("PEP" %in% colnames(myData)) { myData$PEP <- round(myData$PEP, 6) }
  #
  # NB: Excel can take 1048576 rows max!
  # At some point we may need a solution for if we ever have more rows!
  if (nRws > 1048576 - 2) {
    stop("Unhandled extreme case: more rows than the maximum allowed limit (= 1048576) in Excel!")
  }
  #
  if (sheetnm %in% wb_get_sheet_names(WorkBook)) { WorkBook <- wb_remove_worksheet(WorkBook, sheetnm) }
  WorkBook <- wb_add_worksheet(WorkBook, sheetnm)
  #WB <- openxlsx2::wb_remove_worksheet(WB, sheetnm);WB <- openxlsx2::wb_add_worksheet(WB, sheetnm)
  #
  #WB <- openxlsx2::wb_load(fl);WB <- openxlsx2::wb_add_worksheet(WB, sheetnm)
  #
  tblNm <- tolower(gsub(" |\\$|\\.|-", "_", sheetnm))
  dims <- openxlsx2::wb_dims(rows = hdRg[2], cols = 1)
  # We will not write the first row of the data, but a nicely behaved dummy row without any NAs, NaNs or Inf...
  dummyData <- myData[1, , drop = FALSE]
  tst <- sapply(colnames(myData), function(x) { sum(c("numeric", "integer") %in% class(myData[[x]])) }) > 0
  wNum <- setNames(which(tst), NULL)
  wTxt <- setNames(which(!tst), NULL)
  dummyData[wNum] <- 0
  dummyData[wTxt] <- "Hello world!"
  WorkBook <- openxlsx2::wb_add_data_table(WorkBook,
                                           sheetnm,
                                           dummyData,
                                           dims,
                                           table_name = tblNm,
                                           table_style = "TableStyleMedium2",
                                           banded_rows = TRUE,
                                           banded_cols = FALSE)
  #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
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
    dims <- openxlsx2::wb_dims(rows = tblRws[1]#:tblRws[2]
                               , cols = blocks$Columns[[i]])
    WorkBook <- openxlsx2::wb_add_fill(WorkBook, sheetnm, dims, color = openxlsx2::wb_color("yellow"))
  }
  #
  # Write top level header
  if (wrtHeader) {
    for (i in 1:nrow(tstKol)) {
      dims <- openxlsx2::wb_dims(hdRg[1], tstKol$First[i])
      WorkBook <- openxlsx2::wb_add_data(WorkBook, sheetnm, tstKol$Group[i], dims) # Write 1st header level
      if (tstKol$Last[i] > tstKol$First[i]) { # Merge if 1st header level is same for several columns
        dims <- openxlsx2::wb_dims(hdRg[1], tstKol$First[i]:tstKol$Last[i])
        WorkBook <- openxlsx2::wb_merge_cells(WorkBook, sheetnm, dims)
      }
    }
  }
  #openxlsx2::wb_save(WorkBook, repFl);xl_open(repFl)
  #
  # Style the table
  # - vertical align to center
  dims <- openxlsx2::wb_dims(1#:tblRws[2]
                             , colRg)
  WorkBook <- openxlsx2::wb_add_cell_style(WorkBook, sheetnm, dims, vertical = "center")
  #a1 <- WorkBook$worksheets[[3]]$sheet_data$cc; a2 <- WorkBook$worksheets[[3]]$sheet_data$cc; length(unique(a1$c_s)); length(unique(a2$c_s))
  kl <- "id" # Not idKl!
  if (tblMode == "pep") { kl <- c(kl, aacol) }
  w <- which(xlTabs[[sheetnm]] %in% kl) # Peptide ID column, Amino Acid counts
  if (length(w)) { WorkBook <- openxlsx2::wb_set_col_widths(WorkBook, sheetnm, w, 5) }
  # Expression and annotations columns, PEP, some others...
  if (tblMode != "SAINTexpress") { kl <- c(quantcol, "PEP") }
  kl <- c(quantcol, CoreCol[which(CoreCol != "id")])
  if (tblMode != "SAINTexpress") {
    if (Annotate) { kl <- c(kl, annot.col2) }
  }
  w <- unique(c(which(xlTabs[[sheetnm]] %in% kl),
                grep("Regulated - ", xlTabs[[sheetnm]]))) # Specific columns
  if (length(w)) { WorkBook <- openxlsx2::wb_set_col_widths(WorkBook, sheetnm, w, 15) }
  #
  if (tblMode == "PG") {
    w <- which(xlTabs[[sheetnm]] %in% xmlCovCol)
    if (length(w)) {
      WorkBook <- openxlsx2::wb_set_col_widths(WorkBook, sheetnm, w, 100)
      dims <- openxlsx2::wb_dims(tblRws[1]#:tblRws[2]
                                 , w)
      WorkBook <- openxlsx2::wb_add_cell_style(WorkBook, sheetnm,
                                               dims, vertical = "top")
    }
    w <- which(xlTabs[[sheetnm]] %in% covcol[which(!covcol %in% xmlCovCol)])
    if (length(w)) { WorkBook <- openxlsx2::wb_set_col_widths(WorkBook, sheetnm, w, 11) }
  }
  if (tblMode == "SAINTexpress") {
    w <- which(xlTabs[[sheetnm]] %in% boostcol)		
    if (length(w)) {
      WorkBook <- openxlsx2::wb_set_col_widths(WorkBook, sheetnm, w, 50)
    }
  }
  #
  w <- grep("([Ee]vidence|[Pp]eptide) IDs", xlTabs[[sheetnm]]) # Peptide and protein IDs columns
  if (length(w)) { WorkBook <- openxlsx2::wb_set_col_widths(WorkBook, sheetnm, w, 10) }
  #
  w <- which(xlTabs[[sheetnm]] %in% qualFlt) # Quality filter, contaminants
  if (length(w)) { WorkBook <- openxlsx2::wb_set_col_widths(WorkBook, sheetnm, w, 14) }
  WorkBook <- openxlsx2::wb_add_cell_style(WorkBook, sheetnm, "A1", vertical = "top")
  if (wrtHeader) {
    # Style the header
    # - Header borders
    WorkBook <- openxlsx2::wb_add_border(WorkBook, sheetnm,
                                         openxlsx2::wb_dims(hdRg, colRg))
    #
    dims <- openxlsx2::wb_dims(match("Header 1", styleNms), 1)
    stl <- openxlsx2::wb_get_cell_style(WorkBook, "tmp", dims)
    if (stl != "") {
      dims <- openxlsx2::wb_dims(hdRg[2], colRg)
      WorkBook <- openxlsx2::wb_set_cell_style(WorkBook, sheetnm, dims, stl)
    }
    dims <- openxlsx2::wb_dims(match("Header 2", styleNms), 1)
    stl <- openxlsx2::wb_get_cell_style(WorkBook, "tmp", dims)
    if (stl != "") {
      dims <- openxlsx2::wb_dims(hdRg[1], colRg)
      WorkBook <- openxlsx2::wb_set_cell_style(WorkBook, sheetnm, dims, stl)
    }
    #
    # Row heights:
    # This should be done differently:
    # Figure out how many characters can be fitted in each column's width, then adjust column height accordingly.
    # For now this heuristic will do.
    WorkBook <- openxlsx2::wb_set_row_heights(WorkBook,
                                              sheetnm,
                                              hdRg,
                                              c(max(c(30, min(c(ceiling(max(nchar(tstKol$Group))/20)*10, 45)))),
                                                max(c(60, min(c(ceiling(max(nchar(tmpKol2))/10)*20, 180))))))
  }
  # Freeze panes
  m <- match(CoreCol, xlTabs[[sheetnm]])
  WorkBook <- openxlsx2::wb_freeze_pane(WorkBook, sheetnm, first_active_row = c(1, 3)[wrtHeader+1],
                                        firstActiveCol = max(m, na.rm = TRUE)+1)
  #
  #openxlsx2::wb_save(WorkBook, repFl);xl_open(repFl)
  #
  # Conditional formatting
  w1 <- which(!is.na(sapply(xlTabs[[sheetnm]], function(x) { #x <- xlTabs[[sheetnm]][1]
    x <- match(ColumnsTbl$Grp[match(x, ColumnsTbl$Col)], styleNms)
  })))
  w2 <- colRg
  w2 <- w2[which(!w2 %in% w1)]
  if (length(w1)) {
    CS1 <- setNames(lapply(xlTabs[[sheetnm]][w1], function(x) { #x <- xlTabs[[sheetnm]][w1][1]
      m <- match(x, ColumnsTbl$Col)
      res <- list(Group = ColumnsTbl$Grp[m],
                  Class = ColumnsTbl$Class[m])
      res$tmpCell <- openxlsx2::wb_dims(match(res$Group, styleNms), 1)
      res$Style <- openxlsx2::wb_get_cell_style(WorkBook, "tmp", res$tmpCell)
      return(res)
    }), xlTabs[[sheetnm]][w1])
    CS1tst <- aggregate(1:length(CS1), list(sapply(CS1, paste, collapse = "-")), function(x) {
      c(min(x), max(x))
    })
    CS1tst[, c("Min", "Max")] <- as.data.frame(CS1tst$x)
    CS1tst$x <- NULL
    CS1tst <- CS1tst[order(CS1tst$Min),]
    CS1tst$Range <- apply(CS1tst[, c("Min", "Max")], 1, function(x) { list(w1[x[[1]]:x[[2]]]) })
    for (rw in 1:nrow(CS1tst)) {
      mn <- CS1tst$Min[rw]
      mx <- CS1tst$Max[rw]
      rg <- unlist(CS1tst$Range[[rw]])
      dims <- openxlsx2::wb_dims(tblRws[1]#:tblRws[2]
                                 , rg)
      WB <- openxlsx2::wb_set_cell_style(WorkBook, sheetnm, dims, CS1[[mn]]$Style)
    }
  }
  if (length(w2)) {
    CS2 <- setNames(lapply(xlTabs[[sheetnm]][w2], function(x) { #x <- xlTabs[[sheetnm]][w2][1]
      m <- match(x, ColumnsTbl$Col)
      res <- list(Group = ColumnsTbl$Grp[m],
                  Class = ColumnsTbl$Class[m])
      res$Style <- Styles[[res$Group]]
      return(res)
    }), xlTabs[[sheetnm]][w2])
    CS2tst <- aggregate(1:length(CS2), list(sapply(CS2, paste, collapse = "-")), function(x) {
      c(min(x), max(x))
    })
    CS2tst[, c("Min", "Max")] <- as.data.frame(CS2tst$x)
    CS2tst$x <- NULL
    CS2tst <- CS2tst[order(CS2tst$Min),]
    CS2tst$Range <- apply(CS2tst[, c("Min", "Max")], 1, function(x) { list(w2[x[[1]]:x[[2]]]) })
    #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
    #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
    # for (rw in 1:nrow(CS2tst)) { #rw <- 1
    #   print(CS2tst$Group.1[[rw]])
    # }
    for (rw in 1:nrow(CS2tst)) { #rw <- 1
      mn <- CS2tst$Min[rw]
      mx <- CS2tst$Max[rw]
      rg <- unlist(CS2tst$Range[[rw]])
      dims <- openxlsx2::wb_dims(tblRws[1]#:tblRws[2]
                                 , rg)
      if ("Style" %in% names(CS2[[mn]])) {
        if (CS2[[mn]]$Style %in% names(DecoList)) {
          deco <- DecoList[[CS2[[mn]]$Style]]
          WorkBook <- openxlsx2::wb_add_font(WorkBook, sheetnm, dims,
                                             italic = c("", "true")[("italic" %in% deco)+1],
                                             bold = c("", "true")[("bold" %in% deco)+1])
        }
        if (CS2[[mn]]$Style %in% names(SignifList)) {
          signi <- SignifList[[CS2[[mn]]$Style]]
          WorkBook <- openxlsx2::wb_add_numfmt(WorkBook, sheetnm, dims,
                                               paste(c("0.", rep(0, signi)), collapse = ""))
        }
        if (CS2[[mn]]$Style %in% names(ColScaleList)) {
          condFrmt <- ColScaleList[[CS2[[mn]]$Style]]
          WorkBook <- openxlsx2::wb_add_conditional_formatting(WorkBook, sheetnm, dims,
                                                               type = "colorScale", style = condFrmt)
          #WB <- openxlsx2::wb_add_font(WB, sheetnm, dims, color = openxlsx2::wb_color(hex = "FF000000"))
        }
        if (CS2[[mn]]$Style %in% names(ColList)) {
          clr <- ColList[[CS2[[mn]]$Style]]
          WorkBook <- openxlsx2::wb_add_font(WorkBook, sheetnm, dims,
                                             color = openxlsx2::wb_color(hex = "white"))
        }
        if (CS2[[mn]]$Style %in% names(HAlignList)) {
          hal <- HAlignList[[CS2[[mn]]$Style]]
          WorkBook <- openxlsx2::wb_add_cell_style(WorkBook, sheetnm, dims, horizontal = hal)
        }
        if (CS2[[mn]]$Style %in% names(ContainsList)) {
          cntn <- ContainsList[[CS2[[mn]]$Style]]
          stl <- capture.output(cntn$style)
          deco <- tolower(gsub(".* ", "", gsub(" *$", "", grep("Font decoration", stl, value = TRUE))))
          cntnFrt <- gsub(".* ", "", gsub(" *$", "", grep("Font colour", stl, value = TRUE)))
          cntnBgrd <- gsub(".* ", "", gsub(" *$", "", grep("Cell fill background", stl, value = TRUE)))
          cntnHrz <- gsub(".* ", "", gsub(" *$", "", grep("Cell horz. align", stl, value = TRUE)))
          WorkBook <- openxlsx2::wb_add_dxfs_style(WorkBook, CS2[[mn]]$Style,
                                                   font_color = openxlsx2::wb_color(cntnFrt),
                                                   bg_fill = openxlsx2::wb_color(cntnBgrd),
                                                   text_bold = c("", "single")[("bold" %in% deco)+1],
                                                   text_italic = c("", "italic")[("italic" %in% deco)+1])
          WorkBook <- openxlsx2::wb_add_conditional_formatting(WorkBook, sheetnm, dims,
                                                               type = "containsText", style = CS2[[mn]]$Style, rule = cntn$rule)
        }
      }
    }
  }
  #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
  #
  # Now add the rest of the data
  #  - CC = actual data!
  sheetMtch <- match(sheetnm, wb_get_sheet_names(WorkBook))
  cc <- WorkBook$worksheets[[sheetMtch]]$sheet_data$cc
  cc_3 <- cc[which(cc$row_r == "3"),]
  cc_12 <- cc[which(cc$row_r %in% c("1", "2")),]
  rownames(cc_12) <- NULL
  #uNum <- unique(cc_3$typ[wNum]) # Worked until at least 1.10, doesn't work for 1.14 (I don't know when the break occurred)
  #uTxt <- unique(cc_3$typ[wTxt]) # Worked until at least 1.10, doesn't work for 1.14 (I don't know when the break occurred)
  #wNum2 <- which(cc_3$typ == 2) #
  #wTxt2 <- which(cc_3$typ == 4)
  stopifnot(#(length(uNum) == 1)&&(uNum == "2"), # Worked until at least 1.10, doesn't work for 1.14 (I don't know when the break occurred)
            #(length(uTxt) == 1)&&(uTxt == "4"), # Worked until at least 1.10, doesn't work for 1.14 (I don't know when the break occurred)
            length(unique(c(wNum, wTxt))) == nCol,
            sum(wNum %in% wTxt) == 0) # We only cover types 2 and 4
  tmp <- cc_3
  tmp$v <- tmp$is <- ""
  tmp$c_t[wTxt] <- "inlineStr"
  tmp$c_t[wNum] <- ""
  ###
  cc_Rest <- 1:(nRws*nCol) - 1
  cc_Rest <- cc_Rest %% nCol
  cc_Rest <- cc_Rest + 1
  cc_Rest <- tmp[cc_Rest,]
  a <- 1:nrow(cc_Rest)-1
  a <- a - (a %% nCol)
  a <- a/nCol
  a <- a + 3
  cc_Rest$row_r <- a
  cc_Rest$r <- do.call(paste0, c(cc_Rest[, c("c_r", "row_r")]))
  # Numeric values go into column "v", text goes in "is" with some xml formatting
  lngRws <- ((1:nRws)-1)*nCol
  opt <- getOption("scipen") # To avoid scientific notation when writing as text
  options(scipen = 999)
  # Deal with numbers
  if (length(wNum)) {
    for (i in wNum) {
      rg <- lngRws+i
      cc_Rest$v[rg] <- myData[[i]]
    }
  }
  options(scipen = opt)
  # Deal with text
  if (length(wTxt)) {
    # We need to escape some characters!
    scpChar <- data.frame(origChar = c("\"", "'", "<", ">", "&"),
                          scpdChar = c("&quot;", "&apos;", "&lt;", "&gt;", "&amp;"))
    for (i in wTxt) {
      rg <- lngRws+i
      dat <- strsplit(myData[[i]], "")
      for (j in 1:5) { #j <- 1 #j <- 2 #j <- 3 #j <- 4 #j <- 5
        tst <- lapply(dat, function(x) { which(x == scpChar$origChar[j]) })
        w <- which(vapply(tst, length, 1) > 0)
        if (length(w)) {
          tmpScpd <- vapply(w, function(x) { #x <- w[1]
            y <- dat[[x]]
            y[tst[[x]]] <- scpChar$scpdChar[j]
            y <- paste(y, collapse = "")
          }, "")
          myData[[i]][w] <- tmpScpd
        }
      }
      g <- grepl(" ", myData[[i]])+1
      cc_Rest$is[rg] <- paste0("<is><t", c("", " xml:space=\"preserve\"")[g], ">", myData[[i]], "</t></is>")
    }
  }
  # Deal with infinites
  w <- which(cc_Rest$v %in% c("-Inf", "Inf"))
  if (length(w)) {
    cc_Rest$v[w] <- "#NUM!" 
    cc_Rest$c_t[w] <- "e"
  }
  w <- which((cc_Rest$is %in% c("<is><t>-Inf</t></is>", "<is><t>Inf</t></is>",
                                "<is><t>-Inf</t xml:space=\"preserve\"></is>", "<is><t xml:space=\"preserve\">Inf</t></is>"))|(is.infinite(cc_Rest$is)))
  if (length(w)) {
    cc_Rest$v[w] <- "#NUM!" 
    cc_Rest$c_t[w] <- "e"
  }
  # Deal with NaNs
  w <- which(cc_Rest$v == "NaN")
  if (length(w)) {
    cc_Rest$v[w] <- "#VALUE!" 
    cc_Rest$c_t[w] <- "e"
  }
  w <- which((cc_Rest$is %in% c("<is><t>NaN</t></is>", "<is><t xml:space=\"preserve\">NaN</t></is>"))|(is.nan(cc_Rest$is)))
  if (length(w)) {
    cc_Rest$v[w] <- "#VALUE!" 
    cc_Rest$c_t[w] <- "e"
    cc_Rest$is[w] <- ""
  }
  # Deal with NAs
  w <- which((is.na(cc_Rest$v))|(cc_Rest$v == "NA"))
  if (length(w)) {
    cc_Rest$v[w] <- "#N/A"
    cc_Rest$c_t[w] <- "e"
  }
  w <- which((cc_Rest$is %in% c("<is><t>NA</t></is>", "<is><t xml:space=\"preserve\">NA</t></is>"))|(is.na(cc_Rest$is)))
  if (length(w)) {
    cc_Rest$v[w] <- "#N/A"
    cc_Rest$c_t[w] <- "e"
    cc_Rest$is[w] <- ""
  }
  #
  rownames(cc_Rest) <- NULL
  cc <- rbind(cc_12, cc_Rest)
  rownames(cc) <- as.character(1:nrow(cc))
  #
  # Fix range of conditional formatting
  cf <- names(WorkBook$worksheets[[sheetMtch]]$conditionalFormatting)
  if (length(cf)) {
    cfNms0 <- cfNms <- names(WorkBook$worksheets[[sheetMtch]]$conditionalFormatting)
    w <- which(gsub("[A-Z]", "", cfNms) == "3:3")
    if (length(w)) {
      cfNms <- gsub("3$", nRws+2, cfNms[w])
      names(WorkBook$worksheets[[sheetMtch]]$conditionalFormatting)[w] <- cfNms
    }
    #
    # Now, when I compare a cc table created with the classic method with this one, it seems that one key difference is that
    # in the former conditional formatting styles are only written in the first cell of each column.
    # Not the first cell of a rectangular range, but each first cell of each column affected.
    # I am not sure that this is important to fix, but since I am trying to get as close as possible to the real thing...
    # cfStrt <- gsub(":.*", "", cfNms0)
    # cfEnd0 <- gsub(".*:", "", cfNms0)
    # cfEnd <- gsub(".*:", "", cfNms)
    # cfEndRw <- gsub("[A-Z]+", "", cfEnd)
    # cfRg <- lapply(1:length(cfStrt), function(x) {
    #   cc_3$r[match(cfStrt[x], cc_3$r):match(cfEnd0[x], cc_3$r)]
    # })
    # cfRg2 <- lapply(1:length(cfStrt), function(x) {
    #   rg <- cfRg[[x]]
    #   kol <- gsub("[0-9]+", "", rg)
    #   lapply(1:length(rg), function(y) {
    #     cc_ <- cc[which(cc$c_r == kol[y]),]
    #     rs <- cc_$r[which(cc_$r == rg[y]):which(cc_$row_r == cfEndRw[x])]
    #     rs <- rs[which(rs != rg[y])]
    #   })
    # })
    # cfRg2 <- unlist(cfRg2)
    # cc$c_s[match(cfRg2, cc$r)] <- ""
  }
  #
  WorkBook$worksheets[[sheetMtch]]$sheet_data$cc <- cc
  #
  #  - Row attributes!
  ra <- WorkBook$worksheets[[sheetMtch]]$sheet_data$row_attr
  tmp <- rep(3, nRws-1)
  tmp <- ra[tmp, ]
  tmp$r <- (2:nRws)+2
  ra <- rbind(ra, tmp)
  rownames(ra) <- as.character(ra$r)
  WorkBook$worksheets[[sheetMtch]]$sheet_data$row_attr <- ra
  # Also xml, etc...
  try(WorkBook$worksheets[[sheetMtch]]$sheet_data$cc_out <- NULL, silent = TRUE) # Worked until at least 1.10, doesn't work for 1.14 (I don't know when the break occurred)
  # I don't think this was necessary anyway...
  fullDims <- wb_dims(rows = 1:(nRws+2), cols = 1:ncol(myData))
  tblDims <- wb_dims(rows = 2:(nRws+2), cols = 1:ncol(myData))
  WorkBook$worksheets[[sheetMtch]]$dimension <- paste0("<dimension ref=\"", fullDims, "\"/>")
  WorkBook$tables$tab_ref[match(sheetMtch, WorkBook$tables$tab_sheet)] <- tblDims
  xml <- WorkBook$tables$tab_xml[length(WorkBook$tables$tab_xml)]
  xml <- as.character(xml) # Fix weird character encoding shenanigans (should I?)
  xmlID <- gsub("\".*", "", unlist(strsplit(xml, " id=\""))[2])
  # Code to check that xml
  #writeClipboard(xml)
  #writeClipboard(gsub("/><.*", "/>", xml))
  #writeClipboard(gsub("/><.*", "/>", tstXML))
  #a <- unlist(strsplit(gsub("/><", "/>___CUTHERE___<", xml), "___CUTHERE___"));writeClipboard(a)
  #tst <- Isapply(strsplit(gsub("/><", "/>___CUTHERE___<", tstXML), "___CUTHERE___"), unlist)
  #l <- apply(tst, 2, function(x) { length(unique(x)) })
  #(all but the first should be same)
  #tst2 <- tst[, which(l > 1)]
  #tst2
  #grep("ref=", a, value = TRUE) # Only the range in the first need be encoded
  xml <- gsub("ref=\"[A-Z]+[0-9]+:[A-Z]+[0-9]+\"", paste0("ref=\"", tblDims, "\""), xml)
  pat <- paste0(" id=\"", xmlID,"\" name=\"", tblNm, "\" displayName=\"", tblNm, "\" ")
  xml <- paste0(unlist(strsplit(xml, pat)))
  xml <- paste0(xml[1], " id=\"", xmlID, "\" name=\"", tblNm, "\" displayName=\"", tblNm, "\" ", xml[2])
  WorkBook$tables$tab_xml[length(WorkBook$tables$tab_xml)] <- xml
  #
  #writeClipboard(xml)
  #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
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
  filt <- c("General Peptides information",
            unique(ColumnsTbl$Class[match(intCols[[rev(names(intCols))[1]]], ColumnsTbl$Col)]),
            "QC filters")
  if (MakeRatios) {
    filt <- c(filt,
              unique(ColumnsTbl$Class[match(ratCols[[rev(names(ratCols))[1]]], ColumnsTbl$Col)]),
              "Regulated")
  }
  tst <- tst[which(!tst$Group %in% filt),]
  if (nrow(tst)) {
    tst$Min <- sapply(tst$Col, min)
    tst$Max <- sapply(tst$Col, max)
    tst$Max <- tst$Max+1 # For openxlsx2
    tst <- tst[order(tst$Min, decreasing = FALSE),]
    for (i in 1:nrow(tst)) {
      WorkBook <- openxlsx2::wb_group_cols(WorkBook,
                                           sheetnm,
                                           tst$Min[i]:tst$Max[i],
                                           tst$Hide[i])
    }
  }
  #
  if (tblMode == "PG") {
    # - Coverage columns:
    #   - Xml encoded coverage
    kols <- xmlCovCol
    kols <- kols[which(kols %in% xlTabs[[sheetnm]])]
    if (length(kols)) {
      for (kol in kols) {
        dims <- openxlsx2::wb_dims(tblRws[1]:tblRws[2], match(kol, xlTabs[[sheetnm]]))
        kol2 <- gsub(" - ", " ", cleanNms(kol))
        WorkBook <- openxlsx2::wb_add_data(WorkBook, sheetnm, myData[, kol2], dims,
                                           col_names = FALSE)
      }
    }
    #   - Data bar summaries
    kols <- grep("\\[%\\]", covcol, value = TRUE)
    kols <- kols[which(kols %in% xlTabs[[sheetnm]])]
    if (length(kols)) {
      for (kol in kols) {
        dims <- openxlsx2::wb_dims(tblRws[1]:tblRws[2], match(kol,
                                                              xlTabs[[sheetnm]]))
        WorkBook <- openxlsx2::wb_add_conditional_formatting(WorkBook, sheetnm,
                                                             dims, type = "dataBar",
                                                             style = c("FF63C384",
                                                                       "0063C384"),
                                                             rule = c(0, 100))
      }
    }
  }
  # Below: this is done now with the "In list" column, so as not to allow unformating of columns
  # Unfortunately, openxlsx2 does not allow for easy font color changes because of how it deals with styles.
  #
  # if (prot.list.Cond) {
  #   kol <- c(paste0("Leading protein", c("", " ID"), "s"), "Protein")[match(tblMode, c("pep", "PG", "SAINTexpress"))]
  #   #print(kol %in% colnames(myData))
  #   if (tblMode == "SAINTexpress") { w <- which(myData[[kol]] %in% prot.list)+2 } else {
  #     w <- proteoCraft::grsep2(prot.list, myData[[kol]])+2  
  #   }
  #   if (length(w)) {
  #     dims <- openxlsx2::wb_dims(match("Protein list", styleNms), 1)
  #     stl <- openxlsx2::wb_get_cell_style(WorkBook, "tmp", dims)
  #     if (stl != "") {
  #       for (i in w) {
  #         dims <- openxlsx2::wb_dims(i, colRg)
  #         WorkBook <- openxlsx2::wb_set_cell_style(WorkBook, sheetnm, dims, stl)
  #         if (tblMode == "PG") {
  #           WorkBook <- openxlsx2::wb_set_row_heights(WorkBook, sheetnm, dims, 50)
  #         }
  #       }
  #     }
  #   }
  # }
  #
  #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
}
