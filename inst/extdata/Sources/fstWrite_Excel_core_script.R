
# Filter for cases where number of characters is more than the Excel limit:
## Should only ever be character columns (lists of spectra IDs, sequences...)
w <- which(sapply(colnames(tempData), \(x) { class(tempData[[x]]) }) == "character")
if (length(w)) {
  for (i in w) { #i <- w[1L]
    w1 <- which(nchar(tempData[[colnames(tempData)[i]]]) > ExcelMax)
    if (length(w1)) {
      tempData[[colnames(tempData)[i]]][w1] <- paste0(substr(tempData[[colnames(tempData)[i]]][w1], 1L, ExcelMax-3L), "...")
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
  tmp <- setdiff(tmp, covcol)
  xlTabs[[TbNm]] <- tmp
}
# Write tab, apply formatting etc...
for (l in LS) { if (!exists(l)) { assign(l, list()) } }
sheetnmsA <- names(xlTabs)
sheetnmsB <- setdiff(sheetnmsA, c("Description", "Quality control"))
cat(paste0("    - building ", tblMode2, " tab...\n"))
if (tblMode == "pep") {
  nms <- sub(" ///NL///$", "", intNms(names(intRf), TRUE, "pep"))
  if (MakeRatios) {
    nms <- c(nms, ratNms(names(ratRf), TRUE))
  }
  datCol <- ColumnsTbl$Col[unique(which(ColumnsTbl$Class %in% nms))]
}
if (tblMode == "PG") {
  nms <- sub(" ///NL///$", "", intNms(names(intColsTbl), TRUE))
  if ((MakeRatios)&&(!is.null(ratRf))) {
    nms <- c(nms, ratNms(names(ratColsTbl), TRUE))
  }
  datCol <- ColumnsTbl$Col[unique(c(which(ColumnsTbl$Class %in% nms),
                                    grep("est\\. copies/cell", ColumnsTbl$Class)))]
}
#
if (tblMode == "SAINTexpress") {
  datCol <- ColumnsTbl$Col[which(ColumnsTbl$Class == "log2(rat.), avg.")]
}
stopifnot(length(datCol) > 0L)
wrtHeader <- TRUE
for (sheetnm in sheetnmsB) { #sheetnm <- sheetnmsB[1L] #sheetnm <- sheetnmsB[2L]
  # Very important early filter, to avoid discrepancies between header and columns!
  # Should be done at least before the header is written
  xlTabs[[sheetnm]] <- intersect(xlTabs[[sheetnm]], colnames(tempData))
  xlTabs[[sheetnm]] <- intersect(xlTabs[[sheetnm]], ColumnsTbl$Col) # The order matters!
  #
  if (!sheetnm %in% names(data_filt)) { data_filt[[sheetnm]] <- 1L:nrow(tempData) }
  if (!sheetnm %in% names(data_order)) { data_order[[sheetnm]] <- 1L:length(data_filt[[sheetnm]]) }
  #WorkBook <- wb_add_worksheet(WorkBook, sheetnm, grid_lines = FALSE)
  #saveFun(WorkBook, file = "WorkBook_bckp.RData")
  #
  # For testing:
  # data_filt[[sheetnm]] <- data_order[[sheetnm]] <- 1L:10L
  #loadFun("WorkBook_bckp.RData")
  wtst <- data_filt[[sheetnm]]
  sheetrows <- length(wtst)+2L # delete me? not used anywhere I think
  ord <- data_order[[sheetnm]]
  colRg <- 1L:length(xlTabs[[sheetnm]])
  #
  # Write the data
  # Styles list which applies to this tab and is in the right order (important for group processing/grouping!)
  ColumnsTbl2 <- data.frame(Col = xlTabs[[sheetnm]], Grp = "", Class = "")
  w <- which(ColumnsTbl2$Col %in% ColumnsTbl$Col)
  ColumnsTbl2[w, c("Grp", "Class")] <- ColumnsTbl[match(xlTabs[[sheetnm]][w], ColumnsTbl$Col), c("Grp", "Class")] 
  # Test for series of consecutive columns from the same class
  tstKol <- lapply(1L:nrow(ColumnsTbl2), \(x) { #x <- 1L
    res <- x
    w <- which(ColumnsTbl2$Class == ColumnsTbl2$Class[x])
    y <- x-1L
    while (y %in% w) {
      res <- c(y, res)
      y <- y-1L
    }
    y <- x+1L
    while (y %in% w) {
      res <- c(res, y)
      y <- y+1L
    }
    return(paste(c(min(res), max(res)), collapse = ";"))
  })
  tstKol <- as.data.frame(t(sapply(strsplit(unique(unlist(tstKol)), ";"), as.numeric)))
  colnames(tstKol) <- c("First", "Last")
  tstKol$Group <- ColumnsTbl2$Class[tstKol$First]
  tstKol <- tstKol[which(tstKol$Group != ""), , drop = FALSE]
  whLevel[[sheetnm]] <- rep(1L, length(xlTabs[[sheetnm]])) # For now all columns are to write in the 1st header row; this may change
  #
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
  tblRws <- c(1L, nRws)
  hdRg <- 1L:2L
  tblRws <- tblRws + 2L
  # Round PEPs, we don't need all of those zeros!!!
  if ("PEP" %in% colnames(myData)) { myData$PEP <- round(myData$PEP, 6L) }
  #
  # NB: Excel can take 1048576 rows max!
  # At some point we may need a solution for if we ever have more rows!
  if (nRws > 1048576L - 2L) {
    stop("Unhandled extreme case: more rows than the maximum allowed limit (= 1048576) in Excel!")
    # We'll cross that bridge when we hit it...
  }
  #
  if (sheetnm %in% wb_get_sheet_names(WorkBook)) { WorkBook <- wb_remove_worksheet(WorkBook, sheetnm) }
  WorkBook <- wb_add_worksheet(WorkBook,
                               sheetnm,
                               grid_lines = FALSE)
  #WorkBook <- openxlsx2::wb_remove_worksheet(WorkBook, sheetnm);WorkBook <- openxlsx2::wb_add_worksheet(WorkBook, sheetnm, grid_lines = FALSE)
  #
  #WorkBook <- openxlsx2::wb_load(fl);WorkBook <- openxlsx2::wb_add_worksheet(WorkBook, sheetnm, grid_lines = FALSE)
  #
  tblNm <- tolower(gsub(" |\\$|\\.|-", "_", sheetnm))
  dims <- openxlsx2::wb_dims(rows = hdRg[2L], cols = 1L)
  # We will not write the first row of the data, but a nicely behaved dummy row without any NAs, NaNs or Inf...
  dummyData <- myData[1L, , drop = FALSE]
  # Check storage type before writing data to avoid issues - e.g. writing numerics as text
  for (k in colnames(myData)) {
    if (is.factor(myData[[k]])) { myData[[k]] <- as.character(myData[[k]]) }
    if (is.list(myData[[k]])) { myData[[k]] <- vapply(myData[[k]], paste, "", collapse = ";") }
    if (is.matrix(myData[[k]]) || is.data.frame(myData[[k]])) { myData[[k]] <- drop(myData[[k]]) }
  }
  typeTst <- setNames(vapply(colnames(myData), \(x) {
    #inherits(myData[[x]], c("numeric", "integer"))
    is.numeric(myData[[x]])
  }, TRUE), colnames(myData))
  wNum <- setNames(which(typeTst), NULL)
  wTxt <- setNames(which(!typeTst), NULL)
  dummyData[wNum] <- 0L
  dummyData[wTxt] <- "Hello world!"
  if (tblNm %in% WorkBook$tables$tab_name) { # Let's see if that works...
    WorkBook$tables <- WorkBook$tables[which(WorkBook$tables$tab_name != tblNm),]
  }
  WorkBook <- openxlsx2::wb_add_data_table(WorkBook,
                                           sheetnm,
                                           dummyData,
                                           dims,
                                           table_name = tblNm,
                                           table_style = "TableStyleMedium2",
                                           banded_rows = TRUE,
                                           banded_cols = FALSE)
  WorkBook <- openxlsx2::wb_add_cell_style(WorkBook,
                                           sheetnm,
                                           wb_dims(rows = hdRg[2L], cols = 1L:ncol(dummyData)),
                                           wrap_text = TRUE)
  #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
  #
  # Attempt to override the banding in the data columns,
  # because #N/A values are not affected by conditional formatting
  # and seeing them in alternate colours looks... weird.
  tst <- data.frame(IsQuant = xlTabs[[sheetnm]] %in% datCol)
  tst$IsChange <- c(FALSE, sapply(2L:nrow(tst), \(x) {
    as.numeric(tst$IsQuant[x] != tst$IsQuant[x-1L])
  }))
  tst$Block <- cumsum(tst$IsChange)
  blocks <- aggregate(1L:nrow(tst), list(tst$Block), unique)
  colnames(blocks) <- c("Block", "Columns")
  blocks$IsQuant <- tst$IsQuant[match(blocks$Block, tst$Block)]
  w <- which(blocks$IsQuant)
  for (i in w) {
    dims <- openxlsx2::wb_dims(rows = tblRws[1L]#:tblRws[2L]
                               , cols = blocks$Columns[[i]])
    WorkBook <- openxlsx2::wb_add_fill(WorkBook, sheetnm, dims, color = openxlsx2::wb_color("yellow"))
  }
  #
  # Write top level header
  if (wrtHeader) {
    for (i in 1L:nrow(tstKol)) {
      dims <- openxlsx2::wb_dims(hdRg[1L], tstKol$First[i])
      WorkBook <- openxlsx2::wb_add_data(WorkBook, sheetnm, tstKol$Group[i], dims) # Write 1st header level
      if (tstKol$Last[i] > tstKol$First[i]) { # Merge if 1st header level is same for several columns
        dims <- openxlsx2::wb_dims(hdRg[1L], tstKol$First[i]:tstKol$Last[i])
        WorkBook <- openxlsx2::wb_merge_cells(WorkBook, sheetnm, dims)
      }
    }
  }
  #openxlsx2::wb_save(WorkBook, repFl);xl_open(repFl)
  #
  # Style the table
  # - vertical align to center
  dims <- openxlsx2::wb_dims(1L#:tblRws[2L]
                             , colRg)
  WorkBook <- openxlsx2::wb_add_cell_style(WorkBook, sheetnm, dims, vertical = "center")
  #a1 <- WorkBook$worksheets[[3L]]$sheet_data$cc; a2 <- WorkBook$worksheets[[3L]]$sheet_data$cc; length(unique(a1$c_s)); length(unique(a2$c_s))
  colWdth <- setNames(vapply(strsplit(colnames(myData), "///((NL)|(VS))///"), \(x) { #x <- strsplit(colnames(myData)[35L], "///((NL)|(VS))///")
    x <- unlist(x)
    x <- max(vapply(strsplit(x, ""), \(y) {
      length(y) + sum(y %in% LETTERS)
    }, 1L)) + 4L
    max(c(5L, min(c(35L, x))))
  }, 1L), colnames(myData))
  colHgth <- (max(vapply(strsplit(colnames(myData), "///((NL)|(VS))///"), length, 1L))+0.5)*15.75
  WorkBook <- openxlsx2::wb_set_col_widths(WorkBook, sheetnm, 1L:length(colWdth), colWdth)
  #
  w <- which(xlTabs[[sheetnm]] %in% xmlCovCol)
  if (length(w)) {
    dims <- openxlsx2::wb_dims(tblRws[1L], w)
    WorkBook <- openxlsx2::wb_add_cell_style(WorkBook, sheetnm,
                                             dims, vertical = "top")
  }
  #
  if (wrtHeader) {
    # Style the header
    # - Header borders
    WorkBook <- openxlsx2::wb_add_border(WorkBook, sheetnm,
                                         openxlsx2::wb_dims(hdRg, colRg))
    #
    dims <- openxlsx2::wb_dims(match("Header 1", styleNms), 1L)
    stl <- openxlsx2::wb_get_cell_style(WorkBook, "tmp", dims)
    if (stl != "") {
      dims <- openxlsx2::wb_dims(hdRg[2L], colRg)
      WorkBook <- openxlsx2::wb_set_cell_style(WorkBook, sheetnm, dims, stl)
    }
    dims <- openxlsx2::wb_dims(match("Header 2", styleNms), 1L)
    stl <- openxlsx2::wb_get_cell_style(WorkBook, "tmp", dims)
    if (stl != "") {
      dims <- openxlsx2::wb_dims(hdRg[1L], colRg)
      WorkBook <- openxlsx2::wb_set_cell_style(WorkBook, sheetnm, dims, stl)
    }
    #
    # Row heights:
    WorkBook <- openxlsx2::wb_set_row_heights(WorkBook,
                                              sheetnm,
                                              hdRg,
                                              c(max(c(30L, min(45L, c(ceiling(max(nchar(tstKol$Group))/20)*10L)))),
                                                colHgth))
  }
  # Freeze panes
  m <- match(CoreCol, xlTabs[[sheetnm]])
  WorkBook <- openxlsx2::wb_freeze_pane(WorkBook,
                                        sheetnm,
                                        first_active_row = c(1L, 3L)[wrtHeader+1L],
                                        firstActiveCol = max(m, na.rm = TRUE)+1L)
  #
  #openxlsx2::wb_save(WorkBook, repFl);xl_open(repFl)
  #
  # Conditional formatting
  w1 <- which(!is.na(vapply(xlTabs[[sheetnm]], \(x) { #x <- xlTabs[[sheetnm]][1L]
    match(ColumnsTbl$Grp[match(x, ColumnsTbl$Col)], styleNms)
  }, 1L)))
  w2 <- colRg
  w2 <- setdiff(w2, w1)
  #View(ColumnsTbl[w1,])
  #View(ColumnsTbl[w2,])
  if (length(w1)) {
    CS1 <- setNames(lapply(xlTabs[[sheetnm]][w1], \(x) { #x <- xlTabs[[sheetnm]][w1][1L]
      m <- match(x, ColumnsTbl$Col)
      res <- list(Group = ColumnsTbl$Grp[m],
                  Class = ColumnsTbl$Class[m])
      res$tmpCell <- openxlsx2::wb_dims(match(res$Group, styleNms), 1L)
      res$Style <- openxlsx2::wb_get_cell_style(WorkBook, "tmp", res$tmpCell)
      return(res)
    }), xlTabs[[sheetnm]][w1])
    CS1tst <- vapply(CS1, paste, "", collapse = "-")
    CS1tst <- aggregate(1L:length(CS1), list(CS1tst), \(x) {
      c(min(x), max(x))
    })
    CS1tst[, c("Min", "Max")] <- as.data.frame(CS1tst$x)
    CS1tst$x <- NULL
    CS1tst <- CS1tst[order(CS1tst$Min),]
    CS1tst$Ranges <- apply(CS1tst[, c("Min", "Max")], 1L, \(x) {
      x <- w1[x[[1L]]:x[[2L]]]
      split(x, cumsum(c(TRUE, diff(x) != 1L)))
    })
    for (rw in 1L:nrow(CS1tst)) {
      mn <- CS1tst$Min[rw]
      #mx <- CS1tst$Max[rw] # not used, since it's used to access CS1 and the values for mn and mx are the same
      colRngs <- CS1tst$Ranges[[rw]]
      for (colRng in colRngs) { #colRng <- colRngs[[1L]]
        dims <- openxlsx2::wb_dims(tblRws[1L]#:tblRws[2L]
                                   , colRng)
        WorkBook <- openxlsx2::wb_set_cell_style(WorkBook, sheetnm, dims, CS1[[mn]]$Style)
      }
    }
  }
  if (length(w2)) {
    CS2 <- setNames(lapply(xlTabs[[sheetnm]][w2], \(x) { #x <- xlTabs[[sheetnm]][w2][1L] #x <- xlTabs[[sheetnm]][w2][5L]
      m <- match(x, ColumnsTbl$Col)
      res <- list(Group = ColumnsTbl$Grp[m],
                  Class = ColumnsTbl$Class[m])
      res$Style <- Styles[[res$Group]]
      return(res)
    }), xlTabs[[sheetnm]][w2])
    CS2tst <- vapply(CS2, paste, "", collapse = "-")
    CS2tst <- aggregate(1L:length(CS2), list(CS2tst), \(x) {
      c(min(x), max(x))
    })
    CS2tst[, c("Min", "Max")] <- as.data.frame(CS2tst$x)
    CS2tst$x <- NULL
    CS2tst <- CS2tst[order(CS2tst$Min),]
    CS2tst$Ranges <- apply(CS2tst[, c("Min", "Max")], 1L, \(x) {
      x <- w2[x[[1L]]:x[[2L]]]
      split(x, cumsum(c(TRUE, diff(x) != 1L)))
    })
    #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
    #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
    # for (rw in 1L:nrow(CS2tst)) { #rw <- 1L
    #   print(CS2tst$Group.1[[rw]])
    # }
    for (rw in 1L:nrow(CS2tst)) { #rw <- 1L #rw <- 3L #rw <- 4L
      mn <- CS2tst$Min[rw]
      #mx <- CS2tst$Max[rw] # not used, since it's used to access CS2 and the values for mn and mx are the same
      if ("Style" %in% names(CS2[[mn]])) {
        colRngs <- CS2tst$Ranges[[rw]]
        for (colRng in colRngs) { #colRng <- colRngs[[1L]]
          dims <- openxlsx2::wb_dims(tblRws[1L]#:tblRws[2L]
                                     , colRng)
          dims2 <- openxlsx2::wb_dims(tblRws[1L]:(tblRws[1L]+1L), colRng)
          stlNm <- CS2[[mn]]$Style
          if (stlNm %in% names(DecoList)) {
            deco <- DecoList[[stlNm]]
            WorkBook <- openxlsx2::wb_add_font(WorkBook, sheetnm, dims,
                                               italic = c("", "true")[("italic" %in% deco)+1L],
                                               bold = c("", "true")[("bold" %in% deco)+1L])
          }
          if (stlNm %in% names(SignifList)) {
            signi <- SignifList[[stlNm]]
            WorkBook <- openxlsx2::wb_add_numfmt(WorkBook, sheetnm, dims,
                                                 paste(c("0.", rep("0", signi)), collapse = ""))
          }
          if (stlNm %in% names(ColScaleList)) {
            condFrmt <- ColScaleList[[stlNm]]
            WorkBook <- openxlsx2::wb_add_conditional_formatting(WorkBook, sheetnm, dims,
                                                                 type = "colorScale", style = condFrmt)
            #WorkBook <- openxlsx2::wb_add_font(WorkBook, sheetnm, dims, color = openxlsx2::wb_color(hex = "FF000000"))
          }
          if (stlNm %in% names(ColList)) {
            clr <- ColList[[stlNm]]
            WorkBook <- openxlsx2::wb_add_font(WorkBook, sheetnm, dims,
                                               color = openxlsx2::wb_color(hex = "white"))
          }
          if (stlNm %in% names(HAlignList)) {
            hal <- HAlignList[[stlNm]]
            WorkBook <- openxlsx2::wb_add_cell_style(WorkBook, sheetnm, dims, horizontal = hal)
          }
          if (stlNm %in% names(ContainsList)) {
            cntn <- ContainsList[[stlNm]]
            stl <- capture.output(cntn$style)
            deco <- tolower(gsub(".* ", "", gsub(" *$", "", grep("Font decoration", stl, value = TRUE))))
            cntnFrt <- gsub(".* ", "", gsub(" *$", "", grep("Font colour", stl, value = TRUE)))
            cntnBgrd <- gsub(".* ", "", gsub(" *$", "", grep("Cell fill background", stl, value = TRUE)))
            cntnHrz <- gsub(".* ", "", gsub(" *$", "", grep("Cell horz. align", stl, value = TRUE)))
            if (!stlNm %in% dxfs_Style_nms) {
              WorkBook <- openxlsx2::wb_add_dxfs_style(WorkBook, stlNm,
                                                       font_color = openxlsx2::wb_color(cntnFrt),
                                                       bg_fill = openxlsx2::wb_color(cntnBgrd),
                                                       text_bold = c("", "single")[("bold" %in% deco)+1L],
                                                       text_italic = c("", "italic")[("italic" %in% deco)+1L])
              dxfs_Style_nms <- union(dxfs_Style_nms, stlNm)
            }
            WorkBook <- openxlsx2::wb_add_conditional_formatting(WorkBook, sheetnm, dims,
                                                                 type = "containsText", style = stlNm, rule = cntn$rule)
          }
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
    sum(wNum %in% wTxt) == 0L) # We only cover types 2 and 4
  tmp <- cc_3
  tmp$v <- tmp$is <- ""
  tmp$c_t[wTxt] <- "inlineStr"
  tmp$c_t[wNum] <- ""
  ###
  cc_Rest <- 1L:(nRws*nCol) - 1L
  cc_Rest <- cc_Rest %% nCol
  cc_Rest <- cc_Rest + 1L
  cc_Rest <- tmp[cc_Rest,]
  a <- 1L:nrow(cc_Rest) - 1L
  a <- a - (a %% nCol)
  a <- a/nCol
  a <- a + 3L
  cc_Rest$row_r <- a
  cc_Rest$r <- do.call(paste0, c(cc_Rest[, c("c_r", "row_r")]))
  # Numeric values go into column "v", text goes in "is" with some xml formatting
  lngRws <- ((1L:nRws)-1L)*nCol
  opt <- getOption("scipen") # To avoid scientific notation when writing as text
  options(scipen = 999L)
  # Deal with numbers
  if (length(wNum)) {
    for (i in wNum) {
      rg <- lngRws+i
      cc_Rest$v[rg] <- myData[[i]]
    }
  }
  options(scipen = opt)
  # Deal with text
  a <- if (MakeRatios) {
    if (tblMode2 == "SAINTexpress") {
      KolEdit(xmlCovCol, ratColsTbl)
    } else {
      KolEdit(xmlCovCol, intColsTbl, ratColsTbl)
    }
  } else { KolEdit(xmlCovCol, intColsTbl) }
  wTxt2 <- wTxt[which(!colnames(myData[wTxt]) %in% a)]
  if (length(wTxt2)) {
    rg <- unlist(lapply(wTxt2, \(i) { lngRws+i }))
    dat <- unlist(lapply(wTxt2, \(i) { myData[[i]] }))
    g <- grepl(" ", dat)+1L
    dat <- htmltools::htmlEscape(dat)
    cc_Rest$is[rg] <- paste0("<is><t", c("", " xml:space=\"preserve\"")[g], ">", dat, "</t></is>")
  }
  # Deal with infinites
  w <- which(cc_Rest$v %in% c("-Inf", "Inf"))
  if (length(w)) {
    cc_Rest$v[w] <- "#NUM!" 
    cc_Rest$c_t[w] <- "e"
  }
  w <- which((cc_Rest$is %in% c("<is><t>-Inf</t></is>", "<is><t>Inf</t></is>",
                                "<is><t>-Inf</t xml:space=\"preserve\"></is>",
                                "<is><t xml:space=\"preserve\">Inf</t></is>"))
             |(is.infinite(cc_Rest$is)))
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
  w <- which((cc_Rest$is %in% c("<is><t>NaN</t></is>",
                                "<is><t xml:space=\"preserve\">NaN</t></is>"))
             |(is.nan(cc_Rest$is)))
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
  w <- which((cc_Rest$is %in% c("<is><t>NA</t></is>",
                                "<is><t xml:space=\"preserve\">NA</t></is>"))
             |(is.na(cc_Rest$is)))
  if (length(w)) {
    cc_Rest$v[w] <- "#N/A"
    cc_Rest$c_t[w] <- "e"
    cc_Rest$is[w] <- ""
  }
  #
  rownames(cc_Rest) <- NULL
  cc <- rbind(cc_12, cc_Rest)
  rownames(cc) <- as.character(1L:nrow(cc))
  #
  # Fix range of conditional formatting
  cf <- WorkBook$worksheets[[sheetMtch]]$conditionalFormatting
  # - method 1: works in some versions of openxlsx2
  if (is.data.frame(cf) && nrow(cf)) {
    g <- grep("[A-Z]+3$", cf$sqref)
    cf$sqref[g] <- gsub("3$", as.character(nRws+2L), cf$sqref[g])
    WorkBook$worksheets[[sheetMtch]]$conditionalFormatting <- cf
  }
  # - method 2: works in some others
  if (is.character(cf) && length(cf)) {
    cfNms <- names(WorkBook$worksheets[[sheetMtch]]$conditionalFormatting)
    w <- which(gsub("[A-Z]", "", cfNms) == "3:3")
    if (length(w)) {
      cfNms[w] <- gsub("3$", nRws+2L, cfNms[w])
      names(WorkBook$worksheets[[sheetMtch]]$conditionalFormatting)[w] <- cfNms
    }
    
    WorkBook$worksheets[[sheetMtch]]$conditionalFormatting[w] <- gsub("type=\"percentile\" val=\"50\"",
                                                                      "type=\"num\" val=\"0\"",
                                                                      WorkBook$worksheets[[sheetMtch]]$conditionalFormatting[w],
                                                                      fixed = TRUE)
    
  }
  #
  WorkBook$worksheets[[sheetMtch]]$sheet_data$cc <- cc
  #
  #  - Row attributes!
  ra <- WorkBook$worksheets[[sheetMtch]]$sheet_data$row_attr
  tmp <- rep(3L, nRws-1L)
  tmp <- ra[tmp, ]
  tmp$r <- (2L:nRws)+2L
  ra <- rbind(ra, tmp)
  rownames(ra) <- as.character(ra$r)
  WorkBook$worksheets[[sheetMtch]]$sheet_data$row_attr <- ra
  # Also xml, etc...
  try(WorkBook$worksheets[[sheetMtch]]$sheet_data$cc_out <- NULL, silent = TRUE) # Worked until at least 1.10, doesn't work for 1.14 (I don't know when the break occurred)
  # I don't think this was necessary anyway...
  fullDims <- wb_dims(rows = 1L:(nRws+2L), cols = 1L:ncol(myData))
  tblDims <- wb_dims(rows = 2L:(nRws+2L), cols = 1L:ncol(myData))
  WorkBook$worksheets[[sheetMtch]]$dimension <- paste0("<dimension ref=\"", fullDims, "\"/>")
  WorkBook$tables$tab_ref[match(sheetMtch, WorkBook$tables$tab_sheet)] <- tblDims
  xml <- WorkBook$tables$tab_xml[length(WorkBook$tables$tab_xml)]
  xml <- as.character(xml) # Fix weird character encoding shenanigans (should I?)
  xmlID <- gsub("\".*", "", unlist(strsplit(xml, " id=\""))[2L])
  # Code to check that xml
  #writeClipboard(xml)
  #writeClipboard(gsub("/><.*", "/>", xml))
  #writeClipboard(gsub("/><.*", "/>", tstXML))
  #a <- unlist(strsplit(gsub("/><", "/>___CUTHERE___<", xml), "___CUTHERE___"));writeClipboard(a)
  #tst <- Isapply(strsplit(gsub("/><", "/>___CUTHERE___<", tstXML), "___CUTHERE___"), unlist)
  #l <- apply(tst, 2L, \(x) { length(unique(x)) })
  #(all but the first should be same)
  #tst2 <- tst[, which(l > 1L)]
  #tst2
  #grep("ref=", a, value = TRUE) # Only the range in the first need be encoded
  xml <- gsub("ref=\"[A-Z]+[0-9]+:[A-Z]+[0-9]+\"", paste0("ref=\"", tblDims, "\""), xml)
  pat <- paste0(" id=\"", xmlID,"\" name=\"", tblNm, "\" displayName=\"", tblNm, "\" ")
  xml <- paste0(unlist(strsplit(xml, pat)))
  xml <- paste0(xml[1L], " id=\"", xmlID, "\" name=\"", tblNm, "\" displayName=\"", tblNm, "\" ", xml[2L])
  WorkBook$tables$tab_xml[length(WorkBook$tables$tab_xml)] <- xml
  #
  #writeClipboard(xml)
  #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
  #
  # - Groupings
  tst <- aggregate(ColumnsTbl$Hide, list(ColumnsTbl$Class), unique)
  stopifnot(is.logical(tst$x))
  colnames(tst) <- c("Group", "Hide")
  tst$Col <- lapply(tst$Group, \(x) {
    m <- which(ColumnsTbl$Class == x)
    return(which(xlTabs[[sheetnm]] %in% ColumnsTbl$Col[m]))
  })
  tst <- tst[which(lengths(tst$Col) > 0L),]
  filt <- c("General Peptides information",
            unique(ColumnsTbl$Class[match(intCols[[rev(names(intCols))[1L]]], ColumnsTbl$Col)]),
            "QC filters")
  if ((MakeRatios)&&(!is.null(ratRf))) {
    filt <- c(filt,
              unique(ColumnsTbl$Class[match(ratCols[[rev(names(ratCols))[1L]]], ColumnsTbl$Col)]),
              "Regulated")
  }
  tst <- tst[which(!tst$Group %in% filt),]
  if (nrow(tst)) {
    tst$Min <- sapply(tst$Col, min)
    tst$Max <- sapply(tst$Col, max)
    tst$Max <- tst$Max+1L # For openxlsx2
    tst <- tst[order(tst$Min, decreasing = FALSE),]
    for (i in 1L:nrow(tst)) {
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
    kols <- intersect(xmlCovCol, xlTabs[[sheetnm]])
    if (length(kols)) {
      for (kol in kols) {
        dims <- openxlsx2::wb_dims(tblRws[1L]:tblRws[2L], match(kol, xlTabs[[sheetnm]]))
        kol2 <- gsub(" - ", " ", cleanNms(kol))
        WorkBook <- openxlsx2::wb_add_data(WorkBook, sheetnm, myData[, kol2], dims,
                                           col_names = FALSE)
      }
    }
    #   - Data bar summaries
    kols <- intersect(grep("\\[%\\]", covcol, value = TRUE), xlTabs[[sheetnm]])
    if (length(kols)) {
      for (kol in kols) {
        dims <- openxlsx2::wb_dims(tblRws[1L]:tblRws[2L], match(kol,
                                                                xlTabs[[sheetnm]]))
        WorkBook <- openxlsx2::wb_add_conditional_formatting(WorkBook, sheetnm,
                                                             dims, type = "dataBar",
                                                             style = c("FF63C384",
                                                                       "0063C384"),
                                                             rule = c(0L, 100L))
      }
    }
  }
  # Below: this is done now with the "In list" column, so as not to allow unformating of columns
  # Unfortunately, openxlsx2 does not allow for easy font color changes because of how it deals with styles.
  #
  # if (prot.list.Cond) {
  #   kol <- c(paste0("Leading protein", c("", " ID"), "s"), "Protein")[match(tblMode, c("pep", "PG", "SAINTexpress"))]
  #   #print(kol %in% colnames(myData))
  #   w <- if (tblMode == "SAINTexpress") { which(myData[[kol]] %in% prot.list)+2L } else {
  #     grsep2(prot.list, myData[[kol]])+2L
  #   }
  #   if (length(w)) {
  #     dims <- openxlsx2::wb_dims(match("Protein list", styleNms), 1L)
  #     stl <- openxlsx2::wb_get_cell_style(WorkBook, "tmp", dims)
  #     if (stl != "") {
  #       for (i in w) {
  #         dims <- openxlsx2::wb_dims(i, colRg)
  #         WorkBook <- openxlsx2::wb_set_cell_style(WorkBook, sheetnm, dims, stl)
  #         if (tblMode == "PG") {
  #           WorkBook <- openxlsx2::wb_set_row_heights(WorkBook, sheetnm, dims, 50L)
  #         }
  #       }
  #     }
  #   }
  # }
  #
  #wb_save(WorkBook, paste0(wd, "/tst.xlsx"));xl_open(paste0(wd, "/tst.xlsx"))
}
