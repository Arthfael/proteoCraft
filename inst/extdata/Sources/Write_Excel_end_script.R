# Special quality control tab
sheetnm <- "Quality control"
sheetnmsA <- unique(c(sheetnmsA, sheetnm))
if (sheetnm %in% wb_get_sheet_names(WorkBook)) { WorkBook <- wb_remove_worksheet(WorkBook, sheetnm) }
WorkBook <- wb_add_worksheet(WorkBook, sheetnm, grid_lines = FALSE)
dimsA <- wb_dims(2L, 2L)
WorkBook <- wb_add_data(WorkBook, sheetnm, "Summary plots and tables", dimsA)
stl <- wb_get_cell_style(WorkBook, "tmp", wb_dims(match("Header 1", styleNms), 1L))
if (stl != "") { WorkBook <- wb_set_cell_style(WorkBook, sheetnm, dimsA, stl) }
WorkBook <- wb_add_cell_style(WorkBook, sheetnm, dimsA, wrap_text = FALSE)
WorkBook <- wb_set_col_widths(WorkBook, sheetnm, 1L, 5L)
tmp <- ceiling(nchar("Summary plots and tables")/8.43)*8.43+5
WorkBook <- wb_set_col_widths(WorkBook, sheetnm, 2L, tmp)
WorkBook <- wb_set_row_heights(WorkBook, sheetnm, 1L, 15L)
colWdths <- rep(8.43, 1000L) # Should be enough, what say you?
colWdths[1L] <- 5L
colWdths[2L] <- tmp
fls <- c(list.files(paste0(wd, "/Summary plots"), "\\.svg$", full.names = TRUE, recursive = TRUE),
         grep(" VS ", list.files(paste0(wd, "/Workflow control"), "\\.svg$", full.names = TRUE, recursive = TRUE), value = TRUE, invert = TRUE))
nImgs <- length(fls)
o <- 24L
if (nImgs) {
  fls2 <- sub("\\.svg$", ".png", fls)
  invisible(lapply(1L:nImgs, \(i) { rsvg::rsvg_png(fls[i], file = fls2[i]) }))
  flsTbl <- data.frame(File = fls2,
                       x = ceiling(1L:nImgs/2),
                       y = (((1L:nImgs)+1L) %% 2L) + 1L,
                       Width = 0L,
                       Height = 0L)
  flsTblList <- sapply(1L:nImgs, list)
  for (i in 1L:nImgs) {
    flsTblList[[i]] <- png::readPNG(flsTbl$File[i])
    flsTbl[i, c("Height", "Width")] <- dim(flsTblList[[i]])[1L:2L]
  }
  #View(flsTbl[, c("File", "Height", "Width")])
  # We want a 3000*3000 image to fit into a rough square of 20 rows and 6 columns
  # A row should be 20 pixels high
  # A column should be 64 pixels wide
  # 20*20
  # 64*6
  flsTbl$Max <- apply(flsTbl[, c("Width", "Height")], 1L, max)
  flsTbl$Width_Xl <- floor(4L*flsTbl$Width/flsTbl$Max)
  flsTbl$Height_Xl <- floor(4L*flsTbl$Height/flsTbl$Max)
  flsTbl$Col <- vapply(1L:nrow(flsTbl), \(i) {
    w <- which((flsTbl$x < flsTbl$x[i]) & (flsTbl$y == flsTbl$y[i]))
    if (length(w)) {
      res <- (sum(flsTbl$Width_Xl[w]+1L)*8.43+8.43)*1.5
      res <- which(cumsum(colWdths[2L:length(colWdths)]) > res)[1L]
    } else { res <- 2L }
    return(res)
  }, 1L)
  flsTbl$Row <- vapply(1L:nrow(flsTbl), \(i) {
    w <- which((flsTbl$y < flsTbl$y[i]) & (flsTbl$x == flsTbl$x[i]))
    if (length(w)) {
      res <- as.integer(ceiling(sum(flsTbl$Height_Xl[w]*6L+1L)))
    } else { res <- 0L }
    return(res+3L)
  }, 1L)
  #View(flsTbl[, c("File", "Height", "Width", "Row", "Col")])
  for (i in 1L:nImgs) {
    WorkBook <- wb_add_image(WorkBook, sheetnm, wb_dims(flsTbl$Row[i],
                                                        flsTbl$Col[i]),
                             flsTbl$File[i],
                             flsTbl$Width_Xl[i],
                             flsTbl$Height_Xl[i])
  }
  o <- o + max(flsTbl$Row)
}
# QC tables
dms <- wb_dims(2L + o, 2L)
if (exists("Exp_summary")) {
  WorkBook <- wb_add_data_table(WorkBook, sheetnm, Exp_summary,
                                dms, col_names = TRUE, row_names = FALSE,
                                table_name = "Experiment_overview",
                                first_column = TRUE,
                                banded_rows = TRUE)
  XpSum_OS <- nrow(Exp_summary) + 2L 
} else { XpSum_OS <- 0L }
if (exists("Modifs")) {
  temp <- Modifs[, c("Full name", "Mark", "Type", "AA")]
  w <- which(vapply(colnames(temp), \(x) { inherits(temp[[x]], "list") }, TRUE))
  if (length(w)) { for (k in colnames(temp)[w]) { temp[[k]] <- vapply(temp[[k]], paste, "", collapse = ", ") } }
  dms <- wb_dims(2L + XpSum_OS + o, 2L)
  WorkBook <- wb_add_data_table(WorkBook, sheetnm, temp,
                                dms, col_names = TRUE, row_names = FALSE,
                                table_name ="Amino_acid_compositional_biases",
                                first_column = TRUE, banded_rows = TRUE)
  Mods_OS <- nrow(Modifs) + 2L
} else { Mods_OS <- 0L }
dms <- wb_dims(2L + XpSum_OS + Mods_OS + o, 2L)
WorkBook <- wb_add_data_table(WorkBook, sheetnm, AA_biases,
                              dms, col_names = TRUE, row_names = FALSE,
                              table_name ="Modifications", first_column = TRUE, banded_rows = TRUE)
#wb_save(WorkBook, repFl);xl_open(repFl)
#
if ("tmp" %in% wb_get_sheet_names(WorkBook)) {
  WorkBook <- wb_remove_worksheet(WorkBook, "tmp")
}
tmp <- setNames(wb_get_order(WorkBook), wb_get_sheet_names(WorkBook))
mdpptbs <- grep("-mod. pept.$", names(tmp), value = TRUE)
nms <- c("Protein groups",
         "SAINTexpress",
         "Coverage",
         "All peptidoforms",
         mdpptbs,
         "Description",
         "Quality control")
nms <- intersect(nms, names(tmp))
tmp <- tmp[nms]
tmp <- tmp[which(!is.na(tmp))]
names(tmp) <- NULL
WorkBook <- wb_set_order(WorkBook, tmp)
dflt <- c("Protein groups", mdpptbs, "All peptidoforms")
dflt <- intersect(dflt, nms)[1L]
# This bit isn't currently working (report asap to Jan!)
m <- match(dflt, nms)
WorkBook <- wb_set_selected(WorkBook, m)
WorkBook <- wb_set_active_sheet(WorkBook, m)
WorkBook <- wb_set_base_font(WorkBook, 11L, font_name = "Calibri")
#
#
cat("    ---> writing table...\n")
wb_save(WorkBook, repFl)
#xl_open(repFl)
if (nImgs) { unlink(fls2) }
#
# Edit .xlsx (for the bits which openxlsx2 cannot handle well at the moment - or rather which I haven't yet figured out how to make it do!)
# - Unzip
cat("         final edits...\n")
dr <- paste0(wd, "/_unzipped")
unzip(repFl, exdir = dr)
# - Introduce new lines in column headers + fix selected and default tabs
library(xml2)
main <- paste0(dr, "/xl/workbook.xml")
doc <- xml2::read_xml(main)
ns <- xml2::xml_ns(doc)
sheets <- xml2::xml_find_all(doc, ".//d1:sheets/d1:sheet", ns)
nms <- xml2::xml_attr(sheets, "name")
sheetVis <- as.character(nms) == dflt
xml2::xml_attr(sheets, "tabSelected") <- NULL # Remove tabSelected from every sheet
idx <- match(dflt, nms) # Select desired sheet
xml2::xml_attr(sheets[[idx]], "tabSelected") <- "1"
# Set activeTab
view <- xml2::xml_find_first(doc, ".//d1:workbookView", ns)
xml2::xml_attr(view, "activeTab") <- as.character(idx - 1L)
# Write back
xml2::write_xml(doc, main)
#cat(tmp)
#
xmlFls <- list.files(dr, "\\.xml$", recursive = TRUE, full.names = TRUE)
xmlDat <- setNames(lapply(xmlFls, readLines, encoding = "UTF-8"),
                   xmlFls)
w1 <- which(vapply(xmlFls, \(fl) {
  any(grepl("///((VS)|(NL))///", xmlDat[[fl]]))
}, TRUE))
w2 <- which(xmlFls %in% paste0(dr, "/xl/worksheets/sheet", as.character(which(!sheetVis)), ".xml"))
w <- union(w1, w2)
nChars <- setNames(vapply(xmlFls[w], \(fl) { nchar(xmlDat[[fl]]) }, 1L), xmlFls[w])
chunk_size <- 5e7
# Some of these files are very, very... VERY large, so we want to process by chunks... but we also do not want to miss anything!
# So we will process by chunks:
pats1 <- c("///NL///", "///VS/// ")
rpls1 <- c("&#10;", "&#10;")
rplFun1 <- \(x) { gsub(pats1[1L], rpls1[1L], gsub(pats1[2L], rpls1[2L], x)) }
pats2 <- c("tabSelected=\"1\"")
rpls2 <- c("tabSelected=\"0\"")
rplFun2 <- \(x) { gsub(pats2[1L], rpls2[1L], x) }
pat <- paste(union(pats1, pats2), collapse = "|")
maxL <- max(nchar(union(pats1, pats2)))
for (fl in xmlFls[w]) { #fl <- xmlFls[w][1L] #fl <- names(nChars)[which.max(nChars)]
  tmp <- xmlDat[[fl]]
  nc <- nChars[fl]
  n <- ceiling(nc/chunk_size)
  rg <- round(seq_len(n)/n*nc) # NB: here the order of division/multiplication is important here to prevent integer overflow!!!
  if (nc <= chunk_size) {
    tst <- data.frame(start = 1L,
                      end = nc)
  } else {
    # If more than one chunk, let's verify that there is no overlaps between the area surrounding breaks and pattern matches!
    brks <- rg[seq_len(n-1L)]
    brksTst <- vapply(brks, \(x) {
      grepl(pat, substr(tmp, x-maxL, x+maxL+1L))
    }, TRUE)
    rg <- rg[setdiff(seq_len(n), which(brksTst))]
    n <- length(rg)
    tst <- data.frame(end = rg)
    tst$start <- c(1L, tst$end[seq_len(n-1L)]+1L)
  }
  tst$dat <- vapply(seq_len(nrow(tst)), \(i) {
    substr(tmp, tst$start[i], tst$end[i])
  }, "")
  if (fl %in% xmlFls[w1]) {
    tst$dat <- rplFun1(tst$dat)
  }
  if (fl %in% xmlFls[w2]) {
    tst$dat <- rplFun2(tst$dat)
  }
  con <- base::file(fl, "wb")
  for (x in tst$dat) {
    writeChar(enc2utf8(x), con, eos = NULL, useBytes = TRUE)
  }
  close(con)
}
# Note: it would be faster to concatenate the chunks from all files and process the whole thing together once.
# Parallelization could also be considered carefully here (e.g. serialize chunks to disk then parLapply over indices using readr::read_rds() to write each the current node)
#
gc()
# - Save final report
setwd(dr)
fls <- list.files(".", recursive = TRUE, all.files = TRUE)
zip(zipfile = repFl,
    files = fls)
setwd(wd)
xl_open(repFl)
cat("        Done!\n")
# shell(paste0("RMDIR /S /Q \"", dr, "\""), mustWork = FALSE)
