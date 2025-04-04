# Special quality control tab
sheetnm <- "Quality control"
sheetnmsA <- unique(c(sheetnmsA, sheetnm))
if (sheetnm %in% wb_get_sheet_names(WorkBook)) { WorkBook <- wb_remove_worksheet(WorkBook, sheetnm) }
WorkBook <- wb_add_worksheet(WorkBook, sheetnm)
WorkBook <- wb_add_data(WorkBook, sheetnm, "Summary plots and tables", wb_dims(2, 2))
stl <- wb_get_cell_style(WorkBook, "tmp", wb_dims(match("Header 1", styleNms), 1))
if (stl != "") { WorkBook <- wb_set_cell_style(WorkBook, sheetnm, wb_dims(2, 2), stl) }
WorkBook <- wb_add_cell_style(WorkBook, sheetnm, wb_dims(2, 2), wrap_text = FALSE)
WorkBook <- wb_set_col_widths(WorkBook, sheetnm, 1, 5)
tmp <- ceiling(nchar("Summary plots and tables")/8.43)*8.43+5
WorkBook <- wb_set_col_widths(WorkBook, sheetnm, 2, tmp)
WorkBook <- wb_set_row_heights(WorkBook, sheetnm, 1, 15)
colWdths <- rep(8.43, 1000) # Should be enough, what say you?
colWdths[1] <- 5
colWdths[2] <- tmp
fls <- c(list.files(paste0(wd, "/Summary plots"), ".jpeg", full.names = TRUE, recursive = TRUE),
         list.files(paste0(wd, "/Workflow control"), ".jpeg", full.names = TRUE, recursive = TRUE))
fls <- fls[which(file.exists(fls))]
nImgs <- length(fls)
o <- 24
if (nImgs) {
  flsTbl <- data.frame(File = fls,
                       x = ceiling(1:nImgs/2),
                       y = (((1:nImgs)+1) %% 2) + 1,
                       Width = 0,
                       Height = 0)
  flsTblList <- sapply(1:nImgs, list)
  for (i in 1:nImgs) {
    flsTblList[[i]] <- readJPEG(flsTbl$File[i])
    flsTbl[i, c("Height", "Width")] <- dim(flsTblList[[i]])[1:2]
  }
  # We want a 3000*3000 image to fit into a rough square of 20 rows and 6 columns
  # A row should be 20 pixels high
  # A column should be 64 pixels wide
  # 20*20
  # 64*6
  flsTbl$Max <- apply(flsTbl[, c("Width", "Height")], 1, max)
  flsTbl$Width_Xl <- 4*flsTbl$Width/flsTbl$Max
  flsTbl$Height_Xl <- 4*flsTbl$Height/flsTbl$Max
  flsTbl$Col <- sapply(1:nrow(flsTbl), function(i) {
    w <- which((flsTbl$x < flsTbl$x[i])&(flsTbl$y == flsTbl$y[i]))
    if (length(w)) {
      res <- (sum(flsTbl$Width_Xl[w]+1)*8.43+8.43)*1.5
      res <- which(cumsum(colWdths[2:length(colWdths)]) > res)[1]
    } else { res <- 2 }
    return(res)
  })
  flsTbl$Row <- sapply(1:nrow(flsTbl), function(i) {
    w <- which((flsTbl$y < flsTbl$y[i])&(flsTbl$x == flsTbl$x[i]))
    if (length(w)) {
      res <- ceiling(sum(flsTbl$Height_Xl[w]*24/4+1))
    } else { res <- 0 }
    return(res+3)
  })
  for (i in 1:nImgs) {
    WorkBook <- wb_add_image(WorkBook, sheetnm, wb_dims(flsTbl$Row[i], flsTbl$Col[i]), flsTbl$File[i],
                             flsTbl$Width_Xl[i]*1.2, flsTbl$Height_Xl[i]*1.2)
  }
  o <- o+max(flsTbl$Row)
}
# QC tables
dms <- wb_dims(2+o, 2)
if (exists("Exp_summary")) {
  WorkBook <- wb_add_data_table(WorkBook, sheetnm, Exp_summary,
                                dms, col_names = TRUE, row_names = FALSE,
                                table_name = "Experiment_overview",
                                first_column = TRUE,
                                banded_rows = TRUE)
  XpSum_OS <- nrow(Exp_summary) + 2 
} else { XpSum_OS <- 0 }
if (exists("Modifs")) {
  temp <- Modifs[, c("Full name", "Mark", "Type", "AA")]
  w <- which(sapply(colnames(temp), function(x) { class(temp[[x]]) }) == "list")
  if (length(w)) { for (k in colnames(temp)[w]) { temp[[k]] <- sapply(temp[[k]], paste, collapse = ", ") } }
  dms <- wb_dims(2 + XpSum_OS + o, 2)
  WorkBook <- wb_add_data_table(WorkBook, sheetnm, temp,
                                dms, col_names = TRUE, row_names = FALSE,
                                table_name ="Amino_acid_compositional_biases",
                                first_column = TRUE, banded_rows = TRUE)
  Mods_OS <- nrow(Modifs) + 2
} else { Mods_OS <- 0 }
dms <- wb_dims(2 + XpSum_OS + Mods_OS + o, 2)
WorkBook <- wb_add_data_table(WorkBook, sheetnm, AA_biases,
                              dms, col_names = TRUE, row_names = FALSE,
                              table_name ="Modifications", first_column = TRUE, banded_rows = TRUE)
#wb_save(WorkBook, repFl);xl_open(repFl)
#
tmp <- setNames(wb_get_order(WorkBook), wb_get_sheet_names(WorkBook))
mdpptbs <- grep("-mod. pept.$", names(tmp), value = TRUE)
nms <- c("Protein groups",
         "SAINTexpress",
         "Coverage",
         "All peptidoforms",
         mdpptbs,
         "Description",
         "Quality control")
if ("tmp" %in% names(tmp)) { nms <- c(nms, "tmp") } # For when re-running/debugging
nms <- nms[which(nms %in% names(tmp))]
tmp <- tmp[nms]
tmp <- tmp[which(!is.na(tmp))]
names(tmp) <- NULL
WorkBook <- wb_set_order(WorkBook, tmp)
dflt <- c("Protein groups", mdpptbs, "All peptidoforms")
dflt <- dflt[which(dflt %in% nms)][1]
nms <- wb_get_sheet_names(WorkBook)
WorkBook <- wb_set_selected(WorkBook, match(dflt, nms))
WorkBook <- wb_set_active_sheet(WorkBook, match(dflt, nms))
WorkBook <- wb_set_base_font(WorkBook, 11, font_name = "Calibri")
#
if ("tmp" %in% wb_get_sheet_names(WorkBook)) { WorkBook <- wb_remove_worksheet(WorkBook, "tmp") }
#
cat("Writing table...\n")
wb_save(WorkBook, repFl)
xl_open(repFl)
