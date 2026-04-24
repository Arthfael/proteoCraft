#VP_list
#subDr
#wd
#insrt
#
subDr2 <- subDr
if (!grepl(topattern(paste0(wd, "/")), subDr2)) {
  subDr2 <- paste0(wd, "/", subDr2)
}
stopifnot(dir.exists(subDr2))
fl <- paste0(subDr2, "/Thresholds", insrt,".xlsx")
thresh <- lapply(names(VP_list$Thresholds$Absolute), function(x) {
  y <- VP_list$Thresholds$Absolute[[x]]
  x <- data.frame(Group = rep(cleanNms(x), nrow(y)))
  return(cbind(x, y))
})
thresh <- plyr::rbind.fill(thresh)
thresh$Name <- NULL
if ("Root" %in% colnames(thresh)) { thresh$Root <- gsub(" - $", "", thresh$Root) }
thresh$Value <- thresh$Text.value
thresh$Text.value <- NULL
wb <- wb_workbook()
wb <- wb_set_creators(wb, "Me")
wb <- wb_add_worksheet(wb, "Thresholds")
dms <- wb_dims(2, 1)
wb <- wb_add_data(wb, "Thresholds", "Absolute thresholds", dms)
wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                  italic = "true", underline = "single")
dms <- wb_dims(3, 2)
wb <- wb_add_data_table(wb, "Thresholds", thresh, dms,
                        col_names = TRUE, table_style = "TableStyleMedium2",
                        banded_rows = TRUE, banded_cols = FALSE)
tmp1 <- rbind(colnames(thresh), thresh)
colnames(tmp1) <- paste0("V", 1:ncol(tmp1))
if ("FDR" %in% names(VP_list$Thresholds)) {
  dms <- wb_dims(3+nrow(thresh)+3, 1)
  wb <- wb_add_data(wb, "Thresholds", "FDR thresholds", dms)
  wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                    italic = "true", underline = "single")
  fdrThresh <- VP_list$Thresholds$FDR
  fdrThresh$Group <- cleanNms(fdrThresh$Sample)
  fdrThresh$Sample <- NULL
  colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.up")] <- "Colour (up)"
  colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.down")] <- "Colour (down)"
  colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.line")] <- "Colour (line)"
  fdrThresh <- fdrThresh[, c("Group", colnames(fdrThresh)[which(colnames(fdrThresh) != "Group")])]
  dms <- wb_dims(3+nrow(thresh)+4, 2)
  wb <- wb_add_data_table(wb, "Thresholds", fdrThresh, dms,
                          col_names = TRUE, table_style = "TableStyleMedium2",
                          banded_rows = TRUE, banded_cols = FALSE)
  wb <- wb_set_col_widths(wb, "Thresholds", 1, 3)
  tmp2 <- rbind(colnames(fdrThresh), fdrThresh)
  colnames(tmp2) <- paste0("V", 1:ncol(tmp2))
  tst <- plyr::rbind.fill(tmp1, tmp2)
} else {
  tst <- tmp1
}
tst <- setNames(apply(tst, 2, function(x) { max(nchar(x), na.rm = TRUE)*1.5 }), NULL)
wb <- wb_set_col_widths(wb, "Thresholds", 1:(length(tst)+1), c(6, tst))
wb_save(wb, fl)
#xl_open(fl)
#
