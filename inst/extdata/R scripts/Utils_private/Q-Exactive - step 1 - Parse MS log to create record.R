logsDir <- "C:/Xcalibur/system/Exactive/log" # Should contain all of the instruments's logs... for ever!
stopifnot(dir.exists(logsDir))
#logsDir <- "D:/logs" 
logFls <- list.files(logsDir, pattern = "Thermo Exactive--")
PrcsDir <- "B:/group/mspecgrp/Prices"
stopifnot(dir.exists(PrcsDir))

Res <- lapply(logFls, function(logFl) { #logFl <- logFls[1]
  txt <- readLines(paste0(logsDir, "/", logFl))
  g1 <- grep("Preparing next acquisition with\\: ", txt)
  if (!length(g1)) { return() }
  #g2 <- grep("MS acquisition end\\. Accumulated MS size\\: ", txt)
  g1pr <- length(txt)
  if (length(g1) >= 2) { g1pr <- c(g1[2:length(g1)]-1, g1pr) }
  g2 <- sapply(1:length(g1), function(x) {
    #x <- 1
    g1[x]+grep("MS acquisition end\\. Accumulated MS size\\: ", txt[(g1[x]+1):g1pr[x]])[1]
  })
  g1 <- g1[which(!is.na(g2))]
  g2 <- g2[which(!is.na(g2))]
  tst <- grepl("MS acquisition end\\. Accumulated MS size\\: ", txt[g2])
  stopifnot(sum(!tst) == 0)
  flnms <- gsub("[^>]+>.+", "", gsub("<[^>]+>.+", "", gsub(".+<InstrumentD ", "", txt[g1])))
  flnms[which(flnms == "")] <- "... "
  flnms <- gsub("[^ ]+\\.meth", "... ", flnms)
  sqnm <- sapply(strsplit(txt[g1], "</?SequenceFileName>"), function(x) { x[[2]] })
  init <- gsub("\\+[0-9]{2}:[0-9]{2}$", "", gsub("^\\[Time=|\\].+", "", txt[g1]))
  if (!length(init)) { return() }
  finit <- gsub("\\+[0-9]{2}:[0-9]{2}$", "", gsub("^\\[Time=|\\].+", "", txt[g2]))
  diff <- round(as.POSIXct(finit)-as.POSIXct(init), 2)
  if (!is.finite(min(diff))) { stop() }
  stopifnot(min(diff) > 0)
  date <- sapply(strsplit(gsub(" .+", "", init), "-"), function(x) { paste(rev(unlist(x)), collapse = "/") })
  sz <- as.numeric(gsub(".+MS acquisition end\\. Accumulated MS size\\: |KB$", "", txt[g2]))
  vial <- sapply(strsplit(txt[g1], "</?Vial>"), function(x) { x[[2]] })
  tmp <- data.frame("start" = as.POSIXct(init),
                    "end" = as.POSIXct(finit),
                    "filename" = flnms,
                    "date" = date,
                    "time min" = diff,
                    "size kb" = sz,
                    "sequence path" = sqnm,
                    "vial" = vial,
                    "blank" = c("", "+")[grepl("^[BGR][0-9]$", vial)+1],
                    check.names = FALSE)
  return(tmp)
})
Res <- plyr::rbind.fill(Res)

Res <- Res[order(Res$start, decreasing = FALSE),]
Yrs <- as.numeric(gsub("-.*", "", Res$start))
Mnths <- as.numeric(gsub("-.*", "", gsub("^[0-9]+-", "", Res$start)))
Mnths <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[Mnths]
require(openxlsx)
wb <- createWorkbook()
styleDflt <- createStyle(halign = "left")
styleTim <- createStyle(halign = "left", numFmt = "LONGDATE")
wTm <- which(colnames(Res) %in% c("start", "end"))
wNm <- which(colnames(Res) == "filename")
wSq <- which(colnames(Res) == "sequence path")
for (Yr in unique(Yrs)) {
  Yr <- as.character(Yr)
  print(Yr)
  w <- grep(paste0("^", Yr, "-"), Res$start)
  sheet  <- addWorksheet(wb, Yr)
  writeDataTable(wb, Yr, Res[w,], tableStyle = "TableStyleMedium2")
  addStyle(wb, Yr, styleDflt, (1:nrow(Res))+1, 1:ncol(Res), TRUE, TRUE)
  addStyle(wb, Yr, styleTim, (1:nrow(Res))+1, wTm, TRUE, TRUE)
  setColWidths(wb, Yr, 1:ncol(Res), 12)
  setColWidths(wb, Yr, c(wTm, wNm), 20)
  setColWidths(wb, Yr, wSq, 100)
}
pthRoot <- paste0(PrcsDir, "/runtimes-", Yrs[1], "-", Yrs[length(Yrs)], Mnths[length(Mnths)])
saveWorkbook(wb, paste0(pthRoot, ".xlsx"), overwrite = TRUE)
openxlsx::openXL(paste0(pthRoot, ".xlsx"))
save(Res, file = paste0(pthRoot, ".RData"))
#load(paste0(pthRoot, ".RData"))
