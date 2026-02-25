# XML coverage columns
require(openxlsx2)
source(parSrc)
xmlCovCol %<o% c()
frstProt <- gsub(";.*", "", PG$"Leading protein IDs")
tmpLs <- list()
kolRt <- "Peptide IDs"
if (scrptType == "withReps") { Grps <- VPAL$values }
if (scrptType == "noReps") { Grps <- c("", Exp) }
for (grp in Grps) { #grp <- VPAL$values[1]
  if (scrptType == "withReps") { pepidkol <- paste0(kolRt, " - ", grp) }
  if (scrptType == "noReps") {
    if (grp == "") { pepidkol <- kolRt } else { pepidkol <- paste0(kolRt, " - ", grp) }
  }
  temp1 <- data.frame(Accession = frstProt,
                      `Peptide IDs` = PG[[pepidkol]],
                      check.names = FALSE)
  temp1$Sequence <- db$Sequence[match(temp1$Accession, db$"Protein ID")]
  temp1$"Peptide IDs" <- strsplit(temp1$"Peptide IDs", ";")
  tmp <- listMelt(temp1$`Peptide IDs`, 1:nrow(temp1))
  tmp <- data.table(Row = as.integer(tmp$L1), pepID = tmp$value)
  tmp$Seq <- pep$Sequence[match(tmp$pepID, pep$id)]
  tmp <- tmp[, list(Seq = list(Seq)), keyby = list(Row)]
  tmp <- as.data.frame(tmp)
  temp1$PepSeq <- tmp$Seq[match(1:nrow(temp1), tmp$Row)]
  tmpLs[[grp]] <- temp1
}
tmpLs <- plyr::rbind.fill(tmpLs)
tmpLs <- tmpLs[, c("Sequence", "PepSeq")]
require(openxlsx2)
#clusterExport(parClust, "Coverage", envir = environment())
clusterCall(parClust, function() library(openxlsx2))
res <- parApply(parClust, tmpLs[, c("Sequence", "PepSeq")], 1, function(x) { #x <- tmpLs[1, c("Sequence", "PepSeq")]
  prsq <- x[[1]]
  ppsq <- x[[2]]
  if (!is.na(prsq)) {
    rs <- proteoCraft::Coverage(prsq, ppsq, "XML", colour = "#006666")
  } else {
    rs <- list(openxlsx2::fmt_txt("", color = openxlsx2::wb_color(hex = "grey")))
  }
  return(rs)
})
res <- sapply(res, function(x) { x[[1]] })
n <- nrow(PG)
covRt <- "1st ID cov."
for (grp in Grps) { #grp <- VPAL$values[1]
  if (scrptType == "withReps") { kol <- paste0(covRt, " - ", grp) }
  if (scrptType == "noReps") {
    if (grp == "") { kol <- covRt } else { kol <- paste0(covRt, " - ", grp) }
  }
  xmlCovCol <- unique(c(xmlCovCol, kol))
  PG[[kol]] <- res[(1:n)+(n*(match(grp, Grps)-1))]
}
# This chunk is because of Titin. You know who you are, you bad Excel-breaking protein! 
tst <- setNames(lapply(xmlCovCol, function(x) {
  w <- setNames(which(parSapply(parClust, substr(PG[[x]], 1, ExcelMax), function(y) {
    "try-error" %in% class(try(openxlsx2::as_xml(y), silent = TRUE))
  })), NULL)
}), cleanNms(xmlCovCol, start = FALSE))
sapply(tst, length)
tmp <- fmt_txt("This sequence is too long to display in an Excel cell!",
               bold = TRUE, italic = TRUE, color = wb_colour("red"))
for (kol in xmlCovCol) { #kol <- xmlCovCol[1]
  w <- tst[[cleanNms(kol, start = FALSE)]]
  if (length(w)) { PG[w, kol] <- tmp }
}
