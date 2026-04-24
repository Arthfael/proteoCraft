dbOrd %<o% 1:nrow(db)
protDeflt %<o% NULL
if (!exists("prot.list")) {
  prot.list %<o% c()
}
wM <- c()
if (length(prot.list)) {
  wM <- which(db$`Protein ID` %in% prot.list)
  if (length(wM)) {
    dbOrd <- c(wM,
               which(!db$`Protein ID` %in% prot.list))
  }
}
protHeads %<o% setNames(gsub("^>", "", db$Header[dbOrd]), db$`Protein ID`[dbOrd])
nc <- nchar(protHeads)
w <- which(nc > 70)
if (length(w)) { protHeads[w] <- paste0(substr(protHeads[w], 1, 70), "...") }
protHeads2 %<o% protHeads
names(protHeads) <- NULL # Otherwise names can cause shenanigans (be used in shiny instead of the value)
lM <- length(wM)
if (lM) {
  protDeflt <- protHeads[1:lM]
}
nrm2PrtSlctDflt <- c()
if ((exists("Param"))&&
    ("Norma.Prot.Ratio.to.proteins" %in% colnames(Param))&&
    ("character" %in% class(Param$Norma.Prot.Ratio.to.proteins))&&
    (nchar(Param$Norma.Prot.Ratio.to.proteins))) {
  nrm2PrtSlctDflt <- unlist(strsplit(Param$Norma.Prot.Ratio.to.proteins, ";"))
  wM2 <- which(db$`Protein ID`[dbOrd] %in% nrm2PrtSlctDflt)
  if (!length(wM2)) { nrm2PrtSlctDflt <- NULL } else { nrm2PrtSlctDflt <- protHeads[wM2] }
} else { nrm2PrtSlctDflt <- NULL }
names(nrm2PrtSlctDflt) <- NULL
