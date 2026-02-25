w <- which(!prot.list %in% db$`Protein ID`)
if (length(w)) {
  urls <- paste0("https://rest.uniprot.org/uniprotkb/", prot.list[w], ".fasta")
  tst <- setNames(lapply(urls, function(url) {
    rs <- try(readLines(url), silent = TRUE)
    if ("try-error" %in% class(rs)) { rs <- list(Outcome = FALSE) } else {
      rs <- list(Outcome = TRUE, Res = rs)
    }
    return(rs)
  }), prot.list[w])
  pY <- names(tst)[which(sapply(tst, function(x) { x$Outcome }))]
  pN <- names(tst)[which(sapply(tst, function(x) { !x$Outcome }))]
  if (length(pN)) { prot.list <- prot.list[which(!prot.list %in% pN)] }
  if (length(pY)) {
    tst <- lapply(tst[pY], function(x) { Format.DB(x$Res, in.env = TRUE) })
    tst <- plyr::rbind.fill(tst)
    tst$Source <- "UniprotKB"
    tst$"Potential contaminant" <- ""
    tst$"Protein of interest" <- TRUE
    db <- plyr::rbind.fill(tst, db)
  }
}
