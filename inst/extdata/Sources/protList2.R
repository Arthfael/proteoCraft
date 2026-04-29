w <- which(!prot.list %in% db$`Protein ID`)
if (length(w)) {
  urls <- paste0("https://rest.uniprot.org/uniprotkb/", prot.list[w], ".fasta")
  tst <- setNames(lapply(urls, \(url) {
    rs <- try(readLines(url), silent = TRUE)
    rs <- if (inherits(rs, "try-error")) { list(Outcome = FALSE) } else {
      list(Outcome = TRUE, Res = rs)
    }
    return(rs)
  }), prot.list[w])
  pY <- names(tst)[which(vapply(tst, \(x) { x$Outcome }, TRUE))]
  pN <- names(tst)[which(vapply(tst, \(x) { !x$Outcome }, TRUE))]
  if (length(pN)) { prot.list <- prot.list[which(!prot.list %in% pN)] }
  if (length(pY)) {
    tst <- lapply(tst[pY], \(x) { Format.DB(x$Res, in.env = TRUE) })
    tst <- plyr::rbind.fill(tst)
    tst$Source <- "UniprotKB"
    tst$"Potential contaminant" <- ""
    tst$"Protein of interest" <- TRUE
    db <- plyr::rbind.fill(tst, db)
  }
}
