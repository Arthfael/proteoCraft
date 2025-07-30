# Update peptide-to-protein mappings
# FragPipe only reports one protein per PSM, and among other search software, MaxQuant at least is also not fully exhaustive.
# While this may be correct from their point of view, I think I should report all matches.
setwd(wd)
#ev$Proteins <- gsub(";CON_", ";", gsub("^CON_", "", gsub(";CON__", ";", gsub("^CON__", "", ev$Proteins))))
if (Update_Prot_matches) {
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              paste0("body_add_fpar(Report, fpar(ftext(\" - Checking ", names(SearchSoft), "'s peptide sequence to protein assignments:\", prop = WrdFrmt$Body_text), fp_p = WrdFrmt$just))"))
  msg <- "ing peptide-to-protein matches...\n"
  if (exists("Reuse_Prot_matches")) {
    if ((!is.logical(Reuse_Prot_matches))||(is.na(Reuse_Prot_matches))) {
      Reuse_Prot_matches <- FALSE
    }
  }
  if (Reuse_Prot_matches) {
    msg <- paste0("Re-load", msg)
    cat(msg)
    if (!exists("evmatch")) { loadFun(paste0(wd, "/evmatch.RData")) }
    Pep2Prot <- evmatch
  } else {
    msg <- paste0("Check", msg)
    cat(msg)
    #
    source(parSrc, local = FALSE) # Always a good idea to check your cluster before such a big function...
    #
    Seq <- unique(ev$Sequence)
    DB <- db
    mtchSrc <- paste0(libPath, "/extdata/R scripts/Sources/ProtMatch2.R")
    source(mtchSrc)
    Pep2Prot <- evmatch
  }
  wh1 <- which(ev$Sequence %in% Pep2Prot$Sequence)
  wh2 <- which(Pep2Prot$Sequence %in% ev$Sequence)
  mtch1 <- match(ev$Sequence[wh1], Pep2Prot$Sequence)
  if (!"Proteins" %in% colnames(ev)) { ev$Proteins <- "" } else {
    source(parSrc, local = FALSE)
    tmpPs <- unique(Pep2Prot$Sequence[wh2])
    tmpE <- ev$Proteins[match(tmpPs, ev$Sequence)]
    tmpP <- Pep2Prot$Proteins[match(tmpPs, Pep2Prot$Sequence)]
    tmpE <- strsplit(tmpE, ";")
    tmpP <- strsplit(tmpP, ";")
    f0 <- function(x) { paste(sort(x), collapse = ";") }
    tst1 <- parSapply(parClust, tmpE, f0)
    tst2 <- parSapply(parClust, tmpP, f0)
    wN <- which(tst1 != tst2)
    ViewTst <- FALSE
    if ((length(wN))&&(ViewTst)) {
      # Below some code to help with investigations...
      tst <- data.frame(Seq = tmpPs[wN],
                        Orig = tst1[wN],
                        Corr = tst2[wN])
      tst$Orig <- strsplit(tst$Orig, ";")
      #View(tst)
      tst$Corr <- strsplit(tst$Corr, ";")
      f0 <- function(x, y) { x[which(!x %in% y)] }
      tst$In_Orig_only <- lapply(1:nrow(tst), function(x) { f0(unlist(tst$Orig[[x]]),
                                                               unlist(tst$Corr[[x]])) })
      tst$In_Corr_only <- lapply(1:nrow(tst), function(x) { f0(unlist(tst$Corr[[x]]),
                                                               unlist(tst$Orig[[x]])) })
      tst$In_Orig_only_in_db <- vapply(tst$In_Orig_only, function(x) { sum(x %in% db$`Protein ID`) }, 1)
      sum(unlist(tst$In_Orig_only_in_db))
      sum(!unlist(tst$In_Orig_only_in_db))
      w <- which(tst$In_Orig_only_in_db > 0)
      View(tst[w,])
      # lapply(1:nrow(tst), function(x) { f0(unlist(tst$Orig[[x]]),
      #                                      unlist(tst$Corr[[x]])) })
    }
    tst <- sum(tst1 != tst2)
    if (tst) {
      msg <- paste0("Corrected ", tst, " out of ", nrow(Pep2Prot), " assignments (~", round(tst/nrow(Pep2Prot), 2), "%)!
Note that we did not take into account retention time or ion mobility!
Discrepancies with the original search engine matches can have several causes:
 1) Protein present in original column but not corrected results:
   a) If redundant entries were present in the search fasta, then they would have usually been filtered by this workflow when it loads and parses the fasta database.
   b) More rarely, the accession is present in the db... but for every peptide we have checked so far the protein sequence was not compatible with it (even allowing for I/L ambiguity)!!!
 2) Present in corrected but not original:
   a) Some search engines seem to miss some protein matches (at least old versions of MaxQuant) or to report only one out of several possible protein matches (some versions of FragPipe).
   b) In the case of I/L ambiguity, for newer search engines using modern spectrum prediction, retention time and/or ion mobility prediction models, this may be actually correct and reflect incompatibility of the latter characteristics.
")
    } else { msg <- "All assignments validated.\n" }
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
  }
  ev$Proteins[wh1] <- Pep2Prot$Proteins[mtch1]
  if (!Reuse_Prot_matches) {
    # Save in case you need to reload:
    evmatch <- ev[, c("Sequence", "Proteins")]
    saveFun(evmatch, file = paste0(wd, "/evmatch.RData"))
  }
} else {
  kol <- c("Leading proteins", "Proteins")
  kol <- kol[which(kol %in% colnames(ev))]
  tmp <- ev[, kol, drop = FALSE]
  for (k in kol) { tmp[[k]] <- strsplit(tmp[[k]], ";") }
  ev$Proteins <- parApply(parClust, tmp, 1, function(x) {
    paste(unique(unlist(x)), collapse = ";")
  })
}
tst <- unique(unlist(strsplit(ev$Proteins, ";")))
if ("NA" %in% tst) { stop("\"NA\" is not an accepted protein accession!") }
#View(ev[, c("Sequence", "Proteins")])
