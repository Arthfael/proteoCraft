# Update peptide-to-protein mappings
# FragPipe only reports one protein per PSM, and among other search software, MaxQuant at least is also not fully exhaustive.
# While this may be correct from their point of view, I think I should report all matches.
I_eq_L %<o% TRUE # Could be made a global parameter

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
    #rstudioapi::documentOpen(mtchSrc)
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
      msg <- paste0("Corrected ", tst, " out of ", length(tst1), " assignments (~", round(100*tst/length(tst1), 2), "%)!
Note that we do not take into account retention time or ion mobility!
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
#test <- c()
#if (!is.null(prot.list)) {
#  a <- paste0(";", ev$Proteins, ";")
#  b <- prot.list
#  a1 <- paste0(";", b, ";")
#  test <- unique(unlist(lapply(a1, function(x) { ev$id[grep(x, a)] })))
#}
#kol <- which(toupper(colnames(ev)) %in% c("CONTAMINANT", "POTENTIAL CONTAMINANT"))
#ev <- ev[which((is.na(ev[[kol]]))|(ev[[kol]] == "")|(ev$id %in% test)),]
# Test if there are still any evidences without any matching proteins from the database
ev$"Tryptic peptide?" <- TRUE
w <- which(ev$Proteins == "")
l <- length(w)
if (l) {
  tst <- (l>1)+1
  msg <- paste0("There ", c("is", "are")[tst], " ", length(w), " PSM", c("", "s")[tst],
                " (", round(100*l/nrow(ev)), "%) without any matching protein from the database, probably from ",
                c("a ", "")[tst], "non-tryptic peptide", c("", "s")[tst], ".")
  ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
  ev$"Tryptic peptide?"[w] <- FALSE
  tmpEV <- data.frame(Seq = unique(ev$Sequence[w]))
  tmpEV$Seq2 <- tmpEV$Seq
  tmpDB <- db[, c("Protein ID", "Sequence")]
  if (I_eq_L) {
    tmpEV$Seq2 <- gsub("I", "L", tmpEV$Seq2)
    tmpDB$Sequence <- gsub("I", "L", tmpDB$Sequence)
  }
  source(parSrc) # Check cluster for corruption!
  saveRDS(tmpEV, "tmpEV.RDS")
  saveRDS(tmpDB, "tmpDB.RDS")
  invisible(clusterCall(parClust, function() {
    tmpEV <<- readRDS("tmpEV.RDS")
    tmpDB <<- readRDS("tmpDB.RDS")
    return()
  }))
  tmpEV$Proteins <- parSapply(parClust, 1:nrow(tmpEV), function(x) {
    paste(tmpDB$"Protein ID"[grep(tmpEV$Seq2[x], tmpDB$Sequence)], collapse = ";")
  })
  invisible(clusterCall(parClust, function() {
    rm(tmpDB)
    return()
  }))
  unlink("tmpEV.RDS")
  unlink("tmpDB.RDS")
  ev$Proteins[w] <- tmpEV$Proteins[match(ev$Sequence[w], tmpEV$Seq)]
  ev <- ev[which(ev$Proteins != ""),]
}
# Also remove those protein columns we will re-create later
ev$"Leading proteins" <- NULL
ev$"Leading razor protein" <- NULL
