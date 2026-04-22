# Peptides: create pep object
require(data.table)
source(parSrc, local = FALSE)
if (scrptType == "withReps") {
  tmp_EM <- Exp.map
  colnames(tmp_EM)[which(colnames(tmp_EM) == "Ref.Sample.Aggregate")] <- "Parent_sample"
}
if (scrptType == "noReps") {
  ev$MQ.Exp <- FracMap$MQ.Exp[match(ev$"Raw file path", FracMap$"Raw file")]
  tmp_EM <- aggregate(FracMap$MQ.Exp, list(FracMap$"Parent sample"), unique)
  colnames(tmp_EM) <- c("Parent_sample", "MQ.Exp")
  tmp_EM$Use <- FracMap$Use[match(tmp_EM$Parent_sample, FracMap$"Parent sample")]
}
tmp <- as.character(ev$id)
tmp2 <- data.table(id = tmp, mqxp = ev$MQ.Exp, mod = ev$"Modified sequence")
tmp2 <- tmp2[order(ev$id, decreasing = FALSE),]
tmpFl1 <- tempfile(fileext = ".rds")
tmpFl2 <- tempfile(fileext = ".rds")
clusterExport(parClust, list("tmpFl1", "tmpFl2", "scrptType"), envir = environment())
readr::write_rds(tmp_EM, tmpFl1)
readr::write_rds(tmp2, tmpFl2)
invisible(clusterCall(parClust, \() {
  library(data.table)
  assign("tmp_EM", readr::read_rds(tmpFl1), envir = .GlobalEnv)
  assign("tmp2", readr::read_rds(tmpFl2), envir = .GlobalEnv)
  return()
}))
unlink(tmpFl1)
unlink(tmpFl2)
myGrps <- "ALLMYSAMPLESTUDUDUDUMMMDADA"
if (scrptType == "withReps") {
  myGrps <- c(myGrps, RSA$values)
}
if (scrptType == "noReps") {
  myGrps <- c(myGrps, Exp)
}
tmp4 <- setNames(parLapply(parClust, myGrps, \(i) { #i <- myGrps[1L] #i <- myGrps[2L]
  tmp3 <- data.table::copy(tmp2) # Because of how data.tables work! No idea how that plays with clusters...
  if (i != "ALLMYSAMPLESTUDUDUDUMMMDADA") {
    mqe <- unique(unlist(tmp_EM$MQ.Exp[which(tmp_EM$Parent_sample == i)]))
    tmp3 <- tmp3[which(tmp3$mqxp %in% mqe), c("id", "mod")]
  }
  tmp3 <- as.data.frame(tmp3[, list(IDs = paste(id, collapse = ";")),
                             keyby = list(ModSeq = mod)])
  return(tmp3)
}), myGrps)
pep %<o% magrittr::set_colnames(tmp4[["ALLMYSAMPLESTUDUDUDUMMMDADA"]],
                                c("Modified sequence", "Evidence IDs"))
for (i in myGrps[2L:length(myGrps)]) {
  tmp <- tmp4[[i]]
  ki <- paste0("Evidence IDs - ", i)
  pep[[ki]] <- tmp$IDs[match(pep$"Modified sequence", tmp$ModSeq)]
  pep[which(is.na(pep[[ki]])), ki] <- ""
}
pep$id <- 1L:nrow(pep)
rvmtch2 <- match(pep$"Modified sequence", ev$"Modified sequence")
if ("Modified sequence_verbose" %in% colnames(ev)) {
  pep$"Modified sequence_verbose" <- ev$"Modified sequence_verbose"[rvmtch2]
}
pep$Sequence <- ev$Sequence[rvmtch2]
pep$Proteins <- ev$Proteins[rvmtch2]
tmp <- data.table(mod = ev$"Modified sequence", PEP = ev$PEP)
tmp <- as.data.frame(tmp[, list(PEP = min(PEP, na.rm = TRUE)), by = list(mod)])
pep$PEP <- tmp$PEP[match(pep$"Modified sequence", tmp$mod)]
mtch <- match(ev$Sequence, pep$Sequence)
mtch2 <- match(ev$"Modified sequence", pep$"Modified sequence")
# Amino Acid counts
sq <- pep$Sequence
clusterExport(parClust, "sq", envir = environment())
tmp <- parSapply(parClust, proteoCraft::AA, \(aa) {
  nchar(sq) - nchar(gsub(aa, "", sq))
})
colnames(tmp) <- paste0(colnames(tmp), " Count")
pep[, colnames(tmp)] <- tmp
ev[, paste0(AA, " Count")] <- pep[mtch, paste0(AA, " Count")]
for (aa in c("O", "U")) { # Only keep the selenocysteine and pyrrolysine amino acid columns if they are non-empty (NB: pyrrolysine should really only be in some bacteria)
  if (sum(ev[[paste0(aa, " Count")]]) == 0L) {
    ev[[paste0(aa, " Count")]] <- NULL
    pep[[paste0(aa, " Count")]] <- NULL
  }
}
#
pep$Length <- nchar(pep$Sequence)
ev$Length <- pep$Length[mtch]
tmp <- data.table(mod = ev$"Modified sequence", Intensity = ev$Intensity)
tmp$Intensity[which(!is.all.good(tmp$Intensity, 2L))] <- NA
w2 <- which(ev$MQ.Exp %in% unique(unlist(tmp_EM$MQ.Exp[which(tmp_EM$Use)])))
tmp2 <- copy(tmp)
tmp2 <- tmp2[w2, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(mod)]
pep$Intensity <- tmp2$Intensity[match(pep$"Modified sequence", tmp2$mod)]
if ((sum(!tmp_EM$Use))||(length(unique(ev$MQ.Exp)) > length(unique(unlist(tmp_EM$MQ.Exp))))) {
  tmp3 <- copy(tmp)
  tmp3 <- tmp3[, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(mod)]
  pep$"Intensity (all samples including unused)" <- tmp3$Intensity[match(pep$"Modified sequence", tmp3$mod)]
}
ev$"Peptide ID" <- pep$id[mtch2]
#
# Peptide name
tmp <- pep[, c("Proteins", "Modified sequence")]
w <- which(nchar(tmp$Proteins) > 15L)
tmp$Proteins[w] <- paste0(substr(tmp$Proteins[w], 1L, 12L), "...")
tmp$"Modified sequence" <- gsub("_", "", tmp$"Modified sequence")
pep$"Peptide name" <- do.call(paste, c(tmp, sep = "\n"))
ev$"Peptide name" <- pep$"Peptide name"[mtch2]

# Modification site probabilities, score differences, site IDs, etc...
if (exists("Modifs")) {
  wPr <- which(paste0(Modifs$"Full name", " Probabilities") %in% colnames(ev))
  if (length(wPr)) {
    for (i in wPr) { #i <- wPr[1L]
      for (j in c(" Probabilities", " Score Diffs")) { #j <- " Probabilities"
        j1 <- paste0(Modifs$"Full name"[i], j)
        temp <- ev[, c("Modified sequence", j1, "PEP")]
        temp <- temp[which(temp[[j1]] != ""),]
        a0 <- unique(temp$"Modified sequence")
        clusterExport(parClust, list("temp", "j1"), envir = environment())
        b1 <- parLapply(parClust, a0, \(x) { temp[which(temp$"Modified sequence" == x), j1] })
        b2 <- parLapply(parClust, a0, \(x) { temp$PEP[which(temp$"Modified sequence" == x)] })
        l <- vapply(b1, \(x) { length(unique(x)) }, 1L)
        wb <- which(l == 1L)
        if (length(wb)) { b1[wb] <- sapply(b1[wb], \(x) { unique(x) }) } 
        wb <- which(l > 1L)
        if (length(wb)) {
          temp1 <- cbind(b1[wb], b2[wb])
          clusterExport(parClust, list("temp1", "annot_to_tabl", "Isapply", "is.all.good", "insertElems"), envir = environment())
          b1[wb] <- parApply(parClust, temp1, 1L, \(x) {
            p <- 1L-x[[2L]]
            p <- p/sum(p)
            x <- x[[1L]]
            s <- unlist(strsplit(unique(gsub("\\([^\\)]+\\)", "", x)), ""))
            an <- Isapply(x, \(y) {
              as.numeric(annot_to_tabl(y, Nterm = FALSE, Cterm = FALSE, numeric_data = TRUE)[[1L]]$Annotations)
            })
            an <- sweep(an, 1L, p, "*")
            n <- which(apply(an, 2L, \(y) { length(is.all.good(y)) }) == 0L)
            an <- colSums(an, na.rm = TRUE)
            an[n] <- NA
            n <- which(!is.na(an))
            an <- paste0("(", round(an[n], 3L), ")")
            x <- paste(insertElems(s, n-1L, an), collapse = "")
            return(x)
          })
        }
        b2 <- unlist(b1)
        if (length(b2) != length(b1)) { stop("Something went awry!") }
        pep[[j1]] <- ""
        wp <- which(pep$"Modified sequence" %in% a)
        pep[wp,j1] <- b2[match(pep$"Modified sequence"[wp], a)]
      }
      j1 <- paste0(j, " site IDs")
      pep[[j1]] <- NA
      w <- which(!is.na(ev[[j1]]))
      if (length(w)) {
        e <- ev[w,]
        temp <- aggregate(as.character(e[[j1]]), list(e$"Modified sequence"), \(x) {
          paste(sort(unlist(strsplit(x, ";"))), collapse = ";")
        })
        w1 <- which(pep$"Modified sequence" %in% e$"Modified sequence")
        pep[w1, j1] <- temp$x[match(pep$"Modified sequence"[w1], temp$Group.1)]
      }
    }
  }
}

pep[, c("Reverse", "Potential contaminant")] <- ev[rvmtch2, c("Reverse", "Potential contaminant")]
if (scrptType == "withReps") {
  if ((Param$Norma.Pep.Intens)||(Param$Norma.Pep.Ratio)) {
    if (Param$Norma.Ev.Intens) {
      tst <- aggregate(ev$"Normalisation group", list(ev$"Modified sequence"), unique)
      if (is.character(tst$x)) { # It could be that the same sequence is in different normalization groups - if so we do not want to re-use PSM-level normalisation groups!
        pep$"Normalisation group" <- ev$"Normalisation group"[rvmtch2]
      } else { pep$"Normalisation group" <- "Standard" }
    } else { pep$"Normalisation group" <- "Standard" }
    pep$MQ.Exp <- ev$MQ.Exp[match(pep$"Modified sequence", ev$"Modified sequence")]
    nms <- Norm.Groups$names
    w <- which(!nms %in% colnames(pep))
    if (length(w)) { pep[, nms[w]] <- tmp_EM[match(pep$MQ.Exp, tmp_EM$MQ.Exp), nms[w]] }
    pep$"Normalisation group" <- do.call(paste, c(pep[, c(nms, "Normalisation group")], sep = "_"))
  }
}
if ("Quantity Quality" %in% colnames(ev)) { # DiaNN specific:
  tmp <- data.table(Qual = ev$"Quantity Quality", ModSeq = ev$"Modified sequence")
  tmp <- tmp[, list(Qual= mean(Qual)), by = list(ModSeq)]
  tmp <- as.data.frame(tmp)
  pep$"Quantity Quality" <- tmp$Qual[match(pep$"Modified sequence", tmp$ModSeq)]
}
# Interesting note by Vadim on the balance between using more noisy peptides for quant vs fewer high quality ones:
#   https://github.com/vdemichev/DiaNN/issues/1102
# As I expected, they find out that using more is better than filtering.
