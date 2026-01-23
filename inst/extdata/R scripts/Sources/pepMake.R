# Peptides: create pep object
require(data.table)
source(parSrc, local = FALSE)
clusterExport(parClust, list("wd", "scrptType"), envir = environment())
if (scrptType == "withReps") {
  tmp_EM <- Exp.map
  tmp_EM[which(colnames(tmp_EM) == "Ref.Sample.Aggregate")] <- "Parent_sample"
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
readr::write_rds(tmp_EM, paste0(wd, "/tmp1.RDS"))
readr::write_rds(tmp2, paste0(wd, "/tmp2.RDS"))
invisible(clusterCall(parClust, function() {
  library(data.table)
  tmp_EM <<- readr::read_rds(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readr::read_rds(paste0(wd, "/tmp2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
myGrps <- "ALLMYSAMPLESTUDUDUDUMMMDADA"
if (scrptType == "withReps") {
  myGrps <- c(myGrps, RSA$values)
}
if (scrptType == "noReps") {
  myGrps <- c(myGrps, Exp)
}
tmp4 <- setNames(parLapply(parClust, myGrps, function(i) { #i <- myGrps[1] #i <- myGrps[2]
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
for (i in myGrps[2:length(myGrps)]) {
  tmp <- tmp4[[i]]
  ki <- paste0("Evidence IDs - ", i)
  pep[[ki]] <- tmp$IDs[match(pep$"Modified sequence", tmp$ModSeq)]
  pep[which(is.na(pep[[ki]])), ki] <- ""
}
pep$id <- 1:nrow(pep)
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
tmp <- parSapply(parClust, proteoCraft::AA, function(aa) {
  nchar(sq) - nchar(gsub(aa, "", sq))
})
colnames(tmp) <- paste0(colnames(tmp), " Count")
pep[, colnames(tmp)] <- tmp
ev[, paste0(AA, " Count")] <- pep[mtch, paste0(AA, " Count")]
for (aa in c("O", "U")) { # Only keep the selenocysteine and pyrrolysine amino acid columns if they are non-empty (NB: pyrrolysine should really only be in some bacteria)
  if (sum(ev[[paste0(aa, " Count")]]) == 0) {
    ev[[paste0(aa, " Count")]] <- NULL
    pep[[paste0(aa, " Count")]] <- NULL
  }
}
#
pep$Length <- nchar(pep$Sequence)
ev$Length <- pep$Length[mtch]
tmp <- data.table(mod = ev$"Modified sequence", Intensity = ev$Intensity)
tmp$Intensity[which(!is.all.good(tmp$Intensity, 2))] <- NA
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
w <- which(nchar(tmp$Proteins) > 13)
tmp$Proteins[w] <- paste0(substr(tmp$Proteins[w], 1, 12), "...")
tmp$"Modified sequence" <- gsub("_", "", tmp$"Modified sequence")
pep$"Peptide name" <- do.call(paste, c(tmp, sep = "\n"))
ev$"Peptide name" <- pep$"Peptide name"[mtch2]

# Modification site probabilities, score differences, site IDs, etc...
if (exists("Modifs")) {
  wPr <- which(paste0(Modifs$"Full name", " Probabilities") %in% colnames(ev))
  if (length(wPr)) {
    for (i in wPr) { #i <- wPr[1]
      for (j in c(" Probabilities", " Score Diffs")) { #j <- " Probabilities"
        j1 <- paste0(Modifs$"Full name"[i], j)
        temp <- ev[, c("Modified sequence", j1, "PEP")]
        temp <- temp[which(temp[[j1]] != ""),]
        a0 <- unique(temp$"Modified sequence")
        clusterExport(parClust, list("temp", "j1"), envir = environment())
        b1 <- parLapply(parClust, a0, function(x) { temp[which(temp$"Modified sequence" == x), j1] })
        b2 <- parLapply(parClust, a0, function(x) { temp$PEP[which(temp$"Modified sequence" == x)] })
        l <- vapply(b1, function(x) { length(unique(x)) }, 1)
        wb <- which(l == 1)
        if (length(wb)) { b1[wb] <- sapply(b1[wb], function(x) { unique(x) }) } 
        wb <- which(l > 1)
        if (length(wb)) {
          temp1 <- cbind(b1[wb], b2[wb])
          clusterExport(parClust, "temp1", envir = environment())
          b1[wb] <- parApply(parClust, temp1, 1, function(x) {
            p <- 1-x[[2]]
            p <- p/sum(p)
            x <- x[[1]]
            s <- unlist(strsplit(unique(gsub("\\([^\\)]+\\)", "", x)), ""))
            an <- proteoCraft::Isapply(x, function(y) {
              as.numeric(proteoCraft::annot_to_tabl(y, Nterm = FALSE, Cterm = FALSE, numeric_data = TRUE)[[1]]$Annotations)
            })
            an <- sweep(an, 1, p, "*")
            n <- which(apply(an, 2, function(y) { length(proteoCraft::is.all.good(y)) }) == 0)
            an <- colSums(an, na.rm = TRUE)
            an[n] <- NA
            n <- which(!is.na(an))
            an <- paste0("(", round(an[n], 3), ")")
            x <- paste(proteoCraft::insertElems(s, n-1, an), collapse = "")
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
        temp <- aggregate(as.character(e[[j1]]), list(e$"Modified sequence"), function(x) {
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
      if ("character" %in% class(tst$x)) { # It could be that the same sequence is in different normalization groups - if so we do not want to re-use PSM-level normalisation groups!
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
if ("Quantity Quality" %in% colnames(ev)) { # Dia-NN specific:
  tmp <- data.table(Qual = ev$"Quantity Quality", ModSeq = ev$"Modified sequence")
  tmp <- tmp[, list(Qual= mean(Qual)), by = list(ModSeq)]
  tmp <- as.data.frame(tmp)
  pep$"Quantity Quality" <- tmp$Qual[match(pep$"Modified sequence", tmp$ModSeq)]
}
# Interesting note by Vadim on the balance between using more noisy peptides for quant vs fewer high quality ones:
#   https://github.com/vdemichev/DiaNN/issues/1102
# As I expected, they find out that using more is better than filtering.
