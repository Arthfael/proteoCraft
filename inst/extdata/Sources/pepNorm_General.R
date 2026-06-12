
# Classic normalisations which just shift each vector to match a given level
#
# if (normSequence[[nrmStp]]$Method == "mode") {
#   pack <- "modeest"
#   cran_req <- unique(c(cran_req, pack))
#   if (!require(pack, character.only = TRUE, quietly = TRUE)) { pak::pak(pack) }
#   require(pack, character.only = TRUE)
# }
#
# Create averaging function
if ("funCall" %in% names(normSequence[[prevStp]])) { eval(parse(text = normSequence[[prevStp]]$funCall)) }
#
# Initial values for outputs
tmpDat2 <- NA
wAG2 <- wAG1
Outcome <- TRUE
txt2 <- ""
#
currSamples <- intersect(allSamples, colnames(tmpDat1))
#View(tmpDat1[wAG1, currSamples])
#
# NB:
# To preserve original global and between groups data scale, we calculate:
# - m1: original global mean
# - m2: original group mean
# - m3: post-normalisation group mean
# - m4: post-normalisation global mean
# and apply:
# - after group normalisation: + (m2 - m3)
# - globally: + (m1 - m4)
# Those 4 values used to be calculated with the median, however using the mean here is actually better mathematically.
#
m1 <- mean(unlist(tmpDat1[wAG1, currSamples]), na.rm = TRUE) # Original global scale
tmpDat2 <- tmpDat1[, currSamples]*NA
normFlt <- tmpDat1$id
if (normSequence[[nrmStp]]$Method == "GO terms") {
  Norma.Prot.Ratio.to.GO %<o% unlist(strsplit(Param$Norma.Prot.Ratio.to.GO, ";"))
  if (length(Norma.Prot.Ratio.to.GO)) {
    tmp <- listMelt(strsplit(tmpDat1$Proteins, ";"), tmpDat1$id)
    tmp$"GO-ID" <- db$"GO-ID"[match(tmp$value, db$`Protein ID`)]
    tmp <- tmp[which(!is.na(tmp$"GO-ID")),]
    tmp <- listMelt(strsplit(tmp$"GO-ID", ";"), tmp$L1)
    tmp <- as.data.table(tmp)
    tmp <- tmp[which(tmp$value %in% Norma.Prot.Ratio.to.GO),]
    normFlt <- tmpDat1$id[unique(tmp$L1)]
  } else {
    Outcome <- FALSE
  }
}
if (normSequence[[nrmStp]]$Method == "proteins") {
  Prot.Ratio.ref.Acc %<o% unique(unlist(strsplit(Param$Norma.Prot.Ratio.to.proteins, ";")))
  Prot.Ratio.ref.Acc <- Prot.Ratio.ref.Acc[which(Prot.Ratio.ref.Acc %in% db$"Protein ID")]
  if (length(Prot.Ratio.ref.Acc)) {
    tmp <- listMelt(strsplit(tmpDat1$Proteins, ";"), tmpDat1$id)
    tmp <- tmp[which(tmp$value %in% Prot.Ratio.ref.Acc),]
    normFlt <- tmpDat1$id[unique(tmp$L1)]
  } else {
    Outcome <- FALSE
  }
}
if (normSequence[[nrmStp]]$Method == "biotinylated proteins") {
  wbiot %<o% grep("biot", Modifs$"Full name", ignore.case = TRUE)
  l <- length(wbiot)
  if (length(wbiot)) {
    tmp <- if (l == 1L) { Modifs$"Full name"[wbiot] } else {
      paste0(paste(Modifs$"Full name"[wbiot[1L:(l-1L)]], collapse = "\", \""), "\" and \"", Modifs$"Full name"[wbiot[l]])
    }
    warning(paste0("Modifications \"", tmp, "\" were detected as biotinylations, check that this is correct!"))
    normFlt <- tmpDat1$id[grep(topattern(Modifs$Mark[wbiot], start = FALSE), tmpDat1$"Modified sequence")]
  } else {
    Outcome <- FALSE
  }
}
if (normSequence[[nrmStp]]$Method == "histones") {
  
  stop("Work in progress, check back later!")
  
}
Outcome <- length(normFlt) > 0L
if (Outcome) {
  tstNorm <- try({
    for (lGrp in NormGrps$Group) { #lGrp <- NormGrps$Group[1L] # Longitudinal group (peptide class)
      grpMtch <- match(NormGrps$IDs[[match(lGrp, NormGrps$Group)]],
                       tmpDat1$id[wAG1])
      grpMtch <- grpMtch[which(!is.na(grpMtch))]
      grpMtch2 <- grpMtch[which(tmpDat1$id[grpMtch] %in% normFlt)]
      stopifnot(length(grpMtch2) > 0L)
      for (wGrp in RG$values) { #wGrp <- RG$values[1L] # Transversal group (ratios group, i.e. comparison group)
        smpls <- unique(Exp.map$Ref.Sample.Aggregate[which(Exp.map[[RG$column]] == wGrp)])
        smpls <- smpls[which(smpls %in% colnames(tmpDat1))]
        if (length(smpls)) {
          m2 <- unlist(tmpDat1[grpMtch, smpls])
          m2 <- mean(m2[which(is.finite(m2))]) # Original group scale
          if (normSequence[[nrmStp]]$Method == "Levenberg-Marquardt") {
            tmp_nrm <- AdvNorm.IL(tmpDat1[grpMtch2, ], "Modified sequence", smpls, TRUE, 5L)
            colnames(tmp_nrm) <- gsub(topattern("AdvNorm."), "", colnames(tmp_nrm))
            m <- tmpDat1[grpMtch2, smpls] - tmp_nrm[, smpls]
            m <- colMeans(m, na.rm = TRUE)
          } else {
            m <- vapply(smpls, \(smpl) {
              x <- tmpDat1[grpMtch2, smpl]
              normFun(x[which(is.finite(x))])
            }, 1)
          }
          #m <- m-mean(m)
          tmpDat2[grpMtch, smpls] <- sweep(tmpDat1[grpMtch, smpls], 2L, m, "-")
          m3 <- unlist(tmpDat2[grpMtch, smpls])
          m3 <- mean(m3[which(is.finite(m3))]) # Posterior group scale
          tmpDat2[grpMtch, smpls] <- tmpDat2[grpMtch, smpls] + (m2 - m3) # -> Preserve group scale
        }
      }
    }
    m4 <- mean(unlist(tmpDat2[wAG1, ]), na.rm = TRUE) # Posterior global scale
    tmpDat2 <- tmpDat2 + (m1 - m4) # -> Preserve global scale
    txt2 <- if (normSequence[[nrmStp]]$Method == "Levenberg-Marquardt") {
      paste0("normalized using the ", normSequence[[nrmStp]]$Method, " procedure")
    } else {
      paste0("normalized to the ", normSequence[[nrmStp]]$Method)
    }
  }, silent = TRUE)
  Outcome <- !inherits(tstNorm, "try-error")
}

