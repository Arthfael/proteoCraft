
# Classic normalisations which just shift each vector to match a given level
#
if (normSequence[[nrmStp]]$Method == "mode") {
  pack <- "modeest"
  cran_req <- unique(c(cran_req, pack))
  if (!require(pack, character.only = TRUE, quietly = TRUE)) { pak::pkg_install(pack) }
  require(pack, character.only = TRUE)
}
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
m1 <- normFun(unlist(tmpDat1[wAG1, allSamples])) # Keep original global scale
tmpDat2 <- tmpDat1[, allSamples]*NA
normFlt <- tmpDat1$id
if (normSequence[[nrmStp]]$Method == "GO terms") {
  Norma.Prot.Ratio.to.GO %<o% unlist(strsplit(Param$Norma.Prot.Ratio.to.GO, ";"))
  if (length(Norma.Prot.Ratio.to.GO)) {
    tmp <- proteoCraft::listMelt(strsplit(tmpDat1$Proteins, ";"), tmpDat1$id)
    tmp$"GO-ID" <- db$"GO-ID"[match(tmp$value, db$`Protein ID`)]
    tmp <- tmp[which(!is.na(tmp$"GO-ID")),]
    tmp <- proteoCraft::listMelt(strsplit(tmp$"GO-ID", ";"), tmp$L1)
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
    tmp <- proteoCraft::listMelt(strsplit(tmpDat1$Proteins, ";"), tmpDat1$id)
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
    if (l == 1) { tmp <- Modifs$"Full name"[wbiot] } else {
      tmp <- paste0(paste(Modifs$"Full name"[wbiot[1:(l-1)]], collapse = "\", \""), "\" and \"", Modifs$"Full name"[wbiot[l]])
    }
    warning(paste0("Modifications \"", tmp, "\" were detected as biotinylations, check that this is correct!"))
    normFlt <- tmpDat1$id[grep(topattern(Modifs$Mark[wbiot], start = FALSE), tmpDat1$"Modified sequence")]
  } else {
    Outcome <- FALSE
  }
}
Outcome <- length(normFlt) > 0
if (Outcome) {
  tstNorm <- try({
    for (lGrp in NormGrps$Group) { #lGrp <- NormGrps$Group[1] # Longitudinal group (peptide class)
      grpMtch <- match(NormGrps$IDs[[match(lGrp, NormGrps$Group)]],
                       tmpDat1$id[wAG1])
      grpMtch <- grpMtch[which(!is.na(grpMtch))]
      grpMtch2 <- grpMtch[which(tmpDat1$id[grpMtch] %in% normFlt)]
      stopifnot(length(grpMtch2) > 0)
      for (wGrp in RG$values) { #wGrp <- RG$values[1] # Transversal group (ratios group, i.e. comparison group)
        smpls <- unique(Exp.map$Ref.Sample.Aggregate[which(Exp.map[[RG$column]] == wGrp)])
        smpls <- smpls[which(smpls %in% colnames(tmpDat1))]
        if (length(smpls)) {
          m2 <- normFun(proteoCraft::is.all.good(unlist(tmpDat1[grpMtch, smpls]))) # Keep group's original scale
          if (normSequence[[nrmStp]]$Method == "Levenberg-Marquardt") {
            tmp_nrm <- proteoCraft::AdvNorm.IL(tmpDat1[grpMtch2, ], "Modified sequence", smpls, TRUE, 5)
            colnames(tmp_nrm) <- gsub(topattern("AdvNorm."), "", colnames(tmp_nrm))
            m <- tmpDat1[grpMtch, smpls] - tmp_nrm[, smpls]
            m <- colMeans(m, na.rm = TRUE)
          } else {
            m <- sapply(smpls, function(smpl) { normFun(proteoCraft::is.all.good(tmpDat1[grpMtch2, smpl])) })
          }
          #m <- m-mean(m)
          tmpDat2[grpMtch, smpls] <- sweep(tmpDat1[grpMtch, smpls], 1, m, "-")
          m3 <- normFun(proteoCraft::is.all.good(unlist(tmpDat2[grpMtch, smpls]))) # Keep group's original scale
          tmpDat2[grpMtch, smpls] <- tmpDat2[grpMtch, smpls] + (m2 - m3)
        }
      }
    }
    m4 <- normFun(unlist(tmpDat2[wAG1, ])) # Keep original global scale
    tmpDat2 <- tmpDat2 + (m1 - m4)
    if (normSequence[[nrmStp]]$Method == "Levenberg-Marquardt") {
      txt2 <- paste0("normalized using the ", normSequence[[nrmStp]]$Method, " procedure")
    } else {
      txt2 <- paste0("normalized to the ", normSequence[[nrmStp]]$Method)
    }
  }, silent = TRUE)
  Outcome <- !"try-error" %in% class(tstNorm)
}

