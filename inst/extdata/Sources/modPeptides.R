#### Code chunk - Modified peptides analysis
if (scrptTypeFull == "withReps_PG_and_PTMs") {
  myIDcol <- "Leading proteins"
}
if (scrptTypeFull == "withReps_PTMs_only") {
  myIDcol <- "Proteins"
  PTM_normalize %<o% FALSE
}
#
ptms.PVal <- pvalue.col[which(pvalue.use)]
if (names(ptms.PVal) == "DEqMS") {
  warning("Currently DEqMS doesn't work for PTM-modified peptides, using standard limma instead")
  ptms.PVal <- pvalue.col["Moderated"]
}
ptms.ratios.ref %<o% setNames("log2(Ratio) - ",
                              "Original")
PepLabKol %<o% setNames(c("Mod. sequence", "Proteins", "Common protein names", "PEP"),
                        c("Mod. sequence", "Protein(s)", "Protein name(s)", "PEP"))
pep$"Mod. sequence" <- gsub("^_|_$", "", pep$"Modified sequence")
if (("Phospho.analysis" %in% colnames(Param))&&(Param$Phospho.analysis)) {
  warning("This argument is deprecated - use the \"PTM.analysis\" argument now.\nI understand you want to perform phosphopeptides analysis though so will be doing the editing for you... this time!\nCarry on please...")
  if ("PTM.analysis" %in% colnames(Param)) {
    a <- toupper(unlist(strsplit(Param$PTM.analysis, ";")))
    if (!"PHOSPHO" %in% a) { Param$PTM.analysis <- paste(c(a, "PHOSPHO"), collapse = ";") }
  } else { Param$PTM.analysis <- "PHOSPHO" }
}
PTMs %<o% if ("PTM.analysis" %in% colnames(Param)) { PTMs <- unlist(strsplit(Param$PTM.analysis, ";")) } else { c() }
if (length(PTMs)) {
  msg <- "Modified peptides analysis"
  ReportCalls <- AddMsg2Report(Space = FALSE)
  PTMs_ref.ratios %<o% list()
  PTMs_FDR.thresholds %<o% list()
  PTMs_pep %<o% list()
  PTMs_Reg_filters %<o% list()
  PTMs_int.ref %<o% list()
  PTMs_rat.ref %<o% list()
  PTMs_SAM_thresh %<o% list()
  if (F.test) {
    PTMs_F_test_data %<o% list()
    #PTMs_F_test_ref_ratios %<o% list() # Not needed
  }
  if (("PTM.analysis_Norm" %in% colnames(Param))&&(Param$PTM.analysis_Norm != "")) {
    tmp <- as.logical(unlist(strsplit(as.character(Param$PTM.analysis_Norm), ";")))
    if ((!length(tmp) %in% c(1L, length(PTMs)))||(NA %in% tmp)) {
      warning("I cannot make sense of parameter PTM.analysis_Norm, defaulting to TRUE")
      tmp <- TRUE
    }
    if ((length(tmp) == 1L)&&(length(PTMs) > 1L)) {
      warning(paste0("Recycling PTM.analysis_Norm (value = ", as.character(tmp), ") over PTMs..."))
      tmp <- rep(tmp, length(PTMs))
    }
    PTM_normalize %<o% setNames(tmp, PTMs)
  } else { PTM_normalize %<o% setNames(rep(TRUE, length(PTMs)), PTMs) }
  if (enrichGO) {
    PTMs_GO_Plots %<o% list()
    PTMs_Reg_GO_terms %<o% list()
    PTMs_GO_enrich.dat %<o% list()
    PTMs_GO_enrich.FCRt %<o% list()
    PTMs_GO_enrich.tbl %<o% list()
    PTMs_GO_Plots %<o% list()
    PTMs_Reg_GO_terms %<o% list()
    ### Check that Cytoscape is installed and can run, then launch it.
    Src <- paste0(libPath, "/extdata/Sources/Cytoscape_init.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    # Initialize ClueGO
    Src <- paste0(libPath, "/extdata/Sources/ClueGO_init.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
  }
  #
  # Edit materials and methods text
  tmp <- PTMs
  lT <- length(tmp)
  if (lT > 1L) { tmp <- paste0(paste(tmp[seq_len(lT-1L)], collapse = ", "), " and ", tmp[lT]) }
  tmp2 <- setNames(sapply(PTMs, \(Ptm) {
    PTM_normalize[[Ptm]]
  }), PTMs)
  tmp2u <- unique(tmp2)
  l2u <- length(tmp2u)
  if ((l2u == 1L)&&(tmp2u)) {
    tmp2 <- ", re-normalizing values to account for average parent protein group(s) fold change."
  } else {
    if (l2u > 1L) {
      tmp2 <- names(tmp2)[which(tmp2)]; l2 <- length(tmp2)
      if (l2 > 1L) { tmp2 <- paste0(paste(tmp2[seq_len(lT-1L)], collapse = ", "), " and ", tmp2[lT], "") }
      tmp2 <- ", normalizing values to correct for parent protein group(s) fold change."
    } else { tmp2 <- "." }
  }
  l <- length(DatAnalysisTxt)
  DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l], " Statistical analysis was performed for ",
                              c("", "PTMs ")[(lT > 1L)+1L], tmp,
                              "-modified peptides as for protein groups", tmp2)
  #
  ptmNms <- setNames(vapply(PTMs, \(ptm) { #ptm <- PTMs[1L]
    a <- unlist(strsplit(gsub("\\)$", "", ptm), "\\("))
    Ptm <- if (length(a) > 1L) {
      paste0(toupper(substr(a[1L], 1L, 1L)), substr(a[1L], 2L, nchar(ptm)), "(", a[2L], ")")
    } else {
      paste0(toupper(substr(a[1L], 1L, 1L)), substr(a[1L], 2L, nchar(ptm)))
    }
    return(Ptm)
  }, ""), PTMs)
  ptmsTst <- setNames(lapply(PTMs, \(ptm) { #ptm <- PTMs[1L]
    Ptm <- ptmNms[[ptm]]
    w <- which(Modifs$"Full name" == Ptm)
    if (length(w) != 1L) { return(0) }
    p <- Modifs$Mark[w]
    ppat <- paste0("\\(", p, "\\)|\\(", p, ",|,", p, "\\)|,", p, ",") # Pattern to catch all instances of the mod
    tmp <- grepl(ppat, pep$"Modified sequence")
    g <- which(tmp)
    return(g)
  }), PTMs)
  ptmsTst2 <- lengths(ptmsTst)
  w0 <- which(ptmsTst2 == 0L)
  lW0 <- length(w0)
  if (lW0) {
    msg <- ptmNms[w0]
    if (lW0 > 2L) { msg <- c(paste0(msg[1L:(lW0-1L)], collapse = ", "), msg) }
    msg <- paste0(" !!! No peptides found for ", paste0(msg, collapse = " and "), ", skipping th", c("is", "ese")[(lW0 > 1L)+1L],
                  " modification", c("", "s")[(lW0 > 1L)+1L])
    warning(msg)
    PTMs <- PTMs[w0]
    ptmNms <- ptmNms[PTMs]
    ptmsTst <- ptmsTst[PTMs]
  }
  for (ptm in PTMs) { #ptm <- PTMs[1L]
    pep[[ptm]] <- 1L:nrow(pep) %in% ptmsTst[[ptm]]
  }
  ptms.ref %<o% pep.ref[length(pep.ref)]
  # ptms.ratios.ref %<o% setNames(pep.ratios.ref[length(pep.ratios.ref)],
  #                               "Original")
  for (ptm in PTMs) { #ptm <- PTMs[1L]
    source(parSrc, local = FALSE)
    ReportCalls <- AddMsg2Report(Msg = paste0(" - ", ptm), Space = FALSE)
    modDirs <- c("", "/t-tests")
    if (F.test) { modDirs <- c(modDirs, "/F-tests") }
    modDirs <-  paste0(wd, "/Reg. analysis/", ptm, modDirs)
    for (dr in modDirs) { if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }}
    dirlist <- union(dirlist, modDirs)
    #
    Ptm <- ptmNms[ptm]
    g <- ptmsTst[[ptm]]
    w <- which(Modifs$"Full name" == Ptm)
    p <- Modifs$Mark[w]
    ppat <- paste0("\\(", p, "\\)|\\(", p, ",|,", p, "\\)|,", p, ",") # Pattern to catch all instances of the mod
    ptmsh <- substr(p, 1L, 1L)
    ptmpep <- pep[ptmsTst[[ptm]], ]
    #pep[[paste0(Ptm, " ID")]] <- ""
    #ptmpep <- pep[g,]
    #pep[g, paste0(Ptm, " ID")] <- ptmpep$ModPep_ID <- seq_len(nrow(ptmpep))
    temp <- ptmpep[, c("Modified sequence", myIDcol)]
    temp[[myIDcol]] <- strsplit(temp[[myIDcol]], ";")
    temp$"Modified sequence" <- gsub(paste0("[^A-Z", ptmsh, "]"), "",
                                     gsub(ppat, ptmsh, temp$"Modified sequence"))
    #ptmpep[, c("Match(es)", paste0(Ptm, "-site(s)"))] <- ""
    dbsmall <- db[which(db$"Protein ID" %in% unique(unlist(temp[[myIDcol]]))), c("Protein ID", "Sequence")]
    # On I/L ambiguity remaining even with newer DIA methods taking into account RT, IM and fragments intensity, see https://github.com/vdemichev/DiaNN/discussions/1631
    dbsmall$"Seq*" <- gsub("I", "L", dbsmall$Sequence)
    temp$"ModSeq*" <- gsub("I", "L", temp$"Modified sequence")
    #
    kol <- c(myIDcol, "ModSeq*")
    temp$`ModSeq*` <- strsplit(temp$`ModSeq*`, "")
    temp2 <- temp[, kol]
    clusterExport(parClust, list("temp2", "dbsmall", "ptmsh", "listMelt"), envir = environment())
    temp3 <- parApply(parClust, temp2, 1L, \(x) {
      #x <- temp[1, kol]
      #A <- sapply(seq_len(nrow(temp)), \(i) {
      #print(i)
      #x <- temp[i, kol]
      m <- unlist(x[[2L]])
      m <- data.frame(Seq = m, Mod.seq = m, Test = FALSE)
      w1 <- which(m$Seq == ptmsh)
      w2 <- which(m$Seq != ptmsh)
      m$Mod.seq[w1-1L] <- paste0(ptmsh, m$Mod.seq[w1-1L])
      m$Test[w1-1L] <- TRUE
      m <- m[w2,]
      l <- nrow(m)
      m$Offset <- 0L:(l-1L)
      q <- unlist(x[[1L]])
      mtch <- match(q, dbsmall$"Protein ID")
      wN <- which(!is.na(mtch))
      mtch <- mtch[wN]
      q <- dbsmall$"Protein ID"[mtch[wN]]
      if (length(mtch)) {
        seq <- strsplit(dbsmall$"Seq*"[mtch], "")
        matches <- lapply(seq, \(S) { #S <- seq[1L]
          S <- unlist(S)
          lS <- length(S)
          m1 <- m
          m1$Match <- apply(m1[, c("Seq", "Offset")], 1L, \(y) { which(S == y[1L]) - as.numeric(y[2L]) })
          M <- unlist(m1$Match)
          M <- M[which(M > 0L)]
          M <- aggregate(M, list(M), length)
          M <- M[order(-M$x),]
          M <- M$Group.1[which(M$x == l)]
          # Check that peptides are tryptic:
          #test <- sapply(M, \(y) {
          #  # r1: on the N-terminal end, is the peptide preceded by K, R or (if starting at position 2, M)?
          #  if (y > 1L) {
          #    r1 <- if (y == 2L) { S[y-1L] %in% c("K", "R", "M") } else { S[y-1L] %in% c("K", "R") }
          #  } else { r1 <- TRUE }
          #  # r2: on the C-terminal end, is this a tryptic peptide or the last peptide in the protein?
          #  r2 <- (m1$Seq[l] %in% c("K", "R"))|(y+l-1L == lS)
          #  return(r1+r2 == 2L)
          #})
          #M <- M[which(test)]
          return(M)
        })
        names(matches) <- q
        matches <- matches[which(lengths(matches) > 0L)]
        if (length(matches)) {
          matches <- listMelt(matches, ColNames = c("Match", "Protein"))
          matches <- aggregate(matches$Protein, list(matches$Match), paste, collapse = ";")
          colnames(matches) <- c("Match", "Proteins")
          w <- which(m$Test)
          matches$Sites <- sapply(matches$Match, \(y) {
            y <- paste(sapply(w, \(z) { paste0(m$Mod.seq[z], y+m$Offset[z]) }), collapse = "-")
          })
          matches$Match <- apply(matches[, c("Match", "Proteins")], 1L, paste, collapse = " ")
          matches$Sites <- apply(matches[, c("Sites", "Proteins")], 1L, paste, collapse = " ")
          matches <- apply(matches[, c("Match", "Sites")], 2L, paste, collapse = "/")
        } else { matches <- c(NA, NA) }
      } else { matches <- c(NA, NA) }
      return(matches)
      #})
    })
    ptmpep[, c("Match(es)", paste0(Ptm, "-site(s)"))] <- as.data.frame(t(temp3))
    ptmpep[[paste0(Ptm, "-site")]] <- gsub(" .+", "", ptmpep[[paste0(Ptm, "-site(s)")]])
    ptmpep <- ptmpep[which(!is.na(ptmpep$`Match(es)`)),]
    ptmpep$tmp1 <- gsub("^_|_$", "", ptmpep$"Modified sequence")
    ptmpep$tmp2 <- ptmpep[[paste0(Ptm, "-site(s)")]]
    nc <- nchar(ptmpep$tmp2)
    w <- which(nc > 25L)
    ptmpep$tmp2[w] <- paste0(substr(ptmpep$tmp2[w], 1L, 22L), "...")
    ptmpep$Code <- apply(ptmpep[, paste0("tmp", 1L:2L)], 1L, paste, collapse = "\n")
    ptmpep$tmp1 <- NULL
    ptmpep$tmp2 <- NULL
    ptmpep$Name <- ""
    w <- which(ptmpep[[myIDcol]] != "")
    temp2 <- as.data.frame(t(sapply(strsplit(gsub("[/,;].+$", "", ptmpep[w, paste0(Ptm, "-site(s)")]), " "), unlist)))
    colnames(temp2) <- c("Site", "Protein")
    temp2$Protein <- db$"Common Name"[match(temp2$Protein, db$"Protein ID")]
    temp2$ModSeq <- gsub("^_|_$", "", ptmpep$`Modified sequence`[w])
    ptmpep$Name[w] <- do.call(paste, c(temp2[, c("Site", "ModSeq", "Protein")], sep = "\n"))
    ptmpep$Name[which(ptmpep$Name == "")] <- paste0("Unknown source ", ptm, "-modified peptide #", seq_along(which(ptmpep$Name == "")))
    #View(ptmpep[,c("Match(es)", "Modified sequence", "Code", paste0(Ptm, "-site(s)"))])
    if (grepl("^[P,p]hospho( \\([A-Z]+\\))?$", ptm)) {
      p_col <- paste0(gsub(" |\\(|\\)", ".", ptm), ".Probabilities")
      scd_col <- paste0(gsub(" |\\(|\\)", ".", ptm), ".Score.Diffs")
      if (sum(c(p_col, scd_col) %in% colnames(ptmpep)) == 2L) {
        temp <- try(phos_QC(ptmpep, P_col = p_col, ScD_col = scd_col), silent = TRUE)
        if (inherits(temp, "try-error")) {
          warning("No phospho QC performed, check colnames or the phos_QC function!")
        } else { ptmpep$High_Quality_ptmpep <- temp }
      }
    }
    # Convert peptides intensities to log10
    # (I've had it! I should've done this way earlier)
    ref <- ptms.ref[1L]
    ref2 <- paste0("log10(", gsub(" - $", ") - ", gsub("Evidence intensities", "PSMs. int.", ref)))
    kols <- grep(topattern(ref), colnames(ptmpep), value = TRUE)
    for (kol in kols) {
      smpl <- gsub(topattern(ref), "", kol)
      ptmpep[[paste0(ref2, smpl)]] <- log10(ptmpep[[kol]])
      ptmpep[[kol]] <- NULL
    }
    ptms.ref <- c("Original" = ref2)
    #
    # Optional: normalize to parent protein group
    if ((scrptTypeFull == "withReps_PG_and_PTMs")&&(PTM_normalize[[Ptm]])) {
      # Step 1: normalize ratios:
      # We essentially want to correct the fold change of each modified peptide by that of the parent protein
      #
      #
      #
      ################################
      #         IMPORTANT!!!         #
      ################################
      #
      #  This can only work if any protein groups-level re-normalisation is back-propagated onto peptides!!!
      #  This is currently the case. Make sure it stays so!
      #
      #
      #
      #
      # Complete re-write on 16-04-2026 after the switch to a contrasts-centric approach (where any group may be a)
      # from a ctrl vs .
      kolPp1 <- paste0(ptms.ref, RSA$values)
      kolPG <- paste0(Prot.Expr.Root, RSA$values)
      wXst <- which((kolPp1 %in% colnames(ptmpep))&(kolPG %in% colnames(PG)))
      kolPp1 <- kolPp1[wXst]
      kolPG <- kolPG[wXst]
      grps <- Exp.map[match(RSA$values[wXst], Exp.map$Ref.Sample.Aggregate), VPAL$column]
      pep2prot <- listMelt(strsplit(ptmpep$`Protein group IDs`, ";"), 1L:nrow(ptmpep), c("id", "row"))
      pep2prot$id <- as.integer(pep2prot$id)
      prtQuntDat <- PG[, c("id", kolPG)] # log-base = 10...
      if (NAsReplMeth == 1L) { # For method 1 (imputation) we should do it early...
        tmp <- Data_Impute2(prtQuntDat, grps)
        prtQuntDat[, kolPG] <- tmp$Imputed_data[, kolPG]
      }
      prtQuntDat[, kolPG] <- sweep(prtQuntDat[, kolPG], 1L, rowMeans(prtQuntDat[, kolPG], na.rm = TRUE), "-")
      prtQuntDat <- prtQuntDat[match(pep2prot$id, prtQuntDat$id), kolPG]
      if (NAsReplMeth == 2L) { # ... for method 2 (median) late is better (after sweeping by row means, i.e. converting as log10FC)
        tmp <- as.matrix(prtQuntDat)
        w <- which(is.na(tmp))
        tmp[w] <- median(tmp, na.rm = TRUE)
        prtQuntDat <- tmp
      }
      pep2prot[, kolPG] <- prtQuntDat[, kolPG]
      pep2prot <- data.table(pep2prot[, c("row", kolPG)])
      pep2prot <- pep2prot[, lapply(.SD, mean, na.rm = TRUE), by = row, .SDcols = kolPG]
      pep2prot <- as.data.frame(pep2prot)
      # maximum PG fold change should not be too high in most experiments
      #summary(10^unlist(pep2prot[, kolPG]))   # de-logged
      #summary(unlist(pep2prot[, kolPG]))      # as log10
      # (refer to PG-level volcano plots to see whether these make sense...)
      ptms.ref["ReNorm."] <- paste0("ReNorm. ", ptms.ref["Original"])
      tmp <- ptmpep[, kolPp1] - pep2prot[, kolPG]
      kolPp2 <- paste0(ptms.ref["ReNorm."], RSA$values[wXst])
      ptmpep[, kolPp2] <- tmp
      #
      df1 <- ptmpep[, kolPp1]
      df2 <- ptmpep[, kolPp2]
      tst <- length(is.all.good(unlist(df1))) == length(is.all.good(unlist(df2)))
      if (!tst) {
        warning(paste0("Are you expecting re-normalisation of ", ptm,
                       "-modified peptides to results in losses of valid intensity values? If not, investigate, because it happened!"))
      }
      pepPlotFun(df1,
                 df2,
                 "Re-normalisation (intensities)",
                 modDirs[1L])
      if (length(prot.list)) {
        source(parSrc, local = FALSE)
        dr <- paste0(modDirs[1L], "/Heatmaps")
        pepHtmp(prot.list,
                ptmpep,
                ptms.ref["Original"],
                dr,
                paste0(ptm, "-mod. pept. log2 heatmap"),
                "original",
                is.log = TRUE,
                cl = parClust)
        pepHtmp(prot.list,
                ptmpep,
                ptms.ref["ReNorm."],
                dr,
                paste0(ptm, "-mod. pept. log2 heatmap"),
                "re-normalised",
                is.log = TRUE,
                cl = parClust)
        #openwd(dr)
      }
    }
    pepRf <- ptms.ref[length(ptms.ref)]
    # Now calculate ratios
    tmp <- make_Rat2(ptmpep,
                     experiment.map = Exp.map,
                     int.root = pepRf,
                     rat.root = ptms.ratios.ref)
    ptmpep[, colnames(tmp)] <- tmp
    #
    dataType <- "modPeptides" # used for PCA and stats
    #
    Src <- paste0(libPath, "/extdata/Sources/dimRed_plots.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    #View(ptmpep[, grep(topattern(pepRf), colnames(ptmpep), value = TRUE)])
    #
    # Calculate average intensities and ratios, as well as Welch's t-test and moderated P-values;
    # For unpaired replicates a permutations t-test is also performed.
    samDir <- paste0(modDirs[2L], "/SAM")
    ebamDir <- paste0(modDirs[2L], "/EBAM")
    #
    Src <- paste0(libPath, "/extdata/Sources/Stat_tests.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    # Also mean expression over whole dataset
    kls <- grep(topattern(pepRf), colnames(ptmpep), value = TRUE)
    ptmpep$"Mean Expr." <- apply(ptmpep[, kls], 1L, \(x) { mean(is.all.good(unlist(x))) })
    # Create list of control ratio values for the purpose of identifying thresholds for plots:
    #
    if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
      ref.rat <- NULL
    }
    if (Param$Ratios.Thresholds == threshMsg) {
      stop("This option is deprecated!")
      # PTMs_ref.ratios[[Ptm]] <- ref.rat <- setNames(lapply(VPAL$values, \(x) { #x <- VPAL$values[1L]
      #   if (RatConGrps == "Ratio groups") {
      #     x1 <- unique(Exp.map[which(Exp.map[[VPAL$column]] == x), RG$column])
      #   }
      #   if (RatConGrps == "Experiments") {
      #     x1 <- unique(Exp.map$Experiment[which(Exp.map[[VPAL$column]] == x)])
      #     x1 <- unique(Exp.map[which(Exp.map$Experiment == x1), RG$column])
      #   }
      #   if (RatConGrps == "Whole dataset") {
      #     x1 <- unique(Exp.map[[RG$column]])
      #   }
      #   x <- grep(paste0(topattern(paste0(ptms.ratios.ref[length(ptms.ratios.ref)], x1, "_REF.to.REF_")), "[0-9]+"),
      #             colnames(ptmpep), value = TRUE)
      #   x <- if (length(x)) { is.all.good(as.numeric(unlist(ptmpep[, x]))) } else { NULL }
      #   return(x)
      # }), VPAL$values)
    }
    #
    # Estimate P-value significance for a set of accepted FDRs:
    ## NB: For graphical reasons (volcano plots), there is only support for 4 different FDR values. This should suffice anyway.
    temp_thrsh <- c()
    A <- myContrasts$Contrast
    test <- sapply(A, \(x) { #x <- A[6L]
      x <- paste0(ptms.PVal, x)
      r <- x %in% colnames(ptmpep)
      if (r) { r <- length(is.all.good(as.numeric(ptmpep[[x]]))) > 0L }
      return(r)
    })
    A <- A[which(test)]
    stopifnot(length(A) > 0L)
    for (a in A) { #a <- A[1L]
      temp <- FDR(data = ptmpep,
                  aggregate = a,
                  pvalue_root = ptms.PVal,
                  fdr = BH.FDR,
                  returns = c(TRUE, TRUE, FALSE),
                  inputType = "log")
      ptmpep[, colnames(temp$`Significance vector`)] <- temp$`Significance vector`
      temp_thrsh <- c(temp_thrsh, temp$Thresholds)
      PTMs_FDR.thresholds[[Ptm]] <- temp_thrsh
    }
    ptmpep$"1-PEP" <- 1 - ptmpep$PEP
    ptmpep$"log10(1-PEP)" <- log10(ptmpep$"1-PEP")
    a <- grep(topattern(pepRf), colnames(ptmpep), value = TRUE)
    a <- a[which(!grepl("\\.REF$", a))]
    ptmpep$"Av. log10 abundance" <- apply(ptmpep[, a], 1L, \(x) { mean(is.all.good(unlist(x))) })
    ptmpep$"Rel. av. log10 abundance" <- ptmpep$"Av. log10 abundance"/max(is.all.good(ptmpep$"Av. log10 abundance"))
    # Arbitrary thresholds
    P <- Param
    P$Plot.labels <- "Name"
    P$Plot.metrics <- paste0("X:Ratio_Mean.log2;Y:", ptms.PVal)
    ReportCalls <- AddMsg2Report(Msg = " -> ", ptm, " t-tests volcano plots", Space = FALSE)
    #k1 <- grep(topattern(paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)])), colnames(ptmpep), value = TRUE)
    #df1 <- ptmpep[, k1]
    #subDr <- gsub(topattern(wd), "", modDirs[2L])
    subDr <- modDirs[2L]
    if (useSAM) {
      # In this case, we bypass the original decision and base it off SAM even though we plot Student's P-values
      for (i in names(PTMs_SAM_thresh[[Ptm]])) { #i <- names(PTMs_SAM_thresh[[Ptm]])[1L]
        dec <- PTMs_SAM_thresh[[Ptm]][[i]]$decision
        mKol <- rev(colnames(dec))[1L]
        FCkol <- paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)], i)
        stopifnot(FCkol %in% names(ptmpep))
        regKol <- paste0("Regulated - ", i)
        ptmpep[[regKol]] <- "non significant"
        fdrs <- as.numeric(gsub("FDR$", "", colnames(dec)[which(colnames(dec) != mKol)]))
        fdrs <- sort(fdrs, decreasing = TRUE)
        for (f in fdrs) { #f <- fdrs[1L]
          w <- which(ptmpep[[mKol]] %in% dec[which(dec[[paste0(f, "FDR")]] == "+"), mKol])
          if (length(w)) {
            ptmpep[which(ptmpep[w, FCkol] > 0), regKol] <- paste0("up, FDR = ", f*100, "%")
            ptmpep[which(ptmpep[w, FCkol] < 0), regKol] <- paste0("down, FDR = ", f*100, "%")
          }
        }
      }
    }
    stopCluster(parClust)
    source(parSrc)
    volcPlot_args2 <- volcPlot_args
    volcPlot_args2$Prot <- ptmpep
    volcPlot_args2$X.root <- ptms.ratios.ref
    volcPlot_args2$Y.root <- ptms.PVal
    volcPlot_args2$parameters <- P
    volcPlot_args2$FDR.thresh <- PTMs_FDR.thresholds[[Ptm]]
    volcPlot_args2$arbitrary.lines <- arbitrary.thr
    volcPlot_args2$IDs.col <- "Code"
    volcPlot_args2$Proteins.col <- "Proteins"
    volcPlot_args2$title <- paste0(Ptm, " volcano plot ")
    volcPlot_args2$subfolder <- subDr
    volcPlot_args2$plotly_labels <- c(PepLabKol, paste0(Ptm, "-site"))
    volcPlot_args2$curved_Thresh <- PTMs_SAM_thresh[[Ptm]]
    volcPlot_args2$cl <- parClust
    tempVPptm <- do.call(Volcano.plot, volcPlot_args2)
    #k2 <- grep(topattern(paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)])), colnames(ptmpep), value = TRUE)
    #df2 <- ptmpep[, k2]
    #pepPlotFun(df1, df2, "Before VS after volc. plot", modDirs[1L], FALSE)
    stopCluster(parClust)
    source(parSrc)
    #
    # Save plotly plots
    dr <- subDr
    myPlotLys <- tempVPptm$`Plotly plots`
    Src <- paste0(libPath, "/extdata/Sources/save_Plotlys.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    VP_list <- tempVPptm
    insrt <- ""
    Src <- paste0(libPath, "/extdata/Sources/thresholds_Excel.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    #
    ptmpep <- tempVPptm$Protein_groups_file
    volcano.plots[[Ptm]] <- tempVPptm$Plots
    n2 <- names(volcano.plots[[Ptm]]$Labelled)
    for (ttl in n2) {
      plot <- volcano.plots[[Ptm]]$Labelled[[ttl]]
      ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
    }
    if ((create_plotly)&&(!create_plotly_local)) {
      plot_ly[[paste0(Ptm, "_Volcano plots (t-tests)")]] <- tempVPptm$"Plotly plots"
    }
    # Specificity mark for untested proteins
    # We can assume that proteins with PSMs only in the specific pull-down samples are actually specifically enriched!
    # -> Label them as such!
    # NB: This used to be for pull-downs only but I am extending it to the whole workflow for now.
    #     This could also be parameter-controlled.
    # NB: For the F-test, this is done within the function
    #if (IsPullDown) {
    kolTR <- paste0("Regulated - ", myContrasts$Contrast)
    w <- which(kolTR %in% colnames(ptmpep))
    if (length(w)) {
      kolTR <- kolTR[w]
      for (i in w) { #i <- w[1L]
        koleA <- paste0("Evidence IDs - ", myContrasts$A_samples[[i]])
        koleB <- paste0("Evidence IDs - ", myContrasts$B_samples[[i]])
        tmpA <- do.call(cbind, lapply(koleA, \(x) { lengths(strsplit(ptmpep[[x]], ";")) }))
        tmpB <- do.call(cbind, lapply(koleB, \(x) { lengths(strsplit(ptmpep[[x]], ";")) }))
        tstA <- rowSums(tmpA > 0L) == length(koleA)  # = all non-missing
        tstB <- rowSums(tmpB == 0L) == length(koleB) # = all missing
        wA <- which(tstA & tstB)
        if (length(wA)) {
          evcount <- rowSums(tmpA[wA,])
          txtup <- paste0("Specific: at least ", evcount, " PSMs/sample")
          ptmpep[wA, kolTR[i]] <- txtup
        }
      }
    } else { warning("I would expect \"Regulated ...\" columns in the modified peptides table by this stage!") }
    #
    # Create t-test filters:
    ## These can then be used for further steps down the line, such as volcano plots, etc...
    g <- grep("^Regulated - ", colnames(ptmpep), value = TRUE)
    g1 <- gsub("^Regulated - ", "", g)
    up <- grep("^up|^Specific", unique(unlist(ptmpep[, g])), value = TRUE)
    down <- grep("^down|^Anti-specific", unique(unlist(ptmpep[, g])), value = TRUE) # The "anti-specific" part will only become relevant if in future I disconnect symmetry and/or adding specific tags from pull-down experiments
    PTMs_Reg_filters[[Ptm]] <- list()
    PTMs_Reg_filters[[Ptm]]$"t-tests" <- list()
    if ("con" %in% filter_types) {
      PTMs_Reg_filters[[Ptm]]$"t-tests"$"By condition" <- list()
      for (i in seq_along(g)) {
        PTMs_Reg_filters[[Ptm]]$"t-tests"$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                                          Filter_up = sort(which(ptmpep[[g[i]]] %in% up)),
                                                                          Filter_down = sort(which(ptmpep[[g[i]]] %in% down)),
                                                                          Filter = sort(which(ptmpep[[g[i]]] %in% c(up, down))))
      }
    }
    if ("ref" %in% filter_types) {
      PTMs_Reg_filters[[Ptm]]$"t-tests"$"By reference" <- list()
      g2 <- as.data.frame(t(as.data.frame(strsplit(g1, "___"))))
      colnames(g2) <- VPAL$names
      tst <- apply(Exp.map[, RRG$names, drop = FALSE], 1L, paste, collapse = "___")
      tmp <- apply(g2[,RRG$names, drop = FALSE], 1L, paste, collapse = "___")
      g2$Ref <- sapply(tmp, \(x) { y <- unique(tst[which((Exp.map$Reference)&(tst == x))]) })
      for (i in unique(g2$Ref)) {
        w <- which(g2$Ref == i)
        PTMs_Reg_filters[[Ptm]]$"t-tests"$"By reference"[[i]] <- list(Columns = g[w],
                                                                      Filter_up = sort(which(apply(ptmpep[, g[w], drop = FALSE], 1L, \(x) {
                                                                        length(which(x %in% up))
                                                                      }) > 0L)),
                                                                      Filter_down = sort(which(apply(ptmpep[, g[w], drop = FALSE], 1L, \(x) {
                                                                        length(which(x %in% down))
                                                                      }) > 0L)),
                                                                      Filter = sort(which(apply(ptmpep[, g[w], drop = FALSE], 1L, \(x) {
                                                                        length(which(x %in% c(up, down)))
                                                                      }) > 0L)))
      }
    }
    if (sum(c("dat", "dat2") %in% filter_types)) {
      PTMs_Reg_filters[[Ptm]]$"t-tests"$"Whole dataset" <- list(Columns = g,
                                                                Filter_up = sort(which(apply(ptmpep[, g, drop = FALSE], 1L, \(x) {
                                                                  length(which(x %in% up))
                                                                }) > 0L)),
                                                                Filter_down = sort(which(apply(ptmpep[, g, drop = FALSE], 1L, \(x) {
                                                                  length(which(x %in% down))
                                                                }) > 0L)),
                                                                Filter = sort(which(apply(ptmpep[, g, drop = FALSE], 1L, \(x) {
                                                                  length(which(x %in% c(up, down)))
                                                                }) > 0L)))
    }
    if (F.test) {
      #kol <- grep(topattern(pepRf), colnames(ptmpep), value = TRUE)
      #View(ptmpep[, kol])
      ReportCalls <- AddMsg2Report(Msg = " -> ", ptm, " F-tests volcano plots", Space = FALSE)
      # NB: id.col below should be unique!
      stopifnot(length(unique(ptmpep$Name)) == nrow(ptmpep))
      #
      dataType <- "modPeptides"
      #
      FSrc %<o% paste0(libPath, "/extdata/Sources/run_F_test.R")
      #rstudioapi::documentOpen(FSrc)
      tstFtst <- try(source(FSrc, local = FALSE), silent = TRUE)
      #
      if (!inherits(tstFtst, "try-error")) {
        #F_test_ref_ratios %<o% F_volc$`Reference ratios` # Not needed
        volcano.plots[[Ptm]]$"F-tests_Unlabelled" <- F_volc$Plots$"Unlabelled"
        volcano.plots[[Ptm]]$"F-tests_Labelled" <- F_volc$Plots$"Labelled"
        n2 <- names(volcano.plots[[Ptm]]$"F-tests_Labelled")
        dir <- modDirs[3L]
        for (ttl in n2) {
          plot <- volcano.plots[[Ptm]]$"F-tests_Labelled"[[ttl]]
          ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
        }
        # Legacy code for web-hosted plotly plots:
        if ((create_plotly)&&(!create_plotly_local)) { plot_ly[[paste0(Ptm, "_Volcano plots (F-tests)")]] <- F_volc$"Plotly plots" }
        # Create F-test filters:
        g <- grep("^mod\\. F-test Regulated - ", colnames(PTMs_F_test_data[[Ptm]]), value = TRUE)
        g1 <- gsub("^mod\\. F-test Regulated - ", "", g)
        PTMs_Reg_filters[[Ptm]]$"F-tests" <- list()
        if ("con" %in% filter_types) {
          PTMs_Reg_filters[[Ptm]]$"F-tests"$"By condition" <- list()
          for (i in seq_along(g)) {
            PTMs_Reg_filters[[Ptm]]$"F-tests"$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                                              Filter_up = sort(which(PTMs_F_test_data[[Ptm]][[g[i]]] %in% up)),
                                                                              Filter_down = sort(which(PTMs_F_test_data[[Ptm]][[g[i]]] %in% down)),
                                                                              Filter = sort(which(PTMs_F_test_data[[Ptm]][[g[i]]] %in% c(up, down))))
          }
        }
        if ("ref" %in% filter_types) {
          PTMs_Reg_filters[[Ptm]]$"F-tests"$"By reference" <- list()
          g2 <- as.data.frame(t(as.data.frame(strsplit(g1, " - "))))
          row.names(g2) <- NULL
          colnames(g2) <- c("Group", "Analysis", "Condition")
          g2$Ref <- ""
          w <- grep(" VS ", g2$Condition)
          if (length(w)) {
            g2$tmp[w] <- sapply(strsplit(gsub("\\)$", "", g2$Condition[w]), split = " VS "), \(x) { x[[2L]] })
            g2$Ref[w] <- apply(g2[w, c("Group", "Analysis", "tmp")], 1L, \(x) {
              paste0(x[[1L]], " - ", x[[2L]], ", ref: ", x[[3L]])
            })
          }
          w <- which((grepl(" vs ", g2$Condition))&(!grepl(" VS ", g2$Condition)))
          if (length(w)) {
            g2$tmp[w] <- sapply(strsplit(g2$Condition[w], split = " vs "), \(x) { x[[2L]] })
            g2$Ref[w] <- apply(g2[w, c("Group", "Analysis", "tmp")], 1L, \(x) {
              paste0(x[[1L]], ", ref: ", x[[2L]])
            })
          }
          w <- which(!grepl(" vs ", g2$Condition))
          if (length(w)) { g2$Ref[w] <- g1[w] }
          for (i in unique(g2$Ref)) {
            w <- which(g2$Ref == i)
            PTMs_Reg_filters[[Ptm]]$"F-tests"$"By reference"[[i]] <- list(Columns = g[w],
                                                                          Filter_up = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1L, \(x) {
                                                                            length(which(x %in% up))
                                                                          }) > 0L)),
                                                                          Filter_down = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1L, \(x) {
                                                                            length(which(x %in% down))
                                                                          }) > 0L)),
                                                                          Filter = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1L, \(x) {
                                                                            length(which(x %in% c(up, down)))
                                                                          }) > 0L)))
          }
        }
        if ("dat" %in% filter_types) {
          PTMs_Reg_filters[[Ptm]]$"F-tests"$"By analysis" <- list()
          g2 <- as.data.frame(t(as.data.frame(strsplit(g1, " - "))))
          colnames(g2) <- c("Group", "Analysis", "Condition")
          g2$Group_Analysis <- apply(g2[,c("Group", "Analysis")], 1L, paste, collapse = "_")
          for (i in unique(g2$Group_Analysis)) {
            w <- which(g2$Group_Analysis == i)
            PTMs_Reg_filters[[Ptm]]$"F-tests"$"By analysis"[[i]] <- list(Columns = g[w],
                                                                         Filter_up = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1L, \(x) {
                                                                           length(which(x %in% up))
                                                                         }) > 0L)),
                                                                         Filter_down = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1L, \(x) {
                                                                           length(which(x %in% down))
                                                                         }) > 0L)),
                                                                         Filter = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1L, \(x) {
                                                                           length(which(x %in% c(up, down)))
                                                                         }) > 0L)))
          }
        }
        if ("dat2" %in% filter_types) {
          PTMs_Reg_filters[[Ptm]]$"F-tests"$"Whole dataset" <- list(Columns = g,
                                                                    Filter_up = sort(which(apply(PTMs_F_test_data[[Ptm]][, g, drop = FALSE], 1L, \(x) {
                                                                      length(which(x %in% up))
                                                                    }) > 0L)),
                                                                    Filter_down = sort(which(apply(PTMs_F_test_data[[Ptm]][, g, drop = FALSE], 1L, \(x) {
                                                                      length(which(x %in% down))
                                                                    }) > 0L)),
                                                                    Filter = sort(which(apply(PTMs_F_test_data[[Ptm]][, g, drop = FALSE], 1L, \(x) {
                                                                      length(which(x %in% c(up, down)))
                                                                    }) > 0L)))
        }
      } else { warning("The F-test did not generate any plots! Investigate!") }
    }
    #
    # Gene-Set Enrichment Analysis (GSEA)
    if (runGSEA) {
      dataType <- "modPeptides"
      GSEAmode <- "standard"
      Src <- paste0(libPath, "/extdata/Sources/GSEA.R")
      #rstudioapi::documentOpen(Src)
      source(Src, local = FALSE)
    }
    #
    # Heatmap
    g <- paste0(pepRf, RSA$values)
    grps <- Exp.map[match(RSA$values, Exp.map$Ref.Sample.Aggregate), VPAL$column]
    w <- which(g %in% colnames(ptmpep))
    g <- g[w]; grps <- grps[w]
    temp <- ptmpep[, g]
    rownames(temp) <- gsub("\n", " ", ptmpep$Name)
    kol <- gsub(topattern(pepRf), "", g)
    colnames(temp) <- gsub(topattern(pepRf), "", colnames(temp))
    kol <- cleanNms(kol)
    colnames(temp) <- cleanNms(colnames(temp))
    w <- which(!is.finite(as.matrix(temp)), arr.ind = TRUE)
    temp[w] <- NA
    rwMns <- rowMeans(temp[, kol], na.rm = TRUE)
    w <- which(!is.na(rwMns))
    temp <- temp[w,]
    rwMns <- rwMns[w]
    temp[, kol] <- sweep(temp[, kol], 1L, rwMns, "-")
    temp2 <- as.matrix(Data_Impute2(temp, grps)$Imputed_data)
    hcluster <- hclust(dist(t(temp2)))
    hdendro <- as.dendrogram(hcluster)
    hord <- order.dendrogram(hdendro)
    hord <- colnames(temp)[hord]
    vcluster <- hclust(dist(temp2))
    vdendro <- as.dendrogram(vcluster)
    vord <- order.dendrogram(vdendro)
    vord <- rownames(temp)[vord]
    require(ggdendro)
    hdendro.plot <- ggdendrogram(data = hdendro) + theme(axis.text.y = element_text(size = 0.1),
                                                         plot.margin = margin(0L, 0L, 0L, 0L, "cm"))
    vdendro.plot <- ggdendrogram(data = vdendro, rotate = TRUE, labels = FALSE, leaf_labels = FALSE) +
      theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
            panel.background = element_rect(fill = "transparent", colour = NA), 
            plot.background = element_rect(fill = "transparent", colour = NA),
            plot.margin = margin(0L, 0L, 0L, 0L, "cm"))
    # Data wrangling
    temp2 <- set_colnames(reshape::melt.data.frame(temp), c("Sample", "value"))
    temp2$Label <- rownames(temp)
    temp2$Sample <- as.character(temp2$Sample)
    temp2$Colour <- "grey"
    temp2$Size <- 1L
    # Extract the order of the tips in the dendrograms
    # Order the levels according to their position in the clusters
    temp2$Xmin <- match(temp2$Sample, hord)-1L
    temp2$Ymin <- match(temp2$Label, vord)-1L
    Xscale <- length(unique(temp2$Sample))
    temp2$Label2 <- temp2$Label
    w <- which(nchar(temp2$Label2) > 25L)
    temp2$Label2[w] <- paste0(substr(temp2$Label2[w], 1L, 22L), "...")
    # Create heatmap plot
    w1 <- which(temp2$Colour == "green")
    w2 <- which((temp2$Xmin == max(temp2$Xmin))&(temp2$Colour == "green"))
    nm <- paste0("Heatmap\n", Ptm, "-modified peptides")
    heatmap.plot <- ggplot(temp2) +
      geom_rect(aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
      geom_text(data = temp2, aes(y = Ymin+0.5, label = Label2),
                x = length(unique(temp2$Sample)), colour = "grey", hjust = 0, vjust = 0.5, size = 2L) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.margin = margin(0L, 0L, 0L, 0L, "cm")) +
      scale_fill_gradient2(low = "darkblue", mid = "lightgrey", high = "darkred") +
      xlab(NULL) + ylab(NULL)
    leg <- get_legend(heatmap.plot)
    heatmap.plot <- heatmap.plot + theme(legend.position = "none")
    htmp <- arrangeGrob(grobs = list(hdendro.plot, leg, heatmap.plot, vdendro.plot),
                        widths = c(1L, 2L*Xscale, 1L, 3L), heights = c(3L, 5L),
                        layout_matrix = rbind(c(NA, 1L, NA, 2L), c(3L, 3L, 3L, 4L)),
                        padding = 5L)
    nm <- gsub("\n", " - ", nm)
    suppressMessages({
      ggsave(paste0(modDirs[1L], "/", nm, ".jpeg"), htmp, width = 20L, height = 20L, units = "in", dpi = 600L)
      ggsave(paste0(modDirs[1L], "/", nm, ".pdf"), htmp, width = 20L, height = 20L, units = "in", dpi = 600L)
    })
    ReportCalls <- AddPlot2Report(Title = nm, Dir = modDirs[1L])
    #system(paste0("open \"", modDirs[1L], "/", nm, ".jpeg", "\""))
    #system(paste0("open \"", modDirs[1L], "/", nm, ".pdf", "\""))
    #
    # Gene Ontology terms enrichment analysis
    if (enrichGO) {
      if ((!exists("GO_mappings"))&&(file.exists("GO_mappings.RData"))) {
        loadFun("GO_mappings.RData")
      }
      if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) {
        loadFun("GO_terms.RData")
      }
      p <- strsplit(ptmpep$Proteins, ";")
      test <- sapply(annot.col, \(x) { x %in% colnames(db) })
      annot.col2 <- annot.col[which(test)]
      tmpDB <- db[which(db$Observed), c("Protein ID", annot.col2)]
      clusterExport(parClust, c("tmpDB", "annot.col2"), envir = environment())
      temp <- parLapply(parClust, p, \(x) {
        m <- match(x, tmpDB$"Protein ID")
        y <- tmpDB[m, annot.col2]
        if (length(m) > 1L) {
          y <- apply(y, 2L, \(z) {
            z <- z[which(!is.na(z))]
            paste(sort(unique(unlist(strsplit(as.character(z), ";")))), collapse = ";")
          })
        }
        y <- setNames(as.character(unlist(y)), annot.col2)
        return(y)
      })
      temp <- as.data.frame((do.call(rbind, temp)))
      for (i in annot.col2) {
        w <- which(is.na(temp[[i]]))
        if (length(w)) { temp[w, i] <- "" }
      }
      ptmpep[, annot.col2] <- temp
      Kol <- unique(c("Modified sequence", "Evidence IDs", "Sequence",
                      myIDcol, "Proteins", "Leading razor proteins",
                      "PEP", "id", "Protein group IDs", "Razor protein group ID", 
                      "Leading razor protein", "Gene names", "Protein names", "Gene names (all)",
                      "Protein names (all)", "Weights", "Match(es)", "Code", "Name", annot.col2,
                      "Potential contaminant"))
      Kol <- Kol[which(Kol %in% colnames(ptmpep))]
      PTMs_GO_enrich.dat[[Ptm]] <- list()
      PTMs_GO_enrich.FCRt[[Ptm]] <- list()
      PTMs_GO_enrich.tbl[[Ptm]] <- list()
      PTMs_GO_Plots[[Ptm]] <- list()
      PTMs_Reg_GO_terms[[Ptm]] <- list()
      Tsts <- c("t-tests", "F-tests")
      WhTsts <- which(Tsts %in% names(PTMs_Reg_filters[[Ptm]]))
      for (tt in WhTsts) { #tt <- 1 #tt <- 2
        tstrt <- Tsts[tt]
        stopifnot(!is.na(tstrt))
        dir <- paste0(modDirs[1L], "/GO enrich/", tstrt)
        modDirs <- c(modDirs, dir)
        if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
        dirlist <- union(dirlist, dir)
        filt <- PTMs_Reg_filters[[Ptm]][[tstrt]]
        #By <- c("By condition", "By reference", "By analysis", "Whole dataset")
        By <- "By condition"
        By <- By[which(By %in% names(filt))]
        if (length(By)) {
          for (bee in By) { #bee <- By[1L]
            flt <- filt[[bee]]
            if (bee == "Whole dataset") { flt <- list("Whole dataset" = flt) }
            tstbee <- paste0(tstrt, "_", tolower(bee))
            if (length(flt)) {
              if (tt == 1L) { tmpdat <- ptmpep }
              if (tt == 2L) { tmpdat <- PTMs_F_test_data[[Ptm]] }
              flt <- if (tt %in% 1L:2L) {
                flt[myContrasts$Contrast[which(myContrasts$Secondary == "")]]
              } else {
                flt[order(names(flt))]
              }
              for (kol in Kol) { tmpdat[, kol] <- ptmpep[, kol] }
              flt <- flt[order(names(flt))]
              reg <- setNames(lapply(flt, \(x) { list(x$Columns) }), names(flt))
              reg <- set_colnames(reshape2::melt(reg), c("Name", "Bleh", "For"))
              reg$Bleh <- NULL
              # PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]] <- paste0("Mean ", c(ptms.ratios.ref[length(ptms.ratios.ref)],
              #                                                           "log2(Ratio) - ")[tt])
              PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]] <- ptms.ratios.ref
              reg$ParentFC <- gsub(".*Regulated - ", PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]], reg$Name)
              reg$FCname <- paste0(PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]], reg$For)
              tmpdat[, Kol] <- ptmpep[, Kol]
              #tmpdat <- get(c("ptmpep", "PTMs_F_test_data[[Ptm]]", "ptmpep", "PTMs_allSAINTs")[tt]) # PTMs_allSAINTs doesn't exist
              UF <- unique(reg$For)
              temPTM <- as.data.frame(do.call(cbind, lapply(UF, \(x) { #x <- UF[1L]
                x <- reg$ParentFC[which(reg$For == x)]
                x <- if (length(x) > 1L) { apply(tmpdat[, x], 1L, log_ratio_av) } else { tmpdat[[x]] }
                return(x)
              })))
              colnames(temPTM) <- paste0(PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]], UF)
              temPTM[, Kol] <- ptmpep[, Kol]
              #temPTM$"First protein" <- sapply(strsplit(temPTM[[myIDcol]], ";"), \(x) { unlist(x)[1L] })
              if (scrptTypeFull == "withReps_PG_and_PTMs") {
                temp <- listMelt(strsplit(temPTM$`Protein group IDs`, ";"), temPTM$id)
                temp$Genes <- PG$Genes[match(temp$value, PG$id)]
                temp <- aggregate(temp$Genes, list(temp$L1), \(x) {
                  paste(sort(unique(unlist(strsplit(x, ";")))), collapse = ";")
                })
              }
              if (scrptTypeFull == "withReps_PTMs_only") {
                temp <- listMelt(strsplit(temPTM$Proteins, ";"), temPTM$id)
                temp$Genes <- db$Gene[match(temp$value, db$`Protein ID`)]
                temp <- aggregate(temp$Genes, list(temp$L1), \(x) {
                  paste(sort(unique(unlist(strsplit(x, ";")))), collapse = ";")
                })
              }
              temPTM$Genes <- temp$x[match(temPTM$id, temp$Group.1)]
              PTMs_GO_enrich.dat[[Ptm]][[tstbee]] <- temPTM
              PTMs_GO_enrich.tbl[[Ptm]][[tstbee]] <- reg
              flt <- setNames(lapply(UF, \(x) { flt[[x]]$Filter }), UF)
              ttr <- btr <- ""
              if (length(By) > 1L) { ttr <- btr <- paste0(tolower(bee), "_") }
              Pep.Ref.Filt <- tmpFilt <- setNames(lapply(names(flt), \(x) { 1L:nrow(PTMs_GO_enrich.dat[[Ptm]][[tstbee]]) }), names(flt))
              if (GO.enrich.MultiRefs) {
                Pep.Ref.Filt <- try(setNames(lapply(names(flt), \(x) { #x <- names(flt)[1L]
                  m <- match(x, myContrasts$Contrast)
                  A_ <- myContrasts$A_samples[[m]]
                  B_ <- myContrasts$B_samples[[m]]
                  y <- Exp.map[which(Exp.map$Ref.Sample.Aggregate %in% c(A_, B_)), GO.enrichment.Ref.Aggr$column]
                  z <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[GO.enrichment.Ref.Aggr$column]] %in% y)]
                  w1 <- which(apply(ptmpep[, paste0(pepRf, z)], 1L, \(x) {
                    sum(is.all.good(x, 2L))
                  }) > 0L)
                  w2 <- flt[[x]] # Required for if we are imputing missing values
                  return(sort(unique(c(w1, w2))))
                }), names(flt)), silent = TRUE)
                if (inherits(Pep.Ref.Filt, "try-error")) {
                  warning("Invalid \"GO.enrichment.Ref.Aggr\" argument: multiple references for GO enrichment are only feasible if each enrichment filter maps to a single reference! Skipping...")
                  GO.enrich.MultiRefs <- FALSE
                  Pep.Ref.Filt <- tmpFilt
                }
              }
              if ((length(Pep.Ref.Filt) > 1L)||(!is.na(Pep.Ref.Filt))) {
                flt <- setNames(lapply(names(flt), \(x) { flt[[x]][which(flt[[x]] %in% Pep.Ref.Filt[[x]])] }), names(flt))
              }
              # Also save the reference filters 
              nms <- names(PTMs_Reg_filters[[Ptm]][[tstrt]][[bee]])
              for (nm in nms) {
                PTMs_Reg_filters[[Ptm]][[tstrt]][[bee]][[nm]]$Background_filter <- Pep.Ref.Filt[[nm]]
              }
              #
              ReportCalls <- AddMsg2Report(Msg = " - GO terms enrichment analysis", Space = FALSE)
              #
              Mode <- "regulated"
              dataType <- "modPeptides"
              if ((!exists("GO_mappings"))&&(file.exists("GO_mappings.RData"))) {
                loadFun("GO_mappings.RData")
              }
              if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) {
                loadFun("GO_terms.RData")
              }
              #
              Src <- paste0(libPath, "/extdata/Sources/GO_enrich.R")
              #rstudioapi::documentOpen(Src)
              source(Src, local = FALSE)
              #
              clueGO_outDir <- dir
              clueGO_type <- "Enrichment (Right-sided hypergeometric test)"
              Src <- paste0(libPath, "/extdata/Sources/ClueGO_enrich.R")
              #rstudioapi::documentOpen(Src)
              source(Src, local = FALSE)
              #
              # Cleanup - do it now, not within sources!
              try(rm(list = allArgs), silent = TRUE)
              #
              PTMs_GO_Plots[[Ptm]][[tstbee]] <- goRES
              #
              if (!is.null(names(PTMs_GO_Plots[[Ptm]][[tstbee]]))) {
                if ("All_GO_terms" %in% names(PTMs_GO_Plots[[Ptm]][[tstbee]])) {
                  GO_terms <- PTMs_GO_Plots[[Ptm]][[tstbee]]$All_GO_terms
                  PTMs_GO_Plots[[Ptm]][[tstbee]]$All_GO_terms <- NULL  
                }
                n2 <- names(PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_plots)
                for (ttl in n2) { #ttl <- n2[1L]
                  plot <- PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_plots[[ttl]]
                  ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
                }
                if ((create_plotly)&&(!create_plotly_local)) {
                  plot_ly[[paste0(Ptm, "_GO plots - Regulated vs Observed - ", tstbee)]] <- PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_plot_ly
                }
                temp <- PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms
                temp$Mapping <- NULL
                if ("Offspring" %in% colnames(temp)) {
                  temp$Offspring <- parSapply(parClust, temp$Offspring, paste, collapse = ";")
                }
                temp$`Protein table row(s)` <- NULL
                colnames(temp) <- cleanNms(colnames(temp), start = FALSE)
                gn <- grep("^Genes", colnames(temp), value = TRUE)
                pr <- grep("^Proteins", colnames(temp), value = TRUE)
                pp <- grep("^Peptide IDs", colnames(temp), value = TRUE)
                kn <- grep("^Count", colnames(temp), value = TRUE)
                pv <- grep("^Pvalue", colnames(temp), value = TRUE)
                zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
                lf <- grep("^logFC", colnames(temp), value = TRUE)
                si <- grep("^Significance", colnames(temp), value = TRUE)
                #lp <- grep("^Leading protein IDs", colnames(temp), value = TRUE)
                #kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pp, lp))]
                kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pp, pr))]
                #temp <- temp[, c(kl, si, gn, pp, lp, kn, pv, zs, lf)]
                temp <- temp[, c(kl, si, gn, pp, pr, kn, pv, zs, lf)]
                w <- apply(temp[, pv, drop = FALSE], 1L, \(x) { sum(!is.na(x)) }) > 0L
                temp <- temp[w,]
                tst <- apply(temp[, kn, drop = FALSE], 1L, \(x) { sum(x[which(!is.na(x))]) })
                temp <- temp[order(tst, decreasing = TRUE),]
                temp <- temp[order(temp$Ontology, decreasing = FALSE),]
                PTMs_Reg_GO_terms[[Ptm]][[tstbee]] <- temp
                write.csv(temp, file = paste0(dir, "/", Ptm, " GO terms - ", tstbee, ".csv"), row.names = FALSE)
                w <- which(vapply(colnames(temp), \(x) { is.character(temp[[x]]) }, TRUE))
                if (length(w)) {
                  for (i in w) { #i <- w[1L]
                    w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
                    if (length(w1)) {
                      temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1L, ExcelMax-3L), "...")
                    }
                  }
                }
                require(openxlsx)
                HdrStl <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                                      wrapText = TRUE, numFmt = "TEXT", fontSize = 12L)
                wb <- createWorkbook()
                kount <- 0L
                for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1L]
                  w <- which(temp$Ontology == Ontologies[ont])
                  if (length(w)) {
                    kount <- kount + 1L
                    addWorksheet(wb, ont)
                    writeData(wb, ont, temp[w,])
                    setRowHeights(wb, ont, 1L, 60L)
                    setColWidths(wb, ont, 1L:ncol(temp), 12L)
                    setColWidths(wb, ont, which(colnames(temp) == "Term"), 45L)
                    setColWidths(wb, ont, which(colnames(temp) %in% gn), 20L)
                    setColWidths(wb, ont, which(colnames(temp) %in% pr), 20L)
                    setColWidths(wb, ont, which(colnames(temp) %in% pp), 20L)
                    addStyle(wb, ont, HdrStl, 1L, 1L:ncol(temp))
                    addStyle(wb, ont, createStyle(numFmt = "0"), 2L:(length(w)+1L),
                             which(colnames(temp) %in% kn), gridExpand = TRUE)
                    addStyle(wb, ont, createStyle(numFmt = "0.000"), 2L:(length(w)+1L),
                             which(colnames(temp) %in% c(zs, lf, pv)), gridExpand = TRUE)
                  }
                }
                if (kount) {
                  saveWorkbook(wb, paste0(dir, "/", Ptm, " GO terms - ", tstbee, ".xlsx"), overwrite = TRUE)
                  if (tt == 1L) {
                    Kol2 <- paste0("Significance - ", cleanNms(VPAL$values), " ", max(BH.FDR)*100, "%")
                    Kol2 <- Kol2[which(Kol2 %in% colnames(PTMs_Reg_GO_terms[[Ptm]][[tstbee]]))]
                  }
                  if (tt == 2L) {
                    Kol2 <- paste0("Significance - ", names(PTMs_Reg_filters[[Ptm]][[tstrt]][[bee]]), " ", max(BH.FDR)*100, "%")
                    Kol2 <- Kol2[which(Kol2 %in% colnames(PTMs_Reg_GO_terms[[Ptm]][[tstbee]]))]
                  }
                  if (length(Kol2)) {
                    w <- if (length(Kol2) > 1L) {
                      which(apply(PTMs_Reg_GO_terms[[Ptm]][[tstbee]][, Kol2], 1L, \(x) {"+" %in% x}))
                    } else { which(sapply(PTMs_Reg_GO_terms[[Ptm]][[tstbee]][,Kol2], \(x) {"+" %in% x})) }
                    write.csv(PTMs_Reg_GO_terms[[Ptm]][[tstbee]][w,],
                              file = paste0(dir, "/", Ptm, " regulated GO terms - ", tstbee, ".csv"),
                              row.names = FALSE)
                  }
                  # Summary table and heatmap of number of co-regulated GO terms
                  Kol3 <- Kol2
                  N <- length(Kol3)
                  if (N > 1L) {
                    myRng <- 2L:(N+1L)
                    temp <- as.data.frame(matrix(rep("", (N+1L)^2L), ncol = N+1L))
                    W <- lapply(Kol3, \(x) { which(PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms[[x]] == "+") })
                    names(W) <- gsub(paste0("^Significance - | ", max(BH.FDR)*100, "%$"), "", Kol3)
                    temp[myRng, 1L] <- temp[1L, myRng] <- names(W)
                    for (i in myRng) { #i <- 2
                      temp[i, myRng] <- sapply(Kol3, \(x) {
                        sum((PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms[[x]] == "+")&(PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms[[Kol3[i]]] == "+"),
                            na.rm = TRUE)
                      })
                    }
                    names(W) <- cleanNms(gsub(" [0-9]+%$", "", names(W)))
                    tst <- sapply(strsplit(names(W), " - "), length)
                    tst <- (min(tst) > 1L)&(length(unique(tst)) == 1L)
                    if (tst) {
                      tst <- as.data.frame(t(sapply(strsplit(names(W), " - "), unlist)))
                      l <- apply(tst, 2L, \(x) { length(unique(x)) })
                      tst <- tst[, which(l > 1L), drop = FALSE]
                      names(W) <- apply(tst, 1L, paste, collapse = " - ")
                    }
                    temp[myRng, 1L] <- temp[1L, myRng] <- names(W)
                    nm <- paste0("N. of co-regulated GO terms\n", tstrt, "\n(", tolower(bee), ")")
                    write.csv(temp, file = paste0(dir, "/", gsub("\n", " - ", gsub("\n\\(", " (", nm)), ".csv"), row.names = FALSE)
                    temp2 <- temp[myRng, myRng]
                    colnames(temp2) <- temp[1L, myRng]
                    rownames(temp2) <-  temp[myRng, 1L]
                    for (i in 1L:nrow(temp2)) { temp2[[i]] <- as.numeric(temp2[[i]]) }
                    if (max(is.all.good(unlist(temp2)))) {
                      temp2 <- as.matrix(temp2)
                      basic.heatmap(temp2, "N. of co-regulated GO terms", paste0(tstrt, "\n(", tolower(bee), ")"),
                                    save = c("pdf", "jpeg"), folder = dir)
                    } else { warning(paste0(tstrt, " co-regulated GO terms heatmap: nothing to plot")) }
                  } else { warning(paste0(tstrt, " performed for only one condition, skipping.")) }
                }
              }
            } else { warning(paste0("Filter ", tstbee, " has length 0, skipping.")) }
          }
        } else { warning(paste0("No filters available for ", tstrt, ", skipping.")) }
      }
    }
    PTMs_pep[[Ptm]] <- ptmpep
    ReportCalls <- AddSpace2Report()
    PTMs_int.ref[[Ptm]] <- ptms.ref
    #PTMs_rat.ref[[Ptm]] <- ptms.ratios.ref
    #write.csv(ptmpep, paste0("Tables/", Ptm, "-modified peptides.csv"), row.names = FALSE)
    #w <- grsep2(prot.list, ptmpep$Proteins)
    #if (length(w)) {
    #  write.csv(ptmpep[w,], paste0("Tables/", Ptm, "-modified peptide - proteins in list.csv"), row.names = FALSE)
    #}
  }
} else {
  if (scrptTypeFull == "withReps_PTMs_only") {
    stop("Really? There is no PTM-modified class of peptides to analyze? Why did you run this workflow then?")
  }
}

# Once we are done with using GO stuff, backup the final versions
if (Annotate) {
  # (Testing for existence in case we are rerunning part of the script)
  if (exists("GO_terms")) {
    saveFun(GO_terms, file = "GO_terms.RData")
    rm(GO_terms)
  }
  if (exists("GO_mappings")) {
    saveFun(GO_mappings, file = "GO_mappings.RData")
    rm(GO_mappings)
  }
  # otherwise they are also unused further down.
  .obj <- .obj[which(!.obj %in% c("GO_terms", "GO_mappings"))]
}
if (CytoScape) {
  # Then close Cytoscape:
  tst <- try({
    invisible({
      RCy3::closeSession(save.before.closing = FALSE)
      cmd <- paste0("taskkill/im \"", gsub(".+/", "", CytoScVrs), "\" /f")
      #cat(cmd)
      shell(cmd)
    })
  })
  if (inherits(tst, "try-error")) {
    cat("(This error can be ignored as it does not interrupt script execution.\nIt looks like you are trying to close Cytoscape although it is already open (maybe running this script by bits after an interruption?)\n\n")
  }
}
