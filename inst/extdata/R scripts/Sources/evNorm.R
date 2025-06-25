# Optional - Normalize evidence MS1 intensities, then, if applicable, MS2 reporter (Isobaric labelling) or fragment (DIA) intensities
source(parSrc, local = FALSE)
if (Param$Norma.Ev.Intens) {
  msg <- paste0(evNm, "s-level normalisations:\n------------------------------------\n")
  ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
  # Define groups - this will ensure that, if phospho (or other) -enrichment took place, these peptides will be normalized separately
  ev$"Normalisation group" <- "Standard"
  if ("PTM-enriched" %in% colnames(Frac.map)) {
    # (In column "PTM enriched", use NA to indicate no enrichment!!!)
    if (!"PTM-enriched" %in% colnames(Frac.map)) { Frac.map$"PTM-enriched" <- NA }
    ptmChck <- unique(Frac.map$"PTM-enriched")
    ptmChck <- ptmChck[which(!is.na(ptmChck))]
    if (sum(!ptmChck %in% Modifs$`Full name`)) { stop("Some of the modifications in column \"PTM-enriched\" of Fractions map are invalid!") }
    # Here is what we want to do for those modifications:
    # - For enriched samples, keep only peptides with the target modification
    # - For other samples, remove all peptides with the modification which were also found in enriched samples
    if (length(ptmChck)) {
      for (ptm in ptmChck) { #ptm <- ptmChck[1]
        # Below "modified" means "modified with ptm" and "enriched" means "enriched for ptm"
        mrk <- Modifs$Mark[match(ptm, Modifs$`Full name`)]
        rw1 <- Frac.map$"Raw file"[which(Frac.map$"PTM-enriched" == ptm)]
        rw0 <- Frac.map$"Raw file"[which((is.na(Frac.map$"PTM-enriched"))|(Frac.map$"PTM-enriched" != ptm))]
        w1 <- ev$id[which(ev$"Raw file path" %in% rw1)] # PSMs from enriched runs
        ev$"Normalisation group"[match(w1, ev$id)] <- ptm
        w0 <- ev$id[which(ev$"Raw file path" %in% rw0)] # PSMs from non-enriched runs
        w2 <- ev$id[which(!ev$"Raw file path" %in% c(rw0, rw1))] # Any others
        w1a <- w1[grep(mrk, ev$"Modified sequence"[match(w1, ev$id)])] # Modified PSMs from enriched samples (i.e. what we were trying to enrich!)
        if (length(w1a)) {
          w0a <- w0[which(!ev$"Modified sequence"[match(w0, ev$id)] %in% unique(ev$"Modified sequence"[match(w1a, ev$id)]))] # Un-modified PSMs from non-enriched runs
          l1 <- length(w1)-length(w1a) # This is the number of un-modified PSMs we are removing from enriched runs
          l0 <- length(w0)-length(w0a) # This is the number of modified PSMs we are removing from non-enriched runs
          if (l1) {
            msg <- paste0("Removing ", l1, " peptide evidences without the ", ptm, " modification from ", ptm,
                          "-enriched raw files (", round(100*l1/length(w1), 2), "%)!")
            ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
          }
          if (l0) {
            msg <- paste0("Removing ", l0, " ", ptm, "-modified peptide evidences from un-enriched samples!")
            ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
          }
          ev <- ev[which(ev$id %in% c(w0a, w1a, w2)),]
        } else {
          msg <- paste0("Not a single ", ptm, "-modified evidence found in ", ptm, "-enriched raw files, investigate!")
          ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
        }
      }
    }
  }
  #
  #
  #
  # Now, after much thought, I have decided to exclude Norm.Groups from PSMs-level pre-normalisations.
  # Why?
  #  - Norm.Groups can better be used later at peptides and protein level.
  #  - At evidences level, there is sometimes no unique relationship between experimental factors and rows of the table (TMT).
  #  - When dealing with fractionated data, this may result in too many groups and no cross-fractions normalisation.
  #    In these cases, I feel that normalisation of fractions is more important for unbiased quantitation.
  #
  # Below the code is commented... but be AWARE that uncommenting it is insufficient.
  # You would also need to check the subsequent code, and test on both an LFQ and a fractionated TMT (CaBeAHlavata2) dataset.
  #
  # Start of commented chunk
  #
  # nms <- Norm.Groups$names
  # tmp <- listMelt(Exp.map$MQ.Exp, 1:nrow(Exp.map), c("MQ.Exp", "row"))
  # tmp[, nms] <- Exp.map[tmp$row, nms]
  # tmp <- aggregate(tmp[, nms], list(tmp$MQ.Exp), unique)
  # colnames(tmp) <- c("MQ.Exp", nms)
  # for (nm in nms) { #nm <- nms[1]
  #   w2 <- which(vapply(tmp[[nm]], length, 1) > 1)
  #   tmp[w2, nm] <- "Mixed_values!"
  #   tmp[[nm]] <- sapply(tmp[[nm]], unlist)
  # }
  # ev[, nms] <- tmp[match(ev$MQ.Exp, tmp$MQ.Exp), nms]
  # ev$"Normalisation group" <- do.call(paste, c(ev[, c(nms, "Normalisation group")], sep = "_"))
  #
  # End of commented chunk
  #
  #
  # Check that no PSM is assigned NA as normalisation group
  #aggregate(ev$"Normalisation group", list(ev$MQ.Exp), function(x) { sum(is.na(x)) == 0 })
  # Check normalisation groups
  #aggregate(ev$MQ.Exp, list(ev$"Normalisation group"), function(x) { length(unique(x)) }) # Just one...
  # ... except when we have fractions:
  # Indeed, whilst those groups are fine, we also want to normalize per fraction at this stage,
  # so that each series of equivalent fractions from different samples get normalized to each other.
  #
  mrmgrps <- unique(ev$"Normalisation group")
  fr <- unique(ev$Fraction)
  tmp <- data.frame(Grp = as.character(sapply(mrmgrps, function(x) { rep(x, length(fr)) })),
                    Frac = as.character(rep(fr, length(mrmgrps))))
  tmp$Nm <- do.call(paste, c(tmp, sep = "_"))
  ev$"Normalisation group + Fraction" <- NA
  for (i in 1:nrow(tmp)) {
    w <- which((ev$"Normalisation group" == tmp$Grp[i])&(ev$Fraction == tmp$Frac[i]))
    ev$"Normalisation group + Fraction"[w] <- tmp$Nm[i]
  }
  #aggregate(ev$MQ.Exp, list(ev$"Normalisation group + Fraction"), function(x) { length(unique(x)) })
  Norma.Ev.Intens.Groups %<o% set_colnames(aggregate(ev$"Normalisation group + Fraction", list(ev$"Raw file path"), unique),
                                           c("Raw file", "Groups"))
  tst <- aggregate(Norma.Ev.Intens.Groups$"Raw file", list(Norma.Ev.Intens.Groups$Groups), length)
  w <- which(tst$x > 1)
  #tst <- aggregate(Norma.Ev.Intens.Groups$"Raw file", list(Norma.Ev.Intens.Groups$Groups), list)
  if (length(w)) {
    msg <- " - Normalizing MS1-level PSM intensities"
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
    if (length(w) > 1) { cat(" (per fraction/PTM enrichment group)\n") }
    msg <- "   - Classic normalisation to the median"
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
    Norma.Ev.Intens.Groups <- Norma.Ev.Intens.Groups[which(Norma.Ev.Intens.Groups$Groups %in% tst$Group.1[w]),]
    # (Per fractions X PTM enrichment group)
    # Step 1a:
    ev.col["Normalisation"] <- paste0("norm. ", ev.col["Original"])
    ev[[ev.col["Normalisation"]]] <- NA
    Norm.Ev %<o% data.frame(Group = unique(Norma.Ev.Intens.Groups$Groups))
    Grps2 <- MQ.Exp
    Grps2Kol <- "MQ.Exp"
    if (LabelType == "Isobaric") {
      Grps2 <- Iso
      Grps2Kol <- "Isobaric.set"
    }
    for (grp2 in Grps2) { Norm.Ev[[paste0("Grp", grp2)]] <- 1 }
    for (grp in Norm.Ev$Group) { #grp <- Norm.Ev$Group[1]
      r <- Norma.Ev.Intens.Groups$"Raw file"[which(Norma.Ev.Intens.Groups$Groups == grp)]
      wg <- which(ev$"Raw file path" %in% r)
      M <- 10^median(is.all.good(log10(unlist(ev[wg, ev.col["Original"]])))) # For preserving original scale
      #M <- 10^mlv(is.all.good(log10(unlist(w[wg, ev.col["Original"]]))), method = "Parzen")[1]
      for (grp2 in Grps2) { #grp2 <- Grps2[1]
        w2 <- which(ev[wg, Grps2Kol] == grp2)
        if (length(w2)) {
          m <- 10^median(is.all.good(log10(ev[wg[w2], ev.col["Original"]])))
          #m <- 10^mlv(is.all.good(log10(ev[wg[w2], ev.col["Original"]])), method = "Parzen")[1]
          ev[wg[w2], ev.col["Normalisation"]] <- ev[wg[w2], ev.col["Original"]]*M/m
          Norm.Ev[match(grp, Norm.Ev$Group), paste0("Grp", grp2)] <- m/M
        }
      }
    }
    if (("Adv.Norma.Ev.Intens" %in% colnames(Param))&&(Param$Adv.Norma.Ev.Intens != FALSE)) {
      msg <- c("   - Levenberg-Marquardt normalisation..."#,
               #"     (this can take some time)",
               #"     ..."
      )
      ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Print = FALSE)
      cat(paste0(msg, "\n", collapse = "\n"))
      txtAdv <- "Evidence MS1 intensities"
      ev.col["Advanced normalisation"] <- paste0("AdvNorm. ", ev.col["Original"])
      ev[[ev.col["Advanced normalisation"]]] <- NA
      AdvNorm.Ev %<o% data.frame(Group = unique(Norma.Ev.Intens.Groups$Groups))
      for (grp2 in Grps2) { AdvNorm.Ev[[paste0("Grp", grp2)]] <- 1 }
      for (grp in Norm.Ev$Group) { #grp <- Norm.Ev$Group[1]
        w <- which(Norma.Ev.Intens.Groups$Groups == grp)
        r <- Norma.Ev.Intens.Groups$"Raw file"[which(Norma.Ev.Intens.Groups$Groups == grp)]
        wg <- which(ev$"Raw file path" %in% r)
        tmp <- data.table(Uniq = ev$`Unique State`[wg], Grp = ev[wg, Grps2Kol], Int = ev[wg, ev.col["Normalisation"]])
        tmp <- tmp[, list(x = sum(Int, na.rm = TRUE)), keyby = list(Group.1 = Uniq, Group.2 = Grp)]
        tmp <- as.data.frame(tmp)
        tmp <- suppressMessages(reshape2::dcast(tmp, Group.1~Group.2))
        colnames(tmp) <- c("Unique State", paste0("Grp", colnames(tmp)[2:ncol(tmp)]))
        kol <- paste0("Grp", Grps2)
        w <- which(kol %in% colnames(tmp))
        kol <- kol[w] 
        tmp2 <- AdvNorm.IL(tmp, "Unique State", kol, FALSE, 5)
        m2 <- setNames(vapply(Grps2[w], function(x) {
          mean(tmp[[paste0("Grp", x)]]/tmp2[[paste0("AdvNorm.Grp", x)]], na.rm = TRUE)
        }, 1), kol)
        for (grp2 in Grps2) {
          w2 <- which(ev[wg, Grps2Kol] == grp2)
          ev[wg[w2], ev.col["Advanced normalisation"]] <- ev[wg[w2], ev.col["Normalisation"]]/m2[paste0("Grp", grp2)]
        }
        AdvNorm.Ev[match(grp, AdvNorm.Ev$Group), kol] <- 1/m2
      }
      cat("     Done!\n")
    }
    cat("\n")
  } else {
    cat(" - Skipping MS1-level PSM intensity normalisation since all groups contain exactly 1 MS file.\n\n")
  }
  if (LabelType == "Isobaric") {
    msg <- " - Normalizing Reporter intensities"
    if (length(Iso) > 1) { msg <- c(msg, paste0("   (per ", IsobarLab, " sample)\n")) }
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
    msg <- " - Classic normalisation to the median"
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
    # Per combined sample
    er1 <- ev.ref["Normalisation"] <- paste0("Norm. ", ev.ref["Original"])
    er0 <- ev.ref[match("Normalisation", names(ev.ref))-1]
    k0 <- paste0(er0, get(IsobarLab))
    k1 <- paste0(er1, get(IsobarLab))
    ev[, k1] <- NA
    Norm.Ev.RepIntens %<o% data.frame(Group = Iso)
    for (ch in get(IsobarLab)) { Norm.Ev.RepIntens[[paste0("Channel_", ch)]] <- 1 }
    for (i in Iso) { #i <- Iso[1]
      wg <- which(ev$Isobaric.set == i)
      M3 <- 10^median(is.all.good(log10(unlist(ev[wg, k0])))) # For preserving original scale
      #M <- 10^mlv(is.all.good(log10(unlist(w[wg, ev.col["Original"]]))), method = "Parzen")[1]
      m3 <- vapply(get(IsobarLab), function(ch) { 10^median(is.all.good(log10(ev[wg, paste0(er0, ch)]))) }, 1)
      ev[wg, k1] <- sweep(ev[wg, k0], 2, M3/m3, "*")
      Norm.Ev.RepIntens[match(i, Norm.Ev.RepIntens$Group), paste0("Channel_", get(IsobarLab))] <- m3/M3
    }
    tstAdvNrm <- FALSE
    if (("Adv.Norma.Ev.Intens" %in% colnames(Param))&&(Param$Adv.Norma.Ev.Intens != FALSE)) {
      if (Param$Adv.Norma.Ev.Intens.Type == "C") {
        msg <- c("   - Levenberg-Marquardt normalisation...", "     (please wait)", "     ...")
        ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Print = FALSE)
        cat(paste0(msg, "\n", collapse = "\n"))
        # Per combined sample
        # Because this is computationally expensive, we are doing it per Fraction then averaging:
        # the values should be the same across fractions
        er1 <- ev.ref["Advanced normalisation"] <- paste0("AdvNorm. ", ev.ref["Original"])
        er0 <- ev.ref["Normalisation"]
        k0 <- paste0(er0, get(IsobarLab))
        k1 <- paste0(er1, get(IsobarLab))
        ev[, k1] <- NA
        AdvNorm.Ev.RepIntens %<o% data.frame(Group = Iso)
        tmpEv <- ev[, c("Isobaric.set", "Unique State", k0)]
        if ("Fraction" %in% colnames(ev)) { tmpEv$Fraction <- ev$Fraction } else { tmpEV$Fraction <- 1 }
        clusterCall(parClust, function() library(proteoCraft))
        m4 <- tstRI <- list()
        msg <- "     Estimating normalisation factors within..."
        ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Print = FALSE)
        for (i in Iso) { #i <- Iso[1]
          msg <- paste0("     ...", IsobarLab, " sample ", i, "...")
          ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Print = FALSE)
          wi <- which(tmpEv$Isobaric.set == i)
          #View(tmpEv[wi, k0])
          tmp <- as.data.table(tmpEv[wi, c("Unique State", "Fraction", k0)])
          tmp <- tmp[, lapply(.SD, sum, na.rm = TRUE), keyby = list(`Unique State` = `Unique State`, Fraction = Fraction)]
          tmp <- as.data.frame(tmp)
          Fr <- unique(tmp$Fraction)
          clusterExport(parClust, list("tmp", "k0", "Fr"), envir = environment())
          # Create normalized data
          tmp2 <- setNames(parLapply(parClust, Fr, function(fr) { #fr <- Fr[1]
            dat <- tmp[which(tmp$Fraction == fr),]
            proteoCraft::AdvNorm.IL(dat, "Unique State", k0, FALSE, 5)
          }), paste0("Fr. ", Fr))
          # Compute normalisation factors
          tmp2F <- as.data.frame(sapply(Fr, function(fr) {
            isbrLb <- get(IsobarLab)
            sapply(1:length(isbrLb), function(x) {
              setNames(mean(tmp[which(tmp$Fraction == fr), k0[x]]/(tmp2[[paste0("Fr. ", fr)]][, paste0("AdvNorm.", k0[x])]),
                            na.rm = TRUE), paste0("Ch. ", isbrLb[x]))
            })
          }))
          colnames(tmp2F) <- paste0("Fr. ", Fr)
          stopifnot(!sum(rownames(tmp2F) != paste0("Ch. ", get(IsobarLab))))
          tstRI[[i]] <- tmp2F
          tmp2F$Label <- get(IsobarLab)
          # Number of valid values
          tmp2K <- as.data.frame(sapply(Fr, function(fr) {
            vapply(1:length(get(IsobarLab)), function(x) {
              length(is.all.good(log10(tmp[which(tmp$Fraction == fr), k0[x]])))
            }, 1)
          }))
          colnames(tmp2K) <- paste0("Fr. ", Fr)
          rownames(tmp2F) <- paste0("Ch. ", get(IsobarLab))
          tmp2K$Label <- get(IsobarLab)
          # Calculate weighted mean
          tmp2NormFact <- setNames(vapply(get(IsobarLab), function(x) {
            weighted.mean(tmp2F[match(x, tmp2F$Label), paste0("Fr. ", Fr)],
                          tmp2K[match(x, tmp2K$Label), paste0("Fr. ", Fr)])
          }, 1), paste0("Ch ", get(IsobarLab)))
          tmp2NormFact[which(is.na(tmp2NormFact))] <- 1 # Better not normalize than corrupt data!
          m4[[i]] <- tmp2NormFact
        }
        # Apply results
        m4 <- as.data.frame(t(sapply(m4, unlist)))
        for (i in Iso) {
          wi <- which(tmpEv$Isobaric.set == i)
          ev[wi, k1] <- sweep(ev[wi, k0], 2, unlist(m4[which(rownames(m4) == i), paste0("Ch ", get(IsobarLab))]), "/")
        }
        AdvNorm.Ev.RepIntens[, paste0("Channel_", get(IsobarLab))] <- 1/m4[match(AdvNorm.Ev.RepIntens$Group, rownames(m4)), paste0("Ch ", get(IsobarLab))]
        tstAdvNrm <- TRUE
        cat("     Done!\n")
      } else {
        stop("Not implemented yet, I need to update the current AdvNorm function as it takes too long or crashes.")
      }
    }
    # Combine results of all reporter intensity normalisation steps
    Norm.Ev.RepIntens.All %<o% Norm.Ev.RepIntens
    rownames(Norm.Ev.RepIntens.All) <- Norm.Ev.RepIntens.All$Group; Norm.Ev.RepIntens.All$Group <- NULL
    Norm.Ev.RepIntens.All <- suppressMessages(reshape2::melt(Norm.Ev.RepIntens.All))
    Norm.Ev.RepIntens.All$Iso <- Iso
    Norm.Ev.RepIntens.All$Channel <- as.numeric(gsub("^Channel_", "", Norm.Ev.RepIntens.All$variable))
    Norm.Ev.RepIntens.All$Fraction <- "All"
    if (tstAdvNrm) {
      tmp2tst <- suppressMessages(melt(tstRI))
      colnames(tmp2tst)[which(colnames(tmp2tst) == "L1")] <- "Iso"
      tmp2tst$Channel <- get(IsobarLab)
      tmp2tst$Fraction <- as.numeric(gsub("^Fr\\. ", "", tmp2tst$variable))
      for (i in Iso) {
        for (ch in get(IsobarLab)) {
          w1 <- which((Norm.Ev.RepIntens.All$Iso == i)&(Norm.Ev.RepIntens.All$Channel == ch))
          w2 <- which((tmp2tst$Iso == i)&(tmp2tst$Channel == ch))
          tmp2tst$value[w2] <- tmp2tst$value[w2]*Norm.Ev.RepIntens.All$value[w1]
        }
      }
      tmp2tst$Fraction <- factor(tmp2tst$Fraction, levels = 1:max(tmp2tst$Fraction))
      Norm.Ev.RepIntens.All <- tmp2tst
    }
    Norm.Ev.RepIntens.All$Channel <- factor(Norm.Ev.RepIntens.All$Channel, levels = 0:max(Norm.Ev.RepIntens.All$Channel))
    Norm.Ev.RepIntens.All$Iso <- factor(paste0(IsobarLab, " sample ", Norm.Ev.RepIntens.All$Iso), levels = paste0(IsobarLab, " sample ", Iso))
    w <- (("Adv.Norma.Ev.Intens" %in% colnames(Param))&(Param$Adv.Norma.Ev.Intens != FALSE))+1
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " MS1 intensities and ", IsobarLab, " reporter intensities were re-normalized ",
                             c("to the median", "using the Levenberg-Marquardt procedure to minimize sample-to-sample differences")[w], ".")
    # Visualize results
    dir <- paste0(wd, "/Workflow control/", evNm, "s/Reporter intensities")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    ttl <- "Reporter intensities trend VS fractions"
    plot <- ggplot(Norm.Ev.RepIntens.All) +
      geom_tile(aes(x = Fraction, y = Channel, fill = value), linewidth = 1, width = 1) +
      scale_fill_viridis(begin = 0.25) +
      coord_fixed() + theme_bw() + ylab(paste0(IsobarLab, " channel")) + ggtitle(ttl)
    if (length(Iso) > 1) { plot <- plot + facet_wrap(~Iso) }
    print(plot) # This type of QC plot does not need to pop up, the side panel is fine
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
  } else {
    w <- (("Adv.Norma.Ev.Intens" %in% colnames(Param))&(Param$Adv.Norma.Ev.Intens != FALSE))+1
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " MS1 intensities were re-normalized ",
                             c("to the median", "using the Levenberg-Marquardt procedure to minimize sample-to-sample differences")[w], ".")
  }
  if ((LabelType == "LFQ")&&(Param$Label == "DIA")&&("MS2_intensities" %in% colnames(ev))) {
    # If isobaric, re-scale MS2 intensities to apply MS1 normalisation factors
    # Important: this bit must remain after the normalisation of MS1 intensities
    DatAnalysisTxt <- gsub("\\.$", ", then individual MS2 intensities were re-scaled using MS1 intensities.", DatAnalysisTxt)
    kol <- gsub("Intensity", "MS2 intensities", ev.col[length(ev.col)])
    stopifnot(grepl("MS2 intensities", kol))
    ev[[kol]] <- ev$MS2_intensities # (Let's keep this as a numeric list)
    for (smpl in RSA$values) {
      mqe <- unlist(Exp.map$MQ.Exp[match(smpl, Exp.map$Ref.Sample.Aggregate)])
      w <- which(ev$MQ.Exp %in% mqe)
      mRt <- median(is.all.good(ev[w, ev.col[length(ev.col)]]/ev[w, ev.col["Original"]]))
      clusterExport(parClust, "mRt", envir = environment())
      ev[[kol]][w] <- parLapply(parClust, ev[[kol]][w], function(x) { x*mRt })
    }
    if (Param$Norma.Ev.Intens&&Param$Norma.Ev.Intens.show) {
      kol2 <- unique(c("id", "MQ.Exp", "MS2_intensities", kol))
      tst <- ev[, kol2]
      # Here it is easier to sum per row (otherwise this makes for very slow processing, creates a very huge table and plot, with little added value)
      kolz <- colnames(tst)[which(!colnames(tst) %in% c("id", "MQ.Exp"))]
      for (kl in kolz) {
        if (class(tst[[kl]]) != "list") { tst[[kl]] <- vapply(strsplit(tst[[kl]], ";"), as.numeric, 1) }
        tst[[kl]] <- parSapply(parClust, tst[[kl]], sum)
      }
      tst <- reshape2::melt(tst, id.vars = c("id", "MQ.Exp"))
      tst$value <- log10(tst$value)
      tst$variable <- as.character(tst$variable)
      tst$Norm <- "Original"
      tst$Norm[which(tst$variable == kol)] <- "Normalised"
      tst$Norm <- factor(tst$Norm, levels = c("Original", "Normalised"))
      ttl <- paste0(evNm, "s intensity normalisation")
      dir <- paste0(wd, "/Workflow control/", evNm, "s/Normalisation")
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      plot <- ggplot(tst) +
        geom_violin(aes(x = MQ.Exp, y = value, color = Norm, fill = Norm), alpha = 0.25) +
        geom_boxplot(aes(x = MQ.Exp, y = value, color = Norm, fill = Norm), alpha = 0.5) +
        scale_color_viridis_d(begin = 0.25) +
        scale_fill_viridis_d(begin = 0.25) +
        facet_wrap(~Norm, scales = "free") + theme_bw() + ggtitle(ttl) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      print(plot) # This type of QC plot does not need to pop up, the side panel is fine
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ReportCalls <- AddPlot2Report()
    }
  }
}
# If isobaric, re-scale reporter intensities to total evidence intensities:
# Important: this bit must remain after the normalization of MS1 intensities
if (LabelType == "Isobaric") {
  DatAnalysisTxt <- gsub("\\.$", ", then reporter intensities were re-scaled using MS1 intensities.", DatAnalysisTxt)
  ev.ref["Adjusted"] <- paste0("adj. ", ev.ref["Original"])
  k0 <- paste0(ev.ref[match("Adjusted", names(ev.ref))-1], get(IsobarLab))
  k1 <- paste0(ev.ref["Adjusted"], get(IsobarLab))
  temp <- rowSums(ev[, k0], na.rm = TRUE)
  ev[, k1] <- sweep(ev[, k0], 1, ev[, ev.col[length(ev.col)]]/temp, "*")
  if (Param$Norma.Ev.Intens&&Param$Norma.Ev.Intens.show) {
    er0 <- ev.ref[match("Normalisation", names(ev.ref))-1]
    er1 <- ev.ref["Adjusted"]
    a0 <- paste0(er0, get(IsobarLab))
    a1 <- paste0(er1, get(IsobarLab))
    test <- as.data.table(ev[, c("MQ.Exp", a0, a1)])
    test <- data.table::melt(test, id.vars = "MQ.Exp")
    test <- as.data.frame(test)
    test$Norm <- NA
    test$Norm[which(test$variable %in% a0)] <- "Original"
    test$Norm[which(test$variable %in% a1)] <- "Normalised"
    test$Norm <- factor(test$Norm, levels = c("Original", "Normalised"))
    test$value <- log10(test$value)
    #aggregate(temp, list(ev$MQ.Exp), function(x) { sum(!is.na(x)) })
    #aggregate(test$value, list(test$MQ.Exp, test$Norm), function(x) { sum(!is.na(x)) })
    ttl <- paste0(evNm, "s reporter intensity re-scaling")
    dir <- paste0(wd, "/Workflow control/", evNm, "s/Normalisation")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    test$Channel <- as.numeric(gsub(topattern(c(er0, er1)), "", as.character(test$variable)))
    test$Channel <- factor(test$Channel, levels = sort(unique(test$Channel)))
    plot <- ggplot(test) +
      geom_violin(aes(x = Channel, y = value, color = Channel, fill= Channel), alpha = 0.25) +
      geom_boxplot(aes(x = Channel, y = value, color = Channel, fill = Channel), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(MQ.Exp ~ Norm, scales = "free", space = "free") + theme_bw() + ggtitle(ttl) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(plot) # This type of QC plot does not need to pop up, the side panel is fine
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ReportCalls <- AddPlot2Report()
  }
  #Isobaric data: valid values
  kol <- grep(topattern(ev.ref["Original"]), colnames(ev), value = TRUE)
  kol <- grep(" count ", kol, value = TRUE, invert = TRUE)
  tst1 <- sapply(MQ.Exp, function(x) {
    vapply(gsub(topattern(ev.ref["Original"]), "", kol), function(y) {
      w <- which(vapply(Exp.map$MQ.Exp, function(z) { x %in% unlist(z) }, TRUE)&(Exp.map$"Isobaric label" == y))
      stopifnot(length(w) <= 1)
      if (length(w)) { res <- cleanNms(Exp.map$Ref.Sample.Aggregate[w]) } else { res <- "" }
      return(res)
    }, "")
  })
  tst2 <- set_rownames(sapply(MQ.Exp, function(x) {
    vapply(kol, function(y) {
      sum(ev[which(ev$MQ.Exp == x), y] > 0)
    }, 1)
  }), rownames(tst1))
  tst3 <- set_rownames(sapply(MQ.Exp, function(x) {
    vapply(kol, function(y) {
      round(median(is.all.good(log10(ev[which(ev$MQ.Exp == x), y]))), 2)
    }, 1)
  }), rownames(tst1))
  tst <- rbind(rep("", length(MQ.Exp)),
               c("Samples", rep("", length(MQ.Exp)-1)),
               colnames(tst1),
               tst1,
               rep("", length(MQ.Exp)),
               c("Number of valid values", rep("", length(MQ.Exp)-1)), 
               colnames(tst2),
               tst2,
               rep("", length(MQ.Exp)),
               c("Median log10", rep("", length(MQ.Exp)-1)), 
               colnames(tst3),
               tst3,
               rep("", length(MQ.Exp)))
  colnames(tst) <- NULL
  data.table::fwrite(tst, paste0(wd, "/Workflow control/Valid values test.csv"), quote = FALSE, sep = ",", col.names = FALSE, na = "NA")
  #system(paste0("open \"", wd, "/Workflow control/Valid values test.csv\""))
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              "body_add_fpar(Report, fpar(ftext(\"Number of valid values per sample and channel:\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
  ReportCalls$Objects$Valid_values <- as.data.frame(tst)
  ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_table(Report, ReportCalls$Objects$Valid_values)")
  ReportCalls <- AddSpace2Report()
  ## To do here:
  ## - Format table
  ## - Add fraction-specificity
}
