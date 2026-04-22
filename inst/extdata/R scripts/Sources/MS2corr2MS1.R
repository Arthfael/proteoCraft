# DIA only MS2-based correction of MS1-based quantitative values.
# Only used if QuantUMS was off for DiaNN
# Use MS2-based measurements to improve estimates of precursor abundance.
# In the future (and also for Prot.Quant) we should take the data's noisiness into account, i.e.,
# "what is the confidence we have on the averaged MS2-based profile?"
MS2_based_Correction %<o% TRUE
# A column named MS2_intensities with individual MS2 fragment intensities should be present!
if ((LabelType == "LFQ")&&(sum(isDIA))&&("MS2_intensities" %in% colnames(ev))) { # We only run if we are in DIA mode and...
  # ... we are either not using DiaNN or we are but we did not run QuantUMS
  if ((MS2_based_Correction)&&((SearchSoft != "DIANN")||(!QuantUMS))) {
    msg <- "Refining MS1-level measurements using MS2 data..."
    ReportCalls <- AddMsg2Report(Space = FALSE)
    dir <- paste0(wd, "/Workflow control/", evNm, "s/MS2-based MS1 correction")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    if (exists("dirlist")) { dirlist <- unique(c(dirlist, dir)) }
    #
    if (scrptType == "withReps") {
      ref <- if ("MS2-corrected" %in% names(ev.col)) { ev.col[match("MS2-corrected", names(ev.col))-1L] } else {
        ev.col[length(ev.col)]
      }
      nuRef <- ev.col["MS2-corrected"] <- paste0("MS2-corr. ", ev.col["Original"])
    }
    if (scrptType == "noReps") {
      ref <- if ("MS2-corrected" %in% names(int.cols)) { int.cols[match("MS2-corrected", names(int.cols))-1L] } else {
        int.cols[length(int.cols)]
      }
      nuRef <- int.col <- int.cols["MS2-corrected"] <- paste0("MS2-corr. ", int.cols["Original"])
    }
    #
    if (!"Search_ID" %in% colnames(ev)) {
      stopifnot(length(inDirs) == 1L) # If length(inDirs) > 1 then column Search_ID should be always present!!!
      ev$Search_ID <- inDirs
      names(isDIA) <- inDirs
    }
    GrpsVct <- ev$Search_ID
    Grps <- names(isDIA)[which(isDIA)]
    Grps <- Grps[which(Grps %in% GrpsVct)]
    chckMS2Corr <- FALSE
    ev[[nuRef]] <- ev[[ref]]
    for (grp in Grps[]) { #grp <- Grps[1L]
      grpMtch <- match(grp, Grps)
      w <- which(GrpsVct == grp)
      MS2Tbl <- data.table(IDs = ev$id[w], MS1 = ev[w, ref], MS2 = ev$MS2_intensities[w],
                           mod = ev$"Modified sequence"[w], Z = ev$Charge[w])
      MS2Tbl <- MS2Tbl[, list(IDs = list(IDs), MS2 = list(MS2), MS1 = list(MS1),
                              MS1_Av = mean(MS1, na.rm = TRUE)),
                       by = list(`Modified sequence` = MS2Tbl$mod, Charge = MS2Tbl$Z)]
      MS2Tbl <- as.data.frame(MS2Tbl)
      source(parSrc, local = FALSE)
      clusterExport(parClust, "LFQ.lm", envir = environment())
      tst1 <- parSapply(parClust, MS2Tbl$IDs, length)
      wMult <- which(tst1 > 1L)
      MS2Tbl <- MS2Tbl[wMult,]
      MS2Tbl$FiltMS2 <- parLapply(parClust, MS2Tbl$MS2, \(x) {
        #x <- MS2Tbl$MS2[[1L]]
        # if (length(x) == 1L) {
        #   x <- magrittr::set_colnames(data.frame(A = unlist(x)), NULL)
        # } else {
        L <- min(sapply(x, length))
        x <- as.data.frame(sapply(x, \(y) { y[1L:L] })) # Crop if necessary (shouldn't be the case)
        # Re-order fragments
        tst <- rowMeans(x, na.rm = TRUE)
        x <- x[order(tst, decreasing = TRUE),]
        # }
        return(x)
      })
      # MS2Tbl$L <- parSapply(parClust, MS2Tbl$MS2, \(x) {
      #   min(sapply(x, length))
      # })
      MS2Tbl$CorrVal <- parLapply(parClust, MS2Tbl$FiltMS2, \(x) { #x <- MS2Tbl$FiltMS2[[1L]]
        Dat <- x
        # if (ncol(x) > 1L) {
        k <- colnames(Dat)
        Dat$id <- 1L:nrow(Dat)
        Dat$Weights <- rowSums(Dat[, k], na.rm = TRUE)
        #DefArg(LFQ.lm); ids = Dat$id; InputTabl = Dat; IntensCol = k; Summary.method = "weighted.mean";Summary.weights = "Weights"; Min.N = 1L; Max.N = 6L;Is.log = FALSE
        x <- LFQ.lm(Dat$id,
                    InputTabl = Dat,
                    IntensCol = k,
                    Summary.method = "weighted.mean",
                    Summary.weights = "Weights",
                    Min.N = 1L,
                    Max.N = 6L,
                    Is.log = FALSE)
        # } else { x <- Dat[1L, 1L] }
        return(x)
      })
      MS2Tbl$Av <- parSapply(parClust, MS2Tbl$CorrVal, mean)
      MS2Tbl$CorrVal_Scl <- parApply(parClust, MS2Tbl[, c("CorrVal", "Av", "MS1_Av")], 1L, \(x) {
        x[[1L]]*x[[3L]]/x[[2L]]
      })
      MS2Tbl$CorrVal2 <- parApply(parClust, MS2Tbl[, c("IDs", "MS1", "CorrVal_Scl")], 1L, \(x) {
        data.frame(ID = x[[1L]], OrigVal = x[[2L]], CorrVal = x[[3L]])
      })
      tst <- plyr::rbind.fill(MS2Tbl$CorrVal2)
      tst$Ratio1 <- tst$CorrVal/tst$OrigVal
      m <- match(tst$ID, ev$id)
      tst[, c("PEP", "Sample", "ModSeq", "Z", "Proteins", "Seq")] <- ev[m, c("PEP", "Raw file path", "Modified sequence",
                                                                             "Charge", "Proteins", "Sequence")]
      tst$Weights <- -log10(tst$PEP)/5
      tst$Weights[which(!is.all.good(tst$Weights, 2L))] <- 3
      tst$Weights[which(tst$Weights > 3)] <- 3
      tst$CorrVal2 <- (tst$OrigVal+tst$CorrVal*tst$Weights)/(1+tst$Weights)
      tst$Ratio2 <- tst$CorrVal2/tst$OrigVal
      # Add corrected data to ev
      w <- which(ev$id %in% tst$ID)
      ev[w, nuRef] <- tst$CorrVal2[match(ev$id[w], tst$ID)]
      # Additional stuff for plotting
      tst$PepL <- nchar(tst$Seq)
      tmp <- listMelt(MS2Tbl$IDs, MS2Tbl$L)
      tst$N_of_fragments <- tmp$L1[match(tst$ID, tmp$value)]
      tst$N_of_fragments <- factor(tst$N_of_fragments, levels = 1L:max(tst$N_of_fragments))
      tst$MQ.Exp <- Frac.map$MQ.Exp[match(tst$Sample, Frac.map$"s name")]
      tmp <- listMelt(Exp.map$MQ.Exp, Exp.map$Ref.Sample.Aggregate)
      tmp$L1 <- cleanNms(tmp$L1)
      tst$Sample_name <- tmp$L1[match(tst$MQ.Exp, tmp$value)]
      tmp <- listMelt(Exp.map$MQ.Exp, Exp.map[[VPAL$column]])
      tmp$L1 <- cleanNms(tmp$L1)
      tst$Samples_group <- tmp$L1[match(tst$MQ.Exp, tmp$value)]
      ttl <- paste0("Group ", grpMtch, " - Original MS1 vs MS2-corrected MS1 (PEP)")
      scl <- 2L*(ceiling(sqrt(length(unique(tst$Sample))))+1L)
      plot <- ggplot(tst) + geom_point(aes(x = log10(OrigVal), y = log10(CorrVal), colour = -log10(PEP)),
                                       size = 0.01, alpha = 0.1, shape = 16L) + coord_fixed() +
        theme_bw() + ggtitle(ttl, subtitle = grp) + facet_wrap(~Sample)
      #poplot(plot)
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150L, width = scl, height = scl, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150L, width = scl, height = scl, units = "in")
      # ttl <- paste0("Group ", grpMtch, " - Original MS1 vs MS2-corrected MS1 (Nb. of fragments)")
      # plot <- ggplot(tst) + geom_point(aes(x = log10(OrigVal), y = log10(CorrVal), colour = N_of_fragments),
      #                                  size = 0.01, alpha = 0.1, shape = 16) + coord_fixed() +
      #   theme_bw() + ggtitle(ttl, subtitle = grp) + facet_grid(N_of_fragments~Sample)
      # poplot(plot)
      # ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300)
      # ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)
      ttl <- paste0("Group ", grpMtch, " - Original MS1 vs final corrected MS1 (PEP)")
      plot <- ggplot(tst) +
        geom_scattermore(aes(x = log10(OrigVal), y = log10(CorrVal2), colour = -log10(PEP)),
                         size = 0.01, alpha = 0.1, shape = 16L) + coord_fixed() +
        theme_bw() + ggtitle(ttl, subtitle = grp) + facet_wrap(~Sample)
      #poplot(plot)
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150L, width = scl, height = scl, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150L, width = scl, height = scl, units = "in")
      if (chckMS2Corr) {
        tst0 <- data.table(Samples_group = tst$Samples_group, Sample = tst$Sample_name,
                           ModSeq = tst$ModSeq, Charge = tst$Z,
                           Original = tst$OrigVal, Corrected = tst$CorrVal2)
        tst1 <- dcast(tst0, ModSeq+Charge~Sample, value.var = "Original", fun.aggregate = sum)
        tst2 <- dcast(tst0, ModSeq+Charge~Sample, value.var = "Corrected",  fun.aggregate = sum)
        tst1 <- as.data.frame(tst1)
        tst2 <- as.data.frame(tst2)
        Smpls <- unique(tst0$Sample)
        smplGrps <- setNames(tst0$Samples_group[match(Smpls, tst0$Sample)], Smpls)
        Smpls2Grps <- setNames(aggregate(Smpls, list(smplGrps), c), c("Group", "Sample"))
        # Quick median normalisation
        M1 <- median(as.matrix(tst1[, Smpls]), na.rm = TRUE)
        M2 <- median(as.matrix(tst2[, Smpls]), na.rm = TRUE)
        m1 <- sapply(Smpls, \(smpl) { median(tst1[[smpl]], na.rm = TRUE)})
        m2 <- sapply(Smpls, \(smpl) { median(tst2[[smpl]], na.rm = TRUE)})
        tst1[, Smpls] <- sweep(tst1[, Smpls], 2L, M1/m1, "*")
        tst2[, Smpls] <- sweep(tst2[, Smpls], 2L, M2/m2, "*")
        # Additional row normalisation
        tst1[, Smpls] <- sweep(tst1[, Smpls], 1L, rowMeans(tst1[, Smpls], na.rm = TRUE), "/")
        tst2[, Smpls] <- sweep(tst2[, Smpls], 1L, rowMeans(tst2[, Smpls], na.rm = TRUE), "/")
        #
        #
        intraCVs_Orig <- setNames(sapply(Smpls2Grps$Sample, \(smpls) {
          x <- unlist(tst1[, smpls])
          sd(x)/mean(x)
        }), Smpls2Grps$Group)
        intraCVs_Corr <- setNames(sapply(Smpls2Grps$Sample, \(smpls) {
          x <- unlist(tst2[, smpls])
          sd(x)/mean(x)
        }), Smpls2Grps$Group)
        tmp <- unlist(tst1[, Smpls])
        globalCV_Orig <- sd(tmp)/mean(tmp)
        tmp <- unlist(tst2[, Smpls])
        globalCV_Corr <- sd(tmp)/mean(tmp)
      }
      # 
    }
  }
}
