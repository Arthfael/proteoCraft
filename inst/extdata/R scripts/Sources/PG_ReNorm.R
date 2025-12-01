#### Re-normalize protein group expression values
Norma.Prot.Ratio.classic %<o% FALSE
Norma.Prot.Ratio.to.Biot %<o% FALSE
Norma.Prot.Ratio.to.proteins %<o% FALSE
if (("Norma.Prot.Ratio" %in% colnames(Param))&&(Param$Norma.Prot.Ratio)) {
  quant.data.norm <- quant.data
  ExpKol <- paste0(Prot.Expr.Root, Exp.map$Ref.Sample.Aggregate)
  RatKol <- paste0(Prot.Rat.Root, Exp.map$Ref.Sample.Aggregate)
  whichRat <- which(RatKol %in% colnames(quant.data.norm))
  RatKol <- RatKol[whichRat]
  RefKol <- paste0(Prot.Expr.Root, RRG$values, ".REF")
  Norma.Prot.Ratio.classic <- TRUE # The default
  if ((IsBioID2)&&("Norma.Prot.Ratio.to.Biot" %in% colnames(Param))&&(is.logical(Param$Norma.Prot.Ratio.to.Biot))&&(Param$Norma.Prot.Ratio.to.Biot)) {
    nrmFlt <- which(PG$"Biot. peptides count" > 0)
    if (length(w)) {
      Norma.Prot.Ratio.classic <- Norma.Prot.Ratio.to.proteins <- FALSE
      Norma.Prot.Ratio.to.Biot <- TRUE
      txTmp <- "biotinylated protein groups"
    }
  }
  if (("Norma.Prot.Ratio.to.proteins" %in% colnames(Param))&&
      (is.character(Param$Norma.Prot.Ratio.to.proteins))&&
      (!Param$Norma.Prot.Ratio.to.proteins %in% c("", "NA"))&&
      (!Norma.Prot.Ratio.to.Biot)) {
    Prot.Ratio.ref.Acc %<o% unique(unlist(strsplit(Param$Norma.Prot.Ratio.to.proteins, ";")))
    Prot.Ratio.ref.Acc <- Prot.Ratio.ref.Acc[which(Prot.Ratio.ref.Acc %in% db$"Protein ID")]
    if (length(Prot.Ratio.ref.Acc)) {
      Prot.Ratio.ref.Acc <- Prot.Ratio.ref.Acc[which(Prot.Ratio.ref.Acc %in% unique(unlist(strsplit(PG$"Protein IDs", ";"))))] # Having "Leading protein IDs" here was too stringent
      Norma.Prot.Ratio.classic <- FALSE
      l <- length(Prot.Ratio.ref.Acc)
      if (l) {
        temp <- strsplit(PG$"Protein IDs", ";")
        temp <- listMelt(temp, PG$id)
        temp$Norm <- temp$value %in% Prot.Ratio.ref.Acc
        temp <- temp$L1[which(temp$Norm)]
        nrmFlt <- which(PG$id %in% temp)
        if (length(nrmFlt)) {
          Norma.Prot.Ratio.classic <- Norma.Prot.Ratio.to.Biot <- FALSE
          Norma.Prot.Ratio.to.proteins <- TRUE
          txTmp <- paste0("the following proteins: ",
                          paste(db$`Full Name`[match(Prot.Ratio.ref.Acc, db$`Protein ID`)],
                                collapse = "-"))
        }
      }
    }
  }
  if (Norma.Prot.Ratio.classic) {
    nrmFlt <- 1:nrow(quant.data)
    Norma.Prot.Ratio.to.Biot <- Norma.Prot.Ratio.to.proteins <- FALSE
    txTmp <- "all protein groups"
  }
  #
  # nrmFlt is the rows filter we will use to ensure that we normalize to the correct proteins:
  # - Norma.Prot.Ratio.classic: the filter includes all rows
  # - Norma.Prot.Ratio.to.Biot: the filter includes biotinylated proteins
  # - Norma.Prot.Ratio.to.proteins: the filter includes specific user-defined protein(s)
  #
  tst <- sum(c(Norma.Prot.Ratio.to.Biot, Norma.Prot.Ratio.to.proteins, Norma.Prot.Ratio.classic))
  if (tst) {
    stopifnot(tst == 1)
    globMed1 <- median(is.all.good(unlist(quant.data.norm[#nrmFlt # For global scales, we use all
      , c(ExpKol, RefKol)])))
    # First normalize within groups
    for (grp in RG$values) { #grp <- RG$values[1]
      smpls <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[RG$column]] == grp)]
      xpKol <- paste0(Prot.Expr.Root, smpls)
      rtKol <- paste0(Prot.Rat.Root, smpls)
      wRat <- which(rtKol %in% colnames(quant.data.norm))
      rtKol <- rtKol[wRat]
      rfgrps <- unique(Exp.map[which(Exp.map[[RG$column]] == grp), RRG$column])
      rfKol <- paste0(Prot.Expr.Root, rfgrps, ".REF")
      rtMed <- apply(quant.data.norm[nrmFlt, rtKol], 2, function(x) { median(is.all.good(x)) })
      xpMed1 <- median(is.all.good(unlist(quant.data.norm[nrmFlt, xpKol]))) # Before, for re-scaling
      # 1a - apply normalization at ratios level...
      quant.data.norm[, rtKol] <- sweep(quant.data.norm[, rtKol, drop = FALSE],
                                        2,
                                        rtMed,
                                        "-")
      # 1b - ... and probagate it to the corresponding columns
      quant.data.norm[, xpKol[wRat]] <- sweep(quant.data.norm[, xpKol[wRat], drop = FALSE],
                                              2,
                                              rtMed/log2(10), # Convert from base 2 to 10 
                                              "-")
      # 2 - restore original scaling
      xpMed2 <- median(is.all.good(unlist(quant.data.norm[nrmFlt, xpKol]))) # After, for re-scaling
      quant.data.norm[, c(xpKol, rfKol)] <- sweep(quant.data.norm[, c(xpKol, rfKol)],
                                                  2,
                                                  xpMed1-xpMed2,
                                                  "+")
    }
    # ... then re-normalize globally.
    globMed2 <- median(is.all.good(unlist(quant.data.norm[#nrmFlt
      , c(ExpKol, RefKol)])))
    quant.data.norm[, c(ExpKol, RefKol)] <- sweep(quant.data.norm[, c(ExpKol, RefKol)],
                                                  2,
                                                  globMed1-globMed2,
                                                  "+")
    # Sometimes, if you are recycling to e.g. a single protein, you may get randomly singular values in one row
    # This breaks stats so we need to fix that.
    #  - Expression:
    tst <- apply(quant.data.norm[, c(ExpKol, RefKol)], 1, function(x) {
      x <- x[which(!is.na(x))]
      length(unique(x))
    })
    if (1 %in% tst) {
      #aggregate(tst, list(tst), length)
      # Add a small value
      w <- which(tst > 0)
      quant.data.norm[w, c(ExpKol, RefKol)] <- quant.data.norm[w, c(ExpKol, RefKol)] + runif(length(w)*length(c(ExpKol, RefKol)), min = 0, max = 1e-10)
    }
    #  - Ratios:
    tst <- apply(quant.data.norm[, RatKol], 1, function(x) {
      x <- x[which(!is.na(x))]
      length(unique(x))
    })
    if (1 %in% tst) {
      #aggregate(tst, list(tst), length)
      # Add a small value
      w <- which(tst > 0)
      quant.data.norm[w, RatKol] <- quant.data.norm[w, RatKol] + runif(length(w)*length(RatKol), min = 0, max = 1e-10)
    }
    #
    ttl <- paste0("re-normalisation (median of ",
                  c("Biot. PGs", "specific proteins", "all PGs")[which(c(Norma.Prot.Ratio.to.Biot,
                                                                         Norma.Prot.Ratio.to.proteins,
                                                                         Norma.Prot.Ratio.classic))],
                  ")")
    DatAnalysisTxt <- gsub("\\.\\.\\.$", paste0(" and re-normalized to the median value of ", txTmp, "."), DatAnalysisTxt)
    ReportCalls$Calls <- AddTxt2Report(paste0("Re-normalizing Protein groups-level expression values  the median value of ",
                                       txTmp))
  }
  if (Norma.Prot.Ratio.classic) {
    if (!"Adv.Norma.Prot.Intens" %in% colnames(Param)) {
      ObjNm <- "Norma.Prot.Ratio.Adv"
      if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
        msg <- "Do you want to apply Levenberg-Marquardt (slow) re-normalisation at protein groups level?"
        tmp <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
        if (is.na(tmp)) { tmp <- FALSE }
        ObjNm %<c% tmp
        tmp <- AllAnsw[1,]
        tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
        tmp$Value <- list(get(ObjNm))
        m <- match(ObjNm, AllAnsw$Parameter)
        if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
      }
    } else { Norma.Prot.Ratio.Adv <- Param$Adv.Norma.Prot.Intens }
    if (Norma.Prot.Ratio.Adv) {
      ReportCalls$Calls <- AddTxt2Report("Applying Levenberg-Marquardt normalisation at Protein groups-level...")
      ids <- PG$id[nrmFlt]
      # Let's first normalize per ratio group...
      clusterExport(parClust, list("quant.data.norm", "Exp.map", "RG", "RRG", "nrmFlt", "is.all.good",
                                   "Prot.Expr.Root", "Prot.Rat.Root", "ids"), envir = environment())
      tmpNrm <- parLapply(parClust, RG$values, function(grp) { #grp <- RG$values[1]
        smpls <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[RG$column]] == grp)]
        xpKol <- paste0(Prot.Expr.Root, smpls)
        rfgrps <- unique(Exp.map[which(Exp.map[[RG$column]] == grp), RRG$column])
        rfKol <- paste0(Prot.Expr.Root, rfgrps, ".REF")
        grpMed1 <- median(is.all.good(unlist(quant.data.norm[#nrmFlt
          , c(xpKol, rfKol)])))
        tmp <- quant.data.norm[nrmFlt, xpKol]
        tmp$id <- ids
        xpMed1 <- apply(quant.data.norm[nrmFlt, xpKol], 2, median, na.rm = TRUE)
        tmp2 <- proteoCraft::AdvNorm.IL(tmp, "id", xpKol, TRUE, 5)
        tmp2$id <- NULL
        colnames(tmp2) <- xpKol
        xpMed2 <- apply(tmp2, 2, median, na.rm = TRUE)
        tmp2 <- sweep(quant.data.norm[, xpKol],
                      2, 
                      xpMed2-xpMed1,
                      "+")
        # References
        rrg <- RRG$values[which(vapply(RRG$values, function(x) {
          unique(Exp.map[which(Exp.map[[RRG$column]] == x), RG$column] == grp)
        }, TRUE))]
        #xpMed2-xpMed1
        rfMed1 <- apply(quant.data.norm[nrmFlt, rfKol, drop = FALSE], 2, median, na.rm = TRUE)
        rfMed2 <- vapply(rrg, function(x) { #x <- rrg[1]
          x <- Exp.map$Ref.Sample.Aggregate[which((Exp.map[[RRG$column]] == x)&(Exp.map$Reference))]
          median(rowMeans(tmp2[, paste0(Prot.Expr.Root, x), drop = FALSE], na.rm = TRUE), na.rm = TRUE)
        }, 1)
        #rfMed2-rfMed1
        tmp3 <- sweep(quant.data.norm[, rfKol, drop = FALSE],
                      2, 
                      rfMed2-rfMed1,
                      "+")
        resXp <- cbind(tmp2, tmp3)
        # Re-apply original scaling
        grpMed2 <- median(is.all.good(unlist(resXp[#nrmFlt
          , c(xpKol, rfKol)])))
        resXp <- sweep(resXp,
                       2, 
                       grpMed1-grpMed2,
                       "+")
        return(resXp)
      })
      tmpNrm2 <- tmpNrm
      tmpNrm <- do.call(cbind, tmpNrm)
      quant.data.norm[, colnames(tmpNrm)] <- tmpNrm
      if (length(RG$values) > 1) {
        # ... then between ratio groups...
        tmpNrm2 <- as.data.frame(sapply(tmpNrm2, function(x) { #x <- tmpNrm[[1]]
          k <- grep("\\.REF$", colnames(x), invert = TRUE, value = TRUE)
          rowMeans(x[nrmFlt, k], na.rm = TRUE)
        }))
        colnames(tmpNrm2) <- RG$values
        globMed1 <- median(is.all.good(unlist(tmpNrm2)))
        xpMed1 <- apply(tmpNrm2, 2, median, na.rm = TRUE)
        tmpNrm2$id <- ids
        tmp2 <- proteoCraft::AdvNorm.IL(tmpNrm2, "id", RG$values, TRUE, 5)
        tmp2$id <- NULL
        colnames(tmp2) <- RG$values
        xpMed2 <- apply(tmp2, 2, median, na.rm = TRUE)
        globMed2 <- median(is.all.good(unlist(tmp2)))
        xpMed1 <- c(sapply(Exp.map[[RG$column]], function(x) { xpMed1[x] }),
                    sapply(RRG$values, function(x) {
                      xpMed1[Exp.map[match(x, Exp.map[[RRG$column]]) # Only works because RRGs are always contained in RGs
                                     , RG$column]]
                    }))
        xpMed2 <- c(sapply(Exp.map[[RG$column]], function(x) { xpMed2[x] }),
                    sapply(RRG$values, function(x) {
                      xpMed2[Exp.map[match(x, Exp.map[[RRG$column]]) # Only works because RRGs are always contained in RGs
                                     , RG$column]]
                    }))
        quant.data.norm[, c(ExpKol, RefKol)] <- sweep(quant.data.norm[, c(ExpKol, RefKol)],
                                                      2, xpMed2+globMed1-xpMed1-globMed2,
                                                      "+")
      }
      # Re-calculate ratios
      Rat <- as.data.frame(sapply(Exp.map$Ref.Sample.Aggregate[whichRat], function(x) {
        #x <- Exp.map$Ref.Sample.Aggregate[whichRat][1]
        xp <- paste0(Prot.Expr.Root, x)
        rf <- paste0(Prot.Expr.Root, Exp.map[match(x, Exp.map$Ref.Sample.Aggregate), RRG$column], ".REF")
        x1 <- quant.data.norm[[xp]]
        x0 <- quant.data.norm[[rf]]
        return((x1-x0)/log10(2))
      }))
      colnames(Rat) <- RatKol
      quant.data.norm[, RatKol] <- Rat
      #
      ttl <- "re-normalisation (Levenberg-Marquardt)"
      DatAnalysisTxt <- gsub(" and re-normalized to the median value of .+",
                             paste0(" and re-normalized using the Levenberg-Marquardt procedure."), DatAnalysisTxt)
      ReportCalls$Calls <- AddTxt2Report("Done!")
    }
  }
  if (Norma.Prot.Ratio.classic||Norma.Prot.Ratio.to.Biot||Norma.Prot.Ratio.to.proteins) {
    # Once this is done, update ref-to-ref ratios
    res <- make_RefRat(data = quant.data.norm,
                       int.root = Prot.Expr.Root,
                       rat.root = Prot.Rat.Root,
                       logInt = 10)
    quant.data.norm[, colnames(res)] <- res
    #... ok, it's done, but what about peptides?
    # If we want to normalize, e.g., peptide fold changes to parent protein fold changes,
    # then we must first propagate the normalisations to the peptides
    kPpE <- c(paste0(pep.ref[length(pep.ref)], Exp.map$Ref.Sample.Aggregate),
              paste0(pep.ref[length(pep.ref)], unique(Exp.map[[RRG$column]]), ".REF"))
    kPpR <- paste0(pep.ratios.ref[length(pep.ratios.ref)], Exp.map$Ref.Sample.Aggregate)
    kPrE <- c(paste0(Prot.Expr.Root, Exp.map$Ref.Sample.Aggregate),
              paste0(Prot.Expr.Root, unique(Exp.map[[RRG$column]]), ".REF"))
    kPrR <-  paste0(Prot.Rat.Root, Exp.map$Ref.Sample.Aggregate)
    kPpR <- kPpR[which(kPpR %in% colnames(pep))]
    kPrR <- kPrR[which(kPrR %in% colnames(quant.data))]
    pep.quant.data.norm <- pep[, c(kPpE, kPpR)]
    diffE <- 10^vapply(kPrE, function(x) {
      mean(is.all.good(unlist(quant.data.norm[x] - quant.data[x])))
    }, 1)
    diffR <- vapply(kPrR, function(x) {
      mean(is.all.good(unlist(quant.data.norm[x] - quant.data[x])))
    }, 1)
    pep.quant.data.norm[, kPpE] <- sweep(pep.quant.data.norm[, kPpE], 2, diffE, "*")
    pep.quant.data.norm[, kPpR] <- sweep(pep.quant.data.norm[, kPpR], 2, diffR, "+")
    # Also pep.ref.ratios
    if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
      pep.Ref.Ratios.norm %<o% NULL
    }
    if (Param$Ratios.Thresholds == threshMsg) {
      pep.Ref.Ratios.norm %<o% make_RefRat(data = pep.quant.data.norm)
    }
    #
    # And now check and visualize
    dirPG <- paste0(wd, "/Workflow control/Protein groups")
    if (!dir.exists(dirPG)) { dir.create(dirPG, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dirPG))
    dirPep <- paste0(wd, "/Workflow control/Peptides")
    if (!dir.exists(dirPep)) { dir.create(dirPep, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dirPep))
    nrmPlots <- list()
    # - PG expression values
    Samples <- cleanNms(gsub(topattern(Prot.Expr.Root), "", ExpKol))
    temp1 <- data.table(quant.data[, ExpKol])
    colnames(temp1) <- Samples
    temp1 <- melt(temp1, measure.vars = Samples)
    temp2 <- data.table(quant.data.norm[, ExpKol])
    colnames(temp2) <- Samples
    temp2 <- melt(temp2, measure.vars = Samples)
    temp1$Norm <- "Original"
    temp2$Norm <- "Normalised"
    temp <- rbind(temp1, temp2)
    rm(temp1, temp2)
    temp$Sample <- factor(temp$variable, levels = Samples)
    temp$variable <- NULL
    temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
    temp <- temp[which(is.all.good(temp$value, 2)),]
    ttlI1 <- paste0("Protein groups ", ttl, ", expression")
    intPlot1 <- ggplot(temp) +
      geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(Norm~.) + ggtitle(ttlI1) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    #print(intPlot1)
    nrmPlots[["PG_int"]] <- list(Path = paste0(dirPG, "/", ttlI1),
                                 Plot = proteoCraft::plotEval(intPlot1),
                                 Ext = "jpeg")
    #suppressMessages({
    #  ggsave(paste0(dirPG, "/", ttlI1, ".jpeg"), intPlot1, dpi = 300)
    #  ggsave(paste0(dirPG, "/", ttlI1, ".pdf"), intPlot1, dpi = 300)
    #})
    #ReportCalls <- AddPlot2Report(Plot = intPlot1, Title = ttlI1)
    #
    # - PG ratios
    Samples <- cleanNms(gsub(topattern(Prot.Rat.Root), "", RatKol))
    temp1 <- data.table(quant.data[, RatKol])
    colnames(temp1) <- Samples
    temp1 <- melt(temp1, measure.vars = Samples)
    temp2 <- data.table(quant.data.norm[, RatKol])
    colnames(temp2) <- Samples
    temp2 <- melt(temp2, measure.vars = Samples)
    temp1$Norm <- "Original"
    temp2$Norm <- "Normalised"
    temp <- rbind(temp1, temp2)
    rm(temp1, temp2)
    temp$Sample <- factor(temp$variable, levels = Samples)
    temp$variable <- NULL
    temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
    temp <- temp[which(is.all.good(temp$value, 2)),]
    ttlR1 <- paste0("Protein groups ", ttl, ", ratios")
    ratPlot1 <- ggplot(temp) +
      geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(Norm~.) + ggtitle(ttlR1) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    #print(ratPlot1)
    nrmPlots[["PG_rat"]] <- list(Path = paste0(dirPG, "/", ttlR1),
                                 Plot = proteoCraft::plotEval(ratPlot1),
                                 Ext = "jpeg")
    #suppressMessages({
    #  ggsave(paste0(dirPG, "/", ttlR1, ".jpeg"), ratPlot1, dpi = 300)
    #  ggsave(paste0(dirPG, "/", ttlR1, ".pdf"), ratPlot1, dpi = 300)
    #})
    #ReportCalls <- AddPlot2Report(Plot = ratPlot1, Title = ttlR1)
    #
    # - Peptide intensities
    kl <- grep("\\.REF$", kPpE, value = TRUE, invert = TRUE)
    Samples <- cleanNms(gsub(topattern(pep.ref[length(pep.ref)]), "", kl))
    temp1 <- data.table(pep[, kl])
    colnames(temp1) <- Samples
    temp1 <- melt(temp1, measure.vars = Samples)
    temp2 <- data.table(pep.quant.data.norm[, kl])
    colnames(temp2) <- Samples
    temp2 <- melt(temp2, measure.vars = Samples)
    temp1$Norm <- "Original"
    temp2$Norm <- "Normalised"
    temp <- rbind(temp1, temp2)
    rm(temp1, temp2)
    temp$Sample <- factor(temp$variable, levels = Samples)
    temp$variable <- NULL
    temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
    temp$value <- suppressWarnings(log10(temp$value)) # Peptide intensities are not log-transformed!
    temp <- temp[which(is.all.good(temp$value, 2)),]
    ttlI2 <- paste0("Peptides ", ttl, " from PGs, intensity")
    intPlot2 <- ggplot(temp) +
      geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(Norm~.) + ggtitle(ttlI2) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    #print(intPlot2)
    nrmPlots[["Pep_int"]] <- list(Path = paste0(dirPep, "/", ttlI2),
                                  Plot = proteoCraft::plotEval(intPlot2),
                                  Ext = "jpeg")
    #suppressMessages({
    #  ggsave(paste0(dirPep, "/", ttlI2, ".jpeg"), intPlot2, dpi = 300)
    #  ggsave(paste0(dirPep, "/", ttlI2, ".pdf"), intPlot2, dpi = 300)
    #})
    #ReportCalls <- AddPlot2Report(Plot = intPlot2, Title = ttlI2)
    #
    # - Peptides ratios
    Samples <- cleanNms(gsub(topattern(pep.ratios.ref[length(pep.ratios.ref)]), "", kPpR))
    temp1 <- data.table(pep[, kPpR])
    colnames(temp1) <- Samples
    temp1 <- melt(temp1, measure.vars = Samples)
    temp2 <- data.table(pep.quant.data.norm[, kPpR])
    colnames(temp2) <- Samples
    temp2 <- melt(temp2, measure.vars = Samples)
    temp1$Norm <- "Original"
    temp2$Norm <- "Normalised"
    temp <- rbind(temp1, temp2)
    rm(temp1, temp2)
    temp$Sample <- factor(temp$variable, levels = Samples)
    temp$variable <- NULL
    temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
    temp <- temp[which(is.all.good(temp$value, 2)),]
    ttlR2 <- paste0("Peptides ", ttl, " from PGs, ratios")
    ratPlot2 <- ggplot(temp) +
      geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
      scale_color_viridis_d(begin = 0.25) +
      scale_fill_viridis_d(begin = 0.25) +
      facet_grid(Norm~.) + ggtitle(ttlR2) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    #print(ratPlot2)
    nrmPlots[["Pep_rat"]] <- list(Path = paste0(dirPep, "/", ttlR2),
                                  Plot = proteoCraft::plotEval(ratPlot2),
                                  Ext = "jpeg")
    #suppressMessages({
    #  ggsave(paste0(dirPep, "/", ttlR2, ".jpeg"), ratPlot2, dpi = 300)
    #  ggsave(paste0(dirPep, "/", ttlR2, ".pdf"), ratPlot2, dpi = 300)
    #})
    #ReportCalls <- AddPlot2Report(Plot = ratPlot2, Title = ttlR2)
    #
    nrmPlots2 <- lapply(nrmPlots, function(x) {
      x$Ext <- "pdf"
      return(x)
    })
    nrmPlots <- c(nrmPlots, nrmPlots2)
    # Save plots
    tst <- parLapply(parClust, nrmPlots, function(x) {
      suppressMessages({
        ggplot2::ggsave(paste0(x$Path, ".", x$Ext), x$Plot, dpi = 300)
      })
    })
    #
    if (exists("acceptPGNorm")) {
      acceptPGNorm <- as.logical(acceptPGNorm)
      if ((length(acceptPGNorm) != 1)||(is.na(acceptPGNorm))) { acceptPGNorm <- TRUE }
    } else { acceptPGNorm %<o% TRUE }
    msg <- "Accept re-normalisation results?"
    appNm <- paste0(dtstNm, " - PG re-normalisation")
    ui <- fluidPage(
      useShinyjs(),
      setBackgroundColor( # Doesn't work
        color = c(#"#F8F8FF",
          "#E1F7E7"),
        gradient = "linear",
        direction = "bottom"
      ),
      extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
      titlePanel(tag("u", "Protein Groups (and peptides) re-normalisation"),
                 appNm),
      br(),
      checkboxInput("acceptPGNorm", msg, acceptPGNorm),
      actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
      # fluidRow(column(4, withSpinner(plotOutput("pgIntensities"))),
      #          column(4, withSpinner(plotOutput("pgRatios")))),
      #br(),
      # fluidRow(column(4, withSpinner(plotOutput("pepIntensities"))),
      #          column(4, withSpinner(plotOutput("pepRatios")))),
      fluidRow(column(6, withSpinner(imageOutput("pgIntensities", width = "100%", height = "100%"))),
               column(6, withSpinner(imageOutput("pgRatios", width = "100%", height = "100%")))),
      br(),
      br(),
      br(),
      fluidRow(column(6, withSpinner(imageOutput("pepIntensities", width = "100%", height = "100%"))),
               column(6, withSpinner(imageOutput("pepRatios", width = "100%", height = "100%")))),
      br(),
      br()
    )
    IMGs <- c(paste0(dirPG, "/", c(ttlI1, ttlR1), ".jpeg"),
              paste0(dirPep, "/", c(ttlI2, ttlR2), ".jpeg"))
    IMGsDims <- as.data.frame(t(parSapply(parClust, IMGs, function(x) { #x <- IMGs[1]
      a <- jpeg::readJPEG(x)
      setNames(dim(a)[1:2], c("height", "width"))
    })))
    IMGsDims$height <- screenRes$width*0.35*IMGsDims$height/max(IMGsDims$height)
    IMGsDims$width <- screenRes$width*0.35*IMGsDims$width/max(IMGsDims$width)
    server <- function(input, output, session) {
      # output$pgIntensities <- renderPlot(intPlot1)
      # output$pgRatios <- renderPlot(ratPlot1)
      # output$pepIntensities <- renderPlot(intPlot2)
      # output$pepRatios <- renderPlot(ratPlot2)
      output$pgIntensities <- renderImage({
        list(src = IMGs[1], height = IMGsDims$height[1], width = IMGsDims$width[1])
      }, deleteFile = FALSE)
      output$pgRatios <- renderImage({
        list(src = IMGs[2], height = IMGsDims$height[2], width = IMGsDims$width[2])
      }, deleteFile = FALSE)
      output$pepIntensities <- renderImage({
        list(src = IMGs[3], height = IMGsDims$height[3], width = IMGsDims$width[3])
      }, deleteFile = FALSE)
      output$pepRatios <- renderImage({
        list(src = IMGs[4], height = IMGsDims$height[4], width = IMGsDims$width[4])
      }, deleteFile = FALSE)
      #
      observeEvent(input$acceptPGNorm, { assign("acceptPGNorm", as.logical(input[["acceptPGNorm"]]), envir = .GlobalEnv) })
      observeEvent(input$saveBtn, { stopApp() })
      #observeEvent(input$cancel, { stopApp() })
      session$onSessionEnded(function() { stopApp() })
    }
    eval(parse(text = runApp), envir = .GlobalEnv)
    #
    if (acceptPGNorm) {
      quant.data <- quant.data.norm
      pep[, colnames(pep.quant.data.norm)] <- pep.quant.data.norm
      if (Param$Ratios.Thresholds == threshMsg) {
        pep.Ref.Ratios <- pep.Ref.Ratios.norm
        pep[, colnames(pep.Ref.Ratios)] <- pep.Ref.Ratios
      }
    }
  }
  ReportCalls <- AddSpace2Report()
}
