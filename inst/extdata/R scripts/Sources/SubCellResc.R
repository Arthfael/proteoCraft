# For Sub-cellular localisation analysis
# (or any other analysis where we want to apply some additional a priori known normalisation factor):
# Re-scale to old sub-cellular fraction relative levels, as measured:
#  - Sample specific evidence intensity range prior to evidence normalisation,
#  - The "Proportion" columns in Exp.map, which is the ratio of loaded to total available sample amount.
# Because we would still want to have normalized data, if either of the following arguments is TRUE,
# the re-scaling is done within the VPAL groups:
# - "Norma.Ev.Intens"
# - "Norma.Pep.Intens"
# - "Norma.Pep.Ratio"
# - "Norma.Prot.Intens"
# - "Norma.Prot.Ratio"
if (LocAnalysis) {
  stop("Check this code before running it!")
  if (!"Proportion" %in% colnames(Exp.map)) {
    warning("The \"Proportion\" column is absent from the Experiment map, we will assume that the same proportion of each fraction was processed and analyzed.")
    Exp.map$Proportion <- 1 # In case the column was omitted
  } else { Exp.map$Proportion <- as.numeric(Exp.map$Proportion) }
  #
  Prot.Expr.Root2 %<o% setNames(gsub("Expr\\.", "resc. Expr.", Prot.Expr.Root), "SubCell. Profile")
  pep.ref2 %<o% setNames(paste0("resc. ", pep.ref[1L]), "SubCell. Profile")
  # Best way to do this:
  # - Compare final values at peptidoforms level
  # - Get original, uncorrected evidences values and aggregate them into peptidoforms level values
  # - Apply the corresponded ratios to protein groups and peptidoforms-level data to recreate original scale
  #
  # Also, take into account the proportion of each fraction which was actually processed
  prevRef <- pep.ref[length(pep.ref)]
  kol1 <- paste0(prevRef, RSA$values)
  WhInColNms <- which(kol1 %in% colnames(pep))
  kol1 <- kol1[WhInColNms]
  tmp1 <- pep[, c("Modified sequence", kol1)]
  kol2 <- kol2a <- c("Modified sequence", "Raw file path", "MQ.Exp", "Experiment", ev.col["Original"])
  if (LabelType == "Isobaric") {
    kol2b <- grep(paste0(topattern(ev.ref["Original"]), "[0-9]+$"), colnames(ev), value = TRUE)
    chan <- gsub(topattern(ev.ref["Original"]), "", kol2b)
    kol2 <- c(kol2a, kol2b)
    tmp2 <- ev[, kol2]
    tst <- apply(tmp2[, kol2b], 1L, sum)
    tmp2[, kol2b] <- sweep(tmp2[, kol2b], 1L, tmp2[[ev.col["Original"]]]/tst, "*")
    tmp2a <- aggregate(tmp2[, kol2b], list(tmp2$"Modified sequence", tmp2$MQ.Exp), sum)
    tmp2 <- data.frame(`Modified sequence` = pep$"Modified sequence", check.names = FALSE)
    for (mqexp in MQ.Exp) { #mqexp <- MQ.Exp[1L]
      em <- Exp.map[which(Exp.map$MQ.Exp == mqexp),]
      m <- match(chan, em$"Isobaric label")
      w <- which(!is.na(m))
      kol2c <- paste0(pep.ref["Original"], em$Ref.Sample.Aggregate[m[w]])
      tmp2[, kol2c] <- 0
      wa <- which(tmp2a$Group.2 == mqexp)
      wb <- which(tmp2$"Modified sequence" %in% tmp2a$Group.1[wa])
      tmp2[wb, kol2c] <- tmp2a[wa[match(tmp2$"Modified sequence"[wb], tmp2a$Group.1[wa])], kol2b[w]]
    }
  } else {
    tmp2a <- as.data.table(ev[, kol2])
    colnames(tmp2a)[which(colnames(tmp2a) == ev.col["Original"])] <- "Int"
    tmp2a <- tmp2a[, list(x = sum(Int)),
                   by = list(`Modified sequence` = `Modified sequence`,
                             `Raw file` = `Raw file path`,
                             MQ.Exp = MQ.Exp,
                             Experiment = Experiment)]
    tmp3 <- listMelt(Exp.map$MQ.Exp, Exp.map$Ref.Sample.Aggregate)
    tmp2a$RSA <- tmp3$L1[match(tmp2a$MQ.Exp, tmp3$value)]
    tmp2 <- data.frame(`Modified sequence` = pep$"Modified sequence", check.names = FALSE)
    for (rsa in tmp3$L1) { #rsa <- tmp3$L1[1L]
      kol2c <- paste0(pep.ref["Original"], rsa)
      tmp2[[kol2c]] <- 0
      w <- which(tmp2a$RSA == rsa)
      tmp2b <- copy(tmp2a)
      tmp2b <- tmp2b[w, list(x = sum(x)), by = list(`Modified sequence`)]
      tmp2b <- as.data.frame(tmp2b)
      w <- which(tmp2$"Modified sequence" %in% tmp2b$`Modified sequence`)
      tmp2[w, kol2c] <- tmp2b$x[match(tmp2$"Modified sequence"[w], tmp2b$`Modified sequence`)]
    }
  }
  # re-order
  kol2c <- paste0(pep.ref["Original"], RSA$values)
  kol2c <- kol2c[WhInColNms]
  tmp2 <- tmp2[, c("Modified sequence", kol2c)]
  # Calculate ratio before/after, per samples group
  BefAft <- tmp2[, kol2c]/tmp1[, kol1]
  colnames(BefAft) <- RSA$values[WhInColNms]
  BefAft <- vapply(VPAL$values, \(x) {
    x <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)]
    x <- x[which(x %in% RSA$values[WhInColNms])]
    x <- median(unlist(BefAft[, x]), na.rm = TRUE)
    return(x)
  }, 1)
  # This is the ratio between values at the start of the workflow and final peptidoforms values 
  #
  # Next, we can apply proportion of total material loaded
  # Now...
  # At this stage we should have nicely normalized, i.e., "aligned", data
  # We should thus work the same way within groups of replicates.
  # Thus, if there are different loaded amounts per individual sample,
  # when we apply corrections aimed at restoring original samples' relative intensity scales,
  # we should average at this stage over replicates of the same condition.
  Props <- Exp.map[match(RSA$values[WhInColNms], Exp.map$Ref.Sample.Aggregate),
                   c("Proportion", VPAL$column)]
  Props <- aggregate(Props$Proportion, list(Props[[VPAL$column]]), mean)
  Props <- Props$x[match(VPAL$values, Props$Group.1)]
  BefAft <- BefAft/Props
  # Note: this is done currently at sample group level. It may make sense to do it at subcellular fraction level in the future...
  # but there are also risks, e.g. when the perturbation studied changes cell morphology dramatically.
  # => make it an option?
  #
  # Apply correction:
  # We want to calculate the average change within sample groups (compartments/fractions x treatment group)
  # so as to still profit from the normalisation did on peptides/proteins.
  # (Presumably it is ok to normalize within those groups, as samples should be replicates)
  # For now the covariates aggregate used is VPAL, but we could map it to a custom one using a new Param
  for (grp in VPAL$values) { #grp <- VPAL$values[1L]
    smpls <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == grp)]
    PepKol1 <- paste0(prevRef, smpls)
    PepKol2 <- paste0(pep.ref2, smpls) 
    w <- which(PepKol1 %in% colnames(pep))
    pep[, PepKol2[w]] <- pep[, PepKol1[w]]*BefAft[grp]
    PGKol1 <- paste0(Prot.Expr.Root, smpls)
    PGKol2 <- paste0(Prot.Expr.Root2, smpls) 
    w <- which(PGKol1 %in% colnames(quantData))
    quantData[, PGKol2[w]] <- quantData[, PGKol1[w]] + log10(BefAft[grp])
    # NB: Ratios are not re-scaled...
    # and anyway, the ratios should be done within groups which should correspond to the correct proportions 
  }
  l <- length(DatAnalysisTxt)
  DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                              " Expression values were re-scaled per fraction to reflect total original protein amount.")
  # Visualize
  dir <- paste0(wd, "/Workflow control/Re-scaling/")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  for (i in 1L:2L) {
    if (i == 1L) {
      rt1 <- prevRef
      rt2 <- pep.ref2
      ttl <- "Peptide intensities - effect of post-normalisation rescaling"
      temp <- pep
    }
    if (i == 2L) {
      rt1 <- Prot.Expr.Root
      rt2 <- Prot.Expr.Root2
      ttl <- "PGs expression values - effect of post-normalisation rescaling"
      temp <- quantData
    }
    kol1 <- paste0(rt1, RSA$values)
    kol2 <- paste0(rt2, RSA$values)
    w <- which((kol1 %in% colnames(temp))&(kol2 %in% colnames(temp)))
    tst <- temp[, c(kol1[w], kol2[w])]
    tst <- reshape2::melt(tst, measure.vars = c(kol1[w], kol2[w]))
    tst$value <- suppressWarnings(log10(tst$value))
    tst <- tst[which(is.all.good(tst$value, 2L)),]
    tst$variable <- as.character(tst$variable)
    tst2 <- data.frame(variable = unique(tst$variable))
    tst2[, c("Type", "Sample")] <- as.data.frame(t(sapply(strsplit(tst2$variable, " - "), unlist)))
    tst2[, RSA$names] <- as.data.frame(t(sapply(strsplit(tst2$Sample, "___"), unlist)))
    tst2$Group <- Exp.map[match(tst2$Sample, Exp.map$Ref.Sample.Aggregate), Volcano.plots.Aggregate.Level$aggregate]
    tst2$Sample <- cleanNms(tst2$Sample)
    tst2$Group <- cleanNms(tst2$Group)
    tst2$Type <- factor(tst2$Type, levels = gsub(" - $", "", c(rt1, rt2)))
    tst[, colnames(tst2)] <- tst2[match(tst$variable, tst2$variable), colnames(tst2)]
    plot <- ggplot(tst) +
      geom_violin(aes(x = Sample, y = value, colour = Group, fill = Group), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, colour = Group, fill = Group), alpha = 0.5) +
      scale_color_viridis(begin = 0.25, discrete = TRUE, option = "C") +
      scale_fill_viridis(begin = 0.25, discrete = TRUE, option = "C") +
      facet_grid(Type~.) + theme_bw() + ggtitle(ttl) +
      theme(axis.text.x = element_text(angle = 90L, vjust = 0.5, hjust = 1))
    if (i == 2L) {
      print(plot) # This type of QC plot does not need to pop up, the side panel is fine
    }
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    })
    ReportCalls <- AddPlot2Report()
  }
}
