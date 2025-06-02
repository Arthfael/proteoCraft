# Run peptides-level normalisations
# - Can I rewrite this to handle also proteins level re-normalisation?
lNorm <- length(normSequence) #lNorm <- 1
if (lNorm) {
  # Initial values
  pepNorm %<o% list()
  NormGrps %<o% setNames(aggregate(pep$id, list(pep$"Normalisation group"), list),
                         c("Group", "IDs"))
  nrmDr <- paste0(wd, "/Workflow control/Peptides/Intensities")
  #
  # Preliminary step: log10 transformation
  #i <- grep("^Norm[0-9]+", names(pep.ref))
  #if (length(i)) { i <- min(i) - 1 } else { i <- length(pep.ref) }
  strt <- c("Original", "Imputed")
  strt <- rev(strt[which(strt %in% names(pep.ref))])[1]
  rf <- pep.ref[strt]
  kol <- grep(topattern(rf), colnames(pep), value = TRUE)
  tmpDat1 <- do.call(cbind, lapply(kol, function(k) { log10(pep[[k]]) }))
  tmpDat1 <- as.data.frame(tmpDat1)
  colnames(tmpDat1) <- allSamples <- gsub(topattern(rf), "", kol)
  addKol <- setNames(c("id", "Group", "Modified sequence", "Proteins"),
                     c("id", "Normalisation group", "Modified sequence", "Proteins"))
  tmpDat1[, addKol] <- pep[, names(addKol)]
  wAG1 <- which(parApply(parClust, tmpDat1[, allSamples], 1, function(x) { length(proteoCraft::is.all.good(x)) }) > 0)
  pepNorm[[1]] <- list(Data = tmpDat1,
                       Filter = wAG1,
                       Pass = TRUE)
  #
  tmp <- Data_Impute2(tmpDat1[wAG1, allSamples],
                      Exp.map[match(allSamples, Exp.map$Ref.Sample.Aggregate), VPAL$column])
  tmp <- as.matrix(tmp$Imputed_data)
  tst <- pcaBatchPlots(tmp,
                       strt,
                       VPAL$column,
                       Exp.map,
                       "Ref.Sample.Aggregate",
                       "",
                       nrmDr,
                       "PCA plot - before norm.")
  #
  for (nrmStp in 1:lNorm) { #nrmStp <- 1 #nrmStp <- nrmStp+1
    cat(" ->", nrmStp, normSequence[[nrmStp]]$Method, "\n")
    rg <- 1:nrmStp # (and not 1:(nrmStp-1): the first in the list is pre-norm data -> there is an offset)
    prevStp <- max(which(sapply(rg, function(i) { pepNorm[[i]]$Pass })))
    tmpDat1 <- pepNorm[[prevStp]]$Data # The first in the list is pre-norm data -> offset
    #View(tmpDat1)
    #normSequence[[nrmStp]]$Method
    w <- which(!addKol %in% colnames(tmpDat1))
    if (length(w)) { tmpDat1[, addKol[w]] <- pep[, names(addKol)[w]] }
    wAG1 <- pepNorm[[prevStp]]$Filter
    nrmSrc <- paste0(libPath, "/extdata/R scripts/Sources/", normSequence[[nrmStp]]$Source)
    source(nrmSrc, local = FALSE)
    #rstudioapi::documentOpen(nrmSrc)
    # These shources should all take as inputs:
    #   - tmpDat1 -> data from previous step
    #   - wAG1 -> filter (though they are all currently the same values... as should be expected)
    # ... and create the following outputs:
    #   - tmpDat2 -> normalized output
    #   - Normalisation -> details of the normalisation applied
    #   - wAG2 -> filter
    #   - txt2 -> text for materials and methods
    #   - Outcome -> did it work, or was it accepted by the user?
    #View(tmpDat2)
    pepNorm[[nrmStp+1]] <- list(Data = tmpDat2,
                                Normalisation = normSequence[[nrmStp]],
                                Filter = wAG2,
                                Text = txt2,
                                Pass = Outcome)
  }
  wNorm <- which(sapply(pepNorm, function(x) { x$Pass }))
  wNorm <- wNorm[which(wNorm != 1)]
  if (length(wNorm)) {
    # Visualisations
    dat <- lapply(c(1, wNorm), function(i) {
      x <- as.data.frame(pepNorm[[i]]$Data[, allSamples])
      colnames(x) <- proteoCraft::cleanNms(colnames(x))
      x[, addKol] <- pep[, names(addKol)]
      x <- x[pepNorm[[i]]$Filter,]
      x <- melt(x, id.vars = addKol)
      colnames(x) <- c(addKol, "Sample", "value")
      tmp <- proteoCraft::cleanNms(Exp.map[[VPAL$column]])
      m <- match(x$Sample, proteoCraft::cleanNms(Exp.map$Ref.Sample.Aggregate))
      x$"Samples group" <- tmp[m]
      x$Replicate <- Exp.map$Replicate[m]
      if (i == 1) { tmp <- "Original" } else { tmp <- normSequence[[i-1]]$Method }
      x$Norm <- paste0(i, " - ", tmp)
      return(x)
    })
    dat <- do.call(rbind, dat)
    unique(dat$Norm)
    dat$Norm <- factor(dat$Norm, levels = c("1 - Original",
                                            paste0(wNorm, " - ", sapply(normSequence[wNorm-1], function(x) { x$Method }))))
    ttl <- "Peptides intensity normalisation"
    kolz <- "."
    if (length(unique(pep$Normalisation_group)) > 1) { kolz <- c(kolz, "Normalisation_group") }
    form <- as.formula(paste0("Norm~", paste0(kolz, collapse = "+")))
    plot <- ggplot(dat) +
      geom_violin(aes(x = Sample, y = value, color = `Samples group`, fill = `Samples group`), alpha = 0.25) +
      geom_boxplot(aes(x = Sample, y = value, color = `Samples group`, fill = `Samples group`), alpha = 0.5) +
      scale_color_viridis(begin = 0.25, discrete = TRUE, option = "C") +
      scale_fill_viridis(begin = 0.25, discrete = TRUE, option = "C") +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      geom_vline(xintercept = 0) + ggtitle(ttl) + facet_grid(form) +
      theme(strip.text.y = element_text(angle = 0)) +
      coord_fixed(0.5)
    poplot(plot, 12, 22)
    ggsave(paste0(nrmDr, "/", ttl, ".jpeg"), plot, dpi = 300, height = 10, units = "in")
    ggsave(paste0(nrmDr, "/", ttl, ".pdf"), plot, dpi = 300, height = 10, units = "in")
    ReportCalls <- AddPlot2Report(Dir = nrmDr)
    #
    finNorm <- max(wNorm)
    newDat <- pepNorm[[finNorm]]$Data[, allSamples]
    tmp <- as.data.frame(newDat)
    tmp <- proteoCraft::Data_Impute2(tmp,
                                     Exp.map[match(allSamples, Exp.map$Ref.Sample.Aggregate), VPAL$column])
    tmp <- as.matrix(tmp$Imputed_data)
    tst <- pcaBatchPlots(tmp,
                         "Normalisation",
                         VPAL$column,
                         Exp.map,
                         "Ref.Sample.Aggregate",
                         "",
                         nrmDr,
                         "PCA plot - final norm.")
    # De-log
    newDatLin <- newDat
    for (smpl in RSA$values) {
      newDatLin[[smpl]] <- 10^newDatLin[[smpl]]
    }
    # Assign results to pep
    pep.ref["Normalisation"] <- paste0("norm. ", pep.ref["Original"])
    pep[, paste0(pep.ref["Normalisation"], RSA$values)] <- newDatLin[, RSA$values]
    saveFun(pepNorm, paste0(nrmDr, "/pep_intens_norm.RData"))
    #
    # MatMet
    TxtSteps <- unlist(lapply(pepNorm[wNorm], function(x) { x$Text }))
    TxtSteps <- TxtSteps[which(nchar(TxtSteps) > 0)]
    l <- length(TxtSteps)
    for (wrd in c("normalized", "corrected")) {
      pat <- paste0("^", wrd, " ")
      g <- grep(pat, TxtSteps)
      g <- g[which(g %in% (g+1))]
      if (length(g)) { TxtSteps[g] <- gsub(pat, "", TxtSteps[g]) }
    }
    TxtSteps <- paste0(paste(TxtSteps[1:(l-1)], collapse = ", "), ", then ", TxtSteps[l])
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " Peptide intensities were ", TxtSteps, ". Peptidoform-level ratios were then calculated.")
    #
  }
}
