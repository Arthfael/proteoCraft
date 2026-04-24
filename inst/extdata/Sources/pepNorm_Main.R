# Run peptides-level normalisations
# - Can I rewrite this to handle also proteins level re-normalisation?
lNorm <- length(normSequence) #lNorm <- 1L
if (lNorm) {
  # Initial values
  pepNorm %<o% list()
  NormGrps %<o% setNames(aggregate(pep$id, list(pep$"Normalisation group"), list),
                         c("Group", "IDs"))
  nrmDr <- paste0(wd, "/Workflow control/Peptides/Intensities")
  #
  # Preliminary step: log10 transformation
  #i <- grep("^Norm[0-9]+", names(pep.ref))
  #i <- if (length(i)) { min(i) - 1L } else { length(pep.ref) }
  strt <- c("Original", "Imputation")[Impute+1L]
  rf <- pep.ref[strt]
  allSamples <- Exp.map$Ref.Sample.Aggregate
  kol <- paste0(rf, allSamples)
  w <- which(kol %in% colnames(pep))
  Exp.map <- Exp.map[w,]
  allSamples <- allSamples[w]
  kol <- kol[w]
  tmpDat1 <- do.call(cbind, lapply(kol, \(k) {
    val <- log10(pep[[k]])
    w <- which(!is.all.good(val, 2L))
    val[w] <- NA # We shouldn't have to deal with other types of invalid values! They break the rest.
    return(val)
  }))
  tmpDat1 <- as.data.frame(tmpDat1)
  colnames(tmpDat1) <- allSamples <- gsub(topattern(rf), "", kol)
  addKol <- setNames(c("id", "Group", "Modified sequence", "Proteins"),
                     c("id", "Normalisation group", "Modified sequence", "Proteins"))
  tmpDat1[, addKol] <- pep[, names(addKol)]
  clusterExport(parClust, "is.all.good", envir = environment())
  wAG1 <- which(parApply(parClust, tmpDat1[, allSamples], 1L, \(x) { length(is.all.good(x)) }) > 0L)
  pepNorm[[1L]] <- list(Data = tmpDat1,
                       Filter = wAG1,
                       Pass = TRUE)
  #View(tmpDat1[wAG1, allSamples])
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
  for (nrmStp in 1L:lNorm) { #nrmStp <- 1L #nrmStp <- nrmStp+1L
    rg <- 1L:nrmStp # (and not 1L:(nrmStp-1): the first in the list is pre-norm data -> there is an offset)
    prevStp <- max(which(vapply(rg, \(i) { pepNorm[[i]]$Pass }, TRUE)))
    cat("\n +++ Step", nrmStp, normSequence[[nrmStp]]$Method, "\n    input = step", prevStp-1L, "\n")
    tmpDat1 <- pepNorm[[prevStp]]$Data # The first in the list is pre-norm data -> offset
    #View(tmpDat1)
    #normSequence[[nrmStp]]$Method
    w <- which(!addKol %in% colnames(tmpDat1))
    if (length(w)) { tmpDat1[, addKol[w]] <- pep[, names(addKol)[w]] }
    wAG1 <- pepNorm[[prevStp]]$Filter
    nrmSrc <- paste0(libPath, "/extdata/Sources/", normSequence[[nrmStp]]$Source)
    #rstudioapi::documentOpen(nrmSrc)
    source(nrmSrc, local = FALSE)
    # These sources should all take as inputs:
    #   - tmpDat1 -> data from previous step
    #   - wAG1 -> filter (though they are all currently the same values... as should be expected)
    # ... and create the following outputs:
    #   - tmpDat2 -> normalized output
    #   - Normalisation -> details of the normalisation applied
    #   - wAG2 -> filter
    #   - txt2 -> text for materials and methods
    #   - Outcome -> did it work, or was it accepted by the user?
    #View(tmpDat2)
    pepNorm[[nrmStp+1L]] <- list(Data = tmpDat2,
                                Normalisation = normSequence[[nrmStp]],
                                Filter = wAG2,
                                Text = txt2,
                                Pass = Outcome)
  }
  wNorm <- which(vapply(pepNorm, \(x) { x$Pass }, TRUE))
  wNorm <- wNorm[which(wNorm != 1L)]
  if (length(wNorm)) {
    # Visualisations
    dat <- lapply(c(1L, wNorm), \(i) {
      df <- as.data.frame(pepNorm[[i]]$Data)
      currSamples <- allSamples[which(allSamples %in% colnames(df))]
      x <- df[, currSamples]
      colnames(x) <- cleanNms(colnames(x))
      x[, addKol] <- pep[, names(addKol)]
      x <- x[pepNorm[[i]]$Filter,]
      x <- melt(x, id.vars = addKol)
      colnames(x) <- c(addKol, "Sample", "value")
      tmp <- cleanNms(Exp.map[[VPAL$column]])
      m <- match(x$Sample, cleanNms(Exp.map$Ref.Sample.Aggregate))
      x$"Samples group" <- tmp[m]
      x$Replicate <- Exp.map$Replicate[m]
      tmp <- if (i == 1L) { "Original" } else { normSequence[[i-1L]]$Method }
      x$Norm <- paste0(i, " - ", tmp)
      return(x)
    })
    dat <- do.call(rbind, dat)
    unique(dat$Norm)
    dat$Norm <- factor(dat$Norm, levels = c("1 - Original",
                                            paste0(wNorm, " - ", vapply(normSequence[wNorm-1L], \(x) { x$Method }, ""))))
    ttl <- "Peptides intensity normalisation"
    kolz <- "."
    if (length(unique(pep$Normalisation_group)) > 1L) { kolz <- c(kolz, "Normalisation_group") }
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
    poplot(plot, 12L, 22L)
    ggsave(paste0(nrmDr, "/", ttl, ".jpeg"), plot, dpi = 150L, height = 10L, units = "in")
    ggsave(paste0(nrmDr, "/", ttl, ".pdf"), plot, dpi = 150L, height = 10L, units = "in")
    ReportCalls <- AddPlot2Report(Dir = nrmDr)
    #
    finNorm <- max(wNorm)
    newDat <- as.data.frame(pepNorm[[finNorm]]$Data)
    currSamples <- allSamples[which(allSamples %in% colnames(newDat))]
    newDat <- newDat[, currSamples]
    tmp <- newDat
    tmp <- Data_Impute2(tmp,
                        Exp.map[match(currSamples, Exp.map$Ref.Sample.Aggregate), VPAL$column])
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
    wHere <- which(RSA$values %in% colnames(newDatLin))
    Ref.Sample.Aggregate$values <- RSA$values <- RSA$values[wHere]
    Exp.map <- Exp.map[which(Exp.map$Ref.Sample.Aggregate %in% RSA$values),]
    Volcano.plots.Aggregate.Level$values <- VPAL$values <- unique(Exp.map[[VPAL$column]])
    Ratios.Ref.Groups$values <- RRG$values <- unique(Exp.map[[RRG$column]])
    Ratios.Groups$values <- RG$values <- unique(Exp.map[[RG$column]])
    for (smpl in RSA$values) {
      newDatLin[[smpl]] <- 10^newDatLin[[smpl]]
    }
    # Assign results to pep
    pep.ref["Normalisation"] <- paste0("norm. ", pep.ref["Original"])
    pep[, paste0(pep.ref["Normalisation"], RSA$values)] <- newDatLin[, RSA$values]
    saveFun(pepNorm, paste0(nrmDr, "/pep_intens_norm.RData"))
    #
    # MatMet
    TxtSteps <- unlist(lapply(pepNorm[wNorm], \(x) { x$Text }))
    TxtSteps <- TxtSteps[which(nchar(TxtSteps) > 0L)]
    l <- length(TxtSteps)
    for (wrd in c("normalized", "corrected")) {
      pat <- paste0("^", wrd, " ")
      g <- grep(pat, TxtSteps)
      g <- g[which(g %in% (g+1L))]
      if (length(g)) { TxtSteps[g] <- gsub(pat, "", TxtSteps[g]) }
    }
    TxtSteps <- paste0(paste(TxtSteps[1L:(l-1L)], collapse = ", "), ", then ", TxtSteps[l])
    l <- length(DatAnalysisTxt)
    DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l], " Peptide intensities were ", TxtSteps,
                                ". Peptidoform-level ratios were then calculated.")
    #
  }
}
