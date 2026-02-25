# Calculate peptide ratios
# We are simplifying by getting rid of those
makePepRat %<o% FALSE
if (makePepRat) {
  pep.ratios.root %<o% "Ratio"
  pep.ratios.ref %<o% c(Ratios = paste0("log2(", pep.ratios.root, ") - "))
  tmpRt <- make_Rat2(pep,
                     myContrasts,
                     RRG,
                     FALSE,
                     2,
                     Exp.map,
                     int.root = pep.ref[length(pep.ref)],
                     rat.root = pep.ratios.ref)
  pep[ colnames(tmpRt)] <- tmpRt
  #
  # Visualize:
  kol <- paste0(pep.ratios.ref[1], RSA$values)
  kol <- intersect(kol, colnames(pep))
  temp <- pep[, kol]
  colnames(temp) <- sub(topattern(pep.ratios.ref[1]), "", colnames(temp))
  temp <- reshape::melt(temp, measure.vars = colnames(temp))
  temp <- temp[which(is.all.good(temp$value, 2)),]
  temp[, RSA$names] <- Isapply(strsplit(as.character(temp$variable), "___"), unlist)
  temp$variable <- as.factor(cleanNms(as.character(temp$variable)))
  aggr <- RSA$names
  tst <- vapply(aggr, function(a) { length(get(substr(a, 1, 3))) }, 1)
  aggr <- aggr[which(tst > 1)]
  facets <- (length(aggr) > 0)
  if (facets) {
    if (length(aggr) > 1) { form <- paste0(aggr[1], "~", paste(aggr[2:length(aggr)], collapse = "+")) } else {
      form <- paste0(".~", aggr[1])
    }
    form <- as.formula(form)
  }
  dir <- paste0(wd, "/Workflow control/Peptides/Ratios")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  ttl <- "Peptide ratios distribution"
  plot <- ggplot(temp) +
    geom_histogram(aes(x = value, fill = variable), bins = ceiling((max(temp$value)-min(temp$value))/0.05)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_d(begin = 0.25) +
    ggtitle(ttl) + theme_bw() + theme(legend.position = "none")
  if (facets) { plot <- plot + facet_grid(form) }
  #poplot(plot, 12, 22)
  print(plot)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
  })
  ReportCalls <- AddPlot2Report()
  #
  #### Code chunk - Normalize peptide ratios
  # Legacy, I really do not think that this is a good idea at this stage if intensities have been well normalized.
  # Even if we are dealing with SILAC (not supported yet) where ratios are better measured and should be more precise,
  # I would simply re-calculate updated intensities to reflect the ratios at an early stage, then focus on intensities from there.
  # Currently breaks immediately: before running this, the old code should be checked and the parameters shiny app updated to support the necessary parameters.
  if (Param$Norma.Pep.Ratio) {
    Src <- paste0(libPath, "/extdata/R scripts/Sources/rep_pepRat_ReNorm.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
}
