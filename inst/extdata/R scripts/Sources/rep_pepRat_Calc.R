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
                     2L,
                     Exp.map,
                     int.root = pep.ref[length(pep.ref)],
                     rat.root = pep.ratios.ref)
  kol <- colnames(tmpRt)
  pep[, kol] <- tmpRt
  #
  # Visualize:
  temp <- pep[, kol]
  colnames(temp) <- sub(topattern(pep.ratios.ref[1L]), "", colnames(temp))
  temp <- reshape::melt(temp, measure.vars = colnames(temp))
  colnames(temp) <- c("contrast", "ratio")
  temp$contrast <- factor(sub("\\) - \\(", ")\n- (", temp$contrast),
                          levels = sub("\\) - \\(", ")\n- (", myContrasts$Contrast))
  temp <- temp[which(is.all.good(temp$ratio, 2L)),]
  dir <- paste0(wd, "/Workflow control/Peptides/Ratios")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  ttl <- "Peptide ratios distribution"
  plot <- ggplot(temp) +
    geom_histogram(aes(x = ratio, fill = contrast), bins = ceiling((max(temp$ratio)-min(temp$ratio))/0.05)) +
    scale_y_continuous(expand = c(0L, 0L)) +
    scale_fill_viridis_d(begin = 0.25) +
    ggtitle(ttl) + theme_bw() + theme(legend.position = "none") + facet_wrap(~contrast)
  #poplot(plot, 12L, 22L)
  print(plot)
  suppressMessages({
    ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
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
