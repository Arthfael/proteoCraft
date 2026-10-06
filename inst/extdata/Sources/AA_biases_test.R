# Test for amino acid biases:
AA_biases %<o% AA_bias(Ev = ev, DB = db)
#View(AA_biases)
write.csv(AA_biases, paste0(wd, "/Workflow control/AA_biases.csv"), row.names = FALSE)
qcDir <- paste0(wd, "/Summary plots")
if (!dir.exists(qcDir)) { dir.create(qcDir, recursive = TRUE) }
if (scrptType == "withReps") { dirlist <- union(dirlist, qcDir) }
ttl <- "Amino acid observational biases"
plot <- ggplot2::ggplot(AA_biases) +
  ggplot2::geom_col(ggplot2::aes(x = AA, y = log2(Ratio), fill = AA)) +
  ggplot2::geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  ggplot2::xlab("Amino acid") + ggplot2::ylab("log2 ratio(freq. obs. dataset / freq. parent proteome)") +
  ggplot2::scale_y_continuous(expand = c(0L, 0L)) +
  ggplot2::scale_fill_viridis_d(begin = 0.25) +
  ggplot2::ggtitle(ttl, subtitle = "(observed dataset VS parent proteome)") + ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none")
print(plot)
suppressMessages({
  ggplot2::ggsave(paste0(qcDir, "/", ttl, ".svg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
})
plotLy <- plotly::ggplotly(plot, tooltip = c("x", "y"))
plotLy$x$layout$xaxis$autorange <- TRUE
plotLy$x$layout$yaxis$autorange <- TRUE
plotLy <- htmlwidgets::onRender(plotLy, global_autorange)
plotLy <- plotly::config(plotLy,
                         modeBarButtonsToRemove = c("select2d", "lasso2d"))
plotLy <- plotly::plotly_build(plotLy)
#plotLy <- plotly::partial_bundle(plotLy)
setwd(qcDir)
htmlwidgets::saveWidget(plotly::partial_bundle(plotLy), paste0(qcDir, "/", ttl, ".html"), selfcontained = TRUE)
#htmlwidgets::saveWidget(plotLy, paste0(qcDir, "/", ttl, ".html"), selfcontained = TRUE)
setwd(wd)
if ((!exists("QC_plotLys")) && file.exists(qcBckUpFl)) { loadFun(qcBckUpFl) }
if (!exists("QC_plotLys")) { QC_plotLys <- list() }
QC_plotLys[[ttl]] <- plotLy
saveFun(QC_plotLys, qcBckUpFl)
