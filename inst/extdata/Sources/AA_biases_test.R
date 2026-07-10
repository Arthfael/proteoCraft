# Test for amino acid biases:
AA_biases %<o% AA_bias(Ev = ev, DB = db)
#View(AA_biases)
write.csv(AA_biases, paste0(wd, "/Workflow control/AA_biases.csv"), row.names = FALSE)
qcDir <- paste0(wd, "/Summary plots")
if (!dir.exists(qcDir)) { dir.create(qcDir, recursive = TRUE) }
if (scrptType == "withReps") { dirlist <- union(dirlist, qcDir) }
ttl <- "Amino acid observational biases"
plot <- ggplot(AA_biases) +
  geom_col(aes(x = AA, y = log2(Ratio), fill = AA)) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  xlab("Amino acid") + ylab("log2 ratio(freq. obs. dataset / freq. parent proteome)") +
  scale_y_continuous(expand = c(0L, 0L)) +
  scale_fill_viridis_d(begin = 0.25) +
  ggtitle(ttl, subtitle = "(observed dataset VS parent proteome)") + theme_bw() +
  theme(legend.position = "none")
print(plot)
suppressMessages({
  ggsave(paste0(qcDir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
  ggsave(paste0(qcDir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
})
plotLy <- ggplotly(plot, tooltip = c("x", "y"))
if (!exists("QC_plotLys")) { QC_plotLys %<o% list() }
setwd(qcDir)
saveWidget(plotLy, paste0(qcDir, "/", ttl, ".html"), selfcontained = TRUE)
setwd(wd)
QC_plotLys[[ttl]] <- plotLy
