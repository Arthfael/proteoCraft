# Test for amino acid biases:
AA_biases %<o% AA_bias(Ev = ev, DB = db)
#View(AA_biases)
write.csv(AA_biases, paste0(wd, "/Workflow control/AA_biases.csv"), row.names = FALSE)
dir <- paste0(wd, "/Summary plots")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
ttl <- "Amino acid observational biases"
plot <- ggplot(AA_biases) +
  geom_col(aes(x = AA, y = Ratio, fill = AA)) +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed") +
  xlab("Amino acid") + ylab("Ratio of frequency in observed dataset / parent proteome") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d(begin = 0.25) +
  ggtitle(ttl, subtitle = "(observed dataset VS parent proteome)") + theme_bw() +
  theme(legend.position = "none")
print(plot)
ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\"Amino acids frequency biases:\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_table(Report, AA_biases)")
ReportCalls <- AddPlot2Report()
ReportCalls <- AddSpace2Report()
