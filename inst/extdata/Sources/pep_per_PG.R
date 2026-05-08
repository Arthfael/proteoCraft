# Peptidoforms per Protein Group
tmp <- aggregate(PG$id, list(PG$`Peptides count`), length) # Faster than data.table here
colnames(tmp) <- c("Peptides count", "Protein groups")
pal <- colorRampPalette(c("brown", "yellow"))(max(tmp$"Peptides count")-1L)
tmp$Colour <- c("blue", pal)[tmp$`Peptides count`]
tmp2 <- summary(PG$`Peptides count`)
tmp2 <- data.frame(Variable = c(names(tmp2), "", "Protein groups", "Protein groups with 2+ peptidoforms"),
                   Value = c(as.character(signif(as.numeric(tmp2), 3L)),
                             "",
                             as.character(c(nrow(PG), sum(PG$"Peptides count" >= 2L)))))
tmp2$Txt <- apply(tmp2[, c("Variable", "Value")], 1L, \(x) {
  x <- x[which(x != "")]
  x <- if (length(x)) { paste(x, collapse = ": ") } else { "" }
  return(x)
})
tmp2$X <- max(as.numeric(tmp2$Value[match("Max.", tmp2$Variable)]))*0.98
tmp2$Y <- max(tmp$"Protein groups")*(0.98-(0L:(nrow(tmp2) - 1L))*0.02)
ttl <- "Peptidoforms per PG"
plot <- ggplot(tmp) +
  geom_col(aes(x = `Peptides count`, y = `Protein groups`, fill = Colour),
           colour = NA) +
  theme_bw() + ggtitle(ttl) +
  geom_text(data = tmp2, aes(x = X, y = Y, label = Txt), hjust = 1, size = 3L) +
  scale_fill_identity() +
  scale_x_continuous(breaks = seq(10L, floor(max(tmp$`Peptides count`)/10L)*10L, by = 10L), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
#  coord_trans(x = "log10", y = "log10")
poplot(plot, 12L, 22L) # This type of QC plot does not need to pop up, the side panel is fine
dir <- paste0(wd, "/Summary plots")
if (scrptType == "withReps") { dirlist<- unique(c(dirlist, dir)) }
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
suppressMessages({
  ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300L)
  ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L)
})
