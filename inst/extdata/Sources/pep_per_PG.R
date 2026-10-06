# Peptidoforms per Protein Group
tmp <- aggregate(PG$id, list(PG$`Peptides count`), length) # Faster than data.table here
colnames(tmp) <- c("Peptides count", "Protein groups")
pal <- colorRampPalette(c("brown", "yellow"))(max(tmp$"Peptides count")-1L)
tmp$Colour <- c("blue", pal)[tmp$`Peptides count`]
tmp2 <- summary(PG$`Peptides count`)
tmp2 <- data.frame(Variable = c("Protein groups", "Protein groups with 2+ peptidoforms", "", names(tmp2)),
                   Value = c(as.character(c(nrow(PG), sum(PG$"Peptides count" >= 2L))),
                             "",
                             as.character(signif(as.numeric(tmp2), 3L))))
tmp2$Txt <- apply(tmp2[, c("Variable", "Value")], 1L, \(x) {
  x <- setdiff(x, "")
  x <- if (length(x)) { paste(x, collapse = ": ") } else { "" }
  return(x)
})
tmpTxt <- paste(tmp2$Txt, collapse = "\n")
# tmp2$X <- max(as.numeric(tmp2$Value[match("Max.", tmp2$Variable)]))*0.98
# tmp2$Y <- max(tmp$"Protein groups")*(0.98-(0L:(nrow(tmp2) - 1L))*0.02)
ttl <- "Peptidoforms per PG"
plot <- ggplot2::ggplot(tmp) +
  ggplot2::geom_col(ggplot2::aes(x = `Peptides count`, y = `Protein groups`, fill = Colour),
                    colour = NA) +
  ggplot2::theme_bw() + ggplot2::ggtitle(ttl) +
  #ggplot2::geom_text(data = tmp2, ggplot2::aes(x = X, y = Y, label = Txt), hjust = 1, size = 3L) +
  ggplot2::geom_text(label = tmpTxt,
            x = max(as.numeric(tmp2$Value[match("Max.", tmp2$Variable)]))*0.8,
            y = max(tmp$"Protein groups")*0.8,
            hjust = 0.5,
            vjust = 0.5,
            size = 3L) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_x_continuous(breaks = seq(10L, floor(max(tmp$`Peptides count`)/10L)*10L, by = 10L),
                              expand = expansion(mult = c(0, 0.05))) +
  ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
# ggplot2::coord_trans(x = "log10", y = "log10")
poplot(plot, 12L, 22L) # This type of QC plot does not need to pop up, the side panel is fine
qcDir <- paste0(wd, "/Summary plots")
if (scrptType == "withReps") { dirlist <- union(dirlist, qcDir) }
if (!dir.exists(qcDir)) { dir.create(qcDir, recursive = TRUE) }
suppressMessages({
  ggplot2::ggsave(paste0(qcDir, "/", ttl, ".svg"), plot, dpi = 300L)
})
plotLy <- plotly::ggplotly(plot, tooltip = c("x", "y"))
plotLy <- plotly::config(plotLy,
                         modeBarButtonsToRemove = c("select2d", "lasso2d"))
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
