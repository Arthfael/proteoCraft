# - Check variance/intensity dependency before or after normalisation
ttl <- "Peptides Variance vs Intensity dependency"
kol <- grep(topattern(pep.ref[rfnm]), colnames(pep), value = TRUE)
temp <- pep[, kol]
w <- which(apply(temp, 1, function(x) { length(is.all.good(x)) }) > 0)
temp <- temp[w,]
Aggr <- VPAL
if (LocAnalysis) { Aggr <- parse.Param.aggreg("Exp;Com") }
tst <- setNames(lapply(Aggr$values, function(x) { #x <- Aggr$values[2]
  em <- Exp.map[which(Exp.map[[Aggr$column]] == x),]
  if (nrow(em)) {
    kl <- paste0(pep.ref[rfnm], em$Ref.Sample.Aggregate)
    return(rowMeans(temp[, kl, drop = FALSE], na.rm = TRUE))
  } else {
    return(NULL)
  }
}), Aggr$values)
tst <- tst[which(!sapply(tst, is.null))]
tst <- as.data.frame(do.call(cbind, tst))
kol <- colnames(tst) <- cleanNms(colnames(tst))
tst <- apply(tst, 1, function(x) { kol[which(x == max(x))][1] })
tst2 <- as.data.table(ev[, c("Charge", "Modified sequence")])
tst2 <- tst2[, list(x = round(mean(Charge))), by = list(Group.1 = `Modified sequence`)]
tst2 <- as.data.frame(tst2)
temp2 <- data.frame("log10(Mean)" = log10(rowMeans(temp, na.rm = TRUE)),
                    "log10(Variance)" = log10(apply(temp, 1, var, na.rm = TRUE)),
                    "Strongest in..." = tst,
                    check.names = FALSE)
w <- which(!is.na(temp2$`Strongest in...`))
temp2 <- temp2[w,]
temp2$"Main charge" <- paste0("Z = ", tst2$x[match(pep$`Modified sequence`[w], tst2$Group.1)])
temp2$"Main charge" <- factor(temp2$"Main charge", levels = paste0("Z = ", as.character(1:8)))
dir <- paste0(wd, "/Workflow control")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
plot <- ggplot(temp2) +
  geom_scattermore(aes(x = `log10(Mean)`, y = `log10(Variance)`, colour = `Strongest in...`),
                   size = 1, alpha = 1) + ggtitle(ttl, subtitle = rfnm) + coord_fixed(0.3) + theme_bw() +
  scale_colour_viridis_d(begin = 0.25) +
  facet_grid(`Strongest in...` ~ `Main charge`) + theme(strip.text.y.right = element_text(angle = 0))
#poplot(plot, 12, 22)
ggsave(paste0(dir, "/", ttl, " - ", rfnm, ".jpeg"), plot, width = 10, height = 10, units = "in", dpi = 300)
ggsave(paste0(dir, "/", ttl, " - ", rfnm, ".pdf"), plot, width = 10, height = 10, units = "in", dpi = 300)
ReportCalls <- AddPlot2Report(Title = paste0(ttl, " - ", rfnm))
#

