options(stringsAsFactors = FALSE)
library(proteoCraft)
library(ggplot2)
library(rstudioapi)

wd <- selectDirectory(path = "D:\\groups_temp")
setwd(wd)
MQ.load()

temp <- cor_mod_seq(ev)
Modifs <- temp$PTMs
ev <- temp$Peptides
ev <- ev[which(is.all.good(ev$Intensity, 2)),]

Counts <- matrix(rep(0, 6), ncol = 2)
colnames(Counts) <- c("TMT", "Unlabelled")
rownames(Counts) <- c("Reverse", "Contaminants", "Sample")
t1 <- !is.na(ev$Reverse) & ev$Reverse == "+"
g1 <- grepl("tm|kt", ev$Modified.sequence[which(t1)])
t2 <- !is.na(ev$Potential.contaminant[which(!t1)]) & ev$Potential.contaminant[which(!t1)] == "+"
g2 <- grepl("tm|kt", ev$Modified.sequence[which(!t1)][which(t2)])
g3 <- grepl("tm|kt", ev$Modified.sequence[which(!t1)][which(!t2)])
Counts[1, 1] <- sum(g1)
Counts[1, 2] <- sum(!g1)
Counts[2, 1] <- sum(g2)
Counts[2, 2] <- sum(!g2)
Counts[3, 1] <- sum(g3)
Counts[3, 2] <- sum(!g3)
Counts[,1]/Counts[,2]

Comp <- data.frame(Sequence = unique(ev$Sequence))
g <- grepl("tm|kt", ev$Modified.sequence)
temp <- aggregate(ev$Intensity[which(g)], list(ev$Sequence[which(g)]), sum, na.rm = TRUE)
Comp$TMT <- NA
w <- which(Comp$Sequence %in% temp$Group.1)
Comp$TMT[w] <- temp$x[match(Comp$Sequence[w], temp$Group.1)]
temp <- aggregate(ev$Intensity[which(!g)], list(ev$Sequence[which(!g)]), sum, na.rm = TRUE)
Comp$Unlabelled <- NA
w <- which(Comp$Sequence %in% temp$Group.1)
Comp$Unlabelled[w] <- temp$x[match(Comp$Sequence[w], temp$Group.1)]
Comp$"TMT (log10)" <- log10(Comp$TMT)
Comp$"Unlabelled (log10)" <- log10(Comp$Unlabelled)
w <- which(apply(Comp[, c("TMT (log10)", "Unlabelled (log10)")], 1, function(x) { length(is.all.good(x)) }) == 2)

linearMod <- lm(`TMT (log10)` ~ `Unlabelled (log10)`, data =  Comp)
print(linearMod)

plot <- ggplot(Comp[w,]) + geom_point(aes(x = log10(Unlabelled), y = log10(TMT))) +
  geom_abline(slope = linearMod$coefficients[[2]], intercept = linearMod$coefficients[[1]], colour = "red") +
  theme_bw()
poplot(plot)

