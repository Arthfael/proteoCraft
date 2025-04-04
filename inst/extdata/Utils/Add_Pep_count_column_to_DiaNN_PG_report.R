require(reshape2)
require(rstudioapi)
require(ggplot2)

if ((exists("DIANNdir"))&&(!is.null(DIANNdir))) { dflt <- DIANNdir } else { dflt <- "...Search_Folder/*" } 
DIANNdir <- selectDirectory("Select DiaNN output folder", path = "...Search_Folder/*")
setwd(DIANNdir)

fls <- list.files()
tsvs <- grep("\\.tsv$", fls, value = TRUE)
pgsrepFl <- grep("\\.pg_matrix\\.tsv$", tsvs, value = TRUE)
pgsrepFl <- grep("-first-pass\\.pg_matrix\\.tsv$", pgsrepFl, invert = TRUE, value = TRUE)
reportFl <- gsub("\\.pg_matrix\\.tsv$", "\\.tsv", grep("\\.pg_matrix\\.tsv$", pgsrepFl, value = TRUE))

report <- read.delim(reportFl)
pgsrep <- read.delim(pgsrepFl)

# Annotate with peptide numbers
tmp <- aggregate(report$Protein.Group, list(report$Modified.Sequence), unique)
test <- sapply(tmp$x, length) # Should always be 1, but isn't for now; simple fix:
#View(tmp[which(test > 1),])
if (max(test) > 1) {
  tmp2 <- tmp[which(test>1),]
  tmp <- tmp[which(test==1),]
  tmp2 <- melt(setNames(tmp2$x, tmp2$Group.1))
  colnames(tmp2) <- c("x", "Group.1")
  tmp <- rbind(tmp2, tmp)
}
# Unlist, now those do not need to be lists
tmp$x <- unlist(tmp$x)
tmp$Group.1 <- unlist(tmp$Group.1)
#
tmp <- aggregate(tmp$Group.1, list(tmp$x), length)
pgsrep$"Number of peptides" <- tmp$x[match(pgsrep$Protein.Group, tmp$Group.1)]
sum(pgsrep$`Number of peptides` > 1)
summary(pgsrep$`Number of peptides`)
pgsrep$Protein.Names[which(pgsrep$`Number of peptides` == max(pgsrep$`Number of peptides`))]
write.table(pgsrep, gsub("\\.tsv$", "_PepNb.tsv", pgsrepFl), sep = "\t", quote = FALSE, row.names = FALSE)

#
plot <- ggplot(pgsrep) + geom_histogram(aes(x = `Number of peptides`, fill = `Number of peptides`, group = `Number of peptides`), binwidth = 1) + theme_bw() +
  scale_fill_viridis_c()
proteoCraft::poplot(plot)
