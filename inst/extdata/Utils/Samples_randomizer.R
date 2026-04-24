# Samples randomizer
if (!require(svDialogs)) { install.packages("svDialogs") }
require(svDialogs)

if (exists("SampleGrps")) { dflt <- paste(SampleGrps, collapse = ";") } else { dflt <- "Treatment;Control"}
SampleGrps <- unique(unlist(strsplit(dlg_input("Enter names of sample groups (semicolon-separated)", "Treatment;Control")$res, ";")))
# Example:
#Samples <- as.character(sapply(c("Striatum", "Cortex"), function(x) { paste0(x, "_", c("WT", "KO")) }))
NRep <- as.numeric(dlg_input("Enter the number of replicates per group", 3)$res)
#NRep <- 3
Samples <- as.character(sapply(SampleGrps, function(x) { paste0(x, "_", 1:NRep) }))
l <- length(Samples)

Samples <- sample(Samples, l)
print(Samples)
writeClipboard(Samples)
