# Prepare data for clustering and dimensionality reduction plots
Src <- paste0(libPath, "/extdata/Sources/cluster_Heatmap_Prep.R")
#rstudioapi::documentOpen(Src)
dataType <- "PG"
source(Src)
if (scrptType == "withReps") {
  dataType <- "peptides"
  source(Src)
}
