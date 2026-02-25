###############################################################
#                                                             #
#      Example script showing how to load a parquet file      #
#                                                             #
###############################################################

# 1/ Make sure arrow is installed
# Modify this if you do want to use pak as a package manager... (Though why would you? It's vastly superior to anything base R has to offer!)
if (!require(pak)) { install.packages("pak") }
if (!require(arrow)) { pak::pkg_install("arrow") }
library(arrow)

# 2/ Select target directory
# (change path if necessary!)
wd <- rstudioapi::selectDirectory(path = dirname(rstudioapi::getActiveDocumentContext()$path))

# 3/ PSMs file
parquetFl <- paste0(wd, "/report.parquet")
stopifnot(file.exists(parquetFl))

# 4/ Load data as data.frame
data <- arrow::read_parquet(parquetFl)

# 5/ Optional: save as .tsv
tsvFl <- paste0(wd, "/report.tsv")
if (!file.exists(tsvFl)) {
  if (!require(data.table)) { pak::pkg_install("data.table") }
  data.table::fwrite(data, tsvFl, sep = "\t", row.names = FALSE, na = "NA")
}

# Done!
