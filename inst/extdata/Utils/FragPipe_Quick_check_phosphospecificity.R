wd <- "MY FRAGPIPE OUTPUT DIRECTORY"
setwd(wd)
require(data.table)
tmp <- data.table::fread(paste0(wd, "/psm.tsv"), integer64 = "numeric", check.names = FALSE, data.table = FALSE)

g <- grep("\\(79\\.9663\\)", tmp$`Assigned Modifications`)
spec <- length(g)/nrow(tmp)

cat(paste0("Phospho-specificity: ", round(spec*100, 1), "%"))
