options(stringsAsFactors = FALSE)
wd <- "D:\\groups_temp\\"
wd <- "E:\\"
setwd(wd)
require(tools)

pat <- "Regulation analysis - state of the art\\.R$"
pat <- "Regulation analysis - master script\\.R$"
pat <- "Regulation analysis - detailed script\\.R$"
regscripts <- list.files(wd, pat, recursive = TRUE)

regscripts <- data.frame(path = regscripts, dir = dirname(regscripts))
regscripts$"Last.modified_(char)" <- sapply(regscripts$path, function(x) { #x <- regscripts$path[1]
  as.character(file.info(x)$mtime)
})
regscripts$Last.modified <- sapply(regscripts$path, function(x) { #x <- regscripts$path[1]
  file.info(x)$mtime
})
regscripts$Size <- sapply(regscripts$path, function(x) { #x <- regscripts$path[1]
  file.info(x)$size
})
regscripts$MD5sums <- md5sum(regscripts$path)
regscripts <- regscripts[order(as.numeric(regscripts$Last.modified), decreasing = TRUE),]

refscript <- "H:\\R scripts\\Regulation analysis - latest\\Regulation analysis - state of the art.R"
refscriptmod <- file.info(refscript)$mtime
as.numeric(refscriptmod) >= max(as.numeric(regscripts$Last.modified))
as.numeric(refscriptmod) > max(as.numeric(regscripts$Last.modified))
refscriptmd5 <- md5sum(refscript)
w <- which(regscripts$MD5sums == refscriptmd5)
if (length(w)) { for (i in w) { print(regscripts$path[i]) } }
