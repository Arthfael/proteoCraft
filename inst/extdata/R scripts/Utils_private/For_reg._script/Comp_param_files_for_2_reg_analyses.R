options(stringsAsFactors = FALSE)
require(proteoCraft)

RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
dfltLocsFl <- paste0(homePath, "/Default_locations.xlsx")
dfltLocs <- openxlsx2::read_xlsx(dfltLocsFl)
searchDir <- dfltLocs$Path[match("Search folder", dfltLocs$Folder)]

wd1 <- rstudioapi::selectDirectory(path = searchDir)
wd2 <- rstudioapi::selectDirectory(path = wd1)

setwd(wd1)
param1 <- read.csv("parameters.csv", header = FALSE)
setwd(wd2)
param2 <- read.csv("parameters.csv", header = FALSE)

param <- data.frame(V1 = unique(c(param1$V1, param2$V1)))
for (i in 1:2) { # i <- 1
  param[[paste0("Param", i)]] <- ""
  temp <- get(paste0("param", i))
  w <- which(param$V1 %in% temp$V1)
  param[w, paste0("Param", i)] <- temp$V2[match(param$V1[w], temp$V1)]
}
w <- which(param$Param1 != param$Param2)
if (length(w) > 0) { View(param[w,]) }

kol <- 2
w <- which(param1[[kol]] != param2[[kol]])
print(w)
View(cbind(param1[[1]], param1[[kol]], param2[[kol]])[w,])

kol <- 3
w <- which(param1[[kol]] != param2[[kol]])
if (length(w) > 0) {
  print(w)
  View(cbind(param1[[kol]], param2[[kol]])[w,])
}

openwd()
