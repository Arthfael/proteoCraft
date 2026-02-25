options(stringsAsFactors = FALSE)
require(proteoCraft)
require(ggplot2)

RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
dfltLocsFl <- paste0(homePath, "/Default_locations.xlsx")
dfltLocs <- openxlsx2::read_xlsx(dfltLocsFl)
searchDir <- dfltLocs$Path[match("Search folder", dfltLocs$Folder)]

wd <- rstudioapi::selectDirectory(path = searchDir)
setwd(wd)
MQ.load(pep = FALSE, prot = FALSE)
plot <- ggplot(ev) + geom_density(stat = "density", aes(Retention.length*60))
poplot(plot)
scale <- (max(ev$Retention.time)-min(ev$Retention.time))/10
midrange <- (min(ev$Retention.time)+max(ev$Retention.time))/2
w <- which((ev$Retention.time <= midrange-scale*2)&(ev$Retention.time <= midrange+scale*2))
cpw <- round(mean(ev$Retention.length[w])*6)*10
message(paste0("Values to use:\n - Chromatographic peak width: ", cpw, ",\n - Dynamic exclusion window: ", cpw/2))

msmsScans <- read.delim("msmsScans.txt")
msms <- read.delim("msms.txt")

test <- aggregate(ev$id, list(ev$MS.MS.scan.number), list)
test <- test[order(sapply(test$x, length), decreasing = TRUE),]
colnames(test) <- c("MS.MS.scan", "Evidence IDs")

max(msms$Scan.number)
max(msmsScans$MS.MS.IDs)
max(msmsScans$)
