options(stringsAsFactors = FALSE)
require(proteoCraft)
require(ggplot2)
wd <- choose.dir("...Search_Folder/")
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
