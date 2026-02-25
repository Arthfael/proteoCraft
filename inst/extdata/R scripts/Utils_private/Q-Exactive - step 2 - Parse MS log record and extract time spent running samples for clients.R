require(aRmel)
require(svDialogs)
#require(openxlsx)
PrcsDir <- "B:/group/mspecgrp/Prices"
stopifnot(dir.exists(PrcsDir))

dflt <- rev(unlist(strsplit(date(), " ")))[1]
msg <- "Which year do you want values for?"
RunTimesFl <- normalizePath(choose.files(paste0(PrcsDir, "/*.RData"), multi = FALSE), winslash = "/")
load(RunTimesFl)
RunTimes <- Res
rm(Res)
RunTimes$date <- as.POSIXct(RunTimes$date)
# Year
RunTimes$Year <- gsub("-.+", "", RunTimes$start)
Years <- unique(RunTimes$Year)
Year <- as.integer(dlg_list(Years, as.integer(gsub("-.+", "", Sys.Date()))-1, title = "Which year do you want to extract? (escape to process all years)")$res)
wd <- paste0(PrcsDir, "/", Year, "/")
if (!dir.exists(wd)) { dir.create(wd, recursive = TRUE) }
setwd(wd)

RunTimes$"sequence path" <- gsub("\\\\", "/", RunTimes$"sequence path")
RunTimes$ForClient <- grepl("^[A-Z]:/Data/(([^/]+gr(ou)?p)|(Solgate)|(Valanx)|(Kondrashov))/", RunTimes$"sequence path")
if ((length(Year) == 1)&&(!is.na(Year))) {
  RunTimes <- RunTimes[grep(paste0("^", Year, "-"), RunTimes$start),]
}
# Inspect both classes of folders to check that the regex works well:
ClientProjects <- sort(unique(gsub("/[^/]+\\.sld$", "", gsub("^[A-Z]:/Data/", "", RunTimes$"sequence path"[which(RunTimes$ForClient)]))))
ClientFolders <- sort(unique(gsub("/.+", "", ClientProjects)))
OtherFolders <- sort(unique(gsub("/.+", "", gsub("^[A-Z]:/Data/", "", RunTimes$"sequence path"[which(!RunTimes$ForClient)]))))
#RunTimes$Sequence.exists <- file.exists(RunTimes$"sequence path")
# Having sld files would allow parsing them and fixing truncated file names
# I need to update the raw file archiving script to also transfer those files

ClientFolders
OtherFolders
Blanks <- RunTimes$blank == "+"
AllClientSamples <- sum(RunTimes$ForClient)
ClientBlanks <- sum((RunTimes$ForClient)&(RunTimes$blank == "+"))
ClientSamples <- AllClientSamples-ClientBlanks
AvClienRunTime <- sum(RunTimes$"time min"*RunTimes$ForClient*(!Blanks))/(ClientSamples)
BatchSize <- (ClientSamples)/length(ClientProjects)

Days <- max(RunTimes$end) - min(RunTimes$start)
print(paste0("Time covered by the record: ", round(Days, 2), " days"))
Days <- as.numeric(Days)
ttl <- paste0("Q-Exactive HF instrument time use report - ", Sys.Date())
if ((length(Year) == 1)&&(!is.na(Year))) {
  ttl <- paste0(Year, " ", ttl)
}
FracRunTime <- (sum(RunTimes$"time min")/(24*60))/Days
FracClients2TotTime <- (sum(RunTimes$ForClient*RunTimes$"time min")/(24*60))/Days
FracClients2RunTime <- sum(RunTimes$ForClient*as.numeric(RunTimes$"time min"))/sum(as.numeric(RunTimes$"time min"))
Report <- c(ttl,
            paste(rep("-", nchar(ttl)), collapse = ""),
            "",
            paste0("Tot. number of client samples (including blanks): ", AllClientSamples),
            paste0("Tot. number of client samples (not including blanks): ", ClientSamples),
            paste0("Ratio client samples/blanks: ", round(ClientSamples/ClientBlanks, 2)),
            paste0("Av. client run time: ", round(AvClienRunTime/60, 2), " h"),
            paste0("Batches in ", Year, ": ", length(ClientProjects)),
            paste0("Av. batch (project) size: ", round(BatchSize, 0), " samples"),
            paste0("Time spent running samples: ", as.character(round(sum(RunTimes$"time min")/60)), " h = ", round(100*FracRunTime, 1), "% of record scope"),
            paste0("Time spent running samples for clients: ", as.character(round(sum(RunTimes$"time min"*RunTimes$ForClient)/60)), " h = ", round(100*FracClients2TotTime, 1), "% of record scope"),
            paste0(round(100*FracClients2RunTime, 1), "% of actual run time was spent running samples for clients"))
Report
Report <- c(Report, "", "Projects:", ClientProjects, "")
print(Report)
write(Report, paste0(ttl, ".txt"))
writeClipboard(Report)
#openwd()

require(ggplot2)
RunTimes$Month <- gsub("-[0-9]+.+", "", gsub("^[0-9]+-", "", RunTimes$start))
w <- which(RunTimes$ForClient)
temp <- aggregate(RunTimes$"time min"[w], list(RunTimes$Month[w]), sum)
colnames(temp) <- c("Month", "Time (min)")
MonthLength <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
temp$Group <- "Ratio"
temp$"Ratio: used for clients / available" <- temp$"Time (min)"/(60*24*MonthLength[as.numeric(gsub("^[0-9]+-", "", temp$Month))])
ttl <- "Ratio - MS time used for client projects VS available"
lab <- "Average"
if ((length(Year) == 1)&&(!is.na(Year))) {
  ttl <- paste0(Year, " ", ttl)
  lab <- paste0(Year, " average")
}
temp$`Time (min)` <- as.numeric(temp$`Time (min)`)
temp$`Ratio: used for clients / available` <- as.numeric(temp$`Ratio: used for clients / available`)
plot <- ggplot(temp) +
  geom_line(aes(x = Month, y = `Ratio: used for clients / available`, group = Group)) +
  geom_hline(yintercept = as.numeric(FracClients2TotTime), colour = "red", linetype = "dashed") +
  annotate("text", y = as.numeric(FracClients2TotTime)+0.02, x = unique(temp$Month[10]), colour = "red", label = lab, hjust = 0.5) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0, 1) + ggtitle(ttl)
grDevices::windows(width = 20, height = 12); print(plot)
ggsave(paste0(ttl, ".jpeg"), plot, dpi = 150)
