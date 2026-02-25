########################
# Check MS performance #
########################
#
# This script will read and parse all instrument logs and plot all relevant values for a time interval of the user's choice.
# Each time interval gets saved into its own sub-folder.
# Each plot is created in two versions:
# - jpeg
# - interactive html
options(stringsAsFactors = FALSE)

# Install required packages (RTools should be installed)
packs <- c("assertthat", "ggplot2", "svDialogs", "openxlsx", "plotly", "htmlwidgets")
for (pack in packs) { if (!require(pack, character.only = TRUE)) { install.packages(pack, type = "win.binary", dependencies = TRUE) } }
for (pack in packs) { require(pack, character.only = TRUE) }

# Set local work directory, in which pictures will be saved
wd <- "C:/Users/Thermo/Desktop/MS performance plots"
if (!dir.exists(wd)) { dir.create(wd, recursive = TRUE) }
setwd(wd)

# Identify all log files
LogDir <- "C:/Xcalibur/system/Exactive/log" # (their directory)
LogFls <- grep("^InstrumentTemperature--[0-9]{4}-[0-9]{2}-[0-9]{2}\\.log$", list.files(LogDir), value = TRUE)
LogFls <- data.frame(File = LogFls,
                     Date = gsub("^InstrumentTemperature--|\\.log$", "", LogFls),
                     Year = gsub("^InstrumentTemperature--|-[0-9]{2}-[0-9]{2}\\.log$", "", LogFls),
                     Month = gsub("^InstrumentTemperature--[0-9]{4}-|-[0-9]{2}\\.log$", "", LogFls),
                     Day = gsub("^InstrumentTemperature--[0-9]{4}-[0-9]{2}-|\\.log$", "", LogFls))
# (Do not convert to integers except for sorting, otherwise left-hand side zeros could be lost and cause issues!)
LogFls$DDMMYYYY <- apply(LogFls[, c("Day", "Month", "Year")], 1, paste, collapse = "")
LogFls <- LogFls[order(as.integer(LogFls$Year), as.integer(LogFls$Month), as.integer(LogFls$Day), decreasing = FALSE),]

# Define limits of analysis
mStart <- mEnd <- NA
StartMsg <- paste0("Enter a start date (DDMMYYYY format, not earlier than ", LogFls$DDMMYYYY[1],")")
EndMsg <- paste0("Enter an end date (DDMMYYYY format, not later than ", LogFls$DDMMYYYY[nrow(LogFls)], ")")
Start <- dlg_input(StartMsg, LogFls$DDMMYYYY[1])$res
if (!is.na(Start)) { mStart <- match(Start, LogFls$DDMMYYYY) }
while ((is.na(Start))||(nchar(Start) != 8)||(is.na(mStart))) {
  Start <- dlg_input(gsub("^E", "Again: e", StartMsg), LogFls$DDMMYYYY[1])$res
  if (!is.na(Start)) { mStart <- match(Start, LogFls$DDMMYYYY) }
}
End <- dlg_input(EndMsg, LogFls$DDMMYYYY[nrow(LogFls)])$res
if (!is.na(End)) { mEnd <- match(End, LogFls$DDMMYYYY) }
while ((is.na(End))||(nchar(End) != 8)||(is.na(mEnd))) {
  End <- dlg_input(gsub("^E", "Again: e", EndMsg), LogFls$DDMMYYYY[nrow(LogFls)])$res
  if (!is.na(End)) { mEnd <- match(End, LogFls$DDMMYYYY) }
}
while (mEnd < mStart) {
  StartMsg <- "You chose a start date before the end date!"
  Start <- dlg_input(gsub("^E", "You chose a start date before the end date!\nAgain: e", StartMsg), LogFls$DDMMYYYY[1])$res
  if (!is.na(Start)) { mStart <- match(Start, LogFls$DDMMYYYY) }
  while ((is.na(Start))||(nchar(Start) != 8)||(is.na(mStart))) {
    Start <- dlg_input(gsub("^E", "Again: e", StartMsg), LogFls$DDMMYYYY[1])$res
    if (!is.na(Start)) { mStart <- match(Start, LogFls$DDMMYYYY) }
  }
  End <- dlg_input(gsub("^E", "Now, again: e", EndMsg), LogFls$DDMMYYYY[nrow(LogFls)])$res
  if (!is.na(End)) { mEnd <- match(End, LogFls$DDMMYYYY) }
  while ((is.na(End))||(nchar(End) != 8)||(is.na(mEnd))) {
    End <- dlg_input(gsub("^E", "Again: e", EndMsg), LogFls$DDMMYYYY[nrow(LogFls)])$res
    if (!is.na(End)) { mEnd <- match(End, LogFls$DDMMYYYY) }
  }
}

# Parse optional annotations table
# This is a table of events to overlay on each plot. You may want to include:
# - Any abnormal thing you noticed on the instrument
# - Times when the instrument is not operational
# - Any intervention by the engineer to fix and issue or perform preventive maintenance
# - Possibly mass calibrations
# Currently each event is plotted as a single blue line (jpeg) or arrow (html),
# but it would be easy to adapt the code to draw segments from a start to an end date.
AnnotDir <- "T:/Instrument log analysis"
Annotate <- "Events - Instrument 1.xlsx" %in% list.files(AnnotDir)
if (Annotate) {
  Annotations <- read.xlsx(paste0(AnnotDir, "/Events - Instrument 1.xlsx"), check.names = FALSE)
  Annotations$date <- as.POSIXct(strptime(apply(Annotations[, c("Day", "Month", "Year")], 1, paste, collapse = "/"),
                                          "%d/%m/%Y",
                                          Sys.timezone()), Sys.timezone())
  Annotations$event <- gsub("\\. ", ".\n", gsub(" - ", " -\n", gsub(": ", ":\n", gsub(", ", ",\n", Annotations$Event))))
  Annotations$Time <- as.integer(Annotations$date)
}

# Process log files
fls <- paste0(LogDir, "/", LogFls$File[mStart:mEnd])
tst <- unlist(lapply(fls, readLines))  
grhead <- grep("^Time \\[sec\\]\t", tst)
for (i in grhead) { # Do not ask me why it has to be done so... super inefficient
  if (i == rev(grhead)[1]) { j <- length(tst) } else { j <- c(1:length(tst))[grhead[match(i, grhead) + 1]-1] }
  write(tst[i:j], "temp.txt")
  temp <- read.delim("temp.txt", check.names = FALSE)
  temp <- temp[, which(colnames(temp) != "")]
  if (i == grhead[1]) { Logs <- temp } else {
    kol <- colnames(Logs)
    kol2 <- colnames(temp)
    kol2n <- kol2[which(!kol2 %in% kol)]
    koln <- kol[which(!kol %in% kol2)]
    Logs[, kol2n] <- NA
    temp[, koln] <- NA
    Logs <- rbind(Logs, temp)
  }
}
if (file.exists("temp.txt")) { unlink("temp.txt") }
Logs$date <- as.POSIXlt(Logs$Date)
Logs$"Instrument status" <- "Other"
Logs$"Instrument status"[which(round(Logs$`Capillary Temperature [?C]`) %in% c(310:330))] <- "HESI source"
Logs$"Instrument status"[which(round(Logs$`Capillary Temperature [?C]`) %in% c(240:280))] <- "Easy Spray source"
Logs$"Instrument status"[which(round(Logs$`Capillary Temperature [?C]`) %in% c(90:239))] <- "Exchanging capillary"
Logs$"Instrument status"[which(round(Logs$`Capillary Temperature [?C]`) < 90)] <- "Venting"
corecol <- c("Time [sec]", "Date", "Up-Time [days]", "date", "Instrument status")
othercol <- colnames(Logs)[which(!colnames(Logs) %in% corecol)]

# Create plots
Months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
Range <- data.frame(Day = substr(c(Start, End), 1, 2),
                    Month = Months[as.numeric(substr(c(Start, End), 3, 4))],
                    Year = substr(c(Start, End), 5, 8))
if (Range$Year[2] != Range$Year[1]) {
  Tag <- paste(apply(Range, 1, paste, collapse = "."), collapse = " - ")
} else {
  if (Range$Month[2] != Range$Month[1]) {
    Tag <- paste0(paste(apply(Range[, c("Day", "Month")], 1, paste, collapse = "."), collapse = "-"), ".", Range$Year[1])
  } else {
    if (Range$Day[2] != Range$Day[1]) {
      Tag <- paste0(paste(Range$Day, collapse = "-"), ".", paste(Range[1, c("Month", "Year")], collapse = "."))
    } else { Tag <- paste(Range[1,], collapse = ".") }
  }
}
TargDir <- paste0(wd, "/", Tag)
if (!dir.exists(TargDir)) { dir.create(TargDir, recursive = TRUE) }
Range$Time <- as.integer(ISOdatetime(as.integer(Range$Year), match(Range$Month, Months), as.integer(Range$Day), c(0, 23), c(0, 0), c(0, 59)))
wRg <- which((Logs$`Time [sec]` >= Range$Time[1])&(Logs$`Time [sec]` <= Range$Time[2]))
if (length(wRg)) {
  for (col in othercol) { #col <- othercol[1]
    ylab <- col
    M <- max(Logs[[col]], na.rm = TRUE)
    m <- min(Logs[[col]], na.rm = TRUE)
    temp <- data.frame(date = Logs$date[wRg], Y = Logs[[col]][wRg],
                       Group = col, "Instrument status" = Logs$"Instrument status"[wRg], check.names = FALSE)
    ttl <- gsub("/", "-",  gsub(" +", " ", gsub("(\\(|\\[).+(\\)|\\])", "", col)))
    if ((M != m)&&(M/m > 1000)&&(m > 0)) {
      temp$Y <- log10(temp$Y)
      ylab <- paste0("log10(", col, ")")
      M <- log10(M)
      m <- log10(m)
    }
    yScl <- M-m
    # ggplot jpeg version
    plot1 <- ggplot(temp, aes(x = date, y = Y)) + geom_point(aes(group = Group, color = `Instrument status`), size = 0.1) +
      ylab(ylab) + ggtitle(ttl, subtitle = Tag) + theme_bw()
    if (Annotate) {
      wRg2 <- which((Annotations$Time >= Range$Time[1])&(Annotations$Time <= Range$Time[2]))
      if (length(wRg2)) {
        # plotly html version
        plot2 <- suppressWarnings(plot1 +
                                    geom_point(data = Annotations[wRg2,], aes(x = date, text = event), y = m-yScl*0.025, fill = "blue", color = "blue", size = 3, shape = 17) +
                                    #geom_segment(data = Annotations, aes(x = date, xend = date, text = event),
                                    #             y = m, yend = M+yScl*0.05, colour = "blue", linetype = "dotted", alpha = 0.5) +
                                    ylim(m-yScl*0.05, M+yScl*0.05))
        plot2 <- ggplotly(plot2, tooltip = "text")
        plot1 <- plot1 +
          geom_vline(data = Annotations[wRg2,], aes(xintercept = date), colour = "blue", linetype = "dashed", alpha = 0.5) +
          geom_text(data = Annotations, aes(x = date+1, label = Event),
                    y = m+yScl, colour = "blue", angle = 60, size = 2, hjust = 0, vjust = 1) +
          ylim(m, M+yScl)
      }
    } else { plot2 <- ggplotly(plot1) }
    ggsave(paste0(TargDir, "/", ttl, " - ", Tag, ".jpeg"), plot1, dpi = 150, width = 20, height = 5)
    #windows(20, 12); print(plot1)
    # plotly
    saveWidget(partial_bundle(plot2), paste0(TargDir, "/", ttl, " - ", Tag, ".html"))
    #system(paste0("open \"", TargDir, "/", ttl, " - ", Tag, ".html"))
  }
} else { warning("No data to plot within selected time range!") }

