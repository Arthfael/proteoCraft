# Script to extract timsTOF HT total times per year...
# This will attempt to extract the information from both the archive and the timsTOF PC
thisYear <- as.integer(gsub("-.*", "", Sys.Date()))
wd <- paste0("B:/group/mspecgrp/Prices/", thisYear)
if (!dir.exists(wd)) { dir.create(wd, recursive = TRUE) }
setwd(wd)

topattern <- function (x, start = TRUE, end = FALSE, collapse = "|") {
  x <- gsub("\\\\", "\\\\\\\\", as.character(x))
  x <- gsub("\\.", "\\\\.", x)
  x <- gsub("\\*", "\\\\*", x)
  x <- gsub("\\$", "\\\\$", x)
  x <- gsub("\\^", "\\\\^", x)
  x <- gsub("\\+", "\\\\+", x)
  x <- gsub("\\?", "\\\\?", x)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  x <- gsub("\\[", "\\\\[", x)
  x <- gsub("\\]", "\\\\]", x)
  x <- gsub("\\(", "\\\\(", x)
  x <- gsub("\\)", "\\\\)", x)
  x <- gsub("\\|", "\\\\|", x)
  if (start) {
    x <- paste0("^", x)
  }
  if (end) {
    x <- paste0(x, "$")
  }
  if ((length(x) > 1) && (collapse != FALSE)) {
    x <- paste(x, collapse = collapse)
  }
  return(x)
}

# This script should be run on the timsTOF's control PC, so that it can access both
# - the MS archive
# - any local files not yet transferred to the archive
require(parallel)
require(XML)
parClust <- makeCluster(detectCores()-1)

targDirs <- c("B:/archive/mspecgrp/MS/Acquired_data",
              "B:/group/lsfgrp/Mass_Spec/Acquired_data_v2",
              "B:/group/mspecgrp/_Archive/Raw",
              "D:/Data/Projects"
              )
xsTst <- dir.exists(targDirs)
stopifnot(sum(!xsTst) == 0)
dDrsTblsLst  <- setNames(parLapply(parClust, targDirs, function(dr) {
  #if (dir.exists(dr)) {
  tmp <- list.dirs(dr, TRUE, TRUE)
  Ds <- grep("[0-9]+\\.d$", tmp, value = TRUE)
  return(data.frame(Raw.file = Ds,
                    Source = dr,
                    check.names = FALSE))
}), targDirs)
vapply(targDirs, function(x) {
  length(dDrsTblsLst[[x]]$Raw.file)
}, 1)
dDrsTbl <- plyr::rbind.fill(dDrsTblsLst )
dDrsTbl$Run_ID <- as.integer(gsub(".*_|\\.d$", "", dDrsTbl$Raw.file))
w <- which(is.na(dDrsTbl$Run_ID))
if (length(w)) {
  print(dDrsTbl$Raw.file[w])
  dDrsTbl$Run_ID[w] <- max(dDrsTbl$Run_ID, na.rm = TRUE)*100 + 1:length(w)
}
dDrsTbl <- dDrsTbl[order(dDrsTbl$Run_ID),] # May also make the parLapply below more efficient
dDrsTbl$FlInfo <- paste0(dDrsTbl$Raw.file, "/SampleInfo.xml")
dDrsTbl$FlInfo.exists <- file.exists(dDrsTbl$FlInfo)
wNoInfo <- which(!dDrsTbl$FlInfo.exists)
if (length(wNoInfo)) {
  warning("Some .d folders do not seem to contain a \"SampleInfo.xm;\" file... (probably failed injections)")
  View(dDrsTbl[wNoInfo,])
  #proteoCraft::openwd(dDrsTbl$Raw.file[wNoInfo[1]])
}
dDrsTbl <- dDrsTbl[which(dDrsTbl$FlInfo.exists),]
dDrsTbl$LC_method_fl <- parLapply(parClust, dDrsTbl$FlInfo, function(fl) {
  res <- ""
  try({
    #fl <- dDrsTbl$FlInfo[w][1]
    x <- XML::xmlToList(fl)
    if (".attrs" %in% names(x$Sample)) { res <- x$Sample$.attrs["Method"] } else {
      res <- x$Sample["Method"]
    }
    res <- gsub("\\?.+", "", gsub("\\\\", "/", res))
  }, silent = TRUE)
  return(res)
})

dDrsTbl$LC_method <- paste0(dDrsTbl$Raw.file, "/", dDrsTbl$Run_ID, ".m/hystar.method")
dDrsTbl$LC_method.exists <- file.exists(dDrsTbl$LC_method)
wN <- which(!dDrsTbl$LC_method.exists)
if (length(wN)) {
  warning("Some .d folders do not seem to contain an LC method! Investigate")
  View(dDrsTbl[wN,])
  #proteoCraft::openwd(dDrsTbl$Raw.file[wN[1]])
}
dDrsTbl <- dDrsTbl[which(dDrsTbl$LC_method.exists),]
#
uDrs <- aggregate(dDrsTbl$Source, list(dDrsTbl$Run_ID), length)
w <- which(uDrs$x == 2)
w <- which((!dDrsTbl$Run_ID %in% uDrs$Group.1[w])|((dDrsTbl$Run_ID %in% uDrs$Group.1[w])&(dDrsTbl$Source == targDirs[1])))
dDrsTbl <- dDrsTbl[w,]
#
dDrsTbl$LC_meth <- parLapply(parClust, dDrsTbl$LC_method, function(fl) { #fl <- LC_method[1]
  tst <- FALSE
  lc <- NA
  x <- try(XML::xmlToList(fl), silent = TRUE)
  if (!"try-error" %in% class(x)) {
    y <- try(x$LCMethodData$ModuleMethods$ModuleMethodData$text, silent = TRUE)
    if (!"try-error" %in% class(y)) {
      lc <- try(XML::xmlToList(y), silent = TRUE)
      tst <- TRUE
    }
  }
  return(list(Success = tst, Method = lc))
})
dDrsTbl <- dDrsTbl[which(sapply(dDrsTbl$LC_meth, function(x) { x$Success })),]
dDrsTbl$LC_meth <- lapply(dDrsTbl$LC_meth, function(x) { #x <- dDrsTbl$LC_meth[[1]]
  x$Method
})
dDrsTbl$LC_method_length <- vapply(dDrsTbl$LC_meth, function(x) {
  #x <- dDrsTbl$LC_meth[[1]]
  x <- try(x$HyStarMethodData$NoStandardMethodData, silent = TRUE)
  if (!"try-error" %in% class(x)) {
    if (length(x) > 1) {
      x <- x$TotalRunTime
    }
    return(as.numeric(x))
  } else { return(NA) }
}, 1)
tst <- vapply(dDrsTbl$LC_method_length, length, 1)
dDrsTbl$LC_method_length[which(tst == 3)[1]]

dDrsTbl <- dDrsTbl[which(!is.na(dDrsTbl$LC_method_length)),]
tmp <- parLapply(parClust, dDrsTbl$Raw.file, function(dFl) { #dFl <- dDrsTbl$Raw.file[1] #dFl <- rev(dDrsTbl$Raw.file)[1]
  tst <- try({
    sbdr <- grep("\\.m$", list.dirs(dFl), invert = TRUE, value = TRUE)
    sbdr <- sbdr[which(sbdr != dFl)]
    lgFl <- grep("execution-log", list.files(sbdr, full.names = TRUE), value = TRUE)
    if (length(lgFl) != 1) { stop(dFl) }
    lg <- readLines(lgFl)
    #system(paste0("open \"", lgFl, "\""))
    strt <- grep("STARTED", lg, value = TRUE)
    strt <- unlist(strsplit(gsub("\t.*", "", strt), " - "))[1]
    stp <- grep("COMPLETED", lg, value = TRUE)
    if (!length(stp)) { stp <- grep("ABORTED", lg, value = TRUE) }
    if (!length(stp)) { stop(dFl) }
    end <- unlist(strsplit(gsub("\t.*", "", stp), " - "))[1]
    stp <- unlist(strsplit(gsub("\t.*", "", stp), " - "))[2]
    if (!grep("^000\\.", stp)) { stop("This doesn't (yet) support methods this long!") }
    tm <- unlist(strsplit(gsub("^000\\.", "", stp), ":"))
    strt <- gsub("/", "-", strt)
    strt_day <- unlist(strsplit(gsub(" .*", "", strt), "-"))
    w <- which(nchar(strt_day) == 4)
    if (w == 1) {
      end <- gsub("/", "-", end)
      end_day <- unlist(strsplit(gsub(" .*", "", end), "-"))
      strt_day <- paste(rev(strt_day), collapse = "-")
      end_day <- paste(rev(end_day), collapse = "-")
      strt <- paste0(strt_day, " ", gsub("^[^ ]+ ", " ", strt))
      end <- paste0(end_day, " ", gsub("^[^ ]+ ", " ", end))
    }
  }, silent = TRUE)
  if (!"try-error" %in% class(tst)) {
    res <- list(Success = TRUE, Date = strt, End = end, Length = sum(as.numeric(tm) * c(1, 1/60, 1/3600)))
  } else { res <- list(Success = FALSE) }
  return(res)
})
w <- which(vapply(tmp, function(x) { x$Success }, TRUE))
dDrsTbl <- dDrsTbl[w,]
tmp <- tmp[w]
dDrsTbl$Date <- sapply(tmp, function(x) { x$Date })
dDrsTbl$End <- sapply(tmp, function(x) { x$End })
# Unfortunately, Month and Days are mixed up:
# Some %$%^U$$$^$  - possibly me - changed from the horrible US date format to a rational one at some point.
dDrsTbl[, c("Month", "Day", "Year")] <- as.data.frame(t(sapply(strsplit(gsub(" .*", "", dDrsTbl$Date), "-"),
                                                               function(x) { as.integer(unlist(x))[1:3] })))
tst <- aggregate(dDrsTbl$Raw.file, list(dDrsTbl$Source, dDrsTbl$Year), function(x) {
  length(unique(x))
})
tst <- reshape::cast(tst, Group.1~Group.2, value = "x")
colnames(tst)[1] <- "Location"
tst
Months <- c("January", "February", "March", "April", "May", "June",
            "July", "August", "September", "October", "November", "December")
dDrsTbl$Month <- Months[as.integer(dDrsTbl$Month)]
dDrsTbl$byMonth <- do.call(paste, c(dDrsTbl[, c("Month", "Year")], sep = "-"))
lev <- as.character(sapply(sort(unique(dDrsTbl$Year)), function(x) {
  paste0(Months, "-", x)
}))
dDrsTbl$byMonth <- factor(dDrsTbl$byMonth, levels = lev)
dDrsTbl$Folder <- gsub("/[^/]+$", "", dDrsTbl$Raw.file)
dDrsTbl$Dataset <- gsub("^.*/", "", dDrsTbl$Folder)
dDrsTbl$"Run length (h)" <- sapply(tmp, function(x) { x$Length })
dDrsTbl$"Run length (min)" <- dDrsTbl$"Run length (h)"*60
scale <- summary(c(dDrsTbl$LC_method_length, dDrsTbl$"Run length (min)"))
scale <- scale[c("Min.", "Max.")]
require(ggplot2)
require(viridis)
require(plotly)
require(htmlwidgets)
plot <- ggplot(dDrsTbl) + geom_point(aes(x = LC_method_length, y = `Run length (min)`, color = Run_ID,
                                         text = Date), size = 0.3) + theme_bw() +
  geom_abline(slope = 1, intercept = 0) + coord_fixed(xlim = scale, ylim = scale) +
  xlab("From the LC method") + ylab("From the run's log") +
  ggtitle("Methods for estimating instrument time") +
  facet_wrap(~byMonth) + scale_color_viridis(option = "A")
#proteoCraft::poplot(plot)
#print(plot)
plotLy <- ggplotly(plot)
saveWidget(plotLy, paste0(wd, "/Bruker run times.html"))
system(paste0("open \"", wd, "/Bruker run times.html"))

# Handle times
require(lubridate)
# Now, at some point, some sod decided to change the date format!!!
dDrsTbl$Date_POSIXct <- parse_date_time(dDrsTbl$Date, "%m-%d-%Y %I:%M:%S %p")
dDrsTbl$End_POSIXct <- parse_date_time(dDrsTbl$End, "%m-%d-%Y %I:%M:%S %p")
w <- which(is.na(dDrsTbl$Date_POSIXct)|is.na(dDrsTbl$End_POSIXct))
if (length(w)) {
  dDrsTbl$Date_POSIXct[w] <- parse_date_time(dDrsTbl$Date[w], "%d-%m-%Y %I:%M:%S %p")
  dDrsTbl$End_POSIXct[w] <- parse_date_time(dDrsTbl$End[w], "%d-%m-%Y %I:%M:%S %p")
}

# Also use this script to check which folders in the timsTOF still contain .d files
# (useful for cleaning up)

# Filter by whether these were runs for clients or not
dDrsTbl$File_name <- gsub(".*/", "", dDrsTbl$Raw.file)
dDrsTbl$blank <- c("", "+")[grepl("blank|zig|zag|blnk", dDrsTbl$Raw.file)+1]
pat1 <- paste0("^(", paste(paste0("(", vapply(targDirs,
                                              topattern, "", start = FALSE), ")"), collapse = "|"), ")/")
pat2 <- paste0(pat1, "(([^/]+ ?gr(ou)?p)|(ext_)|(Solgate)|(Valanx)|(LTuriak)|(PCF)|(Vicoso))")

dDrsTbl$ForClient <- grepl(pat2, dDrsTbl$Raw.file)

Years <- sort(unique(dDrsTbl$Year))
Years <- Years[which(!is.na(Years))]
dflt <- c(thisYear-1, max(Years))
dflt <- dflt[which(dflt %in% Years)][1]
Year <- as.integer(dlg_list(Years, dflt, title = "Year to extract")$res)
aggregate(dDrsTbl$Date_POSIXct, list(dDrsTbl$Source), max, na.rm = TRUE)
dDrsTbl2 <- dDrsTbl[which(dDrsTbl$Year == thisYear-1),]
unique(dDrsTbl2$Year)
ClientProjects <- sort(unique(gsub("/[^/]+\\.d$", "",
                                   gsub(pat1, "", dDrsTbl2$Raw.file[which(dDrsTbl2$ForClient)]))))
ClientFolders <- sort(unique(gsub("/.+", "", ClientProjects)))
OtherFolders <- sort(unique(gsub("/[^/]+\\.d$", "",
                                 gsub(pat1, "", dDrsTbl2$Raw.file[which(!dDrsTbl2$ForClient)]))))

ClientProjects
ClientFolders
# To check that we do not miss any of the client projects:
OtherFolders
Blanks <- dDrsTbl2$blank == "+"
AllClientSamples <- sum(dDrsTbl2$ForClient)
ClientBlanks <- sum((dDrsTbl2$ForClient)&(dDrsTbl2$blank == "+"))
ClientSamples <- AllClientSamples-ClientBlanks
AvClienRunTime <- sum(dDrsTbl2$`Run length (min)`*dDrsTbl2$ForClient*(!Blanks))/(ClientSamples)
BatchSize <- (ClientSamples)/length(ClientProjects)

recStart <- min(dDrsTbl2$End_POSIXct, na.rm = TRUE)
recEnd <- max(dDrsTbl2$End_POSIXct, na.rm = TRUE)
recStart
recEnd
Days <- recEnd - recStart
print(paste0("Time covered by the record: ", round(Days, 2), " days"))
Days <- as.numeric(Days)
ttl <- paste0("timsTOF HT instrument time use report - ", Sys.Date())
if ((length(Year) == 1)&&(!is.na(Year))) {
  ttl <- paste0(Year, " ", ttl)
}
FracRunTime <- (sum(dDrsTbl2$`Run length (min)`)/(24*60))/Days
FracClients2TotTime <- (sum(dDrsTbl2$ForClient*dDrsTbl2$`Run length (min)`)/(24*60))/Days
FracClients2RunTime <- sum(dDrsTbl2$ForClient*as.numeric(dDrsTbl2$`Run length (min)`))/sum(as.numeric(dDrsTbl2$`Run length (min)`))
Report <- c(ttl,
            paste(rep("-", nchar(ttl)), collapse = ""),
            "",
            paste0("Tot. number of client samples (including blanks): ", AllClientSamples),
            paste0("Tot. number of client samples (not including blanks): ", ClientSamples),
            paste0("Ratio client samples/blanks: ", round(ClientSamples/ClientBlanks, 2)),
            paste0("Av. client run time: ", round(AvClienRunTime/60, 2), " h"),
            paste0("Batches in ", Year, ": ", length(ClientProjects)),
            paste0("Av. batch (project) size: ", round(BatchSize, 0), " samples"),
            paste0("Time spent running samples: ", as.character(round(sum(dDrsTbl2$`Run length (min)`)/60)), " h = ", round(100*FracRunTime, 1), "% of record scope"),
            paste0("Time spent running samples for clients: ", as.character(round(sum(dDrsTbl2$`Run length (min)`*dDrsTbl2$ForClient)/60)), " h = ", round(100*FracClients2TotTime, 1), "% of record scope"),
            paste0(round(100*FracClients2RunTime, 1), "% of actual run time was spent running samples for clients"))
Report
Report <- c(Report, "", "Projects:", ClientProjects, "")
cat(paste0(Report, "\n"))
write(Report, paste0(wd, "/", ttl, ".txt"))
#openwd()
