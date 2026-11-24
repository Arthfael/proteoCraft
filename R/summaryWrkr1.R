#' .summaryWrkr1
#' 
#' @description
#' Worker function used by MQ.summary() to extract a raw file's TIC or BPC.\cr
#' Can now process:\cr
#'  - Thermo .raw files\cr
#'  - Bruker .d folders
#'
#' @param x File index (a row of the table which is passed to the "raw" argument).
#' @param raw Table of MS raw files. Must contain a column named "Path" comtaining a absolute raw file paths.
#' @param type Type of chromatogram to export. Can be "tic" or "bpc" (MS1-level), or, for Bruker files, "pasef_tic" (also just "pasef") or "pasef_bpc".

.summaryWrkr1 <- function(x,
                          raw,
                          type) {
  xt <- tolower(sub(".*\\.", "", raw$Path[x]))
  type <- tolower(type)
  type[type == "pasef"] <- "pasef_tic"
  isPasef <- grepl("^pasef", type)   # logical
  type <- sub("^pasef_", "", type)   # = tic or bpc
  if (!type %in% c("bpc", "tic")) {
    warning("Unrecognized \"type\" argument!")
    return()
  }
  if ((xt == "raw") && isPasef) { return() }
  rawFun <- \(fl,
              tp = type) {
    if (!file.exists(fl)) { return() }
    res <- rawrr::readChromatogram(fl,
                                   type = tp)
    return(data.frame("Retention time" = as.numeric(res$times),
                      "Intensity" = as.numeric(res$intensities),
                      check.names = FALSE))
  }
  dFun <- \(dr,
            tp = type,
            pasef = isPasef) {
    analysisTDF_fl <- paste0(dr, "/analysis.tdf")
    if (!file.exists(analysisTDF_fl)) { return() }
    con <- RSQLite::dbConnect(RSQLite::SQLite(), analysisTDF_fl)
    frames <- RSQLite::dbGetQuery(con, "
  SELECT Id, Time, MsMsType, SummedIntensities, MaxIntensity
  FROM Frames")
    RSQLite::dbDisconnect(con)
    # According to Claude.AI, possible MsMsType values are:
    #   Value   Meaning
    #   0       MS1 (survey scan)
    #   2       MRM
    #   8	      ddaPASEF MS/MS
    #   9	      diaPASEF MS/MS
    #   10	    prm-PASEF MS/
    msType <- list(standard = 0L,
                   pasef = c(8L, 9L))[[pasef+1L]]
    msCol <- c(tic = "SummedIntensities",
               bpc = "MaxIntensity")[tp]
    res <- subset(frames,
                  MsMsType %in% msType)
    return(data.frame(# Time is in seconds, we want minutes
      "Retention time" = res$Time/60,
      # Intensities are stored as integer - which could lead to overflow if we perform any computations on them
      # (e.g. if we normalize and multiply before dividing)
      # -> convert to numeric!
      "Intensity" = as.numeric(res[[msCol]]),
      check.names = FALSE))
  }
  RES <- switch(xt,
                raw = rawFun(raw$Path[x]),
                d = dFun(raw$Path[x]))
  if (!is.data.frame(RES)) { return() }
  RES$"Raw file path" <- raw$Path[x]
  RES$"Raw file" <- raw$`Raw file`[x]
  return(RES)
}
