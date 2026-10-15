#' load_Context
#'
#' @description
#' This function takes a dataset details record file (as created when we run the Start_analysis.R app), then it reloads in environment:\cr
#'  - dataset name\cr
#'  - input and output directory\cr
#'  - temporary work directory\cr
#' and if the latter exists, sets work directory to it.
#' 
#' @param record The dataset details record .txt file to parse.
#' @param startDir In case record is not provided, you can specify here the directory in which the selection window will start. Ignored otherwise.
#' @param clean If set to TRUE, will first remove all objects existing in the current environment. 
#'
#' @export

load_Context <- function(record,
                         startDir,
                         clean = FALSE) {
  TESTING <- FALSE
  #DefArg(load_Bckp); TESTING <- TRUE
  if (clean) { rm(list = setdiff(ls(), c("TESTING", "clean")), envir = .GlobalEnv) }
  misFun <- if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    \(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  # Cleanup workspace here
  #
  if (misFun(record)) {
    if (misFun(startDir)) {
      if (exists("wd") && is.character(wd) && (length(wd) == 1L) && dir.exists(wd)) {
        defltdir <- wd
      } else {
        homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
        dlft <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
        defltdir <- dlft$Path[match("Temporary folder", dlft$Folder)]
      }
    } else { defltdir <- startDir }
    #record <- normalizePath(choose.files(paste0(defltdir, "/*.txt"), multi = FALSE), winslash = "/")
    record <- rstudioapi::selectFile("Select dataset details record .txt file to parse...",
                                     path = paste0(defltdir, "/*.txt"),
                                     filter = "Text file (*.txt)")
  }
  record <- readr::read_lines(record)
  assign("dtstNm", sub("^ *-> *", "", record[grep("^Dataset(/project)? name:", record)+1L]), envir = .GlobalEnv)
  gI1 <- grep("^Input director(y|(ies)):", record)+1L
  gI2 <- grep("^ *-> *", record, invert = TRUE)
  gI2 <- min(gI2[gI2 > gI1])-1L
  assign("inDirs", sub("^ *-> *", "", record[gI1:gI2]), envir = .GlobalEnv)
  gO <- grep("^Final output directory:", record)
  if (length(gO)) {
    assign("outdir", sub("^ *-> *", "", record[gO+1L]), envir = .GlobalEnv)
  }
  assign("wd", sub("^ *-> *", "", record[grep("^Temporary work directory:", record)+1L]), envir = .GlobalEnv)
  if (!dir.exists(wd)) { dir.create(wd, recursive = TRUE) }
  lapply(inDirs, \(indir) {
    if (!dir.exists(indir)) { warning(paste0("Input directory \"", indir, "\" does not exist!")) }
  })
  setwd(wd)
  assign(".obj", c(".obj", intersect(c("dtstNm", "inDirs", "outdir", "wd"), ls(envir = .GlobalEnv))),
         envir = .GlobalEnv)
}
