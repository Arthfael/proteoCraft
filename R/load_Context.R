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
  
  if (clean) { rm(list = ls(), envir = .GlobalEnv) }
  library(proteoCraft)
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::load_Bckp)
  #TESTING <- TRUE
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  # Cleanup workspace here
  #
  if (misFun(record)) {
    if (misFun(startDir)) {
      if ((exists("wd"))&&("character" %in% class(wd))&&(dir.exists(wd))) {
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
  record <- readLines(record)
  dtstNm <<- gsub("^ *-> *", "", record[grep("^Dataset name:", record)+1])
  inDirs <<- gsub("^ *-> *", "", record[grep("^Input directory:", record)+1])
  outdir <<- gsub("^ *-> *", "", record[grep("^Final output directory:", record)+1])
  wd <<- gsub("^ *-> *", "", record[grep("^Temporary work directory:", record)+1])
  if (!dir.exists(wd)) { dir.create(wd, recursive = TRUE) }
  lapply(inDirs, function(indir) {
    if (!dir.exists(indir)) { stop(paste0("Input directory \"", indir, "\" does not exist!")) }
  })
  setwd(wd)
  .obj <<- c(".obj", "dtstNm", "inDirs", "outdir", "wd")
}
