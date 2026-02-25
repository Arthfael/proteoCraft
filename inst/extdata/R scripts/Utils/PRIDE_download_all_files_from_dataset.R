################################################################################
# A script to automate downloading all of the (raw) files from a PRIDE dataset #
################################################################################
# Initialization
options(stringsAsFactors = FALSE)
packs <- c("devtools", "RCurl", "svDialogs", "jsonlite")
for (pack in packs) {
  if (!require(pack, character.only = TRUE)) { install.packages(pack) }
  require(pack, character.only = TRUE)
}
suppressMessages(install_github("PRIDE-R/prideR"))
library(prideR)
# Modified prideR::get.ProjectSummary function

get.ProjectSummary2 <- function (accession) {
  fromJSON(file = paste0("http://www.ebi.ac.uk/pride/ws/archive/projects/", accession), method = "C")
}

# Choose dataset
dataset <- dlg_input("Enter the PRIDE accession of the target dataset", "PXD######")$res
#dataset <- "PXD028312"

# Choose file types of interest
FileTypes <- c("all", "raw", "d", "mzML", "mzXML", "txt", "tsv", "csv", "zip", "rar")
filetypes <- dlg_list(FileTypes, 1, TRUE, "Select file type(s) of interest:")$res
if ("all" %in% filetypes) { filetypes <- "all" }
if ("d" %in% filetypes) { filetypes <- c(filetypes, "d.zip", "d.rar") }
if (exists("dataset")) {
  tst <- suppressWarnings(try(prideR::get.ProjectSummary(dataset), silent = TRUE))
  if (class(tst) == "try-error") {
    cat("Sounds like the EBI's own function still cannot access their own website...
Sigh...
There, lemme fix that for you...
")
    tst <- try(get.ProjectSummary2(dataset), silent = TRUE)
  }
  if (class(tst) != "try-error") {
    dir <- paste0("...Destination/", dataset)
    if (!dir.exists(dir)){ dir.create(dir, recursive = TRUE) }
    setwd(dir)
    ftp <- gsub("^ftp://", "https://", tst[["_links"]]$datasetFtpUrl$href)
    tst2 <- readLines(url(ftp))
    files <- grep("<tr><td valign=\"top\"><img src=\"/icons/[^\"]+\\.gif\" alt=\"\\[.+\\]\"></td><td><a href=\"", tst2, value = TRUE)
    files <- gsub("<.+", "", gsub(".+<a href[^>]+>", "", files))
    if ((length(filetypes) > 1)||(filetypes != "all")) {
      files <- unlist(sapply(filetypes, function(filetype) { grep(paste0(gsub("\\.", "\\.", paste0(".", filetype)), "$"), files, value = TRUE) }))
    }
    # TO DO: also parse tst2 for combined file sizes, and compare vs available disc space
    FlsPths <- paste0(ftp, "/", files)
    if (length(FlsPths)) {
      cat(paste0("Downloading ", length(FlsPths), " files..."))
      for (FlPth in FlsPths) { #FlPth <- FlsPths[1]
        destFl <- paste0(dir, "/", basename(FlPth))
        f <- CFILE(destFl, mode = "wb")
        curlPerform(url = FlPth, writedata = f@ref)
        close(f)
      }
      cat(paste0("All your files have been successfully downloaded to \"", dir, "\"\n"))
    } else { warning("No files found, nothing to download!") }
  } else { cat("...
and sounds like my fix too is broken... sorry!") }
}
unload("prideR") # Critical because this package currently corrupts the data.frame method for View(),
# meaning if left loaded it is now impossible to open data.frames in this session!

