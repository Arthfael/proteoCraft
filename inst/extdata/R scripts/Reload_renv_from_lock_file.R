# Re-load a project from a renv lock file
#
# This is meant for analysis reproducibility, to restore an older R environment. Before running this script:
#  - Check that you have the correct version of R selected in RStudio
#    (comparing against the version recorded in the sessionInfo.txt file).
#  - Create the folder you will want to work from.
# Then run this script and select the folder you just created.
# If this was not already a project folder, a project will be initiated and the R session will restart.
# Browse to the lock file you want to reload and select it, the old renv will be reloaded.

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#                                       |        |                                       #
#                                       |        |                                       #
#                         This script will restart the R session,                        #
#                     which will interrupt the sequence of commands.                     #
#                                       |        |                                       #
#                                       |        |                                       #
#                                      \          /                                      #
#                                       \        /                                       #
#                                        \      /                                        #
#                                         \    /                                         #
#                                          \  /                                          #
#                                           \/                                           #
#                                  Run it line by line!                                  #
#                                                                                        #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

stopifnot(interactive())
library(renv)
wd <- rstudioapi::selectDirectory("Select project directory where you will run the analysis (containing or not the renv.lock file to restore)")
#if (!dir.exists(wd)) { dir.create(wd, recursive = TRUE) }
setwd(wd)
#proteoCraft::openwd()
nuEnv <- (!file.exists("renv.lock"))
if (nuEnv) { renv::init() }
lockFl <- rstudioapi::selectFile("Select renv lock file to restore", path = wd, filter = "renv lock files (*.lock)")
docsPath <- normalizePath(Sys.getenv("HOME"), winslash = "/")
if (grepl("^~", lockFl)) { lockFl <- gsub("^~", docsPath, lockFl) }
#
# Check that the R version you are running is the same as that in the lock file!
tmp1 <- readLines(lockFl)
tmp1 <- tmp1[grep("\"R\": *\\{", tmp1)+1]
tmp1 <- gsub(" *\"Version\": *\"|\" *,.*", "", tmp1)
tmp2 <- version
tmp2 <- paste0(tmp2$major, ".", tmp2$minor)
if (tmp1 != tmp2) {
  stop("This renv lock file uses ", tmp1, " but you are currently running ", tmp2,
       "!\nYou first need to switch to the correct version of R then restart RStudio before moving ahead!")
}
#
# Note on exclusions:
# - pak should be an exception if we use it
# - rawrr can cause some issues and it's my fault: at some point, I used 1.11.14 from an older one buried deep inside Bioconductor;
#   pak cannot find it, so we instead install directly from there.
pakOpt <- getOption("renv.config.pak.enabled")
excl <- c()
if (pakOpt) { excl <- union("pak", excl) }
library(evaluate)
kount <- 0
renvXprs <- "kount <- kount + 1; cat(\" - Attempt #\", kount, \"\n\"); renv::restore(lockfile = lockFl, prompt = FALSE, exclude = excl)"
#eval(parse(text = renvXprs))
myEval <- evaluate::evaluate(renvXprs, stop_on_error = 0L)
tstEval <- function(x = myEval) { sum(vapply(x, function(y) { sum(c("error", "simpleError") %in% class(y)) > 0 }, TRUE)) }
whError <- function(x = myEval) { which(vapply(x, function(y) { sum(c("error", "simpleError") %in% class(y)) > 0 }, TRUE)) }
if ((tstEval())&&(pakOpt)) {
  options(renv.config.pak.enabled = FALSE)
  excl <- setdiff(excl, "pak")
  myEval <- evaluate::evaluate(renvXprs, stop_on_error = 0L)
}
# Below not tested
if (tstEval()) {
  tmp <- unlist(strsplit(as.character(myEval[[whError()]]), "\'"))[[2]] # Not tested
  if (tmp == "rawrr 1.11.14") {
    url <- "https://bioconductor.org/packages/3.19/bioc/src/contrib/Archive/rawrr/rawrr_1.11.14.tar.gz"
    require(curl)
    dstfl <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/rawrr_1.11.14.tar.gz")
    curl_download(url, dstfl)
    install.packages(dstfl)
    unlink(dstfl)
    excl <- union("rawrr", excl)
    myEval <- evaluate::evaluate(renvXprs, stop_on_error = 0L)
  }
}
# If this fails, you may have it an issue with the moving Bioconductor support window...
# Basically, likely the failing package version was installed from a Bioconductor version older than R's current one.
if (tstEval()) {
  out <- Filter(is.character, myEval)
  out <- unlist(strsplit(unlist(out), "\n"))
  tmp <- unique(grep("ERROR", out, ignore.case = TRUE, value = TRUE))
  tmp <- gsub(" *- *Downloading *| *from.*", "", grep(" *- *Downloading *", tmp, value = TRUE))
  tmp <- setdiff(tmp, excl)
  while (tstEval()&&(length(tmp))) {
    excl <- union(tmp, excl)
    myEval <- evaluate::evaluate(renvXprs, stop_on_error = 0L)
    out <- Filter(is.character, myEval)
    out <- unlist(strsplit(unlist(out), "\n"))
    tmp <- unique(grep("ERROR", out, ignore.case = TRUE, value = TRUE))
    tmp <- gsub(" *- *Downloading *| *from.*", "", grep(" *- *Downloading *", tmp, value = TRUE))
    tmp <- setdiff(tmp, excl)
  }
}
if (!tstEval()) {
  if (length(excl)) {
    warning(paste0("The following packages were not reloaded correctly and will need a manual install (or someone to hack the lock file):\n",
                   paste(paste0(" - ", excl), collapse = "\n"), "\n"))
  }
} else {
  stop("Your lock file could not be reloaded automatically, probably because of Bioconductor version shenanigans...\nWhat to do next? Open the renv file, check which packages and versions are missing, locate them (good luck) then install them manually.\n\nHAVE\n...\nFUN!")
}
options(renv.config.pak.enabled = pakOpt)
renv::status()
renv::project()
require(proteoCraft)
renv %<o% TRUE
