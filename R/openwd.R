#' openwd
#'
#' @description 
#' A function to open your a directory (work directory by default). Only works on PCs.
#' @param dir Directory to open. If missing, defaults to current work directory.
#'
#' @examples
#' openwd() # Open current work directory

#' 
#' @export

openwd <- function(dir) {
  syst <- Sys.info()
  syst <- tolower(syst[match("sysname", names(syst))])
  if (missing("dir")) { dir <- getwd() }
  dir <- suppressWarnings(normalizePath(dir))
  if (!dir.exists(dir)) {
    cat(paste0("Directory does not exist:\n", dir, "\n"))
  } else {
    options(warn = -1)
    if (syst == "windows") { cmd <- paste0("explorer \"", dir, "\"") } else {
      cmd <- paste0("xdg-open ", dir) # Not tested
    }
    tst <- try(shell(cmd, intern = TRUE), silent = TRUE)
    if (length(tst)) {
      if (syst == "windows") {
        cat("Hmmm... something went wrong, I couldn't open this directory...\n")
      } else {
        cat("Sorry, this function currently only works on a PC!\n")
        # (There used to be a poor joke about the "Steve Jobs cult" here, but... well he's been dead for years now, so I think I can let this go.)
      }
    }
    options(warn = 0)
  }
}
