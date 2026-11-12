#' .safe_listFls
#' 
#' @description 
#' A safer function to list all files (or directories) from a path.\cr
#' \cr
#' NB:\cr
#' ===\cr
#' This has a bit of a history...\cr
#' At some stage, our workflows read MS methods from raw spectra files. Sometimes, however, these, or the search folder, have been moved, so the files must be located.\cr
#' We could prompt the user for their location (we still do if all else fails), but we aim to interrupt the workflow to ask for user input as rarely as possible, especially since this is nice-to-have (used for materials and methods), not essential. So, the scripts first attempt to find the files automatically by scanning expected locations...\cr
#' And it worked well enough, until SOMEONE decided to run a MeroX search on a couple of Bruker .d files, a search which got stuck at 7% but still managed to generate 50 odd gigabytes of tens of millions of small files...\cr
#' ...\cr
#'\cr
#' Stress test welcome >:D\cr
#'\cr
#' So.\cr
#' We needed a safe way to look for files/folders, but without spending a lifetime looking.\cr
#' Enter, this function:\cr
#' (courtesy of chatGPT, because that's really not the type of commands I am familiar with)\cr
#' \cr
#' 
#' @param path Path from which to list files.
#' @param type Type of things tp list. Default = "files" but could also be "dir" (directories).
#' @param timeout Default = 10
#' 
#' @returns
#' Returns the .
#' 
#' @export

.safe_listFls <- function(path,
                          type = "files",
                          timeout = 10L) {
  type <- substring(tolower(type), 1L, 3L)
  stopifnot(type %in% c("fil", "dir"))
  sysNm <- Sys.info()[['sysname']]
  start <- Sys.time()
  out <- character()
  if (type == "fil") {
    if (sysNm == "Windows") {
      cmd <- sprintf('Get-ChildItem -File -Recurse -Path %s | ForEach-Object { $_.FullName }', shQuote(path))
      p <- processx::process$new("powershell",
                                 c("-Command", cmd),
                                 stdout = "|",
                                 stderr = "|")
    }
    if (sysNm %in% c("Linux", "Darwin")) { # Version for Unix-like OS, not tested!
      p <- processx::process$new("find",
                                 c(path, "-type", "f"),
                                 stdout = "|", stderr = "|")
    }
  }
  if (type == "dir") {
    if (sysNm == "Windows") {
      cmd <- sprintf('Get-ChildItem -Directory -Recurse -Path %s | ForEach-Object { $_.FullName }', shQuote(path))
      p <- processx::process$new("powershell",
                                 c("-Command", cmd),
                                 stdout = "|",
                                 stderr = "|")
    }
    if (sysNm %in% c("Linux", "Darwin")) { # Version for Unix-like OS, not tested!
      p <- processx::process$new("find",
                                 c(path, "-type", "d"),
                                 stdout = "|",
                                 stderr = "|")
    }
  }
  repeat {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    if (elapsed >= timeout) {
      p$kill()
      break
    }
    poll <- p$poll_io(200)
    if (!is.null(poll[["output"]]) && poll[["output"]] == "ready") {
      chunk <- p$read_output_lines()
      out <- c(out, chunk)
    }
    if (!p$is_alive()) {
      # drain remaining output
      # out <- c(out, p$read_output_lines())
      break
    }
  }
  out <- gsub("\\\\", "/", out)
  return(out)
}
