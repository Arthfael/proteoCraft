#' .pyConfig
#' 
#' A function to automatically check our python configuration, upgrade pip, install sdrf-pipelines and check that python scripts is in path.
#' 
#' @param pathWarn Logical, default = FALSE. If TRUE, and it finds in the output of Sys.getenv("PATH") some paths which do not exist, this will issue a warning. However, this may trigger false alerts - it seems that some versions of R/RStudio have paths there which have since slightly changed - without breaking functionality.
#' 
#' @export

.pyConfig <- function(pathWarn = FALSE) {
  tst <- try({
    rs <- system2("py", args = "list", stdout = TRUE)
    pyVersions <- sub("[ [].*", "", rs[2L:length(rs)])
    k <- 0
    tst <- FALSE
    while ((!tst) && (k < length(pyVersions))) {
      k <- k+1L
      pyVers <- pyVersions[k]
      tst <- try({
        # Update pip
        cmd <- paste0("py -", pyVers, " -m pip install --upgrade pip")
        rs <- system(cmd, show.output.on.console = FALSE)
        if (rs != 0L) { return(FALSE) }
        # Install sdrf-pipelines
        cmd <- paste0("py -", pyVers, " -m pip install sdrf-pipelines")
        rs <- system(cmd, show.output.on.console = FALSE)
        if (rs != 0L) { stop("Could not install sdrf-pipelines, investigate!!!") }
        # Install sdrf-pipelines with ontologies
        cmd <- paste0(cmd, "[ontology]")
        rs <- system(cmd, show.output.on.console = FALSE)
        if (rs != 0L) { stop("Could not install sdrf-pipelines[ontology], investigate!!!") }
        # NB: The SDRF pipelines are necessary for checking that edited SDRFs files are correct.
        # If installation fails, try instead the following commands:
        #      "py -3.14 -m pip install sdrf-pipelines"
        #      "py -3.14 -m pip install sdrf-pipelines[ontology]"
        # (replacing 3.14 with the desired version of Python.)
        # This allows you to bypass issues with a too recent version of Python not having yet an update for
        #
      }, silent = TRUE)
      tst <- !inherits(tst, "try-error")
    }
    if (!tst) { stop() }
    #
    # Configure PATH
    py_path <- normalizePath(system("where python", intern = TRUE), winslash = "/")
    py_path <- grep("/Python/", py_path, value = TRUE)
    py_list <- c()
    checkSDRF <- FALSE
    if (length(py_path)) {
      py_path <- sub("/Python/.*", "/Python/", py_path)
      py_list <- lapply(py_path, \(x) { #x <- py_path[1L]
        x <- list.dirs(x, full.names = TRUE, recursive = TRUE)
        x <- grep("/Scripts$", x, value = TRUE)
        x <- grep(paste0("-", pyVers, "[-/]"), x, value = TRUE)
        return(x)
      })
      py_list <- py_list[lengths(py_list) >= 1L]
    }
    if (!length(py_list)) { stop() }
    all_pyScrptPaths <- unlist(py_list)
    if (!length(all_pyScrptPaths)) {
      warning("Could not find .../Scripts/ directory!")
      return(FALSE)
    }
    tmp <- unique(unlist(strsplit(Sys.getenv("PATH"), ";"))) # Note: this is the value created when R started and may include additional stuff! Not exactly the same as actually reading the PATH environment variable!
    tmp <- base::union(tmp, normalizePath(all_pyScrptPaths, winslash = "\\"))
    w <- which(!file.exists(tmp))
    lW <- length(w)
    if (lW) {
      warning(paste0("The following path", c("", "s")[(lW>1L)+1L],
                     " in your PATH environmental variable do", c("es", "")[(lW>1L)+1L], " not exist:\n",
                     paste0(" -> ", tmp[w], "\n", collapse = "")))
    }
    Sys.setenv(PATH = paste(tmp, collapse = ";"))
    #
    pyParser <- unlist(lapply(all_pyScrptPaths, list.files, pattern = "parse_sdrf\\.exe$", recursive = TRUE, full.names = TRUE))
    checkSDRF <- length(pyParser) > 0L
  }, silent = TRUE)
  res <- if (inherits(tst, "try-error")) {
    list(FALSE)
  } else {
    list(out = TRUE,
         vers = pyVers,
         parse_sdrf = pyParser)
  }
  return(res)
}
