#' .pyConfig
#' 
#' A function to automatically check our python configuration, upgrade pip, install sdrf-pipelines and check that python scripts is in path.
#' 
#' @export

.pyConfig <- function() {
  # Update pip
  cmd <- "pip install --upgrade pip"
  rs <- system(cmd, show.output.on.console = FALSE)
  if (rs == 1L) { return(FALSE) }
  # Install sdrf-pipelines
  cmd <- "pip install sdrf-pipelines"
  rs <- system(cmd, show.output.on.console = FALSE)
  # Install sdrf-pipelines with ontologies
  cmd <- "pip install sdrf-pipelines[ontology]"
  rs <- system(cmd, show.output.on.console = FALSE)
  #
  # Configure PATH
  py_path <- normalizePath(system("where python", intern = TRUE), winslash = "/")
  py_path <- grep("/Python/", py_path, value = TRUE)
  py_list <- c()
  checkSDRF <- FALSE
  if (length(py_path)) {
    py_path <- gsub("/Python/.*", "/Python/", py_path)
    py_list <- lapply(py_path, \(x) {
      x <- list.dirs(x, full.names = TRUE, recursive = TRUE)
      x <- grep("/Scripts$", x, value = TRUE)
      return(x)
    })
    py_list <- py_list[which(lengths(py_list) >= 1L)]
  }
  if (!length(py_list)) { return(FALSE) }
  all_pyScrptPaths <- unlist(py_list)
  if (!length(all_pyScrptPaths)) {
    warning("Could not find .../Scripts/ directory!")
    return(FALSE)
  }
  tmp <- unique(unlist(strsplit(Sys.getenv("PATH"), ";")))
  tmp <- union(tmp, normalizePath(all_pyScrptPaths, winslash = "\\"))
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
  return(TRUE)
}
