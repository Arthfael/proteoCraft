#' .pyConfig
#'
#' Finds a usable Python (>= minVersion), creates a dedicated virtual environment for it (under tools::R_user_dir("proteoCraft", "data")), installs sdrf-pipelines in it, and checks that the parse_sdrf command works.\cr
#' Nothing is installed into the user's own Python.
#' 
#' This new version is 99% the work of Claude.AI and should be a vast improvement over my previous version.
#'
#' @param minVersion,maxVersion Character, allowed Python versions (maxVersion = NULL: no limit).
#' @param venvRoot Folder in which to create the python virtual environment.
#' @param venvPrefix Default = "py", prefix for the name of the virtual environment. Will be supplemented with the python version, e.g. "py3.14"
#' @param ontology Logical, default = TRUE. Install the sdrf-pipelines[ontology] extra.
#' @param addToPath Logical, default = TRUE. Append the venv's script directory to PATH for this R session.\cr
#'  (The returned `parse_sdrf` full path works regardless.)
#' @param reinstall Logical, default = FALSE. If TRUE, rebuild the venv even if a working one exists.
#' @param pathWarn Logical, default = FALSE. Warn about non-existent entries in PATH (may give false alerts).
#' @param verbose Logical, default = TRUE.
#'
#' @return A list with elements out (logical), vers, python (venv interpreter), venv,
#'  parse_sdrf (full path of the executable) and log (what happened). On failure out = FALSE.
#'
#' @export

.pyConfig <- function(minVersion = "3.9",
                      maxVersion = NULL,
                      venvRoot = NULL,
                      venvPrefix = "py",
                      ontology = TRUE,
                      addToPath = TRUE,
                      reinstall = FALSE,
                      pathWarn = FALSE,
                      verbose = TRUE) {
  isWin <- .Platform$OS.type == "windows"
  log <- character()
  say <- \(...) {
    m <- paste0(...)
    log <<- c(log, m)
    if (verbose) { message("[.pyConfig] ", m) }
  }
  #
  # Run a command, never throw, always return list(status, output)
  run <- \(cmd,
           args = character(),
           timeout = 0L) {
    out <- suppressWarnings(try(system2(cmd, args, stdout = TRUE, stderr = TRUE,
                                        timeout = timeout), silent = TRUE))
    if (inherits(out, "try-error")) { return(list(status = 127L, output = as.character(out))) }
    st <- attr(out, "status")
    list(status = if (is.null(st)) { 0L } else { as.integer(st) },
         output = as.character(out))
  }
  tail3 <- \(x) { paste(utils::tail(x, 3L), collapse = " | ") }
  #
  # Venv layout helpers
  venvBin <- \(v) { gsub("\\\\", "/", file.path(v, if (isWin) { "Scripts" } else { "bin" })) }
  venvPy <- \(v) { gsub("\\\\", "/", file.path(venvBin(v), if (isWin) { "python.exe" } else { "python" })) }
  venvParser <- \(v) { gsub("\\\\", "/", file.path(venvBin(v), if (isWin) { "parse_sdrf.exe" } else { "parse_sdrf" })) }
  parserOK <- \(v) {
    file.exists(venvParser(v)) && (run(venvParser(v), "--help", timeout = 120L)$status == 0L)
  }
  #
  finish <- \(v, vers) {
    if (addToPath) {
      sep <- .Platform$path.sep
      cur <- strsplit(Sys.getenv("PATH"), sep, fixed = TRUE)[[1L]]
      bin <- venvBin(v)
      nrm <- \(x) normalizePath(x, winslash = "/", mustWork = FALSE)
      if (pathWarn) {
        w <- which(!file.exists(cur) & nzchar(cur))
        if (length(w)) {
          warning("The following PATH entries do not exist:\n", paste0(" -> ", cur[w], "\n", collapse = ""))
        }
      }
      if (!(nrm(bin) %in% nrm(cur))) { Sys.setenv(PATH = paste(c(cur, bin), collapse = sep)) }
    }
    list(out = TRUE,
         vers = vers,
         python = venvPy(v),
         venv = gsub("\\\\", "/", v),
         parse_sdrf = venvParser(v),
         log = log)
  }
  fail <- \(why) {
    say(why)
    warning(why, call. = FALSE)
    list(out = FALSE, vers = NA_character_, python = NA_character_, venv = NA_character_,
         parse_sdrf = character(0L), log = log)
  }
  #
  okRoot <- (!missing(venvRoot)) && (length(venvRoot) == 1L) && (!is.na(venvRoot)) && is.character(venvRoot)
  if (!okRoot) {
    okRoot <- FALSE
    homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
    locDirsFl2 <- paste0(homePath, "/Default_locations.xlsx")
    if (file.exists(locDirsFl2)) {
      locDirs <- openxlsx2::read_xlsx(locDirsFl2)
      venvRoot <- locDirs$Path[match("Temporary folder", locDirs$Folder)]
      okRoot <- dir.exists(venvRoot)
    }
  } 
  if (!okRoot) {
    venvRoot <- tools::R_user_dir("proteoCraft", "data")
  }
  if (!dir.exists(venvRoot)) {
    dir.create(venvRoot, recursive = TRUE, showWarnings = FALSE)
  }
  #
  # 1) Reuse an existing working venv
  if ((!reinstall) && dir.exists(venvRoot)) {
    vs <- list.dirs(venvRoot, recursive = FALSE, full.names = TRUE)
    vs <- vs[grepl(paste0(topattern(venvPrefix), "[0-9]+\\.[0-9]+$"), basename(vs))]
    if (length(vs)) {
      vers <- sub(topattern(venvPrefix), "", basename(vs))
      vs <- vs[order(numeric_version(vers), decreasing = TRUE)]
      for (v in vs) {
        if (parserOK(v)) {
          say("Reusing existing environment: ", gsub("\\\\", "/", v))
          return(finish(v, sub(topattern(venvPrefix), "", basename(v))))
        }
      }
    }
  }
  #
  # 2) Collect candidate interpreters
  cands <- character()
  ev <- Sys.getenv("RETICULATE_PYTHON")
  if (nzchar(ev)) { cands <- c(cands, ev) }
  if (isWin) {
    # Works with both the legacy py.exe launcher and the new Python install manager
    r <- run("py", "-0p", timeout = 30)
    if (r$status == 0L) {
      m <- regmatches(r$output, regexpr("[A-Za-z]:[\\\\/].*\\.exe", r$output))
      cands <- c(cands, trimws(m))
    }
  }
  cands <- c(cands, unname(Sys.which(c("python3", paste0("python3.", 9:15), "python"))))
  if (isWin) {
    # a) py launcher (legacy py.exe and the new Python install manager both understand -0p)
    r <- run("py", "-0p", timeout = 30)
    if (r$status == 0L) {
      m <- regmatches(r$output, regexpr("[A-Za-z]:[\\\\/].*\\.exe", r$output))
      cands <- c(cands, trimws(m))
    }
    # b) Windows registry (PEP 514): every standard installer registers itself here,
    #    independently of whether the py launcher exists
    for (hive in c("HCU", "HLM")) {
      reg <- try(utils::readRegistry("Software\\Python", hive = hive, maxdepth = 4L), silent = TRUE)
      if (inherits(reg, "try-error") || !length(reg)) { next }
      u <- unlist(reg)
      cands <- c(cands,
                 u[grepl("ExecutablePath$", names(u))],
                 file.path(u[grepl("InstallPath\\.\\(Default\\)$", names(u))], "python.exe"))
    }
    # c) Everything `where` can see on PATH
    r <- run("where", c("python", "python3"), timeout = 30)
    if (r$status == 0L) { cands <- c(cands, trimws(r$output)) }
    # d) Common install folders, in case the above all miss
    la <- Sys.getenv("LOCALAPPDATA")
    pf <- unique(c(Sys.getenv("ProgramFiles"), Sys.getenv("ProgramFiles(x86)"), "C:"))
    cands <- c(cands, Sys.glob(c(
      file.path(la, "Programs", "Python", "Python3*", "python.exe"), # python.org per-user
      file.path(la, "Python", "*", "python.exe"),                   # new install manager
      file.path(pf, "Python3*", "python.exe"),
      file.path(pf, "Python", "Python3*", "python.exe")
    )))
  } else {
    cands <- c(cands, Sys.glob(c(
      "/usr/bin/python3", "/usr/bin/python3.*",
      "/usr/local/bin/python3", "/usr/local/bin/python3.*",
      "/opt/homebrew/bin/python3", "/opt/homebrew/bin/python3.*", # macOS Homebrew (Apple Silicon)
      "~/.pyenv/versions/*/bin/python"
    )))
  }
  cands <- c(cands, unname(Sys.which(c("python3", paste0("python3.", 9:15), "python"))))
  cands <- cands[nzchar(cands)]
  cands <- unique(normalizePath(cands[file.exists(cands)], winslash = "/", mustWork = FALSE))
  if (!length(cands)) {
    return(fail("No Python interpreter found. Please install Python (>= " , minVersion, ")."))
  }
  #
  # 3) Probe them: must run, have venv + ensurepip, and be in the allowed version range
  probe <- \(py) { #py <- cands[1L]
    r <- run(py, c("-c", shQuote("import sys, venv, ensurepip; print(*sys.version_info[:3]); print(sys.executable)")),
             timeout = 30)
    if (r$status != 0L) { return(NULL) }
    o <- r$output
    i <- grep("^[0-9]+ [0-9]+ [0-9]+$", o)
    if (!length(i)) { return(NULL) }
    v <- as.integer(strsplit(o[i[1L]], " ")[[1L]])
    list(path = py, major_minor = paste(v[1L:2L], collapse = "."),
         ver = numeric_version(paste(v, collapse = ".")),
         real = if (length(o) > i[1L]) o[i[1L]+1L] else py)
  }
  info <- Filter(Negate(is.null), lapply(cands, probe))
  if (!length(info)) {
    return(fail("Python was found but none of the interpreters works (or they lack the 'venv' module - on Debian/Ubuntu: apt install python3-venv)."))
  }
  ok <- vapply(info, \(x) {
    (x$ver >= numeric_version(minVersion)) &&
      (is.null(maxVersion) || (x$ver < numeric_version(maxVersion)))
  }, TRUE)
  info <- info[ok]
  if (!length(info)) {
    return(fail("No Python interpreter in the allowed version range (>= ", minVersion, ")."))
  }
  info <- info[order(vapply(info, \(x) { as.character(x$ver) }, ""), decreasing = TRUE)]
  info <- info[order(do.call(c, lapply(info, \(x) { x$ver })), decreasing = TRUE)]
  info <- info[!duplicated(vapply(info, \(x) { normalizePath(x$real, winslash = "/", mustWork = FALSE) }, ""))]
  #
  # 4) Try each candidate, newest first; the first that yields a working parse_sdrf wins
  spec <- if (ontology) { "sdrf-pipelines[ontology]" } else { "sdrf-pipelines" }
  for (x in info) {
    v <- file.path(venvRoot, paste0(venvPrefix, x$major_minor))
    say("Trying Python ", as.character(x$ver), " (", x$path, ")")
    unlink(v, recursive = TRUE, force = TRUE)
    r <- run(x$path, c("-m", "venv", shQuote(v)))
    if ((r$status != 0L) || !file.exists(venvPy(v))) {
      say("  could not create venv: ", tail3(r$output)); unlink(v, recursive = TRUE, force = TRUE); next
    }
    # Upgrading pip is nice to have, not essential
    r <- run(venvPy(v), c("-m", "pip", "install", "--upgrade", "pip"))
    if (r$status != 0L) { say("  pip upgrade failed (continuing): ", tail3(r$output)) }
    r <- run(venvPy(v), c("-m", "pip", "install", shQuote(spec)))
    if (r$status != 0L) {
      say("  installing ", spec, " failed: ", tail3(r$output))
      unlink(v, recursive = TRUE, force = TRUE); next
    }
    if (!parserOK(v)) {
      say("  installed, but parse_sdrf does not run.")
      unlink(v, recursive = TRUE, force = TRUE); next
    }
    say("Success with Python ", x$major_minor)
    return(finish(v, x$major_minor))
  }
  fail(paste0("Could not install ", spec, " with any of the Python versions found (",
              paste(vapply(info, \(x) { x$major_minor }, ""), collapse = ", "),
              "). Details in the $log element or run with verbose = TRUE."))
}
