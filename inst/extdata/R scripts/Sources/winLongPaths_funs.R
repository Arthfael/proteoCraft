################################
# Windows long paths functions #
################################
#
# Safe(r) series of functions for file manipulations in Windows:
# These are designed to be more resilient against Windows long path issues,
# when copying data to a server resource which does not have long paths enabled.
#
lp <- function(p) {
  p <- normalizePath(p, winslash = "\\", mustWork = FALSE)
  paste0("\\\\?\\", p)
}
safe_file_exists <- function(files) {
  res <- file.exists(files)
  #res <- as.logical(file.access(files) + 1) # Alternative
  w <- which(!res)
  if (length(w)) {
    # That's the GODDAM Windows long path issue!!!
    # Luckily, chatGPT has a great solution!
    p <- normalizePath(files[w], winslash = "\\")
    p2 <- paste0("\\\\?\\", p)
    res[w] <- file.exists(p2)
  }
  return(res)
}
safe_dir_exists <- function(dirs) { # Not tested!
  res <- dir.exists(dirs)
  w <- which(!res)
  if (length(w)) {
    # That's the GODDAM Windows long path issue!!!
    # Luckily, chatGPT has a great solution!
    res[w] <- dir.exists(lp(dirs[w]))
  }
  return(res)
}
safe_file_size <- function(files) {
  res <- file.size(files)
  w <- which(is.na(res))
  if (length(w)) {
    res[w] <- file.size(lp(files[w]))
  }
  return(res)
}
is_system_or_hidden <- function(file) {
  out <- suppressWarnings(system2("attrib", shQuote(file), stdout = TRUE, stderr = TRUE))
  if (!length(out)) { return(FALSE) }
  out <- gsub(proteoCraft::topattern(file, start = FALSE), "", gsub("\\\\", "/", out))
  return(grepl("[SH]", out))
}
safe_file_copy <- function(from, to, copy.date = TRUE, overwrite = TRUE) {
  ok <- file.copy(lp(from), lp(to), copy.date = copy.date, overwrite = overwrite)
  if (!ok) { ok <- is_system_or_hidden(from) }
  return(ok)
}
safe_dir_create <- function(paths) { # Recursive by default, which is good
  paths <- lp(paths)
  paths <- shQuote(paths)
  for (p in paths) {
    system2("cmd", c("/c", "mkdir", p))
  }
}