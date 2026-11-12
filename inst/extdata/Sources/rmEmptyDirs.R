# Remove empty directories:
dirlist <- base::setdiff(list.dirs(wd, full.names = TRUE, recursive = TRUE),
                         wd)
dirlist <- dirlist[order(nchar(dirlist), decreasing = TRUE)]
for (dir in dirlist) { #d <- dirlist[1L]
  if (!length(list.files(dir, recursive = TRUE))) {
    unlink(dir, recursive = TRUE)
    dirlist <- base::setdiff(dirlist, dir)
  }
}
