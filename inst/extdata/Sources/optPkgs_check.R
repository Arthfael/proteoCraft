# This script asks the user whether to install optional packages and, if necessary, update them and their dependencies,
# so that optional analyses may be run.
optPkgs_opt <- listMelt(.optPkgs)
optPkgs_opt <- aggregate(optPkgs_opt$L1, list(optPkgs_opt$value), .commAnd, sort = TRUE, unique = TRUE)
colnames(optPkgs_opt) <- c("Package", "Analysis")
optPkgs_opt$nc <- nchar(optPkgs_opt$Package)
mxCh <- max(nchar(optPkgs_opt$Package)) + 10L
optPkgs_opt$Pack_ <- optPkgs_opt$Package
optPkgs_opt$Pack_ <- vapply(1L:nrow(optPkgs_opt), \(i) {
  x <- optPkgs_opt$Package[i]
  nc <- optPkgs_opt$nc[i]
  if (nc < mxCh) { x <- paste(c(x, rep(" ", mxCh-nc)), collapse = "") }
  return(x)
}, "")
optPkgs_opt$opt <- do.call(paste, c(optPkgs_opt[, c("Pack_", "Analysis")], sep = " ---[needed for]--->           "))
allOpt <- "install all optional packages"
depOpt <- " + only update dependencies if explicitly required"
opt <- c(allOpt,
         depOpt,
         optPkgs_opt$opt)
myScriptMyChoice <- svDialogs::dlg_list(opt,
                                        allOpt,
                                        TRUE,
                                        "Select optional package(s) to install/update...")$res
updt <- !depOpt %in% myScriptMyChoice
m <- if (allOpt %in% myScriptMyChoice) {
  1L:nrow(optPkgs_opt)
} else {
  match(myScriptMyChoice, optPkgs_opt$opt)
}
pkgs2install <- optPkgs_opt$Package[m]
if (length(pkgs2install)) {
  pak::pak(pkgs2install, upgrade = updt)
}

pkgsCheck %<o% (vapply(.optPkgs, \(pkgs) {
  sum(!vapply(pkgs, \(pkg) { requireNamespace(pkg) }, TRUE))
}, 1L) == 0L)
