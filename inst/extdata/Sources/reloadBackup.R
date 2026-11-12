# Re-load a backup?
#
load_a_Bckp %<o% c(TRUE, FALSE)[match(svDialogs::dlg_message("Re-load a backup?", "yesno")$res, c("yes", "no"))]
#ls0 <- ls(.GlobalEnv, all.names = TRUE)
if (load_a_Bckp) {
  tst <- try({
    locDirs %<o% openxlsx2::read_xlsx(locDirs_fl)
    load_Bckp(startDir = locDirs$Path[match("Temporary folder", locDirs$Folder)])
  }, silent = TRUE)
  #
  # Check for conflicts with our functions
  ls1 <- ls(.GlobalEnv, all.names = TRUE)
  #ls1 <- base::setdiff(ls1, ls0)
  ls2 <- ls(getNamespace("proteoCraft"), all.names = TRUE)
  inter <- base::intersect(ls1, ls2)
  if (length(inter)) {
    # rmv <- c(TRUE, FALSE)[match(svDialogs::dlg_message("Some reloaded objects overlap with objects from proteoCraft. Remove them?", "yesno")$res,
    #                             c("yes", "no"))]
    # if (rmv) {
    #   rm(list = inter)
    #   .obj <- base::setdiff(.obj, inter)
    # }
    inter2 <- vapply(inter, \(x) { paste(c(x, rep(" ", 250L-nchar(x))), collapse = "") }, "")
    rmv <- svDialogs::dlg_list(inter2,
                               inter2,
                               TRUE,
                               "These reloaded objects overlap with objects from proteoCraft. Remove them?")$res
    if (length(rmv)) {
      rmv <- sub(" *$", "", rmv)
      rm(list = rmv)
      .obj <- base::setdiff(.obj, rmv)
    }
  }
  #
  # Update values!
  RPath %<o% as.data.frame(library()$results)
  RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
  libPath %<o% paste0(RPath, "/proteoCraft")
  homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
  xplorSrc %<o% paste0(libPath, "/extdata/Sources/xplorData.R")
  locDirs_fl %<o% paste0(homePath, "/Default_locations.xlsx")
  locDirs %<o% openxlsx2::read_xlsx(locDirs_fl)
  if (exists("N.clust")) {
    tmp <- max(c(round(parallel::detectCores()*0.95)-1L, 1L))
    if (N.clust > tmp) {
      warning("This backup was clearly run on a computer with more CPUs, reloading with fewer cores")
    }
  }
}
