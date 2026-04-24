### Check that Cytoscape is installed and can run, then launch it.
# This will be used to run a ClueGO analysis and to draw protein interaction networks 
if ((!exists("CytoScape_check"))||(length(CytoScape_check) != 1L)||(!is.logical(CytoScape_check))||(is.na(CytoScape_check))) {
  CytoScape_check <- FALSE
}
CytoScape_check %<o% CytoScape_check
if (!CytoScape_check) {
  dlg_message("Check that no other RStudio session (or other process) is using Cytoscape before clicking \"ok\" to continue...", "ok")
  CytoScape_check <- TRUE
}
invisible(suppressMessages({
  if (CytoScape) {
    tst <- try(RCy3::cytoscapePing(), silent = TRUE)
    kount <- 0L
    while ((inherits(tst, "try-error"))&&(kount < 12L)) { # Give it 1 min max!
      kount <- kount + 1L
      if (kount == 1L) {
        #dlg_message("Check that no other RStudio session (or other process) is using Cytoscape before clicking \"ok\" to continue...", "ok")
        shell(paste0("open \"", CytoScExe[1L], "\""))
      }
      Sys.sleep(5L)
      tst <- try(RCy3::cytoscapePing(), silent = TRUE)
    }
    if (inherits(tst, "try-error")) {
      msg <- "Could not connect to Cytoscape! Skipping creation of network .cx files..."
      ReportCalls <- AddMsg2Report()
      CytoScape <- FALSE
    }
  }
}))
