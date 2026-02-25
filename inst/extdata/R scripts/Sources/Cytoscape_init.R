### Check that Cytoscape is installed and can run, then launch it.
# This will be used to run a ClueGO analysis and to draw protein interaction networks 
if ((!exists("CytoScape_check"))||(length(CytoScape_check) != 1)||(!is.logical(CytoScape_check))||(is.na(CytoScape_check))) {
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
    kount <- 0
    while (("try-error" %in% class(tst))&&(kount < 12)) { # Give it 1 min max!
      kount <- kount + 1
      if (kount == 1) {
        #dlg_message("Check that no other RStudio session (or other process) is using Cytoscape before clicking \"ok\" to continue...", "ok")
        shell(paste0("open \"", CytoScExe[1], "\""))
      }
      Sys.sleep(5)
      tst <- try(RCy3::cytoscapePing(), silent = TRUE)
    }
    if ("try-error" %in% class(tst)) {
      msg <- "Could not connect to Cytoscape! Skipping creation of network .cx files..."
      ReportCalls <- AddMsg2Report()
      CytoScape <- FALSE
    }
  }
}))
