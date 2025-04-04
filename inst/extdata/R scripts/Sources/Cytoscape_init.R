### Check that CytoScape is installed and can run, then launch it.
# This will be used to run a ClueGO analysis and to draw protein interaction networks 
if (CytoScape) {
  #suppressMessages(tst <- try(cytoscapePing(), silent = TRUE))
  #tst <- "tst"
  #class(tst) <- "try-error" # Create dummy tst
  dlg_message("Check that no other RStudio session (or other process) is using Cytoscape before clicking \"ok\" to continue...", "ok")
  shell(paste0("open \"", CytoScVrs, "\""))
  suppressMessages(tst <- try(cytoscapePing(), silent = TRUE))
  kount <- 0
  while (("try-error" %in% class(tst))&&(kount < 12)) { # Give it 1 min max!
    kount <- kount + 1
    Sys.sleep(5)
    suppressMessages(tst <- try(cytoscapePing(), silent = TRUE))
  }
  if ("try-error" %in% class(tst)) {
    msg <- "Could not connect to Cytoscape! Skipping creation of network .cx files..."
    ReportCalls <- AddMsg2Report()
    CytoScape <- FALSE
  }
}
