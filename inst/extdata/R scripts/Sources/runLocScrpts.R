# Run local startup scripts
locScriptsDir <- gsub("/[^/]+$", "/proteoCraft_localScripts", locDirs$Path[match("Temporary folder", locDirs$Folder)])
if ((length(locScriptsDir) == 1L)&&(dir.exists(locScriptsDir))) {
  # Load local startup scripts
  strtScripts <- list.files(locScriptsDir, ".R$", full.names = TRUE)
  # Also handle links
  strtScriptsLnks <- list.files(locScriptsDir, ".R\\.lnk$", full.names = TRUE)
  if (length(strtScriptsLnks)) {
    if (!require(RDCOMClient)) { pak::pak("omegahat/RDCOMClient") }
    library(RDCOMClient)
    get_shortcut_target <- \(lnk_path) { # Thanks istaGPT...
      #lnk_path <- strtScriptsLnks[1L]
      shell <- RDCOMClient::COMCreate("WScript.Shell")
      shortcut <- shell$CreateShortcut(normalizePath(lnk_path, mustWork = TRUE))
      shortcut$TargetPath()
    }
    strtScriptsLnks <- vapply(strtScriptsLnks, get_shortcut_target, "")
    strtScripts <- union(strtScripts, strtScriptsLnks)
  }
  # Source all startup scripts
  if (length(strtScripts)) {
    strtScripts <- strtScripts[order(strtScripts)]
    for (fl in strtScripts) { source(fl) }
  }
}
