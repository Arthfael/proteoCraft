# Run local startup scripts
locScriptsDir <- sub("/[^/]+$", "/proteoCraft_localScripts", locDirs$Path[match("Temporary folder", locDirs$Folder)])
if ((length(locScriptsDir) == 1L) && dir.exists(locScriptsDir)) {
  # Load local startup scripts
  strtScripts <- list.files(locScriptsDir, ".R$", full.names = TRUE)
  # Also handle links
  strtScriptsLnks <- list.files(locScriptsDir, ".R\\.lnk$", full.names = TRUE)
  if (length(strtScriptsLnks)) {
    # Resolve Windows .lnk targets using PowerShell
    get_shortcut_target <- \(lnk_path) {
      lnk_path <- normalizePath(lnk_path, mustWork = TRUE)
      result <- system2("powershell.exe",
                        args = c("-NoProfile",
                                 "-NonInteractive",
                                 "-Command",
                                 paste0("$shell = New-Object -ComObject WScript.Shell; ",
                                        "$shell.CreateShortcut($args[0]).TargetPath"),
                                 lnk_path),
                        stdout = TRUE,
                        stderr = TRUE)
      if ((!length(result)) || (!nzchar(result[1L]))) {
        stop("Could not resolve Windows shortcut: ", lnk_path)
      }
      return(result[1L])
    }
    strtScriptsLnks <- vapply(strtScriptsLnks, get_shortcut_target, "")
    strtScripts <- union(strtScripts, strtScriptsLnks)
  }
  # Source all startup scripts
  if (length(strtScripts)) {
    strtScripts <- normalizePath(strtScripts, winslash = "/")
    strtScripts <- strtScripts[order(strtScripts)]
    for (fl in strtScripts) { source(fl) }
  }
}
