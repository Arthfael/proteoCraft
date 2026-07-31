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
      lnk_path <- gsub("'", "''", lnk_path, fixed = TRUE)
      result <- system2("powershell.exe",
                        args = c("-NoProfile",
                                 "-Command",
                                 paste0("$s = New-Object -ComObject WScript.Shell; ",
                                        "$s.CreateShortcut('", lnk_path, "').TargetPath")),
                        stdout = TRUE,
                        stderr = TRUE)
      if ((!length(result)) || (!nzchar(result[1L]))) {
        stop("Could not resolve shortcut: ", lnk_path)
      }
      return(normalizePath(result[1L], winslash = "/"))
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
