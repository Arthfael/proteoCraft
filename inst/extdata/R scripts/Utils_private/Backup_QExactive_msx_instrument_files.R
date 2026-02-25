# Run this through a scheduled Windows task, e.g.:
#
# ...> "C:\PROGRA~1\R\R-4.4.0\bin\Rscript.exe" "B:\group\lsfgrp\Mass_Spec\R scripts\Backup_QExactive_msx_instrument_files.R"
#
# or 
#

#
#
isInteractive <- interactive()
targDir <- "C:/Xcalibur/system/Exactive/instrument/msx_instrument_files"
if (!dir.exists(targDir)) {
  msg <- "Local target folder not found! Calibration backup script failed, investigate!"
  if (isInteractive) {
    svDialogs::dlg_message(msg, "ok")
  }
} else {
  dstDir <- "B:/group/mspecgrp/MassSpecgrp_Updated/Instruments/Thermo/Thermo Q-Exactive HF/Master Calibration File backup - QE-HF1"
  if (!dir.exists(dstDir)) {
    if (isInteractive) {
      svDialogs::dlg_message("Destination folder not found! Calibration backup script failed, investigate!", "ok")
    }
  } else {
    Fls <- list.files(targDir, recursive = TRUE, full.names = TRUE)
    if (length(Fls)) {
      today <- Sys.Date()
      nuDir <- paste0(dstDir, "/Bckp_", today)
      if (!dir.exists(nuDir)) { dir.create(nuDir, recursive = TRUE) }
      nuFls <- paste0(nuDir, "/", basename(Fls))
      w <- which((!file.exists(nuFls))|(file.size(nuFls) != file.size(Fls)))
      fs::file_copy(Fls[w], nuFls[w], overwrite = TRUE)
    }
  }
}
