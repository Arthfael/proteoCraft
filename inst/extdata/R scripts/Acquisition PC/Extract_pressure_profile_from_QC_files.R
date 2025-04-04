# Script to iterate through all QC files and, for each, it it doesn't find the pressure profiles:
# - convert it to mzML containing pressure information
# - run a python script to save the pressure profiles
# - remove the mzML
# For now the script is not scheduled to run automatically in the background, but this could be done in the future.

TargQCDir <- "//archive3.ist.local/archive-lsfgrp/MS/MS Proteomics service - Archive/Standard runs/QCloud/processed"
ThRwFlPrsr <- "C:/PROGRA~1/ThermoRawFileParser/ThermoRawFileParser.exe"
Pr2CSV <- "C:/Users/Thermo/Downloads/mzML_pressure_to_csv.py"


if ((dir.exists(TargQCDir))&&(file.exists(ThRwFlPrsr))&&(file.exists(Pr2CSV))) {
  fls <- list.files(TargQCDir, recursive = TRUE)
  raws <- paste0(TargQCDir, "/", grep("\\.raw$", fls, value = TRUE, ignore.case = TRUE))
  LdPmpPrss <- gsub("\\.raw", "_AD_1_LoadingPump_Pressure_0.csv", raws)
  NCPmpPrss <- gsub("\\.raw", "_AD_2_NC_Pump_Pressure_0.csv", raws)
  mzMLs <- gsub("\\.raw$", ".mzML", raws)
  tst1 <- file.exists(LdPmpPrss)
  tst2 <- file.exists(NCPmpPrss)
  w <- which(tst1+tst2 < 2) 
  if (length(w)) {
    # Use for loop so it can more easily be stopped then restarted
    for (i in w) { #i <- w[1]
      cat(paste0(raws[i], "\n"))
      tst <- file.exists(mzMLs[i])
      if (!tst) {
        cat("Converting file to mzML...\n")
        # Convert to mzML
        cmd <- paste0(ThRwFlPrsr, " -i=\"", raws[i], "\" -b=\"", mzMLs[i], "\" -z --allDetectors")
        #cat(cmd)
        tst <- !shell(cmd)
        # Interesting behaviour here:
        # - Surrounding the ThermoRawFileParser path with quotes does not work (and this is not because of the "PROGRA~1 issue")
        # - Calling "system" instead of "shell" fails on this PC, not on msmonster01
      }
      if (tst) {
        # Extract pressure profile
        cat("Extracting pressure profile...\n")
        #cmd <- paste0(Python, " ", Pr2CSV, " \"", mzMLs, "\"")
        cmd <- paste0("python ", Pr2CSV, " \"", mzMLs[i], "\"")
        #cat(cmd)
        tst <- !shell(cmd)
      }
      if ((tst)&&(file.exists(LdPmpPrss[i]))&&(file.exists(LdPmpPrss[i]))) {
        # Remove mzML
        unlink(mzMLs[i])
      }
      cat("\n")
    }
  }
}
  
