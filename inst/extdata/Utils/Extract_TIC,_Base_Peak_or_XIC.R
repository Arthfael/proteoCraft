if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
if (!require(rawrr)) { BiocManager::install("rawrr", force = TRUE) }
require(rawrr)
#installRawFileReaderDLLs()
#installRawrrExe()

Raws <- choose.files(".raw", "Select one or a few \".raw\" files", multi = TRUE)
#Raws <- "Q:\\MS\\Acquired data\\Standard runs\\QCloud_QC01_20201118180416.raw"
masses <- c(524, 670)
tol <- 100
# Extract all TICs/Base Peaks/the XICs for masses of interest, from all raw files
TICs <- lapply(Raws, function(raw) { readChromatogram(raw, type = "tic") })
BasePeaks <- lapply(Raws, function(raw) { readChromatogram(raw, type = "bpc") })
XICs <- lapply(Raws, function(raw) { readChromatogram(raw, masses, tol) })

# Plots a single XIC
readChromatogram(Raws[1], masses, tol) |> plot()


