require(rawrr)
#rawrr::installRawFileReaderDLLs()
#rawrr::installRawrrExe()

Raws <- choose.files(".raw", "Select one or a few \".raw\" files", multi = TRUE)

masses <- c(524, 670)
tol <- 100
# Extract all TICs/Base Peaks/the XICs for masses of interest, from all raw files
TICs <- lapply(Raws, \(raw) { rawrr::readChromatogram(raw, type = "tic") })
BasePeaks <- lapply(Raws, \(raw) { rawrr::readChromatogram(raw, type = "bpc") })
XICs <- lapply(Raws, \(raw) { rawrr::readChromatogram(raw, masses, tol) })

# Plots a single XIC
rawrr::readChromatogram(Raws[1L], masses, tol) |> plot()


