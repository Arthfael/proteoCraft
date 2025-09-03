alphaDIA_fl <- "D:/groups_temp/novagrp/LCMS_GaNoDMarano3_p/alphaDIA/precursors.tsv"
alphaDIA <- data.table::fread(alphaDIA_fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
# As I read it here https://alphadia.readthedocs.io/en/latest/methods/output-format.html
# the output speclib.mbr.hdf contains all fragment-level identifications...
alphaDIA_LibFl <- gsub("/[^/]+$", "/speclib.mbr.hdf", alphaDIA_fl)

if (!require(rhdf5)) { pak::pkg_install("rhdf5") }
library(rhdf5)
readHDF5lib <- function(file) {
  #suppressWarnings(myLib1 <- rhdf5::h5read(file, "library", bit64conversion = "double", read.attributes = TRUE))
  # NB: for now bit64conversion = "double" is hard-coded here:
  # I get warnings but those only suggest (mild) loss of precision, not that some low values go to NA.
  # Trying different bit64conversion parameters, in "precursor_df" the data is fundamentally corrupted to NA during import for some values...
  # myLib2 <- rhdf5::h5read(file, "library") # ... either using the default behavior...
  # myLib3 <- rhdf5::h5read(file, "library", bit64conversion = "bit64") # ... or explicit "bit64".
  # identical(myLib1$fragment_intensity_df, myLib2$fragment_intensity_df)
  # identical(myLib1$fragment_intensity_df, myLib3$fragment_intensity_df)
  # identical(myLib1$fragment_mz_df, myLib2$fragment_mz_df)
  # identical(myLib1$fragment_mz_df, myLib3$fragment_mz_df)
  # identical(myLib1$precursor_df, myLib2$precursor_df)
  # identical(myLib1$precursor_df, myLib3$precursor_df)
  # # Let's look at precursor_df:
  # pDF1 <- as.data.frame(myLib1$precursor_df)
  # pDF2 <- as.data.frame(myLib2$precursor_df)
  # pDF3 <- as.data.frame(myLib3$precursor_df)
  # i1 <- which(vapply(colnames(pDF1), function(x) { sum(c("integer", "numeric", "double") %in% typeof(pDF1[[x]])) }, 1) > 0)
  # i2 <- which(vapply(colnames(pDF2), function(x) { sum(c("integer", "numeric", "double") %in% typeof(pDF2[[x]])) }, 1) > 0)
  # i3 <- which(vapply(colnames(pDF3), function(x) { sum(c("integer", "numeric", "double") %in% typeof(pDF3[[x]])) }, 1) > 0)
  # # (Don't use class above are some get "integer64" depending on the bit64conversion parameter)
  # colnames(pDF3[i1[which(!i1 %in% i3)]])
  # w1 <- which(is.na(as.matrix(pDF1[, i1])), arr.ind = TRUE)
  # w2 <- which(is.na(as.matrix(pDF2[, i2])), arr.ind = TRUE)
  # w3 <- which(is.na(as.matrix(pDF3[, i3])), arr.ind = TRUE)
  #
  #
  # Run this line to check the internal structure of the HDF5 file:
  #rhdf5::h5ls(file)
  #
  fragment_intensity_df <- as.data.frame(rhdf5::h5read(file, "library/fragment_intensity_df", bit64conversion = "double"))
  fragment_mz_df <- as.data.frame(rhdf5::h5read(file, "library/fragment_mz_df", bit64conversion = "double"))
  precursor_df <- as.data.frame(rhdf5::h5read(file, "library/precursor_df", bit64conversion = "double"))
  # Information linking precursor_df to corresponding fragments is in columns "frag_start_idx" "frag_stop_idx"
  fullLib <- apply(precursor_df[, c("frag_start_idx", "frag_stop_idx")], 1, function(x) { x[[1]]:x[[2]] })
  fullLib <- setNames(fullLib, 1:nrow(precursor_df))
  fullLib <- utils::stack(fullLib)
  colnames(fullLib) <- c("libFragment", "libPrecursor")
  #min(fullLib$libFragment)
  # RRRRRRRAAHAHHHHAHHHHH DO I HATE IT WHEN INDEXING STARTS AT 0!!!
  # Ahem...
  # Anyhow...
  fullLib$libFragment <- fullLib$libFragment+1
  k1 <- colnames(precursor_df)
  k2 <- colnames(fragment_mz_df)
  k3 <- colnames(fragment_intensity_df)
  fullLib[, k1] <- precursor_df[fullLib$libPrecursor, k1]
  fullLib[, paste0("m/z_", k2)] <- fragment_mz_df[fullLib$libFragment, k2]
  fullLib[, paste0("Intensity_", k3)] <- fragment_intensity_df[fullLib$libFragment, k3]
  fullLib$libFragment <- NULL
  fullLib$libPrecursor <- NULL
  return(fullLib)
}
alphaDIA_Lib <- readHDF5lib(alphaDIA_LibFl)





#DIANN_fl <- "D:/groups_temp/novagrp/LCMS_GaNoDMarano3_p/REANAL-FP_with_buryryl/report.parquet"
#DIANN <- arrow::read_parquet(DIANN_fl1)
DIANN_fl <- "D:/groups_temp/novagrp/LCMS_GaNoDMarano2_p/diaNN/diaNN_w_comblib/report.tsv"
DIANN <- data.table::fread(DIANN_fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
list.files(dirname(DIANN_fl))

DIANN_LibFl <- gsub("/report\\.tsv$", "/lib.tsv", DIANN_fl)
file.exists(DIANN_LibFl)
DIANN_Lib <- data.table::fread(DIANN_LibFl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)

View(DIANN)
tst <- vapply(strsplit(DIANN$Fragment.Quant.Raw, ";"), length, 1)
summary(tst)
#tst <- vapply(strsplit(DIANN$Fragment.Quant.Corrected, ";"), length, 1) # Does not always exist!!!
#summary(tst)

DIANN$Fragment.Correlations
DIANN$Fragment.Quant.Corrected
colnames(DIANN)

# No MS2 info in precursors file, I should work from hdf file:


