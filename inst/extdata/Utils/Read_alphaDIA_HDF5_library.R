# Read alphaDIA HDF5 format spectra library

file_path <- ".../alphaDIA/library/speclib.mbr.hdf"

myLib <- rhdf5::h5read(file_path, "library", bit64conversion = "double")
names(myLib)

fragment_intensity_df <- as.data.frame(myLib$fragment_intensity_df)
fragment_mz_df <- as.data.frame(myLib$fragment_mz_df)
precursor_df <- as.data.frame(myLib$precursor_df)
# Information linking precursor_df to corresponding fragments is in columns "frag_start_idx" "frag_stop_idx"
fullLib <- apply(precursor_df[, c("frag_start_idx", "frag_stop_idx")], 1L, \(x) { x[[1L]]:x[[2L]] })
fullLib <- setNames(fullLib, 1L:nrow(precursor_df))
fullLib <- utils::stack(fullLib)
colnames(fullLib) <- c("libFragment", "libPrecursor")
min(fullLib$libFragment)
# RRRRRRRAAHAHHHHAHHHHH DO I HATE IT WHEN INDEXING STARTS AT 0!!!
# Ahem...
# Anyhow...
fullLib$libFragment <- fullLib$libFragment+1L
k1 <- colnames(precursor_df)
k2 <- colnames(fragment_mz_df)
k3 <- colnames(fragment_intensity_df)
fullLib[, k1] <- precursor_df[fullLib$libPrecursor, k1]
fullLib[, paste0("m/z_", k2)] <- fragment_mz_df[fullLib$libFragment, k2]
fullLib[, paste0("Intensity_", k3)] <- fragment_intensity_df[fullLib$libFragment, k3]
fullLib$libFragment <- NULL
fullLib$libPrecursor <- NULL
