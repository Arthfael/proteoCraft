HypTun <- readLines(HypTunFl)
HypTun <- gsub("'tuning_results/' \\+ str\\(args.output_file\\) \\+ '/Gridsearch_Results_'",
               "str(args.output_file) + '/Gridsearch_Results_'",
               HypTun)
HypTun <- gsub("'n_conv_layers'", "'n_conv_layers_'", HypTun)
N <- 1
if (Step == 1) { kol <- paramKol }
if (Step == 2) { kol <- paramKol[which(paramKol %in% colnames(currParam))] }
for (param in kol) { #param <- kol[1]
  lstnm <- paste0(param, "_list")
  g <- grep(paste0("^", lstnm, " = "), HypTun)
  if (!length(g)) {
    lstnm <- paste0(param, "_rate_list")
    g <- grep(paste0("^", lstnm, " = "), HypTun)
  }
  if (!length(g)) {
    lstnm <- paste0(param, "_size_list")
    g <- grep(paste0("^", lstnm, " = "), HypTun)
  }
  stopifnot(length(g) == 1)
  if (Step == 1) { tmp <- unlist(strsplit(gsub(".*\\[|\\].*", "", HypTun[g]), ", ?")) }
  if (Step == 2) { tmp <- sort(unique(currParam[[param]])) }
  N <- N*length(tmp)
  if (Step == 2) { HypTun[[g]] <- paste0(lstnm, " = [", paste(tmp, collapse = ", "), "]") }
}
# Bug fix
w <- which(HypTun == "test_runs = np.array([[2*i, 2*i+1] for i in test_runs]).astype(int).ravel()")
HypTun <- c(HypTun[1:w],
            "#################################################################################################",
            "### Insertion to remove cases where the above lines generate run indices outside of the range ###",
            "#################################################################################################",
            "train_runs = train_runs[train_runs+1 <= q_df.shape[1]]",
            "test_runs = test_runs[test_runs+1 <= q_df.shape[1]]",
            "#################################################################################################",
            "### Insertion ends                                                                            ###",
            "#################################################################################################",
            HypTun[(w+1):length(HypTun)])
write(HypTun, paste0(tuneDir, "/HyperTuning.py"))

GrdSrchFls <- list.files(tuneDir, "Gridsearch_Results_", full.names = TRUE)
if (Step == 1) {
  Step1Tries <- min(c(100, N))
  if (length(GrdSrchFls)) {
    GrdSrchFls <- GrdSrchFls[order(file.info(GrdSrchFls)$ctime, decreasing = FALSE)]
    GrdSrchFls <- GrdSrchFls[1:min(c(length(GrdSrchFls), Step1Tries))]
  }
}
if (Step == 2) { Step2Tries <- min(c(50, N)) }
cmd <- paste0("\"", pyPath, "\" \"", tuneDir, "/HyperTuning.py\" --peptide_file \"", OHfl, "\" --n_runs \"", length(Grps),
              "\" --seq_length \"60\" --output_file \"", tuneDir, "\"")
# This command will only work without reticulate because it doesn't edit the script!!!
# Reticulate cannot provide arguments to a python script.
# If trying to run it in reticulate, edit the script (see example code above) and save it locally
#cat(cmd)
