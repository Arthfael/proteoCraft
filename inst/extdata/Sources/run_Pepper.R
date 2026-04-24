# Run Pepper correction
if (runPepper) {
  PepScrptsDir %<o% paste0(RPath, "/proteoCraft/extdata/Pepper")
  stopifnot(dir.exists(PepScrptsDir))
  pyInit <- paste0(PepScrptsDir, "/Python_start.R")
  source(pyInit)
  #
  # Placeholder for when this gets ported to the main script
  # if (scrptType == "withReps") {
  #   if ("Pepper" %in% names(ev.col)) { ref <- ev.col[match("Pepper", names(ev.col))-1] } else {
  #     ref <- ev.col[length(ev.col)]
  #   }
  #   nuRef <- ev.col["Pepper"] <- paste0("Pepper ", ev.col["Original"])
  # }
  if (scrptType == "noReps") {
    if ("Pepper" %in% names(int.cols)) { ref <- int.cols[match("Pepper", names(int.cols))-1L] } else {
      ref <- int.cols[length(int.cols)]
    }
    nuRef <- int.col <- int.cols["Pepper"] <- "Pepper-adj. intensity"
  }
  #
  # Detect recommended training parameters for all datasets
  trainParam <- suppressWarnings(readLines(paste0(PepScrptsDir, "/train_models.sh")))
  args <- c("peptide_file", "n_runs", "seq_length", "output_file", "filter_size", "n_filters", "n_layers", "n_nodes",
            "dropout", "learning_rate", "batch", "random_run")
  g <- grep("^python3 ", trainParam)
  l <- length(trainParam)
  trainParam <- sapply(args, \(arg) {
    sapply(seq_along(g), \(x) {
      h <- grep(paste0("--", arg), trainParam)
      h <- gsub(paste0(".*'--", arg, "' +"), "", trainParam[h[which(h %in% g[x]:c(g, l)[x+1L])]])
      pat1 <- paste0("^", substr(h, 1L, 1L))
      pat2 <- paste0(substr(h, 1L, 1L), ".*$")
      h <- gsub(pat2, "", gsub(pat1, "", h))
      return(h)
    })
  })
  trainParam <- as.data.frame(trainParam)
  rownames(trainParam) <- gsub(".+/|\\.tsv$", "", trainParam$peptide_file)
  #
  # Hyper-parameter tuning
  # a) Prepare dataset for one-hot-encoding
  smplKol <- "Experiment"
  cat("Samples (= Experiments):\n -", paste0(unique(SamplesMap[[smplKol]]), collapse = "\n - "), "\n")
  w <- which(ev[1L,] %in% SamplesMap[[smplKol]])
  stopifnot(length(w) > 0L)
  pepDir <- paste0(wd, "/Pepper")
  if (!dir.exists(pepDir)) { dir.create(pepDir, recursive = TRUE) }
  tempEv <- Pepper_TrainingData(ev, Modifs, SamplesMap, smplKol, smplKol, path = paste0(pepDir, "/ModPep4Pepper.tsv"), intCol = ref)
  Modifs <- tempEv$PTMs
  tempEv <- tempEv$Data
  #
  # Edit and run python script for One-Hot-Encoding
  #################################################
  Grps <- colnames(tempEv)[(max(grep("^Charge [0-9]+$", colnames(tempEv)))+1L):length(colnames(tempEv))]
  lGrps <- length(Grps)
  cat(paste0("Dataset = ", dtstNm, ", number of samples = ", lGrps, "\n"))
  #
  OHfl <- paste0(pepDir, "/OneHotEncodedPepQuants.tsv")
  OHE <- readLines(paste0(PepScrptsDir, "/onehot_encode_peptide_sequences.py"))
  arg1 <- which(OHE == "input_filename = sys.argv[1L]")
  arg2 <- which(OHE == "runs = int(sys.argv[2L])")
  arg3 <- which(OHE == "output_filename = sys.argv[3L]")
  OHE[arg1] <- paste0("input_filename = \"", paste0(pepDir, "/ModPep4Pepper.tsv"), "\"")
  OHE[arg2] <- paste0("runs = int(", length(Grps), ")")
  OHE[arg3] <- paste0("output_filename = \"", OHfl, "\"")
  write(OHE, paste0(pepDir, "/OneHotEncode.py"))
  #
  cmd <- paste0("python \"", pepDir, "/OneHotEncode.py\"")
  #cat(cmd)
  system(cmd)
  stopifnot(file.exists(OHfl))
  #
  # Hyperparameter tuning:
  tuneDir <- paste0(pepDir, "/Hyperparameter tuning")
  if (!dir.exists(tuneDir)) { dir.create(tuneDir, recursive = TRUE) }
  CleanUp <- FALSE
  GrdSrchFls <- list.files(tuneDir, "Gridsearch_Results_", full.names = TRUE)
  if (length(GrdSrchFls)) {
    opt <- setNames(c(FALSE, TRUE), c("No                                                                                                                                              ",
                                      "Yes                                                                                                                                             "))
    CleanUp <- opt[dlg_list(names(opt), names(opt)[1L],
                            title = paste0(dtstNm, ": shall we start from scratch hyperparameter tuning?"))$res]
    names(CleanUp) <- NULL
  }
  # NB: This code is written in a way that it will always sample Step1Tries + Step2Tries parameter combinations.
  # 1 - randomly
  Step <- 1L
  paramKol <- c("n_conv_layers", "n_filters", "filter_size", "n_layers", "n_nodes", "dropout", "batch", "learning_rate")
  if (CleanUp) { for (fl in GrdSrchFls) { unlink(fl) } } # Cleanup for fresh tuning
  HypTunFl <- paste0(PepScrptsDir, "/gridseach_parameters_neural_network.py")
  stopifnot(file.exists(HypTunFl))
  HypTunSrc <- paste0(PepScrptsDir, "/HyperTun.R")
  source(HypTunSrc)
  L <- length(GrdSrchFls)
  N <- Step1Tries - L
  kount <- 0L
  while ((N > 0L)&&(kount < N+1L)) {
    clusterExport(parClust, "cmd", envir = environment())
    # Below, we are never running less tests than the number of cores, so we do not waste a round on a few tests.
    tst <- parSapply(parClust, 1L:max(c(N.clust, N)), \(x) {
      try(system(cmd), silent = TRUE)
    })
    GrdSrchFls <- list.files(tuneDir, "Gridsearch_Results_", full.names = TRUE)
    L <- length(GrdSrchFls)
    N <- Step1Tries - L
    kount <- kount + 1L
  }
  Step1Tries <- L # Update it if we got any bonus tests
  GrdSrch <- as.data.frame(t(sapply(strsplit(gsub("\\.tsv$", "", GrdSrchFls), "_"), \(x) {
    x <- unlist(x)
    x <- suppressWarnings(as.numeric(x))
    x <- x[which(!is.na(x))]
    return(x)
  })))
  colnames(GrdSrch) <- paramKol
  tst <- apply(GrdSrch, 2L, \(x) { length(unique(x)) })
  GrdSrch <- GrdSrch[, which(tst > 1L), drop = FALSE]
  paramKol2 <- paramKol[which(tst > 1L)]
  rs <- plyr::rbind.fill(lapply(GrdSrchFls, read.delim))
  rs$X <- NULL
  GrdSrch[, gsub("\\.", "_", colnames(rs))] <- rs
  GrdSrch <- GrdSrch[order(GrdSrch$Val_percent_improvement, decreasing = TRUE),]
  bstImprov <- max(GrdSrch$Val_percent_improvement)
  # 2 - explore immediate neighborhood of best 10 results
  Step <- 2L
  currParam <- GrdSrch[1L:10L, c(paramKol2, "No_of_epochs")]
  source(HypTunSrc)
  L <- length(GrdSrchFls)
  N <- Step1Tries + Step2Tries - L
  kount <- 0L
  while ((N)&&(kount < N+1L)) {
    clusterExport(parClust, "cmd", envir = environment())
    # Below, we are never running less tests than the number of cores, so we do not waste a round on a few tests.
    tst <- parSapply(parClust, 1L:max(c(N.clust, N)), \(x) {
      try(system(cmd), silent = TRUE)
    })
    GrdSrchFls <- list.files(tuneDir, "Gridsearch_Results_", full.names = TRUE)
    L <- length(GrdSrchFls)
    N <- Step1Tries + Step2Tries - L
    kount <- kount + 1L
  }
  GrdSrch <- as.data.frame(t(sapply(strsplit(gsub("\\.tsv$", "", GrdSrchFls), "_"), \(x) {
    x <- unlist(x)
    x <- suppressWarnings(as.numeric(x))
    x <- x[which(!is.na(x))]
    return(x)
  })))
  colnames(GrdSrch) <- paramKol
  tst <- apply(GrdSrch, 2L, \(x) { length(unique(x)) })
  GrdSrch <- GrdSrch[, which(tst > 1L), drop = FALSE]
  paramKol2 <- paramKol[which(tst > 1L)]
  rs <- plyr::rbind.fill(lapply(GrdSrchFls, read.delim))
  rs$X <- NULL
  GrdSrch[, gsub("\\.", "_", colnames(rs))] <- rs
  GrdSrch <- GrdSrch[order(GrdSrch$Val_percent_improvement, decreasing = TRUE),]
  bstImprov2 <- max(GrdSrch$Val_percent_improvement)
  print("")
  print(paste0(" - Step 1: best improvement = ", round(bstImprov, 2L), "%"))
  print("")
  print(paste0(" - Step 2: best improvement = ", round(bstImprov2, 2L), "%"))
  print("")
  #View(GrdSrch)
  #
  # Get final, best parameters
  bstParam <- GrdSrch[1L, c(paramKol2, "No_of_epochs")]
  print(bstParam)
  bstParam$peptide_file <- OHfl
  bstParam$n_runs <- as.character(length(Grps))
  bstParam$seq_length <- "60"
  bstParam$random_run <- "0" # Obviously: we will do a proper training now!
  bstParam$output_file <- dtstNm
  #
  # Create model, using our own dataset and the above hypertuned parameters
  TrCoeff <- readLines(paste0(PepScrptsDir, "/peptide_coefficient_predictor.py"))
  
  # Parameters from train_models.sh"
  w <- which(TrCoeff %in% c("parser = argparse.ArgumentParser()", "args = parser.parse_args()", "print(args)"))
  TrCoeff[w] <- ""
  basePat <- paste0("^ *", topattern("parser.add_argument(", start = FALSE), " *")
  for (arg in colnames(trainParam)) { #arg <- colnames(trainParam)[1L]
    pat <- paste0(basePat, "\'--", arg)
    g <- grep(pat, TrCoeff)
    stopifnot(length(g) == 1L)
    val <- bstParam[[arg]]
    tst1 <- suppressWarnings(!is.na(as.numeric(val)))
    tst2 <- suppressWarnings(!is.na(as.integer(val)))
    if (tst1) {
      if (tst2) { TrCoeff[g] <- paste0(arg, " = float(", val, ")") } else {
        TrCoeff[g] <- paste0(arg, " = int(", val, ")")
      }
    } else { TrCoeff[g] <- paste0(arg, " = '", val, "'") }
  }
  # Default parameters
  G <- grep(basePat, TrCoeff)
  if (length(G)) {
    for (g in G) { #g <- G[1L]
      tmp <- TrCoeff[g]
      arg <- gsub("\'.*", "", gsub(paste0(basePat, "\'--"), "", tmp))
      stopifnot(nchar(arg) > 0L)
      val <- unlist(strsplit(tmp, "[,\\(] *default *= *"))
      if (!length(val) == 2L) { stop(paste0("No default for argument ", arg, "! Give me a hand here, human!")) }
      val <- gsub(" *[,\\)].*", "", val[2L])
      if (tst1) {
        if (tst2) { TrCoeff[g] <- paste0(arg, " = float(", val, ")") } else {
          TrCoeff[g] <- paste0(arg, " = int(", val, ")")
        }
      } else { TrCoeff[g] <- paste0(arg, " = '", val, "'") }
    }
  }
  TrCoeff <- gsub(" args\\.", " ", TrCoeff)
  TrCoeff <- gsub("\\(args\\.", "(", TrCoeff)
  arg1a <- which(TrCoeff == "data_df = pd.read_csv(peptide_file, sep = '\\t', index_col = 0)")
  TrCoeff[arg1a] <- "data_df = pd.read_csv(peptide_file, sep = '\t', index_col = 0)"
  #TrCoeff[arg4a] <- paste0("predicted_coefficients.to_csv(\"", output_file + \".tsv\", sep = '\\t')")
  g <- grep("trained_models/", TrCoeff) #TrCoeff[g]
  TrCoeff[g] <- gsub("trained_models/", paste0(pepDir, "/"), TrCoeff[g])
  # Small modification to adjust for a deprecation warning
  w <- which(TrCoeff == "init = tf.compat.v1.initialize_all_variables()")
  TrCoeff[w] <- "init = tf.compat.v1.global_variables_initializer()"
  # Bug fix
  w <- which(TrCoeff == "test_runs = np.array([[2*i, 2*i+1] for i in test_runs]).astype(int).ravel()")
  TrCoeff <- c(TrCoeff[1L:w],
               "#################################################################################################",
               "### Insertion to remove cases where the above lines generate run indices outside of the range ###",
               "#################################################################################################",
               "train_runs = train_runs[train_runs+1 <= q_df.shape[1L]]",
               "test_runs = test_runs[test_runs+1 <= q_df.shape[1L]]",
               "#################################################################################################",
               "### Insertion ends                                                                            ###",
               "#################################################################################################",
               TrCoeff[(w+1L):length(TrCoeff)])
  # Add code to save final layer dimensions
  w <- which(TrCoeff == "print(\"Saved model to disk\")")
  TrCoeff <- c(TrCoeff[1L:w],
               "#################################################################################################",
               "### Insertion: this code saves the dimensions of the last weights layer to a small local file ###",
               "#################################################################################################",
               "weights = model.get_weights()",
               "dim = np.shape(weights[len(weights)-1])",
               "dim = np.matrix(dim)",
               #dim = np.reshape(dim, [2,1])",
               paste0("np.savetxt(\"",  pepDir, "/final_layer_dims.tsv\", dim, delimiter = '\t', fmt='%i')"),
               "#################################################################################################",
               "### Insertion ends                                                                            ###",
               "#################################################################################################",
               TrCoeff[(w+1L):length(TrCoeff)])
  write(TrCoeff, paste0(pepDir, "/TrainCoeff.py"))
  
  py_clear_last_error()
  cat(" - Running coefficients training script...\n")
  cmd <- paste0("python \"", pepDir, "/TrainCoeff.py\"")
  #cat(cmd)
  tst <- try(system(cmd), silent = TRUE)
  stopifnot(tst == 0L)
  #
  modlFl <- paste0(pepDir, "/", dtstNm, "_Coefficient_Predictor_Model.h5")
  if (!tst) { Outcome <- file.exists(modlFl) }
  if (Outcome) { cat("Success!!!\n") } else { cat("Failure!!!!!\n") }
  py_clear_last_error()
  #
  # Transfer to whole dataset
  # a) Prepare dataset for one-hot-encoding
  smplKol <- "Experiment"
  w <- which(ev[1L,] %in% SamplesMap[[smplKol]])
  stopifnot(length(w) > 0L)
  tempEv2 <- Pepper_ProcessData(ev, Modifs, SamplesMap, smplKol, smplKol, path = paste0(pepDir, "/ModPep4Pepper.tsv"), intCol = ref,
                                filter = FALSE)
  tempEv2 <- tempEv2$Data
  tmp <- data.frame(Peptides = nrow(tempEv2), Proteins = length(unique(tempEv2$Protein)), Samples = length(SamplesMap$Experiment))
  write.table(tmp, paste0(pepDir, "/Dataset dimensions.tsv"), sep = "\t", row.names = FALSE)
  Grps2 <- colnames(tempEv2)[(max(grep("^Charge [0-9]+$", colnames(tempEv2)))+1L):length(colnames(tempEv2))]
  lGrps2 <- length(Grps2)
  #
  # Edit and run python script for One-Hot-Encoding
  #################################################
  #
  OHfl <- paste0(pepDir, "/OneHotEncodedPepQuants_Full.tsv")
  OHE <- readLines(paste0(PepScrptsDir, "/onehot_encode_peptide_sequences.py"))
  arg1 <- which(OHE == "input_filename = sys.argv[1L]")
  arg2 <- which(OHE == "runs = int(sys.argv[2L])")
  arg3 <- which(OHE == "output_filename = sys.argv[3L]")
  OHE[arg1] <- paste0("input_filename = \"", paste0(pepDir, "/ModPep4Pepper.tsv"), "\"")
  OHE[arg2] <- paste0("runs = int(", length(Grps2), ")")
  OHE[arg3] <- paste0("output_filename = \"", OHfl, "\"")
  write(OHE, paste0(pepDir, "/OneHotEncode.py"))
  #
  cmd <- paste0("python \"", pepDir, "/OneHotEncode.py\"")
  #cat(cmd)
  system(cmd)
  #
  # Edit and run coefficients predictor transfer script
  #####################################################
  TrsfrCoeff <- readLines(paste0(PepScrptsDir, "/peptide_coefficient_predictor_transfer.py"))
  # Parameters from train_models.sh"
  w <- which(TrsfrCoeff %in% c("parser = argparse.ArgumentParser()", "args = parser.parse_args()", "print(args)"))
  TrsfrCoeff[w] <- ""
  basePat <- paste0("^ *", topattern("parser.add_argument(", start = FALSE), " *")
  bstParam2 <- bstParam
  tmp <- read.delim(paste0(pepDir, "/Dataset dimensions.tsv"))
  bstParam2$n_runs <- tmp$Samples
  bstParam2$output_file <- pepDir
  for (arg in colnames(trainParam)) { #arg <- colnames(trainParam)[1L]
    pat <- paste0(basePat, "\'--", arg)
    g <- grep(pat, TrsfrCoeff)
    stopifnot(length(g) == 1L)
    if (arg == "peptide_file") { val <- OHfl } else { val <- bstParam2[[arg]] }
    tst1 <- suppressWarnings(!is.na(as.numeric(val)))
    tst2 <- suppressWarnings(!is.na(as.integer(val)))
    if (tst1) {
      if (tst2) { TrsfrCoeff[g] <- paste0(arg, " = float(", val, ")") } else {
        TrsfrCoeff[g] <- paste0(arg, " = int(", val, ")")
      }
    } else { TrsfrCoeff[g] <- paste0(arg, " = '", val, "'") }
  }
  # Default parameters
  G <- grep(basePat, TrsfrCoeff)
  if (length(G)) {
    for (g in G) { #g <- G[1L]
      tmp <- TrsfrCoeff[g]
      arg <- gsub("\'.*", "", gsub(paste0(basePat, "\'--"), "", tmp))
      stopifnot(nchar(arg) > 0L)
      val <- unlist(strsplit(tmp, "[,\\(] *default *= *"))
      if (!length(val) == 2L) { stop(paste0("No default for argument ", arg, "! Give me a hand here, human!")) }
      val <- gsub(" *[,\\)].*", "", val[2L])
      if (tst1) {
        if (tst2) { TrsfrCoeff[g] <- paste0(arg, " = float(", val, ")") } else {
          TrsfrCoeff[g] <- paste0(arg, " = int(", val, ")")
        }
      } else { TrsfrCoeff[g] <- paste0(arg, " = '", val, "'") }
    }
  }
  TrsfrCoeff <- gsub(" args\\.", " ", TrsfrCoeff)
  TrsfrCoeff <- gsub("\\(args\\.", "(", TrsfrCoeff)
  arg1a <- which(TrsfrCoeff == "data_df = pd.read_csv(peptide_file, sep = '\\t', index_col = 0)")
  TrsfrCoeff[arg1a] <- "data_df = pd.read_csv(peptide_file, sep = '\t', index_col = 0)"
  TrsfrCoeff <- gsub("\"trained_models/\" +\\+ +output_file +\\+ +", "output_file + \"/\" + ", TrsfrCoeff)
  # Small modification to adjust for a deprecation warning
  w <- which(TrsfrCoeff == "init = tf.compat.v1.initialize_all_variables()")
  TrsfrCoeff[w] <- "init = tf.compat.v1.global_variables_initializer()"
  #
  # The most important:
  # The model path is hard-coded into the transfer script. Change the path to the one we want.
  w <- grep("^target_model\\.load_weights\\(", TrsfrCoeff)
  TrsfrCoeff[w] <- paste0("target_model.load_weights('", pepDir, "/", dtstNm, "_Coefficient_Predictor_Model.h5') ")
  #
  # Similarly, the dataset dimensions are hard coded too
  w1 <- which(TrsfrCoeff == "model = define_model(2438, 28)")
  tmp <- read.delim(paste0(pepDir, "/Dataset dimensions.tsv"))
  TrsfrCoeff[w1] <- paste0("model = define_model(", tmp$Proteins, ", len(train_runs))")
  w2 <- which(TrsfrCoeff == "target_model = define_model(5668, 612)")
  tmp <- unlist(read.table(paste0(pepDir, "/final_layer_dims.tsv"), col.names = FALSE))
  TrsfrCoeff[w2] <- paste0("target_model = define_model(", tmp[1L], ", ", tmp[2L], ")")
  w3 <- which(TrsfrCoeff == "            if run_count == 612:")
  TrsfrCoeff <- c(TrsfrCoeff[1L:(w3-1L)],
                  "            self.alphas = tf.Variable(np.random.rand(protein_count, run_count), trainable = True, dtype = 'float32')",
                  TrsfrCoeff[(w3+11L):length(TrsfrCoeff)])
  #
  write(TrsfrCoeff, paste0(pepDir, "/TransferCoeff.py"))
  py_clear_last_error()
  cat(" - Running coefficients transfer script...\n")
  #saveImgFun(paste0(pepDir, "/Pepper_bckp.RData"))
  #loadFun(paste0(dtst, "/Pepper_bckp.RData"))
  #
  cmd <- paste0("python \"", pepDir, "/TransferCoeff.py\"")
  #cat(cmd)
  tst <- try(system(cmd), silent = TRUE)
  coeffFl <- paste0(pepDir, "/_transfer_inferred_coefficients_run0.tsv")
  if (!tst) { Outcome <- file.exists(modlFl) }
  if (Outcome) { cat("Success!!!\n") } else { cat("Failure!!!!!\n") }
  py_clear_last_error()
  #
  # Apply correction
  # Adj. int. = int. / coeff!!!
  tmp <- read.delim(coeffFl)
  colnames(tmp) <- c("ID", "Pepper coefficient")
  tmp$ID <- as.integer(gsub("^id_", "", tmp$ID))
  tmp$ModSeq <- ev$"Modified sequence"[match(tmp$ID, ev$id)]
  #
  ttl <- "Pepper - distribution of coefficients"
  plot <- ggplot(tmp) + geom_histogram(aes(x = log10(`Pepper coefficient`)), fill = "red", bins = 250L) + 
    scale_y_continuous(expand = c(0L, 0L)) +
    theme_bw() + ggtitle(ttl)
  poplot(plot)
  suppressMessages({
    ggsave(paste0(pepDir, "/", ttl, ".jpeg"), plot, dpi = 150L)
  })
  #
  tmp2 <- ev[[ref]]/tmp$"Pepper coefficient"[match(ev$"Modified sequence", tmp$ModSeq)]
  tmp2 <- tmp2*median(ev[[ref]])/median(tmp2)
  ev[[nuRef]] <- tmp2
  #
  ttl <- "Pepper - corrected intensities distribution"
  tmp <- rbind(data.frame("Modified sequence"= ev$"Modified sequence",
                          "Intensity" = ev[[ref]],
                          "Type" = "Original",
                          check.names = FALSE),
               data.frame("Modified sequence"= ev$"Modified sequence",
                          "Intensity" = ev[[nuRef]],
                          "Type" = "Pepper-corrected",
                          check.names = FALSE))
  plot <- ggplot(tmp) + geom_density(stat = "density", aes(x = log10(`Intensity`), fill = Type), alpha = 0.3) +
    scale_y_continuous(expand = c(0L, 0L)) +
    theme_bw() + ggtitle(ttl)
  poplot(plot)
  suppressMessages({
    ggsave(paste0(pepDir, "/", ttl, ".jpeg"), plot, dpi = 150L)
  })
  #
  # Test results ourselves
  g <- grep(";", ev$Proteins, invert = TRUE)
  tmp <- ev[g, c("Modified sequence", "Proteins", "Experiment", ref, nuRef)]
  colnames(tmp) <- c("ModSeq", "Proteins", "Experiment", "Original", "Pepper")
  tmp <- as.data.table(tmp)
  obj <- setNames(c("Original", "Pepper"), c("Before", "After"))
  call <- paste0(vapply(Exp, \(x) { paste0(x, " = sd(", x, ")") }, ""), collapse = ", ") # start building our call
  for (i in 1L:2L) {
    tmp2 <- copy(tmp)
    tmp2 <- dcast(tmp2, ModSeq + Proteins ~ Experiment, value.var = obj[i], fun.aggregate = sum, na.rm = TRUE)
    tmp2 <- as.data.frame(tmp2)
    w <- which(tmp2[, Exp] == 0, arr.ind = TRUE)
    tmp2[, Exp][w] <- NA
    tmp2 <- as.data.table(tmp2)
    call2 <- paste0("tmp2 <- tmp2[, list(", call, "), by = list(Proteins = Proteins)]")
    eval(parse(text = call2), envir = .GlobalEnv)
    tmp2 <- as.data.frame(tmp2)
    assign(names(obj)[i], tmp2)
    tmp2 <- colSums(tmp2[, Exp], na.rm = TRUE)
    assign(paste0(names(obj)[i], "_sums"), tmp2)
  }
  msg <- paste0("Pepper:\n\n#######\n   Final improvement per sample:\n\n\t -> Sum of intra-protein SDs, ratio before/after\n\t\t",
                paste(Exp, collapse = "\t"), "\n\t\t",
                paste(round(Before_sums/After_sums, 1L), collapse = "\t"), "\n\n\t -> Overall: ",
                round(sum(is.all.good(unlist(Before[, Exp])))/sum(is.all.good(unlist(After[, Exp]))), 1L), "\n\n")
  cat(msg)
  write(msg, paste0(pepDir, "/Final outcome.txt"))
}
