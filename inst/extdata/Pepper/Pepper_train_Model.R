#
if (!require(svDialogs)) { install.packages("svDialogs") }
library(svDialogs)
if (!require(qs)) { install.packages("qs") }
library(qs)
if (!require(data.table)) { install.packages("data.table") }
library(data.table)
if (!require(parallel)) { install.packages("parallel") }
library(parallel)
nClst <- detectCores()-1

# Create parallel processing cluster
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(stopCluster(parClust), silent = TRUE)
  parClust <- makeCluster(nClst, type = "SOCK")
}
try(setDTthreads(threads = nClst), silent = TRUE)

# Peptide directory
homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
PepScrptsDir <- paste0(RPath, "/proteoCraft/extdata/Pepper")
stopifnot(dir.exists(PepScrptsDir))
pyInit <- paste0(PepScrptsDir, "/Python_start.R")
source(pyInit)



# Detect recommended training parameters for all datasets
trainParam <- readLines(paste0(PepScrptsDir, "/train_models.sh"))
args <- c("peptide_file", "n_runs", "seq_length", "output_file", "filter_size", "n_filters", "n_layers", "n_nodes",
          "dropout", "learning_rate", "batch", "random_run")
g <- grep("^python3 ", trainParam)
l <- length(trainParam)
trainParam <- sapply(args, function(arg) {
  sapply(1:length(g), function(x) {
    h <- grep(paste0("--", arg), trainParam)
    h <- gsub(paste0(".*'--", arg, "' +"), "", trainParam[h[which(h %in% g[x]:c(g, l)[x+1])]])
    pat1 <- paste0("^", substr(h, 1, 1))
    pat2 <- paste0(substr(h, 1, 1), ".*$")
    h <- gsub(pat2, "", gsub(pat1, "", h))
    return(h)
  })
})
trainParam <- as.data.frame(trainParam)
rownames(trainParam) <- gsub(".+/|\\.tsv$", "", trainParam$peptide_file)

# Hyper-parameter tuning on one of our own datasets
# a) Prepare one of our own datasets for one-hot-encoding
#evFl <- normalizePath(choose.files(paste0(homePath, "/evidence.tsv")), winslash = "/")
evFl <- rstudioapi::selectFile("Select PSMs tsv file, as saved by one of our data analysis workflows in subfolder .../Tables/",
                               path = paste0(homePath, "/*.tsv"),
                               filter = "PSMs tsv file (*.tsv)")
dtst <- gsub("/Tables$", "", dirname(evFl))
dtstnm <- gsub(".*/", "", dtst)
dtstnm <- dlg_input("Enter a name for this dataset", dtstnm)$res
msg <- "Select experiment/samples .csv map"
#mapFl <- normalizePath(choose.files(paste0(dtst, "/Experiment map.csv"), msg), winslash = "/")
mapFl <- rstudioapi::selectFile(msg,
                                path = paste0(dtst, "/Experiment map.csv"),
                                filter = "csv file (*.csv)")
msg <- "Select Modifications.csv file"
#modFl <- normalizePath(choose.files(paste0(dtst, "/Workflow control/Modifications.csv"), msg), winslash = "/")
modFl <- rstudioapi::selectFile(msg,
                                path = paste0(dtst, "/Workflow control/Modifications.csv"),
                                filter = "csv file (*.csv)")
ev <- data.table::fread(evFl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
smplsMap <- read.csv(mapFl, check.names = FALSE)
if ("Ref.Sample.Aggregate" %in% colnames(smplsMap)) {
  smplKol <- "RSA"
  smplKol2 <- "Ref.Sample.Aggregate"
  if ("RSA" %in% colnames(ev)) {
    smplKol <- "RSA"
  } else {
    tmp <- listMelt(smplsMap$MQ.Exp, smplsMap$Ref.Sample.Aggregate)
    ev$RSA <- tmp$L1[match(ev$MQ.Exp, tmp$value)]
  }
} else {
  smplKol2 <- dlg_list(colnames(smplsMap), title = "Select parent samples column")$res
  cat(paste0("Samples:\n - ", paste0(unique(smplsMap[[smplKol2]]), collapse = "\n - "), "\n"))
  w <- which(ev[1,] %in% smplsMap[[smplKol2]])
  stopifnot(length(w) > 0)
  smplKol <- smplKol2
  if (length(w) > 1) {
    if (smplKol2 %in% colnames(ev)[w]) {
      smplKol <- smplKol2
    } else {
      smplKol <- dlg_list(colnames(ev)[w], title = "I'm not sure which column to use for samples in the evidence file, could you clarify for me?")$res
    }
  }
}
Modifs <- read.csv(modFl, check.names = FALSE)
Modifs$AA <- strsplit(Modifs$AA, ", *")
#
pepDir <- paste0(dtst, "/Pepper")
if (!dir.exists(pepDir)) { dir.create(pepDir, recursive = TRUE) }
tempEv <- Pepper_TrainingData(ev, Modifs, smplsMap, smplKol, smplKol2, path = paste0(pepDir, "/ModPep4Pepper.tsv"))
Modifs <- tempEv$PTMs
tempEv <- tempEv$Data
#View(tempEv)
#
# Edit and run python script for One-Hot-Encoding
#################################################
Grps <- colnames(tempEv)[(max(grep("^Charge [0-9]+$", colnames(tempEv)))+1):length(colnames(tempEv))]
lGrps <- length(Grps)
tmp <- data.frame(Peptides = nrow(tempEv), Proteins = length(unique(tempEv$Protein)), Samples = lGrps)
write.table(tmp, paste0(pepDir, "/Dataset dimensions.tsv"), sep = "\t", row.names = FALSE)
cat(paste0("Dataset = ", dtstnm, ", number of samples = ", lGrps, "\n"))
#
OHfl <- paste0(pepDir, "/OneHotEncodedPepQuants.tsv")
OHE <- readLines(paste0(PepScrptsDir, "/onehot_encode_peptide_sequences.py"))
arg1 <- which(OHE == "input_filename = sys.argv[1]")
arg2 <- which(OHE == "runs = int(sys.argv[2])")
arg3 <- which(OHE == "output_filename = sys.argv[3]")
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
  CleanUp <- opt[dlg_list(names(opt), names(opt)[1],
                          title = paste0(dtstnm, ": shall we start from scratch hyperparameter tuning?"))$res]
  names(CleanUp) <- NULL
}
# NB: This code is written in a way that it will always sample Step1Tries + Step2Tries parameter combinations.
# 1 - randomly
Step <- 1
paramKol <- c("n_conv_layers", "n_filters", "filter_size", "n_layers", "n_nodes", "dropout", "batch", "learning_rate")
if (CleanUp) { for (fl in GrdSrchFls) { unlink(fl) } } # Cleanup for fresh tuning
HypTunFl <- paste0(PepScrptsDir, "/gridseach_parameters_neural_network.py")
stopifnot(file.exists(HypTunFl))
HypTunSrc <- paste0(PepScrptsDir, "/HyperTun.R")
source(HypTunSrc)
L <- length(GrdSrchFls)
N <- Step1Tries - L
kount <- 0
while ((N > 0)&&(kount < N+1)) {
  clusterExport(parClust, "cmd", envir = environment())
  # Below, we are never running less tests than the number of cores, so we do not waste a round on a few tests.
  tst <- parSapply(parClust, 1:max(c(nClst, N)), function(x) {
    try(system(cmd), silent = TRUE)
  })
  GrdSrchFls <- list.files(tuneDir, "Gridsearch_Results_", full.names = TRUE)
  L <- length(GrdSrchFls)
  N <- Step1Tries - L
  kount <- kount + 1
}
Step1Tries <- L # Update it if we got any bonus tests
GrdSrch <- as.data.frame(t(sapply(strsplit(gsub("\\.tsv$", "", GrdSrchFls), "_"), function(x) {
  x <- unlist(x)
  x <- suppressWarnings(as.numeric(x))
  x <- x[which(!is.na(x))]
  return(x)
})))
colnames(GrdSrch) <- paramKol
tst <- apply(GrdSrch, 2, function(x) { length(unique(x)) })
GrdSrch <- GrdSrch[, which(tst > 1), drop = FALSE]
paramKol2 <- paramKol[which(tst > 1)]
rs <- plyr::rbind.fill(lapply(GrdSrchFls, read.delim))
rs$X <- NULL
GrdSrch[, gsub("\\.", "_", colnames(rs))] <- rs
GrdSrch <- GrdSrch[order(GrdSrch$Val_percent_improvement, decreasing = TRUE),]
bstImprov <- max(GrdSrch$Val_percent_improvement)
# 2 - explore immediate neighborhood of best 10 results
Step <- 2
currParam <- GrdSrch[1:10, c(paramKol2, "No_of_epochs")]
source(HypTunSrc)
L <- length(GrdSrchFls)
N <- Step1Tries + Step2Tries - L
kount <- 0
while ((N > 0)&&(kount < N+1)) {
  clusterExport(parClust, "cmd", envir = environment())
  # Below, we are never running less tests than the number of cores, so we do not waste a round on a few tests.
  tst <- parSapply(parClust, 1:max(c(nClst, N)), function(x) {
    try(system(cmd), silent = TRUE)
  })
  GrdSrchFls <- list.files(tuneDir, "Gridsearch_Results_", full.names = TRUE)
  L <- length(GrdSrchFls)
  N <- Step1Tries + Step2Tries - L
  kount <- kount + 1
}
GrdSrch <- as.data.frame(t(sapply(strsplit(gsub("\\.tsv$", "", GrdSrchFls), "_"), function(x) {
  x <- unlist(x)
  x <- suppressWarnings(as.numeric(x))
  x <- x[which(!is.na(x))]
  return(x)
})))
colnames(GrdSrch) <- paramKol
tst <- apply(GrdSrch, 2, function(x) { length(unique(x)) })
GrdSrch <- GrdSrch[, which(tst > 1), drop = FALSE]
paramKol2 <- paramKol[which(tst > 1)]
rs <- plyr::rbind.fill(lapply(GrdSrchFls, read.delim))
rs$X <- NULL
GrdSrch[, gsub("\\.", "_", colnames(rs))] <- rs
GrdSrch <- GrdSrch[order(GrdSrch$Val_percent_improvement, decreasing = TRUE),]
bstImprov2 <- max(GrdSrch$Val_percent_improvement)
print("")
print(paste0(" - Step 1: best improvement = ", round(bstImprov, 2), "%"))
print("")
print(paste0(" - Step 2: best improvement = ", round(bstImprov2, 2), "%"))
print("")
#View(GrdSrch)
#
# Get final, best parameters
bstParam <- GrdSrch[1, c(paramKol2, "No_of_epochs")]
print(bstParam)
bstParam$peptide_file <- OHfl
bstParam$n_runs <- as.character(length(Grps))
bstParam$seq_length <- "60"
bstParam$random_run <- "0" # Obviously: we will do a proper training now!
bstParam$output_file <- dtstnm
# Create model, using our own dataset and the above hypertuned parameters
TrCoeff <- readLines(paste0(PepScrptsDir, "/peptide_coefficient_predictor.py"))

# Parameters from train_models.sh"
w <- which(TrCoeff %in% c("parser = argparse.ArgumentParser()", "args = parser.parse_args()", "print(args)"))
TrCoeff[w] <- ""
basePat <- paste0("^ *", topattern("parser.add_argument(", start = FALSE), " *")
for (arg in colnames(trainParam)) { #arg <- colnames(trainParam)[1]
  pat <- paste0(basePat, "\'--", arg)
  g <- grep(pat, TrCoeff)
  stopifnot(length(g) == 1)
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
  for (g in G) { #g <- G[1]
    tmp <- TrCoeff[g]
    arg <- gsub("\'.*", "", gsub(paste0(basePat, "\'--"), "", tmp))
    stopifnot(nchar(arg) > 0)
    val <- unlist(strsplit(tmp, "[,\\(] *default *= *"))
    if (!length(val) == 2) { stop(paste0("No default for argument ", arg, "! Give me a hand here, human!")) }
    val <- gsub(" *[,\\)].*", "", val[2])
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
TrCoeff <- c(TrCoeff[1:w],
             "#################################################################################################",
             "### Insertion to remove cases where the above lines generate run indices outside of the range ###",
             "#################################################################################################",
             "train_runs = train_runs[train_runs+1 <= q_df.shape[1]]",
             "test_runs = test_runs[test_runs+1 <= q_df.shape[1]]",
             "#################################################################################################",
             "### Insertion ends                                                                            ###",
             "#################################################################################################",
             TrCoeff[(w+1):length(TrCoeff)])
# Add code to save final layer dimensions
w <- which(TrCoeff == "print(\"Saved model to disk\")")
TrCoeff <- c(TrCoeff[1:w],
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
             TrCoeff[(w+1):length(TrCoeff)])
write(TrCoeff, paste0(pepDir, "/TrainCoeff.py"))

py_clear_last_error()
cat(" - Running coefficients training script...\n")
cmd <- paste0("python \"", pepDir, "/TrainCoeff.py\"")
#cat(cmd)
tst <- try(system(cmd), silent = TRUE)
stopifnot(tst == 0)
#
modlFl <- paste0(pepDir, "/", dtstnm, "_Coefficient_Predictor_Model.h5")
if (tst == 0) { Outcome <- file.exists(modlFl) }
if (Outcome) { cat("Success!!!\n") } else { cat("Failure!!!!!\n") }
py_clear_last_error()

# Store for later use
saveModel <- c(TRUE, FALSE)[match(dlg_message("Should we store this model for later use?", "yesno")$res, c("yes", "no"))]
if (saveModel) {
  if (!require(shiny)) { install.packages("shiny") }
  library(shiny)
  if (!require(DT)) { install.packages("DT") }
  library(DT)
  modDir <- paste0(homePath, "/Pepper_models")
  if (!dir.exists(modDir)) { dir.create(modDir, recursive = TRUE) }
  file.copy(modlFl, modDir)
  modDatFl <- paste0(modDir, "/Models_data.csv")
  kol <- c("Species", "Column", "Sample prep", "Acquisition method", "Search software")
  new <- !file.exists(modDatFl)
  if (!new) {
    modDat <- read.csv(modDatFl, check.names = FALSE)
    tst <- sum(!c("Name", kol) %in% colnames(modDat))
    if (tst) { new <- TRUE }
  }
  if (new) {
    modlFls <- list.files(modDir, "\\.h5$")
    modDat <- data.frame("Name" = modlFls,
                         "Species" = NA,
                         "Column" = NA,
                         "Sample prep" = NA,
                         "Acquisition method" = NA,
                         "Search software" = NA,
                         check.names = FALSE)
  }
}
