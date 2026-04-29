############################
#                          #
# Run an F-test with limma #
#                          #
############################

if (!require(limma)) { pak::pak("limma") }
library(limma)

# Check our parent cluster
source(parSrc, local = FALSE)

# Argument names
allArgs <- c("myData",
             "intRef",
             "ratRef",
             "namesCol",
             "namesRoot",
             "runLoc",
             "ohDeer",
             "idCol",
             "protCol",
             "mtchCol",
             "plotlyLab",
             "Alpha",
             "refRat")
w <- which(allArgs %in% .obj)
if (length(w)) {
  warning(paste0("Source-specific objects ", paste(allArgs[w], collapse = "/"),
                 " overlap with reserved object names in this workflow and should be renamed!"))
}
limpaMode <- FALSE
labCol <- c("Label", "Labels")
if (dataType == "modPeptides") {
  intRef <- ptms.ref[length(ptms.ref)]
  ratRef <- ptms.ratios.ref
  namesCol <- "Name"
  namesRoot <- "Pep"
  runLoc <- LocAnalysis
  ohDeer <- paste0(wd, "/Reg. analysis/", ptm, "/ANOVA")
  idCol <- "Name"
  protCol <- "Proteins"
  mtchCol <- "Modified sequence"
  plotlyLab = c(PepLabKol, paste0(Ptm, "-site"))
  Size <- "Rel. av. log10 abundance"
  Alpha <- "-log10(PEP)"
  labCol <- intersect(labCol, colnames(ptmpep))
  allKol <- unique(c("id", idCol, namesCol, "Genes", Size, Alpha, plotlyLab,
                     "Potential contaminant", labCol, mtchCol))
  allKol <- intersect(allKol, colnames(ptmpep))
  g <- grep(topattern(intRef), colnames(ptmpep), value = TRUE)
  myData <- ptmpep[, c(allKol, g)]
  ptmpep[, g] <- log2(ptmpep[, g]) # Because limma expects log2 expression data and for peptide input data isn't log-transformed!
  myData$"-log10(PEP)" <- -log10(myData$PEP)
}
if (dataType == "PG") {
  intRef <- Prot.Expr.Root
  ratRef <- Prot.Rat.Root
  idCol <- protCol <- "Protein IDs"
  namesCol <- "Leading protein IDs"
  mtchCol <- "Protein IDs"
  plotlyLab <- PrLabKol
  Size <- "Rel. av. log10 abundance"
  Alpha <- "Rel. log10(Peptides count)"
  labCol <- intersect(labCol, colnames(PG))
  allKol <- unique(c("id", idCol, namesCol, "Genes", Size, Alpha, plotlyLab,
                     "Potential contaminant", labCol, mtchCol))
  allKol <- intersect(allKol, colnames(PG))
  if (quantAlgo == "limpa") {
    limpaMode <- TRUE
    myData <- quantData_list$EList_obj
    w <- which(rownames(myData$genes) %in% PG[[mtchCol]])
    m <- match(rownames(myData$genes)[w], PG[[mtchCol]])
    myData$genes[, allKol] <- NA
    myData$genes[w, allKol] <- PG[m, allKol]
    colnames(myData) <- sub(".* - ", intRef, colnames(myData))
  } else {
    g <- grep(topattern(intRef), colnames(PG), value = TRUE)
    myData <- PG[, c(allKol, g)]
    myData[, g] <- myData[, g]/log10(2L) # Because limma expects log2 expression data
  }
  namesRoot <- "PG"
  runLoc <- FALSE
  ohDeer <- paste0(wd, "/Reg. analysis/F-tests")
}
refRat <- NULL
BH.FDR_F <- sort(BH.FDR, decreasing = TRUE)
#
if (!dir.exists(ohDeer)) { dir.create(ohDeer, recursive = TRUE) }

Param$Ratios.Contamination.Rates <- abs(as.numeric(Param$Ratios.Contamination.Rates))
if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
  lfcThr <- Param$Ratios.Contamination.Rates
}
if (Param$Ratios.Thresholds == threshMsg) {
  stop("This option is deprecated!")
}
if (lfcThr != 0L) {
  warning("Please note that limma allows applying a non-null lfc threshold, but does not recommend to do it (see ?decideTests for a discussion).")
}

# Step 1 - run F-test
# -------------------
# voomaLmFit() isn't compatible with duplicateCorrelation()
# For an ANOVA, I think controlling for variance is more important than blocks, so I will go for voomaLmFit() instead of lmFit() with optional duplicateCorrelation()
if (limpaMode) {
  voomFit <- voomaLmFit(myData, designMatr, keep.EList = TRUE)
  NA_Filt <- 1L:nrow(myData)
} else {
  kol <- paste0(intRef, rownames(expMap))
  kol <- intersect(kol, colnames(myData))
  tmpVal <- myData[, kol]
  colnames(tmpVal) <- rownames(designMatr)[match(rownames(designMatr), gsub("___", "_", expMap[[RSA$limmaCol]]))]
  NA_Filt <- as.data.frame(sapply(1L:ncol(designMatr), \(x) { #x <- 1L
    kols <- row.names(designMatr)[which(designMatr[, x] == 1L)]
    apply(tmpVal[, kols, drop = FALSE], 1L, \(y) { length(is.all.good(y)) }) >= 2L
  }))
  NA_Filt <- which(rowSums(NA_Filt) == ncol(NA_Filt)) # We will apply this result to decideTests
  voomFit <- voomaLmFit(tmpVal[NA_Filt, ], designMatr, keep.EList = TRUE)
  voomFit$genes <- myData[NA_Filt, namesCol] # For convenience
}
voomFit <- eBayes(voomFit) # Run F-test
fit_postHoc <- contrasts.fit(voomFit, contrMatr)
fit_postHoc <- eBayes(fit_postHoc)
if (limpaMode) {
  my_F_Data <- myData$genes
} else {
  my_F_Data <- myData[, unique(c(namesCol, idCol, protCol, mtchCol))]
  my_F_Data[[F_Root]] <- NA
}
my_F_Data[NA_Filt, F_Root] <- -log10(voomFit$F.p.value)
# Do NOT:
#  - use contrasts.fit() before eBayes() for the F-test!
#  - make limma tests artificially one-sided; instead:
#     - apply FDR correction on two-sided limma tests
#     - filter for logFCs using only biologically relevant side of test
#
# Step 2 - run post-hoc tests (global moderated t-tests with contrasts)
# ---------------------------------------------------------------------
F_PVal_postHoc <- paste0(F_Root, " - ", myContrasts$Contrast)
my_F_Data[, F_PVal_postHoc] <- NA
my_F_Data[NA_Filt, F_PVal_postHoc] <- -log10(fit_postHoc$p.value)

# Step 3 - make decision for each FDR value
# -----------------------------------------
# The decision is based on:
#  - The global F-test
#  - Individual t-tests
#  - Individual fold changes
# 
# a) Significance columns for the F-test
# ......................................
fdrKol <- setNames(paste0("mod. F-test Significant-FDR=", BH.FDR_F*100, "%"),
                   BH.FDR_F*100)
F_fdr <- list()
F_fdr$F_test <- FDR(my_F_Data,
                    pvalue_col = F_Root,
                    fdr = BH.FDR_F,
                    returns = c(TRUE, TRUE, FALSE))
my_F_Data[, fdrKol] <- F_fdr$F_test$`Significance vector`
#
# b) Significance columns for each post-hoc test
# ..............................................
# Note that this could also be done globally. But usually here people consider each test individually.
# The decision as to whether to calculate FDR thresholds for each post-hoc test indiviually,
# or globally for all, could be parameter controlled in the future.
tmp <- my_F_Data[, F_PVal_postHoc, drop = FALSE]
globalFDR <- FALSE
if (globalFDR) {
  stop("Check that bit first!")
  tmp2 <- melt(tmp)
  tmp3 <- FDR(tmp2,
              pvalue_col = "value",
              fdr = BH.FDR_F,
              returns = c(TRUE, TRUE, FALSE))
  tmp3$`Significance vector`[1:10,]
  n <- length(F_PVal_postHoc)
  m <- nrow(tmp3$`Significance vector`)/n
  F_fdr[myContrasts$Contrast] <- lapply(1L:n, \(i) {
    x <- tmp3
    x$`Significance vector` <- x$`Significance vector`[(1L:m) + (i-1L)*m,]
    return(x)
  })
} else {
  tmpFl <- tempfile(fileext = ".rds")
  readr::write_rds(tmp, tmpFl)
  clusterExport(parClust, list("F_Root", "BH.FDR_F", "tmpFl", "FDR", "is.all.good"))
  invisible(clusterCall(parClust, \() {
    tmp <- readr::read_rds(tmpFl)
    assign("tmp", tmp, .GlobalEnv)
    return()
  }))
  unlink(tmpFl)
  F_fdr[myContrasts$Contrast] <- parLapply(parClust, myContrasts$Contrast, \(i) { #i <- myContrasts$Contrast[1L]
    F_pv_i <- paste0(F_Root, " - ", i)
    return(FDR(tmp,
               pvalue_col = F_pv_i,
               fdr = BH.FDR_F,
               returns = c(TRUE, TRUE, FALSE)))
  })
}
fdrKol_contr <- c()
for (i in myContrasts$Contrast) {
  fdrKol_i <- paste0(fdrKol, " - ", i)
  fdrKol_contr <- c(fdrKol_contr, fdrKol_i)
  my_F_Data[, fdrKol_i] <- F_fdr[[i]]$`Significance vector`
}
#
# c) Final decision, based on the prior Significance columns + logFC
# We define a single logFC threshold for the whole test
# First get log2FC columns
xCol <- paste0(ratRef, myContrasts$Contrast)
my_F_Data[NA_Filt, xCol] <- fit_postHoc$coefficients[, myContrasts$Contrast]
regRoot_F <- "mod. F-test Regulated - "
regKol <- paste0(regRoot_F, myContrasts$Contrast)
my_F_Data[, regKol] <- "non significant"
for (fdr in BH.FDR_F) { #fdr <- BH.FDR_F[1L]
  fdrkl1 <- fdrKol[as.character(fdr*100)]
  tst1 <- my_F_Data[[fdrkl1]] == "+"
  for (i in 1L:nrow(myContrasts)) { #i <- 1L
    contr <- myContrasts$Contrast[[i]]
    upOnly <- myContrasts$`Up-only`[[i]]
    rgkl <- paste0(regRoot_F, contr)
    fdrkl2 <- paste0(fdrKol[as.character(fdr*100)], " - ", contr)
    tst2 <- my_F_Data[[fdrkl2]] == "+"
    fckl <- paste0(ratRef, contr)
    tst3 <- my_F_Data[[fckl]] >= 0
    tst4 <- abs(my_F_Data[[fckl]]) < lfcThr
    wUp <- which(tst1&tst2&tst3)
    my_F_Data[wUp, rgkl] <- paste0("up, FDR = ", fdr*100, "%")
    if (!upOnly) {
      wDown <- which(tst1&tst2&!tst3)
      my_F_Data[wDown, rgkl] <- paste0("down, FDR = ", fdr*100, "%")
    }
    wNS <- which(tst1&tst2&tst4)
    my_F_Data[wNS, rgkl] <- "too small FC"
  }
}
#View(my_F_Data)
kol <- setdiff(allKol, colnames(my_F_Data))
if (length(kol)) { my_F_Data[, kol] <- myData[, kol] }
#
# Volcano plots
my_F_Data[[Param$Plot.labels]] <- my_F_Data[[labCol[1L]]]
volcPlot_args2 <- volcPlot_args
volcPlot_args2$Prot <- my_F_Data
volcPlot_args2$X.root <- ratRef
volcPlot_args2$Y.root <- paste0(F_Root, " - ")
volcPlot_args2$FDR.root <- "mod. F-test Significant-FDR="
volcPlot_args2$IDs.col <- idCol
volcPlot_args2$Proteins.col <- protCol
volcPlot_args2$subfolder <- ohDeer
volcPlot_args2$reg.root <- regRoot_F
volcPlot_args2$title <- "F_test volcano plot "
volcPlot_args2$cl <- parClust
# For testing:
#DefArg(Volcano.plot)
#invisible(lapply(names(volcPlot_args2), \(x) { assign(x, volcPlot_args2[[x]], envir = .GlobalEnv); return() }))
#TESTING <- TRUE
F_volc <- do.call(Volcano.plot, volcPlot_args2)
stopCluster(parClust)
source(parSrc)
#
# Save plotly plots
dr <- ohDeer
myPlotLys <- F_volc$`Plotly plots`
Src <- paste0(libPath, "/extdata/Sources/save_Plotlys.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
F_thresh <- F_volc$Thresholds
thresh <- lapply(names(F_thresh$Absolute), \(x) {
  y <- F_thresh$Absolute[[x]]
  x <- data.frame(Group = rep(cleanNms(x), nrow(y)))
  return(cbind(x, y))
})
thresh <- plyr::rbind.fill(thresh)
thresh$Name <- NULL
thresh$Root <- gsub(" - $", "", thresh$Root)
thresh$Value <- thresh$Text.value
thresh$Text.value <- NULL
# FDR thresholds
fdrThresh <- F_thresh$FDR
if (!is.null(fdrThresh)) {
  colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.up")] <- "Colour (up)"
  colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.down")] <- "Colour (down)"
  colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.line")] <- "Colour (line)"
  if ("Test" %in% colnames(fdrThresh)) {
    fdrThresh <- fdrThresh[, c("Test", colnames(fdrThresh)[which(colnames(fdrThresh) != "Test")])]
  }
}
fl <- paste0(ohDeer, "/Thresholds.xlsx")
wb <- wb_workbook()
wb <- wb_set_creators(wb, "Me")
wb <- wb_add_worksheet(wb, "Thresholds")
dms <- wb_dims(2L, 1L)
wb <- wb_add_data(wb, "Thresholds", "Absolute thresholds", dms)
wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                  italic = "true", underline = "single")
dms <- wb_dims(3L, 2L)
wb <- wb_add_data_table(wb, "Thresholds", thresh, dms,
                        col_names = TRUE, table_style = "TableStyleMedium2",
                        banded_rows = TRUE, banded_cols = FALSE)
dms <- wb_dims(nrow(thresh)+6L, 1L)
wb <- wb_add_data(wb, "Thresholds", "FDR thresholds", dms)
wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                  italic = "true", underline = "single")
dms <- wb_dims(nrow(thresh)+7L, 2L)
if (!is.null(fdrThresh)) {
  wb <- wb_add_data_table(wb, "Thresholds", fdrThresh, dms,
                          col_names = TRUE, table_style = "TableStyleMedium2",
                          banded_rows = TRUE, banded_cols = FALSE)
  wb <- wb_set_col_widths(wb, "Thresholds", 1L, 3L)
  tmp1 <- rbind(colnames(thresh), thresh)
  colnames(tmp1) <- paste0("V", 1L:ncol(tmp1))
  tmp2 <- rbind(colnames(fdrThresh), fdrThresh)
  colnames(tmp2) <- paste0("V", 1L:ncol(tmp2))
  tst <- plyr::rbind.fill(tmp1, tmp2)
  tst <- setNames(apply(tst, 2L, \(x) { max(nchar(x), na.rm = TRUE) }), NULL)
  wb <- wb_set_col_widths(wb, "Thresholds", 1L:(length(tst)+1L), c(3L, tst))
}
wb_save(wb, fl)
#xl_open(fl)
#
# Final data to be exported to PG or pep
if (limpaMode) {
  myData <- myData$genes
}
myData[, colnames(my_F_Data)] <- my_F_Data
# Remove log2-transformed data from myData to avoid overwriting quant data in PG or pep
kol <- grep(topattern(intRef), colnames(myData), value = TRUE, invert = TRUE)
#
F_kols <- c(grep(topattern(F_Root), colnames(myData), value = TRUE),
            grep(topattern(regRoot_F), colnames(myData), value = TRUE))
if (dataType == "modPeptides") {
  w <- which(ptmpep[[mtchCol]] %in% myData[[mtchCol]])
  #which(!ptmpep[[mtchCol]] %in% myData[[mtchCol]])
  m <- match(ptmpep[w, mtchCol], myData[[mtchCol]])
  ptmpep[, F_kols] <- NA
  ptmpep[w, F_kols] <- myData[m, F_kols]
  PTMs_F_test_data[[Ptm]] <- my_F_Data
}
if (dataType == "PG") {
  w <- which(PG[[mtchCol]] %in% myData[[mtchCol]])
  #which(!PG[[mtchCol]] %in% myData[[mtchCol]])
  m <- match(PG[w, mtchCol], myData[[mtchCol]])
  PG[, F_kols] <- NA
  PG[w, F_kols] <- myData[m, F_kols]
  F_test_data %<o% my_F_Data
}
#
# Also a posteriori F-test P-values histogram:
nbin <- ceiling(max(c(20L, 10L*round(nrow(myData)/1000L))))
bd <- (0L:nbin)/nbin
if (F_Root %in% colnames(myData)) {
  temp <- data.frame(value = is.all.good(10L^(-myData[[F_Root]])))
  ttl <- "Histogram: F-test moderated Pvalue"
  plot <- ggplot(temp, aes(x = value)) +
    geom_histogram(bins = nbin, colour = "black", alpha = 0.25, fill = "green") +
    guides(fill = "none") + theme_bw() + ggtitle(ttl)
  #poplot(plot)
  ttla <- gsub(": ?", " - ", ttl)
  suppressMessages({
    ggsave(paste0(ohDeer, "/", ttla, ".jpeg"), plot, dpi = 300L)
    ggsave(paste0(ohDeer, "/", ttla, ".pdf"), plot, dpi = 300L)
  })
  ReportCalls <- AddPlot2Report(Title = ttla)
}
# Cleanup
for (i in allArgs) { try(rm(i), silent = TRUE) }
