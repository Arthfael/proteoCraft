############################
#                          #
# Run an F-test with limma #
#                          #
############################

if (!require(limma)) { pak::pak("limma") }
library(limma)

# Check our parent cluster
source(parSrc, local = FALSE)
parallel::clusterExport(parClust, list("AltHyp", "Nested"), envir = environment())
invisible(parallel::clusterCall(parClust, function() { library(siggenes); return() }))

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

if (dataType == "modPeptides") {
  myData <- ptmpep
  intRef <- pepRf
  ratRef <- pepRatRf
  namesCol <- "Name"
  namesRoot <- "Pep"
  runLoc <- LocAnalysis
  ohDeer <- paste0(wd, "/Reg. analysis/", ptm, "/F-tests")
  idCol <- "Name"
  protCol <- "Proteins"
  mtchCol <- "Modified sequence"
  plotlyLab = c(PepLabKol, paste0(Ptm, "-site"))
  myData$"-log10(PEP)" <- -log10(myData$PEP)
  Size <- "Rel. av. log10 abundance"
  Alpha <- "-log10(PEP)"
  refRat <- ref.rat
}
if (dataType == "PG") {
  myData <- PG
  intRef <- Prot.Expr.Root
  ratRef <- Prot.Rat.Root
  namesCol <- mtchCol <- "Leading protein IDs"
  namesRoot <- "PG"
  runLoc <- FALSE
  ohDeer <- paste0(wd, "/Reg. analysis/F-tests")
  idCol <- protCol <- "Protein IDs"
  plotlyLab <- PrLabKol
  Size <- "Rel. av. log10 abundance"
  Alpha <- "Rel. log10(Peptides count)"
  refRat <- Ref.Ratios
}
labKol <- c("Label", "Labels")
labKol <- labKol[which(labKol %in% colnames(myData))]
if ((length(labKol))&&(!Param$Plot.labels %in% colnames(myData))) {
  myData[[Param$Plot.labels]] <- myData[[labKol[1]]]
}
if (!dir.exists(ohDeer)) { dir.create(ohDeer, recursive = TRUE) }
#
# First check that we have means per group
expMap_F <- Exp.map[match(rownames(designMatr), Exp.map$Ref.Sample.Aggregate),]
# Replace hyphens by dots to avoid issues with evaluating contrasts
for (nuCoeff in Coefficients) {
  Coeff <- gsub("___$", "", nuCoeff)
  stopifnot(Coeff %in% colnames(Exp.map))
  l <- length(grep("-", Exp.map[[Coeff]])) # Not expMap_F in case we are re-running a small chunk
  if (l) {
    stopifnot(!nuCoeff %in% colnames(Exp.map)) # Not expMap_F in case we are re-running a small chunk
    expMap_F[[nuCoeff]] <- gsub("-", ".", expMap_F[[Coeff]])
  }
}
Group_ <- do.call(paste, c(expMap_F[, Coefficients, drop = FALSE], sep = "_"))
Group_ <- as.factor(Group_)
expMap_F$Group_ <- Group_
#
expContr_F <- expContrasts_F
expContr_F$Map <- lapply(expContr_F$All, function(x) { #x <- expContr_F$All[[1]]
  em <- expMap_F[which(expMap_F$Group_ %in% unlist(x)),]
  em$Expression_Column <- paste0(intRef, em$Ref.Sample.Aggregate)
  em$Ratios_Column <- paste0(ratRef, em$Ref.Sample.Aggregate)
  em[which(em$Expression_Column %in% colnames(myData)),]
  em <- em[c(which(em$Reference),
             which(!em$Reference)),]
  return(em)
})
contrCall_F <- paste0("contrMatr_F %<o% makeContrasts(",
                      paste(expContr_F$Contrasts, collapse = ", "),
                      ", levels = designMatr)")
#cat(contrCall_F, "\n")
eval(parse(text = contrCall_F), envir = .GlobalEnv)

# Average LFCs
# (Calculate average ratios)
whSimple <- which(expContr_F$Type == "Simple") # One way contrasts
tmp <- lapply(whSimple, function(x) { #x <- whSimple[1]
  em <- expContr_F$Map[x][[1]]
  nm <- expContr_F$name[x]
  w0 <- which(em$Reference)
  w1 <- which(!em$Reference)
  w1r <- which((!em$Reference)&(em$Ratios_Column %in% colnames(myData)))
  kr1 <- paste0("Mean ", ratRef, nm) # Eventually we can rename those, either to contrast name or to full contrast (X-Y)
  RES <- data.frame(V1 = rep(NA, nrow(myData)))
  colnames(RES) <- kr1
  if (!length(w1r)) {
    warning("This contrast does not have ratio values at this stage. This calls for a rewrite of how ratios are defined!")
  }
  if (!length(w0)) {
    stop("Invalid contrast! Investigate!")
  }
  #
  # Ratios have two different origins/roles in proteomics data analyses:
  #  - Some software (MaxQuant) can in some cases (SILAC) measure ratios more precisely than intensities.
  #    Thus sometimes software-provided ratios are quantitatively more informative than intensities!
  #    In these cases, this intensity-based workflow should as early as feasible propagate ratios-level information onto quantities:
  #          - if the software provides Ia, Ib and R(a/b), then the new intensities (J) will be defined by Ja+Jb = Ia+Ib and Ja/Jb = R(a/b)
  #          - after that, ratios can be ignored again.
  #    For historical reasons (this evolved from a SILAC workflow), this workflow also calculates ratios at different stages.
  #    This should be removed! Only provide average LFC at the end!
  #    Note though that we will still have to measure ref-to-ref or intra-group ratios... unless we also rewrite all of the statistics.
  #  - We also plot a log2 average ratio/fold change (LFCs) in volcano plots. This is the only thing we need.
  #
  # Yet one should never forget that ratios should not be the main focus of analysis. Certainly they are not normal so should not be tested statistically! 
  #
  # As part of rewriting this to focus on intensity data, we will calculate ratios here directly from expression values.
  # We will then be able to progressively remove the redundant ratios from the earlier parts of the workflow.
  #
  if (Nested) {
    rps <- unique(em$Replicate)
    rats <- sapply(rps, function(r) {
      k0 <- em$Expression_Column[which((em$Replicate == r)&(em$Reference))]
      k1 <- em$Expression_Column[which((em$Replicate == r)&(!em$Reference))]
      stopifnot(length(k0) == 1, length(k1) == 1)
      return((myData[[k1]]-myData[[k0]])/log10(2))
    })
    if (length(rps) == 1) {
      RES[[kr1]] <- rats
    } else {
      rats <- cbind(rats)
      rs <- as.data.frame(t(parallel::parApply(parClust, rats, 1, Av_SE_fun)))
      RES[[kr1]] <- rs[, 1]
    }
  } else {
    k1 <- em$Expression_Column[w1]
    k0m <- paste0("Mean ", intRef, unique(em[w0, VPAL$column]))
    rats <- sweep(myData[, k1], 1, myData[[k0m]], "-")/log10(2)
    if (length(w1) == 1) {
      RES[[kr1]] <- rats
    } else {
      rs <- as.data.frame(t(parallel::parApply(parClust, rats, 1, Av_SE_fun)))
      RES[[kr1]] <- rs[, 1]
    }
  }
  return(RES)
})
tmp <- do.call(cbind, tmp)
xCol <- colnames(tmp)
w <- which(!xCol %in% colnames(myData))
if (length(w)) {
  warning("Strange...")
  myData[, xCol[w]] <- tmp[, w]
}
whDouble <- which(expContr_F$Type == "Double") # Two way contrasts
if (length(whDouble)) {
  tmp2 <- lapply(whDouble, function(x) { #x <- whDouble[1]
    contr <- unlist(strsplit(gsub("^\\(|\\)$", "", expContr_F$Contrasts[x]), "\\) - \\("))
    m <- match(contr, expContr_F$Contrasts[whSimple])
    res <- tmp[[m[1]]] - tmp[[m[2]]]
  })
  tmp2 <- do.call(cbind, tmp2)
  colnames(tmp2) <- paste0("Mean ", ratRef, expContr_F$name[whDouble])
  myData[, colnames(tmp2)] <- tmp2
  xCol <- c(xCol, colnames(tmp2))
}
# Ref ratios
if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
  refRat_F <- NULL
}
if (Param$Ratios.Thresholds == threshMsg) {
  refRat_F <- setNames(lapply(1:nrow(expContr_F), function(x) { #x <- 1 #x <- 4
    x <- unique(c(unlist(strsplit(expContr_F$x1[x], " - ")),
                  unlist(strsplit(expContr_F$x0[x], " - "))))
    x <- gsub("^Group_", "", x)
    x <- unique(expMap_F[which(expMap_F$Group_ %in% x), VPAL$column])
    x <- is.all.good(as.numeric(unlist(refRat[x])))
    return(x)
  }), expContr_F$name)
}
#
kol <- paste0(intRef, rownames(designMatr))
wY <- which(kol %in% colnames(myData))
kol <- kol[wY]
tmpVal <- myData[, kol]
colnames(tmpVal) <- rownames(designMatr)[wY]
parallel::clusterExport(parClust, list("tmpVal", "Nested", "AltHyp"), envir = environment())

# - Moderated F-test (limma):
# Do not use duplicateCorrelation(), it is for duplicate rows!!!
## Test for heteroskedasticity - if the est is significant, we use a more robust call 
varMatr <- apply(tmpVal, 1, var, na.rm = TRUE)
meanMatr <- rowMeans(tmpVal, na.rm = TRUE)
varMeanMatr <- data.frame("variance" = varMatr,
                          "mean" = meanMatr,
                          row.names = names(varMatr))
lmMod <- lm(variance ~ mean, data = varMeanMatr)
heteroSkedTst <- lmtest::bptest(lmMod)
## Filter
NA_Filt <- as.data.frame(sapply(1:ncol(designMatr), function(x) { #x <- 1
  kols <- row.names(designMatr)[which(designMatr[,x] == 1)]
  apply(tmpVal[, kols, drop = FALSE], 1, function(y) { length(proteoCraft::is.all.good(y)) }) >= 2
}))
NA_Filt <- which(rowSums(NA_Filt) == ncol(NA_Filt)) # We will apply this result to decideTests
## Fit model
if (heteroSkedTst$p.value < 0.01) {
  warning("Data may not be heteroskedastic, using trend = TRUE and robust = TRUE")
}
fit <- voomaLmFit(tmpVal[NA_Filt,]
                  , designMatr)
fit$genes <- myData[NA_Filt
  , namesCol]
fit <- contrasts.fit(fit, contrMatr_F)
fit <- eBayes(fit)
# Don't use limma.one.sided, it's for the t-test P values!
my_F_Data <- myData[, c(namesCol, idCol, protCol, mtchCol, paste0("Mean ", ratRef, expContr_F$name))]
my_F_Data[[F_Root]] <- NA
my_F_Data[NA_Filt
  , F_Root] <- -log10(fit$F.p.value)
#topTable - not used, using decideTests instead
# tpTbl <- topTable(fit, NULL, nrow(my_F_Data)#length(NA_Filt)
#                   , adjust.method = "none",
#                   sort.by = "F", resort.by = "none",
#                   confint = 0.95)
regRoot_F <- "mod. F-test Regulated - "
regKol <- paste0(regRoot_F, expContr_F$name)
my_F_Data[, regKol] <- "non significant"
# First level: get decision for each FDR value
fdrKol <- paste0("mod. F-test Significant-FDR=", rev(BH.FDR)*100, "%")
my_F_Data[, fdrKol] <- ""
for (fdr in rev(BH.FDR)) { #fdr <- rev(BH.FDR)[1]
  dTsts <- decideTests(fit, c("nestedF", "global")[Nested+1], "none", fdr)
  dTsts <- as.data.frame(dTsts)
  m <- match(rownames(my_F_Data), rownames(dTsts))
  wMtch <- which(!is.na(m))
  m <- m[wMtch]
  dTsts <- dTsts[m,]
  fdrkl <- paste0("mod. F-test Significant-FDR=", fdr*100, "%")
  for (i in colnames(dTsts)) { #i <- colnames(dTsts)[1]
    ctrst <- expContr_F$name[match(i, expContr_F$Contrasts)]
    rgkl <- paste0(regRoot_F, ctrst)
    wUp <- which(dTsts[[i]] == 1)
    wDown <- which(dTsts[[i]] == -1)
    my_F_Data[c(wUp, wDown), fdrkl] <- "+"
    my_F_Data[wMtch[wUp], rgkl] <- paste0("up, FDR = ", fdr*100, "%")
    my_F_Data[wMtch[wDown], rgkl] <- paste0("down, FDR = ", fdr*100, "%")
  }
}
# Second level: apply lfc filtering
if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
  lfcThr <- Param$Ratios.Contamination.Rates
}
if (Param$Ratios.Thresholds == threshMsg) {
  d <- abs(unlist(refRat_F))
  d <- sort(d, decreasing = TRUE)
  if (Param$Ratios.Contamination.Rates > 0) {
    d <- d[floor(length(d)*as.numeric(Param$Ratios.Contamination.Rates))]
  } else {
    d <- max(d)+0.000001
  }
  lfcThr <- d
}
if (lfcThr != 0) {
  warning("Please note that limma allows applying a non-null lfc threshold, but does not recommend to do it (see ?decideTests for a discussion).")
}
for (i in 1:nrow(expContr_F)) { #i <- 1
  fckl <- paste0("Mean ", ratRef, expContr_F$name[i])
  rgkl <- paste0(regRoot_F, expContr_F$name[i])
  w <- which((my_F_Data[[rgkl]] != "non significant")&(abs(my_F_Data[[fckl]]) < lfcThr))
  my_F_Data[w, rgkl] <- "too small FC"
}
myData[, c(F_Root, fdrKol, regKol)] <- my_F_Data[, c(F_Root, fdrKol, regKol)]
#View(my_F_Data)

# Mulcom - essentially a "moderated" Dunnett's test
# Was a good idea... but it is reproducibly causing R to crash...
# I will give it up for now...
if ((!exists("Mulcom"))||(length(Mulcom) != 1)||(!is.logical(Mulcom))||(is.na(Mulcom))) { Mulcom <- FALSE }
Mulcom %<o% Mulcom
if (Mulcom) {
  if(!require(Biobase)) { pak::pak("Biobase") }
  library(Biobase)
  if(!require(Mulcom)) { pak::pak("Mulcom") }
  library(Mulcom)
  # See https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/Mulcom/inst/doc/MulcomVignette.pdf
  m <- match(colnames(tmpVal), rownames(designMatr))
  grps <- apply(designMatr[m,], 1, function(x) { colnames(designMatr)[which(x == 1)] })
  grps <- gsub("^Group_", "", grps)
  grps2 <- match(grps, unique(grps))-1
  tmpVal2 <- proteoCraft::Data_Impute2(tmpVal, grps)
  tmpVal2 <- as.matrix(tmpVal2$Imputed_data)
  # rownames(tmpVal2) <- myData[[idCol]]
  # colnames(tmpVal2) <- gsub("___", "_", colnames(tmpVal2))
  # phDat <- data.frame(Sample = colnames(tmpVal2),
  #                     Group = grps2)
  # rownames(phDat) <- colnames(tmpVal2)
  # phDat <- Biobase::AnnotatedDataFrame(data = phDat)
  # ftDat <- myData[, unique(c(idCol, protCol)), drop = FALSE]
  # rownames(ftDat) <- rownames(tmpVal2)
  # ftDat <- Biobase::AnnotatedDataFrame(data = ftDat)
  # tmpVal2 <- Biobase::ExpressionSet(assayData = tmpVal2,
  #                                   phenoData = phDat,
  #                                   featureData = ftDat)
  # mulcom_scores <- mulScores(tmpVal2, grps2)
  #data(benchVign)
  #eset <- exprs(Affy); index = Affy$Groups
  #mulcom_scores <- mulScores(Affy, Affy$Groups)
  #
  # A slightly rewritten function because the original crashes the R session during the C call... though the reason why is not obvious...
  # It is also rewritten more concisely for clarity.
  # Unfortunately... it still crashes!
  mulScores2 <- function(eset, index) {
    #eset <- tmpVal2; index = grps2
    data(benchVign)
    eset <- exprs(Affy); index = Affy$Groups
    
    if ((!is.vector(index))||(sum(c("ExpressionSet", "matrix", "data.frame") %in% class(eset)) == 0)) {
      stop("error in input files", call. = FALSE)
    }
    mulcom <- new("MULCOM")
    if ("ExpressionSet" %in% class(eset)) {
      data <- as.vector(as.matrix(exprs(eset)))
      rwNms <- featureNames(eset)
    } else {
      if ("data.frame" %in% class(eset)) { eset <- as.matrix(eset) }
      data <- as.vector(eset)
      rwNms <- rownames(eset)
    } 
    Grps <- index
    nGrps <- length(levels(factor(Grps)))
    nGrps_1 <- nGrps - 1
    nRws <- nrow(eset)
    nCol <- ncol(eset)
    reference <- c(0)
    means <- mse <- c(seq(0, 0, length = nRws*nGrps_1))
    SS <- c(seq(0, 0, length = nRws*nGrps))
    harmonic_means <- c(seq(0, 0, length = nGrps_1))
    sss2 <- c(seq(0, 0, length = nRws))
    out <- .Call("Single_SimulationC", as.double(data), 
                 as.double(means), as.double(harmonic_means), 
                 as.double(SS), as.double(sss2), as.double(mse), 
                 as.integer(nRws), as.integer(nCol), as.integer(Grps), 
                 as.integer(nGrps), as.integer(reference),
                 PACKAGE = "Mulcom")
    mulcom@FC <- t(data.frame(matrix(out[[2]],
                                     ncol = nGrps_1,
                                     byrow = TRUE),
                              row.names = rwNms))
    mulcom@HM <- matrix(out[[3]],
                        ncol = nGrps_1,
                        byrow = TRUE)
    mse <- data.frame(matrix(out[[6]],
                             ncol = nGrps_1,
                             byrow = TRUE),
                      row.names = rwNms)
    mulcom@MSE_Corrected <- t(mse)
    return(mulcom)
  }
  mulcomScore <- mulScores2(tmpVal2, grps2)
}

#
# Plot F-test results
# Q-Q plot
ttl <- "mod. F-test QQ plot"
fl <- paste0(ohDeer, "/", ttl)
jpeg(file = paste0(fl, ".jpeg"), width = 400, height = 350)
qqt(as.data.frame(fit$F), df = fit$df.prior + fit$df.residual, pch = 16, cex = 0.2)
abline(0,1)
dev.off()
#
kol <- unique(c(idCol, protCol, plotlyLab, "Potential contaminant", Param$Plot.labels, xCol, Alpha, fdrKol, regKol))
kol <- kol[which(kol %in% colnames(myData))]
my_F_Data[, kol] <- myData[, kol]
my_F_Data[["Rel. av. log10 abundance"]] <- myData[[Size]]
tmp <- lapply(1:length(xCol), function(x) { myData[, F_Root, drop = FALSE] })
tmp <- do.call(cbind, tmp)
dummyPVCol <- gsub(topattern(paste0("Mean ", ratRef)), paste0(F_Root, " - "), xCol)
my_F_Data[, dummyPVCol] <- tmp                
contr <- data.frame(Contrast = expContr_F$name)
if ("Target" %in% colnames(Exp.map)) {
  tmp <- gsub("^\\(|\\) - \\(.*", "", contr$Contrast)
  contr$Target <- Exp.map$Target[match(tmp, Exp.map[[VPAL$column]])]
}
aggr_dummy <- data.frame(Aggregate.Name = "Contrast", Characteristics = "Contrast")
aggr_list_dummy <- list(Contrast = unlist(contr))
#
F_volc <- Volcano.plot(Prot = my_F_Data, mode = "custom", experiments.map = contr,
                       X.root = paste0("Mean ", ratRef),
                       Y.root = paste0(F_Root, " - "),
                       aggregate.map = aggr_dummy, aggregate.list = aggr_list_dummy,
                       aggregate.name = "Contrast",
                       parameters = Param,
                       save = c("jpeg", "pdf"),
                       FDR.root = "mod. F-test Significant-FDR=",
                       Ref.Ratio.values = refRat_F,
                       ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                       arbitrary.lines = arbitrary.thr,
                       proteins = prot.list, proteins_split = protsplit, IDs.col = idCol,
                       Proteins.col = protCol,
                       return = FALSE, return.plot = TRUE, title = "F-test volcano plot ",
                       subfolder = ohDeer, subfolderpertype = FALSE,
                       Symmetrical = TRUE,
                       Alpha = Alpha, Size = "Rel. av. log10 abundance", Size.max = 2,
                       plotly = create_plotly, plotly_local = create_plotly_local,
                       plotly_labels = plotlyLab,
                       Ref.Ratio.method = paste0("obs", RefRat_Mode),
                       cl = parClust,
                       reg.root = regRoot_F)
#
# Save plotly plots
dr <- ohDeer
myPlotLys <- F_volc
Src <- paste0(libPath, "/extdata/R scripts/Sources/save_Volcano_plotlys.R")
#rstudioapi::documentOpen(Src)
source(Src, local = FALSE)
#
tmp_Fdat <- F_volc$Protein_groups_file
if (!is.null(tmp_Fdat)) { # Legacy code: since 6.3.1.7, we do not get the Regulated columns from the Volcano.plot function anymore!
  stop("Uh... why was the condition TRUE?!?!")
  tmp_Fdat[[F_Root]] <- myData[[F_Root]]
  tmp_Fdat[[mtchCol]] <- myData[[mtchCol]]
  kol <- colnames(tmp_Fdat)
  g <- grep("^Regulated - ", kol)
  colnames(tmp_Fdat)[g] <- paste0("mod. F-test ", colnames(tmp_Fdat)[g])
  kol <- colnames(tmp_Fdat)
  kol <- kol[which(!kol %in% dummyPVCol)]
  tmp_Fdat <- tmp_Fdat[, kol]
  w <- which(!kol %in% colnames(myData))
  #kol[w]
  myData[, kol[w]] <- tmp_Fdat[, kol[w]]
}
#
F_thresh <- F_volc$Thresholds
thresh <- lapply(names(F_thresh$Absolute), function(x) {
  y <- F_thresh$Absolute[[x]]
  x <- data.frame(Group = rep(cleanNms(x), nrow(y)))
  return(cbind(x, y))
})
thresh <- plyr::rbind.fill(thresh)
thresh$Name <- NULL
thresh$Root <- gsub(" - $", "", thresh$Root)
thresh$Value <- thresh$Text.value
thresh$Text.value <- NULL
# Legacy code: since aRmel 6.3.1.7, we do not get the Regulated columns from the Volcano.plot function anymore!
# fdrThresh <- F_thresh$FDR
# if (!is.null(fdrThresh)) {
#   #stop("Uh... why was the condition TRUE?!?!")
#   colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.up")] <- "Colour (up)"
#   colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.down")] <- "Colour (down)"
#   colnames(fdrThresh)[which(colnames(fdrThresh) == "fdr.col.line")] <- "Colour (line)"
#   if ("Test" %in% colnames(fdrThresh)) {
#     fdrThresh <- fdrThresh[, c("Test", colnames(fdrThresh)[which(colnames(fdrThresh) != "Test")])]
#   }
# }
fl <- paste0(ohDeer, "/Thresholds.xlsx")
wb <- wb_workbook()
wb <- wb_set_creators(wb, "Me")
wb <- wb_add_worksheet(wb, "Thresholds")
dms <- wb_dims(2, 1)
wb <- wb_add_data(wb, "Thresholds", "Absolute thresholds", dms)
wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                  italic = "true", underline = "single")
dms <- wb_dims(3, 2)
wb <- wb_add_data_table(wb, "Thresholds", thresh, dms,
                        col_names = TRUE, table_style = "TableStyleMedium2",
                        banded_rows = TRUE, banded_cols = FALSE)
dms <- wb_dims(3+nrow(thresh)+3, 1)
wb <- wb_add_data(wb, "Thresholds", "FDR thresholds", dms)
wb <- wb_add_font(wb, "Thresholds", dms, "Calibri", wb_color(hex = "FF000000"), bold = "true",
                  italic = "true", underline = "single")
dms <- wb_dims(3+nrow(thresh)+4, 2)
# if (!is.null(fdrThresh)) {
#   wb <- wb_add_data_table(wb, "Thresholds", fdrThresh, dms,
#                           col_names = TRUE, table_style = "TableStyleMedium2",
#                           banded_rows = TRUE, banded_cols = FALSE)
#   wb <- wb_set_col_widths(wb, "Thresholds", 1, 3)
#   tmp1 <- rbind(colnames(thresh), thresh)
#   colnames(tmp1) <- paste0("V", 1:ncol(tmp1))
#   tmp2 <- rbind(colnames(fdrThresh), fdrThresh)
#   colnames(tmp2) <- paste0("V", 1:ncol(tmp2))
#   tst <- plyr::rbind.fill(tmp1, tmp2)
#   tst <- setNames(apply(tst, 2, function(x) { max(nchar(x), na.rm = TRUE) }), NULL)
#   wb <- wb_set_col_widths(wb, "Thresholds", 1:(length(tst)+1), c(3, tst))
# }
wb_save(wb, fl)
#xl_open(fl)
#
if (dataType == "modPeptides") {
  ptmpep <- myData
  #PTMs_F_test_data[[Ptm]] <- tmp_Fdat
  PTMs_F_test_data[[Ptm]] <- my_F_Data
}
if (dataType == "PG") {
  PG <- myData
  #F_test_data %<o% tmp_Fdat
  F_test_data %<o% my_F_Data
}
#
# Cleanup
for (i in allArgs) { try(rm(i), silent = TRUE) }
