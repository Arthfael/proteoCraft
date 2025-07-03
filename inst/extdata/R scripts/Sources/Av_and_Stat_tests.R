if (!dir.exists(samDir)) { dir.create(samDir, recursive = TRUE) }
if (!dir.exists(ebamDir)) { dir.create(ebamDir, recursive = TRUE) }
source(parSrc, local = FALSE)
clusterExport(parClust, list("AltHyp", "Nested", "samDir", "ebamDir"), envir = environment())
clusterCall(parClust, function() library(siggenes))
if (dataType == "modPeptides") {
  #reNorm <- FALSE
  myData <- ptmpep
  intRef <- pepRf
  ratRef <- pepRatRf
  namesCol <- "Name"
  namesRoot <- "Pep"
  runLoc <- LocAnalysis
  ohDeer <- paste0(wd, "/Reg. analysis/", ptm, "/t-tests")
}
if (dataType == "PG") {
  #reNorm <- Norma.Prot.Ratio.classic
  myData <- PG
  intRef <- Prot.Expr.Root
  ratRef <- Prot.Rat.Root
  namesCol <- "Leading protein IDs"
  namesRoot <- "PG"
  runLoc <- FALSE
  ohDeer <- paste0(wd, "/Reg. analysis/t-tests")
}
expContr <- expContrasts
if (!dir.exists(ohDeer)) { dir.create(ohDeer, recursive = TRUE) }
#
# First calculate mean per group
cat(" - Mean sample group intensities\n")
for (smplGrp in VPAL$values) { #smplGrp <- VPAL$values[1] #smplGrp <- VPAL$values[2]
  w <- which(Exp.map[[VPAL$column]] == smplGrp)
  e1 <- Exp.map[w,]
  stopifnot(length(unique(e1$Ref.Sample.Aggregate)) == length(e1$Ref.Sample.Aggregate))
  smpls1 <- e1$Ref.Sample.Aggregate
  col1 <- paste0(intRef, smpls1)
  #col2 <- paste0(ratRef, smpls1)
  w1 <- which(col1 %in% colnames(myData))
  #w2 <- which(col2 %in% colnames(myData))
  if (length(w1)) { # Calculate average expression
    e1 <- e1[w1,]; col1 <- col1[w1]
    e1 <- e1[order(e1$Replicate),]
    ke1 <- paste0("Mean ", intRef, smplGrp)
    ke1s <- paste0(gsub(" - $", ": SE - ", intRef), smplGrp)
    tempVal <- myData[, col1, drop = FALSE]
    row.names(tempVal) <- myData[[namesCol]]
    if (length(w1) == 1) {
      cat(paste0("Warning: Only 1 replicate for ", i1a, "!\n"))
      myData[[ke1]] <- tempVal[[1]]
      if (runLoc) { myData[[kl1]] <- tempLoc[[1]] }
    } else {
      myData[, c(ke1, ke1s)] <- as.data.frame(t(parApply(parClust, tempVal, 1, Av_SE_fun)))
      if (runLoc) {
        myData[, c(kl1, kl1s)] <- as.data.frame(t(parApply(parClust, tempLoc, 1, Av_SE_fun)))
      }
    }
    if (runLoc) {
      tempLoc <- myData[, paste0(Prot.Expr.Root2, smpls1), drop = FALSE]
      row.names(tempLoc) <- myData[[namesCol]]
      kl1 <- paste0("Mean ", Prot.Expr.Root2, smplGrp)
      kl1s <- paste0(gsub(" - $", ": SE - ", Prot.Expr.Root2), smplGrp)
    }
  }
}

expContr$Map <- lapply(expContr$All, function(x) { #x <- expContr$All[[1]]
  em <- expMap[which(expMap$Group_ %in% unlist(x)),]
  em$Expression_Column <- paste0(intRef, em$Ref.Sample.Aggregate)
  em$Ratios_Column <- paste0(ratRef, em$Ref.Sample.Aggregate)
  em[which(em$Expression_Column %in% colnames(myData)),]
  em <- em[c(which(em$Reference),
             which(!em$Reference)),]
  return(em)
})

# Average LFCs
# (Calculate average ratios)
tmp <- lapply(1:nrow(expContr), function(x) { #x <- 1 #x <- 2
  #print(x)
  em <- expContr$Map[x][[1]]
  nm <- expContr$name[x]
  w0 <- which(em$Reference)
  w1 <- which(!em$Reference)
  w1r <- which((!em$Reference)&(em$Ratios_Column %in% colnames(myData)))
  kr1 <- paste0("Mean ", ratRef, nm) # Eventually we can rename those, either to contrast name or to full contrast (X-Y)
  kr1s <- paste0(gsub(" - $", ": SE - ", ratRef), nm) # Eventually we can rename those, either to contrast name or to full contrast (X-Y)
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
  #  - Some software (MaxQuant) can in some cases (SILAC) measure ratios more precisely than intensity.
  #    Thus sometimes software-provided ratios are quantitatively more informative than intensities!
  #    In these cases, this intensity-based workflow should as early as feasible propagate ratios-level information onto quantities:
  #          - if the software provides Ia, Ib and R(a/b), then the new intensities (J) will be defined by Ja+Jb = Ia+Ib and Ja/Jb = R(a/b)
  #          - after that, ratios can be ignored again.
  #    For historical reasons (it evolved from a SILAC workflow), this workflow also calculates ratios at different stages.
  #    This should be removed! Only provide average LFC at the end!
  #    Note though that we will still have to measure ref-to-ref or intra-group ratios... unless we also rewrite all of the statistics.
  #  - We also plot a log2 average ratio/fold change (LFCs) in volcano plots. This is the only thing we need.
  #
  # As part of rewriting this to focus on intensity data, we will calculate ratios here directly from expression.
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
      RES[, c(kr1, kr1s)] <- as.data.frame(t(parApply(parClust, rats, 1, Av_SE_fun)))
    }
  } else {
    k1 <- em$Expression_Column[w1]
    k0m <- paste0("Mean ", intRef, unique(em[w0, VPAL$column]))
    rats <- sweep(myData[, k1], 1, myData[[k0m]], "-")/log10(2)
    if (length(w1) == 1) {
      RES[[kr1]] <- rats
    } else {
      RES[, c(kr1, kr1s)] <- as.data.frame(t(parApply(parClust, rats, 1, Av_SE_fun)))
    }
  }
  return(RES)
})
tmp <- do.call(cbind, tmp)
myData[, colnames(tmp)] <- tmp

# Statistical tests:
#
kol <- paste0(intRef, rownames(designMatr))
kol <- kol[which(kol %in% colnames(myData))]
tmpVal <- myData[, kol]
clusterExport(parClust, list("tmpVal", "Nested", "AltHyp"), envir = environment())
#
# - Student's and Welch's t-tests:
for (ii in 1:2) {
  tmp <- lapply(1:nrow(expContr), function(x) { #x <- 1
    em <- expContr$Map[x][[1]]
    k0 <- em$Expression_Column[which(em$Reference)]
    k1 <- em$Expression_Column[which(!em$Reference)]
    clusterExport(parClust, list("k1", "k0", "ii"), envir = environment())
    RES <- parSapply(parClust, 1:nrow(tmpVal), function(y) {
      v0 <- as.numeric(tmpVal[y, k0])
      v1 <- as.numeric(tmpVal[y, k1])
      if (Nested) {
        w <- which((proteoCraft::is.all.good(v0, 2))&(proteoCraft::is.all.good(v1, 2)))
        v0 <- v0[w]; v1 <- v1[w]
      } else {
        v0 <- proteoCraft::is.all.good(v0)
        v1 <- proteoCraft::is.all.good(v1)
      }
      if ((length(unique(v0)) > 1)&&(length(unique(v1)) > 1)) {
        res <- -log10(t.test(x = v1, y = v0, paired = Nested, alternative = AltHyp,
                             var.equal = (ii == 1))$p.value)
      } else { res <- NA }
      return(res)
    })
    return(RES)
  })
  tmp <- as.data.frame(do.call(cbind, tmp))
  colnames(tmp) <- paste0(c(StudentRoot, WelchRoot)[ii], expContr$name)
  myData[, colnames(tmp)] <- tmp
}
#
# - Permutations t-test:
#   (We must obviously ignore nesting here!!!)
tmp <- lapply(1:nrow(expContr), function(x) { #x <- 1
  RES <- rep(NA, nrow(tmpVal))
  em <- expContr$Map[x][[1]]
  em <- em[order(em$Reference, em$Replicate),] # Redundant, but good safety check
  k0 <- em$Expression_Column[which(em$Reference)]
  k1 <- em$Expression_Column[which(!em$Reference)]
  test0a <- parApply(parClust, tmpVal[, k0], 1, function(x) {
    length(unique(proteoCraft::is.all.good(x)))
  }) > 1
  test1a <- parApply(parClust, tmpVal[, k1], 1, function(x) {
    length(unique(proteoCraft::is.all.good(x)))
  }) > 1
  w <- which((test0a)&(test1a))
  if (length(w)) {
    map <- data.frame(Sample = em[[RSA$column]],
                      Role_ = em[[VPAL$column]])
    clusterExport(parClust, list("k1", "k0", "map", "AltHyp"), envir = environment())
    RES[w] <- parApply(parClust, tmpVal[w, c(k0, k1)], 1, function(x) {
      #x <- tmpVal[w[1], c(k0, k1)]
      x <- data.frame(value = as.numeric(x),
                      Sample = as.factor(map$Sample),
                      Group = as.factor(map$Role_))
      x <- x[which(proteoCraft::is.all.good(x$value, 2)),]
      res <- -log10(coin::pvalue(coin::independence_test(value ~ Group,
                                                         data = x, alternative = AltHyp)))
      return(res)
    })
  }
  return(RES)
})
tmp <- set_colnames(do.call(cbind, tmp), paste0(permRoot, expContr$name))
myData[, colnames(tmp)] <- tmp

# - Moderated t-test (limma) and modified DEqMS version:
#     We will run two versions:
#      - a classic limma moderated t-test
#      - the DEqMS-corrected version
#
# Do not use duplicateCorrelation(), it is for duplicate rows in probe arrays!!!
TESTs <- c("limma", "DEqMS")
for (TEST in TESTs) { #TEST <- TESTs[1] #TEST <- TESTs[2]
  tmpVal2 <- tmpVal
  if (TEST == "limma") {
    pRoot <- modRoot
    insrt <- "limma mod. t-test"
  }
  if (TEST == "DEqMS") {
    grpTst <- aggregate(1:nrow(expMap), list(expMap$Group_), list)
    clusterExport(parClust, list("grpTst", "tmpVal2"), envir = environment())
    tst <- t(parSapply(parClust, 1:nrow(tmpVal2), function(x) {
      sapply(grpTst$x, function(y) { length(proteoCraft::is.all.good(tmpVal2[x, unlist(y)])) })
    }))
    tst <- apply(tst, 1, min)
    wOK <- which(tst >= 2)
    tmpVal2 <- tmpVal2[wOK,]
    nmRoot <- deqmsRoot
    insrt <- "DEqMS mod. t-test"
  }
  #fit <- lmFit(tmpVal2, designMatr, trend = TRUE, robust = TRUE)
  fit <- lmFit(tmpVal2, designMatr)
  fit$genes <- myData[[namesCol]]
  fit <- contrasts.fit(fit, contrMatr)
  fit <- eBayes(fit) # Note: the default proportion of 1% in eBayes is the Bayesian prior, eBayes is robust to deviations from this value.
  if (AltHyp != "two.sided") {
    fit$p.value <- limma.one.sided(fit, c(FALSE, TRUE)[match(AltHyp, c("greater", "lower"))])
  }
  # Do not use topTable!
  # - It adjusts P values (with "BH") by default...
  #   If using it (e.g. rewriting to include systematic P-value adjustment), set adjust.method = "none"!!!
  # - The way it tests a hypothesis seems better suited to F.tests than t-tests!
  # Here you could use decideTests() if wanting to switch to directly using decisions from limma
  if (TEST == "limma") {
    myData[, paste0(modRoot, expContr$name)] <- -log10(fit$p.value)
  }
  # Plot moderated t-test results
  for (contr in colnames(contrMatr)) { #contr <- colnames(contrMatr)[1]
    nm <- expContr$name[match(contr, colnames(contrMatr))]
    nm <- cleanNms(nm)
    # Q-Q plot
    ttl <- paste0("mod. t-test QQ plot - ", nm)
    fl <- paste0(ohDeer, "/", ttl)
    jpeg(file = paste0(fl, ".jpeg"), width = 400, height = 350)
    qqt(as.data.frame(fit$t)[[contr]], df = fit$df.prior + fit$df.residual, pch = 16, cex = 0.2)
    abline(0,1)
    dev.off()
    # MA plot
    ttl <- paste0(insrt, " MA plot - ", nm)
    fl <- paste0(ohDeer, "/", ttl)
    jpeg(file = paste0(fl, ".jpeg"), width = 400, height = 350)
    fit$genes <- NULL
    plotMD(fit, column = match(contr, colnames(fit$p.value)))
    dev.off()
  }
  if (TEST == "DEqMS") {
    # DEqMS:
    # This is meant as en extension to limma to account "for variance dependence on the number of quantified peptides or PSMs
    # for statistical testing of differential protein expression."
    # https://www.bioconductor.org/packages/devel/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html#overview-of-deqms
    bioc_req <- unique(c(bioc_req, "DEqMS"))
    biocInstall("DEqMS")
    cran_req <- unique(c(cran_req, "matrixStats"))
    if (!require(matrixStats, quietly = TRUE)) { install.packages("matrixStats") }
    library(matrixStats)
    if (dataType == "PG") {
      countCol <- grep("^Evidences count - ", colnames(myData), value = TRUE)
    }
    if (dataType == "modPeptides") {
      idCol <- grep("^Evidence IDs - ", colnames(myData), value = TRUE)
      countCol <- gsub("^Evidence IDs - ", "Evidences count - ", idCol)
      myData[, countCol] <- sapply(idCol, function(k) { sapply(strsplit(myData[[k]], ";"), length) })
    }
    psm.counts <- data.frame(count = rowMins(as.matrix(myData[wOK, countCol])),
                             row.names = row.names(myData)[wOK])
    psm.counts <- psm.counts+1 # Pseudo-counts, to allow for PGs with 0 PSMs
    # # We need a least 1 and will treat 0 and 1 together
    # w <- which(psm.counts$count == 0)
    # psm.counts$count[w] <- 1
    fit$count <- psm.counts[rownames(fit$coefficients), "count"]
    fit2 <- spectraCounteBayes(fit)
    # Variance plot
    tst <- try({
      ttl <- paste0(insrt, " variance plot")
      fl <- paste0(ohDeer, "/", ttl)
      jpeg(file = paste0(fl, ".jpeg"), width = 400, height = 350)
      VarianceBoxplot(fit2, n = max(psm.counts), main = ttl, xlab = "PSM count")
      dev.off()
    }, silent = TRUE)
    #
    for (i in 1:ncol(contrMatr)) {
      contrNm <- expContrasts$name[match(colnames(contrMatr)[i], expContrasts$Contrasts)]
      DEqMS.results <- outputResult(fit2, coef_col = i)
      kol <- paste0(deqmsRoot, contrNm)
      myData[[kol]] <- NA
      myData[wOK, kol] <- -log10(DEqMS.results$sca.P.Value)
    }
  }
}

# PolySTest
usePolySTest <- FALSE
if (usePolySTest) {
  packs <- c("SummarizedExperiment", "PolySTest")
  bioc_req <- unique(c(bioc_req, packs))
  for (pack in packs) {
    if (!require(pack, character.only = TRUE)) { biocInstall(pack) }
  }
  ohDeer2 <- paste0(wd, "/Reg. analysis/PolySTests")
  #
  kol <- paste0(intRef, rownames(designMatr))
  kol <- kol[which(kol %in% colnames(myData))]
  tmpVal3 <- as.matrix(myData[, kol])
  w <- which(is.infinite(tmpVal3)&is.na(tmpVal3), arr.ind = TRUE)
  tmpVal3[w] <- NA
  tmpVal3 <- as.data.frame(tmpVal3)
  row.names(tmpVal3) <- 1:nrow(myData)
  kolst <- setNames(lapply(VPAL$values, function(x) { #x <- VPAL$values[1]
    x <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)]
    x <- paste0(intRef, x)
    x[which(x %in% colnames(tmpVal3))]
  }), VPAL$values)
  tst <- vapply(1:nrow(tmpVal3), function(x) {
    min(vapply(VPAL$values, function(y) { length(is.all.good(unlist(tmpVal3[x, kolst[[y]]]))) }, 1))
  }, 1)
  w <- which(tst > 2)
  tmpVal3 <- tmpVal3[w,]
  smplsDat <- data.frame(Sample = gsub(topattern(intRef), "", kol))
  smplsDat[, c("Condition", "Replicate", "Reference")] <- Exp.map[match(smplsDat$Sample, Exp.map$Ref.Sample.Aggregate), c(VPAL$column, "Replicate", "Reference")]
  cond <- c(unique(smplsDat$Condition[which(smplsDat$Reference)]),
            unique(smplsDat$Condition[which(!smplsDat$Reference)]))
  dat <- SummarizedExperiment(assays = list(quant = tmpVal3),
                              colData = smplsDat,
                              metadata = list(NumCond = length(cond),
                                              NumReps = length(unique(smplsDat$Replicate))))
  allComps <- create_pairwise_comparisons(cond, 1)
  # fixed functions
  rp_unpairedMod <- function (tData, trefData, cl = parClust) {
    if (exists("NumRPPairs")) { NumRPPairs <- NumRPPairs } else { NumRPPairs <- 100 }
    if (ncol(tData) != ncol(trefData)) {
      stop("The number of columns in the two datasets must be the same")
    }
    NumReps <- ncol(tData)
    tpRPvalues <- matrix(NA, ncol = NumRPPairs, nrow = nrow(tData), dimnames = list(rows = rownames(tData), cols = seq_len(NumRPPairs)))
    parallel::clusterExport(cl, list("NumReps", "tData", "trefData", "RPStats", "rowMins"), environment())
    RPparOut <- parallel::parLapply(cl, seq_len(NumRPPairs), function(x) {
      tRPMAData <- tData[, sample(seq_len(NumReps)), drop = FALSE] - 
        trefData[, sample(seq_len(NumReps)), drop = FALSE]
      RPMAUp_pvalues <- RPStats(tRPMAData, NumReps)
      RPMADown_pvalues <- RPStats(-tRPMAData, NumReps)
      ttt <- matrixStats::rowMins(cbind(RPMAUp_pvalues, 
                                        RPMADown_pvalues), na.rm = TRUE) * 2
      ttt[ttt > 1] <- 1
      names(ttt) <- names(RPMAUp_pvalues)
      return(ttt)
    })
    save(RPparOut, tData, trefData, file = paste0(ohDeer2, "/tmp.csv"))
    for (p in seq_len(NumRPPairs)) {
      names(RPparOut[[p]]) <- rownames(tData)
      tpRPvalues[names(RPparOut[[p]]), p] <- RPparOut[[p]]
    }
    tpRPvalues[!is.finite(tpRPvalues)] <- NA
    pRPvalues <- rowMeans(tpRPvalues, na.rm = TRUE)
    qRPvalues <- rep(NA, length(pRPvalues))
    names(qRPvalues) <- names(pRPvalues)
    tqs <- p.adjust(na.omit(pRPvalues), method = "BH")
    qRPvalues[names(tqs)] <- tqs
    return(list(pRPvalues = pRPvalues, qRPvalues = qRPvalues))
  }
  PolySTest_unpaired <- function (fulldata, allComps, statTests = c("limma", "Miss_Test", "t_test", "rank_products", "permutation_test")) {
    check_for_polystest(fulldata)
    Data <- assay(fulldata, "quant")
    NumReps <- metadata(fulldata)$NumReps
    NumCond <- metadata(fulldata)$NumCond
    Reps <- rep(seq_len(NumCond), NumReps)
    conditions <- unique(colData(fulldata)$Condition)
    NumComps <- nrow(allComps)
    RRCateg <- matrix(NA, ncol = nrow(allComps), nrow = 2)
    RRCateg[1, ] <- match(allComps[, 1], conditions)
    RRCateg[2, ] <- match(allComps[, 2], conditions)
    Data <- Data - rowMeans(Data, na.rm = TRUE)
    tests <- statTests
    if (any(!tests %in% c("limma", "Miss_Test", "t_test", "rank_products", 
                          "permutation_test"))) {
      stop("Invalid test name(s) specified. They should one or more of:\n         'limma', 'Miss_Test', 't_test', 'rank_products', 'permutation_test'")
    }
    p_values <- q_values <- matrix(NA, nrow = nrow(Data), ncol = length(tests) * 
                                     NumComps)
    rownames(p_values) <- rownames(q_values) <- rownames(Data)
    colnames(p_values) <- paste0("p_values_", rep(tests, each = NumComps), 
                                 "_", rep(seq_len(NumComps), length(tests)))
    colnames(q_values) <- paste0("q_values_", rep(tests, each = NumComps), 
                                 "_", rep(seq_len(NumComps), length(tests)))
    test_funcs <- list(
      limma = function(Data, RRCateg) {
        lm_out <- limma_unpaired(Data, NumCond, NumReps, RRCateg)
        list(pvals = lm_out$plvalues, qvals = lm_out$qlvalues, 
             Sds = lm_out$Sds)
      },
      Miss_Test = function(Data, RRCateg) {
        MissingStats <- MissingStatsDesign(Data, RRCateg, NumCond, 
                                           NumReps)
        list(pvals = MissingStats$pNAvalues, qvals = MissingStats$qNAvalues)
      },
      t_test = function(tData, trefData) {
        ttest_out <- ttest_unpaired(tData, trefData)
        list(pvals = ttest_out$ptvalues, qvals = ttest_out$qtvalues)
      },
      rank_products = function(tData, trefData) {
        rp_out <- rp_unpairedMod(tData, trefData)
        list(pvals = rp_out$pRPvalues, qvals = rp_out$qRPvalues)
      },
      permutation_test = function(tData, trefData) {
        perm_out <- perm_unpaired(tData, trefData)
        list(pvals = perm_out$pPermutvalues, qvals = perm_out$qPermutvalues)
      })
    Sds <- NULL
    for (test in statTests[statTests %in% names(test_funcs)]) {
      if (test %in% c("limma", "Miss_Test")) {
        message("Running ", test, " test")
        res <- test_funcs[[test]](Data, RRCateg)
        p_values[, grep(paste0("p_values_", test), colnames(p_values))] <- res$pvals
        q_values[, grep(paste0("q_values_", test), colnames(q_values))] <- res$qvals
        if (test == "limma") 
          Sds <- res$Sds
        message(test, " completed")
      }
    }
    lratios <- NULL
    pb <- txtProgressBar(0.9, NumComps)
    for (vs in seq_len(NumComps)) {
      if (!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::setProgress(0.1 + 0.3/(NumComps) * vs, detail = paste("tests for comparison", 
                                                                     vs, "of", NumComps))
      }
      tData <- Data[, Reps == RRCateg[1, vs]]
      trefData <- Data[, Reps == RRCateg[2, vs]]
      for (test in statTests[statTests %in% c("t_test", "rank_products", 
                                              "permutation_test")]) {
        res <- test_funcs[[test]](tData, trefData)
        p_values[, grep(paste0("p_values_", test), colnames(p_values))[vs]] <- res$pvals
        q_values[, grep(paste0("q_values_", test), colnames(q_values))[vs]] <- res$qvals
      }
      lratios <- cbind(lratios,
                       rowMeans(Data[, Reps == RRCateg[1, vs], drop = FALSE],
                                na.rm = TRUE) - rowMeans(Data[, Reps == RRCateg[2, vs],
                                                              drop = FALSE], na.rm = TRUE))
      setTxtProgressBar(pb, vs)
    }
    close(pb)
    message("tests completed")
    fulldata <- PolySTest:::prepare_output_data(fulldata, p_values, q_values, 
                                                lratios, tests, allComps)
    return(fulldata)
  }
  permtest_pairedMod <- function (tMAData, cl = parClust) {
    NumReps <- ncol(tMAData)
    if (exists("NumPermCols")) { NumPermCols <- NumPermCols } else { NumPermCols <- 7 }
    if (exists("NumTests")) { NTests <- NumTests } else { NTests <- 1000 }
    if (ncol(tMAData) < NumPermCols) {
      AddDat <- matrix(sample(as.vector(tMAData),
                              (NumPermCols - ncol(tMAData)) * nrow(tMAData), replace = TRUE),
                       nrow = nrow(tMAData))
      PermMAData <- cbind(tMAData, AddDat)
    } else {
      PermMAData <- tMAData
    }
    #
    tMAData <- as.matrix(tMAData)
    RealStats <- PolySTest::StatsForPermutTest(tMAData, Paired = TRUE)
    parallel::clusterExport(cl, list("NumReps", "PermMAData", "RPStats"), environment())
    PermutOut <- parallel::parLapply(cl, seq_len(NTests), function(x) {
      indat <- apply(PermMAData, 1, function(y) {
        y <- unlist(y)
        sample(y, NumReps) * sample(c(1, -1), NumReps, replace = TRUE)
      })
      PolySTest::StatsForPermutTest(t(indat), TRUE)
    })
    PermutOut <- matrix(unlist(PermutOut), nrow = nrow(tMAData))
    pPermutvalues <- apply(cbind(RealStats, PermutOut), 1, function(x) {
      ifelse(is.na(x[1]), NA, (1 + sum(x[1] < x[-1], na.rm = TRUE))/(sum(!is.na(x))))
    })
    qPermutvalues <- rep(NA, length(pPermutvalues))
    names(qPermutvalues) <- names(pPermutvalues)
    tqs <- p.adjust(na.omit(pPermutvalues), method = "BH")
    qPermutvalues[names(tqs)] <- tqs
    return(list(pPermutvalues = pPermutvalues, qPermutvalues = qPermutvalues))
  }
  PolySTest_paired <- function(fulldata,
                               allComps,
                               statTests = c("limma", "Miss_Test", "t_test", "rank_products", "permutation_test"),
                               cl = parClust) {
    check_for_polystest(fulldata)
    Data <- assay(fulldata, "quant")
    NumReps <- metadata(fulldata)$NumReps
    NumCond <- metadata(fulldata)$NumCond
    if (!all(statTests %in% c("limma", "Miss_Test", "t_test", "rank_products", "permutation_test"))) {
      stop("Invalid statistical test specified. Available tests are:\n         limma, Miss_Test, t_test, rank_products, permutation_test")
    }
    tests <- statTests
    MAData <- PolySTest:::create_ratio_matrix(fulldata, allComps)
    MAReps <- rep(seq_len(nrow(allComps)), NumReps)
    conditions <- unique(colData(fulldata)$Condition)
    NumComps <- nrow(allComps)
    RRCateg <- matrix(NA, ncol = nrow(allComps), nrow = 2)
    RRCateg[1, ] <- match(allComps[, 1], conditions)
    RRCateg[2, ] <- match(allComps[, 2], conditions)
    p_values <- q_values <- matrix(NA, nrow = nrow(MAData), ncol = length(statTests)*NumComps)
    rownames(p_values) <- rownames(q_values) <- rownames(MAData)
    colnames(p_values) <- paste0("p_values_", rep(statTests, 
                                                  each = NumComps), "_", rep(seq_len(NumComps), length(statTests)))
    colnames(q_values) <- paste0("q_values_", rep(statTests, 
                                                  each = NumComps), "_", rep(seq_len(NumComps), length(statTests)))
    test_funcs <- list(
      limma = function(tMAData) {
        limma_out <- limma_paired(tMAData, NumComps, NumReps)
        list(pvals = limma_out$plvalues, qvals = limma_out$qlvalues, 
             Sds = limma_out$Sds)
      },
      Miss_Test = function(tMAData) {
        MissingStats <- MissingStatsDesign(Data, RRCateg, NumCond, 
                                           NumReps)
        list(pvals = MissingStats$pNAvalues, qvals = MissingStats$qNAvalues)
      },
      t_test = function(tMAData) {
        ttest_out <- ttest_paired(tMAData)
        list(pvals = ttest_out$ptvalues, qvals = ttest_out$qtvalues)
      },
      rank_products = function(tMAData) {
        RPMAUp_pvalues <- RPStats(tMAData, NumReps)
        RPMADown_pvalues <- RPStats(-tMAData, NumReps)
        ttt <- rowMins(cbind(RPMAUp_pvalues, RPMADown_pvalues), 
                       na.rm = TRUE) * 2
        ttt[ttt > 1] <- 1
        tqs <- rep(NA, length(ttt))
        tqs[!is.na(ttt)] <- p.adjust(na.omit(ttt), method = "BH")
        list(pvals = ttt, qvals = tqs)
      },
      permutation_test = function(tMAData) {
        perm_out <- permtest_pairedMod(tMAData, cl = cl)
        list(pvals = perm_out$pPermutvalues, qvals = perm_out$qPermutvalues)
      })
    Sds <- NULL
    for (test in statTests[statTests %in% names(test_funcs)]) {
      if (test %in% c("limma", "Miss_Test")) {
        message("Running ", test, " test")
        res <- test_funcs[[test]](MAData)
        p_values[, grep(paste0("p_values_", test), colnames(p_values))] <- res$pvals
        q_values[, grep(paste0("q_values_", test), colnames(q_values))] <- res$qvals
        if (test == "limma") 
          Sds <- res$Sds
        message(test, " completed")
      }
    }
    lratios <- NULL
    pb <- txtProgressBar(0.9, NumCond)
    for (vs in seq_len(NumComps)) {
      if (!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::setProgress(0.1 + 0.3/NumComps * vs, detail = paste("tests for comparison", vs, "of", NumComps))
      }
      tMAData <- MAData[, MAReps == vs, drop = FALSE]
      for (test in statTests[statTests %in% c("t_test", "rank_products", "permutation_test")]) {
        res <- test_funcs[[test]](tMAData)
        p_values[, grep(paste0("p_values_", test), colnames(p_values))[vs]] <- res$pvals
        q_values[, grep(paste0("q_values_", test), colnames(q_values))[vs]] <- res$qvals
      }
      lratios <- cbind(lratios, rowMeans(MAData[, MAReps == vs], na.rm = TRUE))
      setTxtProgressBar(pb, vs)
    }
    message("tests completed")
    close(pb)
    fulldata <- PolySTest:::prepare_output_data(fulldata, p_values, q_values, lratios, tests, allComps)
    return(fulldata)
  }
  #
  if (Nested) {
    dat <- PolySTest_paired(dat, allComps)
  } else {
    dat <- PolySTest_unpaired(dat, allComps)
  }
  compNames <- metadata(dat)$compNames
  if (!dir.exists(ohDeer2)) { dir.create(ohDeer2, recursive = TRUE) }
  fl <- paste0(ohDeer2, "/All comparisons.jpeg")
  jpeg(filename = fl)
  plotPvalueDistr(dat, compNames)
  dev.off()
  for (comp in compNames) { #comp <- compNames[1]
    comp2 <- gsub("_vs_", " vs ", comp)
    kol <- paste0(PolySTestRoot, comp2)
    ttl <- paste(proteoCraft::cleanNms(unlist(strsplit(comp, "_vs_"))), collapse = " vs ")
    fl <- paste0(ohDeer2, "/", ttl, ".jpeg")
    jpeg(filename = fl)
    plotPvalueDistr(dat, compNames)
    dev.off()
    myData[[kol]] <- ""
    tmp <- as.data.frame(sapply(BH.FDR, function(Sphinx) { #Sphinx <- BH.FDR[1]
      rs <- rep("", nrow(myData))
      w <- as.integer(row.names(dat)[which(rowData(dat)[, paste0("FDR_PolySTest_", comp)] <= Sphinx)])
      if (length(w)) { rs[w] <- "+" }
      return(rs)
    }))
    tmp <- do.call(paste, c(tmp, sep = "/"))
    myData[[kol]] <- tmp
  }
}

# SAM and EBAM
# For SAM we get a P-value, which can be plotted - although this is not the main objective. It should not be too different
# from the ones we get from t- or F-tests (although calculation details may vary, and siggenes does not appear to support one-tail tests),
# but the SAM procedure itself can be used to identify regulated genes.
#
tmpVal2 <- as.matrix(tmpVal)
tmpVal2[which(!is.finite(tmpVal2), arr.ind = TRUE)] <- NA
tmpVal2[which(is.nan(tmpVal2), arr.ind = TRUE)] <- NA
tmpVal2 <- as.data.frame(tmpVal2)
# Careful: siggenes truncates names at 50 characters!!!
# So we will just create artificial short names
rownames(tmpVal2) <- paste0(namesRoot, as.character(1:nrow(tmpVal2)))
require(siggenes)
clusterCall(parClust, function() library(siggenes))
if (exists("samThresh")) { rm(samThresh) }
for (taest in c("SAM", "EBAM")) { #taest <- "SAM" #taest <- "EBAM"
  tstTst <- try({
    tmp <- lapply(1:nrow(expContr), function(x) { #x <- 1
      nm <- expContr$name[x]
      nm2 <- cleanNms(nm)
      em <- expContr$Map[x][[1]]
      em <- em[order(em$Reference, em$Replicate),] # Redundant, but good safety check
      # Important because siggenes assumes replicates start at 1!
      u <- unique(em$Replicate)
      em$Replicate <- match(em$Replicate, u)
      #
      k0 <- em$Expression_Column[which(em$Reference)]
      k1 <- em$Expression_Column[which(!em$Reference)]
      tmpVal3 <- tmpVal2[, c(k1, k0)]
      tst <- parApply(parClust, tmpVal3, 1, function(x) { length(proteoCraft::is.all.good(x)) })
      # Create filter for missing values appropriate for the current test:
      if (taest == "SAM") {
        wh <- which(tst > 1)
      }
      if (taest == "EBAM") {
        wh <- which(tst == ncol(tmpVal3))
      }
      # Test
      if (length(wh)) {
        if (Nested) {
          clss <- c(-em$Replicate[which(em$Reference)],
                    em$Replicate[which(!em$Reference)])
        } else {
          clss <- ifelse(em$Reference, 0, 1)
        }
        if (taest == "SAM") {
          samK <- paste0(samRoots, nm)
          tmpSIGGENES <- siggenes::sam(tmpVal3[wh,], clss, method = "d.stat", s.alpha = seq(0, 1, 0.01)#, gene.names = rownames(tmpVal3)[wh]
                                       , control = samControl(delta = (1:1000)/1000) # Because I like overkill
                                       , rand = mySeed)
          #print(tmpSIGGENES)
          #summary(tmpSIGGENES)
          # NB: The current siggenes vignette has a typo: obviously FDR isn't equal to p0*False/FDR but to p0*False/Called
          if ("s0" %in% slotNames(tmpSIGGENES)) { s0 <- slot(tmpSIGGENES, "s0") } else { s0 <- numeric(0) }
          #si <- median(apply(tmpVal3[wh,], 1, sd, na.rm = TRUE))
          pVal <- -log10(tmpSIGGENES@p.value)
          RES <- data.frame(A = rep(NA, nrow(myData)),
                            B = rep(NA, nrow(myData)))
          colnames(RES) <- samK
          RES[wh, samK[1]] <- pVal
        }
        if (taest == "EBAM") {
          ebamK <- paste0(ebamRoot, nm)
          a0 <- suppressWarnings(find.a0(tmpVal3[wh,], clss))
          tmpSIGGENES <- siggenes::ebam(a0, delta = 1-BH.FDR, rand = mySeed)
          #print(tmpSIGGENES)
          #summary(tmpSIGGENES)
        }
        clusterExport(parClust, c("nm2", "tmpSIGGENES", "taest"), envir = environment())
        tmpSig <- setNames(parLapply(parClust, bhFDRs, function(f) { #f <- bhFDRs[1]
          d <- siggenes::findDelta(tmpSIGGENES, fdr = f, prec = 10)
          if (!is.null(d)) {
            if ("numeric" %in% class(d)) {
              d <- data.frame(Delta = d["Delta"],
                              Called = d["Called"],
                              FDR = d["FDR"])
            } else {
              d <- as.data.frame(d)
            }
            dr <- c(samDir, ebamDir)[match(taest, c("SAM", "EBAM"))]
            XLfun <- eval(parse(text = paste0("siggenes::", tolower(taest), "2excel")), envir = .GlobalEnv)
            dlt <- max(d$Delta[which(d$FDR <= f)])
            if (!length(dlt)) {
              warning(paste0(taest, ": poor Delta estimate for group ", nm2, " at ", 100*f, "% FDR"))
              dlt <- 0
            }
            XLfun(tmpSIGGENES, dlt, paste0(dr, "/", taest, " - ", nm2, ", ", f, " FDR.csv"), what = "both")
            jpeg(filename = paste0(dr, "/Delta plot - ", nm2, ", ", f, " FDR.jpeg"))
            plot(tmpSIGGENES, dlt)
            dev.off()
            res <- siggenes::list.siggenes(tmpSIGGENES, dlt)
            res <- list(Delta.est = d,
                        Delta = dlt,
                        Result = res)
          } else {
            res <- list(Delta.est = NA,
                        Delta = NA,
                        Result = c())
          }
          return(res)
        }), paste0(as.character(bhFDRs), "FDR"))
        tmpSig2 <- as.data.frame(sapply(names(tmpSig), function(f) { #f <- names(tmpSig)[1]
          x <- tmpSig[[f]]$Result
          y <- rep("", nrow(myData))
          if (length(x)) {
            m <- match(x, rownames(tmpVal3))  
            y[m] <- "+"
          }
          return(y)
        }))
        tmpSig3 <- do.call(paste, c(tmpSig2, sep = "/"))
        tmpSig2[[namesCol]] <- myData[[namesCol]]
        D <- sapply(names(tmpSig), function(f) { tmpSig[[f]]$Delta })
        F <- sapply(names(tmpSig), function(f) { tmpSig[[f]]$Delta.FDR })
        D <- data.frame(d = setNames(D, NULL),
                        FDR = as.numeric(gsub("FDR$", "", names(D))))
        D <- aggregate(D$FDR, list(D$d), max)
        colnames(D) <- c("D", "FDR")
        #unique(tmpSig3)
        if (taest == "SAM") {
          RES[[samK[2]]] <- tmpSig3
          RES <- list(Val = RES,
                      siggenesOut = tmpSIGGENES,
                      decision = tmpSig2,
                      S0 = s0,
                      #Si = si,
                      d = D,
                      degFr = length(c(k0, k1))-2)
        }
        if (taest == "EBAM") {
          RES <- data.frame(tmpSig3)
          colnames(RES) <- ebamK
          RES <- list(Val = RES)
        }
      }
      return(RES)
    })
  }, silent = TRUE)
  if (!"try-error" %in% class(tstTst)) {
    tmp <- cbind.data.frame(lapply(tstTst, function(x) { x$Val }))
    myData[, colnames(tmp)] <- tmp
    if (taest == "SAM") {
      samThresh <- setNames(lapply(tstTst, function(x) { x[c("siggenesOut", "decision", "S0", #"Si",
                                                             "d", "degFr")] }), expContr$name)
    }
  }
}

# Edge ODP and LRT
bioc_req <- unique(c(bioc_req, "edge"))
biocInstall("edge")
#
Coeff2 <- Coefficients
if (length(Coefficients) > 1) {
  comb <- gtools::combinations(length(Coefficients), 2, Coefficients)
  synTst <- apply(comb, 1, function(x) {
    x1 <- expMap[[x[[1]]]]
    x2 <- expMap[[x[[2]]]]
    tst <- aggregate(x2, list(x1), function(x) { length(unique(x)) })
    return((max(tst$x) == 1)&(min(tst$x) == 1))
  })
  Coeff2Rmov <- comb[which(synTst), 2]
  Coeff2 <- Coefficients[which(!Coefficients %in% Coeff2Rmov)]
} else { Coeff2Rmov <- c() }
#comb <- comb[which(!synTst),]
wAdj <- which(!Coefficients %in% c("Replicate", VPAL$names, paste0(c("Replicate", VPAL$names), "___"), Coeff2Rmov))
tstAdj <- length(wAdj) > 0
vpals <- aggregate(expMap$Reference, list(expMap[[VPAL$column]]), unique)
vpals <- vpals$Group.1[which(!vpals$x)]
kol <- paste0(intRef, expMap$Ref.Sample.Aggregate)
tmpDat <- myData[, kol]
nr <- nrow(myData)
# clusterExport(parClust, list("tmpDat", "expMap", "RG", "VPAL", "Coefficients", "Coeff2", "intRef", "tstAdj", "wAdj", "WorkFlow", "nr",
#                              "odpRoot", "lrtRoot", "Exp"),
#               envir = environment())
pVals <- try(lapply(vpals, function(vpal) { #vpal <- vpals[1] #vpal <- vpals[2] #vpal <- vpals[3]
  # Does not seem to work when run as parLapply, I could not figure out the issue
  grp <- unique(expMap[which(expMap[[VPAL$column]] == vpal), RG$column])
  em2 <- expMap[c(which((expMap[[RG$column]] == grp)&(expMap$Reference)),
                  which(expMap[[VPAL$column]] == vpal)),]
  k1 <- paste0(intRef, em2$Ref.Sample.Aggregate[which(!em2$Reference)])
  k0 <- paste0(intRef, em2$Ref.Sample.Aggregate[which(em2$Reference)])
  k1 <- k1[which(k1 %in% colnames(tmpDat))]
  k0 <- k0[which(k0 %in% colnames(tmpDat))]
  tmpVal <- tmpDat[, c(k0, k1)]
  tmpVal2 <- as.matrix(tmpVal)
  tmpVal2[which(!is.finite(tmpVal2), arr.ind = TRUE)] <- NA
  tmpVal2[which(is.nan(tmpVal2), arr.ind = TRUE)] <- NA
  tst1 <- apply(tmpVal2[, k1], 1, function(x) { length(proteoCraft::is.all.good(x)) > 1 })
  tst0 <- apply(tmpVal2[, k0], 1, function(x) { length(proteoCraft::is.all.good(x)) > 1 })
  wOK <- which(tst1&tst0)
  tmpVal2 <- tmpVal2[wOK,]
  colnames(tmpVal2) <- gsub(proteoCraft::topattern(intRef), "", colnames(tmpVal2))
  colnames(tmpVal2) <- proteoCraft::cleanNms(colnames(tmpVal2))
  #
  bioVar <- magrittr::set_rownames(as.matrix(em2[, Coeff2[which(Coeff2 != "Time.point")], drop = FALSE]), NULL)
  bioVar <- bioVar[, which(apply(bioVar, 2, function(x) { length(unique(x)) }) > 1), drop = FALSE]
  stopifnot(ncol(bioVar) > 0)
  bsCall <- "de_obj <- edge::build_study(data = tmpVal2, bio.var = bioVar)"
  if (tstAdj) {
    adjVar <- magrittr::set_rownames(as.matrix(em2[, Coefficients[wAdj], drop = FALSE]), NULL)
    adjVar <- adjVar[, which(apply(adjVar, 2, function(x) { length(unique(x)) }) > 1), drop = FALSE]
    if (ncol(adjVar)) {
      bsCall <- gsub("\\)$", ", adj.var = adjVar)", bsCall)
      covCall <- paste0("cov <- data.frame(", paste0(Coefficients[wAdj], " = ", Coefficients[wAdj], collapse = ", "), ")")
      eval(parse(text = covCall), envir = .GlobalEnv)
    }
  }
  if (WorkFlow == "TIMECOURSE") {
    bsCall <- gsub("\\)$", ", tme = em2$Time.point, sampling = \"timecourse\", basis.df = 4)", bsCall)
  }
  #cat(bsCall)
  eval(parse(text = bsCall), envir = .GlobalEnv)
  # Not used, good for checking:
  full_model <- edge::fullModel(de_obj) # Not used, good for checking
  null_model <- edge::nullModel(de_obj) # Not used, good for checking
  full_matrix <- fullMatrix(de_obj)
  null_matrix <- nullMatrix(de_obj)
  #
  #pData(de_obj)
  #getMethod("fit_models", signature = "deSet")
  res <- as.data.frame(matrix(rep(NA, nr*2), ncol = 2))
  #
  # optimal discovery procedure
  odpTmp <- try({
    # Fit
    odp_fit <- edge::fit_models(de_obj, stat.type = "odp")
    # Currently this is often failing and I do not know why!!!
    #de_odp <- try(edge::odp(de_obj, bs.its = 50, seed = mySeed), silent = TRUE)
    de_odp <- try(edge::odp(de_obj, odp_fit, bs.its = 50, n.mods = 50, verbose = TRUE, seed = mySeed), silent = TRUE)
    if ("try-error"%in% class(de_odp)) {
      de_odp <- edge::odp(de_obj, odp_fit, bs.its = 50, n.mods = 50, verbose = TRUE, seed = mySeed, lambda = 0)
    }
    #summary(de_odp)
    # lambda = 0 is from troubleshooting error:
    # "Error in smooth.spline(lambda, pi0, df = smooth.df) : missing or infinite values in inputs are not allowed",
    # cf. https://support.bioconductor.org/p/105623/
  }, silent = TRUE)
  odpTst <- (!"try-error" %in% class(odpTmp))
  if (odpTst) {
    #qvals1 <- de_odp@qvalueObj$qvalues
    pvals1 <- de_odp@qvalueObj$pvalues
    #lfdr1 <- de_odp@qvalueObj$lfdr
    #pi01 <- de_odp@qvalueObj$pi0
    res[[1]][wOK] <- -log10(pvals1)
  } else { cat(odpTmp[1], "\n") }
  #
  # likelihood ratio test
  lrtTmp <- try({
    # Fit
    lrt_fit <- edge::fit_models(de_obj, stat.type = "lrt")
    de_lrt <- edge::lrt(de_obj, lrt_fit, seed = mySeed)
  }, silent = TRUE)
  lrtTst <- (!"try-error" %in% class(lrtTmp))
  if (lrtTst) {
    #qvals2 <- de_lrt@qvalueObj$qvalues
    pvals2 <- de_lrt@qvalueObj$pvalues
    #lfdr2 <- de_lrt@qvalueObj$lfdr
    #pi02 <- de_lrt@qvalueObj$pi0
    res[[2]][wOK] <- -log10(pvals2)
  } else { cat(lrtTmp[1], "\n") }
  colnames(res) <- paste0(c(odpRoot, lrtRoot), vpal)
  res <- res[, which(c(odpTst, lrtTst)), drop = FALSE]
  return(res)
}), silent = TRUE)
if (!"try-error" %in% class(pVals)) {
  pVals <- do.call(cbind, pVals)
  if (ncol(pVals)) {
    myData[, colnames(pVals)] <- pVals
  } else {
    warning("The ODP and LRT procedures failed, investigate...")
  }
}

# Make sure that the data is numeric and not a matrix!!!
meanCol <- paste0("Mean ", intRef, VPAL$values)
welchCol <- paste0("Welch's t-test -log10(Pvalue) - ", VPAL$values)
permCol <- paste0("Permutations t-test -log10(Pvalue) - ", VPAL$values)
modCol <- paste0("Moderated t-test -log10(Pvalue) - ", VPAL$values)
samCol <- paste0("SAM -log10(Pvalue) - ", VPAL$values)
meanCol <- meanCol[which(meanCol %in% colnames(myData))]
welchCol <- welchCol[which(welchCol %in% colnames(myData))]
permCol <- permCol[which(permCol %in% colnames(myData))]
modCol <- modCol[which(modCol %in% colnames(myData))]
samCol <- samCol[which(samCol %in% colnames(myData))]
for (k in c(meanCol, welchCol, permCol, modCol, samCol)) { myData[[k]] <- as.numeric( myData[[k]]) }
#
# Assign results
if (dataType == "modPeptides") {
  ptmpep <- myData
  if (!exists("PTMs_SAM_thresh")) { PTMs_SAM_thresh %<o% list() }
  if (exists("samThresh")) {
    PTMs_SAM_thresh[[Ptm]] <- samThresh
  }
}
if (dataType == "PG") {
  PG <- myData
  if (exists("samThresh")) {
    SAM_thresh %<o% samThresh
  }
}
