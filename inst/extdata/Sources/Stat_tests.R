# Perform statistical tests
# (except ANOVA and SAINTexpress, which have their own sources)

bhFDRs %<o% sort(BH.FDR, decreasing = FALSE)
samRoots %<o% c(sub(" -log10\\(", " ", sub("\\) - $", " - ", samRoot)),
                samRoot,
                paste0("SAM regulated-FDR=", paste(100*bhFDRs, collapse = "/"), "% FDR - "))
samSubDir %<o% "Reg. analysis/SAM"
ebamSubDir %<o% "Reg. analysis/EBAM"
ebamRoot %<o% paste0("EBAM regulated-FDR=", paste(100*bhFDRs, collapse = "/"), "% FDR - ")
#
samDir <- paste0(wd, "/", samSubDir)
ebamDir <- paste0(wd, "/", ebamSubDir)
if (!dir.exists(samDir)) { dir.create(samDir, recursive = TRUE) }
if (!dir.exists(ebamDir)) { dir.create(ebamDir, recursive = TRUE) }
if (!exists("limmaFits")) { limmaFits <- list() }
limmaFits %<o% limmaFits
if ((!dataType %in% names(limmaFits)) || (!inherits(limmaFits[[dataType]], "list"))) {
  limmaFits[[dataType]] <- list()
}
#
source(parSrc, local = FALSE)
#dataType <- "modPeptides" #dataType <- "PG"
if (dataType == "modPeptides") {
  myData <- ptmpep
  intRef <- pepRf
  namesCol <- "Name"
  namesRoot <- "Pep"
  ohDeer <- paste0(wd, "/Reg. analysis/", ptm, "/t-tests")
}
if (dataType == "peptides") {
  myData <- pep
  intRef <- pep.ref[length(pep.ref)]
  namesCol <- "Modified sequence"
  namesRoot <- "Pep"
  ohDeer <- paste0(wd, "/Reg. analysis/t-tests")
  # pep is not log-transformed at this stage, we need to log(10)-transform it (same as modPeptides and PG)!
  # We must then make sure that this modified quant does not get pushed back to the original table.
  # (Note: for limma the log base is adjusted to 2)
  kol1 <- intersect(paste0(intRef, RSA$values), colnames(myData))
  intRef <- paste0("log10 ", intRef) 
  kol2 <- paste0(intRef, RSA$values)
  myData[, kol2] <- log10(myData[, kol1])
  myData <- myData[, which(!colnames(myData) %in% kol1)]
}
if (dataType == "PG") {
  myData <- PG
  intRef <- #if (exists("prExpr_roots")) { prExpr_roots["Quantitation"] } else {
    Prot.Expr.Root #}
  namesCol <- "Leading protein IDs"
  namesRoot <- "PG"
  ohDeer <- paste0(wd, "/Reg. analysis/t-tests")
}
quantCol <- paste0(intRef, RSA$values)
if (!dir.exists(ohDeer)) { dir.create(ohDeer, recursive = TRUE) }

# Statistical tests:
cat(" -> Statistical tests\n")
#
# - Moderated t-test (limma) and modified DEqMS version:
#     We will run two versions:
#      - a classic limma moderated t-test
#      - the DEqMS-corrected version
#
bioc_req <- unique(c(bioc_req, "DEqMS"))
biocInstall("DEqMS")
cran_req <- unique(c(cran_req, "matrixStats"))
if (!require(matrixStats, quietly = TRUE)) { pak::pak("matrixStats") }
library(matrixStats)
TESTs <- c("limma", "DEqMS")[1L:((dataType == "PG") + 1L)]
limpaMode <- (dataType != "modPeptides") && (quantAlgo == "limpa")
for (TEST in TESTs) { #TEST <- TESTs[1L] #TEST <- TESTs[2L]
  if ((!TEST %in% names(limmaFits[[dataType]])) || (!inherits(limmaFits[[dataType]][[TEST]], "list"))) {
    limmaFits[[dataType]][[TEST]] <- list()
  }
  cat(paste0("   - ", TEST, " moderated t-test\n"))
  if (TEST == "limma") {
    pRoot <- modRoot
    insrt <- "limma mod. t-test"
  }
  if (TEST == "DEqMS") {
    pRoot <- deqmsRoot
    insrt <- "DEqMS mod. t-test"
  }
  if (limpaMode) {
    tmpData <- quantData_list$EList_obj
  } else {
    # Note: limpa can also be used for peptides... but that is for another day!
    wOK <- 1L:nrow(myData)
    if (TEST == "DEqMS") {
      grpTst <- aggregate(1L:nrow(expMap), list(expMap[[VPAL$limmaCol]]), list)
      tmpFl <- tempfile(fileext = ".rds")
      readr::write_rds(myData[, quantCol], tmpFl)
      clusterExport(parClust, list("grpTst", "is.all.good", "quantCol", "tmpFl"), envir = environment())
      invisible(clusterCall(parClust, \() {
        assign("tmp", readr::read_rds(tmpFl), envir = .GlobalEnv)
      }))
      tst <- t(parSapply(parClust, 1L:nrow(myData), \(x) {
        vapply(grpTst$x, \(y) { length(is.all.good(tmp[x, unlist(y)])) }, 1L)
      }))
      tst <- apply(tst, 1L, min)
      wOK <- which(tst >= 2L)
    }
    #
    tmpData <- myData[wOK, quantCol]/log10(2L) # For limma, use log2 data!!!
    rownames(tmpData) <- myData[wOK, namesCol]
  }
  if (Nested) {
    if (!RSA$limmaCol %in% colnames(expMap)) {
      expMap[[RSA$limmaCol]] <- do.call(paste, c(expMap[, paste0(RSA$names, "_._"), drop = FALSE], sep = "___"))
    }
    Block <- expMap[match(rownames(designMatr), gsub("___", "_", as.character(expMap[[RSA$limmaCol]]))),
                    Blocking.factors$limmaCol]
    corfit <- tryCatch(
      duplicateCorrelation(tmpData,
                           designMatr,
                           block = Block),
      warning = \(w) {
        message(w$message)
        return(list(consensus = 0))
      }
    )
    fit <- if (limpaMode) {
      limpa::dpcDE(tmpData, designMatr, block = Block)
    } else {
      lmFit(tmpData, designMatr, block = Block, correlation = corfit$consensus)
    }
  } else {
    fit <- if (limpaMode) {
      limpa::dpcDE(tmpData, designMatr)
    } else {
      lmFit(tmpData, designMatr)
    }
  }
  #NB: voomaLmFit doesn't work as a replacement for lmFit here unless we either:
  # - impute,
  # - use limpa (but then we should use dpcDE instead!)
  # - filter out any row with even a single NA!
  fit$genes <- rownames(tmpData)
  limmaFits[[dataType]][[TEST]]$fitted_values <- fitted(fit)
  fit <- contrasts.fit(fit, contrMatr)
  fit2 <- eBayes(fit) # Note: eBayes() performs empirical Bayes moderation of variances; the default settings are appropriate, and e.g. do not assume 1% of proteins are differentially expressed.
  fit2$genes <- fit$genes
  # Do not use topTable!
  # - It adjusts P values (with "BH") by default...
  #   If using it (e.g. rewriting to include systematic P-value adjustment), set adjust.method = "none"!!!
  # - The way it tests a hypothesis seems better suited to F.tests than t-tests!
  # Here you could use decideTests() if wanting to switch to directly using decisions from limma
  if (TEST == "limma") {
    limmaFits[[dataType]][[TEST]]$fit <- fit2
    kols1 <- paste0(sub(" -log10\\(", " ", sub("\\) - $", " - ", deqmsRoot)), colnames(contrMatr))
    kols2 <- paste0(deqmsRoot, colnames(contrMatr))
    myData[, c(kols1, kols2)] <- NA_real_
    w <- which(myData[[namesCol]] %in% fit2$genes)
    m <- match(myData[w, namesCol], fit2$genes)
    myData[w, kols1] <- fit2$p.value[m,]
    myData[w, kols2] <- -log10(fit2$p.value[m,])
    #View(myData[w, c(kols1, kols2)])
    # Plot moderated t-test results
    for (contr in colnames(contrMatr)) { #contr <- colnames(contrMatr)[1L]
      # Q-Q plot
      ttl <- paste0("mod. t-test QQ plot - ", contr)
      fl <- paste0(ohDeer, "/", ttl)
      jpeg(file = paste0(fl, ".jpeg"), width = 400L, height = 350L)
      qqt(as.data.frame(fit2$t)[[contr]], df = fit2$df.prior + fit2$df.residual, pch = 16L, cex = 0.2)
      abline(0, 1)
      dev.off()
      # MA plot
      ttl <- paste0(insrt, " MA plot - ", contr)
      fl <- paste0(ohDeer, "/", ttl)
      jpeg(file = paste0(fl, ".jpeg"), width = 400L, height = 350L)
      fit2$genes <- NULL
      plotMD(fit2, column = match(contr, colnames(fit2$p.value)))
      dev.off()
    }
  }
  if (TEST == "DEqMS") {
    # DEqMS:
    # This is meant as en extension to limma to account "for variance dependence on the number of quantified peptides or PSMs
    # for statistical testing of differential protein expression."
    # https://www.bioconductor.org/packages/devel/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html#overview-of-deqms
    if (dataType == "PG") {
      countCol <- grep("^Evidences count - ", colnames(myData), value = TRUE)
    }
    if (dataType == "modPeptides") {
      idCol <- grep("^Evidence IDs - ", colnames(myData), value = TRUE)
      countCol <- gsub("^Evidence IDs - ", "Evidences count - ", idCol)
      myData[, countCol] <- sapply(idCol, \(k) { lengths(strsplit(myData[[k]], ";")) })
    }
    psm.counts <- data.frame(count = rowMins(as.matrix(myData[, countCol])))
    rownames(psm.counts) <- myData[[namesCol]]
    psm.counts$count[which(is.na(psm.counts$count))] <- 0L
    psm.counts$count <- psm.counts$count+1L # Pseudo-counts, to allow for PGs with 0 PSMs
    # # We need a least 1 and will treat 0 and 1 together
    # w <- which(psm.counts$count == 0L)
    # psm.counts$count[w] <- 1L
    fit2$count <- rep(1L, nrow(fit2))
    w <- which(rownames(fit2$coefficients) %in% row.names(psm.counts))
    m <- match(rownames(fit2$coefficients)[w], row.names(psm.counts))
    fit2$count[w] <- psm.counts$count[m]
    fit3 <- try(spectraCounteBayes(fit2), silent = TRUE)
    if (!inherits(fit3, "try-error")) {
      # Variance plot
      tst <- try({
        ttl <- paste0(insrt, " variance plot")
        fl <- paste0(ohDeer, "/", ttl)
        jpeg(file = paste0(fl, ".jpeg"), width = 400L, height = 350L)
        VarianceBoxplot(fit3, n = max(psm.counts), main = ttl, xlab = "PSM count")
        dev.off()
      }, silent = TRUE)
      #
      limmaFits[[dataType]][[TEST]]$fit <- fit3
      #
      for (i in 1L:ncol(contrMatr)) { #i <- 1L
        contrNm <- colnames(contrMatr)[i]
        DEqMS.results <- outputResult(fit3, coef_col = i)
        kol1 <- paste0(sub(" -log10\\(", " ", sub("\\) - $", " - ", deqmsRoot)), contrNm)
        kol2 <- paste0(deqmsRoot, contrNm)
        myData[[kol1]] <- myData[[kol2]] <- NA_real_
        w <- which(myData[[namesCol]] %in% row.names(DEqMS.results))
        m <- match(myData[[namesCol]][w], row.names(DEqMS.results))
        myData[w, kol1] <- DEqMS.results$sca.P.Value[m]
        myData[w, kol1] <- -log10(DEqMS.results$sca.P.Value[m])
      }
    }
  }
}

# msqrob2
if ((dataType == "PG") && ("QFeatures_obj" %in% names(quantData_list))) {
  if (!require(msqrob2)) { pak::pak("msqrob2") }
  library(msqrob2)
  if (!require(ExploreModelMatrix)) { pak::pak("ExploreModelMatrix") }
  library(ExploreModelMatrix)
  cat("   - MSqRob test\n")
  #
  # Supplement experiment structure in the object
  #msqrobForm2 <- as.formula(limmaForm)
  msqrobForm <- paste0("~ 0 + ", VPAL$limmaCol)
  if (hasBatch) {
    msqrobForm <- paste0(msqrobForm, " + ", Batch.effect$limmaCol)
  }
  msqrobForm2 <- as.formula(msqrobForm)
  inDat <- quantData_list$QFeatures_obj
  colDat <- inDat@colData
  colDat$Sample <- sub(".* - ", "", rownames(inDat@colData))
  colDat[, Factors] <- Exp.map[match(colDat$Sample, Exp.map$Ref.Sample.Aggregate), Factors]
  for (fct in Factors) { colDat[[fct]] <- factor(colDat[[fct]], levels = FactorsLevels[[fct]]) }
  colDat[[RSA$limmaCol]] <- expMap[match(colDat$Sample, rownames(expMap)), RSA$limmaCol]
  colDat[[VPAL$limmaCol]] <- expMap[match(colDat$Sample, rownames(expMap)), VPAL$limmaCol]
  if (exists("Batch.effect") && length(Batch.effect$values)) {
    colDat[[Batch.effect$limmaCol]] <- expMap[match(colDat$Sample, rownames(expMap)), Batch.effect$limmaCol]
  }
  inDat@colData <- colDat
  # Model
  qf <- msqrob(inDat,
               i = "PG",
               formula = msqrobForm2,
               robust = TRUE)
  #
  models <- rowData(qf[["PG"]])[["msqrobModels"]]
  #models[[1L]]@params
  # vd <- ExploreModelMatrix::VisualizeDesign(sampleData = colData(qf),
  #                                           designFormula = msqrobForm2,
  #                                           textSizeFitted = 4)
  # vd$plotlist
  #
  designMsqrob <- model.matrix(msqrobForm2, colData(qf))
  whSingle <- which(myContrasts$Secondary == "")
  whDouble <- which(myContrasts$Secondary != "")
  Lit <- c("A", "B", "C", "D")
  w <- which(Lit %in% colnames(myContrasts))
  tmp <- myContrasts[, Lit[w]]
  nms <- VPAL$names
  if (length(Exp) == 1L) { nms <- setdiff(nms, "Experiment") }
  fct <- paste(paste0(nms, "_._"), collapse = "_")
  for (lit in Lit[w]) {
    tmp[[lit]] <- paste0(fct, tmp[[lit]])
  }
  myContrasts$msqrob2 <- ""
  myContrasts$msqrob2[whSingle] <- do.call(paste, c(tmp[whSingle, c("A", "B")], sep = " - "))
  if (length(whDouble)) {
    tmp1 <- do.call(paste, c(tmp[whDouble, c("A", "B")], sep = " - "))
    tmp2 <- do.call(paste, c(tmp[whDouble, c("C", "D")], sep = " - "))
    myContrasts$msqrob2[whDouble] <- paste0("(", tmp1, ") - (", tmp2, ")")
  }
  contrMsqrob <- limma::makeContrasts(contrasts = myContrasts$msqrob2,
                                      levels = designMsqrob)
  qf <- hypothesisTest(qf,
                       i = "PG",
                       contrast = contrMsqrob)
  #
  MSqRob_infer %<o% lapply(colnames(contrMsqrob), \(k) {
    x <- rowData(qf[["PG"]])[[k]]
    colnames(x) <- paste0(colnames(x), " - ", k)
    return(x)
  })
  MSqRob_infer <- do.call(cbind, MSqRob_infer)
  kol1 <- colnames(MSqRob_infer) <- paste0("MSqRob ", colnames(MSqRob_infer))
  tmp <- sub(" - .*", "", kol1)
  kol2 <- vapply(strsplit(kol1, " - "), \(x) {
    paste(x[2L:length(x)], collapse = " - ")
  }, "")
  kol2 <- myContrasts$Contrast[match(kol2, myContrasts$msqrob2)]
  kol2 <- paste0(tmp, " - ", kol2)
  colnames(MSqRob_infer) <- kol2
  #View(MSqRob_infer)
  w <- which(myData[[namesCol]] %in% rownames(MSqRob_infer))
  m <- match(myData[w, namesCol], rownames(MSqRob_infer))
  myData[, kol2] <- NA
  myData[w, kol2] <- MSqRob_infer[m, kol2]
  pval <- paste0("MSqRob pval - ", myContrasts$Contrast)
  myData[, sub(topattern("MSqRob pval - "),
               sub(" -log10\\(", " - ", sub("\\) - $", " - ", msqrobRoot)), pval)] <- myData[, pval]
  myData[, sub(topattern("MSqRob pval - "), msqrobRoot, pval)] <- -log10(myData[, pval])
  pvalue.col["MSqRob"] <- msqrobRoot
}

# Others (classic Student's and Welch's t-test, permutations test) and ROTC
# In general any test which is not compatible with complex experiment designs (deals only with 2 groups) should be run here
pairwise_coin_test <- \(data,
                        alternative = "two.sided",
                        skipBlocks = TRUE # Otherwise in most cases we cannot run the test, because blocks = replicates and N = 1 -> no permutations!
                        ) {
  formTxt <- "values ~ group"
  data$group <- as.factor(data$group)
  tst <- c("batch", "block") %in% colnames(data)
  if (skipBlocks) { tst[2L] <- FALSE }
  if (sum(tst)) {
    if (tst[1L]) { data$batch <- as.factor(data$batch) }
    if (tst[2L]) { data$block <- as.factor(data$block) }
    formTxt <- paste0(formTxt, " | ", paste(c("batch", "block")[which(tst)], collapse = " + "))
  }
  form <- as.formula(formTxt)
  # permutation‑based Wilcoxon/Mann–Whitney–type test (allows for covariates, and optionally for blocking)
  tst <- try({ coin::pvalue(coin::independence_test(form, data = data, alternative = alternative)) }, silent = TRUE)
  if (inherits(tst, "try-error")) { return(NA) }
  return(tst)
}
# Export data to cluster
tmpFl <- tempfile(fileext = ".rds")
exports <- list("tmpFl", "Nested", "hasBatch", "pairwise_coin_test", "expMap")
if (Nested) {
  parse.Param.aggreg.2("Blocking.factors")
  blockCol <- Blocking.factors$limmaCol
  exports <- append(exports, "blockCol")
}
hasBatch %<o% (("Batch.effect" %in% colnames(Param))&&(nchar(Param$Batch.effect)))
if (hasBatch) {
  parse.Param.aggreg.2("Batch.effect")
  batchCol <- Batch.effect$limmaCol
  exports <- append(exports, "batchCol")
}
clusterExport(parClust, exports, envir = environment())
tmpData <- myData[, grep(topattern(intRef), colnames(myData), value = TRUE)]
rownames(tmpData) <- myData[[namesCol]]
readr::write_rds(tmpData, tmpFl)
invisible(clusterCall(parClust, \() {
  tmpData <- readr::read_rds(tmpFl)
  assign("tmpData", tmpData, envir = .GlobalEnv)
  return()
}))
unlink(tmpFl)
whSingle <- which(!myContrasts$isDouble)
if (length(whSingle)) {
  #
  # - Student's and Welch's t-tests and permutations test
  cat("   - t-tests and permutation tests\n")
  logRoots <- c(StudentRoot, WelchRoot, permRoot)
  rawRoots <- sub(" -log10\\(", " ", sub("\\) - $", " - ", logRoots))
  clusterExport(parClust, "rawRoots", envir = environment())
  tmpTsts <- lapply(whSingle, \(i) { #i <- 1L #i <- 2L #i <- 3L
    # Get two groups from the contrast
    A <- myContrasts$A[[i]]
    B <- myContrasts$B[[i]]
    # Extract sub-map
    em <- expMap[which(expMap[[VPAL$limmaCol]] %in% c(A, B)),]
    # Check that columns are in input data
    em <- em[which(paste0(intRef, rownames(em)) %in% colnames(tmpData)),]
    if (!nrow(em)) { return(matrix(rep(NA_real_, 3L*nrow(tmpData)), ncol = 3L)) }
    if (Nested) {
      uBlck <- unique(em[[blockCol]])
      tmp <- lapply(uBlck, \(x) {
        #list( # No: each block should contain exactly 1 sample per group!
        c(rownames(em)[which((em[[VPAL$limmaCol]] == A)&(em[[blockCol]] == x))],
          rownames(em)[which((em[[VPAL$limmaCol]] == B)&(em[[blockCol]] == x))])
      })
      #w <- which(vapply(tmp, \(x) { sum(lengths(x) >= 1L) }, 1L) == 2L)
      w <- which(lengths(tmp) == 2L) # Each block should contain exactly 1 sample per group!
      tmp <- tmp[w]
      tmp <- do.call(rbind, tmp)
      A_ <- tmp[, 1L]
      B_ <- tmp[, 2L]
      em <- em[match(c(A_, B_), rownames(em)),]
    } else {
      A_ <- rownames(em)[which(em[[VPAL$limmaCol]] == A)]
      B_ <- rownames(em)[which(em[[VPAL$limmaCol]] == B)]
    }
    if ((length(A_) < 2L)||(length(B_) < 2L)) { return(rep(NA, 3L)) }
    A_k <- paste0(intRef, A_)
    B_k <- paste0(intRef, B_)
    nms <- paste0(rawRoots, myContrasts$Contrast[[i]])
    altHyp <- c("two.sided", "greater")[myContrasts[i, "Up-only"]+1L]
    dt <- data.frame(sample = c(A_, B_),
                     group = c(rep(A, length(A_)), rep(B, length(B_))))
    if (Nested) { dt$block <- em[[blockCol]] }
    if (hasBatch) { dt$batch <- em[[batchCol]] }
    clusterExport(parClust, list("em", "A", "B", "A_k", "B_k", "altHyp", "dt", "Nested"), envir = environment())
    RES <- parLapply(parClust, 1L:nrow(tmpData), \(ii) { #ii <- 1L
      vA <- as.numeric(tmpData[ii, A_k])
      vB <- as.numeric(tmpData[ii, B_k])
      dt$values <- c(vA, vB)
      if (Nested) {
        w <- which((is.finite(vA))&(is.finite(vB)))
        w <- c(w, w+length(vA))
      } else {
        w <- which(is.finite(dt$values))
      }
      dt <- dt[w,]
      vA <- dt$values[which(dt$group %in% A)]
      vB <- dt$values[which(dt$group %in% B)]
      if ((length(unique(vA)) < 2L) || (length(unique(vB)) < 2L)) { return(rep(NA, 3L)) }
      tst1 <- try(t.test(x = vB, y = vA, paired = Nested, alternative = altHyp, var.equal = TRUE)$p.value, silent = TRUE)
      tst2 <- try(t.test(x = vB, y = vA, paired = Nested, alternative = altHyp, var.equal = FALSE)$p.value, silent = TRUE)
      tst3 <- try(pairwise_coin_test(dt, altHyp), silent = TRUE)
      if (inherits(tst1, "try-error")) { tst1 <- NA_real_ }
      if (inherits(tst2, "try-error")) { tst2 <- NA_real_ }
      if (inherits(tst3, "try-error")) { tst3 <- NA_real_ }
      res <- c(tst1,
               tst2,
               tst3)
      return(res)
    })
    RES <- do.call(rbind, RES)
    colnames(RES) <- nms
    return(RES)
  })
  tmpTsts <- as.data.frame(do.call(cbind, tmpTsts))
  kol <- kol12raw <- c(paste0(rawRoots[1L], myContrasts$Contrast[whSingle]),
                       paste0(rawRoots[2L], myContrasts$Contrast[whSingle]))
  kolLog <- kol12log <- c(paste0(logRoots[1L], myContrasts$Contrast[whSingle]),
                          paste0(logRoots[2L], myContrasts$Contrast[whSingle]))
  tmpTsts_ <- tmpTsts[, kol12raw]
  kol3raw <- paste0(rawRoots[3L], myContrasts$Contrast[whSingle])
  kol3log <- paste0(logRoots[3L], myContrasts$Contrast[whSingle])
  tstPerm <- vapply(kol3raw, \(x) { sum(is.finite(tmpTsts[, x])) }, 1L)
  if (!max(tstPerm)) {
    warning("The permutations tests failed.")
  } else {
    tmpTsts_[, kol3raw] <- tmpTsts[, kol3raw]
    kolRaw <- c(kol12raw, kol3raw)
    kolLog <- c(kol12log, kol3log)
  }
  # To check results
  checkTTests <- FALSE
  if (checkTTests) {
    tmp <- dfMelt(tmpTsts_)
    tmp <- tmp[which(is.finite(tmp$value)),]
    tmp[, c("test", "contrast")] <- do.call(rbind, strsplit(as.character(tmp$variable), " Pvalue - "))
    plot <- ggplot(tmp) + geom_density(stat = "density", aes(x = value, color = test)) +
      facet_grid(test~contrast) + theme_bw() + ggtitle("Classic tests: dist. of P-values")
    poplot(plot)
  }
  #
  myData[, kolRaw] <- tmpTsts_[, kolRaw]
  myData[, kolLog] <- -log10(tmpTsts_[, kolRaw])
  #View(tmpTsts_)
  #summary(tmpTsts_)
  #
  # - Reproducibility-Optimized Test Statistic (ROTS)
  cat("   - ROTS tests\n")
  bioc_req <- union(bioc_req, "ROTS")
  if (!require(ROTS)) { pak::pak("ROTS") }
  library(ROTS)
  ROTS_res %<o% setNames(lapply(whSingle, \(i) { #i <- 1L #i <- 2L #i <- 3L
    # Get two groups from the contrast
    A_ <- myContrasts$A_samples[[i]]
    B_ <- myContrasts$B_samples[[i]]
    groups <- c(rep(0L, length(A_)), rep(1L, length(B_)))
    A_k <- paste0(intRef, A_)
    B_k <- paste0(intRef, B_)
    tst <- as.data.frame(sapply(list(A = A_k,
                                     B = B_k), \(kol) {
                                       apply(tmpData[, kol], 1L, \(x) { sum(is.finite(x)) }) >= 2L
                                     }))
    w <- which(rowSums(tst) == 2)
    return(ROTS::ROTS(data = tmpData[w, c(A_k, B_k)]/log10(2L),
                      groups = groups,
                      B = 100,
                      K = 500,
                      seed = mySeed))
  }), myContrasts$Contrast[whSingle])
  #View(ROTS_res[[myContrasts$Contrast[whSingle[1L]]]])
  root <- sub(" -log10\\(", " ", sub("\\) - $", " - ", rotsRoot))
  myData[, paste0(root, myContrasts$Contrast[whSingle])] <- NA_real_
  myData[, paste0(rotsRoot, myContrasts$Contrast[whSingle])] <- NA_real_
  for (contr in myContrasts$Contrast[whSingle]) { #contr <- myContrasts$Contrast[whSingle[1L]]
    w <- which(myData[[namesCol]] %in% rownames(ROTS_res[[contr]]$data))
    myData[w, paste0(root, contr)] <- ROTS_res[[contr]]$pvalue[match(myData[w, namesCol],
                                                                     rownames(ROTS_res[[contr]]$data))]
    myData[w, paste0(rotsRoot, contr)] <- -log10(myData[w, paste0(root, contr)])
  }
  #
  #
  # SAM and EBAM
  # For SAM we get a P-value, which can be plotted.
  # It should not be too different from the ones we get from t- or F-tests
  # (although calculation details may vary, and siggenes does not appear to support one-tail tests).
  #
  tmpData2 <- as.matrix(tmpData)
  tmpData2[which(!is.finite(tmpData2), arr.ind = TRUE)] <- NA_real_
  tmpData2[which(is.nan(tmpData2), arr.ind = TRUE)] <- NA_real_
  tmpData2 <- as.data.frame(tmpData2)
  # Careful: siggenes truncates names at 50 characters!!!
  # So we will just create artificial short names
  rownames(tmpData2) <- paste0(namesRoot, as.character(1L:nrow(tmpData2)))
  require(siggenes)
  clusterExport(parClust, list("samDir", "ebamDir"), envir = environment())
  invisible(clusterCall(parClust, \() {
    library(siggenes)
    return()
  }))
  if (exists("samThresh")) { rm(samThresh) }
  for (taest in c("SAM", "EBAM")) { #taest <- "SAM" #taest <- "EBAM"
    cat(paste0("   - ", taest, " test\n"))
    tstTst <- try({
      tmp <- lapply(whSingle, \(x) { #x <- whSingle[1L]
        nm <- myContrasts$Contrast[x]
        em <- expMap
        em <- em[order(em$Replicate_._),]
        em <- em[c(which(em[[VPAL$limmaCol]] %in% myContrasts$B[[x]]),
                   which(em[[VPAL$limmaCol]] %in% myContrasts$A[[x]])),]
        # Important because siggenes assumes replicates start at 1!
        u <- unique(em$Replicate_._)
        em$Replicate <- match(em$Replicate_._, u)
        #
        wA <- which(em[[VPAL$limmaCol]] %in% myContrasts$A[[x]])
        wB <- which(em[[VPAL$limmaCol]] %in% myContrasts$B[[x]])
        A_ <- rownames(em)[wA]
        B_ <- rownames(em)[wB]
        A_k <- paste0(intRef, A_)
        B_k <- paste0(intRef, B_)
        tmpData3 <- tmpData2[, c(B_k, A_k)]
        tst <- parApply(parClust, tmpData3, 1L, \(x) { sum(is.finite(x)) })
        # Create filter for missing values appropriate for the current test:
        if (taest == "SAM") {
          wh <- which(tst > 1L)
        }
        if (taest == "EBAM") {
          wh <- which(tst == ncol(tmpData3))
        }
        # Test
        if (length(wh)) {
          if (Nested) {
            clss <- c(-em$Replicate[wB],
                      em$Replicate[wA])
          } else {
            clss <- ifelse(1L:nrow(em) %in% wB, 0L, 1L)
          }
          if (taest == "SAM") {
            samK <- paste0(samRoots, nm)
            msg <- capture.output({
              tmpSIGGENES <- siggenes::sam(tmpData3[wh,], clss, method = "d.stat", s.alpha = seq(0, 1, 0.01)#, gene.names = rownames(tmpData3)[wh]
                                           , control = samControl(delta = (1:1000)/1000) # Because I like overkill
                                           , rand = mySeed)
            })
            #print(tmpSIGGENES)
            #summary(tmpSIGGENES)
            # NB: The current siggenes vignette has a typo: obviously FDR isn't equal to p0*False/FDR but to p0*False/Called
            s0 <- if ("s0" %in% slotNames(tmpSIGGENES)) { slot(tmpSIGGENES, "s0") } else { numeric(0) }
            #si <- median(apply(tmpData3[wh,], 1L, sd, na.rm = TRUE))
            nr <- nrow(myData)
            RES <- data.frame(A = rep(NA, nr),
                              B = rep(NA, nr),
                              C = rep(NA, nr),
                              D = rep(NA, nr))
            colnames(RES) <- samK
            RES[wh, samK[1L]] <- tmpSIGGENES@p.value
            RES[wh, samK[2L]] <- -log10(tmpSIGGENES@p.value)
          }
          if (taest == "EBAM") {
            ebamK <- paste0(ebamRoot, nm)
            msg <- capture.output({ a0 <- find.a0(tmpData3[wh,], clss) })
            tmpSIGGENES <- siggenes::ebam(a0, delta = 1-BH.FDR, rand = mySeed)
            #print(tmpSIGGENES)
            #summary(tmpSIGGENES)
          }
          clusterExport(parClust, c("nm", "tmpSIGGENES", "taest"), envir = environment())
          tmpSig <- setNames(parLapply(parClust, BH.FDR, \(f) { #f <- BH.FDR[1L]
            d <- siggenes::findDelta(tmpSIGGENES, fdr = f, prec = 10)
            if (!is.null(d)) {
              if (is.numeric(d)) {
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
                warning(paste0(taest, ": poor Delta estimate for group ", nm, " at ", 100*f, "% FDR"))
                dlt <- 0
              }
              XLfun(tmpSIGGENES, dlt, paste0(dr, "/", taest, " - ", nm, ", ", f, " FDR.csv"), what = "both")
              jpeg(filename = paste0(dr, "/Delta plot - ", nm, ", ", f, " FDR.jpeg"))
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
          }), paste0(as.character(BH.FDR), "FDR"))
          tmpSig2 <- as.data.frame(sapply(names(tmpSig), \(f) { #f <- names(tmpSig)[1L]
            x <- tmpSig[[f]]$Result
            y <- rep("", nrow(myData))
            if (length(x)) {
              m <- match(x, rownames(tmpData3))  
              y[m] <- "+"
            }
            return(y)
          }))
          tmpSig3 <- do.call(paste, c(tmpSig2, sep = "/"))
          tmpSig2[[namesCol]] <- myData[[namesCol]]
          D <- sapply(names(tmpSig), \(f) { tmpSig[[f]]$Delta })
          dF <- sapply(names(tmpSig), \(f) { tmpSig[[f]]$Delta.FDR })
          D <- data.frame(d = setNames(D, NULL),
                          FDR = as.numeric(gsub("FDR$", "", names(D))))
          D <- aggregate(D$FDR, list(D$d), max)
          colnames(D) <- c("D", "FDR")
          #unique(tmpSig3)
          if (taest == "SAM") {
            RES[[samK[3L]]] <- tmpSig3
            RES <- list(Val = RES,
                        siggenesOut = tmpSIGGENES,
                        decision = tmpSig2,
                        S0 = s0,
                        #Si = si,
                        d = D,
                        degFr = length(c(B_, A_))-2L)
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
    if (!inherits(tstTst, "try-error")) {
      tmp <- cbind.data.frame(lapply(tstTst, \(x) { x$Val }))
      myData[, colnames(tmp)] <- tmp
      if (taest == "SAM") {
        samThresh <- setNames(lapply(tstTst, \(x) { x[c("siggenesOut", "decision", "S0", #"Si",
                                                        "d", "degFr")] }), myContrasts$Contrast)
      }
    }
  }
}

# Make sure that the data is numeric and not a matrix!!!
kol <- grep("-log10\\(Pvalue\\) - ", colnames(myData), value = TRUE)
for (k in kol) { myData[[k]] <- as.numeric( myData[[k]]) }
#
# Assign results
if (!exists("samThresh")) { samThresh <- list() }
#
# 
myData <- myData[, grep(topattern(intRef), colnames(myData), value = TRUE, invert = TRUE)] # because in some cases we have changed this data!
#
if (dataType == "modPeptides") {
  ptmpep[, colnames(myData)] <- myData
  if (!exists("PTMs_SAM_thresh")) { PTMs_SAM_thresh <- list() }
  PTMs_SAM_thresh %<o% PTMs_SAM_thresh
  PTMs_SAM_thresh[[Ptm]] <- samThresh
}
if (dataType == "peptides") {
  pep[, colnames(myData)] <- myData
  SAM_thresh %<o% samThresh
}
if (dataType == "PG") {
  PG[, colnames(myData)] <- myData
  SAM_thresh %<o% samThresh
}
