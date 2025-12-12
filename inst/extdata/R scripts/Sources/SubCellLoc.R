#### Sub-Cellular localisation analysis
if (Annotate&&LocAnalysis) {
  msg <- "Sub-Cellular localisation analysis:"
  ReportCalls <- AddMsg2Report(Space = FALSE)
  dir <- paste0(wd, "/pRoloc")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  packs <- c("pRoloc", "pRolocGUI") # We will not use pRolocGUI as it seems to be using an outdated version of shinydashboardPlus (specifically a retired function)
  for (pack in packs) { #pack <- packs[1]
    bioc_req <- unique(c(bioc_req, pack))
    biocInstall(pack)
  }
  SubCellFracAggr %<o% parse.Param.aggreg(Param_filter(Param$Ratios.Groups.Ref.Aggregate.Level, "Com"))
  # Impute missing values
  AllKol <- paste0(prtRfRoot, RSA$values)
  w <- which(AllKol %in% colnames(PG))
  grps <- Exp.map[match(RSA$values[w], Exp.map$Ref.Sample.Aggregate), VPAL$column] 
  AllKol <- AllKol[w]
  wNC <- which(PG$`Potential contaminant` != "+")
  tempDat <- PG[wNC, AllKol]
  tempDat <- Data_Impute2(tempDat, grps)$Imputed_data
  #
  pRolocData <- SVMparams <- list()
  # For sub-cellular localisation analysis, we need to define a series of compartment markers to predict protein location
  # We can use the very granular markers already defined above (pRoloc or built-in), or define new ones here
  ObjNm <- "CompGOTerms2"
  .obj <- unique(c(ObjNm, .obj))
  if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
    msg <- "Enter a list of GO Cell Compartment (GO CC) terms for compartments of interest (semicolon-separated).
You may also include:
 - 1 for a low-resolution markers list (nucleoplasm - chromatin - cytoplasm)
 - 2 for a more exhaustive list.

Example: \"GO:0031012;2\"
"
    tmp <- unlist(strsplit(dlg_input(msg, 2)$res, "[;,] ?"))
    if ("1" %in% tmp) { tmp <- unique(c(tmp, "GO:0005654", "GO:0000785", "GO:0005737")) }
    if ("2" %in% tmp) { tmp <- unique(c(tmp, CompGOTerms)) }
    tmp <- tmp[which(tmp %in% GO_terms$ID[which(GO_terms$Ontology == "CC")])] # (Also neatly removes "1" and "2"...)
    ObjNm %<c% tmp
    AllAnsw <- AllAnsw[which(AllAnsw$Parameter != ObjNm),]
    tmp <- AllAnsw[1,]
    tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
    tmp$Value <- list(get(ObjNm))
    m <- match(ObjNm, AllAnsw$Parameter)
    if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
  }
  if (length(CompGOTerms2)) {
    CompGOTerms2 <- data.frame(Term = CompGOTerms2)
    CompGOTerms2$Offspring <- lapply(CompGOTerms2$Term, function(x) { GOCCOFFSPRING[[x]] })
    allOffspr <- unlist(CompGOTerms2$Offspring)
    allOffspr <- aggregate(allOffspr, list(allOffspr), length)
    CompGOTerms2$Offspring <- lapply(CompGOTerms2$Offspring, function(x) {
      x[which(!x %in% c(CompGOTerms2$Term, allOffspr$Group.1[which(allOffspr$x > 1)]))]
    })
    CompGOTerms2$Name <- gsub(" \\[GO:[0-9]{7}\\]$", "", GO_terms$Term[match(CompGOTerms2$Term, GO_terms$ID)])
    CompGOTerms2$All <- apply(CompGOTerms2[, c("Term", "Offspring")], 1, function(x) { unique(unlist(x)) })
    tst <- setNames(lapply(CompGOTerms2$All, function(x) { grep(paste(x, collapse = "|"), PG$"GO-ID") }),
                    CompGOTerms2$Name)
    tst2 <- unique(unlist(tst))
    tst2 <- tst2[which(vapply(tst2, function(x) { sum(vapply(tst, function(y) { x %in% y }, TRUE)) }, 1) == 1)]
    tst <- lapply(tst, function(x) { x[which(x %in% tst2)] })
    tst <- tst[which(vapply(tst, length, 1) > 0)]
    SubCellMark2 %<o% listMelt(tst) # Overwrite former value
    SubCellMark2$value <- PG$Label[SubCellMark2$value]
    SubCellMark2 <- setNames(SubCellMark2$L1, SubCellMark2$value)
    if (length(SubCellMark2)) {
      msg <- " - Performing per-sample pRoloc analysis"
      ReportCalls <- AddMsg2Report(Space = FALSE)
      ttls <- c()
      pRolocVisMeth <- "t-SNE"
      tmpPG <- PG[, c("id", "Leading protein IDs", "Protein IDs", "Protein names", "Genes", "Label")]
      exports <- list("SubCellMark2", "Exp.map", "SubCellFracAggr", "prtRfRoot", "tempDat", "VPAL", "tmpPG", "Aggregates",
                      "prot.list", "pRolocVisMeth", "Exp", "dir", "wNC")
      clusterExport(parClust, exports, envir = environment())
      invisible(clusterCall(parClust, function() {
        library(graphics)
        library(Biobase)
        library(MSnbase)
        library(pRoloc)
        library(proteoCraft)
        return()
      }))
      f0 <- function(grp) { #grp <- SubCellFracAggr$values[1]
        ttl_s <- c()
        grp1 <- cleanNms(grp)
        w <- which(Exp.map[[SubCellFracAggr$column]] == grp)
        RES <- list(Success = FALSE)
        if (length(w)) {
          em <- Exp.map[w,]
          kol <- paste0(prtRfRoot, unique(em$Ref.Sample.Aggregate))
          kol <- kol[which(kol %in% colnames(tempDat))]
          temp1 <- tempDat[, kol]
          # Add a small noise to avoid non-unicity issue
          sd <- sd(unlist(temp1))
          noise <- rnorm(length(unlist(temp1)), 0, sd*10^-6)
          temp1 <- temp1 + noise
          #
          temp1 <- 10^temp1
          tst <- rowSums(temp1)
          wAG <- which((!is.na(tst))&(is.finite(tst))&(tst > 0))
          if (length(wAG)) {
            #cat(paste0(grp1, ":\n", paste(rep("-", nchar(grp1)+1), collapse = ""), "\n"))
            m <- match(unique(em[[VPAL$column]]), em[[VPAL$column]])
            temp1 <- sweep(temp1[wAG,], 1, tst[wAG], "/")
            rownames(temp1) <- tmpPG$Label[wAG]
            colnames(temp1) <- cleanNms(unique(em[[VPAL$column]]))
            temp1 <- as.matrix(temp1)
            temp2 <- em[m, Aggregates]
            rownames(temp2) <- cleanNms(unique(em[[VPAL$column]]))
            temp2 <- AnnotatedDataFrame(temp2)
            temp3 <- AnnotatedDataFrame(tmpPG[wAG, c("id", "Leading protein IDs", "Protein IDs", "Protein names", "Genes")])
            rownames(temp3) <- tmpPG$Label[wAG]
            MSnData <- MSnSet(exprs = temp1,
                              pData = temp2,
                              fData = temp3)
            MSnData <- addMarkers(MSnData, SubCellMark2)
            #pRolocData[[grp]] <<- MSnData
            #getMarkers(MSnData)
            # Define features of interest
            plotPLmtchs <- FALSE
            if (length(prot.list)) {
              PLmtchs <- grsep2(prot.list, tmpPG$"Leading protein IDs"[wAG])
              if (length(PLmtchs)) {
                plotPLmtchs <- TRUE
                foi1 <- FeaturesOfInterest(description = "Proteins of interest",
                                           fnames = featureNames(MSnData)[PLmtchs])
                description(foi1)
                foi(foi1)
              }
            }
            # the plotDist function is unsatisfactory
            #plotDist(MSnData, featureNames(MSnData))
            #
            # Markers hierarchical dendrogram
            ttl <- paste0("Markers hier. clust. - ", grp1)
            ttl2 <- paste0("Markers hierarchical clustering\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            mrkHClust(MSnData, main = ttl2)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #
            # Average markers class profile plots
            hc <- mrkHClust(MSnData, plot = FALSE)
            ## order of markers according to histogram
            mm <- getMarkerClasses(MSnData)
            m_order <- levels(factor(mm))[order.dendrogram(hc)]
            ## average marker profile
            fmat <- mrkConsProfiles(MSnData)
            ttl <- paste0("Marker classes av. profiles - ", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plotConsProfiles(fmat, order = m_order)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #
            # Dimensionality reduction
            pRolocVisMeth2 <- pRolocVisMeth
            if (pRolocVisMeth == "sLDA") { pRolocVisMeth2 <- "lda" }
            ttl <- paste0(pRolocVisMeth, " - ", grp1)
            ttl2 <- gsub(" - ", "\n", ttl)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plot2D(MSnData, fcol = "markers", main = ttl2, method = pRolocVisMeth2)
            addLegend(MSnData, fcol = "markers", cex = 0.7, where = "bottomright", ncol = 2)
            if (plotPLmtchs) {
              highlightOnPlot(MSnData, foi1, col = "black", lwd = 2)
              legend("topright", c("Proteins of interest"),
                     bty = "n", col = c("purple"),
                     pch = 1)
            }
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #plot3D(MSnData, fcol = "markers", main = ttl2) # Much worse - but also much older - than plotly
            #
            #pRolocVis(MSnData) # Not used: buggy
            #
            # Assess the resolution of the fractionation for the different compartments
            hlq <- QSep(MSnData)
            ttl <- paste0("Subcell. res. heatmap - ", grp1)
            ttl2 <- paste0("Sub-cellular resolution heatmap\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            levelPlot(hlq, main = ttl2)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            ttl <- paste0("Subcell. res. boxplot - ", grp1)
            ttl2 <- paste0("Sub-cellular resolution boxplot\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plot(hlq, main = ttl2)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #
            # Prediction of compartment assignment
            ## ... Unsupervised (kept as an example here, but doesn't look very useful)
            #kcl <- MLInterfaces::MLearn( ~ ., MSnData,  kmeansI, centers = length(Com))
            #plot(kcl, exprs(MSnData))
            #hcl <- MLInterfaces::MLearn( ~ ., MSnData,
            #                             MLInterfaces::hclustI(distFun = dist,
            #                       cutParm = list(k = length(Com))))
            #plot(hcl, labels = FALSE)
            #pcl <- MLearn( ~ ., MSnData,  pamI(dist), k = length(Com))
            #plot(pcl, data = exprs(MSnData))
            #
            ## ... Supervised
            # This is really, reaaally slow!
            cat(" Predicting compartment localisation using SVM (100 iterations), please wait...\n")
            wghts <- classWeights(MSnData, fcol = "markers")
            wghts <- wghts[which(wghts < 1)]
            params <- suppressMessages(suppressWarnings(svmOptimisation(MSnData, "markers", times = 10, xval = 5,
                                                                        class.weights = wghts, verbose = TRUE)))
            # The line above takes a millenium and a half
            # Unfortunately it cannot be parallelized or optimized easily without rewriting pRoloc
            # But this loop could be parallelized at least
            #
            #SVMparams[[grp]] <- params
            ttl <- paste0("SVM opt. boxplot - ", grp1)
            ttl2 <- paste0("SVM optimisation boxplot\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plot(params, main = ttl2)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            ttl <- paste0("SVM opt. heatmap - ", grp1)
            ttl2 <- paste0("SVM optimisation heatmap\n", grp1)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            levelPlot(params)
            dev.off()
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            #f1Count(params)
            params2 <- getParams(params)
            MSnData <- svmClassification(MSnData, params,
                                         fcol = "markers", class.weights = wghts)
            #processingData(MSnData)
            p1 <- getPredictions(MSnData, fcol = "svm")
            p1 <- fData(p1)$svm.pred
            minprob <- median(fData(MSnData)$svm.scores)
            p2 <- getPredictions(MSnData, fcol = "svm", t = minprob)
            p2 <- fData(p2)$svm.pred
            #table(p1, p2)
            ptsze <- exp(fData(MSnData)$svm.scores) - 1
            ttl <- paste0(pRolocVisMeth, " with predictions - ", grp1)
            ttl2 <- gsub(" - ", "\n", ttl)
            ttl_s <- c(ttl_s, ttl)
            grDevices::pdf(paste0(dir, "/", ttl, ".pdf"), width = 10, height = 10); par(cex.main = 1)
            plot2D(MSnData, fcol = "markers", main = ttl2, method = pRolocVisMeth2)
            addLegend(MSnData, fcol = "markers", cex = 0.7, where = "bottomright", ncol = 2)
            if (plotPLmtchs) {
              highlightOnPlot(MSnData, foi1, col = "black", lwd = 2)
              legend("topright", c("Proteins of interest"),
                     bty = "n", col = c("purple"),
                     pch = 1)
              dev.off()
            }
            #system(paste0("open \"", dir, "/", ttl, ".pdf", "\""))
            lokol <- paste0("Localisation - ", grp)
            #PG[[lokol]] <<- ""
            w <- which(p2 != "unknown")
            #PG[wNC[wAG[w]], lokol] <<- p2[w]
            #pRolocData[[grp]] <<- MSnData
            #cat("\n")
            RES <- list(Success = TRUE,
                        Titles = ttl_s,
                        Which = wNC[wAG[w]],
                        Column = lokol,
                        Values = p2[w],
                        pRolocData = MSnData,
                        SVMparams = params)
          }
        }
        return(RES)
      }
      tmpClass <- parLapply(parClust, SubCellFracAggr$values, f0)
      names(tmpClass) <- SubCellFracAggr$values
      grps <- names(tmpClass)[which(vapply(tmpClass, function(x) { x$Success }, TRUE))]
      for (grp in grps) {
        ttls <- c(ttls, tmpClass[[grp]]$Titles)
        lokol <- tmpClass[[grp]]$Column
        w <- tmpClass[[grp]]$Which
        PG[[lokol]] <- ""
        PG[w, lokol] <- tmpClass[[grp]]$Values
        pRolocData[[grp]] <- tmpClass[[grp]]$MSnData
        SVMparams[[grp]] <- tmpClass[[grp]]$SVMparams
      }
      ReportCalls <- AddSpace2Report()
      if (length(ttls)) {
        SilentPDF2JPEG <- function(ttl) {
          suppressMessages(
            suppressWarnings(
              pdf_convert(paste0(ttl, ".pdf"), "jpeg", filenames = paste0(ttl, ".jpeg"), dpi = 600)
            )
          )
        }
        clusterExport(parClust, list("dir", "SilentPDF2JPEG", "pdf_convert"), envir = environment())
        parSapply(parClust, ttls, function(ttl) { try(SilentPDF2JPEG(paste0(dir, "/", ttl)), silent = TRUE) })
      }
      w <- which(!PG$Label %in% names(SubCellMark2))
      #View(PG[w, grep("^Localisation - ", colnames(PG), value = TRUE)])
    }
  }
  #
  # Re-localisation analysis
  SubCellFracAggr2 %<o% parse.Param.aggreg(Param_filter(Param$Volcano.plots.Aggregate.Level, "Com"))
  SSD.Root %<o% "log10(SSD) - "
  SSD.Pval.Root %<o% "Welch's t-test on SSDs -log10(Pvalue) - "
  WhRef <- lapply(SubCellFracAggr2$values, function(x) {
    unique(Exp.map$Reference[which(Exp.map[[SubCellFracAggr2$column]] == x)])
  })
  tst <- vapply(WhRef, length, 1)
  if (max(tst) == 1) {
    tempDat2 <- 10^tempDat
    wh0 <- which(WhRef)
    wh1 <- which(!WhRef)
    if ((length(wh0) == 1)&&(length(wh1))) {
      msg <- " - Performing re-localisation analysis"
      ReportCalls <- AddMsg2Report(Space = FALSE)
      LocAnalysis2 <- TRUE
      EM0 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == SubCellFracAggr2$values[wh0]),]
      EM0$Replicate <- as.numeric(EM0$Replicate)
      EM0 <- EM0[order(EM0$Replicate, EM0$Compartment),]
      grp0 <- SubCellFracAggr2$values[wh0]
      # Currently, this only supports one type of grouping defined by SubCellFracAggr2,
      # itself equivalent to removing "Compartment" from Sample groups ("VPAL")
      comb <- gtools::combinations(max(EM0$Replicate), 2, unique(EM0$Replicate))
      # SSDs should normalize for expression level, otherwise we will pick up proteins whose expression changes!
      # Thus NormSSDs is always TRUE for now.
      NormSSDs <- TRUE
      if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
        RefSSDs %<o% NULL
      }
      if (Param$Ratios.Thresholds == threshMsg) {
        RefSSDs %<o% as.numeric(apply(comb, 1, function(i) { #i <- comb[1,]
          # NB: differs now from the way Ref.Ratios is written - but for a good reason.
          # We operate using different groupings for subcellular re-localisation analysis.
          A <- tempDat2[, paste0(prtRfRoot, EM0$Ref.Sample.Aggregate[which(EM0$Replicate == i[[1]])])]
          B <- tempDat2[, paste0(prtRfRoot, EM0$Ref.Sample.Aggregate[which(EM0$Replicate == i[[2]])])]
          if (NormSSDs) {
            A <- sweep(A, 1, rowSums(A), "/") 
            B <- sweep(B, 1, rowSums(B), "/") 
          }
          res <- log10(rowSums((A-B)^2))
          return(res)
        }))
        RefSSDs %<o% setNames(lapply(SubCellFracAggr2$values, function(x) { RefSSDs }), SubCellFracAggr2$values)
      }
      SSDs %<o% list()
      SSD.FDR.thresh %<o% c()
      for (wh in wh1) { #wh <- wh1[1]
        grp <- SubCellFracAggr2$values[wh]
        EM1 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == grp),]
        EM1$Replicate <- as.numeric(EM1$Replicate)
        EM1 <- EM1[order(EM1$Replicate, EM1$Compartment),]
        # Calculate sum of squared differences
        if (grepl("Rep", SubCellFracAggr$aggregate)) {
          rps <- unique(EM1$Replicate)
          SSDs[[grp]] <- as.data.frame(sapply(rps, function(rp) {
            P0 <- tempDat2[, paste0(prtRfRoot, EM0$Ref.Sample.Aggregate[which(EM0$Replicate == rp)])]
            P1 <- tempDat2[, paste0(prtRfRoot, EM1$Ref.Sample.Aggregate[which(EM1$Replicate == rp)])]
            if (NormSSDs) {
              P0 <- sweep(P0, 1, rowSums(P0), "/") 
              P1 <- sweep(P1, 1, rowSums(P1), "/") 
            }
            res <- P0-P1
            res <- log10(rowSums((P0-P1)^2))
            return(res)
          }))
          colnames(SSDs[[grp]]) <- SSDkols <- paste0(SSD.Root, EM1[match(rps, EM1$Replicate), SubCellFracAggr$column])
        } else {
          comb <- gtools::permutations(max(as.numeric(Rep)), 2, as.numeric(Rep), repeats.allowed = TRUE)
          SSDs[[grp]] <- as.data.frame(apply(comb, 1, function(i) {
            P0 <- tempDat2[, paste0(prtRfRoot, EM0$Ref.Sample.Aggregate[which(EM0$Replicate == i[[1]])])]
            P1 <- tempDat2[, paste0(prtRfRoot, EM1$Ref.Sample.Aggregate[which(EM1$Replicate == i[[2]])])]
            if (NormSSDs) {
              P0 <- sweep(P0, 1, rowSums(P0), "/") 
              P1 <- sweep(P1, 1, rowSums(P1), "/") 
            }
            res <- log10(rowSums((P0-P1)^2))
            return(res)
          }))
          colnames(SSDs[[grp]]) <- SSDkols <- paste0(SSD.Root, apply(comb, 1, function(i) {
            paste0(EM1[match(i[[2]], EM1$Replicate), SubCellFracAggr$column], "_vs_", EM0[match(i[[1]], EM0$Replicate), SubCellFracAggr$column])
          }))
        }
        SSDs[[grp]][[paste0("Mean ", SSD.Root, grp)]] <- apply(SSDs[[grp]][, SSDkols, drop = FALSE], 1, function(x) { log10(mean(10^x)) })
        PG[, colnames(SSDs[[grp]])] <- NA
        PG[wNC, colnames(SSDs[[grp]])] <- SSDs[[grp]]
      }
      # Check distribution
      # Test expression values:
      g <- c(grep(topattern("log10(SSD) - "), colnames(PG), value = TRUE),
             grep(topattern("Mean log10(SSD) - "), colnames(PG), value = TRUE))
      test <- PG[, g]
      colnames(test) <- gsub(topattern("log10(SSD) - ", start = FALSE), "", colnames(test))
      w <- grep("^Mean ", colnames(test))
      colnames(test)[w] <- paste0(gsub("^Mean ", "", colnames(test)[w]), "___Mean")
      test <- test[which(apply(test, 1, function(x) { length(is.all.good(x)) }) > 0),]
      test <- suppressMessages(melt.data.frame(test))
      test$variable <- as.character(test$variable)
      test[, SubCellFracAggr$names] <- ""
      w <- rep(FALSE, nrow(test))
      test[which(!w), SubCellFracAggr$names] <- Isapply(strsplit(test$variable[which(!w)], "___"), unlist)
      a <- SubCellFracAggr$names
      w <- which(vapply(a, function(x) { length(unique(test[[x]])) }, 1) > 1)
      if (length(w)) { a <- a[w] }
      test[[a[1]]] <- factor(test[[a[1]]], levels = sort(unique(test[[a[1]]])))
      test <- test[which(is.all.good(test$value, 2)),]
      test2 <- set_colnames(aggregate(test$value, list(test$variable), median), c("variable", "value"))
      test2[, a] <- test[match(test2$variable, test$variable), a]
      MinMax <- c(min(test$value), max(test$value))
      nbinz <- ceiling((MinMax[2]-MinMax[1])/0.1)
      binz <- c(0:nbinz)/nbinz
      binz <- binz*(MinMax[2] - MinMax[1]) + MinMax[1]
      binz[1] <- binz[1] - 0.000001
      testI <- data.frame(Intensity = (binz[2:(nbinz+1)]+binz[seq_len(nbinz)])/2)
      for (v in unique(test$variable)) {
        wv <- which(test$variable == v)
        testI[[v]] <- vapply(seq_len(nbinz), function(x) {
          sum((test$value[wv] > binz[x])&(test$value[wv] <= binz[x+1]))
        }, 1)
      }
      testI <- reshape2::melt(testI, id.vars = "Intensity")
      testI[, a] <- test[match(testI$variable, test$variable), a]
      testI$variable <- cleanNms(testI$variable)
      dir <- paste0(wd, "/Reg. analysis/Localisation")
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      ttl <- "Distribution of log10(SSD) values"
      plot <- ggplot(testI) +
        geom_area(aes(x = Intensity, y = value, fill = variable, group = variable,
                      colour = variable), alpha = 0.25) +
        scale_color_viridis_d(begin = 0.25) +
        scale_fill_viridis_d(begin = 0.25) +
        geom_vline(data = test2, aes(xintercept = value), linetype = "dashed", color = "grey") +
        ggtitle(ttl) + theme_bw() + theme(legend.position = "none", strip.text.y = element_text(angle = 0)) +
        scale_y_continuous(limits = c(0, max(testI$value)*1.1), expand = c(0, 0))
      if (length(a) == 1) { plot <- plot + facet_wrap(as.formula(paste0("~", a))) } else {
        plot <- plot + facet_grid(as.formula(paste0(a[1], "~", paste(a[2:length(a)], collapse = "+"))))
      }
      poplot(plot, 12, 22)
      suppressMessages({
        ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
        ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")
      })
      ReportCalls <- AddPlot2Report()
      # Statistical test
      M <- median(unlist(PG[wNC, grep(topattern(SSD.Root), colnames(PG), value = TRUE)]))
      for (wh in wh1) { #wh <- wh1[1]
        grp <- SubCellFracAggr2$values[wh]
        SSDs[[grp]][[paste0(SSD.Pval.Root, grp)]] <- apply(SSDs[[grp]][, SSDkols], 1, function(x) {
          -log10(t.test(x, alternative = "greater", mu = M)$p.value)
        })
        temp <- FDR(SSDs[[grp]], grp, SSD.Pval.Root, fdr = BH.FDR, returns = c(TRUE, TRUE), method = "BH")
        SSDs[[grp]][, gsub("^Significant-", "Signif. SSDs-", colnames(temp$`Significance vector`))] <- temp$`Significance vector`
        SSD.FDR.thresh <- c(SSD.FDR.thresh, temp$Thresholds)
        PG[, colnames(SSDs[[grp]])] <- NA
        PG[wNC, colnames(SSDs[[grp]])] <- SSDs[[grp]]
      }
      temp <- PG[, grep("^Regulated - ", colnames(PG), value = TRUE, invert = TRUE)]
      subDr <- "Reg. analysis/Localisation"
      tempVP3 <- Volcano.plot(Prot = temp,
                              mode = "custom",
                              experiments.map = Exp.map,
                              X.root = paste0("Mean ", SSD.Root),
                              Y.root = SSD.Pval.Root,
                              aggregate.map = Aggregate.map,
                              aggregate.name = SubCellFracAggr2$aggregate,
                              aggregate.list = Aggregate.list, parameters = Param,
                              save = c("jpeg", "pdf"), labels = "FDR",
                              Ref.Ratio.values = RefSSDs,
                              Ref.Ratio.method = paste0("obs", RefRat_Mode),
                              ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                              FDR.thresh = SSD.FDR.thresh, FDR.root = "Signif. SSDs-FDR=",
                              arbitrary.lines = arbitrary.thr,
                              proteins = prot.list, proteins_split = protsplit,
                              return = TRUE, return.plot = TRUE,
                              title = "SSDs volcano plot ",
                              subfolder = subDr,
                              subfolderpertype = FALSE, Symmetrical = FALSE,
                              Alpha = "Rel. log10(Peptides count)",
                              Size = "Av. log10 abundance", Size.max = 2,
                              plotly = create_plotly, plotly_local = create_plotly_local,
                              plotly_labels = PrLabKol,
                              X.normalized = FALSE,
                              cl = parClust)
      if (!class(tempVP3) %in% c("try-error", "character")) {
        #
        # Save plotly plots
        dr <- paste0(wd, "/", subDr)
        myPlotLys <- tempVP3
        Src <- paste0(libPath, "/extdata/R scripts/Sources/save_Volcano_plotlys.R")
        #rstudioapi::documentOpen(Src)
        source(Src, local = FALSE)
        #
        VP_list <- tempVP3
        insrt <- ""
        Src <- paste0(libPath, "/extdata/R scripts/Sources/thresholds_Excel.R")
        #rstudioapi::documentOpen(Src)
        source(Src, local = FALSE)
        #
        g <- grep("Regulated - ", colnames(tempVP3$Protein_groups_file), value = TRUE)
        PG[, gsub("^Regulated - ", "Re-localized - ", g)] <- tempVP3$Protein_groups_file[,g]
        volcano.plots$Localisation_Unlabelled <- tempVP3$Plots$Unlabelled
        volcano.plots$Localisation_Labelled <- tempVP3$Plots$Labelled
        n2 <- names(volcano.plots$Localisation_Labelled)
        dir <- paste0(wd, "/Reg. analysis/Localisation")
        for (ttl in n2) {
          plot <- volcano.plots$Localisation_Labelled[[ttl]]
          ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
        }
        if ((create_plotly)&&(!create_plotly_local)) { plot_ly$"Localisation" <- tempVP3$"Plotly plots" }
        # Edit wording + create filters
        g <- grep("^Re-localized - ", colnames(PG), value = TRUE)
        for (gi in g) { #gi <- g[1]
          PG[which(PG[[gi]] == "too small FC"), gi] <- "unchanged distr."
          wUp <- grep("^up, FDR = ", PG[[gi]])
          if (length(wUp)) {
            PrtWidth <- 5
            dir <- paste0(wd, "/Reg. analysis/Localisation/Relocalised proteins")
            if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
            dirlist <- unique(c(dirlist, dir))
            plotNorm <- FALSE
            grp <- gsub("^Re-localized - ", "", gi)
            EM1 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == grp),]
            EM1$Replicate <- as.numeric(EM1$Replicate)
            EM1 <- EM1[order(EM1$Replicate, EM1$Compartment),]
            for (wup in wUp) { #wup <- wUp[1]
              # Protein heatmap
              Seq <- db$Sequence[match(unlist(strsplit(PG$`Leading protein IDs`[wup], ";"))[1], db$`Protein ID`)]
              Prfl <- data.frame(Sample = c(EM0$Ref.Sample.Aggregate, EM1$Ref.Sample.Aggregate))
              if (plotNorm) {
                PGRoot <- Prot.Expr.Root
                PepRoot <- pep.ref[length(pep.ref)]
              } else {
                PGRoot <- Prot.Expr.Root2
                PepRoot <- pep.ref2
              }
              Prfl$value <- 10^as.numeric(PG[wup, paste0(PGRoot, Prfl$Sample)])
              Prfl[, RSA$names] <- Isapply(strsplit(Prfl$Sample, "___"), unlist)
              Prfl$Replicate <- as.integer(Prfl$Replicate)
              Prfl[, c(SubCellFracAggr$aggregate, SubCellFracAggr2$aggregate)] <- Exp.map[match(Prfl$Sample, Exp.map$Ref.Sample.Aggregate),
                                                                                          c(SubCellFracAggr$column, SubCellFracAggr2$column)]
              Prfl$x <- 0
              Prfl$xend <- Prfl$x + Prfl$value/max(Prfl$value)*PrtWidth
              Prfl$y <- 0:(nrow(Prfl)-1)
              Prfl$"SubCell. Frac." <- factor(Prfl$Compartment, levels = Com)
              Prfl$Entity <- "Protein"
              Prfl$Angle <- 0
              Prfl$Size <- 3
              kol <- paste0(pep.ref[length(pep.ref)], Prfl$Sample)
              kol <- paste0(PepRoot, Prfl$Sample)
              pepPrfl <- pep[match(as.integer(unlist(strsplit(PG$`Peptide IDs`[wup], ";"))), pep$id),
                             c("Modified sequence", "Sequence", kol)]
              Match <- vapply(pepPrfl$Sequence, function(x) { nchar(unlist(strsplit(Seq, x))[1]) }, 1)
              pepPrfl <- pepPrfl[order(Match),]
              pepPrfl[, kol] <- sweep(pepPrfl[, kol], 1, rowMax(as.matrix(pepPrfl[, kol])), "/")
              mSeq <- unique(pepPrfl$"Modified sequence") 
              colnames(pepPrfl) <- gsub(topattern(PepRoot), "", colnames(pepPrfl))
              pepPrfl$Sequence <- NULL
              pepPrfl <- reshape2::melt(pepPrfl, id.vars = "Modified sequence")
              colnames(pepPrfl) <- c("Entity", "Sample", "value")
              pepPrfl$Sample <- as.character(pepPrfl$Sample)              
              pepPrfl[, RSA$names] <- Isapply(strsplit(pepPrfl$Sample, "___"), unlist)
              pepPrfl$Replicate <- as.integer(pepPrfl$Replicate)
              pepPrfl[, c(SubCellFracAggr$aggregate, SubCellFracAggr2$aggregate)] <- Exp.map[match(pepPrfl$Sample, Exp.map$Ref.Sample.Aggregate),
                                                                                             c(SubCellFracAggr$column, SubCellFracAggr2$column)]
              pepPrfl$x <- PrtWidth + match(pepPrfl$Entity, mSeq)
              pepPrfl$xend <- pepPrfl$x + pepPrfl$value
              pepPrfl$y <- Prfl$y[match(pepPrfl$Sample, Prfl$Sample)]
              pepPrfl$"SubCell. Frac." <- factor(pepPrfl$Compartment, levels = Com)
              pepPrfl$Angle <- 60
              pepPrfl$Size <- 2
              pepPrfl <- pepPrfl[which(pepPrfl$value > 0),]
              temp <- rbind(Prfl, pepPrfl)
              temp2 <- aggregate(temp[, c("x", "Angle", "Size")], list(temp$Entity), unique)
              colnames(temp2) <- c("Entity", "X", "Angle", "Size")
              temp3 <- aggregate(temp$y, list(temp[[SubCellFracAggr$aggregate]]), function(x) { list(Min = min(x),
                                                                                                     Max = max(x)+1) })
              colnames(temp3) <- c("Samples group", "x")
              temp3[, c("Min", "Max")] <- apply(temp3$x, 2, as.numeric)
              temp3$Mean <- apply(temp3[, c("Min", "Max")], 1, mean)
              temp3[, SubCellFracAggr$names] <- Isapply(strsplit(as.character(temp3$`Samples group`), "___"), unlist)
              temp3$`Sample group` <- cleanNms(temp3$`Samples group`)
              temp4 <- aggregate(temp$y, list(temp[[SubCellFracAggr2$aggregate]]), function(x) { list(Min = min(x),
                                                                                                      Max = max(x)+1) })
              colnames(temp4) <- c("Samples group 2", "x")
              temp4[, c("Min", "Max")] <- apply(temp4$x, 2, as.numeric)
              temp4$Mean <- apply(temp4[, c("Min", "Max")], 1, mean)
              temp4$"Samples group 2" <- cleanNms(temp4$"Samples group 2")
              temp5 <- aggregate(temp$y, list(temp$Sample), function(x) { list(Min = min(x), Max = max(x)+1) })
              colnames(temp5) <- c("Sample", "x")
              temp5[, c("Min", "Max")] <- apply(temp5$x, 2, as.numeric)
              temp5$Mean <- apply(temp5[, c("Min", "Max")], 1, mean)
              temp5$"SubCell. frac." <- Exp.map$Compartment[match(temp5$Sample, Exp.map$Ref.Sample.Aggregate)]
              NCSc <- max(nchar(temp5$`SubCell. frac.`))
              ttl <- paste0("Subcell. fract. distr. - ", PG$`Leading protein IDs`[wup])
              plot <- ggplot(temp) +
                geom_rect(aes(xmin = x, xmax = xend, ymin = y, ymax = y+1, fill = `SubCell. Frac.`)) +
                scale_fill_viridis_d(begin = 0.25) +
                geom_segment(aes(x = x, xend = x), y = -2, yend = max(temp$y + 1.5), colour = "lightgrey") +
                geom_segment(x = max(temp$x)+1, xend =  max(temp$x)+1, y = -2, yend = max(temp$y + 1.5), colour = "lightgrey") +
                geom_segment(x = PrtWidth + 0.5, xend = PrtWidth + 0.5, y = -2, yend = max(temp$y + 1.5)) +
                geom_hline(yintercept = max(temp$y + 1.5)) +
                geom_text(data = temp2, aes(x = X, label = Entity, angle = Angle, size = Size), y = max(temp$y + 2.5), hjust = 0) +
                geom_text(data = temp4, aes(y = Mean, label = `Samples group 2`), x = -3-NCSc, hjust = 1, vjust = 0.5, size = 2.5) +
                geom_segment(data = temp4, aes(y = Min+0.1, yend = Max-0.1), x = -2.8-NCSc, xend = -2.8-NCSc) +
                geom_text(data = temp3, aes(y = Mean, label = Replicate), x = -2-NCSc, hjust = 1, vjust = 0.5, size = 1.8) +
                geom_segment(data = temp3, aes(y = Min+0.1, yend = Max-0.1), x = -1.8-NCSc, xend = -1.8-NCSc) +
                geom_text(data = temp5, aes(y = Mean, label = `SubCell. frac.`), x = -1, hjust = 1, vjust = 0.5, size = 1.8) +
                coord_fixed() + theme_bw() + scale_size_identity() +
                ggtitle("Distribution across subcellular fractions", subtitle = PG$Label[wup]) +
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(),  axis.line = element_line(colour = "black"),
                      plot.margin = margin(0, 0, 0, 0, "cm"), legend.position = "none") +
                xlab(NULL) + ylab(NULL) + #theme(legend.position = "none") +
                xlim(-5-NCSc, max(temp$xend+1)) + ylim(0, max(temp$y+20))
              #poplot(plot, 12, 22)
              suppressMessages({
                ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 150, width = 20, height = 12, units = "in")
                ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150, width = 20, height = 12, units = "in")
              })
            }
          }
        }
        g1 <- gsub("^Re-localized - ", "", g)
        up <- grep("^reloc\\., FDR = ", unique(unlist(PG[, g])), value = TRUE)
        Reg_filters$Localisation <- list()
        if ("con" %in% filter_types) {
          Reg_filters$Localisation$"By condition" <- list()
          rat <- paste0("Mean ", SSD.Root, g1)
          for (i in seq_along(g)) { #i <- 1
            Reg_filters$Localisation$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                                     Filter_up = sort(which(PG[[g[i]]] %in% up)),
                                                                     Filter_down = c(),
                                                                     Filter = sort(which(PG[[g[i]]] %in% up)),
                                                                     Ratios = PG[[rat[i]]],
                                                                     Background_filter = seq_len(nrow(PG)))
            # To do here: check that Ratios values above are correct
          }
        }
        #
        # Z-scored clustering heatmaps of re-localized proteins
        clustMode <- "re-localisation"
        Src <- paste0(libPath, "/extdata/R scripts/Sources/cluster_Heatmap_Main.R")
        #rstudioapi::documentOpen(Src)
        source(Src, local = FALSE)
        #
      } else { warning("No localisation volcano plots created, investigate!") }
    } else {
      if (!length(wh0)) {
        warning("Skipping re-localisation analysis: no reference group...")
      }
      if (length(wh0) > 0) {
        warning("Skipping re-localisation analysis: too many reference groups, only 1 allowed...")
      }
      if (!length(wh1)) {
        warning("Skipping re-localisation analysis: no non-reference group...")
      }
    }
  } else {
    warning("Skipping re-localisation analysis: some subcellular fraction groups contain both reference and non-reference samples...")
  }
  ReportCalls <- AddSpace2Report()
  rm(list = ls()[which(!ls() %in% .obj)])
  Script <- readLines(ScriptPath)
  gc()
  invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
  saveImgFun(BckUpFl)
  #loadFun(BckUpFl)
  source(parSrc, local = FALSE)
}
