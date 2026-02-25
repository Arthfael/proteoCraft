# Isobaric label purity correction
if ((LabelType == "Isobaric")&&
    ("Label.Purities.file" %in% colnames(Param))&&
    (!Param$Label.Purities.file %in% c("", " ", "NA", NA))&&
    (file.exists(Param$Label.Purities.file))) {
  if (!require("matlib", quietly = TRUE)) { install.packages("matlib") }
  cran_req <- c(cran_req, "matlib")
  Iso.purity <- read.csv(Param$Label.Purities.file)
  msg <- paste0("Correcting ", evNm, " reporter intensities for labels purity.")
  ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE)
  DatAnalysisTxt <- paste0(DatAnalysisTxt, " Reporter intensities were corrected for ", IsobarLab, " label purity factors.")
  if (colnames(Iso.purity)[1] == "New.format") {
    Iso.purity <- read.csv(Param$Label.Purities.file, check.names = FALSE)
    Iso.purity <- Iso.purity[, 2:ncol(Iso.purity)]
    Iso.purity$Isobaric.set <- lapply(strsplit(as.character(Iso.purity$Isobaric.set), ";"), as.integer)
    test <- sort(unique(unlist(Iso.purity$Isobaric.set)))
    if (exists("Iso")) {
      if (!sum(!Iso %in% test)) {
        test <- data.frame(Isobaric.set = test, check.names = FALSE)
        test$"Purity table" <- lapply(test$Isobaric.set, function(x) {
          Iso.purity[which(vapply(Iso.purity$Isobaric.set, function(y) { x %in% y }, TRUE)),
                     which(colnames(Iso.purity) != Isobaric.set)]
        })
        kol <- paste0("corr. ", ev.ref["Original"], get(IsobarLab))
        kol <- kol[which(paste0(ev.ref["Original"], get(IsobarLab)) %in% colnames(ev))]
        ev[, kol] <- NA
        for (i in Iso) { #i <- Iso[1]
          wI <- which(ev$"Raw file path" %in% Frac.map$"Raw file"[which(Frac.map$Isobaric.set == i)])
          e <- ev[wI,]
          pt <- test$"Purity table"[[which(test$Isobaric.set == i)]]
          lb <- Exp.map$"Isobaric label"[which(Exp.map$Isobaric.set == i)]
          lb2 <- Exp.map$"Isobaric label details"[which(Exp.map$Isobaric.set == i)]
          kol <- paste0(ev.ref[length(ev.ref)], lb)
          w <- which(kol %in% colnames(e))
          kol <- kol[w] ; lb <- lb[w] ; lb2 <- lb2[w]
          o <- order(as.numeric(lb))
          kol <- kol[o] ; lb <- lb[o] ; lb2 <- lb2[o]
          w <- match(lb2, pt$"Isobaric label details")
          pt <- pt[w,]
          w <- which(colnames(pt) != "Isobaric label details")
          A <- matrix(rep(0, length(lb)*(length(lb)+length(w)-1)), ncol = length(lb))
          for (j in 1:length(lb2)) {
            tmp <- as.numeric(pt[j, w])
            tmp[which(is.na(tmp))] <- 0
            A[j:(j+length(w)-1), j] <- tmp
          }
          w <- which(vapply(1:length(w), function(x) {
            sum(vapply(1:length(lb), function(y) { A[y+x-1, y] == 100 }, TRUE))
          }, 1) == length(lb))
          A <- A[(1:length(lb))+w-1,]
          source(parSrc, local = FALSE)
          exports <- list("A", "e", "kol")
          clusterExport(parClust, exports, envir = environment())
          clusterCall(parClust, function() library(matlib))
          clusterCall(parClust, function() library(proteoCraft))
          temp <- as.data.frame(t(parApply(parClust, e[,kol], 1, function(x) {
            b <- as.numeric(x)
            b[which(!is.all.good(b, 2))] <- 0
            sb <- sum(b)
            if (sb > 0) {
              #showEqn(round(A, 3), b)
              #res <- as.numeric(gsub(".+= +", "", Solve(A, b))) # I stopped using this since it seems to return approximations sometimes
              res <- solve(A, b)
              # There are cases where we will get negative values which we will have to truncate:
              res[which(res < 0)] <- 0
              res <- #round(
                res*sb/sum(res)
              #, 0)# I used to round here, but do not think it's necessary
              # print(res/b)
            } else { res <- rep(NA, length(kol)) }
            res <- setNames(res, kol)
            return(res)
          })))
          e[, gsub(topattern(ev.ref[length(ev.ref)]), paste0("corr. ", ev.ref["Original"]), kol)] <- temp[, kol]
          ev[wI,] <- e
        }
        ev.ref["Corrected"] <- paste0("corr. ", ev.ref["Original"])
        cat("Reporter Intensity correction done!\n")
      } else { warning("Some Isobaric sets are not defined in the labels purity table, skipping!") }
    } else { stop("There should be an \"Iso\" (Isobaric sets) object by now!") }
  } else {
    Iso.purity$Isobaric.set <- as.character(Iso.purity$Isobaric.set)
    test <- sort(unique(as.integer(unlist(strsplit(Iso.purity$Isobaric.set, ";")))))
    if (exists("Iso")) {
      if (!sum(!Iso %in% test)) {
        test <- test[which(test %in% Iso)]
        test <- data.frame(Isobaric_set = test)
        test$Purity.table <- lapply(test$Isobaric_set, function(x) {
          Iso.purity[which(vapply(lapply(strsplit(Iso.purity$Isobaric.set, ";"), as.integer), function(y) { x %in% unlist(y) }, TRUE)),
                     which(colnames(Iso.purity) != "Isobaric.set")]
        })
        kol <- paste0("corr. ", ev.ref["Original"], get(IsobarLab))
        kol <- kol[which(paste0(ev.ref["Original"], get(IsobarLab)) %in% colnames(ev))]
        ev[, kol] <- NA
        for (i in Iso) { #i <- Iso[1]
          wI <- which(ev$"Raw file path" %in% Frac.map$"Raw file"[which(Frac.map$Isobaric.set == i)])
          e <- ev[wI,]
          pt <- test$Purity.table[[which(test$Isobaric_set == i)]]
          lb <- Exp.map$"Isobaric label"[which(Exp.map$Isobaric.set == i)]
          lb2 <- Exp.map$"Isobaric label details"[which(Exp.map$Isobaric.set == i)]
          kol <- paste0(ev.ref[length(ev.ref)], lb)
          w <- which(kol %in% colnames(e))
          kol <- kol[w] ; lb <- lb[w] ; lb2 <- lb2[w]
          o <- order(as.numeric(lb))
          kol <- kol[o] ; lb <- lb[o] ; lb2 <- lb2[o]
          pt <- pt[which(pt$Isobaric.label.details %in% lb2),]
          kol2 <- c("MI_minus.2", "MI_minus.1", "Monoisotopic", "MI_plus.1", "MI_plus.2")
          kol3 <- c("Who_minus.2", "Who_minus.1", "Who_plus.1", "Who_plus.2")
          pt[, kol2] <- sweep(pt[, kol2], 1, rowSums(pt[, kol2]), "/")
          A <- as.data.frame(matrix(rep(0, length(lb)*(length(lb)+4)), nrow = length(lb)))
          colnames(A) <- c("-2", "-1", lb2, "+1", "+2")
          kount <- 0
          for (l in lb2) { #l <- lb2[2]
            kount <- kount+1
            tmp <- unlist(pt[which(pt$Isobaric.label.details == l), kol3])
            wc <- wc2 <- which(tmp != "")
            wc2[which(wc2 > 2)] <- wc2[which(wc2 > 2)]+1
            fact <- pt[which(pt$Isobaric.label.details == l), kol2, drop = FALSE]
            fact <- fact[, sort(c(3, wc2))]
            colnames(fact) <- c(l, tmp[wc])[order(c(3, wc2))]
            fact <- fact[, which(colnames(fact) %in% colnames(A))]
            if (length(fact)) {
              A[kount, colnames(fact)] <- fact
            }
          }
          A <- as.matrix(A[,3:(length(lb)+2)])
          exports <- list("A", "e", "kol")
          clusterExport(parClust, exports, envir = environment())
          clusterCall(parClust, function() library(matlib))
          clusterCall(parClust, function() library(proteoCraft))
          temp <- as.data.frame(t(parApply(parClust, e[, kol], 1, function(x) {
            b <- as.numeric(x)
            b[which(!is.all.good(b, 2))] <- 0
            sb <- sum(b)
            if (sb > 0) {
              #showEqn(round(A, 3), b)
              #res <- as.numeric(gsub(".+= +", "", Solve(A, b))) # I stopped using this since it seems to return approximations sometimes
              res <- solve(A, b)
              # There are cases where we will get negative values which we will have to truncate:
              res[which(res < 0)] <- 0
              res <- #round(
                res*sb/sum(res)
              #, 0)# I used to round here, but do not think it's necessary
              # print(res/b)
            } else { res <- rep(NA, length(kol)) }
            return(res)
          })))
          e[, gsub(topattern(ev.ref[length(ev.ref)]), paste0("corr. ", ev.ref["Original"]), kol)] <- temp
          ev[wI,] <- e
        }
        ev.ref["Corrected"] <- paste0("corr. ", ev.ref["Original"])
        cat("Reporter Intensity correction done!\n")
      } else { warning("Some Isobaric sets are not defined in the labels purity table, skipping!") }
    } else { stop("There should be an \"Iso\" (Isobaric sets) object by now!") }
  }
}
