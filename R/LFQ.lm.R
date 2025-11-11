#' LFQ.lm
#'
#' @description 
#' A function which summarizes individual peptide/fragment quantitative values into a single quantitative vector, using the Levenberg-Marquardt procedure.
#' Used at two levels:
#' - In Prot_Quant, to summarize peptidoform-level information into a protein group-levels quantitative vector.
#' - In the rarely used MS2corr2MS1.R source, to summarize a peptidoform's MS1 (precursor) and MS2 (fragment) measurements into a single quantitative vector.
#' The output is log10-transformed.
#' 
#' @param ids Integer, a vector/list of ids in InputTabl to summarize.
#' @param InputTabl Table containing individual entities (peptides or fragments) to summarize.
#' @param id Name of the IDs column in InputTabl. Default = "id"
#' @param IntensCol Names of the quantitative columns in InputTabl
#' @param Summary.method The summary method used for ratios (the Levenberg-Marquardt algorithm is used for intensities). One of "mean", "median", or "weighted.mean".
#' @param Summary.weights If a "weighted.mean" summary method is chosen, then a vector of weights must be provided. This should be the name of a column of the peptides table
#' @param Intensity.weights default = FALSE; if TRUE, will take into account individual peptide intensities when calculating average profile if "Summary.method" is either mean "mean" or "weighted.mean" (thus, it will actually be a weighted mean regardless). 
#' @param Min.N How many peptides should at least be present? Should be at the very least 1 (default = 2).
#' @param Max.N How many peptides can we use at most for the Levenberg-Marquardt procedure? Using too many peptides can be an issue, e.g. with huge proteins like Titin. Default = 50. The most intense peptides will be selected.
#' @param Is.log Are input intensities log-transformed? TRUE by default.
#' @param reNorm Integer. If set to:\cr
#' - 1 (default), will re-normalize the profile's median intensity to that of the highest median intensity peptide. 
#' - 2 (former default), will re-normalize the profile to the highest original single intensity value. 
#' - 0, will re-normalize the profile to its own median intensity.
#' 
#' @examples
#' lfq <- LFQ.lm(ids,
#'               InputTabl = Pep,
#'               IntensCol = Pep.Intens.Nms,
#'               Summary.method = Summary.method,
#'               Summary.weights = Summary.weights,
#'               Min.N = Min.N,
#'               Max.N = Max.N,
#'               Is.log = TRUE)
#' # Used within Prot.Quant
#' 
#' @export

LFQ.lm <- function(ids,
                   InputTabl,
                   id = "id",
                   IntensCol,
                   Summary.method = "median",
                   Summary.weights,
                   Min.N = 2,
                   Max.N = 50,
                   Is.log = TRUE,
                   reNorm = 1) {
  Viz = FALSE
  if ((missing(reNorm))||(!is.numeric(reNorm))||(length(reNorm) != 1)||(is.na(reNorm))||(!reNorm %in% 0:2)) {
    reNorm <- 1
  }
  #proteoCraft::DefArg(proteoCraft::LFQ.lm)
  #ids = IDsInputTabl = MSAll;IntensCol = Samples;Summary.method = "median";Min.N = 2;Max.N = 50;Is.log = FALSE
  #ids <- temp.ids;InputTabl = tmpPep;IntensCol = paste0(Pep.Intens.root, Samples);Summary.method = "median";Min.N = 2;Max.N = 50;Is.log = TRUE
  #ids <- temp.ids[prot.list[1]]; Viz <- TRUE
  summaryFun <- get(Summary.method)
  rs <- setNames(rep(NA, length(IntensCol)), IntensCol) # This is the default vector of NAs to replace if the data is complete enough
  # Check that we have at least enough peptidoforms:
  # NB: Currently done at peptidoforms level, i.e. accepts two peptidoforms of the same primary sequence as different.
  # If I ever decide to change it, this would be the place.
  if (missing("Summary.weights")) {
    Summary.weights <- "Summary.weights"
    InputTabl[[Summary.weights]] <- 1
  }
  mtch <- match(unlist(ids), InputTabl[[id]])
  if (length(mtch) >= Min.N) { # Are there enough values?
    temp1 <- InputTabl[mtch, IntensCol, drop = FALSE]
    if (!Is.log) { temp1 <- log10(temp1) } # The rest assumes log-transformed data
    wights <- InputTabl[mtch, Summary.weights]
    w1 <- which(apply(temp1[, IntensCol, drop = FALSE], 1, function(y) {
      length(proteoCraft::is.all.good(y))
    }) > 0)
    if (length(w1)) { # Are there enough valid values?
      # Remove peptides with only non-valid or missing values
      temp1 <- temp1[w1, , drop = FALSE]; wights <- wights[w1]
      tst2 <- vapply(IntensCol, function(y) {
        length(proteoCraft::is.all.good(temp1[[y]]))
      }, 1)
      if (max(tst2) > 1) { # Are there columns with at least 2 valid values?
        if (length(w1) > Max.N) { # (No need to re-calculate tst1)
          # Remove excess peptides, for cases where we have much more than we need and including all would slow down the calculations
          # (I'm looking at you, Titin!!!)
          Ranks <- nrow(temp1) + 1 - rank(apply(temp1[, IntensCol, drop = FALSE], 1, function(x) {
            sum(proteoCraft::is.all.good(x))
          }))
          wR <- which(Ranks <= Max.N)
          temp1 <- temp1[wR, IntensCol, drop = FALSE]
          wights <- wights[wR]
        }
        # Columns with at least 1 valid value
        wNN <- which(vapply(IntensCol, function(y) {
          length(proteoCraft::is.all.good(temp1[[y]]))
        }, 1) > 0)
        av <- apply(temp1[, wNN, drop = FALSE], 1, function(y) {
          median(proteoCraft::is.all.good(y))
        })
        temp2 <- sweep(temp1[, wNN, drop = FALSE], 1, av, "-") # Normalized profiles row-wise to the median
        f <- rep(0, nrow(temp2) - 1)
        diff.log.v <- function(...) {
          p <- list(...)
          res <- proteoCraft::diff.log(p, dat = temp2)
          return(res)
        }
        LM <- minpack.lm::nls.lm(par = f,
                                 fn = diff.log.v,
                                 lower = unlist(f)-1,
                                 upper = unlist(f)+1) # Align
        temp3 <- sweep(temp2, 1, c(0, LM$par), "-") # Fine LM row-wise normalization
        if (Viz) {
          tmp1 <- as.matrix(temp1)
          tmp1[which(!is.finite(tmp1), arr.ind = TRUE)] <- NA
          tmp2 <- as.matrix(temp2)
          tmp2[which(!is.finite(tmp2), arr.ind = TRUE)] <- NA
          tmp3m <- as.matrix(temp3)
          tmp3m[which(!is.finite(tmp3m), arr.ind = TRUE)] <- NA
          grDevices::windows(width = 10, height = 10)
          par(cex.main = 0.3)
          gplots::heatmap.2(tmp1, Colv = NULL, Rowv = NULL,
                            main = "Original", xlab = NULL, ylab = NULL,
                            key = TRUE, keysize = 1,
                            trace = "none", density.info = c("none"),
                            na.color = "black",
                            sepcolor = "blue", dendrogram = "none",
                            cexRow = 0.7,
                            cexCol = 0.7)
          grDevices::windows(width = 10, height = 10)
          par(cex.main = 0.3)
          gplots::heatmap.2(tmp2, Colv = NULL, Rowv = NULL,
                            main = "Row normalized", xlab = NULL, ylab = NULL,
                            key = TRUE, keysize = 1,
                            trace = "none", density.info = c("none"),
                            na.color = "black",
                            sepcolor = "blue", dendrogram = "none",
                            cexRow = 0.7,
                            cexCol = 0.7)
          grDevices::windows(width = 10, height = 10)
          par(cex.main = 0.3)
          gplots::heatmap.2(tmp3m, Colv = NULL, Rowv = NULL,
                            main = "Aligned", xlab = NULL, ylab = NULL,
                            key = TRUE, keysize = 1,
                            trace = "none", density.info = c("none"),
                            na.color = "black",
                            sepcolor = "blue", dendrogram = "none",
                            cexRow = 0.7,
                            cexCol = 0.7)
          #dev.off()
        }
        # Average aligned profiles
        temp3 <- apply(temp3, 2, function(y) {
          wAG <- which(proteoCraft::is.all.good(y, 2))
          y <- y[wAG]
          wghts <- wights[wAG] # These are all 1 if method is not weighted mean
          if (length(y)) {
            if (Summary.method == "weighted.mean") { y <- summaryFun(y, wghts) } else { y <- summaryFun(y) }
          } else {
            y <- NA # In previous versions was 0 here
            # However this cannot be justified:
            # If there is no valid observation then there is no reason to have "no change" when in reality what we have is "nothing"
          }
          return(y)
        })
        # Apply best-flyer hypothesis logic for estimating absolute quant level
        if (reNorm %in% 0:1) {
          m0 <- median(temp3[wNN])
          temp3 <- temp3 - m0
          if (reNorm == 1) {
            m1 <- median(proteoCraft::is.all.good(apply(temp1, 1, function(x) { median(proteoCraft::is.all.good(x)) })))
            temp3 <- temp3 + m1
          }
          rs[wNN] <- temp3
        }
        if (reNorm == 2) {
          m <- max(proteoCraft::is.all.good(unlist(temp1[, wNN])))
          m <- as.data.frame(which(temp1[, wNN, drop = FALSE] == m, arr.ind = TRUE))
          rs[wNN] <- as.numeric(temp3 + temp1[m$row[1], wNN[m$col[1]]] -  temp3[m$col[1]])
        }
       
      } else {
        rs <- setNames(apply(temp1, 2, function(y) {
          y <- proteoCraft::is.all.good(y)
          if (!length(y)) { y <- NA }
          return(y)
        }), IntensCol)
      }
      if (!Is.log) { rs <- 10^rs }
    }
  }
  return(rs)
}
