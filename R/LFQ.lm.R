#' LFQ.lm
#'
#' @description 
#' A function which summarizes individual peptide/fragment quantitative values into a single quantitative vector, using the Levenberg-Marquardt procedure.
#' Used at two levels:
#' - For DIA-NN input, to summarize a peptidoform's MS1 (precursor) and MS2 (fragment) measurements into a single quantitative vector.
#' - In Prot.Quant, to summarize peptidoform-level information into a protein group-levels quantitative vector.
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
#' @param Viz Should we plot the results? For testing only, default = FALSE.
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
                Viz = FALSE,
                Is.log = TRUE) {
  #proteoCraft::DefArg(proteoCraft::LFQ.lm)
  #ids = IDsInputTabl = MSAll;IntensCol = Samples;Summary.method = "median";Min.N = 2;Max.N = 50;Is.log = FALSE
  #ids <- temp.ids;InputTabl = tmpPep;IntensCol = paste0(Pep.Intens.root, Samples);Summary.method = "median";Min.N = 2;Max.N = 50;Is.log = TRUE
  #ids <- temp.ids[prot.list[1]]; Viz <- TRUE
  sum.func <- get(Summary.method)
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
    temp2 <- InputTabl[mtch, IntensCol, drop = FALSE]
    if (!Is.log) { temp2 <- log10(temp2) } # The rest assumes log-transformed data
    wights <- InputTabl[mtch, Summary.weights]
    w1 <- which(apply(temp2[, IntensCol, drop = FALSE], 1, function(y) {
      length(proteoCraft::is.all.good(y))
    }) > 0)
    if (length(w1)) { # Are there enough valid values?
      # Remove peptides with only non-valid or missing values
      temp2 <- temp2[w1, , drop = FALSE]; wights <- wights[w1]
      tst2 <- sapply(IntensCol, function(y) {
        length(proteoCraft::is.all.good(temp2[[y]] ))
      })
      if (max(tst2) > 1) { # Are there columns with at least 2 valid values?
        if (length(w1) > Max.N) { # (No need to re-calculate tst1)
          # Remove excess peptides, for cases where we have much more than we need and including all would slow down the calculations
          # (I'm looking at you, Titin!!!)
          Ranks <- nrow(temp2) + 1 - rank(apply(temp2[, IntensCol, drop = FALSE], 1, function(x) {
            sum(proteoCraft::is.all.good(x))
          }))
          wR <- which(Ranks <= Max.N)
          temp2 <- temp2[wR, IntensCol, drop = FALSE]
          wights <- wights[wR]
        }
        # Columns with at least 1 valid value
        wNN <- which(sapply(IntensCol, function(y) {
          length(proteoCraft::is.all.good(temp2[[y]]))
        }) > 0)
        av <- apply(temp2[, wNN, drop = FALSE], 1, function(y) {
          median(proteoCraft::is.all.good(y))
        })
        temp3 <- sweep(temp2[, wNN, drop = FALSE], 1, av, "-") # Normalized profiles row-wise to the median
        f <- rep(0, nrow(temp3) - 1)
        diff.log.v <- function(...) {
          p <- list(...)
          res <- proteoCraft::diff.log(p, dat = temp3)
          return(res)
        }
        LM <- minpack.lm::nls.lm(par = f,
                                 fn = diff.log.v,
                                 lower = unlist(f)-1,
                                 upper = unlist(f)+1) # Align
        temp4 <- sweep(temp3, 1, c(0, LM$par), "-") # Fine LM row-wise normalization
        if (Viz) {
          tmp2 <- as.matrix(temp2)
          tmp2[which(!is.finite(tmp2), arr.ind = TRUE)] <- NA
          tmp3 <- as.matrix(temp3)
          tmp3[which(!is.finite(tmp3), arr.ind = TRUE)] <- NA
          tmp4m <- as.matrix(temp4)
          tmp4m[which(!is.finite(tmp4m), arr.ind = TRUE)] <- NA
          grDevices::windows(width = 10, height = 10)
          par(cex.main = 0.3)
          gplots::heatmap.2(tmp2, Colv = NULL, Rowv = NULL,
                            main = "Original", xlab = NULL, ylab = NULL,
                            key = TRUE, keysize = 1,
                            trace = "none", density.info = c("none"),
                            na.color = "black",
                            sepcolor = "blue", dendrogram = "none",
                            cexRow = 0.7,
                            cexCol = 0.7)
          grDevices::windows(width = 10, height = 10)
          par(cex.main = 0.3)
          gplots::heatmap.2(tmp3, Colv = NULL, Rowv = NULL,
                            main = "Row normalized", xlab = NULL, ylab = NULL,
                            key = TRUE, keysize = 1,
                            trace = "none", density.info = c("none"),
                            na.color = "black",
                            sepcolor = "blue", dendrogram = "none",
                            cexRow = 0.7,
                            cexCol = 0.7)
          grDevices::windows(width = 10, height = 10)
          par(cex.main = 0.3)
          gplots::heatmap.2(tmp4m, Colv = NULL, Rowv = NULL,
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
        temp4 <- apply(temp4, 2, function(y) {
          wAG <- which(proteoCraft::is.all.good(y, 2))
          y <- y[wAG]
          wghts <- wights[wAG] # These are all 1 if method is not weighted mean
          if (length(y)) {
            if (Summary.method == "weighted.mean") { y <- sum.func(y, wghts) } else { y <- sum.func(y) }
          } else {
            y <- NA # In previous versions was 0 here
            # However this cannot be justified:
            # If there is no valid observation then there is no reason to have "no change" when in reality what we have is "nothing"
          }
          return(y)
        })
        # Apply best-flyer hypothesis logic for estimating absolute quant level
        m <- max(proteoCraft::is.all.good(unlist(temp2[,wNN])))
        m <- as.data.frame(which(temp2[, wNN, drop = FALSE] == m, arr.ind = TRUE))
        rs[wNN] <- as.numeric(temp4 + temp2[m$row[1], wNN[m$col[1]]] -  temp4[m$col[1]])
      } else {
        rs <- setNames(apply(temp2, 2, function(y) {
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
