#' FDR
#'
#' @description 
#' A function to calculate FDR thresholds and assess P-value significance.
#' 
#' @param data The dataframe that contains the values.
#' @param aggregate The aggregate of experimental factor levels (e.g. "Exp1___KO___Treated") to be pasted to the root of P-value column names to create the column name.
#' @param pvalue_root The root of the name of the P-values column. If the tag "-log10" is found in that name, P-values will assume to be -log10 transformed and will be reverse transformed accordingly.
#' @param pvalue_col Alternatively to the method of constructing a P-values column name used above, you may simply provide its name with this argument.
#' @param fdr Vector of acceptable False Discovery Rate values. Should be in the 0-1 range - but if any value is larger than 1 the function will automatically divide values by 100 (and throw an error if values are not now in the 0-1 range). Only required if returning a vector of significance or threshold values, see returns argument.
#' @param returns Logical vector of length 3. Default = c(TRUE, FALSE, FALSE).\cr
#'  - Set the first value to TRUE to output a vector of significance for each value in fdr (hence, this requires the fdr argument).\cr
#'  - Set the second to TRUE to output threshold values (NB: these are not-log transformed) for each value in fdr (hence, this requires the fdr argument).\cr
#'  - Set the third to TRUE to get the adjusted p-values.\cr
#' @param method Which method to use. Default = "BH" for Benjamini-Hochberg. Alternatively use "BY" for Benjamini-Yekutieli. 
#' @param SIMPLIFY Logical. Simplify the output from a list into a single object, if the output is a list of length 1? Default = TRUE although this will possibly move to FALSE in the future.
#'
#' @details
#' The default method is the Benjamini-Hochberg procedure.
#' Can also perform the Benjamini-Yekutieli procedure - when under positive dependence assumptions.
#' 
#' @returns
#' A list with 3 elements ("Significance vector", "Thresholds" and "Adj. P-values") or, if SIMPLIFY = TRUE and sum(returns) == 1, the single element of the list to return.
#' 
#' @examples
#' temp <- FDR(data = PG,
#'             aggregate,
#'             pvalue_root = "-log10 Pvalue.",
#'             fdr = BH.FDR)
#'
#' @export

FDR <- function(data,
                aggregate,
                pvalue_root,
                pvalue_col,
                fdr,
                returns = c(TRUE, FALSE, FALSE),
                method = "BH",
                SIMPLIFY = TRUE) {
  #proteoCraft::DefArg(proteoCraft::FDR)
  #
  # Check arguments
  # - returns
  stopifnot(length(returns) <= 3,
            "logical" %in% class(returns))
  if (length(returns) == 2) { returns <- c(returns, FALSE) }
  if (length(returns) == 1) { returns <- c(returns, FALSE, FALSE) }
  w <- which(is.na(returns))
  if (length(w)) {
    returns[w] <- c(TRUE, FALSE, FALSE)[w]
  }
  if (!sum(returns)) { warning("You haven't selected any type of output so nothing will be returned!") }
  # - method
  method <- gsub(" |-|_|\\.", "", toupper(method))
  stopifnot(method %in% c("BH", "BY"))
  # - fdr
  reqFDR <- sum(returns[1:2])
  if (reqFDR) {
    suppressWarnings(fdr <- as.numeric(fdr) )
    fdr <- fdr[which(!is.na(fdr))]
    lFDR <- length(fdr)
    stopifnot(lFDR > 0)
    stopifnot(sum(fdr <= 0) == 0)
    if (sum(fdr > 1)) { fdr <- fdr/100 }
    stopifnot(sum(fdr <= 0) == 0,
              sum(fdr > 1) == 0)
  }
  # - pvalue_col
  if (missing(pvalue_col)) {
    pvalue_col <- paste0(pvalue_root, aggregate)
  }
  # - SIMPLIFY
  if ((missing(SIMPLIFY))||(length(SIMPLIFY) != 1)||(!is.logical(SIMPLIFY))||(is.na(SIMPLIFY))) {
    SIMPLIFY <- TRUE
  }
  #
  P <- data.frame(Pvalues = data[[pvalue_col]])
  if (grepl("log10", pvalue_col)) { # If P-values are logged (auto-detect), then...
    P$Pvalues <- 10^(-P$Pvalues) #... de-log them.
  }
  P$Order <- 1:nrow(P)
  P$Rank <- NA
  P$All.Good <- proteoCraft::is.all.good(P$Pvalues, 2)
  w1 <- which(P$All.Good)
  N <- length(w1)
  if (N) {
    P$Rank[w1] <- rank(P$Pvalues[w1])
    P <- P[order(P$Rank, decreasing = FALSE),]
    w2 <- which(P$All.Good)
    rg2 <- 1:N
    if (method == "BH") { CN <- 1 }
    if (method == "BY") { CN <- sum(1/rg2) }
    BH.rank <- c() # Vector of FDR rank
    if (reqFDR) {
      thresh <- signifKols <- crtKols <- fdrKols <- c()
    }
    adjKol <- paste0("Adj. P-values")
    P[[adjKol]] <- NA
    P[w2, adjKol] <- P$Pvalues[w2]*N*CN/P$Rank[w2] # Basic calculation
    P[w2, adjKol] <- vapply(w2, function(x) { # (ensure monotonicity)
      min(c(1, P[w2[x:N], adjKol]))
    }, 1)
    #P$pAdjTst <- p.adjust(P$Pvalues, "BH") # You can verify that we get the same value using p.adjust()
    if (reqFDR) {
      for (j in 1:lFDR) { #j <- 3
        crtKols[j] <- crtKol <- paste0("Critical_", fdr[j])
        P[[crtKol]] <- NA
        P[w2, crtKol] <- P$Rank[w2]*fdr[j]/(N*CN)
        tt <- which(P$Pvalues[w2] <= P[w2, crtKol])
        fdrKols[j] <- fdrKol <- as.character(fdr[j])
        if (length(tt)) {
          M <- max(tt)
          thresh[j] <- P$Pvalues[w2][M]
          BH.rank[[fdrKol]] <- P$Rank[max(tt)]
        } else {
          thresh[j] <- 0
          BH.rank[[fdrKol]] <- 0
        }
        end <- "%"
        if (!missing(aggregate)) { end <- paste0(end, " - ", aggregate) }
        names(thresh)[j] <- paste0("Threshold-FDR=", fdr[j]*100, end)
        signifKols[j] <- signifKol <- paste0("Significant-FDR=", fdr[j]*100, end)
        P[[signifKol]] <- ""
        P[which(P$Rank <= BH.rank[[fdrKol]]), signifKol] <- "+"
      }
    }
    P <- P[order(P$Order, decreasing = FALSE),]
    res <- list()
    if (returns[1]) { res$"Significance vector" <- P[, signifKols, drop = FALSE] }
    if (returns[2]) { res$"Thresholds" <- thresh }
    if (returns[3]) { res$"Adj. P-values" <- P[, adjKol, drop = FALSE] }
    if ((SIMPLIFY)&&(length(res) == 1)) { res <- res[[1]] } # Simplify - probably a bad idea...
    return(res)
  } else {
    warning("This column contained no valid P-values to process!")
    return()
  }
}
