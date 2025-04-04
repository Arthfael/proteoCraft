#' FDR
#'
#' @description 
#' A function to calculate FDR thresholds and assess P-value significance.
#' Default is the Benjamini-Hochberg procedure.
#' Can also do the Benjamini-Yekutieli procedure under positive dependence assumptions.
#' 
#' @param data The dataframe that contains the values.
#' @param aggregate The current aggregate to be pasted to the root of P-value column names to create the column name.
#' @param pvalue_root The root of the name of the P-values column. Default = "Ratio_Minus.log10.Pvalue." If the tag "log10" is found in that name, P-values will assume to be -log10 transformed and will be reverse transformed accordingly.
#' @param pvalue_col Alternatively to the method of constructing a P-values column name used above, you may simply provide its name with this argument.
#' @param fdr The acceptable False Discovery Rate values.
#' @param returns Default = c(TRUE, FALSE) Set first to TRUE to output vector of significance; set second to TRUE to output threshold values (NB: these are not-log transformed).
#' @param method The method to use. Default = "BH" for Benjamini-Hochberg. Alternatively use "BY" for Benjamini-Yekutieli. 
#' 
#' @examples
#' temp <- FDR(data = PG,
#'             aggregate,
#'             pvalue_root = "Ratio_Minus.log10.Pvalue.",
#'             fdr = BH.FDR)
#'
#' @export

FDR <- function(data,
                aggregate,
                pvalue_root,
                pvalue_col,
                fdr,
                returns = c(TRUE, FALSE),
                method = "BH") {
  #proteoCraft::DefArg(proteoCraft::FDR)
  stopifnot(length(returns) == 2, class(returns) == "logical")
  method <- gsub(" |-|_|\\.", "", toupper(method))
  stopifnot(method %in% c("BH", "BY"))
  if (sum(returns) == 0) { warning("You haven't selected any type of output so nothing will be returned!") }
  if (sum(fdr > 1) > 0) { fdr <- fdr/100 }
  stopifnot(sum(fdr <= 0) == 0)
  if (missing(pvalue_col)) { 
    pvalue_col <- paste0(pvalue_root, aggregate)
  }
  P <- data.frame(Pvalues = data[[pvalue_col]])
  if (grepl("log10", pvalue_col)) { # If P-values are logged (auto-detect), then...
    P$Pvalues <- 10^(-P$Pvalues) #... de-log them.
  }
  P$Order <- c(1:nrow(P))
  P$Rank <- NA
  P$All.Good <- proteoCraft::is.all.good(P$Pvalues, 2)
  w1 <- which(P$All.Good)
  P$Rank[w1] <- rank(P$Pvalues[w1])
  P <- P[order(P$Rank),]
  w2 <- which(P$All.Good)
  N <- length(w2)
  if (N > 0) {
    BH.rank <- c() # Vector of FDR rank
    k <- c()
    thresh <- c()
    for (j in 1:length(fdr)) { #j <- 3
      P[[paste0("Critical_", fdr[j])]] <- NA
      if (method == "BH") { CN <- 1 }
      if (method == "BY") { CN <- sum(1/c(1:N)) }
      P[w2, paste0("Critical_", fdr[j])] <- P$Rank[w2]*fdr[j]/(N*CN)
      tt <- which(P$Pvalues[w2] <= P[w2, paste0("Critical_", fdr[j])])
      if (length(tt) > 0) {
        M <- max(tt)
        thresh[j] <- P$Pvalues[w2][M]
        BH.rank[[as.character(fdr[j])]] <- P$Rank[max(tt)]
      } else {
        thresh[j] <- 0
        BH.rank[[as.character(fdr[j])]] <- 0
      }
      end <- "%"
      if (!missing(aggregate)) { end <- paste0(end, " - ", aggregate) }
      names(thresh)[j] <- paste0("Threshold-FDR=", fdr[j]*100, end)
      k[j] <- paste0("Significant-FDR=", fdr[j]*100, end)
      P[[k[j]]] <- ""
      P[which(P$Rank <= BH.rank[[as.character(fdr[j])]]), k[j]] <- "+"
    }
    P <- P[order(P$Order),]
    res <- list()
    if (returns[1]) { res$"Significance vector" <- P[, k, drop = FALSE] }
    if (returns[2]) { res$"Thresholds" <- thresh }
    if (length(res) == 1) { res <- res[[1]] }
    return(res)
  } else { return("This column contained no valid P-values to process!") }
}
