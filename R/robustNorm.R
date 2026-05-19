#' robustNorm
#' 
#' @description
#' This is meant as an alternative to AdvNorm.IL() for applying a more robust method to aligning the samples.
#' I have only tested, not written this function. ChatGPT and Perplexity did, in a merry-go-round, where I bounced the code from one to the other. Yes. We are there.
#' The idea for a more robust loss function was floated to my by chatGPT as  was discussing normalization with 'it'.
#' And because I was too busy with other stuff and unfamiliar with the robuts loss function landscape, I asked it to make it for me.
#' 
#' 
#' @param df Data frame containing expression values.
#' @param myColumns Name of the columns containing expression data. THIS SHOULD BE LOG-TRANSFORMED!
#' @param loss Loss function. Options are "l2" (= sum of squared errors), "huber", "tukey" or "cauchy".
#' @param k Tuning constants for the method, allowing to override the defaults.
#' @param max_iter Integer, default = 50L.
#' @param tol Tolerance, default = 1e-6
#' 
#' @returns
#' A list:
#'  - data = data.frame with normalized values
#'  - sample_offset = optimized offsets
#'  - row_mean = row_mean,
#'  - weights
#'  - iterations = last iteration before convergence 
#'  - loss = loss function used
#'  - sigma = scaled Median Absolute Deviation of the residuals
#' 
#' @examples
#' normData <- robustNorm(pep, paste0("log10(Int.) - ", mySamples))
#' 
#' @export

robustNorm <- function(df,
                       myColumns,
                       loss = c("l2", "huber", "tukey", "cauchy"),
                       k = NULL,
                       max_iter = 50L,
                       tol = 1e-6) {
  DefArg(robustNorm)
  loss <- match.arg(loss)
  #
  X <- as.matrix(df[, myColumns])
  p <- nrow(X)
  s <- ncol(X)
  # Tuning constants
  if (is.null(k)) {
    k <- switch(loss,
                l2 = Inf,
                huber = 1.345,
                tukey = 4.685,
                cauchy = 2.385)
  }
  # Initialize
  sample_offset <- rep(0, s)
  # Initial protein means
  row_mean <- rowMeans(X, na.rm = TRUE)
  # Helper: weighted mean
  wmean <- \(x, w) {
    ok <- is.finite(x) & is.finite(w) & (w > 0)
    if (!any(ok)) { return(NA_real_) }
    return(sum(x[ok]*w[ok])/sum(w[ok]))
  }
  # Helper: robust weights
  robust_weights <- \(r, loss, k) {
    u <- r / k
    w <- switch(loss,
                l2 = rep(1, length(r)),
                huber = ifelse(abs(u) <= 1, 1, 1/abs(u)),
                tukey = {
                  w <- (1-u^2L)^2L
                  w[abs(u) >= 1] <- 0
                  w
                },
                cauchy = 1/(1+u^2L))
    w[which(!is.finite(w))] <- 0
    return(w)
  }
  # Main IRLS loop
  for (iter in 1L:max_iter) {
    sample_offset_old <- sample_offset
    # Current fitted values
    fitted <- outer(row_mean, sample_offset, "+")
    residuals <- X - fitted
    # Robust scale estimate (MAD)
    sigma <- mad(as.vector(residuals), na.rm = TRUE)
    if ((!is.finite(sigma)) || (sigma == 0L)) { sigma <- 1 }
    # Standardized residuals
    r_std <- residuals/sigma
    # Robust weights
    W <- robust_weights(r_std, loss, k)
    W[which(is.na(X))] <- NA
    # Update protein means
    X_adj <- sweep(X, 2L, sample_offset, "-")
    row_mean <- vapply(seq_len(p), \(i) { wmean(X_adj[i,], W[i,]) }, 1)
    # Update sample offsets
    X_centered <- sweep(X, 1L, row_mean, "-")
    sample_offset <- vapply(seq_len(s), \(j) { wmean(X_centered[, j], W[, j]) }, 1)
    # Identifiability constraint
    sample_offset <- sample_offset - mean(sample_offset, na.rm = TRUE)
    # Convergence
    delta <- max(abs(sample_offset - sample_offset_old), na.rm = TRUE)
    if (delta < tol) {
      message("Converged at iteration ", iter)
      break
    }
  }
  # Normalized matrix
  X_norm <- sweep(X, 2L, sample_offset, "-")
  df[, myColumns] <- X_norm
  return(list(data = df,
              sample_offset = sample_offset,
              row_mean = row_mean,
              weights = W,
              iterations = iter,
              loss = loss,
              sigma = sigma))
}
