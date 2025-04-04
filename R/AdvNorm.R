#' AdvNorm
#'
#' @description
#' A function that performs advanced normalisation on dataframe containing a column of fractions names, peptide IDs, and 2 or more expression columns.  
#' 
#' @param df The input data frame, containing unique peptide identifiers, fraction numbers, and expression data.
#' @param pep.ids.col The name of the column in df containing peptide identifiers. Default is "Unique State", a column I usually create which contains modified sequence pasted to charge.
#' @param frac.col The name of the column in df which contains fraction numbers or raw file names.
#' @param exprs.col The name of the columns in df which contain expression values.
#' @param exprs.log Default = FALSE, meaning the input data is not log-transformed. Set instead to the base of the log if input data is log-transformed. If set to TRUE, assumes base 10. The data returned will be transformed or not as per the input.
#' 
#' @examples
#' adv.norm.data <- AdvNorm(df = data,
#'                          pep.ids.col = "Unique State",
#'                          frac.col = "Raw.file",
#'                          exprs.col = paste0("Reporter intensity ", c(0:9)),
#'                          exprs.log = FALSE)
#' require(ggplot2)
#' require(reshape2)
#' test <- melt(adv.norm.data)
#' test$value <- log10(test$value)
#' plot <- ggplot(test) + geom_density(stat = "density", aes(x = value, colour = variable))
#' grDevices::windows(width = 10, height = 10)
#' print(plot)
#' @export

AdvNorm <- function (df,
                     pep.ids.col = "Unique State",
                     frac.col,
                     exprs.col,
                     exprs.log = FALSE) {
  warning("Deprecation warning:\nThis function hasn't been used in, like, forever!\nCheck its code and update it (using the much newer AdvNorm.IL as inspiration) before using it! Also, add a pinch of data.tables!")
  if (exprs.log != FALSE) {
    if ((is.logical(exprs.log))&&(exprs.log)) { exprs.log <- 10 }
    for (i in exprs.col) { df[[i]] <- exprs.log^df[[i]] }
  }
  n.exprs <- length(exprs.col)
  arg.list <- list()
  for (i in 1:(length(exprs.col)-1)) { arg.list[[paste0("exprs.", i)]] <- 1 }
  for (i in 2:length(unique(df[[frac.col]]))) { arg.list[[paste0("frac.", i)]] <- 1 }
  names <- c(exprs.col, sort(unique(df[[frac.col]])))
  # Create function to calculate all residuals
  diff.log <- function(nms, cut, dat, p) {
    # There is also proteoCraft::diff.log! Can I replace one by the other?
    p <- unlist(p)
    stopifnot(length(p) + 2 == length(nms))
    x1 <- p[1:(cut-1)]
    x2 <- p[(cut):length(p)]
    x3 <- nms[2:cut]
    x4 <- nms[(cut+2):length(nms)]
    for (i in 1:length(x4)) {
      dat[which(dat[[frac.col]] == x4[i]), nms[1:cut]] <- dat[which(dat[[frac.col]] == x4[i]), nms[1:cut]]/x2[i]
    }
    dat[, x3] <- sweep(dat[, x3], 2, x1, "/")
    dat2 <- aggregate(dat[, nms[1:cut]],
                      list(dat[[pep.ids.col]]), sum)
    dat3 <- as.matrix(dat2[, nms[1:cut]])
    dat3 <- as.data.frame(t(apply(dat3, 1, function(z) {
      z1 <- c(1:(length(z)-1))
      z1 <- sapply(z1, function(w) {
        w <- proteoCraft::is.all.good(log10(z[w]) - log10(z[(w+1):length(z)]))
      })
      return(proteoCraft::is.all.good(unlist(z1)))
    })))
    dat3 <- unlist(dat3)
    return(dat3)
  }
  diff.log.2 <- function(...) {
    p <- list(...)
    res <- diff.log(nms = names, cut = n.exprs, dat = df, p)
    return(res)
  }
  # Perform optimisation of normalisation parameters:
  LM <- minpack.lm::nls.lm(par = arg.list,
                           fn = diff.log.2)
  # Apply to expression data
  res <- df[, exprs.col]
  a <- names[(n.exprs+2):length(names)]
  b <- as.numeric(LM$par)
  for (i in 1:length(a)) {
    w <- which(res[[frac.col]] == a[i])
    res[w,] <- res[w,]/b[n.exprs:length(b)][i]
  }
  res[, 2:ncol(res)] <- sweep(res[, 2:ncol(res)], 2, b[1:(n.exprs-1)], "/")
  if (exprs.log != FALSE) {
    for (i in exprs.col) { res[[i]] <- log(res[[i]], exprs.log) }
  }
  colnames(res) <- paste0("AdvNorm.", colnames(res))
  return(res)
}
