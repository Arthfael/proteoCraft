#' norm.mod
#'
#' @description 
#' Create a normal distribution model for normally distributed data.
#' The function assumes the data is normal, and will thus NOT check for normality by itself.
#' 
#' 
#' @param values The vector of normally distributed data
#' @param bins How many bins should we make? Default = 100.
#' @param adjust.bins TRUE by default. If set, will adjust bin number so there is at least a ratio of points to bins of 5, and at least 3 bins; failing that it will throw an error.
#' @param mu The mu parameter of the normal distribution. Usually I will be using this function on data where I want mu to be 0 (assumed if input mu is FALSE). If TRUE, the estimate will return a mu value.
#' @param use.frac What fraction of the data to use (from the center)? It will exclude outliers to avoid spreading the bins too far to both sides, so we get good sampling in the centre. Default is 0.9.
#' 
#' @examples
#' norm.mod(values = x)
#' 
#' @export

norm.mod <- function(values, bins = 100, mu = FALSE, adjust.bins = TRUE, use.frac = 0.9) {
  require(nls2)
  x1 <- proteoCraft::is.all.good(values)
  l <- length(x1)
  if (l) {
    x1 <- sort(x1)
    if (!mu) { d <- abs(x1) } else { d <- abs(x1 - median(x1)) }
    x1 <- x1[order(d)]
    x1 <- x1[1:round(l*use.frac)]
    l <- length(x1)
    if (adjust.bins) {
      if (floor(l/5) < 3) { stop("There is not enough data to work on!") }
      bins <- min(floor(l/5), bins)
    }
    if (!mu) {
      Max <- max(abs(x1))
      Min <- -Max   
    } else {
      Max <- max(x1)
      Min <- min(x1)
    }
    s <- (Max - Min)/bins
    binz <- data.frame(begin = Min + c(0:(bins-1))*s, end = Min + c(1:bins)*s)
    centrez <- rowMeans(binz)
    binz[1,1] <- binz[1,1]*1.1
    binz[bins,2] <- binz[bins,2]*1.1
    contentz <- apply(binz, 1, function(y) { length(which((x1 >= y[1])&(x1 <= y[2]))) })
    contentz <- contentz/sum(contentz)
    contentz <- data.frame(X = centrez, Y = contentz)
    if (!mu) {
      NLS <- nls2::nls2(Y ~ k*exp(-1/2*(X)^2/sigma^2),
                        algorithm = "brute-force",
                        start = c(sigma = sd(x1), k = 1),
                        data = contentz)
    } else {
      NLS <- nls2::nls2(Y ~ k*exp(-1/2*(X - mu)^2/sigma^2),
                        algorithm = "brute-force",
                        start = c(mu = median(x1), sigma = sd(x1), k = 1),
                        data = contentz)
    }
    NLS <- summary(NLS)$parameters[,"Estimate"]
    return(NLS) 
  } else { stop("There is no data to work on!") }
}
