#' norm.thresh
#'
#' @description 
#' Find thresholds to include only a given proportion of data.
#' 
#' @param values The vector of normally distributed values for which we want a threshold.
#' @param ref The reference. If left as NULL, proportions are from the values vector. If provided as a vector of values, the function aims to find thresholds including a given proportion (based on modelled normal distributions) p of reference to values.
#' @param p The relevant proportion to include (either absolute, or contamination by a reference).
#' @param shapiro.threshold Which threshold can we accept for the Shapiro-Wilk normality test's P-value? Default = 0.05
#' 
#' @examples
#' norm.thresh(values = test$"log2.Ratio.Exp1 Astrocyte ApoE4/4",
#'               ref = test$"log2.Ratio.Exp1 Astrocyte ApoE3/3")
#' 
#' @export

norm.thresh <- function(values, ref = NULL, p = 0.01, shapiro.threshold = 0.05) {
  stopifnot(p > 0, p < 1)
  v1 <- proteoCraft::is.all.good(values)
  modV <- proteoCraft::norm.mod2(values = v1)
  if (!"Mean" %in% names(modV)) { stop("The values may not be normally distributed! If you know they are, try again or use a lower threshold for the Shapiro-Wilk test threshold.") }
  if (modV$Mean/modV$SD > 0.2) { warning("The data may not be normalised! For this to work, the data must be normalised.") }
  fV <- function(x) {dnorm(x, mean = 0, sd = modV$SD)}
  a <- modV$SD*5
  a <- c(1:100)*a*1.5/100
  if (!is.null(ref)) {
    r1 <- proteoCraft::is.all.good(ref)
    modR <- proteoCraft::norm.mod2(values = r1)
    if (!"Mean" %in% names(modR)) { stop("The reference values' may not be normally distributed! If you know they are, try again or use a lower threshold for the Shapiro-Wilk threshold.") }
    if (modR$Mean/modR$SD > 0.2) { warning("The data may not be normalised! For this to work, the data must be normalised.") }
    fR <- function(x) {dnorm(x, mean = 0, sd = modR$SD)}
    f3 <- function(x) {fR(x)/fV(x) - p}
  } else { f3 <- function(x) {fV(x) - p} }
  if ((!is.null(ref))&&(modR$SD >= modV$SD)) { res <- Inf } else {
    t <- f3(a)
    n1 <- max(which(t > 0))
    n2 <- min(which(t < 0))
    # This is an approximation - but for all intents and purposes good enough here:
    # it means the threshold is beyond the range of data of interest
    x1 <- a[n1]
    x2 <- a[n2]
    y1 <- t[n1]
    y2 <- t[n2]
    # Approximation of linearity:
    res <- x2 + (x1 - x2) * y2/(y2 - y1)
    # If it does not work, this loop can get us a good approximation:
    if (abs(f3(res)) > 0.0001) {
      while(abs(f3(res)) > 0.0001) {
        x <- res*c(1.5, 1, 1/1.5)
        y <- f3(x)
        s <- sign(y)
        su <- unique(sign(y))
        if (length(su) == 1) {
          y <- abs(y)
          res <- x[which(y == min(y))]
        } else {
          bx <- sapply(su, function(t) {
            t <- which(s == t)
            x1 <- x[t]
            y1 <- abs(y[t])
            x1 <- x1[which(y1 == min(y1))]
            return(x1)
          })
          by <- f3(bx)
          res <- bx[2]+(bx[1]-bx[2])*by[2]/(by[2]-by[1])
        }
      }
    }
  }
  return(res)
}
