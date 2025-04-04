#' Thresholds
#'
#' @description 
#' Find thresholds to include only a given proportion of data.
#' The data has to be normalized as the function assumes symmetry of the problem relative to 0!!!
#' Unused, old, deprecated, candidate for deletion.
#' 
#' @param values The vector of distributed values for which we want a threshold.
#' @param ref The reference. If left as NULL, proportions are from the values vector. If provided as a vector of values, the function aims to find thresholds including a given proportion (based on modelled distributions) p of reference to values.
#' @param p The relevant proportion to include (either absolute, or contamination by a reference).
#' @param Pvalue_threshold Which threshold can we accept for the distribution test's P-value? Default = 0.05
#' @param Distribution Which Distribution does the data follow. Allowed values are "Normal" and "Cauchy" (default).
#' 
#' @examples
#' Thresholds(values = test$Treated,
#'            ref = test$Control)
#' 
#' @export

Thresholds <- function(values, ref = NULL, p = 0.01, Pvalue_threshold = 0.05, Distribution = "Cauchy") {
  stopifnot(p > 0, p < 1, Distribution %in% c("Cauchy", "Normal"))
  if (Distribution == "Cauchy") { Location <- "Location"; Scale <- "Scale"; FUN <- cauchy.mod }
  if (Distribution == "Normal") { Location <- "Mean"; Scale <- "SD"; FUN <- proteoCraft::norm.mod2 }
  v1 <- proteoCraft::is.all.good(values)
  modV <- FUN(values = v1)
  if (!Location %in% names(modV)) { stop("The chosen distribution may be incorrect! If you know it is correct, try again or use a lower threshold for the distribution test's threshold.") }
  LocationV <- modV[[Location]]; ScaleV <- modV[[Scale]]
  if (Distribution == "Cauchy") { fV <- function(x) { dcauchy(x, location = 0, scale = ScaleV) } }
  if (Distribution == "Normal") { fV <- function(x) {dnorm(x, mean = 0, sd = ScaleV) } }
  a <- 5*ScaleV
  a <- c(1:100)*a*1.5/100
  if (abs(LocationV/ScaleV) > 0.5) { warning("The data may not be normalised! For this to work, the data must be normalised.") }
  if (!is.null(ref)) {
    r1 <- proteoCraft::is.all.good(ref)
    modR <- FUN(values = r1)
    if (!Location %in% names(modR)) { stop("The chosen distribution may be incorrect! If you know it is correct, try again or use a lower threshold for the distribution test's threshold.") }
    LocationR <- modR[[Location]]; ScaleR <- modR[[Scale]]
    if (Distribution == "Cauchy") { fR <- function(x) { dcauchy(x, location = 0, scale = ScaleR) } }
    if (Distribution == "Normal") { fR <- function(x) { dnorm(x, mean = 0, sd = ScaleR) } }
    f3 <- function(x) {fR(x)/fV(x) - p}
    if (abs(LocationR/ScaleR) > 0.5) { warning("The data may not be normalised! For this to work, the data must be normalised.") }
  } else { f3 <- function(x) {fV(x) - p} }
  if ((!is.null(ref))&&(modR[[Scale]] >= modV[[Scale]])) {
    res <- Inf
  } else {
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
          bx <- sapply(su, function(tt) {
            tt <- which(s == tt)
            x1 <- x[tt]
            y1 <- abs(y[tt])
            return(x1[which(y1 == min(y1))])
          })
          by <- f3(bx)
          res <- bx[2]+(bx[1]-bx[2])*by[2]/(by[2]-by[1])
        }
      }
    }
  }
  return(res)
}
