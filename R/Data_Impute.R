#' Data_Impute
#'
#' @description
#' A function to replace missing values in a data table with imputed values.
#' Useful for when you have to apply a method which requires a table full of finite numeric data.
#' The function calculates the SD and Mean per row and fills blanks with random samplings of a normal (default) distribution with the same SD and Mean.
#' The output is a list with 2 elements:
#' - The data frame with added imputations.
#' - The positions that were imputed. This is so after whichever method is to be applied on the data has been applied
#' 
#' @param data A matrix or data.frame with more than 1 column, containing missing or invalid (e.g. infinite).
#' @param method By default, "normal". Alternative: "cauchy". Note that for Cauchy we are using Mean as Location and SD as Scale.
#' @param seed Use a fixed seed for reproducibility purposes.
#' 
#' @examples
#' New_Data <- Data_Impute(data)
#' 
#' @export

Data_Impute <- function(data, method = "normal", seed) {
  if (missing(seed)) { if (exists("mySeed")) { seed <- mySeed } else { seed <- 1234567 } }
  set.seed(seed)
  norm::rngseed(seed)
  #
  Klass <- class(data)
  if (!Klass %in% c("data.frame", "matrix")) { stop("The \"data\" must be a matrix or data.frame!") }
  method <- tolower(method)
  if(!method %in% c("normal", "cauchy")) {
    warning("I can only deal with methods \"normal\" or \"cauchy\", defaulting to \"normal\"!")
  }
  if (method == "normal") { samfun <- rnorm }
  if (method == "cauchy") { samfun <- rcauchy }
  kol <- colnames(data)
  if (is.null(kol)) { kol <- paste0("C", c(1:ncol(data))) }
  l <- length(kol)
  misses <- as.data.frame(sapply(kol, function(x) { !proteoCraft::is.all.good(data[[x]], 2) }))
  test <- as.data.frame(t(apply(data, 1, function(x) {
    x <- proteoCraft::is.all.good(x)
    return(c(mean(x), sd(x)))
  })))
  colnames(test) <- c("mean", "sd")
  w <- which(proteoCraft::is.all.good(test$sd, 2))
  a <- loess(sd ~ mean, data = test[w,], surface = "direct")
  test$loess <- predict(a, newdata = test)
  w1 <- which(proteoCraft::is.all.good(test$loess, 2))
  w2 <- which((!proteoCraft::is.all.good(test$sd, 2))|(test$sd == 0))
  test$sd[w2] <- sapply(test$mean[w2], function(x) {
    y1 <- which(test$mean[w] > x); y2 <- which(test$mean[w] < x)
    if (length(y1) > 0) {
      y1 <- min(test$mean[w][y1]); y1 <- mean(test$loess[which(test$mean[w1] == y1)])
    } else { y1 <- as.numeric(NA) }
    if (length(y2) > 0) {
      y2 <- max(test$mean[w][y2]); y2 <- mean(test$loess[which(test$mean[w1] == y2)])
    } else { y2 <- as.numeric(NA) }
    x <- mean(proteoCraft::is.all.good(c(y1, y2)))
    return(x)
  })
  data[,c("Mean", "SD")] <- test[,c("mean", "sd")]
  data$which <- apply(data[,kol], 1, function(x) { list(which(!proteoCraft::is.all.good(unlist(x), 2))) })
  W <- which(sapply(data$which, function(x) { length(unlist(x)) > 0 }))
  data2 <- as.data.frame(t(apply(data[W,c(kol, "Mean", "SD", "which")], 1, function(x) {
    x <- unlist(x)
    x1 <- x[1:l]; m <- x[l+1]; sd <- x[l+2]; x2 <- x[(l+3):length(x)]
    l2 <- length(x2)
    x1[x2] <- samfun(l2, m, sd)
    return(x1)
  })))
  data <- data[,kol]
  data[W,] <- data2
  if (class(data) != Klass) { data <- get(paste0("as.", Klass))(data) }
  res <- list("Imputed_data" = data,
              "Positions_Imputed" = misses)
  return(res)
}
