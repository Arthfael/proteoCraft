#' cauchy.mod
#'
#' @description
#' Create a Cauchy distribution model of the data.
#' 
#' @param values The vector of Cauchy-distributed data.
#' @param print.Density Should we print the density plot? Default = FALSE.
#' @param print.QQ Should we print the Q-Q plot? Default = FALSE.
#' @param Khmaladze.threshold Which threshold can we accept for the Khmaladze test? Values accepted are 0.01, 0.025, 0.05 (default) and 0.1.
#' 
#' @examples
#' cauchy.mod(values = x)
#' 

cauchy.mod <- function(values,
                       Khmaladze.threshold = 0.05,
                       print.Density = FALSE,
                       print.QQ = FALSE) {
  warning("Deprecation warning:\nThis function hasn't been used in, like, forever!\nCheck its code and update it before using it!")
  plotEval <- function(plot) { ggplotify::as.ggplot(ggplotify::as.grob(plot)) }
  if (!Khmaladze.threshold %in% c(0.01, 0.025, 0.05, 0.1)) {
    warning("Khmaladze.threshold must be one of 0.01, 0.025, 0.05 or 0.1! Reverting to default = 0.05")
    Khmaladze.threshold <- 0.05
  }
  v <- proteoCraft::is.all.good(values)
  Khmaladze <- GofKmt::KhmaladzeTrans(v, F0 = "Cauchy")
  CritValue <- as.data.frame(t(sapply(strsplit(Khmaladze$CritValue, split = "%: "), function(x) { as.numeric(unlist(x)) })))
  colnames(CritValue) <- c("Level (%)", "Critical value")
  CritValue$Level <- as.numeric(CritValue$`Level (%)`)/100
  RES <- list()
  RES[["Khmaladze-Wilk value"]] <- Khmaladze$TestStat
  RES[["Khmaladze-Wilk critical value"]] <- CritValue$`Critical value`[which(CritValue$Level == Khmaladze.threshold)]
  RES[["Khmaladze-Wilk P-value"]] <- Khmaladze.threshold
  RES[["Follows a Cauchy distribution"]] <- Khmaladze$TestStat < RES[["Khmaladze-Wilk critical value"]]
  if (RES[["Follows a Cauchy distribution"]]) {
    tmp <- isobar::fitCauchy(v)
    RES[["Location"]] <- tmp@param@location
    RES[["Scale"]] <- tmp@param@scale
    tmp2 <- as.data.frame(v)
    colnames(tmp2) <- c("values")
    f <- function(x) {dcauchy(x, location = tmp@param@location, scale = tmp@param@scale)}
    plot1 <- ggplot2::ggplot(tmp2) +
      ggplot2::geom_density(stat = "density", ggplot2::aes(x = values)) +
      ggplot2::stat_function(fun = f, color = "red")
    RES[["Density plot"]] <- plotEval(plot1)
    #####################################################################################
    # Q-Q plot in ggplot2
    # taken from http://stackoverflow.com/questions/4357031/qqnorm-and-qqline-in-ggplot2
    #####################################################################################
    y <- quantile(v, c(0.25, 0.75))
    x <- qcauchy(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    d <- data.frame(resids = v)
    plot2 <- ggplot2::ggplot(d, ggplot2::aes(sample = resids)) +
      ggplot2::stat_qq() +
      ggplot2::geom_abline(slope = slope, intercept = int)
    RES[["Q-Q plot"]] <- plotEval(plot2)
    if (print.Density) {
      windows(height = 10, width = 10)
      print(plot1)
    }
    if (print.QQ) {
      windows(height = 10, width = 10)
      print(plot2)
    }
  } else {
    warning("The Khmaladze test failed! The data may not fit a Cauchy distribution.")
  }
  return(RES)
}
