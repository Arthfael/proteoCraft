#' norm.mod2
#'
#' @description 
#' Create a normal distribution model for normally distributed data.
#' 
#' 
#' @param values The vector of normally distributed data
#' @param shapiro.threshold Which threshold can we accept for the Shapiro-Wilk normality test's P-value? Default = 0.05
#' @param print.Density Should we print the density plot? Default = FALSE.
#' @param print.QQ Should we print the Q-Q plot? Default = FALSE.
#' 
#' @examples
#' norm.mod2(values = x)
#' 
#' @export

norm.mod2 <- function(values, shapiro.threshold = 0.05, print.Density = FALSE, print.QQ = FALSE) {
  v <- proteoCraft::is.all.good(values)
  if (length(v) <= 5000) { v1 <- v } else { v1 <- sample(v, size = 5000, replace = FALSE, prob = NULL) }
  shapiro <- shapiro.test(v1)
  RES <- list()
  RES[["Shapiro-Wilk P-value"]] <- shapiro$p.value
  RES[["Shapiro-Wilk P-value threshold"]] <- shapiro.threshold
  RES[["Is normal"]] <- RES[["Shapiro-Wilk P-value"]] >= RES[["Shapiro-Wilk P-value threshold"]]
  if (RES[["Is normal"]]) {
    fit <- MASS::fitdistr(x = v, densfun = "normal")
    M <- fit$estimate[["mean"]]
    RES[["Mean"]] <- M
    SD <- fit$estimate[["sd"]]
    RES[["SD"]] <- SD
    require(ggplot2)
    tmp <- as.data.frame(v)
    colnames(tmp) <- c("values")
    f <- function(x) { dnorm(x, mean = M, sd = SD) }
    plot1 <- ggplot2::ggplot(tmp) +
      ggplot2::geom_density(stat = "density",
                            ggplot2::aes(x = values)) +
      ggplot2::stat_function(fun = f, colour = "red") +
      ggplot2::geom_text(label = paste0("Shapiro test P-value = ", shapiro$p.value),
                         x = -2.975*SD, y = f(0)*1.095, hjust = 0) +
      ggplot2::xlim(-3*SD, 3*SD) + ggplot2::ylim(0, f(0)*1.1) +
      ggplot2::theme_bw()
    RES[["Density plot"]] <- plot1
    #####################################################################################
    # Q-Q plot in ggplot2
    # taken from http://stackoverflow.com/questions/4357031/qqnorm-and-qqline-in-ggplot2
    #####################################################################################
    y <- quantile(v, c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    d <- data.frame(resids = v)
    plot2 <- ggplot2::ggplot(d,
                             ggplot2::aes(sample = resids)) + ggplot2::stat_qq() +
      ggplot2::geom_abline(slope = slope, intercept = int)
    RES[["Q-Q plot"]] <- plot2
    if (print.Density) {
      windows(height = 10, width = 10)
      print(plot1)
    }
    if (print.QQ) {
      windows(height = 10, width = 10)
      print(plot2)
    }
  } else { warning("The data is not normal!!!") }
  return(RES)
}
