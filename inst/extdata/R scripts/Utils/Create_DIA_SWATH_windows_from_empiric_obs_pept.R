# Header
options(stringsAsFactors = FALSE)
library(ggplot2)
library(reshape2)
library(proteoCraft)
library(modeest)
library(MASS)
library(truncdist)
require(minpack.lm)

# Note:
# This appears to be useless on a QE-HF, as variable window size is not supported.
# Despite what Moritz told me, it does not appear that creating one DIA experiment per DIA window works:
# The software then complains about the lack of inclusion list, and the list must be global.

# Choose working directory:
wd <- rstudioapi::selectDirectory()
setwd(wd)

# Choose destination directory:
destd <- rstudioapi::selectDirectory()

# Choose MaxQuant evidence.txt file to process:
ev <- rstudioapi::selectFile()
ev <- read.delim(ev)
rm(list = ls()[which(!ls() %in% c("wd", "ev", "destd"))])

ev$temp <- apply(ev[,c("Modified.sequence", "Charge")], 1, function(x) {paste(x, collapse = "_")})
test <- aggregate(ev$m.z, list(ev$temp), unique); colnames(test) <- c("Unique_ID", "m.z")
test$Charge <- as.factor(gsub("^[^0-9]+", "", test$Unique_ID))
test <- test[order(test$m.z),]
test$"log10(m/z)" <- log10(test$m.z)

N.windows <- 40
Upper <- 1400
Lower <- 400

# Estimate distribution:
bin.scale <- max(50, round(nrow(test)/100)) # (at least 50 bins...)
bins <- c(0:bin.scale)*(Upper-Lower)/bin.scale + Lower
halfblur <- 10
kounts <- Isapply(c(1:bin.scale), function(x) {
  x1 <- length(which((test$m.z >= bins[x])&(test$m.z < c(bins, bins[length(bins)+1])[x+1])))
  x2 <- length(which((test$m.z >= bins[max(1,x-halfblur)])&(test$m.z < bins[min(bin.scale, x-1+halfblur)])))/(2*halfblur)
  return(c(x1, x2))
})
bins <- data.frame(m.z = (bins[1:bin.scale]+bins[2:(bin.scale+1)])/2,
                   counts = kounts$V1,
                   blurred.counts = kounts$V2)
for (i in c("", "blurred.")) {
  bins[[paste0(i, "density")]] <- bins[[paste0(i, "counts")]]/sum(bins[[paste0(i, "counts")]])
}
bins$"log10(m/z)" <- log10(bins$m.z)
loki <- mlv(test$"log10(m/z)", method = "mfv")
skael <- sd(test$"log10(m/z)")
fit <- fitdistr(x = test$"log10(m/z)", "normal")
phun <- function(par, x) {
  x$DTrunc <- dtrunc(x = x$"log10(m/z)", spec = "norm", a = log10(Lower), b = log10(Upper),
                     mean = par[1], sd = par[2])
  x$DTrunc <- x$DTrunc/sum(x$DTrunc)
  return(x$density - x$DTrunc)
}
LM <- nls.lm(fn = phun, par = fit$estimate, x = bins)
model <- LM$par

## Calculate windows:
# Empirically:
temp <- c(1:(N.windows-1))*nrow(test)/N.windows
temp <- sapply(temp, function(x) {
  mean(c(max(test$m.z[floor(x)]), min(test$m.z[ceiling(x)])))
  # I guess there are cases where I found such complicated calculations necessary...
})
temp <- c(Lower, temp, Upper)
temp <- round(temp, 1)
SWATH.windows <- as.data.frame(t(sapply(2:length(temp), function(x) {
  return(c(temp[x-1]-0.5, temp[x]+0.5))
})))
colnames(SWATH.windows) <- c("Start", "End")
SWATH.windows$Y <- c(1:nrow(SWATH.windows))

# From estimated normal distribution:
left <- pnorm(q = log10(Lower), model[["mean"]], model[["sd"]])
right <- pnorm(q = log10(Upper), model[["mean"]], model[["sd"]])
span <- right - left
qnorm(p = left, model[["mean"]], model[["sd"]])
breaks <- round(10^sapply(c(0:N.windows), function(x) {qnorm(p = left + span*x/40, model[["mean"]], model[["sd"]])}), 1)
SWATH.windows$Model_Start <- breaks[1:N.windows]-0.5
SWATH.windows$Model_End <- breaks[2:(N.windows+1)]+0.5
SWATH.windows$Model_Centre <- round((SWATH.windows$Model_Start+SWATH.windows$Model_End)/2, 1)
SWATH.windows$Model_Width <- round(SWATH.windows$Model_End-SWATH.windows$Model_Start, 1)

# Plots
plot2 <- ggplot(SWATH.windows) + geom_segment(aes(x = Start, xend = End, y = Y, yend = Y), colour = "green") +
  geom_segment(aes(x = Model_Start, xend = Model_End, y = Y-0.1, yend = Y-0.1), colour = "red") +
  theme_bw()
poplot(plot2)
phun <- function(x) { dnorm(x, mean = model[["mean"]], sd = model[["sd"]]) }
scaling <- mean(bins$blurred.density[c(max(which(bins$"log10(m/z)" <= loki)), min(which(bins$"log10(m/z)" >= loki)))])
scphun <- function(x) { phun(x)*scaling/phun(loki) }
temp <- melt(bins[,c("log10(m/z)", "density", "blurred.density")], id.vars = "log10(m/z)")
colnames(temp)[which(colnames(temp) == "value")] <- "density"
plot <- ggplot(temp) + geom_point(aes(x = `log10(m/z)`, y = density, colour = variable)) +
  stat_function(fun = scphun, color = "red") +
  geom_vline(data = SWATH.windows, aes(xintercept = log10(Model_End)),
             colour = "red", linetype = "dotted") + theme_bw()
poplot(plot, height = 10, width = 10)

SWATH.windows2 <- SWATH.windows[,c("Model_Start", "Model_End")]
colnames(SWATH.windows2) <- gsub("^Model_", "", colnames(SWATH.windows2))
write.table(SWATH.windows2, file = paste0(destd, "\\Non-linear_isolation_windows.txt"),
            sep = "\t", row.names = FALSE)

SWATH.windows3 <- SWATH.windows[,c("Model_Centre", "Model_Width")]
colnames(SWATH.windows3) <- gsub("^Model_", "", colnames(SWATH.windows3))
write.table(SWATH.windows3, file = paste0(destd, "\\Non-linear_isolation_windows_centre_and_width.txt"),
            sep = "\t", row.names = FALSE)
shell(paste0("explorer \"", destd, "\""), intern = TRUE)
openwd()
