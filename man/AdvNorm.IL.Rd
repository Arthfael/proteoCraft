% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AdvNorm.IL.R
\name{AdvNorm.IL}
\alias{AdvNorm.IL}
\title{AdvNorm.IL}
\usage{
AdvNorm.IL(df, ids.col = "Unique State", exprs.col, exprs.log = FALSE, K = 5)
}
\arguments{
\item{df}{The input data frame, containing unique peptide identifiers, fraction numbers, and expression data.}

\item{ids.col}{The name of the column in df containing peptide identifiers. Default is "Unique State", a column I usually create which contains modified sequence pasted to charge.}

\item{exprs.col}{The name of the columns in df which contain expression values.}

\item{exprs.log}{Default = FALSE, meaning the input data is not log-transformed. Set instead to the base of the log if input data is log-transformed. If set to TRUE, the data is assumed to be log10 base. The data returned will be transformed or not as per the input.}

\item{K}{A fold factor defining the minimum (1/K) and maximum (*K) range of accepted values for correction factors.}
}
\description{
A function that performs simplified advanced normalisation on the columns of a data.frame.
The ".IL" suffix in the name is because at first it was written with a TMT or iTRAQ reporter intensities data frame in mind.
However, this can work for any data frame of quantitative values.
Input must contain a column of fractions names, IDs, and 2 or more expression columns.

This takes long but is better than median only.

11/01/2024 - following major rewriting of what was a very poorly written function,
and in parts thanks to the data.table package, processing has been sped up ~60 times!
}
\examples{
adv.norm.data <- AdvNorm.IL(data, "Unique State", paste0("Reporter intensity ", c(0:9)), FALSE, K = 5)
require(ggplot2)
require(reshape2)
test <- melt(adv.norm.data)
test$value <- log10(test$value)
plot <- ggplot(test) + geom_density(stat = "density", aes(x = value, colour = variable))
grDevices::windows(width = 10, height = 10)
print(plot)

}
