% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_aggr.R
\name{list_aggr}
\alias{list_aggr}
\title{list_aggr}
\usage{
list_aggr(
  origin,
  dest,
  list.col,
  other.col,
  dest.match.col,
  insep = ";",
  outsep = ";"
)
}
\arguments{
\item{origin}{The origin data frame.}

\item{dest}{The destination data frame.}

\item{list.col}{The (character) column of IDs to use for the melting process. See the "insep" argument.}

\item{other.col}{The other column(s) of interest to be melted.}

\item{dest.match.col}{The name of the column in the destination data frame whose IDs will be matched to those listed in list.col.}

\item{insep}{The separator in list.col; default = ";". Change to the correct separator (e.g. ",") if needed.}

\item{outsep}{The separator to use in the output.}
}
\description{
A function to create melt lists.
}
\examples{
test <- data.frame("values" = c(1:3))
test$ids <- list("a;b", "c", "d;e")
test2 <- data.frame(ids = c("a", "b", "c", "d", "e", "f"))
list_aggr(test, test2, "ids", "values", "ids", insep = ";", outsep = ";")

}
