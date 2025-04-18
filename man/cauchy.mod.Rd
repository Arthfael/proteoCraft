% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cauchy.mod.R
\name{cauchy.mod}
\alias{cauchy.mod}
\title{cauchy.mod}
\usage{
cauchy.mod(
  values,
  Khmaladze.threshold = 0.05,
  print.Density = FALSE,
  print.QQ = FALSE
)
}
\arguments{
\item{values}{The vector of Cauchy-distributed data.}

\item{Khmaladze.threshold}{Which threshold can we accept for the Khmaladze test? Values accepted are 0.01, 0.025, 0.05 (default) and 0.1.}

\item{print.Density}{Should we print the density plot? Default = FALSE.}

\item{print.QQ}{Should we print the Q-Q plot? Default = FALSE.}
}
\description{
Create a Cauchy distribution model of the data.
}
\examples{
cauchy.mod(values = x)

}
