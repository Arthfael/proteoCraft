% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data_Impute.R
\name{Data_Impute}
\alias{Data_Impute}
\title{Data_Impute}
\usage{
Data_Impute(data, method = "normal", seed)
}
\arguments{
\item{data}{A matrix or data.frame with more than 1 column, containing missing or invalid (e.g. infinite).}

\item{method}{By default, "normal". Alternative: "cauchy". Note that for Cauchy we are using Mean as Location and SD as Scale.}

\item{seed}{Use a fixed seed for reproducibility purposes.}
}
\description{
A function to replace missing values in a data table with imputed values.
Useful for when you have to apply a method which requires a table full of finite numeric data.
The function calculates the SD and Mean per row and fills blanks with random samplings of a normal (default) distribution with the same SD and Mean.
The output is a list with 2 elements:
- The data frame with added imputations.
- The positions that were imputed. This is so after whichever method is to be applied on the data has been applied
}
\examples{
New_Data <- Data_Impute(data)

}
