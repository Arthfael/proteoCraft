% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data_Impute2.R
\name{Data_Impute2}
\alias{Data_Impute2}
\title{Data_Impute2}
\usage{
Data_Impute2(quant_data, groups, is.log = TRUE, seed)
}
\arguments{
\item{quant_data}{A matrix or data.frame with more than 1 column, containing missing or invalid (e.g. infinite).}

\item{groups}{A vector of groups of length equal to the number of columns of the input data frame.}

\item{is.log}{Default = TRUE. Are the input values log transformed? If set to FALSE, the data will be log-transformed prior to processing, then de-logged.}

\item{seed}{Use a fixed seed for reproducibility purposes.}
}
\description{
Better function to impute missing values than my former version, "Data_Impute". Uses the imputeLCMD package.
Works in two steps:
 - Step 1: within user-defined groups of user-defined related samples (replicates from a treatment),
           on which it performs MAR imputation using KNN.
 - Step 2: globally on all remaining missing values, which it imputed as MNAR using QRILC.
The output is a list of the two arrays: new quantitative values and positions which were imputed.
}
\examples{
New_Data <- Data_Impute2(data)

}
