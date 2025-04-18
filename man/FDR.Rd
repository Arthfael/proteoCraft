% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDR.R
\name{FDR}
\alias{FDR}
\title{FDR}
\usage{
FDR(
  data,
  aggregate,
  pvalue_root,
  pvalue_col,
  fdr,
  returns = c(TRUE, FALSE),
  method = "BH"
)
}
\arguments{
\item{data}{The dataframe that contains the values.}

\item{aggregate}{The current aggregate to be pasted to the root of P-value column names to create the column name.}

\item{pvalue_root}{The root of the name of the P-values column. Default = "Ratio_Minus.log10.Pvalue." If the tag "log10" is found in that name, P-values will assume to be -log10 transformed and will be reverse transformed accordingly.}

\item{pvalue_col}{Alternatively to the method of constructing a P-values column name used above, you may simply provide its name with this argument.}

\item{fdr}{The acceptable False Discovery Rate values.}

\item{returns}{Default = c(TRUE, FALSE) Set first to TRUE to output vector of significance; set second to TRUE to output threshold values (NB: these are not-log transformed).}

\item{method}{The method to use. Default = "BH" for Benjamini-Hochberg. Alternatively use "BY" for Benjamini-Yekutieli.}
}
\description{
A function to calculate FDR thresholds and assess P-value significance.
Default is the Benjamini-Hochberg procedure.
Can also do the Benjamini-Yekutieli procedure under positive dependence assumptions.
}
\examples{
temp <- FDR(data = PG,
            aggregate,
            pvalue_root = "Ratio_Minus.log10.Pvalue.",
            fdr = BH.FDR)

}
