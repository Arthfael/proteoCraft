% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanNms.R
\name{cleanNms}
\alias{cleanNms}
\title{cleanNms}
\usage{
cleanNms(
  Samples,
  Experiment = Exp,
  sep = "___",
  rep = " ",
  start = TRUE,
  end = FALSE
)
}
\arguments{
\item{Samples}{The vector of sample names, or of column names containing sample names. Each name should be made of "___"-aggregated factor levels.}

\item{Experiment}{The Experiment object.}

\item{sep}{Expected separator. Default = "___". NB: avoid too simple separators which may be used as also parts of factor level names, otherwise weird things may happen.}

\item{rep}{Replacement, user-friendlier separator. Default = " "}

\item{start}{Default = TRUE, which means the pattern will require a match to the start of the string.}

\item{end}{Default = FALSE. Set to TRUE for the pattern to require a match to the end of the string.}
}
\description{
Meant as an internal function to clean up "___"-aggregated names for nicer plot annotations.
(NB: should work for other separators, we are just using "___" here).
For instance, it will turn "Exp1___Treated___KO___1" and "Exp1___Treated___KO" into the much user friendlier "Exp1 Treated KO 1" and "Exp1 Treated KO",
or, if - as is usually the case - there is only one Experiment, into "Treated KO 1" and "Treated KO". 

(Of course, ideally this should not even be necessary. But when I started what would become these workflows nearly a decade ago I did not even know how to use lists in R...
At some point it would be nice to rewrite everything using Exp.map to match individual sample names to their values for each factor and their corresponding column names for each test etc...
One more thing to do...)
}
