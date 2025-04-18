% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listMelt.R
\name{listMelt}
\alias{listMelt}
\title{listMelt}
\usage{
listMelt(List, Names = NULL, ColNames = c("value", "L1"))
}
\arguments{
\item{List}{A list.}

\item{Names}{Names to apply to each element, if missing. Can be convenient in case you want to preserve Name's original non-character data type (e.g. numeric).}

\item{ColNames}{Column names of the output data frame. Default = c("value", "L1"); note the reverse order compared with classic melt!}
}
\description{
Wrapper function for utils::stack(), allowing much faster and efficient melting of simple lists, i.e. lists where each element is a 1-dimension vector.
This won't work on more complex lists!!!

Allows the following:
 - Much faster melting than reshape::melt.list() for large lists
 - Avoids using equally slow reshape::melt(), which is apparently deprecated and causing issues with data.tables
 - Simplifies the syntax of the rewrites
}
\examples{
temp <- listMelt(strsplit(PG$"Protein IDs", ";"), PG$id, c("Protein accession", "PG id"))

}
