% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse.Param.aggreg.2.R
\name{parse.Param.aggreg.2}
\alias{parse.Param.aggreg.2}
\title{parse.Param.aggreg.2}
\usage{
parse.Param.aggreg.2(
  param.nm,
  parameters = Param,
  aggregates = Aggregates,
  map = Aggregate.map,
  aggr.list = Aggregate.list,
  param.list.nm = "Param.aggreg",
  parsed.param.nm,
  obj.list.nm = ".obj"
)
}
\arguments{
\item{param.nm}{The name of the parameter to parse.}

\item{parameters}{The parameters file. Default = Param}

\item{aggregates}{The annotations (basic aggregate building blocks) as a named vector. Default: Aggregates}

\item{map}{The reference aggregates map. All existing aggregates in this map must exist. Default: Aggregate.map}

\item{aggr.list}{The list of aggregates. Default: Aggregate.list}

\item{param.list.nm}{The name of the vector of parsed parameters. Default = "Param.aggreg"}

\item{parsed.param.nm}{The name of the parsed parameter variable to create. If not specified, defaults to param.nm}

\item{obj.list.nm}{The name of the list of objects to never purge from the global environment. Default = ".obj"}
}
\description{
A new function to parse an aggregates-type parameter. Extends the functionalities of the former version.
The old function is still valid for parsing to temporary variables.
}
\examples{
parse.Param.aggreg.2("Ratios.Plot.split")

}
