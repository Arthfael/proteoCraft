% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Isapply.R
\name{Isapply}
\alias{Isapply}
\title{Isapply}
\usage{
Isapply(
  x,
  FUN,
  ...,
  simplify = TRUE,
  USE.NAMES = TRUE,
  col.names,
  convert.to = "df"
)
}
\arguments{
\item{x}{a list, typically returned from \code{lapply()}.}

\item{FUN}{the function to be applied to each element of \code{X}:
    see \sQuote{Details}.  In the case of functions like
    \code{+}, \code{\%*\%}, the function name must be backquoted or quoted.}

\item{...}{optional arguments to \code{FUN}.}

\item{simplify}{logical or character string; should the result be
    simplified to a vector, matrix or higher dimensional array if
    possible?  For \code{sapply} it must be named and not abbreviated.
    The default value, \code{TRUE}, returns a vector or matrix if appropriate,
    whereas if \code{simplify = "array"} the result may be an
    \code{\link[base]{array}} of \dQuote{rank}
    (\eqn{=}\code{length(dim(.))}) one higher than the result
    of \code{FUN(X[[i]])}.}

\item{USE.NAMES}{logical; if \code{TRUE} and if \code{X} is character,
    use \code{X} as \code{\link[base]{names}} for the result unless it had names
    already.  Since this argument follows \code{\dots} its name cannot
    be abbreviated.}

\item{col.names}{Optional column names to provide to the output, if it will be 2D. Ignored (with a warning) if output is 1D.}

\item{convert.to}{By default ("df" or 1), converts 2D output to a data.frame. If set as "m" or "2", converts it to a matrix.}
}
\description{
I love the sapply function, but its output can be difficult to deal with:
It can require different treatment depending on whether results per x are expected to be length 1 or more.
This is a wrapper so I don't have to rewrite the same stuff all the time, hopefully it will help.
It does not always work currently...
}
