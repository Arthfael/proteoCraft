% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chooseFrom.R
\name{chooseFrom}
\alias{chooseFrom}
\title{chooseFrom}
\usage{
chooseFrom(title, choices, default = 1, multiple = FALSE)
}
\arguments{
\item{title}{The caption displayed at the top of the choices.}

\item{choices}{Available items to choose from.}

\item{default}{Integer, the choice by default (default = 1).}

\item{multiple}{Are multiple choices allowed? Default = FALSE}
}
\description{
A wrapper function to allow easy selection from a range of options.
Completely useless really... what was I thinking.
It does in a more complicated way what the wrapped function did...
Candidate for deletion...
}
\examples{
temp <- chooseFrom("What do you want?", c("Peace", "War"), 2) # (Please choose 1!)

}
