% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyNumInput.R
\name{shinyNumInput}
\alias{shinyNumInput}
\title{shinyNumInput}
\usage{
shinyNumInput(
  startVal,
  min = NA,
  max = NA,
  step = NA,
  naRep = NA,
  width = "300px",
  root
)
}
\arguments{
\item{startVal}{Starting values of each numeric value in the column.}

\item{min}{Minimum allowed value}

\item{max}{Maximum allowed value}

\item{step}{Interval to use when stepping between min and max}

\item{naRep}{What to replace NA values with? Default = NA; set e.g. to 1 to use 1 as default replacement for incorrect initial inputs.}

\item{width}{The width of the input, e.g. "400px", or "100%"; see htmltools::validateCssUnit().}

\item{root}{Root of the created input names. Names generated will be of the form "root___i", where i is the 1:n row number.}
}
\description{
A function to create numeric input columns for use in DT data table columns.
}
\examples{
# Use within a shiny server as below:

server <- function(input, output, session) {

  Table$Counts <- shinyNumInput(Table$Counts, 0, Inf, 1, root = "Counts")
  Table$Intensity <- shinyNumInput(Table$Counts, root = "Intensity")

}

}
