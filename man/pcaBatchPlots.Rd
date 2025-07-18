% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcaBatchPlots.R
\name{pcaBatchPlots}
\alias{pcaBatchPlots}
\title{pcaBatchPlots}
\usage{
pcaBatchPlots(
  dat,
  root,
  batches = myBatches,
  map = Exp.map,
  SamplesCol = "Ref.Sample.Aggregate",
  intRoot = "log10 Intensity - ",
  dir = btchDir,
  ttl = "PCA plot - SVA batch corr.",
  openMe = FALSE,
  isRef,
  make_Avg = FALSE
)
}
\arguments{
\item{dat}{Input data matrix, expected to be log-transformed!}

\item{root}{Just a name to stick to the}

\item{batches}{Relevant batches, default = myBatches; should be valid column names of map!}

\item{map}{Map connecting samples to batches. Note that this may produce incorrect results if dat has repeated column names!}

\item{SamplesCol}{Name of the sample colum in the map, default = "Ref.Sample.Aggregate"}

\item{intRoot}{Root of the expression columns in the file, default = "log10 Intensity - "; if a character of length 0, all columns are assumed to be expression columns!}

\item{dir}{Directory where to save the plots, default = btchDir}

\item{ttl}{Plot title root, default = "PCA plot - SVA batch corr."}

\item{openMe}{Should we open the plotly plot? (default = FALSE)}

\item{isRef}{Whether a column in dat is an IRS reference column (mixed sample).}

\item{make_Avg}{In case isRef is missing, should we plot also the average of each batch? FALSE by default.}
}
\description{
A function to make PCA (ggplot and plotly) plots for the purpose of monitoring the effects of batch corrections.
}
\examples{
  scoresLst <- PCAlyLst <- PCsLst <- list()
  tmp <- pcaBatchPlots(imputedExpdata,
                       "original",
                       map = smplsMap,
                       SamplesCol = "New")
  PCAlyLst[["original"]] <- tmp$PlotLy
  scoresLst[["original"]] <- tmp$Scores
  PCsLst[["original"]] <- tmp$PCs
  tmp <- pcaBatchPlots(imputedExpdata2,
                       "corrected",
                       map = smplsMap,
                       SamplesCol = "New")
  PCAlyLst[["corrected"]] <- tmp$PlotLy
  scoresLst[["corrected"]] <- tmp$Scores
  PCsLst[["corrected"]] <- tmp$PCs

}
