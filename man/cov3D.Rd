% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov3D.R
\name{cov3D}
\alias{cov3D}
\title{cov3D}
\usage{
cov3D(
  pdb,
  peptides,
  colscale,
  path,
  intensities = NULL,
  display = TRUE,
  asRatios = FALSE,
  I_eq_L = TRUE
)
}
\arguments{
\item{pdb}{A PDB file, either already read into R or a path to an existing such file.\cr The first model in the file will be detected, parsed, and the collection of peptides of interest will be mapped onto the model.}

\item{peptides}{The collection of peptide sequences to map on the protein sequence(s).}

\item{colscale}{Name of a color scale acceptable by plotly. The default is "viridis" if asRatios = FALSE and "plasma" if asRatios = TRUE}

\item{path}{Where to save our httml file.}

\item{intensities}{Values to map to the "colscale" argument so that peptides can be printed with different colors, e.g. mapped to abundance.}

\item{display}{Logical: should we open the plot? TRUE by default.}

\item{asRatios}{Is the data fed to the intensities argument actually ratios data (default = FALSE)? If TRUE, this changes the default color scale from "viridis" to "plasma"}

\item{I_eq_L}{Should we consider I and L identical? Currently, by default, TRUE for both DIA and DDA: see https://github.com/vdemichev/DiaNN/discussions/1631}
}
\description{
This function can be used to create a 3d coverage map of identified pepides overlaid over a protein model in PDB format.
}
\examples{
cov3D(pdb, peptides)

}
