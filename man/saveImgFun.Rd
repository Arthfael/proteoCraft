% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveImgFun.R
\name{saveImgFun}
\alias{saveImgFun}
\title{saveImgFun}
\usage{
saveImgFun(file)
}
\arguments{
\item{file}{The path/file to save to.}
}
\description{
Function to save all objects in the environment quickly.
The backend used is qs2::qs_savem under Windows and fastSave::save.lbzip2 under Linux.
}
\examples{
proteoCraft::saveImgFun("Backup.qs2")
proteoCraft::loadFun("Backup.qs2")

}
