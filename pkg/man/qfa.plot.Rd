\name{qfa.plot}
\alias{qfa.plot}
\title{Plots fitted model and data for all the colonies in results of qfa.fit}
\description{
Produces a pdf of plate by plate arrays of timecourse plots showing the fitted curve from qfa.fit and the data for each colony. Useful way to inspect model fit & data.
}
\usage{
qfa.plot(file,results,d,fmt="\%Y-\%m-\%d_\%H-\%M-\%S",barcodes=c(),
master.plates=c(),treatments=c(),screen.names=c(),backgrounds=c())}
\arguments{
  \item{file}{The file to output the plots to.}
  \item{results}{The output of qfa.fit which contains the fitted curve parameters of colony growth you wish to plot.}
  \item{d}{The original data.frame fed to qfa.fit containing all of the timecourse data}
  \item{fmt}{The format in which Date.Time of measurement and inoculation time are stored}
  \item{barcodes}{Plot only for the plates with barcodes in this character vector; all by default.}
  \item{master.plates}{Plot only for the plates from master.plates in this character vector; all by default.}
  \item{treatments}{Plot only for the plates with treatments in this character vector; all by default.}
  \item{screen.names}{Plot only for the plates with screen.names in this character vector; all by default.}
  \item{backgrounds}{Plot only for the plates with genetic background in this character vector; all by default.}
}
\examples{
}
\keyword{qfa}