\name{qfa.plot}
\alias{qfa.plot}
\title{Plots fitted model and data for all the colonies in results of qfa.fit}
\description{
Produces a multipage pdf of growth curves.  Each page corresponds to a single plate and growth curves are arrayed on the page according to their position on the plate.  Both observations and the fitted growth curve from qfa.fit are shown for each culture.  Where available, various fitness estimates are displayed for each culture, together with culture genotype.  These .pdfs are useful for visually checking quality of model fit & data.
}
\usage{
qfa.plot(file,results,d,fmt="\%Y-\%m-\%d_\%H-\%M-\%S",barcodes=c(),
master.plates=c(),treatments=c(),screen.names=c(),backgrounds=c(),
maxg=0,maxt=0,logify=FALSE)}
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
  \item{maxg}{Upper cell density (y-axis limit) for all growth curve plots.  Default value is zero.  If maxg=0, then upper fitnesses are chosen automatically to show all datapoints.}
  \item{maxt}{Growth curve is plotted from time t = 0 to maxt.}
  \item{logify}{Boolean indicating whether growth curve plots should be on a semilog scale.  Cell density (y-axis) only.}
}
\examples{
}
\keyword{qfa}