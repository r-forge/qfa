\name{qfa.acf}
\alias{qfa.acf}
\title{Check autocorrelations in qfa.simulate posterior samples}
\description{Calculates and plots autocorrelation functions for model parameter posteriors.}
\usage{qfa.acf(posterior,lag.max=100,plot=TRUE,tol=0.2,printout=FALSE)}
\arguments{
  \item{posterior}{qfa.posterior output with posteriors for the variables you wish to calculate ACFs for.}
  \item{lag.max}{The maximum lag that the ACFs should be calculated to.}
  \item{plot}{If TRUE, plots the ACF for every posterior}
  \item{tol}{This determines what magnitude of autocorrelation constitutes a problematic autocorrelation, and horizontal lines are drawn on the ACF plots at (+/-tol)}
  \item{printout}{If TRUE, this prints out the number, and proportion, of autocorrelations with magnitude greater than tol for each posterior.}
}
\value{
R list with elements ACF (list of autocorrelation functions for all the parameters), nAC (Vector of number of abs(autocorrelations)>tol for each parameter), pAC (Vector of % of autocorrelations out of those checked with magnitude greater than tol).
}
\examples{
}
\keyword{qfa}