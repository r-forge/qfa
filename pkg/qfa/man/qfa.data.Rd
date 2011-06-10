\name{qfa.data}
\alias{qfa.data}
\title{Create data object for bayesian inference using rjags}
\description{
Given timecourse data from rod.read, or alternatively a data.frame including Barcode, Row, Col, Growth, Date.Time, Inoc.Time, where each row corresponds to a single measurement of a single colony at that Date.Time,this creates the data object required to fit the Bayesian model. This can take quite a while for large numbers of colonies.
}
\usage{
qfa.data(d,edgestrip=TRUE,split=c(),fmt="\%Y-\%m-\%d_\%H-\%M-\%S",orflist=c())
}
\arguments{
  \item{d}{The data.frame containing the timecourse data for each colony.}
  \item{edgestrip}{If TRUE, the model is not fitted to colonies on the edges of plates (max and min of Row and Col).}
  \item{split}{A vector of length 2, whose first element gives the number of parts to split the ORFs in to, and whose second element gives which part to create the data object for. This allows one to split a big job up by ORFs. No splitting by default.}
  \item{fmt}{The date.time format that the inoculation times (Inoc.Time) and measurement times (Date.Time) are stored in.}
  \item{orflist}{A character vector containing the ORFs to be included in the model: all in the data by default. Split will be performed on these if specified.}
}
\value{
An R list with elements: og (growth data array), time (time data array), E (# ORFs), M (number of reps for each ORF), N (number of measurements for each ORF-repeat), ORF (list of ORFs in order, useful as ORF names not stored in JAGS).
}
\keyword{qfa}