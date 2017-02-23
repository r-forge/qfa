\name{numericalfitness}
\alias{numericalfitness}
\title{Numerical fitness estimates from timecourse data}
\description{
Generates numerical area under curve (nAUC) and numerical single time point (nSTP) measures of fitness.  This function generates a piecewise linear interpolation between the available datapoints in the data frame obsdat.  Interpolation is continued beyond the observed range of times assuming a density equal to the nearest observed time.  

nAUC is the area under the interpolating function, and is estimated by integrating between limits 0 and AUCLim.  nSTP is the interpolated cell density estimate at t=STP.    
}
\usage{
numericalfitness(obsdat,AUCLim,STP,nrate=FALSE)
}
\arguments{
  \item{obsdat}{A data frame containing a numeric Expt.Time column (time since culture inoculation, typically in days) and a numeric Growth column (typically cell density estimates estimated from a timeseries of culture photographs, quantified by software such as Colonyzer: http://research.ncl.ac.uk/colonyzer/)}
  \item{AUCLim}{Upper limit of integration for calculating area under growth curve as a measure of culture fitness}
  \item{STP}{Single time point at which to estimate cell density (measure of fitness, as used in SGA for example: http://en.wikipedia.org/wiki/Synthetic_genetic_array)}
  \item{nrate}{Calculate numerical growth rate phenotypes and add these to the output results.}
}
\value{
A named vector specifying nAUC and nSTP.
}
\keyword{qfa}