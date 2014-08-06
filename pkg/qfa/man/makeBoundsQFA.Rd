\name{makeBoundsQFA}
\alias{makeBoundsQFA}
\title{Generate generalised logistic parameter bounds for QFA}
\description{
Generates constraints on allowable values of generalised logistic model parameters appropriate for studying microbial growth curves using Quantitative Fitness Analysis (QFA): http://research.ncl.ac.uk/qfa/.  Biologically plausible growth rates limit parameter r and scaling of cell density (between values of 0.0 and 1.0) limit parameters g and K.  The lower limit for K is set to about the best guess for the inoculum density (inocguess). To ensure that fitted generalised logistic models do not deviate drastically from Verhulst's original logistic model, v is constrained to lie between 0.1 and 10.0.  Finally, the inoculum density parameter g is required to be non-negative and less than or equal to the upper allowable limit for K (since real microbial population densities must not decrease from their inoculum density).
}
\usage{
bounds=makeBoundsQFA(inocguess,d,minK,fixG,globalOpt)
}
\arguments{
  \item{inocguess}{Either a numeric value representing an initial guess at the inoculum density (cell density estimate at t=0) or NULL.  If inocguess is set to NULL, then the initial guess for the inoculum density will be estimated from the observed dataset where possible. }
  \item{d}{A data frame containing a numeric Expt.Time column (time since culture inoculation, typically in days) and a numeric Growth column (typically cell density estimates estimated from a timeseries of culture photographs, quantified by software such as Colonyzer: http://research.ncl.ac.uk/colonyzer/).  
  
  Note that, unlike several other curve-fitting functions in this package, this dataframe should be filtered to only contain "good" observations.  For example it should not contain cell density estimates which are classified as being below the limit of detection for the imaging system used.  Ideally, cell densities should not be affected by plate condensation.}
  \item{minK}{Lowest allowable value of K (stationary phase cell density) allowable for a culture which is not classed as dead or missing.}
  \item{fixG}{Boolean specifying whether inoculum density parameter g should be constrained to lie within a narrow region around inocguess.}
  \item{globalOpt}{Boolean specifying whether curve fitting algorithm should search for inoculum density parameter g across a wide search space.}
}
\value{
A named list of vectors (each of length 2) specifying the lower and upper bounds for the space in which each of the generalised logistic model parameters K, r, g and v will be searched for, together with an updated version of the best inoculum density guess (inocguess).
}
\keyword{qfa}