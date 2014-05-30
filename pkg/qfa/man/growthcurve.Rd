\name{growthcurve}
\alias{growthcurve}
\title{Wrapper function for generating QFA generalised logistic model fits to observed timecourse data}
\description{
Checks that there are enough observations above the detection threshold to fit the model (if not, returns parameters equivalent to a dead colony).  Then, generates bounds on parameter values appropriate for QFA (including Colonyzer image analysis).  Carries out least squares fitting of model to data, switching between global optimisation (de.fit) or local optimisation from a sensible starting point (data.fit).  If the experiment has been stopped before the culture reaches stationary phase, a second round of optimisation is carried out, with a more tightly bound K parameter, in order to improve the estimation of r.  Regardless of whether local or global optimisation were selected, if a suspicious combination of high r and low K is selected for a culture, global optimisation in a search space corresponding to a sick culture is attempted.  If a lower squared error is generated, the corresponding parameter values are proposed for the current culture.  As one final check, if a completely dead culture (parameter values K=g, r=0) fits the data better than the best found so far, that set of parameters are returned for the current culture.
}
\usage{
pars=growthcurve(obsdat,inocguess,fixG,globalOpt,detectThresh,minK,logTransform)
}
\arguments{
  \item{obsdat}{A data frame containing a numeric Expt.Time column (time since culture inoculation, typically in days) and a numeric Growth column (typically cell density estimates estimated from a timeseries of culture photographs, quantified by software such as Colonyzer: http://research.ncl.ac.uk/colonyzer/).} 
  \item{inocguess}{Either a numeric value representing an initial guess at the inoculum density (cell density estimate at t=0) or NULL.  If inocguess is set to NULL, then the initial guess for the inoculum density will be estimated from the observed dataset where possible.}
  \item{fixG}{Boolean specifying whether inoculum density parameter g should be constrained to lie within a narrow region around inocguess.  Default TRUE.}
  \item{globalOpt}{Boolean specifying whether curve fitting algorithm should search for inoculum density parameter g across a wide search space. Default FALSE.}
  \item{detectThresh}{Minimimum cell density that is reliably detectable by this combination of plate, imager and image analysis software.  Default 0.}
   \item{minK}{Lowest allowable value of K (stationary phase cell density) allowable for a culture which is not classed as dead or missing.  Default 0.}
\value{
A list of generalised logistic model parameters (K, r, g, v) which give the best fit to the observed data in obsdat, together with objval, the objective value (squared error).
}
\keyword{qfa}