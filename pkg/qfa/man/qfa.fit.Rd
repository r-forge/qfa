\name{qfa.fit}
\alias{qfa.fit}
\title{Growth curve modelling}
\description{
Given a series of culture density observations from colonyzer.read, this function will fit the generalised logistic growth model to timecourse observations for all colonies by least squares using either the L-BFGS-B algorithm in R's optim function, or the differential evolution, stochastic global optimisation package DEoptim.  It will also calculate a numerical Area Under Curve (nAUC) fitness measure by integrating under a loess smooothed version of the dataset if there are sufficient observations or under a linear interpolation between observations if observations are too infrequent. 
}
\usage{
qfa.fit(d,inocguess,ORF2gene="ORF2GENE.txt",fmt="\%Y-\%m-\%d_\%H-\%M-\%S",minK=0.025,
detectThresh=0.0005,globalOpt=FALSE,logTransform=FALSE,fixG=TRUE,AUCLim=5,STP=20,
nCores=1,glog=TRUE,...)
}
\arguments{
  \item{d}{The data.frame containing the timecourse data for each colony (returned from colonyzer.read).}
  \item{inocguess}{The best guess for starting density of viable cells in each colony.  This is the g parameter in the generalised logistic model.  Typically, for dilute inoculum 384 format spotted cultures, this value cannot be observed directly by photography.  inocguess should be in the same units as the values in the Growth column in d.  If fixG=TRUE, only values of g within the range 0.9*inocguess and 1.1*inocguess will be assessed during optimisation.  Otherwise values within 0.01*inocguess and 100.0*inocguess will be tried.  Without a sensible independent estimate for inoculum density, the best we can do is to estimate it based on observed data.  Estimating inoculum density will only work well if the inoculum density is high enough to be measurable (e.g. pinned cultures or conc. spotted) and is clearly observed.  Clearly observed means: no condensation on plates immediately after they are placed in incubator for example.  If we are making an independent estimate of inoculum density, then we should also reset the time at which the experiment "begins".  This experiment start time should be the time at which the inoculum density is observed.}
  \item{ORF2gene}{The location of the text file whose first column is of the relevant ORF names and whose second column is of corresponding gene names. If human readable gene names are not important and unique strain identifiers will suffice, set to FALSE.}
  \item{fmt}{The date.time format that the inoculation time (Inoc.Time) and measurement times (Date.Time) are stored in}
  \item{minK}{The minimum value of K above which a strain is said to be alive.  Strains with K optimised to lie below this value will be classified as dead, by setting r to be zero.}
  \item{detectThresh}{The minimum detectable cell density (or Growth value) which reliably identifies the presence of cells.  Cell densities below this value are classified as noise and discarded.}
  \item{globalOpt}{Flag indicating whether qfa.fit should use the slower, but more robust DEoptim global optimisation functions to fit the generalised logistic model to the data, or the quicker optim function.}
  \item{logTransform}{Experimental flag signalling use of different objective function for optimisation.  You should probably ignore this or set it to FALSE}
  \item{fixG}{Flag indicating whether to allow g parameter to vary over a wide or narrow range during optimisation.  fixG=TRUE corresponds to narrow constraints on g.}
  \item{AUCLim}{Numerical AUC (nAUC) is calculated as the integral of an approximation of the growth curve between time 0 and AUCLim}
  \item{STP}{Time to use for \dQuote{Single Time Point} fitness estimate.  Defaults to 20 days (very late in growth curve) which is like carrying capacity.}

  \item{nCores}{Can attempt to split model fitting load across multiple parallel cores.  Experimental, probably best to leave this value set to default (1)}
  \item{glog}{When set to TRUE, carry out generalised (asymmetric) logistic model fit to growth curve data.  When set to FALSE, carry out simpler logistic model fit (as in Addinall et al. 2011)}
  \item{...}{Extra arguments passed to optim} 
}
\value{
R data.frame, similar to that returned by the colonyzer.read function.  The major difference is that instead of a row for every cell density observation for every culture, this object summarises all timecourse density observations for each culture with fitted generalised logistic parameters and nAUC.

\itemize{
\item Barcode - See colonyzer.read
\item Row - See colonyzer.read
\item Col - See colonyzer.read
\item Background - See colonyzer.read
\item Treatment - See colonyzer.read
\item Medium - See colonyzer.read
\item ORF - See colonyzer.read
\item K - Generalised logistic model carrying capacity
\item r - Generalised logistic model rate parameter
\item g - Generalised logistic model inoculum density (referred to in vignette as $g_0$)
\item v - Generalised logistic model shape parameter (set to 1 to recover logistic model)
\item obj - Objective function value at selected optimum
\item t0 - Time of first detectable cell density observation (i.e. above detectThresh)
\item tshift - Shift applied to observation times before fitting logistic model (need to apply same shift before overlaying curve on expt. obs.).  Default is zero (expt. starts at inoculation time specified in experimental description file), but if qfa.fit function is called with inocguess=NULL, then the start of experiment is redefined as the time of the first reliable density observation.
\item nAUC - Numerical Area Under Curve.  This is a model-free fitness estimate.
\item nSTP - Single Time Point fitness.  Cell density at time STP, as estimated with approximating function.  This is a model-free fitness estimate.
\item d0 - Normalised cell density of first observation (be careful about condensation on plates when using this).  Note this is not necessarily the density at t0.
\item Screen.Name - See colonyzer.read
\item Library.Name - See colonyzer.read
\item MasterPlate.Number - See colonyzer.read
\item Timeseries.order - See colonyzer.read
\item Inoc.Time - See colonyzer.read
\item TileX - See colonyzer.read
\item TileY - See colonyzer.read
\item XOffset - See colonyzer.read
\item YOffset - See colonyzer.read
\item Threshold - See colonyzer.read
\item EdgeLength - See colonyzer.read
\item EdgePixels - See colonyzer.read
\item RepQuad - See colonyzer.read
}

}

\keyword{qfa.fit}
