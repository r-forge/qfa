\name{qfa.fit}
\alias{qfa.fit}
\title{Growth curve modelling}
\description{
Given a series of culture density observations from colonyzer.read, this function will fit the generalised logistic growth model to timecourse observations for all colonies by least squares using either the L-BFGS-B algorithm in R's optim function, or the differential evolution, stochastic global optimisation package DEoptim.  It will also calculate a numerical Area Under Curve (nAUC) fitness measure by integrating under a loess smooothed version of the dataset if there are sufficient observations or under a linear interpolation between observations if observations are too infrequent. 
}
\usage{
qfa.fit(d,inocguess,ORF2gene="ORF2GENE.txt",fmt="\%Y-\%m-\%d_\%H-\%M-\%S",minK=0.025,
detectThresh=0.0005,globalOpt=FALSE,logTransform=FALSE,fixG=TRUE,AUCLim=5,STP=20,
nCores=1,glog=TRUE,modelFit=TRUE...)
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
  \item{glog}{Boolean (TRUE or FALSE) specifying whether to carry out generalised (asymmetric) logistic model fit to growth curve data.  When set to FALSE, carry out simpler logistic model fit (as in Addinall et al. 2011)}
  \item{modelFit}{Boolean (TRUE or FALSE) specifying whether to carry out any model fitting at all.  When set to FALSE, only numerical fitness estimates such as nr, nMDP, nAUC are generated}
  \item{checkSlow}{Boolean (TRUE or FALSE) specifying whether to re-optimise curve-fitting for slow-growing strains.  If TRUE, slow-growing or dead strains are identified heuristically and a second round of curve fitting using global (but slower) optimisation is carried out.  Heuristic identification of slow-growing strains is currently experimental, it seems we have over-tuned these to datasets we capture at Newcastle.  If you notice a banding pattern in your MDR or r fitnesses, please set checkSlow to FALSE.}
  \item{...}{Extra arguments passed to optim} 
}
\value{
R data.frame, similar to that returned by the colonyzer.read function.  The major difference is that instead of a row for every cell density observation for every culture, this object summarises all timecourse density observations for each culture with fitted generalised logistic parameters and numerical fitness estimates.

\itemize{
\item Barcode - Unique plate identifier
\item Row - Row number (counting from top of image) of culture in rectangular gridded array
\item Col - Column number (counting from left of image) of culture in rectangular gridded array
\item ScreenID - Unique identifier for this QFA screen
\item Treatment  - Conditions applied externally to plates (e.g. temperature(s) at which cultures were grown, UV irradiation applied, etc.)
\item Medium - Nutrients/drugs in plate agar
\item ORF - Systematic, unique identifier for genotype in this position in arrayed library
\item Screen.Name - Name of screen (identifies biological repeats, and experiment)
\item Library.Name	- Name of library, specifying particular culture location
\item MasterPlate Number - Library plate identifier
\item Timeseries order - Sequential photograph number
\item Inoc.Time - User specified date and time of inoculation (specified in ExptDescription.txt file)
\item TileX - Culture tile width (pixels)
\item TileY - Culture tile height (pixels)
\item XOffset - x-coordinate of top left corner of rectangular tile bounding culture (number of pixels from left of image)
\item YOffset - y-coordinate of top left corner of rectangular tile bounding culture (number of pixels from top of image)
\item Threshold - Global pixel intensity threshold used for image segmentation (after lighting correction)
\item EdgeLength - Number of culture pixels classified as being microcolony edge pixels (useful for classifying contaminants in cultures grown from dilute inoculum)
\item EdgePixels - Number of pixels classified as culture on edge of square tile
\item RepQuad - Integer identifying which of the quadrants of a 1536 plate were used to inoculate the current 384 plate (set equal to 1 for all cultures for 1536 format for example)
\item K - Generalised logistic model carrying capacity
\item r - Generalised logistic model rate parameter
\item g - Generalised logistic model inoculum density (referred to in vignette as $g_0$)
\item v - Generalised logistic model shape parameter (set to 1 to recover logistic model)
\item objval - Objective function value at selected optimum
\item tshift - Shift applied to observation times before fitting logistic model (need to apply same shift before overlaying curve on expt. obs.).  Default is zero (expt. starts at inoculation time specified in experimental description file), but if qfa.fit function is called with inocguess=NULL, then the start of experiment is redefined as the time of the first reliable density observation.
\item t0 - Time of first detectable cell density observation (i.e. above detectThresh)
\item d0 - Normalised cell density of first observation (be careful about condensation on plates when using this).  Note this is not necessarily the density at t0.
\item nAUC - Numerical Area Under Curve.  This is a model-free fitness estimate.
\item nSTP - Single Time Point fitness.  Cell density at time STP, as estimated with approximating function.  This is a model-free fitness estimate.
\item nr - Numerical estimate of intrinsic growth rate.  Growth rate estimated by fitting smoothing function to log of data, calculating numerical slope estimate across range of data and selecting the maximum estimate (should occur during exponential phase).
\item nr_t - Time at which maximum slope of log observations occurs
\item maxslp - Numerical estimate of maximum slope of growth curve.  Slope estimated by fitting smoothing function to untransformed data and calculating numerical slope estimate of smoothed version of data and selecting the maximum estimate (should occur approximately half way through growth).  This fitness measure will be affected by both rate of growth and final colony size.  Final colony size is expected to be strongly affected by competition between cultures.
\item maxslp_t - Time at which maximum slope of observations occurs
\item Client - Client for whom screen was carried out
\item ExptDate - A representative/approximate date for the experiment (note that genome-wide QFA screens typically take weeks to complete)
\item User - Person who actually carried out screen
\item PI - Principal investigator leading project that screen is part of
\item Condition - The most important defining characteristic of screen, as specified by user (e.g. the temperature screen was carried out at if screen is part of multi-temperature set of screens, or the query mutation if part of a set of screens comparing query mutations, or the drugs present in the medium if part of a set of drug screens)
\item Inoc - Qualitative identifier of inoculation type (e.g. "DIL" for dilute inoculum, "CONC" for concentrated).  Used to distinguish between experiments carried out with different methods of inoculation.
\item Gene - Identifier for genotype at a particular location on an agar plate.  Typically prefer unambiguous, systematic gene names here.
\item TrtMed - Combination of treatment and medium identifiers, specifying the environment in which the cells have grown
}

}

\keyword{qfa.fit}
