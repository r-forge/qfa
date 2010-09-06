\name{qfa.epi}
\alias{qfa.epi}
\title{Finds genetic interaction strengths and p-values}
\description{
Fits a genetic independence model between control strains and double mutant strains, either using rjags and a Bayesian linear regression model, or lm and maximum likelihood. For each ORF, the probability that it is a false discovery of a suppressor or enhancer is calculated. These probabilities are then fdr corrected and returned along with genetic interaction scores.
}
\usage{
qfa.epi(double,control,qthresh,rjags=FALSE,orfdict="ORF2GENE.txt",
GISthresh=0.5,edgestrip=TRUE,plot=TRUE,modcheck=TRUE,fitfunct=mdrmdp)
}
\arguments{
  \item{double}{Either a qfa.posterior or the results of qfa.fit for the double mutants}
  \item{control}{Either a qfa.posterior or the results of qfa.fit for the control strains}
  \item{qthresh}{The FDR corrected cut off}
  \item{rjags}{If TRUE, fits genetic independence model using rjags; if FALSE, fits genetic indepdence model using lm.}
  \item{orfdict}{If using Bayesian posteriors, or if gene names not present in double or control, location of file giving a column of ORFs first and a column of corresponding gene names second - so gene names can be plotted}
  \item{GISthresh}{When returning interaction hitlists, this variable determines the cutoff for strength of genetic interaction.}
  \item{edgestrip}{If TRUE, this does not include colonies from the edges of the plates in calculations.}
  \item{plot}{If TRUE, then an epistasis plot is made}
  \item{modcheck}{If TRUE and rjags is FALSE, then diagnostic residual plots are output to \dQuote{ModelCheck.pdf}}
  \item{fitfunct}{The name of a fitness function whose arguments are, in order, (K,r,g) (carrying capacity, rate and initial size of logistic growth model).}
}
\value{
Returns a R list, with elements Results - which displays the p-values and genetic interaction scores for all ORFs in an R data.frame - and Suppressors and Enhancers, which are R data.frames with only those ORFs of that type (enhancer/suppressor), with abs(GIS)>GISthresh and fdr adjust p-value<qthresh.
}
\examples{
}
\keyword{qfa}