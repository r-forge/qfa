\name{qfa.simulate}
\alias{qfa.simulate}
\title{Simulate the posterior of a qfa.model}
\description{
Does initial burn in and update of qfa.model, followed by sampling of the posterior.
}
\usage{
qfa.simulate(model,initial.update=5*10^4,samples=5*10^4,thins=10,orfs=c(),
params=c("k1","r1","g","k11","k12","k21","k22",
"r11","r12","r21","r22","s1","s2"),...)
}
\arguments{
  \item{model}{qfa.model you wish to simulate the posterior of}
  \item{initial.update}{The number of updates that should be done before sampling.}
  \item{samples}{The size of the sample to be drawn once the initial updating is finished.}
  \item{thins}{The thinning to be applied to the samples.}
  \item{orfs}{The names of the ORFs in the model, in the correct order, usually passed from $ORF in previous qfa.data output.}
  \item{params}{The names of the parameters that should be sampled; by default, this samples the ORF level parameters relevant to ORF level growth dynamics (\dQuote{k1},\dQuote{r1},\dQuote{g}), as well as the highest level parameters of the model to aid in autocorrelation tests.}
  \item{...}{Additional arguments passed to \code{\link{coda.samples}}}
}
\value{
An R list with an element for each parameter giving the posterior sample. If the parameter is present for every orf (i.e. (\dQuote{k1},\dQuote{r1})) the element in the posterior list will be a matrix with a column for each ORF. This can push the limits of R memory for large posterior samples for a lot of ORFs.
}
\examples{
}
\keyword{qfa}
