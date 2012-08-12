\name{qfa.model}
\alias{qfa.model}
\title{Initialize JAGS model}
\description{
Builds and returns a hierarchical Bayesian model in rjags. The model is based on the logistic growth model for the ORFs' timecourses, and the ORFs should come from the same genetic background and treatment group.
}
\usage{
qfa.model(jagd,Kloc,gloc,Ksd="(2/3)*Kloc",g.range=10,...)
}
\arguments{
  \item{jagd}{The data object, output from qfa.data, containing all the timecourse data for each ORF}
  \item{Kloc}{What you expect a typical final colony size to be in this background & treatment; important for sensible location parameters of hyperpriors.}
  \item{gloc}{Best guess at initial colony size}
  \item{Ksd}{A positive number, or an expression involving Kloc, giving the standard error associated with your guess at typical final colony size, Kloc. The more certain about Kloc, the smaller Ksd should be in relation to Kloc. By default, Ksd is 2/3 of Kloc, which is somewhat uninformative about typical final colony size.}
  \item{g.range}{A percentage which determines the range over which initial colony size is allowed to vary: [(1-g.range/100)*gloc,(1+g.range/100)*gloc].} 
  \item{...}{Further arguments passed to \code{\link{jags.model}}}
}
\value{
Returns a hierarchical Bayesian QFA model in rjags. Kloc and Ksd are the mean and sd of the gamma prior for mean ORF mean final colony size across all the ORFs in this model. The inital colony size, g, has a uniform prior bound by [(1-g.range/100)*gloc,(1+g.range/100)*gloc].
}
\examples{
}
\keyword{qfa}
