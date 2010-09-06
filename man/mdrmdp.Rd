\name{mdrmdp}
\alias{mdrmdp}
\title{Fitness Function: MDR*MDP}
\description{
Calculates the phenotype of (maximum doubling rate)*(maximum doubling potential) for a colony. 
}
\usage{
mdrmdp(K,r,g)
}
\arguments{
  \item{K}{The carrying capacity of the colony (final colony size)}
  \item{r}{The colony growth rate parameter}
  \item{g}{The initial colony size}
}
\examples{
# Calculate mdrmdp with final size 6*10^5, rate 5.2 and initial size 116
mdrmdp(6*10^5,5.2,116)
}
\keyword{qfa}
