\name{qfa.epiplot}
\alias{qfa.epiplot}
\title{Makes an epistasis plot from the full results of qfa.epi}
\description{
Plots control fitnesses on the x-axis and double mutant fitnesses on the y-axis, with those deemed to be hits (by fdr adjust p-value) coloured.
}
\usage{
qfa.epiplot(results,qthresh,fitratio=FALSE,ref.orf="YOR202W")
}
\arguments{
  \item{results}{The full results (for all ORFs) of qfa.epi.}
  \item{qthresh}{The fdr adjusted cutoff point for determining hits.}
  \item{fitratio}{The ratio of background mutant fitness to wildtype fitness, from the genetic indepdence model.
If FALSE, this is estimated from results using linear regression.}
  \item{ref.orf}{A reference ORF, for which horizontal and vertical lines will be drawn for its fitnesses. (HIS3 default)}
}
\examples{
}
\keyword{qfa}