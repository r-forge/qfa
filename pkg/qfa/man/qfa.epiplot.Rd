\name{qfa.epiplot}
\alias{qfa.epiplot}
\title{Makes an epistasis plot from the full results of qfa.epi}
\description{
Plots control fitnesses on the x-axis and double mutant fitnesses on the y-axis, with those deemed to be hits (by fdr adjust p-value) coloured.
}
\usage{
qfa.epiplot(results,qthresh,fitratio=FALSE,ref.orf="YOR202W",xxlab="Control Fitness",
yylab="Query Fitness",mmain="Epistasis Plot",fmax=0)
}
\arguments{
  \item{results}{The results of interaction analysis returned by the qfa.epi function.}
  \item{qthresh}{The fdr adjusted cutoff point for determining hits.}
  \item{fitratio}{The ratio of background mutant fitness to wildtype fitness, from the genetic indepdence model.
If FALSE, this is estimated from results using linear regression.}
  \item{ref.orf}{ORF for a reference strain (typically wild-type or a surrogate), whose fitness will be marked on the control and query axes of the interaction plot (horizontal and vertical blue lines). HIS3 is the default strain.}
  \item{xxlab}{x axis label}
  \item{yylab}{y axis label}
  \item{mmain}{Plot label}
  \item{fmax}{Maximum fitness range for both x-axis (control axis) and y-axis (query axis).  If 0, axis ranges are automatically chosen to include all data points.}
}
\examples{
}
\keyword{qfa}