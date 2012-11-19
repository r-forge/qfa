\name{SHM_main}
\alias{SHM_main}
\title{Runs SHM C code}
\description{
For use with SHM post processed Colonyzer output.
Runs SHM C code
}
\usage{
SHM_main(burn,iters,thin,CAPL,QFA.I,QFA.y,QFA.x,QFA.NoORF,QFA.NoTIME,PRIORS)
}
\arguments{
  \item{burn}{Burn in period}
  \item{iters}{Number of iterations in final sample}
  \item{thin}{Amount of thinning}
  \item{CAPL}{Cutoff number of orfs}
  \item{QFA.I}{QFA.I}
  \item{QFA.y}{QFA.y}
  \item{QFA.x}{QFA.x}
  \item{QFA.NoORF}{QFA.NoORF}
  \item{QFA.NoTIME}{QFA.NoTIME}
  \item{PRIORS}{PRIORS}
}
\value{
SHM output
}
\keyword{qfaBayes}
