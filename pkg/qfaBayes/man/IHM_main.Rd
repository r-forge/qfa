\name{IHM_main}
\alias{IHM_main}
\title{Runs IHM C code}
\description{
For use with control and query SHM output.
Runs IHM C code
}
\usage{
IHM_main(burn,iters,thin,QFA.IA,QFA.yA,QFA.NoORFA,QFA.IB,QFA.yB,QFA.NoORFB,PRIORS)
}
\arguments{
  \item{burn}{Burn in period}
  \item{iters}{Number of iterations in final sample}
  \item{thin}{Amount of thinning}
  \item{QFA.IA}{QFA.IA}
  \item{QFA.yA}{QFA.yA}
  \item{QFA.NoORFA}{QFA.NoORFA}
  \item{QFA.IB}{QFA.IB}
  \item{QFA.yB}{QFA.yB}
  \item{QFA.NoORFB}{QFA.NoORFB}
  \item{PRIORS}{PRIORS}
}
\value{
IHM output
}
\keyword{qfaBayes}
