\name{JHM_main}
\alias{JHM_main}
\title{Runs JHM C code}
\description{
For use with JHM post processed Colonyzer output.
Runs JHM C code
}
\usage{
JHM_main(burn,iters,thin,QFA.IA,QFA.yA,QFA.xA,QFA.NoORFA,QFA.NoTIMEA,QFA.IB,QFA.yB,QFA.xB,QFA.NoORFB,QFA.NoTIMEB,PRIORS)
}
\arguments{
  \item{burn}{Burn in period}
  \item{iters}{Number of iterations in final sample}
  \item{thin}{Amount of thinning}
  \item{QFA.IA}{QFA.IA}
  \item{QFA.yA}{QFA.yA}
  \item{QFA.xA}{QFA.xA}
  \item{QFA.NoORFA}{QFA.NoORFA}
  \item{QFA.NoTIMEA}{QFA.NoTIMEA}
  \item{QFA.IB}{QFA.IB}
  \item{QFA.yB}{QFA.yB}
  \item{QFA.xB}{QFA.xB}
  \item{QFA.NoORFB}{QFA.NoORFB}
  \item{QFA.NoTIMEB}{QFA.NoTIMEB}
  \item{PRIORS}{PRIORS}
}
\value{
JHM output
}
\keyword{qfaBayes}
