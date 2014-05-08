\name{JHM_main}
\alias{JHM_main}
\title{Runs JHM C code}
\description{
For use with JHM post processed Colonyzer output.
Runs JHM C code
}
\usage{
JHM_main(burn,iters,thin,adaptive_phase,adaptive_period,QFA.IA,QFA.yA,
  QFA.xA,QFA.NoORFA,QFA.NoTIMEA,QFA.IB,QFA.yB,QFA.xB,QFA.NoORFB,
  QFA.NoTIMEB,PRIORS,TUNING)}
\arguments{
  \item{burn}{Burn-in period}
  \item{iters}{Number of iterations in final sample}
  \item{thin}{Amount of thinning (after burn-in)}
  \item{adaptive_phase}{Total length of adaptive phase within burn in period}
  \item{adaptive_period}{Period to calculate adaptive statistic over during adaptive phase}
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
  \item{TUNING}{TUNING}
}
\value{
JHM output
}
\keyword{qfaBayes}
