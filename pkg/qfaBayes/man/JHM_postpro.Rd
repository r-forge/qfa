\name{JHM_postpro}
\alias{JHM_postpro}
\title{Proforms post processing for JHM}
\description{
For use with Colonyzer output.
Proforms post processing for JHM.
}
\usage{
JHM_postpro(a,TreatA,Screen_a,MPlate_a,b,TreatB,Screen_b,MPlate_b)
}
\arguments{
\item{a}{Control Colonyzer output}
  \item{TreatA}{Control list of treatments}
  \item{Screen_a}{Control list of Screens}
  \item{MPlate_a}{Control list of Master Plates}
  \item{remove_row_a}{plate rows to remove for control}
  \item{remove_col_a}{plate columns to remove for control}
  \item{b}{Query Colonyzer output}
  \item{TreatB}{Query list of treatments}
  \item{Screen_b}{Query list of Screens}
  \item{MPlate_b}{Query list of Master Plates}
  \item{remove_row_b}{plate rows to remove for query}
  \item{remove_col_b}{plate columns to remove for query}
}
\value{
Post processed Colonyzer output for JHM.
}
\keyword{qfaBayes}
