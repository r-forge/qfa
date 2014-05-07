\name{JHM_postpro}
\alias{JHM_postpro}
\title{Performs post processing for JHM}
\description{
For use with Colonyzer output.
Performs post-processing for JHM.
}
\usage{
JHM_postpro(a,Treatment_a,Screen_a,MPlate_a,remove_row_a,remove_col_a,
  b,Treatment_b,Screen_b,MPlate_b,remove_row_b,remove_col_b)
}
\arguments{
\item{a}{Control Colonyzer output}
  \item{Treatment_a}{Control list of treatments}
  \item{Screen_a}{Control list of Screens}
  \item{MPlate_a}{Control list of Master Plates}
  \item{remove_row_a}{plate rows to remove for control}
  \item{remove_col_a}{plate columns to remove for control}
  \item{b}{Query Colonyzer output}
  \item{Treatment_b}{Query list of treatments}
  \item{Screen_b}{Query list of Screens}
  \item{MPlate_b}{Query list of Master Plates}
  \item{remove_row_b}{plate rows to remove for query}
  \item{remove_col_b}{plate columns to remove for query}
}
\value{
Post-processed Colonyzer output for JHM.
}
\keyword{qfaBayes}
