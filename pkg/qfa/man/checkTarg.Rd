\name{checkTarg}
\alias{checkTarg}
\title{checkTarg}
\description{
Converts gene name to systematic feature name if necessary.  Searches in FName, Alias and Gene columns of a data.frame returned from readSGD for the target gene.  Returns FName (ORF) if appropriate.
}
\usage{
checkTarg(targ,sgd)
}
\arguments{
  \item{targ}{Systematic gene or feature name, standard gene name or alias}
  \item{sgd}{Data frame representing features.tab file from SGD.  As returned from readSGD.}
}
\value{
String representing systematic name of gene or feature.  Returns "MISSING" if target not found in sgd object.
}
\keyword{qfa}