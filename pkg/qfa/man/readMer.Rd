\name{readMer}
\alias{readMer}
\title{readMer}
\description{
Parse the .mer version of the Lydall lab QFA_database.  .mer files are text dumps of the filemaker database used by the lab to capture experimental metadata.  Ensures that SGA.Number column contains an integer.
}
\usage{
readMer(fname)
}
\arguments{
  \item{fname}{Path to .mer file describing QFA_Database.}
}
\value{
Subset of rows from sgd dataframe containing information about genes/features within specified search range of genes in gvec.  FName column contains systematic name of gene/features.  Gene column contains standard name (where available).
}
\keyword{qfa}