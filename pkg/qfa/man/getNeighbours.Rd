\name{getNeighbours}
\alias{getNeighbours}
\title{getNeighbours}
\description{
Finds all genes within kb kilobases of genes listed in gvec according to data in sgd object, returned by readSGD.
}
\usage{
getNeighbours(gvec,kb,sgd)
}
\arguments{
  \item{gvec}{Vector of gene names (systematic or standard gene names, or if you're lucky, aliases).}
  \item{kb}{Search range: search for genes whose centres are kb kilobases from centres of genes in gvec}
  \item{sgd}{Data frame representing features.tab file from SGD.  As returned from readSGD.}
}
\value{
Subset of rows from sgd dataframe containing information about genes/features within specified search range of genes in gvec.  FName column contains systematic name of gene/features.  Gene column contains standard name (where available).
}
\keyword{qfa}