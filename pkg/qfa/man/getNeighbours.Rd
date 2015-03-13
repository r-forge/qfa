\name{getNeighbours}
\alias{getNeighbours}
\title{getNeighbours}
\description{
Finds all genes within kb kilobases of genes listed in gvec according to data in sgd object, returned by readSGD.
}
\usage{
getNeighbours(gvec,kb,sgd,geneCoord="Mid",geneDist="Mid",geneDir="All")
}
\arguments{
  \item{gvec}{Vector of gene names (systematic or standard gene names, or if you're lucky, aliases).}
  \item{kb}{Search range: search for genes whose centres are kb kilobases from centres of genes in gvec}
  \item{sgd}{Data frame representing features.tab file from SGD.  As returned from readSGD.}
  \item{geneCoord}{Definition of position of target genes (from gvec).  Options are: "Start", "Stop" and Mid" for the start, end and middle of the target gene respectively.  Default value is "Mid".}
  \item{geneDist}{Gene feature used to define distance between target gene position and other genes on the same chromosome.  Options are: "Start", "Stop", "Mid" (as for geneCoord) and "Min" which is the nearest out of "Start" and "Stop".  Default value is "Mid".}
  \item{geneDir}{Direction from target in which to search for neighbouring genes.  Options are: "Up", "Down" and "All" for just upstream, just downstream and both respectively.  Default value is "All".}
}
\value{
Subset of rows from sgd dataframe containing information about genes/features within specified search range of genes in gvec.  FName column contains systematic name of gene/features.  Gene column contains standard name (where available).
}
\keyword{qfa}