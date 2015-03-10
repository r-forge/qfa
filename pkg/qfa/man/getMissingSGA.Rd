\name{getMissingSGA}
\alias{getMissingSGA}
\title{getMissingSGA}
\description{
For a given list of SGA numbers, finds all of the genotypes which we can classify as being missing in the first stage of each of the SGAs.  We can assume that genotypes which are missing in the first stage of an SGA (where stains are not under marker-specific selection) are simply missing from the library (rather than synthetic lethal interaction with background mutation for example).  Cultures with a specific genotype are classified as being missing if their fitness is less than thresfrac times the median fitness of all genotypes observed in that screen.
}
\usage{
getMissingSGA(QFAs,threshfrac)
}
\arguments{
  \item{QFAs}{List of QFA numbers to examine.}
  \item{threshfrac}{Genotypes whose median fitness is less than thresfrac times the genome-wide median will be classified as dead.}
}
\value{
Subset of rows from fitness reports containing information about genotypes classified as missing.  ORF column contains systematic gene names.  Gene column contains standard gene names.
}
\keyword{qfa}