\name{qfa.fit}
\alias{qfa.fit}
\title{Optimization procedure for logistic growth modelling}
\description{
Given timecourse data from rod.read, or alternatively a data.frame including Barcode, Row, Col, Growth, Date.Time, Inoc.Time, Tile.Dimensions.X, Tile.Dimensions.Y, where each row corresponds to a single measurement of a single colony at that Date.Time, this function will fit the logistic growth model to all colonies using the L-BFGS-B algorithm via R's optim function.
}
\usage{
qfa.fit(d,inocguess,ORF2gene="ORF2GENE.txt",
fmt="\%Y-\%m-\%d_\%H-\%M-\%S",edgestrip=TRUE,...)
}
\arguments{
  \item{d}{The data.frame containing the timecourse data for each colony.}
  \item{inocguess}{The best guess at initial colony size: optimization on this parameter is restricted to within 10 percent of inocguess.}
  \item{ORF2gene}{The location of the text file whose first column is of the relevant ORF names and whose second column is of corresponding gene names. If gene names are not wanted, set to FALSE.}
  \item{fmt}{The date.time format that the inoculation time (Inoc.Time) and measurement times (Date.Time) are stored in}
  \item{edgestrip}{If TRUE, the model is not fitted to colonies on the edge of plates (max and min of Row and Col).}
  \item{...}{Extra arguments passed to optim} 
}
\value{
R data.frame with a row for each colony including the fitted model parameters for that colony.
}
\keyword{qfa}
