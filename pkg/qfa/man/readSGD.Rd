\name{readSGD}
\alias{readSGD}
\title{readSGD}
\description{
Parses the SGD_features.tab file, assigning useful column names and handling examples of EOF found in quotes in middle of file.  Adds gene coordinate "Mid" (mean of "Start" and "Stop").  If no filename argument passed to this function, attempts to download latest version from SGD.
}
\usage{
readSGD(fname)
}
\arguments{
  \item{fname}{Path to local copy of SGD_features.tab.  Omit this argument to get latest version from SGD}
}
\value{
Data frame containing gene feature information from SGD_features.tab
}
\keyword{qfa}