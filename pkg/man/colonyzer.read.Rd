\name{colonyzer.read}
\alias{colonyzer.read}
\title{Reading of raw cell density timecourse data from Colonyzer output files}
\description{
Reads in and binds together all of the Colonyzer output files in a directory, so they are checked and ready for bayesian or likelihood inference.
Colonyzer is open source image analysis software for quantifying cell densities on agar plates: http://research.ncl.ac.uk/colonyzer/
}
\usage{
colonyzer.read(path=".",files=c(),experiment="ExptDescription.txt",ORF2gene="",libraries="LibraryDescriptions.csv",background="")
}
\arguments{
  \item{path}{The path to the folder containing the Colonyzer .dat files to be read: working directory by default.}
  \item{files}{Character vector giving locations of Colonyzer .dat files to be read - overrides path}
  \item{experiment}{A text file describing the inoculation times, library and plate number for unique plates.  MUST ADD DETAILED FORMAT SPEC. Taken relative to path if specified.}
  \item{ORF2gene}{Path to a file linking unique ORF numbers with human readable gene names (e.g. standard names from SGD).}
  \item{libraries}{CSV file describing the array strains present in each row-column coordinate of each plate in a series of libararies. MUST ADD DETAILED FORMAT SPEC.}
  \item{background}{Experiment identifier (e.g. query mutation)}
}
\value{
An R data.frame where each row corresponds to a single observation on a single colony, with the value of the growth measurement in 'Growth', and the date and time of the measurement in 'Date.Time'. Other information about the observation is stored in the other columns.
}
\examples{
}
\keyword{qfa}