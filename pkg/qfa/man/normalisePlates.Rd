\name{normalisePlates}
\alias{normalisePlates}
\title{Normalising culture fitness by plate}
\description{
Sometimes estimated culture fitnesses vary systematically depending on the plate on which they are inoculated.  Agar in individual plates could come from different batches, and therefore have slightly different levels of nutrients or water.  Plates could be inoculated at different times, and stored at slightly different temperatures for example.  Depending on inoculation method, inoculation time specified may be less accurate for individual plates.  Any of these issues could effect simulated fitness slightly.  This function allows us to normalise culture fitnessses across plates to eliminate such effects.  It should only really be used for small differences.  In order to preserve real signal, the experimental sources of larger differences should be corrected before analysis instead of by normalising them away.

Starting with a data frame describing the output from the qfa.fit function (with optional added columns from the makeFitness function) normalisePlates finds all unique groups in that data frame, calculates a median value from the indicated column for all plates in a given group and then normalises the fitnesses of each culture on each plate so that the median fitness on each plate is equal to the median fitness for all plates in the group.  The function returns a vector which can be added to the original data frame or used to over-write the original, raw data.

Typically a "group" will be a description of treatment (e.g. temperature) or of growth medium (e.g. drug added to solid agar).
}
\usage{
normalisePlates(d,column)}
\arguments{
  \item{d}{Dataframe (output from qfa.fit) for normalisation.}
  \item{column}{String name of column to normalise (typically the fitness measure of interest).  See qfa.fit and makeFitness help files for descriptions of available culture fitness measures.}
  \item{groupcol}{String name of column labelling group membership of each plate/barcode.  Typically this will be the Treatment or Medium column (or some combination of both if there are more than one treatment and more than one medium within the dataframe d.}
}
\keyword{normalisePlates}