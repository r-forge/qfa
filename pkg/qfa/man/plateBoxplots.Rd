\name{plateBoxplots}
\alias{plateBoxplots}
\title{Plate Boxplots}
\description{
This function generates quality control boxplots summarising fitness patterns between plates (useful in a multi-plate screen, useful for identifying potential issues with media used or timing of inoculation) as well as within plates (useful for identifying Row, Column or Quadrant effects).


The first several pages of output contain boxplots summarising the dependence of fitness on library plate number.  All plate fitnesses are grouped by the groupcol column and by the RepQuad column before boxplots are drawn.

The remaining pages are invididual plate summaries, showing fitness distributions (check for bimodality), fitnesses by row (check for up-down pattern of left-right trend), by column (check for up-down pattern or left-right trend) and by quadrant (check for any single outlier). 

Ideally there should be no trend or pattern in any of these boxplots. 

Boxplots with notches can be used as an initial guide to the significance of any observed differences (boxes with notches that do not overlap may be significantly different).
}
\usage{
plateBoxplots(dataframe,outputfile,fitmax=185,groupcol="Treatment")
}
\arguments{
  \item{dataframe}{Dataframe containing fitnesses to be summarised.  Typically output from qfa.fit function.  Must add a "fit" column before passing to this function.  For example: dataframe$fit=dataframe$MDRMDP}
  \item{outputfile}{Output file name.}
  \item{fitmax}{Upper limit for y-axis in plots (maximum observable fitness in this combination of experiment and fitness definition. }
  \item{groupcol}{String name of column labelling group membership of each plate/barcode.  Typically this will be the Treatment or Medium column (or some combination of both if there are more than one treatment and more than one medium within the dataframe d).  Plates are grouped by unique values of groupcol before being summarised with a boxplot.}
}
\value{
Generates a multi-page .pdf report demonstrating within-plate culture variance for all plates, media, treatments and replicates in a QFA experiment.}
\keyword{qfa}