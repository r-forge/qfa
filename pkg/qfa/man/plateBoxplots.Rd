\name{plateBoxplots}
\alias{plateBoxplots}
\title{Plate Boxplots}
\description{
This function generates quality control boxplots summarising fitness patterns between plates (useful in a multi-plate screen, for identifying potential issues with media used or timing of inoculation) as well as within plates (useful for identifying Row, Column or Quadrant effects which might arise due to pinning or spotting issues).

The first several pages of output contain boxplots summarising the dependence of fitness on library plate number.  All plate fitnesses are grouped by the groupcol column and by the RepQuad column before boxplots are drawn.  You should examine these pages for evidence of any outlier plates or groups of plates (suggests there is a batch effect in the media used) or for any trend with plate number (suggests there is an issue with inoculation timing, or the time between inoculation and moving the plate into the incubator). 

The remaining pages are invididual plate summaries, showing fitness distributions (check for bimodality), fitnesses by row (check for up-down pattern or a left-right trend), by column (check for up-down pattern or a left-right trend) and by quadrant (check for any single outlier quadrant).  Row or column trends might indicate a lighting issue, a problem with agar thickness or temperature gradients across the plates for example.  Outliers in the quadrant plots indicate a problem with spotting (different numbers of cells being inoculated from each of the source plates (typically 4)).  Bimodality in the fitness distribution on a plate is also an indication of spotting problems. 

Ideally there should be no trend or pattern in any of these boxplots. 

Boxplots with notches can be used as an initial indication of whether any observed differences are significant: pairs of boxes with notches that do not overlap may be significantly different.  For example if the notch for quadrant 4 of a plate does not overlap with the notch for quadrant 1 of the same plate, that would suggest that the fitnesses of spots in quadrant 1 are significantly different from those in quadrant 4.
}
\usage{
plateBoxplots(dataframe,outputfile,fitmax=185,groupcol="Treatment", fdef="fit")
}
\arguments{
  \item{dataframe}{Dataframe containing fitnesses to be summarised.  Typically output from qfa.fit function.  Must add a "fit" column before passing to this function.  For example: dataframe$fit=dataframe$MDRMDP}
  \item{outputfile}{Output file name.}
  \item{fitmax}{Upper limit for y-axis in plots (maximum observable fitness in this combination of experiment and fitness definition. }
  \item{groupcol}{String name of column labelling group membership of each plate/barcode.  Typically this will be the Treatment or Medium column (or some combination of both if there are more than one treatment and more than one medium within the dataframe d).  Plates are grouped by unique values of groupcol before being summarised with a boxplot.}
  \item{fdef}{Fitness definition chosen (e.g. "r","nAUC","MDR","MDRMDP").}
}
\value{
Generates a multi-page .pdf report demonstrating within-plate and between-plate culture variance for all plates, media, treatments and replicates in a QFA experiment.}
\keyword{qfa}