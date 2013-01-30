\name{visTool}
\alias{visTool}
\title{Interactive fitness plots}
\description{
Creates an interactive, searchable version of a selection of fitness plots
}
\usage{
visTool(complexes,orf2gene,GISfiles)
}
\arguments{
  \item{complexes}{Data frame describing the names and members of functionally related complexes.  Largely taken from Benschopp et al. 2010 (http://dx.doi.org/10.1016/j.molcel.2010.06.002).  Contains at least three columns: "X..Complex.Name", "CompList" and "comments" describing the name of a protein complex, the systematic gene names (Y numbers) for members of the complex and any comments.}
  \item{orf2gene}{Data frame decribing a one-to-one relationship between standard gene names and systematic gene names.  Contains two columns: "ORF" and "Gene".  This association is used to convert between standard and systematic gene names during searching.  Note that only one standard gene name is associated with each systematic name, and so, during any search for standard gene names, if you happen to choose a version which is not included in this list, your target gene will not be found.}
  \item{GISfiles}{A list of objects output from the qfa package function report.epi, each of which compares the fitnesses of a library of strains grown as part of a control experiment and as part of a query experiment.}
}
\value{
Returns an interactive plot for comparing control and query fitnesses and visualising interaction strengths for each of the report.epi objects in the GISfiles argument.

Instructions for use:

Windows mouse
~~~~~~~~~~~~~~~
L click: Highlight gene/Rotate text position
R click: SGD (or press 'w' on keyboard)
M click: Remove last gene (or press 'd' on keyboard)

Mac mouse
~~~~~~~~~~~~~~~
Click: Highlight gene/Rotate text position

Keyboard
~~~~~~~~~~~~~~~
Left/Right arrow: change plot
Up/Down arrow: change functional complex highlighted
u: add new genes to list of functional complexes
z: select tool (toggle on and off)
s: add selection 
c: clear selection
w: open last gene highlighted in SGD
d: unhighlight last gene highlighted
t: toggle colours indicating positive and negative interaction
r: begin zoom (now click on top left and bottom right of new zoomed plot)
p: print current plot to QFAVisualisation.pdf
q: quit
}
\examples{
\dontrun{visToolDemo()}
}
\keyword{qfa}