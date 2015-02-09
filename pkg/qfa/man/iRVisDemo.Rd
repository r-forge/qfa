\name{iRVisDemo}
\alias{iRVisDemo}
\title{Interactive fitness plots from Addinal et al. 2011 PLoS Genetics}
\description{
Creates an interactive, searchable version of a selection of fitness plots from Addinall et al. 2011 PLoS Genetics (http://dx.doi.org/10.1371/journal.pgen.1001362).
}
\usage{
iRVisDemo(groupFun=buildComplexes, fitmax=0)
}
\arguments{
\item{groupFun}{Specify the name of a function to build a data frame containing related groups of genes to be highlighted together.  Current options are buildBenschop (default) and buildGO.  buildBenschop allows highlighting of functionally related genes from Benschop et al. (2010): http://dx.doi.org/10.1016/j.molcel.2010.06.002.  buildGO allows highlighting of functionally related genes from the Gene Ontology database: http://www.geneontology.org/}
\item{fitmax}{Upper limit for fitness values on axes.  Specifying a value here is useful if you need to make multiple plots.  Fixing the limits of your axes will allow fitness plots to be compared more easily.  This is partiuclarly important when making a multiple panel figure for a manuscript, for example.}
}
\value{
Returns an interactive plot for comparing control and query fitnesses and visualising interaction strengths for data presented in Addinall et al. 2011 PLoS Genetics (http://dx.doi.org/10.1371/journal.pgen.1001362).

More detailed description and instructions can be found here: http://qfa.r-forge.r-project.org/visTool

Controls: Windows mouse
~~~~~~~~~~~~~~~
Left click: Highlight gene/Rotate text position
Right click: Open SGD page for gene (alternatively, press 'w' on keyboard)
Middle click: Remove last gene (alternatively, press 'd' on keyboard)

Controls: OSX mouse
~~~~~~~~~~~~~~~
Click: Highlight gene/Rotate text position

Controls: Keyboard
~~~~~~~~~~~~~~~
Left/Right arrow: Switch to next/previous fitness plot (comparison between different pair of experiments).
Up/Down arrow: Change group of genes currently highlighted (in purple).
b: Print experimental metadata report to console window
c: Clear highlighting from currently selected genes.
d: Remove highlighting from last gene highlighted.
i: In rank plot mode, toggle between x-axis rank order.
k: In rank plot mode, toggle between fitness ratio, GIS, query fitness and control fitness on y-axis.
l: Toggle plot style between fitness plot (default) and rank plot mode.
p: Save current plot as vector graphic to QFAVisualisation.ps.  Other filetypes can also be generated - 'n': .png and 'm': .pdf
q: Quit tool and print gene names currently selected to console window.
r: Enter zoom mode.  Click on top left and bottom right of area to inspect (experimental, may have to quit to zoom out...).
s: Highlight genes encircled using select tool ('z').
t: Toggle between colours (red/green) indicating positive & negative interaction & no colouring.
u: Add new group of genes to list of highlightable groups.
v: Print list of currently selected genes to console and copy standard gene names to clipboard
w: Open SGD web-page for last gene highlighted.
x: Print list of currently selected genes to console and copy systematic gene names to clipboard
z: Toggle select tool on/off. Select genes for highlighting by drawing a polygon on plot.  Press 's' to add selection when finished.
}
\examples{
\dontrun{visToolDemo()}
}
\keyword{qfa}