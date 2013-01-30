\name{visToolDemo}
\alias{visToolDemo}
\title{Interactive fitness plots from Addinal et al. 2011 PLoS Genetics}
\description{
Creates an interactive, searchable version of a selection of fitness plots from Addinall et al. 2011 PLoS Genetics (http://dx.doi.org/10.1371/journal.pgen.1001362).
}
\usage{
visToolDemo()
}
\arguments{
}
\value{
Returns an interactive plot for comparing control and query fitnesses and visualising interaction strengths for data presented in Addinall et al. 2011 PLoS Genetics (http://dx.doi.org/10.1371/journal.pgen.1001362).

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