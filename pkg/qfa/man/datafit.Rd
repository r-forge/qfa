\name{data.fit}
\alias{data.fit}
\title{Fitting generalised logistic model to growth data by least squares (using optim function)}
\description{
Uses the optim function to carry out least squares fit of logistic model to growth data.  
}
\usage{
data.fit(tim,growth,inocguess,xybounds,inits=list(),logTransform=FALSE,verbose=FALSE)
}
\arguments{
\item{tim}{Vector of culture size observation times}
\item{growth}{Vector of observed culture sizes}
\item{inocguess}{Inocluation density estimate (same physical units as growth vector)}
\item{xybounds}{List of upper and lower bounds for model parameters.  NEED TO ADD EXAMPLE OF CONTENTS.}
\item{inits}{List of initial guesses for each parameter.  If empty list, generate initial guesses from mean of bounds}
\item{logTransform}{Boolean specifying whether we shoudl carry out model fit on log scale.}
\item{verbose}{Boolean specifying level of reporting during model fitting}
}

\value{
A named vector of best parameter estimates along with the objective function (square distance between data and model) at the solution.
}
\keyword{qfa}