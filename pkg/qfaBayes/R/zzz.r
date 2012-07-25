.First.lib <- function(lib, pkg)
{
library.dynam("qfaBayes", pkg, lib)
require(MASS)
}

.Last.lib <- function(lib, pkg)
{
library.dynam.unload("qfaBayes", pkg, lib)
pos <- match("package:MASS", search())
if(! is.na(pos)){
detach(pos = pos)
}
}


.onLoad <- function(lib, pkg){
library.dynam("qfaBayes", pkg, lib)
require(MASS)
     packageStartupMessage( "Bayesian QFA" )
}
