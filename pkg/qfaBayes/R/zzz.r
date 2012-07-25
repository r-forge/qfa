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
     # do whatever needs to be done when the package is loaded
     # some people use it to bombard users with 
     # messages using 
library.dynam("qfaBayes", pkg, lib)
require(MASS)
     packageStartupMessage( "my package is so cool" )
     packageStartupMessage( "so I will print these lines each time you load it")
}
