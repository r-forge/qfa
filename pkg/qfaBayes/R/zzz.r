.onLoad <- function(lib, pkg){
library.dynam("qfaBayes", pkg, lib)
     packageStartupMessage( "Bayesian QFA" )
}
