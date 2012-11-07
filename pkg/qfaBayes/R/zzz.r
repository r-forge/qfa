.onLoad <- function(lib, pkg){
library.dynam("qfaBayes", pkg, lib)
     packageStartupMessage( "
Bayesian QFA 
(independent C code avalible at https://github.com/jhncl/Model-in-C/tree/master/CCode_FINAL)" )
}
