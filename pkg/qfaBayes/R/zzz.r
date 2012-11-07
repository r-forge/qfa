".onLoad" <- function(lib, pkg)
{
  library.dynam(pkg, pkg, lib)
     packageStartupMessage("Bayesian QFA")
}


