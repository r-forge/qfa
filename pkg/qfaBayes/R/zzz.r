".onLoad" <- function(lib, pkg){
  library.dynam(pkg, pkg, lib)
  packageStartupMessage("\nBayesian QFA", 
    "\nStart with demos SHM_simple_C, IHM_simple_C and JHM_simple_C.\n")
}


