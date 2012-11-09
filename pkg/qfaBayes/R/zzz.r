".onLoad" <- function(lib, pkg)
{
  library.dynam(pkg, pkg, lib)
     packageStartupMessage("\nBayesian QFA" 
"\nStart with demos SHM_C, IHM_C and JHM_C and then edit these scripts for results with suitable sample size, thinning and burn in.",
"\nC code outside R avalible at https://github.com/jhncl/Model-in-C/tree/master/CCode_FINAL\n")
}


