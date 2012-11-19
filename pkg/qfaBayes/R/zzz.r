".onLoad" <- function(lib, pkg)
{
  library.dynam(pkg, pkg, lib)
     packageStartupMessage("\nBayesian QFA", 
"\nStart with demos SHM_simple_C, SHM_C, IHM_C and JHM_C and then edit these scripts for results with suitable sample size, thinning and burn in.",
"\nThe demos within this package make heavy use of local storage. It is therefore recommended you move to a directory with large storage capability", 
"\nTo change tuning paramaters use C code outside R avalible at https://github.com/jhncl/Model-in-C/tree/master/CCode_FINAL\n")
}


