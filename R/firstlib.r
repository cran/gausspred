.First.lib <- function(lib,pkg)
{
   library.dynam("gausspred",pkg,lib)
   cat("gausspred 1.0.0 loaded\n", 
#        "COPY RIGHT 2008 (c) Longhai Li,http://math.usask.ca/~longhai\n",
       "Type ?begin.gaussianpred for help\n",sep="")
}

