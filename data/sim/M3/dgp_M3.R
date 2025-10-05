#################################################
########  Import packages and functions  ########  
#################################################

if (Sys.info()["sysname"] == "Windows") {
  Sys.setenv(LANG = "en")
} else{
  Sys.setlocale("LC_MESSAGES", "en_US.utf8")
}

library(spatstat) 
library(parallel)
source("../../../func/func_dgp.R")



#########################################
########  Generate data from M3  ########  
#########################################

A <- c(10, 20, 40) # Sidelength
itr <- 500 # Number of replicates
kappa <- c(0.25, 0.75, 0.2) # Intensities of the parent processes
bandwidth <- c(0.6, 0.3, 1) # Standard deviation of Gaussian kernel
clusters <- rep("Thomas", length(kappa))
xi <- matrix(c(0, -0.7, 0, -0.9, 0, 0, 0.3, 0.1, 0), 3, 3, byrow = T) # Interaction parameters
mc.cores <- if (Sys.info()["sysname"] == "Windows") 1L else parallel::detectCores()
M3.list <- vector("list", length(A))


for (a in seq_along(A)){
  win <- owin(xrange = c(-A[a]/2, A[a]/2), yrange = c(-A[a]/2, A[a]/2))
  # Intensity function
  rho <- list(expression( 3*exp( -2*( (x/A[a])^2 + (y/A[a])^2 ) ) ),
              expression( 2*exp( -2*( (x/A[a])^2 - (y/A[a])^2 ) ) ),
              0.75)
  M3.list[[a]] = rmkppm(rho=rho, kappa=kappa, bandwidth=bandwidth, xi=xi, clusters=clusters,
                        win=win, nsim=itr, nu.ker=NULL, epsth=1e-03,
                        dimyx = c(300, 300), mc.cores=mc.cores)
}



saveRDS(M3.list, paste0("M3_A", paste(A, collapse = "-"), ".rds"))