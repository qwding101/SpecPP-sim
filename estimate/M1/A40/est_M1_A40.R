if (Sys.info()["sysname"] == "Windows") {
  Sys.setenv(LANG = "en")
} else{
  Sys.setlocale("LC_MESSAGES", "en_US.utf8")
}

##############################################
########## Import data and functions #########
##############################################

library(spatstat)
library(parallel)
library(doParallel)
source("../../../func/func_sim.R", chdir = T)

### Simulation condition
model = "M1"
sidelength = A1 = A2 = 40 # Correspond to pp.data[[3]]
kappa = c(0.25, 0.75, 0.2) # Parent intensity
sigma = c(0.6, 0.3, 1) # Standard deviation in normal kernel
nu = NULL
xi = matrix(c(0, 0.7, 0, 0.9, 0, 0, 0.3, 0.1, 0), 3, 3, byrow = T) # Interaction parameters
kernels = rep("Thomas", length(kappa))

df.itr = data.frame(c(1,2,1), c(1,2,2))
pp.data = readRDS("../../../data/sim/M1/M1_A10-20-40.rds")
pp.data = pp.data[[3]] # Only select A = 40

# Remove the 3rd process
for (i in seq_along(pp.data)){
  pp.data[[i]] =  pp.data[[i]][marks(pp.data[[i]]) != 3][drop=3]
}
pp.data = as.solist(pp.data)






##########################################################
### Raw periodogram estimator with estimated intensity ###
##########################################################

result.11 = parallel.period.fun(i=1, j=1, pp.data, inten.fun.sep="~1", endpt = 1.5)
result.22 = parallel.period.fun(i=2, j=2, pp.data, inten.fun.sep="~1", endpt = 1.5)
result.12 = Map(function(x,y) Re(x*Conj(y)), result.11$DFT.list, result.22$DFT.list)
period.est.I.hat = list(I11 = result.11$period.list, I22 = result.22$period.list, I12 = result.12)

saveRDS(period.est.I.hat, paste0(model, "_A", sidelength,"_period.rds"))





#######################################################################
###### Kernel spectral density estimator with optimal bandwidth #######
#######################################################################

### Bandwidth
band.fixed = rep(area(pp.data[[1]])^(-1/6), length(pp.data))
band.fixed[1]

### Estimation
period.est.smooth.fixed = vector("list", nrow(df.itr))
for (r in 1:nrow(df.itr)){
  period.est.smooth.fixed[[r]] = parallel.period.2Dsmooth.fun(pp.data,
                                                              i=df.itr[r,1], j=df.itr[r,2],
                                                              inten.fun.sep="~1",
                                                              bandwidth = band.fixed)
}

saveRDS(period.est.smooth.fixed, paste0(model, "_A", sidelength,"_ksde_opt.rds"))








##############################################################################
###### Kernel spectral density estimator with cross-validated bandwidth ######
##############################################################################

### Bandwidth selection
band.range = seq(0.12, 0.55, 0.025)
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
band.cv = foreach(q = 1:length(pp.data),
                  .combine = 'c',
                  .packages = 'spatstat') %dopar%{
                    cv = bandwidth.cv.fun.sim2(ppp = pp.data[[q]],
                                               inten.fun.sep = "~1",
                                               band.range = band.range)
                    cv$`Optimal bandwidth`
                  }
stopCluster(cl)
table(band.cv)
saveRDS(band.cv, paste0(model, "_A", sidelength,"_optbandwidth.rds"))
end_time = Sys.time()
end_time - start_time


### Estimation
period.est.smooth = vector("list", nrow(df.itr))
for (r in 1:nrow(df.itr)){
  period.est.smooth[[r]] = parallel.period.2Dsmooth.fun(pp.data,
                                                        i=df.itr[r,1], j=df.itr[r,2],
                                                        inten.fun.sep = "~1",
                                                        bandwidth = band.cv)
}

saveRDS(period.est.smooth, paste0(model, "_A", sidelength,"_ksde_cv.rds"))