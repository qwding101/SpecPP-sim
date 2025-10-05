###########################################################
###### Import functions to calculate IBIAS and IMSE #######
###########################################################

source("../../func/func_sim.R", chdir = T)



###################################
###### Simulation parameters ###### 
###################################

model = "M1"
a = 0.025 # Data taper
Ac = c("A10","A20","A40")
df.itr = data.frame(c(1,2,1), c(1,2,2))
kappa = c(0.25, 0.75, 0.2) # Parent intensity
sigma = c(0.6, 0.3, 1) # Standard deviation in normal kernel
nu = NULL
xi = matrix(c(0, 0.7, 0, 0.9, 0, 0, 0.3, 0.1, 0), 3, 3, byrow = T) # Interaction parameters
kernels = rep("Thomas", length(kappa))
inten.fun.sep = list(X1 = 0.5, X2 = 1.5) # Intensity




########################################
######### Compute true spectra ######### 
########################################

true.spec = list(A10 = list(), A20 = list(), A40 = list())
freq.li = list(A10 = list(omega1 = freq.grid.fun(10, 10)$omega1,
                          omega2 = freq.grid.fun(10, 10)$omega2),
               A20 = list(omega1 = freq.grid.fun(20, 20)$omega1,
                          omega2 = freq.grid.fun(20, 20)$omega2),
               A40 = list(omega1 = freq.grid.fun(40, 40)$omega1,
                          omega2 = freq.grid.fun(40, 40)$omega2))

for (a.idx in seq_along(true.spec)){
  for (r in 1:nrow(df.itr)){
    true.spec[[a.idx]][[r]] = outer(freq.li[[a.idx]]$omega2, freq.li[[a.idx]]$omega1,
                                    FUN = function(w2, w1, i, j, inten.fun.sep, a, kappa, xi,
                                                   kernels, bandwidth, nu.ker) 
                                      f.fun(sqrt(w1^2+w2^2), i, j, inten.fun.sep, a, kappa, xi,
                                            kernels, bandwidth, nu.ker),
                                    i = df.itr[r,1], j = df.itr[r,2], inten.fun.sep=inten.fun.sep, a = a,
                                    kappa = kappa, xi = xi, kernels = kernels, bandwidth = sigma, nu.ker = nu)
    attr(true.spec[[a.idx]][[r]], "omega1") = freq.li[[a.idx]]$omega1 # Col frequency
    attr(true.spec[[a.idx]][[r]], "omega2") = freq.li[[a.idx]]$omega2 # Row frequency
  }
}



########################################
###### Import spectrum estimators ###### 
########################################

I.hat = list(A10 = readRDS(paste0("../../estimate/", model, "/", Ac[1], "/", model, "_", Ac[1],"_period.rds")),
             A20 = readRDS(paste0("../../estimate/", model, "/", Ac[2], "/", model, "_", Ac[2],"_period.rds")),
             A40 = readRDS(paste0("../../estimate/", model, "/", Ac[3], "/", model, "_", Ac[3],"_period.rds")))
f.opt = list(A10 = readRDS(paste0("../../estimate/", model, "/", Ac[1], "/", model, "_", Ac[1],"_ksde_opt.rds")),
             A20 = readRDS(paste0("../../estimate/", model, "/", Ac[2], "/", model, "_", Ac[2],"_ksde_opt.rds")),
             A40 = readRDS(paste0("../../estimate/", model, "/", Ac[3], "/", model, "_", Ac[3],"_ksde_opt.rds")))
f.cv = list(A10 = readRDS(paste0("../../estimate/", model, "/", Ac[1], "/", model, "_", Ac[1],"_ksde_cv.rds")),
            A20 = readRDS(paste0("../../estimate/", model, "/", Ac[2], "/", model, "_", Ac[2],"_ksde_cv.rds")),
            A40 = readRDS(paste0("../../estimate/", model, "/", Ac[3], "/", model, "_", Ac[3],"_ksde_opt.rds")))





##################################
###### Compute IBIAS & IMSE ######
##################################

summary.list = list(A10 = matrix(NA, ncol = nrow(df.itr), nrow = 6),
                    A20 = matrix(NA, ncol = nrow(df.itr), nrow = 6),
                    A40 = matrix(NA, ncol = nrow(df.itr), nrow = 6))
freq.remove = c(-pi/10, pi/10) # Exclude low-frequency values (∥ω∥∞ < 0.1π) when evaluating IBIAS and IMSE

for (a.idx in seq_along(summary.list)){
  colnames(summary.list[[a.idx]]) = c("11", "22", "12")
  rownames(summary.list[[a.idx]]) = c("IBIAS (I.hat)", "IMSE (I.hat)",
                                      "IBIAS (f.hat.opt)", "IMSE (f.hat.opt)",
                                      "IBIAS (f.hat.cv)", "IMSE (f.hat.cv)")
  for (r in 1:nrow(df.itr)){
    summary.list[[a.idx]][1, r] = IBIAS.fun(I.hat[[a.idx]][[r]], true.spec[[a.idx]][[r]], remove = freq.remove)
    summary.list[[a.idx]][2, r] = IMSE.fun(I.hat[[a.idx]][[r]], true.spec[[a.idx]][[r]], remove = freq.remove)
    summary.list[[a.idx]][3, r] = IBIAS.fun(f.opt[[a.idx]][[r]], true.spec[[a.idx]][[r]], remove = freq.remove)
    summary.list[[a.idx]][4, r] = IMSE.fun(f.opt[[a.idx]][[r]], true.spec[[a.idx]][[r]], remove = freq.remove)
    summary.list[[a.idx]][5, r] = IBIAS.fun(f.cv[[a.idx]][[r]], true.spec[[a.idx]][[r]], remove = freq.remove)
    summary.list[[a.idx]][6, r] = IMSE.fun(f.cv[[a.idx]][[r]], true.spec[[a.idx]][[r]], remove = freq.remove)
  }
}

summary.list

saveRDS(summary.list, paste0(model, "_summary.rds"))
