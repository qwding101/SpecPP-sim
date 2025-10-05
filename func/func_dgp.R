#########################################
###########  Import packages  ###########
#########################################
library(spatstat) 
library(parallel)


#########################################
###########  Import functions  ##########
#########################################

# Jalilian et al.'s (2015) codes for generating multivariate product-shot-noise Cox process
# Modified from https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12339

### Kernel functions
bkernels <- list()
# Gaussian kernel
bkernels$Thomas <- function(r, bandwidth, ...){ 
  exp(- r^2/(2 * bandwidth^2)) / (2 * pi * bandwidth^2)
}
# Variance-Gamma (Bessel) kernel with bandwidth and shape parameter nu.ker
bkernels$VarGamma <- function(r, bandwidth, nu.ker){
  stopifnot(nu.ker > -1/2)
  const <- 1 / (4 * pi * nu.ker * bandwidth^2)
  u <- r/bandwidth
  u <- ifelse(u > 0, (u^nu.ker) * besselK(u, nu.ker) / (2^(nu.ker - 1) * gamma(nu.ker)), 1)
  return(abs(const * u))
}
# Cauchy kernel
bkernels$Cauchy <- function(r, bandwith, ...){ 
  ((1 + (r / bandwith)^2)^(-1.5)) / (2 * pi * bandwith^2)
}


### c_l function in Jalilian et al. (2015) p. 1024
c.fun <- list()
c.fun$Thomas <- function(r, bandwith, kappa, ...){
  (kappa/(bkernels[["Thomas"]](0, bandwith)^2)) * exp(-r^2/(4*bandwith^2))/(4*pi*bandwith^2*kappa)
}
c.fun$VarGamma <- function(r, bandwith, kappa, nu.ker){
  (kappa/(bkernels[["VarGamma"]](0, bandwith, nu.ker))^2) * besselK(r/bandwith, nu.ker) * (r/bandwith)^(2*nu.ker+1) / (pi*bandwith^2*2^(2*nu.ker+2)*gamma(2*nu.ker+2)*kappa)
}
c.fun$Cauchy <- function(r, bandwith, ...){
  (kappa/(bkernels[["Cauchy"]](0, bandwith))^2) * (1+r^2/(4*bandwith^2))^(-1.5) / (8*pi*bandwith^2*kappa)
}


### Pair correlation function
g.fun <- function(r, i, j, kappa, xi, kernels, bandwith, nu.ker=NULL){
  
  num.parent = length(kappa)
  
  if (i > num.parent || j > num.parent){
    stop("Index out of bound.", call. = F)
  }
  if (identical(dim(xi), rep(num.parent, 2)) == F){
    stop(paste("xi should be a", num.parent,"by", num.parent, "matrix."), call. = F)
  }
  if (length(kernels) != num.parent || length(bandwith) != num.parent){
    stop("Please specify a kernel function, and its parameter(s), for each parent process.", 
         call. = F)
  }
  
  
  exp.term <- function(i, j, num.parent, kappa, xi, kernels, bandwith, nu.ker=NULL){
    sum.term <- 0
    for (l in (1:num.parent)[-c(i,j)]){
      sum.term <- sum.term + kappa[l]*xi[l,i]*xi[l,j]*c.fun[[kernels[l]]](r, bandwith[l], kappa[l], nu.ker[l])
    }
    return(exp(sum.term))
  }
  
  # Return pair correlation function if i = j
  if (i == j){ 
    first.term <- 1 + (bkernels[[kernels[i]]](0, bandwith[i], nu.ker[i])^2 * 
                        c.fun[[kernels[i]]](r, bandwith[i], kappa[i], nu.ker[i]))/kappa[i]
  }
  # Return cross-pair correlation function if i =/= j
  else{
    first.term <- (1+xi[i,j]*bkernels[[kernels[i]]](0, bandwith[i], nu.ker[i])*c.fun[[kernels[i]]](r, bandwith[i], kappa[i], nu.ker[i])) *
      (1+xi[j,i]*bkernels[[kernels[j]]](0, bandwith[j], nu.ker[j])*c.fun[[kernels[j]]](r, bandwith[j], kappa[j], nu.ker[j]))
  }
  
  return(first.term * exp.term(i, j, num.parent, kappa, xi, kernels, bandwith, nu.ker))
  
}






### Simulation function: general cases
# kappa: intensity of parent Poisson process
# rho: intensity function of individual Cox process
# bandwidth: kernel bandwidth
# nu.ker: shape parameter in variance gamma kernel
# xi: matrix of interaction parameters
# cluster: a string vector indicating the type of cluster kernel
# dimyx: dimension of discretization 

rmkppm <- function(rho, kappa, bandwidth, xi, win=owin(), clusters=NULL, nsim=1,
                   nu.ker=NULL, ij=NULL, eps = NULL, dimyx = NULL, xy = NULL, epsth=0.001, mc.cores=1L)
{
  m <- length(rho)
  if ((length(kappa) != m) || length(bandwidth) != m ) 
    stop("kappa and bandwidth paramters must be of the same size.")
  if (!all(dim(xi) == c(m, m)))
    stop("xi paramter is not a matrix of correct dimensions.")
  if (is.null(clusters))
    clusters <- rep("Thomas", m)
  else if(length(clusters) != m)
    stop("clusters must be a vector of the size of the number of components.")
  if (is.null(nu.ker))
    nu.ker <- rep(-1/4, m)
  diag(xi) <- 0
  if (all(xi == 0))
    return(rmkppm0(rho=rho, kappa=kappa, bandwidth=bandwidth, win=win, clusters=clusters, nsim=nsim,
                   nu.ker=nu.ker, ij=ij, eps = eps, dimyx = dimyx, xy = xy, epsth=epsth, mc.cores=mc.cores))
  if (is.list(rho) == F){
    rho <- as.list(rho)
  }
  frame <- boundingbox(win)
  dframe <- diameter(frame)
  W <- as.mask(win, eps = eps, dimyx = dimyx, xy = xy)
  wx <- as.vector(raster.x(W))
  wy <- as.vector(raster.y(W))
  
  sigma <- rmax <- numeric(m)
  for (i in 1:m)
  {
    if(is.im(rho[[i]])){
      rho[[i]] <- as.im(rho[[i]], dimyx=W$dim)
    }
    if(is.expression(rho[[i]])){
      rho[[i]] <- as.im(rho[[i]], W=win, dimyx=W$dim)
    }
    
    keri <- function(r){ bkernels[[clusters[i]]](r, bandwidth[i], nu.ker[i]) }
    keri0 <- keri(0)
    sigma[i] <- kappa[i] / keri0
    kerithresh <- function(r){ keri(r) / keri0 - epsth}
    rmax[i] <- uniroot(kerithresh, lower = bandwidth[i], upper = 5 * dframe)$root # 4 * bandwidth[i] #
  }
  dilated <- grow.rectangle(frame, max(rmax))
  
  corefun <- function(idumm){
    library(spatstat) # Add for parLapply() to work
    Phi <- lapply(kappa, rpoispp, win=dilated) # Simulate parent process with intensity kappa
    fr <- vector("list", length=m)
    for (i in 1:m){
      keri <- function(r){ bkernels[[clusters[i]]](r, bandwidth[i], nu.ker[i]) }
      keri0 <- keri(0)
      fr[[i]] <- keri(crossdist.default(wx, wy, Phi[[i]]$x, Phi[[i]]$y)) / keri0
    }
    
    if (is.null(ij))
      ij <- 1:m
    xp <- yp <- mp <- NULL
    for (i in 1:m){
      Si <- rowSums(fr[[i]])  / sigma[i] # Shot-noise field
      E <- matrix(1, nrow=length(wx), ncol=m)
      for (j in (1:m)[-i]){
        E[, j] <- apply(1 + xi[j, i] * fr[[j]], 1, prod) * exp(-xi[j, i] * sigma[j]) # Product fields
      }
      values <-  Si * apply(E, 1, prod) # Shot-noise * compound fields
      Lam <- im(values, xcol=W$xcol, yrow=W$yrow, unitname = unitname(W))
      rhoi <- rho[[i]]
      
      #Xi <- rpoispp(eval.im(rhoi * Lam)) # Eq (8), which has bug when some intensity values are negative
      Xi <- rpoispp(eval.im(ifelse(rhoi * Lam < 0, 0, rhoi * Lam))) # Eq (8) but force negative values to be zero

      xp <- c(xp, Xi$x)
      yp <- c(yp, Xi$y)
      mp <- c(mp, rep(ij[i], Xi$n))
    }
    
    simout <- ppp(xp, yp, window=win, marks=as.factor(mp))
    # attr(simout, "parents") <- Phi
  }
  outlist <- if (Sys.info()["sysname"] == "Windows") {
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl=cl, c('bkernels'))
    out = parLapply(cl, 1:nsim, corefun)
    stopCluster(cl)
    out
    # lapply(1:nsim, corefun)
  } 
  else mclapply(1:nsim, corefun, mc.cores=mc.cores) 
  if (nsim == 1){
    return(outlist[[1]])
  } 
  names(outlist) <- paste("Simulation", 1:nsim)
  return(outlist)
}


### Simulation function: null model of independent shot-noise components
rmkppm0 <- function(rho, kappa, bandwidth, win=owin(), clusters=NULL, nsim=1,
                    nu.ker=NULL, ij=NULL, eps = NULL, dimyx = NULL, xy = NULL, epsth=0.001, mc.cores=1L)
{
  m <- length(rho)
  if ((length(kappa) != m) || length(bandwidth) != m ) 
    stop("kappa and bandwidth paramters must be of the same size.")
  if (is.null(clusters))
    clusters <- rep("Thomas", m)
  else if(length(clusters) != m)
    stop("clusters must be a vector of the size of the number of components.")
  if (is.null(nu.ker))
    nu.ker <- rep(-1/4, m)
  rho <- as.list(rho)
  if (is.null(ij))
    ij <- 1:m
  corefun0 <- function(dumm)
  {
    xp <- yp <- mp <- NULL
    for (i in 1:m)
    {
      rhoi <- rho[[i]]
      if (is.numeric(rhoi))
        mui <- rhoi / kappa[i]
      else if (is.im(rhoi)) 
        mui <- eval.im(rhoi / kappa[i])
      Xi <- switch(clusters[i], 
                   Thomas = rThomas(kappa[i], bandwidth[i], mui, win=win),
                   Cauchy =  rCauchy(kappa[i], bandwidth[i], mui, win=win, eps=epsth),
                   VarGamma = rVarGamma(kappa[i], nu.ker=nu.ker[i], bandwidth[i], mui, win=win, eps=epsth, nu.pcf=NULL))
      xp <- c(xp, Xi$x)
      yp <- c(yp, Xi$y)
      mp <- c(mp, rep(ij[i], Xi$n))
    }
    
    out <- ppp(xp, yp, window=win, marks=as.factor(mp))
    #  attr(out, "parents") <- Phi
  }
  outlist <- if (mc.cores == 1) lapply(1:nsim, corefun0) 
  else mclapply(1:nsim, corefun0, mc.cores=mc.cores)
  if (nsim == 1) 
    return(outlist[[1]])
  names(outlist) <- paste("Simulation", 1:nsim)
  return(outlist)
}




