###################################################
###### Import functions and model parameters ######
###################################################

source("../../func/func_sim.R", chdir = T)

kappa = c(0.25, 0.75, 0.2) # Parent intensity
a = 0.025 # Data taper
inten.fun.sep = list()
df.itr = data.frame(c(1,2,1), c(1,2,2))
bandwidth.li = list(sigma = c(0.6, 0.3, 1), sigma = c(0.6, 0.3, 1), sigma = c(0.6, 0.3, 1))
nu.li = list(NULL, NULL, NULL)
xi.1 = matrix(c(0, 0.7, 0, 0.9, 0, 0, 0.3, 0.1, 0), 3, 3, byrow = T) 
xi.2 = matrix(c(0, -0.7, 0, -0.9, 0, 0, 0.3, 0.1, 0), 3, 3, byrow = T)
xi.li = list(xi.1, xi.1, xi.2)
kernels.li = list(rep("Thomas", length(kappa)),
                  rep("Thomas", length(kappa)),
                  rep("Thomas", length(kappa)))




########################################################
###### Plot the pseudo-spectra for M1, M2, and M3 ######
########################################################

linetype = c("solid", "dashed", "dotdash")
ylim = list(c(0,0.4), c(0,0.8), c(-0.3,0.8))


cairo_pdf("fig_1b.pdf", width = 8, height = 2.5, bg = "transparent")

par(mfrow = c(1, 3),  mar = c(4, 3, 2, 0.5))
w = seq(0, sqrt((1.5*pi)^2 + (1.5*pi)^2), .01) # Frequency sequence
for (i in seq_along(bandwidth.li[1:3])){
  if (i == 1){ # The intensity of M1
    inten.fun.sep$X1 = function(x, x_i){ 
      return(0.5)
    }
    inten.fun.sep$X2 = function(x, x_i){ 
      return(1.5)
    }
  }else{ # The intensity functions of M2 and M3
    inten.fun.sep$X1 = function(x, x_i){
      if (x_i == 1) return(3*exp(-2*x^2))
      if (x_i == 2) return(exp(-2*x^2))
      # The product of above two terms is 3*exp( -2*( (x1)^2 + (x2)^2 ) )
    }
    inten.fun.sep$X2 = function(x, x_i){
      if (x_i == 1) return(2*exp(-2*x^2))
      if (x_i == 2) return(exp(2*x^2))
      # The product of above two terms is 2*exp( -2*( (x1)^2 - (x2)^2 ) )
    }
  }
  
  for (r in 1:nrow(df.itr)){
    # Compute the pseudo-spectra
    f.h = f.fun(w, i=df.itr[r,1], j=df.itr[r,2], inten.fun.sep=inten.fun.sep,
                a=a, kappa=kappa, xi=xi.li[[i]], kernels=kernels.li[[i]],
                bandwidth=bandwidth.li[[i]], nu.ker=nu.li[[i]])
    if (r == 1){
      plot(w, f.h, type = "l", las = 1, bty='l', ylim = ylim[[i]],
           lty = linetype[r], lwd = 2,
           xlab = expression(paste("||", bold("\u03c9"), "||")), ylab = "",
           main = paste0("M",i),
           cex.main = 2, cex.axis = 1.5, cex.lab = 1.75)
      if (i == max(seq_along(bandwidth.li[1:3]))){
        legend('topright',
               legend = c(expression(italic(F[h]^{(paste(1,",",1))})),
                          expression(italic(F[h]^{(paste(2,",",2))})),
                          expression(italic(F[h]^{(paste(1,",",2))}))),
               lty = linetype, #col = color,
               lwd = 2, xpd = T, bty = 'n', seg.len = 4)
      }
    }else{
      lines(w, f.h, lty = linetype[r], lwd = 2)
    }
  }
  abline(h=0, col = rgb(0, 0, 0, .3), lty = "dotted")
}

dev.off()
