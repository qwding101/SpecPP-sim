##############################################################################
#### This version yields both real and imaginary parts of the estimators. ####
#### Everything else is the same as func_sim.R.                           ####
##############################################################################


library(spatstat)
library(ggplot2)
library(foreach)
library(doParallel)
library(parallel)
library(cowplot)
library(igraph)
library(ggraph)
library(dplyr)

fixpd.fun = function(mat){
  min.eigen = min(eigen(mat, only.values=T)$values)
  if (min.eigen <= 0){
    mat = mat + (1e-3 - min.eigen)*diag(1, nrow(mat))
  }
  return(mat)
}

between.fun = function(x, lower, upper) return(x >= lower & x <= upper)

cate.comb.fun = function(cate){
  pos1 = pos2 = c()
  for (i in seq_along(cate)){
    pos1 = append(pos1, rep(cate[i], length(cate)-i+1))
    pos2 = append(pos2, cate[i:length(cate)])
  }
  return(cbind(pos1, pos2))
}

ha.fun = function(x, a){
  # x is vectorized, while a is a number
  if (a == 0){
    result = rep(1, length(x))
  }else{
    result = rep(NA, length(x))
    idx = x >= -0.5 & x <= a - 0.5
    result[idx] = ((x[idx]+0.5)/a) - ((1/(2*pi))*sin(2*pi*(x[idx]+0.5)/a))
    result[x > a-0.5 & x < 0.5-a] = 1
    idx = x >= 0.5-a & x <= 0.5
    result[idx] = ((0.5-x[idx])/a) - ((1/(2*pi))*sin(2*pi*(0.5-x[idx])/a))
  }
  return(result)
}

Bartlett.fun = function(v, b){
  return(ifelse(abs(v/b) > 1, 0, (1-abs(v/b))/b))
}

KK.fun = function(x1, x2, omega1, omega2, b1, b2=b1){
  return(Bartlett.fun(omega1 - x1, b1)*Bartlett.fun(omega2 - x2, b2))
}

freq.grid.fun = function(A1, A2, ext.factor = NULL, return.comb = F, endpt = 1.5){

  freq = list()
  # Generate frequency grid
  stopifnot("The value of side lengths should be larger than zero." = A1*A2 > 0)
  A1 = round(A1); A2 = round(A2)
  k1 = -A1:A1; k2 = -A2:A2
  omega1 = endpt*pi*k1/A1; omega2 = endpt*pi*k2/A2

  # Extend frequency gird to handle boundary effect when conducting kernel smoothing
  # E.g., original frequency: -1.5π ~ 1.5π
  #    => extended frequency: -1.5π*ext.factor ~ 1.5π*ext.factor
  if (!is.null(ext.factor)){
    dx = omega1[2] - omega1[1]
    dy = omega2[2] - omega2[1]
    ext.x = seq(max(omega1), max(omega1)*ext.factor, by = dx)[-1]
    ext.y = seq(max(omega2), max(omega2)*ext.factor, by = dy)[-1]
    omega1 = c(sort(-ext.x), omega1, ext.x)
    omega2 = c(sort(-ext.y), omega2, ext.y)
  }
  freq[[1]] = omega1
  freq[[2]] = omega2
  names(freq) = c("omega1", "omega2")

  if (return.comb){
    # Below code is basically `as.matrix(expand.grid(list(omega1, omega2, NA)))` but much faster
    freq[[3]] = cbind(omega1 = rep(omega1, length(omega2)),
                      omega2 = rep(omega2, each = length(omega1)),
                      Density = NA)
    names(freq)[3] = "omega.comb"
  }

  return(freq)
}

smooth.fun = function(w, period.mat, w.k1, w.k2, b1=1, b2=b1, loo = F){

  dy = w.k2[2] - w.k2[1]

  # Narrow the range of discretized values to iterate,
  # which is related to the kernel bandwidth b1 and b2
  idx.w1 = w.k1 >= w[1] - b1 & w.k1 <= w[1] + b1
  idx.w2 = w.k2 >= w[2] - b2 & w.k2 <= w[2] + b2
  w1.range = w.k1[idx.w1]
  w2.range = w.k2[idx.w2]
  w1.pt2smooth.actual = length(w1.range) # Number of points to smooth (including w) for x-direction in reality
  w2.pt2smooth.actual = length(w2.range) # Number of points to smooth (including w) for y-direction in reality

  # Detect small bandwidth problem: there is no any neighbor to smooth for at least one direction
  if (w1.pt2smooth.actual == 1 || w2.pt2smooth.actual == 1){
    if (loo){
      stop(paste0("Either the minimal bandwidth is too small or the frequency grid is too coarse. ",
                 "Given current frequency grid, the bandwidth should be larger than ", round(dy,4), "."),
           call. = F)
    }else{
      stop(paste0("Either the bandwidth is too small or the frequency grid is too coarse. ",
                 "Given current frequency grid, the bandwidth should be larger than ", round(dy,4), "."),
           call. = F)
    }
  }

  # Given a (w[c], w[r]) pair, calculate the kernel product matrix
  # for the range of discretized values x1.range and x2.range.
  kernel.product = t(outer(w1.range, w2.range, KK.fun,
                           omega1 = w[1], omega2 = w[2], b1 = b1, b2 = b2))

  if (loo){ # For LOOCV, remove the kernel density for central point
    kernel.product[which(w2.range == w[2]), which(w1.range == w[1])] = 0
  }

  # Volume = height of endpoint * subrectangle area, where area = dx*dy,
  # which is cancelled out in below ratio
  numerator = sum((kernel.product*period.mat[idx.w2, idx.w1]))
  denominator = sum(kernel.product)

  return(ifelse(denominator == 0, 0, numerator/denominator))
}

scale2unitarea = function(im, length.x, length.y, return.fun = T){
  if (!is.im(im)){
    im = as.im.funxy(im)
  }
  im.new = im
  im.new$xrange = im.new$yrange = c(-0.5, 0.5)
  im.new$xstep = im$xstep/length.x
  im.new$ystep = im$ystep/length.y
  im.new$xcol = im$xcol/length.x
  im.new$yrow = im$yrow/length.y

  if (return.fun){
    return(as.function.im(im.new))
  } else{
    return(im.new)
  }
}

plot.period.fun = function(period, freq.list = NULL, remove=NULL, title=NULL,
                           palette="Spectral", legend.range=NULL){

  if (is.matrix(period)){
    if (is.null(freq.list)){
      omega.comb = attr(period, "freq.list")$omega.comb
    }else{
      omega.comb = freq.list$omega.comb
    }
    period.long = as.data.frame(omega.comb)
    period.long$dist = sqrt(period.long$omega1^2 + period.long$omega2^2)
    period.long$Density = as.vector(t(period))
  }
  else if (is.data.frame(period)){ # Seems useless now because no function exports long format?
    period.long = period
  }
  else{
    stop("The periodogram should be a matrix or long-format data.frame.", call. = F)
  }

  if (!is.null(remove)){
    period.long$Density[between.fun(period.long$omega1, -remove[1], remove[1]) &
                        between.fun(period.long$omega2, -remove[2], remove[2])] = NA
  }

  ggplot(period.long, aes(omega1, omega2)) +
    geom_raster(aes(fill = Density)) +
    scale_fill_distiller(palette = palette, name = "", direction = 1,
                         limit = legend.range, na.value = "transparent") +
    labs(title = title, x = expression(omega[1]), y = expression(omega[2])) +
    theme_classic()
}

plot.pairs.fun = function(period.smooth.mat.list, ppp, xnorm = T, freq.list = NULL,
                          shared.legend = T, remove = NULL){

  cate = levels(marks(ppp))
  plot.layout = matrix(NA, ncol = length(cate), nrow = length(cate))
  plot.layout[lower.tri(plot.layout, diag = T)] = seq_along(period.smooth.mat.list)

  titles = rep("", length(plot.layout))
  titles[diag(plot.layout)] = cate

  legend.range = switch(shared.legend + 1, # Basically `ifelse()` but ifelse() cannot return NULL
                    NULL,
                    c(min(sapply(period.smooth.mat.list, min)),
                      max(sapply(period.smooth.mat.list, max)))
  )

  period.plot.li = vector("list", length(period.smooth.mat.list))


  if (xnorm){
    if (is.null(freq.list)){
      omega.comb = attr(period.smooth.mat.list, "freq.list")$omega.comb
    }else{
      omega.comb = freq.list$omega.comb
    }
    period.long = as.data.frame(omega.comb)
    period.long$dist = sqrt(period.long$omega1^2 + period.long$omega2^2)
for (r in seq_along(period.smooth.mat.list)){

  period.long$Density = as.vector(t(period.smooth.mat.list[[r]]))

  period.plot.li[[r]] = group_by(period.long, dist) %>%
    summarise(avg.density = mean(Density)) %>%
    ggplot(aes(dist, avg.density)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    theme_classic() +
    labs(title = titles[r], x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(ylim = legend.range)
}

result = plot_grid(plotlist = period.plot.li[as.vector(plot.layout)],
                   nrow = length(cate), ncol = length(cate), byrow = F,
                   scale = 0.85)
result = result +
  draw_label(expression(paste("||", bold("\u03c9"), "||")), x=0.5, y=0, vjust=-0.5, size=12, angle=0) +
  draw_label("Radially averaged density", x=0, y=0.5, vjust= 1.5, size=12, angle=90)

  }else{
for (r in seq_along(period.smooth.mat.list)){
  period.plot.li[[r]] = plot.period.fun(period.smooth.mat.list[[r]],
                                        freq.list = attr(period.smooth.mat.list, "freq.list"),
                                        remove = remove,
                                        title = names(period.smooth.mat.list)[r],
                                        legend.range = legend.range) +
    labs(title = titles[r], x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5))
}
if (shared.legend){
  period.plot.li = lapply(period.plot.li, function(x) x + guides(fill = "none"))
  common.legend = get_legend(plot.period.fun(period.smooth.mat.list[[1]],
                                             legend.range = legend.range))
  result = plot_grid(plotlist = period.plot.li[as.vector(plot.layout)],
                     nrow = length(cate), ncol = length(cate), byrow = F, scale = 0.85)
  result = plot_grid(result, common.legend, ncol = 2,
                     rel_widths = c(1, .1))

}else{
  result = plot_grid(plotlist = period.plot.li[as.vector(plot.layout)],
                     nrow = length(cate), ncol = length(cate), byrow = F,
                     scale = 0.85)
}

result = result +
  draw_label(expression(omega[1]), x=0.5, y=0, vjust=-0.5, angle=0) +
  draw_label(expression(omega[2]), x=0, y=0.5, vjust= 1.5, angle=90)
  }

  return(result)
}

plot.coher.fun = function(sp.est, coh.max.mat, partial.coh.max.mat,
                          xnorm = T, ylim = NULL){

  # sp.est: a list of kernel spectral estimate matrices of the original point process
  # coh.max.mat: a m by m matrix storing maximum coherences
  # partial.coh.max.mat: a m by m matrix storing maximum partial coherences

  cate = attr(sp.est, "cate")
  cate.comb = data.frame(attr(sp.est, "cate.comb"))
  coh.allw = data.frame(cbind(
    attr(coh.max.mat, "CohTable"),
    attr(partial.coh.max.mat, "CohTable")[, -(1:2)])
  )
  coh.allw$dist = sqrt(coh.allw$omega1^2 + coh.allw$omega2^2)
  colnames(coh.allw) = c("omega1", "omega2",
                         paste0("C", 1:(nrow(cate.comb)-length(cate))),
                         paste0("PC", 1:(nrow(cate.comb)-length(cate))), "dist")

  plot.order = diag(rep(999, length(cate)))
  plot.order[upper.tri(plot.order)] = paste0("C", 1:(nrow(cate.comb)-length(cate)))
  plot.order[lower.tri(plot.order)] = t(plot.order)[lower.tri(plot.order)]
  plot.order[lower.tri(plot.order)] = paste0("P", plot.order[lower.tri(plot.order)])
  diag(plot.order) = cate
  plot.order = as.vector(plot.order)


  plot.list = vector("list", length = length(plot.order))

  if (xnorm){

    for (k in seq_along(plot.order)){
      if (plot.order[k] %in% cate){
        plot.list[[k]] = ggplot() +
          annotate("text", x = 1, y = 1, size = 5, label = plot.order[k]) +
          theme_void()
      }else{
        df = coh.allw[, c("omega1","omega2", "dist", plot.order[k])]
        colnames(df)[4] = "Value"
        plot.list[[k]] = group_by(df, dist) %>%
          summarise(avg.value = mean(Value)) %>%
          ggplot(aes(dist, avg.value)) +
          geom_line() +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
          theme_classic() +
          labs(x = NULL, y = NULL) +
          theme(plot.title = element_text(hjust = 0.5),
                panel.background = element_rect(fill='transparent'),
                plot.background = element_rect(fill='transparent', color=NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          coord_cartesian(ylim = ylim)
      }
    }
    common.title = ggdraw() +
      draw_label("Radially averaged (partial) coherence",
                 fontface = 'bold', x = 0, hjust = 0)
    result = plot_grid(plotlist = plot.list,
                       nrow = length(cate), ncol = length(cate),
                       byrow = F, scale = 0.8)
    result = result +
      draw_label(expression(paste("||", bold("\u03c9"), "||")), x=0.5, y=0, vjust=-0.5, size=12, angle=0) +
      draw_label("Coherence", x=1, y=0.5, vjust=1.5, angle=-90) +
      draw_label("Partial coherence", x=0, y=0.5, vjust= 1.5, angle=90)
    result = plot_grid(common.title, result, ncol = 1, rel_heights = c(0.05, 1))

  }else{

    for (k in seq_along(plot.order)){
      if (plot.order[k] %in% cate){
        plot.list[[k]] = ggplot() +
          annotate("text", x = 1, y = 1, size = 5, label = plot.order[k]) +
          theme_void()
      }else{
        df = coh.allw[, c("omega1","omega2", plot.order[k])]
        colnames(df)[3] = "Density"
        plot.list[[k]] = plot.period.fun(df, palette = "Spectral") +
          labs(x = NULL, y = NULL)
      }
    }
    result = plot_grid(plotlist = plot.list,
                       nrow = length(cate), ncol = length(cate), byrow = F, scale = 0.8)
    result = result +
      draw_label("Coherence", x=1, y=0.5, vjust=1.5, angle=-90) +
      draw_label("Partial coherence", x=0, y=0.5, vjust= 1.5, angle=90)
  }

  return(result)
}

period.fun.sim = function(i, j, ppp, inten.fun.sep="est", a = 0.025,
                          return.DFT = F, A1 = NULL, A2 = A1, ext.factor = NULL, endpt = 1.5){

  # Shift the observation window to the origin
  ppp = shift(ppp, origin = "centroid")

  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = freq.grid.fun(A1, A2, ext.factor = ext.factor, endpt = endpt)

  exp.term = function(w, x, a, A){
    return(ha.fun(x/A, a)*exp(-1i*x*w))
  }

  int.prod = function(i, w_q, w_p, a, A1, A2, inten.fun.sep){
    int.fun = Vectorize( # Vectorize the argument w
      function(i, w, a, A, x_i){
        AA = ifelse(attr(inten.fun.sep, "scale") == "scaled", 1, A)
        real.part = integrate(function(x, a, A, x_i, i, w) Re(ha.fun(x,a)*inten.fun.sep[[i]](x*AA,x_i)*exp(-1i*x*A*w)),
                         lower=-1/2, upper=1/2, a=a, A=A, x_i=x_i, i=i, w=w)$value
        imag.part = integrate(function(x, a, A, x_i, i, w) Im(ha.fun(x,a)*inten.fun.sep[[i]](x*AA,x_i)*exp(-1i*x*A*w)),
                         lower=-1/2, upper=1/2, a=a, A=A, x_i=x_i, i=i, w=w)$value
        return(real.part + 1i*imag.part)
      }
    )
    return(area(ppp.li[[i]])*int.fun(i, w_p, a, A=A1, x_i=1)*int.fun(i, w_q, a, A=A2, x_i=2))
  }

  Jh.fun = function(i){
    if (is.list(inten.fun.sep)){ # Use population intensity (scaled to [-0.5,0.5]^2 input)
      if (is.numeric(inten.fun.sep[[i]])){
        lambda = inten.fun.sep[[i]]
        stopifnot("The population intensity should be larger than zero."= lambda > 0)
        inten.fun.sep[[i]] = function(x, x_i){ 
          if (x_i == 1) return(lambda) 
          if (x_i == 2) return(1) 
        }
      }else if (!is.function(inten.fun.sep[[i]])){
        stop("The element in the list should be a positive value or function.", call. = F)
      }
      attr(inten.fun.sep, "scale") = "scaled"
    }else if (inten.fun.sep == "est"){ # Use plug-in estimator of the intensity function, if inten.fun.sep == "est"
      inten.fun.sep = list()
      coef.est = ppm(ppp.li[[i]] ~ I(x^2) + I(y^2))$coef
      inten.fun.sep[[i]] = function(x, x_i){ 
        if (x_i == 1) return(exp(coef.est[1] + coef.est[2]*x)) 
        if (x_i == 2) return(exp(coef.est[3]*x)) 
      }
      attr(inten.fun.sep, "scale") = "original"
    }else if (inten.fun.sep == "~1"){
      inten.fun.sep = list()
      inten.fun.sep[[i]] = function(x, x_i){
        if (x_i == 1) return(intensity.ppp(ppp.li[[i]]))
        if (x_i == 2) return(1)
      }
      attr(inten.fun.sep, "scale") = "original"
    }

    const = ((2*pi)^(-2/2))*(area(ppp.li[[i]])^(-1/2))/(1+a*(5/(4*pi^2) - 4/3))
    V1 = outer(freq.list$omega1, ppp.li[[i]]$x, exp.term, a=a, A=A1)
    V2 = outer(freq.list$omega2, ppp.li[[i]]$y, exp.term, a=a, A=A2)
    J_h.woconst = V2 %*% t(V1)

    mat.left = outer(freq.list$omega2, freq.list$omega1[freq.list$omega1 < 0], int.prod,
                     i=i, a=a, A1=A1, A2=A2, inten.fun.sep=inten.fun.sep)
    mat.center = outer(freq.list$omega2, 0, int.prod,
                       i=i, a=a, A1=A1, A2=A2, inten.fun.sep=inten.fun.sep)
    mat.right = Conj(matrix(rev(as.vector(mat.left)), ncol = ncol(mat.left), nrow = nrow(mat.left)))
    C_h.woconst = cbind(mat.left, mat.center, mat.right)

    return(const*(J_h.woconst-C_h.woconst))
  }


  if (is.null(marks(ppp))){
    marks(ppp) = 1
    ppp.li = list()
    ppp.li[[1]] = ppp
  }else{
    ppp.li = split(ppp)
  }


  if (i == j){

    DFT = Jh.fun(i) # Centered DFT
    colnames(DFT) = freq.list$omega1
    rownames(DFT) = freq.list$omega2
    period.mat = Mod(DFT)^2 # Periodogram for ith process

    if (return.DFT){
      # Return both DFT and periodogram
      # DFT of both ith and jth processes can be used to compute the cross-spectrum
      return(list(Periodogram = period.mat, DFT = DFT))
    }else{
      # Only return the periodogram
      return(period.mat)
    }
  }
  else{
    # Calculate the cross-spectrum directly
    period.mat = Jh.fun(i)*Conj(Jh.fun(j)) # This part is different from func_sim.R
    colnames(period.mat) = freq.list$omega1
    rownames(period.mat) = freq.list$omega2
    return(period.mat)
  }
}

parallel.period.fun = function(i, j=i, ppp.list, inten.fun.sep="est", a=0.025, endpt=1.5){
  if (Sys.info()["sysname"] == "Windows") {
    cl = makeCluster(detectCores() - 1)
    clusterExport(cl=cl, c('ha.fun', 'freq.grid.fun'))
    clusterEvalQ(cl=cl, library(spatstat))
    result = parLapply(cl, ppp.list,
                       period.fun.sim, i=i, j=j, inten.fun.sep=inten.fun.sep,
                       a=a, return.DFT = T, endpt = endpt)
    stopCluster(cl)
  }else{
    result = mclapply(ppp.list,
                      period.fun.sim, i=i, j=j, inten.fun.sep=inten.fun.sep,
                      a=a, return.DFT = T, endpt = endpt,
                      mc.cores=parallel::detectCores())
  }

  if (i == j){
    return(list(period.list = lapply(result, function(x) x$Periodogram),
                DFT.list = lapply(result, function(x) x$DFT)))
  }else{
    return(result)
  }
}

period.fun = function(i, j, ppp,
                      inten.formula = "~1", data.covariate = NULL,
                      a = 0.025, return.DFT = F, A1 = NULL, A2 = A1,
                      ext.factor = NULL, endpt = 1.5){
  # Shift the observation window to the origin
  ppp = shift(ppp, origin = "centroid")
  
  if (!is.null(data.covariate)){
    data.covariate = lapply(data.covariate, shift, origin = "centroid")
  }
  
  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = freq.grid.fun(A1, A2, ext.factor = ext.factor, endpt = endpt, return.comb = T)
  
  exp.term = function(w, x, a, A){
    return(ha.fun(x/A, a)*exp(-1i*x*w))
  }
  
  H.h.lambda.1 = Vectorize(function(k, w2, w1, a, A1, A2, inten.fitted){
    surface = as.im(X = function(x,y,a,A1,A2,w1,w2) ha.fun(x,a)*ha.fun(y,a)*inten.fitted(x*A1,y*A2)*exp(-1i*(A1*x*w1+A2*y*w2)),
                    W = owin(xrange = c(-.5,.5), yrange = c(-.5,.5)),
                    a = a, A1 = A1, A2 = A2, w1 = w1, w2 = w2)
    return(area(ppp.li[[k]])*integral.im(surface))
  }, vectorize.args = c("w2","w1") # Vectorize the argument w1 & w2 for `outer()` to use
  )
  
  Jh.fun = function(k){
    pppk = ppp[marks(ppp) == k] # Select data from kth process
    
    if (is.null(inten.formula)){ # Estimate intensity function by kernel smoothing
      inten.fitted = as.function(density.ppp(unmark(pppk)))
    } else if (is.list(inten.formula)){ # Use population intensity functions
      stopifnot("The list should include the population intensity functions for all univariate processes."= length(inten.formula) == length(levels(marks(ppp))))
      names(inten.formula) = levels(marks(ppp))
      inten.fitted = inten.formula[[k]]
    } else if (is.character(inten.formula)){
      # Estimate intensity function by log-linear model
      inten.ppm = ppm(Q = as.formula(paste("unmark(pppk)", inten.formula)),
                      data = data.covariate)
      inten.fitted = as.function(predict(inten.ppm))
    }
    
    # Compute DFT
    const = ((2*pi)^(-2/2))*(area(pppk)^(-1/2))/(1+a*(5/(4*pi^2) - 4/3))
    V1 = outer(freq.list$omega1, pppk$x, exp.term, a=a, A=A1)
    V2 = outer(freq.list$omega2, pppk$y, exp.term, a=a, A=A2)
    J_h.woconst = V2 %*% t(V1)
    mat.left = outer(freq.list$omega2, freq.list$omega1[freq.list$omega1 < 0], H.h.lambda.1,
                     k=k, a=a, A1=A1, A2=A2, inten.fitted=inten.fitted)
    mat.center = outer(freq.list$omega2, 0, H.h.lambda.1,
                       k=k, a=a, A1=A1, A2=A2, inten.fitted=inten.fitted)
    mat.right = Conj(matrix(rev(as.vector(mat.left)), ncol = ncol(mat.left), nrow = nrow(mat.left)))
    C_h.woconst = cbind(mat.left, mat.center, mat.right)
    
    DFT = const*(J_h.woconst-C_h.woconst)
    attr(DFT, "inten.fitted") = inten.fitted
    return(DFT)
  }
  
  
  if (is.null(marks(ppp))){
    marks(ppp) = 1
    ppp.li = list()
    ppp.li[[i]] = ppp
  }else{
    ppp.li = split(ppp)
  }
  
  if (i == j){
    
    DFT = Jh.fun(i) # Centered DFT
    colnames(DFT) = freq.list$omega1
    rownames(DFT) = freq.list$omega2
    period.mat = Mod(DFT)^2 # Periodogram for ith process
    attr(period.mat, "freq.list") = freq.list
    
    if (return.DFT){
      # Return both DFT and periodogram
      # DFT of both ith and jth processes can be used to compute the cross-spectrum
      period.DFT = list(Periodogram = period.mat, DFT = DFT)
      attr(period.DFT, "freq.list") = freq.list
      return(period.DFT)
    }else{
      # Only return the periodogram
      return(period.mat)
    }
  }
  else{
    # Calculate the cross-spectrum directly
    period.mat = Jh.fun(i)*Conj(Jh.fun(j))
    colnames(period.mat) = freq.list$omega1
    rownames(period.mat) = freq.list$omega2
    attr(period.mat, "freq.list") = freq.list
    return(period.mat)
  }
}

bandwidth.cv.fun.sim = function(ppp, inten.fun.sep = "est", a = 0.025,
                                band.range, correct = T, A1 = NULL, A2 = A1, endpt = 1.5){

  # Crate frequency index array (freq.idx.df) to sum in likelihood calculation
  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = freq.grid.fun(A1, A2, ext.factor = NULL, return.comb = T, endpt = endpt)
  freq.idx.df = cbind(omega1.ind = rep(seq_along(freq.list$omega1), length(freq.list$omega2)),
                      omega2.ind = rep(seq_along(freq.list$omega2), each = length(freq.list$omega1)),
                      val = NA)

  # Bandwidth to evaluate
  band.mat = cbind(b1 = (A2/A1)*band.range, b2 = band.range)

  # Matrix to save the LOOCV result
  cv = vector("list", 3)
  cv[[3]] = matrix(c(band.range, rep(NA, length(band.range))),
                  ncol = length(band.range), nrow = 2, byrow = T)
  rownames(cv[[3]]) = c("Bandwidth", "Whittle likelihood")
  # Enumerate all possible combination of spectra to compute
  if (is.null(marks(ppp))){ # Univariate point pattern
    cate = 1
  }else{ # Multivariate point pattern
    cate = levels(marks(ppp))
  }
  cate.comb = cate.comb.fun(cate)
  period.mat.list = period.mat.ext.list = vector("list", nrow(cate.comb))

  # Function to compute the Whittle likelihood without the summation part
  ind.Whittle.fun = function(w, period.est, period.est.loo, cate.comb){

    I.mat = F.mat = matrix(NA, length(cate), length(cate))
    I.val = sapply(period.est, function(x){x[w[2], w[1]]})
    F.val = sapply(period.est.loo, function(x){x[w[2], w[1]]})

    I.mat = as.matrix(xtabs(val ~ ., cbind(data.frame(cate.comb), val = I.val)))
    F.mat = as.matrix(xtabs(val ~ ., cbind(data.frame(cate.comb), val = F.val)))
    if (nrow(cate.comb) > 1){
      # For univariate case, just skip below code because it will be still the same
      I.mat = I.mat + Conj(t(I.mat)) - diag(diag(I.mat))
      F.mat = F.mat + Conj(t(F.mat)) - diag(diag(F.mat))
    }


    attr(F.mat, "class") = NULL
    attr(F.mat, "call") = NULL
    F.mat = fixpd.fun(F.mat)

    return(Re(sum(diag(I.mat %*% solve(F.mat)))) +
             log(prod(Re(eigen(F.mat, only.values = T)$values)))) # Changed
  }

  # Function to compute the Whittle likelihood for a specific bandwidth
  sum.Whittle.fun = function(band.mat.row, cate.comb,
                             freq.list, freq.ext.list, freq.idx.df,
                             period.raw, period.ext){

    omega.comb = freq.list$omega.comb
    period.smooth.mat.loo.list = vector("list", nrow(cate.comb))

    for (r in 1:nrow(cate.comb)){
      # Compute the leave-one-out smoothed periodogram estimator: F matrix
      f.ij = apply(omega.comb, MARGIN = 1, smooth.fun,
                              period.mat = period.ext[[r]],
                              w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                              b1 = band.mat.row[1], b2 = band.mat.row[2],
                              loo = T)
      period.smooth.mat.loo.list[[r]] = matrix(f.ij,
                                               ncol=length(freq.list$omega1),
                                               nrow=length(freq.list$omega2), byrow=T)
    }

    # Compute the Whittle likelihood (use na.rm = T to treat NaN as zero)
    freq.idx.df[, 3] = apply(freq.idx.df[,1:2], MARGIN = 1, ind.Whittle.fun,
                             period.est=period.raw, period.est.loo=period.smooth.mat.loo.list,
                             cate.comb=cate.comb)
    return(sum(freq.idx.df[, 3], na.rm = T))
  }


  ### Whether to correct the edge effect in kernel smoothing by extending frequency grid
  if (correct){
    freq.ext.list = freq.grid.fun(A1, A2, ext.factor = 2, return.comb = T, endpt = endpt)
  }else{
    freq.ext.list = freq.list
  }


  ### Compute I matrix for each frequency
  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  DFT.list = foreach(k = seq_along(cate),
                     .export = c('period.fun.sim', 'freq.grid.fun', 'ha.fun'),
                     .packages = 'spatstat') %dopar%{
                       period.fun.sim(i = cate[k], j = cate[k], ppp = ppp,
                                  inten.fun.sep=inten.fun.sep, a = a, A1 = A1, A2 = A2,
                                  ext.factor = 2, return.DFT = T, endpt = endpt)$DFT
                     }

  names(DFT.list) = cate


  for (r in 1:nrow(cate.comb)){
    # Calculate the periodogram estimators for both extended and original frequency grids
    period.mat.ext.list[[r]] = DFT.list[[cate.comb[r,1]]]*Conj(DFT.list[[cate.comb[r,2]]]) # Changed
    period.mat.list[[r]] = period.mat.ext.list[[r]][between.fun(freq.ext.list$omega2, min(freq.list$omega2), max(freq.list$omega2)),
                                                    between.fun(freq.ext.list$omega1, min(freq.list$omega1), max(freq.list$omega1))]
  }

  ### For each bandwidth, conduct LOOCV
  cv[[3]][2,] = foreach(j = 1:nrow(band.mat),
                        .combine = 'c',
                        .export = c("Bartlett.fun","KK.fun","smooth.fun","between.fun","fixpd.fun")
                        ) %dopar% {
                           sum.Whittle.fun(band.mat.row = band.mat[j,],
                                           cate.comb = cate.comb,
                                           freq.list = freq.list,
                                           freq.ext.list = freq.ext.list,
                                           freq.idx.df = freq.idx.df,
                                           period.raw = period.mat.list,
                                           period.ext = period.mat.ext.list)
                        }
  stopCluster(cl)

  cv[[1]] = cv[[3]][1, which.min(cv[[3]][2,])]
  cv[[2]] = min(cv[[3]][2,])
  names(cv) = c("Optimal bandwidth", "Likelihood", "Result")

  if (which.min(cv[[3]][2,]) %in% c(1, length(band.range))){
   warning(paste0("The optimal bandwidth lies on the endpoint of the `band.range`. ",
                  "We suggest extending the range of `band.range` and executing the bandwidth search again."),
           call. = F)
  }
  return(cv)
}

bandwidth.cv.fun.sim2 = function(ppp, inten.fun.sep = "est", a = 0.025,
                                band.range, correct = T, A1 = NULL, A2 = A1, endpt = 1.5){

  # Crate frequency index array (freq.idx.df) to sum in likelihood calculation
  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = freq.grid.fun(A1, A2, ext.factor = NULL, return.comb = T, endpt = endpt)
  freq.idx.df = cbind(omega1.ind = rep(seq_along(freq.list$omega1), length(freq.list$omega2)),
                      omega2.ind = rep(seq_along(freq.list$omega2), each = length(freq.list$omega1)),
                      val = NA)

  # Bandwidth to evaluate
  band.mat = cbind(b1 = (A2/A1)*band.range, b2 = band.range)

  # Matrix to save the LOOCV result
  cv = vector("list", 3)
  cv[[3]] = matrix(c(band.range, rep(NA, length(band.range))),
                   ncol = length(band.range), nrow = 2, byrow = T)
  rownames(cv[[3]]) = c("Bandwidth", "Whittle likelihood")
  # Enumerate all possible combination of spectra to compute
  if (is.null(marks(ppp))){ # Univariate point pattern
    cate = 1
  }else{ # Multivariate point pattern
    cate = levels(marks(ppp))
  }
  cate.comb = cate.comb.fun(cate)
  period.mat.list = period.mat.ext.list = vector("list", nrow(cate.comb))

  # Function to compute the Whittle likelihood without the summation part
  ind.Whittle.fun = function(w, period.est, period.est.loo, cate.comb){

    I.mat = F.mat = matrix(NA, length(cate), length(cate))
    I.val = sapply(period.est, function(x){x[w[2], w[1]]})
    F.val = sapply(period.est.loo, function(x){x[w[2], w[1]]})

    I.mat = as.matrix(xtabs(val ~ ., cbind(data.frame(cate.comb), val = I.val)))
    F.mat = as.matrix(xtabs(val ~ ., cbind(data.frame(cate.comb), val = F.val)))
    if (nrow(cate.comb) > 1){
      # For univariate case, just skip below code because it will be still the same
      I.mat = I.mat + Conj(t(I.mat)) - diag(diag(I.mat))
      F.mat = F.mat + Conj(t(F.mat)) - diag(diag(F.mat))
    }

    attr(F.mat, "class") = NULL
    attr(F.mat, "call") = NULL
    F.mat = fixpd.fun(F.mat)

    return(Re(sum(diag(I.mat %*% solve(F.mat)))) +
             log(prod(Re(eigen(F.mat, only.values = T)$values)))) # Changed
  }

  # Function to compute the Whittle likelihood for a specific bandwidth
  sum.Whittle.fun = function(band.mat.row, cate.comb,
                             freq.list, freq.ext.list, freq.idx.df,
                             period.raw, period.ext){

    omega.comb = freq.list$omega.comb
    period.smooth.mat.loo.list = vector("list", nrow(cate.comb))

    for (r in 1:nrow(cate.comb)){
      # Compute the leave-one-out smoothed periodogram estimator: F matrix
      f.ij = apply(omega.comb, MARGIN = 1, smooth.fun,
                              period.mat = period.ext[[r]],
                              w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                              b1 = band.mat.row[1], b2 = band.mat.row[2],
                              loo = T)
      period.smooth.mat.loo.list[[r]] = matrix(f.ij,
                                               ncol=length(freq.list$omega1),
                                               nrow=length(freq.list$omega2), byrow=T)
    }

    # Compute the Whittle likelihood (use na.rm = T to treat NaN as zero)
    freq.idx.df[, 3] = apply(freq.idx.df[,1:2], MARGIN = 1, ind.Whittle.fun,
                             period.est=period.raw, period.est.loo=period.smooth.mat.loo.list,
                             cate.comb=cate.comb)
    return(sum(freq.idx.df[, 3], na.rm = T))
  }


  ### Whether to correct the edge effect in kernel smoothing by extending frequency grid
  if (correct){
    freq.ext.list = freq.grid.fun(A1, A2, ext.factor = 2, return.comb = T, endpt = endpt)
  }else{
    freq.ext.list = freq.list
  }


  ### Compute I matrix for each frequency
  DFT.list = vector("list", length(cate))
  for (k in seq_along(cate)){
    DFT.list[[k]] = period.fun.sim(i = cate[k], j = cate[k], ppp = ppp,
                                   inten.fun.sep=inten.fun.sep, a = a, A1 = A1, A2 = A2,
                                   ext.factor = 2, return.DFT = T, endpt = endpt)$DFT
  }
  names(DFT.list) = cate


  for (r in 1:nrow(cate.comb)){
    # Calculate the periodogram estimators for both extended and original frequency grids
    period.mat.ext.list[[r]] = DFT.list[[cate.comb[r,1]]]*Conj(DFT.list[[cate.comb[r,2]]]) # Changed
    period.mat.list[[r]] = period.mat.ext.list[[r]][between.fun(freq.ext.list$omega2, min(freq.list$omega2), max(freq.list$omega2)),
                                                    between.fun(freq.ext.list$omega1, min(freq.list$omega1), max(freq.list$omega1))]
  }

  ### For each bandwidth, conduct LOOCV
  for (j in 1:nrow(band.mat)){
    cv[[3]][2,j] = sum.Whittle.fun(band.mat.row = band.mat[j,],
                                   cate.comb = cate.comb,
                                   freq.list = freq.list,
                                   freq.ext.list = freq.ext.list,
                                   freq.idx.df = freq.idx.df,
                                   period.raw = period.mat.list,
                                   period.ext = period.mat.ext.list)
  }

  cv[[1]] = cv[[3]][1, which.min(cv[[3]][2,])]
  cv[[2]] = min(cv[[3]][2,])
  names(cv) = c("Optimal bandwidth", "Likelihood", "Result")

  if (which.min(cv[[3]][2,]) %in% c(1, length(band.range))){
    warning(paste0("The optimal bandwidth lies on the endpoint of the `band.range`. ",
                   "We suggest extending the range of `band.range` and executing the bandwidth search again."),
            call. = F)
  }
  return(cv)
}


period.2Dsmooth.fun.sim = function(ppp, i = NULL , j = i, inten.fun.sep = "est",
                                   bandwidth, correct = T, a = 0.025, A1 = NULL, A2 = A1, endpt = 1.5){
  ### Construct frequency grids
  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = freq.grid.fun(A1, A2, ext.factor = NULL, return.comb = T, endpt = endpt) # Original grid
  omega.comb = freq.list$omega.comb
  if (correct){
    freq.ext.list = freq.grid.fun(A1, A2, ext.factor = 2, return.comb = T, endpt = endpt) # Grid for edge-correction
  }else{
    freq.ext.list = freq.list
  }

  b = c(A2/A1, 1)*bandwidth



  ### List all combinations
  if (!(is.null(i) | is.null(j))){ # If both i and j are specified, compute the spectra for [i,j]
    cate.comb = cbind(i, j)
    cate = 1
  }else{ # Otherwise, enumerate all combinations of spectra to compute
    if (!is.multitype(ppp)){ # Univariate point pattern
      cate = 1
    }else{ # Multivariate point pattern
      cate = levels(marks(ppp))
    }
    cate.comb = cate.comb.fun(cate)
  }

  ### Spectra estimation
  if (length(cate) == 1){ # For univariate process, no need to use parallel computing
    # Calculate the periodogram
    period.mat.ext = period.fun.sim(i = cate.comb[1,1], j = cate.comb[1,2], ppp = ppp,
                                    inten.fun.sep = inten.fun.sep, a = a,
                                    A1 = A1, A2 = A2, ext.factor = 2, endpt = endpt)
    # Calculate the kernel smoothed spectral estimator
    val = apply(omega.comb, MARGIN = 1, smooth.fun,
                           period.mat = period.mat.ext,
                           w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                           b1 = b[1], b2 = b[2])
    # Return the kernel smoothed spectral estimator
    period.smooth.mat.list = list()
    period.smooth.mat.list[[1]] = matrix(val,
                                         ncol = length(freq.list$omega1),
                                         nrow = length(freq.list$omega2), byrow = T)
    colnames(period.smooth.mat.list[[1]]) = freq.list$omega1
    rownames(period.smooth.mat.list[[1]]) = freq.list$omega2
  }else{ # For multivariate process, use parallel computing

    cl = makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    DFT.list = foreach(k = seq_along(cate),
                       .export = c('period.fun.sim', 'freq.grid.fun', 'ha.fun'),
                       .packages = 'spatstat') %dopar%{
                         period.fun.sim(i = cate[k], j = cate[k], ppp = ppp,
                                        inten.fun.sep = inten.fun.sep, a = a,
                                        A1 = A1, A2 = A2, ext.factor = 2, return.DFT = T, endpt = endpt)$DFT
                       }
    names(DFT.list) = cate
    period.smooth.mat.list = foreach(r = 1:nrow(cate.comb),
                                     .export = c("Bartlett.fun","KK.fun","smooth.fun")) %dopar% {
                                       # Calculate the periodogram
                                       period.mat.ext = DFT.list[[cate.comb[r,1]]]*Conj(DFT.list[[cate.comb[r,2]]]) # Changed
                                       # Calculate the kernel smoothed spectral estimator
                                       val = apply(omega.comb, MARGIN = 1, smooth.fun,
                                                              period.mat = period.mat.ext,
                                                              w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                                                              b1 = b[1], b2 = b[2])
                                       # Return the kernel smoothed spectral estimator
                                       period.smooth.mat = matrix(val,
                                                                  ncol = length(freq.list$omega1),
                                                                  nrow = length(freq.list$omega2), byrow = T)
                                       colnames(period.smooth.mat) = freq.list$omega1
                                       rownames(period.smooth.mat) = freq.list$omega2
                                       period.smooth.mat
                                     }
    stopCluster(cl)
  }

  names(period.smooth.mat.list) = paste0(cate.comb[,1], ", ", cate.comb[,2])
  if (nrow(cate.comb) == 1){
    return(period.smooth.mat.list[[1]])
  }else{
    return(period.smooth.mat.list)
  }

}


# Now bandwidth is a vector
parallel.period.2Dsmooth.fun = function(ppp.list, i, j, inten.fun.sep = "est",
                                        bandwidth, correct = T, a = .025, endpt = 1.5){

  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  result = foreach(q = 1:length(ppp.list),
                   .export = c("Bartlett.fun","ha.fun","KK.fun","smooth.fun",
                               "freq.grid.fun","period.fun.sim","period.2Dsmooth.fun.sim"),
                   .packages = 'spatstat') %dopar%{
                      period.2Dsmooth.fun.sim(ppp = ppp.list[[q]], i = i, j = j,
                                              inten.fun.sep = inten.fun.sep,
                                              bandwidth = bandwidth[q],
                                              correct = correct, a = a, endpt = endpt)
                       }
  stopCluster(cl)
  return(result)
}

period.2Dsmooth.fun = function(ppp, i = NULL, j = i,
                               inten.formula = "~1", data.covariate = NULL,
                               bandwidth, correct = T, a = 0.025,
                               A1 = NULL, A2 = A1, endpt = 1.5, equal = F){
  ### Construct frequency grids
  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = freq.grid.fun(A1, A2, ext.factor = NULL, return.comb = T, endpt = endpt) # Original grid
  omega.comb = freq.list$omega.comb
  if (correct){
    freq.ext.list = freq.grid.fun(A1, A2, ext.factor = 2, return.comb = T, endpt = endpt) # Grid for edge-correction
  }else{
    freq.ext.list = freq.list
  }

  if (equal){
    b = rep(bandwidth, 2)
  }else{
    b = c(A2/A1, 1)*bandwidth
  }
  names(b) = c("b1", "b2")

  ### List all combinations
  if (!(is.null(i) | is.null(j))){ # If both i and j are specified, compute the spectra for [i,j]
    cate.comb = cbind(i, j)
    cate = 1
  }else{ # Otherwise, enumerate all combinations of spectra to compute
    if (is.null(marks(ppp))){ # Univariate point pattern
      cate = 1
    }else{ # Multivariate point pattern
      cate = levels(marks(ppp))
    }
    cate.comb = cate.comb.fun(cate)
  }

  ### Spectra estimation
  if (length(cate) == 1){ # For univariate process, no need to use parallel computing
    # Calculate the periodogram
    period.mat.ext = period.fun(i = cate.comb[1,1], j = cate.comb[1,2], ppp = ppp,
                                inten.formula = inten.formula, data.covariate = data.covariate,
                                a = a, A1 = A1, A2 = A2, ext.factor = 2)
    # Calculate the kernel smoothed periodogram
    val = apply(omega.comb, MARGIN = 1, smooth.fun,
                period.mat = period.mat.ext,
                w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                b1 = b[1], b2 = b[2])
    # Return the kernel smoothed periodogram
    period.smooth.mat.list = list()
    period.smooth.mat.list[[1]] = matrix(val,
                                         ncol = length(freq.list$omega1),
                                         nrow = length(freq.list$omega2), byrow = T)
    colnames(period.smooth.mat.list[[1]]) = freq.list$omega1
    rownames(period.smooth.mat.list[[1]]) = freq.list$omega2

    inten.fitted.list = list(period.mat.ext)

  }else{ # For multivariate process, use parallel computing
    cl = makeCluster(detectCores() - 1)
    registerDoParallel(cl)

    DFT.list = foreach(k = seq_along(cate),
                       .export = c('period.fun', 'freq.grid.fun', 'ha.fun'),
                       .packages = 'spatstat') %dopar%{
                         period.fun(i = cate[k], j = cate[k], ppp = ppp,
                                    inten.formula = inten.formula, data.covariate = data.covariate,
                                    a = a, A1 = A1, A2 = A2, ext.factor = 2, return.DFT = T, endpt = endpt)$DFT
                       }
    names(DFT.list) = cate



    period.smooth.mat.list = foreach(r = 1:nrow(cate.comb),
                                     .export = c("Bartlett.fun","KK.fun","smooth.fun")) %dopar% {
                                       # Calculate the periodogram
                                       period.mat.ext = DFT.list[[cate.comb[r,1]]]*Conj(DFT.list[[cate.comb[r,2]]])
                                       if (cate.comb[r,1] == cate.comb[r,2]){
                                         period.mat.ext = Re(period.mat.ext)
                                       }
                                       # Calculate the kernel smoothed periodogram
                                       val = apply(omega.comb, MARGIN = 1, smooth.fun,
                                                   period.mat = period.mat.ext,
                                                   w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                                                   b1 = b[1], b2 = b[2])
                                       # Return the kernel smoothed periodogram
                                       period.smooth.mat = matrix(val,
                                                                  ncol = length(freq.list$omega1),
                                                                  nrow = length(freq.list$omega2), byrow = T)
                                       colnames(period.smooth.mat) = freq.list$omega1
                                       rownames(period.smooth.mat) = freq.list$omega2
                                       period.smooth.mat
                                     }
    stopCluster(cl)

    inten.fitted.list = DFT.list
  }



  names(period.smooth.mat.list) = paste0(cate.comb[,1], ", ", cate.comb[,2])
  if (nrow(cate.comb) == 1){
    out = period.smooth.mat.list[[1]]
  }else{
    out = period.smooth.mat.list
  }

  inten.fitted.list = lapply(inten.fitted.list, attr, which = "inten.fitted")
  names(inten.fitted.list) = cate
  inten.fitted.scaled = lapply(inten.fitted.list, scale2unitarea, length.x = A1, length.y = A2)
  attr(out, "inten.fitted.scaled") = inten.fitted.scaled
  attr(out, "a") = a
  attr(out, "A1") = A1
  attr(out, "A2") = A2
  attr(out, "freq.list") = freq.list
  attr(out, "cate") = cate
  attr(out, "cate.comb") = cate.comb
  attr(out, "bandwidth") = b

  return(out)
}

bandwidth.cv.fun = function(ppp, inten.formula = NULL, data.covariate = NULL,
                            a = 0.025, band.range, correct = T, A1 = NULL, A2 = A1,
                            endpt = 1.5, equal = F){
  # Crate frequency index array (freq.idx.df) to sum in likelihood calculation
  if (is.null(A1) | is.null(A2)){
    A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  }
  freq.list = freq.grid.fun(A1, A2, ext.factor = NULL, return.comb = T, endpt = endpt)

  freq.idx.df = cbind(omega1.ind = rep(seq_along(freq.list$omega1), length(freq.list$omega2)),
                      omega2.ind = rep(seq_along(freq.list$omega2), each = length(freq.list$omega1)),
                      val = NA)
  # Bandwidth to evaluate
  if (equal){
    band.mat = cbind(b1 = band.range, b2 = band.range)
  }else{
    band.mat = cbind(b1 = (A2/A1)*band.range, b2 = band.range)
  }

  # Matrix to save the LOOCV result
  cv = vector("list", 3)
  cv[[3]] = matrix(c(band.range, rep(NA, length(band.range))),
                   ncol = length(band.range), nrow = 2, byrow = T)
  rownames(cv[[3]]) = c("Bandwidth", "Whittle likelihood")

  # Enumerate all possible combination of spectra to compute
  if (is.null(marks(ppp))){ # Univariate point pattern
    cate = 1
  }else{ # Multivariate point pattern
    cate = levels(marks(ppp))
  }
  cate.comb = cate.comb.fun(cate)
  period.mat.list = period.mat.ext.list = vector("list", nrow(cate.comb))


  # Function to compute the Whittle likelihood without the summation part
  ind.Whittle.fun = function(w, period.est, period.est.loo, cate.comb){

    I.val = sapply(period.est, function(x){x[w[2], w[1]]})
    F.val = sapply(period.est.loo, function(x){x[w[2], w[1]]})

    I.mat = as.matrix(xtabs(val ~ ., cbind(data.frame(cate.comb), val = I.val)))
    F.mat = as.matrix(xtabs(val ~ ., cbind(data.frame(cate.comb), val = F.val)))
    if (nrow(cate.comb) > 1){
      # For univariate case, just skip below code because it will be still the same
      I.mat = I.mat + Conj(t(I.mat)) - diag(diag(I.mat))
      F.mat = F.mat + Conj(t(F.mat)) - diag(diag(F.mat))
    }

    attr(F.mat, "class") = NULL
    attr(F.mat, "call") = NULL
    F.mat = fixpd.fun(F.mat)

    return(Re(sum(diag(I.mat %*% solve(F.mat)))) +
             log(prod(Re(eigen(F.mat, only.values=T)$values))))
  }

  # Function to compute the Whittle likelihood across different bandwidth
  sum.Whittle.fun = function(band.mat.row, cate.comb,
                             freq.list, freq.ext.list, freq.idx.df,
                             period.raw, period.ext){
    period.smooth.mat.loo.list = vector("list", nrow(cate.comb))
    for (r in 1:nrow(cate.comb)){
      # Compute the leave-one-out smoothed periodogram estimator: F matrix
      val.smoothed = apply(freq.list$omega.comb, MARGIN = 1, smooth.fun,
                           period.mat = period.ext[[r]],
                           w.k1 = freq.ext.list$omega1, w.k2 = freq.ext.list$omega2,
                           b1 = band.mat.row[1], b2 = band.mat.row[2],
                           loo = T)
      # F matrix
      period.smooth.mat.loo.list[[r]] = matrix(val.smoothed,
                                               ncol=length(freq.list$omega1),
                                               nrow=length(freq.list$omega2), byrow=T)
    }

    # Compute the Whittle likelihood (use na.rm = T to treat NaN as zero)
    freq.idx.df[, 3] = apply(freq.idx.df[, 1:2], MARGIN = 1, ind.Whittle.fun,
                             period.est=period.raw, period.est.loo=period.smooth.mat.loo.list,
                             cate.comb=cate.comb)
    return(sum(freq.idx.df[, 3], na.rm = T))
  }


  ### Whether to correct the edge effect in kernel smoothing by extending frequency grid
  if (correct){
    freq.ext.list = freq.grid.fun(A1, A2, ext.factor = 2, return.comb = T, endpt = endpt)
  }else{
    freq.ext.list = freq.list
  }


  ### Compute I matrix for each frequency
  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)

  DFT.list = foreach(k = seq_along(cate),
                     .export = c('period.fun', 'freq.grid.fun', 'ha.fun'),
                     .packages = 'spatstat') %dopar%{
                       period.fun(i = cate[k], j = cate[k], ppp = ppp,
                                  inten.formula = inten.formula, data.covariate = data.covariate,
                                  a = a, A1 = A1, A2 = A2, ext.factor = 2, return.DFT = T, endpt = endpt)$DFT
                     }
  names(DFT.list) = cate

  for (r in 1:nrow(cate.comb)){
    # Calculate the periodogram for both extended and original frequency grids
    period.mat.ext.list[[r]] = DFT.list[[cate.comb[r,1]]]*Conj(DFT.list[[cate.comb[r,2]]])
    if (cate.comb[r,1] == cate.comb[r,2]){
      period.mat.ext.list[[r]] = Re(period.mat.ext.list[[r]])
    }
    period.mat.list[[r]] = period.mat.ext.list[[r]][between.fun(freq.ext.list$omega2, min(freq.list$omega2), max(freq.list$omega2)),
                                                    between.fun(freq.ext.list$omega1, min(freq.list$omega1), max(freq.list$omega1))]
  }

  ### For each bandwidth, conduct LOOCV
  cv[[3]][2,] = foreach(j = 1:nrow(band.mat),
                        .combine = 'c',
                        .export = c("Bartlett.fun","KK.fun","smooth.fun","between.fun","fixpd.fun")
                        ) %dopar% {
                           sum.Whittle.fun(band.mat.row = band.mat[j,],
                                           cate.comb = cate.comb,
                                           freq.list = freq.list,
                                           freq.ext.list = freq.ext.list,
                                           freq.idx.df = freq.idx.df,
                                           period.raw = period.mat.list,
                                           period.ext = period.mat.ext.list)
                        }
  stopCluster(cl)

  cv[[1]] = cv[[3]][1, which.min(cv[[3]][2,])]
  cv[[2]] = min(cv[[3]][2,])
  names(cv) = c("Optimal bandwidth", "Likelihood", "Result")

  if (which.min(cv[[3]][2,]) %in% c(1, length(band.range))){
   warning(paste0("The optimal bandwidth lies on the endpoint of the `band.range`. ",
                  "We suggest extending the range of `band.range` and executing the bandwidth search again."),
           call. = F)
  }
  return(cv)
}



###### This function is used to compute H matrices inside the inverse Fourier transform of L2
Hmatrix = function(sp.est, ppp){
  # sp.est: a list of kernel spectral estimate matrices of the original point process (ppp)
  # ppp: point pattern

  inten.fitted.scaled = attr(sp.est,"inten.fitted.scaled")
  a = attr(sp.est,"a")
  cate = attr(sp.est,"cate")
  cate.comb = attr(sp.est,"cate.comb")

  ppp = shift(ppp, origin = "centroid")
  A1 = diff(ppp$window$xrange); A2 = diff(ppp$window$yrange)
  ppp.li = split.ppp(ppp)

  H_h2 = (a*(5/(4*(pi^2)) - 4/3) + 1)^2
  H_h2l.val = matrix(0, nrow = length(cate), ncol = length(cate), dimnames = list(cate, cate))
  H_h2ll.val = matrix(NA, nrow = length(cate), ncol = length(cate), dimnames = list(cate, cate))

  for (k in seq_along(cate)){
    H_h2l.val[k,k] = sum((ha.fun(coords.ppp(ppp.li[[cate[k]]])$x/A1, a)^2)*(ha.fun(coords.ppp(ppp.li[[cate[k]]])$y/A2, a)^2))/area(ppp)
  }

  for (r in 1:nrow(cate.comb)){
    Xi = cate.comb[r, 1]
    Xj = cate.comb[r, 2]
    Xi.x = coords.ppp(ppp.li[[Xi]])$x/A1
    Xi.y = coords.ppp(ppp.li[[Xi]])$y/A2
    Xj.x = coords.ppp(ppp.li[[Xj]])$x/A1
    Xj.y = coords.ppp(ppp.li[[Xj]])$y/A2

    H_h2ll.val[Xi, Xj] = 0.5*(sum((ha.fun(Xi.x, a)^2)*(ha.fun(Xi.y, a)^2)*inten.fitted.scaled[[Xj]](x = Xi.x, y = Xi.y)) +
                                sum((ha.fun(Xj.x, a)^2)*(ha.fun(Xj.y, a)^2)*inten.fitted.scaled[[Xi]](x = Xj.x, y = Xj.y)))/area(ppp)
    H_h2ll.val[Xj, Xi] = H_h2ll.val[Xi, Xj]
  }

  return(list(H_h2 = H_h2, H_h2l.val = H_h2l.val, H_h2ll.val = H_h2ll.val))
}





###### Function to estimate the spectrum of intensity reweighted process
IR_spec = function(w1, w2, sp.est, i = NULL, j = NULL, H.list = NULL, ppp){
  # w1, w2: frequency vector (only allow frequency values evaluated in `sp.est`)
  # sp.est: a list of kernel spectral estimate matrices of the original point process
  # i, j: pair index (optional)
  # H.list: the result from `Hmatrix()` (optional)
  # ppp: point pattern (this argument is only required if you don't specify `H.list`)

  sp.est.full = sp.est
  A1 = attr(sp.est,"A1")
  A2 = attr(sp.est,"A2")
  cate = attr(sp.est,"cate")
  cate.comb = attr(sp.est,"cate.comb")
  omega1 = attr(sp.est,"freq.list")$omega1
  omega2 = attr(sp.est,"freq.list")$omega2
  inten.fitted.scaled = attr(sp.est,"inten.fitted.scaled")

  if (is.null(i) + is.null(j) == 2){ # Return a m by m matrix, computing all pairwise local spectra
    i = j = seq_along(cate)
  }else if (is.null(i) + is.null(j) == 1){
    stop("Please specify the index for the other process, or leave both i and j unspecified.", call. = F)
  }else{ # Only compute the local spectrum of the (i, j)-th pair
    cate = unique(c(i, j))
    subset.comb.idx = apply(cate.comb, 1, function(x) all(x %in% c(i,j)))
    cate.comb = cate.comb[subset.comb.idx, , drop = F]
    sp.est = sp.est[subset.comb.idx]
  }

  ### Estimate the inverse Fourier transform of L2
  ## F.hat matrix given a frequency pair (w1, w2)
  w.idx = c(which(omega1 == w1), which(omega2 == w2))
  F.val = sapply(sp.est, function(x){x[w.idx[2], w.idx[1]]})
  F.mat = as.matrix(xtabs(val ~ ., cbind(data.frame(cate.comb), val = F.val)))
  F.mat = F.mat + Conj(t(F.mat)) - diag(diag(F.mat))
  attr(F.mat, "class") = NULL; attr(F.mat, "call") = NULL
  ## Obtain the subcomponents in L2 if they were not supplied
  if (is.null(H.list)){
    H.list = Hmatrix(sp.est = sp.est.full, ppp = ppp)
  }
  ## Inverse Fourier transform of L2
  invFT.L2.est = (H.list$H_h2*F.mat[i, j] - (2*pi)^(-2)*H.list$H_h2l.val[i, j])/H.list$H_h2ll.val[i, j]


  ### Return the result, which is a function of w1, w2
  unit.mat = diag(1, length(inten.fitted.scaled))
  rownames(unit.mat) = colnames(unit.mat) = names(inten.fitted.scaled)
  out = (2*pi)^(-2)*unit.mat[i,j] + invFT.L2.est

  attr(out, "cate") = names(inten.fitted.scaled)
  attr(out, "eval.freq") = c(w1, w2)

  return(out)
}



###### Function to calculate the (partial) coherence matrix for a given frequency
coh.by.freq = function(w1, w2, sp.est, type = "partial", i = NULL, j = NULL,
                       sp.IR = NULL, H.list = NULL, ppp = NULL){
  # w1, w2: frequency vector (only allow frequency values evaluated in `sp.est`)
  # sp.est: a list of kernel spectral estimate matrices of the point process
  # partial: if TRUE, compute partial coherence; otherwise, compute coherence
  # i, j: pair index (optional)
  # sp.IR: spectrum estimate of the IR process (if this is provided, then w1, w2, H.list, ppp are not required)
  # H.list: the result from `Hmatrix()` (this argument is only required if you don't specify `sp.IR`)
  # ppp: point pattern (this argument is only required if you don't specify `H.list`)


  if (is.null(sp.IR)){
    sp.IR = IR_spec(w1, w2, sp.est, H.list = H.list, ppp = ppp)
    sp.IR = fixpd.fun(sp.IR)
    nom = switch(type,
                 partial = solve(sp.IR),
                 normal = sp.IR,
                 stop("`type` should be 'partial' or 'normal'.", call. = F))
    cate = attr(sp.est, "cate")
  }else if (is.matrix(sp.IR)){
    sp.IR = fixpd.fun(sp.IR)
    nom = switch(type,
                 partial = solve(sp.IR),
                 normal = sp.IR,
                 stop("`type` should be 'partial' or 'normal'.", call. = F))
    cate = attr(sp.IR, "cate")
  }else{
    stop("`sp.IR` should be a matrix.", call. = F)
  }

  D = Mod(diag(1/diag(nom)))
  coh.mat = D %*% Mod(nom)^2 %*% D

  colnames(coh.mat) = cate
  row.names(coh.mat) = cate

  if (is.null(i) + is.null(j) == 2){ # Return a m by m matrix
    return(coh.mat)
  }else if (is.null(i) + is.null(j) == 1){
    stop("Please specify the index for the other process, or leave both i and j unspecified.",
         call. = F)
  }else{ # Only compute the coherence of the (i, j)-th pair
    return(coh.mat[i, j])
  }
}




# H.list = Hmatrix(period.smooth.mat, X)
# cate = attr(period.smooth.mat,"cate")
# omega1 = attr(period.smooth.mat,"freq.list")$omega1
# omega2 = attr(period.smooth.mat,"freq.list")$omega2
#
# sp.IR = IR_spec(w1 = omega1[1], w2 = omega2[1], sp.est = period.smooth.mat, H.list = H.list)
# coh.partial = coh.by.freq(sp.IR = sp.IR, H.list = H.list)
# coh.normal = coh.by.freq(sp.IR = sp.IR, H.list = H.list, type = "normal")




###### Function to calculate the maximum (partial) coherence matrix across frequency
coh.max.fun = function(sp.est, ppp, type = "partial", all = F){
  # sp.est: a list of kernel spectral estimate matrices of the point process
  # ppp: point pattern
  # all: if TRUE, extract the maximum (partial) coherence across all frequencies, which is not recommended

  H.list = Hmatrix(sp.est, ppp = ppp)
  cate = attr(sp.est, "cate")
  omega.comb = attr(sp.est, "freq.list")$omega.comb[, 1:2]

  if (all){ # Compute (partial) coherence for all frequencies.
    cat(paste0("Use all", nrow(omega.comb), "frequencies to pick the maximum",
               ifelse(type == "partial", " partial", ""),
               " coherence."))
  }else{ # Compute (partial) coherence only for frequencies which share no overlapping smoothing neighbors
    omega1 = attr(sp.est, "freq.list")$omega1
    omega2 = attr(sp.est, "freq.list")$omega2
    b = attr(sp.est, "bandwidth")
    origin.idx = c((length(omega1)+1)/2, (length(omega2)+1)/2)
    omega.diff = c(omega1[2]-omega1[1], omega2[2]-omega2[1])
    idx.diff = ceiling(2*b/omega.diff)
    omega1 = omega1[c(rev(seq(from = origin.idx[1], to = 1, by = -idx.diff[1])[-1]),
                      seq(from = origin.idx[1], to = length(omega1), by = idx.diff[1]))]
    omega2 = omega2[c(rev(seq(from = origin.idx[2], to = 1, by = -idx.diff[2])[-1]),
                      seq(from = origin.idx[2], to = length(omega2), by = idx.diff[2]))]
    omega.comb = cbind(omega1 = rep(omega1, length(omega2)),
                       omega2 = rep(omega2, each = length(omega1)))

    cat(paste0("Number of frequencies to pick the maximum",
               ifelse(type == "partial", " partial", ""),
               " coherence:"),
        nrow(omega.comb),
        "(this value should not be too small).")
  }


  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  out = foreach(r = 1:nrow(omega.comb), .combine = "rbind",
                .export = c("IR_spec", "coh.by.freq", "fixpd.fun"),
                .inorder = T) %dopar% {
                  sp.IR = IR_spec(w1 = omega.comb[r,1], w2 = omega.comb[r,2],
                                  sp.est = sp.est, H.list = H.list)
                  coh.mat = coh.by.freq(sp.IR = sp.IR, H.list = H.list, type = type)
                  return(coh.mat[upper.tri(coh.mat, diag = F)])
                }
  stopCluster(cl)


  coh.w = cbind(omega.comb, out)
  max.coh = apply(coh.w[, -c(1:2), drop = F], MARGIN = 2, max)
  max.coh.mat = diag(1, length(cate))
  max.coh.mat[upper.tri(max.coh.mat, diag = F)] = max.coh
  max.coh.mat = max.coh.mat + t(max.coh.mat) - diag(diag(max.coh.mat))
  row.names(max.coh.mat) = colnames(max.coh.mat) = cate

  attr(max.coh.mat, 'CohTable') = coh.w
  class(max.coh.mat) = "SuppAttr" # Define a new class for customized print method which does not print attribute

  return(max.coh.mat)
}

print.SuppAttr = function(x, ...){
  print(matrix(as.numeric(x),
               nrow = attributes(x)$dim[1],
               ncol = attributes(x)$dim[2],
               dimnames = list(attributes(x)$dimnames[[1]],
                               attributes(x)$dimnames[[2]])),
        ...)
}


### Kernel functions
bkernels = list()
# Gaussian kernel
bkernels$Thomas = function(r, bandwidth, ...){
  exp(- r^2/(2 * bandwidth^2)) / (2 * pi * bandwidth^2)
}
# Variance-Gamma (Bessel) kernel with bandwidth and shape parameter nu.ker
bkernels$VarGamma = function(r, bandwidth, nu.ker){
  stopifnot(nu.ker > -1/2)
  const = 1 / (4 * pi * nu.ker * bandwidth^2)
  u = r/bandwidth
  u = ifelse(u > 0, (u^nu.ker) * besselK(u, nu.ker) / (2^(nu.ker - 1) * gamma(nu.ker)), 1)
  return(abs(const * u))
}
# Cauchy kernel
bkernels$Cauchy = function(r, bandwidth, ...){
  ((1 + (r / bandwidth)^2)^(-1.5)) / (2 * pi * bandwidth^2)
}

### c_l function in Jalilian et al. (2015) p. 1024
c.fun = list()
c.fun$Thomas = function(r, bandwidth, kappa, ...){
  (kappa/(bkernels[["Thomas"]](0, bandwidth)^2)) * exp(-r^2/(4*bandwidth^2))/(4*pi*bandwidth^2*kappa)
}
c.fun$VarGamma = function(r, bandwidth, kappa, nu.ker){
  (kappa/(bkernels[["VarGamma"]](0, bandwidth, nu.ker))^2) * besselK(r/bandwidth, 2*nu.ker+1) * (r/bandwidth)^(2*nu.ker+1) / (pi*bandwidth^2*2^(2*nu.ker+2)*gamma(2*nu.ker+2)*kappa)
}
c.fun$Cauchy = function(r, bandwidth, ...){
  (kappa/(bkernels[["Cauchy"]](0, bandwidth))^2) * (1+r^2/(4*bandwidth^2))^(-1.5) / (8*pi*bandwidth^2*kappa)
}

### Pair correlation function of multivariate product-shot-noise Cox process
g.fun = function(r, i, j, kappa, xi, kernels, bandwidth, nu.ker=NULL){

  num.parent = length(kappa)

  if (i > num.parent || j > num.parent){
    stop("Index out of bound.", call. = F)
  }
  if (identical(dim(xi), rep(num.parent, 2)) == F){
    stop(paste("xi should be a", num.parent,"by", num.parent, "matrix."), call. = F)
  }
  if (length(kernels) != num.parent || length(bandwidth) != num.parent){
    stop("Please specify a kernel function, and its parameter(s), for each parent process.",
         call. = F)
  }


  exp.term = function(i, j, num.parent, kappa, xi, kernels, bandwidth, nu.ker=NULL){
    sum.term = 0
    for (l in (1:num.parent)[-c(i,j)]){
      sum.term = sum.term + kappa[l]*xi[l,i]*xi[l,j]*c.fun[[kernels[l]]](r, bandwidth[l], kappa[l], nu.ker[l])
    }
    return(exp(sum.term))
  }

  # Return pair correlation function if i = j
  if (i == j){
    first.term = 1 + (bkernels[[kernels[i]]](0, bandwidth[i], nu.ker[i])^2 *
                        c.fun[[kernels[i]]](r, bandwidth[i], kappa[i], nu.ker[i]))/kappa[i]
  }
  # Return cross-pair correlation function if i =/= j
  else{
    first.term = (1+xi[i,j]*bkernels[[kernels[i]]](0, bandwidth[i], nu.ker[i])*c.fun[[kernels[i]]](r, bandwidth[i], kappa[i], nu.ker[i])) *
  (1+xi[j,i]*bkernels[[kernels[j]]](0, bandwidth[j], nu.ker[j])*c.fun[[kernels[j]]](r, bandwidth[j], kappa[j], nu.ker[j]))
  }

  return(first.term * exp.term(i, j, num.parent, kappa, xi, kernels, bandwidth, nu.ker))

}

### Theoretical pseudo-spectra of multivariate product-shot-noise Cox process
f.fun = function(w, i, j, inten.fun.sep, a, kappa, xi, kernels, bandwidth, nu.ker=NULL){

  # Function of the Hankel transform
  Hankel.int = Vectorize( # Vectorize with respect to w
    function(w, i, j, kappa, xi, kernels, bandwidth, nu.ker=NULL){
      integrate(function(x, w, i, j, kappa, xi, kernels, bandwidth, nu.ker=NULL)
                x*besselJ(w*x, 0)* (g.fun(x, i=i, j=j, kappa=kappa, xi=xi,
                                          kernels=kernels, bandwidth=bandwidth,
                                          nu.ker=nu.ker) - 1),
                lower=0, upper=Inf, w=w, i=i, j=j, kappa=kappa, xi=xi,
                kernels=kernels, bandwidth=bandwidth, nu.ker=nu.ker)$value
    },
    vectorize.args = c("w")
  )
  Hankel = Hankel.int(w, i, j, kappa, xi, kernels, bandwidth, nu.ker=nu.ker)

  H_h2 = (a*(5/(4*(pi^2)) - 4/3)+1)^2
  # H_h2 = integrate(function(x, a) ha.fun(x,a)^2, lower = -0.5, upper = 0.5, a=a)$value^2

  if (all(sapply(inten.fun.sep, is.numeric))){ # If the model is homogeneous (constant intensities for all species)

    if (i == j){ # If i = j, return marginal pseudo-spectrum
      return((inten.fun.sep[[i]]*(1/(2*pi) + inten.fun.sep[[i]]*Hankel))/(2*pi))
    }else{ # If i =/= j, return cross pseudo-spectrum
      return(inten.fun.sep[[i]]*inten.fun.sep[[j]]*Hankel/(2*pi))
    }

  }else{
    H_h2l2 = integrate(function(x, inten.fun.sep, a, x_i) ha.fun(x,a)^2*inten.fun.sep[[i]](x,x_i)*inten.fun.sep[[j]](x,x_i),
                     lower=-0.5, upper=0.5, inten.fun.sep=inten.fun.sep, a=a, x_i=1)$value *
    integrate(function(x, inten.fun.sep, a, x_i) ha.fun(x,a)^2*inten.fun.sep[[i]](x,x_i)*inten.fun.sep[[j]](x,x_i),
              lower=-0.5, upper=0.5, inten.fun.sep=inten.fun.sep, a=a, x_i=2)$value

    if (i == j){
      H_h2l1 = integrate(function(x, inten.fun.sep, a, x_i) ha.fun(x,a)^2*inten.fun.sep[[i]](x,x_i),
                         lower=-0.5, upper=0.5, inten.fun.sep=inten.fun.sep, a=a, x_i=1)$value *
        integrate(function(x, inten.fun.sep, a, x_i) ha.fun(x,a)^2*inten.fun.sep[[i]](x,x_i),
                  lower=-0.5, upper=0.5, inten.fun.sep=inten.fun.sep, a=a, x_i=2)$value
      return((2*pi)^(-2)*H_h2l1/H_h2 + (H_h2l2*Hankel)/(H_h2*2*pi))
    }
    else{
      return((H_h2l2*Hankel)/(H_h2*2*pi))
    }
  }
}

# For below two functions:
# est.spec.li: a list of matrices, each matrix is a replication of spectrum estimate
# true.spec: the true spectrum matrix

IBIAS.fun = function(est.spec.li, true.spec, remove = NULL){
  numer = Reduce("+", est.spec.li) / length(est.spec.li)
  if (!is.null(remove)){
    stopifnot("The first value should be smaller than the second one." = remove[1] < remove[2])
    omega1 = attr(true.spec, "omega1")
    omega2 = attr(true.spec, "omega2")
    remove.w1.bol = between.fun(omega1, remove[1], remove[2])
    remove.w2.bol = between.fun(omega2, remove[1], remove[2])
    numer[remove.w1.bol, remove.w2.bol] = NA
  }
  num.remove = ifelse(is.null(remove), 0, sum(remove.w1.bol)*sum(remove.w2.bol))
  return(sum((numer/true.spec - 1)^2, na.rm = T)/(prod(dim(numer)) - num.remove))
}

IMSE.fun = function(est.spec.li, true.spec, remove = NULL){
  RQE = lapply(est.spec.li, function(x) (x/true.spec - 1)^2)
  RQE.sum = Reduce("+", RQE)
  if (!is.null(remove)){
    stopifnot("The first value should be smaller than the second one." = remove[1] < remove[2])
    omega1 = attr(true.spec, "omega1")
    omega2 = attr(true.spec, "omega2")
    remove.w1.bol = between.fun(omega1, remove[1], remove[2])
    remove.w2.bol = between.fun(omega2, remove[1], remove[2])
    RQE.sum[remove.w1.bol, remove.w2.bol] = NA
  }
  num.remove = ifelse(is.null(remove), 0, sum(remove.w1.bol)*sum(remove.w2.bol))
  return(sum(RQE.sum, na.rm = T)/(length(est.spec.li)*(prod(dim(RQE.sum)) - num.remove)))
}
