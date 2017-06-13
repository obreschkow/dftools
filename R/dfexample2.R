#' Example of using mftools
#'
#' This function generates n mock galaxies with observing errors, drawn from a predefined mass function and effective observing volume. It then uses \code{\link{dffit}} to recover the input function
#'
#' @param n Number of gtalaxies
#'
#' @examples
#' dfexample2()
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfexample2 = function(p.true = c(-2,10,-1.3,2/3,7,0.3), seed = 3,
                      sigma.x = 0.25, sigma.y = 0.1, corr.xy = 0.7) {
  
  # initialize
  set.seed(seed)
  mrange = c(8,12)
  jrange = c(1,4)
  cov = rbind(c(sigma.x^2,corr.xy*sigma.x*sigma.y),
              c(corr.xy*sigma.x*sigma.y,sigma.y^2))
  
  # define 2D distribution function in the (x1,x2)-plane, where
  # x1 = log10(baryon mass [Msun])
  # x2 = log10(baryon specific angular momentum [kpc km/s])
  phi = function(x,p) {
    mu = 10^(x[,1]-p[2])
    phi_m = log(10)*10^p[1]*mu^(p[3]+1)*exp(-mu)
    dj = x[,2]-p[4]*(x[,1]-p[5]) # offset from mean M-j relation
    phi_j = exp(-dj^2/2/p[6]^2)
    return(phi_m*phi_j)
  }
  
  # define 2D effective volume function
  library('pracma')
  veff = function(x) {
    return(2e-10*10^(1.5*x[,1])*(pracma::erf((x[,2]-3)*1)+1))
  }
  
  # expected source count density function
  sc = function(x) {
    return(phi(x,p.true)*veff(x))
  }
  
  # compute expected number of galaxies
  sc.array = function(x,y) {
    s = array(NA,length(x))
    for (i in seq(length(x))) {
      s[i] = sc(cbind(x[i],y[i]))
    }
    return(array(s,dim(x)))
  }
  n = simpson2d(sc.array,6,14,-1,6)
  if (n>1) {n = round(n)}
  cat(sprintf('Expected number of galaxies = %d\n',n))
  if (n<10) stop('Expected number of galaxies cannot be smaller than 10.')
  if (n>1e5) stop('Expected number of galaxies cannot be larger than 10^5.')
  
  # check if data lies nicely within boundaries
  if (simpson2d(sc.array,mrange[1],mrange[2],jrange[1],jrange[2])<0.9*n) {
    stop('With this parameter choice, more than 10% of the expected sources lie outside the integration range.')
  }
  
  # make observing errors
  x.err = aperm(array(cov,c(2,2,n)),c(3,2,1))
  
  # find maximum of sc(x)
  fminfct = function(x) {
    x = rbind(x)
    return(-sc(x))
  }
  opt = optim(c(10,3),fminfct,method='L-BFGS-B',lower=c(mrange[1],jrange[1]),upper=c(mrange[2],jrange[2]))
  xmax = rbind(opt$par)
  scmax = -opt$value
  
  # draw random sample
  x = array(NA,c(n,2))
  count = 0
  while (count<n) {
    xtry = rbind(c(runif(1,mrange[1],mrange[2]),runif(1,jrange[1],jrange[2])))
    if (sc(xtry)>runif(1,0,scmax)) {
      count = count+1
      x[count,] = MASS::mvrnorm(1,xtry,x.err[count,,])
    }
  }
  
  # fit model to sample
  cat('Fit six parameters:\n')
  if (corr.xy==0) {x.err = cbind(rep(sigma.x,n),rep(sigma.y,n))} # just to try if everything works well in this mode, too
  df = dffit(x, veff, x.err, phi = phi, p.initial = p.true, x.grid = list(seq(mrange[1],mrange[2],0.1),
                                                                 seq(jrange[1],jrange[2],0.1)))
  
  # plot parameter covariances and initial parameters
  dfplotcov(df, p = p.true)
  
  # plot DF
  dfplot2(df, p.ref = p.true, xlab = 'log10(Mass/Msun)', ylab = 'log10(j/[kpc km/s])')
  legend(1.3e8,8e3,c('Mock observations','Observational uncertainties','Input distribution function (DF)','Fitted DF','Source count model from fitted DF'),
         pch=c(20,NA,NA,15,NA),col=c('black','grey','red','blue','purple'),lty=c(NA,1,2,2,1),bty='n')
  
  # console output
  cat(sprintf('Total computation time       = %.2fs\n',df$fit$status$walltime.total))
  cat(sprintf('Computation time for fitting = %.2fs\n',df$fit$status$walltime.fitting))
  
  invisible(df)
  
}