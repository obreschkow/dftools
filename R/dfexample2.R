#' Example of fitting a galaxy mass function in the presence of large-scale structure
#' 
#' This function generates n mock galaxies with masses, mass errors and distances. The mock sample is drawn from a custom Schechter function, predefined selection function and predefined large-scale structure density \code{f(r)}. The function then runs \code{\link{dffit}} to recover the inpu Schechter function from these noisy and LSS-biased data.
#'
#' @importFrom graphics legend text
#'  
#' @param seed A positive integer used as seed for the random number generator.
#' @param sigma Gaussian observing error in log10(mass).
#'
#' @examples
#' dfexample2()
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export


dfexample2 = function(seed = 1, sigma = 0.3) {

  # survey model
  rmax = 100
  p = c(-2,10,-1.3)
  dmax = function(x) 1e-3*sqrt(10^x)
  f = function(x,r) pracma::erf((dmax(x)-r)/dmax(x)*20)*0.5+0.5
  dVdr = function(r) 0.2*r^2
  g = function(r) 1+0.9*sin((r/100)^0.6*2*pi)
  
  # make data
  cat('Generate mock data with observing errors and large-scale structure (LSS):\n')
  dat = dfmockdata(seed = seed, sigma = sigma, p = p, f = f, g = g, dVdr = dVdr, rmax = rmax)
  dat <<- dat
  cat(sprintf('=> %d galaxies\n',dat$n))
  
  # fit data
  selection = list(f,dVdr,0,rmax)
  cat('Fit mock data without any bias correction:\n')
  survey1 = dffit(dat$x,selection,NULL,r=dat$r,correct.lss.bias = FALSE)
  cat('Fit mock data while correcting for observational errors (Eddington bias):\n')
  survey2 = dffit(dat$x,selection,dat$x.err,r=dat$r,correct.lss.bias = FALSE)
  cat('Fit mock data while correcting for observational errors and LSS:\n')
  survey3 = dffit(dat$x,selection,dat$x.err,r=dat$r,correct.lss.bias = TRUE)
  
  # plot effective volumes
  x = seq(6,12,0.01)
  dfplotveff(survey3,xlab='log10(M/Msun)',legend=FALSE)
  lines(x,dat$veff.lss(x),type='l',ylim=c(0,6),col='black')
  legend('bottomright',c('Recovered model from data used for fit','Input model with LSS used to generate the data','Input model without without LSS'),
         lwd=c(2,2,2),col=c('blue','black','red'),bty='n')
  
  # plot MFs
  mfplot(survey1,xlim=c(1e7,2e11),ylim=c(2e-5,5),nbins=20,bin.xmin=7,bin.xmax=11,col.fit='orange',col.data.input='orange')
  mfplot(survey2,xlim=c(1e7,2e11),ylim=c(2e-5,5),nbins=20,bin.xmin=7,bin.xmax=11,col.fit='red',col.data.posterior='red',show.input.data=FALSE,add=TRUE)
  mfplot(survey3,xlim=c(1e7,2e11),ylim=c(2e-5,5),nbins=20,bin.xmin=7,bin.xmax=11,col.fit='blue',col.data.posterior='blue',show.input.data=FALSE,add=TRUE)
  lines(10^x,survey1$model$gdf(x,p),lty=2)
  legend('topright',c('Fit without bias corrections','Correcting observational errors',
                      'Correcting observational errors + LSS','True input model'),
         lwd=c(2,2,2,1),lty=c(1,1,1,2),col=c('orange','red','blue','black'),bty='n')
  x = 8; y = 2e-3
  points(10^x,y,pch=20)
  segments(10^(x-sigma),y,10^(x+sigma),y,lwd=1)
  text(10^x,y*1.6,expression('Assumed observing error'~sigma),cex=0.8)
  
}