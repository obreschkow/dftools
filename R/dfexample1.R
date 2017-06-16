#' Example of using dftools
#'
#' This function generates n mock galaxies with observing errors, drawn from a Schechter function with predefined parameters, using a predefined effective volume function. It then uses \code{\link{dffit}} to recover the input function.
#'
#' @param n Number of gtalaxies
#' @param seed Random number seed
#' @param sigma Observational uncertainty, defined as standard-deviation in log10(Mass)
#' @param p.true Input parameters for the Schechter function
#' @param include.mrp If \code{TRUE}, an MRP function will also be fitted
#'
#' @examples
#' mfexample()
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfexample1 <- function(n = 5000, seed = 3, sigma = 0.5, p.true = c(-2,9,-1.3), include.mrp = FALSE) {

  # user parameters
  set.seed(seed)
  veff.scale = function(x) {1e-9*10^(0.8*x)}

  # make rescaled Veff(x)
  xmin = 2
  xmax = 13
  dx = 0.01
  x = seq(xmin,xmax,dx)
  phi.true = dfmodel(x, p.true)
  f = n/sum(phi.true*veff.scale(x)*dx)
  veff.fn <- function(x) {veff.scale(x)*f}
  veff = veff.fn(x)

  # resample
  cat(sprintf('Generate %d galaxies, with observing Gaussian errors\n',n))
  rho.true = phi.true*dx*veff
  cum = cumsum(rho.true/sum(rho.true))
  x.obs = array(1e99,n) # = log10(mass)
  xmin = min(x)
  xmax = max(x)
  r = runif(n,min=0.0)
  for (i in seq(n)) {
    while (x.obs[i]<xmin | x.obs[i]>xmax) {
      index = which.min(abs(cum-r[i]))
      x.obs[i] = x[index]+rnorm(1)*sigma
    }
  }

  # fit
  cat('Fit a Schechter function to the observed data (grey points)\n')
  bundle = dffit(x.obs, veff.fn, rep(sigma,n), write.fit = T, xmin = xmin, xmax = xmax, dx = dx, p.initial = p.true)
  p.fit = bundle$fit$parameters$p.optimal

  # make posterior masses
  cat('Determine posterior masses\n')
  bundle = dfposteriors(bundle)
  
  # MRP
  if (include.mrp) {
    cat('Fit a MRP function (4 parameter model) to the observed data\n')
    bundle2 = dffit(x.obs, veff.fn, rep(sigma,n), write.fit = T, xmin = xmin, xmax = xmax, dx = dx, phi = 'MRP', p.initial = c(p.true,1))
  }
  
  # plot covariances
  dfplotcov(bundle,p.true,title='Schechter parameter covariances')
  if (include.mrp) dfplotcov(bundle2,c(p.true,1),line.color='red',point.color='#ffbbbb',title='MRP parameter covariances')

  # plot main plot
  nbins = max(4,round(sqrt(n)/2))
  mfplot(bundle, xlim=c(1e4,2e10), ylim=c(1e-4,2), nbins = nbins, bin.xmin=4.8, bin.xmax=10, col.data = 'grey')
  mfplot(bundle, nbins = nbins, bin.xmin=4.8, bin.xmax=10, bin.type = 3, add = TRUE, show.uncertainties = FALSE, col.data='black')
  lines(10^x,pmax(1e-10,dfmodel(x,p.true)),col='black',lty=2)
  lines(10^x,veff/max(veff)*200,col='orange')
  if (include.mrp) {lines(10^x,bundle2$fit$functions$phi.fit(x),col='red')}

  # legend
  y = 1.5e-2
  x = 5.5
  points(10^x,y,pch=20)
  segments(10^(x-sigma),y,10^(x+sigma),y,lwd=1)
  text(10^x,y*1.6,expression('Assumed observing error'~sigma),cex=0.8)
  legend(2e4,6e-3,c('Input Veff (arbitrary vertical units)',
                    'Input Schechter function',
                    'Fitted Schechter function with 68%-uncertainties',
                    'Binned input data',
                    'Binned posterior data',
                    'Note: Horizontal bars are bin widths.',
                    'Vertical bars are Poission errors (correlated for purple points).'),
         bty='n',lwd=1.5,col=c('orange','black','blue','grey','black',NA,NA),
         pch=c(NA,NA,NA,20,20),
         lty=c(1,2,1,1,1,NA,NA),cex=0.8)

  # console output
  cat(sprintf('Total computation time       = %.2fs\n',bundle$fit$status$walltime.total))
  cat(sprintf('Computation time for fitting = %.2fs\n',bundle$fit$status$walltime.fitting))
  
  invisible(bundle)
}
