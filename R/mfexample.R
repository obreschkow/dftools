#' Example of using mftools
#'
#' This function generates n mock galaxies with observing errors, drawn from a predefined mass function and effective observing volume. It then uses \code{\link{fit.df}} to recover the input function
#'
#' @param n Number of gtalaxies
#'
#' @examples
#' example()
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

mfexample <- function(n = 25000, seed = 0, sigma = 0.5, p.true = c(-2,9,-1.3)) {

  # user parameters
  set.seed(seed)
  veff.scale = function(x) {1e-9*(10^x)}

  # make rescaled Veff(x)
  dx = 0.01
  x = seq(2,12.5,dx)
  phi.true = mfmodel(x, p.true)
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
  mf = fit.df(x.obs, veff.fn, rep(sigma,n), log = TRUE, write.fit = T, log.integration.range = x, p.initial = p.true)
  p.fit = mf$fit$parameters$p.optimal

  # make posterior masses
  cat('Determine posterior masses\n')
  mf = mfposterior(mf)

  # plot
  nbins = max(4,round(sqrt(n)/2))
  mfplot(mf, xlim=c(1e4,2e10), ylim=c(1e-4,2), nbins = nbins, bin.xmin=4.8, bin.xmax=10, col.bin='grey')
  mfplot(mf, nbins = nbins, bin.xmin=4.8, bin.xmax=10, bin.type = 2, add = TRUE, show.uncertainties = FALSE)
  mfplot(mf, nbins = nbins, bin.xmin=4.8, bin.xmax=10, bin.type = 3, add = TRUE, show.uncertainties = FALSE, col.bin='purple')
  lines(10^x,pmax(1e-10,mfmodel(x,p.true)),col='black',lty=2)
  lines(10^x,veff/max(veff)*200,col='orange')

  # MRP
  cat('Fit a MRP function (4 parameter model) to the observed data\n')
  mf2 = fit.df(x.obs, veff.fn, rep(sigma,n), log = TRUE, write.fit = T, log.integration.range = x, mass.function = 'MRP', p.initial = c(p.true,1))
  lines(10^x,mf2$fit$functions$mass.function(x),col='red')

  # legend
  y = 1.5e-2
  x = 5.5
  points(10^x,y,pch=20)
  segments(10^(x-sigma),y,10^(x+sigma),y,lwd=1)
  text(10^x,y*1.6,expression('Assumed observing error'~sigma),cex=0.8)
  legend(2e4,6e-3,c('Input Veff (arbitrary vertical units)',
                    'Input Schechter function',
                    'Binned mock observations with errors',
                    'Fitted Schechter function with 68%-uncertainties',
                    'Fitted MRP function',
                    'MF from randomly sampled posterior mass PDFs',
                    'MF from full posterior mass PDFs',
                    'Note: Horizontal bars are bin widths.',
                    'Vertical bars are Poission errors (correlated for purple points).'),
         bty='n',lwd=1.5,col=c('orange','black','grey','blue','red','black','purple',NA,NA),
         pch=c(NA,NA,20,NA,NA,20,20,NA,NA),
         lty=c(1,2,1,1,1,1,1,NA,NA),cex=0.8)

  # plot covariances
  mfplotcovariance(mf,p.true,title='Schechter parameter covariances')
  mfplotcovariance(mf2,c(p.true,1),line.color='red',point.color='#ffbbbb',title='MRP parameter covariances')
}
