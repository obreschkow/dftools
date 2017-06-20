#' Example of using dftools
#'
#' This function generates n mock galaxies with masses and observing errors, drawn from a custom Schechter function and predefined effective volume function. It then uses \code{\link{dffit}} to recover the inpu Schechter function from the noisy data.
#'
#' @param n Number of gtalaxies
#' @param seed Random number seed
#' @param sigma Observational uncertainty, defined as standard-deviation in log10(Mass)
#' @param p.true Input parameters for the Schechter function
#' @param include.mrp If \code{TRUE}, an MRP function will also be fitted
#'
#' @examples
#' dfexample1()
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfexample1 <- function(n = 1e4, seed = 1, sigma = 0.5, p.true = c(-2,10,-1.3), include.mrp = FALSE) {

  # user parameters
  set.seed(seed)
  veff.scale = function(x) {1e-9*10^(0.8*x)}

  # make rescaled Veff(x)
  dat = dfmockdata(n = n, seed = seed, sigma = sigma, p = p.true)

  # fit
  cat('Fit a Schechter function to the observed data (grey points)\n')
  survey = dffit(dat$x, dat$veff, dat$x.err, write.fit = T,
                 xmin = dat$xmin, xmax = dat$xmax, p.initial = dat$p.true)
  p.fit = survey$fit$parameters$p.optimal

  # make posterior masses
  cat('Determine posterior masses\n')
  survey = dfposteriors(survey)
  
  # MRP
  if (include.mrp) {
    cat('Fit a MRP function (4 parameter model) to the observed data\n')
    survey2 = dffit(dat$x, dat$veff, dat$x.err, write.fit = T, xmin = dat$xmin, xmax = dat$xmax, phi = 'MRP', p.initial = c(p.true,1))
  }
  
  # plot covariances
  dfplotcov(survey,p.true,title='Schechter parameter covariances')
  if (include.mrp) dfplotcov(survey2,c(p.true,1),line.color='red',point.color='#ffbbbb',title='MRP parameter covariances')

  # plot main plot
  x = seq(2,12,0.01)
  nbins = max(4,round(sqrt(n)/2))
  mfplot(survey, xlim=c(1e6,2e11), ylim=c(1e-4,2), nbins = nbins, bin.xmin=7, bin.xmax=11, col.data = 'grey')
  mfplot(survey, nbins = nbins, bin.xmin=6, bin.xmax=11, bin.type = 3, add = TRUE, show.uncertainties = FALSE, col.data='black')
  lines(10^x,pmax(1e-10,dfmodel(x,p.true)),col='black',lty=2)
  lines(10^x,dat$veff(x)/max(dat$veff.values),col='orange')
  if (include.mrp) {lines(10^x,survey2$fit$functions$phi.fit(x),col='red')}

  # legend
  y = 1.5e-2
  x = 7
  points(10^x,y,pch=20)
  segments(10^(x-sigma),y,10^(x+sigma),y,lwd=1)
  text(10^x,y*1.6,expression('Assumed observing error'~sigma),cex=0.8)
  legend(2e6,3e-3,c('Input Veff (arbitrary vertical units)',
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
  cat(sprintf('Total computation time       = %.2fs\n',survey$fit$status$walltime.total))
  cat(sprintf('Computation time for fitting = %.2fs\n',survey$fit$status$walltime.fitting))
  
  invisible(survey)
}
