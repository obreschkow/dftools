#' Example of fitting a galaxy mass function
#'
#' This function generates n mock galaxies with masses and observing errors, drawn from a custom Schechter function and predefined effective volume function. It then uses \code{\link{dffit}} to recover the inpu Schechter function from the noisy data.
#' 
#' @importFrom graphics legend text segments
#' 
#' @param seed A positive integer used as seed for the random number generator.
#' @param n Number of gtalaxies
#' @param sigma Observational uncertainty, defined as standard-deviation in log10(Mass)
#' @param p.true Input parameters for the Schechter function
#' @param include.mrp If \code{TRUE}, an MRP function will also be fitted
#'
#' @examples
#' # Run basic example
#' dfexample1()
#' 
#' # Run basic example and also fit an MRP function
#' dfexample1(include.mrp = TRUE)
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfexample1 <- function(n = 1e3, seed = 1, sigma = 0.5, p.true = c(-2,10,-1.3), include.mrp = FALSE) {

  # make rescaled Veff(x)
  cat('Generate mock data with observing errors:\n')
  dat = dfmockdata(n = n, seed = seed, sigma = sigma, p = p.true)
  cat(sprintf('=> %d galaxies\n',dat$n))

  # fit
  cat('Fit a Schechter function to the mock data:\n')
  survey = dffit(dat$x, dat$veff, dat$x.err, write.fit = T,
                 xmin = dat$xmin, xmax = dat$xmax, p.initial = dat$p)
  p.fit = survey$fit$parameters$p.optimal
  
  # MRP
  if (include.mrp) {
    cat('Fit a MRP function (4 parameter model) to the mock data:\n')
    survey2 = dffit(dat$x, dat$veff, dat$x.err, write.fit = T, xmin = dat$xmin, xmax = dat$xmax, gdf = 'MRP', p.initial = c(p.true,1))
  }
  
  # plot covariances
  dfplotcov(survey,expectation2 =p.true,title='Schechter parameter covariances',model.col='blue',model2.cross.lty = 2)
  if (include.mrp) dfplotcov(survey2,expectation2=c(p.true,1),model.col='red',title='MRP parameter covariances',model2.cross.lty = 2)

  # plot main plot
  nbins = max(4,round(sqrt(n)/2))
  mfplot(survey, xlim=c(1e6,2e11), ylim=c(1e-4,2), p = p.true, nbins = nbins, bin.xmin=7, bin.xmax=11)
  x = seq(2,10.65,0.01)
  if (include.mrp) {lines(10^x,survey2$fit$gdf(x),col='red')}

  # legend
  list = rep(TRUE,5)
  list[3] = include.mrp
  x = 7; y = 1.5e-2
  points(10^x,y,pch=20)
  segments(10^(x-sigma),y,10^(x+sigma),y,lwd=1)
  text(10^x,y*1.6,expression('Observing error'~sigma),cex=0.8)
  legend(2e6,6e-3,c('Input Schechter function',
                    'Fitted Schechter function with 68%-uncertainties',
                    'Fitted MRP function',
                    'Binned input data',
                    'Binned posterior data')[list],
         bty='n',lwd=c(1,1.5,1,1,1)[list],
         col=c('black','blue','red','purple','black')[list],
         pch=c(NA,NA,NA,20,20)[list],
         lty=c(2,1,1,1,1)[list],cex=0.8)

  # console output
  cat(sprintf('Total computation time       = %.2fs\n',survey$fit$status$walltime.total))
  cat(sprintf('Computation time for fitting = %.2fs\n',survey$fit$status$walltime.fitting))
  
  invisible(survey)
}
