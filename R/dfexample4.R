#' Example of fitting a step-wise mass function
#'
#' This function generates n mock galaxies with masses and observing errors, drawn from a custom Schechter function and predefined effective volume function. It then uses \code{\link{dffit}} to recover the inpu Schechter function from the noisy data.
#' 
#' @importFrom graphics legend text segments
#' 
#' @param seed A positive integer used as seed for the random number generator.
#' @param n Number of gtalaxies
#' @param sigma Observational uncertainty, defined as standard-deviation in log10(Mass)
#'
#' @examples
#' # Run basic example
#' dfexample4()
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfexample4 <- function(n = 1e4, seed = 1, sigma = 0) {
  
  # make analytic non-Schechter gdf
  gdf.analytic = function(x,p) {
    component1 = dfmodel(x,p[1:3])
    component2 = p[4]*exp(-(x-9.3)^2*2)
    return(component1+component2)
  }
  
  # draw random data
  p.true = c(-2,10,-1.3,1)
  dat = dfmockdata(n = n, seed = seed, sigma = sigma, gdf = gdf.analytic, p = p.true,
                   xmax = 13)
  
  # make step-wise linear fitting function
  xpoints = seq(7,11,1)
  ss = dfssmodel(xpoints,method = 'linear')
  
  # fitt step-wise linear function
  survey = dffit(dat$x, dat$veff, dat$x.err, gdf = ss$gdf, p.initial = ss$p.initial,
                 xmin = dat$xmin, xmax = dat$xmax)
  
  # plot
  mfplot(survey,xlim=c(1e6,1e11))
  x = seq(4,12,0.01)
  lines(10^x,gdf.analytic(x,p.true))
  
  invisible(survey)
}
survey = dfexample4()
