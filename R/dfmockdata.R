#' Example data
#'
#' This function produces an example of log-masses \code{x} with Gaussian uncertainties \code{x.err} and a selection \code{selection}. This selection is a list of an array with galaxy-by-galaxy effective volumes and a function to complete the effective volume beyond the observed mass range. The data is taken from a publication by Westmeier et al. (2017) and represents HI-Masses in the Sculptor structure.
#' 
#' @importFrom pracma erf
#' 
#' @examples
#' data <- dfdata()
#' print(data$x)
#' print(data$x.err)
#' print(data$selection[[1]])
#'
#' # These data can then be used to fit a MF, e.g.
#' df = dffit(data$x, data$selection, data$x.err)
#'
#' @seealso \code{\link{dffit}}
#'
#' @author Danail Obreschkow
#'
#' @export

dfmockdata <- function(n = 1000,
                       seed = 0,
                       f = function(x,r) {mass = 10^x; rmax = sqrt(mass); f = pracma::erf((rmax-r)/(rmax+1))*0.5+0.5},
                       dVdr.scale = function(r) {r^2},
                       gdf = function(x,p) dfmodel(x, p, type = 'Schechter'),
                       p = c(-2,10,-1.3),
                       sigma = 0.0,
                       rmin = 0, rmax = 20,
                       xmin = 0, xmax = 12#, dx = 0.01 # must be a large enough range to contain all possible data
                       ) {
  
  set.seed(seed)
  
  # rescale observing volume to match the requested number of galaxies
  veff.scale.elemental = function(x) {
    fct = function(r) f(x,r)*dVdr.scale(r)
    return(integrate(fct,rmin,rmax)$value)
  }
  veff.scale = function(x) sapply(x,veff.scale.elemental)
  scd.scale = function(x) veff.scale(x)*gdf(x,p)
  n.scale = integrate(scd.scale,xmin,xmax)$value
  rescaling.factor = n/n.scale
  dVdr = function(r) dVdr.scale(r)*rescaling.factor
  
  # make veff.function
  veff.function.elemental = function(x) {
    fct = function(r) f(x,r)*dVdr(r)
    return(integrate(fct,rmin,rmax)$value)
  }
  veff.function = Vectorize(veff.function.elemental)
  
  # make source count density function
  scd = function(x) veff.function(x)*gdf(x,p)
  
  # sample data
  # cum = function(x) integrate(scd,xmin,x)$value
  # rand = runif(n,0,n)
  # x = array(NA,c(n,1))
  # for (i in seq(n)) {
  #   fct = function(x) cum(x)-rand[i]
  #   x[i] = uniroot(fct,c(xmin,xmax))$root
  # }
    
  # sample data
  dx = 0.001
  xgrid = seq(xmin,xmax,dx)
  m = length(xgrid)
  x = array(NA,c(n,1))
  cum = cumsum(scd(xgrid))
  cum = cum/cum[m]
  rand = runif(n)
  for (i in seq(n)) {
    index = which.min(abs(cum-rand[i]))
    if (index<=2) stop('There are sampled values below xmin. Choose a lower value for xmin.')
    if (index>=(m-1)) stop('There are sampled values above xmax. Choose a larger value for xmax.')
    x[i] = xgrid[index]
  }
  
  # add observing errors
  if (!is.null(sigma)) {
    x = x+rnorm(n)*sigma
    x.err = rep(sigma,n)
  } else {
    x.err = NULL
  }
  
  # add distances
  r = NULL
  
  # make effective volumes for each observation, that an observer would assign, not knowning the observational error in x
  veff.values = veff.function(x)
  
  # return
  return(list(x = x, x.err = x.err, r = r,
              f = f, dVdr = dVdr, veff.function = veff.function, veff.values = veff.values,
              rescaling.factor = rescaling.factor))
  
}