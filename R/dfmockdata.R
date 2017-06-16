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

dfmockdata <- function(n = 100,
                       seed = 0,
                       f = function(x,r) {mass = 10^x; rmax = sqrt(mass); f = pracma::erf((rmax-r)/(rmax+1))*0.5+0.5},
                       dVdr.scale = function(r) {r^2},
                       gdf = function(x,p) dfmodel(x, p, type = 'Schechter'),
                       p = c(-2,10,-1.3),
                       sigma = 0.1,
                       rmin = 0, rmax = 20,
                       xmin = 0, xmax = 12, dx = 0.01 # must be a large enough range to contain all possible data
                       ) {
  
  set.seed(seed)
  
  # make x-grid
  xgrid = seq(xmin,xmax,dx)
  m = length(xgrid)
  
  # rescale observing volume to match the requested number of galaxies
  n.expected = 0
  for (i in seq(m)) {
    intfct = function(r) f(xgrid[i],r)*dVdr.scale(r)
    veff.scale = integrate(intfct,rmin,rmax)$value
    n.expected = n.expected+veff.scale*gdf(xgrid[i],p)
  }
  rescaling.factor = n/n.expected
  dVdr = function(r) dVdr.scale(r)*rescaling.factor
  
  # make veff.function
  veff.function.elemental = function(x) {
    intfct = function(r) f(x,r)*dVdr.scale(r)
    return(integrate(intfct,rmin,rmax)$value)
  }
  veff.function = function(x) {
    sapply(x,veff.function.elemental,USE.NAMES = FALSE)
  }
  
  # make source count density function
  scd = function(x) {
    veff.function(x)*gdf(x,p)
  }
    
  # sample data
  x = array(NA,c(n,1))
  cum = cumsum(scd(xgrid))
  cum = cum/cum[m]
  rand = runif(n)
  for (i in seq(n)) {
    index = which.min(abs(cum-rand[i]))
    if (index==1) stop('There are sampled values below xmin. Choose a lower value for xmin.')
    if (index==m) stop('There are sampled values above xmax. Choose a larger value for xmax.')
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