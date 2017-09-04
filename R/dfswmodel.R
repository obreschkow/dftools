#' Make step-wise linear distribution function
#'
#' This function converts a continuous one-dimensional distribution function (DF), such as a galaxy mass function, into a step-wise DF on an arbitrary number of bins.
#' 
#' @importFrom pracma grad
#' @importFrom stats approx splinefun
#'
#' @param gdf is the one-dimensional generative DF to be converted into a step-wise DF.
#' @param xedges is a vector of observables x defining the edges of the bins for the step-wise DF.
#' @param method	specifies the interpolation method to be used. Choices are "constant", "linear" and "spline". In the latter two cases the generative distribution function is linearly interpolated beyond the range of \code{xedges}.
#' @param extrapolate is a logical flag. If \code{TRUE}, the step-wise DF is linearly extraplated outside the range of \code{xedges}. This extrapolation is useful when fitting the stepw-wise DF to data.
#' 
#' @return \code{dfswmodel} returns a list with two components \code{gdf} and {p.initial}. The first component is a function of log-masses \code{x} and a parameter vector \code{p}. The second component is the initial parameter-vector producing a step-wise mass function that closely matches the Schechter function specified by \code{p}.
#'
#' @examples
#' # For an overall illustration run
#' dfexample(3)
#' 
#' # The follwoing is a step-by-step example.
#' # First, define and visulize a smooth Schechter function (a special type of galaxy mass function)
#' p.true = c(-2,11,-1.4) # Schechter function parameters
#' gdf = function(x) dfmodel(x, p.true, type='Schechter') # Schechter function
#' x = seq(5.8,12,0.01) # plotting range
#' plot(10^x, gdf(x), type='l', log='xy', ylim=c(1e-5,2))
#' 
#' # Define the edges of non-equal bins
#' xbin = c(6,7,7.5,seq(8,9,0.2),seq(9.5,11,0.5),11.8) # edges of bins
#' 
#' # Make and plot step-wise constant function
#' swconstant = dfswmodel(gdf, xbin, 'constant') # step-wise constant function
#' lines(10^x, swconstant$gdf(x, swconstant$p.initial), col='blue', lwd=3)
#' 
#' # Make and plot step-wise linear function
#' swlinear = dfswmodel(gdf, xbin, 'linear') # step-wise linear function
#' lines(10^x, swlinear$gdf(x, swlinear$p.initial), col='orange', lwd=2)
#' 
#' # Make and plot step-wise spline function
#' swspline = dfswmodel(gdf, xbin, 'spline') # step-wise linear function
#' lines(10^x, swspline$gdf(x, swspline$p.initial), col='red')
#' 
#' # Plot mid-points of bins
#' points(10^swconstant$xmid, 10^swconstant$p.initial, pch=20)
#'
#' @seealso \code{\link{dffit}}
#'
#' @author Danail Obreschkow
#'
#' @export

dfswmodel <- function(gdf = function(x) dfmodel(x, dfmodel(output = 'initial'), 'density'),
                      xedges = seq(6,12,0.5),
                      method = 'linear',
                      extrapolate = TRUE) {
  
  # input handling
  nedges = length(xedges)
  if (nedges<2) stop('xedges needs at least two elements')
  nbins = nedges-1
  if (any(xedges[2:nedges]<=xedges[1:nbins])) stop('xedges must be a monotonically increasing vector.')
  xmid = (xedges[1:nbins]+xedges[2:nedges])/2
  ymid = log10(gdf(xmid))
  if (any(!is.finite(ymid))) stop('cannot evaluate Schechter function at all values of xpoint.')
  
  # interpolation function
  if (method=='constant') {
    gdf = function(x,p) {
      y = approx(xedges,c(p,0),xout=x,method='constant',f=0,rule=1)$y
      if (!extrapolate) {
        y[x<=xedges[1]] = NA
        y[x>=xedges[nedges]] = NA
      } else {
        slope.left = (p[2]-p[1])/(xmid[2]-xmid[1])
        slope.right = (p[nbins]-p[nbins-1])/(xmid[nbins]-xmid[nbins-1])
        list = x<=xedges[1]
        y[list] = (x[list]-xmid[1])*slope.left+p[1]
        list = x>=xedges[nedges]
        y[list] = (x[list]-xmid[nbins])*slope.right+p[nbins]
      }
      return(10^y)
    }
  } else if (method=='linear') {
    gdf = function(x,p) {
      y = approx(xmid,p,xout=x,method='linear',rule=1)$y
      slope.left = (p[2]-p[1])/(xmid[2]-xmid[1])
      slope.right = (p[nbins]-p[nbins-1])/(xmid[nbins]-xmid[nbins-1])
      list = x<=xmid[1]
      y[list] = (x[list]-xmid[1])*slope.left+p[1]
      list = x>=xmid[nbins]
      y[list] = (x[list]-xmid[nbins])*slope.right+p[nbins]
      if (!extrapolate) {
        y[x<=xedges[1]] = NA
        y[x>=xedges[nedges]] = NA
      }
      return(10^y)
    }
  } else if (method=='spline') {
    gdf = function(x,p) {
      f = splinefun(xmid,p)
      y = f(x)
      slope.left = pracma::grad(f,xmid[1])
      slope.right = pracma::grad(f,xmid[nbins])
      list = x<xmid[1]
      y[list] = (x[list]-xmid[1])*slope.left+p[1]
      list = x>xmid[nbins]
      y[list] = (x[list]-xmid[nbins])*slope.right+p[nbins]
      if (!extrapolate) {
        y[x<=xedges[1]] = NA
        y[x>=xedges[nedges]] = NA
      }
      return(10^y)
    }
  } else {
    stop('unknown interpolation method')
  }
  return(list(gdf = gdf, p.initial = ymid, xmid = xmid, xmin = xedges[1], xmax = xedges[nedges]))
}