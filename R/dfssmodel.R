#' Make step-wise linearmass function
#'
#' This function converts a Schechter function into a step-wise linear mass function of an arbitrary number of bins, to fit a quasi non-parametric mass functions.
#' 
#' @importFrom pracma grad
#'
#' @param xpoints Vector of log-masses used as interpolation points.
#' @param p Schechter function parameters used to determine the initial values of the step-wise paraemters.
#' 
#' @return \code{dfssmodel} returns a list with two components \code{gdf} and {p.initial}. The first component is a function of log-masses \code{x} and a parameter vector \code{p}. The second component is the initial parameter-vector producing a step-wise mass function that closely matches the Schechter function specified by \code{p}.
#'
#' @examples
#' # Basic example of a stepwise linear MF
#' xpoints = seq(6.5,12,0.5)
#' model = dfssmodel(xpoints)
#' plot(10^xpoints,10^model$p.initial,log='xy',ylim=c(1e-5,1),pch=20)
#' lines(10^seq(4,12,0.01),model$gdf(seq(4,12,0.01),model$p.initial))
#' 
#' # Example of fitting a stepwise linear MF
#' dfexample4()
#'
#' @seealso \code{\link{dffit}}
#'
#' @author Danail Obreschkow
#'
#' @export

dfssmodel <- function(xpoints,
                      method = 'linear',
                      p.Schechter = dfmodel(output = 'initial')) {
  
  # input handling
  npoints = length(xpoints)
  if (npoints<2) stop('xpoints needs at least two elements')
  if (any(xpoints[2:npoints]<=xpoints[1:npoints-1])) stop('xpoints must be a monotonically increasing vector.')
  ypoints = log10(dfmodel(xpoints, p = p.Schechter))
  if (any(!is.finite(ypoints))) stop('cannot evaluate Schechter function at all values of xpoint.')
  
  # interpolation function
  gdf = function(x,p) {
    if (method=='constant') {
      y = approx(xpoints,p,xout=x,method='constant',f=0.5,rule=1)$y
      interpolate = FALSE
    } else if (method=='linear') {
      y = approx(xpoints,p,xout=x,method='linear',rule=1)$y
      slope.left = (p[2]-p[1])/(xpoints[2]-xpoints[1])
      slope.right = (p[npoints]-p[npoints-1])/(xpoints[npoints]-xpoints[npoints-1])
      interpolate = TRUE
    } else if (method=='spline') {
      f = splinefun(xpoints,p)
      y = f(x)
      slope.left = pracma::grad(f,xpoints[1])
      slope.right = pracma::grad(f,xpoints[npoints])
      interpolate = TRUE
    } else {
      stop('unknown interpolation method')
    }
    if (interpolate) {
      list = x<xpoints[1]
      y[list] = (x[list]-xpoints[1])*slope.left+p[1]
      list = x>xpoints[npoints]
      y[list] = (x[list]-xpoints[npoints])*slope.right+p[npoints]
    } else {
      y[x<=xpoints[1]] = NA
      y[x>=xpoints[npoints]] = NA
    }
    return(10^y)
  }
  
  return(list(gdf = gdf, p.initial = ypoints))
}