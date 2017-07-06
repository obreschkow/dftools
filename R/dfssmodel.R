#' Make step-wise mass function
#'
#' This function converts a Schechter function into a step-wise mass function of an arbitrary number of bins, to fit a quasi non-parametric mass functions.
#'
#' @param xpoints Vector of log-masses used as interpolation points.
#' @param p Schechter function parameters used to determine the initial values of the step-wise paraemters.
#' @param method specifies the interpolation method to be used. Choices are "linear" or "constant".
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
                      p.Schechter = dfmodel(output = 'initial'),
                      method = 'linear') {
  
  # input handling
  npoints = length(xpoints)
  if (npoints<2) stop('xpoints needs at least two elements')
  if (any(xpoints[2:npoints]<=xpoints[1:npoints-1])) stop('xpoints must be a monotonically increasing vector.')
  ypoints = log10(dfmodel(xpoints, p = p.Schechter))
  if (any(!is.finite(ypoints))) stop('cannot evaluate Schechter function at all values of xpoint.')
  
  # interpolation function
  gdf = function(x,p) {
    y = approx(xpoints,p,x,rule=2,method=method,f=0.5)$y
    if (method == 'linear') {
      list = x<xpoints[1]
      slope = (p[2]-p[1])/(xpoints[2]-xpoints[1])
      y[list] = (x[list]-xpoints[1])*slope+p[1]
      list = x>xpoints[npoints]
      slope = (p[npoints]-p[npoints-1])/(xpoints[npoints]-xpoints[npoints-1])
      y[list] = (x[list]-xpoints[npoints])*slope+p[npoints]
    }
    return(10^y)
  }
  
  return(list(gdf = gdf, p.initial = ypoints))
}