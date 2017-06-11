#' Evaluate posterior masses
#'
#' Given a fitted df, this function computes posterior masses.
#'
#' @param df List produced by \code{\link{fit.df}}
#'
#' @return Returns a structured list, which, in addition to the input argument, also contains the sublist \code{posteriors}. Note that the posteriors can be visualized using \code{\link{dfplot}} with the argument \code{bin.use.posteriors = TRUE}.
#'
#' @seealso See examples in \code{\link{fit.df}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfposteriors <- function(df) {
  
  # Input handling
  xr = df$input$options$x.mesh
  dx = xr[2]-xr[1]
  n = length(xr)
  n.data = dim(df$input$data$x)[1]
  n.dim = dim(df$input$data$x)[2]
  
  # Make prior
  prior = df$fit$functions$source.count(xr)
  
  # Make posterior
  rho.unbiased = array(0,length(xr))
  m0 = m1 = md = array(NA,n.data)
  for (i in seq(n.data)) {
    
    rho.observed = exp(-(df$input$data$x[i]-xr)^2/2/df$input$data$x.err[i]^2)
    rho.corrected = rho.observed*prior
    
    # total source counts
    s = sum(rho.corrected)
    rho.unbiased = rho.unbiased+rho.corrected/(s*dx)
    
    # mean and standard deviation
    m0[i] = sum(xr*rho.corrected)/s
    m1[i] = sqrt(sum((xr-m0[i])^2*rho.corrected)/s)
    
    # find maximum point (mode) using parabolic fit
    i0 = which.max(rho.corrected)
    y0 = rho.corrected[i0]
    yp = rho.corrected[min(n,i0+1)]
    ym = rho.corrected[max(1,i0-1)]
    xpeak = -(yp-ym)/(2*(yp+ym-2*y0))
    md[i] = xr[i0]+xpeak*dx
  }
  
  df$posterior$source.count.density = rho.unbiased
  df$posterior$x.mean = m0
  df$posterior$x.stdev = m1
  df$posterior$x.mode = md
  df$posterior$x.mode.correction = md-df$input$mass
  df$posterior$x.random = m0+m1*rnorm(n.data)
  
  return(df)
}