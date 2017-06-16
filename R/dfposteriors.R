#' Evaluate posterior masses
#'
#' Given a fitted df, this function computes posterior masses.
#'
#' @param df List produced by \code{\link{dffit}}
#'
#' @return Returns a structured list, which, in addition to the input argument, also contains the sublist \code{posteriors}. Note that the posteriors can be visualized using \code{\link{dfplot}} with the argument \code{bin.use.posteriors = TRUE}.
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfposteriors <- function(bundle) {
  
  if (is.null(bundle$data$x.err)) stop('Posterior data PDFs can only be produced if the data is uncertain, i.e. if x.err is given.')
  
  # Input handling
  x = bundle$data$x
  x.err = bundle$data$x.err
  x.mesh = bundle$grid$x
  x.mesh.dv = bundle$grid$dvolume
  n.data = dim(x)[1]
  n.dim = dim(x)[2]
  
  # Make prior
  prior = bundle$grid$scd
  prior[!is.finite(prior)] = 0
  prior = pmax(0,prior)
  
  # Make inverse covariances
  invC = array(NA,c(n.data,n.dim,n.dim))
  if (length(dim(x.err))==2) {
    if (!(dim(x.err)[1]==n.data & dim(x.err)[2]==n.dim)) {
      stop('Unknown format for x.err in .corefit.')
    }
    if (n.dim==1) {
      for (i in seq(n.data)) {
        invC[i,,] = 1/x.err[i,]^2
      }
    } else {
      for (i in seq(n.data)) {
        invC[i,,] = diag(1/x.err[i,]^2)
      }
    }
  } else if (length(dim(x.err))==3) {
    if (n.dim==1) stop('Unknown format for x.err in .corefit.')
    if (!(dim(x.err)[1]==n.data & dim(x.err)[2]==n.dim & dim(x.err)[3]==n.dim)) {
      stop('Unknown format for x.err in .corefit.')
    }
    for (i in seq(n.data)) {
      invC[i,,] = solve(x.err[i,,])
    }
  } else {
    stop('Unknown format for x.err in .corefit.')
  }
    
  # produce posteriors
  m0 = m1 = md = array(NA,c(n.data,n.dim))
  rho.unbiased = array(0,dim(x.mesh))
  for (i in seq(n.data)) {
    
    # make prior PDF for data point i
    d = x[i,]-t(x.mesh)
    rho.observed = exp(-colSums(d*(invC[i,,]%*%d))/2)
    
    # make posterior PDF for data point i
    rho.corrected = rho.observed*prior
    s = sum(rho.corrected)
    rho.unbiased = rho.unbiased+rho.corrected/(s*x.mesh.dv)
    
    # mean, standard deviation and mode
    for (j in seq(n.dim)) {
      m0[i,j] = sum(x.mesh[,j]*rho.corrected)/s
      m1[i,j] = sqrt(sum((x.mesh[,j]-m0[i,j])^2*rho.corrected)/s)
    }
    md[i,] = x.mesh[which.max(rho.corrected),]
  }
  
  bundle$posterior = list(x.mean = m0, x.stdev = m1, x.mode = md, x.random = m0+m1*array(rnorm(n.data*n.dim),c(n.data,n.dim)),
                          scd = rho.unbiased)
  
  invisible(bundle)
}