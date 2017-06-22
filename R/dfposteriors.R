#' Evaluate posterior masses
#'
#' This function is called after fitting a generative distribution function usnig \code{\link{dffit}}. It uses the output list of \code{dffit} to deterine the posterior probability distributions of the uncertain data, given the fitted generative distribution function (such as a galaxy mass function). \code{dfposteriors} works for models of any dimension.
#'
#' @param survey List produced by \code{\link{dffit}}
#'
#' @return Returns a structured list containing the input list \code{survey} with an additional sublist \code{posteriors}. This sublist contains several measures of the posterior data: means, standard deviations and modes, as well as a random value drawn from the actual posterior distribution. This random value can be used to plot unbiased distribution functions, such as mass functions. The output also contains the array \code{$grid$scd.posterior} giving the posterior source count density of the input data evaluated on the grid \code{$grid$x}.
#' 
#' @examples
#' # make a mock survey drawn from a Schechter function, where masses have 0.3 dex observing errors
#' dat = dfmockdata(sigma = 0.3)
#' 
#' # fit a Schechter function to the noisy data
#' survey = dffit(dat$x,dat$veff,dat$x.err)
#' 
#' # make posterior masses
#' survey = dfposteriors(survey)
#' 
#' # plot mass function (MF) of raw data and input model as dashed line
#' mfplot(survey,bin.type=1,col.data='grey',xlim=c(1e6,1e12))
#' lines(10^survey$grid$x, survey$model$gdf(survey$grid$x,c(-2,10,-1.3)),lty=2)
#' 
#' # add MF using random values from posterior data PDFs:
#' mfplot(survey,bin.type=2,col.data='orange',add=TRUE)
#' 
#' # add MF using the full posterior PDFs of all data points:
#' mfplot(survey,bin.type=3,col.data='black',add=TRUE)
#'
#' @seealso Note that the posteriors can be visualized using \code{\link{dfplot}} with the argument \code{bin.use.posteriors = TRUE}. See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfposteriors <- function(survey) {
  
  if (is.null(survey$data$x.err)) stop('Posterior data PDFs can only be produced if the data is uncertain, i.e. if x.err is given.')
  
  # Input handling
  x = survey$data$x
  x.err = survey$data$x.err
  x.mesh = survey$grid$x
  x.mesh.dv = survey$grid$dvolume
  n.data = dim(x)[1]
  n.dim = dim(x)[2]
  
  # Make prior
  prior = survey$grid$scd
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
  rho.unbiased = rho.unbiased.sqr = array(0,dim(x.mesh))
  for (i in seq(n.data)) {
    
    # make prior PDF for data point i
    d = x[i,]-t(x.mesh)
    rho.observed = exp(-colSums(d*(invC[i,,]%*%d))/2)
    
    # make posterior PDF for data point i
    rho.corrected = rho.observed*prior
    s = sum(rho.corrected)
    rho.unbiased = rho.unbiased+rho.corrected/(s*x.mesh.dv)
    rho.unbiased.sqr = rho.unbiased.sqr+(rho.corrected/(s*x.mesh.dv))^2
    
    # mean, standard deviation and mode
    for (j in seq(n.dim)) {
      m0[i,j] = sum(x.mesh[,j]*rho.corrected)/s
      m1[i,j] = sqrt(sum((x.mesh[,j]-m0[i,j])^2*rho.corrected)/s)
    }
    md[i,] = x.mesh[which.max(rho.corrected),]
  }
  
  survey$posterior = list(x.mean = m0, x.stdev = m1, x.mode = md, x.random = m0+m1*array(rnorm(n.data*n.dim),c(n.data,n.dim)))
  survey$grid$scd.posterior = rho.unbiased
  survey$grid$effective.counts = rho.unbiased^2/rho.unbiased.sqr # this equation gives the effective number of sources per bin
  survey$grid$effective.counts[!is.finite(survey$grid$effective.counts)] = 0
  
  invisible(survey)
}