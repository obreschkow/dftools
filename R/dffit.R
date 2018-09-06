#' Fit a generative distribution function, such as a galaxy mass function
#'
#' This function finds the most likely P-dimensional model parameters of a D-dimensional distribution function (DF) generating an observed set of N objects with D-dimensional observables x, accounting for measurement uncertainties and a user-defined selection function. For instance, if the objects are galaxies, \code{dffit} can fit a mass function (D=1), a mass-size distribution (D=2) or the mass-spin-morphology distribution (D=3). A full description of the algorithm can be found in Obreschkow et al. (2017).
#'
#' @importFrom akima interp
#' @importFrom stats optim rpois quantile approxfun cov var
#'
#' @param x is an N-by-D matrix (or N-element vector if D=1) containing the observed D-dimensional properties of N objects. For instance, x can be N-element vector containing the logarithmig masses of N galaxies, to which a mass function should be fitted.
#' 
#' @param selection Specifies the effective volume \code{V(xval)} in which an object of (D-dimensional) true property (e.g. true log-mass) \code{xval} can be observed. This volume can account for complex selection criteria, as long as they can be expressed as a function of the \emph{true} properties, i.e. before measurement errors occur. This volume can be specified in five ways:\cr\cr
#' 
#' (1) If \code{selection} can be a single positive number. This number will be interpreted as a constant volume, \code{V(xval)=selection}, in which all objects are fully observable. \code{V(xval)=0} is assumed outside the "observed domain". This domain is defined as \code{min(x)<=xval<=max(x)} for a scalar observable (D=1), or as \code{min(x[,j])<=xval[j]<=max(x[,j])} for all j=1,...,D if D>1. This mode can be used for volume-complete surveys or for simulated galaxies in a box.\cr\cr
#' 
#' (2) \code{selection} can be an N-element vector. The elements will be interpreted as the volumes of each galaxy. \code{V(xval)} is interpolated (linearly in \code{1/V}) for other values \code{xval}. \code{V(xval)=0} is assumed outside the observed domain, except if D=1 (as in the case of fitting mass functions), where \code{V(xval)=0} is assumed only if \code{xval<min(x)}, whereas for \code{xval>max(x)}, \code{V(xval)} is set equal to the maximum effective volume, i.e. the maximum of \code{selection}. \cr\cr
#' 
#' (3) \code{selection} can be a function of D variables, which directly specifies the effective volume for any \code{xval}, i.e. \code{Veff(xval)=selection(xval)}.\cr\cr
#' 
#' (4) \code{selection} can be a list (\code{selection = list(veff.values, veff.userfct)}) of an \code{N}-element vector \code{veff.values} and a \code{D}-dimensional function \code{veff.userfct}. In this case, the effective volume is computed using a hybrid scheme of modes (2) and (3): \code{V(xval)} will be interpolated from the N values of \code{veff.values} inside the observed domain, but set equal to \code{veff.userfct} outside this domain.\cr\cr
#' 
#' (5) \code{selection} can be a list of two functions and one 2-element vector: \code{selection = list(f, dVdr, rmin, rmax)}, where \code{f = function(xval,r)} is the isotropic selection function and \code{dVdr = function(r)} is the derivative of the total survey volume as a function of comoving distance \code{r}. The scalars \code{rmin} and \code{rmax} (can be \code{0} and \code{Inf}) are the minimum and maximum comoving distance limits of the survey. Outside these limits \code{V(xval)=0} will be assumed.\cr\cr
#' 
#' @param x.err specifies the observational uncertainties of the N data points. These uncertainties can be specified in four ways:\cr\cr
#' 
#' (1) If \code{x.err = NULL}, the measurements \code{x} will be considered exact.\cr\cr
#' 
#' (2) If \code{x.err} is a \code{N-by-D} matrix, the scalars \code{x.err[i,j]} are interpreted as the standard deviations of Gaussian uncertainties on \code{x[i,j]}.\cr\cr
#' 
#' (3) If \code{x.err} is a \code{N-by-D-by-D} array, the \code{D-by-D} matrices \code{x.err[i,,]} are interpreted as the covariance matrices of the D observed values \code{x[i,]}.\cr\cr
#' 
#' (4) \code{x.err} can also be a function of a D-vector \code{x} and an integer \code{i}, such that \code{x.err(x,i)} is the prior probability distribution function of the data point \code{i} to have the true value \code{x}. This function must be vectorized in the first argument, such that calling \code{x.err(x,i)} with \code{x} being a N-by-D matrix returns an N-element vector.\cr\cr
#' 
#' @param r Optional N-element vector specifying the comoving distances of the N objects (e.g. galaxies). This vector is only needed if \code{correct.lss.bias = TRUE}.
#' @param gdf Either a string or a function specifying the DF to be fitted. A string is interpreted as the name of a predefined mass function (i.e. functions of one obervable, \code{D=1}). Available options are \code{'Schechter'} for Schechter function (3 parameters), \code{'PL'} for a power law (2 parameters), or \code{'MRP'} for an MRP function (4 parameters). Alternatively, \code{gdf = function(xval,p)} can be any function of the \code{P} observable(s) \code{xval} and a list of parameters \code{p}. The function \code{gdf(xval,p)} must be fully vectorized in \code{xval}, i.e. it must output a vector of \code{N} elements if \code{xval} is an \code{N-by-P} array (such as the input argument \code{x}). Note that if \code{gdf} is given as a function, the argument \code{p.initial} is mandatory.
#' @param p.initial is a P-vector specifying the initial model parameters for fitting the DF.
#' @param prior is an optional function specifying the priors on the P model parameters. This function must take a P-dimensional vector \code{p} as input argument and return a scalar proportional to the natural logarithm of the prior probability of the P-dimensional model parameter \code{p}. In other words, \code{prior(p)} is an additative term for the log-likelihood. Note that \code{prior(p)} must be smooth and finite for the full parameter space. If \code{prior=NULL} (default), uniform priors are assumed.
#' @param obs.selection is an optional selection function of a D-vector \code{x}, which specifies the fraction (between 0 and 1) of data with \emph{observed} values \code{x}, whose true value passes the selection function \code{selection} can be observed. Normally this fraction is assumed to be 1 and all the selections are assumed to be specified relative to the true value in \code{selection}. However, if the data were selected, for example, with a sharp cut in the observed values, this can be specified in \code{obs.selection}. The function \code{obs.selection(x)} must be vectorized and return a vector of \code{N} elements, if \code{x} is an \code{N-by-D} matrix (such as the input argument \code{x}).
#' @param obs.sel.cov is an optional \code{D-by-D} matrix, only used if \code{obs.selection} is set. It specifies the mean covariance matrix of the data to estimate how much data has been scattered outside the observing range. If \code{obs.selection} is set, but \code{obs.sel.cov} is not specified, the code estimates \code{obs.sel.cov} from the mean errors of the data (only works for 1D data).
#' @param n.iterations Maximum number of iterations in the repeated fit-and-debias algorithm to evaluate the maximum likelihood.
#' @param correct.lss.bias If \code{TRUE} the \code{distance} values are used to correct for the observational bias due to galaxy clustering (large-scale structure). The overall normalization of the effective folume is chosen such that the expected mass contained in the survey volume is the same as for the uncorrected effective volume.
#' @param lss.weight If \code{correct.lss.bias==TRUE}, this optional function of a \code{P}-vector is the weight-function used for the mass normalization of the effective volume. For instance, to preserve the number of galaxies, choose \code{lss.weight = function(x) 1}, or to perserve the total mass, choose \code{lss.weight = function(x) 10^x} (if the data \code{x} are log10-masses).
#' @param lss.errors is a logical flag specifying whether uncertainties computed via resampling should include errors due to the uncerainty of large-scale structure (LSS). If \code{TRUE} the parameter uncerainties are estimated by refitting the LSS correction at each resampling iteration. This argument is only considered if \code{correct.lss.bias=TRUE} and \code{n.bootstrap>0}.
#' @param n.bootstrap If \code{n.bootstrap} is an integer larger than one, the data is resampled \code{n.bootstrap} times using a non-parametric bootstrapping method to produce more accurate covariances. These covariances are given as matrix and as parameter quantiles in the output list. If \code{n.bootstrap = NULL}, no resampling is performed.
#' @param n.jackknife If \code{n.jackknife} is an integer larger than one, the data is jackknife-resampled \code{n.jackknife} times, removing exactly one data point from the observed set at each iteration. This resampling adds model parameters, maximum likelihood estimator (MLE) bias corrected parameter estimates (corrected to order 1/N). If \code{n.jackknife} is larger than the number of data points N, it is automatically reduced to N.  If \code{n.jackknife = NULL}, no such parameters are deterimed.
#' @param xmin,xmax,dx are \code{P}-element vectors (i.e. scalars for 1-dimensional DF) specifying the points (\code{seq(xmin[i],xmax[i],by=dx[i])}) used for some numerical integrations.
#' @param keep.eddington.bias If \code{TRUE}, the data is not corrected for Eddington bias. In this case no fit-and-debias iterations are performed and the argument \code{n.iterations} will be ignored.
#' @param write.fit If \code{TRUE}, the best-fitting parameters are displayed in the console.
#' @param add.gaussian.errors If \code{TRUE}, Gaussian estimates of the 16 and 84 percentiles of the fitted generative distribution function (gdf) are included in the sublist \code{grid} of the output.
#' @param make.posteriors If \code{TRUE}, posterior probability distributions of the observed data are evaluated from the best fitting model.
#' 
#' @details
#' For a detailed description of the method, please refer to the peer-reviewed publication by Obreschkow et al. 2017 (in prep.).
#'
#' @return The routine \code{dffit} returns a structured list, which can be interpreted by other routines, such as \code{\link{dfwrite}}, \code{\link{dfplot}}, \code{\link{dfplotcov}}, \code{\link{dfplotveff}}. The list contains the following sublists:
#' 
#' \item{data}{contains the input arguments \code{x}, \code{x.err}, \code{r} and \code{lss.weight}, as well as the integers \code{n.data} and \code{n.dim} specifying the number of objects (N) to be fitted and the dimension (D) of the observables. In other words, \code{n.data} and \code{n.dim} are the number of rows and columns of \code{x}, respectively.}
#' 
#' \item{selection}{is a list describing the selection function of the data used when fitting the generative DF. The most important entries in this list are:\cr\cr
#' \code{veff} is a function of a D-dimensional vector, specifying the effective volume associated with an object of properties xval. If LSS is corrected for, i.e. if the argument \code{correct.lss.bias} is set to \code{TRUE}, this function is the final effecive volume, including the effect of LSS.\cr\cr
#' \code{veff.no.lss} is a function of a D-dimensional vector, specifying the effecive volume associate with an object, if LSS were not accounted for. This function is indentical to \code{veff}, if LSS is not accounted for in the fit, i.e. if the argument \code{correct.lss.bias} is set to \code{FALSE}.\cr\cr
#' \code{mode} is an integer specifying the format of the input argument \code{selection}.}
#' 
#' \item{fit}{is a list describing the fitted generative distribution function. Its most important entries are:\cr\cr
#' \code{p.best} is a P-vector giving the most likely model parameters according to the MML method.\cr\cr
#' \code{p.sigma} and \code{p.covariance} are the standard deviations and covariance matrices of the best-fitting parameters in the Gaussian approximation from the Hessian matrix of the modified likelihood function.\cr\cr
#'  \code{gdf} is a function of a D-dimensional vector, which is the generative DF, evaluated at the parameters \code{p.best}.\cr\cr
#'  \code{scd} is a function of a D-dimensional vector, which gives the predicted source counts of the most likely model, i.e. scd(x)=gdf(x)*veff(x).}
#' 
#' \item{posteriors}{is a list specifying the posterior PDFs of the observed data, given the best-fitting model. It contains the following entries:\cr\cr
#' \code{x.mean} is an N-by-D dimensional array giving the D-dimensional means of the posterior PDFs of the N objects.\cr\cr
#' \code{x.mode} is an N-by-D dimensional array giving the D-dimensional modes of the posterior PDFs of the N objects.\cr\cr
#' \code{x.stdev} is an N-by-D dimensional array giving the D-dimensional standard deviations of the posterior PDFs of the N objects.\cr\cr
#' \code{x.random} is an N-by-D dimensional array giving one random D-dimensional value drawn from the posterior PDFs of each of the N objects.}
#' 
#' \item{model}{is a list describing the generative DF used to model the data. The main entries of this list are:\cr\cr
#' \code{gdf(xval,p)} is the generative DF to be fitted, written as phi(x|theta) in the reference publication.\cr\cr
#' \code{gdf.equation} is a text-string representing the analytical equation of \code{gdf}.\cr\cr
#' \code{parameter.names} is a P-vector of expressions (see \code{\link{expression}}), specifying the names of the parameters.\cr\cr
#' \code{n.para} is an integer specifying the number P of model parameters.}
#' 
#' \item{grid}{is a list of arrays with numerical evaluations of different functions on a grid in the D-dimensional observable space. This grid is used for numerical integrations and graphical representations. The most important list entries are:\cr\cr
#' \code{x} is a M-by-D array of M points, defining a regular cartesian grid in the D-dimensional observable space.\cr\cr
#' \code{xmin} and \code{xmax} are D-vectors specifying the lower and upper boundary of the grid in the D-dimensional observable space.\cr\cr
#' \code{dx} is a D-vector specifying the steps between grid points.\cr\cr
#' \code{dvolume} is a number specifying the D-dimensional volume associated with each grid point.\cr\cr
#' \code{n.points} is the number M of grid points.\cr\cr
#' \code{gdf} is an M-vector of values of the best-fitting generative DF at each grid point. The additional entries \code{gdf.error.neg} and \code{gdf.error.pos} define the 68\%-confidence range in the Hessian approximation of the parameter covariances. The optional entries \code{gdf.quantile.#} are different quantiles, generated if the input argument \code{n.bootstrap} is set.\cr\cr
#' \code{veff} is an M-vector of effective volumes at each grid point.\cr\cr
#' \code{scd} is an M-vector of predicted source counts according to the best-fitting model.\cr\cr
#' \code{scd.posterior} is an M-vector of observed source counts derived from the posterior PDFs of each object.\cr\cr
#' \code{effective.counts} is an M-vector of fractional source counts derived from the posterior PDFs of each object.}
#' 
#' \item{options}{is a list of various optional input arguments of \code{dffit}.}
#'
#' @keywords schechter function
#' @keywords mass function
#' @keywords fit
#'
#' @examples
#' # For a quick overview of some key functionalities run
#' dfexample()
#' # with varying integer arguments 1, 2, 3, 4.
#' 
#' # The following examples introduce the basics of dftools step-by step.
#' # First, generate a mock sample of 1000 galaxies with 0.5dex mass errors, drawn from
#' # a Schechter function with the default parameters (-2,11,-1.3):
#' dat = dfmockdata(n=1000, sigma=0.5)
#' 
#' # show the observed and true log-masses (x and x.true) as a function of true distance r
#' plot(dat$r,dat$x,col='grey'); points(dat$r,dat$x.true,pch=20)
#' 
#' # fit a Schechter function to the mock sample without accounting for errors
#' survey1 = dffit(dat$x, dat$veff)
#' 
#' # plot fit and add a black dashed line showing the input MF
#' ptrue = dfmodel(output='initial')
#' mfplot(survey1, xlim=c(1e6,2e12), ylim=c(2e-4,2),
#' show.data.histogram = TRUE, p = ptrue, col.fit = 'purple')
#' 
#' # now, do the same again, while accountting for measurement errors in the fit
#' # this time, the posterior data, corrected for Eddington bias, is shown as black points
#' survey2 = dffit(dat$x, dat$veff, dat$x.err)
#' mfplot(survey2, show.data.histogram = NA, add = TRUE)
#'
#' # show fitted parameter PDFs and covariances with true input parameters as black points
#' dfplotcov(list(survey2,survey1,ptrue),pch=c(20,20,3),col=c('blue','purple','black'),nstd=15)
#'
#' # show effective volume function
#' dfplotveff(survey2)
#'
#' # now create a smaller survey of only 30 galaxies with 0.5dex mass errors
#' dat = dfmockdata(n=30, sigma=0.5)
#' 
#' # fit a Schechter function and determine uncertainties by resampling the data
#' survey = dffit(dat$x, dat$veff, dat$x.err, n.bootstrap = 30)
#' 
#' # show best fit with 68% Gaussian uncertainties from Hessian and posterior data as black points
#' mfplot(survey, show.data.histogram = TRUE, uncertainty.type = 1)
#' 
#' # show best fit with 68% and 95% resampling uncertainties and posterior data as black points
#' mfplot(survey, show.data.histogram = TRUE, uncertainty.type = 3)
#' 
#' # add input model as dashed lines
#' lines(10^survey$grid$x, survey$model$gdf(survey$grid$x,ptrue), lty=2)
#'
#' @author Danail Obreschkow
#'
#' @export

dffit <- function(x,
                  selection = NULL,
                  x.err = NULL,
                  r = NULL,
                  gdf = 'Schechter',
                  p.initial = NULL,
                  prior = NULL,
                  obs.selection = NULL,
                  obs.sel.cov = NULL,
                  n.iterations = 100,
                  correct.lss.bias = FALSE,
                  lss.weight = NULL,
                  lss.errors = TRUE,
                  n.bootstrap = NULL,
                  n.jackknife = NULL,
                  xmin = 5,
                  xmax = 13,
                  dx = 0.01,
                  keep.eddington.bias = FALSE,
                  write.fit = FALSE,
                  add.gaussian.errors = TRUE,
                  make.posteriors = TRUE) {

  # Set timer
  tStart = Sys.time()
  
  # Initialize main dataframe
  survey = list(data = list(x = x, x.err = x.err, r = r, lss.weight = lss.weight),
                selection = list(),
                model = list(),
                grid = list(xmin = xmin, xmax = xmax, dx = dx),
                fit = list(),
                options = list(p.initial = p.initial, prior = prior,
                               n.iterations = n.iterations,
                               n.bootstrap = n.bootstrap, n.jackknife = n.jackknife,
                               keep.eddington.bias = keep.eddington.bias,
                               correct.lss.bias = correct.lss.bias,
                               lss.weight = lss.weight, lss.errors = lss.errors),
                tmp = list(selection = selection, gdf = gdf, obs.selection = obs.selection, obs.sel.cov = obs.sel.cov))
  
  # Check and pre-process input arguments
  survey = .handle.input(survey)
  survey = .make.veff(survey)
  survey = .make.grid(survey)
  survey = .make.obs.filter(survey)
  survey$tmp$rho.observed = .make.prior.pdfs(survey)
  
  # Find most likely generative model
  fit = .corefit(survey)
  if (survey$options$correct.lss.bias) survey$selection$veff = .get.veff.lss(survey,fit$p.best)
  survey$fit = list(p.best = fit$p.best, p.sigma = sqrt(diag(fit$p.covariance)), p.covariance = fit$p.covariance,
                    gdf = function(x) survey$model$gdf(x,fit$p.best),
                    scd = function(x) survey$model$gdf(x,fit$p.best)*survey$selection$veff(x),
                    lnL = fit$lnL,
                    ln.evidence = fit$ln.evidence,
                    status = fit$status)
  survey$grid$gdf = c(survey$fit$gdf(survey$grid$x))
  survey$grid$veff = c(survey$selection$veff(survey$grid$x))
  survey$grid$scd = survey$grid$gdf*survey$grid$veff
  
  # Determine Gaussian uncertainties
  if (survey$fit$status$converged & add.gaussian.errors) survey = .add.Gaussian.errors(survey)
  
  # Resample to determine more accurate uncertainties with quantiles
  if (survey$fit$status$converged & !is.null(survey$options$n.bootstrap)) survey = .resample(survey)
  
  # Resample to determine more accurate uncertainties with quantiles
  if (survey$fit$status$converged & !is.null(survey$options$n.jackknife)) survey = .jackknife(survey)
  
  # Make posteriors
  if (survey$fit$status$converged & !is.null(survey$data$x.err) & make.posteriors) survey = .dfposteriors(survey)
  
  # Write best fitting parameters
  if (write.fit) dfwrite(survey)

  # Finalize
  survey$fit$status$walltime = as.double(Sys.time())-as.double(tStart)
  survey[which(names(survey)=='tmp')] = NULL
  invisible(survey)

}

.handle.input <- function(survey) {
  
  # Handle x
  if (length(survey$data$x)<1) stop('Give at least one data point.')
  if (is.null(dim(survey$data$x))) {
    survey$data$x = cbind(as.vector(survey$data$x)) # make col-vector
  } else if (length(dim(survey$data$x))==1) {
    survey$data$x = cbind(survey$data$x)
  } else if (length(dim(survey$data$x))!=2) {
    stop('x cannot have more than two dimensions.')
  }
  survey$data$n.data = dim(survey$data$x)[1]
  survey$data$n.dim = dim(survey$data$x)[2]
  for (i in seq(survey$data$n.dim)) {
    if (any(survey$data$x[,i]<survey$grid$xmin[i])) stop('xmin cannot be larger than smallest observed value x.')
    if (any(survey$data$x[,i]>survey$grid$xmax[i])) stop('xmax cannot be smaller than largest observed value x.')
  }
  
  # Handle x.err
  if (!is.null(survey$data$x.err)) {
    if (is.function(survey$data$x.err)) {
      # ok
    } else {
      if (all(survey$data$x.err==0)) {
        survey$data$x.err = NULL
      } else {
        if (is.null(dim(survey$data$x.err))) {
          survey$data$x.err = cbind(as.vector(survey$data$x.err)) # make col-vector
        } else if (length(dim(survey$data$x.err))==1) {
          survey$data$x.err = cbind(survey$data$x.err)
        } else  if (length(dim(survey$data$x.err))==2) {
          if (any(dim(survey$data$x)!=dim(survey$data$x.err))) stop('Size of x.err not compatible with size of x.')
        } else if (length(dim(survey$data$x.err))==3) {
          if (survey$data$n.dim==1) stop('For one-dimensional distribution function x.err cannot have 3 dimensions.')
          if (!(dim(survey$data$x.err)[1]==survey$data$n.data & dim(survey$data$x.err)[2]==survey$data$n.dim & dim(survey$data$x.err)[3]==survey$data$n.dim)) {
            stop('Size of x.err not compatible with size of x.')
          }
        } else {
          stop('x.err cannot have more than three dimensions.')
        }
        if (min(survey$data$x.err)<0) stop('All values of x.err must be positive.')
        if (length(dim(survey$data$x.err))==2) {
          for (i in seq(survey$data$n.dim)) {
            if (any(survey$data$x[,i]-survey$data$x.err[,i]<survey$grid$xmin[i])) stop('xmin cannot be larger than smallest observed value x-x.err.')
            if (any(survey$data$x[,i]+survey$data$x.err[,i]>survey$grid$xmax[i])) stop('xmax cannot be smaller than largest observed value x+x.err.')
          }
        } else if (length(dim(survey$data$x.err))==3) {
          for (i in seq(survey$data$n.dim)) {
            if (any(survey$data$x[,i]-sqrt(survey$data$x.err[,i,i])<survey$grid$xmin[i])) stop('xmin cannot be larger than smallest observed value x-x.err.')
            if (any(survey$data$x[,i]+sqrt(survey$data$x.err[,i,i])>survey$grid$xmax[i])) stop('xmax cannot be smaller than largest observed value x+x.err.')
          }
        }
      }
    }
  }
  
  # Handle distance
  if (!is.null(survey$data$r)) {
    survey$data$r = as.vector(survey$data$r)
    if (length(survey$data$r)!=survey$data$n.data) stop('The number of distance values must be equal to the number of data points (=number of columns of x).')
    if (min(survey$data$r)<=0) stop('All distance values must be positive.')
  }
  
  # Handle correct.lss.bias
  if (survey$options$correct.lss.bias) {
    if (is.null(survey$data$r)) stop('Distances must be given of correct.lss.bias = TRUE.')
  }
  
  # Handle priors
  if (is.null(survey$options$prior)) {
    survey$options$prior = function(p) 0
  }
  
  # Handle gdf
  if (is.function(survey$tmp$gdf)) {
    if (is.null(survey$options$p.initial)) stop('For user-defined distribution functions initial parameters must be given.')
    survey$model$gdf = survey$tmp$gdf
    survey$model$gdf.equation = NULL
    survey$model$parameter.names = NULL
  } else {
    if (is.null(survey$options$p.initial)) survey$options$p.initial = dfmodel(output = 'initial', type = survey$tmp$gdf)
    survey$model$gdf = function(x,p) {dfmodel(x, p, type = survey$tmp$gdf)}
    survey$model$gdf.equation = dfmodel(output = 'equation', type = survey$tmp$gdf)
    survey$model$parameter.names = dfmodel(output = 'names', type = survey$tmp$gdf)
  }
  survey$model$n.para = length(survey$options$p.initial)
  
  # Handle n.iterations
  if (is.null(survey$options$n.iterations)) stop('n.iterations must be a positive integer.')
  if (survey$options$n.iterations<1) stop('n.iterations must be a positive integer.')
  
  # Handle n.bootstrap
  if (!is.null(survey$options$n.bootstrap)) {
    if (survey$options$n.bootstrap<2) stop('n.bootstrap must be 2 or larger.')
  }
  
  # Handle n.jackknife
  if (!is.null(survey$options$n.jackknife)) {
    if (survey$options$n.jackknife<2) stop('n.jackknife must be 2 or larger.')
  }
  
  invisible(survey)
}

#' @export
.make.veff = function(survey) {
  
  # Generates the function veff.function(xval) from various selection function types.
  
  # Handle selection
  if (is.null(survey$tmp$selection)) stop('A selection function must be given.')
  s = survey$tmp$selection
  x = survey$data$x
  r = survey$data$r
  n.dim = survey$data$n.dim
  n.data = survey$data$n.data
  xmin = apply(x,2,min)
  xmax = apply(x,2,max)
  mode = NULL
  veff.values = NULL
  veff.userfct = NULL
  veff.function = NULL
  f.function = NULL
  dVdr = NULL
  rmin = NULL
  rmax = NULL
  
  # Mode 1: Constant effective volume inside observed domain
  if (is.double(s)) {
    if (length(s)==1) {
      if (s<=0) stop('s = Vconstant mustt be positive.')
      mode = 1
      veff.function.elemental = function(xval) {
        if (any(xval<xmin) | any(xval>xmax)) {
          return(0)
        } else {
          return(s)
        }
      }
    }
  }
  
  # Mode 2: Interpolated effective volume inside observed domain
  if (is.double(s)) {
    if (length(s)==n.data) {
      if (min(s)<=0) stop('All values of selection (=Veff) must be positive.')
      mode = 2
      veff.values = s
      if (n.dim==1) {
        vapprox = function(xval) {
          f = approxfun(x[,1],1/veff.values,rule=2)
          return(1/f(xval))
        }
        veff.max = max(veff.values)
        veff.function.elemental = function(xval) {
          if (xval<xmin) {
            return(0)
          } else if (xval>xmax) {
            return(veff.max)
          } else {
            return(vapprox(xval))
          }
        }
      } else if (n.dim==2) {
        vapprox = function(xval) {
          z = 1/(akima::interp(x[,1],x[,2],1/veff.values,xval[1],xval[2],duplicate='mean'))$z
          if (is.na(z)) {return(0)} else {return(z)}
        }
        veff.function.elemental = function(xval) {
          if (any(xval<xmin) | any(xval>xmax)) {
            return(0)
          } else {
            return(vapprox(xval))
          }
        }
      } else {
        stop('Linear interpolation of Veff not implemented for DF with more than 2 dimensions. Use a different selection type.')
      }
    }
  }
  
  # Mode 3: Effective volume given directly
  if (is.function(s)) {
    mode = 3
    test = try(s(rbind(x)))
    if (is.double(test) & length(test)==n.data) {
      veff.function = s
    } else {
      veff.function.elemental = s
    }
    veff.userfct = s
  }
  
  # Mode 4: Hybrid of 2 and 3
  if (is.list(s)) {
    if (length(s)==2) {
      if (is.double(s[[1]]) & is.function(s[[2]])) {
        veff.values = s[[1]]
        test = try(s[[2]](rbind(x)))
        if (is.double(test) & length(test)==n.data) {
          veff.userfct = function(xval) s[[2]](rbind(xval))
        } else {
          veff.userfct = s[[2]]
        }
        if (min(veff.values)<=0) stop('All values of selection (=Veff) must be positive.')
        mode = 4
        if (n.dim==1) {
          vapprox = function(xval) {
            f = approxfun(x[,1],1/veff.values,rule=2)
            return(1/f(xval))
          }
        } else if (n.dim==2) {
          vapprox = function(xval) {
            z = 1/(akima::interp(x[,1],x[,2],1/veff.values,xval[1],xval[2],duplicate='mean'))$z
            if (is.na(z)) {return(veff.userfct(xval))} else {return(z)}
          }
        } else {
          stop('Linear interpolation of Veff not implemented for DF with more than 2 dimensions. Use a different selection type.')
        }
        veff.function.elemental = function(xval) {
          if (any(xval<xmin) | any(xval>xmax)) {
            return(veff.userfct(xval))
          } else {
            return(vapprox(xval))
          }
        }
      }
    }
  }
  
  # Mode 5: Selection given as {f(x,r), dVdr(r), rmin, rmax}
  if (is.list(s)) {
    if (length(s)==4) {
      if (is.function(s[[1]]) & is.function(s[[2]]) & is.double(s[[3]]) & is.double(s[[4]])) {
        mode = 5
        rmin = s[[3]]
        rmax = s[[4]]
        if (!is.null(r)) {
          if (rmin>min(r)) stop('rmin cannot in selection = list(f, dVdr, rmin, rmax) cannot be larger than minimal observed r')
          if (rmax<max(r)) stop('rmax cannot in selection = list(f, dVdr, rmin, rmax) cannot be smaller than maximal observed r')
        }
        test = try(s[[1]](NA,NA)*s[[2]](NA),silent=TRUE)
        if (!is.double(test)) stop('In the argument selection = list(f, dVdr, rmin, rmax), the functions f(xval,r) and dVdr(r) must work if r is a vector.')
        f.function = s[[1]]
        dVdr = s[[2]]
        veff.function.elemental = function(xval) {
          f = function(r) {s[[1]](xval,r)*s[[2]](r)}
          return(integrate(f,rmin,rmax,stop.on.error=FALSE)$value)
        }
      }
    }
  }
  
  if (is.null(mode)) stop('Unknown selection format.')
  
  # apply to all elements
  if (is.null(veff.function)) {
    veff.function = function(xval) {
      if (length(dim(xval))==2) {
        return(apply(xval,1,veff.function.elemental))
      } else if (length(dim(xval))==0) {
        if (n.dim==1) {
          return(sapply(xval,veff.function.elemental))
        } else {
          if (length(xval)!=n.dim) stop('Incorrect argument.')
          return(veff.function.elemental(xval))
        }
      } else {
        stop('Incorrect argument.')
      }
    }
  }
  
  survey$selection$veff = veff.function
  survey$selection$veff.no.lss = veff.function
  survey$selection$veff.input.values = veff.values
  survey$selection$veff.input.function = veff.userfct
  survey$selection$f = f.function
  survey$selection$dVdr = dVdr
  survey$selection$rmin = rmin
  survey$selection$rmax = rmax
  survey$selection$mode = mode
  
  invisible(survey)
}

.get.veff.lss = function(survey,p) {
  
  if (survey$selection$mode!=5) {
    cat('To correct LSS bias, the selection must be specified in the format\n')
    cat('selection = list(f, dVdr, rmin, rmax)\n')
    stop()
  }
  
  # initialize input
  x = survey$data$x
  r = survey$data$r
  n.dim = survey$data$n.dim
  n.data = survey$data$n.data
  rmin = survey$selection$rmin
  rmax = survey$selection$rmax
  simpson.integration = n.dim==1
  
  # evaluate integrals
  integrand.lss = function(x,r) {survey$selection$f(x,r)*survey$model$gdf(x,p)}
  integral = array(NA,n.data)
  if (simpson.integration) {
    for (i in seq(n.data)) {
      integral[i] = integrate(integrand.lss,survey$grid$xmin,survey$grid$xmax,r[i],stop.on.error=FALSE)$value
    }
  } else {
    for (i in seq(n.data)) {
      integral[i] = sum(integrand.lss(survey$grid$x,r[i]))*survey$grid$dvolume
    }
  }
  
  # make veff.lss function
  veff.lss.function.elemental = function(xval) {
    list = survey$selection$f(xval,r)>0
    return(sum(survey$selection$f(xval,r[list])/integral[list]))
  }
  veff.lss.scale = Vectorize(veff.lss.function.elemental)
  
  # renormalize veff to conserve total
  if (is.null(survey$options$lss.weight)) {
    weight = function(x) 1
  } else {
    weight = survey$options$lss.weight
  }
  int.ref = function(x) survey$selection$veff.no.lss(x)*survey$model$gdf(x,p)*weight(x)
  int.exp = function(x) veff.lss.scale(x)*survey$model$gdf(x,p)*weight(x)
  if (simpson.integration) {
    reference = integrate(int.ref,survey$grid$xmin,survey$grid$xmax,stop.on.error=FALSE)$value
    expectation = integrate(int.exp,survey$grid$xmin,survey$grid$xmax,stop.on.error=FALSE)$value
  } else {
    reference = sum(int.ref(survey$grid$x))*survey$grid$dvolume
    expectation = sum(int.exp(survey$grid$x))*survey$grid$dvolume
  }
  normalization.factor = reference/expectation
  veff.lss = function(x) veff.lss.scale(x)*normalization.factor
  
  # return output
  return(veff.lss)
}

.make.grid = function(survey) {
  
  n.dim = survey$data$n.dim
  if (length(survey$grid$xmin)!=n.dim) stop('xmin must be a P-element vector, where P is the number of columns of x.')
  if (length(survey$grid$xmax)!=n.dim) stop('xmax must be a P-element vector, where P is the number of columns of x.')
  if (length(survey$grid$dx)!=n.dim) stop('dx must be a P-element vector, where P is the number of columns of x.')
  x.grid = list()
  nx = array(NA,n.dim)
  for (i in seq(n.dim)) {
    if (survey$grid$xmax[i]<survey$grid$xmin[i]+survey$grid$dx[i]) stop('xmax cannot be smaller than xmin+dx.')
    x.grid[[i]] = seq(survey$grid$xmin[i],survey$grid$xmax[i],survey$grid$dx[i])
    nx[i] = length(x.grid[[i]])
  }
  x.mesh.dv = prod(survey$grid$dx)
  x.mesh = array(NA,c(prod(nx),n.dim))
  k1 = 1
  k2 = prod(nx)
  for (i in seq(n.dim)) {
    k2 = k2/nx[i]
    x.mesh[,i] = rep(rep(x.grid[[i]],each=k1),k2)
    k1 = k1*nx[i]
  }
  survey$grid$x = x.mesh
  survey$grid$dvolume = x.mesh.dv
  survey$grid$veff = survey$selection$veff(x.mesh)
  survey$grid$n.points = dim(x.mesh)[1]
  invisible(survey)
}

.make.obs.filter = function(survey) {
  
  # Make filter function that corrects the density of true x-values for the fact that a fraction of these values have been scattered outside the observed range of x
  # (important, for example, if there is a sharp sample cut in the observed input values of x)
  
  # Handle observational selection
  survey$selection$obs.selection = survey$tmp$obs.selection
  if (!is.null(survey$tmp$obs.selection)) {
    if (!is.function(survey$tmp$obs.selection)) {
      stop('obs.selection must be a function of a D-dimensional vector')
    }
  }
  if (is.null(survey$selection$obs.selection)) {
    survey$selection$obs.sel.cov = NULL
  } else {
    if (is.null(survey$tmp$obs.sel.cov)) {
      survey$selection$obs.sel.cov = mean(survey$data$x.err^2)
    } else {
      survey$selection$obs.sel.cov = survey$tmp$obs.sel.cov
    }
  }
  
  # put filter on grid
  n.mesh = dim(survey$grid$x)[1]
  n.dim = dim(survey$data$x)[2]
  survey$grid$obs.filter = rep(1,n.mesh)
  if (!is.null(survey$selection$obs.selection)) {
    inv.covariance = solve(survey$selection$obs.sel.cov)
    gauss = rep(NA,n.mesh)
    for (i in seq(n.mesh)) {
      mu = t(array(survey$grid$x[i,],c(n.dim,n.mesh)))
      dx = survey$grid$x-mu
      for (j in seq(n.mesh)) {
        gauss[j] = exp(-dx[j,]%*%inv.covariance%*%dx[j,]/2)
      }
      survey$grid$obs.filter[i] = sum(gauss*survey$selection$obs.selection(survey$grid$x))/sum(gauss)
    }
  }
  
  invisible(survey)
  
}

#' @export
.corefit <- function(survey, supress.warning = FALSE) {
  
  # simplify input
  x = survey$data$x
  gdf = survey$model$gdf
  p.initial = survey$options$p.initial
  x.mesh = survey$grid$x
  x.mesh.dv = survey$grid$dvolume
  keep.eddington.bias = survey$options$keep.eddington.bias
  parameter.prior = survey$options$prior
  obs.filter = survey$grid$obs.filter
  
  # get array sizes
  n.data = dim(x)[1]
  n.dim = dim(x)[2]
  n.mesh = dim(x.mesh)[1]
  
  # Input handling
  if (!survey$options$correct.lss.bias) veff.mesh = c(survey$grid$veff)
  if (is.null(survey$data$x.err) & !survey$options$correct.lss.bias) {
    n.iterations = 1
  } else {
    n.iterations = survey$options$n.iterations
  }
  
  # Iterative algorithm
  running = TRUE
  k = 0
  offset = 0
  chain = array(NA,c(n.iterations,length(p.initial)+1))
  
  # definie debiasing function
  debias = function(p) {
    rho.unbiased = array(0,n.mesh)
    if (is.null(survey$data$x.err) | keep.eddington.bias) {
      for (i in seq(n.data)) {
        rho.unbiased = rho.unbiased+survey$tmp$rho.observed[[i]]
      }
    } else {
      # predicted source counts (up to a factor x.mesh.dv)
      prior = gdf(x.mesh,p)*veff.mesh
      prior[!is.finite(prior)] = 0
      prior = pmax(0,prior)
      for (i in seq(n.data)) {
        rho.corrected = survey$tmp$rho.observed[[i]]*prior
        rho.corrected = rho.corrected/(sum(rho.corrected)*x.mesh.dv)
        rho.unbiased = rho.unbiased+rho.corrected
      }
    }
    return(rho.unbiased)
  }
  
  while (running) {
    
    k = k+1
    
    # determine veff LSS
    if (survey$options$correct.lss.bias) {
      veff.lss = .get.veff.lss(survey,p.initial)
      veff.mesh = c(veff.lss(x.mesh))
    }
    
    # make unbiased source density function
    rho.unbiased = debias(p.initial)
    
    # make -ln(L)
    neglogL = function(p) {
      phi = gdf(x.mesh,p)*obs.filter
      # safety operations (adding about 50% computation time)
      phi[!is.finite(phi)] = 0
      phi = pmax(.Machine$double.xmin,phi) # also vectorizes the array
      # end safety operations
      return(sum(phi*veff.mesh-log(phi)*rho.unbiased)*x.mesh.dv-offset-parameter.prior(p))
    }
    
    # test
    if (all(!is.finite(gdf(x.mesh,p.initial)))) stop('cannot evaluate GDF at initial parameters provided')
    test = try(neglogL(p.initial),silent=TRUE)
    if (!is.finite(test)) stop('cannot evaluate likelihood at initial parameters provided')
    
    # maximize ln(L)
    opt = optim(p.initial,neglogL,hessian=TRUE,control=list(maxit=1e5)) # reltol=1e-20,abstol=1e-20,
    offset = opt$value+offset
    chain[k,] = c(opt$par,opt$value)
    
    # assess convergence
    if ((is.null(survey$data$x.err) & !survey$options$correct.lss.bias) | keep.eddington.bias) { # exist without extra iterations
      converged = opt$convergence==0
      running = FALSE
    } else {
      # asses convergence
      if (k==1) {
        converged = FALSE
        d.old = Inf
      } else {
        d = abs(opt$value-value.old)
        converged = d>=d.old
        d.old = d
      }
      value.old = opt$value
      
      if (converged) {
        running = FALSE
      } else if (k == n.iterations) {
        converged = FALSE
        running = FALSE
        if (!supress.warning) {
          cat('WARNING: Maximum number of iteration reached. Consider increasing n.iterations and/or providing better initial parameters.\n')
        }
      }
      
      # prepare initial values for next iteration
      p.initial = opt$par
    }
  }
  
  # make output
  if (det(opt$hessian)<1e-12) {
    converged = FALSE
    cov = NULL
    if (!supress.warning) {
      cat('WARNING: Fit ill-conditionned. Consider providing better initial parameters or selection arguments.\n')
    }
    ln.evidence = FALSE
  } else {
    cov = solve(opt$hessian)
    n.para = length(opt$par)
    ln.evidence = -offset+0.5*n.para*log(2*pi)+0.5*log(det(cov))
  }
  
  # finalize output
  fit = list(p.best = opt$par, p.covariance = cov, lnL = function(p) -neglogL(p),
             status = list(n.iterations = k, converged = converged, chain = chain[1:k,]), ln.evidence = ln.evidence)
  return(fit)
  
}

.make.prior.pdfs = function(survey) {
  
  # simplify input
  x = survey$data$x
  x.err = survey$data$x.err
  x.mesh = survey$grid$x
  x.mesh.dv = survey$grid$dvolume
  
  # get array sizes
  n.data = dim(x)[1]
  n.dim = dim(x)[2]
  n.mesh = dim(x.mesh)[1]
  
  rho.observed = list()

  if (is.null(x.err)) {
  
    for (i in seq(n.data)) {
      difference = colSums(abs(t(x.mesh)-x[i,]))
      index = which.min(difference)[1]
      rho.observed[[i]] = array(0,n.mesh)
      rho.observed[[i]][index] = rho.observed[[i]][index]+1/x.mesh.dv
    }
    
  } else {
    
    if (is.function(x.err)) {
      
      for (i in seq(n.data)) {
        rho.observed[[i]] = c(x.err(c(x.mesh),i))
      }
      
    } else {
    
      # Make inverse covariances
      invC = array(NA,c(n.data,n.dim,n.dim))
      if (length(dim(x.err))==2) {
        if (!(dim(x.err)[1]==n.data & dim(x.err)[2]==n.dim)) {
          stop('Unknown format for x.err in .corefit.')
        }
        if (n.dim==1) {
          for (i in seq(n.data)) {
            invC[i,,] = 1/max(survey$grid$dx/2.35,x.err[i,])^2
          }
        } else {
          for (i in seq(n.data)) {
            invC[i,,] = diag(1/max(survey$grid$dx/2.35,x.err[i,])^2)
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
      
      # make a priori observed PDFs
      rho.observed = list()
      for (i in seq(n.data)) {
        d = x[i,]-t(x.mesh)
        rho.observed[[i]] = exp(-colSums(d*(invC[i,,]%*%d))/2)
      }
    }
    
    # normalize all prior PDFs
    for (i in seq(n.data)) {
      rho.observed[[i]] = rho.observed[[i]]/sum(rho.observed[[i]])/x.mesh.dv
    }
  }
  
  return(rho.observed)
  
}

.add.Gaussian.errors <- function(survey) {

  cov = survey$fit$p.covariance
  eig = eigen(cov)
  np = survey$model$n.para
  nx = survey$grid$n.points
  
  # sample surface of covariance ellipsoid
  nsteps = 500
  y.new = array(NA,c(nx,nsteps+2*np))
  for (i in seq(nsteps+2*np)) {
    if (i<=nsteps) {
      e = rnorm(np)
      e = e/sqrt(sum(e^2))
    } else {
      e = array(0,np)
      if (i>nsteps+np) {
        e[i-nsteps-np] = 1
      } else {
        e[i-nsteps] = -1
      }
    }
    v = c(eig$vectors%*%(sqrt(eig$values)*e))
    p.new = survey$fit$p.best+v
    y.new[,i] = survey$model$gdf(survey$grid$x,p.new)
  }

  survey$grid$gdf.error.neg = array(NA,nx)
  survey$grid$gdf.error.pos = array(NA,nx)
  for (i in seq(nx)) {
    survey$grid$gdf.error.neg[i] = survey$grid$gdf[i]-min(c(y.new[i,],Inf),na.rm=TRUE)
    survey$grid$gdf.error.pos[i] = max(c(y.new[i,],-Inf),na.rm=TRUE)-survey$grid$gdf[i]
  }

  invisible(survey)
}

.jackknife = function(survey) {
  
  # input handling
  n.data = survey$data$n.data
  np = survey$model$n.para
  if (n.data<=np) stop('Jackknifing requires more data points than model parameters.')
  n.jn = min(n.data,survey$options$n.jackknife)
  if (survey$options$n.jackknife<2) stop('n.jackknife must be at least 2.')
  
  # set up jackknife survey
  jn = survey
  jn$options$p.initial = survey$fit$p.best
  jn$grid$veff = survey$grid$veff*(n.data-1)/n.data
  jn$data$n.data = n.data-1
  jn$options$n.iterations = 1
  
  # run jackknifing resampling
  rejected = sample(seq(n.data),size=n.jn,replace=FALSE)
  p.new = array(NA,c(n.jn,np))
  for (i in seq(n.jn)) {
    cat(sprintf('\rJackknifing %4.2f%%',100*i/n.jn))
    list = setdiff(seq(n.data),rejected[i])
    jn$data$x = as.matrix(survey$data$x[list,])
    jn$tmp$rho.observed = survey$tmp$rho.observed
    jn$tmp$rho.observed[[i]] = NULL
    p.new[i,] = .corefit(jn, supress.warning = TRUE)$p.best
  }
  cat('\r')
  
  # correct estimator bias
  p.reduced = apply(p.new, 2, mean, na.rm = T)
  survey$fit$p.best.mle.bias.corrected = n.data*survey$fit$p.best-(n.data-1)*p.reduced
  survey$fit$gdf.mle.bias.corrected = function(x) survey$model$gdf(x,survey$fit$p.best.mle.bias.corrected)
  survey$fit$scd.mle.bias.corrected = function(x) survey$fit$gdf.mle.bias.corrected(x)*survey$selection$veff(x)
  survey$grid$gdf.mle.bias.corrected = survey$fit$gdf.mle.bias.corrected(survey$grid$x)
  survey$grid$scd.mle.bias.corrected = survey$fit$scd.mle.bias.corrected(survey$grid$x)
  
  invisible(survey)
  
}

.resample = function(survey) {
  
  # input handling
  n.data = survey$data$n.data
  np = survey$model$n.para
  if (n.data<2*np) stop('Bootstrapping requires at least twice as many data points as model parameters.')
  if (survey$options$n.bootstrap<2) stop('n.bootstrap must be at least 2.')
  
  # set up resample survey
  b = survey
  b$options$p.initial = survey$fit$p.best
  b$selection$veff = b$selection$veff.no.lss = survey$selection$veff
  b$options$correct.lss.bias = survey$options$correct.lss.bias & survey$options$lss.errors
  
  # randomly resample and refit the DF
  p.new = array(NA,c(survey$options$n.bootstrap,np))
  for (iteration in seq(survey$options$n.bootstrap)) {
    cat(sprintf('\rResampling: %4.2f%%',100*iteration/survey$options$n.bootstrap))
    b$data$n.data = max(np,rpois(1,n.data))
    s = sample.int(n.data,b$data$n.data,replace=TRUE)
    b$data$x = as.matrix(survey$data$x[s,])
    b$tmp$rho.observed = list()
    for (i in seq(b$data$n.data)) {
      b$tmp$rho.observed[[i]] = survey$tmp$rho.observed[[s[i]]]
    }
    if (!is.null(survey$data$r)) b$data$r = survey$data$r[s]
    p.new[iteration,] = .corefit(b, supress.warning = TRUE)$p.best
  }
  cat('\r')
  
  # compute covariance
  survey$fit$p.covariance.resample = array(NA,c(np,np))
  for (i in seq(np)) {
    for (j in seq(np)) {
      survey$fit$p.covariance.resample[i,j] = var(p.new[1:iteration,i],p.new[1:iteration,j],na.rm=TRUE)
    }
  }
  
  # make parameter quantiles
  q = c(0.02,0.16,0.84,0.98)
  p.quant = array(NA,c(length(q),np))
  for (i in seq(np)) {
    p.quant[,i] = quantile(p.new[,i],q,names=FALSE,na.rm=TRUE)
  }
  survey$fit$p.quantile.02 = p.quant[1,]
  survey$fit$p.quantile.16 = p.quant[2,]
  survey$fit$p.quantile.84 = p.quant[3,]
  survey$fit$p.quantile.98 = p.quant[4,]

  # make DF quantiles
  s = array(NA,c(survey$options$n.bootstrap,survey$grid$n.points))
  for (i in seq(survey$options$n.bootstrap)) {
    s[i,] = survey$model$gdf(survey$grid$x,p.new[i,])
  }
  y.quant = array(NA,c(4,survey$grid$n.points))
  for (i in seq(survey$grid$n.points)) {
    list = !is.na(s[,i]) & is.finite(s[,i]) & (s[,i]>0)
    y.quant[,i] = quantile(s[list,i],q,names=FALSE)
  }
  survey$grid$gdf.quantile.02 = y.quant[1,]
  survey$grid$gdf.quantile.16 = y.quant[2,]
  survey$grid$gdf.quantile.84 = y.quant[3,]
  survey$grid$gdf.quantile.98 = y.quant[4,]
  
  invisible(survey)

}

.dfposteriors <- function(survey) {
  
  if (is.null(survey$data$x.err)) stop('Posterior data PDFs can only be produced if the data is uncertain, i.e. if x.err is given.')
  
  # Input handling
  x = survey$data$x
  x.err = survey$data$x.err
  x.mesh = survey$grid$x
  x.mesh.dv = survey$grid$dvolume
  n.data = dim(x)[1]
  n.dim = dim(x)[2]
  n.mesh = dim(x.mesh)[1]
  
  # Make priors
  prior = survey$grid$scd
  prior[!is.finite(prior)] = 0
  prior = pmax(0,prior)
  
  # produce posteriors
  m0 = m1 = md = array(NA,c(n.data,n.dim))
  rho.unbiased = rho.unbiased.sqr = array(0,dim(x.mesh)[1])
  for (i in seq(n.data)) {
    
    # make posterior PDF for data point i
    rho.corrected = survey$tmp$rho.observed[[i]]*prior
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