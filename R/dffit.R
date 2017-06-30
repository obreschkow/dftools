#' Fit a generative distribution function, such as a galaxy mass function
#'
#' This function fits galaxy mass function (MF) to a discrete set of \code{N} galaxies with noisy data. More generally, \code{dffit} finds the most likely \code{P}-dimensional distribution function (DF) generating \code{N} objects \code{i=1,...,N} with uncertain measurements \code{P} observables. For instance, if the objects are galaxies, it can fit a MF (\code{P=1}), a mass-size distribution (\code{P=2}) or the mass-spin-morphology distribution (\code{P=3}). A full description of the algorithm can be found in Obreschkow et al. (2017).
#'
#' @importFrom akima interp
#' @importFrom stats optim rpois quantile approxfun cov
#'
#' @param x Normally \code{x} is a \code{N}-element vector, representing the log-masses (log10(M/Msun)) of \code{N} galaxies. More generally, \code{x} can be either a vector of \code{N} elements or a matrix of \code{N-by-P} elements, containing the values of one or \code{P} observables of \code{N} objects, respectively.
#' @param selection Specifies the effective volume \code{Veff(xval)} in which a galaxy of log-mass \code{xval} can be observed; or, more generally, the volume in which an object of observed values \code{xval[1:P]} can be observed. This volume can be specified in five ways: (1) If \code{selection} is a single positive number, it will be interpreted as a constant volume, \code{Veff(xval)=selection}, in which all galaxies are fully observable. \code{Veff(xval)=0} is assumed outside the "observed domain". This domain is defined as \code{min(x)<=xval<=max(x)} for one observable (\code{P=1}), or as \code{min(x[,j])<=xval[j]<=max(x[,j])} for all \code{j=1,...,P} if \code{P>1}. This mode can be used for volume-complete surveys or for simulated galaxies in a box. (2) If \code{selection} is a vector of \code{N} elements, they will be interpreted as the effective volumes for each of the \code{N} galaxies. \code{Veff(xval)} is interpolated (linearly in \code{1/Veff}) for other values \code{xval}. \code{Veff(xval)=0} is assumed outside the observed domain. (3) \code{selection} can be a function of \code{P} variables, which directly specifies the effective volume for any \code{xval}, i.e. \code{Veff(xval)=selection(xval)}. (4) \code{selection} can also be a list (\code{selection = list(veff.values, veff.userfct)}) of an \code{N}-element vector \code{Veff.values} and a \code{P}-dimensional function \code{veff.userfct}. In this case, the effective volume is computed using a hybrid scheme of modes (2) and (3): \code{Veff(xval)} will be interpolated from the \code{N} values of \code{Veff.values} inside the observed domain, but set equal to \code{veff.userfct} outside this domain. (5) Finally, \code{selection} can be a list of two functions and one optional 2-element vector: \code{selection = list(f, dVdr, rmin, rmax)}, where \code{f = function(xval,r)} is the isotropic selection function and \code{dVdr = function(r)} is the derivative of the total survey volume as a function of comoving distance \code{r}. The scalars \code{rmin} and \code{rmax} (can be \code{0} and \code{Inf}) are the minimum and maximum comoving distance limits of the survey. Outside these limits \code{Veff=0} will be assumed.
#' @param x.err Optional vector or array specifying the observational errors of \code{x}. If \code{x} is a vector then \code{x.err} must also be a vector of same length. Its elements are interpreted as the standard deviations of Gaussian uncertainties in \code{x}. If \code{x} is a \code{N-by-P} matrix representing \code{N} objects with \code{P} observables, then \code{x.err} must be either a \code{N-by-P} matrix or a \code{N-by-P-by-P} array. In the first case, the elements \code{x.err[i,]} are interpreted as the standard deviations of Gaussian uncertainties on \code{x[i,]}. In the second case, the \code{P-by-P} matrices \code{x.err[i,,]} are interpreted as the covariance matrices of the \code{P} observed values \code{x[i,]}.
#' @param r Optional vector of \code{N} elements specifying the comoving distances of the \code{N} galaxies. This vector is only needed if \code{correct.lss.bias = TRUE}.
#' @param gdf Either a string or a function specifying the DF to be fitted. A string is interpreted as the name of a predefined mass function (i.e. functions of one obervable, \code{P=1}). Available options are \code{'Schechter'} for Schechter function (3 parameters), \code{'PL'} for a power law (2 parameters), or \code{'MRP'} for an MRP function (4 parameters). Alternatively, \code{gdf = function(xval,p)} can be any function of the \code{P} observable(s) \code{xval} and a list of parameters \code{p}. IMPORTANT: The function \code{gdf(xval,p)} must be fully vectorized in \code{xval}, i.e. it must output a vector of \code{N} elements if \code{xval} is an \code{N-by-P} array (such as \code{x}). Note that if \code{gdf} is given as a function, the argument \code{p.initial} is mandatory.
#' @param p.initial Initial model parameters for fitting the DF.
#' @param n.iterations Maximum number of iterations in the repeated fit-and-debias algorithm to evaluate the maximum likelihood.
#' @param correct.lss.bias If \code{TRUE} the \code{distance} values are used to correct for the observational bias due to galaxy clustering (large-scale structure). The overall normalization of the effective folume is chosen such that the expected mass contained in the survey volume is the same as for the uncorrected effective volume.
#' @param n.resampling If \code{n.resampling} is an integer larger than one, the data is resampled \code{n.resampling} times, removing exactly one data points from the set at each iteration. This reampling adds realistic parameter uncertainties with quantiles. If \code{n.resampling = NULL}, no resampling is performed.
#' @param n.jackknife If \code{n.jackknife} is an integer larger than one, the data is jackknife-resampled \code{n.jackknife} times, removing exactly one data point from the observed set at each iteration. This resampling adds model parameters, maximum likelihood estimator (MLE) bias corrected parameter estimates (corrected to order 1/N). If \code{n.jackknife} is larger than the number of data points N, it is automatically reduced to the number of data points.  If \code{n.jackknife = NULL}, no sucm parameters are deterimed.
#' @param write.fit If \code{TRUE}, the best-fitting parameters are displayed in the console.
#' @param xmin,xmax,dx are \code{P}-element vectors (i.e. scalars for 1-dimensional DF) specifying the points (\code{seq(xmin[i],xmax[i],by=dx[i])}) used for some numerical integrations.
#' @param keep.eddington.bias If \code{TRUE}, the data is not corrected for Eddington bias. In this case no fit-and-debias iterations are performed and the argument \code{n.iterations} will be ignored.
#' @param make.posteriors If \code{TRUE}, posterior probability distributions of the observed data are evaluated from the best fitting model.
#' 
#' @details
#' For a detailed description of the method, please refer to the peer-reviewed publication by Obreschkow et al. 2017 (in prep.).
#'
#' @return \code{dffit} returns a structured list, which can be interpreted by other functions, such as \code{\link{dfwrite}}, \code{\link{dfplot}}, \code{\link{dfplotcov}}, \code{\link{dfplotveff}}. The list contains the following sublists:
#' \item{data}{is a list containing the input data, that is the array of observations \code{x}, their Gaussian uncertainties \code{x.err} and the distances if the objects \code{r}.}
#' \item{selection}{is a list describing the selection function underlying the observations, namely the function \code{veff(x)}, which is derived from the input argument \code{selection}.}
#' \item{model}{is a list describing the generative distribution function used to model. The main entry of this list is the function \code{gdf(xval,p)} (often written as phi(x|theta) in the literature).}
#' \item{grid}{is a list of arrays representing a grid in the observables used for numerical integrations. Most importantly, the N-by-P array \code{x} contains the grid points, the N-element vector \code{gdf} gives the corresponding values of the fitted generative distribution function and the N-element vector \code{scd} representing the source count density given no measurement errors.}
#' \item{fit}{is a list describing the fitted generative distribution function. It contains the array \code{p.best} giving the most likely model parameters, as well as their Gaussian uncertainties \code{p.sigma} and covariance matrix \code{p.covariance}. The list also contains the function \code{gdf(x)}, which is the model function evaluated at the parameters \code{p.best}; and the function \code{scd(x)=gdf(x)*veff(x)} representing the expected source count density.}
#' \item{posteriors}{is a list of arrays specifying the posterior PDFs of the observed data, given the best-fitting model. The posterior PDFs are given via their means, standard deviations and modes, as well as a random value drawn from the posterior PDFs. This random value can be used to plot unbiased distribution functions, such as mass functions.}
#' \item{options}{is a list of various optional input arguments of \code{dffit}.}
#'
#' @keywords schechter function
#' @keywords mass function
#' @keywords fit
#'
#' @examples
#' # For a quick overview of some key functionalities of dftools, run the following three examples:
#' # 1) Example of fitting a mass function under strong Eddington bias
#' dfexample1()
#' 
#' # 2) Example of fitting a mass function in the presence of large-scale structure
#' dfexample2()
#' 
#' # 3) Example of fitting a 2D distribution function:
#' dfexample3()
#' 
#' # The following examples introduce the basics of dftools step-by step.
#' # First, generate a mock sample of 1000 galaxies with 0.5dex mass errors
#' dat = dfmockdata(sigma=0.5)
#' 
#' # show the observed and true log-masses (x and x.true) as a function of distance r
#' plot(dat$r,dat$x,col='grey'); points(dat$r,dat$x.true,pch=20)
#' 
#' # fit a Schechter function to the mock sample without accounting for errors
#' survey = dffit(dat$x, dat$veff)
#' 
#' # plot fit and add a black dashed line showing the input MF
#' mfplot(survey, xlim=c(1e6,2e11), ylim=c(2e-4,2), show.data.histogram = TRUE)
#' lines(10^survey$grid$x, pmax(2e-4,survey$model$gdf(survey$grid$x,c(-2,10,-1.3))),lty=2)
#' 
#' # do the same again, while accountting for measurement errors in the fit
#' survey = dffit(dat$x, dat$veff, dat$x.err)
#' mfplot(survey, xlim=c(1e6,2e11), ylim=c(2e-4,2), show.data.histogram = TRUE)
#' lines(10^survey$grid$x, pmax(2e-4,survey$model$gdf(survey$grid$x,c(-2,10,-1.3))),lty=2)
#'
#' # show fitted parameter PDFs and covariances with true input parameters as black points
#' dfplotcov(survey, reference = c(-2,10,-1.3))
#'
#' # show effective volume function
#' dfplotveff(survey)
#'
#' # Now create a smaller survey of only 30 galaxies with 0.5dex mass errors
#' dat = dfmockdata(30, sigma=0.5)
#' 
#' # fit a Schechter function and determine uncertainties by resampling the best fit
#' survey = dffit(dat$x, dat$veff, dat$x.err, n.resampling = 30)
#' 
#' # show best fit with 68% Gaussian uncertainties from Hessian and posterior data
#' mfplot(survey, show.data.histogram = TRUE, uncertainty.type = 1)
#' 
#' # show best fit with 68% and 95% resampling uncertainties and posterior data
#' mfplot(survey, show.data.histogram = TRUE, uncertainty.type = 3)
#' 
#' # add input model as dashed lines
#' lines(10^survey$grid$x, pmax(2e-4,survey$model$gdf(survey$grid$x,c(-2,10,-1.3))),lty=2)
#'
#' @author Danail Obreschkow
#'
#' @export

dffit <- function(x, # normally log-mass, but can be multi-dimensional
                  selection = NULL,
                  x.err = NULL, # Gaussian standard deviations of x
                  r = NULL, # Distance (or redshift)
                  gdf = 'Schechter',
                  p.initial = NULL,
                  n.iterations = 100,
                  correct.lss.bias = FALSE,
                  n.resampling = NULL,
                  n.jackknife = NULL,
                  write.fit = TRUE,
                  xmin = 4,
                  xmax = 12,
                  dx = 0.01,
                  keep.eddington.bias = FALSE,
                  make.posteriors = TRUE) {

  # Set timer
  tStart = Sys.time()
  
  # Initialize main dataframe
  survey = list(data = list(x = x, x.err = x.err, r = r),
                selection = list(),
                model = list(),
                grid = list(xmin = xmin, xmax = xmax, dx = dx),
                fit = list(),
                options = list(p.initial = p.initial, n.iterations = n.iterations,
                               n.resampling = n.resampling, n.jackknife = n.jackknife,
                               keep.eddington.bias = keep.eddington.bias, correct.lss.bias = correct.lss.bias),
                tmp = list(selection = selection, gdf = gdf))
  
  # Check and pre-process input arguments
  survey = .handle.input(survey)
  survey = .make.veff(survey)
  survey = .make.grid(survey)
  survey[which(names(survey)=='tmp')] = NULL
  
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
  if (survey$fit$status$converged) survey = .add.Gaussian.errors(survey)
  
  # Resample to determine more accurate uncertainties with quantiles
  if (survey$fit$status$converged & !is.null(survey$options$n.resampling)) survey = .resample(survey)
  
  # Resample to determine more accurate uncertainties with quantiles
  if (survey$fit$status$converged & !is.null(survey$options$n.jackknife)) survey = .jackknife(survey)
  
  # Make posteriors
  if (survey$fit$status$converged & !is.null(survey$data$x.err) & make.posteriors) survey = .dfposteriors(survey)
  
  # Write best fitting parameters
  if (write.fit) dfwrite(survey)

  # Finalize
  survey$fit$status$walltime = as.double(Sys.time())-as.double(tStart)
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
      if (min(survey$data$x.err)<=0) stop('All values of x.err must be positive.')
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
  
  # Handle gdf
  if (is.function(survey$tmp$gdf)) {
    survey$model$gdf = survey$tmp$gdf
    survey$model$gdf.equation = NA
    if (is.null(survey$options$p.initial)) stop('For user-defined distribution functions initial parameters must be given.')
  } else {
    survey$model$gdf = function(x,p) {dfmodel(x, p, type = survey$tmp$gdf)}
    survey$model$gdf.equation = dfmodel(output = 'equation', type = survey$tmp$gdf)
    if (is.null(survey$options$p.initial)) survey$options$p.initial = dfmodel(output = 'initial', type = survey$tmp$gdf)
  }
  survey$model$n.para = length(survey$options$p.initial)
  
  # Handle n.iterations
  if (is.null(survey$options$n.iterations)) stop('n.iterations must be a positive integer.')
  if (survey$options$n.iterations<1) stop('n.iterations must be a positive integer.')
  
  # Handle n.resampling
  if (!is.null(survey$options$n.resampling)) {
    if (survey$options$n.resampling<2) stop('n.resampling must be 2 or larger.')
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
      } else if (n.dim==2) {
        vapprox = function(xval) {
          z = 1/(akima::interp(x[,1],x[,2],1/veff.values,xval[1],xval[2],duplicate='mean'))$z
          if (is.na(z)) {return(0)} else {return(z)}
        }
      } else {
        stop('Linear interpolation of Veff not implemented for DF with more than 2 dimensions. Use a different selection type.')
      }
      veff.function.elemental = function(xval) {
        if (any(xval<xmin) | any(xval>xmax)) {
          return(0)
        } else {
          return(vapprox(xval))
        }
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
  
  survey$selection = list(veff = veff.function, veff.no.lss = veff.function,
                          veff.input.values = veff.values,
                          veff.input.function = veff.userfct,
                          f = f.function,
                          rmin = rmin, rmax = rmax,
                          mode = mode)
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
  
  # renormalize veff
  int.ref = function(x) survey$selection$veff.no.lss(x)*survey$model$gdf(x,p)*10^x
  int.exp = function(x) veff.lss.scale(x)*survey$model$gdf(x,p)*10^x
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

#' @export
.corefit <- function(survey, supress.warning = FALSE) {
  
  # simplify input
  x = survey$data$x
  x.err = survey$data$x.err
  gdf = survey$model$gdf
  p.initial = survey$options$p.initial
  x.mesh = survey$grid$x
  x.mesh.dv = survey$grid$dvolume
  n.iterations = survey$options$n.iterations
  keep.eddington.bias = survey$options$keep.eddington.bias
  
  # get array sizes
  n.data = dim(x)[1]
  n.dim = dim(x)[2]
  n.mesh = dim(x.mesh)[1]
  
  # Input handling
  if (!survey$options$correct.lss.bias) veff.mesh = c(survey$grid$veff)
  if (is.null(x.err) & !survey$options$correct.lss.bias) n.iterations = 1
  
  if (!is.null(x.err)) {
    
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
    
    # make a priori observed PDFs
    rho.observed = list()
    for (i in seq(n.data)) {
      d = x[i,]-t(x.mesh)
      rho.observed[[i]] = exp(-colSums(d*(invC[i,,]%*%d))/2)
      if (sum(rho.observed[[i]])<0.01) {
        index = which.min(abs(colSums(d)))[1]
        rho.observed[[i]][index] = 1
      }
    }
  }
  
  # Iterative algorithm
  running = TRUE
  k = 0
  offset = 0
  chain = array(NA,c(n.iterations,length(p.initial)+1))
  
  # definie debiasing function
  debias = function(p) {
    rho.unbiased = array(0,n.mesh)
    if (is.null(x.err)) {
      for (i in seq(n.data)) {
        difference = colSums(abs(t(x.mesh)-x[i,]))
        index = which.min(difference)[1]
        rho.unbiased[index] = rho.unbiased[index]+1/x.mesh.dv
      }
    } else {
      if (keep.eddington.bias) {
        prior = array(1,n.mesh)
      } else {
        # predicted source counts (up to a factor x.mesh.dv)
        prior = gdf(x.mesh,p)*veff.mesh
        prior[!is.finite(prior)] = 0
        prior = pmax(0,prior)
      }
      for (i in seq(n.data)) {
        rho.corrected = rho.observed[[i]]*prior
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
      phi = gdf(x.mesh,p)
      # safety operations (adding about 50% computation time)
      phi[!is.finite(phi)] = 0
      phi = pmax(.Machine$double.xmin,phi) # also vectorizes the array
      # end safety operations
      return(sum(phi*veff.mesh-log(phi)*rho.unbiased)*x.mesh.dv-offset)
      #lambda = pmax(.Machine$double.xmin,phi*veff.mesh)
      #return(sum(lambda-log(lambda)*rho.unbiased)*x.mesh.dv-offset)
    }
    
    # maximize ln(L)
    opt = optim(p.initial,neglogL,hessian=TRUE)#,control=list(reltol=1e-20,abstol=1e-20,maxit=1e4))
    offset = opt$value+offset
    chain[k,] = c(opt$par,opt$value)
    
    # assess convergence
    if ((is.null(x.err) & !survey$options$correct.lss.bias)| keep.eddington.bias) { # exist without extra iterations
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
    ln.evidence = 0.5*log(det(cov))-offset
  }
  
  # finalize output
  fit = list(p.best = opt$par, p.covariance = cov, lnL = function(p) -neglogL(p),
             status = list(n.iterations = k, converged = converged, chain = chain[1:k,]), ln.evidence = ln.evidence)
  return(fit)
  
}

.add.Gaussian.errors <- function(survey) {

  cov = survey$fit$p.covariance

  # make Gaussian uncertainties of DF
  eig = eigen(cov)
  np = survey$model$n.para
  nx = survey$grid$n.points
  index = 0
  nsteps.tot.max = 64
  nsteps = max(2,round(nsteps.tot.max^(1/np))) # a larger number of steps leads to a more accurate sampling of the covariance ellipsoid
  y.new = array(NA,c(nx,nsteps^np))
  k = atan(seq(-pi/2,pi/2,length=nsteps))
  kv = array(NA,c(nsteps^np,np))
  for (i in seq(np)) {
    kv[,i] = rep(k,nsteps^(i-1),each=nsteps^(np-i))
  }
  
  # sample surface of covariance ellipsoid
  for (i in seq(nsteps^np)) {
    e = kv[i,]#/max(1e-20,sqrt(sum(kv[i,]^2)))
    v = c(eig$vectors%*%(sqrt(eig$values)*e))
    p.new = survey$fit$p.best+v
    y.new[,i] = survey$model$gdf(survey$grid$x,p.new)
    if (any(!is.finite(y.new[,i])) | any(y.new[,i]<0)) y.new[,i] = NA
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
  if (n.data<2) stop('Jackknifing requires at least two objects.')
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
    if (is.null(survey$data$x.err)) {
      jn$data$x.err = NULL
    } else if (length(dim(survey$data$x.err))==2) {
      jn$data$x.err = as.matrix(survey$data$x.err[list,])
    } else if (length(dim(survey$data$x.err))==3) {
      jn$data$x.err = as.matrix(survey$data$x.err[list,,])
    }
    p.new[i,] = .corefit(jn, supress.warning = TRUE)$p.best
  }
  cat('\r')
  
  # estimate covariance
  ok = !is.na(rowSums(p.new))
  cov.jn = cov(p.new[ok,])*(n.data-1)
  
  # compute poisson covariance
  jn = survey
  jn$options$p.initial = survey$fit$p.best
  jn$options$n.iterations = 1
  q = c(0.16,0.5,0.84)
  p.pois = array(NA,c(3,np))
  for (i in seq(3)) {
    n.new = qpois(q[i],n.data) 
    jn$grid$veff = survey$grid$veff*n.new/n.data
    p.pois[i,] = .corefit(jn, supress.warning = TRUE)$p.best
  }
  cov.pois = cov(p.pois)
  
  # estimate combined covariance
  if (is.na(cov.pois[1,1])) {
    survey$fit$p.covariance.jackknife = cov.jn
  } else {
    survey$fit$p.covariance.jackknife = cov.jn+cov.pois
  }
  
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
  if (n.data<2) stop('Jackknifing requires at least two objects.')
  if (survey$options$n.resampling<2) stop('n.resampling must be at least 2.')
  
  # set up resample survey
  b = survey
  b$options$p.initial = survey$fit$p.best
  x.err.mean = mean(b$data$x.err)
  
  # randomly resample and refit the DF
  cum = cumsum(survey$grid$scd/sum(survey$grid$scd))
  p.new = array(NA,c(survey$options$n.resampling,np))
  for (iteration in seq(survey$options$n.resampling)) {
    cat(sprintf('\rResampling: %4.2f%%',100*iteration/survey$options$n.resampling))
    n.new = max(2,rpois(1,n.data))
    x.obs = array(NA,n.new)
    r = runif(n.new)
    for (i in seq(n.new)) {
      index = which.min(abs(cum-r[i]))
      x.obs[i] = survey$grid$x[index]
    }
    x.obs = x.obs+rnorm(n.new)*x.err.mean
    b$data$n.data = n.new
    b$data$x = cbind(x.obs)
    b$data$x.err = cbind(rep(x.err.mean,n.new))
    p.new[iteration,] = .corefit(b, supress.warning = TRUE)$p.best
  }
  cat('\r')
  
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
  s = array(NA,c(survey$options$n.resampling,survey$grid$n.points))
  for (i in seq(survey$options$n.resampling)) {
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