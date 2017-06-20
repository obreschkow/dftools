#' Fit a generative distribution function, such as a galaxy mass function
#'
#' This function fits galaxy mass function (MF) to a discrete set of \code{N} galaxies with noisy data. More generally, \code{dffit} finds the most likely \code{P}-dimensional distribution function (DF) generating \code{N} objects \code{i=1,...,N} with uncertain measurements \code{P} observables. For instance, if the objects are galaxies, it can fit a MF (\code{P=1}), a mass-size distribution (\code{P=2}) or the mass-spin-morphology distribution (\code{P=3}). A full description of the algorithm can be found in Obreschkow et al. (2017).
#'
#' @importFrom akima interp
#'
#' @param x Normally \code{x} is a \code{N}-element vector, representing the log-masses (log10(M/Msun)) of \code{N} galaxies. More generally, \code{x} can be either a vector of \code{N} elements or a matrix of \code{N-by-P} elements, containing the values of one or \code{P} observables of \code{N} objects, respectively.
#' @param selection Specifies the effective volume \code{Veff(xval)} in which a galaxy of log-mass \code{xval} can be observed; or, more generally, the volume in which an object of observed values \code{xval[1:P]} can be observed. This volume can be specified in five ways: (1) If \code{selection} is a single positive number, it will be interpreted as a constant volume, \code{Veff(xval)=selection}, in which all galaxies are fully observable. \code{Veff(xval)=0} is assumed outside the "observed domain". This domain is defined as \code{min(x)<=xval<=max(x)} for one observable (\code{P=1}), or as \code{min(x[,j])<=xval[j]<=max(x[,j])} for all \code{j=1,...,P} if \code{P>1}. This mode can be used for volume-complete surveys or for simulated galaxies in a box. (2) If \code{selection} is a vector of \code{N} elements, they will be interpreted as the effective volumes for each of the \code{N} galaxies. \code{Veff(xval)} is interpolated (linearly in \code{1/Veff}) for other values \code{xval}. \code{Veff(xval)=0} is assumed outside the observed domain. (3) \code{selection} can be a function of \code{P} variables, which directly specifies the effective volume for any \code{xval}, i.e. \code{Veff(xval)=selection(xval)}. (4) \code{selection} can also be a list (\code{selection = list(veff.values, veff.userfct)}) of an \code{N}-element vector \code{Veff.values} and a \code{P}-dimensional function \code{veff.userfct}. In this case, the effective volume is computed using a hybrid scheme of modes (2) and (3): \code{Veff(xval)} will be interpolated from the \code{N} values of \code{Veff.values} inside the observed domain, but set equal to \code{veff.userfct} outside this domain. (5) Finally, \code{selection} can be a list of two functions and one optional 2-element vector: \code{selection = list(f, dVdr, rrange)}, where \code{f = function(xval,r)} is the isotropic selection function and \code{dVdr = function(r)} is the derivative of the total survey volume as a function of comoving distance \code{r}. The optional vector \code{rrange} (with default \code{rrange=c(0,Inf)}) gives the minimum and maximum comoving distance limits of the survey. Outside these limits \code{Veff=0} will be assumed.
#' @param x.err Optional vector or array specifying the observational errors of \code{x}. If \code{x} is a vector then \code{x.err} must also be a vector of same length. Its elements are interpreted as the standard deviations of Gaussian uncertainties in \code{x}. If \code{x} is a \code{N-by-P} matrix representing \code{N} objects with \code{P} observables, then \code{x.err} must be either a \code{N-by-P} matrix or a \code{N-by-P-by-P} array. In the first case, the elements \code{x.err[i,]} are interpreted as the standard deviations of Gaussian uncertainties on \code{x[i,]}. In the second case, the \code{P-by-P} matrices \code{x.err[i,,]} are interpreted as the covariance matrices of the \code{P} observed values \code{x[i,]}.
#' @param distance Optional vector of \code{N} elements specifying the comoving distances of the \code{N} galaxies. This vector is only needed if \code{correct.lss.bias = TRUE}.
#' @param gdf Either a string or a function specifying the DF to be fitted. A string is interpreted as the name of a predefined mass function (i.e. functions of one obervable, \code{P=1}). Available options are \code{'Schechter'} for Schechter function (3 parameters), \code{'PL'} for a power law (2 parameters), or \code{'MRP'} for an MRP function (4 parameters). Alternatively, \code{gdf = function(xval,p)} can be any function of the \code{P} observable(s) \code{xval} and a list of parameters \code{p}. IMPORTANT: The function \code{gdf(xval,p)} must be fully vectorized in \code{xval}, i.e. it must output a vector of \code{N} elements if \code{xval} is an \code{N-by-P} array (such as \code{x}). Note that if \code{gdf} is given as a function, the argument \code{p.initial} is mandatory.
#' @param p.initial Initial model parameters for fitting the DF.
#' @param n.iterations Maximum number of iterations in the repeated fit-and-debias algorithm to evaluate the maximum likelihood.
#' @param n.resampling Integer (>0) specifying the number of iterations for the resampling of the most likely DF used to evaluate realistic parameter uncertainties with quantiles. If \code{n.resampling = NULL}, no resampling is performed.
#' @param correct.mle.bias The maximum likelihood estimator (MLE) of a finite dataset can be biased – a general property of the ML approach. If \code{TRUE}, \code{dffit} also outputs the parameters, where this estimator bias has been corrected, to first order in 1/N, using jackknifing.
#' @param correct.lss.bias If \code{TRUE} the \code{distance} values are used to correct for the observational bias due to galaxy clustering (large-scale structure).
#' @param lss.normalization Integer value, determining the type of normalization used when accounting for cosmic large scale structure. Use \code{lss.normalization=1} to conserve the mass (sum of 10^x) in the survey; \code{lss.normalization=2} to conserve the sum of exp(x) in the survey; \code{lss.normalization=3} to conserve the number of objects (galaxies) in the survey.
#' @param write.fit If \code{TRUE}, the best-fitting parameters are displayed in the console.
#' @param xmin,xmax,dx are \code{P}-element vectors (i.e. scalars for 1-dimensional DF) specifying the points (\code{seq(xmin[i],xmax[i],by=dx[i])}) used for some numerical integrations.
#' @param keep.eddington.bias If \code{TRUE}, the data is not corrected for Eddington bias. In this case no fit-and-debias iterations are performed and the argument \code{n.iterations} will be ignored.
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
#' \item{options}{is a list of various optional input arguments of \code{dffit}.}
#'
#' @keywords schechter function
#' @keywords mass function
#' @keywords fit
#'
#' @examples
#' # full example in 1D (subject to non-trivial Eddington bias)
#' dfexample1()
#' 
#' # full example in 2D:
#' dfexample2()
#' 
#' # basic example of fitting data without including measurement errors
#' dat = dfmockdata(100)
#' survey = dffit(dat$x, dat$veff)
#' mfplot(survey, xlim=c(1e6,2e11), ylim=c(2e-4,2), show.data.histogram = TRUE)
#' 
#' # add a black dashed line showing input MF used to generate the data
#' lines(10^survey$grid$x, survey$model$gdf(survey$grid$x,c(-2,10,-1.3)),lty=2)
#' 
#' # now, include measurement errors in fit
#' survey = dffit(dat$x, dat$veff, dat$x.err)
#' mfplot(survey, xlim=c(1e6,2e11), ylim=c(2e-4,2), show.data.histogram = TRUE)
#' lines(10^survey$grid$x, survey$model$gdf(survey$grid$x,c(-2,10,-1.3)),lty=2)
#'
#' # show fitted parameter PDFs and covariances with true input parameters as black points
#' dfplotcov(survey, p = c(-2,10,-1.3))
#'
#' # show fitted effective volume function
#' dfplotveff(survey)
#'
#' # determine Schechter function uncertainties from resampling and
#' # evaluate bias-corrected MLE (plotted in red)
#' dat = dfmockdata(30)
#' survey = dffit(dat$x, dat$veff, dat$x.err, n.resampling = 1e2, correct.mle.bias = TRUE)
#' mfplot(survey, uncertainty.type=3, nbins=10, bin.xmin=6.5, bin.xmax=9.5, xlim=c(1e6,2e11), ylim=c(2e-4,2))
#' 
#' # add a red dashed line with the bias corrected ML estimator
#' lines(10^survey$grid$x,survey$grid$gdf.mle.bias.corrected,col='red',lty=2)
#' 
#' # add a black dashed line showing the input model
#' lines(10^survey$grid$x, survey$model$gdf(survey$grid$x,c(-2,10,-1.3)),lty=2)
#'
#' # evaluate posteriors of the mass measurements and
#' # visualize the change in the mass mode between observation and posterior
#' # i.e. the Eddington bias correction
#' survey = dfposteriors(survey)
#' plot(survey$dat$x,survey$posterior$x.mode)
#' lines(survey$dat$x,survey$dat$x)
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
                  n.resampling = NULL,
                  correct.mle.bias = FALSE,
                  correct.lss.bias = FALSE,
                  lss.normalization = 1,
                  write.fit = TRUE,
                  xmin = 4,
                  xmax = 12,
                  dx = 0.01,
                  keep.eddington.bias = FALSE) {

  # Set timer
  tStart = Sys.time()
  
  # Initialize main dataframe
  survey = list(data = list(x = x, x.err = x.err, r = r),
                selection = list(),
                model = list(),
                grid = list(xmin = xmin, xmax = xmax, dx = dx),
                fit = list(),
                options = list(p.initial = p.initial, n.iterations = n.iterations, n.resampling = n.resampling,
                               keep.eddington.bias = keep.eddington.bias,  correct.lss.bias = correct.lss.bias,
                               lss.normalization = lss.normalization, correct.mle.bias = correct.mle.bias),
                tmp = list(selection = selection, gdf = gdf))
  
  # Check and pre-process input arguments
  survey = .handle.input(survey)
  survey = .make.veff(survey)
  survey = .make.grid(survey)
  survey[which(names(survey)=='tmp')] = NULL
  
  # Find most likely generative model
  survey = .corefit(survey)
  if (write.fit) dfwrite(survey)
  
  # Determine Gaussian uncertainties
  survey = .add.Gaussian.errors(survey)
  
  # Resample to determine more accurate uncertainties with quantiles
  if (survey$fit$status$converged & !is.null(survey$options$n.resampling)) {survey = .resample(survey)}
  
  # Estimator bias correction
  if (survey$fit$status$converged & survey$options$correct.mle.bias) {survey = .correct.mle.bias(survey)}

  # Finalize
  survey$fit$status$walltime.total = as.double(Sys.time())-as.double(tStart)
  invisible(survey)

}

.correct.lss.bias = function(survey) {
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
  
  # Handle x.err
  if (!is.null(survey$data$x.err)) {
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
    if (min(survey$data$x)<=0) stop('All values of x.err must be positive.')
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
    if (is.null(survey$options$lss.normalization)) stop('lss.normalization must be 1, 2 or 3.')
    if (survey$options$lss.normalization<1 | survey$options$lss.normalization>3) stop('lss.normalization must be 1, 2 or 3.')
    if (survey$data$n.dim>1) stop('Currently correct.lss.bias can only be TRUE for one-dimensional DFs.')
  }
  
  # Handle gdf
  if (is.function(survey$tmp$gdf)) {
    survey$model$gdf = survey$tmp$gdf
    survey$model$gdf.equation = NA
    if (is.null(survey$options$p.initial)) stop('For user-defined distribution functions initial parameters must be given.')
  } else {
    survey$model$gdf <- function(x,p) {dfmodel(x, p, type = survey$tmp$gdf)}
    survey$model$gdf.equation = dfmodel(output = 'equation', type = survey$tmp$gdf)
    if (is.null(survey$options$p.initial)) survey$options$p.initial = dfmodel(output = 'initial', type = survey$tmp$gdf)
  }
  survey$model$n.para = length(survey$options$p.initial)
  
  # Handle n.iterations
  if (is.null(survey$options$n.iterations)) stop('n.iterations must be a positive integer.')
  if (survey$options$n.iterations<1) stop('n.iterations must be a positive integer.')
  
  # Handle n.resampling
  if (!is.null(survey$options$n.resampling)) {
    if (survey$options$n.resampling<10) stop('n.resampling must be 10 or larger.')
  }
  
  # Handle correct.mle.bias
  if (survey$options$correct.mle.bias) {
    if (survey$data$n.data<2) stop('At least two data points must be given if correct.mle.bias = TRUE.')
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
  
  # Mode 5: Selection given as {f(x,r), dVdr(r), range}
  if (is.list(s)) {
    if (length(s)>=2 & length(s)<=3) {
      if (is.function(s[[1]]) & is.function(s[[2]])) {
        if (length(s)==3) {
          if (is.double(s[[3]])) {
            if (length(s[[3]])==2) {
              rmin = s[[3]][1]
              rmax = s[[3]][2]
            } else {
              stop('Unknown selection format.')
            }
          } else {
            stop('Unknown selection format.')
          }
        } else {
          rmin = 0
          rmax = Inf
        }
        mode = 5
        test = try(s[[1]](NA,NA)*s[[2]](NA),silent=TRUE)
        if (!is.double(test)) stop('In the argument selection = list(f, dVdr, ...), the functions f(xval,r) and dVdr(r) must work if r is a vector.')
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
  
  survey$selection = list(veff = veff.function,
                          veff.input.values = veff.values, veff.input.function = veff.userfct,
                          rmin = rmin, rmax = rmax,
                          mode = mode)
  invisible(survey)
}

.make.veff.lss = function(survey,p) {
  
  s = survey$tmp$selection
  x = survey$data$x
  r = survey$data$r
  n.dim = survey$data$n.dim
  if (n.dim>1) stop('Currently LSS can only be corrected in one-dimensional distribution functions.')
  n.data = survey$data$n.data
  rmin = survey$selection$rmin
  rmax = survey$selection$rmax
  
  # decided not to implement this
  # if (survey$selection$mode==5) {
  #   integrand.lss = function(x,r,p) {s[[1]](x,r)*survey$model$gdf(x,p)}
  #   veff.lss.function.elemental = function(xval) {
  #     sm = 0
  #     i = 1
  #     for (i in seq(n.data)) {
  #       if (r[i]>=rmin & r[i]<=rmax & s[[1]](xval,r[i])>0) {
  #         sm = sm+s[[1]](xval,r[i])/integrate(integrand.lss,survey$grid$xmin,survey$grid$xmax,r[i],p,stop.on.error=FALSE)$value
  #       }
  #     }
  #     return(sm)
  #   }
  # }
  
  # determine minimum log-mass, which would be detected at distance r[i] for each i
  if (is.null(survey$selection$xlim)) {
    survey$selection$xmin = array(NA,n.data)
    for (i in seq(n.data)) {
      survey$selection$xmin[i] = min(x[r>=r[i]])
    }
  }
  #survey$selection$xmin = 2*log10(survey$data$r)+6
  
  # evaluate integrals
  integral = array(NA,n.data)
  for (i in seq(n.data)) {
    integral[i] = integrate(survey$model$gdf,survey$selection$xmin[i],survey$grid$xmax,p,stop.on.error=FALSE)$value
  }
  
  # make veff.lss function
  veff.lss.function.elemental = function(xval) sum(1/integral[xval>=survey$selection$xmin])
  veff.lss.scale = Vectorize(veff.lss.function.elemental)
  
  # renormalize veff
  #veff.lss.scale = veff.lss.theory
  if (survey$options$lss.normalization == 1) {
    integrand = function(x) veff.lss.scale(x)*survey$model$gdf(x,p)*10^x
    expectation = integrate(integrand,survey$grid$xmin,survey$grid$xmax,stop.on.error=FALSE)$value
    normalization.factor = sum(10^x)/expectation
  } else if (survey$options$lss.normalization == 2) {
    integrand = function(x) veff.lss.scale(x)*survey$model$gdf(x,p)*exp(x)
    expectation = integrate(integrand,survey$grid$xmin,survey$grid$xmax,stop.on.error=FALSE)$value
    normalization.factor = sum(exp(x))/expectation
  } else if (survey$options$lss.normalization == 3) {
    integrand = function(x) veff.lss.scale(x)*survey$model$gdf(x,p)
    expectation = integrate(integrand,survey$grid$xmin,survey$grid$xmax,stop.on.error=FALSE)$value
    normalization.factor = n.data/expectation
  } else {
    stop('unknown LSS normalization')
  }
  veff.lss = function(x) veff.lss.scale(x)*normalization.factor
  
  # return output
  survey$selection$veff = veff.lss
  invisible(survey)
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
  
  # Start timer
  tStart = Sys.time()
  
  # simplify input
  x = survey$data$x
  x.err = survey$data$x.err
  gdf = survey$model$gdf
  p.initial = survey$options$p.initial
  x.mesh = survey$grid$x
  x.mesh.dv = survey$grid$dvolume
  n.iterations = survey$options$n.iterations
  keep.eddington.bias = survey$options$keep.eddington.bias
  n.data = survey$data$n.data
  n.dim = survey$data$n.dim

  # Input handling
  n.mesh = dim(x.mesh)[1]
  if (!survey$options$correct.lss.bias) veff.mesh = c(survey$grid$veff)
  if (is.null(x.err)) n.iterations = 1

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
    }
  }
 
  # Iterative algorithm
  running = TRUE
  k = 0
  chain = array(NA,c(n.iterations,length(p.initial)+1))
  
  while (running) {

    k = k+1
    
    # determine veff LSS
    if (survey$options$correct.lss.bias) {
      survey = .make.veff.lss(survey,p.initial)
      veff.mesh = c(survey$selection$veff(survey$grid$x))
    }

    # make unbiased source density function
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
        prior = gdf(x.mesh,p.initial)*veff.mesh
        prior[!is.finite(prior)] = 0
        prior = pmax(0,prior)
      }
      for (i in seq(n.data)) {
        rho.corrected = rho.observed[[i]]*prior
        rho.corrected = rho.corrected/(sum(rho.corrected)*x.mesh.dv)
        rho.unbiased = rho.unbiased+rho.corrected
      }
    }

    # make -ln(L)
    neglogL = function(p) {
      phi = gdf(x.mesh,p)
      # safety operations (adding about 50% computation time)
      phi[!is.finite(phi)] = 0
      phi = pmax(.Machine$double.xmin,phi) # also vectorizes the array
      # end safety operations
      return(sum(phi*veff.mesh-log(phi)*rho.unbiased)*x.mesh.dv)
    }

    # maximize ln(L)
    opt = optim(p.initial,neglogL,hessian=TRUE,control=list(parscale=rep(1,length(p.initial)),reltol=1e-10,abstol=1e-10,maxit=1e4))
    chain[k,] = c(opt$par,opt$value)
    
    # assess convergence
    if (is.null(x.err) | keep.eddington.bias) { # exist without extra iterations
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
  
  # Timer
  if (is.null(survey$fit$status$walltime.fitting)) {survey$fit$status$walltime.fitting = 0}
  dt = as.double(Sys.time())-as.double(tStart)
  survey$fit$status$walltime.fitting = survey$fit$status$walltime.fitting+dt

  # make output
  if (det(opt$hessian)<1e-12) {
    converged = FALSE
    cov = NULL
    if (!supress.warning) {
      cat('WARNING: Fit ill-conditionned. Consider providing better initial parameters or selection arguments.\n')
    }
  } else {
    cov = solve(opt$hessian)
  }
  
  # finalize output
  if (survey$options$correct.lss.bias) {survey = .make.veff.lss(survey,p.initial)}
  survey$fit$p.best = opt$par
  survey$fit$p.sigma = sqrt(diag(cov))
  survey$fit$p.covariance = cov
  survey$fit$gdf = function(x) survey$model$gdf(x,opt$par)
  survey$fit$scd = function(x) survey$fit$gdf(x)*survey$selection$veff(x)
  survey$fit$logL = function(p) -neglogL(p)
  survey$fit$status = list(n.fit.and.debias = k, converged = converged, chain = chain[1:k,])
  survey$grid$gdf = c(survey$fit$gdf(survey$grid$x))
  survey$grid$scd = c(survey$fit$scd(survey$grid$x))
  
  invisible(survey)

}

.add.Gaussian.errors <- function(survey) {

  cov = survey$fit$p.covariance

  # make Gaussian uncertainties of DF
  eig = eigen(cov)
  np = survey$model$n.para
  nx = survey$grid$n.points
  index = 0
  nsteps.tot.max = 216
  nsteps = max(2,round(nsteps.tot.max^(1/np))) # a larger number of steps leads to a more accurate sampling of the covariance ellipsoid
  y.new = array(NA,c(nx,nsteps^np))
  k = seq(-1,1,length=nsteps)
  step = array(1,np)
  step[1] = 0
  
  # sample surface of covariance ellipsoid
  for (i in seq(nsteps^np)) {

    # next step
    overflow = TRUE
    j = 1
    while (overflow) {
      step[j] = step[j]+1
      if (step[j]>nsteps) {
        step[j] = 1
        j = j+1
      } else {
        overflow = FALSE
      }
    }

    # evaluate DF
    if (any(k[step]!=0)) {
      e = k[step]/sqrt(sum(k[step]^2))
      v = c(eig$vectors%*%(sqrt(eig$values)*e))
      p.new = survey$fit$p.best+v
      y.new[,i] = survey$model$gdf(survey$grid$x,p.new)
      if (any(!is.finite(y.new[,i])) | any(y.new[,i]<0)) y.new[,i] = NA
    }
  }

  survey$grid$gdf.error.neg = array(NA,nx)
  survey$grid$gdf.error.pos = array(NA,nx)
  for (i in seq(nx)) {
    survey$grid$gdf.error.neg[i] = survey$grid$gdf[i]-min(c(y.new[i,],Inf),na.rm=TRUE)
    survey$grid$gdf.error.pos[i] = max(c(y.new[i,],-Inf),na.rm=TRUE)-survey$grid$gdf[i]
  }

  invisible(survey)
}

.correct.mle.bias <- function(survey) {
  n = dim(survey$data$x)[1]
  np = survey$model$n.para
  if (n<2) {
    stop('Bias correction requires at least two objects.')
  } else if (n>=1e3) {
    cat('WARNING: bias correction normally not relevant for more than 1000 objects.\n')
  }
  p.new = array(NA,c(np,n))
  b = survey
  b$options$p.initial = survey$fit$p.best
  b$options$n.iterations = 1
  b$grid$veff = survey$grid$veff*(n-1)/n
  b$data$n.data = n-1
  for (i in seq(n)) {
    list = setdiff(seq(n),i)
    if (is.null(survey$data$x.err)) {
      x.err = NULL
    } else if (length(dim(survey$data$x.err))==2) {
      x.err = as.matrix(survey$data$x.err[list,])
    } else if (length(dim(survey$data$x.err))==3) {
      x.err = as.matrix(survey$data$x.err[list,,])
    }
    b$data$x = as.matrix(survey$data$x[list,])
    b$data$x.err = x.err
    b = .corefit(b, supress.warning = TRUE)
    p.new[,i] = b$fit$p.best
  }
  p.reduced = apply(p.new, 1, mean, na.rm = T)
  
  survey$fit$p.best.mle.bias.corrected = n*survey$fit$p.best-(n-1)*p.reduced
  survey$fit$gdf.mle.bias.corrected = function(x) survey$model$gdf(x,survey$fit$p.best.mle.bias.corrected)
  survey$fit$scd.mle.bias.corrected = function(x) survey$fit$gdf.mle.bias.corrected(x)*survey$selection$veff(x)
  survey$grid$gdf.mle.bias.corrected = survey$fit$gdf.mle.bias.corrected(survey$grid$x)
  survey$grid$scd.mle.bias.corrected = survey$fit$scd.mle.bias.corrected(survey$grid$x)
  
  invisible(survey)
}

.resample <- function(survey, seed = 1) {

  # randomly resample and refit the DF
  if (dim(survey$data$x)[2]>1) stop('resampling only available for one-dimensional distribution functions.')
  set.seed(seed)
  np = survey$model$n.para
  #x = seq(survey$grid$xmin,survey$grid$xmax,survey$grid$dx)
  density = survey$grid$scd
  cum = cumsum(density/sum(density))
  n.data = dim(survey$data$x)[1]
  b = survey
  b$options$n.iterations = 1
  b$data$x.err = NULL
  p.new = array(NA,c(survey$options$n.resampling,np))
  for (iteration in seq(survey$options$n.resampling)) {
    n.new = max(2,rpois(1,n.data))
    x.obs = array(NA,n.new)
    r = runif(n.new)
    for (i in seq(n.new)) {
      index = which.min(abs(cum-r[i]))
      x.obs[i] = survey$grid$x[index]
    }
    b$data$n.data = n.new
    b$data$x = cbind(x.obs)
    p.new[iteration,] = .corefit(b, supress.warning = TRUE)$fit$p.best
  }

  # make parameter quantiles
  q = c(0.02,0.16,0.84,0.98)
  p.quant = array(NA,c(length(q),np))
  for (i in seq(np)) {
    p.quant[,i] = quantile(p.new[,i],q,names=FALSE)
  }
  survey$fit$p.quantile.02 = p.quant[1,]
  survey$fit$p.quantile.16 = p.quant[2,]
  survey$fit$p.quantile.84 = p.quant[3,]
  survey$fit$p.quantile.98 = p.quant[4,]

  # make DF quantiles
  s = array(NA,c(survey$options$n.resampling,survey$grid$n.points))
  for (iteration in seq(survey$options$n.resampling)) {
    s[iteration,] = survey$model$gdf(survey$grid$x,p.new[iteration,])
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
