#' Fit a generative distribution function, such as a galaxy mass function
#'
#' This function fits galaxy mass function (MF) to a discrete set of \code{N} galaxies with noisy data. More generally, \code{dffit} finds the most likely \code{P}-dimensional distribution function (DF) generating \code{N} objects \code{i=1,...,N} with uncertain measurements \code{P} observables. For instance, if the objects are galaxies, it can fit a MF (\code{P=1}), a mass-size distribution (\code{P=2}) or the mass-spin-morphology distribution (\code{P=3}). A full description of the algorithm can be found in Obreschkow et al. (2017).
#'
#' @importFrom akima interpp
#'
#' @param x Normally \code{x} is a \code{N}-element vector, representing the log-masses (log10(M/Msun)) of \code{N} galaxies. More generally, \code{x} can be either a vector of \code{N} elements or a matrix of \code{N-by-P} elements, containing the values of one or \code{P} observables of \code{N} objects, respectively.
#' @param selection Specifies the effective volume \code{Veff(xval)} in which a galaxy of log-mass \code{xval} can be observed; or, more generally, the volume in which an object of observed values \code{xval[1:P]} can be observed. This volume can be specified in five ways: (1) If \code{selection} is a single positive number, it will be interpreted as a constant volume, \code{Veff(xval)=selection}, in which all galaxies are fully observable. \code{Veff(xval)=0} is assumed outside the "observed domain". This domain is defined as \code{min(x)<=xval<=max(x)} for one observable (\code{P=1}), or as \code{min(x[,j])<=xval[j]<=max(x[,j])} for all \code{j=1,...,P} if \code{P>1}. This mode can be used for volume-complete surveys or for simulated galaxies in a box. (2) If \code{selection} is a vector of \code{N} elements, they will be interpreted as the effective volumes for each of the \code{N} galaxies. \code{Veff(xval)} is interpolated (linearly in \code{1/Veff}) for other values \code{xval}. \code{Veff(xval)=0} is assumed outside the observed domain. (3) \code{selection} can be a function of \code{P} variables, which directly specifies the effective volume for any \code{xval}, i.e. \code{Veff(xval)=selection(xval)}. (4) \code{selection} can also be a list (\code{selection = list(Veff.values, Veff.fct)}) of an \code{N}-element vector \code{Veff.values} and a \code{P}-dimensional function \code{Veff.fct}. In this case, the effective volume is computed using a hybrid scheme of modes (2) and (3): \code{Veff(xval)} will be interpolated from the \code{N} values of \code{Veff.values} inside the observed domain, but set equal to \code{Veff.fct} outside this domain. (5) Finally, \code{selection} can be a list of two functions and one optional 2-element vector: \code{selection = list(f, dVdr, rrange)}, where \code{f = function(xval,r)} is the isotropic selection function and \code{dVdr = function(r)} is the derivative of the total survey volume as a function of comoving distance \code{r}. The optional vector \code{rrange} (with default \code{rrange=c(0,Inf)}) gives the minimum and maximum comoving distance limits of the survey. Outside these limits \code{Veff=0} will be assumed.
#' @param x.err Optional vector or array specifying the observational errors of \code{x}. If \code{x} is a vector then \code{x.err} must also be a vector of same length. Its elements are interpreted as the standard deviations of Gaussian uncertainties in \code{x}. If \code{x} is a \code{N-by-P} matrix representing \code{N} objects with \code{P} observables, then \code{x.err} must be either a \code{N-by-P} matrix or a \code{N-by-P-by-P} array. In the first case, the elements \code{x.err[i,]} are interpreted as the standard deviations of Gaussian uncertainties on \code{x[i,]}. In the second case, the \code{P-by-P} matrices \code{x.err[i,,]} are interpreted as the covariance matrices of the \code{P} observed values \code{x[i,]}.
#' @param distance Optional vector of \code{N} elements specifying the comoving distances of the \code{N} galaxies. This vector is only needed if \code{correct.lss.bias = TRUE}.
#' @param phi Either a string or a function specifying the DF to be fitted. A string is interpreted as the name of a predefined mass function (i.e. functions of one obervable, \code{P=1}). Available options are \code{'Schechter'} for Schechter function (3 parameters), \code{'PL'} for a power law (2 parameters), or \code{'MRP'} for an MRP function (4 parameters). Alternatively, \code{phi = function(xval,p)} can be any function of the \code{P} observable(s) \code{xval} and a list of parameters \code{p}. IMPORTANT: The function \code{phi(xval,p)} must be fully vectorized in \code{xval}, i.e. it must output a vector of \code{N} elements if \code{xval} is an \code{N-by-P} array (such as \code{x}). Note that if \code{phi} is given as a function, the argument \code{p.initial} is mandatory.
#' @param p.initial Initial model parameters for fitting the DF.
#' @param n.iterations Maximum number of iterations in the repeated fit-and-debias algorithm to evaluate the maximum likelihood.
#' @param n.resampling Integer (>0) specifying the number of iterations for the resampling of the most likely DF used to evaluate realistic parameter uncertainties with quantiles. If \code{n.resampling = NULL}, no resampling is performed.
#' @param correct.mle.bias The maximum likelihood estimator (MLE) of a finite dataset can be biased – a general property of the ML approach. If \code{TRUE}, \code{dffit} also outputs the parameters, where this estimator bias has been corrected, to first order in \code{1/N}, using jackknifing.
#' @param correct.lss.bias If \code{TRUE} the \code{distance} values are used to correct for the observational bias due to galaxy clustering (large-scale structure).
#' @param lss.sigma Cosmic variance of survey volume. Explicitly, \code{lss.sigma} is interpreted as the expected relative RMS on the total number of galaxies in the survey volume. This value must be determined from a cosmological model.
#' @param write.fit If \code{TRUE}, the best-fitting parameters are displayed in the console.
#' @param x.grid sets the grid on which the numerical integration is performed. This grid must be specified as \code{x.grid = list(x1, ..., xp)}, where \code{xi} is an equally spaced vector with the grid values for the i-th observable. In the case of a MF (a 1-dimensional DF), \code{x.grid = list(seq(xmin,xmax,by=dx))}, where \code{xmin} and \code{xmax} are the minimum/maximum values of the log-mass considered in the numerical integratino and \code{dx} is the grid spacing.
#'
#' @return Returns a structured list. The sublist \code{input} contains all relevant input arguments. The sublist \code{fit} contains all the output arguments of the MLE algorithm. The output can be visualized using \code{\link{mfplot}}, \code{\link{plot.df}} and \code{\link{dfwrite}}.
#'
#' @keywords schechter function
#' @keywords mass function
#' @keywords fit
#'
#' @examples
#' # basic example
#' data = mfdata()
#' df = dffit(data$x, data$selection)
#' mfplot(df, xlim=c(2e6,5e10))
#' 
#' # include measurement errors in fit and also plot histogram of source counts
#' df = dffit(data$x, data$selection, data$x.err)
#' mfplot(df, xlim=c(2e6,5e10), show.data.histogram = TRUE)
#'
#' # show fitted parameter PDFs and covariances
#' dfplotcov(df)
#'
#' # show fitted effective volume function
#' dfplotveff(df)
#'
#' determine Schechter function uncertainties from resampling and evaluate bias-corrected MLE (plotted in red)
#' df = dffit(data$x, data$selection, data$x.err, n.resampling = 1e2, correct.mle.bias = TRUE)
#' mfplot(df, uncertainty.type=3, nbins=10, bin.xmin=6.5, bin.xmax=9.5, xlim=c(2e6,5e10), ylim=c(2e-3,1.5))
#' lines(10^df$fit$evaluation$x,df$fit$evaluation$y.bias.corrected,col='red',lty=2)
#'
#' # evaluate posteriors of the mass measurements and visualize the change in the mass mode between observation and posterior
#' # i.e. the Eddington bias correction
#' df = posterior.data(df)
#' plot(data$x,df$posterior$x.mode.correction,pch=20,xlab='log10(Mass)',ylab='Eddington bias correction = posterior-observed log-mass mode')
#' abline(h=1)
#'
#' @author Danail Obreschkow
#'
#' @export

dffit <- function(x, # normally log-mass, but can be multi-dimensional
                  selection = NULL,
                  x.err = NULL, # Gaussian standard deviations of x
                  distance = NULL,
                  phi = 'Schechter',
                  p.initial = NULL,
                  n.iterations = 100,
                  n.resampling = NULL,
                  correct.mle.bias = FALSE,
                  correct.lss.bias = FALSE,
                  lss.sigma = NULL,
                  write.fit = TRUE,
                  x.grid = list(seq(4,12,0.01))) {

  tStart.global = Sys.time() # timer
  wt.fitting <<- 0
  
  # Input handling ########################################################################

  # Handle x
  if (length(x)<1) stop('Give at least one data point.')
  if (is.null(dim(x))) {
    x = cbind(as.vector(x)) # make col-vector
  } else if (length(dim(x))==1) {
    x = cbind(x)
  } else if (length(dim(x))!=2) {
    stop('x cannot have more than two dimensions.')
  }
  n.data = dim(x)[1]
  n.dim = dim(x)[2]

  # Handle x.err
  if (!is.null(x.err)) {
    if (is.null(dim(x.err))) {
      x.err = cbind(as.vector(x.err)) # make col-vector
    } else if (length(dim(x.err))==1) {
      x.err = cbind(x.err)
    } else  if (length(dim(x.err))==2) {
      if (any(dim(x)!=dim(x.err))) stop('Size of x.err not compatible with size of x.')
    } else if (length(dim(x.err))==3) {
      if (n.dim==1) stop('For one-dimensional distribution function x.err cannot have 3 dimensions.')
      if (!(dim(x.err)[1]==n.data & dim(x.err)[2]==n.dim & dim(x.err)[3]==n.dim)) {
        stop('Size of x.err not compatible with size of x.')
      }
    } else {
      stop('x.err cannot have more than three dimensions.')
    }
    if (min(x)<=0) stop('All values of x.err must be positive.')
  }

  # Handle distance
  if (!is.null(distance)) {
    distance = as.vector(distance)
    if (length(distance)!=n.data) stop('The number of distance values must be equal to the number of data points (=number of columns of x).')
    if (min(distance)<=0) stop('All distance values must be positive.')
  }

  # Handle selection
  if (is.null(selection)) stop('A selection function must be given.')
  veff.list = .make.veff(selection,x)
  veff.function = veff.list$veff.function

  # Handle phi
  if (is.function(phi)) {
    phi.function = phi
    phi.equation = NA
    if (is.null(p.initial)) stop('For user-defined distribution functions initial parameters must be given.')
  } else {
    phi.function <- function(x,p) {dfmodel(x, p, type = phi)}
    phi.equation = dfmodel(output = 'equation', type = phi)
    if (is.null(p.initial)) p.initial = dfmodel(output = 'initial', type = phi)
  }

  # Handle n.iterations
  if (is.null(n.iterations)) stop('n.iterations must be a positive integer.')
  if (n.iterations<1) stop('n.iterations must be a positive integer.')

  # Handle n.resampling
  if (!is.null(n.resampling)) {
    if (n.resampling<10) stop('n.resampling must be 10 or larger.')
  }

  # Handle correct.lss.bias
  if (correct.lss.bias) {
    if (is.null(distance)) stop('Distances must be given of correct.lss.bias = TRUE.')
    if (is.null(lss.sigma)) stop('lss.sigma must be given of correct.lss.bias = TRUE.')
    if (lss.sigma<0) stop('lss.sigma cannot be negative.')
    veff.function.lss = 0 #xxx
  }

  # Handle correct.mle.bias
  if (correct.mle.bias) {
    if (n.data<2) stop('At least two data points must be given if correct.mle.bias = TRUE.')
  }

  # Handle x.grid
  if (!is.list(x.grid)) stop('x.grid must be a list of p vectors, where p is the number of columns of x.')
  if (length(x.grid)!=n.dim) stop('x.grid must be a list of p vectors, where p is the number of columns of x.')
  dx = nx = array(NA,n.dim)
  for (i in seq(n.dim)) {
    dx[i] = x.grid[[i]][2]-x.grid[[i]][1]
    if (dx[i]<=0) stop('x.grid must be made of vectors with monotonically increasing elements.')
    nx[i] = length(x.grid[[i]])
    if (nx[i]<2) stop('x.grid must be made of vectors with at least 2 elements, each.')
    dtmp = x.grid[[i]][3:nx[i]]-x.grid[[i]][2:(nx[i]-1)]
    if (any(abs(dtmp-dx[i])>dx[i]*1e-10)) stop('x.grid must be made of vectors with equally spaced elements.')
  }
  x.mesh.dv = prod(dx)
  x.mesh = array(NA,c(prod(nx),n.dim))
  k1 = 1
  k2 = prod(nx)
  for (i in seq(n.dim)) {
    k2 = k2/nx[i]
    x.mesh[,i] = rep(rep(x.grid[[i]],each=k1),k2)
    k1 = k1*nx[i]
  }
  veff.mesh = veff.function(x.mesh)

  # Find most likely generative model ######################################################
  best.fit = .corefit(x, x.err,
                      veff.mesh,
                      phi.function, p.initial,
                      x.mesh, x.mesh.dv,
                      n.iterations = n.iterations)

  # Initialize output parameters ###########################################################
  phi.fit <- function(x) {phi.function(x,best.fit$p.optimal)}
  source.count <- function(x) {phi.fit(x)*veff.function(x)}
  input = list(data = list(x = x, x.err = x.err, distance = distance),
               selection = list(veff.function = veff.function,
                                veff.values = veff.list$veff.values,
                                veff.fct = veff.list$veff.fct,
                                veff.mesh = veff.mesh),
               distribution.function = list(phi = phi.function,
                                            phi.equation = phi.equation),
               options = list(p.initial = p.initial,
                              n.iterations = n.iterations, n.resampling = n.resampling,
                              correct.lss.bias = correct.lss.bias, lss.sigma = lss.sigma,
                              correct.mle.bias = correct.mle.bias, x.grid = x.grid,
                              x.mesh = x.mesh, x.mesh.dv = x.mesh.dv))
  fit = list(parameters = list(p.optimal = best.fit$p.optimal,
                               p.covariance = best.fit$covariance),
             functions = list(phi.fit = phi.fit,
                              source.count = source.count,
                              logL = best.fit$logL),
             evaluation = list(x = x.mesh,
                               y = apply(x.mesh,1,phi.fit)),
             status = list(n.iterations = best.fit$n.fit.and.debias,
                           converged = best.fit$converged,
                           chain = best.fit$chain))
  df = list(input = input, fit = fit)
  
  # Additional processing for converged solutions ###############################################
  if (df$fit$status$converged) {

    # Bias correction
    if (correct.mle.bias) {df = .correct.mle.bias(df)}

    # Determine uncertainties
    df = .add.Gaussian.errors(df)
    if (!is.null(n.resampling)) {df = .resample(df)}

  }

  # Finalize ###########################################################
  if (write.fit) dfwrite(df)
  df$fit$status$walltime.fitting = wt.fitting
  df$fit$status$walltime.total = as.double(Sys.time())-as.double(tStart.global)
  
  invisible(df)

}

#' @export
#'
.corefit <- function(x, x.err,
                     veff.mesh,
                     phi.function, p.initial,
                     x.mesh, x.mesh.dv,
                     n.iterations, debiasing.first.iteration = TRUE,
                     supress.warning = FALSE) {
  
  tStart.fitting = Sys.time()

  # Input handling
  n.data = dim(x)[1]
  n.dim = dim(x)[2]
  n.mesh = dim(x.mesh)[1]
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
  value.old = Inf
  chain = array(NA,c(n.iterations,length(p.initial)+1))
  
  while (running) {

    k = k+1

    # make unbiased source density function
    rho.unbiased = array(0,n.mesh)
    if (is.null(x.err)) {
      for (i in seq(n.data)) {
        difference = colSums(abs(t(x.mesh)-x[i,]))
        index = which.min(difference)[1]
        rho.unbiased[index] = rho.unbiased[index]+1/x.mesh.dv
      }
    } else {
      if (k==1 & !debiasing.first.iteration) {
        prior = array(1,n.mesh)
      } else {
        # predicted source counts (up to a factor x.mesh.dv)
        prior = phi.function(c(x.mesh),p.initial)*veff.mesh
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
      phi = phi.function(x.mesh,p)
      # safety operations (adding about 50% computation time)
      phi[!is.finite(phi)] = 0
      phi = pmax(.Machine$double.xmin,phi) # also vectorizes the array
      # end safety operations
      return(sum(phi*veff.mesh-log(phi)*rho.unbiased)*x.mesh.dv)
    }

    # maximize ln(L)
    tStart.fitting = Sys.time()
    opt = optim(p.initial,neglogL,hessian=TRUE,control=list(parscale=rep(1,length(p.initial)),reltol=1e-10,abstol=1e-10,maxit=1e3))
    chain[k,] = c(opt$par,opt$value)
    wt.fitting <<- wt.fitting + as.double(Sys.time())-as.double(tStart.fitting)
    
    # assess convergence
    if (is.null(x.err)) { # exist without extra iterations
      converged = opt$convergence==0
      running = FALSE
    } else {
      # asses convergence
      if (k==1) {
        d = 1e98
        d.old = 1e99
      } else {
        d = abs(opt$value-value.old)/n.data
      }
      converged = d>=d.old
      value.old = opt$value
      d.old = d

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
  } else {
    cov = solve(opt$hessian)
  }
  
  return(list(p.optimal = opt$par, covariance = cov, n.fit.and.debias = k,
              converged = converged,
              chain = chain[1:k,],
              logL = function(p) -neglogL(p)))

}

#' @export
#'
.make.veff = function(selection,x) {

  # Generates the function veff.function(xval) from various selection function types.
  # The argument xval must be a M-by-P array. veff.function(xval) returns a vector of
  # M elements.

  if (length(dim(x))==0) stop('x must be a N-by-P array.')
  n.data = dim(x)[1]
  n.dim = dim(x)[2]
  xmin = apply(x,2,min)
  xmax = apply(x,2,max)
  mode = NULL
  veff.values = NULL
  veff.fct = NULL

  # Mode 1: Constant effective volume inside observed domain
  if (is.double(selection)) {
    if (length(selection)==1) {
      if (selection<=0) stop('selection = Vconstant mustt be positive.')
      mode = 1
      veff.function.elemental = function(xval) {
        if (any(xval<xmin) | any(xval>xmax)) {
          return(0)
        } else {
          return(selection)
        }
      }
    }
  }

  # Mode 2: Interpolated effective volume inside observed domain
  if (is.double(selection)) {
    if (length(selection)==n.data) {
      if (min(selection)<=0) stop('All values of selection (=Veff) must be positive.')
      mode = 2
      veff.values = selection
      if (n.dim==1) {
        vapprox = function(xval) {
          f = approxfun(x[,1],1/veff.values,rule=2)
          return(1/f(xval))
        }
      } else if (n.dim==2) {
        vapprox = function(xval) {
          z = 1/(akima::interpp(x[,1],x[,2],1/veff.values,xval[1],xval[2],duplicate='mean'))$z
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
  if (is.function(selection)) {
    mode = 3
    veff.function.elemental = selection
    veff.fct = selection
  }

  # Mode 4: Hybrid of 2 and 3
  if (is.list(selection)) {
    if (length(selection)==2) {
      if (is.double(selection[[1]]) & is.function(selection[[2]])) {
        veff.values = selection[[1]]
        veff.fct = selection[[2]]
        if (min(veff.values)<=0) stop('All values of selection (=Veff) must be positive.')
        mode = 4
        if (n.dim==1) {
          vapprox = function(xval) {
            f = approxfun(x[,1],1/veff.values,rule=2)
            return(1/f(xval))
          }
        } else if (n.dim==2) {
          vapprox = function(xval) {
            z = 1/(akima::interpp(x[,1],x[,2],1/veff.values,xval[1],xval[2],duplicate='mean'))$z
            if (is.na(z)) {return(veff.fct(xval))} else {return(z)}
          }
        } else {
          stop('Linear interpolation of Veff not implemented for DF with more than 2 dimensions. Use a different selection type.')
        }
        veff.function.elemental = function(xval) {
          if (any(xval<xmin) | any(xval>xmax)) {
            return(veff.fct(xval))
          } else {
            return(vapprox(xval))
          }
        }
      }
    }
  }

  # Mode 5: selection given as {f(x,r), dVdr(r), range}
  if (is.list(selection)) {
    if (length(selection)>=2 & length(selection)<=3) {
      if (is.function(selection[[1]]) & is.function(selection[[2]])) {
        if (length(selection)==3) {
          if (is.double(selection[[3]])) {
            if (length(selection[[3]])==2) {
              mode = 5.1
              rmin = selection[[3]][1]
              rmax = selection[[3]][2]
            }
          }
        } else {
          mode = 5
          rmin = 0
          rmax = Inf
        }
        test = try(selection[[1]](NA,NA)*selection[[2]](NA),silent=TRUE)
        if (!is.double(test)) stop('In the argument selection = list(f, dVdr, ...), the functions f(xval,r) and dVdr(r) must not contain if-statements and work for r being a vector.')
        veff.function.elemental = function(xval) {
          f = function(r) {selection[[1]](xval,r)*selection[[2]](r)}
          return(integrate(f,rmin,rmax)$value)
        }
      }
    }
  }

  if (is.null(mode)) stop('Unknown selection.')

  # apply to all elements
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

  return(list(veff.function = veff.function,
              veff.values = veff.values,
              veff.fct = veff.fct))
}

.add.Gaussian.errors <- function(df) {

  # make Gaussian uncertainties of parameters
  cov = df$fit$parameters$p.covariance
  df$fit$parameters$p.sigma = sqrt(diag(cov))

  # make Gaussian uncertainties of DF
  eig = eigen(cov)
  np = length(df$fit$parameters$p.optimal)
  nx = length(df$fit$evaluation$x)
  index = 0
  nsteps = 6 # a larger number of steps leads to a more accurate sampling of the covariance ellipsoid
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
      p.new = df$fit$parameters$p.optimal+v
      y.new[,i] = df$input$distribution.function$phi(df$fit$evaluation$x,p.new)
      if (any(!is.finite(y.new[,i])) | any(y.new[,i]<0)) y.new[,i] = NA
    }
  }

  df$fit$evaluation$y.error.neg = array(NA,nx)
  df$fit$evaluation$y.error.pos = array(NA,nx)
  for (i in seq(nx)) {
    df$fit$evaluation$y.error.neg[i] = df$fit$evaluation$y[i]-min(c(y.new[i,],Inf),na.rm=TRUE)
    df$fit$evaluation$y.error.pos[i] = max(c(y.new[i,],-Inf),na.rm=TRUE)-df$fit$evaluation$y[i]
  }

  return(df)
}

.correct.mle.bias <- function(df) {
  n = dim(df$input$data$x)[1]
  np = length(df$fit$parameters$p.optimal)
  if (n<2) {
    stop('Bias correction requires at least two objects.')
  } else if (n>=1e3) {
    cat('WARNING: bias correction normally not relevant for more than 1000 objects.\n')
  }
  p.new = array(NA,c(np,n))
  for (i in seq(n)) {
    list = setdiff(seq(n),i)
    if (is.null(df$input$data$x.err)) {
      x.err = NULL
    } else if (length(dim(df$input$data$x.err))==2) {
      x.err = as.matrix(df$input$data$x.err[list,])
    } else if (length(dim(df$input$data$x.err))==3) {
      x.err = as.matrix(df$input$data$x.err[list,,])
    }
    cf = .corefit(as.matrix(df$input$data$x[list,]), x.err,
                 df$input$selection$veff.mesh*(n-1)/n,
                 df$input$distribution.function$phi,
                 df$fit$parameters$p.optimal,
                 df$input$options$x.mesh, df$input$options$x.mesh.dv,
                 n.iterations = 1, debiasing.first.iteration = TRUE,
                 supress.warning = TRUE)
    p.new[,i] = cf$p.optimal
  }
  p.reduced = apply(p.new, 1, mean, na.rm = T)
  df$fit$parameters$p.optimal.bias.corrected = n*df$fit$parameters$p.optimal-(n-1)*p.reduced
  df$fit$evaluation$y.bias.corrected = df$input$distribution.function$phi(df$fit$evaluation$x,df$fit$parameters$p.optimal.bias.corrected)
  return(df)
}

.resample <- function(df, seed = 1) {

  # randomly resample and refit the DF
  if (dim(df$input$data$x)[2]>1) stop('resampling only available for one-dimensional distribution functions.')
  set.seed(seed)
  np = length(df$fit$parameters$p.optimal)
  x = df$input$options$x.grid[[1]]
  density = pmax(0,df$fit$evaluation$y*df$input$selection$veff.function(x))
  cum = cumsum(density/sum(density))
  n.data = dim(df$input$data$x)[1]
  p.new = array(NA,c(df$input$options$n.resampling,np))
  for (iteration in seq(df$input$options$n.resampling)) {
    n.new = max(2,rpois(1,n.data))
    x.obs = array(NA,n.new)
    r = runif(n.new)
    for (i in seq(n.new)) {
      index = which.min(abs(cum-r[i]))
      x.obs[i] = x[index]
    }
    p.new[iteration,] = .corefit(cbind(x.obs), NULL, df$input$selection$veff.mesh,
                                 df$input$distribution.function$phi,
                                 df$fit$parameters$p.optimal,
                                 df$input$options$x.mesh, df$input$options$x.mesh.dv,
                                 n.iterations = 1)$p.optimal
  }

  # make parameter quantiles
  q = c(0.02,0.16,0.84,0.98)
  p.quant = array(NA,c(length(q),np))
  for (i in seq(np)) {
    p.quant[,i] = quantile(p.new[,i],q,names=FALSE)
  }
  p.quantile = list(p.quantile.02 = p.quant[1,], p.quantile.16 = p.quant[2,],
                    p.quantile.84 = p.quant[3,], p.quantile.98 = p.quant[4,])

  # make DF quantiles
  s = array(NA,c(df$input$options$n.resampling,length(x)))
  for (iteration in seq(df$input$options$n.resampling)) {
    s[iteration,] = df$input$distribution.function$phi(x,p.new[iteration,])
  }
  y.quant = array(NA,c(4,length(x)))
  for (i in seq(length(x))) {
    list = !is.na(s[,i]) & is.finite(s[,i]) & (s[,i]>0)
    y.quant[,i] = quantile(s[list,i],q,names=FALSE)
  }
  y.quantile = list(y.quantile.02 = y.quant[1,], y.quantile.16 = y.quant[2,],
                    y.quantile.84 = y.quant[3,], y.quantile.98 = y.quant[4,])

  # output parameters
  df$fit$parameters = append(df$fit$parameters, p.quantile)
  df$fit$evaluation = append(df$fit$evaluation, y.quantile)
  return(df)

}
