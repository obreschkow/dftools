#' Example of fitting a step-wise mass function
#'
#' This function generates n mock galaxies with masses and observing errors, drawn from a custom Schechter function and predefined effective volume function. It then uses \code{\link{dffit}} to recover the inpu Schechter function from the noisy data.
#' 
#' @importFrom graphics legend text segments
#' @importFrom pracma erf simpson2d
#' @importFrom stats optim
#' 
#' @param case is an integer from 1 to 4, specifying the example. Choose 1 for an example of fitting a Schechter function mock data drawn from a Schechter function. Choose 2 to fit a Schechter function to mock data drawn from a Schechter function in the presence of large-scale structure. Choose 3 to fit a quasi non-parametric mass function to mock data drawn from a modified Scheckter function. Choose 4 to fit a 2D analytic distribution function to mock data in the the mass-specific angular momentum plane.
#' @param seed A positive integer used as seed for the random number generator.
#' 
#' @examples
#' # Run basic example
#' dfexample(1)
#'
#' @seealso See \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfexample <- function(case = 1, seed = 1) {
  
  if (case == 1) {
    survey = .dfexample_basic(seed)
  } else if (case == 2) {
    survey = .dfexample_lss(seed)
  } else if (case == 3) {
    survey = .dfexample_nonparametric(seed)
  } else if (case == 4) {
    survey = .dfexample_mj(seed)
  }
  
  invisible(survey)
  
}

.dfexample_basic <- function(seed) {
  
  # use parameters
  n = 1e3
  sigma = 0.5
  p.true = c(-2,10,-1.3)
  
  # make rescaled Veff(x)
  cat('Generate mock data with observing errors:\n')
  dat = dfmockdata(n = n, seed = seed, sigma = sigma, p = p.true)
  cat(sprintf('=> %d galaxies\n',dat$n))
  
  # fit
  cat('Fit a Schechter function to the mock data:\n')
  survey = dffit(dat$x, dat$veff, dat$x.err, write.fit = T,
                 xmin = dat$xmin, xmax = dat$xmax, p.initial = dat$p)
  p.fit = survey$fit$parameters$p.optimal
  
  # plot covariances
  dfplotcov(survey,expectation2 =p.true,title='Parameter covariance',model.col='blue',model2.cross.lty = 2)
  
  # plot main plot
  nbins = max(4,round(sqrt(n)/2))
  mfplot(survey, xlim=c(1e6,2e11), ylim=c(1e-4,2), p = p.true, nbins = nbins, bin.xmin=7, bin.xmax=11)
  x = seq(2,10.65,0.01)
  
  # legend
  x = 7; y = 1.5e-2
  points(10^x,y,pch=20)
  segments(10^(x-sigma),y,10^(x+sigma),y,lwd=1)
  text(10^x,y*1.8,expression('Observing error'~sigma),cex=0.85)
  legend(2e6,6e-3,c('Input Schechter function',
                    'Fitted Schechter function with 68%-uncertainties',
                    'Binned input data',
                    'Binned posterior data'),
         bty='n',lwd=c(1,1.5,1,1),
         col=c('black','blue','purple','black'),
         pch=c(NA,NA,20,20),
         lty=c(2,1,1,1),cex=0.85)
  
  # console output
  cat(sprintf('Total computation time       = %.2fs\n',survey$fit$status$walltime.total))
  cat(sprintf('Computation time for fitting = %.2fs\n',survey$fit$status$walltime.fitting))
  
  invisible(survey)
}

.dfexample_lss = function(seed) {
  
  # user parameters
  sigma = 0.3
  
  # survey model
  rmax = 100
  p = c(-2,10,-1.3)
  dmax = function(x) 1e-3*sqrt(10^x)
  f = function(x,r) pracma::erf((dmax(x)-r)/dmax(x)*20)*0.5+0.5
  dVdr = function(r) 0.2*r^2
  g = function(r) 1+0.9*sin((r/100)^0.6*2*pi)
  
  # make data
  cat('Generate mock data with observing errors and large-scale structure (LSS):\n')
  dat = dfmockdata(seed = seed, sigma = sigma, p = p, f = f, g = g, dVdr = dVdr, rmax = rmax)
  dat <<- dat
  cat(sprintf('=> %d galaxies\n',dat$n))
  
  # fit data
  selection = list(f,dVdr,0,rmax)
  cat('Fit mock data without any bias correction:\n')
  survey1 = dffit(dat$x,selection,NULL,r=dat$r,correct.lss.bias = FALSE)
  cat('Fit mock data while correcting for observational errors (Eddington bias):\n')
  survey2 = dffit(dat$x,selection,dat$x.err,r=dat$r,correct.lss.bias = FALSE)
  cat('Fit mock data while correcting for observational errors and LSS:\n')
  survey3 = dffit(dat$x,selection,dat$x.err,r=dat$r,correct.lss.bias = TRUE, lss.weight = function(x) 10^x)
  
  # plot effective volumes
  x = seq(6,12,0.01)
  dfplotveff(survey3,xlab='log10(M/Msun)',legend=FALSE)
  lines(x,dat$veff.lss(x),type='l',ylim=c(0,6),col='black')
  legend('bottomright',c('Recovered model from data used for fit','Input model with LSS used to generate the data','Input model without LSS'),
         lwd=c(2,2,2),col=c('blue','black','red'),bty='n',cex=0.85)
  
  # plot MFs
  mfplot(survey1,xlim=c(1e7,2e11),ylim=c(1e-5,5),nbins=20,bin.xmin=7,bin.xmax=11,col.fit='orange',col.data.input='orange',show.data.histogram = FALSE)
  mfplot(survey2,xlim=c(1e7,2e11),ylim=c(1e-5,5),nbins=20,bin.xmin=7,bin.xmax=11,col.fit='red',col.data.posterior='red',show.input.data=FALSE,add=TRUE,show.data.histogram = FALSE)
  mfplot(survey3,xlim=c(1e7,2e11),ylim=c(1e-5,5),nbins=20,bin.xmin=7,bin.xmax=11,col.fit='blue',col.data.posterior='blue',show.input.data=FALSE,add=TRUE,show.data.histogram = FALSE)
  lines(10^x,survey1$model$gdf(x,p),lty=2)
  legend('bottomleft',c('Uncorrected fit','Including mass errors',
                        'Including errors + LSS','Input model'),
         lwd=c(2,2,2,1),lty=c(1,1,1,2),col=c('orange','red','blue','black'),bty='n',
         cex = 0.85)
  x = 8; y = 1.7e-3
  points(10^x,y,pch=20)
  segments(10^(x-sigma),y,10^(x+sigma),y,lwd=1)
  text(10^x,y*1.8,expression('Observing error'~sigma),cex=0.85)
  
}

.dfexample_nonparametric = function(seed) {

  sigma = 0.5
  n = 1e3
  
  # make analytic non-Schechter gdf
  gdf.analytic = function(x,p) {
    component1 = dfmodel(x,p[1:3])
    component2 = p[4]*exp(-(x-9.4)^2*2)
    return(component1+component2)
  }
  
  # draw random data
  p.true = c(-2,10,-1.3,1)
  dat = dfmockdata(n = n, seed = seed, sigma = sigma, gdf = gdf.analytic, p = p.true,
                   xmax = 13)
  
  # make step-wise linear fitting function
  xpoints = seq(6.5,11,length=10)
  ssc = dfssmodel(xpoints,method='constant')
  ssl = dfssmodel(xpoints,method='linear')
  sss = dfssmodel(xpoints,method='spline')
  
  # fitt step-wise linear function
  surveyc = dffit(dat$x, dat$veff, dat$x.err, gdf = ssc$gdf,
                  p.initial = log10(gdf.analytic(xpoints,p.true)),
                  xmin = dat$xmin, xmax = dat$xmax, write.fit = FALSE)
  surveyl = dffit(dat$x, dat$veff, dat$x.err, gdf = ssl$gdf,
                  p.initial = log10(gdf.analytic(xpoints,p.true)),
                  xmin = dat$xmin, xmax = dat$xmax, write.fit = FALSE)
  surveys = dffit(dat$x, dat$veff, dat$x.err, gdf = sss$gdf,
                  p.initial = log10(gdf.analytic(xpoints,p.true)),
                  xmin = dat$xmin, xmax = dat$xmax, write.fit = FALSE)
  
  # plot
  mfplot(surveyc,xlim=c(1e6,4e11),ylim=c(1e-3,3),show.posterior.data=FALSE,col.data.input='#00000030',show.data.histogram = FALSE)
  mfplot(surveys,xlim=c(1e6,4e11),ylim=c(1e-3,3),show.posterior.data=FALSE,show.input.data=FALSE,add=TRUE,show.data.histogram = FALSE,col.fit='red')
  mfplot(surveyl,xlim=c(1e6,4e11),ylim=c(1e-3,3),show.posterior.data=FALSE,show.input.data=FALSE,add=TRUE,show.data.histogram = FALSE,col.fit='orange')
  x = seq(4,12,0.01)
  lines(10^x,gdf.analytic(x,p.true),lty=2)
  
  x = 9.3; y = 2e-2
  points(10^x,y,pch=20)
  segments(10^(x-sigma),y,10^(x+sigma),y,lwd=1)
  text(10^x,y*1.8,expression('Observing error'~sigma),cex=0.85)
  
  legend('bottomleft',c('Input MF','Data with observing errors drawn from the input MF',
                        'Fitted MF using uniform bins',
                        'Fitted MF using linear elements',
                        'Fitted MF using splines'),
         lwd=c(1,1,2,2,2),lty=c(2,1,1,1,1),col=c('black','grey','blue','orange','red'),bty='n',
         pch=c(NA,20,NA,NA,NA),
         cex = 0.85)
  
  invisible(surveys)
}

.dfexample_mj = function(seed) {
  
  # user parameters
  p.true = c(-2,10,-1.3,2/3,7,0.25)
  sigma.x = 0.25
  sigma.y = 0.1
  corr.xy = 0.6
  
  # initialize
  set.seed(seed+1)
  mrange = c(8,12)
  jrange = c(1,4)
  cov = rbind(c(sigma.x^2,corr.xy*sigma.x*sigma.y),
              c(corr.xy*sigma.x*sigma.y,sigma.y^2))
  
  # define 2D distribution function in the (x1,x2)-plane, where
  # x1 = log10(baryon mass [Msun])
  # x2 = log10(baryon specific angular momentum [kpc km/s])
  gdf = function(x,p) {
    mu = 10^(x[,1]-p[2])
    phi_m = log(10)*10^p[1]*mu^(p[3]+1)*exp(-mu)
    dj = x[,2]-p[4]*(x[,1]-p[5]) # offset from mean M-j relation
    phi_j = exp(-dj^2/2/p[6]^2)
    return(phi_m*phi_j)
  }
  
  # define 2D effective volume function
  veff = function(x) {
    return(2e-10*10^(1.5*x[,1])*(pracma::erf((x[,2]-3)*1)+1))
  }
  
  # expected source count density function
  sc = function(x) {
    return(gdf(x,p.true)*veff(x))
  }
  
  # compute expected number of galaxies
  sc.array = function(x,y) {
    s = array(NA,length(x))
    for (i in seq(length(x))) {
      s[i] = sc(cbind(x[i],y[i]))
    }
    return(array(s,dim(x)))
  }
  n = simpson2d(sc.array,6,14,-1,6)
  if (n>1) {n = round(n)}
  if (n<10) stop('Expected number of galaxies cannot be smaller than 10.')
  if (n>1e5) stop('Expected number of galaxies cannot be larger than 10^5.')
  
  # check if data lies nicely within boundaries
  if (simpson2d(sc.array,mrange[1],mrange[2],jrange[1],jrange[2])<0.9*n) {
    stop('With this parameter choice, more than 10% of the expected sources lie outside the integration range.')
  }
  
  # make observing errors
  x.err = aperm(array(cov,c(2,2,n)),c(3,2,1))
  
  # find maximum of sc(x)
  fminfct = function(x) {
    x = rbind(x)
    return(-sc(x))
  }
  opt = optim(c(10,3),fminfct,method='L-BFGS-B',lower=c(mrange[1],jrange[1]),upper=c(mrange[2],jrange[2]))
  xmax = rbind(opt$par)
  scmax = -opt$value
  
  # draw random sample
  cat('Generate mock galaxies with mass M, specific angular momenta j and covariant errors:\n')
  x = array(NA,c(n,2))
  count = 0
  while (count<n) {
    xtry = rbind(c(runif(1,mrange[1],mrange[2]),runif(1,jrange[1],jrange[2])))
    if (sc(xtry)>runif(1,0,scmax)) {
      count = count+1
      x[count,] = MASS::mvrnorm(1,xtry,x.err[count,,])
    }
  }
  cat(sprintf('=> %d galaxies\n',n))
  
  # fit model to sample
  cat('Fit six-parameter model to the mock data:\n')
  if (corr.xy==0) {x.err = cbind(rep(sigma.x,n),rep(sigma.y,n))} # just to try if everything works well in this mode, too
  survey = dffit(x, veff, x.err, gdf = gdf, p.initial = p.true, xmin = c(mrange[1],jrange[1]), xmax = c(mrange[2],jrange[2]), dx = c(0.1,0.1))
  
  # plot parameter covariances and initial parameters
  dfplotcov(survey, expectation2 = p.true)
  
  # plot DF
  dfplot2(survey, p.ref = p.true, xlab = 'Mass/Msun', ylab = 'j/[kpc km/s]',
          xlim = c(1e8,5e11),ylim=c(10,5e3))
  legend('topleft',c('Mock observations','Observational uncertainties','Input distribution function (DF)','Fitted DF','Source count model from fitted DF'),
         pch=c(20,NA,NA,15,NA),col=c('black','grey','red','blue','purple'),lty=c(NA,1,2,2,1),bty='n',cex=0.85)
  
  # console output
  cat(sprintf('Total computation time       = %.2fs\n',survey$fit$status$walltime.total))
  cat(sprintf('Computation time for fitting = %.2fs\n',survey$fit$status$walltime.fitting))
  
  invisible(survey)
  
}
