#' Plot effective volume
#'
#' This routine plots the effective survey volume function \code{Veff(x)} assumed when fitting a one-dimensional generative distribution function (GDF), such as a galaxy mass function. Note that this function \code{Veff(x)} is stored as \code{survey$selection$veff} when fitting a GDF using \code{survey=dffit(...)}.
#'
#' @importFrom magicaxis magaxis magplot
#'
#' @param survey List produced by \code{\link{dffit}}
#' @param xlim 2-element vector with x-axis plotting limits
#' @param ylim 2-element vector with y-axis plotting limits
#' @param xlab Label on x-axis (logarithmic mass or luminosity scale).
#' @param ylab Label on y-axis (volume scale).
#' @param legend Logical flag determining whether the legend will be plotted.
#' @param xpower10 If \code{TRUE}, the model argument x is elevated to the power of 10 in the plots.
#'
#' @seealso See examples in \code{\link{dffit}}. To display \code{Veff(x)} of two-dimensional distribution functions use \code{\link{dfplotveff2}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplotveff <- function(survey,
                       xlim = NULL,
                       ylim = NULL,
                       xlab = 'Observable x',
                       ylab = expression('Effective volume'),
                       legend = TRUE,
                       xpower10 = TRUE) {
  
  blue = '#7777ff'
  n.dim = dim(survey$data$x)[2]
  if (n.dim!=1) {
    if (n.dim==2) {
      stop('Use dfplotveff2 for two-dimensional Veff functions.')
    } else {
      stop('dfplotveff only handles one-dimensional Veff functions.')
    }
  }
  
  list = c(TRUE,FALSE,FALSE,FALSE)
  x = survey$data$x
  dx = max(x)-min(x)
  y = survey$selection$veff(x)
  if (xpower10) {
    if (is.null(xlim)) xlim = 10^(range(x)+dx*c(-0.1,0.1))
    log='xy'
  } else {
    if (is.null(xlim)) xlim = range(x)+dx*c(-0.1,0.1)
    log='y'
  }
  
  par(pty = "m")
  if (is.null(ylim)) ylim = c(1e-5,1)*2*max(y)
  plot(1, 1, type='n', xlim=xlim, ylim=ylim,
       xaxs='i', yaxs='i', xaxt='n', yaxt='n',
       xlab = '', ylab = '', log = log)
  xr = sort(c(survey$grid$x,x))
  if (xpower10) {
    xplot = 10^x
    xrplot = 10^xr
  } else {
    xplot = x
    xrplot = xr
  }
  lines(xrplot,pmax(ylim[1]/2,survey$selection$veff(xr)),col=blue,lwd=3)
  if (!is.null(survey$selection$veff.input.function)) {
    lines(xrplot,survey$selection$veff.input.function(xr))
    list[2] = TRUE
  }
  if (!is.null(survey$selection$veff.input.values)) {
    points(xplot,survey$selection$veff.input.values,pch=20)
    list[3] = TRUE
  }
  if (survey$options$correct.lss.bias) {
    lines(xrplot,survey$selection$veff.no.lss(xr),col='black',lty=2)
    list[4] = TRUE
  }
  magicaxis::magaxis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=3,labels=FALSE,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1)
  if (legend) {
    legend('bottomright',c('Model used for DF fit (veff)','Input function (veff.input.function)',
                           'Input values (veff.input.values)','Model without LSS (veff.no.lss)')[list],
           lwd=c(3,1,NA,1)[list],lty=c(1,1,NA,2),pch=c(NA,NA,20,NA)[list],col=c(blue,'black','black','black')[list],bty='n')
  }
}
