#' Plot effective volume
#'
#' This routine plots the effective survey volume function \code{Veff(x)} assumed when fitting a one-dimensional generative distribution function (GDF), such as a galaxy mass function. Note that this function \code{Veff(x)} is stored as \code{survey$selection$veff} when fitting a GDF using \code{survey=dffit(...)}.
#'
#' @importFrom magicaxis magaxis magplot
#'
#' @param df List produced by \code{\link{dffit}}
#' @param xlab Label on x-axis (logarithmic mass or luminosity scale).
#' @param ylab Label on y-axis (volume scale).
#'
#' @seealso See examples in \code{\link{dffit}}. To display \code{Veff(x)} of two-dimensional distribution functions use \code{\link{dfplotveff2}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplotveff <- function(survey,
                       xlab = 'Observable x',
                       ylab = expression('Effective volume')) {
  
  n.dim = dim(survey$data$x)[2]
  if (n.dim!=1) {
    if (n.dim==2) {
      stop('Use dfplotveff2 for two-dimensional Veff functions.')
    } else {
      stop('dfplotveff only handles one-dimensional Veff functions.')
    }
  }
  
  list = c(TRUE,FALSE,FALSE)
  x = survey$data$x
  dx = max(x)-min(x)
  par(pty = "m")
  ylim = c(1e-4,1)*2*max(survey$selection$veff(x))
  plot(1, 1, type='n', xlim=range(x)+dx*c(-0.1,0.1), ylim=ylim,
       xaxs='i', yaxs='i', xaxt='n', yaxt='n',
       xlab = '', ylab = '', log='y')
  xr = sort(c(survey$grid$x,x))
  lines(xr,pmax(ylim[1]/2,survey$selection$veff(xr)),col='blue',lwd=3)
  if (!is.null(survey$selection$veff.input.function)) {
    lines(xr,survey$selection$veff.input.function(xr))
    list[2] = TRUE
  }
  if (!is.null(survey$selection$veff.input.values)) {
    points(x,survey$selection$veff.input.values,pch=20)
    list[3] = TRUE
  }
  magicaxis::magaxis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=3,labels=FALSE,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1)
  legend('bottomright',c('Model used for DF fit (veff.function)','Input model (veff.fct)','Input values (veff.values)')[list],
         lwd=c(3,1,NA)[list],pch=c(NA,NA,20)[list],col=c('blue','black','black')[list],bty='n')
}
