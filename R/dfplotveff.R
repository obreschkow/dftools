#' Plot effective volume
#'
#' This function plots the function \code{Veff(x)} used for the mass function fit, as stored in the sublist \code{fit} of the output produced by \code{\link{dffit}}.
#'
#' @importFrom magicaxis magaxis magplot
#'
#' @param df List produced by \code{\link{dffit}}
#' @param xlab Label on x-axis (logarithmic mass or luminosity scale).
#' @param ylab Label on y-axis (volume scale).
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplotveff <- function(df,
                      xlab = 'Observable x',
                      ylab = expression('Distribution function'~phi)) {
  
  n.dim = dim(df$input$data$x)[2]
  if (n.dim!=1) {
    if (n.dim==2) {
      stop('Use dfplotveff2 for two-dimensional Veff functions.')
    } else {
      stop('dfplotveff only handles one-dimensional Veff functions.')
    }
  }
  
  list = c(TRUE,FALSE,FALSE)
  x = df$input$data$x
  dx = max(x)-min(x)
  par(pty = "m")
  plot(0, 0, type='n', xlim=range(x)+dx*c(-0.1,0.1), ylim=c(0,1.1*max(df$input$selection$veff.function(x))),
       xaxs='i', yaxs='i', xaxt='n', yaxt='n',
       xlab = '', ylab = '')
  xr = sort(c(df$fit$evaluation$x,x))
  lines(xr,df$input$selection$veff.function(xr),col='blue',lwd=3)
  if (!is.null(df$input$selection$veff.fct)) {
    lines(xr,df$input$selection$veff.fct(xr))
    list[2] = TRUE
  }
  if (!is.null(df$input$selection$veff.values)) {
    points(x,df$input$selection$veff.values,pch=20)
    list[3] = TRUE
  }
  magicaxis::magaxis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=3,labels=FALSE,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1)
  legend('bottomright',c('Model used for DF fit (veff.function)','Input model (veff.fct)','Input values (veff.values)')[list],
         lwd=c(3,1,NA)[list],pch=c(NA,NA,20)[list],col=c('blue','black','black')[list],bty='n')
}
