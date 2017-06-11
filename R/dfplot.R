#' Display fitted distribution function
#'
#' This function displays the distribution function fitted using \code{\link{dffit}}.
#'
#' @importFrom magicaxis magaxis magplot
#'
#' @param df List produced by \code{\link{dffit}}
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param xlim 2-element vector with x-axis plotting limits
#' @param ylim 2-element vector with y-axis plotting limits
#' @param show.uncertainties If \code{TRUE}, uncertainties are displayed around the best fit model.
#' @param uncertainty.type \code{1}: plot Gaussian 1-sigma uncertanties propagated from the Hessian matrix of the likelihood. \code{2}: plot 68 percentile region (from 16 to 84 percent). \code{3} plot 68 (16 to 84) and 95 (2 to 98) percentile regions.
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplot <- function(df,
                   xlab = 'Observable x',
                   ylab = expression('Distribution function'~phi),
                   xlim = NULL,
                   ylim = NULL,
                   log = 'y',
                   col = 'blue',
                   show.uncertainties = TRUE,
                   uncertainty.type = NULL) {
  
  r = col2rgb(col)[1]/255
  g = col2rgb(col)[2]/255
  b = col2rgb(col)[3]/255
  
  n.data = dim(df$input$data$x)[1]
  n.dim = dim(df$input$data$x)[2]
  
  if (n.dim==1) {

    # define plot limits
    if (is.null(xlim)) xlim = range(df$fit$evaluation$x)
    if (is.null(ylim)) ylim = c(1e-3*max(df$fit$evaluation$y),2*max(df$fit$evaluation$y))
  
    # open plot
    par(pty = "m")
    plot(1,1,type='n',log=log,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
           xlim = xlim, ylim = ylim, xlab = '', ylab = '',bty='n')
  
    # plot polygons
    if (show.uncertainties & df$fit$status$converged) {
      poly.x = c(df$fit$evaluation$x,rev(df$fit$evaluation$x))
      if (is.null(uncertainty.type)) {
        if (length(df$fit$evaluation$y.quantile.16)>0) {
          uncertainty.type = 2
        } else {
          uncertainty.type = 1
        }
      }
      if ((uncertainty.type>1) & (!length(df$fit$evaluation$y.quantile.16)>0)) stop('Quantiles not available. Use resampling in dffit.')
      if (uncertainty.type == 3) {
        poly.y.95 = pmax(ylim[1],c(df$fit$evaluation$y.quantile.02,rev(df$fit$evaluation$y.quantile.98)))
        polygon(poly.x,poly.y.95,col=rgb(r,g,b,0.2),border=NA)
      }
      if (uncertainty.type >= 2) {
        poly.y.68 = pmax(ylim[1],c(df$fit$evaluation$y.quantile.16,rev(df$fit$evaluation$y.quantile.84)))
        polygon(poly.x,poly.y.68,col=rgb(r,g,b,0.3),border=NA)
      }
      if (uncertainty.type == 1) {
        poly.y.68 = pmax(ylim[1],c(df$fit$evaluation$y-df$fit$evaluation$y.error.neg,
                                   rev(df$fit$evaluation$y+df$fit$evaluation$y.error.pos)))
        polygon(poly.x,poly.y.68,col=rgb(r,g,b,0.3),border=NA)
      }
    }
  
    # plot central fit
    x = df$fit$evaluation$x[df$fit$evaluation$y>0]
    lines(x,df$fit$evaluation$y[df$fit$evaluation$y>0],col=col,lwd=2)
  
    # axes
    magicaxis::magaxis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=3,labels=FALSE,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1)
    box(which = "plot", lty = "solid", lwd = 1)
    
  } else {
    
    stop('dfplot cannot handle distribution functions with more than one dimension.')
    
  }
}
