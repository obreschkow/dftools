#' Display fitted galaxy mass function
#'
#' This function displays the galaxy mass function (MF) fitted using \code{\link{mffit}}.
#'
#' @importFrom magicaxis magaxis magplot
#'
#' @param mf List produced by \code{\link{mffit}}
#' @param nbins Number of bins to be plotted. This is purely for illustrative purposes. The fitting does not use bins. Choose \code{nbins=NULL} (default) to determine the number of bins automatically or \code{nbins=0} to suppress the bins.
#' @param bin.xmin Left edge of first bin
#' @param bin.xmax Right edge of last bin
#' @param bin.type Integer value defining the type of data to be plotted in bins. 1 = plot raw input masses with Veff values evaluted at these masses, subject to Eddington bias. 2 = plot random masses sampled from posterior mass PDF with corresponding values Veff. 3 = posterior source counts using the full posterior mass PDFs of all sources.
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param xlim 2-element vector with x-axis plotting limits
#' @param ylim 2-element vector with y-axis plotting limits
#' @param show.uncertainties If \code{TRUE}, uncertainties are displayed around the best fit model.
#' @param uncertainty.type \code{1}: plot Gaussian 1-sigma uncertanties propagated from the Hessian matrix of the likelihood. \code{2}: plot 68 percentile region (from 16 to 84 percent). \code{3} plot 68 (16 to 84) and 95 (2 to 98) percentile regions.
#' @param add If \code{TRUE}, the lines are overplotted on the currently open plot.
#' @param col Color of ML fit and uncertainty regions
#' @param lwd Line width of ML fit
#' @param lty Line type of ML fit
#' @param col.bias.correction Line properties of bias corrected MF fit
#' @param lwd.bias.correction Line properties of bias corrected MF fit
#' @param lty.bias.correction Line properties of bias corrected MF fit
#' @param margins Margins (bottom,left,top,right)
#'
#' @return If used as \code{mf = mfplot(mf)} (with \code{nbins=NULL} or \code{nbins>0}), the list \code{mf} is appended an additional entry \code{bin}, containing the binned galaxy data points.
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

mfplot <- function(mf,
                   nbins = NULL,
                   bin.xmin = NULL,
                   bin.xmax = NULL,
                   bin.type = 1,
                   xlab = expression('M [M'['sun']*']'),
                   ylab = expression(phi~'[Mpc'^-3~'dex'^-1~']'),
                   xlim = NULL,
                   ylim = NULL,
                   show.uncertainties = TRUE,
                   uncertainty.type = NULL,
                   show.bias.correction = FALSE,
                   add = FALSE,
                   col = 'blue',
                   lwd = 1.5,
                   lty = 1,
                   col.bin = 'black',
                   col.bias.correction = col,
                   lwd.bias.correction = lwd,
                   lty.bias.correction = 2,
                   margins=c(5.1,4.1,4.1,2.1)) {

  r = col2rgb(col)[1]/255
  g = col2rgb(col)[2]/255
  b = col2rgb(col)[3]/255

  if (dim(mf$input$data$x)[2]>1) stop('mfplot can only handle one-dimensional distribution functions.')
  n.data = dim(mf$input$data$x)[1]

  # make binned data
  make.bins = is.null(nbins) || (!is.na(nbins) && (nbins>0))
  if (make.bins) {

    # make bin intervals
    x = mf$input$data$x
    bin = list()
    if (is.null(nbins)) {
      bin$n = min(100,round(sqrt(n.data)))
    } else {
      bin$n = nbins
    }
    if (is.null(bin.xmin)) {
      bin$xmin = min(x)-(max(x)-min(x))/bin$n*0.25
    } else {
      bin$xmin = bin.xmin
    }
    if (is.null(bin.xmax)) {
      bin$xmax = max(x)+(max(x)-min(x))/bin$n*0.25
    } else {
      bin$xmax = bin.xmax
    }
    bin$wx = bin$xmax-bin$xmin
    bin$dx = bin$wx/bin$n
    bin$xcenter = bin$xmin+(seq(bin$n)-0.5)*bin$dx

    # fill data into bins
    if (bin.type <= 2) {

      if (bin.type == 1) {
        xval = x
      } else {
        xval = mf$posterior$x.rand
      }
      v = mf$input$selection$veff.function(xval)

      bin$phi = bin$count = bin$xmean = array(0,bin$n)
      for (i in seq(n.data)) {
        k = floor((xval[i]-bin$xmin)/bin$wx*0.99999999*bin$n)+1
        if (k>=1 & k<=bin$n) {
          bin$phi[k] = bin$phi[k]+1/bin$dx/v[i]
          bin$count[k] = bin$count[k]+1
          bin$xmean[k] = bin$xmean[k]+x[i]
        }
      }
      bin$xmean = bin$xmean/pmax(1,bin$count)
      mf = append(mf,list(bin=bin))

    } else if (bin.type == 3) {

      xg = mf$input$options$x.grid[[1]]
      dx = xg[2]-xg[1]
      for (k in seq(bin$n)) {
        list = floor((xg-bin$xmin)/bin$wx*0.99999999*bin$n)+1==k
        bin$xmean[k] = sum(mf$posterior$source.count.density[list]*xg[list])/sum(mf$posterior$source.count.density[list])
        bin$count[k] = sum(mf$posterior$source.count.density[list])*dx
        bin$phi[k] = sum(mf$posterior$source.count.density[list]/mf$input$selection$veff.function(xg[list]))/sum(list)
      }

    }

  }

  # define plot limits
  if (is.null(xlim)) {
    xlim = 10^range(mf$fit$evaluation$x)
  }
  if (is.null(ylim)) {
    if (make.bins) {
      ylim = c(1e-3*max(mf$fit$evaluation$y),2*max(max(mf$fit$evaluation$y),bin$phi))
    } else {
      ylim = c(1e-3*max(mf$fit$evaluation$y),2*max(mf$fit$evaluation$y))
    }
  }

  # open plot
  if (add==FALSE) {
    par(pty = "m")
    par(mar=margins)
    plot(1,1,type='n',log='xy',xaxs='i',yaxs='i',xaxt='n',yaxt='n',
         xlim = xlim, ylim = ylim, xlab = '', ylab = '',bty='n')
  }

  # plot polygons
  if (show.uncertainties & mf$fit$status$converged) {
    poly.x = 10^c(mf$fit$evaluation$x,rev(mf$fit$evaluation$x))
    if (is.null(uncertainty.type)) {
      if (length(mf$fit$evaluation$y.quantile.16)>0) {
        uncertainty.type = 2
      } else {
        uncertainty.type = 1
      }
    }
    if ((uncertainty.type>1) & (!length(mf$fit$evaluation$y.quantile.16)>0)) stop('Quantiles not available. Use resampling in mffit.')
    if (uncertainty.type == 3) {
      poly.y.95 = pmax(ylim[1],c(mf$fit$evaluation$y.quantile.02,rev(mf$fit$evaluation$y.quantile.98)))
      polygon(poly.x,poly.y.95,col=rgb(r,g,b,0.2),border=NA)
    }
    if (uncertainty.type >= 2) {
      poly.y.68 = pmax(ylim[1],c(mf$fit$evaluation$y.quantile.16,rev(mf$fit$evaluation$y.quantile.84)))
      polygon(poly.x,poly.y.68,col=rgb(r,g,b,0.3),border=NA)
    }
    if (uncertainty.type == 1) {
      poly.y.68 = pmax(ylim[1],c(mf$fit$evaluation$y-mf$fit$evaluation$y.error.neg,
                                 rev(mf$fit$evaluation$y+mf$fit$evaluation$y.error.pos)))
      polygon(poly.x,poly.y.68,col=rgb(r,g,b,0.3),border=NA)
    }
  }

  # plot central fit
  x = mf$fit$evaluation$x[mf$fit$evaluation$y>0]
  lines(10^x,mf$fit$evaluation$y[mf$fit$evaluation$y>0],col=col,lwd=lwd,lty=lty)
  if (show.bias.correction & mf$fit$status$converged) {
    if (length(mf$fit$parameters$p.optimal.bias.corrected)==0) stop('Bias corrected MLE parameters not available. Use bias.correction in mffit.')
    lines(10^x,mf$input$mass.function(x,mf$fit$parameters$p.optimal.bias.corrected),col=col.bias.correction,lwd=lwd.bias.correction,lty=lty.bias.correction)
  }

  # plot binned data
  if (make.bins) {
    list = bin$phi>0
    points(10^bin$xmean[list],bin$phi[list],pch=20,col=col.bin)
    f.16 = pmax(1e-3,qpois(0.16,bin$count[list])/bin$count[list])
    f.84 = qpois(0.84,bin$count[list])/bin$count[list]
    segments(10^bin$xmean[list],bin$phi[list]*f.16,10^bin$xmean[list],bin$phi[list]*f.84,col=col.bin)
    segments(10^(bin$xmin+seq(0,bin$n-1)[list]*bin$dx),bin$phi[list],
             10^(bin$xmin+seq(1,bin$n)[list]*bin$dx),bin$phi[list],col=col.bin)
  }

  # axes
  if (!add) {
    magicaxis::magaxis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=3,labels=FALSE,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1)
  }

  box(which = "plot", lty = "solid", lwd = 1)

  invisible(mf)
}
