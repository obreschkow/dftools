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
#' @param log String specifying the log-axes as in \code{\link{plot}}.
#' @param xpower10 If \code{TRUE}, the model argument x is elevated to the power of 10 in the plots.
#' @param col Color of fit (see \code{\link{plot}})
#' @param lwd Line width of fit (see \code{\link{plot}})
#' @param lty Line type of fit (see \code{\link{plot}})
#' @param show.data.points If \code{TRUE}, binned data points are displayed.
#' @param show.data.histogram If \code{TRUE}, a histogram of source counts is displayed in a bottom panel.
#' @param show.uncertainties If \code{TRUE}, uncertainties are displayed around the best fit model.
#' @param uncertainty.type \code{1}: plot Gaussian 1-sigma uncertanties propagated from the Hessian matrix of the likelihood. \code{2}: plot 68 percentile region (from 16 to 84 percent). \code{3} plot 68 (16 to 84) and 95 (2 to 98) percentile regions.
#' @param show.bias.correction If \code{TRUE}, the bias corrected MLE is shown instead of the native ML parameters.
#' @param add If \code{TRUE}, the lines are overplotted on the currently open plot.
#' @param nbins Number of bins to be plotted; must be larger than 0. Choose \code{nbins=NULL} (default) to determine the number of bins automatically.
#' @param bin.xmin Left edge of first bin
#' @param bin.xmax Right edge of last bin
#' @param bin.type Integer value defining the type of data to be plotted in bins. 1 = plot raw input masses with Veff values evaluted at these masses, subject to Eddington bias. 2 = plot random masses sampled from posterior mass PDF with corresponding values Veff. 3 = posterior source counts using the full posterior mass PDFs of all sources.
#' @param col.data Color of binned data
#' @param cex.data Size of binned data
#' @param lwd.data Line width of binned data
#' @param col.hist Color of source count histogram
#' @param margins Margins (bottom,left,top,right)
#' 
#' @return Returns the input list \code{df} with the additional sub-list \code{df$bin} that contains the binned data.
#' 
#' @seealso For optimized plotting of galaxy mass functions, use the derived function \code{\link{mfplot}}. See examples in \code{\link{dffit}}.
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
                   xpower10 = FALSE,
                   col = 'blue',
                   lwd = 2,
                   lty = 1,
                   show.data.points = FALSE,
                   show.data.histogram = FALSE,
                   show.uncertainties = TRUE,
                   uncertainty.type = NULL,
                   show.bias.correction = FALSE,
                   add = FALSE,
                   nbins = NULL,
                   bin.xmin = NULL,
                   bin.xmax = NULL,
                   bin.type = 1,
                   col.data = 'black',
                   cex.data = 1,
                   lwd.data = 1,
                   col.hist = 'grey',
                   margins=c(5.1,4.1,4.1,2.1)) {
  
  r = col2rgb(col)[1]/255
  g = col2rgb(col)[2]/255
  b = col2rgb(col)[3]/255
  
  n.data = dim(df$input$data$x)[1]
  n.dim = dim(df$input$data$x)[2]
  
  if (n.dim!=1) {
    if (n.dim==2) {
      stop('Use dfplot2 for two-dimensional distribution functions.')
    } else {
      stop('dfplot only handles one-dimensional distribution functions.')
    }
  }

  # define plot limits
  if (is.null(xlim)) {
    xlim = range(df$fit$evaluation$x)
    if (xpower10) xlim = 10^xlim
  }
  if (is.null(ylim)) ylim = c(1e-3*max(df$fit$evaluation$y),2*max(df$fit$evaluation$y))
  
  # bin data
  df = .bin.data(df,nbins,bin.type,bin.xmin,bin.xmax)
  
  # open plot
  if (!add) {
    par(pty = 'm')
    par(mar = margins)
    if (show.data.histogram) {
      plot(0,0,type='n',yaxs='i',xaxt='n',yaxt='n',xlim=xlim,ylim=c(0,1),xlab='',ylab='',bty='n')
      .plotSub(0,1,0.2,1)
    }
    plot(1,1,type='n',log=log,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
           xlim = xlim, ylim = ylim, xlab = '', ylab = '',bty='n')
  }

  # plot uncertainty regions
  if (show.uncertainties & df$fit$status$converged) {
    poly.x = c(df$fit$evaluation$x,rev(df$fit$evaluation$x))
    if (xpower10) poly.x = 10^poly.x
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
  if (xpower10) x = 10^x
  if (show.bias.correction & df$fit$status$converged) {
    if (length(df$fit$parameters$p.optimal.bias.corrected)==0) stop('Bias corrected MLE parameters not available. Use bias.correction in dffit.')
    lines(x,df$input$distribution.function$phi(df$fit$evaluation$x[df$fit$evaluation$y>0],
                                               df$fit$parameters$p.optimal.bias.corrected),col=col,lwd=lwd,lty=lty)
  } else {
    lines(x,df$fit$evaluation$y[df$fit$evaluation$y>0],col=col,lwd=lwd,lty=lty)
  }

  # plot binned data points
  if (show.data.points) {
    list = df$bin$phi>0
    f.16 = pmax(1e-3,qpois(0.16,df$bin$count[list])/df$bin$count[list])
    f.84 = qpois(0.84,df$bin$count[list])/df$bin$count[list]
    if (xpower10) {
      points(10^df$bin$xmean[list],df$bin$phi[list],pch=20,col=col.data,cex=cex.data)
      segments(10^df$bin$xmean[list],df$bin$phi[list]*f.16,10^df$bin$xmean[list],df$bin$phi[list]*f.84,col=col.data,lwd=lwd.data)
      segments(10^(df$bin$xmin+seq(0,df$bin$n-1)[list]*df$bin$dx),df$bin$phi[list],
               10^(df$bin$xmin+seq(1,df$bin$n)[list]*df$bin$dx),df$bin$phi[list],col=col.data,lwd=lwd.data)
    } else {
      points(df$bin$xmean[list],df$bin$phi[list],pch=20,col=col.data,cex=cex.data)
      segments(df$bin$xmean[list],df$bin$phi[list]*f.16,df$bin$xmean[list],df$bin$phi[list]*f.84,col=col.data,lwd=lwd.data)
      segments((df$bin$xmin+seq(0,df$bin$n-1)[list]*df$bin$dx),df$bin$phi[list],
               (df$bin$xmin+seq(1,df$bin$n)[list]*df$bin$dx),df$bin$phi[list],col=col.data,lwd=lwd.data)
    }
  }
  
  # axes
  if (!add) {
    magicaxis::magaxis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1)
  }
  
  # plot binned data histogram
  if (show.data.histogram) {
    .plotSub(0,1,0,0.2)
    ymax = max(df$bin$count)*1.2
    if (length(grep('x',log))==1) {lg='x'} else {lg=''}
    plot(1,1,type='n',log=lg,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
         xlim = xlim, ylim = c(0,ymax), xlab = '', ylab = '',bty='n')
    xbin = rep(df$bin$xmin+seq(0,df$bin$n)*df$bin$dx,each=2)
    if (xpower10) xbin = 10^xbin
    xhist = c(xlim[1],xbin,xlim[2])
    yhist = c(0,0,rep(df$bin$count,each=2),0,0)
    polygon(xhist,yhist,col=col.hist,border = NA)
    par(xpd=TRUE)
    lines(xlim,rep(ymax,2))
    par(xpd=FALSE)
    magicaxis::magaxis(side=2,ylab='Counts',lwd=NA,lwd.ticks=1,minorn=0)
    magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1,minorn=0)
    .plotSubEnd()
  }

  # axes
  if (!add) {
    magicaxis::magaxis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=3,labels=FALSE,lwd=NA,lwd.ticks=1)
    box(which = "plot", lty = "solid", lwd = 1)
    par(pty = "m")
  }
  
  invisible(df)
}

.bin.data = function(df,nbins,bin.type,bin.xmin,bin.xmax) {
  
  # initialize
  x = df$input$data$x
  bin = list(type = bin.type)
  n.data = length(x)
  
  # determine number of bins
  if (is.null(nbins)) {
    bin$n = min(100,round(sqrt(n.data)))
  } else {
    if (nbins<=0) stop('Choose more than 0 bins.')
    bin$n = nbins
  }
  
  # make bin intervals
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
  wx = bin$xmax-bin$xmin
  bin$dx = wx/bin$n
  bin$xcenter = bin$xmin+(seq(bin$n)-0.5)*bin$dx
    
  # fill data into bins
  if (bin$type <= 2) {
    
    if (bin$type == 1) {
      xval = x
    } else {
      xval = df$posterior$x.rand
    }
    v = df$input$selection$veff.function(xval)
    
    bin$phi = bin$count = bin$xmean = array(0,bin$n)
    for (i in seq(n.data)) {
      k = floor((xval[i]-bin$xmin)/wx*0.99999999*bin$n)+1
      if (k>=1 & k<=bin$n) {
        bin$phi[k] = bin$phi[k]+1/bin$dx/v[i]
        bin$count[k] = bin$count[k]+1
        bin$xmean[k] = bin$xmean[k]+xval[i]
      }
    }
    bin$xmean = bin$xmean/pmax(1,bin$count)
    
  } else if (bin$type == 3) {
    
    xg = df$input$options$x.grid[[1]]
    dx = xg[2]-xg[1]
    for (k in seq(bin$n)) {
      list = floor((xg-bin$xmin)/wx*0.99999999*bin$n)+1==k
      bin$xmean[k] = sum(df$posterior$source.count.density[list]*xg[list])/sum(df$posterior$source.count.density[list])
      bin$count[k] = sum(df$posterior$source.count.density[list])*dx
      bin$phi[k] = sum(df$posterior$source.count.density[list]/df$input$selection$veff.function(xg[list]))/sum(list)
    }
  }
  df = append(df,list(bin=bin))
  return(df)  
}

.plotSub <- function(xleft=0.1,xright=0.3,ybottom=0.1,ytop=0.3) {
  pomdglobal <<- par(no.readonly = T)
  par(omd=c(0,1,0,1))
  xmarg = sum(par()$mai[c(2,4)])
  xplot = par()$pin[1]
  ymarg = sum(par()$mai[c(1,3)])
  yplot = par()$pin[2]
  par(new=T,omd=c(xleft*xplot/(xplot+xmarg),(xright*xplot+xmarg)/(xplot+xmarg),ybottom*yplot/(yplot+ymarg),(ytop*yplot+ymarg)/(yplot+ymarg)))
}

.plotSubEnd <- function() {
  par(oma=c(0,0,0,0))
  par(omi=c(0,0,0,0))
  par(omd=c(0,1,0,1))
}

