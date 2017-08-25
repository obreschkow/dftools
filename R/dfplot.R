#' Display fitted generative distribution function
#'
#' This function displays a one-dimensional generative distribution function fitted using \code{\link{dffit}}.
#'
#' @importFrom magicaxis magaxis magplot
#' @importFrom graphics plot lines points polygon box
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats qpois
#'
#' @param survey List produced by \code{\link{dffit}}
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ylab.histogram y-axis label for histogram
#' @param xlim 2-element vector with x-axis plotting limits
#' @param ylim 2-element vector with y-axis plotting limits
#' @param log String specifying the log-axes as in \code{\link{plot}}.
#' @param p Parameters of a reference distribution function do be over-plotted to the fitted function. Choose \code{NULL} to avoid plotting a reference function.
#' @param veff Function of x to be plotted over the histogram, if \code{show.data.histogram = TRUE}.
#' @param xpower10 If \code{TRUE}, the model argument x is elevated to the power of 10 in the plots.
#' @param show.input.data If \code{TRUE}, the input data is shown in bins. Each bin values is simply the sum 1/Veff(x) of the observed x-values in this bin.
#' @param show.posterior.data If \code{TRUE}, the posterior data, constructucted from all the individual posterior PDFs of the observed data, are snown in bins. Note that posterior data only exists of the fitted data is uncertain, e.g. if the argument \code{x.err} in \code{\link{dffit}} was non-zero.
#' @param show.data.histogram If \code{TRUE}, a histogram of source counts, based on the input data, is displayed in a bottom panel. Choose \code{NA}, to leave space for the histogram without drawing it.
#' @param show.uncertainties If \code{TRUE}, uncertainties are displayed around the best fit model.
#' @param uncertainty.type \code{1}: plot Gaussian 1-sigma uncertanties propagated from the Hessian matrix of the likelihood. \code{2}: plot 68 percentile region (from 16 to 84 percent). \code{3} plot 68 (16 to 84) and 95 (2 to 98) percentile regions.
#' @param show.bias.correction If \code{TRUE}, the bias corrected MLE is shown instead of the native ML parameters.
#' @param add If \code{TRUE}, the lines are overplotted on the currently open plot.
#' @param nbins Number of bins to be plotted; must be larger than 0. Choose \code{nbins=NULL} (default) to determine the number of bins automatically.
#' @param bin.xmin Left edge of first bin
#' @param bin.xmax Right edge of last bin
#' @param col.fit Color of fit (see \code{\link{plot}})
#' @param lwd.fit Line width of fit (see \code{\link{plot}})
#' @param lty.fit Line type of fit (see \code{\link{plot}})
#' @param col.data.input Color of binned input data
#' @param cex.data.input Size of binned input data
#' @param lwd.data.input Line width of binned input data
#' @param col.data.posterior Color of binned posterior data
#' @param cex.data.posterior Size of binned posterior data
#' @param lwd.data.posterior Line width of binned posterior data
#' @param col.hist Color of source count histogram
#' @param col.ref Color of reference distribution function
#' @param lwd.ref Line width of reference distribution function
#' @param lty.ref Line type of reference distribution function
#' @param col.veff Color of reference function \code{veff}
#' @param lwd.veff Line width of reference function \code{veff}
#' @param lty.veff Line type of reference function \code{veff}
#' @param margins Margins (bottom,left,top,right)
#' 
#' @return Returns the input list \code{survey} with the additional sub-list \code{survey$bin} that contains the binned data.
#' 
#' @seealso For optimized plotting of galaxy mass functions, use the derived function \code{\link{mfplot}}. As an example run \code{dfexample(1)}. See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplot <- function(survey,
                   xlab = 'Observable x',
                   ylab = expression('Generative distribution function'~phi),
                   ylab.histogram = 'Selection',
                   xlim = NULL,
                   ylim = NULL,
                   log = 'y',
                   p = NULL,
                   veff = NULL,
                   xpower10 = FALSE,
                   show.input.data = TRUE,
                   show.posterior.data = TRUE,
                   show.data.histogram = FALSE,
                   show.uncertainties = TRUE,
                   uncertainty.type = NULL,
                   show.bias.correction = FALSE,
                   add = FALSE,
                   nbins = NULL,
                   bin.xmin = NULL,
                   bin.xmax = NULL,
                   col.fit = 'blue',
                   lwd.fit = 2,
                   lty.fit = 1,
                   col.data.input = 'purple',
                   cex.data.input = 1,
                   lwd.data.input = 1,
                   col.data.posterior = 'black',
                   cex.data.posterior = 1,
                   lwd.data.posterior = 1,
                   col.hist = 'grey',
                   col.ref = 'black',
                   lwd.ref = 1,
                   lty.ref = 2,
                   col.veff = '#666666',
                   lwd.veff = 1.5,
                   lty.veff = 5,
                   margins=c(5.1,4.1,4.1,2.1)) {
  
  r = col2rgb(col.fit)[1]/255
  g = col2rgb(col.fit)[2]/255
  b = col2rgb(col.fit)[3]/255
  
  n.data = dim(survey$data$x)[1]
  n.dim = dim(survey$data$x)[2]
  
  if (n.dim!=1) {
    if (n.dim==2) {
      stop('Use dfplot2 for two-dimensional distribution functions.')
    } else {
      stop('dfplot only handles one-dimensional distribution functions.')
    }
  }

  # define plot limits
  if (is.null(xlim)) {
    xlim = range(survey$data$x)
    xlim = xlim+c(-0.1,0.1)*(xlim[2]-xlim[1])
    if (xpower10) xlim = 10^xlim
  }
  if (is.null(ylim)) {
    ylim = max(survey$fit$gdf(survey$data$x))*c(2e-4,2)
  }
  
  # bin data
  survey = .bin.data(survey,nbins,bin.xmin,bin.xmax)
  
  # open plot
  if (!add) {
    par(pty = 'm')
    par(mar = margins)
    if (is.na(show.data.histogram) | show.data.histogram) {
      plot(0,0,type='n',yaxs='i',xaxt='n',yaxt='n',xlim=xlim,ylim=c(0,1),xlab='',ylab='',bty='n')
      .plotSub(0,1,0.2,1)
    }
    plot(1,1,type='n',log=log,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
           xlim = xlim, ylim = ylim, xlab = '', ylab = '',bty='n')
  } else {
    if (is.na(show.data.histogram) | show.data.histogram) {
      .plotSub(0,1,0.2,1)
      plot(1,1,type='n',log=log,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
           xlim = xlim, ylim = ylim, xlab = '', ylab = '',bty='n')
    }
  }
  
  # plot uncertainty regions
  if (show.uncertainties & survey$fit$status$converged) {
    poly.x = c(survey$grid$x,rev(survey$grid$x))
    if (xpower10) poly.x = 10^poly.x
    if (is.null(uncertainty.type)) {
      if (length(survey$grid$gdf.quantile.16)>0) {
        uncertainty.type = 2
      } else {
        uncertainty.type = 1
      }
    }
    if ((uncertainty.type>1) & (!length(survey$grid$gdf.quantile.16)>0)) stop('Quantiles not available. Use resampling in dffit.')
    if (uncertainty.type == 3) {
      poly.y.95 = pmax(ylim[1],c(survey$grid$gdf.quantile.02,rev(survey$grid$gdf.quantile.98)))
      list = is.finite(poly.x) & is.finite(poly.y.95)
      polygon(poly.x[list],poly.y.95[list],col=rgb(r,g,b,0.15),border=NA)
    }
    if (uncertainty.type >= 2) {
      poly.y.68 = pmax(ylim[1],c(survey$grid$gdf.quantile.16,rev(survey$grid$gdf.quantile.84)))
      list = is.finite(poly.x) & is.finite(poly.y.68)
      polygon(poly.x[list],poly.y.68[list],col=rgb(r,g,b,0.25),border=NA)
    }
    if (uncertainty.type == 1) {
      poly.y.68 = pmax(ylim[1],c(survey$grid$gdf-survey$grid$gdf.error.neg,
                                 rev(survey$grid$gdf+survey$grid$gdf.error.pos)))
      list = is.finite(poly.x) & is.finite(poly.y.68)
      polygon(poly.x[list],poly.y.68[list],col=rgb(r,g,b,0.25),border=NA)
    }
  }
  
  # plot central fit
  x = survey$grid$x[survey$grid$gdf>0]
  if (xpower10) x = 10^x
  if (show.bias.correction & survey$fit$status$converged) {
    if (length(survey$fit$parameters$p.optimal.bias.corrected)==0) stop('Bias corrected MLE parameters not available. Use bias.correction in dffit.')
    lines(x,survey$input$distribution.function$phi(survey$grid$x[survey$grid$gdf>0],
                                               survey$fit$parameters$p.optimal.bias.corrected),col=col.fit,lwd=lwd.fit,lty=lty.fit)
  } else {
    lines(x,survey$grid$gdf[survey$grid$gdf>0],col=col.fit,lwd=lwd.fit,lty=lty.fit)
  }
  
  # plot reference distribution function
  if (!is.null(p)) {
    if (xpower10) {
      lines(x,survey$model$gdf(log10(x),p),col=col.ref,lwd=lwd.ref,lty=lty.ref)
    } else {
      lines(x,survey$model$gdf(x,p),col=col.ref,lwd=lwd.ref,lty=lty.ref)
    }
  }

  # plot binned input data points
  bin = list()
  for (mode in seq(2)) {
    if (mode==1) {
      show = show.input.data
      if (show) {
        bin$count = survey$bin$count.input
        bin$gdf = survey$bin$gdf.input
        bin$xmean = survey$bin$xmean.input
        col.data = col.data.input
        cex.data = cex.data.input
        lwd.data = lwd.data.input
      }
    } else {
      show = show.posterior.data & !is.null(survey$data$x.err) & !is.null(survey$grid$effective.counts)
      if (show) {
        bin$count = survey$bin$count.posterior
        bin$gdf = survey$bin$gdf.posterior
        bin$xmean = survey$bin$xmean.posterior
        col.data = col.data.posterior
        cex.data = cex.data.posterior
        lwd.data = lwd.data.posterior
      }
    }
    if (show) {
      list = bin$gdf>0
      pm = 0.05
      f.16 = qpois(0.16,bin$count)/bin$count
      f.84 = qpois(0.84,bin$count)/bin$count
      upper = f.16<pm
      f.16 = pmax(f.16,pm)
      if (xpower10) {
        points(10^bin$xmean[list],bin$gdf[list],pch=20,col=col.data,cex=cex.data)
        segments(10^bin$xmean[list],bin$gdf[list]*f.16[list],10^bin$xmean[list],bin$gdf[list]*f.84[list],col=col.data,lwd=lwd.data)
        points(10^bin$xmean[upper],bin$gdf[upper]*pm,pch=25,cex=cex.data*0.7,col=col.data,bg=col.data)
        #stop('sdf')
        segments(10^(survey$bin$xmin+seq(0,survey$bin$n-1)[list]*survey$bin$dx),bin$gdf[list],
                 10^(survey$bin$xmin+seq(1,survey$bin$n)[list]*survey$bin$dx),bin$gdf[list],col=col.data,lwd=lwd.data)
      } else {
        points(bin$xmean[list],bin$gdf[list],pch=20,col=col.data,cex=cex.data)
        segments(bin$xmean[list],bin$gdf[list]*f.16[list],bin$xmean[list],bin$gdf[list]*f.84[list],col=col.data,lwd=lwd.data)
        points(bin$xmean[upper],bin$gdf[upper]*pm,pch=25,cex=cex.data,bg=col.data)
        segments((survey$bin$xmin+seq(0,survey$bin$n-1)[list]*survey$bin$dx),bin$gdf[list],
                 (survey$bin$xmin+seq(1,survey$bin$n)[list]*survey$bin$dx),bin$gdf[list],col=col.data,lwd=lwd.data)
      }
    }
  }
  
  # axes
  if (!add) {
    magicaxis::magaxis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1)
  }
  
  # plot binned data histogram
  if (!is.na(show.data.histogram) & show.data.histogram) {
    .plotSub(0,1,0,0.2)
    ymax = max(survey$bin$histogram)*1.2
    if (length(grep('x',log))==1) {lg='x'} else {lg=''}
    plot(1,1,type='n',log=lg,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
         xlim = xlim, ylim = c(0,ymax), xlab = '', ylab = '',bty='n')
    xbin = rep(survey$bin$xmin+seq(0,survey$bin$n)*survey$bin$dx,each=2)
    if (xpower10) xbin = 10^xbin
    xhist = c(xlim[1],xbin,xlim[2])
    yhist = c(0,0,rep(survey$bin$histogram,each=2),0,0)
    polygon(xhist,yhist,col=col.hist,border = NA)
    if (!is.null(veff)) {
      if (xpower10) {
        x = seq(log10(xlim[1]),log10(xlim[2]),length=200)
      } else {
        x = seq(xlim[1],xlim[2],length=200)
      }
      y = veff(x)
      if (xpower10) {x = 10^x}
      lines(x,y/max(y)*ymax*0.85,col=col.veff,lwd=lwd.veff,lty=lty.veff)
    }
    par(xpd=TRUE)
    lines(xlim,rep(ymax,2))
    par(xpd=FALSE)
    magicaxis::magaxis(side=2,ylab=ylab.histogram,lwd=NA,labels=FALSE,lwd.ticks=NA)
    magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=NA)
  }
  
  if (is.na(show.data.histogram) | show.data.histogram) {
    .plotSubEnd()
  }

  # axes
  if (!add) {
    magicaxis::magaxis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
    magicaxis::magaxis(side=3,labels=FALSE,lwd=NA,lwd.ticks=1)
    box(which = "plot", lty = "solid", lwd = 1)
    par(pty = "m")
  }
  
  if (is.na(show.data.histogram) | show.data.histogram) {
    .plotSub(0,1,0,1)
    plot(1,1,type='n',log=log,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
         xlim = xlim, ylim = c(exp(log(ylim[1])-(log(ylim[2])-log(ylim[1]))*0.25),ylim[2]), xlab = '', ylab = '',bty='n')
    clip(xlim[1],xlim[2],ylim[1],ylim[2])
  }
  
  invisible(survey)
}

.bin.data = function(survey,nbins,bin.xmin,bin.xmax) {
  
  # initialize
  x = survey$data$x
  bin = list()
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
  
  # fill input data into bins
  v = survey$selection$veff(x)
  bin$gdf.input = bin$count.input = bin$xmean.input = array(0,bin$n)
  for (i in seq(n.data)) {
    k = floor((x[i]-bin$xmin)/wx*0.99999999*bin$n)+1
    if (k>=1 & k<=bin$n) {
      bin$gdf.input[k] = bin$gdf.input[k]+1/bin$dx/v[i]
      bin$count.input[k] = bin$count.input[k]+1
      bin$xmean.input[k] = bin$xmean.input[k]+x[i]
    }
  }
  bin$xmean.input = bin$xmean.input/pmax(1,bin$count.input)
    
  # fill posterior data into bins
  if (!is.null(survey$data$x.err) & !is.null(survey$grid$effective.counts)) {
    bin$gdf.posterior = bin$count.posterior = bin$xmean.posterior = array(0,bin$n)
    xg = survey$grid$x
    dx = xg[2]-xg[1]
    for (k in seq(bin$n)) {
      list = floor((xg-bin$xmin)/wx*0.99999999*bin$n)+1==k
      bin$xmean.posterior[k] = sum(survey$grid$scd.posterior[list]*xg[list])/sum(survey$grid$scd.posterior[list])
      bin$count.posterior[k] = mean(survey$grid$effective.counts[list])
      bin$gdf.posterior[k] = sum(survey$grid$scd.posterior[list]/survey$grid$veff[list])/sum(list)
    }
  }
  
  # make historgram counts
  bin$histogram = array(0,bin$n)
  for (i in seq(n.data)) {
    k = floor((x[i]-bin$xmin)/wx*0.99999999*bin$n)+1
    if (k>=1 & k<=bin$n) {
      bin$histogram[k] = bin$histogram[k]+1
    }
  }
  survey$bin = bin
  
  invisible(survey)
}

#' @export
.plotSub <- function(xleft=0.1,xright=0.3,ybottom=0.1,ytop=0.3) {
  par(omd=c(0,1,0,1))
  xmarg = sum(par()$mai[c(2,4)])
  xplot = par()$pin[1]
  ymarg = sum(par()$mai[c(1,3)])
  yplot = par()$pin[2]
  par(new=T,omd=c(xleft*xplot/(xplot+xmarg),(xright*xplot+xmarg)/(xplot+xmarg),ybottom*yplot/(yplot+ymarg),(ytop*yplot+ymarg)/(yplot+ymarg)))
}

#' @export
.plotSubEnd <- function() {
  par(oma=c(0,0,0,0))
  par(omi=c(0,0,0,0))
  par(omd=c(0,1,0,1))
}

#' @export
.dfsun <- function(x,y,cex=1) {
  par(xpd=TRUE)
  points(rep(x,2),rep(y,2),pch=c(1,20),cex=c(1.2,0.45)*cex)
  par(xpd=FALSE)
}

