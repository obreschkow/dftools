#' Display fitted distribution function
#'
#' This function displays a one-dimensional distribution function fitted using \code{\link{dffit}}.
#'
#' @importFrom magicaxis magaxis magplot
#'
#' @param df List produced by \code{\link{dffit}}
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param xlim 2-element vector with x-axis plotting limits
#' @param ylim 2-element vector with y-axis plotting limits
#' @param xpower10 If \code{TRUE}, the model argument x is elevated to the power of 10 in the plots.
#' @param ypower10 If \code{TRUE}, the model argument x is elevated to the power of 10 in the plots.
#' @param col Color of fit (see \code{\link{plot}})
#' @param cex Line width of fit (see \code{\link{plot}})
#' @param pch Line type of fit (see \code{\link{plot}})
#' @param show.data If \code{TRUE}, binned data points are displayed.
#' @param show.data.err If \code{TRUE}, the bias corrected MLE is shown instead of the native ML parameters.
#' @param margins Margins (bottom,left,top,right)
#' 
#' @return Returns the input list \code{df} with the additional sub-list \code{bundle$bin} that contains the binned data.
#' 
#' @seealso For optimized plotting of galaxy mass functions, use the derived function \code{\link{mfplot}}. See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplot2 <- function(bundle,
                   xlab = 'Observable x1',
                   ylab = 'Observable x2',
                   xlim = NULL,
                   ylim = NULL,
                   p.ref = NULL,
                   xpower10 = TRUE,
                   ypower10 = TRUE,
                   show.data = TRUE,
                   show.data.err = TRUE,
                   show.phi = TRUE,
                   show.phi.contours = TRUE,
                   show.sc.contours = TRUE,
                   contour.levels = c(0.5,0.2,0.1),
                   gamma = 0.2,
                   col.data = 'black',
                   cex.data = 0.5,
                   pch.data = 20,
                   col.data.err = '#00000030',
                   lwd.data.err = 1,
                   lty.data.err = 1,
                   col.phi = 'blue',
                   col.phi.contours = col.phi,
                   lwd.phi.contours = 1,
                   lty.phi.contours = 2,
                   col.sc.contours = 'purple',
                   lwd.sc.contours = 1,
                   lty.sc.contours = 1,
                   col.ref.contours = 'red',
                   lwd.ref.contours = 1,
                   lty.ref.contours = 2,
                   margins=c(5.1,4.1,4.1,2.1)) {
  
  n.data = dim(bundle$data$x)[1]
  n.dim = dim(bundle$data$x)[2]
  
  if (n.dim!=2) {
    if (n.dim==1) {
      stop('Use dfplot for one-dimensional distribution functions.')
    } else {
      stop('dfplot2 only handles two-dimensional distribution functions.')
    }
  }
  
  # define plot limits
  x.grid = list(seq(bundle$grid$xmin[1],bundle$grid$xmax[1],bundle$grid$dx[1]),
                seq(bundle$grid$xmin[2],bundle$grid$xmax[2],bundle$grid$dx[2]))
  xrange = c(bundle$grid$xmin[1],bundle$grid$xmax[1])
  yrange = c(bundle$grid$xmin[2],bundle$grid$xmax[2])
  log = ''
  if (is.null(xlim)) {
    if (xpower10) {
      xlim = 10^xrange
      log = paste0(log,'x')
    } else {
      xlim = xrange
    }
  }
  if (is.null(ylim)) {
    if (ypower10) {
      ylim = 10^yrange
      log = paste0(log,'y')
    } else {
      ylim = yrange
    }
  }
  if (xpower10) {
    x = 10^bundle$data$x[,1]
    cx = 10^x.grid[[1]]
  } else {
    x = bundle$data$x[,1]
    cx = x.grid[[1]]
  }
  if (ypower10) {
    y = 10^bundle$data$x[,2]
    cy = 10^x.grid[[2]]
  } else {
    y = bundle$data$x[,2]
    cy = x.grid[[2]]
  }
  
  # make model DF field
  r = col2rgb(col.phi)[1]/255
  g = col2rgb(col.phi)[2]/255
  b = col2rgb(col.phi)[3]/255
  nx = length(x.grid[[1]])
  ny = length(x.grid[[2]])
  img = array(bundle$grid$gdf,c(nx,ny))
  img = img/max(img)
  rgb = array(NA,c(nx,ny,3))
  rgb[,,1] = 1-img^gamma*0.7*(1-r)
  rgb[,,2] = 1-img^gamma*0.7*(1-g)
  rgb[,,3] = 1-img^gamma*(1-b)
  
  # make model source count field
  imgsc = array(bundle$grid$scd,c(nx,ny))
  imgsc = imgsc/max(imgsc)
  
  # make reference DF field
  if (!is.null(p.ref)) {
    yref = bundle$model$gdf(bundle$grid$x,p.ref)
    imgref = array(yref,c(nx,ny))
    imgref = imgref/max(imgref)
  }
  
  # open plot
  par(pty = 'm')
  par(mar = margins)
  plot(1,1,type='n',log=log,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
         xlim = xlim, ylim = ylim, xlab = '', ylab = '',bty='n')
  
  # plot central fit
  if (show.phi) {
    rasterImage(aperm(rgb[,ny:1,],c(2,1,3)),xlim[1],ylim[1],xlim[2],ylim[2])
  }
  if (show.phi.contours) {
    contour(cx,cy,img,add=TRUE,lwd=lwd.phi.contours,lty=lty.phi.contours,col=col.phi.contours,drawlabels=FALSE,levels=contour.levels)
  }
  
  # plot reference
  if (!is.null(p.ref)) {
    contour(cx,cy,imgref,add=TRUE,lwd=lwd.ref.contours,lty=lty.ref.contours,col=col.ref.contours,drawlabels=FALSE,levels=contour.levels)
  }
  
  # plot source counts
  if (show.sc.contours) {
    contour(cx,cy,imgsc,add=TRUE,lwd=lwd.sc.contours,lty=lty.sc.contours,col=col.sc.contours,drawlabels=FALSE,levels=contour.levels)
  }
  
  # plot data points
  if (show.data) {
    if (show.data.err) {
      if (length(dim(bundle$data$x.err))==2) {
        for (i in seq(n.data)) {
          px = bundle$data$x[i,1]+c(-1,+1)*bundle$data$x.err[i,1]
          py = bundle$data$x[i,2]*c(1,1)
          if (xpower10) {px = 10^px}
          if (ypower10) {py = 10^py}
          lines(px,py,col=col.data.err,lwd=lwd.data.err,lty=lty.data.err)
          px = bundle$data$x[i,1]*c(1,1)
          py = bundle$data$x[i,2]+c(-1,+1)*bundle$data$x.err[i,2]
          if (xpower10) {px = 10^px}
          if (ypower10) {py = 10^py}
          lines(px,py,col=col.data.err,lwd=lwd.data.err,lty=lty.data.err)
        }
      } else if (length(dim(bundle$data$x.err))==3) {
        for (i in seq(n.data)) {
          pts = ellipse::ellipse(bundle$data$x.err[i,,],centre=bundle$data$x[i,],level=0.68,draw=F)
          if (xpower10) {px = 10^pts[,1]} else {px = pts[,1]}
          if (ypower10) {py = 10^pts[,2]} else {py = pts[,2]}
          lines(px,py,col=col.data.err,lwd=lwd.data.err,lty=lty.data.err)
        }
      }
    }
    points(x,y,pch=pch.data,cex=cex.data,col=col.data)
  }
  
  # axes
  magicaxis::magaxis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=3,labels=FALSE,lwd=NA,lwd.ticks=1)
  magicaxis::magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1)
  box(which = "plot", lty = "solid", lwd = 1)
  par(pty = "m")
  
  invisible(bundle)
}
