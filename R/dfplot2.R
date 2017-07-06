#' Display fitted distribution function
#'
#' This function displays a one-dimensional distribution function fitted using \code{\link{dffit}}.
#'
#' @importFrom magicaxis magaxis magplot
#' @importFrom graphics contour box rasterImage
#'
#' @param survey List produced by \code{\link{dffit}}
#' @param xlab,ylab axis labels
#' @param xlim,ylim 2-element vectors with plotting ranges
#' @param p.ref Optional vector of parameters of a reference mdeol to be displayed.
#' @param xpower10 If \code{TRUE}, the 1st model argument is elevated to the power of 10 in the plots.
#' @param ypower10 If \code{TRUE}, the 2nd model argument is elevated to the power of 10 in the plots.
#' @param show.data Logical flag to show/hide the observed data.
#' @param show.data.err Logical flag to show/hide the measurement errors of the
#' @param show.gdf Logical flag to show/hide the best fitting model of the generative distribution function as gradual shading.
#' @param show.gdf.contours Logical flag to show/hide the best fitting model of the generative distribution function as iso-contours.
#' @param show.scd.contours Logical flag to show/hide the source count density (for data with no measurement errors) predicted by the best fitting model.
#' @param contour.levels Vector specifying the fractions of the integrated densities contained inside the contours of the gdf and scd.
#' @param col.data Color (see \code{\link{plot}}) of data points.
#' @param cex.data Size (see \code{\link{plot}}) of data points.
#' @param pch.data Point type (see \code{\link{plot}}) of data points.
#' @param col.data.err Color (see \code{\link{plot}}) of error ellipses of data.
#' @param lwd.data.err Line width (see \code{\link{plot}}) of error ellipses of data.
#' @param lty.data.err Line type (see \code{\link{plot}}) of error ellipses of data.
#' @param col.gdf Color (see \code{\link{plot}}) of best fitting generative distribution function.
#' @param col.gdf.contours Color (see \code{\link{plot}}) of best fitting generative distribution function contours.
#' @param lwd.gdf.contours Line width (see \code{\link{plot}}) of best fitting generative distribution function contours.
#' @param lty.gdf.contours Line type (see \code{\link{plot}}) of best fitting generative distribution function contours.
#' @param gamma.gdf Scalar value >0 specifying the brightness-scale of the best fitting generative distribution function.
#' @param col.scd.contours Color (see \code{\link{plot}}) of the source count density contours predicted by the best fitting model.
#' @param lwd.scd.contours Line width (see \code{\link{plot}}) of the source count density contours predicted by the best fitting model.
#' @param lty.scd.contours Line type (see \code{\link{plot}}) of the source count density contours predicted by the best fitting model.
#' @param col.ref.contours Color (see \code{\link{plot}}) of the reference model contours specified by \code{p.ref}.
#' @param lwd.ref.contours Line width (see \code{\link{plot}}) of the reference model contours specified by \code{p.ref}.
#' @param lty.ref.contours Line type (see \code{\link{plot}}) of the reference model contours specified by \code{p.ref}.
#' @param margins Margins (bottom,left,top,right)
#' 
#' @return Returns the input list \code{survey}.
#' 
#' @seealso For an full example of \code{dfplot2} run \code{dfexample(4)}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplot2 <- function(survey,
                   xlab = 'Observable x1',
                   ylab = 'Observable x2',
                   xlim = NULL,
                   ylim = NULL,
                   p.ref = NULL,
                   xpower10 = TRUE,
                   ypower10 = TRUE,
                   show.data = TRUE,
                   show.data.err = TRUE,
                   show.gdf = TRUE,
                   show.gdf.contours = TRUE,
                   show.scd.contours = TRUE,
                   contour.levels = c(0.5,0.2,0.1),
                   col.data = 'black',
                   cex.data = 0.5,
                   pch.data = 20,
                   col.data.err = '#00000030',
                   lwd.data.err = 1,
                   lty.data.err = 1,
                   col.gdf = 'blue',
                   col.gdf.contours = col.gdf,
                   lwd.gdf.contours = 1,
                   lty.gdf.contours = 2,
                   gamma.gdf = 0.2,
                   col.scd.contours = 'purple',
                   lwd.scd.contours = 1,
                   lty.scd.contours = 1,
                   col.ref.contours = 'red',
                   lwd.ref.contours = 1,
                   lty.ref.contours = 2,
                   margins=c(5.1,4.1,4.1,2.1)) {
  
  n.data = dim(survey$data$x)[1]
  n.dim = dim(survey$data$x)[2]
  
  if (n.dim!=2) {
    if (n.dim==1) {
      stop('Use dfplot for one-dimensional distribution functions.')
    } else {
      stop('dfplot2 only handles two-dimensional distribution functions.')
    }
  }
  
  # define plot limits
  x.grid = list(seq(survey$grid$xmin[1],survey$grid$xmax[1],survey$grid$dx[1]),
                seq(survey$grid$xmin[2],survey$grid$xmax[2],survey$grid$dx[2]))
  xrange = c(survey$grid$xmin[1],survey$grid$xmax[1])
  yrange = c(survey$grid$xmin[2],survey$grid$xmax[2])
  if (is.null(xlim)) {
    if (xpower10) {
      xlim = 10^xrange
    } else {
      xlim = xrange
    }
  }
  if (is.null(ylim)) {
    if (ypower10) {
      ylim = 10^yrange
    } else {
      ylim = yrange
    }
  }
  log = ''
  if (xpower10) {
    log = paste0(log,'x')
    x = 10^survey$data$x[,1]
    cx = 10^x.grid[[1]]
  } else {
    x = survey$data$x[,1]
    cx = x.grid[[1]]
  }
  if (ypower10) {
    log = paste0(log,'y')
    y = 10^survey$data$x[,2]
    cy = 10^x.grid[[2]]
  } else {
    y = survey$data$x[,2]
    cy = x.grid[[2]]
  }
  
  # make model DF field
  r = col2rgb(col.gdf)[1]/255
  g = col2rgb(col.gdf)[2]/255
  b = col2rgb(col.gdf)[3]/255
  nx = length(x.grid[[1]])
  ny = length(x.grid[[2]])
  img = array(survey$grid$gdf,c(nx,ny))
  img = img/max(img)
  rgb = array(NA,c(nx,ny,3))
  rgb[,,1] = 1-img^gamma.gdf*0.7*(1-r)
  rgb[,,2] = 1-img^gamma.gdf*0.7*(1-g)
  rgb[,,3] = 1-img^gamma.gdf*(1-b)
  
  # make model source count field
  imgsc = array(survey$grid$scd,c(nx,ny))
  imgsc = imgsc/max(imgsc)
  
  # make reference DF field
  if (!is.null(p.ref)) {
    yref = survey$model$gdf(survey$grid$x,p.ref)
    imgref = array(yref,c(nx,ny))
    imgref = imgref/max(imgref)
  }
  
  # open plot
  par(pty = 'm')
  par(mar = margins)
  plot(1,1,type='n',log=log,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
         xlim = xlim, ylim = ylim, xlab = '', ylab = '',bty='n')
  
  # plot central fit
  if (show.gdf) {
    rasterImage(aperm(rgb[,ny:1,],c(2,1,3)),xlim[1],ylim[1],xlim[2],ylim[2])
  }
  if (show.gdf.contours) {
    contour(cx,cy,img,add=TRUE,lwd=lwd.gdf.contours,lty=lty.gdf.contours,col=col.gdf.contours,drawlabels=FALSE,levels=contour.levels)
  }
  
  # plot reference
  if (!is.null(p.ref)) {
    contour(cx,cy,imgref,add=TRUE,lwd=lwd.ref.contours,lty=lty.ref.contours,col=col.ref.contours,drawlabels=FALSE,levels=contour.levels)
  }
  
  # plot source counts
  if (show.scd.contours) {
    contour(cx,cy,imgsc,add=TRUE,lwd=lwd.scd.contours,lty=lty.scd.contours,col=col.scd.contours,drawlabels=FALSE,levels=contour.levels)
  }
  
  # plot data points
  if (show.data) {
    if (show.data.err) {
      if (length(dim(survey$data$x.err))==2) {
        for (i in seq(n.data)) {
          px = survey$data$x[i,1]+c(-1,+1)*survey$data$x.err[i,1]
          py = survey$data$x[i,2]*c(1,1)
          if (xpower10) {px = 10^px}
          if (ypower10) {py = 10^py}
          lines(px,py,col=col.data.err,lwd=lwd.data.err,lty=lty.data.err)
          px = survey$data$x[i,1]*c(1,1)
          py = survey$data$x[i,2]+c(-1,+1)*survey$data$x.err[i,2]
          if (xpower10) {px = 10^px}
          if (ypower10) {py = 10^py}
          lines(px,py,col=col.data.err,lwd=lwd.data.err,lty=lty.data.err)
        }
      } else if (length(dim(survey$data$x.err))==3) {
        for (i in seq(n.data)) {
          pts = ellipse::ellipse(survey$data$x.err[i,,],centre=survey$data$x[i,],level=0.68,draw=F)
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
  
  invisible(survey)
}
