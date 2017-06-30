#' Plot parameter covariance matrix
#'
#' This function displays the parameter covariances  of the distribution function parameters (e.g. the mass function parameters) fitted using \code{\link{dffit}}.
#'
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm rmvnorm
#' @importFrom graphics par rect lines points plot axis hist
#'
#' @param survey List produced by \code{\link{dffit}}
#' @param expectation Vector of expected parameters
#' @param covariance Matrix of parameter covariances
#' @param expectation2 Vector of expected parameters of an optional 2nd model
#' @param covariance2 Matrix of parameter covariances of an optional 2nd model
#' @param chain Matrix, where each row represents a set of model parameters, to be displayed as a point of a cloud.
#' @param reference Reference parameters to be plotted as points and crosses.
#' @param name Optional list of parameter names
#' @param title Optional plot title
#' @param nstd Width of the plots in multiples of the standard deviations of the model, i.e. the square roots of the diagonal elements of \code{covariance}.
#' @param lower Logical flag indicating whether the lower triangle is shown.
#' @param upper Logical flag indicating whether the lower triangle is shown.
#' @param margins Plot margins (bottom,left,top,right)
#' @param model.col Model color
#' @param model.lwd Model line width
#' @param model.lty Model line type
#' @param model.cex Model point size
#' @param model2.col Second model color
#' @param model2.lwd Second model line width
#' @param model2.lty Second model line type
#' @param model2.cex Second model point size
#' @param chain.col Color of point cloud
#' @param chain.lwd Line width of point cloud ellipses
#' @param chain.lty Line type of point cloud ellipses
#' @param chain.cex Point size at the center of point cloud
#' @param chain.point.cex Point size of cloud points
#' @param chain.point.alpha Transparency of cloud points
#' @param reference.col Color of reference point
#' @param reference.lwd Line width of reference lines
#' @param reference.lty Line type of reference lines
#' @param reference.cex Size of reference point
#' @param text.size.labels Text size of parameter names
#' @param text.size.numbers Text size of numbers
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplotcov <- function(survey = NULL,
                      expectation = NULL, covariance = NULL, # model
                      expectation2 = NULL, covariance2 = NULL, # model 2
                      chain = NULL, # chain
                      reference = NULL, # reference
                      name = NULL, title = '',
                      nstd = 10,
                      lower = TRUE, upper = FALSE,
                      margins = c(4,4,0.5,0.5),
                      model.col = 'blue', model.lwd = 2, model.lty = 1, model.cex = 1,
                      model2.col = 'black', model2.lwd = 2, model2.lty = 1, model2.cex = 1,
                      chain.col = 'purple', chain.lwd = 2, chain.lty = 1, chain.cex = 1, chain.point.cex = 0.1,
                      chain.point.alpha = 0.2,
                      reference.col = 'black', reference.lwd = 1, reference.lty = 2, reference.cex = 1,
                      text.size.labels = 1.1, text.size.numbers = 0.8) {
  
  # input handling
  if (is.null(survey)==is.null(expectation)) stop('Either the survey order expectation must be given, but not both.')
  if (is.null(expectation)!=is.null(covariance)) stop('If the expectation is provided, the covariance must be provided, too, and vice versa.')
  if (is.null(expectation2)!=is.null(covariance2)) stop('If the expectation2 is provided, the covariance2 must be provided, too, and vice versa.')
  if (!is.null(survey)) {
    expectation = survey$fit$p.best
    covariance = survey$fit$p.covariance
  }
  n = length(expectation)
  chain.r = col2rgb(chain.col)[1]/255
  chain.g = col2rgb(chain.col)[2]/255
  chain.b = col2rgb(chain.col)[3]/255
  chain.hist.col = rgb(chain.r,chain.g,chain.b,0.4)
  chain.point.col = rgb(chain.r,chain.g,chain.b,chain.point.alpha)
  
  my.hist = function(x,breaks) {
    b = c(-1e99,breaks,1e99)
    counts = hist(x,breaks=b,plot=F)$counts[2:(length(b)-2)]
    center = (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2
    x = array(rbind(array(breaks),array(breaks)))
    y = c(0,array(rbind(array(counts),array(counts))),0)
    return(list(counts=counts,center=center,x=x,y=y))
  }

  # define function for single squares
  make_single_square <- function(i,j) {

    xlength = sqrt(covariance[i,i])*nstd
    ylength = sqrt(covariance[j,j])*nstd

    xmin = expectation[i]-xlength/2
    xmax = expectation[i]+xlength/2
    ymin = expectation[j]-ylength/2
    ymax = expectation[j]+ylength/2

    E = expectation[c(i,j)] # 2-element vector with the center position
    C = covariance[c(i,j),c(i,j)] # 2-by-2 covariance matrix

    xoffset = i-1 # x-position of bottom left corner of the little square in the coordinates of the whole plot
    yoffset = n-j # y-position of bottom left corner of the little square in the coordinates of the whole plot

    # draw rectangle
    par(xpd=T)
    rect(xleft=xoffset,xright=xoffset+1,ybottom=yoffset,ytop=yoffset+1,lwd=1.5)
    par(xpd=F)

    if (i==j) {
      
      # chain (transparent histogram)
      if (!is.null(chain)) {
        nbins = min(100,max(1,sqrt(dim(chain)[1])))
        h = my.hist(chain[,i],seq(xmin,xmax,length=nbins+1))
        xplot = rep(seq(0,1,length=nbins+1),each=2)+xoffset
        yplot = h$y/max(h$y)*0.9+yoffset
        polygon(xplot,yplot,col=chain.hist.col,border=NA)
      }
      
      # model2 (bell curve)
      if (!is.null(expectation2)) {
        x = seq(xmin,xmax,length=200)
        y = 0.9*sqrt(C[1]/covariance2[i,i])*exp(-(x-expectation2[i])^2/2/covariance2[i,i])
        x = (x-xmin)/xlength+xoffset
        y = pmin(1,y)
        y[2:199] = y[1:198]*0.0005+y[2:199]*0.999000001+y[3:200]*0.0005
        y[y>=1] = NA
        lines(x,y+yoffset,col=model2.col,lty=model2.lty,lwd=model2.lwd)
      }

      # model (bell curve)
      sigma = 1/nstd
      x = seq(0,1,length=200)
      y = 0.9*exp(-(x-0.5)^2/2/sigma^2)
      lines(x+xoffset,y+yoffset,col=model.col,lty=model.lty,lwd=model.lwd)
      
      # refernce (vertical line)
      if (!is.null(reference)) {
        xplot = (reference[i]-E[1])/nstd/sqrt(C[1,1])+0.5
        if (xplot>0 & xplot<1) lines(rep(xplot+xoffset,2),c(0,1)+yoffset,lty=reference.lty,lwd=reference.lwd,col=reference.col)
      }
      
    } else {
      
      # chain (point cloud)
      if (!is.null(chain)) {
        selection = (chain[,i]>xmin) & (chain[,i]<xmax) & (chain[,j]>ymin) & (chain[,j]<ymax)
        xplot = (chain[selection,i]-xmin)/xlength+xoffset
        yplot = (chain[selection,j]-ymin)/ylength+yoffset
        points(xplot,yplot,pch=20,cex=chain.point.cex,col=chain.point.col)
      }
      
      # chain (ellipse and point)
      if (!is.null(chain)) {
        xmu = mean(chain[,i]) 
        ymu = mean(chain[,j])
        cxy = array(NA,c(2,2))
        cxy[1,1] = cov(chain[,i],chain[,i])
        cxy[1,2] = cxy[2,1] = cov(chain[,i],chain[,j])
        cxy[2,2] = cov(chain[,j],chain[,j])
        pts = ellipse::ellipse(cxy,centre=c(xmu,ymu),level=0.68,draw=F)
        xplot = (pts[,1]-xmin)/xlength # convert the coordinates into the coordinates of the plot
        yplot = (pts[,2]-ymin)/ylength
        lines(pmin(1,pmax(0,xplot))+xoffset,pmin(1,pmax(0,yplot))+yoffset,col=chain.col,lty=chain.lty,lwd=chain.lwd)
        pts = ellipse::ellipse(cxy,centre=c(xmu,ymu),level=0.95,draw=F)
        xplot = (pts[,1]-xmin)/xlength # convert the coordinates into the coordinates of the plot
        yplot = (pts[,2]-ymin)/ylength
        lines(pmin(1,pmax(0,xplot))+xoffset,pmin(1,pmax(0,yplot))+yoffset,col=chain.col,lty=chain.lty,lwd=chain.lwd/2)
        points((xmu-xmin)/xlength+xoffset,(ymu-ymin)/ylength+yoffset,pch=20,col=chain.col,cex=chain.cex)
      }
      
      # model2 (ellipse and point)
      if (!is.null(expectation2)) {
        pts = ellipse::ellipse(covariance2[c(i,j),c(i,j)],centre=expectation2[c(i,j)],level=0.68,draw=F)
        xplot = (pts[,1]-xmin)/xlength # convert the coordinates into the coordinates of the plot
        yplot = (pts[,2]-ymin)/ylength
        lines(pmin(1,pmax(0,xplot))+xoffset,pmin(1,pmax(0,yplot))+yoffset,col=model2.col,lty=model2.lty,lwd=model2.lwd)
        pts = ellipse::ellipse(covariance2[c(i,j),c(i,j)],centre=expectation2[c(i,j)],level=0.95,draw=F)
        xplot = (pts[,1]-xmin)/xlength # convert the coordinates into the coordinates of the plot
        yplot = (pts[,2]-ymin)/ylength
        lines(pmin(1,pmax(0,xplot))+xoffset,pmin(1,pmax(0,yplot))+yoffset,col=model2.col,lty=model2.lty,lwd=model2.lwd/2)
        points((expectation2[i]-xmin)/xlength+xoffset,(expectation2[j]-ymin)/ylength+yoffset,pch=20,col=model2.col,cex=model2.cex)
      }

      # model (ellipse and point)
      pts = ellipse::ellipse(C,centre=E,level=0.68,draw=F)
      xplot = (pts[,1]-xmin)/xlength # convert the coordinates into the coordinates of the plot
      yplot = (pts[,2]-ymin)/ylength
      lines(pmin(1,pmax(0,xplot))+xoffset,pmin(1,pmax(0,yplot))+yoffset,col=model.col,lty=model.lty,lwd=model.lwd)
      pts = ellipse::ellipse(C,centre=E,level=0.95,draw=F)
      xplot = (pts[,1]-xmin)/xlength # convert the coordinates into the coordinates of the plot
      yplot = (pts[,2]-ymin)/ylength
      lines(pmin(1,pmax(0,xplot))+xoffset,pmin(1,pmax(0,yplot))+yoffset,col=model.col,lty=model.lty,lwd=model.lwd/2)
      points(0.5+xoffset,0.5+yoffset,pch=20,col=model.col,cex=model.cex)

      # reference (lines and point)
      if (!is.null(reference)) {
        xplot = (reference[i]-E[1])/nstd/sqrt(C[1,1])+0.5
        yplot = (reference[j]-E[2])/nstd/sqrt(C[2,2])+0.5
        if (yplot>0 & yplot<1) lines(c(0,1)+xoffset,rep(yplot+yoffset,2),lwd=reference.lwd,lty=reference.lty,col=reference.col)
        if (xplot>0 & xplot<1) lines(rep(xplot+xoffset,2),c(0,1)+yoffset,lwd=reference.lwd,lty=reference.lty,col=reference.col)
        if (xplot>0 & xplot<1 & yplot>0 & yplot<1) points(xplot+xoffset,yplot+yoffset,pch=20,col=reference.col,cex=reference.cex)
      }

    }

    # add axes
    tickpos = c(0.2,0.5,0.8)
    if (i==1 & j>1) {
      axis(2, at = yoffset+tickpos,labels=sprintf('%4.1f',E[2]+ylength*(tickpos-0.5)),tck=0.015,lwd=0,lwd.ticks=1,cex.axis=text.size.numbers,padj = 0.8)
      if (is.null(name)) {
        axis(2, at = yoffset+0.5,pos=-0.1,labels=bquote(p [.(j)]),cex.axis=text.size.labels,tick=F)
      } else {
        axis(2, at = yoffset+0.5,pos=-0.1,labels=name[j],cex.axis=text.size.labels,tick=F)
      }
    }
    if (j==n) {
      axis(1, at = xoffset+tickpos,labels=sprintf('%4.1f',E[1]+xlength*(tickpos-0.5)),tck=0.015,lwd=0,lwd.ticks=1,cex.axis=text.size.numbers,padj = -0.8)
      if (is.null(name)) {
        axis(1, at = xoffset+0.5,pos=-0.1,labels=bquote(p [.(i)]),cex.axis=text.size.labels,tick=F)
      } else {
        axis(1, at = xoffset+0.5,pos=-0.1,labels=name[i],cex.axis=text.size.labels,tick=F)
      }
    }
  }

  # start a new plot
  pty.old = par()$pty
  mar.old = par()$mar
  par(pty='s')
  par(mar=margins)
  plot(0,0,type='n',xlim=c(0,n),ylim=c(0,n),ann=FALSE,xaxs='i',yaxs='i',xaxt='n',yaxt='n',bty='n')

  # iterate over all parameter-pairs
  for (i in seq(n)) {
    for (j in seq(n)) {
      if (i==j | (j>i & lower) | (i>j & upper)) {
        make_single_square(i,j)
      }
    }
  }

  text(1.2,n-0.2,title,offset=0,adj=0,cex=1)
  par(pty=pty.old)
  par(mar=mar.old)

}
