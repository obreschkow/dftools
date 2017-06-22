#' Plot parameter covariance matrix
#'
#' This function displays the parameter covariances  of the distribution function parameters (e.g. the mass function parameters) fitted using \code{\link{dffit}}.
#'
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm rmvnorm
#' @importFrom graphics par rect lines points plot axis
#'
#' @param survey List produced by \code{\link{dffit}}
#' @param p Reference parameters to be plotted as crosses.
#' @param n.points Number of random points drawn from the multivariate Gaussian.
#' @param n.standarddev Number of standard-deviations shown to the left and right of the mean.
#' @param line.color Line color
#' @param line.lwd Line width
#' @param point.color Point color
#' @param point.cex Point size
#' @param text.size.labels Text size of parameter names
#' @param text.size.numbers Text size of parameter values
#' @param margins Margins (bottom,left,top,right)
#' @param title Title text
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplotcov <- function(survey, p = NULL, n.points = 500, n.standarddev = 5,
                            line.color = 'blue', line.lwd = 2,
                            point.color = '#bbbbff', point.cex = 0.1,
                            text.size.labels = 1.1, text.size.numbers = 0.8,
                            margins = c(4,4,0.5,0.5),
                            title = '') {
  covariance = survey$fit$p.covariance
  .covariance.plot(survey$fit$p.best,
                  survey$fit$p.covariance,
                  p = p, npoints = n.points,
                  nstd = n.standarddev*2,
                  line.color = line.color, line.lwd = line.lwd,
                  point.color = point.color, point.cex = point.cex,
                  text.size.labels = text.size.labels, text.size.numbers = text.size.numbers,
                  margins = margins, title = title)
}

#' @export
.covariance.plot = function(expectation,covariance,p=NULL,
                           name=NULL,
                   npoints=300,nstd=6,
                   line.color='blue',line.lwd=2,
                   point.color='grey',point.cex=0.1,
                   text.size.labels=1.1,
                   text.size.numbers=0.8,
                   margins=c(4,4,0.5,0.5),
                   title = '') {

  # expectation : vector of n elements with the central values of n parameters
  # covariance  : n-by-n covariance matrix
  # name        : list of n parameter names
  # npoints     : number of points per square
  # nstd        : side-length of the squares in standard deviations

  n = length(expectation)

  make_single_square <- function(i,j) {

    xlength = sqrt(covariance[i,i])*nstd
    ylength = sqrt(covariance[j,j])*nstd

    xmin    = expectation[i]-xlength/2
    ymin    = expectation[j]-ylength/2

    E = expectation[c(i,j)] # 2-element vector with the center position
    C = covariance[c(i,j),c(i,j)] # 2-by-2 covariance matrix

    xoffset = i-1 # x-position of bottom left corner of the little square in the coordinates of the whole plot
    yoffset = n-j # y-position of bottom left corner of the little square in the coordinates of the whole plot

    # draw rectangle
    par(xpd=T)
    rect(xleft=xoffset,xright=xoffset+1,ybottom=yoffset,ytop=yoffset+1,lwd=1.5)
    par(xpd=F)

    if (i==j) {

      # draw Gaussian
      sigma = 1/nstd
      x = seq(0,1,length=200)
      y = 0.9*exp(-(x-0.5)^2/2/sigma^2)
      lines(x+xoffset,y+yoffset,col=line.color,lty=1,lwd=line.lwd)
      if (!is.null(p)) {
        xplot = (p[i]-E[1])/nstd/sqrt(C[1,1])+0.5
        if (xplot>0 & xplot<1) lines(rep(xplot+xoffset,2),c(0,1)+yoffset,lty=2)
      }
    } else {

      # draw points
      if ((!is.null(npoints)) && (npoints>0)) {
        pts = mvtnorm::rmvnorm(n=npoints,mean=E,sigma=C)
        xplot = (pts[,1]-xmin)/xlength # convert the coordinates into the coordinates of the plot
        yplot = (pts[,2]-ymin)/ylength
        selection = (xplot>0) & (xplot<1) & (yplot>0) & (yplot<1)
        points(xplot[selection]+xoffset,yplot[selection]+yoffset,pch=20,cex=point.cex,col=point.color)
      }

      # draw the ellipse
      pts = ellipse::ellipse(C,centre=E,level=0.68,draw=F)
      xplot = (pts[,1]-xmin)/xlength # convert the coordinates into the coordinates of the plot
      yplot = (pts[,2]-ymin)/ylength
      lines(xplot+xoffset,yplot+yoffset,col=line.color,lty=1,lwd=line.lwd)
      pts = ellipse::ellipse(C,centre=E,level=0.95,draw=F)
      xplot = (pts[,1]-xmin)/xlength # convert the coordinates into the coordinates of the plot
      yplot = (pts[,2]-ymin)/ylength
      lines(xplot+xoffset,yplot+yoffset,col=line.color,lty=1,lwd=line.lwd/2)

      if (!is.null(p)) {
        xplot = (p[i]-E[1])/nstd/sqrt(C[1,1])+0.5
        yplot = (p[j]-E[2])/nstd/sqrt(C[2,2])+0.5
        if (yplot>0 & yplot<1) lines(c(0,1)+xoffset,rep(yplot+yoffset,2),lty=2)
        if (xplot>0 & xplot<1) lines(rep(xplot+xoffset,2),c(0,1)+yoffset,lty=2)
        if (xplot>0 & xplot<1 & yplot>0 & yplot<1) points(xplot+xoffset,yplot+yoffset,pch=20,col='black')
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
      if (j>=i) {
        make_single_square(i,j)
      }
    }
  }

  text(1.2,n-0.2,title,offset=0,adj=0,cex=1)
  par(pty=pty.old)
  par(mar=mar.old)

}
