#' Plot parameter covariance matrix
#'
#' This function displays the parameter covariances of the distribution function parameters (e.g. the mass function parameters) fitted using \code{\link{dffit}}.
#'
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm rmvnorm
#' @importFrom graphics par rect lines points plot axis hist
#' @importFrom pracma erf circshift
#'
#' @param data List of objects to be plotted. Each object in this list must be one of the following four: (1) a vector of parameters; (2) a list of two objects, a vector of parameters and an associated covariance matrix; (3) a matrix, where each row represents a set of model parameters, to be displayed as a point of a cloud; (4) an output list produced by \code{\link{dffit}}. The plot will be centered on the parameters (or average parameters) of the first object in the data list, unless the argument \code{master} is set to a different index.
#' @param master Integer specifying the index of the object in the data list defining the axes scales
#' @param order Vector of integer specifying the order of drawing.
#' @param names Optional list of parameter names
#' @param col Vector of colors of each object in the data-list
#' @param lwd Vector of line width of each object in the data-list
#' @param lty Vector of line types of each object in the data-list
#' @param cex Vector of point sizes for the central parameters each object in the data-list
#' @param pch Vector of point type for the central parameters each object in the data-list
#' @param cloud.alpha Vector of alpha values setting the transparency of the point clouds 
#' @param cloud.nmax Vector specifying the maximum number of points to be plotted in a point cloud
#' @param hist.alpha Vector of alpha values setting the transparency of the histograms
#' @param show.histogram = Logical vector specifying whether to plot histograms
#' @param show.cloud Logical vector specifying whether to plot the point cloud
#' @param show.gaussian Logical vector specifying whether to draw the Gaussian apprixmation of the parameter distribution
#' @param show.expectation Logical vector specifying whether to mark the expected values (= averages if a parameter set is specified)
#' @param show.ellipse.68 Logical vector specifying whether to draw 68\% confidence ellipses
#' @param show.ellipse.95 Logical vector specifying whether to draw 95\% confidence ellipses
#' @param title Optional plot title
#' @param nstd Width of the plots in multiples of the standard deviations of the model, i.e. the square roots of the diagonal elements of \code{covariance}.
#' @param pmin optional P-vector with lower parameter limits to be potted. If given, \code{pmax} must also be given and \code{nstd} is ignored.
#' @param pmax optional P-vector with lower parameter limits to be potted. If given, \code{pmin} must also be given and \code{nstd} is ignored.
#' @param vert.stretch P-vector specifying the vertical stretch of the diagonal panels.
#' @param nbins Optional integer specifying the number of bins used in histograms. If set to \code{NULL}, this number is determined automatically.
#' @param lower Logical flag indicating whether the lower triangle is shown
#' @param upper Logical flag indicating whether the lower triangle is shown
#' @param margins Plot margins (bottom,left,top,right)
#' @param text.size.labels Text size of parameter names
#' @param text.size.numbers Text size of numbers
#' @param text.offset.labels 2-element vector to adjust the position of the vertical and horizontal parameter names
#' @param text.offset.numbers 2-element vector to adjust the position of the vertical and horizontal numbers
#' @param text.format.numbers String specifying the floating point format of the numbers. (see \code{\link{sprintf}})
#' 
#' @examples
#' # generate two random correlated vectors p[1,] and p[2,] of n elements
#' # and display their covariance plot
#' n = 100
#' p = array(NA,c(n,2))
#' p[,1] = rnorm(n)
#' p[,2] = 0.25*p[,1]+rnorm(n,1,0.5)
#' dfplotcov(list(p))
#' 
#' # now produce the same plot, but increase the number of bins
#' # add in red color the theoretical expectation
#' expected_mean = c(0,1)
#' expected_covariance = cbind(c(1,0.5^2),c(0.5^2,0.5^2+0.25^2))
#' dfplotcov(list(p,list(expected_mean,expected_covariance)), nbins=20, col=c('black','red'))
#' 
#' # Fit a Schechter function to a mock survey, plot the best-fitting parameters
#' # with uncertainties in blue and add input parameters as black crosses
#' dat = dfmockdata(n=1000,sigma=0.5)
#' survey = dffit(dat$x, dat$veff, dat$x.err)
#' ptrue = c(-2,11,-1.3)
#' dfplotcov(list(survey,ptrue), col=c('blue','black'), pch=c(20,3))
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplotcov <- function(data,
                      master = 1,
                      order = seq(length(data)),
                      names = NULL,
                      col = c('blue','black','red','green','orange'),
                      lwd = 1,
                      lty = 1,
                      cex = 1.2,
                      pch = 20,
                      cloud.alpha = 0.1,
                      cloud.nmax = NULL,
                      hist.alpha = 0.4,
                      show.histogram = TRUE,
                      show.cloud = TRUE,
                      show.gaussian = TRUE,
                      show.expectation = TRUE,
                      show.ellipse.68 = TRUE,
                      show.ellipse.95 = TRUE,
                      title = '',
                      nstd = 10,
                      pmin = NULL,
                      pmax = NULL,
                      vert.stretch = 1,
                      nbins = NULL,
                      lower = TRUE, upper = FALSE,
                      margins = c(4,4,0.5,0.5),
                      text.size.labels = 1.1, text.size.numbers = 0.8,
                      text.offset.labels = c(0,0), text.offset.numbers = c(0,0),
                      text.format.numbers = '%4.1f') {
  
  # input handling
  m = length(data)
  if (length(col)==1) col = rep(col,m)
  if (length(lwd)==1) lwd = rep(lwd,m)
  if (length(lty)==1) lty = rep(lty,m)
  if (length(cex)==1) cex = rep(cex,m)
  if (length(pch)==1) pch = rep(pch,m)
  if (length(vert.stretch)==1) vert.stretch = rep(vert.stretch,m)
  if (length(nbins)==1) nbins = rep(nbins,m)
  if (length(cloud.alpha)==1) cloud.alpha = rep(cloud.alpha,m)
  if (length(cloud.nmax)==1) cloud.nmax = rep(cloud.nmax,m)
  if (length(hist.alpha)==1) hist.alpha = rep(hist.alpha,m)
  if (length(show.histogram)==1) show.histogram = rep(show.histogram,m)
  if (length(show.cloud)==1) show.cloud = rep(show.cloud,m)
  if (length(show.expectation)==1) show.expectation = rep(show.expectation,m)
  if (length(show.gaussian)==1) show.gaussian = rep(show.gaussian,m)
  if (length(show.ellipse.68)==1) show.ellipse.68 = rep(show.ellipse.68,m)
  if (length(show.ellipse.95)==1) show.ellipse.95 = rep(show.ellipse.95,m)
  
  my.hist = function(x,breaks) {
    b = c(-1e99,breaks,1e99)
    counts = hist(x,breaks=b,plot=F)$counts[2:(length(b)-2)]
    center = (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2
    x = array(rbind(array(breaks),array(breaks)))
    y = c(0,array(rbind(array(counts),array(counts))),0)
    return(list(counts=counts,center=center,x=x,y=y))
  }

  # define function for single squares
  make_single_square <- function(i,j,k,E,C,P) {
    
    # color handling
    cl = col2rgb(col[k])/255
    hist.col  = rgb(cl[1],cl[2],cl[3],hist.alpha[k])
    cloud.col = rgb(cl[1],cl[2],cl[3],cloud.alpha[k])
    
    # draw rectangle
    xoffset = i-1 # x-position of bottom left corner of the little square in the coordinates of the whole plot
    yoffset = n-j # y-position of bottom left corner of the little square in the coordinates of the whole plot
    if (k==1) {
      par(xpd=T)
      rect(xleft=xoffset,xright=xoffset+1,ybottom=yoffset,ytop=yoffset+1,lwd=1)
      par(xpd=F)
    }

    if (i==j) {
      
      # histogram
      if (show.histogram[k] & !is.null(P)) {
        if (is.null(nbins)) {
          nbins.plot = min(100,max(10,round(sqrt(dim(P)[1]))*nstd*0.05))
        } else {
          nbins.plot = nbins[k]
        }
        h = my.hist(P[,i],seq(xmin[i],xmax[i],length=nbins.plot+1))
        dx = h$center[2]-h$center[1]
        f = area[i]/sum(h$counts)/dx
        xplot = rep(seq(0,1,length=nbins.plot+1),each=2)+xoffset
        yplot = pmin(1,h$y*f)+yoffset
        polygon(xplot,yplot,col=hist.col,border=NA)
      }
      
      # gaussian function
      if (show.gaussian[k] & !is.null(C)) {
        if (!is.na(E[i]) & !is.na(C[i,i])) {
          x = seq(xmin[i],xmax[i],length=200)
          dx = x[2]-x[1]
          y = exp(-(x-E[i])^2/2/C[i,i])
          y = y/sum(y)/dx*area[i]
          x = (x-xmin[i])/(xmax[i]-xmin[i])+xoffset
          y = pmin(1,y)
          y[2:199] = y[1:198]*0.0005+y[2:199]*0.999000001+y[3:200]*0.0005
          y[y>=1] = NA
          lines(x,y+yoffset,col=col[k],lty=lty[k],lwd=lwd[k]*1.5)
        }
      }
      
    } else {
      
      # point cloud
      if (show.cloud[k] & !is.null(P)) {
        selection = (P[,i]>xmin[i]) & (P[,i]<xmax[i]) & (P[,j]>xmin[j]) & (P[,j]<xmax[j])
        if (!is.null(cloud.nmax[k])) {
          selection = seq(length(selection))[selection]
          npts = min(length(selection),cloud.nmax[k])
          selection = selection[1:npts]
        }
        xplot = (P[selection,i]-xmin[i])/(xmax[i]-xmin[i])+xoffset
        yplot = (P[selection,j]-xmin[j])/(xmax[j]-xmin[j])+yoffset
        points(xplot,yplot,pch=20,cex=0.15,col=cloud.col)
      }
      
      # central point
      if (show.expectation[k]) {
        if (!is.na(sum(E[c(i,j)]))) {
          x = (E[i]-xmin[i])/(xmax[i]-xmin[i])
          y = (E[j]-xmin[j])/(xmax[j]-xmin[j])
          if (x>=0 & x<=1 & y>=0 & y<=1) {
            points(x+xoffset,y+yoffset,
                   pch=pch[k],col=col[k],cex=cex[k])
          }
        }
      }
      
      # 68% ellipse
      if (show.ellipse.68[k] & !is.null(C)) {
        if (!is.na(sum(E[c(i,j)])) & !is.na(sum(C[c(i,j),c(i,j)]))) {
          pts = ellipse::ellipse(C[c(i,j),c(i,j)],centre=c(E[i],E[j]),level=0.68,draw=F)
          xplot = (pts[,1]-xmin[i])/(xmax[i]-xmin[i]) # convert the coordinates into the coordinates of the plot
          yplot = (pts[,2]-xmin[j])/(xmax[j]-xmin[j])
          lines(pmin(1,pmax(0,xplot))+xoffset,pmin(1,pmax(0,yplot))+yoffset,col=col[k],lty=lty[k],lwd=lwd[k]*2)
        }
      }
      
      # 95% ellipse
      if (show.ellipse.95[k] & !is.null(C)) {
        if (!is.na(sum(E[c(i,j)])) & !is.na(sum(C[c(i,j),c(i,j)]))) {
          pts = ellipse::ellipse(C[c(i,j),c(i,j)],centre=c(E[i],E[j]),level=0.95,draw=F)
          xplot = (pts[,1]-xmin[i])/(xmax[i]-xmin[i]) # convert the coordinates into the coordinates of the plot
          yplot = (pts[,2]-xmin[j])/(xmax[j]-xmin[j])
          lines(pmin(1,pmax(0,xplot))+xoffset,pmin(1,pmax(0,yplot))+yoffset,col=col[k],lty=lty[k],lwd=lwd[k])
        }
      }

    }

    # add axes
    if (k == 1) {
      tickpos = c(0.2,0.5,0.8)
      if (i==1 & j>1) {
        axis(2, at = yoffset+tickpos,labels=sprintf(text.format.numbers,E[j]+(xmax[j]-xmin[j])*(tickpos-0.5)),tck=0.015,lwd=0,lwd.ticks=1,cex.axis=text.size.numbers,padj = 0.8+text.offset.numbers[1])
        axis(2, at = yoffset+0.5,pos=-0.1+text.offset.labels[1],labels=names[[j]],cex.axis=text.size.labels,tick=F)
      }
      if (j==n) {
        axis(1, at = xoffset+tickpos,labels=sprintf(text.format.numbers,E[i]+(xmax[i]-xmin[i])*(tickpos-0.5)),tck=0.015,lwd=0,lwd.ticks=1,cex.axis=text.size.numbers,padj = -0.8+text.offset.numbers[2])
        axis(1, at = xoffset+0.5,pos=-0.1+text.offset.labels[2],labels=names[[i]],cex.axis=text.size.labels,tick=F)
      }
    }
  }

  # iterate over all parameter-pairs
  for (q in seq(0,m)) {
    
    if (q == 0) {
      k = master
    } else {
      k = order[q]
    }
    
    # make expectation E, covariance C, points P
    type = NA
    d = data[[k]]
    if (is.null(d)) {
      type = 0
    } else if (is.list(d)) {
      if (is.null(d$fit)) {
        
        # object type 2
        type = 2
        E = d[[1]]
        C = d[[2]]
        P = NULL
        
      } else {
        
        # object type 4
        type = 4
        E = d$fit$p.best
        C = d$fit$p.covariance
        P = NULL
        
      }
      n = length(E)
    } else {
      if (length(dim(d))<=1) {
        
        # object type 1
        type = 1
        n = length(d)
        E = d
        C = NULL
        P = NULL
        
      } else {
        
        # object type 3
        type = 3
        n = dim(d)[2]
        E = array(NA,n)
        for (i in seq(n)) {
          E[i] = mean(d[,i]) 
        }
        C = cov(d)
        P = d
        
      }
    }
    if (is.na(type)) stop('unknown object type in dfplotcov.')
    
    if (q == 0) {
      
      if (is.null(C)) stop('dfplotcov: the master object cannot be an expectation only.')
      
      # determin parameter names
      if (is.null(names)) {
        names = {}
        for (i in seq(n)) {
          names[[i]] = bquote(p [.(i)])
        }
        if (type == 4) {
          if (!is.null(d$model$parameter.names)) {
            for (i in seq(n)) {
              names[[i]] = d$model$parameter.names[i]
            }
          }
        }
      }
      
      # determine boundaries
      nref = n
      if (is.null(pmin) & is.null(pmax)) {
        xmin = xmax = xavg = array(NA,n)
        for (i in seq(n)) {
          xlength = sqrt(C[i,i])*nstd
          xavg[i] = E[i]
          xmin[i] = E[i]-xlength/2
          xmax[i] = E[i]+xlength/2
        }
      } else {
        if (length(pmin)==n & length(pmax)==n) {
          xmin = pmin
          xmax = pmax
          xavg = (pmin+pmax)/2
        } else {
          stop('dfplotcov: pmin and pmax must both be P-element vectors.')
        }
      }
      
      # determine histogram normalizations
      area = array(NA,n)
      for (i in seq(n)) {
        area[i] = 0.9*vert.stretch[i]*sqrt(2*pi*C[i,i])
      }
      
      # start a new plot
      pty.old = par()$pty
      mar.old = par()$mar
      par(pty='s')
      par(mar=margins)
      plot(0,0,type='n',xlim=c(0,n),ylim=c(0,n),ann=FALSE,xaxs='i',yaxs='i',xaxt='n',yaxt='n',bty='n')
        
    } else {
    
      if (n!=nref) stop('In dfplotcov, all data objects must have the same number of parameters.')
        
      if (type>0) {
        for (i in seq(n)) {
          for (j in seq(n)) {
            if (i==j | (j>i & lower) | (i>j & upper)) {
              make_single_square(i,j,k,E,C,P)
            }
          }
        }
      }
    }
  }

  text(1.2,n-0.2,title,offset=0,adj=0,cex=1)
  par(pty=pty.old)
  par(mar=mar.old)

}
