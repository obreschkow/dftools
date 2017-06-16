#' Plot effective volume in 2D
#'
#' This function plots the 2D function \code{Veff(x)}, where \code{x} is a 2-element vector.
#'
#' @importFrom rgl open3d decorate3d surface3d
#'
#' @param df List produced by \code{\link{dffit}}
#' @param xlab Label on x-axis (e.g. logarithm of mass).
#' @param ylab Label on y-axis (e.g. logarithm of size or angular momentum).
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplotveff2 <- function(bundle,
                       xlab = 'Observable x1',
                       ylab = 'Observable x2',
                       zlab = 'Veff') {
  
  n.dim = dim(bundle$data$x)[2]
  if (n.dim!=2) {
    if (n.dim==1) {
      stop('Use dfplotveff for one-dimensional Veff functions.')
    } else {
      stop('dfplotveff2 only handles two-dimensional Veff functions.')
    }
  }
  
  x = seq(bundle$grid$xmin[1],bundle$grid$xmax[1],bundle$grid$dx[1])
  y = seq(bundle$grid$xmin[2],bundle$grid$xmax[2],bundle$grid$dx[2])
  
  z = bundle$grid$veff
  z = log10(pmax(1e-99,z))
  z = array(z,c(length(x),length(y)))
  zmax = max(z)
  z[z<(zmax-5)] = NA
  
  if (!is.null(bundle$input$selection$veff.input.function)) {
    z2 = bundle$input$selection$veff.input.function(bundle$grid$x)
    z2 = log10(pmax(1e-99,z2))
    z2 = array(z2,c(length(x),length(y)))
    zmax = max(z2)
    z2[z2<(zmax-5)] = NA
  }
  
  rgl::open3d()
  rgl::decorate3d(xlim = range(x), ylim = range(y), zlim = c(-5,0)+zmax,
                  xlab = xlab, ylab = ylab, zlab = zlab, aspect = TRUE)
  rgl::surface3d(x,y,z,col='blue',alpha=0.3, forceClipregion= TRUE)
  if (!is.null(bundle$input$selection$veff.input.function)) {
    rgl::surface3d(x,y,z2,col='black',alpha=0.3, forceClipregion= TRUE)
  }
  
}
