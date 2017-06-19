#' Plot effective volume in 2D
#'
#' This routine plots the effective survey volume function \code{Veff(x)} assumed when fitting a two-dimensional generative distribution function (GDF). Note that this function \code{Veff(x)} is stored as \code{survey$selection$veff} when fitting a GDF using \code{survey=dffit(...)}.
#'
#' @importFrom rgl open3d decorate3d surface3d
#'
#' @param df List produced by \code{\link{dffit}}
#' @param xlab Label on x-axis (e.g. logarithm of mass).
#' @param ylab Label on y-axis (e.g. logarithm of size or angular momentum).
#'
#' @seealso See examples in \code{\link{dffit}}. To display \code{Veff(x)} of one-dimensional distribution functions use \code{\link{dfplotveff}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfplotveff2 <- function(survey,
                       xlab = 'Observable x1',
                       ylab = 'Observable x2',
                       zlab = 'Veff') {
  
  n.dim = dim(survey$data$x)[2]
  if (n.dim!=2) {
    if (n.dim==1) {
      stop('Use dfplotveff for one-dimensional Veff functions.')
    } else {
      stop('dfplotveff2 only handles two-dimensional Veff functions.')
    }
  }
  
  x = seq(survey$grid$xmin[1],survey$grid$xmax[1],survey$grid$dx[1])
  y = seq(survey$grid$xmin[2],survey$grid$xmax[2],survey$grid$dx[2])
  
  z = survey$grid$veff
  z = log10(pmax(1e-99,z))
  z = array(z,c(length(x),length(y)))
  zmax = max(z)
  z[z<(zmax-5)] = NA
  
  if (!is.null(survey$input$selection$veff.input.function)) {
    z2 = survey$input$selection$veff.input.function(survey$grid$x)
    z2 = log10(pmax(1e-99,z2))
    z2 = array(z2,c(length(x),length(y)))
    zmax = max(z2)
    z2[z2<(zmax-5)] = NA
  }
  
  rgl::open3d()
  rgl::decorate3d(xlim = range(x), ylim = range(y), zlim = c(-5,0)+zmax,
                  xlab = xlab, ylab = ylab, zlab = zlab, aspect = TRUE)
  rgl::surface3d(x,y,z,col='blue',alpha=0.3, forceClipregion= TRUE)
  if (!is.null(survey$input$selection$veff.input.function)) {
    rgl::surface3d(x,y,z2,col='black',alpha=0.3, forceClipregion= TRUE)
  }
  
}
