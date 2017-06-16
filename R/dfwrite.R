#' Write best fitting parameters
#'
#' This function displays the best-fitting model parameters of a mass function previously fitted using \code{\link{dffit}}.
#'
#' @param df List produced by \code{\link{dffit}}
#'
#' @examples
#' data = mfdata()
#' df = dffit(data$x, data$selection, write.fit = FALSE)
#' dfwrite(df)
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfwrite <- function(bundle) {

  p = bundle$fit$p.best

  if (!is.na(bundle$model$gdf.equation)) {
    cat(sprintf('%s\n',bundle$model$gdf.equation))
  }

  if (!bundle$fit$status$converged) {
    for (i in seq(length(p))) {
      cat(sprintf('p[%d] = %7.2f (not converged)\n',i,p[i]))
    }
  } else {
    if (length(bundle$fit$p.quantile.16)>0) {
      sigma.84 = bundle$fit$p.quantile.84-bundle$fit$p.best
      sigma.16 = bundle$fit$p.best-bundle$fit$p.quantile.16
      for (i in seq(length(p))) {
        cat(sprintf('p[%d] = %7.2f (+%3.2f -%3.2f)\n',i,p[i],sigma.84[i],sigma.16[i]))
      }
    } else {
      sigma = bundle$fit$p.sigma
      for (i in seq(length(p))) {
        cat(sprintf('p[%d] = %7.2f (+-%3.2f)\n',i,p[i],sigma[i]))
      }
    }
  }

}
