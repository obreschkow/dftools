#' Write best fitting parameters
#'
#' This function displays the best-fitting model parameters of a mass function previously fitted using \code{\link{dffit}}.
#'
#' @param survey List produced by \code{\link{dffit}}
#'
#' @examples
#' dat = dfmockdata()
#' survey = dffit(dat$x, dat$veff, write.fit = FALSE)
#' dfwrite(survey)
#'
#' @seealso See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

dfwrite <- function(survey) {

  p = survey$fit$p.best

  if (!is.na(survey$model$gdf.equation)) {
    cat(sprintf('%s\n',survey$model$gdf.equation))
  }

  if (!survey$fit$status$converged) {
    for (i in seq(length(p))) {
      cat(sprintf('p[%d] = %7.2f (not converged)\n',i,p[i]))
    }
  } else {
    if (length(survey$fit$p.quantile.16)>0) {
      sigma.84 = survey$fit$p.quantile.84-survey$fit$p.best
      sigma.16 = survey$fit$p.best-survey$fit$p.quantile.16
      for (i in seq(length(p))) {
        cat(sprintf('p[%d] = %7.2f (+%3.2f -%3.2f)\n',i,p[i],sigma.84[i],sigma.16[i]))
      }
    } else {
      sigma = survey$fit$p.sigma
      for (i in seq(length(p))) {
        cat(sprintf('p[%d] = %7.2f (+-%3.2f)\n',i,p[i],sigma[i]))
      }
    }
  }

}
