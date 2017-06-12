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

dfwrite <- function(df) {

  p = df$fit$parameters$p.optimal

  if (!is.na(df$input$distribution.function$phi.equation)) {
    cat(sprintf('%s\n',df$input$distribution.function$phi.equation))
  }

  if (!df$fit$status$converged) {
    for (i in seq(length(p))) {
      cat(sprintf('p[%d] = %7.2f (not converged)\n',i,p[i]))
    }
  } else {
    if (length(df$fit$parameters$p.quantile.16)>0) {
      sigma.84 = df$fit$parameters$p.quantile.84-df$fit$parameters$p.optimal
      sigma.16 = df$fit$parameters$p.optimal-df$fit$parameters$p.quantile.16
      for (i in seq(length(p))) {
        cat(sprintf('p[%d] = %7.2f (+%3.2f -%3.2f)\n',i,p[i],sigma.84[i],sigma.16[i]))
      }
    } else {
      sigma = df$fit$parameters$p.sigma
      for (i in seq(length(p))) {
        cat(sprintf('p[%d] = %7.2f (+-%3.2f)\n',i,p[i],sigma[i]))
      }
    }
  }

}
