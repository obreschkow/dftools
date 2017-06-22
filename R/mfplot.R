#' Display fitted galaxy mass function
#'
#' This function displays the galaxy mass function (MF) fitted using \code{\link{dffit}}, by calling \code{\link{dfplot}} with appropriate parameters.
#' 
#' @param survey List produced by \code{\link{dffit}}
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param log String specifying the log-axes as in \code{\link{plot}}.
#' @param xpower10 If \code{TRUE}, the model argument x is elevated to the power of 10 in the plots.
#' @param ... Other of \code{\link{dfplot}}.
#'
#' @seealso As an example run \code{\link{dfexample1}}. See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

mfplot <- function(survey,
                   xlab = expression('M [M'['sun']*']'),
                   ylab = expression(phi~'[Mpc'^-3~'dex'^-1~']'),
                   log = 'xy',
                   xpower10 = TRUE,
                   ...) {
  survey = dfplot(survey, xlab = xlab, ylab = ylab, log = log, xpower10 = xpower10, ...)
  invisible(survey)
}
