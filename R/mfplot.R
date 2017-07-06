#' Display fitted galaxy mass function
#'
#' This function displays the galaxy mass function (MF) fitted using \code{\link{dffit}}, by calling \code{\link{dfplot}} with appropriate parameters.
#' 
#' @param survey List produced by \code{\link{dffit}}
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param log String specifying the log-axes as in \code{\link{plot}}.
#' @param veff Function of x to be plotted over the histogram, if \code{show.data.histogram = TRUE}.
#' @param xpower10 If \code{TRUE}, the model argument x is elevated to the power of 10 in the plots.
#' @param show.data.histogram If \code{TRUE}, a histogram of source counts, based on the input data, is displayed in a bottom panel.
#' @param ... All other arguments of \code{\link{dfplot}}.
#'
#' @seealso As an example run \code{dfexample(1)}. See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

mfplot <- function(survey,
                   xlab = expression('M [M'['sun']*']'),
                   ylab = expression(phi~'[Mpc'^-3~'dex'^-1~']'),
                   log = 'xy',
                   veff = survey$selection$veff,
                   xpower10 = TRUE,
                   show.data.histogram = TRUE,
                   ...) {
  survey = dfplot(survey, xlab = xlab, ylab = ylab, log = log, veff = veff,
                  xpower10 = xpower10, show.data.histogram = show.data.histogram, ...)
  invisible(survey)
}
