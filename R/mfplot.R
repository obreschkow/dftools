#' Display fitted galaxy mass function
#'
#' This function displays the galaxy mass function (MF) fitted using \code{\link{dffit}}, by calling \code{\link{dfplot}} with appropriate parameters.
#'
#' @seealso As an example run \code{\link{dfexample1()}}. See examples in \code{\link{dffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

mfplot <- function(survey,
                   xlab = expression('M [M'['sun']*']'),
                   ylab = expression(phi~'[Mpc'^-3~'dex'^-1~']'),
                   log = 'xy',
                   xpower10 = TRUE,
                   show.data.points = TRUE,
                   ...) {
  survey = dfplot(survey, xlab = xlab, ylab = ylab, log = log, xpower10 = xpower10,
              show.data.points = show.data.points, ...)
  invisible(survey)
}
