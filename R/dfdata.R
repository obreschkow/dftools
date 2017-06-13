#' Example data
#'
#' This function produces an example of log-masses \code{x} with Gaussian uncertainties \code{x.err} and a selection \code{selection}. This selection is a list of an array with galaxy-by-galaxy effective volumes and a function to complete the effective volume beyond the observed mass range. The data is taken from a publication by Westmeier et al. (2017) and represents HI-Masses in the Sculptor structure.
#'
#' @examples
#' data <- dfdata()
#' print(data$x)
#' print(data$x.err)
#' print(data$selection[[1]])
#'
#' # These data can then be used to fit a MF, e.g.
#' df = dffit(data$x, data$selection, data$x.err)
#'
#' @seealso \code{\link{dffit}}
#'
#' @author Danail Obreschkow
#'
#' @export

dfdata <- function() {

  # data (Westmeier 2017, rounded)
  mass = c(1889050000,1908130000,1572290000,10972900,9872600,1844110000,4276510,112822000,1124150000,332224000,241310000,26226300,8919190,5083300,2119200000,88838000,129331000,813553000,2744830000,23753200,245808000,28925700,344132000,1346930000,132612000,23345300,150683000,29690500,52450000,105312000,103877000)

  mass.sigma = c(609331000,558581000,567870000,3824130,1420330,529762000,5054950,34776300,414191000,88406800,64581400,3696160,2529850,4125370,261901000,49264700,19148400,264810000,1049500000,11766300,114450000,13311400,155055000,413508000,54283400,9067970,114012000,11052600,17015600,33073300,31396900)

  veff.values = c(98.47100,98.47100,98.47100,36.33020,32.59020,98.47100,7.48225,98.47100,98.47100,98.47100,98.47100,70.43790,24.81000,10.81120,98.47100,98.47100,98.47100,98.47100,98.47100,56.40330,98.47100,82.76920,98.47100,98.47100,98.47100,73.12730,98.47100,71.04220,89.04720,98.47100,98.47100)

  veff.fct = function(x) {pmax(0,pmin(98.471,(x-6.53)*75))}
  
  return(list(x = log10(mass), x.err = mass.sigma/mass/log(10),
              selection = list(veff.values = veff.values, veff.fct = veff.fct)))
}