\name{dftools-package}
\alias{dftools-package}
\alias{dftools}
\docType{package}
\title{
Fitting distribution functions such as galaxy mass functions
}
\description{
The package can be used to find the most likely P parameters of a D-dimensional distribution function (DF) generating N objects with D observables per object. For instance, if the objects are galaxies, it can fit a mass function (P=1), a mass-size distribution (P=2) or the mass-spin-morphology distribution (P=3). Fits are performed using the modified maximum likelihood (MML) formalism with the fit-and-debias algorithm. Unlike most common fitting techniques, this method accurately accounts for measurement is uncertainties and complex selection functions. In the context of astrophysical applications, such as fitting a galaxy mass function, the routines can also accurately account for cosmic large-scale structure, inherent to any galaxy survey. A full description of the formalism and algorithm can be found in Obreschkow et al. (2017). The package contains the following routines:\cr\cr
\code{\link{dffit}} is the core routine, used to fit model parameters to a set of observations.\cr\cr
\code{\link{dfplot}}, \code{\link{mfplot}}, \code{\link{dfplot2}} are plotting functions that use the output argument of \code{dffit} as primary input argument. These functions can visualize the input data and fitted model in a single plot.\cr\cr
\code{\link{dfplotveff}} and \code{\link{dfplotveff2}} are a second set of plotting functions that use the output argument of \code{dffit} as primary input argument. These functions visualize the effective volume (= weight) as a function of the observables.\cr\cr
\code{\link{dfmodel}} and \code{\link{dfswmodel}} produce analytical distribution functions, such as galaxy mass functions, to be fitted by \code{dffit}.\cr\cr
\code{\link{dfwrite}} takes the output argument of \code{dffit} and writes the best-fitting solution and additional information on the screen.\cr\cr
\code{\link{dfmockdata}} produces a mock survey given a generative DF, such as a galaxy mass function, and a selection function.\cr\cr
\code{\link{dfexample}} allows the user to run a few illustrative examples of fitting galaxies.
}
\details{
\tabular{ll}{
Package: \tab dftools\cr
Type: \tab Package\cr
Version: \tab 0.5\cr
Date: \tab 2017-06-11\cr
License: \tab GPL-3\cr
}
}
\author{
Danail Obreschkow <danail.obreschkow@icrar.org>
}
\references{
Obreschkow et al., 2017
}
\keyword{ package }
