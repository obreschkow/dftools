\name{dftools-package}
\alias{dftools-package}
\alias{dftools}
\docType{package}
\title{
Fitting of generative distribution functions such as galaxy mass functions
}
\description{
This package can find the most likely P-dimensional distribution function (DF) generating N objects with P observables. For instance, if the objects are galaxies, it can fit a MF (P=1), a mass-size distribution (P=2) or the mass-spin-morphology distribution (P=3). Unlike most common fitting approaches, this method accurately accounts for measurement is uncertainties and non-detections and it does not bin the data. A full description of the algorithm can be found in Obreschkow et al. (2017).
}
\details{
\tabular{ll}{
Package: \tab dftools\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
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
\seealso{
\code{\link{dffit}}
}
\examples{
See examplies in \code{\link{dffit}}.
}
