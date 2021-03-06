% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dfplotveff.R
\name{dfplotveff}
\alias{dfplotveff}
\title{Plot effective volume}
\usage{
dfplotveff(
  survey,
  xlim = NULL,
  ylim = NULL,
  xlab = "Mass M",
  ylab = expression("Effective volume"),
  legend = TRUE,
  xpower10 = TRUE
)
}
\arguments{
\item{survey}{List produced by \code{\link{dffit}}}

\item{xlim}{2-element vector with x-axis plotting limits}

\item{ylim}{2-element vector with y-axis plotting limits}

\item{xlab}{Label on x-axis (logarithmic mass or luminosity scale).}

\item{ylab}{Label on y-axis (volume scale).}

\item{legend}{Logical flag determining whether the legend will be plotted.}

\item{xpower10}{If \code{TRUE}, the model argument x is elevated to the power of 10 in the plots.}
}
\description{
This routine plots the effective survey volume function \code{Veff(x)} assumed when fitting a one-dimensional generative distribution function (GDF), such as a galaxy mass function. Note that this function \code{Veff(x)} is stored as \code{survey$selection$veff} when fitting a GDF using \code{survey=dffit(...)}.
}
\seealso{
See examples in \code{\link{dffit}}. To display \code{Veff(x)} of two-dimensional distribution functions use \code{\link{dfplotveff2}}.
}
\author{
Danail Obreschkow
}
