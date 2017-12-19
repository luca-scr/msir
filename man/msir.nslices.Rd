\name{msir.nslices}
\alias{msir.nslices}

\title{Default number of slices}

\description{This function computes a Sturges' type number of slices to be used as default in the \code{\link{msir}} function.}

\usage{
msir.nslices(n, p)
}

\arguments{
  \item{n}{the number of observations in the sample.}
  \item{p}{the number of predictors in the sample.}
}

%\details{}

\value{The function returns a single value, i.e. the number of slices.}

%\references{}
\author{Luca Scrucca \email{luca@stat.unipg.it}}
\seealso{\code{\link{msir}}}
