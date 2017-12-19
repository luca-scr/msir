\name{summary.msir}
%
\alias{summary.msir}
\alias{print.summary.msir}

\title{Summary and print methods for \code{'msir'} objects}

\description{Summary and print methods for \code{'msir'} objects.}

\usage{
\method{summary}{msir}(object, numdir = object$numdir, std = FALSE, verbose = TRUE, \dots)
\method{print}{summary.msir}(x, digits = max(5, getOption("digits") - 3), \dots )
}

\arguments{
  \item{object}{a \code{'msir'} object}
  \item{numdir}{the number of directions to be shown.}
  \item{std}{if \code{TRUE} the coefficients basis are scaled such that all predictors have unit standard deviation.}
  \item{verbose}{if \code{FALSE} the coefficients basis are omitted; by default \code{verbose = TRUE}.}  
  \item{x}{a \code{'summary.msir'} object.}
  \item{digits}{the significant digits to use.}
  \item{\dots}{additional arguments.}
}

%\details{}
%\value{}
%\references{}
\author{Luca Scrucca \email{luca@stat.unipg.it}}
\seealso{\code{\link{msir}}}
%\examples{}
