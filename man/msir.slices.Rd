\name{msir.slices}
\alias{msir.slices}

\title{Slice a vector into slices of approximately equal size}

\description{Function used for slicing a continuous response variable.}

\usage{
msir.slices(y, nslices)
}

\arguments{
  \item{y}{a vector of \eqn{n} values}
  \item{nslices}{the number of slices, no larger than \eqn{n}}
}

%\details{}

\value{Returns a list with components:
  \item{slice.indicator}{an indicator variable for the slices.}
  \item{nslices}{the actual number of slices produced.}
  \item{slice.sizes}{the number of observations in each slice.}
}

%\references{}
\author{Luca Scrucca \email{luca.scrucca@unipg.it}}
\seealso{\code{\link{msir}}}
