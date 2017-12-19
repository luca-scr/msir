\name{msir.dir}
\alias{msir.dir}

\title{Model-based Sliced Inverse Regression directions}

\description{MSIR estimates a set of up to \eqn{p} orthogonal direction vectors each of length \eqn{p}, the first \eqn{d} of which are estimates of the basis of the dimensional reduction subspace.}

\usage{
msir.dir(object, numdir = object$numdir)
}

\arguments{
  \item{object}{a \code{'msir'} object}
  \item{numdir}{the number of basis vectors to return.}
}

%\details{}

\value{The function returns the estimated directions in the original \eqn{n} dimensional space for plotting.}

\references{Scrucca, L. (2011) Model-based SIR for dimension reduction. \emph{Computational Statistics & Data Analysis}, 55(11), 3010-3026.}

\author{Luca Scrucca \email{luca@stat.unipg.it}}

\seealso{\code{{msir}}}

\examples{
# 1-dimensional simple regression
n = 200; p = 5
b = as.matrix(c(1,-1,rep(0,p-2)))
x = matrix(rnorm(n*p), nrow = n, ncol = p)
y = exp(0.5 * x\%*\%b) + 0.1*rnorm(n)
MSIR = msir(x, y)
msir.dir(MSIR, numdir = 1)
}

\keyword{regression}
