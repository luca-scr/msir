\name{predict.msir}
\alias{predict.msir}

\title{Model-based Sliced Inverse Regression directions}

\description{MSIR estimates a set of \eqn{d \le p} orthogonal direction vectors of length \eqn{p} which are estimates of the basis of the dimensional reduction subspace.}

\usage{
  \method{predict}{msir}(object, dim = 1:object$numdir, newdata, \dots)
}

\arguments{
  
  \item{object}{an object of class \code{'msir'} resulting from a call to \code{\link{msir}.}}

  \item{dim}{the dimensions of the reduced subspace used for prediction.}

  \item{newdata}{a data frame or matrix giving the data. If missing the data obtained from the call to \code{\link{msir}} are used.}

  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{}

\value{
The function returns a matrix of points projected on the subspace spanned by the estimated basis vectors.
}

\references{Scrucca, L. (2011) Model-based SIR for dimension reduction. \emph{Computational Statistics & Data Analysis}, 55(11), 3010-3026.}

\author{Luca Scrucca \email{luca.scrucca@unipg.it}}

\seealso{\code{{msir}}}

\examples{
n <- 200; p <- 5
b <- as.matrix(c(1,-1,rep(0,p-2)))
x <- matrix(rnorm(n*p), nrow = n, ncol = p)
y <- exp(0.5 * x\%*\%b) + 0.1*rnorm(n)
pairs(cbind(y,x), gap = 0)

MSIR <- msir(x, y)
summary(MSIR)
plot(MSIR, which = 1, type = "2Dplot")
all.equal(predict(MSIR), MSIR$dir)
predict(MSIR, dim = 1:2)

x0 <- matrix(rnorm(n*p), nrow = n, ncol = p)
y0 <- exp(0.5 * x0\%*\%b) + 0.1*rnorm(n)
plot(predict(MSIR, dim = 1, newdata = x0), y0)
}
