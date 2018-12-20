\name{msir.bic}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{msir.bic}
\alias{bicDimRed}

\title{BIC-type criterion for dimensionality}

\description{BIC-type criterion for selecting the dimensionality of a dimension reduction subspace.}

\usage{
msir.bic(object, type = 1, plot = FALSE)

bicDimRed(M, x, nslices, type = 1, tol = sqrt(.Machine$double.eps))
}

\arguments{
  \item{object}{a \code{'msir'} object}
  \item{plot}{if \code{TRUE} a plot of the criterion is shown.}
  \item{M}{the kernel matrix. See details below.}
  \item{x}{the predictors data matrix. See details below.}
  \item{type}{See details below.}
  \item{nslices}{the number of slices. See details below.}
  \item{tol}{a tolerance value}
}

\details{
This BIC-type criterion for the determination of the structural dimension selects \eqn{d} as the maximizer of
\deqn{G(d) = l(d) - Penalty(p,d,n)}
where \eqn{l(d)} is the log-likelihood for dimensions up to \eqn{d}, \eqn{p} is the number of predictors, and \eqn{n} is the sample size.
The term \eqn{Penalty(p,d,n)} is the type of penalty to be used:
  \itemize{
  \item \code{type = 1}: \eqn{Penalty(p,d,n) = -(p-d) \log(n)}
  \item \code{type = 2}: \eqn{Penalty(p,d,n) = 0.5 C d (2p-d+1)}, where \eqn{C = (0.5 \log(n) + 0.1 n^(1/3))/2  nslices/n}
  \item \code{type = 3}: \eqn{Penalty(p,d,n) = 0.5 C d (2p-d+1)}, where \eqn{C = \log(n) nslices/n}
  \item \code{type = 4} \eqn{Penalty(p,d,n) = 1/2 d \log(n)}
  }
}

\value{
Returns a list with components:
\item{evalues}{eigenvalues}
\item{l}{log-likelihood}
\item{crit}{BIC-type criterion}
\item{d}{selected dimensionality}

The \code{msir.bic} also assign the above information to the corresponding \code{'msir'} object.
}

\references{
Zhu, Miao and Peng (2006) "Sliced Inverse Regression for CDR Space Estimation", JASA.\cr
Zhu, Zhu (2007) "On kernel method for SAVE", Journal of Multivariate Analysis.
}

\author{Luca Scrucca \email{luca.scrucca@unipg.it}}

%\note{}
\seealso{\code{\link{msir}}}
\examples{
# 1-dimensional symmetric response curve
n <- 200
p <- 5
b <- as.matrix(c(1,-1,rep(0,p-2)))
x <- matrix(rnorm(n*p), nrow = n, ncol = p)
y <- (0.5 * x\%*\%b)^2 + 0.1*rnorm(n)
MSIR <- msir(x, y)
msir.bic(MSIR, plot = TRUE)
summary(MSIR)
msir.bic(MSIR, type = 3, plot = TRUE)
summary(MSIR)
}

\keyword{regression}
