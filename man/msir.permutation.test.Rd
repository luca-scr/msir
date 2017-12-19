\name{msir.permutation.test}
\alias{msir.permutation.test}

\title{Permutation test for dimensionality}

\description{Approximates marginal dimension test significance levels by sampling from the permutation distribution.}

\usage{
msir.permutation.test(object, npermute = 99, numdir = object$numdir, verbose = TRUE)
}

\arguments{
  \item{object}{a \code{'msir'} object.}
  \item{npermute}{number of permutations to compute.}
  \item{numdir}{maximum value of the dimension to test.}
  \item{verbose}{if \code{TRUE} a textual progress bar is shown during computation.}
}

\details{The function approximates significance levels of the marginal dimension tests based on a permutation test.}

\value{
The function returns a list with components:
\item{summary}{a table containing the hypotheses, the test statistics, the permutation p-values.}
\item{npermute}{the number of permutations used.}
Furthermore, it also assigns the above information to the corresponding \code{'msir'} object.
}

\references{Scrucca, L. (2011) Model-based SIR for dimension reduction. \emph{Computational Statistics & Data Analysis}, 55(11), 3010-3026.}

\author{Luca Scrucca \email{luca@stat.unipg.it}}

\seealso{\code{\link[dr]{dr}}}

\examples{
\dontrun{
# 1-dimensional simple regression
n = 200; p = 5
b = as.matrix(c(1,-1,rep(0,p-2)))
x = matrix(rnorm(n*p), nrow = n, ncol = p)
y = exp(0.5 * x\%*\%b) + 0.1*rnorm(n)
MSIR = msir(x, y)
msir.permutation.test(MSIR)
summary(MSIR)
}
}

\keyword{regression}
\keyword{htest}
