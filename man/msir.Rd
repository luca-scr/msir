\name{msir}
\alias{msir}
\alias{print.msir}

\title{Model-based Sliced Inverse Regression (MSIR)}

\description{A dimension reduction method based on Gaussian finite mixture models which provides an extension to sliced inverse regression (SIR). The basis of the subspace is estimated by modeling the inverse distribution within slice using Gaussian finite mixtures with number of components and covariance matrix parameterization selected by BIC or defined by the user.}

\usage{
msir(x, y, nslices = msir.nslices, slice.function = msir.slices, 
     modelNames = NULL, G = NULL, cov = c("mle", "regularized"), \dots)
}

\arguments{
  \item{x}{a \eqn{(n \times p)}{(n x p)} design matrix containing the predictors data values.}
  \item{y}{a \eqn{(n \times 1)}{(n x 1)} vector of data values for the response variable. It can be a numeric vector (regression) but also a factor (classification). In the latter case, the levels of the factor define the slices used.}
  \item{nslices}{the number of slices used, unless \code{y} is a factor. By default the value returned by \code{\link{msir.nslices}}.}
  \item{slice.function}{the slice functions to be used, by default \code{\link{msir.slices}}, but the user can provide a different slicing function.}
  \item{modelNames}{a vector of character strings indicating the Gaussian mixture models to be fitted as described in \code{\link[mclust]{mclustModelNames}}. If a  vector of strings is given they are used for all the slices. If a list of vectors is provided then each vector refers to a single slice.}
  \item{G}{an integer vector specifying the numbers of mixture components used in fitting Gaussian mixture models. If a list of vectors is provided then each vector refers to a single slice.}
  \item{cov}{The predictors marginal covariance matrix. Possible choices are:
  \tabular{lcl}{
    \code{"mle"} \tab = \tab maximum likelihood estimate\cr
    \code{"regularized"} \tab = \tab regularized covariance matrix (see  \code{\link{msir.regularizedSigma}})\cr
    \code{R matrix} \tab = \tab a \eqn{(p \times p)}{(p x p)} user defined covariance matrix
  }
}
  \item{\dots}{other arguments passed to \code{msir.compute}.}
}

%\details{}

\value{
Returns an object of class \code{'msir'} with attributes:
  \item{call}{the function call.}
  \item{x}{the design matrix.}
  \item{y}{the response vector.}
  \item{slice.info}{output from slicing function.}
  \item{mixmod}{a list of finite mixture model objects as described in  \code{\link[mclust]{mclustModel}}.}
  \item{loglik}{the log-likelihood for the mixture models.}
  \item{f}{a vector of length equal to the total number of mixture components containing the fraction of observations in each fitted component within slices.}
  \item{mu}{a matrix of component within slices predictors means.}
  \item{sigma}{the marginal predictors covariance matrix.}
  \item{M}{the msir kernel matrix.}
  \item{evalues}{the eigenvalues from the generalized eigen-decomposition of \code{M}.}
  \item{evectors}{the raw eigenvectors from the generalized eigen-decomposition of \code{M} ordered according to the eigenvalues.}
  \item{basis}{the normalized eigenvectors from the generalized eigen-decomposition of \code{M} ordered according to the eigenvalues.}
  \item{std.basis}{standardized basis vectors obtained by multiplying each coefficient of the eigenvectors by the standard deviation of the corresponding predictor. The resulting coefficients are scaled such that all predictors have unit standard deviation.}
  \item{numdir}{the maximal number of directions estimated.}
  \item{dir}{the estimated MSIR directions from mean-centered predictors.}
}

\references{Scrucca, L. (2011) Model-based SIR for dimension reduction. \emph{Computational Statistics & Data Analysis}, 55(11), 3010-3026.}

\author{Luca Scrucca \email{luca.scrucca@unipg.it}}
% \note{}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{\code{\link[msir]{summary.msir}}, \code{\link[msir]{plot.msir}}, \code{\link[dr]{dr}}}

\examples{
# 1-dimensional simple regression
n = 200; p = 5
b = as.matrix(c(1,-1,rep(0,p-2)))
x = matrix(rnorm(n*p), nrow = n, ncol = p)
y = exp(0.5 * x\%*\%b) + 0.1*rnorm(n)
MSIR = msir(x, y)
summary(MSIR)
plot(MSIR, type = "2Dplot")

# 1-dimensional symmetric response curve
n = 200; p = 5
b = as.matrix(c(1,-1,rep(0,p-2)))
x = matrix(rnorm(n*p), nrow = n, ncol = p)
y = (0.5 * x\%*\%b)^2 + 0.1*rnorm(n)
MSIR = msir(x, y)
summary(MSIR)
plot(MSIR, type = "2Dplot")
plot(MSIR, type = "coefficients")

# 2-dimensional response curve
n = 300; p = 5
b1 = c(1, 1, 1, rep(0, p-3))
b2 = c(1,-1,-1, rep(0, p-3))
b = cbind(b1,b2)
x = matrix(rnorm(n*p), nrow = n, ncol = p)
y = x \%*\% b1 + (x \%*\% b1)^3 + 4*(x \%*\% b2)^2 + rnorm(n)
MSIR = msir(x, y)
summary(MSIR)
plot(MSIR, which = 1:2)
\dontrun{plot(MSIR, type = "spinplot")}
plot(MSIR, which = 1, type = "2Dplot", span = 0.7)
plot(MSIR, which = 2, type = "2Dplot", span = 0.7)
}
\keyword{regression}
\keyword{multivariate}
