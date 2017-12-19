\name{msir.regularizedSigma}
\alias{msir.regularizedSigma}

\title{Regularized estimate of predictors covariance matrix.}

\description{This function computes a regularized version of the covariance matrix of the predictors. Among the possible models the one which maximizes BIC is returned.}

\usage{
msir.regularizedSigma(x, inv = FALSE, model = c("XII", "XXI", "XXX"))
}

\arguments{
  \item{x}{the predictors data matrix.}
  \item{inv}{if \code{TRUE} the inverse of the estimated covariance matrix is returned.}
  \item{model}{available models:
  \tabular{lcl}{
  \code{XII} \tab = \tab diagonal equal variances \cr
  \code{XXI} \tab = \tab diagonal unequal variances \cr
  \code{XXX} \tab = \tab full covariance matrix
  }
}
}

%\details{}

\value{A \eqn{(p \times p)}{(p x p)} covariance matrix estimate.}

% \references{}

\author{Luca Scrucca \email{luca.scrucca@unipg.it}}
\seealso{\code{\link{msir}}}

%\examples{}
%\keyword{ ~kwd1 }
