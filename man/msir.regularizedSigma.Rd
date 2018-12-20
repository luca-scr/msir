\name{msir.regularizedSigma}
\alias{msir.regularizedSigma}

\title{Regularized estimate of predictors covariance matrix.}

\description{This function computes a regularized version of the covariance matrix of the predictors. Among the possible models the one which maximizes BIC is returned.}

\usage{
msir.regularizedSigma(x, inv = FALSE, model = c("XII", "XXI", "XXX"))
}

\arguments{
  \item{x}{Ahe predictors data matrix.}
  \item{inv}{A logical specifying what must be returned. If \code{TRUE} the inverse of the estimated covariance matrix is returned, otherwise the estimated covariance matrix (default).}
  \item{model}{A character string specifying the available models:
  \itemize{
  \item \code{XII}: diagonal equal variances
  \item \code{XXI}: diagonal unequal variances
  \item \code{XXX}: full covariance matrix
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
