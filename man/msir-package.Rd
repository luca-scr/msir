\name{msir-package}
\alias{msir-package}
\docType{package}

\title{Model-based Sliced Inverse Regression (MSIR)}

\description{An R package that implements MSIR, a dimension reduction method based on Gaussian finite mixture models. The basis of the subspace is estimated by modeling the inverse distribution within slice using finite mixtures of Gaussians, with number of components and covariance matrix parameterization selected by BIC or defined by the user.
The method provides an extension to sliced inverse regression (SIR) and allows to overcome the main limitation of SIR, i.e., the failure in the presence of regression symmetric relationships, without the need to impose further assumptions.}

% \details{}

\seealso{\code{\link{msir}}}

\author{Luca Scrucca \email{luca.scrucca@unipg.it}}

\references{
Scrucca, L. (2011) Model-based SIR for dimension reduction. 
\emph{Computational Statistics & Data Analysis}, 55(11), 3010-3026.}

% \examples{}

\keyword{package}
