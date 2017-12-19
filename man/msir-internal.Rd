\name{msir-internal}
%
\alias{msir.fit}
\alias{eigen.decomp}
\alias{mvdnorm}
\alias{normalize}
\alias{msir.parameters}
\alias{msir.components}
\alias{msir.componentsSlice}
\alias{msir.recoverdir}
\alias{spinplot}

\title{Internal \code{'msir'} functions}

\description{Internal functions for package \code{'msir'}.}

\usage{
msir.fit(data, labels, G = NULL, modelNames = NULL, 
         control = emControl(itmax = c(.Machine$integer.max, 50)), 
         initialization = NULL, warn = FALSE, verbose = FALSE, \dots) 

eigen.decomp(X1, X2, inv = FALSE, tol = sqrt(.Machine$double.eps))

mvdnorm(x, mean, sigma, log = FALSE, tol = sqrt(.Machine$double.eps))

normalize(x)
    
msir.parameters(object, numdir = object$numdir)

msir.components(object)
msir.componentsSlice(object)

msir.recoverdir(object, data, normalized = TRUE, std = FALSE)

spinplot(data, \dots)

}

% \arguments{}

\details{These functions are not intended to be called by the user.}

\author{Luca Scrucca \email{luca@stat.unipg.it}}

\seealso{\code{\link{msir}}}

\keyword{internal}
\keyword{package}
