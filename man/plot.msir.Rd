\name{plot.msir}
\alias{plot.msir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for \code{'msir'} objects.}

\description{Plots directions and other information from MSIR estimation.}

\usage{
\method{plot}{msir}(x, which, 
     type = c("pairs", "2Dplot", "spinplot", "evalues", "coefficients"), 
     span = NULL, std = TRUE, ylab, xlab, restore.par = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{'msir'} object.}
  \item{which}{a vector of value(s) giving the directions for which the plot should be drawn.}
  \item{type}{the type of plot to be drawn.}
  \item{span}{the span of smoother (only for \code{type = "pairs" | "2Dplot"}).}
  \item{std}{if \code{TRUE} coefficients are standardized (only for \code{type = "coefficients"}).}
  \item{ylab}{a character string for the y-axis label.}
  \item{xlab}{a character string for the x-axis label.}
  \item{restore.par}{if \code{TRUE} the graphical parameters (see \code{\link{par}}) changed are restored to the previous state. If you want to manipulate the resulting plot you should set \code{restore.par = FALSE}.}
  \item{...}{additional arguments.}
}

%\details{}

%\value{}

\references{Scrucca, L. (2011) Model-based SIR for dimension reduction. \emph{Computational Statistics & Data Analysis}, 55(11), 3010-3026.}

\author{Luca Scrucca \email{luca.scrucca@unipg.it}}

\seealso{\code{\link{msir}}}

\examples{
# 2-dimensional response curve
n = 300; p = 5
b1 = c(1, 1, 1, rep(0, p-3))
b2 = c(1,-1,-1, rep(0, p-3))
b = cbind(b1,b2)
x = matrix(rnorm(n*p), nrow = n, ncol = p)
y = x \%*\% b1 + (x \%*\% b1)^3 + 4*(x \%*\% b2)^2 + rnorm(n)
MSIR = msir(x, y)
summary(MSIR)
plot(MSIR)
plot(MSIR, which = 1:2)
plot(MSIR, type = "2Dplot", which = 1, span = 0.7)
plot(MSIR, type = "2Dplot", which = 2, span = 0.7)
\dontrun{plot(MSIR, type = "spinplot")}
plot(MSIR, type = "evalues")
plot(MSIR, type = "coefficients")
}

\keyword{regression}
\keyword{dplot}