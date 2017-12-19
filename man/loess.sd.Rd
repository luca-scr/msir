\name{loess.sd}
\alias{loess.sd}
\alias{panel.loess}

\title{Local Polynomial Regression Fitting with Variability bands}

\description{Nonparametric estimation of mean function with variability bands.}

\usage{
loess.sd(x, y = NULL, nsigma = 1, ...)

panel.loess(x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, 
            col.smooth = "red", span = 2/3, degree = 2, nsigma = 1, ...) 
}

\arguments{
  \item{x}{a vector of values for the predictor variable \eqn{x}.}
  \item{y}{a vector of values for the response variable \eqn{y}.}
  \item{nsigma}{a multiplier for the standard deviation function.}

  \item{col, bg, pch, cex}{numeric or character codes for the color(s), point type and size of points; see also \code{\link{par}}.}
  \item{col.smooth}{color to be used by \code{lines} for drawing the smooths.}
  \item{span}{smoothing parameter for \code{\link{loess}}.}
  \item{degree}{the degree of the polynomials to be used, see \code{\link{loess}}.}

  \item{...}{further argument passed to the function \code{\link{loess}}.}
}

%\details{}

\value{The function \code{loess.sd} computes the loess smooth for the mean function and the mean plus and minus \code{k} times the standard deviation function.
  
The function \code{panel.loess} can be used to add to a scatterplot matrix panel a smoothing of mean function using loess with variability bands at plus and minus \code{nsigmas} times the standard deviation.}

\references{Weisberg, S. (2005) Applied Linear Regression, 3rd ed., Wiley, New York, pp. 275-278.}

\author{Luca Scrucca \email{luca@stat.unipg.it}}

\seealso{\code{\link{loess}}}

\examples{
data(cars)
plot(cars, main = "lowess.sd(cars)")
lines(l <- loess.sd(cars))
lines(l$x, l$upper, lty=2)
lines(l$x, l$lower, lty=2)
}

\keyword{loess}
