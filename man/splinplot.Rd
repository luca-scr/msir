\name{spinplot}
\alias{spinplot}

\title{Rotating three-dimensional plot}

\description{General function to draw a rgl-based rotating 3D scatterplot.}

\usage{
spinplot(x, y, z, 
         scaling = c("abc", "aaa"), 
         rem.lin.trend = FALSE, 
         uncor.vars = FALSE, 
         fit.ols = FALSE, 
         fit.smooth = FALSE, 
         span = 0.75, 
         ngrid = 25, 
         markby, 
         pch.points = 1, 
         col.points = "black", 
         cex.points = 1,
         col.axis   = "gray", 
         col.smooth = "limegreen", 
         col.ols    = "lightsteelblue",
         background = "white", 
         \dots)
}

\arguments{
  \item{x}{a vector of values for the variable in the horizontal (H) screen axis.}
  \item{y}{a vector of values for the variable in the vertical (V) screen axis.}
  \item{z}{a vector of values for the variable in the out-of-screen (O) axis.}
  \item{scaling}{the scaling applied. Two possible values are \code{"abc"} and \code{"aaa"}.}
  \item{rem.lin.trend}{a logical specifying if the linear trend should be remove. If \code{TRUE} then the vertical axis is replaced by e(V|H,O), i.e. the residuals from a linear fit of the vertical axis variable on the others.}
  \item{uncor.vars}{a logical specifying if uncorrelated H and O variables should be used. If \code{TRUE} then the O variable is replaced by e(O|H), i.e. the residuals of the regression of O on H, hence obtaining a pair of uncorrelated variables.}
  \item{fit.ols}{a logical specifying if a fitted OLS plane should be included.}
  \item{fit.smooth}{a logical specifying if a nonparametric smoothing plane should be included.}
  \item{span}{the span used by \code{\link{loess}} to fit the polynomial surface.}
  \item{ngrid}{the number of grid points to use for displaing the fitted plane.}
  \item{markby}{a variable (usually a factor) to be used for marking the points.}
  \item{pch.points}{a vector of symbols for marking the points.}
  \item{col.points}{a vector of colors for marking the points.}
  \item{cex.points}{the cex for points.}
  \item{col.axis}{the color of the axis.}
  \item{col.ols}{the color to be used for drawing the OLS plane.}
  \item{col.smooth}{the color to be used for drawing the smoothing plane.}
  \item{background}{the color of background space.}
  \item{\dots}{catch further unused arguments.}
}

\details{This function is mainly based on the functionality of the \code{spin-plot} function once available in XLisp-Stat software \url{https://en.wikipedia.org/wiki/XLispStat}, and the adds-on introduced by the Arc software \url{http://www.stat.umn.edu/arc/index.html}.}

%\value{}

\references{
Cook R. D., Weisberg S. (1999) \emph{Applied Regression Including Computing and Graphics}, Wiley, Chapter 8
}

\author{Luca Scrucca \email{luca.scrucca@unipg.it}}

% \seealso{}

\examples{
\dontrun{
x1 <- rnorm(100)
x2 <- rnorm(100)
y  <- x1 + 2*x2^2 + 0.5*rnorm(100)
spinplot(x1, y, x2)
spinplot(x1, y, x2, scaling = "aaa")
spinplot(x1, y, x2, rem.lin.trend = "TRUE")
spinplot(x1, y, x2, fit.smooth = TRUE)
spinplot(x1, y, x2, fit.ols = TRUE)

x <- iris[,1:3]
y <- iris[,5]
spinplot(x)
spinplot(x, markby = y)
spinplot(x, markby = y, col.points = c("dodgerblue2", "orange", "green3"))
spinplot(x, markby = y, pch = c(0,3,1), col.points = c("dodgerblue2", "orange", "green3"))

# to save plots use 
# rgl.postscript("plot.pdf", fmt="pdf")
# or
# rgl.snapshot("plot.png")
}
}
