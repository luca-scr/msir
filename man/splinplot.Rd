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
         background = "white", 
         col.points = "black", 
         col.axis   = "gray", 
         cex.points = 1,
         col.smooth = "limegreen", 
         col.ols    = "lightsteelblue",
         \dots)
}

\arguments{
  \item{x}{}
  \item{y}{}
  \item{z}{}
  \item{scaling}{}
  \item{rem.lin.trend}{}
  \item{uncor.vars}{}
  \item{fit.ols}{}
  \item{fit.smooth}{}
  \item{span}{}
  \item{ngrid}{}
  \item{markby}{}
  \item{pch.points}{}
  \item{background}{}
  \item{col.points}{}
  \item{col.axis}{}
  \item{cex.points}{}
  \item{col.smooth}{}
  \item{col.ols}{}
  \item{\dots}{}
}

\details{The function is mainly based on the \code{spin-plot} function once available in XLisp-Stat software.}

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
rgl.postscript("plot.pdf", fmt="pdf")
# or
rgl.snapshot("plot.png")
}
}
