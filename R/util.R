#-------------------------------------------------------------------#
#                                                                   #
#  Miscellaneous functions                                          #
#                                                                   #
#-------------------------------------------------------------------#

mvdnorm <- function(x, mean, sigma, log = FALSE, tol = sqrt(.Machine$double.eps))
{
# Multivariate normal probability density function (pdf)
  if(is.vector(x)) 
    { x <- matrix(x, ncol = length(x)) }
  else
    { x <- as.matrix(x) }
  if(missing(mean)) 
    { mean <- rep(0, length = ncol(x)) }
  if(missing(sigma)) 
    { sigma <- diag(ncol(x)) }
  #    
  SVD <- svd(sigma)
  pos <- (SVD$d > max(tol*SVD$d[1], 0)) # in case of not full rank covar matrix
  inv.sigma <- SVD$v[,pos,drop=FALSE] %*% (1/SVD$d[pos] *
                                          t(SVD$u[,pos,drop=FALSE]))
  md <- mahalanobis(x, center = mean, cov = inv.sigma, inverted = TRUE)
  logdet <- sum(log(SVD$d[pos]))
  logdens <- -0.5*(ncol(x) * log(2 * pi) + logdet + md)
  if(log) return(logdens)
  else    return(exp(logdens))
}

normalize <- function(x)
{
# Normalize the vector x to have unit length
  x <- as.vector(x)
  x <- x/sqrt(crossprod(x)[1,1])
  return(x)
}

catwrap <- function(x, width = getOption("width"), ...)
{
# version of cat with wrapping at specified width
  cat(paste(strwrap(x, width = width, ...), collapse = "\n"), "\n")
}
