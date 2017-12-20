##
## Permutation test ----
##

msir.permutation.test <- function(object, npermute = 99, numdir = object$numdir, verbose = TRUE)
{
  test.statistic <- function(object, nd)
    { length(object$y) * rev(cumsum(rev(object$evalues)))[1:nd] }
  call <- object$call
  n <- length(object$y)
  z <- scale(object$x, center = TRUE, scale = FALSE) %*% object$evectors
  nd <- min(numdir,length(which(abs(object$evalues)>1.e-8))-1)
  nt   <- nd + 1
  obstest <- test.statistic(object, nt)
  count <- val <- rep(0, nt)
  ysl <- object$slice.info$slice.indicator
  nslices <- object$slice.info$nslices
  G <- lapply(object$mixmod, function(m) m$G)
  modelNames <- lapply(object$mixmod, function(m) m$modelName)
  call$G <- G
  call$modelNames <- modelNames
  #
  if(verbose)
    { cat("Calculating, please be patient...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0, max = npermute, style = 3) }
  #
  for(j in 1:npermute) 
     { perm <- sample.int(n)
       for(col in 0:nd)
          { xx <- if(col == 0) 
                    z[perm,] else cbind(z[,(1:col)], z[perm,-(1:col)])
            pmod <- msir(xx, ysl, nslices = nslices, 
                         G = G, modelNames = modelNames)
            val[col + 1] <- test.statistic(pmod, col + 1)[col + 1]
          }
       val[is.na(val)] <- Inf
       count[val > obstest] <- count[val > obstest] + 1
       if(verbose) setTxtProgressBar(pbar, j)
     }
  if(verbose) { close(pbar)
                cat("\n") }
  pval <- (count)/(npermute + 1)
  ans <- as.data.frame(cbind(obstest, pval))
  rownames(ans) <- paste(0:(nt - 1), "D vs >= ", 1:nt, "D", sep = "")
  colnames(ans) <- c("Stat", "p-value")
  ans <- list(summary = ans, npermute = npermute)
  # assign ans to object in parent environment
  objname <- deparse(substitute(object))
  object$permtest <- ans
  assign(objname, value = object, 
         envir = sys.frame(sys.parent()))
  #
  return(ans)
}


##
## BIC-type criterion ----
##

msir.bic <- function(object, type = 1, plot = FALSE)
{
  svd <- svd(object$sigma, nu=0)
  inv.sqrt.sigma <- svd$v %*% diag(1/sqrt(svd$d)) %*% t(svd$v)
  M <- t(inv.sqrt.sigma) %*% object$M %*% inv.sqrt.sigma
  out <- bicDimRed(M, x = object$x, nslices = length(object$f), type = type)
  if(plot)
    { oldpar <- par(mfrow=c(1,2), mar = c(5,4,2,0.5))
      on.exit(par(oldpar))
      plot(1:length(object$evalues), object$evalues, 
           type = "b", xaxt = "n", panel.first = grid(),
           xlab = "Dimension", ylab = "Eigenvalues")
      axis(1, at = 1:length(object$evalues))       
      par(mar = c(5,0.5,2,4))
      plot(0:(length(out$crit)-1), out$crit, 
           type = "b", panel.first = grid(),
           xlab = "Dimension", xaxt = "n", yaxt = "n")
      axis(1, at = 0:(length(out$crit)-1))       
      axis(4)
      mtext("BIC-type criterion", side = 4, line = 2.5)
      points(out$d, max(out$crit), pch = 20)
    }
  # assign ans to object in parent environment
  objname <- deparse(substitute(object))
  object$bic <- out
  assign(objname, value = object, 
         envir = sys.frame(sys.parent()))
  #
  return(out)  
}

bicDimRed <- function(M, x, nslices, type = 1, tol = sqrt(.Machine$double.eps)) 
{
# BIC-type criterion for the determination of the structural dimension of
# a dimension reduction method.
# The method selects d as the maximizer of
#
#   G(d) = logL(d) - Penalty(p,d,n)
#
# Arguments:
# M = the kernel matrix
# x = the matrix of data
# nslices = number of slices used for slicing y
# type = penalty term used:
#        1)    -(p-d)*log(n)
#        2-3)  0.5*C*d*(2*p-d+1))
#            where C is defined as
#            2) (0.5*log(n) + 0.1*n^(1/3))/2 * nslices/n
#            3) log(n) * nslices/n 
#        4)    1/2*d*log(n)
#
# Returns a list with elements:
# d = the selected dimension 
# crit = the values of the criterion
# evalues = eigenvalues of the kernel matrix M
#
# Reference (penalty 2-3):
# Zhu, Miao and Peng (2006) "Sliced Inverse Regression for CDR Space
#   Estimation", JASA.
# Reference (penalty 1):
# Zhu, Zhu (2007) "On kernel method for SAVE", Journal of Multivariate Analysis
#
   M <- as.matrix(M)
   x <- as.matrix(x)
   n <- nrow(x)
   p <- ncol(x)
   h <- nslices
   # eigen decomposition
   ev <- svd(M)$d
   ev <- (ev+abs(ev))/2
   # BIC-type criterion
   ev2 <- svd(M + diag(p))$d
   ev2 <- (ev2+abs(ev2))/2
   # penalty term
   penalty <- function(d, type = 1)
   { 
     pen <- 0
     if(type == 4) pen <- 1/2*d*log(n)
     else if(type == 1) pen <- -(p-d)*log(n)
     else { if(type == 2) C <- (0.5*log(n) + 0.1*n^(1/3))/2 * h/n
            if(type == 3) C <- log(n) * h/n 
            pen <- 0.5*C*d*(2*p-d+1)
          }
     return(pen)      
   }
   k <- sum(ev > tol)
   crit <- l <- rep(0, k)
   for(d in 0:(k-1)) 
      { l[d+1] <- 0.5*n*sum((log(ev2)+1-ev2)[(1+min(k,d)):p])
        crit[d+1] <- l[d+1] - penalty(d, type) }
   dmax <- which.max(crit) - 1
   names(crit) <- paste("d=", 0:(k-1), sep="")
   # return
   ans <-list(evalues = ev[1:k], l = l[1:k], crit = crit, d = dmax)
   return(ans)
}

