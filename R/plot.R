
plot.msir <- function(x, which, 
                      type = c("pairs", "2Dplot", "spinplot", "evalues", "coefficients"), 
                      span = NULL, 
                      std = TRUE, 
                      ylab, xlab, 
                      restore.par = TRUE, ...)
{
  object <- x  # Argh.  Really want to use 'object' anyway
  y <- object$y
  numdir <- object$numdir
  if(missing(which)) which <- 1:numdir
  dir <- object$dir[,which,drop=FALSE]
  type <- match.arg(type)
  if(type != "spinplot" & isTRUE(restore.par))
    { oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar)) }
  switch(type,
   "pairs" = { if(is.numeric(span))
                 { pairs(data.frame(y, dir), gap = 0,
                         upper.panel = 
                           function(...)
                            panel.loess(col.smooth = "green3", span = span, ...),
                         ...) }
              else  pairs(data.frame(y, dir), gap = 0, ...) },
   "2Dplot" = { plot(dir[,1], y, 
                     xlab = ifelse(missing(xlab), paste("MSIR Dir", which[1], sep=""), xlab), 
                     ylab = ifelse(missing(ylab), "Y", ylab), ...) 
                if(is.numeric(span))
                  { lines(l <- loess.sd(dir[,1], y, k = 1, span = span))
                    lines(l$x, l$upper, lty=2)
                    lines(l$x, l$lower, lty=2) } },
   "spinplot" = { data <- data.frame(dir[,1], y, dir[,2])
                  colnames(data) <- c(colnames(dir)[1], "y", colnames(dir)[2])
                  if(!missing(ylab)) colnames(data)[2] <- ylab
                  if(!missing(xlab)) colnames(data)[c(1,3)] <- xlab
                  spinplot(data, fit.smooth = is.numeric(span), span = span, ...)
                },
   "evalues" = { bic <- msir.bic(object, plot = FALSE)
                 oldpar <- par(mfrow=c(1,2), mar = c(5,4,2,0.5))
                 on.exit(par(oldpar))
                 plot(1:length(object$evalues), object$evalues, 
                      type = "b", xaxt = "n", panel.first = grid(),
                      xlab = "Dimension", ylab = "Eigenvalues")
                 axis(1, at = 1:length(object$evalues))
                 if(!is.null(object$permtest))
                   { pval <- object$permtest$summary[,2]
                     pvalcol <- cut(pval, c(0,0.05,0.10,1), 
                                    include.lowest = TRUE,
                                    labels = c("black", "grey70", "white"))
                     points(1:length(object$evalues), object$evalues,
                            pch = 21, bg = paste(pvalcol))
                     legend("topright", inset = 0.01, 
                            title = "p-value",
                            pch = 21, pt.bg = paste(levels(pvalcol)),
                            legend = c("[0-0.05)", "[0.05-0.10)", "[0.10-1]"))
                     }
                 par(mar = c(5,0.5,2,4))
                 plot(0:(length(bic$crit)-1), bic$crit-max(bic$crit), 
                         type = "b", xlab = "Dimension", 
                         xaxt = "n", yaxt = "n",
                         panel.first = abline(h = -c(0,2,6,10), 
                                              v = 0:(length(bic$crit)-1),
                                              lty = 3, col = "grey"))
                 axis(1, at = 0:(length(bic$crit)-1))       
                 axis(4)
                 mtext(expression(paste(Delta, "BIC-type criterion", sep="")),
                       side = 4, line = 2.5)
                 points(bic$d, 0, pch = 20)
               },
    "coefficients" = { if(missing(which)) which <- 1:object$numdir
                       if(std) B <- object$std.basis[,which,drop=FALSE]
                       else    B <- object$basis[,which,drop=FALSE]
                       evalues <- object$evalues[which]
      ev  <- signif(evalues, 4)
      evp <- round(evalues/sum(object$evalues[1:numdir])*100,2)
      if(length(evalues) != ncol(B))
        stop("number of columns of coef must be equal to num. of eigenvalues!")
      w <- max(strwidth(rownames(B), units="inches"))*1.1
      ww <- c(w,rep((par("fin")[1]-w)/ncol(B),ncol(B)))
      ww <- pmax(0.1, ww/sum(ww))
      layout(matrix(1:(ncol(B)+1), 1, ncol(B)+1), 
             widths = ww, heights = rep(1,ncol(B)+1))
      mar <- par()$mar
      par(mar = c(mar[1], 0, mar[3], 1))
      plot(0, 0, type="n", ylim=c(1,nrow(B)), 
           xlab="", ylab="", axes=FALSE)
      mtext(rownames(B), side=4, at=nrow(B):1, las=1, adj=1, cex=0.8)
      par(mar = c(mar[1], 1, mar[3], 1))
      for(j in 1:ncol(B))
         { plot(0, 0, type="n", 
                xlim = range(0,B[,j])*1.1,
                ylim = c(1,nrow(B)),
                xlab ="", ylab="", yaxt="n")
            mtext(colnames(B)[j], side=3, line=1)
            mtext(paste(ev[j], " (", evp[j], "%)", sep=""), 
                  side=1, line=3, cex = 0.8)
            for(i in 1:nrow(B))
               lines(c(0,B[i,j]), rep(nrow(B)+1-i,2), lwd=2)
            abline(v=0) 
         }
     }
  )
  invisible()
}
