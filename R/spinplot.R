##
## Spin-plot ----
##

spinplot <- function(x, y, z, 
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
                     ...)
{
  if(!requireNamespace("rgl", quietly = TRUE))
    return(warning("'rgl' packages is needed by spinplot."))
    
  if(missing(y) & missing(z))
    { X <- as.matrix(x)
      if(ncol(X) != 3)
        stop("x must be a 3-columns matrix or data.frame!") 
      varnames <- colnames(X) 
    }
  else
   { varnames <- c(deparse(substitute(x)),
                   deparse(substitute(y)),
                   deparse(substitute(z)))
     X <- cbind(x, y, z) 
     colnames(X) <- varnames 
   }
  if(is.null(varnames)) 
     varnames <- c("H", "V", "O")
  #
  if(rem.lin.trend)
    { e <- lm(X[,2]~X[,c(1,3)])$residuals
      X[,2] <- e
      varnames[2] <- paste("e(", varnames[2], "|", varnames[1], 
                           ",", varnames[3], ")", sep="")
    }
  #
  if(uncor.vars)
    { e <- lm(X[,3]~X[,1])$residuals
      X[,3] <- e
      varnames[3] <- paste("e(", varnames[3], "|", varnames[1], ")", sep="")
    }
  #
  scaling <- match.arg(scaling)
  if(scaling == "abc")
    { # abc-scaling
      X <- apply(X, 2, function(x)
                       2*(x-sum(range(x))/2)/diff(range(x))) }
  else
    { # aaa-scaling
      a <- min(apply(X, 2, function(x) 2/diff(range(x))))
      X <- apply(X, 2, function(x) a*(x-sum(range(x))/2)) }
  ax <- apply(X, 2, function(x) c(min(x), max(x)))
  #
  n <- nrow(X) 

  if(missing(markby)) markby <- rep(1,n)
  u <- sort(unique(markby))
  nu <- length(u)
  # specify pch.points recycling if necessary
  if(missing(pch.points)) pch.points <- 1:nu
  pch.points <- rep(pch.points, nu/length(pch.points)+1)[1:nu]
  # specify col.points recycling if necessary
  if(missing(col.points)) col.points <- 1:nu 
  if(is.numeric(col.points))
    col.points <- grDevices::palette()[col.points]
  col.points <- rep(col.points, nu/length(col.points)+1)[1:nu]

  # setup rgl graphical window
  if(rgl::rgl.cur() > 0) 
    rgl::rgl.set(rgl::rgl.cur()) else rgl::rgl.open()
  rgl::rgl.bg(sphere = FALSE, fogtype = "none", lit = TRUE, 
              back = "lines", alpha = 1,
              color = rep(background,2))
  rgl::material3d("point_antialias" = TRUE)
  rgl::rgl.pop("lights")
  rgl::rgl.light(ambient = "black", diffuse = "black", specular = "black")

  # draw the spinplot
  rgl::plot3d(X, type = "n",
              xlab = "", ylab = "", zlab = "",
              box = FALSE, axes = FALSE)
  # draw axis
  alen <- 0.05; awid <- 0.01
  bbox <- rgl::par3d("bbox")
  diffs <- bbox[c(2,4,6)]-bbox[c(1,3,5)]
  rgl::lines3d(c(0,ax[2,1]), c(0,0), c(0,0), color = col.axis)
  rgl::triangles3d(ax[2,1]+c(0,1,0)*alen*diffs[1],
                   0+c(-1,0,1)*awid*diffs[2],
                   0+c(-1,0,1)*awid*diffs[3],lit=FALSE,
                   color = col.axis)
  rgl::lines3d(c(0,0), c(0,ax[2,2]), c(0,0), color = col.axis)
  rgl::triangles3d(0+c(-1,0,1)*awid*diffs[1],
                   ax[2,2]+c(0,1,0)*alen*diffs[2],
                   0+c(-1,0,1)*awid*diffs[3],lit=FALSE,
                   color = col.axis)
  rgl::lines3d(c(0,0), c(0,0), c(0,ax[2,3]), color = col.axis)
  rgl::triangles3d(0+c(-1,0,1)*awid*diffs[1],
                   0+c(-1,0,1)*awid*diffs[2],
                   ax[2,3]+c(0,1,0)*alen*diffs[3],lit=FALSE,
                   color = col.axis)
  # axis labels
  coordtext <- function(i) max(ax[,i])+0.1*diff(ax[,i])
  rgl::rgl.texts(c(coordtext(1),0,0), 
                 c(0,coordtext(2),0), 
                 c(0,0,coordtext(3)), 
                 text = varnames,
                 adj = 0.5, 
                 cex = rgl::par3d("cex")*2/3,
                 color = col.axis)
  # draw points
  for(j in 1:nu)
  {
    i <- which(markby == u[j])
    rgl::pch3d(X[i,], pch = pch.points[j], color = col.points[j],
               cex = rgl::par3d("cex")*cex.points*0.25)
  }
  #
  if(fit.ols)
    { 
      mod <- lm(y ~ x + z, 
                data = data.frame(x = X[,1], y = X[,2], z = X[,3]))
      xgrid <- seq(-1,1,length=ngrid)
      zgrid <- seq(-1,1,length=ngrid)
      pred.grid <- expand.grid(x = xgrid, z = zgrid)
      pred <- matrix(predict(mod, pred.grid), ngrid, ngrid)
      rgl::rgl.surface(xgrid, zgrid, pred, 
                       alpha = 0.5, lit = FALSE, 
                       color = col.ols, 
                       front = "lines", back = "lines")
    }
  #
  if(fit.smooth)
    { 
      mod <- loess(y ~ x*z, span = span, 
                   data = data.frame(x = X[,1], y = X[,2], z = X[,3]))
      xgrid <- seq(-1,1,length=ngrid)
      zgrid <- seq(-1,1,length=ngrid)
      pred.grid <- expand.grid(x = xgrid, z = zgrid)
      pred <- predict(mod, pred.grid)
      rgl::rgl.surface(xgrid, zgrid, pred, 
                       alpha = 0.5, lit = FALSE, 
                       color = col.smooth, 
                       front = "lines", back = "lines")
    }

  # set initial view
  rgl::rgl.viewpoint(theta = 10, phi = 15, fov = 1, zoom = 0.5)
  rgl::par3d("windowRect" = c(100,100,500,500))
  rgl::rgl.bringtotop()
  
  invisible()  
}


 
