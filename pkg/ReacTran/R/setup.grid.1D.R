
##==============================================================================
## Creation of a one-dimensional finite difference grid
##==============================================================================

setup.grid.1D <- function(x.up=0,	x.down=NULL, L=NULL,
  N=NULL, dx.1 =NULL, p.dx.1 = rep(1,length(L)),
  max.dx.1 = 10*dx.1, dx.N =NULL, p.dx.N = rep(1,length(L)),
  max.dx.N = 10*dx.N) {

## Check on the input

  if (is.null(x.down[1]) && is.null(L[1]))
    stop ("Cannot create grid: the length or end point is not specified")

  if (is.null(x.down[1])) {
    x.down <- cumsum(L)
#### KS ####
#	  x.down <- vector(length=length(L))
#    x.check <- x.up
#	  for (i in 1:length(L)) {
#      if (L[i] < 0)
#        stop ("Cannot create grid: L < 0")
#      x.down[i] <- x.check + L[i]
#      x.check <- x.down[i]
#    }
  }
  
  if (is.null(L[1])) {
    L <- diff(c(0,x.down))
#### KS ####
#  	L <- vector(length=length(x.down))
#    x.check <- x.up
#    for (i in 1:length(x.down)) {
#   	  if (x.down[i] < x.check)
#         stop ("Cannot create grid: x.down < x.up")
#	    L[i] <- x.down[i] - x.check
#      x.check <- x.down[i]
#    }
  }

## Calculation of the grid cell sizes

  dx <- vector()

## loop over all zones
  for (i in 1:length(L)) {
    if (is.null(N[i]) && is.null(dx.1[i]))
		  stop ("Cannot create grid: neither N (the number of boxes) nor dx.1 (the size of the first grid cell) is specified")

    if ((is.null(dx.1[i])) && (is.null(dx.N[i])))  { # all grid cells are equal

      if (N[i] < 1)
        stop ("Cannot create grid: N < 1")
      A.dx  <- rep(L[i]/N[i],N[i])
    }
    else if (is.null(dx.N[i])) { # dx.1 has a value
  ## KS Ik denk dat dit niet kan - wel eventueel testen of niet = 0?##
      if (is.null(dx.1[i]))
        stop ("Cannot create grid: dx.1 not specified")
      if (dx.1[i] > L[i])
        stop ("Cannot create grid: dx.1 > L")
      # use gradual increase of grid cell size at upper interface alone
      A.dx <- vector()
      pos <- vector()
      A.dx[1] <- dx.1[i]
      pos[1] <- dx.1[i]
      j <- 1

      while ((is.null(N[i])&&(pos[j] < L[i])) || (!is.null(N[i])&&(j<N[i]))) {
        j <- j + 1
        A.dx[j] <- min(max.dx.1[i],A.dx[j-1]*p.dx.1[i])
        pos[j] <- pos[j-1] + A.dx[j]
      }

      A.dx <- A.dx*(L[i]/pos[length(pos)])

    } else {

      if (is.null(dx.1[i]))
        stop ("Cannot create grid: dx.1 not specified")
        # use gradual increase of grid cell size at both interfaces
      A.dx  <- vector()
      pos <- vector()
      A.dx[1] <- dx.1[i]
      pos[1] <- dx.1[i]
      j <- 1
      if (!is.null(N[i]))
        N.half <- as.integer(N[i]/2)

      while ((is.null(N[i])&&(pos[j] < 0.5*L[i])) ||
             (!is.null(N[i])&&(j<N.half))) {
         j <- j + 1
         A.dx[j] <- min(max.dx.1[i],A.dx[j-1]*p.dx.1[i])
         pos[j] <- pos[j-1] + A.dx[j]
      }
      A.dx <- A.dx*(0.5*L[i]/pos[length(pos)])

      B.dx  <- vector()
      pos <- vector()
      B.dx[1] <- dx.N[i]
      pos[1] <- dx.N[i]
      j <- 1
      while ((is.null(N[i])&&(pos[j] < 0.5*L[i])) ||
            (!is.null(N[i])&&(j<(N[i]-N.half)))) {
        j <- j + 1
        B.dx[j] <- min(max.dx.N[i],B.dx[j-1]*p.dx.N[i])
        pos[j] <- pos[j-1] + B.dx[j]
      }
      B.dx <- B.dx*(0.5*L[i]/pos[length(pos)])

      A.dx <- c(A.dx,B.dx[length(B.dx):1])
    }
    dx <- c(dx,A.dx)
  }
  ## positions and distances

  x.int  <- x.up + diffinv(dx)
  x.mid  <- colMeans(rbind(x.int[-1],x.int[-(length(dx)+1)]))
  dx.aux <- c(dx[1]/2,diff(x.mid),dx[length(dx)]/2)

  ## Packaging of results

  Res <- list(x.up = x.up,
	  			    x.down = x.down[length(x.down)],
		  		    x.mid = x.mid,   # position of centre of the grid cells, vector of length N
              x.int  = x.int,  # position of the grid cell interfaces , vector of length N+1
              dx = dx  ,       # thickness of the grid cells , vector length N
              dx.aux = dx.aux, # auxiliary vector with distances between centre of adjacent cells, first and last: half of cell size, vector of length N+1
              N = length(dx))  # total number of grid cells
  class(Res) <- "grid.1D"
  return(Res)
}

##==============================================================================
## S3 method: Plotting of a one-dimensional finite difference grid
##==============================================================================

plot.grid.1D <- function(x,...) {
  mf <- par(mfrow=c(2,1))
  on.exit(par(mf))
    
  plot(x$x.int,main="position of cells",ylab="x",xlab="index",...)
  points(0.5+(1:x$N),x$x.mid,pch=16)
  legend("topleft",c("x.int","x.mid"),pch=c(1,16))
  plot(x$dx.aux,main="box thickness",ylab="dx",xlab="index",...)
  points(0.5+(1:x$N),x$dx,pch=16)
  legend("topleft",c("dx.aux","dx.mid"),pch=c(1,16))
}

