
##==============================================================================
## 2D Transport of a solute in the sediment
##==============================================================================

tran.2D <- function(C, C.x.up=C[1,], C.x.down=C[nrow(C),],
  C.y.up=C[,1],  C.y.down=C[,ncol(C)],
  flux.x.up=NULL, flux.x.down=NULL, flux.y.up=NULL, flux.y.down=NULL,
  a.bl.x.up=NULL, C.bl.x.up=NULL, a.bl.x.down=NULL, C.bl.x.down=NULL,
  a.bl.y.up=NULL, C.bl.y.up=NULL, a.bl.y.down=NULL, C.bl.y.down=NULL,
  D.grid=NULL, D.x=NULL, D.y=D.x, v.grid=NULL, v.x=0, v.y=0, 
  AFDW.grid=NULL, AFDW.x=1, AFDW.y=AFDW.x,
  VF.grid=NULL,VF.x=1, VF.y=VF.x, A.grid=NULL, 
  A.x=1, A.y=1, dx=NULL, dy=NULL, grid=NULL,
  full.check = FALSE, full.output = FALSE)
											
{

  N <- nrow(C)
  M <- ncol(C)

# DEFAULT INFILLING OF GRID PARAMETERS

# infilling of 2D grid
  if (is.null(grid)) {
     DX    <-  if (is.list(dx)) dx$dx else rep(dx,length.out=N)
     DXaux <-  if (is.list(dx)) dx$dx.aux else 0.5*(c(0,rep(dx,length.out=N))+
                                  c(rep(dx,length.out=N),0))
     DY    <-  if (is.list(dy)) dy$dx else rep(dy,length.out=M)
     DYaux <-  if (is.list(dy)) dy$dx.aux else 0.5*(c(0,rep(dy,length.out=M))+
                                  c(rep(dy,length.out=M),0))

    grid <- list(
      dx    = DX,  dx.aux= DXaux,
      dy    = DY,  dy.aux= DYaux
                 )
   }
# infilling of VF.grid
  if (is.null(VF.grid)) {
    VF.grid <- list()
    # infilling of x-matrix    
    if (is.matrix(VF.x)) {
      if (dim(VF.x) != c(N+1,M))
      stop ("error: VF.x matrix does not have correct (N+1)xM dimensions" )
      VF.grid$x.int <- VF.x
      VF.grid$x.mid <- 0.5*(VF.x[1:N,]+VF.x[2:(N+1),])
    } else if (class(VF.x)=="prop.1D") {
      VF.grid$x.int <- matrix(data=VF.x$int,nrow=(N+1),ncol=M)
      VF.grid$x.mid <- matrix(data=VF.x$mid, nrow=N, ncol=M)
    } else if (length(VF.x) == 1) {
      VF.grid$x.int <- matrix(data=VF.x,nrow=(N+1),ncol=M)
      VF.grid$x.mid <- matrix(data=VF.x,nrow=N,ncol=M)
    } else if (length(VF.x) != N+1) {
      stop("error: VF.x should be a vector of length 1 or N+1")
    } else {  # correct length
      VF.grid$x.int <- matrix(data=VF.x,nrow=(N+1),ncol=M)
      VF.grid$x.mid <- matrix(data=0.5*(VF.x[1:N]  +VF.x[2:(N+1)]),
                     nrow=N, ncol=M)
    }
    # infilling of y-matrix    
    if (is.matrix(VF.y)) {
      if (dim(VF.y) != c(N,M+1))
        stop ("error: VF.y matrix does not have the correct Nx(M+1) dimensions")
      VF.grid$y.int <- VF.y
      VF.grid$y.mid <- 0.5*(VF.y[,1:M]+VF.y[,2:(M+1)])
    } else if (class(VF.y)=="prop.1D") {
      VF.grid$y.int <- matrix(data=VF.y$int,nrow=N,ncol=(M+1))
      VF.grid$y.mid <- matrix(data=VF.y$mid, nrow=N, ncol=M)
    } else if (length(VF.y) == 1) {
      VF.grid$y.int <- matrix(data=VF.y,nrow=N,ncol=(M+1))
      VF.grid$y.mid <- matrix(data=VF.y,nrow=N,ncol=M)
    } else if (length(VF.y) != M+1) {
      stop("error: VF.y should be a vector of length 1 or M+1")
    } else {  # correct length
      VF.grid$y.int <- matrix(data=VF.y,nrow=N,ncol=(M+1))
      VF.grid$y.mid <- matrix(data=0.5*(VF.y[1:N]  +VF.y[2:(N+1)]),
                     nrow=N, ncol=M)
    }
  }  

# infilling of A.grid
  if (is.null(A.grid)) {
    A.grid <- list()
    # infilling of x-matrix    
    if (is.matrix(A.x)) {
      if (sum(abs(dim(A.x) - c(N+1,M)))!=0)
        stop ("error: A.x matrix not of correct dimensions" )
      A.grid$x.int <- A.x
      A.grid$x.mid <- 0.5*(A.x[1:N,]+A.x[2:(N+1),])
    } else if (class(A.x)=="prop.1D") {
      A.grid$x.int <- matrix(data=A.x$int,nrow=(N+1),ncol=M)
      A.grid$x.mid <- matrix(data=A.x$mid, nrow=N, ncol=M)
    } else if (length(A.x) == 1) {
      A.grid$x.int <- matrix(data=A.x,nrow=(N+1),ncol=M)
      A.grid$x.mid <- matrix(data=A.x,nrow=N,ncol=M)
    } else if (length(A.x) != N+1) {
      stop("error: A.x should be a vector of length 1 or N+1")
    } else {  # correct length
      A.grid$x.int <- matrix(data=A.x,nrow=(N+1),ncol=M)
      A.grid$x.mid <- matrix(data=0.5*(A.x[1:N]  +A.x[2:(N+1)]),
                      nrow=N, ncol=M)
    }
    # infilling of y-matrix    
    if (is.matrix(A.y)) {
      if (sum(abs(dim(A.y) - c(N,M+1)))!=0)
        stop ("error: A.y matrix not of correct dimensions")
      A.grid$y.int <- A.y
      A.grid$y.mid <- 0.5*(A.y[,1:M]+A.y[,2:(M+1)])
    } else if (class(A.y)=="prop.1D") {
      A.grid$y.int <- matrix(data=A.y$int,nrow=N,ncol=(M+1))
      A.grid$y.mid <- matrix(data=A.y$mid,nrow=N,ncol=M)
    } else if (length(A.y) == 1) {
      A.grid$y.int <- matrix(data=A.y,nrow=N,ncol=(M+1))
      A.grid$y.mid <- matrix(data=A.y,nrow=N,ncol=M)
    } else if (length(A.y) != M+1) {
       stop("error: A.y should be a vector of length 1 or M+1")
    } else {  # correct length
      A.grid$y.int <- matrix(data=A.y,nrow=N,ncol=(M+1))
      A.grid$y.mid <- matrix(data=0.5*(A.y[1:N]  +A.y[2:(N+1)]),
                      nrow=N, ncol=M)
    }
  }
  
# infilling of AFDW.grid
  if (is.null(AFDW.grid)) {
    AFDW.grid <- list()
    if (class(AFDW.x)=="prop.1D") {
      AFDW.grid$x.int <- matrix(data=AFDW.x$int,nrow=(N+1),ncol=M)
    } else if (is.matrix(AFDW.x) || is.vector(AFDW.x)) { 
      AFDW.grid$x.int <- matrix(data=AFDW.x,nrow=(N+1),ncol=M)
    }  
    if (class(AFDW.y)=="prop.1D") {
      AFDW.grid$y.int <- matrix(data=AFDW.y$int,nrow=N,ncol=(M+1))
    } else if (is.matrix(AFDW.x) || is.vector(AFDW.x)) { 
      AFDW.grid$y.int <- matrix(data=AFDW.y,nrow=N,ncol=(M+1))
    }  
  }  

# infilling of D.grid
  if (is.null(D.grid)) {
    D.grid <- list()
    if (class(D.x)=="prop.1D") {
      D.grid$x.int <- matrix(data=D.x$int,nrow=(N+1),ncol=M)
    } else if (is.matrix(D.x) || is.vector(D.x)) { 
      D.grid$x.int <- matrix(data=D.x,nrow=(N+1),ncol=M)
    }  
    if (class(D.y)=="prop.1D") {
      D.grid$y.int <- matrix(data=D.y$int,nrow=N,ncol=(M+1))
    } else if (is.matrix(D.x) || is.vector(D.x)) { 
      D.grid$y.int <- matrix(data=D.y,nrow=N,ncol=(M+1))
    }  
  }  

# infilling of v.grid
  if (is.null(v.grid)) {
    v.grid <- list()
    if (class(v.x)=="prop.1D") {
      v.grid$x.int <- matrix(data=v.x$int,nrow=(N+1),ncol=M)
    } else if (is.matrix(v.x) || is.vector(v.x)) { 
      v.grid$x.int <- matrix(data=v.x,nrow=(N+1),ncol=M)
    }  
    if (class(v.y)=="prop.1D") {
      v.grid$y.int <- matrix(data=v.y$int,nrow=N,ncol=(M+1))
    } else if (is.matrix(v.x) || is.vector(v.x)) { 
      v.grid$y.int <- matrix(data=v.y,nrow=N,ncol=(M+1))
    }  
  }  


#==============================================================================
# INPUT CHECKS  
#==============================================================================


  if (full.check) {

## check dimensions of input concentrations

    if (!is.null(C.x.up)) {
      if (!((length(C.x.up)==1) || (length(C.x.up)==(M))))
        stop("error: C.x.up should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(C.x.down)) {
      if (!((length(C.x.down)==1) || (length(C.x.down)==(M))))
        stop("error: C.x.down should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(C.y.up)) {
      if (!((length(C.y.up)==1) || (length(C.y.up)==(N))))
        stop("error: C.y.up should be a vector of length 1 or nrow(C)")
    }
    if (!is.null(C.y.down)) {
      if (!((length(C.y.down)==1) || (length(C.y.down)==(N))))
        stop("error: C.y.down should be a vector of length 1 or nrow(C)")
    }

    if (!is.null(C.bl.x.up)) {
      if (!((length(C.bl.x.up)==1) || (length(C.bl.x.up)==(M))))
        stop("error: C.bl.x.up should be a vector of length 1 or ncol(C)")
   }

    if (!is.null(C.bl.x.down)) {
      if (!((length(C.bl.x.down)==1) || (length(C.bl.x.down)==(M))))
        stop("error: C.bl.x.down should be a vector of length 1 or ncol(C)")
    }

    if (!is.null(C.bl.y.up)) {
      if (!((length(C.bl.y.up)==1) || (length(C.bl.y.up)==(N))))
        stop("error: C.bl.y.up should be a vector of length 1 or nrow(C)")
    }

    if (!is.null(C.bl.y.down)) {
      if (!((length(C.bl.y.down)==1) || (length(C.bl.y.down)==(N))))
        stop("error: C.bl.y.down should be a vector of length 1 or nrow(C)")
    }

# check dimensions of input fluxes

    if (!is.null(flux.x.up)) {
      if (!((length(flux.x.up)==1) || (length(flux.x.up)==(M))))
        stop("error: flux.x.up should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(flux.x.down)) {
      if (!((length(flux.x.down)==1) || (length(flux.x.down)==(M))))
        stop("error: flux.x.down should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(flux.y.up)) {
      if (!((length(flux.y.up)==1) || (length(flux.y.up)==(N))))
        stop("error: flux.y.up should be a vector of length 1 or nrow(C)")
    }

    if (!is.null(flux.y.down)) {
      if (!((length(flux.y.down)==1) || (length(flux.y.down)==(N))))
        stop("error: flux.y.down should be a vector of length 1 or nrow(C)")
    }


## check input of grid

    if (is.null(dx) && is.null(dy) && is.null(grid))
      stop("error: dx, dy, and grid cannot be NULL at the same time")

    gn <- names(grid)
    if (! "dx" %in% gn)
      stop("error: grid should be a list that contains 'dx' ")
    if (! "dx.aux" %in% gn)
    	stop("error: grid should be a list that contains 'dx.aux' ")
    if (! "dy" %in% gn)
      stop("error: grid should be a list that contains 'dy' ")
    if (! "dy.aux" %in% gn)
    	stop("error: grid should be a list that contains 'dy.aux' ")
    if (is.null(grid$dx) || is.null(grid$dx.aux))
    	stop("error: the grid should be a list with (numeric) values for 'dx' and 'dx.aux' ")
    if (is.null(grid$dy) || is.null(grid$dy.aux))
    	stop("error: the grid should be a list with (numeric) values for 'dy' and 'dy.aux' ")
    if (any(grid$dx <= 0) || any(grid$dx.aux <= 0) )
    	stop("error: the grid distances dx and dx.aux should always be positive")
    if (any(grid$dy <= 0) || any(grid$dy.aux <= 0) )
    	stop("error: the grid distances dy and dy.aux should always be positive")

## check input of AFDW.grid

    if (is.null(AFDW.x) && is.null(AFDW.y) && is.null(AFDW.grid))
      stop("error: AFDW.x, AFDW.y, and AFDW.grid cannot be NULL at the same time")

    gn <- names(AFDW.grid)
    if (! "x.int" %in% gn)
      stop("error: AFDW.grid should be a list that contains 'x.int', the AFDW values at the vertical interfaces of the grid cells")
    if (! "y.int" %in% gn)
      stop("error: AFDW.grid should be a list that contains 'y.int', the AFDW values at the horizontal interfaces of the grid cells")
    if (is.null(AFDW.grid$x.int))
      stop("error: AFDW.grid$x.int should be a list with (numeric) values")
    if (is.null(AFDW.grid$y.int))
      stop("error: AFDW.grid$y.int should be a list with (numeric) values")
    if (any (AFDW.grid$x.int < 0)||any (AFDW.grid$x.int > 1))
    	stop("error: the AFDW should range between 0 and 1")
    if (any (AFDW.grid$y.int < 0)||any (AFDW.grid$y.int > 1))
	    stop("error: the AFDW should range between 0 and 1")

## check input of D.grid

    if (is.null(D.x) && is.null(D.y) && is.null(D.grid))
      stop("error: D.x, D.y, and D.grid cannot be NULL at the same time")

    gn <- names(D.grid)
    if (! "x.int" %in% gn)
      stop("error: D.grid should be a list that contains 'x.int', the D values at the vertical interfaces of the grid cells")
    if (! "y.int" %in% gn)
      stop("error: D.grid should be a list that contains 'y.int', the D values at the horizontal interfaces of the grid cells")
    if (is.null(D.grid$x.int))
      stop("error: D.grid$x.int should be a list with (numeric) values")
    if (is.null(D.grid$y.int))
      stop("error: D.grid$y.int should be a list with (numeric) values")
    if (any (D.grid$x.int < 0)||any (D.grid$y.int < 0))
    	stop("error: the diffusion coefficient should always be positive")

## check input of v.grid

    if (is.null(v.x) && is.null(v.y) && is.null(v.grid))
      stop("error: v.x, v.y, and v.grid cannot be NULL at the same time")

    gn <- names(v.grid)
    if (! "x.int" %in% gn)
      stop("error: v.grid should be a list that contains 'x.int', the velocity values at the vertical interfaces of the grid cells")
    if (! "y.int" %in% gn)
      stop("error: v.grid should be a list that contains 'y.int', the velocity values at the horizontal interfaces of the grid cells")
    if (is.null(v.grid$x.int))
      stop("error: the advective velocity v.grid$x.int should be a list with (numeric) values")
    if (is.null(v.grid$y.int))
      stop("error: the advective velocity v.grid$y.int should be a list with (numeric) values")

## check input of VF.grid

    gn <- names(VF.grid)
    if (! "x.int" %in% gn)
      stop("error: VF.grid should be a list that contains 'x.int'")
    if (! "y.int" %in% gn)
      stop("error: VF.grid should be a list that contains 'y.int'")
    if (! "x.mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'x.mid'")
    if (! "y.mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'y.mid'")
    if (is.null(VF.grid$x.int) || is.null(VF.grid$y.int) || is.null(VF.grid$x.mid) || is.null(VF.grid$y.mid))
     stop("error: VF should contain (numeric) values")
    if (any (VF.grid$x.int < 0) || any (VF.grid$y.int < 0) || any (VF.grid$x.mid < 0) || any (VF.grid$y.mid < 0))
      stop("error: the VF values should always be positive")

## check input of A.grid
    gn <- names(A.grid)
    if (! "x.int" %in% gn)
      stop("error: A.grid should be a list that contains 'x.int'")
    if (! "y.int" %in% gn)
      stop("error: A.grid should be a list that contains 'y.int'")
    if (! "x.mid" %in% gn)
      stop("error: A.grid should be a list that contains 'x.mid'")
    if (! "y.mid" %in% gn)
      stop("error: A.grid should be a list that contains 'y.mid'")
    if (is.null(A.grid$x.int) || is.null(A.grid$y.int) || is.null(A.grid$x.mid) || is.null(A.grid$y.mid))
     stop("error: the VF should contain (numeric) values")
    if (any (A.grid$x.int < 0) || any (A.grid$y.int < 0) || any (A.grid$x.mid < 0) || any (A.grid$y.mid < 0))
      stop("error: the A values should always be positive")

  }
## FUNCTION BODY: CALCULATIONS

## Impose boundary flux at upstream x-boundary when needed
## Default boundary condition is no gradient
  if (! is.null (flux.x.up[1])) {
    nom <- flux.x.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1] +
           (1-AFDW.grid$x.int[1,])*v.grid$x.int[1,])*C[1,]
    denom <- VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1]+
             AFDW.grid$x.int[1,]*v.grid$x.int[1,])
    C.x.up <- nom/denom
  }

## Impose boundary flux at downstream x-boundary when needed
## Default boundary condition is no gradient
  if (! is.null (flux.x.down[1])) {
  	nom <- flux.x.down - VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/
            grid$dx.aux[N+1] + AFDW.grid$x.int[(N+1),]*v.grid$x.int[(N+1),])*C[N,]
    denom <- -VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/grid$dx.aux[N+1]+
            (1-AFDW.grid$x.int[(N+1),])*v.grid$x.int[(N+1),])
    C.x.down <- nom/denom
  }

# Impose boundary flux at upstream y-boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.y.up[1])) {
    nom <- flux.y.up + VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1] +
           (1-AFDW.grid$y.int[,1])*v.grid$y.int[,1])*C[,1]
    denom <- VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1]+
             AFDW.grid$y.int[,1]*v.grid$y.int[,1])
    C.y.up <- nom/denom
  }

# Impose boundary flux at downstream y-boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.y.down[1]))  {
	  nom <- flux.y.down - VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/
           grid$dy.aux[M+1] + AFDW.grid$y.int[,(M+1)]*v.grid$y.int[,(M+1)])*C[,M]
    denom <- -VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/grid$dy.aux[M+1]+
             (1-AFDW.grid$y.int[,(M+1)])*v.grid$y.int[,(M+1)])
    C.y.down <- nom/denom
  }

## when upper boundary layer is present, calculate new C.x.up
  if (!is.null(a.bl.x.up) & !is.null(C.bl.x.up[1])) {
	  nom <- a.bl.x.up*C.bl.x.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/
           grid$dx.aux[1] + (1-AFDW.grid$x.int[1,])*v.grid$x.int[1,])*C[1,]
    denom <- a.bl.x.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1]+
             AFDW.grid$x.int[1,]*v.grid$x.int[1,])
	  C.x.up <- nom/denom
  }

## when lower boundary layer is present, calculate new C.x.down
  if (!is.null(a.bl.x.down) & !is.null(C.bl.x.down[1])) {
	  nom <- a.bl.x.down*C.bl.x.down + VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/
           grid$dx.aux[(N+1)] + (1-AFDW.grid$x.int[(N+1),])*
           v.grid$x.int[(N+1),])*C[N,]
    denom <- a.bl.x.down + VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/
             grid$dx.aux[(N+1)]+ AFDW.grid$x.int[(N+1),]*v.grid$x.int[(N+1),])
	  C.x.down <- nom/denom
  }

## when left boundary layer is present, calculate new C.y.up
  if (!is.null(a.bl.y.up) & !is.null(C.bl.y.up[1])) {
	  nom <- a.bl.y.up*C.bl.y.up + VF.grid$y.int[,1]*(D.grid$y.int[,1]/
           grid$dy.aux[1] + (1-AFDW.grid$y.int[,1])*v.grid$y.int[,1])*C[,1]
    denom <- a.bl.y.up + VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1]+
             AFDW.grid$y.int[,1]*v.grid$y.int[,1])
	  C.y.up <- nom/denom
  }

## when right boundary layer is present, calculate new C.y.down
  if (!is.null(a.bl.y.down) & !is.null(C.bl.y.down[1]))   {
	  nom <- a.bl.y.down*C.bl.y.down + VF.grid$y.int[,(M+1)]*
           (D.grid$y.int[,(M+1)]/grid$dy.aux[(M+1)] +
           (1-AFDW.grid$y.int[,(M+1)])*v.grid$y.int[,(M+1)])*C[,M]
    denom <- a.bl.y.down + VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/
             grid$dy.aux[(M+1)]+ AFDW.grid$y.int[,(M+1)]*v.grid$y.int[,(M+1)])
	  C.y.down <- nom/denom
  }

## Calculate diffusive part of the flux
  x.Dif.flux <- as.matrix(-VF.grid$x.int * D.grid$x.int *
                diff(rbind(C.x.up, C, C.x.down, deparse.level = 0))/
                matrix(data=grid$dx.aux,nrow=(N+1),ncol=M,byrow=FALSE))
  y.Dif.flux <- as.matrix(-VF.grid$y.int * D.grid$y.int *
                t(diff(t(cbind(C.y.up,C,C.y.down,deparse.level = 0))))/
                matrix(data=grid$dy.aux,nrow=N,ncol=(M+1),byrow=TRUE))

## Calculate advective part of the flux - NOG NEGATIEVE FLUXEN
  x.Adv.flux <- 0
  
  if (any(v.grid$x.int >0) ) {
    vv <- v.grid$x.int
    vv[vv<0]<-0
    x.Adv.flux <-  x.Adv.flux + as.matrix(VF.grid$x.int * vv * (
                 (1-AFDW.grid$x.int) * rbind(C.x.up,C,deparse.level = 0)
                 + AFDW.grid$x.int * rbind(C,C.x.down,deparse.level = 0)))
  }
  if (any (v.grid$x.int < 0))  {
    vv <- v.grid$x.int
    vv[vv>0]<-0
    x.Adv.flux <-  x.Adv.flux + as.matrix(VF.grid$x.int * vv * (
                    AFDW.grid$x.int * rbind(C.x.up,C,deparse.level = 0)
                 + (1-AFDW.grid$x.int) * rbind(C,C.x.down,deparse.level = 0)))

  }
  y.Adv.flux <- 0
  if (any(v.grid$y.int >0) ) {
    vv <- v.grid$y.int
    vv[vv<0]<-0
    y.Adv.flux <-  y.Adv.flux + as.matrix(VF.grid$y.int * vv * (
                 (1-AFDW.grid$y.int) * cbind(C.y.up,C,deparse.level = 0)
                 + AFDW.grid$y.int * cbind(C,C.y.down,deparse.level = 0)))
  }
  if (any (v.grid$y.int < 0)) {
    vv <- v.grid$y.int
    vv[vv>0]<-0
    y.Adv.flux <-  y.Adv.flux + as.matrix(VF.grid$y.int * vv * (
                    AFDW.grid$y.int * cbind(C.y.up,C,deparse.level = 0)
                 + (1-AFDW.grid$y.int) * cbind(C,C.y.down,deparse.level = 0)))
  }

  x.flux <- x.Dif.flux + x.Adv.flux
  y.flux <- y.Dif.flux + y.Adv.flux

## Impose boundary fluxes when needed
## Default boundary condition is no gradient
  if (! is.null (flux.x.up[1]))
    x.flux[1,]   <- flux.x.up
  if (! is.null (flux.x.down[1]))
    x.flux[nrow(x.flux),] <- flux.x.down
    
  if (! is.null (flux.y.up[1]))
    y.flux[,1]   <- flux.y.up
  if (! is.null (flux.y.down[1]))
    y.flux[,ncol(y.flux)] <- flux.y.down

## Calculate rate of change = flux gradient
  dFdx <- - (diff(A.grid$x.int*x.flux) / A.grid$x.mid/grid$dx ) / VF.grid$x.mid
  dFdy <- -t(diff(t(A.grid$y.int*y.flux))/t(A.grid$y.mid)/grid$dy) / VF.grid$y.mid


  if (!full.output) {
    return (list (dC = dFdx + dFdy,                    # Rate of change due to advective-diffuisve transport in each grid cell
                  flux.x.up = x.flux[1,],                # flux across lower boundary interface; positive = IN
                  flux.x.down = x.flux[nrow(x.flux),],   # flux across lower boundary interface; positive = OUT
                  flux.y.up = y.flux[,1],                # flux across lower boundary interface; positive = IN
                  flux.y.down = y.flux[,ncol(y.flux)]))  # flux across lower boundary interface; positive = OUT

  } else {
    return (list (dC = dFdx + dFdy,                    # Rate of change in the centre of each grid cells
                  C.x.up = C.x.up,                     # concentration at upper interface
                  C.x.down = C.x.down,                 # concentration at upper interface
                  C.y.up = C.y.up,                     # concentration at upper interface
                  C.y.down = C.y.down,                 # concentration at upper interface
                  x.flux = x.flux,                     # flux across at the interface of each grid cell
                  y.flux = y.flux,                     # flux across at the interface of each grid cell
                  flux.x.up = x.flux[1,],               # flux across lower boundary interface; positive = IN
                  flux.x.down = x.flux[nrow(x.flux),],  # flux across lower boundary interface; positive = OUT
                  flux.y.up = y.flux[,1],               # flux across lower boundary interface; positive = IN
                  flux.y.down = y.flux[,ncol(y.flux)])) # flux across lower boundary interface; positive = OUT
  }
} # end tran.2D

