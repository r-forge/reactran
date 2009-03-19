
##==============================================================================
## 2D Transport of a solute in the sediment
##==============================================================================

tran.2D <- function(C, C.up=C[1,], C.down=C[nrow(C),],
  C.left=C[,1],  C.right=C[,ncol(C)],
  flux.up=NULL, flux.down=NULL, flux.left=NULL, flux.right=NULL,
  a.bl.up=NULL, C.bl.up=NULL, a.bl.down=NULL, C.bl.down=NULL,
  a.bl.left=NULL, C.bl.left=NULL, a.bl.right=NULL, C.bl.right=NULL,
  D.x=NULL, D.y=D.x, D.grid=NULL,
	v.x=0, v.y=v.x, v.grid=NULL,
	AFDW.x=1, AFDW.y=AFDW.x, AFDW.grid=NULL,
  VF.x=1, VF.y=VF.x, VF.grid=NULL,
	dx=NULL, dy=NULL, grid=NULL,
  full.check = FALSE, full.output = FALSE)
											
{

# DEFAULT INFILLING OF GRID PARAMETERS

  if (is.null(grid))
    grid <- list(dx=rep(dx,length.out=nrow(C)),
                        dx.aux=0.5*(c(0,rep(dx,length.out=nrow(C)))+
                                    c(rep(dx,length.out=nrow(C)),0)),
                        dy=rep(dy,length.out=ncol(C)),
                        dy.aux=0.5*(c(0,rep(dy,length.out=ncol(C)))+
                           c(rep(dy,length.out=ncol(C)),0)))
  if (is.null(AFDW.grid))
    AFDW.grid <- list(x.int=matrix(data=AFDW.x,nrow=(nrow(C)+1),ncol=ncol(C)),
                      y.int=matrix(data=AFDW.y,nrow=nrow(C),ncol=(ncol(C)+1)))
  if (is.null(D.grid))
    D.grid <- list(x.int=matrix(data=D.x,nrow=(nrow(C)+1),ncol=ncol(C)),
                   y.int=matrix(data=D.y,nrow=nrow(C),ncol=(ncol(C)+1)))
  if (is.null(v.grid))
    v.grid <- list(x.int=matrix(data=v.x,nrow=(nrow(C)+1),ncol=ncol(C)),
                   y.int=matrix(data=v.y,nrow=nrow(C),ncol=(ncol(C))+1))
  if (is.null(VF.grid))
    VF.grid <- list(x.int=matrix(data=VF.x,nrow=(nrow(C)+1),ncol=ncol(C)),
                    y.int=matrix(data=VF.y,nrow=nrow(C),ncol=(ncol(C)+1)),
                    x.mid=matrix(data=VF.x,nrow=nrow(C),ncol=ncol(C)),
                    y.mid=matrix(data=VF.y,nrow=nrow(C),ncol=ncol(C)))

#==============================================================================
# INPUT CHECKS  
#==============================================================================


  if (full.check) {

## check dimensions of input concentrations

    if (!is.null(C.up)) {
      if (!((length(C.up)==1) || (length(C.up)==(ncol(C)))))
        stop("error: C.up should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(C.down)) {
      if (!((length(C.down)==1) || (length(C.down)==(ncol(C)))))
        stop("error: C.down should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(C.left)) {
      if (!((length(C.left)==1) || (length(C.left)==(nrow(C)))))
        stop("error: C.left should be a vector of length 1 or nrow(C)")
    }
    if (!is.null(C.right)) {
      if (!((length(C.right)==1) || (length(C.right)==(nrow(C)))))
        stop("error: C.right should be a vector of length 1 or nrow(C)")
    }

    if (!is.null(C.bl.up)) {
      if (!((length(C.bl.up)==1) || (length(C.bl.up)==(ncol(C)))))
        stop("error: C.bl.up should be a vector of length 1 or ncol(C)")
   }

    if (!is.null(C.bl.down)) {
      if (!((length(C.bl.down)==1) || (length(C.bl.down)==(ncol(C)))))
        stop("error: C.bl.down should be a vector of length 1 or ncol(C)")
    }

    if (!is.null(C.bl.left)) {
      if (!((length(C.bl.left)==1) || (length(C.bl.left)==(nrow(C)))))
        stop("error: C.bl.left should be a vector of length 1 or nrow(C)")
    }

    if (!is.null(C.bl.right)) {
      if (!((length(C.bl.right)==1) || (length(C.bl.right)==(nrow(C)))))
        stop("error: C.bl.right should be a vector of length 1 or nrow(C)")
    }

# check dimensions of input fluxes

    if (!is.null(flux.up)) {
      if (!((length(flux.up)==1) || (length(flux.up)==(ncol(C)))))
        stop("error: flux.up should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(flux.down)) {
      if (!((length(flux.down)==1) || (length(flux.down)==(ncol(C)))))
        stop("error: flux.down should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(flux.left)) {
      if (!((length(flux.left)==1) || (length(flux.left)==(nrow(C)))))
        stop("error: flux.left should be a vector of length 1 or nrow(C)")
    }

    if (!is.null(flux.right)) {
      if (!((length(flux.right)==1) || (length(flux.right)==(nrow(C)))))
        stop("error: flux.right should be a vector of length 1 or nrow(C)")
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

    if (is.null(VF.x) && is.null(VF.y) && is.null(VF.grid))
      stop("error: VF.x, VF.y, and VF.grid cannot be NULL at the same time")

    gn <- names(VF.grid)
    if (! "x.int" %in% gn)
      stop("error: VF.grid should be a list that contains 'x.int', the values at the vertical interfaces of the grid cells")
    if (! "y.int" %in% gn)
      stop("error: VF.grid should be a list that contains 'y.int', the values at the horizontal interfaces of the grid cells")
    if (! "x.mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'x.mid', the values at the middle of the grid cells")
    if (! "y.mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'y.mid', the values at the middle of the grid cells")
    if (is.null(VF.grid$x.int) || is.null(VF.grid$y.int) || is.null(VF.grid$x.mid) || is.null(VF.grid$y.mid))
      stop("error: the VF.grid should be a list with (numeric) values")
    if (any (VF.grid$x.int < 0) || any (VF.grid$y.int < 0) || any (VF.grid$x.mid < 0) || any (VF.grid$y.mid < 0))
      stop("error: the VF.grid values should always be positive")

  }
## FUNCTION BODY: CALCULATIONS

  N <- nrow(C)
  M <- ncol(C)

## Impose boundary flux at upper boundary when needed
## Default boundary condition is no gradient
  if (! is.null (flux.up[1])) {
    nom <- flux.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1] +
           (1-AFDW.grid$x.int[1,])*v.grid$x.int[1,])*C[1,]
    denom <- VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1]+
             AFDW.grid$x.int[1,]*v.grid$x.int[1,])
    C.up <- nom/denom
  }

## Impose boundary flux at lower boundary when needed
## Default boundary condition is no gradient
  if (! is.null (flux.down[1])) {
  	nom <- flux.down - VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/
            grid$dx.aux[N+1] + AFDW.grid$x.int[(N+1),]*v.grid$x.int[(N+1),])*C[N,]
    denom <- -VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/grid$dx.aux[N+1]+
            (1-AFDW.grid$x.int[(N+1),])*v.grid$x.int[(N+1),])
    C.down <- nom/denom
  }

# Impose boundary flux at upper boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.left[1])) {
    nom <- flux.left + VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1] +
           (1-AFDW.grid$y.int[,1])*v.grid$y.int[,1])*C[,1]
    denom <- VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1]+
             AFDW.grid$y.int[,1]*v.grid$y.int[,1])
    C.left <- nom/denom
  }

# Impose boundary flux at lower boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.right[1]))  {
	  nom <- flux.right - VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/
           grid$dy.aux[M+1] + AFDW.grid$y.int[,(M+1)]*v.grid$y.int[,(M+1)])*C[,M]
    denom <- -VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/grid$dy.aux[M+1]+
             (1-AFDW.grid$y.int[,(M+1)])*v.grid$y.int[,(M+1)])
    C.right <- nom/denom
  }

## when upper boundary layer is present, calculate new C.up
  if (!is.null(a.bl.up) & !is.null(C.bl.up[1])) {
	  nom <- a.bl.up*C.bl.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/
           grid$dx.aux[1] + (1-AFDW.grid$x.int[1,])*v.grid$x.int[1,])*C[1,]
    denom <- a.bl.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1]+
             AFDW.grid$x.int[1,]*v.grid$x.int[1,])
	  C.up <- nom/denom
  }

## when lower boundary layer is present, calculate new C.down
  if (!is.null(a.bl.down) & !is.null(C.bl.down[1])) {
	  nom <- a.bl.down*C.bl.down + VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/
           grid$dx.aux[(N+1)] + (1-AFDW.grid$x.int[(N+1),])*
           v.grid$x.int[(N+1),])*C[N,]
    denom <- a.bl.down + VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/
             grid$dx.aux[(N+1)]+ AFDW.grid$x.int[(N+1),]*v.grid$x.int[(N+1),])
	  C.down <- nom/denom
  }

## when left boundary layer is present, calculate new C.left
  if (!is.null(a.bl.left) & !is.null(C.bl.left[1])) {
	  nom <- a.bl.left*C.bl.left + VF.grid$y.int[,1]*(D.grid$y.int[,1]/
           grid$dy.aux[1] + (1-AFDW.grid$y.int[,1])*v.grid$y.int[,1])*C[,1]
    denom <- a.bl.left + VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1]+
             AFDW.grid$y.int[,1]*v.grid$y.int[,1])
	  C.left <- nom/denom
  }

## when right boundary layer is present, calculate new C.right
  if (!is.null(a.bl.right) & !is.null(C.bl.right[1]))   {
	  nom <- a.bl.right*C.bl.right + VF.grid$y.int[,(M+1)]*
           (D.grid$y.int[,(M+1)]/grid$dy.aux[(M+1)] +
           (1-AFDW.grid$y.int[,(M+1)])*v.grid$y.int[,(M+1)])*C[,M]
    denom <- a.bl.right + VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/
             grid$dy.aux[(M+1)]+ AFDW.grid$y.int[,(M+1)]*v.grid$y.int[,(M+1)])
	  C.right <- nom/denom
  }

## Calculate diffusive part of the flux
  x.Dif.flux <- as.matrix(-VF.grid$x.int*D.grid$x.int*
                diff(rbind(C.up,C,C.down,deparse.level = 0))/
                matrix(data=grid$dx.aux,nrow=(N+1),ncol=M,byrow=FALSE))
  y.Dif.flux <- as.matrix(-VF.grid$y.int*D.grid$y.int*
                t(diff(t(cbind(C.left,C,C.right,deparse.level = 0))))/
                matrix(data=grid$dy.aux,nrow=N,ncol=(M+1),byrow=TRUE))

## Calculate adcevtive part of the flux
  x.Adv.flux <- as.matrix(VF.grid$x.int*v.grid$x.int*((1-AFDW.grid$x.int)*
                rbind(C.up,C,deparse.level = 0)+AFDW.grid$x.int*
                rbind(C,C.down,deparse.level = 0)))
  y.Adv.flux <- as.matrix(VF.grid$y.int*v.grid$y.int*((1-AFDW.grid$y.int)*
                cbind(C.left,C,deparse.level = 0)+AFDW.grid$y.int*
                cbind(C,C.right,deparse.level = 0)))

  x.flux <- x.Dif.flux + x.Adv.flux
  y.flux <- y.Dif.flux + y.Adv.flux

## Impose boundary fluxes when needed
## Default boundary condition is no gradient
  if (! is.null (flux.up[1]))
    x.flux[1,]   <- flux.up
  if (! is.null (flux.down[1]))
    x.flux[nrow(x.flux),] <- flux.down
    
  if (! is.null (flux.left[1]))
    y.flux[,1]   <- flux.left
  if (! is.null (flux.right[1]))
    y.flux[,ncol(y.flux)] <- flux.right

## Calculate rate of change = flux gradient
  dFdx <- -(diff(x.flux)/grid$dx)/VF.grid$x.mid
  dFdy <- -t(diff(t(y.flux))/grid$dy)/VF.grid$y.mid

  if (!full.output) {
    return (list (dC = dFdx + dFdy))                    # Rate of change due to advective-diffuisve transport in each grid cell
  } else {
    return (list (dC = dFdx + dFdy,                    # Rate of change in the centre of each grid cells
                  C.up = C.up,                         # concentration at upper interface
                  C.down = C.down,                     # concentration at upper interface
                  C.left = C.left,                     # concentration at upper interface
                  C.right = C.right,                   # concentration at upper interface
                  x.flux = x.flux,                     # flux across at the interface of each grid cell
                  y.flux = y.flux,                     # flux across at the interface of each grid cell
                  flux.up = x.flux[1,],                # flux across lower boundary interface; positive = IN
                  flux.down = x.flux[nrow(x.flux),],   # flux across lower boundary interface; positive = OUT
                  flux.left = y.flux[,1],              # flux across lower boundary interface; positive = IN
                  flux.right = y.flux[,ncol(y.flux)])) # flux across lower boundary interface; positive = OUT
  }
} # end tran.2D

