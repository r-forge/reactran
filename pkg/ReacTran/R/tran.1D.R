
##==============================================================================
## Transport in a one-dimensional finite difference grid
##==============================================================================

tran.1D <- function(C, C.up=C[1], C.down=C[length(C)],
     flux.up=NULL, flux.down=NULL, a.bl.up=NULL, C.bl.up=NULL,
		 a.bl.down=NULL, C.bl.down=NULL,
     D=NULL, D.grid=NULL, v=0, v.grid=NULL,
  	 AFDW=1, AFDW.grid=NULL, VF=1, VF.grid=NULL, A=1, A.grid=NULL,
     dx=NULL, grid=NULL, full.check = FALSE, full.output = FALSE) {


### INPUT CHECKS
  N <- length(C)
  if (N == 0)
    stop("C should be a vector with numeric values")
  if (is.null(grid) && is.null(dx))
        stop("error: dx or grid should be specified")

  if (is.null(C.up))   C.up <- C[1]
  if (is.null(C.down)) C.down <- C[N]

  ## default infilling of grid parameters

  if (is.null(AFDW.grid))
    AFDW.grid <- list(int=AFDW)
  if (is.null(D.grid))
    D.grid <- list(int=D)
  if (is.null(D.grid$int)) D.grid$int <- 0
  if (is.null(v.grid))
    v.grid <- list(int=v)
  if (is.null(VF.grid))
    VF.grid <- list(int=rep(VF,length.out=N+1),
                    mid=0.5*(rep(VF,length.out=(N+1))[1:N]+
                    rep(VF,length.out=(N+1))[2:(N+1)]))
  if (is.null(A.grid))
    A.grid <- list(int=rep(A,length.out=(N+1)),
                   mid=0.5*(rep(A,length.out=(length(C)+1))[1:N]+
                   rep(A,length.out=(N+1))[2:(N+1)]))

  if (is.null(grid))
    grid <- list(dx=rep(dx,length.out=N),
                 dx.aux=0.5*(c(0,rep(dx,length.out=N))+
                 c(rep(dx,length.out=N),0)))

  ## check dimensions of input arguments

  if (full.check) {

    if (!is.null(AFDW)) {
      if (!((length(AFDW)==1) || (length(AFDW)==(N+1))))
        stop("error: AFDW should be a vector of length 1 or N+1")
    }

    if (!is.null(D)) {
      if (!((length(D)==1) || (length(D)==(N+1))))
        stop("error: D should be a vector of length 1 or N+1")
    }

    if (!is.null(v)) {
      if (!((length(v)==1) || (length(v)==(N+1))))
        stop("error: v should be a vector of length 1 or N+1")
    }

    if (!is.null(VF)) {
      if (!((length(VF)==1) || (length(VF)==(N+1))))
        stop("error: VF should be a vector of length 1 or N+1")
    }

    if (!is.null(A)) {
      if (!((length(A)==1) || (length(A)==(N+1))))
        stop("error: A should be a vector of length 1 or N+1")
    }

    if (!is.null(dx)) {
      if (!((length(dx)==1) || (length(dx)==N)))
        stop("error: dx should be a vector of length 1 or N")
    }

    ## check input of grid
    gn <- names(grid)
    if (! "dx" %in% gn)
      stop("error: grid should be a list that contains 'dx' ")
    if (! "dx.aux" %in% gn)
	    stop("error: grid should be a list that contains 'dx.aux' ")
    if (is.null(grid$dx) || length(grid$dx) == 0 ||
        is.null(grid$dx.aux) || length(grid$dx.aux) == 0)
    	stop("error: the grid should be a list with (numeric) values for 'dx' and 'dx.aux' ")
    if (any(grid$dx <= 0) || any(grid$dx.aux <= 0) )
    	stop("error: the grid distances dx and dx.aux should always be positive")

    ## check input of AFDW.grid
    gn <- names(AFDW.grid)
    if (! "int" %in% gn)
      stop("error: AFDW.grid should be a list that contains 'int', the AFDW values at the interface of the grid cell ")
    if (is.null(AFDW.grid$int) || length(AFDW.grid$int) == 0)
      stop("error: AFDW.grid should be a list with (numeric) values")
    if (any(AFDW.grid$int < 0)||any(AFDW.grid$int > 1))
	    stop("error: the AFDW should always range between 0 and 1")

    ## check input of D.grid
    gn <- names(D.grid)
    if (! "int" %in% gn)
      stop("error: D.grid should be a list that contains 'int', the D values at the interface of the grid cell ")
    if (is.null(D.grid$int) || length(D.grid$int) == 0)
      stop("error: D.grid should be a list with (numeric) values")
    if (any(D.grid$int < 0))
	    stop("error: the mixing coefficient should always be positive")

    ## check input of v.grid
    gn <- names(v.grid)
    if (! "int" %in% gn)
      stop("error: v.grid should be a list that contains 'int', the v values at the interface of the grid cell ")
    if (is.null(v.grid$int)|| length(v.grid$int) == 0)
      stop("error: the advective velocity v should be a list with (numeric) values")
    if (any (v.grid$int < 0) & any (v.grid$int > 0))
	    stop("error: the advective velocity cannot be both positive and negative within the same domain")

    ## check input of VF.grid
    gn <- names(VF.grid)
    if (! "int" %in% gn)
      stop("error: VF.grid should be a list that contains 'int', the area values at the interface of the grid cell ")
    if (! "mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'mid', the area at the middle of the grid cells")
    if (is.null(VF.grid$int)|| length(VF.grid$int) == 0 ||
        is.null(VF.grid$mid)|| length(VF.grid$mid) == 0)
      stop("error: the volume fraction VF should be a list with (numeric) values")
    if (any (VF.grid$int < 0) || any (VF.grid$mid < 0))
      stop("error: the volume fraction should always be positive")

    ## check input of A.grid
    gn <- names(A.grid)
    if (! "int" %in% gn)
      stop("error: A.grid should be a list that contains 'int', the area values at the interface of the grid cell ")
    if (! "mid" %in% gn)
      stop("error: A.grid should be a list that contains 'mid', the area at the middle of the grid cells")
    if (is.null(A.grid$int) || length(A.grid$int) == 0 ||
        is.null(A.grid$mid) || length(A.grid$mid)==0)
      stop("error: the volume fraction VF should be a list with (numeric) values")
    if (any (A.grid$int < 0) || any (A.grid$mid < 0))
      stop("error: the area A should always be positive")

    ## check input of boundary layer parameters
    if (!is.null(a.bl.up) & !is.null(C.bl.up)) {
    	if (a.bl.up < 0)
        stop("error: the boundary layer transfer coefficient should be positive")
    }

    if (!is.null(a.bl.down) & !is.null(C.bl.down)) {
	    if (a.bl.down < 0)
        stop("error: the boundary layer transfer coefficient should be positive")
    }
  } # end full.check


### FUNCTION BODY: CALCULATIONS


  if (full.output) {
    ## Impose boundary flux at upper boundary when needed
    ## Default boundary condition is zero gradient
    if (! is.null (flux.up)) {
      if (any (v.grid$int >= 0)) {
      ## advection directed downwards
        nom <- flux.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
               (1-AFDW.grid$int[1])*v.grid$int[1])*C[1]
        denom <- VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1]+
                 AFDW.grid$int[1]*v.grid$int[1])
	    } else	{
      ## advection directed upwards
        nom <- flux.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
               AFDW.grid$int[1]*v.grid$int[1])*C[1]
        denom <- VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1]+
                (1-AFDW.grid$int[1])*v.grid$int[1])
	    }
      C.up <- nom/denom
    }
    ## Impose boundary flux at lower boundary when needed
    ## Default boundary condition is no gradient
    if (! is.null (flux.down)) {
      if (any (v.grid$int >= 0)) {
        ## advection directed downwards
	      nom <- flux.down - VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] +
               AFDW.grid$int[N+1]*v.grid$int[N+1])*C[N]
        denom <- -VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1]+
                (1-AFDW.grid$int[N+1])*v.grid$int[N+1])
      } else {
        ## advection directed downwards
	      nom <- flux.down - VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] +
               (1-AFDW.grid$int[N+1])*v.grid$int[N+1])*C[N]
        denom <- -VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1]+
                 AFDW.grid$int[N+1]*v.grid$int[N+1])
      }
      C.down <- nom/denom
    }

    ## when upstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.up) & !is.null(C.bl.up)) {
      if (any (v.grid$int >= 0))	{
        ## advection directed downwards
        nom <- a.bl.up*C.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
              (1-AFDW.grid$int[1])*v.grid$int[1])*C[1]
        denom <- a.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
                 AFDW.grid$int[1]*v.grid$int[1])
	    } else	{
        ## advection directed upwards
        nom <- a.bl.up*C.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
               AFDW.grid$int[1]*v.grid$int[1])*C[1]
        denom <- a.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
                (1-AFDW.grid$int[1])*v.grid$int[1])
	    }
      C.up <- nom/denom
    }

    ## when downstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.down) & !is.null(C.bl.down)) {
      if (any (v.grid$int >= 0))	{
        ## advection directed downwards
       nom <- a.bl.down*C.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/
              grid$dx.aux[N+1] + (1-AFDW.grid$int[N+1])*v.grid$int[N+1])*C[N]
       denom <- a.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] +
                AFDW.grid$int[N+1]*v.grid$int[N+1])
     	} else	{
        ## advection directed upwards
        nom <- a.bl.down*C.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/
               grid$dx.aux[N+1] + AFDW.grid$int[N+1]*v.grid$int[N+1])*C[N]
        denom <- a.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] +
                 (1-AFDW.grid$int[N+1])*v.grid$int[N+1])
	    }
      C.down <- nom/denom
    }

    ## Calculate diffusive part of the flux

    dif.flux <- as.vector(-VF.grid$int*D.grid$int*diff(c(C.up,C,C.down))/
                          grid$dx.aux)
    adv.flux <- rep(0,length.out=length(dif.flux))

    ## Add advective part of the flux if needed
    ## KS: CHANGED
#    if (any (v.grid$int > 0)) {   # advection directed downwards
#	     conc <- AFDW.grid$int*c(C.up,C)
#	     if (any (AFDW.grid$int < 1))
#         conc <- conc +(1-AFDW.grid$int)*c(C,C.down)
#	     adv.flux <- as.vector(VF.grid$int*v.grid$int*conc)
#    }
#    if (any (v.grid$int < 0)) {   # advection directed upwards
#	    conc <- AFDW.grid$int*c(C,C.down)
#	    if (any (AFDW.grid$int < 1))
#        conc <- conc +(1-AFDW.grid$int)*c(C.up,C)
#	    adv.flux <- as.vector(VF.grid$int*v.grid$int*conc)
#    }
      if (any (v.grid$int > 0)) {   # advection directed downwards
       v <- v.grid$int
       v[v.grid$int<0] <- 0
	     conc <- AFDW.grid$int*c(C.up,C)
	     if (any (AFDW.grid$int < 1))
         conc <- conc +(1-AFDW.grid$int)*c(C,C.down)
	     adv.flux <- adv.flux + as.vector(VF.grid$int*v*conc)
    }
    if (any (v.grid$int < 0)) {   # advection directed upwards
       v <- v.grid$int
       v[v.grid$int>0] <- 0

      conc <- AFDW.grid$int*c(C,C.down)
	    if (any (AFDW.grid$int < 1))
        conc <- conc +(1-AFDW.grid$int)*c(C.up,C)
	     adv.flux <- adv.flux + as.vector(VF.grid$int*v*conc)
    }

    flux <- dif.flux + adv.flux

  } else { # not full.output

    ## when upstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.up) & !is.null(C.bl.up))  {
      if (any (v.grid$int >= 0)) { # advection directed downwards
        nom <- a.bl.up*C.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
               (1-AFDW.grid$int[1])*v.grid$int[1])*C[1]
        denom <- a.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
               AFDW.grid$int[1]*v.grid$int[1])
	    } else	{  # advection directed upwards
        nom <- a.bl.up*C.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
               AFDW.grid$int[1]*v.grid$int[1])*C[1]
        denom <- a.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] +
                 (1-AFDW.grid$int[1])*v.grid$int[1])
	    }
      C.up <- nom/denom
    }

    ## when upstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.down) & !is.null(C.bl.down)) {
      if (any (v.grid$int >= 0))	{   # advection directed downwards
        nom <- a.bl.down*C.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/
               grid$dx.aux[N+1] + (1-AFDW.grid$int[N+1])*v.grid$int[N+1])*C[N]
        denom <- a.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] +
                 AFDW.grid$int[N+1]*v.grid$int[N+1])
	    } else	{  # advection directed upwards
        nom <- a.bl.down*C.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/
                 grid$dx.aux[N+1] + AFDW.grid$int[N+1]*v.grid$int[N+1])*C[N]
        denom <- a.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] +
                (1-AFDW.grid$int[N+1])*v.grid$int[N+1])
	    }
      C.down <- nom/denom
    }

    ## Calculate diffusive part of the flux
	  flux <- as.vector(-(VF.grid$int)*D.grid$int*
                      diff(c(C.up,C,C.down))/grid$dx.aux)
    ## Add advective part of the flux if needed
    # KS CHANGED....
#	  if (any (v.grid$int > 0)) 	{   # advection directed downwards
#  		conc <- AFDW.grid$int*c(C.up,C)
#	  	if (any (AFDW.grid$int < 1))
#        conc <- conc +(1-AFDW.grid$int)*c(C,C.down)
#		  flux <- flux + as.vector((VF.grid$int)*v.grid$int*conc)
#	  }

#	  if (any (v.grid$int < 0))  {   # advection directed upwards
#  		conc <- AFDW.grid$int*c(C,C.down)
#		  if (any (AFDW.grid$int < 1))
#        conc <- conc +(1-AFDW.grid$int)*c(C.up,C)
#		  flux <- flux + as.vector((VF.grid$int)*v.grid$int*conc)
#	  }
      if (any (v.grid$int > 0)) {   # advection directed downwards
       v <- v.grid$int
       v[v.grid$int<0] <- 0
	     conc <- AFDW.grid$int*c(C.up,C)
	     if (any (AFDW.grid$int < 1))
         conc <- conc +(1-AFDW.grid$int)*c(C,C.down)
	     flux <- flux + as.vector(VF.grid$int*v*conc)
    }
    if (any (v.grid$int < 0)) {   # advection directed upwards
       v <- v.grid$int
       v[v.grid$int>0] <- 0

      conc <- AFDW.grid$int*c(C,C.down)
	    if (any (AFDW.grid$int < 1))
        conc <- conc +(1-AFDW.grid$int)*c(C.up,C)
	     flux <- flux + as.vector(VF.grid$int*v*conc)
    }

  }

  if (! is.null (flux.up))
    flux[1] <- flux.up
  if (! is.null (flux.down))
    flux[N+1] <- flux.down

    
## Calculate rate of change = Flux gradient       ## Karline: /VF.grid$mid
  dC <- -diff(A.grid$int*flux)/A.grid$mid/grid$dx   / VF.grid$mid

  if (!full.output){
    return (list (dC = dC,                   # Rate of change due to advective-diffuisve transport in each grid cell
  					flux.up = flux[1],               # Flux across lower boundary interface; positive = IN
	  				flux.down = flux[length(flux)])) # Flux across lower boundary interface; positive = OUT

  } else {
  	return (list (dC = dC,                   # Rate of change due to advective-diffuisve transport in each grid cell
	  				C.up = C.up,                     # Concentration at upper interface
		  			C.down = C.down,                 # Concentration at lower interface
			  		dif.flux = dif.flux,             # Diffusive flux across at the interface of each grid cell
				  	adv.flux = adv.flux,             # Advective flux across at the interface of each grid cell
					  flux = flux,                     # Flux across at the interface of each grid cell
  					flux.up = flux[1],               # Flux across lower boundary interface; positive = IN
	  				flux.down = flux[length(flux)])) # Flux across lower boundary interface; positive = OUT
  }

} # end tran.1D

