
##==============================================================================
## Transport in a one-dimensional finite difference grid
##==============================================================================

tran.1D <- function(C, C.up=C[1], C.down=C[length(C)],
     flux.up=NULL, flux.down=NULL, a.bl.up=NULL, C.bl.up=NULL,
		 a.bl.down=NULL, C.bl.down=NULL,
     D=0, v=0, AFDW=1, VF=1, A=1,
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

  if (!is.list(AFDW))
    AFDW <- list(int=AFDW)
  if (!is.list(D))
    D <- list(int=D)
  if (is.null(D$int)) D$int <- 0
  if (!is.list(v))
    v <- list(int=v)
  if (!is.list(VF))
    VF <- list(int=rep(VF,length.out=N+1),
                mid=0.5*(rep(VF,length.out=(N+1))[1:N]+
                rep(VF,length.out=(N+1))[2:(N+1)]))
  if (!is.list(A))
    A <- list(int=rep(A,length.out=(N+1)),
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

    ## check input of AFDW
    gn <- names(AFDW)
    if (! "int" %in% gn)
      stop("error: AFDW should be a list that contains 'int', the AFDW values at the interface of the grid cell ")
    if (is.null(AFDW$int) || length(AFDW$int) == 0)
      stop("error: AFDW should contain (numeric) values")
    if (any(AFDW$int < 0)||any(AFDW$int > 1))
	    stop("error: the AFDW should always range between 0 and 1")

    ## check input of D
    gn <- names(D)
    if (! "int" %in% gn)
      stop("error: D should be a list that contains 'int', the D values at the interface of the grid cell ")
    if (is.null(D$int) || length(D$int) == 0)
      stop("error: D should contain (numeric) values")
    if (any(D$int < 0))
	    stop("error: the mixing coefficient should always be positive")

    ## check input of v
    gn <- names(v)
    if (! "int" %in% gn)
      stop("error: v should be a list that contains 'int', the v values at the interface of the grid cell ")
    if (is.null(v$int)|| length(v$int) == 0)
      stop("error: the advective velocity v should contain (numeric) values")
    if (any (v$int < 0) & any (v$int > 0))
	    stop("error: the advective velocity cannot be both positive and negative within the same domain")

    ## check input of VF
    gn <- names(VF)
    if (! "int" %in% gn)
      stop("error: VF should be a list that contains 'int', the area values at the interface of the grid cell ")
    if (! "mid" %in% gn)
      stop("error: VF should be a list that contains 'mid', the area at the middle of the grid cells")
    if (is.null(VF$int)|| length(VF$int) == 0 ||
        is.null(VF$mid)|| length(VF$mid) == 0)
      stop("error: the volume fraction VF should contain (numeric) values")
    if (any (VF$int < 0) || any (VF$mid < 0))
      stop("error: the volume fraction should always be positive")

    ## check input of A
    gn <- names(A)
    if (! "int" %in% gn)
      stop("error: A should be a list that contains 'int', the area values at the interface of the grid cell ")
    if (! "mid" %in% gn)
      stop("error: A should be a list that contains 'mid', the area at the middle of the grid cells")
    if (is.null(A$int) || length(A$int) == 0 ||
        is.null(A$mid) || length(A$mid)==0)
      stop("error: the surface area A should contain (numeric) values")
    if (any (A$int < 0) || any (A$mid < 0))
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
      if (any (v$int >= 0)) {
      ## advection directed downwards
        nom <- flux.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               (1-AFDW$int[1])*v$int[1])*C[1]
        denom <- VF$int[1]*(D$int[1]/grid$dx.aux[1]+
                 AFDW$int[1]*v$int[1])
	    } else	{
      ## advection directed upwards
        nom <- flux.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               AFDW$int[1]*v$int[1])*C[1]
        denom <- VF$int[1]*(D$int[1]/grid$dx.aux[1]+
                (1-AFDW$int[1])*v$int[1])
	    }
      C.up <- nom/denom
    }
    ## Impose boundary flux at lower boundary when needed
    ## Default boundary condition is no gradient
    if (! is.null (flux.down)) {
      if (any (v$int >= 0)) {
        ## advection directed downwards
	      nom <- flux.down - VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
               AFDW$int[N+1]*v$int[N+1])*C[N]
        denom <- -VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1]+
                (1-AFDW$int[N+1])*v$int[N+1])
      } else {
        ## advection directed downwards
	      nom <- flux.down - VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
               (1-AFDW$int[N+1])*v$int[N+1])*C[N]
        denom <- -VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1]+
                 AFDW$int[N+1]*v$int[N+1])
      }
      C.down <- nom/denom
    }

    ## when upstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.up) & !is.null(C.bl.up)) {
      if (any (v$int >= 0))	{
        ## advection directed downwards
        nom <- a.bl.up*C.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
              (1-AFDW$int[1])*v$int[1])*C[1]
        denom <- a.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
                 AFDW$int[1]*v$int[1])
	    } else	{
        ## advection directed upwards
        nom <- a.bl.up*C.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               AFDW$int[1]*v$int[1])*C[1]
        denom <- a.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
                (1-AFDW$int[1])*v$int[1])
	    }
      C.up <- nom/denom
    }

    ## when downstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.down) & !is.null(C.bl.down)) {
      if (any (v$int >= 0))	{
        ## advection directed downwards
       nom <- a.bl.down*C.bl.down + VF$int[N+1]*(D$int[N+1]/
              grid$dx.aux[N+1] + (1-AFDW$int[N+1])*v$int[N+1])*C[N]
       denom <- a.bl.down + VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
                AFDW$int[N+1]*v$int[N+1])
     	} else	{
        ## advection directed upwards
        nom <- a.bl.down*C.bl.down + VF$int[N+1]*(D$int[N+1]/
               grid$dx.aux[N+1] + AFDW$int[N+1]*v$int[N+1])*C[N]
        denom <- a.bl.down + VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
                 (1-AFDW$int[N+1])*v$int[N+1])
	    }
      C.down <- nom/denom
    }

    ## Calculate diffusive part of the flux

    dif.flux <- as.vector(-VF$int*D$int*diff(c(C.up,C,C.down))/
                          grid$dx.aux)
    adv.flux <- rep(0,length.out=length(dif.flux))

    ## Add advective part of the flux if needed
      if (any (v$int > 0)) {   # advection directed downwards
       vv <- v$int
       vv[v$int<0] <- 0
	     conc <- AFDW$int*c(C.up,C)
	     if (any (AFDW$int < 1))
         conc <- conc +(1-AFDW$int)*c(C,C.down)
	     adv.flux <- adv.flux + as.vector(VF$int*vv*conc)
    }
    if (any (v$int < 0)) {   # advection directed upwards
       vv <- v$int
       vv[v$int>0] <- 0

      conc <- AFDW$int*c(C,C.down)
	    if (any (AFDW$int < 1))
        conc <- conc +(1-AFDW$int)*c(C.up,C)
	     adv.flux <- adv.flux + as.vector(VF$int*vv*conc)
    }

    flux <- dif.flux + adv.flux

  } else { # not full.output

    ## when upstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.up) & !is.null(C.bl.up))  {
      if (any (v$int >= 0)) { # advection directed downwards
        nom <- a.bl.up*C.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               (1-AFDW$int[1])*v$int[1])*C[1]
        denom <- a.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               AFDW$int[1]*v$int[1])
	    } else	{  # advection directed upwards
        nom <- a.bl.up*C.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               AFDW$int[1]*v$int[1])*C[1]
        denom <- a.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
                 (1-AFDW$int[1])*v$int[1])
	    }
      C.up <- nom/denom
    }

    ## when upstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.down) & !is.null(C.bl.down)) {
      if (any (v$int >= 0))	{   # advection directed downwards
        nom <- a.bl.down*C.bl.down + VF$int[N+1]*(D$int[N+1]/
               grid$dx.aux[N+1] + (1-AFDW$int[N+1])*v$int[N+1])*C[N]
        denom <- a.bl.down + VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
                 AFDW$int[N+1]*v$int[N+1])
	    } else	{  # advection directed upwards
        nom <- a.bl.down*C.bl.down + VF$int[N+1]*(D$int[N+1]/
                 grid$dx.aux[N+1] + AFDW$int[N+1]*v$int[N+1])*C[N]
        denom <- a.bl.down + VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
                (1-AFDW$int[N+1])*v$int[N+1])
	    }
      C.down <- nom/denom
    }

    ## Calculate diffusive part of the flux
	  flux <- as.vector(-(VF$int)*D$int*
                      diff(c(C.up,C,C.down))/grid$dx.aux)
    ## Add advective part of the flux if needed
      if (any (v$int > 0)) {   # advection directed downwards
       vv <- v$int
       vv[v$int<0] <- 0
	     conc <- AFDW$int*c(C.up,C)
	     if (any (AFDW$int < 1))
         conc <- conc +(1-AFDW$int)*c(C,C.down)
	     flux <- flux + as.vector(VF$int*vv*conc)
    }
    if (any (v$int < 0)) {   # advection directed upwards
       vv <- v$int
       vv[v$int>0] <- 0

      conc <- AFDW$int*c(C,C.down)
	    if (any (AFDW$int < 1))
        conc <- conc +(1-AFDW$int)*c(C.up,C)
	     flux <- flux + as.vector(VF$int*vv*conc)
    }

  }

  if (! is.null (flux.up))
    flux[1] <- flux.up
  if (! is.null (flux.down))
    flux[N+1] <- flux.down

    
## Calculate rate of change = Flux gradient       ## Karline: /VF$mid
  dC <- -diff(A$int*flux)/A$mid/grid$dx   / VF$mid

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

