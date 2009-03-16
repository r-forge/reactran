

##==============================================================================
## Advective-diffusive transport in a river
##==============================================================================

tran.volume.1D <- function(C, C.up=C[1], C.down=C[length(C)],
   C.lat=0, C.lat.grid=list(mid=rep(C.lat,length.out=length(C))),
   F.up=NA, F.down=NA, F.lat=NA, F.lat.grid=list(mid=rep(F.lat,length.out=length(C))),
	 Disp=NULL,	Disp.grid=list(int=rep(Disp,length.out=(length(C)+1))),
	 flow = 0, flow.grid=list(int=rep(flow,length.out=(length(C)+1))),
	 flow.lat=0,flow.lat.grid=list(mid=rep(flow.lat,length.out=length(C))),
	 V=NULL,V.grid=list(mid=rep(V,length.out=length(C))))
                            
{

## INPUT CHECKS
## check input of grid

  gn <- names(V.grid)
  if (! "mid" %in% gn)
    stop("error: V.grid should be a list that contains 'mid' ")
  if (is.null(V.grid$mid))
	  stop("error: the argument V.grid should be a list with (numeric) values for 'mid'")

## check input of Disp.grid

  gn <- names(Disp.grid)
  if (! "int" %in% gn)
    stop("error: Disp.grid should be a list that contains 'int', the dispersion coefficient at the grid cell interfaces")
  if (is.null(Disp.grid$int))
    stop("error: the argument Disp.grid should be a list with (numeric) values")
  if (any (Disp.grid$int < 0))
	  stop("error: the bulk dispersion coefficient Disp should always be positive")

## check input of flow.grid
  gn <- names(flow.grid)
  if (! "int" %in% gn)
    stop("error: flow.grid should be a list that contains 'int'")
  if (is.null(flow.grid$int))
    stop("error: the argument flow.grid should be a list with (numeric) values")
  if (any (flow.grid$int < 0) & any (flow.grid$int > 0))
	  stop("error: the discharge flow cannot be both positive and negative within the same domain")

## check input of flow.lat.grid
  gn <- names(flow.lat.grid)
  if (! "mid" %in% gn)
    stop("error: flow.lat.grid should be a list that contains 'mid'")
  if (is.null(flow.lat.grid$mid))
    stop("error: the argument flow.lat.grid should be a list with (numeric) values")

## check input of flow.grid
  gn <- names(C.lat.grid)
  if (! "mid" %in% gn)
    stop("error: the argument C.lat.grid should be a list that contains 'mid'")

## check input of flow.grid
  gn <- names(F.lat.grid)
  if (! "mid" %in% gn)
    stop("error: the argument F.lat.grid should be a list that contains 'mid'")

  if (is.na (F.lat.grid$mid[1]))
    F.lat.grid$mid <- C.lat.grid$mid*flow.lat.grid$mid

## FUNCTION BODY: CALCULATIONS

  N <- length(C)

## Calculate the discharge in each box (using the water balance)

  for (i in 2:(N+1))
  	flow.grid$int[i] = flow.grid$int[i-1] + flow.lat.grid$mid[i-1]

## Calculate diffusive part of the mass flow F

  Dif.F <- as.vector(-Disp.grid$int*diff(c(C.up,C,C.down)))
  Adv.F <- vector(length=(N+1))

  for (i in 1:(N+1)) {
	  if (flow.grid$int[i] > 0) # advection directed downstream
   		Adv.F[i] <- flow.grid$int[i]*c(C.up,C)[i]

  	if (flow.grid$int[i] < 0)  # advection directed downstream
 	  	Adv.F[i] <- flow.grid$int[i]*c(C,C.down)[i+1]
	}


  F <- as.vector(Dif.F + Adv.F)

## Impose boundary fluxes when needed
  if (! is.na (F.up))
    F[1]   <- F.up
  if (! is.na (F.down))
    F[length(F)] <- F.down
    
## Calculate rate of change = Flux gradient + lateral input

  dC <- -diff(F)/V.grid$mid + F.lat.grid$mid/V.grid$mid

  return (list (dC = dC,                 # Rate of change in the centre of each grid cells
                F = F,                   # Flux across at the interface of each grid cell
                F.up = F[1],             # Flux across upstream boundary ; positive = IN
                F.down = F[length(F)],   # Flux across downstream boundary ; positive = OUT
				  	    F.lat = F.lat.grid$mid)) # Flux lateral ; positive = IN

}
