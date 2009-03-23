

##==============================================================================
## Advective-diffusive transport in a river
##==============================================================================

tran.volume.1D <- function(C, C.up=C[1], C.down=C[length(C)],
   C.lat=0, F.up=NULL, F.down=NULL, F.lat=NULL,
	 Disp=NULL, flow = 0, flow.lat=0,	 V=NULL,
   full.check = FALSE, full.output = FALSE)

{

## INPUT CHECKS

  N <- length(C)


  if (!is.list(F.lat))
   F.lat <- list(mid=rep(F.lat,length.out=N))

  if (!is.list(Disp))
   Disp <- list(int=rep(Disp,length.out=N+1))

  if (!is.list(flow))
    flow<-list(int=rep(flow,length.out=(N+1)))

  if (!is.list(flow.lat))
    flow.lat <-list(mid=rep(flow.lat,length.out=N))

  if (!is.list(V))
    V <- list(mid=rep(V,length.out=N))

  if (!is.list(C.lat))
    C.lat <- list(mid=rep(C.lat,length.out=length(C)))

  if (full.check) {
## check input of grid
    gn <- names(V)

    if (! "mid" %in% gn)
      stop("error: V should be a list that contains 'mid' ")

  if (is.null(V$mid))
	  stop("error: the argument V should contain (numeric) values ")

## check input of Disp

  gn <- names(Disp)
  if (! "int" %in% gn)
    stop("error: Disp should be a list that contains 'int', the dispersion coefficient at the grid cell interfaces")
  if (is.null(Disp$int))
    stop("error: the argument Disp should be contain (numeric) values")
  if (any (Disp$int < 0))
	  stop("error: the bulk dispersion coefficient Disp should always be positive")

## check input of flow
  gn <- names(flow)
  if (! "int" %in% gn)
    stop("error: flow should be a list that contains 'int'")
  if (is.null(flow$int))
    stop("error: the argument flow should contain (numeric) values")
  if (any (flow$int < 0) & any (flow$int > 0))
	  stop("error: the discharge flow cannot be both positive and negative within the same domain")

## check input of flow.lat
  gn <- names(flow.lat)
  if (! "mid" %in% gn)
    stop("error: flow.lat should be a list that contains 'mid'")
  if (is.null(flow.lat$mid))
    stop("error: the argument flow.lat should contain (numeric) values")

## check input of flow
  gn <- names(C.lat)
  if (! "mid" %in% gn)
    stop("error: the argument C.lat should be a list that contains 'mid'")

## check input of flow
  gn <- names(F.lat)
  if (! "mid" %in% gn)
    stop("error: the argument F.lat should be a list that contains 'mid'")
  }
  if (is.null (F.lat$mid[1]))
    F.lat$mid <- C.lat$mid*flow.lat$mid

## FUNCTION BODY: CALCULATIONS

## Calculate the discharge in each box (using the water balance)
  f1 <- flow$int[1]
  flow$int <- c(f1,f1+cumsum(flow.lat$mid))

## Calculate diffusive part of the mass flow F

  Dif.F <- as.vector(-Disp$int*diff(c(C.up,C,C.down)))

# Advection: positive flows first
  v1 <- flow$int

  if ( any (v1 > 0 )) {
    v1[v1<0] <- 0
	  Adv.F <- v1 *c(C.up,C)
  } else Adv.F <- 0
  
# If there are negative flows:
  if ( any (v1 < 0 )) {
    v1 <- flow$int
    v1[v1>0] <- 0
  	Adv.F <- Adv.F + v1 *c(C,C.down)
  }


  F <- as.vector(Dif.F + Adv.F)

## Impose boundary fluxes when needed
  if (! is.null (F.up))
    F[1]   <- F.up
  if (! is.null (F.down))
    F[length(F)] <- F.down
    
## Calculate rate of change = Flux gradient + lateral input

  dC <- -diff(F)/V$mid + F.lat$mid/V$mid

    if (!full.output){
    return (list (dC = dC,               # Rate of change due to advective-diffuisve transport in each grid cell
                F.up = F[1],             # Flux across upstream boundary ; positive = IN
                F.down = F[length(F)]    # Flux across downstream boundary ; positive = OUT
                ))
  } else {
    return (list (dC = dC,                 # Rate of change in the centre of each grid cells
                F = F,                   # Flux across at the interface of each grid cell
                F.up = F[1],             # Flux across upstream boundary ; positive = IN
                F.down = F[length(F)],   # Flux across downstream boundary ; positive = OUT
				  	    F.lat = F.lat$mid)) # Flux lateral ; positive = IN
  }
}
