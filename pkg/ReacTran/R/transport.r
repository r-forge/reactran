# WEIGUP added
###############################################################################
# Transport of a solid substance in the sediment
###############################################################################

tran1D.solid <- function (y,                   # Concentration in porous medium
                          y.up=y[1],           # bottomwater concentration 
                          y.down=y[length(y)], # Concentration at deep interface
                          flux.up=NA,        # Deposition flux on sediment interface
                          flux.down=NA,      # Deep burial flux 
                          disp,              # Sediment bioturbation rate, on box interface
                          v=0,               # Sedimentation rate, pos in axis direction
                          v.up=0,            # velocity rate, against the axis
                          weight.up=1,       # Upstream weighing coefficient for advection and upwelling
                          por.INF=NULL,      # Porosity at infinite depth
                          por.int=por.INF ,  # Porosity at box interfaces
                          por.mid=por.INF ,  # Porosity in centre of boxes
                          porgrid = list(INF=por.INF,int=por.int,mid=por.mid),
                          dx=NULL,           # thickness of boxes
                          dx.int=dx,         # distance from centre to centre of boxes 
                          grid=NULL)
                            
#==============================================================================
# transport of a solid substance (mass/cm3 solid), 
#==============================================================================

  {
# use the general 1D transport routine, Taking into account porosity changes.
# Porosity is ~ surface in the function tran1D.
# However: 
# As the advection u is multiplied with (1-por.INF ) (assumes steady-state consolidation),
# it is not possible to set surf.int=(1-por.int). Hence the complexity.

# original flux definition:
#   Flux    <- -disp*(1-por.int)*diff(c(y.up,y,y.down))/grid$dx.int +  #diffusion
#                weight..up* v *(1-por.INF ) * c(y.up,y)              -  #advection, in direction of axis
#                weight.up* v.up *(1-por.int)* c(y,y.down)              #advection, in direction of axis
# check input
if (is.null(grid))    grid <- list(dx=dx, dx.int=dx.int)
if (is.null(porgrid)) porgrid <- list(INF=por.INF,int=por.int,mid=por.mid)
 gn <- names(porgrid)
 if (! "INF"  %in% gn) 
  stop("error: porgrid should be a list that contains 'INF', the porosity at infinite depth")
 if (! "int"  %in% gn) 
  stop("error: porgrid should be a list that contains 'int', the porosity at compartment interfaces")
 if (! "mid"  %in% gn) 
  stop("error: porgrid should be a list that contains 'mid', the porosity at compartment middle")

 if (is.null(porgrid$INF) || is.null(porgrid$int) || is.null(porgrid$mid)) 
    stop("error: porgrid should be a list with (numeric) values") 

# correct v for porosity changes
    v   =v*(1-porgrid$INF )/(1-porgrid$int)
    N       <- length(y)

    Flux    <- -disp * diff(c(y.up,y,y.down))/grid$dx.int +  #diffusion
                weight.up* v     * c(y.up,y)             -   # advection directed in pos 
                weight.up* v.up  * c(y,y.down)               # advection in upstream direction

    # in case advection does not use upstream weighing....
    if (any (weight.up < 1))
    Flux    <- Flux                                                  + 
               (1-weight.up)* v   * c(y,y.down)       -  #advection, in direction of axis
               (1-weight.up)* v.up* c(y.up,y)            #advection, in direction of axis
    
    # if boundary fluxes (total per BULK surface per time) are imposed
    Flux <- as.vector((1-porgrid$int)*Flux)
    if (! is.na (flux.up))   Flux[1]   <- flux.up
    if (! is.na (flux.down)) Flux[N+1] <- flux.down
    
    # Rate of change= Flux gradient, corrected for porosity
    dy    <- -diff(Flux)/(1-porgrid$mid)/grid$dx

    return (list (dy = dy,                 # rate of change in each layer
                  flux=Flux,               # flux across upstream interface
                  flux.up = Flux[1],       # Flux across upstream interface; positive = IN 
                  flux.down = Flux[N+1]))  # Flux across downstream interface; positive = OUT

  } # end tran1d.solid

###############################################################################
# Transport of a dissolved substance in the sediment, including upwelling
###############################################################################

tran1D.solute <- function (y,                             # concentration in sediment
                           y.up=y[1],                     # bottomwater concentration 
                           y.down=y[length(y)],           # Concentration at deep interface
                           flux.up=NA,                    # Flux IN sediment at sediment-water interface
                           flux.down=NA,                  # Deep burial flux 
                           disp,                          # Molecular diffusion coeff, 
                           v=0,                           # Sedimentation rate, in the direction of the axis
                           v.up=0,                        # upwelling rate, against the axis
                           weight.up=1,                   # Upstream weighing coefficient for advection and upwelling
                           por.INF =NULL,                 # Porosity at infinite depth
                           por.int=por.INF ,              # Porosity at box interfaces
                           por.mid=por.INF ,              # Porosity in centre of boxes
                           porgrid = list(INF=por.INF,int=por.int,mid=por.mid),                           
                           tortuosity=1,                  # sediment tortuosity coefficient
                           dx=NULL,                       # thickness of boxes
                           dx.int=dx,                     # distance from centre to centre of boxes 
                           grid=NULL)

#==============================================================================
# transport of a dissolved substance (mass/cm3 liquid), 
#==============================================================================
# Uses the general 1D transport routine, Taking into account porosity changes.
# Porosity is ~ surface in the function tran1D.
# However: 
# As the advection u is multiplied with (por.INF ) (assumes steady-state consolidation),
# it is not possible to set surf.int=(por.int). Hence the higher complexity.

  {
# original flux definition:
    # effective diffusion coeff, accounting for porosity and tortuosity
#    SedDisp <- disp * por.int^(tortuosity-1) 

#    Flux    <- -SedDisp * por.int*diff(c(y.up,y,y.down))/grid$dx.int +  #diffusion
#                weight..up* v * por.INF  *c(y.up,y)            -        #advection
#                weight.up* v.up* por.int *c(y,y.down)   
if (is.null(grid))    grid <- list(dx=dx, dx.int=dx.int)
if (is.null(porgrid)) porgrid <- list(INF=por.INF,int=por.int,mid=por.mid)  
 gn <- names(porgrid)
 if (! "INF"  %in% gn) 
  stop("error: porgrid should be a list that contains 'INF', the porosity at infinite depth")
 if (! "int"  %in% gn) 
  stop("error: porgrid should be a list that contains 'int', the porosity at compartment interfaces")
 if (! "mid"  %in% gn) 
  stop("error: porgrid should be a list that contains 'mid', the porosity at compartment middle")

 if (is.null(porgrid$INF) || is.null(porgrid$int) || is.null(porgrid$mid)) 
    stop("error: porgrid should be a list with (numeric) values") 



    disp=disp*porgrid$int^(tortuosity-1)   #
    v   =v* porgrid$INF/porgrid$int    
    N       <- length(y)
    
    Flux    <- -disp * diff(c(y.up,y,y.down))/grid$dx.int +  #diffusion
                weight.up* v     * c(y.up,y)             -   # advection directed in pos 
                weight.up* v.up  * c(y,y.down)               # advection in upstream direction

    # in case advection does not use upstream weighing....
    if (any (weight.up < 1))
    Flux    <- Flux                                                  + 
               (1-weight.up)* v   * c(y,y.down)       -  #advection, in direction of axis
               (1-weight.up)* v.up* c(y.up,y)            #advection, in direction of axis
    
    # if boundary fluxes (total per BULK surface per time) are imposed
    Flux <- as.vector(porgrid$int*Flux)
    if (! is.na (flux.up))   Flux[1]   <- flux.up
    if (! is.na (flux.down)) Flux[N+1] <- flux.down
    
    # Rate of change= Flux gradient, corrected for porosity
    dy    <- -diff(Flux)/porgrid$mid/grid$dx

    return (list (dy = dy,                 # rate of change in each layer
                  flux=Flux,               # flux across upstream interface
                  flux.up = Flux[1],       # Flux across upstream interface; positive = IN 
                  flux.down = Flux[N+1]))  # Flux across downstream interface; positive = OUT

  }

###############################################################################
# estuarine advective-diffusive transport              #
###############################################################################

tran1D.volume <- function(y,                       # concentration to be transported
                          y.up=y[1],               # upstream boundary concentration
                          y.down=y[length(y)],     # downstream boundary concentration
                          input.up=NA,             # Upstream flux, IN system
                          input.down=NA,           # Downstream flux OUT of system
                          flow=0,                  # flow rate (L3/T)
                          flow.up=0,               # Upward flow, against the axis (L3/T)
                          weight.up=1,             # Upstream weighing coefficient for advection and upward flow
                          bulkdisp,                # bulk dispersion coefficient (L3/T)
                          volume)                  # box volume (L3)

   {
    N      <- length(y)
    Input  <- weight.up* flow*c(y.up,y)         -         # advection direction of axis
              weight.up* flow.up*c(y,y.down)     -         # advection against axis    
              bulkdisp*diff(c(y.up,y,y.down))    # diffusion: note disp here is bulk diffusion coefficient !

    # in case advection does not use upstream weighing....
    if (any (weight.up < 1))
    Input    <- Input                                                  + 
                (1-weight.up)* flow  * c(y,y.down)       -  #advection, in direction of axis
                (1-weight.up)* flow.up* c(y.up,y)            #advection, in direction of axis

    # if boundary fluxes are imposed
    if (! is.na (input.up))   Input[1]   <- input.up
    if (! is.na (input.down)) Input[N+1] <- input.down
    Input <- as.vector(Input)
                  
    dy      <- -diff(Input)/volume

    list(dy=dy,                      # rate of change
         input       = Input,
         input.up    = Input[1],     # total (volume integrated) import upstream
         input.down  =-Input[N+ 1])  # total (volume integrated) import downstream
  }


###############################################################################
# Advective-diffusive transport across variable shape  #
###############################################################################

tran1D <- function(y,                       # concentration to be transported
                   y.up=y[1],               # upstream boundary concentration
                   y.down=y[length(y)],     # downstream boundary concentration
                   flux.up=NA,              # Upstream flux, mass/L^2/time, IN system
                   flux.down=NA,            # Downstream flux, mass/L^2/time, OUT of system
                   disp,                    # Diffusion coeff, 
                   v =0,                    # Velocity in direction of axis
                   v.up =0,                 # Upward velocity, against the axis
                   weight.up=1,             # Upstream weighing coefficient for advection and upwelling
                   surf.mid=1,              # Surface in centre of boxes
                   surf.int=surf.mid,       # Surface at box interfaces
                   dx=NULL,                 # thickness of boxes
                   dx.int=dx,               # distance from centre to centre of boxes 
                   grid=list(dx=dx,         # thickness of boxes, vector length N
                   dx.int=dx.int))

#==============================================================================
# transport of a substance in a body with variable shape
# shape is set by surface.
#==============================================================================

  {
    gn <- names(grid)
    if (! "dx" %in% gn) 
    stop("error: grid should be a list that contains 'dx', the box thicknesses")
    if (! "dx.int" %in% gn) 
    stop("error: grid should be a list that contains 'dx.int', the interface distances")

    N       <- length(y)

    Flux    <- -disp * diff(c(y.up,y,y.down))/grid$dx.int +  #diffusion
                weight.up* v     * c(y.up,y)             -   # advection directed in pos 
                weight.up* v.up  * c(y,y.down)               # advection in upstream direction

    # in case advection does not use upstream weighing....
    if (any (weight.up < 1))
    Flux    <- Flux                                                  + 
               (1-weight.up)* v   * c(y,y.down)       -  #advection, in direction of axis
               (1-weight.up)* v.up* c(y.up,y)            #advection, in direction of axis
    
    # if boundary fluxes (total per surface per time) are imposed
    Flux <- as.vector(Flux)
    if (! is.na (flux.up))   Flux[1]   <- flux.up
    if (! is.na (flux.down)) Flux[N+1] <- flux.down
    
    # Rate of change= Flux gradient, corrected for porosity
    dy    <- -diff(surf.int*Flux)/surf.mid/grid$dx

    return (list (dy = dy,                 # rate of change in each layer
                  flux=Flux,               # flux across upstream interface
                  flux.up = Flux[1],       # Flux across upstream interface; positive = IN 
                  flux.down = Flux[N+1]))  # Flux across downstream interface; positive = OUT
  }

#==============================================================================
# Subdivides a certain depth into boxes
#==============================================================================

setup.grid <-  function(len=NULL,        # total length to divide
                        N=NULL,          # Number of boxes to divide                        
                        dx.1 =NULL,      # size of first grid cell - used in types "exp" and "sine"
                        layer = NULL,    # matrix with the discretisation specifications, in case type = "layer"
                        type = "ct",     # type of grid, one of "ct","exp","sine", "layer"      
                        dx.int.up=NULL,  # distance from mid box 1 to upstream interface
                        dx.int.down=NULL,
                        x.0 = 0)         # position of upstream interface; defaults to 0 (start of axis)
{

# grid sizes 
if (type != "layer")
{
 if (is.null(len))  stop ("Cannot create grid: len, the total length is not specified") 
 if (type != "ct" && is.null(N)) 
  stop ("Cannot create grid: N, the number of boxes is not specified") 
}

if (type %in% c("exp","sine"))
{
if (dx.1 > len)      stop ("Cannot create grid: dx.1 > len")             
if (is.null(dx.1))   stop ("Cannot create grid: dx.1, the size of first grid cell is not specified ")             
}

if(type=="ct")                         # all grid sizes are equal
 {                   
 if (is.null(N) && is.null(dx.1)) 
   stop ("Cannot create grid: neither N, the number of boxes, nor dx.1 the size of boxes is specified") 
 if (is.null(N)) N <- round(len/dx.1)
 if (N < 1) stop ("Cannot create grid: dx.1 > len")             
 dx    <- rep(len/N,N)
} else if (type == "exp") {           
# smallest compartment size (dx.1) is first one
# size increases according to a power function: 2nd box: dx.1*pow, 3rd box: dx.1*pow*pow etc...
 pow   <- uniroot( f<- function(p) dx.1*(1-p^N)/(1-p) - len, interval=c(0.1,10),tol=1e-15)$root
 dx    <- vector (length = N)
 dx[1] <- dx.1
 for (i in 2:N) dx[i]<-dx[i-1]*pow

} else if (type=="sine")
# sine function to divide depth; size near both ends is imposed (first)
# size increases/decreases towards middle according to a sine function
{

Todivide <- len - dx.1*N
Sine     <- sin((0:(N-1))/(N-1)*pi)
Sum      <- sum(Sine)
dx       <- dx.1+(Sine/Sum*Todivide)
} else if (type=="layer")            # several layers with different specifications
 {
ii    <- which (!is.na (layer[,2]))
if (any(layer[ii,2]>layer[ii,3])) stop ("cannot create grid: size of first box (col 2) and total size (col 3) not compatible")
  layer[,1] <- cumsum(layer[,1])
  layer[,3] <- cumsum(layer[,3])
  N   <- layer[nrow(layer),1]
  dat <- rbind(c(1,0),layer[,c(1,3)])

  layer <- rbind(c(1,NA,NA),layer)
  for (i in 1:length(ii)) dat <- rbind(dat,c(layer[ii,1]+1,layer[ii+1,2]))
  intx  <- spline(dat,n=N+1)$y
  dx    <- diff(intx)
 }

# distances
N     <- length(dx)
intx  <- diffinv(dx) + x.0
x     <- colMeans(rbind(intx[-1],intx[-(N+1)]))

dx.int <- c(dx[1]/2,diff(x),dx[N]/2)
if(!is.null(dx.int.up))   dx.int[1] <- dx.int.up
if(!is.null(dx.int.down)) dx.int[N+1] <- dx.int.down

return(list(x      = x,      # distance from origin to centre of boxes, vector of length N
            x.int  = intx,   # distance from origin to upper box interface, vector of length N+1
            dx     = dx  ,   # thickness of boxes, vector length N
            dx.int = dx.int, # distance between centre of adjacent boxes, first and last: half of box size, vector of length N+1
            N      = N))     # number of boxes
}

#==============================================================================
# Calculates a certain property over a grid
#==============================================================================

setup.prop <-  function(y.0=NULL,        # value at upstream interface
                        y.INF=y.0,       # value at infinite depth
                        y.coeff =0,      # rate of increase with distance
                        L = 0,           # if present: size of upstream layer where y=y0
                        xy = NULL,       # 2-column matrix with x-y values
                        type = "spline", # one of "spline", "linear"
                        grid)            # discretisation grid 
{
# check input
 gn <- names(grid)
 if (! "dx"     %in% gn) stop("error: grid should be a list that contains dx")
 if (! "dx.int" %in% gn) stop("error: grid should be a list that contains dx.int")

# data series input or exponential 
if (! is.null(xy))
{
  if (! is.matrix(xy))stop("cannot setup property: xy should be a 2-columned matrix or NULL")
  
  if (type=="linear")
  {
  # check the extent of input data-series, inr: should span range of grid$x.int, outr
  outr <- range(grid$x.int) #output range
  inr  <- range(xy[,1])     # input range
  if (outr[1]<inr[1]) {xy <- rbind(c(outr[1],xy[which.min(xy[,1]),2]),xy)}
  if (outr[2]>inr[2]) {xy <- rbind(xy,c(outr[2],xy[which.max(xy[,1]),2]))}
  y.int <- approx(x=xy[,1],y=xy[,2],xout = grid$x.int)$y
  y.mid <- approx(x=xy[,1],y=xy[,2],xout = grid$x)$y
  
  } else {
  f <- splinefun(xy[,1],xy[,2])
  y.int <- f(grid$x.int)
  y.mid <- f(grid$x)

  }  
  y.0  <-y.int[1]
  y.INF<-y.int[length(y.int)]

} else  # exponentially decaying function
{
 if (is.null(L)) L <- 0
 cc <- -abs(y.coeff)
 y.mid <- y.INF + (y.0-y.INF)*exp((grid$x-L)*cc)      #  profile, middle of layers
 y.int <- y.INF + (y.0-y.INF)*exp((grid$x.int-L)*cc)  #  profile, upper interface
 
 if (L!=0)
 {
 y.mid[grid$x    <=L]<-y.0
 y.int[grid$x.int<=L]<-y.0
 }
}
return(list(up     = y.0,    
            INF    = y.INF,
            mid    = y.mid,
            int    = y.int)) 
}

#==============================================================================
# The fiadeiro and veronis scheme: upstream weighing cofficients
#==============================================================================

fiadeiro <- function(v,  # advection rate, L/t, either one value or a vector of length N+1
                     disp, # diffusion rate, L2/t, either one value or a vector of length N+1
                     dx.int, # distances between centre of boxes (including two interfaces)
                     grid=list(dx.int=dx.int))
                     {
# the weighing scheme to reduce numerical diffusion, due to advection
# sigma are weighing coefficients of upstream compartment;
# they vary from 0.5 (high additional diffusion) to 1 (absence of diffusion)
# reference: Fiadeiro and Veronis, 1977   ! Tellus v 29, pp 512-522
Pe    <- v* grid$dx.int/disp
sigma <- (1+(1/tanh(Pe)-1/Pe))/2
return(sigma)
}

                       