###############################################################################
# Creation of a one-dimensional finite difference grid
###############################################################################

setup.grid.1D <- function(x.up=0,					# position of the upstream interface	
												  x.down=NULL,    # position of zonation endpoints, one value (position of downstream interface) or a vector of length N
												  L=NULL,         # thickness of zones, one value (model domain = one zone) or a vector of length N
                          N=NULL,         # desired number of grid cells in a zone, one value or a vector of length N 
                          dx.1 =NULL,     # size of first grid cell in a zone, one value or a vector of length N 
												  p.dx.1 = rep(1,length(L)), # grid cell size growth factor in upper half of the zone, one value or a vector of length N  
												  max.dx.1 = 10*dx.1, # maximum grid cell size in upper half of the zone, one value or a vector of length N      
                          dx.N =NULL,     # size of last grid cell in a zone, one value or a vector of length N 
												  p.dx.N = rep(1,length(L)), # grid cell size growth factor in lower half of the zone, one value or a vector of length N  
												  max.dx.N = 10*dx.N) # maximum grid cell size in lower half of the zone, one value or a vector of length N
{

# Check on the input

if (is.null(x.down[1]) && is.null(L[1]))  stop ("Cannot create grid: the length or end point is not specified") 
if (is.null(x.down[1]))
 {
	 x.down <- vector(length=length(L))
   x.check <- x.up 
	 for (i in 1:length(L))
   {
   if (L[i] < 0) stop ("Cannot create grid: L < 0")
   x.down[i] <- x.check + L[i]
   x.check <- x.down[i]
   }
 }
if (is.null(L[1])) 
 {
	 L <- vector(length=length(x.down))
   x.check <- x.up 
   for (i in 1:length(x.down))
   {
	 if (x.down[i] < x.check) stop ("Cannot create grid: x.down < x.up")
	 L[i] <- x.down[i] - x.check
   x.check <- x.down[i]
   }
 }

# Calculation of the grid cell sizes 

dx <- vector()
# loop over all zones
for (i in 1:length(L))
{
 if (is.null(N[i]) && is.null(dx.1[i])) 
		stop ("Cannot create grid: neither N (the number of boxes) nor dx.1 (the size of the first grid cell) is specified") 

 if ((is.null(dx.1[i])) && (is.null(dx.N[i]))) 
 # all grid cells are equal
 {
 if (N[i] < 1) stop ("Cannot create grid: N < 0")             
 A.dx  <- rep(L[i]/N[i],N[i])
 }
 else if (is.null(dx.N[i])){
 if (is.null(dx.1[i])) stop ("Cannot create grid: dx.1 not specified")
 if (dx.1[i] > L[i]) stop ("Cannot create grid: dx.1 > L")             
 # use gradual increase of grid cell size at uper interface alone
 A.dx <- vector()
 pos <- vector()
 A.dx[1] <- dx.1[i] 
 pos[1] <- dx.1[i]
 j <- 1
 while ((is.null(N[i])&&(pos[j] < L[i]))||(!is.null(N[i])&&(j<N[i])))
 {
 j <- j + 1
 A.dx[j] <- min(max.dx.1[i],A.dx[j-1]*p.dx.1[i])
 pos[j] <- pos[j-1] + A.dx[j]
 }
 A.dx <- A.dx*(L[i]/pos[length(pos)]) 
} else { 
 if (is.null(dx.1[i])) stop ("Cannot create grid: dx.1 not specified")
 # use gradual increase of grid cell size at both interfaces
 A.dx  <- vector()
 pos <- vector()
 A.dx[1] <- dx.1[i] 
 pos[1] <- dx.1[i]
 j <- 1
 if (!is.null(N[i])) N.half <- as.integer(N[i]/2)
 while ((is.null(N[i])&&(pos[j] < 0.5*L[i]))||(!is.null(N[i])&&(j<N.half)))
 {
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
 while ((is.null(N[i])&&(pos[j] < 0.5*L[i]))||(!is.null(N[i])&&(j<(N[i]-N.half))))
 {
 j <- j + 1
 B.dx[j] <- min(max.dx.N[i],B.dx[j-1]*p.dx.N[i])
 pos[j] <- pos[j-1] + B.dx[j]
 }
 B.dx <- B.dx*(0.5*L[i]/pos[length(pos)]) 

 A.dx <- c(A.dx,B.dx[length(B.dx):1])
}
dx <- c(dx,A.dx)
}
# positions and distances

x.int  <- x.up + diffinv(dx) 
x.mid  <- colMeans(rbind(x.int[-1],x.int[-(length(dx)+1)]))
dx.aux <- c(dx[1]/2,diff(x.mid),dx[length(dx)]/2)

# Packaging of results 

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

###############################################################################
# Plotting of a one-dimensional finite difference grid
###############################################################################

plot.grid.1D <- function(x,...)
{
   mf <- par(mfrow=c(2,1))
   on.exit(par(mf))
    
   plot(x$x.mid,main="position of cells",ylab="x.mid",xlab="index",...)
   plot(x$dx,main="box thickness",ylab="dx",xlab="index",...)
}

###############################################################################
# Creation of a two-dimensional finite difference grid
###############################################################################

setup.grid.2D <- function(x.grid=NULL,		# one-dimensional grid - horizontal 
												  y.grid=NULL)    # one-dimensional grid - vertical 
{

# check input
gn <- names(x.grid)
if (! "x.up"   %in% gn) stop("error in setup.2Dgrid: x.grid should be a list that contains x.up")
if (! "x.down" %in% gn) stop("error in setup.2Dgrid: x.grid  should be a list that contains x.down")
if (! "x.mid"  %in% gn) stop("error in setup.2Dgrid: x.grid  should be a list that contains x.mid")
if (! "x.int"  %in% gn) stop("error in setup.2Dgrid: x.grid  should be a list that contains x.int")
if (! "dx"     %in% gn) stop("error in setup.2Dgrid: x.grid  should be a list that contains dx")
if (! "dx.aux" %in% gn) stop("error in setup.2Dgrid: x.grid  should be a list that contains dx.aux")
if (! "N"      %in% gn) stop("error in setup.2Dgrid: x.grid  should be a list that contains N")

gn <- names(y.grid)
if (! "x.up"   %in% gn) stop("error in setup.2Dgrid: y.grid should be a list that contains x.up")
if (! "x.down" %in% gn) stop("error in setup.2Dgrid: y.grid  should be a list that contains x.down")
if (! "x.mid"  %in% gn) stop("error in setup.2Dgrid: y.grid  should be a list that contains x.mid")
if (! "x.int"  %in% gn) stop("error in setup.2Dgrid: y.grid  should be a list that contains x.int")
if (! "dx"     %in% gn) stop("error in setup.2Dgrid: y.grid  should be a list that contains dx")
if (! "dx.aux" %in% gn) stop("error in setup.2Dgrid: y.grid  should be a list that contains dx.aux")
if (! "N"      %in% gn) stop("error in setup.2Dgrid: y.grid  should be a list that contains N")

# Packaging of results 

Res <- list(x.up = x.grid$x.up,
				    x.down = x.grid$x.down,
				    x.mid = x.grid$x.mid,   # position of centre of the grid cells, vector of length N
            x.int  = x.grid$x.int,  # position grid cell interfaces , vector of length N+1
            dx = x.grid$dx,         # thickness of grid cells , vector length N
            dx.aux = x.grid$dx.aux, # auxiliary vector with distances between centre of adjacent cells, first and last: half of cell size, vector of length N+1
            x.N = x.grid$N,         # number of vertical grid layers
            y.left = y.grid$x.up,
				    y.right = y.grid$x.down,
				    y.mid = y.grid$x.mid,   # position of centre of the grid cells, vector of length N
            y.int  = y.grid$x.int,  # position grid cell interfaces , vector of length N+1
            dy = y.grid$dx,         # thickness of grid cells , vector length N
            dy.aux = y.grid$dx.aux, # auxiliary vector with distances between centre of adjacent cells, first and last: half of cell size, vector of length N+1
            y.N = y.grid$N)         # number of horizontal grid layers

class(Res) <- "grid.2D"
return(Res)
}

###############################################################################
# Attaches a property to a 1D grid
###############################################################################

setup.prop.1D <- function(func = NULL,     # function 
												  value = NULL,    # constant value
                          xy = NULL,       # 2-column matrix with x-y values
												  interpolate = "spline", # one of "spline", "linear"
												  grid,...)        # grid 
{
# check input
gn <- names(grid)
if (! "x.up"   %in% gn) stop("error in setup.1Dprop: grid should be a list that contains x.up")
if (! "x.down" %in% gn) stop("error in setup.1Dprop: grid should be a list that contains x.down")
if (! "x.mid"  %in% gn) stop("error in setup.1Dprop: grid should be a list that contains x.mid")
if (! "x.int"  %in% gn) stop("error in setup.1Dprop: grid should be a list that contains x.int")
if (! "dx"     %in% gn) stop("error in setup.1Dprop: grid should be a list that contains dx")
if (! "dx.aux" %in% gn) stop("error in setup.1Dprop: grid should be a list that contains dx.aux")
if (! "N"      %in% gn) stop("error in setup.1Dprop: grid should be a list that contains N")


if (is.null(xy) && is.null(value) && is.null(func)) stop("error in setup.prop: function, value and xy cannot be NULL together")

if (! is.null(func))
# profile specification via function 
{
	y.int <- func(grid$x.int,...)
	y.mid <- func(grid$x.mid,...)
}

if (! is.null(value))
# profile specification via value
{
	y.int <- rep(value,times=length(grid$x.int))
	y.mid <- rep(value,times=length(grid$x.mid))
}

if (! is.null(xy))
# profile speficication via data series input 
{
  if (! is.matrix(xy)) stop("error in setup.prop: xy should be a 2-columned matrix or NULL")
	if (!(interpolate %in% c("linear","spline"))) stop("error in setup.prop: <interpolate> not properly specified, should be <linear> or <spline>")
		
  if (interpolate=="linear")
  {
  # check the range of input data-series (inr))
  # this range should span the range of grid$x.int (outr)
  outr <- range(grid$x.int) # output range
  inr  <- range(xy[,1])     # input range
	# extend the xy matrix if necessary 
	if (outr[1]<inr[1]) {xy <- rbind(c(outr[1],xy[which.min(xy[,1]),2]),xy)}
  if (outr[2]>inr[2]) {xy <- rbind(xy,c(outr[2],xy[which.max(xy[,1]),2]))}
  # linear interpolation
	y.int <- approx(x=xy[,1],y=xy[,2],xout = grid$x.int)$y
  y.mid <- approx(x=xy[,1],y=xy[,2],xout = grid$x.mid)$y
  } else {

	# spline interpolation
	f <- splinefun(xy[,1],xy[,2])
  y.int <- f(grid$x.int)
  y.mid <- f(grid$x.mid)
  }  
}  

Res <- list(mid  = y.mid,
            int  = y.int) 
class(Res) <- "prop.1D"
return(Res)
}

###############################################################################
# Attaches a property to a 2D grid
###############################################################################

setup.prop.2D <- function(func = NULL,     # function 
												  value = NULL,
                          grid,...)        # 2D grid 
{
# check input
gn <- names(grid)
if (! "x.mid"  %in% gn) stop("error in setup.2D.prop: grid should be a list that contains x.mid")
if (! "x.int"  %in% gn) stop("error in setup.2D.prop: grid should be a list that contains x.int")
if (! "y.mid"  %in% gn) stop("error in setup.2D.prop: grid should be a list that contains x.mid")
if (! "y.int"  %in% gn) stop("error in setup.2D.prop: grid should be a list that contains x.int")

if (is.null(func) && is.null(value)) stop("error in setup.prop: function and value should not be both NULL")

if (!is.null(value))
# profile specification via value
{
x.int <- matrix(nrow=length(grid$x.int),ncol=length(grid$y.mid),data=value)
y.int <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.int),data=value) 
mid <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.mid),data=value)
}

if (!is.null(func))
# profile specification via function 
{

x.int <- matrix(nrow=length(grid$x.int),ncol=length(grid$y.mid))
y.int <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.int))
mid <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.mid))

for (i in 1:length(grid$x.int)) x.int[i,] <- func(grid$x.int[i],grid$y.mid,...)
for (i in 1:length(grid$x.mid)) y.int[i,] <- func(grid$x.mid[i],grid$y.int,...)
for (i in 1:length(grid$x.mid)) mid[i,] <- func(grid$x.mid[i],grid$y.mid,...)
}

Res <- list(mid  = mid,
            x.int  = x.int,
            y.int =  y.int) 
class(Res) <- "prop.2D"
return(Res)
}

###############################################################################
# The fiadeiro and veronis scheme: upstream weighing cofficients
###############################################################################

fiadeiro <- function(v, # advective velocity, L/t, either one value or a vector of length N+1
                     D, # diffusion coefficient , L2/t, either one value or a vector of length N+1
                     dx.aux, # auxiliary vector containing the distances between the locations where the concentration is defined (i.e. the grid cell centers and the two interfaces), either one value or a vector of length N+1
                     grid=list(dx.aux=dx.aux))
                     {
# the weighing scheme to reduce numerical diffusion, due to advection
# sigma are AFDW weighing coefficients of upstream compartment;
# AFDW vary from 0.5 (high additional diffusion, centred difference) to 1 (absence of diffusion, backward difference)
# reference: Fiadeiro and Veronis, 1977   ! Tellus v 29, pp 512-522
Pe    <- abs(v)*grid$dx.aux/D
sigma <- (1+(1/tanh(Pe)-1/Pe))/2
return(sigma)
}

###############################################################################
# Calculates the advective velocities assuming steady state compaction
###############################################################################

setup.compaction.1D <- function(v.0 = NULL,   # advective velocity at the sediment-water interface
		                   v.inf = NULL, # addvective velocity at infinite depth 
											 por.0,    # porosity at the sediment-water interface 
											 por.inf,  # porosity at infinite depth 
											 por.grid) # porosity profile (grid property) 
{
# check input
gn <- names(por.grid)
if (! "mid"  %in% gn) stop("error in setup.prop: porosity should be a list that contains mid")
if (! "int"  %in% gn) stop("error in setup.prop: porosity should be a list that contains int")

if (is.null(v.0) && is.null(v.inf)) stop("error in setup.advection: either the sedimentation velocity <v.0> or the burial velocity <v.inf> should be specified")

if (is.null(v.inf))
{
# sedimentation velocity is specified 	
	v.factor <-  v.0*(1-por.0)
	v.inf <- v.factor/(1-por.inf) 
	u.inf <- v.inf
	u.factor <-  u.inf*por.inf
} else {
# burial velocity is specified 	
	v.factor <-  v.inf*(1-por.inf)
	v.0 <- v.factor/(1-por.0) 
	u.inf <- v.inf
	u.factor <-  u.inf*por.inf
	}

v.mid <- v.factor/(1-por.grid$mid)
v.int <- v.factor/(1-por.grid$int)

u.mid <- u.factor/por.grid$mid
u.int <- u.factor/por.grid$int

Res <- list(u=list(mid=u.mid,int=u.int),v=list(mid=v.mid,int=v.int))
class(Res) <- "compaction.1D"
return(Res)
}
###############################################################################
# Transport in a one-dimensional finite difference grid 
###############################################################################

tran.1D <- function(C,                   # Concentration in solid phase of porous medium
                    C.up=C[1],           # Concentration at the upper interface
                    C.down=C[length(C)], # Concentration at lower interface
                    flux.up=NA,          # Flux at the upper interface (sediment-water interface)
                    flux.down=NA,        # Flux at the lower interface (deep sediment interface)
									  a.bl.up=NULL,        # transfer coefficient across the upstream boundary layer; Flux = a.bl.up*(C.bl.up-C[1])
									  C.bl.up=NULL,        # concentration at the top of the upstream boundary layer 
									  a.bl.down=NULL,      # transfer coefficient across the downstream boundary layer; Flux = a.bl.down*(C[N]-C.bl.down)
									  C.bl.down=NULL,      # concentration at the top of the downstream boundary layer 
                    D=NULL,              # Diffusion coefficient (positive), defined on the grid cell interfaces, one value or a vector of length N+1 
			              D.grid=NULL,         # Diffusion coefficient D defined as grid property
									  v=0,                 # Advective velocity (positive or negative), defined on the grid cell interfaces, one value or a vector of length N+1
									  v.grid=NULL,         # Advective velocity v defined as grid property
  	                AFDW=1,              # Advective Finite Difference Weight, weighing coefficient in finite difference scheme for advection (1:backward, 0:forward, 0.5:centred) 
									  AFDW.grid=NULL,      # Advective Finite Difference Weight, defined as grid property
                    VF=1,                # Interface area defined at grid cell interfaces, one value or a vector of length N+1
                    VF.grid=NULL,        # volume fraction VF defined as a grid property
                    A=1,                 # Interface area defined at grid cell interfaces, one value or a vector of length N+1
                    A.grid=NULL,         # surface area A defined as a grid property
                    dx=NULL,             # thickness of grid cells, one value or a vector of length N 
                    grid=NULL,           # grid definition
									  full.check = FALSE,  # logical flag enabling a full check of the consistency of the arguments (FALSE speeds up execution by 50%)
									  full.output = FALSE) # logical flag enabling a full return of the output (FALSE speeds up execution by 20%)
											
{
#==============================================================================
# INPUT CHECKS  
#==============================================================================


#==============================================================================
# DEFAULT INFILLING OF GRID PARAMETERS 
#==============================================================================

if (is.null(AFDW.grid)) AFDW.grid <- list(int=AFDW)
if (is.null(D.grid)) D.grid <- list(int=D)
if (is.null(v.grid)) v.grid <- list(int=v)
if (is.null(VF.grid)) VF.grid <- list(int=rep(VF,length.out=(length(C)+1)),mid=0.5*(rep(VF,length.out=(length(C)+1))[1:(length(C))]+rep(VF,length.out=(length(C)+1))[2:(length(C)+1)]))
if (is.null(A.grid)) A.grid <- list(int=rep(A,length.out=(length(C)+1)),mid=0.5*(rep(A,length.out=(length(C)+1))[1:(length(C))]+rep(A,length.out=(length(C)+1))[2:(length(C)+1)]))
if (is.null(grid)) grid <- list(dx=rep(dx,length.out=length(C)),dx.aux=0.5*(c(0,rep(dx,length.out=length(C)))+c(rep(dx,length.out=length(C)),0)))

# check dimensions of input arguments

if (full.check)
{

if (!is.null(AFDW))
{
 if (!((length(AFDW)==1) || (length(AFDW)==(length(C)+1))))
  stop("error: AFDW should be a vector of length 1 or N+1")
} 

if (!is.null(D))
{
 if (!((length(D)==1) || (length(D)==(length(C)+1))))
  stop("error: D should be a vector of length 1 or N+1")
} 

if (!is.null(v))
{
 if (!((length(v)==1) || (length(v)==(length(C)+1))))
  stop("error: v should be a vector of length 1 or N+1")
} 

if (!is.null(VF))
{
 if (!((length(VF)==1) || (length(VF)==(length(C)+1))))
  stop("error: VF should be a vector of length 1 or N+1")
} 

if (!is.null(A))
{
 if (!((length(A)==1) || (length(A)==(length(C)+1))))
  stop("error: A should be a vector of length 1 or N+1")
} 

if (!is.null(dx))
{
 if (!((length(dx)==1) || (length(dx)==(length(C)))))
  stop("error: dx should be a vector of length 1 or N+1")
} 

# check input of grid

gn <- names(grid)
if (! "dx" %in% gn) 
  stop("error: grid should be a list that contains 'dx' ")
if (! "dx.aux" %in% gn) 
	stop("error: grid should be a list that contains 'dx.aux' ")
if (is.null(grid$dx) || is.null(grid$dx.aux)) 
	stop("error: the grid should be a list with (numeric) values for 'dx' and 'dx.aux' ") 
if (any(grid$dx <= 0) || any(grid$dx.aux <= 0) )
	stop("error: the grid distances dx and dx.aux should always be positive") 

# check input of AFDW.grid

gn <- names(AFDW.grid)
if (! "int" %in% gn) 
  stop("error: AFDW.grid should be a list that contains 'int', the AFDW values at the interface of the grid cell ")
if (is.null(AFDW.grid$int)) 
  stop("error: AFDW.grid should be a list with (numeric) values") 
if (any(AFDW.grid$int < 0)||any(AFDW.grid$int > 1))
	stop("error: the AFDW should always range between 0 and 1") 

# check input of D.grid

gn <- names(D.grid)
 if (! "int" %in% gn) 
  stop("error: D.grid should be a list that contains 'int', the D values at the interface of the grid cell ")
if (is.null(D.grid$int)) 
  stop("error: D.grid should be a list with (numeric) values") 
if (any(D.grid$int < 0))
	stop("error: the mixing coefficient should always be positive") 

# check input of v.grid

gn <- names(v.grid)
if (! "int" %in% gn) 
  stop("error: v.grid should be a list that contains 'int', the v values at the interface of the grid cell ")
if (is.null(v.grid$int)) 
  stop("error: the advective velocity v should be a list with (numeric) values") 
if (any (v.grid$int < 0) & any (v.grid$int > 0))
	stop("error: the advective velocity cannot be both positive and negative within the same domain") 

# check input of VF.grid

gn <- names(VF.grid)
if (! "int" %in% gn) 
  stop("error: VF.grid should be a list that contains 'int', the area values at the interface of the grid cell ")
if (! "mid" %in% gn) 
  stop("error: VF.grid should be a list that contains 'mid', the area at the middle of the grid cells")
if (is.null(VF.grid$int) || is.null(VF.grid$mid)) 
  stop("error: the volume fraction VF should be a list with (numeric) values") 
if (any (VF.grid$int < 0) || any (VF.grid$mid < 0))
  stop("error: the volume fraction should always be positive") 

# check input of A.grid

gn <- names(A.grid)
if (! "int" %in% gn) 
  stop("error: A.grid should be a list that contains 'int', the area values at the interface of the grid cell ")
if (! "mid" %in% gn) 
  stop("error: A.grid should be a list that contains 'mid', the area at the middle of the grid cells")
if (is.null(A.grid$int) || is.null(A.grid$mid)) 
  stop("error: the volume fraction VF should be a list with (numeric) values") 
if (any (A.grid$int < 0) || any (A.grid$mid < 0))
  stop("error: the area A should always be positive") 

# check input of boundary layer parameters
if (!is.null(a.bl.up) & !is.null(C.bl.up))
{
	if (a.bl.up < 0) stop("error: the boundary layer transfer coefficient should be positive") 
}

if (!is.null(a.bl.down) & !is.null(C.bl.down))
{
	if (a.bl.down < 0) stop("error: the boundary layer transfer coefficient should be positive") 
}

}

#==============================================================================
# FUNCTION BODY: CALCULATIONS 
#==============================================================================

N <- length(C)

if (full.output)
{
# Impose boundary flux at upper boundary when needed
# Default boundary condition is no gradient
if (! is.na (flux.up))
{
  if (any (v.grid$int >= 0))
	{
# advection directed downwards 
  nom <- flux.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + (1-AFDW.grid$int[1])*v.grid$int[1])*C[1]
  denom <- VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1]+AFDW.grid$int[1]*v.grid$int[1])
	} else	{
# advection directed upwards
  nom <- flux.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + AFDW.grid$int[1]*v.grid$int[1])*C[1]
  denom <- VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1]+(1-AFDW.grid$int[1])*v.grid$int[1])
	}
  C.up <- nom/denom 
}
# Impose boundary flux at lower boundary when needed
# Default boundary condition is no gradient
if (! is.na (flux.down)) 
{
  if (any (v.grid$int >= 0))
  {
# advection directed downwards 
	nom <- flux.down - VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + AFDW.grid$int[N+1]*v.grid$int[N+1])*C[N] 
  denom <- -VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1]+(1-AFDW.grid$int[N+1])*v.grid$int[N+1])
  } else {
# advection directed downwards 
	nom <- flux.down - VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + (1-AFDW.grid$int[N+1])*v.grid$int[N+1])*C[N] 
  denom <- -VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1]+AFDW.grid$int[N+1]*v.grid$int[N+1])
  }	
  C.down <- nom/denom 
}

# when upstream boundary layer is present, calculate new C.up 
if (!is.null(a.bl.up) & !is.null(C.bl.up))
{
  if (any (v.grid$int >= 0))
	{
# advection directed downwards 
  nom <- a.bl.up*C.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + (1-AFDW.grid$int[1])*v.grid$int[1])*C[1]
  denom <- a.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + AFDW.grid$int[1]*v.grid$int[1])
	} else	{
# advection directed upwards
  nom <- a.bl.up*C.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + AFDW.grid$int[1]*v.grid$int[1])*C[1]
  denom <- a.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + (1-AFDW.grid$int[1])*v.grid$int[1])
	}
  C.up <- nom/denom 
}

# when downstream boundary layer is present, calculate new C.up 
if (!is.null(a.bl.down) & !is.null(C.bl.down))
{
  if (any (v.grid$int >= 0))
	{
# advection directed downwards 
  nom <- a.bl.down*C.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + (1-AFDW.grid$int[N+1])*v.grid$int[N+1])*C[N]
  denom <- a.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + AFDW.grid$int[N+1]*v.grid$int[N+1])
	} else	{
# advection directed upwards
  nom <- a.bl.down*C.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + AFDW.grid$int[N+1]*v.grid$int[N+1])*C[N]
  denom <- a.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + (1-AFDW.grid$int[N+1])*v.grid$int[N+1])
	}
  C.down <- nom/denom 
}

# Calculate diffusive part of the flux 

dif.flux <- as.vector(-VF.grid$int*D.grid$int*diff(c(C.up,C,C.down))/grid$dx.aux)  
adv.flux <- rep(0,length.out=length(dif.flux))

# Add advective part of the flux if needed
if (any (v.grid$int > 0))
# advection directed downwards 
{
	conc <- AFDW.grid$int*c(C.up,C)
	if (any (AFDW.grid$int < 1)) conc <- conc +(1-AFDW.grid$int)*c(C,C.down)
	adv.flux <- as.vector(VF.grid$int*v.grid$int*conc) 
}
if (any (v.grid$int < 0))
# advection directed upwards 
{
	conc <- AFDW.grid$int*c(C,C.down)
	if (any (AFDW.grid$int < 1)) conc <- conc +(1-AFDW.grid$int)*c(C.up,C)
	adv.flux <- as.vector(VF.grid$int*v.grid$int*conc) 
}

flux <- dif.flux + adv.flux

} else {

# when upstream boundary layer is present, calculate new C.up 
if (!is.null(a.bl.up) & !is.null(C.bl.up))
{
  if (any (v.grid$int >= 0))
	{
# advection directed downwards 
  nom <- a.bl.up*C.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + (1-AFDW.grid$int[1])*v.grid$int[1])*C[1]
  denom <- a.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + AFDW.grid$int[1]*v.grid$int[1])
	} else	{
# advection directed upwards
  nom <- a.bl.up*C.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + AFDW.grid$int[1]*v.grid$int[1])*C[1]
  denom <- a.bl.up + VF.grid$int[1]*(D.grid$int[1]/grid$dx.aux[1] + (1-AFDW.grid$int[1])*v.grid$int[1])
	}
  C.up <- nom/denom 
}

# when upstream boundary layer is present, calculate new C.up 
if (!is.null(a.bl.down) & !is.null(C.bl.down))
{
  if (any (v.grid$int >= 0))
	{
# advection directed downwards 
  nom <- a.bl.down*C.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + (1-AFDW.grid$int[N+1])*v.grid$int[N+1])*C[N]
  denom <- a.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + AFDW.grid$int[N+1]*v.grid$int[N+1])
	} else	{
# advection directed upwards
  nom <- a.bl.down*C.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + AFDW.grid$int[N+1]*v.grid$int[N+1])*C[N]
  denom <- a.bl.down + VF.grid$int[N+1]*(D.grid$int[N+1]/grid$dx.aux[N+1] + (1-AFDW.grid$int[N+1])*v.grid$int[N+1])
	}
  C.down <- nom/denom 
}

# Calculate diffusive part of the flux 
	flux <- as.vector(-(VF.grid$int)*D.grid$int*diff(c(C.up,C,C.down))/grid$dx.aux)  
	
# Add advective part of the flux if needed
	if (any (v.grid$int > 0))
# advection directed downwards 
	{
		conc <- AFDW.grid$int*c(C.up,C)
		if (any (AFDW.grid$int < 1)) conc <- conc +(1-AFDW.grid$int)*c(C,C.down)
		flux <- flux + as.vector((VF.grid$int)*v.grid$int*conc) 
	}
	if (any (v.grid$int < 0))
# advection directed upwards 
	{
		conc <- AFDW.grid$int*c(C,C.down)
		if (any (AFDW.grid$int < 1)) conc <- conc +(1-AFDW.grid$int)*c(C.up,C)
		flux <- flux + as.vector((VF.grid$int)*v.grid$int*conc) 
	}
	
}

if (! is.na (flux.up)) flux[1] <- flux.up
if (! is.na (flux.down)) flux[N+1] <- flux.down

    
# Calculate rate of change = Flux gradient
dC <- -diff(A.grid$int*flux)/A.grid$mid/grid$dx

if (!full.output){
return (list (dC = dC))                    # Rate of change due to advective-diffuisve transport in each grid cell 
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

###############################################################################
# Advective-diffusive transport in a river
###############################################################################

tran.volume.1D <- function(C,                   # Tracer concentration 
                          C.up=C[1],           # Tracer concentration at the upstream interface, one value 
                          C.down=C[length(C)], # Tracer concentration at downstream interface, one value 
													C.lat=0,             # Tracer concentration of lateral input, one value or a vector of length N 
													C.lat.grid=list(mid=rep(C.lat,length.out=length(C))),  # Tracer concentration in lateral input defined as grid property
													F.up=NA,             # Tracer input at the upstream interface 
                          F.down=NA,           # Tracer output at the downstream interface 
													F.lat=NA,            # Lateral tracer input, one value or a vector of length N 
													F.lat.grid=list(mid=rep(F.lat,length.out=length(C))),  # Lateral tracer input defined as grid property
													Disp=NULL,           # Dispersion coefficient (positive), defined on the grid cell interfaces, one value or a vector of length N+1 
													Disp.grid=list(int=rep(Disp,length.out=(length(C)+1))),  # Dispersion coefficient defined as grid property
											    flow = 0,            # Water flow at grid interfaces,one value or a vector of length N+1  
													flow.grid=list(int=rep(flow,length.out=(length(C)+1))), # Flow defined at grid interfaces 
													flow.lat=0,          # Lateral water input, one value or a vector of length N  
													flow.lat.grid=list(mid=rep(flow.lat,length.out=length(C))),  # Lateral water input defined as grid property
													V=NULL,              # Volume of grid cells, one value or a vector of length N
													V.grid=list(mid=rep(V,length.out=length(C)))) # Volume of grid cells defined as grid property
                            
{
#==============================================================================
# INPUT CHECKS  
#==============================================================================
# check input of grid

gn <- names(V.grid)
if (! "mid" %in% gn) 
  stop("error: V.grid should be a list that contains 'mid' ")
if (is.null(V.grid$mid)) 
	stop("error: the argument V.grid should be a list with (numeric) values for 'mid'") 

# check input of Disp.grid

gn <- names(Disp.grid)
if (! "int" %in% gn) 
  stop("error: Disp.grid should be a list that contains 'int', the dispersion coefficient at the grid cell interfaces")
if (is.null(Disp.grid$int)) 
  stop("error: the argument Disp.grid should be a list with (numeric) values") 
if (any (Disp.grid$int < 0))
	stop("error: the bulk dispersion coefficient Disp should always be positive") 

# check input of flow.grid
gn <- names(flow.grid)
if (! "int" %in% gn) 
  stop("error: flow.grid should be a list that contains 'int'")
if (is.null(flow.grid$int)) 
  stop("error: the argument flow.grid should be a list with (numeric) values") 
if (any (flow.grid$int < 0) & any (flow.grid$int > 0))
	stop("error: the discharge flow cannot be both positive and negative within the same domain") 

# check input of flow.lat.grid
gn <- names(flow.lat.grid)
if (! "mid" %in% gn) 
  stop("error: flow.lat.grid should be a list that contains 'mid'")
if (is.null(flow.lat.grid$mid)) 
  stop("error: the argument flow.lat.grid should be a list with (numeric) values") 

# check input of flow.grid
gn <- names(C.lat.grid)
if (! "mid" %in% gn) 
  stop("error: the argument C.lat.grid should be a list that contains 'mid'")

# check input of flow.grid
gn <- names(F.lat.grid)
if (! "mid" %in% gn) 
  stop("error: the argument F.lat.grid should be a list that contains 'mid'")

if (is.na (F.lat.grid$mid[1])) F.lat.grid$mid <- C.lat.grid$mid*flow.lat.grid$mid

#==============================================================================
# FUNCTION BODY: CALCULATIONS 
#==============================================================================

N <- length(C)

# Calculate the discharge in each box (using the water balance)

for (i in 2:(N+1))
{
	flow.grid$int[i] = flow.grid$int[i-1] + flow.lat.grid$mid[i-1]
}

# Calculate diffusive part of the mass flow F 

Dif.F <- as.vector(-Disp.grid$int*diff(c(C.up,C,C.down)))  
Adv.F <- vector(length=(N+1))

for (i in 1:(N+1))
{
	if (flow.grid$int[i] > 0)
# advection directed downstream
	{
		Adv.F[i] <- flow.grid$int[i]*c(C.up,C)[i]
	}
	
	if (flow.grid$int[i] < 0)
# advection directed downstream
	{
		Adv.F[i] <- flow.grid$int[i]*c(C,C.down)[i+1]
	}
}

F <- as.vector(Dif.F + Adv.F) 

# Impose boundary fluxes when needed
if (! is.na (F.up))   F[1]   <- F.up
if (! is.na (F.down)) F[length(F)] <- F.down
    
# Calculate rate of change = Flux gradient + lateral input

dC <- -diff(F)/V.grid$mid + F.lat.grid$mid/V.grid$mid 

return (list (dC = dC,                      # Rate of change in the centre of each grid cells 
              F = F,                     # Flux across at the interface of each grid cell 
              F.up = F[1],               # Flux across upstream boundary ; positive = IN 
              F.down = F[length(F)],  # Flux across downstream boundary ; positive = OUT 
					    F.lat = F.lat.grid$mid))      # Flux lateral ; positive = IN

} # end tran1D.river

###############################################################################
# 2D Transport of a solute in the sediment
###############################################################################

tran.2D <- function(C,                   # Concentration in solid phase of porous medium
                    C.up=C[1,],          # Concentration at the upper interface
                    C.down=C[nrow(C),],  # Concentration at lower interface
                    C.left=C[,1],        # Concentration at the upper interface
                    C.right=C[,ncol(C)], # Concentration at lower interface
                    flux.up=NA,          # Flux at the upper interface (sediment-water interface)
                    flux.down=NA,        # Flux at the lower interface (deep sediment interface)
                    flux.left=NA,        # Flux at the upper interface (sediment-water interface)
                    flux.right=NA,       # Flux at the lower interface (deep sediment interface)
									  a.bl.up=NULL,        # transfer coefficient across the upper boundary layer; Flux = a.bl.up*(C.bl.up-C[1,])
									  C.bl.up=NULL,        # concentration at the top of the upper boundary layer 
									  a.bl.down=NULL,      # transfer coefficient across the lower boundary layer; Flux = a.bl.down*(C[N,]-C.bl.down)
									  C.bl.down=NULL,      # concentration at the top of the lower boundary layer 
									  a.bl.left=NULL,      # transfer coefficient across the left boundary layer; Flux = a.bl.left*(C.bl.left-C[,1])
									  C.bl.left=NULL,      # concentration at the top of the left boundary layer 
									  a.bl.right=NULL,     # transfer coefficient across the right boundary layer; Flux = a.bl.right*(C[,M]-C.bl.right)
									  C.bl.right=NULL,     # concentration at the top of the right boundary layer 
                    D.x=NULL,            # Diffusion coefficient (positive), defined on the grid cell interfaces, one value or a vector of length N+1 
                    D.y=D.x,             # Diffusion coefficient (positive), defined on the grid cell interfaces, one value or a vector of length N+1 
									  D.grid=NULL,         # Diffusion coefficient D defined as grid property
									  v.x=0,               # Advective velocity (positive or negative), defined on the grid cell interfaces, one value or a vector of length N+1
									  v.y=v.x,             # Advective velocity (positive or negative), defined on the grid cell interfaces, one value or a vector of length N+1
									  v.grid=NULL,         # Advective velocity v defined as grid property
									  AFDW.x=1,            # Advective Finite Difference Weight, weighing coefficient in finite difference scheme for advection (1:backward, 0:forward, 0.5:centred) 
									  AFDW.y=AFDW.x,       # Advective Finite Difference Weight, weighing coefficient in finite difference scheme for advection (1:backward, 0:forward, 0.5:centred) 
									  AFDW.grid=NULL,      # Weighing coefficient in finite difference scheme for advection (1:backward, 0:forward, 0.5:centred) 
                    VF.x=1,              # Volume fraction at the interface of the grid cells, one value 
                    VF.y=VF.x,           # Volume fraction at the interface of the grid cells, one value 
                    VF.grid=NULL,        # Volume fraction defined as a grid property
									  dx=NULL,             # thickness of grid cells 
									  dy=NULL,             # thickness of grid cells 
                    grid=NULL,           # grid definition
									  full.check = FALSE,  # logical flag enabling a full check of the consistency of the arguments (FALSE speeds up execution by 50%)
									  full.output = FALSE) # logical flag enabling a full return of the output (FALSE speeds up execution by 20%)
											
{

#==============================================================================
# DEFAULT INFILLING OF GRID PARAMETERS 
#==============================================================================

if (is.null(grid)) grid <- list(dx=rep(dx,length.out=nrow(C)),dx.aux=0.5*(c(0,rep(dx,length.out=nrow(C)))+c(rep(dx,length.out=nrow(C)),0)),dy=rep(dy,length.out=ncol(C)),dy.aux=0.5*(c(0,rep(dy,length.out=ncol(C)))+c(rep(dy,length.out=ncol(C)),0)))
if (is.null(AFDW.grid)) AFDW.grid <- list(x.int=matrix(data=AFDW.x,nrow=(nrow(C)+1),ncol=ncol(C)),y.int=matrix(data=AFDW.y,nrow=nrow(C),ncol=(ncol(C)+1)))
if (is.null(D.grid)) D.grid <- list(x.int=matrix(data=D.x,nrow=(nrow(C)+1),ncol=ncol(C)),y.int=matrix(data=D.y,nrow=nrow(C),ncol=(ncol(C)+1)))
if (is.null(v.grid)) v.grid <- list(x.int=matrix(data=v.x,nrow=(nrow(C)+1),ncol=ncol(C)),y.int=matrix(data=v.y,nrow=nrow(C),ncol=(ncol(C))+1))
if (is.null(VF.grid)) VF.grid <- list(x.int=matrix(data=VF.x,nrow=(nrow(C)+1),ncol=ncol(C)),y.int=matrix(data=VF.y,nrow=nrow(C),ncol=(ncol(C)+1)),x.mid=matrix(data=VF.x,nrow=nrow(C),ncol=ncol(C)),y.mid=matrix(data=VF.y,nrow=nrow(C),ncol=ncol(C)))

#==============================================================================
# INPUT CHECKS  
#==============================================================================


if (full.check)
{

# check dimensions of input concentrations

if (!is.null(C.up))
{
 if (!((length(C.up)==1) || (length(C.up)==(ncol(C)))))
  stop("error: C.up should be a vector of length 1 or ncol(C)")
} 
if (!is.null(C.down))
{
 if (!((length(C.down)==1) || (length(C.down)==(ncol(C)))))
  stop("error: C.down should be a vector of length 1 or ncol(C)")
} 
if (!is.null(C.left))
{
 if (!((length(C.left)==1) || (length(C.left)==(nrow(C)))))
  stop("error: C.left should be a vector of length 1 or nrow(C)")
} 
if (!is.null(C.right))
{
 if (!((length(C.right)==1) || (length(C.right)==(nrow(C)))))
  stop("error: C.right should be a vector of length 1 or nrow(C)")
} 

if (!is.null(C.bl.up))
{
 if (!((length(C.bl.up)==1) || (length(C.bl.up)==(ncol(C)))))
  stop("error: C.bl.up should be a vector of length 1 or ncol(C)")
} 

if (!is.null(C.bl.down))
{
 if (!((length(C.bl.down)==1) || (length(C.bl.down)==(ncol(C)))))
  stop("error: C.bl.down should be a vector of length 1 or ncol(C)")
} 

if (!is.null(C.bl.left))
{
 if (!((length(C.bl.left)==1) || (length(C.bl.left)==(nrow(C)))))
  stop("error: C.bl.left should be a vector of length 1 or nrow(C)")
} 

if (!is.null(C.bl.right))
{
 if (!((length(C.bl.right)==1) || (length(C.bl.right)==(nrow(C)))))
  stop("error: C.bl.right should be a vector of length 1 or nrow(C)")
} 

# check dimensions of input fluxes

if (!is.null(flux.up))
{
 if (!((length(flux.up)==1) || (length(flux.up)==(ncol(C)))))
  stop("error: flux.up should be a vector of length 1 or ncol(C)")
} 
if (!is.null(flux.down))
{
 if (!((length(flux.down)==1) || (length(flux.down)==(ncol(C)))))
  stop("error: flux.down should be a vector of length 1 or ncol(C)")
} 
if (!is.null(flux.left))
{
 if (!((length(flux.left)==1) || (length(flux.left)==(nrow(C)))))
  stop("error: flux.left should be a vector of length 1 or nrow(C)")
} 
if (!is.null(flux.right))
{
 if (!((length(flux.right)==1) || (length(flux.right)==(nrow(C)))))
  stop("error: flux.right should be a vector of length 1 or nrow(C)")
} 


# check input of grid

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

# check input of AFDW.grid

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

# check input of D.grid

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

# check input of v.grid

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

# check input of VF.grid

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
#==============================================================================
# FUNCTION BODY: CALCULATIONS 
#==============================================================================

N <- nrow(C)
M <- ncol(C)

# Impose boundary flux at upper boundary when needed
# Default boundary condition is no gradient
if (! is.na (flux.up[1]))
{
  nom <- flux.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1] + (1-AFDW.grid$x.int[1,])*v.grid$x.int[1,])*C[1,]
  denom <- VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1]+AFDW.grid$x.int[1,]*v.grid$x.int[1,])
  C.up <- nom/denom 
}

# Impose boundary flux at lower boundary when needed
# Default boundary condition is no gradient
if (! is.na (flux.down[1])) 
{
	nom <- flux.down - VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/grid$dx.aux[N+1] + AFDW.grid$x.int[(N+1),]*v.grid$x.int[(N+1),])*C[N,] 
  denom <- -VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/grid$dx.aux[N+1]+(1-AFDW.grid$x.int[(N+1),])*v.grid$x.int[(N+1),])
  C.down <- nom/denom 
}

# Impose boundary flux at upper boundary when needed
# Default boundary condition is no gradient
if (! is.na (flux.left[1]))
{
  nom <- flux.left + VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1] + (1-AFDW.grid$y.int[,1])*v.grid$y.int[,1])*C[,1]
  denom <- VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1]+AFDW.grid$y.int[,1]*v.grid$y.int[,1])
  C.left <- nom/denom 
}

# Impose boundary flux at lower boundary when needed
# Default boundary condition is no gradient
if (! is.na (flux.right[1])) 
{
	nom <- flux.right - VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/grid$dy.aux[M+1] + AFDW.grid$y.int[,(M+1)]*v.grid$y.int[,(M+1)])*C[,M] 
  denom <- -VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/grid$dy.aux[M+1]+(1-AFDW.grid$y.int[,(M+1)])*v.grid$y.int[,(M+1)])
  C.right <- nom/denom 
}

# when upper boundary layer is present, calculate new C.up 
if (!is.null(a.bl.up) & !is.null(C.bl.up[1]))
{
	nom <- a.bl.up*C.bl.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1] + (1-AFDW.grid$x.int[1,])*v.grid$x.int[1,])*C[1,]
  denom <- a.bl.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1]+ AFDW.grid$x.int[1,]*v.grid$x.int[1,])
	C.up <- nom/denom 
}

# when lower boundary layer is present, calculate new C.down 
if (!is.null(a.bl.down) & !is.null(C.bl.down[1]))
{
	nom <- a.bl.down*C.bl.down + VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/grid$dx.aux[(N+1)] + (1-AFDW.grid$x.int[(N+1),])*v.grid$x.int[(N+1),])*C[N,]
  denom <- a.bl.down + VF.grid$x.int[(N+1),]*(D.grid$x.int[(N+1),]/grid$dx.aux[(N+1)]+ AFDW.grid$x.int[(N+1),]*v.grid$x.int[(N+1),])
	C.down <- nom/denom 
}

# when left boundary layer is present, calculate new C.left 
if (!is.null(a.bl.left) & !is.null(C.bl.left[1]))
{
	nom <- a.bl.left*C.bl.left + VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1] + (1-AFDW.grid$y.int[,1])*v.grid$y.int[,1])*C[,1]
  denom <- a.bl.left + VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[1]+ AFDW.grid$y.int[,1]*v.grid$y.int[,1])
	C.left <- nom/denom 
}

# when right boundary layer is present, calculate new C.right
if (!is.null(a.bl.right) & !is.null(C.bl.right[1]))
{
	nom <- a.bl.right*C.bl.right + VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/grid$dy.aux[(M+1)] + (1-AFDW.grid$y.int[,(M+1)])*v.grid$y.int[,(M+1)])*C[,M]
  denom <- a.bl.right + VF.grid$y.int[,(M+1)]*(D.grid$y.int[,(M+1)]/grid$dy.aux[(M+1)]+ AFDW.grid$y.int[,(M+1)]*v.grid$y.int[,(M+1)])
	C.right <- nom/denom 
}

# Calculate diffusive part of the flux 
x.Dif.flux <- as.matrix(-VF.grid$x.int*D.grid$x.int*diff(rbind(C.up,C,C.down,deparse.level = 0))/matrix(data=grid$dx.aux,nrow=(N+1),ncol=M,byrow=FALSE))  
y.Dif.flux <- as.matrix(-VF.grid$y.int*D.grid$y.int*t(diff(t(cbind(C.left,C,C.right,deparse.level = 0))))/matrix(data=grid$dy.aux,nrow=N,ncol=(M+1),byrow=TRUE))  

# Calculate adcevtive part of the flux 
x.Adv.flux <- as.matrix(VF.grid$x.int*v.grid$x.int*((1-AFDW.grid$x.int)*rbind(C.up,C,deparse.level = 0)+AFDW.grid$x.int*rbind(C,C.down,deparse.level = 0)))  
y.Adv.flux <- as.matrix(VF.grid$y.int*v.grid$y.int*((1-AFDW.grid$y.int)*cbind(C.left,C,deparse.level = 0)+AFDW.grid$y.int*cbind(C,C.right,deparse.level = 0)))  

x.flux <- x.Dif.flux + x.Adv.flux 
y.flux <- y.Dif.flux + y.Adv.flux 

# Impose boundary fluxes when needed
# Default boundary condition is no gradient
if (! is.na (flux.up[1]))   x.flux[1,]   <- flux.up
if (! is.na (flux.down[1])) x.flux[nrow(x.flux),] <- flux.down
    
if (! is.na (flux.left[1]))  y.flux[,1]   <- flux.left
if (! is.na (flux.right[1])) y.flux[,ncol(y.flux)] <- flux.right

# Calculate rate of change = flux gradient
dFdx <- -(diff(x.flux)/grid$dx)/VF.grid$x.mid
dFdy <- -t(diff(t(y.flux))/grid$dy)/VF.grid$y.mid

if (!full.output){
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

