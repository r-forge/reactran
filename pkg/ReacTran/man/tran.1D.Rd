\name{tran.1D}
\alias{tran.1D}
\title{General one-dimensional advective-diffusive transport}
\description{Estimates the transport term (i.e. the rate of change of the concentration due to diffusion and advection) in a one-dimensional model where the interfaces between grid cells can have a variable cross-sectional area.} 
\usage{tran.1D(C, C.up=C[1],C.down=C[length(C)], 
flux.up=NA, flux.down=NA, a.bl.up=NULL, C.bl.up=NULL, a.bl.down=NULL, C.bl.down=NULL,
D=NULL, D.grid=list(int=D), v=0, v.grid=list(int=v), AFDW=1, AFDW.grid=list(int=AFDW), 
VF=1, VF.grid=list(int=rep(VF,length.out=(length(C)+1)),mid=0.5*(rep(VF,length.out=(length(C)+1))[1:(length(C))]+rep(A,length.out=(length(C)+1))[2:(length(C)+1)])),
A=1, A.grid=list(int=rep(A,length.out=(length(C)+1)),mid=0.5*(rep(A,length.out=(length(C)+1))[1:(length(C))]+rep(A,length.out=(length(C)+1))[2:(length(C)+1)])),
dx=NULL,grid=list(dx=rep(dx,length.out=length(C)),dx.aux=0.5*(c(0,rep(dx,length.out=length(C)))+c(rep(dx,length.out=length(C)),0))),
full.check = FALSE, full.output = FALSE)}

\arguments{
  \item{C }{concentration, expressed per unit volume, defined at the centre of each grid cell. A vector of length N [M/L3] }
  \item{C.up }{concentration at upstream boundary. One value [M/L3]}
  \item{C.down }{concentration at downstream boundary. One value [M/L3]}
  \item{flux.up }{flux across the upstream boundary, positive = INTO model domain. One value [M/L2/T]}
  \item{flux.down }{flux across the downstream boundary, positive = OUT of model domain. One value [M/L2/T]}
  \item{a.bl.up }{transfer coefficient across the upstream boundary layer. Flux = a.bl.up*(C.bl.up-C[1]). One value [L/T]}
  \item{C.bl.up }{concentration at the top of the upstream boundary layer. One value [M/L3]}
  \item{a.bl.down }{transfer coefficient across the downstream boundary layer; Flux = a.bl.down*(C[N]-C.bl.down). One value [L/T]}
  \item{C.bl.down }{concentration at the top of the downstream boundary layer. One value [M/L3]}
  \item{D }{diffusion coefficient, defined on grid cell interfaces. One value or a vector of length N+1 [L2/T]}
  \item{D.grid }{diffusion coefficient packaged as a grid list; the list contains at least the element 'int' (see 'setup.1D.prop')  [L2/T]}
  \item{v }{advective velocity in the x-axis direction, defined on grid cell interfaces. Can be positive (downstream flow) or negative (upstream flow). One value or a vector of length N+1 [L/T]}
  \item{v.grid }{advective velocity packaged as a grid list; the list contains at least the element 'int' (see 'setup.1D.prop') [L/T]}
  \item{AFDW }{weight used in the finite difference scheme for advection, defined on grid cell interfaces; backward = 1, centred = 0.5, forward = 0; default is backward. One value or a vector of length N+1 [-]}
  \item{AFDW.grid }{AFDW coefficients packaged as a grid list; the list contains at least the element 'int' (see 'setup.1D.prop') [-]}
  \item{VF }{Volume fraction defined at the grid cell interfaces, one value or a vector of length N+1 [L2]}
  \item{VF.grid }{Volume fraction packaged as a grid list; the list contains at least the elements 'int' and 'mid' (see 'setup.1D.prop') [L2]}
  \item{A }{Interface area defined at the grid cell interfaces, one value or a vector of length N+1 [L2]}
  \item{A.grid }{Interface area packaged as a grid list; the list contains at least the elements 'int' and 'mid' (see 'setup.1D.prop') [L2]}
  \item{dx }{distance between adjacent cell interfaces (thickness of grid cells). One value or vector of length N [L]}
  \item{grid }{discretization grid, a list containing at least elements 'dx' and 'dx.aux' (see 'setup.1D.grid') [L]}
  \item{full.check }{logical flag enabling a full check of the consistency of the arguments (default = FALSE, TRUE slows down execution by 50 percent)}
  \item{full.output }{logical flag enabling a full return of the output (default = FALSE, TRUE slows down execution by 20 percent)}
}
\value{
  when (full.output = FALSE) then a list containing:
  \item{dC }{The rate of change of the concentration C due to transport, defined in the centre of each grid cell [M/L3/T]}
  when (full.output = TRUE) then a list containing:
  \item{dC }{The rate of change of the concentration C due to transport, defined in the centre of each grid cell [M/L3/T]}
  \item{C.up }{Concentration at the upstream interface. One value [M/L3]}
  \item{C.down }{Concentration at the downstream interface. One value [M/L3]}
  \item{adv.flux }{Advective flux across at the interface of each grid cell. A vector of length N+1 [M/L2/T]}
  \item{dif.flux }{Diffusive flux across at the interface of each grid cell. A vector of length N+1 [M/L2/T]}
  \item{flux }{Total flux across at the interface of each grid cell. A vector of length N+1 [M/L2/T]}
  \item{flux.up }{Flux across the upstream boundary, positive = INTo model domain. One value [M/L2/T]}
  \item{flux.down }{Flux across the downstream boundary, positive = OUT of model domain. One value [M/L2/T]}
}
\author{Filip Meysman <f.meysman@nioo.knaw.nl>, Karline Soetaert <k.soetaert@nioo.knaw.nl> }
\examples{
##################################################################
###### EXAMPLE 1: O2 and OC consumption in sediments        ######
##################################################################

# this example uses only the volume fractions 
# in the reactive transport term

#====================#
# Model formulation  #
#====================#

# Monod consumption of oxygen (O2)

O2.model <- function (t=0,O2,pars=NULL)
{
tran <- tran.1D(C=O2,C.up=C.ow.O2,D.grid=D.grid,v.grid=v.grid,VF.grid=por.grid,grid=grid)$dC 
reac <- - R.O2*(O2/(Ks+O2))
return(list(dCdt = tran+reac))
}

# First order consumption of organic carbon (OC)

OC.model <- function (t=0,OC,pars=NULL)
{
tran <- tran.1D(C=OC,flux.up=F.OC,D.grid=Db.grid,v.grid=v.grid,VF.grid=svf.grid,grid=grid)$dC 
reac <- - k*OC
return(list(dCdt = tran + reac))
}

#======================#
# Parameter definition #
#======================#

# Parameter values

F.OC    <- 25              # input flux organic carbon [micromol cm-2 yr-1]
C.ow.O2 <- 0.25            # concentration O2 in overlying water [micromol cm-3]
por     <- 0.8             # porosity
D       <- 400             # diffusion coefficient O2 [cm2 yr-1]
Db      <- 10              # mixing coefficient sediment [cm2 yr-1]
v       <- 1               # advective velocity [cm yr-1]
k       <- 1              # decay constant organic carbon [yr-1]
R.O2    <- 10               # O2 consumption rate [micromol cm-3 yr-1]
Ks      <- 0.005           # O2 consumption saturation constant 

# Grid definition

L <- 10   # depth of sediment domain [cm]
N <- 100  # number of grid layers
grid <- setup.1D.grid(x.up=0,L=L,N=N)

# Volume fractions 

por.grid <- setup.1D.prop(value=por,grid=grid)
svf.grid <- setup.1D.prop(value=(1-por),grid=grid)
D.grid <- setup.1D.prop(value=D,grid=grid)
Db.grid <- setup.1D.prop(value=Db,grid=grid)
v.grid <- setup.1D.prop(value=v,grid=grid)

#====================#
# Model solution     #
#====================#

# Initial conditions + simulation O2

O2 <- rep(0,length.out=N) 
O2 <- steady.band(y=O2, func=O2.model, nspec=1)$y

# Initial conditions + simulation OC

OC <- rep(0,length.out=N) 
OC <- steady.band(y=OC, func=OC.model, nspec=1)$y

# Plotting output

par(mfrow=c(1,2))

matplot(O2,grid$x.mid,pch=16,type="b",ylim=c(L,0), xlim=c(min(0,min(O2)),max(O2)),
xlab="",ylab="depth [cm]",main=expression("O2 concentration"),axes=FALSE)
abline(h = 0)
axis(pos=0, side=2)
axis(pos=0, side=3)

matplot(OC,grid$x.mid,pch=16,type="b",ylim=c(L,0), xlim=c(min(0,min(OC)),max(OC)),
xlab="",ylab="depth [cm]",main=expression("OC concentration"),axes=FALSE)
abline(h = 0)
axis(pos=0, side=2)
axis(pos=0, side=3)

##################################################################
###### EXAMPLE 2: O2 in a cylindrical and spherical organism  ######
##################################################################

# this example uses only the surface areas 
# in the reactive transport term

#====================#
# Model formulation  #
#====================#

# the numerical model - rate of change=transport-consumption
Cylinder.Model <- function(time,O2,pars)
      return (list(tran.1D(C=O2,C.down=BW,D=Da,A=A.cyl,dx=dx)$dC-Q))

Sphere.Model <- function(time,O2,pars)
      return (list(tran.1D(C=O2,C.down=BW,D=Da,A=A.sphere,dx=dx)$dC-Q))

#======================#
# Parameter definition #
#======================#

# parameter values

BW     <- 2         # mmol/m3,  oxygen conc in surrounding water
Da     <- 0.5       # cm2/d     effective diffusion coeff in organism
R      <- 0.0025    # cm        radius of organism
Q      <- 250000    # nM/cm3/d  oxygen consumption rate/ volume / day
L      <- 0.05      # cm        length of organism (if a cylinder)

# the numerical model

N  <- 40                           # layers in the body
dx <- R/N                          # thickness of each layer
x.mid <- seq(dx/2,by=dx,length.out=N) # distance of center to mid-layer
x.int <- seq(0,by=dx,length.out=N+1)  # distance to layer interface

# Cylindrical surfaces
A.cyl   <- 2*pi*x.int*L               # surface at mid-layer depth

# Spherical surfaces
A.sphere <- 4*pi*x.int^2               # surface of sphere, at each mid-layer

#====================#
# Model solution     #
#====================#

# the analytical solution of cylindrical and spherical model
cylinder <- function(Da,Q,BW,R,r)  BW+Q/(4*Da)*(r^2-R^2)
sphere   <- function(Da,Q,BW,R,r)  BW+Q/(6*Da)*(r^2-R^2)

# solve the model numerically for a cylinder
O2.cyl <- steady.1D (runif(N),func=Cylinder.Model,nspec=1,atol=1e-10)$y

# solve the model numerically for a sphere
O2.sphere <- steady.1D (runif(N),func=Sphere.Model,nspec=1,atol=1e-10)$y

#====================#
# Plotting output    #
#====================#

windows()
par(mfrow=c(1,1))

plot(x.mid,O2.cyl,xlab="distance from centre, cm",ylab="mmol/m3",
main="tran.1D",sub="diffusion-reaction in a cylinder and sphere")
lines(x.mid, cylinder(Da,Q,BW,R,x.mid))

points(x.mid, O2.sphere, pch=18,col="red" )
lines(x.mid, sphere(Da,Q,BW,R,x.mid),col="red")

legend ("topleft",lty=c(1,NA),pch=c(NA,1),
        c("analytical solution","numerical approximation"))
legend ("bottomright",pch=c(1,18),lty=1,col=c("black","red"),
        c("cylinder","sphere"))

##################################################################
###### EXAMPLE 3: O2 consumption in a spherical aggregate   ######
##################################################################

# this example uses both the surface areas and the volume fractions
# in the reactive transport term

#====================#
# Model formulation  #
#====================#

Aggregate.Model <- function(time,O2,pars)
{
tran <- tran.1D(C=O2,C.down=C.ow.O2,D.grid=D.grid,A.grid=A.grid,VF.grid=por.grid,grid=grid)$dC 
reac <- - R.O2*(O2/(Ks+O2))*(O2>0)
return(list(dCdt = tran+reac))
}

#======================#
# Parameter definition #
#======================#

# Parameters

C.ow.O2 <- 0.25            # concentration O2 in overlying water [micromol cm-3]
por     <- 0.8             # porosity
D       <- 400             # diffusion coefficient O2 [cm2 yr-1]
v       <- 0               # advective velocity [cm yr-1]
R.O2    <- 1000000         # O2 consumption rate [micromol cm-3 yr-1]
Ks      <- 0.005           # O2 consumption saturation constant [micromol cm-3]

# Grid definition
R <- 0.025           # radius of the agggregate [cm]
N <- 100             # number of grid layers
grid <- setup.1D.grid(x.up=0,L=R,N=N)

# Volume fractions 

por.grid <- setup.1D.prop(value=por,grid=grid)
D.grid <- setup.1D.prop(value=D,grid=grid)

# Surfaces 

A.mid <- 4*pi*grid$x.mid^2               # surface of sphere, at each mid-layer
A.int <- 4*pi*grid$x.int^2               # surface of sphere, at each mid-layer
A.grid=list(int=A.int,mid=A.mid)

#====================#
# Model solution     #
#====================#

# Numerical solution: staedy state 

O2.agg <- steady.1D (runif(N),func=Aggregate.Model,nspec=1,atol=1e-10)$y

#====================#
# Plotting output    #
#====================#

windows()
par(mfrow=c(1,1))

plot(grid$x.mid,O2.agg,xlab="distance from centre, cm",ylab="mmol/m3",
main="Diffusion-reaction of O2 in a spherical aggregate")
legend ("bottomright",pch=c(1,18),lty=1,col=c("black"),
        c("O2 concentration"))

}
\references{Soetaert and Herman, a guide to ecological modelling - using R as
a simulation platform, 2008. Springer}
\details{The \bold {boundary conditions} are either
\itemize{
  \item (1) zero-gradient 
  \item (2) fixed concentration 
  \item (3) convective boundary layer 
  \item (4) fixed flux 
}
  
  The above order also shows the priority. The default condition is the zero gradient. The convective boundary layer condition overrules the fixed concentration. The fixed flux overrules all other specifications. 
  
  \bold{Transport properties}
  The \emph{diffusion coefficient} (D), the \emph{advective velocity} (v), the \emph{volume fraction} (VF), the \emph{interface surface} (A), and the \emph{advective finite difference weight} (AFDW)  can be either one value or a vector.
  When a vector, this vector must be of length N+1, defined at all grid cell interfaces, including upper and lower boundary.

  The \bold{finite difference grid } (grid) is specified either:
\itemize{
  \item as a grid list, as generated by \code{\link{setup.1D.grid}} or 
  \item by the parameter dx representing the thickness of the grid cells (one value or a vector of length N+1)
}
}  
\note{ 
The advective-diffusion equation is not checked for mass conservation. Sometimes, this is not an issue, for instance when v represents a sinking velocity of particles or a swimming velocity of organisms. In others cases however, mass conservation needs to be accounted for. To ensure mass conservation, the advective velocity must obey certain continuity constraints: in essence the product of the volume fraction (VF), interface surface (A) and advective velocity (v) should be constant. In sediments, one can use 'setup.1D.advection' to ensure that the advective velocities for the pore water and solid phase meet these constraints. 
}
\keyword{utilities}

