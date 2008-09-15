\name{tran1D}
\alias{tran1D}
\title{General 1-D transport}
\description{Estimates the rate of change of a substance due to diffusive and advective 1-D transport in a shape. 
  \cr The shape is defined by the surface areas, which do not need to be constant.
  \cr The default considers a constant surface area equal to 1.
  \cr The axis points from upstream (box 1) to downstream (last box)}
\usage{tran1D(y, y.up=y[1], y.down=y[length(y)], 
flux.up=NA, flux.down=NA, disp, v=0, v.up=0, weight.up=1, 
surf.mid=1,surf.int=surf.mid, 
dx=NULL, dx.int=dx, grid=list(dx=dx,dx.int=dx.int))}
\arguments{
  \item{y }{concentration (or other variable), defined in centre of boxes. A vector of length N, [Mass/Length3]}
  \item{y.up }{concentration at upstream boundary interface. One value, [M/L3]}
  \item{y.down }{concentration at downstream boundary interface. One value, [M/L3]}
  \item{flux.up }{flux across the upper boundary interface, positive = IN medium. One value, [M/L2/Time]}
  \item{flux.down }{flux across the lower boundary interface, positive= OUT of medium. One value, [M/L2/T]}
  \item{disp }{dispersion, diffusion coefficients, defined on box interfaces. One value or a vector of length N+1, [L2/T]}
  \item{v }{velocity in the axis direction, defined on box interfaces. One value or a vector of length N+1, [L/T]}
  \item{v.up }{velocity, against the axis direction, defined on box interfaces. One value or a vector of length N+1, [L/T]}
  \item{weight.up }{upstream weighing coefficient for velocity and upward velocity, defined on box interfaces; default is backward differences. One value or a vector of length N+1, [-]}
  \item{surf.mid }{surface area, defined in centre of boxes. One value or a vector of length N, [L2]}
  \item{surf.int }{surface area, defined at box interfaces. One value or a vector of length N+1, [L2]}
  \item{dx }{thickness of boxes. One value or a vector of length N, [L]}
  \item{dx.int }{distance over the interfaces, i.e. from centre to centre of boxes. One value or a vector of length N+1, [L]}
  \item{grid }{discretisation grid, a list containing at least elements dx (box thicknesses) and dx.int (interface distances), e.g. as calculated by setup.grid, [L]
    If present, overrules dx and dx.int}
}
\value{
  a list containing:
  \item{dy }{rate of change of y in each layer due to transport, [M/L3 liquid/T]}
  \item{flux }{Fluxes across each box interface. A vector of length N+1, [M/L2/Time]}
  \item{flux.up }{Flux across the upper boundary interface, positive = IN medium. One value, [M/L2/Time]}
  \item{flux.down }{Flux across the lower boundary interface, positive= OUT of medium. One value, [M/L2/T]}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{##################################################################
###### EXAMPLE: O2 in a cylindrical and spherical organism  ######
##################################################################

#====================#
# Model equations    #
#====================#

# the numerical model - rate of change=transport-consumption
CylinderMod <- function(time,O2,pars)
      return (list(tran1D(O2,y.down=BW,disp=Da,surf.mid=A,
      surf.int=Ai,dx=dx,dx.int=dxi)$dy-Q))

SphereMod <- function(time,O2,pars)
      return (list(tran1D(O2,y.down=BW,disp=Da,surf.mid=As,
      surf.int=Asi,dx=dx,dx.int=dxi)$dy-Q))

#====================#
# Model application  #
#====================#

# parameters

BW     <- 2         # mmol/m3,  oxygen conc in surrounding water
Da     <- 0.5       # cm2/d     effective diffusion coeff in organism
R      <- 0.0025    # cm        radius of organism
Q      <- 250000    # nM/cm3/d  oxygen consumption rate/ volume / day
L      <- 0.05      # cm        length of organism (if a cylinder)

# the analytical cylindrical and spherical model
cylinder <- function(Da,Q,BW,R,r)  BW+Q/(4*Da)*(r^2-R^2)
sphere   <- function(Da,Q,BW,R,r)  BW+Q/(6*Da)*(r^2-R^2)

# the numerical model: parameters
N  <- 40                           # layers in the body
dx <- R/N                          # thickness of each layer
x  <- seq(dx/2,by=dx,length.out=N) # distance of center to mid-layer
xi <- seq(0,by=dx,length.out=N+1)  # distance to layer interface
dxi<- c(rep(dx,N),dx/2)            # dispersion distances

# Cylindrical surfaces
A   <- 2*pi*x *L               # surface at mid-layer depth
Ai  <- 2*pi*xi*L               # surface at layer interface

# Spherical surfaces
As   <- 4*pi*x^2               # surface of sphere, at each mid-layer
Asi  <- 4*pi*xi^2              # surface at layer interface


# solve the model numerically for a cylinder
CONC  <- steady.1D (runif(N),fun=CylinderMod,nspec=1,atol=1e-10)
O2    <- CONC$y

# solve the model numerically for a sphere
CONC2 <- steady.1D (runif(N),fun=SphereMod,nspec=1,atol=1e-10)
O2b   <- CONC2$y

#====================#
# Plotting output    #
#====================#

plot(x,O2,xlab="distance from centre, cm",ylab="mmol/m3",
main="tran1D",sub="diffusion-reaction in a cylinder and sphere")

lines(x, cylinder(Da,Q,BW,R,x))

points(x, O2b, pch=18,col="red" )
lines(x, sphere(Da,Q,BW,R,x),col="red")

legend ("topleft",lty=c(1,NA),pch=c(NA,1),
        c("analytical solution","numerical approximation"))
legend ("bottomright",pch=c(1,18),lty=1,col=c("black","red"),
        c("cylinder","sphere"))}
\references{Soetaert and Herman, a guide to ecological modelling, Book in preparation}
\details{The \bold {boundary conditions} are either
  \item (1) zero-gradient (default)
  \item (2) concentration boundary 
  \item (3) flux boundary 
  
  if the flux boundary is specified, it overrules the other specifications. 
  
  \bold{Transport properties}
  The \emph{dispersion coefficient} (disp), \emph{velocity} (v) and \emph{upstream velocity}, v.up can be either one value or a vector.
  If they are a vector,they must be of length N+1, defined at all box interfaces, including upper and lower boundary.


  The \bold{surface area} needs to be defined:
  \item (1) at the centre of boxes (surf.mid), either one value or a vector of length N
  \item (2) at the box interfaces (surf.int), either one value, or a vector of length N+1

  If surf.int and surf.mid are one value, they should be equal.
  \cr No attempt is done to test whether the surfaces are compatible. 
  
  The \bold{spatial discretisation} is specified either:
  \item by a grid, as generated by \code{\link{setup.grid}} or 
  \item separate values for the thicknesses: dx.int, the "dispersion distance", i.e. the distance from the centre of one box to the centre of the box below. Either one value or a vector of length N+1. 
  Take care at boundaries!, and dx, the thickness of each box, either one value or a vector of length N.

  If dx.int and dx are one value, they should be equal
  \cr No attempt is done to test whether the discretisation is compatible.}
\note{The advection equation is not conservative. 
  
  This is ok if the advection is for instance a sinking or swimming velocity, then mass conservation is not an issue.
  
  However, if v is liquid flow, and it is not constant, then the rate of change should be adapted to conserve mass. 
  There should be inwelling if flux increases, outwelling if it decreases. This extra flux should be added to the rate of change.
  
  (For the same reason you should be careful when specifying advection rates for bodies whose surface is not constant).}
\keyword{utilities}

