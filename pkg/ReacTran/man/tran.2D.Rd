\name{tran.2D}
\alias{tran.2D}
\title{
  General two-dimensional advective-diffusive transport
}

\description{
  Estimates the transport term (i.e. the rate of change of a concentration
  due to diffusion and advection) in a two-dimensional rectangular model domain.
}

\usage{
tran.2D ( C, C.up=C[1,], C.down=C[nrow(C),],
  C.left=C[,1], C.right=C[,ncol(C)],
  flux.up=NULL, flux.down=NULL, flux.left=NULL, flux.right=NULL,
  a.bl.up=NULL, C.bl.up=NULL, a.bl.down=NULL, C.bl.down=NULL,
  a.bl.left=NULL, C.bl.left=NULL, a.bl.right=NULL, C.bl.right=NULL,
  D.x=NULL, D.y=D.x, D.grid=NULL,
  v.x=0, v.y=v.x, v.grid=NULL,
  AFDW.x=1, AFDW.y=AFDW.x, AFDW.grid=NULL,
  VF.x=1, VF.y=VF.x, VF.grid=NULL,
  dx=NULL, dy=NULL, grid=NULL,
  full.check = FALSE, full.output = FALSE)
}

\arguments{
  \item{C }{concentration, expressed per unit volume, defined at the centre
    of each grid cell; NxM matrix [M/L3].
  }
  \item{C.up }{concentration at upstream boundary; vector of length M [M/L3].
  }
  \item{C.down }{concentration at downstream boundary; vector of length
    M [M/L3].
  }
  \item{C.left }{concentration at left boundary; vector of length N [M/L3].
  }
  \item{C.right }{concentration at right boundary; vector of length N [M/L3].
  }
  \item{flux.up }{flux across the upstream boundary, positive = INTO model
    domain; vector of length M [M/L2/T].
  }
  \item{flux.down }{flux across the downstream boundary, positive = OUT of
    model domain; vector of length M [M/L2/T].
  }
  \item{flux.left }{flux across the left boundary, positive = INTO model
    domain; vector of length N [M/L2/T].
  }
  \item{flux.right }{flux across the right boundary, positive = OUT of
    model domain; vector of length N [M/L2/T].
  }
  \item{a.bl.up }{transfer coefficient across the upstream boundary layer.
    \code{Flux=a.bl.up*(C.bl.up-C[1,])}. One value [L/T].
  }
  \item{C.bl.up }{concentration at the upstream boundary layer; vector of
    length M  [M/L3].
  }
  \item{a.bl.down }{transfer coefficient across the downstream boundary
    layer; \code{Flux=a.bl.down*(C[N,]-C.bl.down)}. One value [L/T].
  }
  \item{C.bl.down }{concentration at the downstream boundary layer; vector of
    length M  [M/L3].
  }
  \item{a.bl.left }{transfer coefficient across the left boundary layer.
    \code{Flux=a.bl.left*(C.bl.left-C[,1])}. One value [L/T].
  }
  \item{C.bl.left }{concentration at the left boundary layer; vector of
    length N  [M/L3].
  }
  \item{a.bl.right }{transfer coefficient across the right boundary layer;
    \code{Flux=a.bl.right*(C[,M]-C.bl.right)}. One value [L/T].
  }
  \item{C.bl.right }{concentration at the right boundary layer; vector of
    length N  [M/L3].
  }
  \item{D.x }{diffusion coefficient in x-direction, defined on grid cell
    interfaces. One value or an (N+1)x(M) matrix [L2/T].
  }
  \item{D.y }{diffusion coefficient in y-direction, defined on grid cell
    interfaces. One value or an (N)x(M+1) matrix [L2/T].
  }
  \item{D.grid }{diffusion coefficient packaged as a grid list; the list
    contains at least the elements \code{x.int} and \code{y.int}
    (see \code{\link{setup.prop.2D}})  [L2/T].
  }
  \item{v.x }{advective velocity in the x-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value or an (N+1)x(M) matrix  [L/T].
  }
  \item{v.y }{advective velocity in the y-direction, defined on grid cell
    interfaces. Can be positive (downstream flow) or negative (upstream flow).
    One value or an (N)x(M+1) matrix  [L/T].
  }
  \item{v.grid }{advective velocity packaged as a grid list; the list
    contains at least the elements \code{x.int} and \code{y.int} (see
    \code{\link{setup.prop.2D}}) [L/T].
  }
  \item{AFDW.x }{weight used in the finite difference scheme for advection
    in the x-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward. One value or an
    (N+1)x(M) matrix [-].
  }
  \item{AFDW.y }{weight used in the finite difference scheme for advection
    in the y-direction, defined on grid cell interfaces; backward = 1,
    centred = 0.5, forward = 0; default is backward. One value or an
    (N)x(M+1) matrix [-].
  }
  \item{AFDW.grid }{AFDW coefficients packaged as a grid list; the list
    contains at least the elements \code{x.int} and \code{y.int}
    (see \code{\link{setup.prop.2D}}) [-].
  }
  \item{VF.x }{Volume fraction at the grid cell interfaces in the x-direction.
    One value or an (N+1)x(M) matrix [-].
  }
  \item{VF.y }{Volume fraction at the grid cell interfaces in the y-direction.
    One value or an (N)x(M+1) matrix [-].
  }
  \item{VF.grid }{Volume fraction packaged as a grid list; the list contains
    at least the elements \code{x.int}, \code{y.int}, 'x.mid' and 'y.mid'
    (see \code{\link{setup.prop.2D}}) [-].
  }
  \item{dx }{distance between adjacent cell interfaces in the x-direction
    (thickness of grid cells). One value or vector of length N [L].
  }
  \item{dy }{distance between adjacent cell interfaces in the y-direction
    (thickness of grid cells). One value or vector of length M [L].
  }
  \item{grid }{discretization grid, a list containing at least elements
    \code{dx}, \code{dx.aux}, \code{dy}, \code{dy.aux}
    (see \code{\link{setup.grid.2D}}) [L].
  }
  \item{full.check }{logical flag enabling a full check of the consistency
    of the arguments (default = \code{FALSE}; \code{TRUE} slows down
    execution by 50 percent).
  }
  \item{full.output }{logical flag enabling a full return of the output
    (default = \code{FALSE}; \code{TRUE} slows down execution by 20 percent).
  }
}
\value{
  \item{when (\code{full.output = FALSE}) then a list containing:}
  \item{dC }{the rate of change of the concentration C due to transport,
    defined in the centre of each grid cell [M/L3/T].
  }
  \item{when (\code{full.output = TRUE}) then a list containing:
  }
  \item{dC }{the rate of change of the concentration C due to transport,
    defined in the centre of each grid cell [M/L3/T].
  }
  \item{C.up }{concentration at the upstream interface.
    A vector of length M [M/L3].
  }
  \item{C.down }{concentration at the downstream interface.
    A vector of length M [M/L3].
  }
  \item{C.left }{concentration at the left interface.
    A vector of length N [M/L3].
  }
  \item{C.right }{concentration at the right interface.
    A vector of length N [M/L3].
  }
  \item{x.flux }{flux across the horizontal interface of the grid cells.
    A (N+1)x(M) matrix  [M/L2/T].
  }
  \item{y.flux }{flux across the vertical interfaces of the grid cells.
    A (N)x(M+1) matrix [M/L2/T].
  }
  \item{flux.up }{flux across the upstream boundary, positive = INTO model
    domain. A vector of length M [M/L2/T].
  }
  \item{flux.down }{flux across the downstream boundary, positive = OUT
    of model domain. A vector of length M [M/L2/T].
  }
  \item{flux.left }{flux across the left boundary, positive = INTO model
    domain. A vector of length N [M/L2/T].
  }
  \item{flux.right }{flux across the right boundary, positive = OUT of model
    domain. A vector of length N [M/L2/T].
  }
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
# Parameters
F        <- 100             # input flux [micromol cm-2 yr-1]
por      <- 0.8             # constant porosity
D        <- 400             # mixing coefficient [cm2 yr-1]
v        <- 1               # advective velocity [cm yr-1]

# Grid definition
x.N <- 4   # number of cells in x-direction
y.N <- 6   # number of cells in y-direction
x.L <- 8   # domain size x-direction [cm]
y.L <- 24  # domain size y-direction [cm]
dx       <- x.L/x.N             # cell size x-direction [cm]
dy       <- y.L/y.N             # cell size y-direction [cm]
 
# Intial conditions 
C <- matrix(nrow=x.N,ncol=y.N,data=0,byrow=FALSE)

# Boundary conditions: fixed concentration  
C.up=rep(1,times=y.N)
C.down=rep(0,times=y.N)
C.left=rep(1,times=x.N)
C.right=rep(0,times=x.N)

# Only diffusion 
tran.2D(full.output=TRUE,C=C,D.x=D,D.y=D,v.x=0,
v.y=0,VF.x=por,VF.y=por,dx=dx,dy=dy,C.up=C.up,C.down=C.down,
C.left=C.left,C.right=C.right)  

# Strong advection, backward (default), central and forward 
#finite difference schemes 
tran.2D(C=C,D.x=D,v.x=100*v,VF.x=por,dx=dx,dy=dy,
C.up=C.up,C.down=C.down,C.left=C.left,C.right=C.right)  
tran.2D(AFDW.x=0.5,C=C,D.x=D,v.x=100*v,VF.x=por,dx=dx,dy=dy,
C.up=C.up,C.down=C.down,C.left=C.left,C.right=C.right)  
tran.2D(AFDW.x=0,C=C,D.x=D,v.x=100*v,VF.x=por,dx=dx,dy=dy,
C.up=C.up,C.down=C.down,C.left=C.left,C.right=C.right)  

# Boundary conditions: fixed fluxes 

flux.up=rep(200,times=y.N)
flux.down=rep(-200,times=y.N)
flux.left=rep(200,times=x.N)
flux.right=rep(-200,times=x.N)
tran.2D(C=C,D.x=D,v.x=0,VF.x=por,dx=dx,dy=dy,
flux.up=flux.up,flux.down=flux.down,
flux.left=flux.left,flux.right=flux.right)  

# Boundary conditions: convective boundary layer on all sides

a.bl <- 800   # transfer coefficient
C.bl.up <- rep(1,times=(y.N)) # fixed conc at boundary layer
C.bl.left <- rep(1,times=(x.N)) # fixed conc at boundary layer
tran.2D(full.output=TRUE,C=C,D.x=D,v.x=0,VF.x=por,
dx=dx,dy=dy,C.bl.up=C.bl.up,a.bl.up=a.bl,C.bl.down=C.bl.up,
a.bl.down=a.bl,C.bl.left=C.bl.left,a.bl.left=a.bl,
C.bl.right=C.bl.left,a.bl.right=a.bl)  

# Runtime test with and without argument checking

n.iterate <-1000

test1 <- function()
{
for (i in 1:n.iterate )
ST<-tran.2D(full.check=TRUE,C=C,D.x=D,v.x=0,VF.x=por,
dx=dx,dy=dy,C.bl.up=C.bl.up,a.bl.up=a.bl,C.down=C.down)  
} 
system.time(test1())

test2 <- function()
{
for (i in 1:n.iterate )
ST<-tran.2D(full.output=TRUE,C=C,D.x=D,v.x=0,VF.x=por,
dx=dx,dy=dy,C.bl.up=C.bl.up,a.bl.up=a.bl,C.down=C.down)  
} 
system.time(test2())

test3 <- function()
{
for (i in 1:n.iterate )
ST<-tran.2D(full.output=TRUE,full.check=TRUE,C=C,D.x=D,v.x=0,
VF.x=por,dx=dx,dy=dy,C.bl.up=C.bl.up,a.bl.up=a.bl,C.down=C.down)  
} 
system.time(test3())

}
\references{
  Soetaert and Herman, a practical guide to ecological modelling - using R as
  a simulation platform, 2009. Springer
}

\details{
  \itemize{
    \item {The x-axis is taken in the vertical pointing downwards (N grid cells).}
    \item {The y-axis is taken in the horizontal pointing to the right (M grid cells).}
  }
  The \bold{boundary conditions} are either
  \itemize{
    \item (1) zero-gradient
    \item (2) fixed concentration
    \item (3) convective boundary layer
    \item (4) fixed flux
  }
  This is also the order of priority. The zero gradient is the default,
  the fixed flux overrules all other.
}  

\keyword{utilities}

