\name{tran.2D}
\alias{tran.2D}
\title{General two-dimensional advective-diffusive transport}
\description{Estimates the transport term (i.e. the rate of change of the concentration due to diffusion and advection) in a two-dimensional rectangular model domain where the interfaces between grid cells can have a variable cross-sectional area.} 
\usage{tran.2D(C, C.up=C[1,], C.down=C[nrow(C),], C.left=C[,1], C.right=C[,ncol(C)], 
flux.up=NA, flux.down=NA, flux.left=NA, flux.right=NA,
a.bl.up=NULL, C.bl.up=NULL, a.bl.down=NULL, C.bl.down=NULL, a.bl.left=NULL, C.bl.left=NULL, a.bl.right=NULL, C.bl.right=NULL,
D.x=NULL, D.y=D.x, D.grid=list(x.int=matrix(data=D.x,nrow=(nrow(C)+1),ncol=ncol(C)),y.int=matrix(data=D.y,nrow=nrow(C),ncol=(ncol(C)+1))), 
v.x=0, v.y=v.x, v.grid=list(x.int=matrix(data=v.x,nrow=(nrow(C)+1),ncol=ncol(C)),y.int=matrix(data=v.y,nrow=nrow(C),ncol=(ncol(C))+1)), 
AFDW.x=1, AFDW.y=AFDW.x, AFDW.grid=list(x.int=matrix(data=AFDW.x,nrow=(nrow(C)+1),ncol=ncol(C)),y.int=matrix(data=AFDW.y,nrow=nrow(C),ncol=(ncol(C)+1))), 
VF.x=1, VF.y=VF.x, VF.grid=list(x.int=matrix(data=VF.x,nrow=(nrow(C)+1),ncol=ncol(C)),y.int=matrix(data=VF.y,nrow=nrow(C),ncol=(ncol(C)+1)),x.mid=matrix(data=VF.x,nrow=nrow(C),ncol=ncol(C)),y.mid=matrix(data=VF.y,nrow=nrow(C),ncol=ncol(C))),
dx=NULL, dy=NULL, grid=list(dx=rep(dx,length.out=nrow(C)),dx.aux=0.5*(c(0,rep(dx,length.out=nrow(C)))+c(rep(dx,length.out=nrow(C)),0)),dy=rep(dy,length.out=ncol(C)),dy.aux=0.5*(c(0,rep(dy,length.out=ncol(C)))+c(rep(dy,length.out=ncol(C)),0))),
full.check = FALSE, full.output = FALSE)}

\arguments{
  \item{C }{Concentration, expressed per unit volume, defined at the centre of each grid cell. A NxM matrix [M/L3] }
  \item{C.up }{Concentration at upstream boundary. A vector of length M [M/L3]}
  \item{C.down }{Concentration at downstream boundary. A vector of length M [M/L3]}
  \item{C.left }{Concentration at left boundary. A vector of length N [M/L3]}
  \item{C.right }{Concentration at right boundary. A vector of length N [M/L3]}
  \item{flux.up }{Flux across the upstream boundary, positive = INTO model domain. A vector of length M [M/L2/T]}
  \item{flux.down }{Flux across the downstream boundary, positive = OUT of model domain. A vector of length M [M/L2/T]}
  \item{flux.left }{Flux across the left boundary, positive = INTO model domain. A vector of length N [M/L2/T]}
  \item{flux.right }{Flux across the right boundary, positive = OUT of model domain. A vector of length N [M/L2/T]}
  \item{a.bl.up }{Transfer coefficient across the upstream boundary layer. Flux = a.bl.up*(C.bl.up-C[1,]). One value [L/T]}
  \item{C.bl.up }{Concentration at the top of the upstream boundary layer. A vector of length M  [M/L3]}
  \item{a.bl.down }{Transfer coefficient across the downstream boundary layer; Flux = a.bl.down*(C[N,]-C.bl.down). One value [L/T]}
  \item{C.bl.down }{Concentration at the top of the downstream boundary layer. A vector of length M  [M/L3]}
  \item{a.bl.left }{Transfer coefficient across the left boundary layer. Flux = a.bl.left*(C.bl.left-C[,1]). One value [L/T]}
  \item{C.bl.left }{Concentration at the top of the left boundary layer. A vector of length N  [M/L3]}
  \item{a.bl.right }{Transfer coefficient across the right boundary layer; Flux = a.bl.right*(C[,M]-C.bl.right). One value [L/T]}
  \item{C.bl.right }{Concentration at the top of the right boundary layer. A vector of length N  [M/L3]}
  \item{D.x }{diffusion coefficient in x-direction, defined on grid cell interfaces. One value or an (N+1)x(M) matrix [L2/T]}
  \item{D.y }{diffusion coefficient in y-direction, defined on grid cell interfaces. One value or an (N)x(M+1) matrix [L2/T]}
  \item{D.grid }{diffusion coefficient packaged as a grid list; the list contains at least the elements 'x.int' and 'y.int' (see 'setup.2D.prop')  [L2/T]}
  \item{v.x }{advective velocity in the x-direction, defined on grid cell interfaces. Can be positive (downstream flow) or negative (upstream flow). One value or an (N+1)x(M) matrix  [L/T]}
  \item{v.y }{advective velocity in the y-direction, defined on grid cell interfaces. Can be positive (downstream flow) or negative (upstream flow). One value or an (N)x(M+1) matrix  [L/T]}
  \item{v.grid }{advective velocity packaged as a grid list; the list contains at least the elements 'x.int' and 'y.int' (see 'setup.2D.prop') [L/T]}
  \item{AFDW.x }{weight used in the finite difference scheme for advection in the x-direction, defined on grid cell interfaces; backward = 1, centred = 0.5, forward = 0; default is backward. One value or an (N+1)x(M) matrix [-]}
  \item{AFDW.y }{weight used in the finite difference scheme for advection in the y-direction, defined on grid cell interfaces; backward = 1, centred = 0.5, forward = 0; default is backward. One value or an (N)x(M+1) matrix [-]}
  \item{AFDW.grid }{AFDW coefficients packaged as a grid list; the list contains at least the elements 'x.int' and 'y.int' (see 'setup.2D.prop') [-]}
  \item{VF.x }{Volume fraction at the grid cell interfaces in the x-direction. One value or an (N+1)x(M) matrix [L2]}
  \item{VF.y }{Volume fraction at the grid cell interfaces in the y-direction. One value or an (N)x(M+1) matrix [L2]}
  \item{VF.grid }{Volume fraction packaged as a grid list; the list contains at least the elements 'x.int', 'y.int', 'x.mid' and 'y.mid' (see 'setup.1D.prop') [L2]}
  \item{dx }{distance between adjacent cell interfaces in the x-direction (thickness of grid cells). One value or vector of length N [L]}
  \item{dy }{distance between adjacent cell interfaces in the y-direction (thickness of grid cells). One value or vector of length M [L]}
  \item{grid }{discretization grid, a list containing at least elements 'dx', 'dx.aux', 'dy', 'dy.aux' (see 'setup.2D.grid') [L]}
  \item{full.check }{logical flag enabling a full check of the consistency of the arguments (default = TRUE, FALSE speeds up execution by 50 percent)}
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
  \item{flux.up }{Flux across the upstream boundary, positive = INTO model domain. One value [M/L2/T]}
  \item{flux.down }{Flux across the downstream boundary, positive = OUT of model domain. One value [M/L2/T]}
}
\author{Filip Meysman <f.meysman@nioo.knaw.nl>, Karline Soetaert <k.soetaert@nioo.knaw.nl> }
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
tran.2D(full.output=TRUE,C=C,D.x=D,D.y=D,v.x=0,v.y=0,VF.x=por,VF.y=por,dx=dx,dy=dy,C.up=C.up,C.down=C.down,C.left=C.left,C.right=C.right)  

# Strong advection, backward (default), central and forward finite difference schemes 
tran.2D(C=C,D.x=D,v.x=100*v,VF.x=por,dx=dx,dy=dy,C.up=C.up,C.down=C.down,C.left=C.left,C.right=C.right)  
tran.2D(AFDW.x=0.5,C=C,D.x=D,v.x=100*v,VF.x=por,dx=dx,dy=dy,C.up=C.up,C.down=C.down,C.left=C.left,C.right=C.right)  
tran.2D(AFDW.x=0,C=C,D.x=D,v.x=100*v,VF.x=por,dx=dx,dy=dy,C.up=C.up,C.down=C.down,C.left=C.left,C.right=C.right)  

# Boundary conditions: fixed fluxes 

flux.up=rep(200,times=y.N)
flux.down=rep(-200,times=y.N)
flux.left=rep(200,times=x.N)
flux.right=rep(-200,times=x.N)
tran.2D(C=C,D.x=D,v.x=0,VF.x=por,dx=dx,dy=dy,flux.up=flux.up,flux.down=flux.down,flux.left=flux.left,flux.right=flux.right)  

# Boundary conditions: convective boundary layer on all sides

a.bl <- 800   # transfer coefficient
C.bl.up <- rep(1,times=(y.N)) # fixed concentration on top of boundary layer
C.bl.left <- rep(1,times=(x.N)) # fixed concentration on top of boundary layer
tran.2D(full.output=TRUE,C=C,D.x=D,v.x=0,VF.x=por,dx=dx,dy=dy,C.bl.up=C.bl.up,a.bl.up=a.bl,C.bl.down=C.bl.up,a.bl.down=a.bl,C.bl.left=C.bl.left,a.bl.left=a.bl,C.bl.right=C.bl.left,a.bl.right=a.bl)  

# Runtime test with and without argument checking

n.iterate <-1000

test1 <- function()
{
for (i in 1:n.iterate )
ST<-tran.2D(full.check=TRUE,C=C,D.x=D,v.x=0,VF.x=por,dx=dx,dy=dy,C.bl.up=C.bl.up,a.bl.up=a.bl,C.down=C.down)  
} 
system.time(test1())

test2 <- function()
{
for (i in 1:n.iterate )
ST<-tran.2D(full.output=TRUE,C=C,D.x=D,v.x=0,VF.x=por,dx=dx,dy=dy,C.bl.up=C.bl.up,a.bl.up=a.bl,C.down=C.down)  
} 
system.time(test2())

test3 <- function()
{
for (i in 1:n.iterate )
ST<-tran.2D(full.output=TRUE,full.check=TRUE,C=C,D.x=D,v.x=0,VF.x=por,dx=dx,dy=dy,C.bl.up=C.bl.up,a.bl.up=a.bl,C.down=C.down)  
} 
system.time(test3())

}
\references{Soetaert and Herman, a guide to ecological modelling - using R as
a simulation platform, 2008. Springer}
\details{The \bold{boundary conditions} are either
\itemize{
  \item (1) zero-gradient 
  \item (2) fixed concentration 
  \item (3) convective boundary layer 
  \item (4) fixed flux 
}
  
  The above order also shows the priority. The default condition is the zero gradient. The convective boundary layer condition overrules the fixed concentration. The fixed flux overrules all other specifications. 
  
  \bold{Transport properties}
  The \emph{diffusion coefficient} (D), the \emph{advective velocity} (v), the \emph{volume fraction} (VF), and the \emph{advective finite difference weight} (AFDW)  can be either one value or a vector.
  When a vector, this vector must be of length N+1, defined at all grid cell interfaces, including upper and lower boundary.

  The \bold{finite difference grid } (grid) is specified either:
\itemize{
  \item as a grid list, as generated by \code{\link{setup.1D.grid}} or 
  \item by the parameter dx representing the thickness of the grid cells (one value or a vector of length N+1)
}
}  
\note{ 
The advective-diffusion equation is not necessarily mass conservative. Sometimes, this is not an issue, for instance when v represents a sinking velocity of particles or a swimming velocity of organisms. In others cases however, mass conservation needs to be accounted for. To ensure mass conservation, the advective velocity must obey certain continuity constraints: in essence the product of the interface surface (A) and advective velocity (v) should be constant. In sediments, one can use 'setup.1D.advection' to ensure that the advective velocities for the pore water and solid phase meet these constraints. 
}
\keyword{utilities}

