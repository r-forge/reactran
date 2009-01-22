\name{setup.2D.prop}
\alias{setup.2D.prop}
\title{Attaches a property to a two-dimensional grid}
\description{This routine calculates the value of a given property at the middle of grid cells ('mid') and at the interfaces of the grid cells ('int'). Two possibilities are available: either specifying a mathematical function ('func') that describes the spatial dependency or by interpolation of a data series (data matrix 'xy'). For example, the routine can be used to specify the porosity, the mixing intensity or other parameters over a grid.}
\usage{setup.2D.prop(func = (F <- function(x) return(rep(value,times=length(x)))), value = NULL, grid,...)}
\arguments{
  \item{func }{Function that describes the spatial dependency}
  \item{value }{Constant value given to the property}
  \item{grid }{List specifying the 2D grid characteristics, see 'setup.2D.grid' for details on the structure of this list }
  \item{... }{Additional arguments that are passed on to func }
}
\value{
  A list containing:
  \item{mid }{property value at the middle of the grid cells, NxN matrix}
  \item{x.int }{property value at the horizontal interfaces of the grid cells, (N+1)xN matrix}
  \item{y.int }{property value at the vertical interfaces of the grid cells, Nx(N+1) matrix}
}
\author{Filip Meysman <f.meysman@nioo.knaw.nl>, Karline Soetaert <k.soetaert@nioo.knaw.nl> }
\examples{
# Inverse quadratic function 
inv.quad <- function(x,y,a=NULL,b=NULL)
{
return(1/((x-a)^2+(y-b)^2))
}

# Construction of the 2D grid 
x.grid <- setup.1D.grid(x.up=0,L=10,N=10)
y.grid <- setup.1D.grid(x.up=0,L=10,N=10)
grid2D <- setup.2D.grid(x.grid,y.grid)

# Attaching the inverse quadratic function to the 2D grid 
setup.2D.prop(func=inv.quad,grid=grid2D,a=5,b=5)

}
\keyword{utilities}

