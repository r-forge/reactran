\name{setup.2D.grid}
\alias{setup.2D.grid}
\title{Creation of a two-dimensional finite difference grid}
\description{Creates a finite difference grid over a two-dimensional model domain staring from two separate one-dimensional grids (for example as created by setup.1D.grid).} 
\usage{setup.2D.grid(x.grid=NULL, y.grid=NULL)}
\arguments{
  \item{x.grid }{list containing the one-dimensional grid in the vertical direction - see setup.1D.grid for the structure of the list}
  \item{y.grid }{list containing the one-dimensional grid in the horizontal direction - see setup.1D.grid for the structure of the list}
}
\value{
  a list containing:
  \item{x.up }{vertical position of the upper interface}
  \item{x.down }{vertical position of the lower interface}
  \item{x.mid }{vertical position of the middle of the grid cells, vector of length N}
  \item{x.int }{vertical position of the horizontal interfaces of the grid cells, vector of length N+1}
  \item{dx }{distance between adjacent cell interfaces (thickness of grid cells), vector of length N}
  \item{dx.aux }{auxiliary vector containing the distance between adjacent cell centers; at the upper and lower boundary calculated as (x[1]-x.up) and (x.down-x[N]) respectively; vector of length N+1}
  \item{x.N }{the total number of grid cells in the vertical direction}
  \item{y.left }{horizontal position of the left interface}
  \item{y.right }{horizontal position of the right interface}
  \item{y.mid }{horizontal position of the middle of the grid cells, vector of length N}
  \item{y.int }{horizontal position of the vertical interfaces of the grid cells, vector of length N+1}
  \item{dy }{distance between adjacent cell interfaces (thickness of grid cells), vector of length N}
  \item{dy.aux }{auxiliary vector containing the distance between adjacent cell centers; at the left and right boundary calculated as (y[1]-xyup) and (y.down-y[N]) respectively; vector of length N+1}
  \item{y.N }{the total number of grid cells in the horizontal direction}
}

\author{Filip Meysman <filip.meysman@vub.ac.be>, Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
# test of the setup.2Dgrid functionality
x.grid <- setup.1D.grid(x.up=0,L=10,N=5)
y.grid <- setup.1D.grid(x.up=0,L=20,N=10)
grid2D <- setup.2D.grid(x.grid,y.grid)
}
\keyword{utilities}

