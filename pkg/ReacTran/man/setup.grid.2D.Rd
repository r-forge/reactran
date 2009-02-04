\name{setup.grid.2D}
\alias{setup.grid.2D}
\title{Creation of a two-dimensional finite difference grid}
\description{Creates a finite difference grid over a rectangular two-dimensional model domain starting from two separate one-dimensional grids (as created by \code{\link{setup.grid.1D}}).} 
\usage{setup.grid.2D(x.grid=NULL, y.grid=NULL)}
\arguments{
  \item{x.grid }{list containing the one-dimensional grid in the vertical direction - see \code{\link{setup.grid.1D}} for the structure of the list}
  \item{y.grid }{list containing the one-dimensional grid in the horizontal direction - see \code{\link{setup.grid.1D}} for the structure of the list}
}
\value{
  a list of type \code{grid.2D} containing:
  \item{x.up }{vertical position of the upper interface; one value}
  \item{x.down }{vertical position of the lower interface; one value}
  \item{x.mid }{vertical position of the middle of the grid cells; vector of length \code{x.N}}
  \item{x.int }{vertical position of the horizontal interfaces of the grid cells; vector of length \code{x.N+1}}
  \item{dx }{distance between adjacent cell interfaces (thickness of grid cells); vector of length \code{x.N}}
  \item{dx.aux }{auxiliary vector containing the distance between adjacent cell centers; at the upper and lower boundary calculated as (\code{x[1]-x.up}) and (\code{x.down-x[x.N]}) respectively; vector of length \code{x.N+1}}
  \item{x.N }{total number of grid cells in the vertical direction; one value}
  \item{y.left }{horizontal position of the left interface; one value}
  \item{y.right }{horizontal position of the right interface; one value}
  \item{y.mid }{horizontal position of the middle of the grid cells; vector of length \code{y.N}}
  \item{y.int }{horizontal position of the vertical interfaces of the grid cells; vector of length \code{y.N+1}}
  \item{dy }{distance between adjacent cell interfaces (thickness of grid cells); vector of length \code{y.N}}
  \item{dy.aux }{auxiliary vector containing the distance between adjacent cell centers; at the left and right boundary calculated as (\code{y[1]-y.left}) and (\code{y.right-y[y.N]}) respectively; vector of length \code{y.N+1}}
  \item{y.N }{total number of grid cells in the horizontal direction; one value}
}

\author{Filip Meysman <f.meysman@nioo.knaw.nl>, Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
# test of the setup.grid.2D functionality
x.grid <- setup.grid.1D(x.up=0,L=10,N=5)
y.grid <- setup.grid.1D(x.up=0,L=20,N=10)
(grid2D <- setup.grid.2D(x.grid,y.grid))
}
\keyword{utilities}

