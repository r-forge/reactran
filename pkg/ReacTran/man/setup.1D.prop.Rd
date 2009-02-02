\name{setup.1D.prop}
\alias{setup.1D.prop}
\title{Attaches a property to a one-dimensional grid}
\description{This routine estimates the value of a given property at the middle of grid cells ('mid') and at the interfaces of the grid cells ('int'). Two possibilities are available: either specifying a mathematical function ('func') that describes the spatial dependency or by interpolation of a data series (data matrix 'xy'). For example, the routine can be used to specify the porosity, the mixing intensity or other parameters over a grid.}
\usage{setup.1D.prop(func = (F <- function(x) return(rep(value,times=length(x)))), value=NULL, xy = NULL, type = "spline", grid,...)}
\arguments{
  \item{func }{Function that describes the depth dependency}
  \item{value }{Constant value given to the property}
  \item{xy }{A two-column matrix where the first column ('x') provides the position, and the second column ('y') provides the values that need interpolation over the grid}
  \item{type }{Specifies how the interpolation should be done, one of "spline" or "linear"; only used if xy is present}
  \item{grid }{List specifying the 1D grid characteristics, see 'setup.1D.grid' for details on the structure of this list }
  \item{... }{Additional arguments that are passed on to func }
}
\value{
  A list containing:
  \item{mid }{property value at the middle of the grid cells, vector of length N}
  \item{int }{property value at at the interface of the grid cells, vector of length N+1}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>, Filip Meysman <filip.meysman@vub.ac.be>}
\examples{

# Construction of the 1D grid 

grid <- setup.1D.grid(x.up=0,L=10,N=10)

# Porosity profile via function specification

exp.profile <- function(x,y.0=NULL,y.inf=NULL,x.att=NULL)
{return(y.inf + (y.0-y.inf)*exp(-x/x.att))}

P.func <- setup.1D.prop(func=exp.profile,grid=grid,y.0=0.9,y.inf=0.5,x.att=3)

# Porosity profile via data series interpolation

P.data <- matrix(ncol=2,data=c(0,3,6,10,0.9,0.65,0.55,0.5))
P.spline <- setup.1D.prop(xy=P.data,grid=grid)
P.linear <- setup.1D.prop(xy=P.data,grid=grid,type="linear")

# Plot different profiles 

plot(grid$x.int,P.func$int,type="l",main="setup.prop, function evaluation")  
points(P.data,cex=1.5,pch=16) 
lines(grid$x.int,P.spline$int,lty="dashed")
lines(grid$x.int,P.linear$int,lty="dotdash")
}
\details{
To interpolate, there are two options:
\itemize{
 \item "spline" gives a smooth profile, but sometimes generates strange profiles - always check the result!
 \item "linear" gives a segmented profile
}
}
\keyword{utilities}

