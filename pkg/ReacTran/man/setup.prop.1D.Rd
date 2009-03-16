\name{setup.prop.1D}
\alias{setup.prop.1D}
\title{
  Attaches a property to a one-dimensional grid
}

\description{
  This routine calculates the value of a given property at the middle of
  grid cells (\code{mid}) and at the interfaces of the grid cells (\code{int}).

  Two possibilities are available: either specifying a mathematical function
  (\code{func}) that describes the spatial dependency of the property,
  or the property is derived from interpolation of a data series
  (data matrix \code{xy}).

  For example, in a sediment model, \code{setup.prop.1D} can be used to
  specify the porosity, the mixing intensity or some other parameter over
  the one-dimensional grid.

  Similarly, in a vertical water column model, \code{setup.prop.1D} can be
  used to specify the sinking velocity of particles or any other model parameter
  varying with water depth.
}

\usage{
setup.prop.1D(func=NULL, value=NULL, xy=NULL,
  interpolate="spline", grid, ...)
}

\arguments{
  \item{func }{function that describes the depth dependency.
  }
  \item{value }{constant value given to the property.
  }
  \item{xy }{a two-column matrix where the first column (\code{x}) provides
    the position, and the second column (\code{y}) provides the values that
    need interpolation over the grid.
  }
  \item{interpolate }{specifies how the interpolation should be done, one
    of "spline" or "linear"; only used if \code{xy} is present.
  }
  \item{grid }{list specifying the 1D grid characteristics, see
    \code{\link{setup.grid.1D}} for details on the structure of this list.
  }
  \item{... }{additional arguments that are passed on to \code{func}.
  }
}

\value{
  A list of type \code{prop.1D} containing:
  \item{mid }{property value at the middle of the grid cells; vector of
    length N (where N is the number of grid cells).
  }
  \item{int }{property value at at the interface of the grid cells;
    vector of length N+1.
  }
}

\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>,
  Filip Meysman <f.meysman@nioo.knaw.nl>
}

\examples{

# Construction of the 1D grid 

grid <- setup.grid.1D(x.up=0,L=10,N=10)

# Porosity profile via function specification

exp.profile <- function(x,y.0=NULL,y.inf=NULL,x.att=NULL)
{return(y.inf + (y.0-y.inf)*exp(-x/x.att))}

P.func <- setup.prop.1D(func=exp.profile,grid=grid,y.0=0.9,
y.inf=0.5,x.att=3)

# Porosity profile via data series interpolation

P.data <- matrix(ncol=2,data=c(0,3,6,10,0.9,0.65,0.55,0.5))
P.spline <- setup.prop.1D(xy=P.data,grid=grid)
P.linear <- setup.prop.1D(xy=P.data,grid=grid,interpolate="linear")

# Plot different profiles 

plot(grid$x.int,P.func$int,type="l",
main="setup.prop, function evaluation")  
points(P.data,cex=1.5,pch=16) 
lines(grid$x.int,P.spline$int,lty="dashed")
lines(grid$x.int,P.linear$int,lty="dotdash")
}

\details{
  There are two options to carry out the interpolation:
  \itemize{
     \item "spline" gives a smooth profile, but sometimes generates
       strange profiles - always check the result!
     \item "linear" gives a segmented profile
  }
}

\keyword{utilities}

