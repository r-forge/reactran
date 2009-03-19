\name{setup.prop.2D}
\alias{setup.prop.2D}
\title{
  Attaches a property to a two-dimensional grid
}

\description{
  Calculates the value of a given property at the middle of grid cells
  (\code{mid}) and at the interfaces of the grid cells (\code{int}).

  Two possibilities are available: either specifying a mathematical function
  (\code{func}) that describes the spatial dependency of the property, or
  asssuming a constant value (\code{value}).

  For example, on a sediment model, the routine can be used to specify the
  porosity, the mixing intensity or other parameters over the grid of the
  reactangular sediment domain.
}

\usage{
setup.prop.2D(func = NULL, value = NULL, grid,...)
}

\arguments{
  \item{func }{function that describes the spatial dependency.
  }
  \item{value }{constant value given to the property.
  }
  \item{grid }{list specifying the 2D grid characteristics, see
    \code{\link{setup.grid.2D}} for details on the structure of this list.
  }
  \item{... }{additional arguments that are passed on to \code{func}.
  }
}

\value{
  A list of type \code{prop.2D} containing:
  \item{mid }{property value at the middle of the grid cells; NxM matrix.
  }
  \item{x.int }{property value at the horizontal interfaces of the grid cells;
    (N+1)xM matrix.
  }
  \item{y.int }{property value at the vertical interfaces of the grid cells;
    Nx(M+1) matrix.
  }
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
# Inverse quadratic function 
inv.quad <- function(x,y,a=NULL,b=NULL)
return(1/((x-a)^2+(y-b)^2))


# Construction of the 2D grid 
x.grid <- setup.grid.1D(x.up=0,L=10,N=10)
y.grid <- setup.grid.1D(x.up=0,L=10,N=10)
grid2D <- setup.grid.2D(x.grid,y.grid)

# Attaching the inverse quadratic function to the 2D grid 
(twoD<-setup.prop.2D(func=inv.quad,grid=grid2D,a=5,b=5))
contour(log(twoD$x.int))
}

\keyword{utilities}

