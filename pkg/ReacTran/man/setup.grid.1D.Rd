\name{setup.grid.1D}
\alias{setup.grid.1D}
\alias{plot.grid.1D}
\title{
  Creation of a one-dimensional finite difference grid
}
\description{
  Subdivides the one-dimensional model domain into one or more zones that
  are each sub-divided into grid cells.

  The resulting grid structure can be used in the other \code{ReacTran}
  functions.

  The grid structure is characterized by the position of the middle of
  the grid cells (\code{x.mid}) and the position of the interfaces between
  cells (\code{x.int}).

  Distances are calculated between the interfaces (\code{dx}), i.e. the
  thickness of the grid cells.
  An auxiliary set of distances is calculated (\code{dx.aux}) between the
  points where the concentrations are specified (at the centers of the grid
  cell and the two external interfaces).

  A more complex grid consisting of multiple zones can be constructed when
  using vectors as arguments. In each zone, one can control the grid
  resolution near the upstream and downstream boundary.
  The grid resolution at the upstream interface changes according to the
  formula \code{dx[i+1] = min(max.dx.1,p.dx.1*dx[i])}.

  A similar formula controls the resolution at the downstream interface.
  This allows refinement of the grid near the interfaces.
}

\usage{
setup.grid.1D(x.up=0, x.down=NULL, L=NULL, N=NULL, dx.1=NULL,
  p.dx.1=rep(1,length(L)), max.dx.1=L, dx.N=NULL,
  p.dx.N=rep(1,length(L)), max.dx.N=L)

\method{plot}{grid.1D}(x, \dots)

}

\arguments{
  \item{x.up }{position of the upstream interface; one value.
  }
  \item{x.down }{position of the endpoint of each zone; one value when the
    model domain covers only one zone (\code{x.down} = position of downstream
    interface), or a vector of length M when the model domain is divided into
    M zones (\code{x.down[M]} = position of downstream interface).
  }
  \item{L }{thickness of zones; one value (model domain = one zone) or a
    vector of length M (model domain = M zones).
  }
  \item{N }{number of grid cells within a zone; one value or a vector of
    length M.
  }
  \item{dx.1 }{size of the first grid cell in a zone; one value or a vector
    of length M.
  }
  \item{p.dx.1 }{factor controlling the increase in grid cell size near
    the upstream boundary; one value or a vector of length N. The default
    value is 1 (constant grid cell size).
  }
  \item{max.dx.1 }{maximum grid cell size in the upstream half of the zone;
    one value or a vector of length M.
  }
  \item{dx.N }{size of the last grid cell in a zone; one value or a vector
    of length M.
  }
  \item{p.dx.N }{factor controlling the increase in grid cell size near the
    downstream boundary; one value or a vector of length N. The default value
    is 1 (constant grid cell size).
  }
  \item{max.dx.N }{maximum grid cell size in the downstream half of the zone;
    one value or a vector of length M.
  }
  
  \item{x }{the object of class \code{grid.1D} that needs plotting.
  }
  \item{...}{additional arguments passed to the function \link{plot}.
  }
}

\value{
  a list of type \code{grid.1D} containing:

  \item{N }{the total number of grid cells.
  }
  \item{x.up }{position of the upstream interface; one value.
  }
  \item{x.down }{position of the downstream interface; one value.
  }
  \item{x.mid }{position of the middle of the grid cells;
    vector of length \code{N}.
  }
  \item{x.int }{position of the interfaces of the grid cells;
    vector of length \code{N+1}.
  }
  \item{dx }{distance between adjacent cell interfaces (thickness of grid
    cells), vector of length \code{N}.
  }
  \item{dx.aux }{auxiliary vector containing the distance between adjacent
    cell centers; at the upper and lower boundary calculated as
    (\code{x[1]-x.up}) and (\code{x.down-x[N]}) respectively; vector of
    length \code{N+1}.
  }
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
# 1D Grid: one zone, constant resolution 
(GR <- setup.grid.1D(x.up=0,L=10,N=10))
plot(GR)

# 1D Grid: one zone, constant resolution, origin not zero 
(GR<-setup.grid.1D(x.up=5,x.down=10,N=10))
plot(GR)

# 1D Grid: one zone, higher resolution near the upstream interface 
(GR<-setup.grid.1D(x.up=0,x.down=10,dx.1=0.1,p.dx.1=1.1))
plot(GR)


# 1D Grid: one zone, higher resolution near the upstream
# and downstream interface 
GR<-setup.grid.1D(x.up=0,x.down=10,
dx.1=0.1,p.dx.1=1.1,dx.N=0.1,p.dx.N=1.1)
plot(GR)

# 1D Grid: two zones, higher resolution near the upstream
# and downstream interface 
(GR<-setup.grid.1D(x.up=0,L=c(5,5),dx.1=c(0.2,0.2),p.dx.1=c(1.1,1.1),
                   dx.N=c(0.2,0.2),p.dx.N=c(1.1,1.1)))
plot(GR)


# 1D Grid: two zones, higher resolution near the upstream
# and downstream interface 
# the number of grid cells in each zone is imposed via N  
(GR <- setup.grid.1D(x.up=0,L=c(5,5),N=c(20,10),dx.1=c(0.2,0.2),
             p.dx.1=c(1.1,1.1),dx.N=c(0.2,0.2),p.dx.N=c(1.1,1.1)))
plot(GR)
}
\keyword{utilities}

