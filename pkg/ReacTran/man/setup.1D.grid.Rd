\name{setup.1D.grid}
\alias{setup.1D.grid}
\alias{plot.grid.1D}
\title{Creation of a one-dimensional finite difference grid}
\description{Subdivides the one-dimensional model domain into several zones with each N grid cells. The resulting grid structure can be used in the transport routines. The grid structure is characterized by the position of the middle of the grid cells (x.mid) and the position of the interfaces between cells (x.int). Distances are calculated between the interfaces (dx), i.e. the thickness of the grid cells. An auxiliary set of distances is calculated (dx.aux) between the points where the concentrations are specified (at the centers of the grid cell and the two external interfaces). A more complex grid can be specified consisting of multiple zones when using vectors as arguments. In each zone, one can control the grid resolution near the upstream and downstraem boundary. The grid resolution at the upstream interface changes according to the formula dx[i+1] = min(max.dx.1,pdx.1*dx[i]). A similar formula controls the resoltion at the downstream interface. } 
\usage{setup.1D.grid(x.up=0, x.down=NULL, L=NULL, N=NULL,
                     dx.1=NULL, p.dx.1 = rep(1,length(L)), max.dx.1 = 10*dx.1,
                     dx.N =NULL, p.dx.N = rep(1,length(L)), max.dx.N = 10*dx.N)

\method{plot}{grid.1D}(x, \dots)

                     }
\arguments{
  \item{x.up }{position of the upstream interface}
  \item{x.down }{position of zone endpoints, one value when the model domain covers only one zone (x.down = position of downstream interface) or a vector of length N, when the model domain is divided into several consecutive zones}
  \item{L }{thickness of zones, one value (model domain = one zone) or a vector of length N}
  \item{N }{desired number of grid cells in a zone, one value or a vector of length N}
  \item{dx.1 }{size of first grid cell in a zone, one value or a vector of length N}
  \item{p.dx.1 }{growth factor of the grid cell size in upper half of the zone, one value or a vector of length N. The default value is 1 (constant grid cell size)}
  \item{max.dx.1 }{maximum grid cell size in upper half of the zone, one value or a vector of length N}
  \item{dx.N }{size of last grid cell in a zone, one value or a vector of length N}
  \item{p.dx.N }{growth factor of the grid cell size in lower half of the zone, one value or a vector of length N. The default value is 1 (constant grid cell size)}
  \item{max.dx.N }{maximum grid cell size in lower half of the zone, one value or a vector of length N}

  \item{x }{an object of class \code{1D.grid} that has to be plotted}
  \item{...}{additional arguments passed to function \code{plot}}

}
\value{
  a list of type \code{grid.1D} containing:
  \item{x.up }{position of the upstream interface}
  \item{x.down }{position of the downstream interface}
  \item{x.mid }{position of the middle of the grid cells, vector of length N}
  \item{x.int }{position of the interfaces of the grid cells, vector of length N+1}
  \item{dx }{distance between adjacent cell interfaces (thickness of grid cells), vector of length N}
  \item{dx.aux }{auxiliary vector containing the distance between adjacent cell centers; at the upper and lower boundary calculated as (x[1]-x.up) and (x.down-x[N]) respectively; vector of length N+1}
  \item{N }{the total number of grid cells}
}
\author{Filip Meysman <filip.meysman@vub.ac.be>, Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
# 1D Grid: one zone, constant resolution 
grid1 <- setup.1D.grid(x.up=0,L=10,N=10)
plot(grid1)

# 1D Grid: one zone, constant resolution, origin not zero 
(GR<-setup.1D.grid(x.up=5,x.down=10,N=10))
plot(GR)

# 1D Grid: one zone, higher resolution near the upstream interface 
(GR<-setup.1D.grid(x.up=0,x.down=10,dx.1=0.1,p.dx.1=1.1))
plot(GR)


# 1D Grid: one zone, higher resolution near the upstream and downstream interface 
GR<-setup.1D.grid(x.up=0,x.down=10,dx.1=0.1,p.dx.1=1.1,dx.N=0.1,p.dx.N=1.1)
plot(GR)

# 1D Grid: two zones, higher resolution near the upstream and downstream interface 
(GR<-setup.1D.grid(x.up=0,L=c(5,5),dx.1=c(0.2,0.2),p.dx.1=c(1.1,1.1),
                   dx.N=c(0.2,0.2),p.dx.N=c(1.1,1.1)))
plot(GR)


# 1D Grid: two zones, higher resolution near the upstream and downstream interface 
# the number of grid cells in each zone is imposed via N  
(GR <- setup.1D.grid(x.up=0,L=c(5,5),N=c(20,10),dx.1=c(0.2,0.2),
                     p.dx.1=c(1.1,1.1),dx.N=c(0.2,0.2),p.dx.N=c(1.1,1.1)))
plot(GR)
}
\keyword{utilities}

