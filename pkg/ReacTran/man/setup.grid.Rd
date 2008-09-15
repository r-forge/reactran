\name{setup.grid}
\alias{setup.grid}
\title{Create a numerical grid}
\description{Subdivides a certain size into grid cells, usable for transport routines. Grid cells are characterised by 
  \item thicknesses (dx), and distances from centre to centre (dx.int)
  \item distance from upstream boundary to middle of cell (x) and from upstream boundary to upstream interface of cell (x.int)}
\usage{setup.grid(len=NULL, N=NULL, dx.1=NULL, layer=NULL, type="ct",
dx.int.up=NULL,dx.int.down=NULL,x.0=0)}
\arguments{
  \item{len }{total size to divide; not used if type = "layer"}
  \item{N }{number of grid cells; not used if type = "layer"; can be null if type == "ct" and dx.1 is given a value}
  \item{dx.1 }{size of first grid cell; only used if type = "ct", "exp" or "sine"; in the latter case, it is also the size of the last grid cell.
    If type =="ct" and both dx.1 and N are specified, then N is re-estimated as round(len/dx.1); if N is rounded, then also dx.1 will be re-estimated}
  \item{layer }{3-columned matrix, specifying the number of cells (col 1), the size of the first grid cell (col 2) and the total length (col 3) 
    in a number of discrete layers; usually the size of the first grid cell will only be specified for the first layer, it can be set NA for the other layers.}
  \item{type }{type of grid, one of "ct","exp","sine","layer", for constant, exponential, sine and layer type respectively}
  \item{dx.int.up }{the "interface" distance to the upstream boundary, i.e. the distance from the centre of 1st box to the upper interface; the default is half the size of the first box.}
  \item{dx.int.down }{the "interface" distance to the downstream boundary, i.e. the distance from the centre of last box to the lower interface; the default is half the size of the first box.}
  \item{x.0 }{distance of upstream interface; defaults to 0 (start of x-axis)}
}
\value{
  a list containing:
  \item{x }{distance from origin to centre of boxes, vector of length N}
  \item{x.int }{distance from origin to box upper interface, vector of length N+1}
  \item{dx }{thickness of boxes, vector length N}
  \item{dx.int }{distance over the interface, i.e. between centre of adjacent boxes, vector of length N+1}
  \item{N }{the number of boxes}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# example
Exp    <- setup.grid(len=100,dx.1=1,N=20,type="exp")
Sine   <- setup.grid(len=100,dx.1=1,N=20,type="sine")
Ct     <- setup.grid(len=100,dx.1=1,N=20,type="ct")

# grid composed of two layers; with the following specifications:
# layer 1: 10 boxes, size of 1st box= 1, total size in layer= 40;
# layer 2: 10 boxes, size of 1st box not specified (NA), total size=60 
Twolay <- setup.grid(len=100,layer=matrix(ncol=3,byrow=TRUE,
                      data=c(10,1,40,10,NA,60)),type="layer")

plot(Exp$dx,main="setup.grid",ylab="dx",xlab="box")  
plot(Sine$dx,pch=16,main="setup.grid",ylab="dx",xlab="box")  
plot(Ct$dx,pch=18,main="setup.grid",ylab="dx",xlab="box")  
plot(Twolay$dx,pch=20,main="setup.grid",ylab="dx",xlab="box")}
\details{The interface distance dx.int.up and dx.int.down are by the default is half the size of the first/last box. 
  \cr For sediments for instance, this can be set equal to the size of the diffusive boundary layer (+ 0.5*dx.1)}
\keyword{utilities}

