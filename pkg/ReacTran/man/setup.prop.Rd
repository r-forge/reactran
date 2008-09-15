\name{setup.prop}
\alias{setup.prop}
\title{Creates a profile of properties over a grid}
\description{Estimates the values of certain properties at the centre of boxes and box interfaces, 
  at upstream boundary and at infinite depth. 
  \cr Input is either as a et of coefficient defining and exponential function or by xy-values   
  
  For instance for porous media, the function can be used to specify porosities, or bioturbation coefficients over a numerical grid.}
\usage{setup.prop(y.0=NULL, y.INF=y.0, y.coeff=0, L=0, 
  xy=NULL, type="spline", grid)}
\arguments{
  \item{y.0 }{value of property y at upstream interface (x=0)}
  \item{y.INF }{value of property y at infinite distance (x=inf) }
  \item{y.coeff }{rate of change with distance x}
  \item{L }{if specified, the thickness of the (upper) layer where the property has value y0; by default this layer is not present}
  \item{xy }{if present, a two-columned matrix, with x-y (distance-property) values, which will be inter/extrapolated to the discretisation grid
    overrules all other settings}
  \item{type }{only if xy is present: how the inter/extrapolation should be done, one of "spline" or "linear"}
  \item{grid }{Discretisation grid, a list containing at least elements dx and dx.int, e.g. as calculated by setup.grid, [L]}
}
\value{
  a list containing:
  \item{up }{property value at upstream boundary, one value}
  \item{INF }{property value at infinite distance, one value}
  \item{mid }{property value at middle of boxes, vector of length N}
  \item{int }{property value at upstream box interfaces, vector of length N+1}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# example
Grid    <- setup.grid(len=10,N=50)	
Poro    <- setup.prop(y.0=0.95,y.INF=0.8,y.coeff=1,grid=Grid)
Db      <- setup.prop(y.0=1,y.INF=0,y.coef=1,L=5,grid=Grid)
mat     <- matrix(ncol=2,data=c(1,3,5,9,0.1,0.5,0.6,0.9))
Pspline <- setup.prop(xy=mat,grid=Grid)
Plin    <- setup.prop(xy=mat,grid=Grid,type="linear")

plot(Grid$x,Poro$mid, type="l",main="setup.prop, L=0")  
plot(Grid$x,Db$mid, type="l",main="setup.prop, L=5")  
plot(Grid$x,Pspline$mid, type="l",main="setup.prop, spline interpolation") 
points(mat,cex=1.5,pch=16) 
plot(Grid$x,Plin$mid, type="l",main="setup.prop, linear interpolation")  
points(mat,cex=1.5,pch=16)}
\details{If input is as data, the data should be specified in xy, a 2-columned matrix with x(distances, 1st column) and y(property values, 2nd column).
 These data are then inter and extrapolated to the numerical grid. To interpolate, there are two options: 
 \item "spline" gives a smooth profile, but can generate strange values - always check the result!
 \item "linear" gives a segmened profile
 
 If input is specified as an exponential it is estimated as:
 \deqn{y= yInf + (y0-yInf)*exp(x*abs(ycoeff))} for x>L
  \deqn{y = y0} for x<=L}
\keyword{utilities}

