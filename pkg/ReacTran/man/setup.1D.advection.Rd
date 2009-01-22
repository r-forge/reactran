\name{setup.1D.advection}
\alias{setup.1D.advection}
\title{Calculates the advective velocities of the pore water and the solid phase in a sediment assuming steady state compaction}
\description{This routine calculates the advective velocities of the pore water and the solid phase in a sediment based on the assumption of steady state compaction.  The velocities of the pore water ('u') and the solid phase ('v') are calculated in the middle ('mid') of the grid cells and the interfaces ('int'). One requires the specification of the porosity at the interface and at infinite depth, as well as the advective velocity of teh solid phase at the interface. }
\usage{setup.1D.advection(v.0 = NULL, v.inf = NULL, por.0, por.inf, por.grid)}
\arguments{
  \item{v.0 }{Advective velocity of the solid phase at the sediment-water interface (also referred to as the sedimentation velocity)}
  \item{v.inf }{Advective velocity of the solid phase at infinite depth (also referred to as the burial velocity)}
  \item{por.0 }{Porosity at the sediment-water interface}
  \item{por.inf }{Porosity at infinite depth }
  \item{por.grid }{Porosity profile specified as a 1D grid property (see 'setup.1D.prop' for details on the structure of this list)}
}
\value{
  A list containing:
  \item{u }{list with pore water advective velocities at the middle of the grid cells ('u$mid') and at the interfaces ('u$int')}
  \item{v }{list with solid phase advective velocities at the middle of the grid cells ('u$mid') and at the interfaces ('u$int')}
}
\author{Filip Meysman <f.meysman@nioo.knaw.nl>, Karline Soetaert <k.soetaert@nioo.knaw.nl> }
\examples{

# setup of the 1D grid 
L <- 10
grid <- setup.1D.grid(x.up=0,L=L,N=20)

# attaching an exponential porosity profile to the 1D grid

exp.profile <- function(x,y.0=NULL,y.inf=NULL,x.att=NULL)
{return(y.inf + (y.0-y.inf)*exp(-x/x.att))}

por.grid <- setup.1D.prop(func=exp.profile,grid=grid,y.0=0.9,y.inf=0.5,x.att=3)

# calculate the advective velocities

dummy <- setup.1D.advection(v.0 = 1, por.0=0.9, por.inf=0.5, por.grid = por.grid)
u.grid <-dummy$u
v.grid <-dummy$v

# plotting the results 

par(mfrow=c(2,1),cex.main=1.2,cex.lab=1.2)

matplot(por.grid$int,grid$x.int,pch=19,ylim=c(L,0), xlim=c(0,1),
xlab="",ylab="depth [cm]",main=expression("porosity"),axes=FALSE)
abline(h = 0)
axis(pos=0, side=2)
axis(pos=0, side=3)

matplot(u.grid$int,grid$x.int,type="l",lwd=2,col="blue",ylim=c(L,0), xlim=c(0,max(u.grid$int,v.grid$int)),
xlab="",ylab="depth [cm]",main=expression("advective velocity [cm yr-1]"),axes=FALSE)
abline(h = 0)
axis(pos=0, side=2)
axis(pos=0, side=3)
lines(v.grid$int,grid$x.int,lwd="2",col="red")
legend(x="bottomright", legend=c("pore water","solid phase"),col=c("blue","red"),lwd=c(2,2))

}
\keyword{utilities}

