\name{tran.1D.river}
\alias{tran.1D.river}
\title{1-D advective-diffusive transport in a river}
\description{Estimates the transport term (i.e. the rate of change of the concentration due to diffusion and advection) in a one-dimensional model of a river with whose cross sectional area is not necessarily constant, e.g. estuaries}
\usage{tran.1D.river(C, C.up=C[1], C.down=C[length(C)],
C.lat=0, C.lat.grid=list(mid=rep(C.lat,length.out=length(C))), 
F.up=NA, F.down=NA, F.lat=NA, F.lat.grid=list(mid=rep(F.lat,length.out=length(C))),
Disp=NULL, Disp.grid=list(int=rep(Disp,length.out=(length(C)+1))),
flow = 0, flow.grid=list(int=rep(flow,length.out=(length(C)+1))),
flow.lat=0, flow.lat.grid=list(mid=rep(flow.lat,length.out=length(C))),
V=NULL,V.grid=list(mid=rep(V,length.out=length(C))))}
\arguments{
  \item{C }{Tracer concentration, defined at the centre of the grid cells. A vector of length N [M/L3]}
  \item{C.up }{Tracer concentration at the upstream interface. One value [M/L3]}
  \item{C.down }{Tracer concentration at downstream interface. One value [M/L3]}
  \item{C.lat }{Tracer concentration of lateral input, defined at grid cell centres. One value or a vector of length N [M/L3]}
  \item{C.lat.grid }{Tracer concentration in lateral input defined as grid property}
  \item{F.up }{Total tracer input at the upstream interface. One value [M/T]}
  \item{F.down }{Total tracer input at downstream interface. One value [M/T]}
  \item{F.lat }{Total lateral tracer input, defined at grid cell centres. One value or a vector of length N [M/T]}
  \item{F.lat.grid }{Total lateral tracer input defined as grid property}
  \item{Disp }{BULK dispersion coefficient (=disp*A.int/dx.int), defined on grid cell interfaces. One value or a vector of length N+1, [L3/T]}
  \item{Disp.grid }{BULK dispersion coefficient defined as grid property}
  \item{flow }{Water flow rate, defined on grid cell interfaces. One value or a vector of length N+1, [L3/T]}
  \item{flow.grid }{Water flow rate defined as grid property}
  \item{flow.lat }{Lateral water flow rate into each volume box, defined at grid cell centres. One value or a vector of length N, [L3/T]}
  \item{flow.lat.grid }{Lateral water flow rate defined as grid property}
  \item{V }{Grid cell volume, defined at grid cell centres. One value or a vector of length N [L3]}
  \item{V.grid }{Grid cell volume defined as grid property}
}
\value{
  \item{dC }{The rate of change of the concentration C due to transport, defined in the centre of each grid cell [M/L3/T]}
  \item{F }{Mass flow across at the interface of each grid cell. A vector of length N+1 [M/T]}
  \item{F.up }{Mass flow across the upstream boundary, positive = INTO model domain. One value [M/T]}
  \item{F.down }{Mass flow across the downstream boundary, positive = OUT of model domain. One value [M/T]}
  \item{F.lat }{Lateral mass input per volume box, positive = INTO model domain. A vector of length N [M/T]}
}
\author{Filip Meysman <f.meysman@nioo.knaw.nl>, Karline Soetaert <k.soetaert@nioo.knaw.nl> }
\examples{

##################################################################
#  EXAMPLE : organic carbon (OC) decay in a widening river       #
##################################################################

# Two scenarios are simulated: the baseline includes only input 
# of organic matter upstream. The second scenario simulates the 
# input of an important side river half way the estuary.  

#====================#
# Model formulation  #
#====================#

river.model <- function (t=0,OC,pars=NULL)
{
tran <- tran.1D.river(C=OC,F.up=F.OC,F.lat=F.lat,Disp=Disp,flow=flow,V=Volume)$dC 
reac <- - k*OC
return(list(dCdt = tran + reac))
}

#======================#
# Parameter definition #
#======================#

# Initialising morphology estuary: 

nbox          <- 500     # number of grid cells
lengthEstuary <- 100000  # length of estuary [m]
BoxLength     <- lengthEstuary/nbox # [m]
Distance      <- seq(BoxLength/2, by=BoxLength, len=nbox) # [m]

# Cross sectional area: sigmoid function of estuarine distance [m2]
CrossArea <- 4000 + 72000 * Distance^5 /(Distance^5+50000^5)

# Volume of boxes                          (m3)
Volume  <- CrossArea*BoxLength

# Transport coefficients
Disp    <- 1000   # m3/s, bulk dispersion coefficient
flow    <- 180    # m3/s, mean river flow

F.OC    <- 180               # input organic carbon [mol s-1]
F.lat.0 <- F.OC              # lateral input organic carbon [mol s-1]

k       <- 10/(365*24*3600)  # decay constant organic carbon [s-1]


#====================#
# Model solution     #
#====================#

F.lat <- rep(0,length.out=nbox)
Conc1 <- steady.band(runif(nbox),fun=river.model,nspec=1)$y   
F.lat <- F.lat.0*dnorm(x=Distance/lengthEstuary, mean = Distance[nbox/2]/lengthEstuary, sd = 1/20, log = FALSE)/nbox 
Conc2 <- steady.band(runif(nbox),fun=river.model,nspec=1)$y   

#====================#
# Plotting output    #
#====================#

matplot(Distance/1000,cbind(Conc1,Conc2),lwd=2,
main="Organic carbon decay in an estuary",xlab="distance [km]",ylab="OC Concentration [mM]",
type="l")  
legend ("topright",lty=1,col=c("black","red"),
        c("baseline","with lateral input"))
}
\references{Soetaert and Herman, a guide to ecological modelling - using R as a simulation platform, 2008. Springer.}
\details{The \bold{boundary conditions} are of type
\itemize{
  \item 1. zero-gradient (default)
  \item 2. fixed concentration 
  \item 3. fixed input 
}
  The \emph{bulk dispersion coefficient} (Disp), and \emph{flow rate} (flow) can be either one value or a vector of length N+1, defined at all grid cell interfaces, including upstream and downstream boundary.
  
  The spatial discretisation is given by the volume of each box (V), which can be one value or a vector of length N+1, defined at the centre of each grid cell.}
\seealso{\code{\link{tran.1D}}}
\keyword{utilities}

