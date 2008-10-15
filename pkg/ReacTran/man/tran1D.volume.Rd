\name{tran1D.volume}
\alias{tran1D.volume}
\title{1-D advective-diffusive volumetric transport}
\description{Estimates the rate of change of substances due to dispersive and advective 
  transport in 1-D water bodies whose volume is not necessarily constant, e.g. estuaries}
\usage{tran1D.volume(y, y.up=y[1], y.down=y[length(y)], 
input.up=NA, input.down=NA, flow, flow.up=0, 
weight.up=1, bulkdisp, volume)}
\arguments{
  \item{y }{concentration (or other variable), defined in centre of boxes. A vector of length N, [Mass/Length3]}
  \item{y.up }{concentration (or other variable), at upstream boundary interface. One value, [M/L3]}
  \item{y.down }{concentration (or other variable), at downstream boundary interface. One value, [M/L3]}
  \item{input.up }{total input across the upper boundary interface, positive = IN system. One value, [M/Time]}
  \item{input.down }{total input across the lower boundary interface, positive= IN of system. One value, [M/T]}
  \item{flow }{flow rate, defined on box interfaces. One value or a vector of length N+1, [L3/T]}
  \item{flow.up }{upward flow rate, against the axis, defined on box interfaces. One value or a vector of length N+1, [L3/T]}
  \item{weight.up }{upstream weighing coefficient for flow and upward flow, defined on box interfaces; default is backward differences. One value or a vector of length N+1, [-]}
  \item{bulkdisp }{BULK dispersion coefficient (=disp*A.int/dx.int), defined on box interfaces. One value or a vector of length N+1, [L3/T]}
  \item{volume }{box volume. One value or a vector of length N, [L3]}
}
\value{
  a list containing:
  \item{dy }{rate of change of y in each layer due to transport, [M/L3/T]}
  \item{input }{Total inputs across each box interface. A vector of length N+1, [M/Time]}
  \item{input.up }{Total input across the upper boundary interface, positive = IN system. One value, [M/Time]}
  \item{input.down }{Total input across the lower boundary interface, positive= IN system. One value, [M/T]}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{##################################################################
######   EXAMPLE : Boundary conditions                      ######
##################################################################

#====================#
# Model equations    #
#====================#

# two imposed boundary concentrations and rate
model1 <- function (t,Conc,pars=NULL,y.up,y.down,rate)   
   return (list(tran1D.volume(Conc,y.up=y.up,y.down=y.down,flow=Flow,
                         bulkdisp=disp,volume=Volume)$dy-Conc*rate))

# two imposed boundary concentrations and rate
model2 <- function (t,Conc,pars=NULL,y.up,rate)   
   return (list(tran1D.volume(Conc,y.up=y.up,flow=Flow,
                         bulkdisp=disp,volume=Volume)$dy-Conc*rate))

#====================#
# Model application  #
#====================#

# Initialising morphology: #

nbox       <- 500 

# parameters defining the morphology

lengthEstuary <- 100000                # m total length of estuary
BoxLength     <- lengthEstuary/nbox
Distance      <- seq(BoxLength/2,by=BoxLength, len=nbox)

# Cross section at middle of boxes         (m2)
# cross sectional surface area: sigmoid function of estuarine distance
CrossSurfArea <- 4000 + 72000 * Distance^5 /(Distance^5+50000^5)

# Volume of boxes                          (m3)
Volume  <- CrossSurfArea*BoxLength

# Transport coefficients
disp    <- 50000   # m3/s, bulk dispersion coefficient
Flow    <- 180     # m3/s, mean river flow

# RUNNING the model:   #


Conc1     <- steady.band(runif(nbox),fun=model1,y.up=0, y.down=50,
                         rate=0,nspec=1)   
Conc2     <- steady.band(runif(nbox),fun=model1,y.up=0, y.down=50,
                         rate=0.025/3600/24,nspec=1)   
Conc3     <- steady.band(runif(nbox),fun=model1,y.up=50,y.down=50,
                         rate=0.025/3600/24,nspec=1)  
# model 2                         
Conc4     <- steady.band(runif(nbox),fun=model2,y.up=50,
                         rate=0.025/3600/24,nspec=1)  

#====================#
# Plotting output    #
#====================#

matplot(Distance,cbind(Conc1$y,Conc2$y,Conc3$y,Conc4$y),lwd=2,
main="tran1D.volume",xlab="distance, km",ylab="Concentration",
type="l",sub="reactive transport in estuary")  
legend ("top",c("rate=0/d","rate=0.025/d","rate=0.025/d",
                 "rate=0.025/d, 0-gradient downstream"),
         title="Boundaries",lty=1:5,col=1:5,lwd=2)}
\references{Soetaert and Herman, a guide to ecological modelling - using R as a simulation platform, 2008. Springer.}
\details{The \bold{boundary conditions} are either of type
  \item 1. zero-gradient (default)
  \item 2. concentration boundary 
  \item 3. flux boundary 

  If the flux boundary is specified, it overrules the other specification  
  
  Often both boundaries will consist of a concentration boundary

  \bold{transport properties}
  The \emph{bulk dispersion coefficient} (disp), and \emph{flow rate}, flow can be either one value or a vector.
  
  If they are a vector, they must be of length N+1, defined at all box interfaces, including upstream and downstream boundary.


  The spatial discretisation is given by:
  \item volume, the volume of each box; either one number or one value for each box.}
\seealso{\code{\link{tran1D.solute}}, \code{\link{tran1D.solid}}, \code{\link{tran1D}}}
\keyword{utilities}

