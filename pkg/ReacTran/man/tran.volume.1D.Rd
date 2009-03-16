\name{tran.volume.1D}
\alias{tran.volume.1D}
\title{
  1-D volumetric advective-diffusive transport in an aquatic system
}

\description{
  Estimates the volumetric transport term (i.e. the rate of change of the
  concentration due to diffusion and advection) in a one-dimensional model of
  an aquatic system (river, estuary).

  Volumetric transport implies the use of flows (mass per unit of time) rather
  than fluxes (mass per unit of area per unit of time) as is done in
  \code{\link{tran.1D}}.

  The \code{tran.volume.1D} routine is particularly suited for modelling
  channels (like rivers, estuaries) where the cross-sectional area changes,
  but where this area chnage is not explicitly modelled as such.

  Another difference withy \code{tran.1D} is that the present routine also
  allows the input of lateral water flow.
}

\usage{
tran.volume.1D(C, C.up=C[1], C.down=C[length(C)],
  C.lat=0, C.lat.grid=list(mid=rep(C.lat,length.out=length(C))),
  F.up=NA, F.down=NA,
  F.lat=NA, F.lat.grid=list(mid=rep(F.lat,length.out=length(C))),
  Disp=NULL, Disp.grid=list(int=rep(Disp,length.out=(length(C)+1))),
  flow = 0, flow.grid=list(int=rep(flow,length.out=(length(C)+1))),
  flow.lat=0, flow.lat.grid=list(mid=rep(flow.lat,length.out=length(C))),
  V=NULL,V.grid=list(mid=rep(V,length.out=length(C))))
}

\arguments{
  \item{C }{tracer concentration, defined at the centre of the grid cells.
    A vector of length N [M/L3].
  }
  \item{C.up }{tracer concentration at the upstream interface. One value [M/L3].
  }
  \item{C.down }{tracer concentration at downstream interface. One value [M/L3].
  }
  \item{C.lat }{tracer concentration in the lateral input, defined at
    grid cell centres. One value or a vector of length N [M/L3].
  }
  \item{C.lat.grid }{tracer concentration in lateral input defined as
    grid property.
  }
  \item{F.up }{total tracer input at the upstream interface. One value [M/T].
  }
  \item{F.down }{total tracer input at downstream interface. One value [M/T].
  }
  \item{F.lat }{total lateral tracer input, defined at grid cell centres.
    One value or a vector of length N [M/T].
  }
  \item{F.lat.grid }{total lateral tracer input defined as grid property.
  }
  \item{Disp }{BULK dispersion coefficient, defined on grid cell interfaces.
    One value or a vector of length N+1, [L3/T].
  }
  \item{Disp.grid }{BULK dispersion coefficient defined as grid property.
  }
  \item{flow }{water flow rate, defined on grid cell interfaces. One value
    or a vector of length N+1, [L3/T].
  }
  \item{flow.grid }{water flow rate defined as grid property.
  }
  \item{flow.lat }{lateral water flow rate into each volume box, defined at
    grid cell centres. One value or a vector of length N, [L3/T].
  }
  \item{flow.lat.grid }{lateral water flow rate defined as grid property.
  }
  \item{V }{grid cell volume, defined at grid cell centres. One value or a
    vector of length N [L3].
  }
  \item{V.grid }{grid cell volume defined as grid property.
  }
}

\value{
  \item{dC }{the rate of change of the concentration C due to transport,
    defined in the centre of each grid cell [M/L3/T].
  }
  \item{F }{mass flow across at the interface of each grid cell. A vector
    of length N+1 [M/T].
  }
  \item{F.up }{mass flow across the upstream boundary, positive = INTO
    model domain. One value [M/T].
  }
  \item{F.down }{mass flow across the downstream boundary, positive = OUT
    of model domain. One value [M/T].
  }
  \item{F.lat }{lateral mass input per volume box, positive = INTO model
    domain. A vector of length N [M/T].
  }
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

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
tran <- tran.volume.1D(C=OC,F.up=F.OC,F.lat=F.lat,Disp=Disp,
flow=flow,V=Volume)$dC 
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
F.lat <- F.lat.0*dnorm(x=Distance/lengthEstuary,
mean = Distance[nbox/2]/lengthEstuary, sd = 1/20, log = FALSE)/nbox 
Conc2 <- steady.band(runif(nbox),fun=river.model,nspec=1)$y   

#====================#
# Plotting output    #
#====================#

matplot(Distance/1000,cbind(Conc1,Conc2),lwd=2,
main="Organic carbon decay in an estuary",xlab="distance [km]",
ylab="OC Concentration [mM]",
type="l")  
legend ("topright",lty=1,col=c("black","red"),
        c("baseline","with lateral input"))
}
\references{
  Soetaert and Herman (2009) A practical guide to ecological modelling -
  using R as a simulation platform. Springer.
}

\details{
  The \bold{boundary conditions} are of type
  \itemize{
    \item 1. zero-gradient (default)
    \item 2. fixed concentration
    \item 3. fixed input
  }
  The \emph{bulk dispersion coefficient} (\code{Disp}) and the \emph{flow rate}
  (\code{flow}) can be either one value or a vector of length N+1, defined at
  all grid cell interfaces, including upstream and downstream boundary.

  The spatial discretisation is given by the volume of each box (\code{V}),
  which can be one value or a vector of length N+1, defined at the centre of
  each grid cell.

  The water flow is mass conservative. Over each volume box, the routine
  calculates internally the downstream outflow of water in terms of the
  upstream inflow and the lateral inflow.
}

\seealso{
  \code{\link{tran.1D}}
}

\keyword{utilities}

