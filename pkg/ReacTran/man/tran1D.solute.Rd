\name{tran1D.solute}
\alias{tran1D.solute}
\title{Transport of dissolved substances in a porous medium}
\description{Estimates the rate of change of dissolved substances due to 1-D diffusive and advective transport in porous media}
\usage{tran1D.solute(y, y.up, y.down=y[length(y)], flux.up=NA, 
flux.down=NA, disp, v=0, v.up=0, weight.up=1, 
por.INF , por.int=por.INF , por.mid=por.INF ,
porgrid = list(INF=por.INF,int=por.int,mid=por.mid), 
tortuosity=1,dx=NULL, dx.int=dx, grid=NULL)}
\arguments{
  \item{y }{concentration per unit of solute in the porous medium, defined in centre of boxes. A vector of length N, [Mass/Length3], e.g. [nmol/cm3 liquid]}
  \item{y.up }{concentration at upstream boundary interface. One value, [M/L3 liquid]}
  \item{y.down }{concentration at downstream boundary interface. One value, [M/L3 liquid]}
  \item{flux.up }{flux across the upper boundary interface, positive = IN medium. One value, [M/L2/Time]}
  \item{flux.down }{flux across the lower boundary interface, positive= OUT of medium. One value, [M/L2/T]}
  \item{disp }{molecular diffusion coefficients, defined on box interfaces. One value or a vector of length N+1, [L2/T]}
  \item{v }{velocity in the axis direction (e.g. sedimentation rate), defined on box interfaces. One value or a vector of length N+1, [L/T]}
  \item{v.up }{upwelling rate, against the axis direction, defined on box interfaces. One value or a vector of length N+1, [L/T]}
  \item{weight.up }{upstream weighing coefficient for advection and upward flow, defined on box interfaces; default is backward differences. One value or a vector of length N+1, [-]}
  \item{por.INF }{(volumetric) porosity at Infinite depth. One value, [-]}
  \item{por.int }{(volumetric) porosity at box interfaces. One value or a vector of length N+1, [-]}
  \item{por.mid }{(volumetric) porosity in centre of boxes. One value or a vector of lenght N, [-]}
  \item{porgrid }{porosity grid, a list containing the elements INF, int and mid, e.g. as calculated by setup.prop, [-]}
  \item{tortuosity }{sediment tortuosity coefficient}
  \item{dx }{thickness of boxes. One value or a vector of length N, [L]}
  \item{dx.int }{distance over the interfaces, i.e. from centre to centre of boxes. One value or a vector of length N+1, [L]}
  \item{grid }{discretisation grid, a list containing at least elements dx and dx.int, e.g. as calculated by setup.grid, [L]}
}
\value{
  a list containing:
  \item{dy }{rate of change of y in each layer due to transport, [M/L3 liquid/T]}
  \item{flux }{Fluxes across each box interface. A vector of length N+1, [M/L2/Time]}
  \item{flux.up }{Flux across the upper boundary interface, positive = IN medium. One value, [M/L2 BULK/Time]}
  \item{flux.down }{Flux across the lower boundary interface, positive= OUT of medium. One value, [M/L2 BULK/T]}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{##################################################################
######  EXAMPLE 1: Different boundary conditions            ######
##################################################################

#====================#
# Model equations    #
#====================#
# model of transport and first-order decay
# MODEL 1: two imposed boundary concentrations
model1 <- function (t,Conc,pars)   
   return (list(tran1D.solute(Conc,y.up,y.down,disp=disp,por.INF =poro ,
                         dx=dx)$dy-Conc*rate))

# MODEL 2: upper boundary conc imposed, zero-gradient lower boundary
model2 <- function (t,Conc,pars)   
   return (list(tran1D.solute(Conc,y.up,disp=disp,por.INF =poro ,
                         dx=dx)$dy-Conc*rate))
#====================#
# Model application  #
#====================#

y.down <- y.up <- 100       # mmol/m3
poro     <- 0.9             # constant porosity
disp     <- 1               # cm2/day
dx       <- 0.01            # 0.01 cm thick slices, constant thickness
N        <- 1000            # 1000 layers
rate     <- 0.1             #/day 1st-order consumption rate
Depth    <- seq(from=dx/2,by=dx,len=N) 


# upper boundary conc imposed, zero-gradient lower boundary
model3 <- function (t,Conc,pars)   
   return (list(tran1D.solute(Conc,flux.up=10,disp=disp,
                     por.INF =poro,dx=dx)$dy-Conc*rate))

# steady-state solution of the different models
Conc1     <- steady.band(runif(N),func=model1,nspec=1)  
Conc2     <- steady.band(runif(N),func=model2,nspec=1)  
Conc3     <- steady.band(runif(N),func=model3,nspec=1)  

#====================#
# model output       #
#====================#

matplot(cbind(Conc1$y,Conc2$y,Conc3$y),Depth,ylab="depth, cm",
        ylim=rev(range(Depth)),xlab="Concentration",type="l",
        lwd=2,main="tran1D.solute")  
legend ("right", title="Boundaries",lty=1:3,col=1:3,lwd=2,
        c("upper and lower: conc","upper: conc, lower: 0-gradient",
          "upper:flux, lower: 0-gradient"))


##################################################################
######  EXAMPLE 2: OXYGEN diagenesis                        ######
##################################################################

#====================#
# Model equations    #
#====================#

O2DIA <- function (time=0,      # time, not used here
                   O2,          # oxygen concentration
                   parms=NULL,  # parameter values; not used
                   k=0,         # 1st-order O2 consumption rate
                   Q=0)         # 0-order consumption rate 

#===============================================
# Oxygen is transported and consumed with 
# 1-st order rate k and 0-th order rate Q
#===============================================

{
# Rate of change due to transport
O2tran  <- tran1D.solute(O2,y.up=bwO2,disp=dispO2,
                        porgrid=porGrid,grid=Grid)
return (list(O2tran$dy - k*O2 - Q)           )
}

#====================#
# Model application  #
#====================#
# parameters and grid
# model grid includes diffusive boundary layer
DBL      <- 0.2 
Grid     <- setup.grid(dx.1=0.05,len=10,
                       dx.int.up = DBL+0.05/2)

# sediment parameters
N        <- Grid$N
Depth    <- Grid$x  # depth at middle of each layer

# porosity gradient
porGrid  <- setup.prop(y.0=0.8,y.INF=0.6,y.coeff=2,grid=Grid)

# porosity in diffusive boundary layer = 1
porGrid$int[1] <- 1                         

# biogeochemical parameters
dispO2   <- 0.4   # diffusion coefficient, cm2/d
bwO2     <- 250   # bottom water oxygen concentration, mmol/m3

# 4 model applications: different consumption rate 
# steady-state solution                          
# constant 1-st order decay   
k    <- rep(0.3,N)       

# initial guess: N random numbers between 0,1
sol  <- steady.band (runif(N), fun=O2DIA, pos=TRUE, k=k,nspec=1)
O2   <- c(bwO2,sol$y)    # add bwO2 concentration to result

# hyperbolic declining 1-st order decay
k2   <- 5*(1-Depth/(Depth+0.1)) 
sol  <- steady.band (runif(N), fun=O2DIA, pos=TRUE, k=k2,nspec=1)
O2   <- cbind(O2,c(bwO2,sol$y))

# hyperbolic increasing 1-st order decay
k3   <- 1*(Depth/(Depth+2.5))   # hyperbolic increase
sol  <- steady.band (runif(N), fun=O2DIA, pos=TRUE, k=k3,nspec=1)
O2   <- cbind(O2,c(bwO2,sol$y))

O2cons <- O2[-1,] * cbind(k,k2,k3)  # true o2-consumption rate

# 0-order decay
Q     <- 10
sol   <- steady.band (runif(N), fun=O2DIA, k=0,Q=Q,nspec=1)
O2    <- cbind(O2,c(bwO2,pmax(0,sol$y)))  # negative O2 removed
O2cons<- cbind(O2cons,rep(Q,N))

#====================#
# model output       #
#====================#

pm   <- par(mfrow = c(1,2))
matplot(O2,c(-DBL,Depth),ylim=c(10,-DBL),main="O2",type="l",
xlab="mmolO2/m3 Liquid",ylab="depth, cm",lwd=c(1,2,1,1),col=1:4)
abline(h=0)

matplot(O2cons,Depth,ylim=c(10,-DBL),main="O2 consumption rate",
        type="l",log="x",xlab="mmolO2/m3/d",
        ylab="depth, cm",lwd=c(1,2,1,1),col=1:4)
abline(h=0)

mtext(outer=TRUE,side=3,"O2 diagenesis")
par(mfrow=pm$mfrow)

##################################################################
######  EXAMPLE 3: SEEP diagenesis                          ######
##################################################################

#====================#
# Model equations    #
#====================#

SEEPDIAmodel <- function (time=0,      # time, not used here
                          state,       # concentrations: O2,HS,SO4,CH4
                          parms=NULL)  # parameter values; not used

#=====================================================================
# Estimates the rate of change of Oxygen, hydrogen sulphide, 
# sulphate, and methane under seepage conditions. 
# Reactions are anaerobic oxidation of methane (AOM) and 
# the reoxidation of sulphide (HSox)
#=====================================================================

{

 O2  <- state[1:N]
 HS  <- state[(N+1)  :(2*N)]
 SO4 <- state[(2*N+1):(3*N)]
 CH4 <- state[(3*N+1):(4*N)]

# Rate of change due to transport
 O2tran  <- tran1D.solute(O2,y.up=bwO2,y.down=deepO2,disp=dispO2,
                        v=w,v.up=upwel,porgrid=porGrid,grid=Grid)
                        
 HStran <- tran1D.solute(HS,y.up=bwHS,y.down=deepHS,disp=dispHS,
                        v=w,v.up=upwel,porgrid=porGrid,grid=Grid)

 SO4tran <- tran1D.solute(SO4,y.up=bwSO4,y.down=deepSO4,disp=dispSO4,
                        v=w,v.up=upwel,porgrid=porGrid,grid=Grid)
                        
 CH4tran <- tran1D.solute(CH4,y.up=bwCH4,y.down=deepCH4,disp=dispCH4,
                        v=w,v.up=upwel,porgrid=porGrid,grid=Grid)

# Anaerobic oxidation of methane
 AOM  <- rAOM * SO4 * CH4

# reoxidation of sulphide
 HSox <- rHSox * HS * O2

# Update the rate of change
#dSV  =  transport     + reaction
 dO2  <- O2tran$dy  -2*HSox
 dHS  <- HStran$dy  -  HSox + AOM
 dSO4 <- SO4tran$dy +  HSox - AOM 
 dCH4 <- CH4tran$dy         - AOM       

# returning values

return(list(c(dO2=dO2,dHS=dHS,dSO4=dSO4,dCH4=dCH4), # rates of change
          O2flux =O2tran$flux.up ,       #O2 sediment-water flux
          O2deepflux =O2tran$flux.down,  #O2 deep (burial) flux
          HSflux =HStran$flux.up ,       #sulphide sediment-water flux
          HSdeepflux =HStran$flux.down,  #sulphide deep(burial) flux
          SO4flux=SO4tran$flux.up,       #Sulphate sediment-water flux
          SO4deepflux=SO4tran$flux.down, #Sulphate deep (burial) flux
          CH4flux=CH4tran$flux.up,       #Methane sediment-water flux
          CH4deepflux=CH4tran$flux.down, #Methane deep(burial) flux
          AOM=AOM,          # profile of anoxic oxidation of methane
          HSox=HSox))       # profile of sulphide reoxidation rates

}


#====================#
# Model run          #
#====================#
# sediment parameters
Grid  <- setup.grid(dx.1=0.2,len=100)
Depth <- Grid$x   # depth of each box
N     <- Grid$N

# porosity gradient: surf por=0.9,deep por=0.7)
porGrid <- setup.prop(y.0=0.9,y.INF=0.7,y.coeff=2,
                    grid=Grid)

w       <- 0           # sediment advection rate, cm/hr
upwel   <- 10/365/24   # upwelling rate, cm/hr

# deep concentrations
deepO2  <-         0   # mmol/m3       
deepHS  <-         0   # mmol/m3       
deepSO4 <-         0   # mmol/m3      
deepCH4 <-     68000   # mmol/m3      

# bottom water concentrations 
bwO2    <-      180    # mmol/m3       
bwHS    <-         0   # mmol/m3        
bwSO4   <-     28900   # mmol/m3       
bwCH4   <-         0   # mmol/m3      

# diffusion coefficients   
dispCH4  <- 314/365/24 # cm2/hr  
dispSO4  <- 169/365/24 # cm2/hr  
dispO2   <- 265/365/24 # cm2/hr  
dispHS   <- 346/365/24 # cm2/hr  

# process rates
rAOM    <-  80e-6/365/24  #/(mmol/m3)/hr  rate anaerobic oxid methane
rHSox   <-  176e-6/365/24 #/(mmol/m3)/hr  rate reoxid hydrogen sulphate            

# Solve steady-state; initialise with random numbers between 0,1
# Three applications, upwelling rate varying from 10 -> 1000 cm/year

# require(inverse.R)
# scenario 1
upwel   <- 10/365/24    # cm/hr
CONC10  <- steady.1D (runif(4*N),fun=SEEPDIAmodel,atol=1e-8,nspec=4,
                      positive=TRUE)$y
AOM10   <- SEEPDIAmodel(state=CONC10)$AOM

# scenario 2
upwel   <- 100/365/24   # cm/hr
CONC100 <- steady.1D (CONC10,fun=SEEPDIAmodel,atol=1e-8,nspec=4,
                      positive=TRUE)$y
AOM100  <- SEEPDIAmodel(state=CONC100)$AOM

# scenario 3
upwel   <- 1000/365/24  # cm/hr
CONC1000<- steady.1D (CONC100,fun=SEEPDIAmodel,atol=1e-8,nspec=4,
                      positive=TRUE)$y
AOM1000 <- SEEPDIAmodel(state=CONC1000)$AOM

# all scenarios combined in one matrix
CONC <- cbind(CONC10,CONC100,CONC1000)
AOM  <- cbind(AOM10,AOM100,AOM1000)

# 3 columns for each substance
O2  <- CONC[1:N,]
HS  <- CONC[(N+1)  :(2*N),]
SO4 <- CONC[(2*N+1):(3*N),]
CH4 <- CONC[(3*N+1):(4*N),]


#====================#
# plotting output    #
#====================#

# layout: 1st row: 3 columns merged and small (0.2); 
#         2nd row: higher (0.8) with 3 columns, equally large (2)
nf <- layout(matrix(c(1,1,1,2,3,4), ncol=3, nrow=2, byrow = TRUE),
             width=c(2,2,2),height=c(.25,.75))

cextxt <- 1.2
cextxt2 <- 1.3

## first the triangle plot

plot(0,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,xlab="",ylab="")
polygon(x=c(0,1,1),y=c(0.5,1,0),col="darkblue")
text(1,0.5,"fluid flow",col="white",cex=1.8,adj=1)

mtext(side=1,outer=FALSE,line=3, at=0.3,"methane",cex=cextxt,
      col="darkblue")
mtext(side=1,outer=FALSE,line=5, at=0.3,"sulphate",cex=cextxt,
      col="darkgreen")
mtext(side=1,outer=FALSE,line=4, at=0.5,"mmol /m3",cex=cextxt)

upwel<-c("10 cm /yr","100 cm /yr","1000 cm /yr")

## then the three sediment depth profiles
for (i in 1:3)
{
plot (AOM[,i],Depth,ylim=c(100,0),col="red",type="l",lwd=3,xlab="",
          ylab="Sedimentdepth, cm")
legend("bottom",upwel[i])
par(new=TRUE)
plot(CH4[,i],Depth,ylim=c(100,0),xlim=c(0,70000),
          col="darkblue",type="l",lwd=3,axes=FALSE,xlab="",ylab="")
axis(side=3)
lines(SO4[,i],Depth,col="darkgreen",lwd=3)
}

mtext(side=1,outer=TRUE,line=-2, "mmol/m3/hr",cex=cextxt)
mtext(side=1,at=0.3,outer=TRUE,line=-2, "AOM",cex=cextxt,
      col="red",adj=0)
mtext(outer=TRUE,side=3,line=-2,cex=1.5,"SEEPDIAmodel")}
\references{Berner R.A., 1980. Early Diagenesis- A Theoretical Approach. Princeton Univ. Press
  \cr Boudreau, B.P., 1997. Diagenetic Models and their Implementation. 
  Modelling transport and Reactions in Aquatic Sediments. Springer, Berlin, 414p.
  \cr Soetaert and Herman, a guide to ecological modelling, Book in preparation.
  
  Reference for the SEEPDIA model (example 3):
  \cr van Oevelen, 2007. HERMES report Deliverable 47}
\details{For a specification of the \emph{boundary conditions} and \emph{spatial discretisation} see function tran1D

  See \code{\link{tran1D.solid}} for details about \emph{transport coefficients} and \emph{porosity}.}
\seealso{\code{\link{tran1D.solid}}, \code{\link{tran1D.volume}}, \code{\link{tran1D}}}
\note{This is just a particular application of the general transport routine tran1D, where
*porosity* provides the *surface area* and where the advection rate is corrected for porosity (to account for steady-state compaction).}
\keyword{utilities}

