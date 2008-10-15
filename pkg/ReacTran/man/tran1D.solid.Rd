\name{tran1D.solid}
\alias{tran1D.solid}
\title{Transport of solid substances in a porous medium}
\description{Estimates the rate of change of solid (particulate) substances due to 1-D diffusive and advective transport in porous media}
\usage{tran1D.solid(y, y.up=y[1],y.down=y[length(y)], 
flux.up=NA, flux.down=NA, disp, v=0, v.up=0, weight.up=1, 
por.INF, por.int=por.INF, por.mid=por.INF,
porgrid = list(INF=por.INF,int=por.int,mid=por.mid), 
dx=NULL, dx.int=dx, grid=NULL)}
\arguments{
  \item{y }{concentration per unit of solid in the porous medium, defined in centre of boxes. A vector of length N, [Mass/Length3], e.g. [nmol/cm3 solid]}
  \item{y.up }{concentration at upstream boundary interface. One value, [M/L3 solid]}
  \item{y.down }{concentration at downstream boundary interface. One value, [M/L3 solid]}
  \item{flux.up }{flux across the upper boundary interface, positive = IN medium. One value, [M/L2/Time]}
  \item{flux.down }{flux across the lower boundary interface, positive = OUT of medium. One value, [M/L2/T]}
  \item{disp }{dispersion coefficients (e.g. sediment turbation rate), defined on box interfaces. One value or a vector of length N+1, [L2/T]}
  \item{v }{velocity in the axis direction (e.g. sedimentation rate), defined on box interfaces. One value or a vector of length N+1, [L/T]}
  \item{v.up }{upward velocity, against the axis direction, defined on box interfaces. One value or a vector of length N+1, [L/T]}
  \item{weight.up }{upstream weighing coefficient for advection and upward flow, defined on box interfaces; default is backward differences. One value or a vector of length N+1, [-]}
  \item{por.INF }{(volumetric) porosity at infinite depth. One value, [-]}
  \item{por.int }{(volumetric) porosity at box interfaces. One value or a vector of length N+1, [-]}
  \item{por.mid }{(volumetric) porosity in middle (centre) of boxes. One value or a vector of lenght N, [-]}
  \item{porgrid }{porosity grid, a list containing the elements INF, int and mid, e.g. as calculated by setup.prop, [-]}
  \item{dx }{thickness of boxes. One value or a vector of length N, [L]}
  \item{dx.int }{distance over the interfaces, i.e. from centre to centre of boxes. One value or a vector of length N+1, [L]}
  \item{grid }{discretisation grid, a list containing at least elements dx and dx.int, e.g. as calculated by setup.grid, [L]}
}
\value{
  a list containing:
  \item{dy }{Rate of change of y in each layer due to transport, [M/L3 solid/T]}
  \item{flux }{Fluxes across each box interface. A vector of length N+1, [M/L2 BULK/Time]}
  \item{flux.up }{Flux across the upper boundary interface, positive = IN medium. One value, [M/L2 BULK/Time]}
  \item{flux.down }{Flux across the lower boundary interface, positive= OUT of medium. One value, [M/L2 BULK/T]}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{##################################################################
######  EXAMPLE 1: Boundary conditions                      ######
##################################################################

#====================#
# Model equations    #
#====================#
# model of 1-D transport and 1-st order decay
# upper boundary conc imposed, zero-gradient lower boundary

model <- function (t,Conc,pars,w)   
   return (list(tran1D.solid(Conc,flux.up=10,disp=disp,v=w,
                        por.INF =poro ,dx=dx)$dy-Conc*rate))

#====================#
# Model application  #
#====================#
	
y.down <- y.up <- 100       # mmol/m3
poro     <- 0.9             # constant porosity
disp     <- 0.01            # cm2/day
w        <- 0.01            # cm/day
dx       <- 0.01            # 0.01 cm thick slices, constant thickness
N        <- 1000            # 1000 layers
rate     <- 0.1             #/day 1st-order consumption rate
Depth    <- seq(from=dx/2,by=dx,len=N) 


# Steady-state solution, different values of advection rate
Conc1    <- steady.band(y=runif(N),func=model,nspec=1,w=0.01) 
Conc2    <- steady.band(y=runif(N),func=model,nspec=1,w=0.001)
Conc3    <- steady.band(y=runif(N),func=model,nspec=1,w=0.1)  

#====================#
# Plotting output    #
#====================#

matplot(cbind(Conc1$y,Conc2$y,Conc3$y),Depth,ylim=rev(range(Depth)),
        ylab="depth, cm",xlab="Concentration",type="l",lwd=2,
        main="tran1D.solid")  
legend ("right",c("w=0.01","w=0.001","w=0.1"),
 title="sedimentation rate",lty=1:3,col=1:3,lwd=2)
  
##################################################################
######  EXAMPLE 2: OMEXDIA: C, N, O2 diagenesis             ######
##################################################################

#====================#
# Model equations    #
#====================#

OMEXDIAmod <-function(time=0,     # time, not used here
                      state,      # FDET, SDET, O2, NO3, NH3, ODU
                      parms,      # parameter values
                      Full=FALSE) # if true: returns output vars
{
with (as.list(pars),{

Flux  <- MeanFlux * (1+sin(2*pi*time/365))

FDET  <- state[1:N]
SDET  <- state[(N+1)  :(2*N)]
O2    <- state[(2*N+1):(3*N)]
NO3   <- state[(3*N+1):(4*N)]
NH3   <- state[(4*N+1):(5*N)]
ODU   <- state[(5*N+1):(6*N)]

# Rate of change due to transport
# Solid substances
FDETtran<- tran1D.solid(FDET,flux.up=Flux*pFast   ,disp=Db,v=w,  
                     weight.up=0.5,porgrid=porGrid, grid=Grid)
                        
SDETtran<- tran1D.solid(SDET,flux.up=Flux*(1-pFast),disp=Db,v=w,  
                     weight.up=0.5,porgrid=porGrid, grid=Grid)

# Solute substances
O2tran  <- tran1D.solute(y=O2, y.up=bwO2, disp=DispO2, v=w,
                     porgrid=porGrid, grid=Grid)

NO3tran <- tran1D.solute(y=NO3, y.up=bwNO3,disp=DispNO3,v=w,
                     porgrid=porGrid, grid=Grid)

NH3tran <- tran1D.solute(y=NH3,y.up=bwNH3,disp=DispNH3/(1+NH3Ads),
                     v=w,porgrid=porGrid, grid=Grid)

ODUtran <- tran1D.solute(y=ODU,y.up=bwODU,disp=DispODU,v=w,
                     porgrid=porGrid, grid=Grid)

# production of DIC and DIN, expressed per cm3 LIQUID/day
DICprod_Min <- (rFast*FDET         +rSlow*SDET       )*
               (1.-porGrid$mid)/porGrid$mid
DINprod_Min <- (rFast*FDET*NCrFdet+rSlow*SDET*NCrSdet)*
               (1.-porGrid$mid)/porGrid$mid

# oxic mineralisation, denitrification, anoxic mineralisation
Oxicminlim <- O2/(O2+ksO2oxic)               # limitation terms
Denitrilim <- (1-O2/(O2+kinO2denit))*NO3/(NO3+ksNO3denit)
Anoxiclim  <- (1-O2/(O2+kinO2anox))*(1-NO3/(NO3+kinNO3anox))
Rescale    <- 1/(Oxicminlim+Denitrilim+Anoxiclim)

OxicMin    <- DICprod_Min*Oxicminlim*Rescale  # oxic mineralisation
Denitrific <- DICprod_Min*Denitrilim*Rescale  # Denitrification
AnoxicMin  <- DICprod_Min*Anoxiclim *Rescale  # an0xic mineralisation

# reoxicidation and ODU deposition
Nitri      <- rnit  *NH3*O2/(O2+ksO2nitri)
OduOx      <- rODUox*ODU*O2/(O2+ksO2oduox)

pDepo      <- min(1,0.233*(w*365)^0.336 )
OduDepo    <- AnoxicMin*pDepo 

# Update the rate of change
dFDET <- FDETtran$dy - rFast*FDET
dSDET <- SDETtran$dy - rSlow*SDET
dO2   <- O2tran$dy  -  OxicMin      -2* Nitri -      OduOx
dNH3  <- NH3tran$dy + (DINprod_Min  - Nitri) / (1.+NH3Ads)
dNO3  <- NO3tran$dy - 0.8*Denitrific + Nitri 
dODU  <- ODUtran$dy + AnoxicMin  - OduOx - OduDepo

# 
if (!Full)return(list(c(dFDET,dSDET,dO2,dNO3,dNH3,dODU))) else 

return(list(c(dFDET,dSDET,dO2,dNO3,dNH3,dODU), 
            c(O2flux=O2tran$flux.up,        #O2 sediment-water flux
             O2deepflux=O2tran$flux.down,   #O2 deep(burial) flux
             NO3flux=NO3tran$flux.up,       #NO3 sediment-water flux
             NO3deepflux=NO3tran$flux.down, #NO3 deep(burial) flux
             NH3flux=NH3tran$flux.up,       #NH3 sediment-water flux 
             NH3deepflux=NH3tran$flux.down, #NH3 deep (burial) flux
             ODUflux=ODUtran$flux.up,       #ODU sediment-water flux
             ODUdeepflux=ODUtran$flux.down, #ODU deep(burial) flux
             OxicMin=OxicMin,               #oxic mineralisation 
             Denitrific=Denitrific,         #denitrification rates
             AnoxicMin=AnoxicMin,           #anoxic mineralisation 
             Nitri=Nitri,                   #nitrification rates
             OduOx=OduOx)))                 #ODU oxidation rates

})

}



#====================#
# Model applications #
#====================#

# 50 layers; total length=15 cm, first box=0.1 cm
Grid  <- setup.grid(N=50,dx.1=0.1,len=15,type="exp")
Depth <- Grid$x                  # depth of each box
N     <- Grid$N
 
# porosity gradients
porGrid <- setup.prop(y.0=0.9,y.INF=0.7,y.coeff=2,grid=Grid)

# Bioturbation profile (at box interfaces)
biot <- 1/365          # cm2/d      - bioturbation coefficient
mixL <- 5              # cm         - depth of mixed layer
Db   <- setup.prop(y.0=biot,y.INF=0,y.coeff=1,L=mixL,grid=Grid)$int

# organic matter dynamics  #

MeanFlux <- 20000/12*100/365 # nmol/cm2/d - C deposition: 20gC/m2/yr
rFast    <- 0.01             #/day        - decay rate fast decay det.
rSlow    <- 0.00001          #/day        - decay rate slow decay det.
pFast    <- 0.9              #-           - fraction fast det. in flux
w        <- 0.1/1000/365     # cm/d       - advection rate
NCrFdet  <- 0.16             # molN/molC  - NC ratio fast decay det.
NCrSdet  <- 0.13             # molN/molC  - NC ratio slow decay det.

# oxygen and DIN dynamics  #

# Nutrient bottom water conditions
bwO2            <- 300       #mmol/m3     Oxygen conc in bottom water
bwNO3           <- 10        #mmol/m3
bwNH3           <- 1         #mmol/m3
bwODU           <- 0         #mmol/m3

# Nutrient parameters
NH3Ads          <- 1.3    #-           Adsorption coeff ammonium
rnit            <- 20.    #/d          Max nitrification rate
ksO2nitri       <- 1.     #umolO2/m3   half-sat O2 in nitrification
rODUox          <- 20.    #/d          Max rate oxidation of ODU
ksO2oduox       <- 1.     #mmolO2/m3   half-sat O2 in oxidation of ODU
ksO2oxic        <- 3.     #mmolO2/m3   half-sat O2 in oxic mineralis
ksNO3denit      <- 30.    #mmolNO3/m3  half-sat NO3 in denitrif
kinO2denit      <- 1.     #mmolO2/m3   half-sat O2 inhib denitrif
kinNO3anox      <- 1.     #mmolNO3/m3  half-sat NO3 inhib anoxic min
kinO2anox       <- 1.     #mmolO2/m3   half-sat O2 inhib anoxic min

# Diffusion coefficients, temp = 10dgC
Temp            <- 10                        # temperature 
DispO2          <- 0.955    +Temp*0.0386
DispNO3         <- 0.844992 +Temp*0.0336
DispNH3         <- 0.84672  +Temp*0.0336
DispODU         <- 0.8424   +Temp*0.0242

# parameter vector 
pars<- c(
MeanFlux  =  MeanFlux   ,rFast     =  rFast      ,
rSlow     =  rSlow      ,pFast     =  pFast      ,
w         =  w          ,NCrFdet   =  NCrFdet    ,
NCrSdet   =  NCrSdet    ,bwO2      =  bwO2       ,
bwNO3     =  bwNO3      ,bwNH3     =  bwNH3      ,
bwODU     =  bwODU      ,NH3Ads    =  NH3Ads     ,
rnit      =  rnit       ,ksO2nitri =  ksO2nitri  ,
rODUox    =  rODUox     ,ksO2oduox =  ksO2oduox  ,
ksO2oxic  =  ksO2oxic   ,ksNO3denit=  ksNO3denit ,
kinO2denit=  kinO2denit ,kinNO3anox=  kinNO3anox ,
kinO2anox =  kinO2anox  ,DispO2    =  DispO2     ,
DispNO3   =  DispNO3    ,DispNH3   =  DispNH3    ,
DispODU   =  DispODU  )

# STEADY-STATE CALCULATIONS USING R-function
#-----------------------------------------------
# Three runs, with different flux
pars["MeanFlux"] <- 15000/12*100/365  # nmol/cm2/day - C deposition: 15gC/m2/yr
OC   <- rep(10,6*N)
DIA  <- steady.1D (y=OC, func=OMEXDIAmod,Full=FALSE,
                   parms=pars, maxiter=100,
                   atol=1e-8,nspec=6,positive=TRUE)
CONC <- DIA$y

pars["MeanFlux"] <- 50000/12*100/365  # nmol/cm2/day - C deposition: 50gC/m2/yr
OC   <- rep(10,6*N)
DIA  <- steady.1D (y=OC, func=OMEXDIAmod,Full=FALSE,
                   parms=pars, maxiter=100,
                   atol=1e-8,nspec=6,positive=TRUE)
CONC <- cbind(CONC,DIA$y)

pars["MeanFlux"] <- 2000/12*100/365  # nmol/cm2/day - C deposition: 2gC/m2/yr
OC   <- rep(10,6*N)
DIA  <- steady.1D (y=OC, func=OMEXDIAmod,Full=FALSE,
                   parms=pars, maxiter=100,
                   atol=1e-8,nspec=6,positive=TRUE)
CONC <- cbind(CONC,DIA$y)



# Concentrations in sediment
FDET  <- CONC[1:N,]
SDET  <- CONC[(N+1)  :(2*N),]
O2    <- CONC[(2*N+1):(3*N),]
NO3   <- CONC[(3*N+1):(4*N),]
NH3   <- CONC[(4*N+1):(5*N),]
ODU   <- CONC[(5*N+1):(6*N),]

TOC  <- (FDET+SDET)*1200/10^9/2.5     # % organic carbon (excess)

par(mfrow=c(2,2))
matplot(TOC,Depth,ylim=c(15,0),xlab="procent" ,main="TOC",
        type="l",lwd=2)
matplot(O2,Depth,ylim=c(15,0),xlab="mmol/m3" ,main="O2",
        type="l",lwd=2)
matplot(NO3,Depth,ylim=c(15,0),xlab="mmol/m3" ,main="NO3",
        type="l",lwd=2)
matplot(NH3,Depth,ylim=c(15,0),xlab="mmol/m3" ,main="NH3",
        type="l",lwd=2)

legend ("bottom",col=1:3,lty=1:3,lwd=2,
legend=c("15gC/m2/yr","50gC/m2/yr","2gC/m2/yr"),title="flux")

mtext(outer=TRUE,side=3,line=-2,cex=1.5,"OMEXDIAmodel")

#---------------------------------------
# FOR an application of a DYNAMIC RUN
# SEE the example of OMEXDIAmodel
# (a Fortran implementation of OMEXDIA
#---------------------------------------}
\references{Berner R.A., 1980. Early Diagenesis- A Theoretical Approach. Princeton Univ. Press
  \cr Boudreau, B.P., 1997. Diagenetic Models and their Implementation. 
  Modelling transport and Reactions in Aquatic Sediments. Springer, Berlin, 414p.
  \cr Soetaert and Herman, a guide to ecological modelling - using R as a simulation platform, 2008. Springer.
  
  Reference for the OMEXDIA model (example 2):
  \cr Soetaert K, PMJ Herman and JJ Middelburg, 1996a. 
A model of early diagenetic processes from the shelf to abyssal depths. 
Geochimica Cosmochimica Acta, 60(6):1019-1040.
\cr Soetaert K, PMJ Herman and JJ Middelburg, 1996b.  
Dynamic response of deep-sea sediments to seasonal variation: a model.  
Limnol. Oceanogr. 41(8): 1651-1668.}
\details{For a specification of the \emph{boundary conditions} and \emph{spatial discretisation} see details of function tran1D
  
  The \emph{bioturbation coefficient} (disp), and \emph{sedimentation rate} (v) can be either one value or a vector.
  \cr If they are a vector, they must be of length N+1, defined at all box interfaces, including upper and lower boundary.

  \bold{Porosity} needs to be defined:
  \item 1. at the box interfaces (por.int), either one value, or a vector of length N+1
  \item 2. at the centre of boxes (por.mid), either one value or a vector of length N
  \item 3. and at infinite depth (por.INF ).
  
  Alternatively, porosity can be inputted as a list, e.g. created by \code{\link{setup.prop}}

  If por.int and por.mid are one value, they should be equal to por.INF .
  \cr No attempt is done to test whether the inputted porosities are compatible.}
\seealso{\code{\link{tran1D.solute}}, \code{\link{tran1D.volume}}, \code{\link{tran1D}}}
\note{This is just a particular application of the general transport routine tran1D, where
*1-porosity* provides the *surface area* and where the advection rate is corrected for the porosity gradient 
(to account for steady-state compaction).}
\keyword{utilities}

