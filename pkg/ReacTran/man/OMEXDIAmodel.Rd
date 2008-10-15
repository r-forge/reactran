\name{OMEXDIAmodel}
\alias{OMEXDIAmodel}
\title{application of the OMEXDIA diagenetic model - fortran DLL}
\description{A 1-D model of Carbon, nitrogen and oxygen diagenesis in a marine sediment.
\cr The model describes six state variables, in \bold{100} layers:
\item 2 fractions of organic carbon (FDET,SDET): fast and slow decaying, solid substance
\item Oxygen (O2), dissolved substance
\item Nitrate (NO3), dissolved substance
\item Ammonium (NH3), dissolved substance
\item Oxygen demand unit (ODU), dissolved substance, as lump sum of all reduced substances other than ammonium

See Soetaert et al., 1996 for further details of the model.}
\usage{OMEXDIAmodel(pars, grid, porgrid, Db, Dyna=TRUE)}
\arguments{
  \item{pars }{ The model parameters, a vector such as \code{\link{OMEXDIAparms}}}
  \item{grid }{ The discretisation grid, e.g. generated with \code{\link{setup.grid}}; should contain 100 boxes}
  \item{porgrid }{ The porosity grid, e.g. generated with \code{\link{setup.prop}}; should contain 100 boxes and be consistent with \code{grid}}
  \item{Db }{The bioturbation profile, defined on box interfaces, e.g. containing 101 values}
  \item{Dyna }{logical; if TRUE: also runs the model dynamically}
}
\value{
  a list containing:
  \item{steady }{The steady-state condition of the state variables, a vector containing steady-state concentrations of
FDET(0-100), SDET(101-200), O2 (201-300), NO3 (301-400), NH3 (401-500) and ODU (501-600)}
  \item{precis }{Only if \code{Dyna} is TRUE: the precision of the steady-state solution}
  \item{Solved }{Only if \code{Dyna} is TRUE: a logical, TRUE when steady-state has been reached }
  \item{times }{the time values for which the model output has been produced}
  \item{FDET }{The concentration profile (columns) of Fast decaying organic carbon, one for each time value (rows)}
  \item{SDET }{The concentration profile (columns) of slow decaying organic carbon, one for each time value (rows)}
  \item{O2 }{The concentration profile (columns) of oxygen, one for each time value (rows)}
  \item{NO3 }{The concentration profile (columns) of nitrate, one for each time value (rows)}
  \item{NH3 }{The concentration profile (columns) of ammonium, one for each time value (rows)}
  \item{ODU }{The concentration profile (columns) of ODU, one for each time value (rows)}
  \item{O2flux }{The sediment-water exhange rate of oxygen, one per time value}
  \item{NO3flux }{The sediment-water exhange rate of nitrate, one per time value}
  \item{NH3flux }{The sediment-water exhange rate of ammonium, one per time value}
  \item{ODUflux }{The sediment-water exhange rate of ODU, one per time value}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# 100 layers; total length=50 cm, first box=0.01 cm
N     <- 100
Grid  <- setup.grid(N=N,dx.1=0.01,len=50,type="exp")
Depth <- Grid$x                  # depth of each box

# porosity gradients
porGrid <- setup.prop(y.0=0.9,y.INF=0.7,y.coeff=2,grid=Grid)

# Bioturbation profile (at box interfaces)
biot <- 1/365          # cm2/d      - bioturbation coefficient
mixL <- 5              # cm         - depth of mixed layer
Db   <- setup.prop(y.0=biot,y.INF=0,y.coeff=1,L=mixL,grid=Grid)$int

#-------------------------------------------
# Three dynamic runs, with different flux
#-------------------------------------------

# Application 1:
pars <- OMEXDIAparms
pars["MeanFlux"] <- 50000/12*100/365  # nmol/cm2/day
out1 <- OMEXDIAmodel(pars,Grid,porGrid,Db)

pars["MeanFlux"] <- 15000/12*100/365  # nmol/cm2/day
out2 <- OMEXDIAmodel(pars,Grid,porGrid,Db)

pars["MeanFlux"] <- 2000/12*100/365  # nmol/cm2/day
out3 <- OMEXDIAmodel(pars,Grid,porGrid,Db)

# Steady-state concentrations in sediment
CONC  <- cbind(out1$steady,out2$steady,out3$steady)

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

#====================#
# Dynamic output     #
#====================#

par(mfrow=c(3,3))
par(oma = c(0,0,3,0))

femmecol<- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
times<-out1$times
image(x=times,y=Depth,z=out1$O2,col= femmecol(100),xlab="time, days",
        ylim=c(5,0),ylab= "Depth, cm",main="Oxygen, mmol/m3")
contour(x=times,y=Depth,z=out1$O2, add = TRUE)
image(x=times,y=Depth,z=out2$O2,col= femmecol(100),xlab="time, days",
        ylim=c(5,0),ylab= "Depth, cm",main="Oxygen, mmol/m3")
contour(x=times,y=Depth,z=out2$O2, add = TRUE)
image(x=times,y=Depth,z=out3$O2,col= femmecol(100),xlab="time, days",
        ylim=c(5,0),ylab= "Depth, cm",main="Oxygen, mmol/m3")
contour(x=times,y=Depth,z=out3$O2, add = TRUE)

image(x=times,y=Depth,z=out1$NO3,col= femmecol(100),xlab="time, days",
        ylim=c(15,0),ylab= "Depth, cm",main="Nitrate, mmol/m3")
contour(x=times,y=Depth,z=out1$NO3, add = TRUE)
image(x=times,y=Depth,z=out2$NO3,col= femmecol(100),xlab="time, days",
        ylim=c(15,0),ylab= "Depth, cm",main="Nitrate, mmol/m3")
contour(x=times,y=Depth,z=out2$NO3, add = TRUE)
image(x=times,y=Depth,z=out3$NO3,col= femmecol(100),xlab="time, days",
        ylim=c(15,0),ylab= "Depth, cm",main="Nitrate, mmol/m3")
contour(x=times,y=Depth,z=out3$NO3, add = TRUE)

image(x=times,y=Depth,z=out1$NH3,col= femmecol(100),xlab="time, days",
        ylim=c(15,0),ylab= "Depth, cm",main="Ammonium, mmol/m3")
contour(x=times,y=Depth,z=out1$NH3, add = TRUE)
image(x=times,y=Depth,z=out2$NH3,col= femmecol(100),xlab="time, days",
        ylim=c(15,0),ylab= "Depth, cm",main="Ammonium, mmol/m3")
contour(x=times,y=Depth,z=out2$NH3, add = TRUE)
image(x=times,y=Depth,z=out3$NH3,col= femmecol(100),xlab="time, days",
        ylim=c(15,0),ylab= "Depth, cm",main="Ammonium, mmol/m3")
contour(x=times,y=Depth,z=out3$NH3, add = TRUE)

mtext(outer=TRUE,side=3,line=-1,cex=1.15,at=c(0.2,0.5,0.85),
       c("High flux","Medium flux","low flux"))
mtext(outer=TRUE,side=3,line=1,cex=1.5,"OMEXDIAmodel")}
\references{Soetaert K, PMJ Herman and JJ Middelburg, 1996a.
A model of early diagenetic processes from the shelf to abyssal depths.
Geochimica Cosmochimica Acta, 60(6):1019-1040.

Soetaert K, PMJ Herman and JJ Middelburg, 1996b.
Dynamic response of deep-sea sediments to seasonal variation: a model.
Limnol. Oceanogr. 41(8): 1651-1668.}
\details{The model application starts by estimating the steady-state condition of the model,
after which a dynamic simulation, with sinusoidal forcing of organic carbon deposition flux is run.

First one year of spinup is executed, and the final condition used as initial value for the
subsequent run, again for one year and which is written to the output.

For efficiency reasons, the OMEXDIA diagenetic model was written in Fortran, and this code linked to
the package. The Fortran code can be found in subdirectory \file{dynload} of the packages subdirectory, (\file{omexdia.f}) and 
an R-file showing how to compile this file to a dynamically linked library (DLL), and how to use this in R is also included 
(\file{omexdia.r}).

To solve the model, a steady-state solver from package rootSolve (here we used \code{\link[rootSolve]{steady.band}}) and an integration routines 
in package deSolve (we use \code{\link[deSolve]{ode.band}}) are used. The routines in these packages have been devised for linking with compiled code (here fortran).

The execution speed of this version of the model is similar to the original application where all code was in Fortran; it is 
almost two orders of magnitude faster than the implementation in \R.

An implementation of the same model, in \R can be found as one of the examples of \code{\link{tran1D.solid}}}
\keyword{misc}

