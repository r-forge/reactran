\name{fiadeiro}
\alias{fiadeiro}
\title{Upstream weighing coefficients for advection}
\description{Upstream weighing coefficients used in the finite difference scheme for advection. These AFVW (advective finite difference weights) reduce numerical dispersion while conserving positivity.}
\usage{fiadeiro(v, D, dx.aux, grid=list(dx.aux=dx.aux))}
\arguments{
  \item{v }{advective velocity, either one value or a vector of length N+1 with N the number of grid cells [L/T]}
  \item{D }{diffusion coefficient, either one value or a vector of length N+1 [L2/T]}
  \item{dx.aux }{auxiliary vector containing the distances between the locations where the concentration is defined (i.e. the grid cell centers and the two outer interfaces), either one value or a vector of length N+1 }
  \item{grid }{discretization grid as calculated by setup.grid}
}
\value{
  the weighing coefficient(s) for the upstream compartments, either one value or a vector of length N+1
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>, Filip Meysman <f.meysman@nioo.knaw.nl>}
\examples{
#=============================================================================
# Model formulation (set of differential equations)
#=============================================================================

# Test model to evaluate different weighing coefficients for finite difference
# formulas and their effect on munerical diffusion. The model describes the
# decay of organic carbon (OC) as it settles through the ocean water column.

model <- function (time,OC,pars,AFDW=1)
{
dOC <- tran.1D(OC,flux.up=F_OC,D=D.eddy,v=v.sink,AFDW=AFDW,dx=dx)$dC - k*OC
return(list(dOC))
}
#=============================================================================
# Parameter set
#=============================================================================

L <- 1000         # water depth model domain [m]
x.att <- 200      # attenuation depth [m]
v.sink.0 <- 10   # sinking velocity at the surface [m d-1]
D.eddy <- 10       # eddy diffusion coefficient [m2 d-1]
F_OC <- 10        # particle flux [mol m-2 d-1]
k <- 0.1          # decay coefficient [d-1]

#=============================================================================
# Model solution for a coarse grid (10 grid cells)
#=============================================================================

# Setting up the grid
N <- 10                              # number of grid layers 
dx <- L/N                            # thickness of boxes[m]
dx.aux <- rep(dx,(N+1))              # auxilliary grid vector
x.int <- seq(from=0,to=L,by=dx)      # water depth at box interfaces [m]
x.mid <- seq(from=dx/2,to=L,by=dx)   # water depth at box centres [m]

# Exponetially declining sink velocity
v.sink <- v.sink.0*exp(-x.int/x.att)  # sink velocity [m d-1]
Pe <- v.sink*dx/D.eddy

# Calculate the weighing coefficients
AFDW <- fiadeiro(v=v.sink,D=D.eddy,dx.aux=dx.aux)

par(mfrow=c(2,1),cex.main=1.2,cex.lab=1.2)

matplot(Pe,x.int,log="x",pch=19,ylim=c(L,0),xlim=c(0.1,1000), 
xlab="",ylab="depth [m]",main=expression("Peclet number"),axes=FALSE)
abline(h = 0)
axis(pos=NA, side=2)
axis(pos=0, side=3)

matplot(AFDW,x.int,pch=19,ylim=c(L,0),xlim=c(0.5,1), 
xlab="",ylab="depth [m]",main=expression("AFDW coefficient"),axes=FALSE)
abline(h = 0)
axis(pos=NA, side=2)
axis(pos=0, side=3)

# Three steady-state solutions for a coarse grid based on:
# (1) backward differences (BD)
# (2) centered differences (CD)
# (3) Fiadeiro & Veronis scheme (FV)
BD <- steady.band(y=runif(N), func=model, AFDW=1.0, nspec=1)$y
CD <- steady.band(y=runif(N), func=model, AFDW=0.5, nspec=1)$y
FV <- steady.band(y=runif(N), func=model, AFDW=AFDW, nspec=1)$y
CONC <- cbind(BD,CD,FV)

windows()

par(mfrow=c(1,2))

# Plotting output
matplot(CONC,x.mid,pch=16,type="b",ylim=c(L,0),
xlab="",ylab="depth [m]",main=expression("conc (Low resolution grid)"),axes=FALSE)
abline(h = 0)
axis(pos=0, side=2)
axis(pos=0, side=3)
legend("bottomright",legend=c("backward diff","centred diff","Fiadeiro&Veronis"),
col=c(1:3),lty=c(1:3),pch=c(16,16,16))


#=============================================================================
# Model solution for a fine grid (1000 grid cells)
#=============================================================================

# Setting up the grid
N <- 1000                              # number of grid layers 
dx <- L/N                            # thickness of boxes[m]
dx.aux <- rep(dx,(N+1))              # auxilliary grid vector
x.int <- seq(from=0,to=L,by=dx)      # water depth at box interfaces [m]
x.mid <- seq(from=dx/2,to=L,by=dx)   # water depth at box centres [m]

# Exponetially declining sink velocity
v.sink <- v.sink.0*exp(-x.int/x.att)  # sink velocity [m d-1]
Pe <- v.sink*dx/D.eddy

# Calculate the weighing coefficients
AFDW <- fiadeiro(v=v.sink,D=D.eddy,dx.aux=dx.aux)

# Three steady-state solutions for a coarse grid based on:
# (1) backward differences (BD)
# (2) centered differences (CD)
# (3) Fiadeiro & Veronis scheme (FV)
BD <- steady.band(y=runif(N), func=model, AFDW=1.0, nspec=1)$y
CD <- steady.band(y=runif(N), func=model, AFDW=0.5, nspec=1)$y
FV <- steady.band(y=runif(N), func=model, AFDW=AFDW, nspec=1)$y
HR_CONC <- cbind(BD,CD,FV)

# Plotting output
matplot(HR_CONC,x.mid,pch=16,type="b",ylim=c(L,0),
xlab="",ylab="depth [m]",main=expression("conc (High resolution grid)"),axes=FALSE)
abline(h = 0)
axis(pos=0, side=2)
axis(pos=0, side=3)
legend("bottomright",legend=c("backward diff","centred diff","Fiadeiro&Veronis"),
col=c(1:3),lty=c(1:3),pch=c(16,16,16))

# Results and conclusions:
# - For this fine grid all three solutions are identical
# - For the coarse grid, the BD and FV solutions show numerical dispersion
}
\references{Fiadeiro Me, Veronis G, 1977
\cr Weighted-mean schemes for finite-difference approximation to advection-diffusion equation.
Tellus v 29, pp 512-522}
\details{The Fiadeiro and Veronis (1977) weighing scheme reduces numerical diffusion, due to advection. 
  It is based on the following rationale:
  \itemize{
\item When weighing is by backward differences (weighing coefficient of the upstream box = 1) state variables remain positive, 
but this scheme has the highest numerical (= artificial) diffusion.
\item When weighing is by centered differences (weighing coefficient = 0.5), the numerical diffusion is lower, 
but state variables may become negative.
}
In the Fiadeiro and Veronis (1977) weighing scheme, the weighing coefficients of the upstream boxes depend on 
the magnitude of the additional true diffusion (D): 
\cr the higher the true diffusion, the closer the weighing coefficients are to centered weighing...
\cr The weighing coefficients thus vary from 0.5 (very high true diffusion) to 1 (absence of true diffusion)

Note: if the substance concentrations decline in the direction of the axis, 
then centered differences will never lead to negative concentrations. 
Thus, under these circumstances, centered differences are to be preferred over the fiadeiro scheme.}
\note{A certain amount of numerical diffusion always remains. 
  \cr Reducing grid size may be a more effective way of reducing numerical dispersion.}
\keyword{utilities}

