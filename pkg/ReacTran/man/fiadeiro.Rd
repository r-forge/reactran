\name{fiadeiro}
\alias{fiadeiro}
\title{Upstream weighing coefficients for advection}
\description{Upstream weighing coefficients for advection that reduce numerical dispersion whilst conserving positivity}
\usage{fiadeiro(v, disp, dx.int, grid=list(dx.int=dx.int))}
\arguments{
  \item{v }{velocity, advection rate, [L/t], either one value or a vector of length N+1, with N the number of boxes}
  \item{disp }{diffusion rate, [L2/t], either one value or a vector of length N+1}
  \item{dx.int }{distances, [L], between centre of boxes (including two interfaces), either one value or a vector of length N+1 }
  \item{grid }{discretisation grid, e.g. as calculated by setup.grid}
}
\value{
  the weighing coefficient(s) for the upstream compartments, either one value or a vector of length N+1
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{##################################################################
######   EXAMPLE : advection weighing schemes               ######
##################################################################

#====================#
# Model equations    #
#====================#
	
model <- function (time,OC,pars,decay=0.05,weight.up=1)
{
 # model describing the concentration of particulate organic C (OC) 
 # that sink out of the water and decay
 
 dOC <- tran1D(OC,flux.up=10,disp=disp,v =sink,
                weight.up=weight.up,dx=dx)$dy - decay*OC
 return(list(dOC))                
}

#====================#
# Model application  #
#====================#

dx    <- 100                  # thickness of boxes
Dist  <- seq(0,1000,by=dx)    # water depth at modeled box interface
Depth <- seq(dx/2,1000,by=dx) # water depth at centre of modeled box
N     <- length(Depth)
cc    <- -0.005

# exponentially declining sinking rate, maximal 100 m/day
sink <- 100*exp(cc*Dist)
disp <- 1000 
Weights <-fiadeiro(sink,disp,dx)


# plot fiadeiro weighting coefficients as a function of sinking rate
plot(sink,Weights,type="b",main="Fiadeiro and Veronis scheme",
     xlab="sinking rate,m/day",
     ylab="upstream weighing coefficient")

# estimate the steady-state solution using a coarse grid (a0 boxes) 
# and based on: 
# (1) backward differences (2) centered diff, (3) fiadeiro scheme
ss   <- steady.band(y=runif(N),func=model,nspec=1)$y
ss   <- cbind(ss,steady.band(y=runif(N),func=model,
                             weight.up=0.5,nspec=1)$y)
ss   <- cbind(ss,steady.band(y=runif(N),func=model,
                             weight.up=Weights,nspec=1)$y)

#====================#
# Plotting output    #
#====================#

matplot(ss,Depth,pch=16,type="b",main="numerical diffusion",
      ylab="water depth, m",
      xlab="concentration of sinking particles",
      ylim=c(1000,0))

# now with increased resolution (1000 boxes)
dx    <- 0.1                   # thickness of boxes only 0.1 m
Dist  <- seq(0,1000,by=dx)    
Depth <- seq(dx/2,1000,by=dx) 
N     <- length(Depth)

# exponentially declining sinking rate, maximal 100 m/day
sink <- 100*exp(cc*Dist)
  
# with this resolution, the weighing coeff does not matter
st   <- steady.band(y=runif(N),fun=model,weight.up=1,nspec=1)$y  
lines(st,Depth,col="darkblue",lty=1)
legend("bottomright",title="numerical diffusion",
      legend=c("upstream","centered","fiadeiro","large resolution"),
      col=c(1:3,"darkblue"),lty=c(1:3,1),pch=c(16,16,16,NA))}
\references{Fiadeiro Me, Veronis G, 1977
\cr Weighted-mean schemes for finite-difference approximation to advection-diffusion equation.
Tellus v 29, pp 512-522}
\details{The Fiadeiro and Veronis (1977) weighing scheme reduces numerical diffusion, due to advection. 
  It is based on the following rationale:
\item When weighing is by backward differences (weighing coefficient of the upstream box = 1) state variables remain positive, 
but this scheme has the highest numerical (= artificial) diffusion.
\item When weighing is by centered differences (weighing coefficient = 0.5), the numerical diffusion is lower, 
but state variables may become negative.

In the Fiadeiro and Veronis (1977) weighing scheme, the weighing coefficients of the upstream boxes depend on 
the magnitude of the additional true diffusion (D): 
\cr the higher the true diffusion, the closer the weighing coefficients are to centered weighing...
\cr The weighing coefficients thus vary from 0.5 (very high true diffusion) to 1 (absence of true diffusion)

Note: if the substance concentrations decline in the direction of the axis, 
then centered differences will never lead to negative concentration. 
Thus, under these circumstances, centered differences are to be preferred over the fiadeiro scheme.}
\note{A certain amount of numerical diffusion always remains. 
  \cr Reducing grid size may be a more effective way of reducing numerical dispersion.}
\keyword{utilities}

