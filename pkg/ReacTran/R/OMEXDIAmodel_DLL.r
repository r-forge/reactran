

##################################################################
######          OMEXDIA: C, N, O2 diagenesis                ######
##################################################################

#====================#
# Model equations    #
#====================#

OMEXDIAmodel <- function (pars,
                          grid,
                          porgrid,
                          Db,
                          Dyna = TRUE) 
{
N <- grid$N
if (length(Db) != 101) stop("Db should contain 101 values")
if (length(pars) != 25) stop("pars should contain 25 values")
if (grid$N != 100)stop("grid should contain 100 boxes")
if (length(porgrid$mid) != 100)stop("porosity grid should contain 100 boxes")
if (is.null(grid$dx))stop("'grid$dx' should be present")
if (is.null(grid$dx.int))stop("'grid$dx.int' should be present")
if (is.null(porgrid$mid))stop("'porgrid$mid' should be present")
if (is.null(porgrid$int))stop("'porgrid$int' should be present")

# First the steady-state condition
OC   <- rep(10,6*N)
DIA  <- steady.band(y=OC,fun="omexdiamod",initfun="initomexdia",
                   initpar=c(pars,grid$dx,grid$dx.int,
                   porgrid$mid,porgrid$int,Db),nspec=6,
                   dllname="ReacTran",nout=8,positive=TRUE)
steady <-DIA$y

# state variables rearranged from ordering per slice -> per spec
ii       <- as.vector(t(matrix(ncol=6,1:(6*N))))   

steady[ii] <-steady
if (!Dyna) return(list(steady=steady,precis=attr(DIA,"precis"),Solved=attr(DIA,"steady"))) else
{
# spinup for one year
times <- 0:365
out   <- ode.band (y=DIA$y, times=times, fun="omexdiamod",
                   initfun="initomexdia",method="lsode",
                   parms=c(pars,grid$dx,grid$dx.int,
                   porgrid$mid,porgrid$int,Db),
                   nspec=6,nout=8,dllname="ReacTran")

# final simulation
CONC  <- out[nrow(out),2:(1+6*N)]
out   <- ode.band (y=CONC, times=times, fun="omexdiamod",
                   initfun="initomexdia",method="lsode",
                   parms=c(pars,grid$dx,grid$dx.int,
                   porgrid$mid,porgrid$int,Db),
                   nspec=6,nout=8,dllname="ReacTran")
# remove time
out      <- out [,-1]

# rearrange from ordering per slice -> per spec
out[,ii] <-out[,1:(6*N)]
FDET <- out[,1:100] 
SDET <- out[,101:200] 
O2   <- out[,201:300] 
NO3  <- out[,301:400] 
NH3  <- out[,401:500] 
ODU  <- out[,501:600] 

vars <- out[,601:ncol(out)]

return(list(times=times,
FDET=FDET,SDET=SDET,O2=O2,NO3=NO3,NH3=NH3,ODU=ODU,
O2flux=vars[1], NO3flux=vars[3],    
NH3flux=vars[5],ODUflux=vars[7],
steady=steady))
} 
}
