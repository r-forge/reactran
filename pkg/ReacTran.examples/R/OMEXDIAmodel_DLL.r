

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
                          Dyna = TRUE)  {
  N <- grid$N
  if (length(Db) != 101)
    stop("Db should contain 101 values")
  if (grid$N != 100)
    stop("grid should contain 100 boxes")
  if (length(porgrid$mid) != 100)
    stop("porosity grid should contain 100 boxes")
  if (is.null(grid$dx))
    stop("'grid$dx' should be present")
  if (is.null(grid$dx.aux))
    stop("'grid$dx.aux' should be present")
  if (is.null(porgrid$mid))
    stop("'porgrid$mid' should be present")
  if (is.null(porgrid$int))
    stop("'porgrid$int' should be present")

## The parameters
## the default parameter values

  Parms <- c(
  ## organic matter dynamics  #
    MeanFlux = 20000/12*100/365,  # nmol/cm2/d - Carbon deposition: 20gC/m2/yr
    rFast    = 0.01            ,  #/day        - decay rate fast decay detritus
    rSlow    = 0.00001         ,  #/day        - decay rate slow decay detritus
    pFast    = 0.9             ,  #-           - fraction fast detritus in flux
    w        = 0.1/1000/365    ,  # cm/d       - advection rate
    NCrFdet  = 0.16            ,  # molN/molC  - NC ratio fast decay detritus
    NCrSdet  = 0.13            ,  # molN/molC  - NC ratio slow decay detritus

  ## oxygen and DIN dynamics  #

  ## Nutrient bottom water conditions
    bwO2            = 300      ,    #mmol/m3     Oxygen conc in bottom water
    bwNO3           = 10       ,    #mmol/m3
    bwNH3           = 1        ,    #mmol/m3
    bwODU           = 0        ,    #mmol/m3

  ## Nutrient parameters
    NH3Ads          = 1.3      ,    #-           Adsorption coeff ammonium
    rnit            = 20.      ,    #/d          Max nitrification rate
    ksO2nitri       = 1.       ,    #umolO2/m3   half-sat O2 in nitrification
    rODUox          = 20.      ,    #/d          Max rate oxidation of ODU
    ksO2oduox       = 1.       ,    #mmolO2/m3   half-sat O2 in oxidation of ODU
    ksO2oxic        = 3.       ,    #mmolO2/m3   half-sat O2 in oxic mineralisation
    ksNO3denit      = 30.      ,    #mmolNO3/m3  half-sat NO3 in denitrification
    kinO2denit      = 1.       ,    #mmolO2/m3   half-sat O2 inhib denitrification
    kinNO3anox      = 1.       ,    #mmolNO3/m3  half-sat NO3 inhib anoxic degr
    kinO2anox       = 1.       ,    #mmolO2/m3   half-sat O2 inhib anoxic min

  ## Diffusion coefficients, temp = 10dgC
    #Temp            = 10                     ,   # temperature
    DispO2          = 0.955    +10*0.0386    ,  #cm2/d
    DispNO3         = 0.844992 +10*0.0336    ,
    DispNH3         = 0.84672  +10*0.0336    ,
    DispODU         = 0.8424   +10*0.0242)


## check parameter inputs

  if (length(pars)) {
    nms <- names(Parms)
    Parms[(namc <- names(pars))]<-pars
    if (length(noNms <- namc[!namc %in% nms]) > 0)
      warning("unknown names in pars: ", paste(noNms, collapse = ", "))
  }

## solve the model
## First the steady-state condition
  OC   <- rep(10,6*N)
  DIA  <- steady.band(y=OC,fun="omexdiamod",initfun="initomexdia",
                     initpar=c(Parms,grid$dx,grid$dx.aux,
                     porgrid$mid,porgrid$int,Db),nspec=6,
                     dllname="ReacTran.examples",nout=8,positive=TRUE)
  steady <-DIA$y

## state variables rearranged from ordering per slice -> per spec
  ii       <- as.vector(t(matrix(ncol=6,1:(6*N))))

  steady[ii] <-steady
  Parms <- unlist(Parms)
  if (!Dyna)
    return(list(steady=steady, precis=attr(DIA,"precis"),
            Solved=attr(DIA,"steady"))) else {
    ## spinup for one year
    times <- 0:365
    out   <- ode.band (y=DIA$y, times=times, fun="omexdiamod",
                       initfun="initomexdia",method="lsode",
                       parms=c(Parms,grid$dx,grid$dx.aux,
                       porgrid$mid,porgrid$int,Db),
                       nspec=6,nout=8,dllname="ReacTran.examples")

    ## final simulation
    CONC  <- out[nrow(out),2:(1+6*N)]
    out   <- ode.band (y=CONC, times=times, fun="omexdiamod",
                       initfun="initomexdia",method="lsode",
                       parms=c(Parms,grid$dx,grid$dx.aux,
                       porgrid$mid,porgrid$int,Db),
                       nspec=6,nout=8,dllname="ReacTran.examples")
    ## remove time
    out      <- out [,-1]

    ## rearrange from ordering per slice -> per spec
    out[,ii] <-out[,1:(6*N)]
    FDET <- out[,1:100]
    SDET <- out[,101:200]
    O2   <- out[,201:300]
    NO3  <- out[,301:400]
    NH3  <- out[,401:500]
    ODU  <- out[,501:600]

    vars <- out[,601:ncol(out)]

    return ( list (times=times, FDET=FDET, SDET=SDET, O2=O2,
                   NO3=NO3, NH3=NH3, ODU=ODU, O2flux=vars[1],
                   NO3flux=vars[3], NH3flux=vars[5],
                   ODUflux=vars[7], steady=steady))
  }
}
