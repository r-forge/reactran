
##################################################################
######          OMEXDIA: C, N, O2 diagenesis                ######
##################################################################

OMEXDIAparms<- c(
# organic matter dynamics  #
MeanFlux = 20000/12*100/365,  # nmol/cm2/d - Carbon deposition: 20gC/m2/yr
rFast    = 0.01            ,  #/day        - decay rate fast decay detritus
rSlow    = 0.00001         ,  #/day        - decay rate slow decay detritus
pFast    = 0.9             ,  #-           - fraction fast detritus in flux
w        = 0.1/1000/365    ,  # cm/d       - advection rate
NCrFdet  = 0.16            ,  # molN/molC  - NC ratio fast decay detritus
NCrSdet  = 0.13            ,  # molN/molC  - NC ratio slow decay detritus

# oxygen and DIN dynamics  #

# Nutrient bottom water conditions
bwO2            = 300      ,    #mmol/m3     Oxygen conc in bottom water
bwNO3           = 10       ,    #mmol/m3
bwNH3           = 1        ,    #mmol/m3
bwODU           = 0        ,    #mmol/m3

# Nutrient parameters
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

# Diffusion coefficients, temp = 10dgC
#Temp            = 10                     ,   # temperature 
DispO2          = 0.955    +10*0.0386    ,  #cm2/d
DispNO3         = 0.844992 +10*0.0336    ,
DispNH3         = 0.84672  +10*0.0336    ,
DispODU         = 0.8424   +10*0.0242)

