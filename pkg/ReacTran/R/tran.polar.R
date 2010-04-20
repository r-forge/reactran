
##==============================================================================
## Diffusive transport in polar coordinates (x, phi)
##==============================================================================

tran.polar <- function(C, C.r.up=NULL,    C.r.down=NULL, 
                        C.phi.up=NULL,  C.phi.down=NULL, 
                        
                        flux.r.up=NULL,   flux.r.down=NULL, 
                        flux.phi.up=NULL, flux.phi.down=NULL, 
                        
                        D.r=1, D.phi=D.r,
                        
                        r=NULL, phi=NULL,
                        
                        full.check = FALSE, full.output = FALSE )
											
{
# ------------------------------------------------------------------------------
# Initialisation
# ------------------------------------------------------------------------------

# a function to check the dimensionality of the system
  checkDim <- function (dr, dphi, nx="dr", nphi="dphi") {
    if (length(dr) != 1 )
      stop (nx,", should have length 1")
    if (length(dphi) != 1)
       stop (nphi,", should have length 1")
  }

  checkDim2 <- function (dr,dphi,nx="dr",nphi="dphi") {
    if (length(dr) != dimC[1]+1)
      stop (nx,", should have length equal to dim(C)[1]+1 = ",dimC[1]+1," it is ",length(dr))
 
    if (length(dphi) != dimC[2]+1)
       stop (nphi,", should have length equal to dim(C)[2]+1 = ",dimC[2]+1," it is ",length(dphi))
  }

  dimC <- dim (C) 
  Type <- 2

  if (max(phi)> 2*pi)
     stop ("phi should be < 2pi")
  if (min(phi)< 0)
     stop ("phi should be > 0 ")
  
  checkDim2(r, phi,"r, the grid in x-direction",
                   "phi, the grid in phi-direction")
  # central values
  Nr       <- dimC[1]
  r.c      <- 0.5*(r[-1]+r[-(Nr+1)])
  dr       <- diff(r)
  drint    <- c(r.c[1]-r[1],diff(r.c),r[Nr+1]-r.c[Nr])
  divr     <- 1/r
  divr[is.na(divr)] <- 0
  
  Nphi <- dimC[2]
  phi.c    <- 0.5*(phi[-1]+phi[-(Nphi +1)])
  dphi     <- diff(phi)
  dphiint  <- c(phi.c[1]-phi[1],diff(phi.c),phi[Nphi+1]-phi.c[Nphi])

  checkDim (D.r,D.phi,"D.r, the diffusion in x-direction ",
                      "D.phi, the diffusion in phi-direction ")
  
    
# check the dimensionality of the boundaries of system

    if (! is.null(C.r.up)){ 
     if( length(C.r.up) != 1 && length(C.r.up) != Nphi)
       stop ("'C.r.up' should have length 1 or equal to ", Nphi)
    } else if (! is.null(flux.r.up)) { 
     if( length(flux.r.up) != 1 && length(flux.r.up) != Nphi)
      stop ("'flux.r.up' should have length 1 or equal to ", Nphi)
    } else  
      stop ("'flux.r.up' OR 'C.r.up' should be specified")

    if (! is.null(C.r.down)) {
     if( length(C.r.down) != 1 && length(C.r.down) != Nphi)
      stop ("'C.r.down' should have length 1 or equal to ", Nphi)
    } else if (! is.null(flux.r.down)){
     if(length(flux.r.down) != 1 && length(flux.r.down) != Nphi)
      stop ("'flux.r.down' should have length 1 or equal to ", Nphi)
    } else  
      stop ("'flux.r.down' OR 'C.r.down' should be specified")

    if (! is.null(C.phi.up)) {
     if ( length(C.phi.up) != 1 && length(C.phi.up) != Nr)
      stop ("'C.phi.up' should have length 1 or equal to ", Nr)
    } else if (! is.null(flux.phi.up)){
     if (length(flux.phi.up) != 1 && length(flux.phi.up) != Nr)
       stop ("'flux.phi.up' should have length 1 or equal to ", Nr)
    } else  
      stop ("'flux.phi.up' OR 'C.phi.up' should be specified")

    if (! is.null(C.phi.down)) { 
     if(length(C.phi.down) != 1 && length(C.phi.down) != Nr)
      stop ("'C.phi.down' should have length 1 or equal to ", Nr)
    } else if (! is.null(flux.phi.down)) {
      if( length(flux.phi.down) != 1 && length(flux.phi.down) != Nr)
       stop ("'flux.phi.down' should have length 1 or equal to ", Nr)
    } else  
      stop ("'flux.phi.down' OR 'C.phi.down' should be specified")


# ------------------------------------------------------------------------------
# Function body: calculation
# ------------------------------------------------------------------------------

## Initialise 'fluxes' in all directions
  Flux.r <- 0
  Flux.phi <- 0

## First rewrite boundary values as "concentration"
## then perform diffusive transport


   if (is.null(C.r.up))
     C.r.up <-  flux.r.up / D.r * drint[1]+ C[1,]
   if (is.null(C.r.down))
     C.r.down <- - flux.r.down / D.r * drint[Nr+1] + C[Nr,]
  
   Flux.r <- D.r * (rbind(C.r.up,C, deparse.level = 0) -
                     rbind(C,C.r.down, deparse.level = 0))/drint  

   if (is.null(C.phi.up))
     C.phi.up <- flux.phi.up/D.phi * dphiint[1]*r.c + C[,1]
   if (is.null(C.phi.down))
     C.phi.down <- - flux.phi.down /D.phi * dphiint[Nphi+1]*r.c + C[,Nphi]

   
   Flux.phi <- D.phi * (cbind(C.phi.up,C,deparse.level = 0) - 
                        cbind(C,C.phi.down, deparse.level = 0))/
                matrix(data= dphiint,nrow=Nr,ncol=(Nphi+1),byrow=TRUE) /r.c

## Calculate rate of change = flux gradient 

  dFdphi <- 0.
  dFdr   <- -diff(Flux.r * r)/dr/r.c
  dFdphi <- -t(diff(t(Flux.phi))/dphi)/r.c

 if (full.output)
  return (list (dC = dFdr + dFdphi ,                  # Rate of change due to advective-diffuisve transport in each grid cell
                C.r.up     = C.r.up,
                C.r.down   = C.r.down,
                C.phi.up   = C.phi.up,
                C.phi.down = C.phi.down,
                r.flux   = Flux.r,                 
                phi.flux = Flux.phi,
                flux.r.up = Flux.r[1,],
                flux.r.down = Flux.r[Nr+1,],
                flux.phi.up = Flux.phi[,1],
                flux.phi.down = Flux.phi[,Nphi+1]
                
                ))  
 else
  return (list (dC = dFdr + dFdphi ,                  # Rate of change due to advective-diffuisve transport in each grid cell
                flux.r.up = Flux.r[1,],
                flux.r.down = Flux.r[Nr+1,],
                flux.phi.up = Flux.phi[,1],
                flux.phi.down = Flux.phi[,Nphi+1]
                ))  
 
} # end tran.polar
