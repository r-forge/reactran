## =============================================================================
## Maps a matrix u from (x,y) to (x2,y2)
## =============================================================================
mapxy <- function(x, y, x2, y2, u) {
  Nx <- length(x)
  Ny <- length(y)

  dx   <- c(diff(x),1)  # 1= for last value
  dy <- c(diff(y),1)

  Du <- dim(u)
  if (Du[1] != Nx)
    stop("u and x not compatible")
  if (Du[2] != Ny)
    stop("u and y not compatible")

  Transf <- function (x2,y2,u) {
    u    <- rbind(c(u[1,1],u[1,]) , cbind(u[,1],u))

  # find embracing values : first interval
    ix <- findInterval(x2, x )
    iy <- findInterval(y2, y)

  # next interval
    ixp1 <- pmin(ix+1,Nx)
    iyp1 <- pmin(iy+1,Ny)

  # interpolation factor
    xfac <- (x2-x[ix])/dx[ix]
    yfac <- (y2-y[iy])/dy[iy]

  # interpolate
    (1-yfac)*((1-xfac)*u[cbind(ix,iy)]+xfac*u[cbind(ixp1,iy)]) +
    yfac*((1-xfac)*u[cbind(ix,iyp1)]+xfac*u[cbind(ixp1,iyp1)])

  } # end Transf

outer(x2, y2, FUN = Transf, u = u)
}

#mapxy(1:87, 1:61, seq(1,87,0.25), seq(1,61,0.25),volcano)->VOLCANO


## =============================================================================
## From polar to cartesian coordinates
## =============================================================================

polar2cart <- function(out, r, phi, x = NULL, y = NULL)
{
  classout <- class(out)[1]
  if (!classout %in% c("steady2D", "deSolve"))
    stop ("Class of 'out' should be one of steady2D or deSolve")

  ATTR <- attributes(out)
  if (length(ATTR$dimens)!= 2)
    stop ("input should be 2-dimensional")

  Nr   <- length(r)  -1
  Np   <- length(phi)-1
  maxr   <- max(r)
  maxphi <- max(phi)
  minr   <- min(r)
  minphi <- min(phi)

  # add leftmost boundary
  r    <- c(r[1],   0.5*(r[-1]+r[-(Nr+1)]), r[Nr+1])
  phi  <- c(phi[1], 0.5*(phi[-1]+phi[-(Np+1)]), phi[Np+1])
                                                # check dimensions...
  dr   <- c(diff(r),1)       # 1= for last value: no interpolation
  dphi <- c(diff(phi),1)

  if (is.null (x))
    if (! is.null(y)) x<-y
  if (is.null (y)) {
    mr <- max(abs(r))
    x <- y <- seq(-mr, mr, len=Nr)
  }


  Transf <- function(x, y,  u) {

    Du <- dim(u)
    if (Du[1] != Nr)
      stop("u and r not compatible")
    if (Du[2] != Np)
      stop("u and phi not compatible")
  # augment u with boundary values
    u    <- rbind(c(u[1,1],u[1,], u[1,Np]) , 
                  cbind(u[,1],u, u[,Np]),
                  c(u[Nr,1],u[Nr,], u[Nr,Np]) )

    R   <- sqrt(x^2+y^2)
    Phi <- atan(y/x)
    Phi[x<0] <- Phi[x<0] + pi
    Phi[x>0 & y < 0]   <-   Phi[x > 0 & y < 0] + 2*pi
    Phi[x==0 & y > 0]  <- pi/2
    Phi[x==0 & y < 0]  <- 3*pi/2
    Phi[x==0 & y == 0] <- 0

  # find embracing values : interval location
    iR <- findInterval(R, r )
    iP <- findInterval(Phi,phi)
    iR [iR == 0]   <- NA
    iR [iR > Nr+1] <- NA
    iP [iP == 0]   <- NA
    iP [iP > Np+1] <- NA

  # next interval
    iPp1 <- pmin(iP+1,Np+1)
    iRp1 <- pmin(iR+1,Nr+1)

  # interpolation factor
    iRfac <- (R-r[iR])/dr[iR]
    iPfac <- (Phi-phi[iP])/dphi[iP]

  # interpolate inbetween 4 values
    (1-iPfac)*((1-iRfac)*u[cbind(iR,iP)]+iRfac*u[cbind(iRp1,iP)]) +
    iPfac*((1-iRfac)*u[cbind(iR,iPp1)]+iRfac*u[cbind(iRp1,iPp1)])

  } # end Transf


  Num   <- prod(ATTR$dimens)
  nspec <- ATTR$nspec
  if (classout[1] == "steady2D") {
    MAP   <- list()
    MAP$y <- NULL
    for (i in 1:nspec)  {
      ii <- (1+(i-1)*Num):(i*Num)
      U <- matrix(out$y[ii], nrow = ATTR$dimens[1], ncol= ATTR$dimens[2])
      MAP$y <- c(MAP$y, as.vector(outer(x, y, FUN = Transf, u = U)))
    }
  attributes(MAP) <- ATTR
  } else {               # Dynamic output
    MAP <- NULL
    nr = nrow(out)
    for ( j in 1: nr) {
      mm <- out[j,1]
      for (i in 1:nspec)  {
        ii <- ((1+(i-1)*Num):(i*Num))+1
        U <- matrix(out[j,ii], nrow = ATTR$dimens[1], ncol= ATTR$dimens[2])
        mm <- c(mm, as.vector(outer(x, y, FUN = Transf, u = U)))
      }
      MAP <- rbind(MAP, mm)
    }
  attributes(MAP)[c("istate","rstate","class","type","nspec","dimens")] <- 
    ATTR[c("istate","rstate","class","type","nspec","dimens")]
  }
  attributes(MAP)$dimens     <- c(length(x), length(y))
  attributes(MAP)$xgrid      <- x
  attributes(MAP)$ygrid      <- y
 
  MAP
}
