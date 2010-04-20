\name{tran.polar}
\alias{tran.polar}
\alias{polar2cart}

\title{
  Diffusive Transport in polar (r, phi) coordinates.
}

\description{
  Estimates the transport term (i.e. the rate of change of a concentration
  due to diffusion) in a polar(r, phi)  coordinates
}

\usage{
tran.polar (C, C.r.up = NULL, C.r.down = NULL, C.phi.up = NULL, 
  C.phi.down = NULL, flux.r.up = NULL, flux.r.down = NULL, 
  flux.phi.up = NULL, flux.phi.down = NULL, 
  D.r = 1, D.phi = D.r, r = NULL, phi = NULL,
  full.check = FALSE , full.output = FALSE)

polar2cart (out, r, phi, x = NULL, y = NULL)

}

\arguments{
  \item{C }{concentration, expressed per unit volume, defined at the centre
    of each grid cell; Nr*Nphi matrix [M/L3].
  }
  \item{C.r.up }{concentration at upstream boundary in r(x)-direction;
    vector of length Nphi [M/L3].
  }
  \item{C.r.down }{concentration at downstream boundary in r(x)-direction;
    vector of length Nphi [M/L3].
  }
  \item{C.phi.up }{concentration at upstream boundary in phi-direction;
    vector of length Nr [M/L3].
  }
  \item{C.phi.down }{concentration at downstream boundary in phi-direction;
    vector of length Nr [M/L3].
  }
  \item{flux.r.up }{flux across the upstream boundary in r-direction,
    positive = INTO model domain; vector of length Nphi [M/L2/T].
  }
  \item{flux.r.down }{flux across the downstream boundary in r-direction,
    positive = OUT of model domain; vector of length Nhpi [M/L2/T].
  }
  \item{flux.phi.up }{flux across the upstream boundary in phi-direction,
    positive = INTO model domain; vector of length Nr [M/L2/T].
  }
  \item{flux.phi.down }{flux across the downstream boundary in phi-direction,
    positive = OUT of model domain; vector of length Nr [M/L2/T].
  }
  \item{D.r }{diffusion coefficient in r-direction, defined on grid cell
    interfaces. One value, a vector of length (Nr+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a (Nr+1)* Nphi matrix [L2/T].
  }
  \item{D.phi }{diffusion coefficient in y-direction, defined on grid cell
    interfaces. One value, a vector of length (Ny+1),
    a \code{prop.1D} list created by \code{\link{setup.prop.1D}},
    or a Nr*(Nphi+1) matrix [L2/T].
  }
  \item{r }{position of adjacent cell interfaces in the r-direction. 
    A vector of length Nr+1 [L].
  }
  \item{phi }{position of adjacent cell interfaces in the phi-direction.
    A vector of length Nphi+1 [L].
  }
  \item{full.check }{logical flag enabling a full check of the consistency
    of the arguments (default = \code{FALSE}; \code{TRUE} slows down
    execution by 50 percent).
  }
  \item{full.output }{logical flag enabling a full return of the output
    (default = \code{FALSE}; \code{TRUE} slows down execution by 20 percent).
  }
  \item{out }{output as returned by \code{tran.polar}, and which is to be
    mapped from polar to cartesian coordinates
  }
  \item{x }{The cartesian x-coordinates to whicht the polar coordinates are
    to be mapped
  }
  \item{y }{The cartesian y-coordinates to whicht the polar coordinates are
    to be mapped
  }
}
\value{
  a list containing:
  \item{dC }{the rate of change of the concentration C due to transport,
    defined in the centre of each grid cell, a Nr*Nphi matrix. [M/L3/T].
  }
  \item{C.r.up }{concentration at the upstream interface in r-direction.
    A vector of length Nphi [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{C.r.down }{concentration at the downstream interface in r-direction.
    A vector of length Nphi [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{C.phi.up }{concentration at the the upstream interface in phi-direction.
    A vector of length Nr [M/L3].  Only when \code{full.output = TRUE}.
  }
  \item{C.phi.down }{concentration at the downstream interface in phi-direction.
    A vector of length Nr [M/L3]. Only when \code{full.output = TRUE}.
  }
  \item{r.flux }{flux across the interfaces in x-direction of the grid cells.
    A (Nr+1)*Nphi matrix  [M/L2/T]. Only when \code{full.output = TRUE}.
  }
  \item{phi.flux }{flux across the interfaces in y-direction of the grid cells.
    A Nr*(Nphi+1) matrix [M/L2/T].  Only when \code{full.output = TRUE}.
  }
  \item{flux.r.up }{flux across the upstream boundary in r-direction,
    positive = INTO model domain. A vector of length Nphi [M/L2/T].
  }
  \item{flux.r.down }{flux across the downstream boundary in r-direction,
    positive = OUT of model domain. A vector of length Nphi [M/L2/T].
  }
  \item{flux.phi.up }{flux across the upstream boundary in phi-direction,
    positive = INTO model domain. A vector of length Nr [M/L2/T].
  }
  \item{flux.phi.down }{flux across the downstream boundary in phi-direction,
   positive = OUT of model domain. A vector of length Nr [M/L2/T].
  }

}

\examples{
## =============================================================================
## Testing the functions
## =============================================================================
# Parameters
F        <- 100             # input flux [micromol cm-2 yr-1]
D        <- 400             # mixing coefficient [cm2 yr-1]

# Grid definition
r.N   <- 4   # number of cells in r-direction
phi.N <- 6   # number of cells in phi-direction
r.L <- 8   # domain size x-direction [cm]
r    <- seq(0,r.L,len=r.N+1)             # cell size x-direction [cm]
phi  <- seq(0,2*pi,len=phi.N+1)          # phi-direction - phi: from 0, 2pi
 
# Intial conditions 
C <- matrix(nrow=r.N, ncol=phi.N, data=0, byrow=FALSE)

# Boundary conditions: fixed concentration  
C.r.up <- rep(1, times=phi.N)
C.r.down <- rep(0, times=phi.N)
C.phi.up   <- rep(1, times=r.N)
C.phi.down <- rep(0, times=r.N)

# Concentration boundary conditions
tran.polar(C=C, D.r=D, D.phi=D, 
  C.r.up=C.r.up, C.r.down=C.r.down, r=r, phi=phi,
  C.phi.up=C.phi.up,C.phi.down=C.phi.down)

# Flux boundary conditions
flux.r.up <- rep(200, times=phi.N)
flux.r.down <- rep(-200, times=phi.N)
flux.phi.up <- rep(200, times=r.N)
flux.phi.down <- rep(-200, times=r.N)

tran.polar(C=C, D.r=D, r=r, phi=phi,
  flux.r.up=flux.r.up, flux.r.down=flux.r.down,
  flux.phi.up=flux.phi.up, flux.phi.down=flux.phi.down,
  full.output = TRUE)


## =============================================================================
## A model with diffusion and first-order consumption
## =============================================================================
N     <- 50          # number of grid cells
XX    <- 4          # total size
rr    <- 0.005       # consumption rate
ini   <- 1           # initial value at x=0
D     <- 400

r     <- seq (2, 4, len=N+1)
phi   <- seq(0, 2*pi, len = N+1)
phi.m <- 0.5*(phi[-1]+phi[-(N+1)])

# The model equations

Diffpolar <- function (t,y,parms)  {
  CONC  <- matrix(nrow = N, ncol = N, data = y)
  dCONC <- tran.polar(CONC, D.r = D, D.phi = D, r = r, phi = phi,
        C.r.up = 0,          C.r.down = 1*sin(5*phi.m), 
        C.phi.up = CONC[,N], C.phi.down = CONC[,1]     )$dC 
        - rr * CONC
  return (list(dCONC))
}

# initial condition: 0 everywhere, except in central point
y <- matrix(nr = N, nc = N, data = 0)
N2  <- ceiling(N/2)

y[N2,N2] <- ini  # initial concentration in the central point...

# solve to steady-state; cyclicBnd = 2, because of C.phi.up, C.phi.down
out <- steady.2D (y = y, func = Diffpolar, parms = NULL,
                dim = c(N,N), lrw = 1e6, cyclicBnd = 2)

image(out)

cart <- polar2cart(out, r = r, phi = phi, 
                        x = seq(-4, 4, len = 100), 
                        y = seq(-4, 4, len = 100))
image(cart)

}

\references{
  Soetaert and Herman, 2009. a practical guide to ecological modelling -
  using R as a simulation platform. Springer
}

\details{
  \code{tran.polar} performs (simplified) transport in polar coordinates

  The \bold{boundary conditions} are either
  \itemize{
    \item (1) fixed concentration
    \item (2) fixed flux
  }

  This is also the order of priority. The fixed flux overrules the other.
  
  \code{polar2cart} maps the polar coordinates to cartesian coordinates
  
  If \code{x} and \code{y} is not provided, then it will create an (x,y)
  grid based on \code{r} : \code{seq(-maxr, maxr, length.out=Nr)}, where
  \code{maxr} is the maximum value of \code{r}, and \code{Nr} is the number
  of elements in \code{r}.
}  

\keyword{utilities}

