## =============================================================================
##
## Several useful properties to be used with setup.prop
##
## =============================================================================

## exponential decline

  p.exp.decrease <- function(x, x.0=0, y.0=1, y.inf=0.5, x.att=1)
       return(y.inf + (y.0-y.inf)*exp(-pmax(x-x.0,0)/x.att))

## exponential increase

 p.exp.increase <- function(x, x.0=0, y.0=1, y.inf=0.5, x.att=1)
       return(1 - (y.inf + (y.0-y.inf)*exp(-pmax(x-x.0,0)/x.att)))


# surface and volume of a sphere
  p.sphere.surf <- function (x) 4*pi*x*x
  p.sphere.vol <- function (x) 4/3*pi*x^3
 

# surface and volume of a spheroid
  p.spheroid.surf <- function (x, b=1) {     # b = long/short radius
    bb <- x*b
    te <- acos(b)
    if (b < 1)      #oblate
      2*pi *(x^2 + bb^2 / sin(te)*log((1+sin(te))/cos(te)))
    else if (b > 1)
      2*pi*(x^2 + x*bb*te/sin(te))
    else
      4 * pi*x^2
  }

  p.spheroid.vol <- function (x, b=1) {
    bb <- x*b
    4/3*pi * x^2 * bb
  }

  p.cylinder.surf <- function (x,L=1){
    2*pi*x*L
  }

  p.cylinder.vol <- function(x,L=1){
    pi*x^2*L
  }
