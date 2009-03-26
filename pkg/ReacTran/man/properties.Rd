\name{p.exp.increase}
\alias{p.exp.increase}
\alias{p.exp.decrease}
\alias{p.sphere.surf}
\alias{p.sphere.vol}
\alias{p.spheroid.surf}
\alias{p.spheroid.vol}
\alias{p.cylinder.surf}
\alias{p.cylinder.vol}
\title{
  Properties to be used with setup.prop.1D
}

\description{
  Functions that define properties as a function of x. To be used with
  \link{setup.prop.1D}
  
  \itemize{
    \item \code{p.exp.decrease}: exponentially decreasing function
      \deqn{
         y = y_{\inf} + (y_0-y_{\inf}) \exp(-\max(0,x-x_0)/x_a)
      }{y=y0+(y0-yinf)*exp(-max(0,(x-x0))/xa}

    \item \code{p.exp.increase}: exponentially increasing function
      \deqn{
         y = 1- (y_{\inf} + (y_0-y_{\inf}) \exp(-\max(0,x-x_0)/x_a))
      }{y=1-(y0+(y0-yinf)*exp(-max(0,(x-x0))/xa)}
    \item \code{p.sphere.surf, p.sphere.vol} the surface and volume of a
      sphere
      
    \item \code{p.spheroid.surf, p.spheroid.vol} the surface and volume of a
      spheroid

    \item \code{p.cylinder.surf, p.cylinder.vol} the surface and volume of a
      cylinder; note that the surface ignores the top and bottom.
  }
}

\usage{
p.exp.decrease (x, x.0=0, y.0=1, y.inf=0.5, x.att=1)

p.exp.increase (x, x.0=0, y.0=1, y.inf=0.5, x.att=1)

p.sphere.surf (x)

p.sphere.vol (x)

p.spheroid.surf (x, b=1)

p.spheroid.vol (x, b=1)

p.cylinder.surf (x, L=1)

p.cylinder.vol (x, L=1)
}

\arguments{
  \item{x }{the x-values for which the property has to be calculated.
  }
  \item{x.0 }{the x-offset; for \code{x <= x.0}, the value will be equal
    to \code{y.0}.
  }
  \item{y.0 }{the y-value for \code{x <= x.0}.
  }
  \item{y.inf }{the y-value for \code{x} = infinite
  }
  \item{x.att }{the attenuation coefficient (exponential decrease value)
  }
  \item{b }{the ratio of long/short radius of the spheroid;
    if b<1: the spheroid is oblate.
  }
  \item{L }{the length of the cylinder
  }
}

\value{
  the property value, estimated for each x-value.
}

\author{
  Filip Meysman <f.meysman@nioo.knaw.nl>,
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
 x<- seq(0,20,len=100)
 plot(x, p.exp.increase(x),main="p.exp.increase")
 lines(x, p.exp.increase(x,x.0 = 5))

 x<- seq(0,20,len=100)
 plot(x, p.exp.decrease(x),main="p.exp.decrease")
 lines(x, p.exp.decrease(x,x.0 = 5))

 mf <- par(mfrow=c(3,2))
 plot(x, p.sphere.surf(x),main="p.sphere.surf")
 plot(x, p.sphere.vol(x),main="p.sphere.vol")
 plot(x, p.spheroid.surf(x,b=0.5),main="p.spheroid.surf")
 plot(x, p.spheroid.vol(x,b=0.5),main="p.spheroid.vol")
 plot(x, p.cylinder.surf(x,L=1),main="p.cylinder.surf")
 plot(x, p.cylinder.vol(x,L=1),main="p.cylinder.vol")
 par("mfrow"=mf)
}
\keyword{utilities}


