\name{p.exp}
\alias{p.exp}
\alias{p.lin}
\alias{p.sig}

\title{
  Common depth dependencies for one-dimensional properties to be used with setup.prop.1D
}

\description{
  Functions that define an y-property as a function of the one-dimensional x-coordinate. These routines can be used to specify properties and parameters as a function of depth in the water column or the sediment. They make a transition from an upper (or upstream) with value \code{y.0} to a lower zone with a value \code{y.inf}. The width of the transition zone is \code{x.att} and the depth where the transiztion zone starts is \code{x.L}. Particularly useful in combination with
  \link{setup.prop.1D}
  
  \itemize{
    \item \code{p.exp}: exponentially decreasing transition 
      \deqn{
         y = y_{\inf} + (y_0-y_{\inf}) \exp(-\max(0,x-x_0)/x_a)
      }{y=y0+(y0-yinf)*exp(-max(0,(x-x0))/xa}

    \item \code{p.lin}: linearly decreasing transition 
      \deqn{
         0 <= x <= x_L => y = y_0
         x_L < x <= x_L + x_{att} => y = y_0 - (y_0-y_{inf})*(x-x_L)/x_{att})
         x >= x_L + x.att => y = y_{inf} 
        }
    \item \code{p.sig}: sogmoidal decreasing transition 
      \deqn{
         y = y_{inf} + (y_0-y_{inf})\exp(-(x-(x_L + 0.5 x_{att}))/(0.25 x_{att}))/(1+\exp(-(x-(x_L+0.5*x_{att}))/(0.25 x_{att})))
        }
  }
}

\usage{
p.exp(x, y.0=1, y.inf=0.5, x.L=0, x.att=1)

p.lin(x, y.0=1, y.inf=0.5, x.L=0, x.att=1)

p.sig(x, y.0=1, y.inf=0.5, x.L=0, x.att=1)
}

\arguments{
  \item{x }{the x-values for which the property has to be calculated.
  }
  \item{y.0 }{the y-value at the origin 
  }
  \item{y.inf }{the y-value at infinity 
  }
  \item{x.L }{the x-coordinate where the transition zone starts; for \code{x <= x.0}, the value will be equal
    to \code{y.0}. For \code{x >> x.L + x.att} the value will be tend to \code{y.inf}
  }
  \item{x.att }{attenuation coefficient in exponential decrease, or the size of the transition zone in the linear and sigmoid decrease
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

 x<- seq(0,5,len=100)
 plot(x, p.exp(x,x.L=2),xlab="x.coordinate", ylab="y value",ylim=c(0,1))
 lines(x, p.lin(x,x.L=2),col="blue")
 lines(x, p.sig(x,x.L=2),col="red")
}
 
\keyword{utilities}


