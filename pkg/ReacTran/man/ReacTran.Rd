\name{ReacTran-package}
\alias{ReacTran-package}
\alias{ReacTran}
\docType{package}
\title{
Reactive transport modelling in 1D or 2-D
}
\description{
  
  R-package ReacTran contains routines for developing models that describe reaction and
  advective-diffusive transport in one or two dimensions.

  It can be applied to describe transport in porous media, in estuaries, and in bodies with
  variable shape.

  It contains:
  \itemize{
    \item Functions to transport properties in one or two dimensions
    \item Functions to setup a grid, usable in transport
    \item Utility functions
  }
  
  
}
\details{
  \tabular{ll}{
    Package: \tab ReacTran\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2009-02-02\cr
    License: \tab GNU Public License 2 or above\cr
  }
}
\author{
  Karline Soetaert (Maintainer)
  
  Filip Meysman
}
\references{
  
% to be filled in...

}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(tran.1D)
example(tran.2D)

## open the directory with documents
browseURL(paste(system.file(package="ReacTran"), "/doc", sep=""))

}
}

\keyword{ package }
