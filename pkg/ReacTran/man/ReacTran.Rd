\name{ReacTran-package}
\alias{ReacTran-package}
\alias{ReacTran}
\docType{package}
\title{
  Reactive transport modelling in 1D or 2-D
}
\description{
  
  The R-package ReacTran contains routines that enable the development of
  reactive transport models in aquatic
  systems (rivers, lakes), porous media (floc aggregates, sediments,...) and even idealized organisms (spherical cells, cylindrical worms,...). The geometry of the model domain is either one-dimensional or two-dimensional. 
  
  It contains:
  \itemize{
    \item Functions to setup a finite-difference grid (1D or 2D)
    \item Functions to attach parameters and properties to this grid (1D or 2D)
    \item Functions to calculate the advective-diffusive transport term over
      the grid (1D or 2D)
    \item Utility functions
  }
  
  
}
\details{
  \tabular{ll}{
    Package: \tab ReacTran\cr
    Type: \tab Package\cr
    Version: \tab 1.1\cr
    Date: \tab 2009-03-27\cr
    License: \tab GNU Public License 2 or above\cr
  }
}
\author{
  Karline Soetaert (Maintainer)
  
  Filip Meysman
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
