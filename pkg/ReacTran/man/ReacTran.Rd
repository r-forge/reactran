\name{ReacTran-package}
\alias{ReacTran-package}
\alias{ReacTran}
\docType{package}
\title{
  Reactive transport modelling in 1D or 2-D
}
\description{
  
  R-package ReacTran contains routines that enable the development of
  one-dimensional and two-dimensional reactive transport models in aquatic
  systems (rivers, lakes) as well porous media (floc aggregates, sediments,...).
  
  It contains:
  \itemize{
    \item Functions to setup a finite-difference grid (1D or 2D)
    \item Functions to attach parameters and properties to the grid (1D or 2D)
    \item Functions to calculate the advective-diffusive transport term over
      the grid (1D or 2D)
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
