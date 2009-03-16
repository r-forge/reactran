
##==============================================================================
## Attaches a property to a 2D grid
##==============================================================================

setup.prop.2D <- function(func = NULL, value = NULL, grid, ...) {

  ## check input
  gn <- names(grid)
  if (! "x.mid"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.mid")
  if (! "x.int"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.int")
  if (! "y.mid"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.mid")
  if (! "y.int"  %in% gn)
    stop("error in setup.2D.prop: grid should be a list that contains x.int")

  if (is.null(func) && is.null(value))
    stop("error in setup.prop: function and value should not be both NULL")

  if (!is.null(value)) { # profile specification via value
    x.int <- matrix(nrow=length(grid$x.int),ncol=length(grid$y.mid),data=value)
    y.int <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.int),data=value)
    mid <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.mid),data=value)
  }

  if (!is.null(func)) { # profile specification via function
    x.int <- matrix(nrow=length(grid$x.int),ncol=length(grid$y.mid))
    y.int <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.int))
    mid <- matrix(nrow=length(grid$x.mid),ncol=length(grid$y.mid))

    for (i in 1:length(grid$x.int))
      x.int[i,] <- func(grid$x.int[i],grid$y.mid,...)
    for (i in 1:length(grid$x.mid))
      y.int[i,] <- func(grid$x.mid[i],grid$y.int,...)
    for (i in 1:length(grid$x.mid))
      mid[i,] <- func(grid$x.mid[i],grid$y.mid,...)
  }

  Res <- list(mid  = mid,
              x.int  = x.int,
              y.int =  y.int)
  class(Res) <- "prop.2D"
  return(Res)
}

