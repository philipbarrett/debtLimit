#####################################################################
# sol.R
#
# Contains the functions for solving the model
# 13mar2017
# Philip Barrett, Washington DC
#####################################################################

sol.wrapper <- function( params, init.guess, st.which.sol, qd=c(0), qe=c(0), def=matrix(0) ){
# Wraps and processes the solution output from the nonlinear solver

  if(st.which.sol == 'core.nl'){
    sol <- sol.core.nl(params, init.guess, qd, qe, def)
  }else{
    stop('Solution method not recognized')
  }

  n <- length(sol$x) / 2
  out <- list()
  out$p <- sol$x[1:n]
  out$d <- sol$x[n+1:n]
      # Reorder the solution elements
  out$err <- matrix( sol$fvec, nrow = n, ncol = 2 )
      # The errors
  out$accy <- max( abs( out$err ) ) < 1e-08
      # Boolean for accuracy of the solution
  out$qd <- qd
  out$qe <- qe
  out$def <- def
      # Record inputs
  out$sol <- sol
      # Extra solution info
  ## NEED TO CHECK FOR LOCAL/GLOBAL CONCAVITY
  return(out)
}

sol.core.nl <- function(params, init.guess, qd, qe, def){
# The core solution routine using nleqslv w/o derivatives.  Inital guess is a matrix.

  x0 <- c( init.guess )
      # Reformat as a vector
  n <- nrow(init.guess)
      # The number of states
  fn <- function(x){
      # The function for which we want to find the root
    p <- x[1:n]
    d <- x[n+1:n]
        # Extract the values
    z.2 <- zed_2( p, d, params, qd, qe, def )
        # The vector of: the implied value of p and the gradient wrt p
    return( z.2 - cbind(p,rep(1,n)) )
        # Becuase we want the slope to be unity at the fixed point
  }

  control <- list( maxit=500, allowSingular =TRUE, ftol=1e-10, xtol=1e-12 )
  sol <- nleqslv( x0, fn, control = control )
      # The solution object
  return( sol )
}
