#####################################################################
# sol.R
#
# Contains the functions for solving the model
# 13mar2017
# Philip Barrett, Washington DC
#####################################################################

sol.wrapper <- function( params, init.guess, st.which.sol,
                         An=c(0), Bn=c(0), Cn=c(0), def=matrix(0), i=1 ){
# Wraps and processes the solution output from the nonlinear solver

  n <- length(params$R)
  out <- list()
  test.p <- 0
  min.p <- 1
  maxit <- 5
  it <- 0
  guess <- init.guess

  while( any( abs( test.p - min.p ) > 1e-03 ) && it < maxit ){
    it <- it+1
    if(st.which.sol == 'core.nl'){
      sol <- sol.core.nl(params, guess, An, Bn, Cn, def)
      out$p <- sol$x[1:n]
      out$d <- sol$x[n+1:n]
          # Reorder the solution elements
      test.p <- out$p
      min.p <- p.init.d( params, out$p, out$d, An, Bn, Cn, def )
      guess[,1] <- min.p
          # Check (roughly) for globality
    }else if(st.which.sol == 'core.nl.i'){
      sol <- sol.core.nl.i(params, guess, An, Bn, Cn, def, i)
      out$p <- guess[,1]
      out$d <- guess[,2]
      out$p[i] <- sol$x[1]
      out$d[i] <- sol$x[2]
          # Reorder the solution elements
      test.p <- out$p[i]
      min.p <- p.init.d( params, out$p, out$d, An, Bn, Cn, def )[i]
      guess[i,1] <- min.p
      guess[i,2] <- out$d[i]
          # Check (roughly) for globality in state i
    }else{
      stop('Solution method not recognized')
    }

  }
  if( it > 1 ){
    message( "Finding global solution required ", it, " iterations" )
    message( "  Global solution satisfied up to ", max(abs(test.p-min.p)) )
  }

  out$err <- matrix( sol$fvec, nrow = n, ncol = 2 )
      # The errors
  out$accy <- max( abs( out$err ) ) < 1e-08
      # Boolean for accuracy of the solution
  out$An <- An
  out$Bn <- Bn
  out$Cn <- Cn
  out$def <- def
      # Record inputs
  out$sol <- sol
      # Extra solution info
  ## NEED TO CHECK FOR LOCAL/GLOBAL CONCAVITY
  return(out)
}

sol.core.nl <- function(params, init.guess, An, Bn, Cn, def){
# The core solution routine using nleqslv w/o derivatives.  Inital guess is a matrix.

  x0 <- c( init.guess )
      # Reformat as a vector
  n <- nrow(init.guess)
      # The number of states
  fn <- function(x){
      # The function for which we want to find the root
    p <- pmax( x[1:n], 0 )
    d <- x[n+1:n]
        # Extract the values
    z.2 <- zed_2( p, d, params, An, Bn, Cn, def )
        # The vector of: the implied value of p and the gradient wrt p
    return( c( z.2 - cbind( p, rep(1,n) ) ) )
        # Becuase we want the slope to be unity at the fixed point
  }

  control <- list( maxit=1000, allowSingular =TRUE, ftol=1e-10, xtol=1e-12 )
  sol <- nleqslv( x0, fn, control = control )
      # The solution object
  return( sol )
}

sol.core.nl.i <- function(params, init.guess, An, Bn, Cn, def, i){
# The core solution routine using nleqslv w/o derivatives for state i only.
# Inital guess is a matrix.

  x0 <- init.guess[i,]
      # Extract x0
  n <- nrow(init.guess)
      # The number of states
  fn <- function(x){
  # The function for which we want to find the root
    p <- pmax( init.guess[,1], 0 )
    d <- init.guess[,2]
        # Extract the values for the unchanging states
    p[i] <- max( x[1], 0 )
    d[i] <- x[2]
        # Update dimension i
    z.2 <- zed_2( p, d, params, An, Bn, Cn, def )[i,]
        # The vector of: the implied value of p and the gradient wrt p
        # TO DO: ADD A DIMENSION-i VERSION OF ZED_2 FOR SPEED HERE
    return( z.2 - c( p[i], 1 ) )
        # Becuase we want the slope to be unity at the fixed point
  }

  control <- list( maxit=1000, allowSingular =TRUE, ftol=1e-10, xtol=1e-12 )
  sol <- nleqslv( x0, fn, control = control )
      # The solution object
  return( sol )
}

sol.core.nloptr.i <- function(params, init.guess, An, Bn, Cn, def, i){
# The core solution routine using constrained optimization derivatives for state i only.

}
