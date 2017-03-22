#####################################################################
# sol.R
#
# Contains the functions for solving the model
# 13mar2017
# Philip Barrett, Washington DC
#####################################################################

sol.wrapper <- function( params, init.guess=NULL,
                         An=c(0), Bn=c(0), Cn=c(0), def=matrix(0), i=1 ){
# Wraps and processes the solution output from the nonlinear solver

  n <- length(params$R)
  if(is.null(init.guess)){
    d.init <- rep( min( sol.nonstoch(params) ), n )
    p.init <- p_init_d( params, rep(0,n), d.init, An, Bn, Cn, def )
    init.guess <- cbind( p.init, d.init )
  }
  maxit <- if(is.null(params$it)) 10 else params$it
  tol <- if(is.null(params$tol)) 1e-02 else params$tol
      # Convergence requirements
  it <- 0
  err <- 2 * tol
  sol.i <- list( p=init.guess[,1], d=init.guess[,2] )
      # Initiate the candidate solution and the loop variables
  while( ( err > tol ) && it < maxit ){
    it <- it + 1
    message("it = ", it )
    for( i in 1:n ){
      p.guess <- sol.i$p
      d.guess <- sol.i$d
      p.guess[i] <- p_init_d_i( params, sol.i$p, sol.i$d, An, Bn, Cn, def, i-1 )
      sol.i <- sol.core.global( params, cbind( p.guess, d.guess), 'core.nl.i',
                                An, Bn, Cn, def, i )
    }
    err <- max( abs( zed_2( sol.i$p, sol.i$d, params, An, Bn, Cn, def ) -
                       c( sol.i$p, rep(1,n) ) ) )
    message("   err = ", round(err,4) )
  }

  sol <- sol.core.global( params, cbind( sol.i$p, sol.i$d ), 'core.nl', An, Bn, Cn, def )
      # The nonlinear solution
  return(sol)
}

sol.core.global <- function( params, init.guess, st.which.sol, An, Bn, Cn, def, i=1 ){
# Returns a global solution in either the i-th state or overall

  n <- length(params$R)
  guess <- init.guess
      # Set up
  maxit <- if(is.null(params$gloabl.it)) 10 else params$global.it
  tol <- if(is.null(params$global.tol)) 1e-04 else params$global.tol
      # Convergence requirements
  it <- 0
  cand.p <- 1
  p.global <- 0
      # Initiate the canidate and global minima for p

  while( any( abs( cand.p - p.global ) > tol ) && it < maxit ){
    it <- it + 1
    if(st.which.sol == 'core.nl'){
      sol <- sol.core.nl(params, guess, An, Bn, Cn, def)
      cand.p <- sol$x[1:n]
      d <- sol$x[n+1:n]
          # Reorder the solution elements
      p.global <- p_init_d( params, cand.p, d, An, Bn, Cn, def )
          # Check for approximate globality
      z.global <- zed( p.global, d, params, An, Cn, def )
          # The values of z at the approximate global solutions
      p.global[ z.global - p.global >= sol$fvec[1:n] ] <-
                              cand.p[ z.global - p.global >= sol$fvec[1:n] ]
          # If the approximate global minimum is inferior to the candidate
          # solution, then replace
      guess <- cbind( p.global, d )
          # Update the guess
    }else if(st.which.sol == 'core.nl.i'){
      sol <- sol.core.nl.i(params, guess, An, Bn, Cn, def, i)
      cand.p <- guess[,1]
      d <- guess[,2]
      cand.p[i] <- sol$x[1]
      d[i] <- sol$x[2]
          # Reorder the solution elements, extracting the i'th elements where
          # required
      p.global.i <- p_init_d_i( params, cand.p, d, An, Bn, Cn, def, i-1 )
      p.global <- cand.p
      p.global[i] <- p.global.i
          # The approximate global solution in dimension i
      z.global.i <- zed( p.global, d, params, An, Cn, def )[i]
          # The values of z at the approximate global solutions
      p.global.i <- if( z.global.i - p.global.i >= sol$fvec[1] ) cand.p[i] else p.global.i
          # Replace the glbal minimum with the candidate if it is inferior
      p.global[i] <- p.global.i
          # Update the global trackced in the loop
      guess[i,] <- c( p.global.i, d[i] )
          # Update the guess
    }else{
      stop('Solution method not recognized')
    }
    # if( all( abs(sol$fvec) < 1e-08 ) ) break()
        # Universal escape
  }

  out <- list( p=guess[,1], d=guess[,2], err = matrix( sol$fvec, nrow = n, ncol = 2 ),
               accy <- max( abs( sol$fvec ) ) < 1e-08, An=An, Bn=Bn, Cn=Cn, def=def,
               sol=sol, it=it, diff=abs( cand.p - p.global ) )
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
    if( p[i] < 0 ){
      return( z.2 - c( p[i], 1 ) + p[i] ^ 2 )
    }
    if( p[i] > 1 ){
      return( z.2 - c( p[i], 1 ) - ( p[i] - 1 ) ^ 2 )
    }
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
