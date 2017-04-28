#####################################################################
# sol.R
#
# Contains the functions for solving the model
# 13mar2017
# Philip Barrett, Washington DC
#####################################################################

sol.wrapper <- function( params, init.guess=NULL,
                         An=c(0), Bn=c(0), Cn=c(0), def=matrix(0), i=1, plot.on=FALSE ){
# Wraps and processes the solution output from the nonlinear solver

  nn <- length(params$R)
  if(is.null(init.guess)){
    init.guess <- d.p.init.wrapper( params, An, Bn, Cn, def )
  }
  maxit <- if(is.null(params$it)) 10 else params$it
  tol <- if(is.null(params$tol)) 1e-02 else params$tol
      # Convergence requirements
  method <- if(is.null(params$inner.method)) 'err' else params$inner.method
      # The iteration method.  Can be 'all' or 'err
  it <- 0
  err <- 2 * tol
  sol.i <- list( p=init.guess[,1], d=init.guess[,2] )
  p.guess <- sol.i$p
  d.guess <- sol.i$d
      # Initiate the candidate solution and the loop variables
  while( ( err > tol ) && it < maxit ){
    it <- it + 1
    message("it = ", it )
    if( method == 'err' ){
      err.i <- zed_2( p.guess, d.guess, params, An, Bn, Cn, def ) - cbind( p.guess, rep(1,nn) )
          # Evaluate the error
      i <- which.max(abs(err.i[,1]))
          # Choose the state to improve based on probability error
      message('   improving state ', i )
      if(abs(err.i[i,1]) > 1e-02){
        p.guess[i] <- p_init_d_i( params, p.guess, d.guess, An, Bn, Cn, def, i-1 )
      }
          # Use iterative search if the initial guess is bad
      sol.i <- sol.core.global( params, cbind( p.guess, d.guess), 'core.nl.i',
                                An, Bn, Cn, def, i )
      p.guess <- sol.i$p
      d.guess <- sol.i$d
          # Solve and update guesses
    }
    if( method == 'all' ){
      for( i in 1:nn ){
        err.i <- zed_2( p.guess, d.guess, params, An, Bn, Cn, def ) - cbind( p.guess, rep(1,nn) )
        if(abs(err.i[i,1]) > 1e-02){
          # p.guess[i] <- p_init_d_i( params, p.guess, d.guess, An, Bn, Cn, def, i-1 )
          p.guess[i] <- 0
        }
            # Use iterative search if the initial guess is bad
        # sol.i <- sol.core.global( params, cbind( p.guess, d.guess), 'core.nl.i',
        #                           An, Bn, Cn, def, i )
        sol.i <- sol.core.local( params, cbind( p.guess, d.guess), 'core.nl.i',
                                  An, Bn, Cn, def, i )
        p.guess <- sol.i$p
        d.guess <- sol.i$d
        if(plot.on){
          plot.z(sol.i$p,sol.i$d,params,An, Bn, Cn,
                 xlim=c(0,min(1,p.guess[i]*3)), ylim=c(0,min(1,p.guess[i]*6)))
        }
      }
    }
    if( !( method %in% c('err', 'all') ) ){
      stop('Inner iteration method not recognized')
    }
    err <- max( abs( zed_2( p.guess, d.guess, params, An, Bn, Cn, def ) -
                       c( p.guess, rep(1,nn) ) ) )
    message("   err = ", round(err,4) )
  }

  sol <- sol.core.global( params, cbind( sol.i$p, sol.i$d ), 'core.nl', An, Bn, Cn, def )
      # The nonlinear solution
  return(sol)
}

sol.core.global <- function( params, init.guess, st.which.sol, An, Bn, Cn, def, i=1 ){
# Returns a global solution in either the i-th state or overall

  nn <- length(params$R)
  guess <- init.guess
      # Set up
  maxit <- if(is.null(params$gloabl.it)) 10 else params$global.it
  tol <- if(is.null(params$global.tol)) 1e-05 else params$global.tol
      # Convergence requirements
  it <- 0
  err <- 2 * tol
  cand.p <- 1
  p.global <- 0
      # Initiate the canidate and global minima for p

  while( ( any( abs( cand.p - p.global ) > tol ) || any( abs(err) > tol ) ) && it < maxit ){
    it <- it + 1
    if(st.which.sol == 'core.nl'){
      sol <- sol.core.nl(params, guess, An, Bn, Cn, def)
      cand.p <- sol$x[1:nn]
      d <- sol$x[nn+1:nn]
          # Reorder the solution elements
      p.global <- p_init_d( params, cand.p, d, An, Bn, Cn, def )
          # Check for approximate globality
      z.global <- zed( p.global, d, params, An, Cn, def, 0 )
          # The values of z at the approximate global solutions
      p.global[ z.global - p.global >= sol$fvec[1:nn] ] <-
                              cand.p[ z.global - p.global >= sol$fvec[1:nn] ]
          # If the approximate global minimum is inferior to the candidate
          # solution, then replace
      guess <- cbind( p.global, d )
          # Update the guess
      err <- matrix( sol$fvec, nrow = nn, ncol = 2 )
          # The error
    }else if(st.which.sol == 'core.nl.i'){
      sol <- sol.core.nl.i(params, guess, An, Bn, Cn, def, i)
      if(max(abs(sol$fvec))>1e-04){
        d.init <- d_init_p_i( params, rep(0,nn), rep(130,nn), An, Bn, Cn, def, i, 200 )
            # Find a better guess for d
        guess[i,'d.guess'] <- d.init
        sol <- sol.core.nl.i(params, guess, An, Bn, Cn, def, i)
      }
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
      z.global.i <- zed( p.global, d, params, An, Cn, def, 0 )[i]
          # The values of z at the approximate global solutions
      p.global.i <- if( z.global.i - p.global.i >= sol$fvec[1] ) cand.p[i] else p.global.i
          # Replace the global minimum with the candidate if it is inferior.  WHY?
      p.global[i] <- p.global.i
          # Update the global tracked in the loop
      guess[i,] <- c( p.global.i, d[i] )
          # Update the guess
      err <- zed_2( p.global, d, params, An, Bn, Cn, def )[i,] - c(p.global[i],1)
          # The error
    }else{
      stop('Solution method not recognized')
    }
  }

  if( any( abs(err) > tol ) ){
    # browser()
    warning('  sol.core.global failing with error = ', max(abs(err)))
  }
      # Universal escape

  out <- list( p=guess[,1], d=guess[,2], err = err,
               accy = max( abs( sol$fvec ) ) < 1e-08, An=An, Bn=Bn, Cn=Cn, def=def,
               sol=sol, it=it, diff=abs( cand.p - p.global ), params=params )
  return(out)
}

sol.core.local <- function( params, init.guess, st.which.sol, An, Bn, Cn, def, i=1 ){
# Returns a local solution in either the i-th state or overall

  # browser()
  nn <- length(params$R)
  guess <- init.guess
      # Set up
  maxit <- if(is.null(params$gloabl.it)) 50 else params$global.it
  tol <- if(is.null(params$global.tol)) 1e-05 else params$global.tol
      # Convergence requirements
  it <- 0
  err.opt <- err <- 2 * tol
      # Initiate the canidate and global minima for p

  while( any( abs(err) > tol ) && it < maxit ){
    it <- it + 1
    if( it > 1 ){
      guess[i,1] <- if(it==2) 0 else guess[i,1] + min( 1e-02, 1 / maxit )
    }
    if(st.which.sol == 'core.nl'){
      sol <- sol.core.nl(params, guess, An, Bn, Cn, def)
      guess[,1] <- sol$x[1:nn]
      guess[,1] <- sol$x[nn+1:nn]
          # Reorder the solution elements
      err <- matrix( sol$fvec, nrow = nn, ncol = 2 )
          # The error
    }else if(st.which.sol == 'core.nl.i'){
      sol <- sol.core.nl.i(params, guess, An, Bn, Cn, def, i)
      if(max(abs(sol$fvec))>1e-04){
        d.init <- d_init_p_i( params, pmax(guess[,'p.guess'], 1e-07), guess[,'d.guess'],
                              An, Bn, Cn, def, i-1, max( 200, max(guess[,'d.guess']) ) )
            # Find a better guess for d.  The increment for p is to guarantee
            # solution when p~=0
        guess[i,'d.guess'] <- d.init
        sol <- sol.core.nl.i(params, guess, An, Bn, Cn, def, i)
      }
      guess[i,1] <- sol$x[1]
      guess[i,2] <- sol$x[2]
          # Reorder the solution elements, extracting the i'th elements where
          # required
      err <- matrix( sol$fvec, nrow = nn, ncol = 2 )
          # The error
    }else{
      stop('Solution method not recognized')
    }
    if( all( abs( err ) <= abs( err.opt ) ) || it == 1 ){
      guess.opt <- guess
      err.opt <- err
      sol.opt <- sol
    }
        # Store the best iteration so far
  }

  if( any( abs(err) > tol ) ){
    # browser()
    warning('  sol.core.local failing with error = ', max(abs(err)))
  }
  # Universal escape

  out <- list( p=guess.opt[,1], d=guess.opt[,2], err = err.opt,
               accy = max( abs( sol.opt$fvec ) ) < 1e-08, An=An, Bn=Bn, Cn=Cn, def=def,
               sol=sol.opt, it=it, params=params )
  return(out)
}

sol.core.nl <- function(params, init.guess, An, Bn, Cn, def){
# The core solution routine using nleqslv w/o derivatives.  Inital guess is a matrix.

  x0 <- c( init.guess )
      # Reformat as a vector
  nn <- nrow(init.guess)
      # The number of states
  fn <- function(x){
      # The function for which we want to find the root
    p <- pmax( x[1:nn], 0 )
    d <- x[nn+1:nn]
        # Extract the values
    z.2 <- zed_2( p, d, params, An, Bn, Cn, def )
        # The vector of: the implied value of p and the gradient wrt p
    return( c( z.2 - cbind( p, rep(1,nn) ) ) )
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
  nn <- nrow(init.guess)
      # The number of states
  fn <- function(x){
  # The function for which we want to find the root
    p <- pmax( init.guess[,1], 0 )
    d <- init.guess[,2]
        # Extract the values for the unchanging states
    p[i] <- x[1]
    d[i] <- max( x[2], 1e-02 )
        # Update dimension i
    z.2 <- zed_2( p, d, params, An, Bn, Cn, def )[i,]
        # The vector of: the implied value of p and the gradient wrt p
    if( p[i] < 0 ){
      p[i] <- 0
      z.2 <- zed_2( p, d, params, An, Bn, Cn, def )[i,]
      return( z.2 - c( 0, 1 ) + x[1] ^ 2 )
    }
    if( p[i] > 1 ){
      p[i] <- 1
      z.2 <- zed_2( p, d, params, An, Bn, Cn, def )[i,]
      return( z.2 - c( 1, 1 ) + ( x[1] - 1 ) ^ 2 )
    }
    return( z.2 - c( p[i], 1 ) )
        # Because we want the slope to be unity at the fixed point
  }

  control <- list( maxit=1000, allowSingular =TRUE, ftol=1e-10, xtol=1e-12 )
  sol <- nleqslv( x0, fn, control = control )
      # The solution object
  return( sol )
}

sol.core.nloptr.i <- function(params, init.guess, An, Bn, Cn, def, i){
# The core solution routine using constrained optimization derivatives for state i only.

}
