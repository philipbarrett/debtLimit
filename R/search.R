#####################################################################
# search.R
#
# Contains the functions for solving the model via iterative search metholds
# 27apr2017
# Philip Barrett, Washington DC
#####################################################################

sol.search.i <- function(params, init.guess, i, An, Bn, Cn, def, tol=1e-05, maxit=20,
                         print.level=0, p.tol=1e-06, p.maxit=100 ){
# Solves the state-i problem using search methods.  Turning down the tolerances
# here does not necessarily lead to more accurate solutions.

  it <- 0
  diff <- 2*tol
      # Loop variables
  p <- init.guess[,1]
  d.new <- d <- init.guess[,2]
      # Exctract p & d
  while( it < maxit & diff > tol ){
    it <- it + 1
    if( print.level > 0 ) message('\n** Iteration **', it )
    d.new[i] <- d_init_p_i(params, p, d, An, Bn, Cn, def, i-1, 1.2 * max(d), 0,
                           d_step_0 = - max(d) / 80, print_level = print.level-1 )
    p.new <- p_min_tanget_i( params, p, d.new, 1, i-1, An, Bn, Cn, def, p.maxit, p.tol,
                             print_level = print.level-1 )
        # Update d and p using i-1 for c++ conversion
    diff <- max(abs( p.new - p))
    if( print.level > 0 ) message('      diff = ', round(diff,6) )
        # Compute the difference
    z.2 <- zed_2( p.new, d.new, params, An, Bn, Cn, def )
        # The level and derivative of the function
    err <- max( abs( z.2[i,] - c(0,1) ) )
        # The error
    if( it > 1 ){
      if( err > err.old ) break()
    }
        # Break out if the difference is increasing - usually a bad sign.
    p <- p.new
    d <- d.new
    diff.old <- diff
    err.old <- err
        # Update everything
  }
  z.2 <- zed_2( p, d, params, An, Bn, Cn, def )
  err <- max( abs( z.2[i,] - c(0,1) ) )
      # Recompute these
  return( list( p=p, d=d, err=err, z.2=z.2, it=it ) )
}

sol.search <- function(params, init.guess=NULL,
                       An=c(0), Bn=c(0), Cn=c(0), def=matrix(0), plot.on=FALSE ){
# Computes an approximate solution using search methods
  if( is.null(init.guess) ){
    init.guess <- cbind( 1e-05, sol.nonstoch(params) )
  }

  nn <- length(params$R)
  maxit <- if(is.null(params$it)) 10 else params$it
  tol <- if(is.null(params$tol)) 1e-02 else params$tol

  it <- 0
  err <- 2 * tol
  p.guess <- init.guess[,1]
  d.guess <- init.guess[,2]
      # Initiate the candidate solution and the loop variables
  while( ( err > tol ) && it < maxit ){
    it <- it + 1
    message("it = ", it )
    err <- 0
    for( i in 1:nn ){
      guess <- cbind(p.guess, d.guess)
      search <- sol.search.i( params, guess, i, An, Bn, Cn, def )
      p.guess[i] <- search$p[i]
      d.guess[i] <- search$d[i]
      err <- max( abs( zed_2( p.guess, d.guess, params, An, Bn, Cn, def ) -
                         matrix( c(0,1), ncol=2, nrow=nn, byrow=TRUE ) ) )
          # Compute the error
      if( plot.on )
        plot.z(p.guess, d.guess, params, An, Bn, Cn, def,
               xlim=c(0,2*p.guess[i]), ylim=c(0,2*p.guess[i]) )
      if( err < tol ) i <- nn + 1
          # Break out if the error looks good
    }
    message( 'err = ', err )
  }
  err <- zed_2( p.guess, d.guess, params, An, Bn, Cn, def ) -
                     matrix( c(0,1), ncol=2, nrow=nn, byrow=TRUE )
      # Recompute the error

  return( list ( p=p.guess, d=d.guess, err=err, max.err=err ) )
}
