#####################################################################
# outerSol.R
#
# Contains the functions for solving the outer part of the solution alogorithm
# 26mar2017
# Philip Barrett, Washington DC
#####################################################################

outer.wrapper <- function(inner.sol, params, Q.init=NULL, D.prime.init=NULL,
                          P.init=NULL, d.grid.init=NULL, print_level=1 ){
# Wraps the outer loop for solution.

  # Create the grids
  if(is.null(params$e.pts)) params$e.pts <- 21
  if(is.null(params$d.pts)) params$d.pts <- 15
  if(is.null(params$x.mult)) params$x.mult <- 2
  e.grid <- e_grid_fn(2*params$surp.sd, params$e.pts, params$d.tri )
  d.grid <- d_grid_fn(inner.sol$d, params$surp.sd, x_sd_mult = params$x.mult, n_pts=params$d.pt )
  # Q <- Q_init( d.grid, inner.sol$d, params$R )

  ## Set default prices and starting continuation debts
  if( is.null(Q.init) ){
    Q <- Q_init( d.grid, inner.sol$d, params$R )
  }else{
    Q <- t( apply( Q.init, 1, function(q) approx( d.grid.init, q, d.grid, rule=2 )$y ) )
  }
  if( is.null(P.init) ){
    P <- 0 * Q
  }else{
    P <- t( apply( P.init, 1, function(q) approx( d.grid.init, q, d.grid, rule=2 )$y ) )
  }
  if( is.null( D.prime.init ) ){
    D.prime <- matrix(0)
    D.prime.flag <- FALSE
  }else{
    D.prime <- t( apply( D.prime.init, 1, function(d) approx( d.grid.init, d, d.grid )$y ) )
    D.prime.flag <- TRUE
  }
  d.bar <- inner.sol$d

  ## Set default convergence parameters
  if(is.null(params$tol)) params$tol <- 1e-04
  if(is.null(params$maxit)) params$maxit <- 50
  if(is.null(params$q.tol)) params$q.tol <- 1e-05
  if(is.null(params$q.maxit)) params$q.maxit <- 50
  if(is.null(params$d.tol)) params$d.tol <- 1e-04
  if(is.null(params$d.maxit)) params$d.maxit <- 60

  ## Now solve the model
  ZZ.QQ <- ziter( P, d.bar, Q, d.grid, e.grid, params$tri, D.prime, D.prime.flag,
                  params, print_level, params$tol, params$maxit, params$q.tol,
                  params$q.maxit, params$d.tol, params$d.maxit, Q_out = TRUE )
      # Iterate to find the equilibrium default probabiltiies
  n <- length(params$R)
  ZZ <- ZZ.QQ[1:n,]
  QQ <- ZZ.QQ[n+1:n,]
  # QQ.2 <- q_hat_mat( ZZ, d.bar, QQ, QQ, d.grid, params$R,params$G, params$lambda,
  #                  params$phi, e.grid, params$v.s.coeff, params$tri,
  #                  D.prime, D.prime.flag, params$trans, print_level+1 )
      # The associated price
  DD <- d_prime_mat( d.bar, QQ, QQ, d.grid, params$G, params$s.shift, params$lambda, e.grid, params$trans,
                     params$v.s.coeff, params$tri, D.prime, D.prime.flag, print_level - 1 )
      #  Expected continuation debt
  QE <- qe_mat( d.bar, QQ, QQ, d.grid, params$G, params$s.shift, params$lambda, e.grid, params$trans,
                params$v.s.coeff, params$tri, D.prime, D.prime.flag, print_level - 1 )
      # Expected continuation price
  p.bar <- sapply( 1:length(params$R), function(i) max(ZZ[i,ZZ[i,]<1-1e-06]) )
      # Implied default probability at the boundary
  return( list( d.grid=d.grid, e.grid=e.grid, d.bar=d.bar, P=ZZ, Q=QQ, D.prime=DD, QE=QE,
                p.bar=p.bar, p.init=inner.sol$p ) )
}

outer.err <- function( sol, params ){
# Computes the error on the outer solution
  P.err <- ziter( sol$P, sol$d.bar, sol$Q, sol$d.grid, sol$e.grid, params$tri,
                  matrix(0,1,1), FALSE, params, print_level=0, maxit=1 ) - sol$P
      # Error on the outermost fixed point
  Q.new <- q_hat_mat( sol$P, sol$d.bar, sol$Q, sol$Q, sol$d.grid, params$R, params$G, params$s.shift,
                      params$lambda, params$phi, sol$e.grid, params$v.s.coeff, params$tri,
                      matrix(0,1,1), FALSE, params$trans, print_level=0, maxit=1 )
  Q.err <- Q.new - sol$Q
      # Error on prices
  D.prime.err <- d_prime_mat( sol$d.bar, sol$Q, sol$Q, sol$d.grid, params$G, params$s.shift, params$lambda,
                              sol$e.grid, params$trans, params$v.s.coeff, params$tri,
                              matrix(0,1,1), FALSE, print_level=0 ) - sol$D.prime
      # Error on average debt.  Here we *do* need to allow for iteration becuase
      # we do not retain the starting guess of
  d.bar.l <- sol$d.grid[apply( sol$Q, 1, function(x) min(which(x==0)))-1]
  d.bar.u <- sol$d.grid[apply( sol$Q, 1, function(x) min(which(x==0)))]
      # Upper and lower bounds on the implied debt limit
  return( list( P=P.err, Q=Q.err, D.prime=D.prime.err, d.bar.range=cbind( d.bar.l, d.bar.u ),
                d.bar=sol$d.bar ) )
}

ABC <- function(sol, params){
# Computes An, Bn, Cn, d.bar updates, where:
#   - An = q^e(d.bar)
#   - Bn = d/dp( q^e(d.bar) )  (computed numerically)
#   - Cn = q(d.bar)

  nn <- length(params$R)
      # Number of states
  Q.pos.idx <- apply( sol$Q, 1, function(x) min(which(x==0)) ) - 1
      # Last point at which Q>0 for each state
  d.bar.apx <- sol$d.grid[ Q.pos.idx ]
      # Approximate d.bar
  An <- sapply( 1:nn, function(i) sol$QE[i,Q.pos.idx[i]])
      # The value of the continuation price at the last strictly positive price
  p.inc <- 1e-06
      # Incremental change to p
  Q.p <- q_hat_mat( sol$P - p.inc, sol$d.bar, sol$Q, sol$Q, sol$d.grid,
                    params$R, params$G, params$lambda, params$phi, sol$e.grid,
                    params$v.s.coeff, params$tri, 0*Q, FALSE, params$trans  )
      # The prices with slightly lower default probabilities (everywhere)
  QE.p <- qe_mat( sol$d.bar, Q.p, Q.p, sol$d.grid, params$G, params$lambda, sol$e.grid,
                  params$trans, params$v.s.coeff, params$tri, 0*Q, FALSE, 0 )
      # QE at ths slightly lower default probability
  QE.p.bar <- sapply( 1:nn, function(i) QE.p[ i, Q.pos.idx[i] ] )
      # QE evaluated at the boundary with slightly lower p
  Bn <- ( An - QE.p.bar ) / p.inc
      # The derivative w.r.t. p
  Cn <- sapply( 1:nn, function(i) sol$Q[i,Q.pos.idx[i]])
      # Return the last strictly positive value for the debt price
  return( list( An=An, Bn=Bn, Cn=Cn ) )
}
