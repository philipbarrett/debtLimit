#####################################################################
# outerSol.R
#
# Contains the functions for solving the outer part of the solution alogorithm
# 26mar2017
# Philip Barrett, Washington DC
#####################################################################

outer.wrapper <- function(d.grid, e.grid, d.bar, params, An, Bn, Cn, Q.init=NULL,
                          D.prime.init=NULL, print_level=1 ){
# Wraps the outer loop for solution.

  ## Set default prices and starting continuation debts
  if( is.null(Q.init) ){
    Q <- Q_init( d.grid, sol.w$d, params$R )
  }
  if( is.null( D.prime.init ) ){
    D.prime.init <- matrix(0)
    D.prime.flag <- FALSE
  }else{
    D.prime.flag <- TRUE
  }

  ## Set default convergence parameters
  if(is.null(params$tol)) params$tol <- 1e-04
  if(is.null(params$maxit)) params$maxit <- 50
  if(is.null(params$q.tol)) params$q.tol <- 1e-05
  if(is.null(params$q.maxit)) params$q.maxit <- 50
  if(is.null(params$d.tol)) params$d.tol <- 1e-04
  if(is.null(params$d.maxit)) params$d.maxit <- 30

  ## Now solve the model
  ZZ.QQ <- ziter( 0 * Q, d.bar, Q, d.grid, e.grid, params$tri, D.prime.init, D.prime.flag,
               params, An, Cn, print_level, params$tol, params$maxit, params$q.tol,
               params$q.maxit, params$d.tol, params$d.maxit, Q_out = TRUE )
      # Iterate to find the equilibrium default probabiltiies
  n <- length(params$R)
  ZZ <- ZZ.QQ[1:n,]
  QQ <- ZZ.QQ[n+1:n,]
  QQ.2 <- q_hat_mat( ZZ, d.bar, QQ, QQ, d.grid, params$R,params$G, params$lambda,
                   params$phi, e.grid, params$v.s.coeff, params$tri,
                   D.prime.init, D.prime.flag, params$trans, print_level+1 )
      # The associated price
  DD <- d_prime_mat( d.bar, QQ, QQ, d.grid, params$G, params$lambda, e.grid, params$trans,
                     params$v.s.coeff, params$tri, D.prime.init, D.prime.flag, print_level - 1 )
      #  Expected continuation debt
  QE <- qe_mat( d.bar, QQ, QQ, d.grid, params$G, params$lambda, e.grid, params$trans,
                params$v.s.coeff, params$tri, D.prime.init, D.prime.flag, print_level - 1 )
      # Expected continuation price
  return( list( d.grid=d.grid, e.grid=e.grid, d.bar=d.bar, P=ZZ, Q=QQ, D.prime=DD, QE=QE ) )
}

outer.err <- function( sol, params, An, Cn ){
# Computes the error on the outer solution
  P.err <- ziter( sol$P, sol$d.bar, sol$Q, sol$d.grid, sol$e.grid, params$tri,
                  matrix(0,1,1), FALSE, params, An, Cn, print_level=3, maxit=1 ) - sol$P
      # Error on the outermost fixed point
  Q.new <- q_hat_mat( sol$P, sol$d.bar, sol$Q, sol$Q, sol$d.grid, params$R, params$G,
                      params$lambda, params$phi, sol$e.grid, params$v.s.coeff, params$tri,
                      matrix(0,1,1), FALSE, params$trans, print_level=2, maxit=1 )
  Q.err <- Q.new - sol$Q
      # Error on prices
  D.prime.err <- d_prime_mat( sol$d.bar, sol$Q, sol$Q, sol$d.grid, params$G, params$lambda,
                              sol$e.grid, params$trans, params$v.s.coeff, params$tri,
                              matrix(0,1,1), FALSE, print_level=1 ) - sol$D.prime
      # Error on average debt.  Here we *do* need to allow for iteration becuase
      # we do not retain the starting guess of
  return( list( P=P.err, Q=Q.err, D.prime=D.prime.err ) )
}
