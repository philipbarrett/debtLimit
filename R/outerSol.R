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

  if( is.null(Q.init) ){
    Q <- Q_init( d.grid, sol.w$d, params$R )
  }
  if( is.null( D.prime.init ) ){
    D.prime.init <- matrix(0)
    D.prime.flag <- FALSE
  }else{
    D.prime.flag <- TRUE
  }

  ZZ <- ziter( 0 * Q, d.bar, Q, Q, d.grid, e.grid, params$tri, D.prime.init, D.prime.flag,
               params, An, Cn, print_level )
      # Iterate to find the equilibrium default probabiltiies
  QQ <- q_hat_mat( ZZ, d.bar, Q, Q, d.grid, params$R,params$G, params$lambda,
                   params$phi, e.grid, params$v.s.coeff, params$tri,
                   D.prime.init, D.prime.flag, params$trans, print_level - 1 )
      # The associated price
  DD <- d_prime_mat( d.bar, QQ, QQ, d.grid, params$G, params$lambda, e.grid, params$trans,
                     params$v.s.coeff, params$tri, D.prime.init, D.prime.flag, print_level - 1 )
      #  Expected continuation debt
  QE <- qe_mat( d.bar, QQ, QQ, d.grid, params$G, params$lambda, e.grid, params$trans,
                params$v.s.coeff, params$tri, D.prime.init, D.prime.flag, print_level - 1 )
      # Expected continuation price
  return( list( d.grid=d.grid, d.bar=d.bar, P=ZZ, Q=QQ, D.prime=DD, QE=QE ) )
}

outer.err <- function( sol ){
# Computes the error on the outer solution

}
