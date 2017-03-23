/***********************************************************************************
 * ziter.cpp
 *
 * Computes the grand iteration over choice of P
 *
 * 23mar2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "ziter.hpp"
#include "qiter.hpp"
#include "zed.hpp"

// [[Rcpp::export]]
arma::mat ziter( arma::mat P, arma::vec d_bar, arma::mat QHat, arma::mat Q,
                     arma::vec d_grid, arma::vec e_grid,
                     bool tri, arma::mat D_prime_0, bool D_prime_0_flag,
                     List params, arma::vec An, arma::vec Cn,
                     int print_level=0, double tol=1e-04, int maxit=100,
                     double q_tol=1e-04, int q_maxit=50, double d_tol=1e-05, int d_maxit=20 ){
// Conducts the iteration on the big P until we get convergence
  int n = P.n_rows ;
  int m = P.n_cols ;
  mat def =zeros(1,1) ;
      // Initialization

  mat trans = params["trans"] ;
  vec R = params["R"] ;
  vec G = params["G"] ;
  vec coeff = params["v.s.coeff"] ;
  double phi = params["phi"] ;
  double lambda = params["lambda"] ;
      // Extraction

  double diff = 2 * tol ;
  int it = 0 ;
  mat P_old = P ;
  mat P_new = P ;
  mat QHat_old = QHat ;
      // Initialize loop variables
  while( diff > tol && it < maxit ){
    it++ ;
    mat QHat_new = q_hat_mat( P_old, d_bar, QHat_old, Q, d_grid, R, G, lambda, phi, e_grid,
                              coeff, tri, D_prime_0, D_prime_0_flag, trans, print_level-1,
                              q_tol, q_maxit=50, d_tol, d_maxit );
        // Update Q
    for( int j = 0 ; j < m ; j++ ){
      P_new.col(j) = zed_d( P_old.col(j), d_grid(j)*ones(n), d_bar, params, An, Cn, def ) ;
    }
        // Update P
    diff = abs( P_new - P_old ).max() ;
        // The difference
    if( print_level > 0 ){
      Rcout << "P iteration = " << it << ", diff = " << diff << std::endl ;
    }
    P_old = P_new ;
  }
  return P_new ;
}
