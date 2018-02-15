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
arma::mat ziter( arma::mat P, arma::vec d_bar, arma::mat QHat, //arma::mat Q,
                     arma::vec d_grid, arma::vec e_grid,
                     bool tri, arma::mat D_prime_0, bool D_prime_0_flag,
                     List params, int print_level=0, double tol=1e-04, int maxit=100,
                     double q_tol=1e-04, int q_maxit=50, double d_tol=1e-05, int d_maxit=20,
                     bool Q_out=false ){
// Conducts the iteration on the big P until we get convergence
  int n = P.n_rows ;
  int m = P.n_cols ;
  mat def =zeros(1,1) ;
      // Initialization

  mat trans = params["trans"] ;
  vec R = params["R"] ;
  vec G = params["G"] ;
  vec coeff = params["v.s.coeff"] ;
  vec shift = params["s.shift"] ;
  double phi = params["phi"] ;
  double lambda = params["lambda"] ;
      // Extraction

  double diff = 2 * tol ;
  double q_diff = 2 * tol ;
  double q_1_step_diff = 2 * tol ;
  int it = 0 ;
  uvec d_bar_idx(n) ;
  vec Cn = zeros(n) ;
  mat P_old = P ;
  mat P_new = P ;
  mat QHat_old = QHat ;
  mat QHat_new = QHat ;
  mat QHat_1_step = QHat ;
      // Initialize loop variables
  for( int i = 0 ; i < n ; i ++ ){
    uvec temp = find( abs( d_grid - d_bar(i) ) <= 0  ) ;
    d_bar_idx(i) = temp(temp.n_elem-1) ;
        // Choose Cn to be the debt price *just* below the default threshold
  }
  if( print_level > 0 ){
    // Rcout << "\nd_bar_idx: " << conv_to<rowvec>::from(d_bar_idx) << std::endl << std::endl ;
    Rcout << std::setw(10) << "Iteration" << std::setw(10) << "diff" << std::setw(10) <<
      "q_diff" << std::setw(10) << "q_1_diff" << std::endl ;
  }

  while( ( diff > tol || q_1_step_diff  > tol ) && it < maxit ){
    it++ ;
    QHat_new = .75 * q_hat_mat( P_old, d_bar, QHat_old, QHat_old, d_grid, R, G, shift, lambda, phi, e_grid,
                              coeff, tri, D_prime_0, D_prime_0_flag, trans, print_level-1,
                              q_tol, q_maxit, d_tol, d_maxit ) +
                              .25 *   QHat_old;
    // Now calculate:
    //   - Cn = vector of prices at default boundary
    //   - mC = matrix of QE given the current price matrix
    for( int i = 0 ; i < n ; i ++ ){
    //   // vec temp = conv_to<vec>::from(QHat_new.row(i)) ;
    //   // uvec idx = find( temp > 0 ) ;
    //   // vec temp2 = temp.elem( idx ) ;
      Cn(i) = QHat_new(i,d_bar_idx(i)) ;
          // Make Cn equal to the vector of last non-zero elements
    }

    for( int j = 0 ; j < m ; j++ ){
      vec An = q_e( d_grid(j), d_bar, QHat_new.col(j), QHat_new, d_grid, G, shift, lambda, e_grid,
                    coeff, tri, D_prime_0, D_prime_0_flag, trans, print_level-1, d_tol, d_maxit ) ;
      P_new.col(j) = .75 * zed_d( P_old.col(j), d_grid(j)*ones(n), d_bar, params, An, Cn, def ) +
                        .25 * P_old.col(j);
    }
        // Update P.  Add gain to limit oscillation
    diff = abs( P_new - P_old ).max() / .75 ;
    q_diff = abs( QHat_new - QHat_old ).max() / .75 ;
        // The difference
    QHat_1_step = q_hat_mat( P_old, d_bar, QHat_old, QHat_old, d_grid, R, G, shift, lambda, phi, e_grid,
                          coeff, tri, D_prime_0, D_prime_0_flag, trans, print_level-1,
                          q_tol, 1, d_tol, d_maxit );
    q_1_step_diff = abs( QHat_1_step - QHat_old ).max() ;
                           // ( )The one-step difference - QHat_old ).max() ;
        // The one-step difference
    if( print_level > 0 ){
      Rcout << std::setw(10) <<  it << std::setprecision(3) << std::setw(10) << diff <<
        std::setw(10) << q_diff << std::setw(10) << q_1_step_diff << std::endl ;
    }
    P_old = P_new ;
    QHat_old = QHat_new ;
  }
  if( Q_out ){
    return join_vert( P_new, QHat_new ) ;
  }
  return P_new ;
}
