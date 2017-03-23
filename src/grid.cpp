/***********************************************************************************
 * grid.cpp
 *
 * Defines the grids on which to calculate prices and default probabilities
 *
 * 22mar2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "grid.hpp"

// [[Rcpp::export]]
arma::vec d_grid_fn( arma::vec d, double x_sd, double x_sd_mult=1.1, int n_pts=7 ){
// Creates the vector of grid points for d

  int n = d.n_elem ;
      // Number of states
  int m = 1 + n * (n_pts + 1) ;
      // Number of points for output vector
  vec out = zeros(m) ;
      // Initialize output
  vec spread = zeros(n_pts+1) ;
      // The spread of points around each default boundary
  for( int i = 0 ; i < n_pts ; i++ ){
    double q = ( i + .5 ) / n_pts ;
    spread(i) = R::qnorm5( q, 0, x_sd * x_sd_mult, 1, 0 ) ;
        // The distance to spread out around the default boundaries
  }
  spread(n_pts) = 1e-04 ;
      // Add a point *just* above the default threshold
  for( int i = 0 ; i < n ; i++ ){
    out.subvec(i*(n_pts+1)+1,(i+1)*(n_pts+1)) = d(i) + spread ;
  }
      // Fill in the grid
  vec out_trim = out( find( out <= max(d) + 1 ) ) ;
      // Throw away points far above the maximum possible
  return sort(out_trim) ;
}

// [[Rcpp::export]]
arma::vec e_grid_fn( double x_sd, int n_pts=7 ){
// Vector of points for surplus shock
  vec out = zeros(n_pts) ;
      // Initialize output
  for( int i = 0 ; i < n_pts ; i++ ){
    double q = ( i + .5 ) / n_pts ;
    out(i) = R::qnorm5( q, 0, x_sd, 1, 0 ) ;
    // The distance to spread out around the default boundaries
  }
  return out ;
}

// [[Rcpp::export]]
arma::mat Q_init( arma::vec d_grid, arma::vec d, arma::vec R ){
// Creates the initial guess of the matrix of debt prices, Q
// Could really do this much much better
  int n = d.n_elem ;
  int m = d_grid.n_elem ;
      // Dimensions of state and debt grid
  mat Q = zeros(n,m) ;
      // Initialize the output
  for( int i = 0 ; i < n; i++ ){
    for( int j = 0 ; j < m ; j++ ){
      Q(i,j) = ( d_grid(j) <= d(i) ) ? ( 1 / R(i) ) : 0 ;
    }
  }
  return Q ;
}


