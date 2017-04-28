/***********************************************************************************
 * eds.cpp
 *
 * Code to compute the epsilon-distinguishable sets as described in
 * Maliars & Maliars (2015), Quantitative Economics
 *
 * 14jan2016
 * Philip Barrett, Chicago
 *
 ***********************************************************************************/


#include "eds.hpp"

// [[Rcpp::export]]
arma::mat keep_mat( arma::mat & X, arma::vec keep ){
// Returns a submatrix of X where the rows are selected if the corresponding
// entry of the keep vector is 1, and otherwise not.

  int n_rows = X.n_rows ;
      // The number of rows of X
  int n_z = sum( keep ) ;
      // The number of points to keep
  mat Z = zeros( n_z, X.n_cols ) ;
      // Initialize the return matrix
  int i_z = 0 ;
      // Counter for rows of z
  for( int i= 0 ; i < n_rows ; i++ ){
    // Loops over rows of X
    if( keep(i) == 1 ){
      // If i is a keeper, then keep it!
      Z.row( i_z ) = X.row( i ) ;
          // Add to Z
      i_z++ ;
          // Increment Z
    }
  }

  return Z ;
}

// [[Rcpp::export]]
arma::mat pcd_scaling( arma::mat & X ) {
// Reutrns the mean-zero, spherical unit-variance equivalent of a cloud of
// points via PC decomposition

  int i_select = 250 ;
      // The number of rows to use in the pcd part
  rowvec X_mean = mean( X, 0 ) ;
      // The vector of column means
  mat X_center = ( X - ones( X.n_rows ) * X_mean ) ;
      // Centering the X matrix on its mean
  ivec select = randi<ivec>(i_select, distr_param(0, X.n_rows - 1 ) ) ;
      // The subset of X to use in the PCD rescaling
  mat X_for_pc = zeros( i_select, X.n_cols ) ;
  for( int i = 0 ; i < i_select ; i++ )
    X_for_pc.row(i) = X_center.row( select(i) ) ;
      // Create the random subset of rows for the pcd part
  mat X_pc = princomp( X_for_pc ) ;
      // The principal components of X
  mat X_pcd = X_center * X_pc ;
      // The principal component projection of X
  rowvec X_sd = stddev( X_pcd, 0, 0 ) ;
      // The vector of standard deviations
  mat X_norm = X_pcd / ( ones( X.n_rows ) * X_sd ) ;
      // Rescale to make spherical unit variance

  return X_norm ;

}

// [[Rcpp::export]]
arma::vec eds_keep( arma::mat & X_norm, arma::vec & eps ){
// Computes the vector of rows to keep in the EDS algorithm
  int n_rows = X_norm.n_rows ;
      // The number of rows of X
  vec keep = ones( n_rows ) ;
      // Initiate the loop indicator
  double dist = 0 ;
      // The Euclidean distance
  for ( int i = 0 ; i < n_rows - 1 ; i++ ){
    // Loop over rows of X_norm
    if( keep(i) == 1 ){
      // If the point is not a keeper, no need to look for nearest neighbours
      for( int j = i + 1 ; j < n_rows ; j++ ){
        // Loop over other rows
        dist = norm( X_norm.row(i) - X_norm.row(j), 2 ) ;
            // Compute the Euclidean distance
        if( dist < eps(i) ) keep(j) = 0 ;
            // If the distance is small set the keeping indicator to zero
      }
    }
  }
  return keep ;
}

// [[Rcpp::export]]
arma::mat p_eps( arma::mat & X, arma::vec & eps ) {
// The EDS algorithm for a vector of distances epsilon

  mat X_norm = pcd_scaling( X ) ;
      // The PCD rescaling of X (removes units)

  /** 1. Calculate the rows to retain **/
  vec keep = eds_keep( X_norm, eps ) ;

  /** 2. Construct a return matrix with those rows **/
  mat Z = keep_mat( X, keep ) ;

  return Z ;
}

// [[Rcpp::export]]
arma::mat p_eps_const( arma::mat & X, double eps ) {
// The EDS algorithm for fixed epsilon

  vec v_eps = eps * ones( X.n_rows ) ;
      // Create the constant vector of epsilons
  return p_eps( X, v_eps ) ;
}

// [[Rcpp::export]]
arma::vec normal_kernel_density( arma::mat & X, double h=0 ) {
// Coputes a normal kernal density given a bandwidth h

  int d = X.n_cols ;
      // The dimension of the state space
  int n = X.n_rows ;
      // The number of data points
  if( h == 0 )
    h = pow( n, -1.0 / ( d + 4.0 ) ) ;
      // Standard setting for the bandwidth
  mat D2 = zeros( n, n ) ;
      // The matrix of pairs of square distances
  double dist2 = 0 ;
      // Temporary (square) distance placeholder

  for( int i = 1 ; i < n ; i++ ){
    // Loop over rows
    for( int j = 0 ; j < i ; j++ ){
      // Loop over columns
      dist2 = pow( norm( X.row(i) - X.row(j), 2 ), 2.0 ) ;
          // The square distance
      D2( i, j ) = dist2 ;
      D2( j, i ) = dist2 ;
          // Fill in the matrix
    }
  }

  mat expD2 = exp( - D2 / ( 2.0 * h * h ) ) ;
      // The exponent term
  double scale_factor = ( n * pow( h, d ) * pow( 2.0 * M_PI, d / 2.0 ) ) ;
  vec g = sum( expD2, 1 ) / scale_factor ;
      // The density vector

  return g ;
}

// [[Rcpp::export]]
arma::vec almost_ergodic_indices( arma::mat & X, double delta, double h=0 ){
// Removes the points with cumulative probability mass less than delta

  vec g = normal_kernel_density( X, h=0 ) ;
      // The kernel density
  double sum_g = sum(g) ;
      // The sum of the densities
  uvec indices = sort_index(g);
      // The sorted indices
  vec keep = ones( X.n_rows ) ;
      // The vector of keepers
  double cumprob = 0 ;
      // Initialize cumulative probability
  int i = 0 ;
      // Initialize the row counter
  while( cumprob < delta & i < X.n_rows ){
    cumprob = cumprob + g( indices( i ) ) / sum_g ;
        // The cumulative probability so far
    keep( indices( i ) ) = 0 ;
        // Eliminate this row
    i++ ;
        // Increment counter
  }
  return keep ;
}

// [[Rcpp::export]]
arma::mat almost_ergodic( arma::mat & X, double delta, double h=0 ){
// Removes the points with cumulative probability mass less than delta

  vec keep = almost_ergodic_indices( X, delta, h ) ;
      // The indices of the retained rows
  mat Z = keep_mat( X, keep ) ;
      // Remove the omitted rows
  return Z ;
}

// [[Rcpp::export]]
arma::vec p_eps_cheap_idx( arma::mat & X, arma::vec & eps, double delta, double h=0 ) {
// Indices of the cheap EDS algorithm with removal of low-probability points

  mat X_norm = pcd_scaling( X ) ;
      // The PCD rescaling of X (removes units)

  /** 1. Calculate the rows to retain **/
  vec keep = eds_keep( X_norm, eps ) ;

  /** 2. Construct a PC matrix with those rows **/
  mat Z_norm = keep_mat( X_norm, keep ) ;
  mat Z = keep_mat( X, keep ) ;

  /** 3. Construct a return matrix with those rows **/
  vec erg_indices_Z = almost_ergodic_indices( Z_norm, delta, h ) ;
      // The indices of the normalized returned points

  /** 4. Now loop over the original indices to remove the low-probability points **/
  vec erg_indices = keep ;
      // Initialize the
  int j = 0 ;
      // The Z index
  for( int i = 0 ; i < X.n_rows ; i++ ){
    if( keep(i) == 1 ){
      erg_indices(i) = erg_indices_Z(j) ;
          // Copy from the probability-reduced set
      j++ ;
          // Increment the Z counter
    }
  }

  return erg_indices ;
      // The indices of X
}

// [[Rcpp::export]]
arma::mat p_eps_cheap( arma::mat & X, arma::vec & eps, double delta, double h=0 ){
// The cheap EDS algorithm with removal of low-probability points

  mat X_norm = pcd_scaling( X ) ;
      // The PCD rescaling of X (removes units)

  /** 1. Calculate the rows to retain **/
  vec keep = eds_keep( X_norm, eps ) ;

  /** 2. Construct a PC matrix with those rows **/
  mat Z_norm = keep_mat( X_norm, keep ) ;
  mat Z = keep_mat( X, keep ) ;

  /** 3. Construct a return matrix with those rows **/
  vec erg_indices = almost_ergodic_indices( Z_norm, delta, h ) ;
      // The indices of the normalized returned points
  mat out = keep_mat( Z, erg_indices ) ;
      // Retain only the indices in the almost-ergodic set

  return out ;
}

// [[Rcpp::export]]
arma::vec p_eps_cheap_const_idx( arma::mat & X, double eps, double delta, double h=0 ) {
// Indices wth constant epsilon
  vec v_eps = eps * ones( X.n_rows ) ;
      // Create the constant vector of epsilons
  return p_eps_cheap_idx( X, v_eps, delta, h ) ;
}

// [[Rcpp::export]]
arma::mat p_eps_cheap_const( arma::mat & X, double eps, double delta, double h=0 ) {
// With constant epsilon
  vec v_eps = eps * ones( X.n_rows ) ;
      // Create the constant vector of epsilons
  return p_eps_cheap( X, v_eps, delta, h ) ;
}





