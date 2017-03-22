/***********************************************************************************
 * pre.cpp
 *
 * Defines various pre-solution functions
 *
 * 21feb2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "pre.hpp"
#include "zed.hpp"

// [[Rcpp::export]]
double p_init_d_i( List params, arma::vec p, arma::vec d, arma::vec An,
                   arma::vec Bn, arma::vec Cn, arma::mat def, int i ){
// Computes the initial guess of p given d as the minimizer of z-p (s.t. z'~=1) in dimension i

  double z_p_max = .1 ;                   // The ceiling on z_p for consideration
  double inc = .001 ;                     // The increment
  vec z_2 = zeros(2) ;                    // Container for the z, z_p output
  double this_p = 0 ;                     // Test value of p
  double this_z = 0 ;                     // Resulting value of p
  double this_z_minus_p = 0 ;             // Test value of z-p
  double this_z_p = 0 ;                   // Test value of z_p
  double opt_p = 0 ;                      // Current optimal p
  double opt_z = 1000 ;                   // Current optimal z. High to accept first trial.
  double opt_z_minus_p = opt_z - opt_p ;  // Current optimal z-p
  double opt_z_p_minus_1 = 1000 ;         // Current optimal z_p
  vec p_temp = p ;
  mat z_temp ;

  int it = 0 ;

  while( this_p < 1 && it < 1 / inc ){

    it++ ;
      // Guard against infinite loop
    p_temp[i] = this_p ;
        // Change p
    z_temp = zed_2( p_temp, d, params, An, Bn, Cn, def ) ;
    this_z = z_temp(i,0) ;
    this_z_p = z_temp(i,1) ;
        // Compute the values of z, z_p
    this_z_minus_p = this_z - this_p ;
          // The update for z-p

            // Rcout << "\nit = " << it << std::endl ;
            // Rcout << "this_p = " << this_p << std::endl ;
            // Rcout << "this_z = " << this_z << std::endl ;
            // Rcout << "this_z_minus_p = " << this_z_minus_p << std::endl ;
            // Rcout << "this_z_p = " << this_z_p << std::endl ;
            // Rcout << "opt_p = " << opt_p << std::endl ;
            // Rcout << "opt_z = " << opt_z << std::endl ;
            // Rcout << "opt_z_minus_p  = " << opt_z_minus_p << std::endl ;
            // Rcout << "opt_z_p_minus_1 = " << opt_z_p_minus_1 << std::endl ;

    if( this_z_minus_p < opt_z_minus_p ){
      // if( fabs( this_z_p - 1 ) < opt_z_p_minus_1 && this_z < 1 ){
      if( this_z < 1 ){
        opt_p = this_p ;
        opt_z_minus_p = this_z_minus_p ;
            // Superior point found!
        opt_z_p_minus_1 = fabs( this_z_p - 1 ) ;
      }
      this_p += inc ;
            // Increment p by a little
    }else{
      this_p = std::max( this_p + inc, this_z - opt_z_minus_p ) ;
          // Increment p.  Use the line going through the current optimum as a
          // way to "step ahead" the choice of p.  This relies on the
          // increasingness of z. NB: when the slope is very near 1, we need to
          // nudge this along a bit with the this_p + inc.
    }
    if( this_z == 1 ){
      this_p = 1 ;
    }
  }

  z_temp = zed_2( p, d, params, An, Bn, Cn, def ) ;
  double z_0 = z_temp(i,0) ;
      // The value of z at p=0
  if( opt_z_minus_p > z_0 ){
    opt_p = 0 ;
  }
      // If p=0 gives a superior solution, use that instead.
  return opt_p ;
}

// [[Rcpp::export]]
arma::vec p_init_d( List params, arma::vec p, arma::vec d, arma::vec An,
                    arma::vec Bn, arma::vec Cn, arma::mat def ){
  // Vector version of p_init_d_i

  int n = p.n_elem ;
      // number of states
  vec out = zeros(n) ;
      // Initialize output
  for( int i = 0 ; i < n ; i++ ){
    out[i] = p_init_d_i(params, p, d, An, Bn, Cn, def, i) ;
  }
  return out ;
}

