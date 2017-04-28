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
  double inc = 1e-04 ;                    // The increment
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
    if( this_z > 1 - 1e-06 ){
      this_p = 1 ;
    }
  }

  z_temp = zed_2( p, d, params, An, Bn, Cn, def ) ;
  double z_0 = z_temp(i,0) ;
      // The value of z at p=0
  if( opt_z_minus_p > z_0 ){
    opt_p = 0 ;
      // If p=0 gives a superior solution, use that instead.
  }else{

  }
  return opt_p ;
}

// [[Rcpp::export]]
arma::vec p_fp( List params, arma::vec p, arma::vec d, arma::vec An,
                   arma::vec Bn, arma::vec Cn, arma::mat def,
                   int maxit=200, double tol=1e-06 ){
// Finds a vector of fixed points given d near p
  int n = d.n_elem ;
      // Number of staes
  p = p + 5e-04 ;
  vec p_old = p ;
      // Initialize p
  double p_inc = .1 / maxit ;
      // Increment for p
  int i = 0 ;
      // Itertion counter
  double diff = 2 * tol ;
      // Initialize the difference
  while( diff > tol & i < maxit ){
    i++ ;
        // Increment loop counter
    p = zed( p_old, d, params, An, Cn, def ) ;
        // Update p
    diff = max( abs( p - p_old ) ) ;
        // Update difference
    p_old = p ;
  }
  return p ;
}

// [[Rcpp::export]]
double d_init_p_i( List params, arma::vec p, arma::vec d, arma::vec An,
                   arma::vec Bn, arma::vec Cn, arma::mat def, int i,
                   double d_max, double d_min=0, int maxit=100, int print_level=0,
                   int max_outer=10, int i_outer=0, double d_step_0=-10 ){
// Computes the initial guess of d given p by reducing d until we find a change
// of sign in zed and then iterating to find z~=p
  int n = p.n_elem ;
  double d_step = d_step_0 ;
      // The initial step length of d
  double tol = 1e-08 ;
      // The initial tolerance
  double err = tol*2 ;
  int it = 0 ;
      // Iteration controls
  double d_opt = d_max ;
  double d_guess = d_max ;
  double d_opt_err = 2 ;
  double z_guess = 2 ;
  double z_guess_err = 2 ;
  double z_guess_err_old = 2 ;
  double p_i = p(i) ;
  vec z_temp = zeros(n) ;
      // Measuring the error on zed and recording the best guess

  if( print_level > 0 ){
    Rcout << "\n** i_outer = " << i_outer << " **" << std::endl ;
    Rcout << std::setw(10) << "Iteration" << std::setw(10) << "d_guess" <<
      std::setw(10) << "p_i" << std::setw(10) << "z_guess" << std::setw(10) <<
      "z_guess_err" << std::setw(10) << "d_opt" << std::setw(10) << "d_opt_err" <<
        std::setw(10) << "d_step" << std::endl ;
    Rcout << std::setw(10) <<  it << std::setprecision(4) << std::setw(10) << d_guess <<
      std::setw(10) << p_i << std::setw(10) << z_guess <<
      std::setw(10) << z_guess_err << std::setw(10) << d_opt << std::setw(10) << d_opt_err <<
        std::setw(10) << d_step << std::endl ;
  }

  while( err > tol && it < maxit && d_guess > d_min ){
    it++ ;
        // Increment the counter
    d(i) = d_guess ;
        // Update the vector d
    z_temp = zed( p, d, params, An, Cn, def, print_level - 1 ) ;
        // Evaluate z
    z_guess =z_temp(i) ;
        // Extract i-th element
    z_guess_err = z_guess - p_i ;
        // Computue the error
    if( z_guess_err * z_guess_err_old < 0 ){
      d_step = - d_step / 3 ;
        // If the error sign flips, start searching in the other direction
    }
    if( fabs(z_guess_err) < fabs( d_opt_err ) ){
      d_opt = d_guess ;
      d_opt_err = z_guess_err ;
        // If the error is smaller, record it as optimal
    }
    if( print_level > 0 ){
      Rcout << std::setw(10) <<  it << std::setprecision(4) << std::setw(10) << d_guess <<
        std::setw(10) << p_i << std::setw(10) << z_guess <<
        std::setw(10) << z_guess_err << std::setw(10) << d_opt << std::setw(10) << d_opt_err <<
          std::setw(10) << d_step << std::endl ;
    }
    if( fabs(d_opt_err) < tol ){
      return d_opt ;
        // If error small enough, return
    }
    d_guess = std::max( 0.0, d_guess + d_step ) ;
        // Increment guess
    z_guess_err_old = z_guess_err ;
        // Store old guess
  }
  if( i_outer < max_outer ){
    return d_init_p_i( params, p, d, An, Bn, Cn, def, i, d_opt + fabs(d_step_0), d_opt - fabs(d_step_0),
                       maxit, print_level, max_outer, i_outer+1, d_step_0 / 3 ) ;
        // If no return then all recursively on interval around best point
  }
  return d_opt ;
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

// [[Rcpp::export]]
arma::vec d_init_p( List params, arma::vec p, arma::vec d, arma::vec An,
                    arma::vec Bn, arma::vec Cn, arma::mat def,
                    double d_max, int maxit=100, int print_level=0, int max_outer=10 ){
  // Vector version of p_init_d_i

  int n = p.n_elem ;
      // number of states
  vec out = zeros(n) ;
      // Initialize output
  for( int i = 0 ; i < n ; i++ ){
    out[i] = d_init_p_i(params, p, d, An, Bn, Cn, def, i, d_max, 0, maxit, print_level,
                        max_outer, 0) ;
  }
  return out ;
}
