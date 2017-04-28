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
// Computes the initial guess of p given d as the minimizer of z-p (s.t. z'~=1)
// in dimension i

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
      // Iteration counter
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
arma::vec p_ax_c_i( List params, arma::vec p, arma::vec d, double m, double c, int i,
                    arma::vec An, arma::vec Bn, arma::vec Cn, arma::mat def,
                    int maxit=200, double tol=1e-06, int print_level=0 ){
// Finds the smallest intersection **from above** of zed = m*p + c in direction i.
// If no intersection, retuns 0 if zed always below m*p+c.  Otherwise 1.
// NB: THIS IS SUPER INEFFICIENT. COMPUTES ALL ZED DIRECTIONS NEEDLESSLY
  int n = d.n_elem ;
      // Number of staes
  double z_i = 0;
  double p_i  = 0;
  vec z_new = p ;
  p(i) = p_i ;
      // Initialize p
  int it = 0 ;
      // Iteration counter
  double diff = 2 * tol ;
      // Initialize the difference

  if( print_level > 0 ){
    Rcout << std::setw(10) << "Iteration" << std::setw(10) << "p(i)" <<
      std::setw(10) << "z_i" << std::setw(10) << "p_i" << std::setw(10) <<
        "diff" << std::endl ;
    Rcout << std::setw(10) <<  it << std::setprecision(4) << std::setw(10) << p(i) <<
      std::setw(10) << z_i << std::setw(10) << p_i << std::setw(10) << diff << std::endl ;
  }

  while( diff > tol & it < maxit ){
    it ++ ;
        // Increment loop counter
    z_new = zed( p, d, params, An, Cn, def ) ;
        // Update p
    z_i = z_new(i) ;
    p_i = ( z_i - c ) / m ;
        // Update p_i according to z_i = m * p_i + c
    diff = fabs( p(i) - p_i ) ;
        // Update difference

    if( print_level > 0 ){
      Rcout << std::setw(10) <<  it << std::setprecision(4) << std::setw(10) << p(i) <<
        std::setw(10) << z_i << std::setw(10) << p_i << std::setw(10) << diff << std::endl ;
    }

    p(i) = p_i ;
        // Update p ;

    /* Bounds on p in [0,1] **/
    if( p_i < 0 ){
      p(i) = 0 ;
      diff = 0 ;
    }
    if( p_i > 1 ){
      p(i) = 1 ;
      diff = 0 ;
    }
  }
  return p ;
}

// [[Rcpp::export]]
arma::vec p_min_tanget_i( List params, arma::vec p, arma::vec d, double m, int i,
                          arma::vec An, arma::vec Bn, arma::vec Cn, arma::mat def,
                          int maxit=200, double tol=1e-06, int print_level=0 ){
// Finds min_c s.t. z(p) = p + c s.t. z'(p)=1
// If z(p) concave everywhere then this will return 0

  double p_l = 0 ;                        // Lower bound on p
  double p_u = 0 ;                        // Upper bound on p
  mat z_2 = zeros(p.n_elem,2) ;           // Container for the z, z_p output
  mat z_2_inc = zeros(p.n_elem,2) ;       // Container for the z, z_p output
  double c_u = -1 ;                       // Value of c at upper bound on p
  vec p_temp = p ;                        // Temporary container
  vec z_temp = p ;                        // Temporary container
  double inc = 1e-06 ;                    // Increment for computng convexity/concavity
  double z_d2 = 0 ;                       // Second derivative of z

  p_temp(i) = p_l ;
  vec z_l = zed( p_temp, d, params, An, Cn, def ) ;
      // Compute the value of z at p=0
  double c_l = z_l(i) ;                   // Value of c at lower bound on p
  double c_p = .5 * ( c_l + c_u ) ;       // The candidate value of c
  double c_u_old = c_u ;                  // For printing
  double c_l_old = c_l ;                  // For printing

  double diff = 2 * tol ;                 // Difference between iterations
  int it = 0 ;                            // Iteration counter

  if( print_level > 0 ){

    Rcout<< "maxit = " << maxit << " tol = " << tol <<std::endl ;

    Rcout << std::setw(10) << "Iteration" << std::setw(10) << "c_l_old" <<
      std::setw(10) << "c_u_old" << std::setw(10) << "c_p" << std::setw(10) <<
        "p_temp(i)" << std::setw(10) << "z_temp(i)" << std::setw(10) <<
          "z_2(1)" << std::setw(10) << "z_2_inc(1)" <<
            std::setw(10) << "z_d2" << std::setw(10) << "diff" << std::endl ;
    Rcout << std::setw(10) <<  it << std::setprecision(4) << std::setw(10) << c_l_old <<
      std::setw(10) << c_u_old << std::setw(10) << c_p <<
        std::setw(10) << p_temp(i) << std::setw(10) << z_temp(i)<<
          std::setw(10) << z_2(1) << std::setw(10) << z_2_inc(1) <<
            std::setw(10) << z_d2 << std::setw(10) << diff << std::endl ;
  }

  while( diff > tol & it < maxit ){

    it ++ ;
        // Increment loop counter
    c_l_old = c_l ;
    c_u_old = c_u ;
        // Storage
    p_temp = p_ax_c_i( params, p, d, 1, c_p, i, An, Bn, Cn, def, maxit, tol, print_level - 1 ) ;
        // Can use p here because p(i) is the only thing changing
    z_temp = zed( p_temp, d, params, An, Cn, def ) ;
        // Update zed
    if( z_temp(i) > 1 - 1e-08 ){
      c_u = c_p ;
          // If the resulting p is very close to unity then we make the upper
          // limit for c equal to the candidate
    }else{
      z_2 = zed_2(p_temp, d, params, An, Bn, Cn, def) ;
          // Store the value and the derivative
      p_temp(i) = p_temp(i) + inc ;
      z_2_inc = zed_2(p_temp, d, params, An, Bn, Cn, def) ;
          // The gradient at a slightly higher point
      z_d2 = ( z_2_inc(i,1) - z_2(i,1) ) / inc ;
          // The second derivative
      if( z_d2 >= 0 ){
        c_l = c_p ;
            // If concave, replace lower bound
      }else{
        c_u = c_p ;
            // If convex, replace upper bound
      }
    }

    if( print_level > 0 ){
      Rcout << std::setw(10) <<  it << std::setprecision(4) << std::setw(10) << c_l_old <<
        std::setw(10) << c_u_old << std::setw(10) << c_p <<
          std::setw(10) << p_temp(i) << std::setw(10) << z_temp(i)<<
            std::setw(10) << z_2(1) << std::setw(10) << z_2_inc(1) <<
              std::setw(10) << z_d2 << std::setw(10) << diff << std::endl ;
    }

    diff = fabs( c_l - c_u );
        // Update difference
    c_p = .5 * ( c_l + c_u ) ;
        // Update the candidate guess for c

  }


  if( print_level > 0 ){
    Rcout << "\n\n** Final evaluation **" << std::endl ;
  }
  p_temp = p_ax_c_i( params, p, d, 1, c_l, i, An, Bn, Cn, def,
                     maxit, tol, print_level ) ;
      // Evaluate for output at the lower bound

  return p_temp ;
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
    z_guess_err = z_guess - std::max( p_i, 1.1 * tol ) ;
        // Computue the error.  Add a little to p_i to make sure we avoid
        // hitting exactly zero (where ther are multiple solutions possible)
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
        // If no return then call recursively on interval around best point
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
