/***********************************************************************************
 * zed.cpp
 *
 * Defines the function z(.)
 *
 * 20feb2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "zed.hpp"
#include "surp.hpp"
#include "q.hpp"
#include "tri.hpp"

/** NEED TO CHANGE ALL THE qe AND q to An, Bn, ALSO THINK ABOUT Cn **/

// [[Rcpp::export]]
arma::vec zed( arma::vec p, arma::vec d, List params, arma::vec An,
               arma::vec Cn, arma::mat def ) {
// Computes zed

  /** 1. Extract parameters **/
  mat trans = params["trans"] ;
  vec R = params["R"] ;
  vec G = params["G"] ;
  vec v_s_coeff = params["v.s.coeff"] ;
  double phi = params["phi"] ;
  double surp_sd = params["surp.sd"] ;
  double lambda = params["lambda"] ;
  int n = R.n_elem ;
  bool tri = params["tri"] ;
  bool d_tri = params["d.tri"] ;
  std::string cont_type = params["cont.type"] ;

  /** 2. Create the vector of bond prices consistent with default probabilities **/
  vec q = q_fn(R, p, trans, lambda, phi, n, cont_type, G, An, def ) ;
      // Candidate debt price at the boundary
  vec qprime = ( cont_type == "fix" ) ? Cn : q ;
      // The price of the continuation debt.

  /** 3. Create threshold values for the surplus shock **/
  vec surp = v_surp( d, v_s_coeff, G, tri ) ;
      // The surplus vector.  Depends on growth and debt in the *preceding* period.
  mat H = mat(n,n) ;
  mat p_H = mat(n,n) ;
      // The matrix H of surplus shocks required to hit the debt limit and their
      // corresponding probabilities

  for( int i = 0 ; i < n ; i++ ){
    for( int j = 0 ; j < n ; j++ ){
      H(i,j) = d(i) * ( ( 1 - lambda ) + lambda * qprime(j) ) / ( q(i) * G(j) ) -
                    d(j) - surp(i) ;
          // Market value of new debt in next period is # of old obligations *
          // current market price.
      if( d_tri ){
        double ub = surp_sd * std::sqrt( 6.0 ) ;
            // Var = 18 * (Upper bound) ^ 2 (symmetric triangle case)
        p_H(i,j) = p_triangle( H(i,j), -ub, 0, ub ) ;
      }else{
        p_H(i,j) = R::pnorm( H(i,j), 0, surp_sd, 1, 0 ) ;
      }
          // Fill these in
    }
  }

  /** 4. Create the output vectors **/
  vec z = zeros(n) ;
  for( int i=0 ; i < n ; i++ ){
    z(i) = dot( trans.row(i), p_H.row(i) ) ;
  }
      // Create the output as the expectation of the default probabilities over the states.
  return z ;
}

// [[Rcpp::export]]
arma::mat zed_2( arma::vec p, arma::vec d, List params,
                 arma::vec An, arma::vec Bn, arma::vec Cn, arma::mat def ) {
// Computes zed and zed_p simultaneously

  /** 1. Select the method for differentiation **/
  std::string method = params["diff.method"] ;
      // Can be: "num", "ana"
  if( method == "num" ){
    return zed_2_num( p, d, params, An, Cn, def ) ;
  }else if( method == "ana" ){
    return zed_2_ana( p, d, params, An, Bn, Cn, def ) ;
  }else{
    Rcout << "No recognized differentiation method submitted" << std::endl ;
  }
}

// [[Rcpp::export]]
arma::mat zed_2_num( arma::vec p, arma::vec d, List params,
                     arma::vec Cn, arma::vec An, arma::mat def){
// Uses numerical differentiation to compute z_p
  int n = p.n_elem ;
      // Number of states
  mat out = zeros( n, 2 ) ;
      // The output matrix
  out.col(0) = zed( p, d, params, An, Cn, def ) ;
      // The level
  for( int i = 0 ; i < n ; i++ ){
    out(i,1) = zed_2_num_d_i( p, d, params, i, An, Cn, def ) ;
  }   // Compute the derivatives
  return out ;
}

// [[Rcpp::export]]
double zed_2_num_d_i( arma::vec p, arma::vec d, List params, int i,
                      arma::vec An, arma::vec Cn, arma::mat def ){
// Uses numerical differentiation to compute z_p in the ith dimension
  vec p_l = p ;
  vec p_h = p ;
      // Initiate the upper and lower boundaries of the p vector
  double inc = 1e-06 ;
      // The increment
  p_l[i] = ( p_l[i] > inc ) ? p_l[i] - inc : p[i] ;
  p_h[i] = ( p_h[i] < 1 - inc ) ? p_h[i] + inc : p[i] ;
      // The evaluation vectors
  double dist = p_h[i] - p_l[i] ;
      // The distance
  vec zed_d = ( zed( p_h, d, params, An, Cn, def ) -
                    zed( p_l, d, params, An, Cn, def ) ) / dist ;
      // The derivative of z (computing in all dimensions)
  return zed_d(i) ;
}

// [[Rcpp::export]]
arma::mat zed_2_ana( arma::vec p, arma::vec d, List params, arma::vec An,
                     arma::vec Bn, arma::vec Cn, arma::mat def){
// Computes (z,z.p) analytically

  /** 1. Extract parameters **/
  mat trans = params["trans"] ;
  vec R = params["R"] ;
  vec G = params["G"] ;
  vec v_s_coeff = params["v.s.coeff"] ;
  double phi = params["phi"] ;
  double surp_sd = params["surp.sd"] ;
  double lambda = params["lambda"] ;
  int n = R.n_elem ;
  bool tri = params["tri"] ;
  bool d_tri = params["d.tri"] ;
  std::string cont_type = params["cont.type"] ;

  /** 2. Create the vector of bond prices consistent with default probabilities **/
  vec q = q_fn(R, p, trans, lambda, phi, n, cont_type, G, An, def ) ;
      // Current period debt price
  mat q_dp = q_d_p( R, p, trans, lambda, phi, n, cont_type, "ana", G, An, Bn, def ) ;
      // The derivative of q w.r.t. p
  vec qprime = ( cont_type == "fix" ) ? Cn : q ;
      // Next period debt price

  /** 3. Create threshold values for the surplus shock **/
  vec surp = v_surp( d, v_s_coeff, G, tri ) ;
      // The surplus vector.  Depends on growth and debt in the *preceding* period.
  mat H = zeros(n,n) ;
  mat p_H = zeros(n,n) ;
  mat p_H_dp = zeros(n,n) ;
  mat H_dp = zeros(n,n) ;
      // The matrix H of surplus shocks required to hit the debt limit, their
      // corresponding probabilities, and the derivatives w.r.t. p

  for( int i = 0 ; i < n ; i++ ){
    for( int j = 0 ; j < n ; j++ ){
      H(i,j) = d(i) * ( ( 1 - lambda ) + lambda * qprime(j) ) / ( q(i) * G(j) ) -
        d(j) - surp(i) ;
          // Market value of new debt in next period is # of old obligations *
          // current market price.
      if( cont_type == "fix" ){
        H_dp(i,j) = - d(i) * ( ( 1 - lambda ) + lambda * qprime(j) ) /
                        ( q(i) * G(j) ) * q_dp(i,i) / q(i) ;
      }else{
        H_dp(i,j) = - d(i) * ( ( 1 - lambda ) + lambda * qprime(j) ) /
                        ( q(i) * G(j) ) * q_dp(i,i) / q(i)
                    + d(i) * lambda * q_dp(j,i) / ( q(i) * G(j) ) ;
      }
          // The derivative of H w.r.t. p.  The first conditional version is
          // "correct", whereas the second should only be used as an
          // approximation.
      p_H(i,j) = R::pnorm( H(i,j), 0, surp_sd, 1, 0 ) ;
      p_H_dp(i,j) = R::dnorm( H(i,j), 0, surp_sd, 0 ) ;

      if( d_tri ){
        double ub = surp_sd * std::sqrt( 6.0 ) ;
            // Var = 18 * (Upper bound) ^ 2 (symmetric triangle case)
        p_H(i,j) = p_triangle( H(i,j), -ub, 0, ub ) ;
        p_H_dp(i,j) = d_triangle( H(i,j), -ub, 0, ub ) ;
      }else{
        p_H(i,j) = R::pnorm( H(i,j), 0, surp_sd, 1, 0 ) ;
        p_H_dp(i,j) = R::dnorm( H(i,j), 0, surp_sd, 0 ) ;
      }
          // Fill these in too
    }
  }
  mat H_temp = H_dp % p_H_dp ;

    // mat q_dp2 = q_d_p( R, p, trans, lambda, phi, n, cont_type, "num", G, An, Bn, def ) ;
    // Rcout << "trans:\n" << trans << std::endl ;
    // Rcout << "H:\n" << H << std::endl ;
    // Rcout << "H_temp:\n" << H_temp << std::endl ;
    // Rcout << "p_H:\n" << p_H << std::endl ;
    // Rcout << "H_dp:\n" << H_dp << std::endl ;
    // Rcout << "p_H_dp:\n" << p_H_dp << std::endl ;
    // Rcout << "q_dp:\n" << q_dp << std::endl ;
    // Rcout << "q_dp2:\n" << q_dp2 << std::endl ;
    // MOAR TESTS THEN DELETE

  /** 4. Create the output vectors **/
  vec z = zeros(n) ;
  vec z_p = zeros(n) ;
  for( int i=0 ; i < n ; i++ ){
    z(i) = dot( trans.row(i), p_H.row(i) ) ;
    z_p(i) = dot( trans.row(i), H_temp.row(i) ) ;
  }
      // Create the output as the expectation of the default probabilities over the states.

  mat out = zeros(n,2) ;
      // Initialize the output
  out.col(0) = z ;
  out.col(1) = z_p ;
      // Populate the output
  return out ;
}


// // [[Rcpp::export]]
// arma::mat zed_2_jac_i( double p, int i_eta, double d, arma::vec d_bar, List params ){
//
//   /** 1. Extract parameters **/
//   mat trans = params["trans"] ;
//   vec v_eta = params["v.eta"] ;
//   double phi = params["phi"] ;
//   vec v_s_coeff = params["v.s.coeff"] ;
//   double eps_bar = params["eps.bar"] ;
//
//   /** 2. . Create the first derivatives of threshold values for the surplus shock **/
//   int m = v_s_coeff.n_elem - 1 ;
//       // Order of surplus rule
//   int n = d_bar.n_elem ;
//       // Number of dimensions
//   double surp = surp_triangle( d, v_s_coeff ) ;
//       // The surplus
//   vec H = ( 1 + v_eta + ( 1 - phi ) * p / ( 1 - p ) ) * d - d_bar - surp ;
//       // The surplus shock that sets end-of-next-period debt to d.bar
//   double H_p = ( 1 - phi ) * ( 1 / ( 1 - p ) + p / pow( 1 - p , 2 ) ) * d ;
//       // Derivative of H wrt p (same in all dimensions - so store as double)
//   vec v_i = zeros(n) ;
//   v_i[i_eta-1] = 1.0 ;
//       // The basis vector in direction i
//
//   double surp_d = d_surp_triangle( d, v_s_coeff ) ;
//       // The surplus derivative w.r.t. d
//   vec H_d = ( 1 + v_eta + ( 1 - phi ) * p / ( 1 - p ) ) - ( v_i + surp_d ) ;
//       // Derivative wrt d is a scalar plus an extra part in the ith element
//
//   /** 3. Create the probabilities and first derivatives **/
//   vec G_eps = zeros(n) ;
//   vec G_eps_2 = zeros(n) ;
//   for( int i = 0 ; i < n ; i++ ){
//     G_eps[i] = d_triangle( H[i], - eps_bar, 0, eps_bar ) ;
//     G_eps_2[i] = d_triangle2( H[i], - eps_bar, 0, eps_bar ) ;
//   }
//       // The first and second derivatives of G wrt eps
//   double z_p = as_scalar( trans.row(i_eta-1) * G_eps ) * H_p ;
//       // The derivative of z wrt p
//   double z_d = as_scalar( trans.row(i_eta-1) * ( G_eps % H_d ) ) ;
//       // The derivative of z wrt d
//
//   /** 4. Create the second derivatives of threshold values for the surplus shock **/
//   double H_p_p = ( 1 - phi ) * ( 2 / pow( 1 - p, 2.0 ) + 2 * p / pow( 1 - p, 3.0 ) ) * d ;
//       // The second derivative of H wrt p
//   double H_p_d = H_p / d ;
//       // The second derivative of H wrt p & d
//
//   /** 5. Create the probabilities and second derivatives**/
//   double z_p_p = as_scalar( trans.row(i_eta-1) * ( G_eps_2 * H_p * H_p + G_eps * H_p_p ) ) ;
//       // Derivative of z_p wrt p
//   double z_p_d = as_scalar( trans.row(i_eta-1) * ( ( G_eps_2 % H_d ) * H_p + G_eps * H_p_d ) ) ;
//       // Derivative of z_p wrt d
//
//   /** 6. Construct the output **/
//   mat out = zeros(2,2) ;
//   out(0,0) = z_p ;
//   out(0,1) = z_d ;
//   out(1,0) = z_p_p ;
//   out(1,1) = z_p_d ;
//       // The output
//   return out ;
//       // The jacobian
// }
