/***********************************************************************************
 * sim.cpp
 *
 * Computes the model simulation
 *
 * 14apr2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "sim.hpp"
#include "surp.hpp"
#include "tri.hpp"
#include "qiter.hpp"

// double zed_sim( double p, int i_eta, double d, arma::vec d_bar, List params ) {
// // Computes zed: Local copy in sim.cpp
//
//   /** 1. Extract parameters **/
//   mat trans = params["trans"] ;
//   vec v_eta = params["v.eta"] ;
//   double phi = params["phi"] ;
//   vec v_s_coeff = params["v.s.coeff"] ;
//   double eps_bar = params["eps.bar"] ;
//
//   /** 2. Create threhold values for the surplus shock **/
//   int m = v_s_coeff.n_elem - 1 ;
//       // Order of surplus rule
//   int n = d_bar.n_elem ;
//       // Number of dimensions
//   vec H = ( 1 + v_eta + ( 1 - phi ) * p / ( 1 - p ) ) * d - d_bar - surp_triangle( d, v_s_coeff ) ;
//       // The surplus shock that sets end-of-next-period debt to d.bar
//   vec G = zeros(n) ;
//   for( int i = 0 ; i < n ; i++ ){
//     G[i] = p_triangle( H[i], - eps_bar, 0, eps_bar ) ;
//   }
//       // The state-dependent vector of default probabilities and its derivative
//   double z = as_scalar( trans.row(i_eta-1) * G ) ;
//       // The outputs
//   return z ;
// }
//
// [[Rcpp::export]]
NumericMatrix sim_core( arma::vec i_idx_R, arma::vec d_bar, arma::vec d_grid, arma::mat P,
                        arma::mat Q, List params, double d0, bool s_flag, arma::vec s_in ){
// Computes the simulation of the model
// Outputs: i_idx, R, G, r-g, expected (r-g)', d, d.bar, default indicator, s(d), eps,
// s, p, q, d', p', q'

  /** 0. Extract parameters **/
  const int n = i_idx_R.n_elem ;
  double phi = params["phi"] ;
  double lambda = params["lambda"] ;
  double s_sd = params["surp.sd"] ;
  double eps_bar = std::sqrt(6) * s_sd ;
  vec v_R = params["R"] ;
  vec v_G = params["G"] ;
  vec v_s_coeff = params["v.s.coeff"] ;
  vec shift = params["s.shift"] ;
  vec i_idx_c = i_idx_R - 1 ;                     // c-style counter
  mat trans = params["trans"] ;
  int n_states =v_R.n_elem ;                      // umber of states

  /** 1..  The independent parts **/
  vec eps = zeros( n ) ;
  NumericVector q_eps = runif(n) ;                // Initiate the vector for eps draws
  for ( int i = 0 ; i < n ; i++ ){
    eps(i) = q_triangle( q_eps(i), - eps_bar, 0, eps_bar ) ;
  }   // The vector of eps

  /** 2. The output **/
  mat out = zeros( n, 16 ) ;
  out.col(0) = i_idx_R ;
  out.col(9) = eps ;
      // The exogenous stuff
  vec e_R_prime = trans * v_R ;
  vec e_G_prime = trans * v_G ;
      // The expected value of eta in the next period.
  vec v_d = zeros(1) ;                        // Placeholder for interpolation
  vec p_prime = zeros(1) ;                    // Placeholder for interpolation
  vec q_prime = zeros(1) ;                    // Placeholder for interpolation
  vec e_grid = zeros(1) ;                     // Need to use a grid of continuation e
  mat m_d_prime = zeros(1,n_states) ;         // Because d_prime returns a matrix

  v_d(0) = d0 ;
  vec this_P = conv_to<vec>::from(P.row(v_R(i_idx_c(0)))) ;
  vec this_Q = conv_to<vec>::from(Q.row(v_R(i_idx_c(0)))) ;
      // Generate initial interpolation points
  interp1( d_grid, this_P, v_d, p_prime ) ;
  interp1( d_grid, this_Q, v_d, q_prime ) ;
  out(0,13) = d0 ;
  out(0,14) = p_prime(0) ;
  out(0,15) = q_prime(0) ;
      // Period 0.  Debt, debt price and default probability.

  for( int i = 1 ; i < n ; i ++ ){
    out(i,1) = v_R(i_idx_c(i)) ;              // The realized value of R
    out(i,2) = v_G(i_idx_c(i)) ;              // The realized value of G
    out(i,3) = out(i,1) - out(i,2) ;          // The realized value of R-G
    out(i,4) = e_R_prime(i_idx_c(i)) - e_G_prime(i_idx_c(i)) ;
        // The expected level of R-G in the next period
    out(i,5) = out(i-1,13) ;
        // Debt
    out(i,6) = d_bar(i_idx_c(i)) ;
        // The debt limiit
    out(i,8) = surp_tri( out(i,5), v_s_coeff, shift(i_idx_c(i)) ) ;
    if(s_flag){
      out(i,10) = s_in(i) ;
      out(i,9) = out(i,10) - out(i,8) ;
    }else{
      out(i,10) = out(i,8) + out(i,9) ;
          // The surplus rule and realized surplus
    }
    out(i,11) = out(i-1,14) ;
    out(i,12) = out(i-1,15) ;
        // The inherited probability of default and debt price
    e_grid(0) = out(i,9) ;
    m_d_prime = d_prime( i_idx_c(i), out(i,5), d_bar, out(i,12), Q, d_grid, v_G, shift,
                         lambda, e_grid, v_s_coeff, true, zeros(1,1), false, 0, 1e-05, 20, false ) ;
    out(i,13) = std::max( m_d_prime(0,i_idx_c(i)), 0.0 ) ;
        // Continuation debt
    out(i,7) = ( out(i,5) <= out(i,6) ) ? 0 : 1 ;
        // Default indicator
    if( out(i,7) == 1 ){
    // If default happens, then debt is written down
      out(i,13) = phi * ( 1 - lambda ) * out(i,12) * out(i,5) ;
      out(i,14) = 1 ;
      out(i,15) = 1 ;
          // Claims worth a fraction of the repayment due are written (can
          // assume certain repayment here because only the continuation debt
          // burden in the next period really matters)
    }else{
      v_d(0) = out(i,5) ;
      vec this_P = conv_to<vec>::from(P.row(v_R(i_idx_c(0)))) ;
      vec this_Q = conv_to<vec>::from(Q.row(v_R(i_idx_c(0)))) ;
          // Generate initial interpolation points
      interp1( d_grid, this_P, v_d, p_prime ) ;
      interp1( d_grid, this_Q, v_d, q_prime ) ;
      out( i, 14 ) = p_prime(0) ;
      out( i, 15 ) = q_prime(0) ;
          // The default probability and continuation debt price
    }
  }

  NumericMatrix nm_out = wrap(out) ;
  colnames(nm_out) = CharacterVector::create("idx", "R", "G", "rmg", "e.rmg", "d", "d.bar",
           "def", "s.d", "eps", "s", "p", "q", "d.prime", "p.prime", "q.prime" ) ;

  return nm_out ;

}

//
// // [[Rcpp::export]]
// NumericMatrix sim( List sol, List dc, const int n = 1e06, double d0 = 0 ){
//   // Wrapper for sim_core
//   List params = sol["params"] ;
//   vec d_bar = sol["d.bar"] ;
//   mat dc_d = dc["d"] ;
//   mat dc_p = dc["p"] ;
//   mat dc_rp = dc["prem"] ;
//   NumericMatrix out = wrap(sim_core( n, d_bar, dc_d, dc_p, dc_rp, params, d0 )) ;
//   colnames(out) = CharacterVector::create("i.eta", "eta", "e.eta", "d", "d.bar", "def",
//                   "s.d", "eps", "s", "p", "rp", "eri", "d.prime", "p.prime", "rp.prime" ) ;
//   return out ;
// }

// #include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace arma;
//

