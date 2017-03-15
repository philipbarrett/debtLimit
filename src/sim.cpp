// #include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;
// using namespace arma;
//
// // [[Rcpp::export]]
// double q_triangle( double q, double a, double b, double c ){
// // Triangle distribution inversion
//   if( q <= ( b - a ) / ( c - a ) ){
//     double disc = q * ( b - a ) * ( c - a )  ;
//         // The discriminant
//         // Because x^2 - 2ay + C = 0, where C = a^2 - disc
//     return a + pow( disc, .5 ) ;
//   }
//   return( c - q_triangle( 1 - q, 0, c - b, c - a ) ) ;
//       // If q is below the mode then call recursively
// }
//
// // [[Rcpp::export]]
// arma::vec markov_sim( const int n, const NumericMatrix M, const int s0,
//                           const int n_s ){
// // Fast Markov simulation
//
//   NumericVector out(n) ;
//   out(0) = s0 ;
//   // Initialize output
//   NumericMatrix sum_M( n_s, n_s ) ;
//   // The row sum of M
//   for( int i=0 ; i < n_s ; i++ ){
//     sum_M( i, 0 ) = M( i , 0 ) ;
//     // The first column
//     for( int j=1 ; j < n_s ; j++ ){
//       sum_M( i, j ) = sum_M( i, j - 1 ) + M( i , j ) ;
//     }
//   } // Create the row sum of M
//
//   NumericVector shks = runif( n ) ;
//   // The vector of random shocks on [0,1]
//   int this_s = 0 ;
//
//   for( int i = 1 ; i < n ; i++ ){
//     this_s = 0 ;
//     // Initialize counters
//     while( shks(i) > sum_M( out(i-1), this_s ) ){
//       this_s++ ;
//     }
//     out( i ) = this_s ;
//   }
//
//   vec out_arma = zeros(n) ;
//   for ( int i = 0 ; i < n ; i++ )
//     out_arma(i) = out(i) ;
//         // Convert to arma type
//
//   return out_arma ;
// }
//
// double p_triangle( double x, double a, double b, double c ){
// // Triangle distribution CDF on (a,b,c)
//   if( x <= a )
//     return 0 ;
//   if( x >= c )
//     return 1 ;
//   double z = 2 / ( c - a ) ;
//   // The peak of the triangle
//   if( x <= b )
//     return 0.5 * z / ( b - a ) * ( x * ( x - 2 * a ) + pow( a, 2 ) ) ;
//   if( x > b )
//     return 1 - p_triangle( c - x, 0, c - b, c - a ) ;
// }
//
// double surp_triangle( double d, vec coeff ){
//   // Computes the surplus function which is y1 below x1, and y2 above x2, and
//   // follows a triangle distribution CDF between, with intemediate point x3
//   double x1 = coeff(0) ;
//   double x2 = coeff(1) ;
//   double x3 = coeff(2) ;
//   double y1 = coeff(3) ;
//   double y2 = coeff(4) ;
//   // Extract the coefficients
//   return y1 + ( y2 - y1 ) * p_triangle( d, x1, x2, x3 ) ;
// }
//
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
// // [[Rcpp::export]]
// arma::mat sim_core( const int n, vec d_bar, arma::mat dc_d, arma::mat dc_p,
//                     arma::mat dc_rp, List params, double d0 ){
// // Computes the simulation of the model
// // Outputs: i_eta, r-g, expected (r-g)', d, d.bar, default indicator, s(d), eps,
// // s, p, RP, effective int rate, d', p', RP'
//
//   /** 0. Extract parameters **/
//   mat trans = params["trans"] ;
//   NumericMatrix trans_nm = params["trans"] ;
//   vec v_eta = params["v.eta"] ;
//   int n_eta = v_eta.n_elem ;
//   double phi = params["phi"] ;
//   vec v_s_coeff = params["v.s.coeff"] ;
//   double eps_bar = params["eps.bar"] ;
//
//   /** 1..  The independent parts **/
//   int s0 = n_eta / 2 ;
//     // Approimate centre of distribution
//   vec i_eta_c = markov_sim( n, trans_nm, s0, n_eta ) ;
//   vec i_eta_R = i_eta_c + 1 ;
//     // The indices of i_eta
//   vec eps = zeros( n ) ;
//   NumericVector q_eps = runif(n) ;
//       // Initiate the vector for eps draws
//   for ( int i = 0 ; i < n ; i++ ){
//     eps(i) = q_triangle( q_eps(i), - eps_bar, 0, eps_bar ) ;
//   }   // The vector of eps
//
//   /** 2. The output **/
//   mat out = zeros( n, 15 ) ;
//   out.col(0) = i_eta_R ;
//   out.col(7) = eps ;
//       // The exogenous stuff
//   out(0,11) = d0 ;
//       // Period 0.  First period or so is inconsitent.  Not a big deal.
//   vec e_eta_prime = trans * v_eta ;
//       // The expected value of eta in the next period.
//   vec p_prime = zeros(1) ;
//   vec rp_prime = zeros(1) ;
//   vec v_d = zeros(1) ;
//   vec this_dc_d = dc_d.col(0) ;
//   vec this_dc_p = dc_p.col(0) ;
//   vec this_dc_rp = dc_rp.col(0) ;
//
//   for( int i = 1 ; i < n ; i ++ ){
//     out(i,1) = v_eta(i_eta_c(i)) ;
//         // The realized value of r-g
//     out(i,2) = e_eta_prime(i_eta_c(i)) ;
//         // The expected level of eta in the next period
//     out(i,3) = out(i-1,12) ;
//         // Debt
//     out(i,4) = d_bar(i_eta_c(i)) ;
//         // The debt limiit
//     out(i,5) = ( out(i,3) <= out(i,4) ) ? 0 : 1 ;
//         // Default indicator
//     out(i,6) = surp_triangle( out(i,3), v_s_coeff ) ;
//     out(i,8) = out(i,6) + out(i,7) ;
//         // The surplus rule and realized surplus
//     out(i,9) = out(i,13) ;
//     out(i,10) = out(i,14) ;
//         // The inherited probability of default and associated risk premium
//     out(i,11) = 1 + out(i,1) + out(i,10) ;
//         // The effective growth rate of debt
//     out(i,12) = std::max( out(i,3) * out(i,11) - out(i,8), 0.0 ) ;
//         // Continuation debt
//     v_d(0) = out(i,3) ;
//     this_dc_d = dc_d.col(i_eta_c(i)) ;
//     this_dc_p = dc_p.col(i_eta_c(i)) ;
//     this_dc_rp = dc_rp.col(i_eta_c(i)) ;
//         // Set up linear interpolation
//     if( out(i,3) < this_dc_d(0) ){
//       p_prime(0) = 0 ;
//       rp_prime(0) = 0 ;
//           // If no risk of default in next period
//     }else{
//       interp1( this_dc_d, this_dc_p, v_d, p_prime) ;
//       interp1( this_dc_d, this_dc_rp, v_d, rp_prime ) ;
//           // Interpolation for continuation default probability and risk premium
//     }
//     out( i, 13 ) = p_prime(0) ;
//     out( i, 14 ) = rp_prime(0) ;
//         // The default probability and risk premium
//   }
//
//   return out ;
//
// }
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
