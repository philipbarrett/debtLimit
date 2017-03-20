/***********************************************************************************
* q.cpp
*
* Defines the pricing function
*
* 13mar2017
* Philip Barrett, DC
*
***********************************************************************************/

#include "q.hpp"

// [[Rcpp::export]]
arma::vec q_fn( arma::vec R, arma::vec p, arma::mat trans, double lambda, double phi,
                int n,  std::string cont_type, arma::vec G, arma::vec An,
                arma::mat def ){
// The bond price function, computed from the equation:
//    Rq = (1-lambda)*phi*p*q + (1-p)*(1-lambda) + lambda*qe
//
// cont_type defines the expected continuation price qe.  Can be:
//  - "avg" the continuation price is the conditional average of the other
//          default prices, so qe=(1-p)*sum_x' pi(x'|x)q(x')
//  - "fix" the continuation price is fixed, given by the vector An
//  - "low" the continuation price is the average over the states with lower R-G,
//          on the grounds that default is likely to occur in the "worse" states
//          (requires that the vector G be submitted). Where qe=sum_x' p(x'|x)q(x')
//          and p(x'|x) = pi(x'|x) for R(x')-G(x') < R(x)-G(x) and 0 o/w
//  - "def" same as low but with a matrix of ones and zeroes defining conditional
//          repayment states.

// TODO: Think about a version wherestate-conditional prices are submitted and
// probability of receiving the bond at that price is a function of p.

  vec q = 0 * R ;
  // Initialize the price vector

  if( cont_type == "avg" ){
    // Straight average over default thresholds

    mat A = diagmat(R) - ( 1 - lambda ) * phi * diagmat(p) - lambda * trans ;
        // Valuation matrix
    vec B = ( 1 - lambda ) * ( 1 - p ) ;
        // Period payoff
        // Rcout << "A:\n" << A << std::endl ;
        // Rcout << "B:\n" << B << std::endl ;
    q = solve( A, B ) ;

  }else if( cont_type == "fix" ){
    // The continuation price is known and given by An
    q = ( ( 1 - lambda ) * ( 1 - p ) + lambda * An ) / ( R - ( 1 - lambda ) * phi * p ) ;

  }else{
    // Average over selected default thresholds based on R-G

    if( cont_type == "low" ){
      // Need to create the list of default states
      def = zeros(n,n) ;
      // Reinitiate the matrix of default states (zero for default, 1 for repayment)
      for( int i = 0 ; i < n ; i++ ){
        for( int j = 0 ; j < n ; j++ ){
          def(i,j) = ( R(i) - G(i) >= R(j) - G(j) ) ? 1.0 : 0.0 ;
          // Because if R-G rises then the state is better
        }
      }
    }

    mat trans2 = def % trans ;
        // Initiate the alternative tranisiton matrix
        // Rcout << "def:\n" << def << std::endl ;
        // Rcout << "trans:\n" << trans << std::endl ;
        // Rcout << "trans2:\n" << trans2 << std::endl ;
    mat A = diagmat(R) - ( 1 - lambda ) * phi * diagmat(p) - lambda * trans2 ;
    vec B = ( 1 - lambda ) * ( 1 - p ) ;
        // Rcout << "A:\n" << A << std::endl ;
        // Rcout << "B:\n" << B << std::endl ;
    q = solve( A, B ) ;
        // Solve for the price
  }
  return q ;
}



// [[Rcpp::export]]
arma::mat q_d_p( arma::vec R, arma::vec p, arma::mat trans, double lambda, double phi,
                int n, std::string cont_type, std::string d_type,
                arma::vec G, arma::vec An, arma::vec Bn, arma::mat def ){
// Computes the derivative of q w.r.t. p
//    Rq' = (1-lambda)*phi*(q + p*q') - (1-lambda) + lambda * d/dp( qe )
//  Where d/dp(qe) depends on cont_type:
//      - "avg" trans * result
//      - "fix" Bn

  mat dq = zeros(n,n) ;
      // Initialize the output matrix

  if( d_type == "num" ){
  // Numerical differentiation
    for( int i = 0 ; i < n ; i++ ){
      dq.col(i) = q_d_p_num_i(R, p, trans, lambda, phi, n, i, cont_type, G, An, def ) ;
    }
    return dq ;
  }

  vec q = q_fn(R, p, trans, lambda, phi, n, cont_type, G, An, def) ;
      // The current price: Required in a bunch of the derivatives

  if( cont_type == "avg" ){
  // Straight average over default thresholds
    mat A = diagmat(R) - ( 1 - lambda ) * phi * diagmat(p) - lambda * trans ;
        // Valuation matrix (kinda)
    mat B = ( 1 - lambda ) * diagmat( - 1 + phi * q ) ;
        // Period payoff (sorta).  Diagonal version from the fact that q is
        // directly affected only by p in the *current* state.
    dq = solve( A, B ) ;

  }else if( cont_type == "fix" ){
  // The continuation derivative is known and given by Bn
      // Rcout << "q:\n" << q << std::endl ;
      // Rcout << "p:\n" << p << std::endl ;
      // Rcout << "Bn:\n" << Bn << std::endl ;
      // Rcout << "R:\n" << R << std::endl ;
    dq = diagmat( ( ( 1 - lambda ) * ( - 1 + phi * q ) +
                      ( 1 - p ) * lambda % Bn ) / ( R - ( 1 - lambda ) * phi * p ) ) ;
        // The derivative
  }
  return dq ;
}


// [[Rcpp::export]]
arma::vec q_d_p_num_i( arma::vec R, arma::vec p, arma::mat trans, double lambda, double phi,
                     int n, int i, std::string cont_type, arma::vec G, arma::vec An,
                     arma::mat def ){
// Computes the vector of numerical derivatives in the i-direction

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
  vec q_d = ( q_fn(R, p_h, trans, lambda, phi, n, cont_type, G, An, def ) -
    q_fn(R, p_l, trans, lambda, phi, n, cont_type, G, An, def ) ) / dist ;
      // The derivative of z (computing in all dimensions)
  return q_d ;
}
