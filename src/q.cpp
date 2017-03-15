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
                int n,  std::string cont_type, arma::vec G, arma::vec q_e,
                arma::mat def ){
// The bond price function, computed from the equation:
//    Rq = (1-lambda)*phi*q + (1-p)*(1-lambda) + lambda*q_e
//
// cont_type defines the expected continuation price q_e.  Can be:
//  - "avg" the continuation price is the conditional average of the other
//          default prices, so q_e=(1-p)*sum_x' pi(x'|x)q(x')
//  - "fix" the continuation price is fixed, given by the vector q_e
//  - "low" the continuation price is the average over the states with lower R-G,
//          on the grounds that default is likely to occur in the "worse" states
//          (requires that the vector G be submitted). Where q_e=sum_x' p(x'|x)q(x')
//          and p(x'|x) = pi(x'|x) for R(x')-G(x') < R(x)-G(x) and 0 o/w
//  - "def" same as low but with a matrix of ones and zeroes defining conditional
//          repayment states.

// TODO: Think about a version wherestate-conditional prices are submitted and
// probability of receiving the bond at that price is a function of p.

  vec q = 0 * R ;
  // Initialize the price vector

  if( cont_type == "avg" ){
    // Straight average over default thresholds

    mat A = diagmat(R) - ( 1 - lambda ) * phi * diagmat(p) -
                  lambda * ( eye(n,n) - diagmat(p) ) * trans ;
        // Valuation matrix
    vec B = ( 1 - lambda ) * ( 1 - p ) ;
        // Period payoff
        // Rcout << "A:\n" << A << std::endl ;
        // Rcout << "B:\n" << B << std::endl ;
    q = solve( A, B ) ;

  }else if( cont_type == "fix" ){
    // The continuation price is known and given by q_e
    q = ( ( 1 - lambda ) * ( 1 - p ) + lambda * q_e ) / ( R - ( 1 - lambda ) * phi * p ) ;

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
