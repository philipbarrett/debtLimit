/***********************************************************************************
 * surp.cpp
 *
 * Computes the surplus
 *
 * 20feb2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "surp.hpp"
#include "tri.hpp"

// [[Rcpp::export]]
arma::vec v_surp( arma::vec d, arma::vec coeff, arma::vec shift, bool tri ){
  int n = d.n_elem ;
      // Length of the vector of debts
  vec out = zeros(n) ;
  for ( int i=0 ; i < n ; i++ ){
    out[i] = surp( d[i], coeff, shift(i), tri ) ;
  }
  return out ;
}

// [[Rcpp::export]]
double surp( double d, arma::vec coeff, double shift=0, bool tri=false ){
  return tri ? surp_tri( d, coeff, shift ) : surp_poly( d, coeff, shift ) ;
}

// [[Rcpp::export]]
double surp_poly( double d, arma::vec coeff, double shift=0 ){
// Computes the surplus polynomial function from a vector of coeffficients

  double x = d / coeff[0] ;
      // The initial entry is the scaling factor
  double out = coeff[1] + shift ;
      // The constant
  int n = coeff.n_elem ;
      // Length of the coefficient vector
  for ( int i = 2 ; i < n ; i++ ){
    out += coeff[i] * pow( x, i - 1 ) ;
      // Add the debt polynomial
  }
  // out += ( G - 1 ) * 100 * coeff[n-1] ;
  //     // Add the linear function of G
  return out ;
}


// [[Rcpp::export]]
double surp_tri( double d, arma::vec coeff, double shift=0 ){
  // Computes the surplus function which is y1 below x1, and y2 above x2, and
  // follows a triangle distribution CDF between, with intemediate point x3
  double x1 = coeff(0) ;
  double x2 = coeff(1) ;
  double x3 = coeff(2) ;
  double y1 = coeff(3) + shift ;
  double y2 = coeff(4) + shift ;
  // Extract the coefficients
  return y1 + ( y2 - y1 ) * p_triangle( d, x1, x2, x3 ) ; // + coeff(5) * ( G - 1 ) * 100 ;
}

// [[Rcpp::export]]
double d_surp_tri( double d, arma::vec coeff, double shift=0 ){
  double x1 = coeff(0) ;
  double x2 = coeff(1) ;
  double x3 = coeff(2) ;
  double y1 = coeff(3) + shift ;
  double y2 = coeff(4) + shift ;
  // Extract the coefficients
  return ( y2 - y1 ) * d_triangle( d, x1, x2, x3 ) ;
}
