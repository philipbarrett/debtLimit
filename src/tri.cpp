/***********************************************************************************
 * tri.cpp
 *
 * Defines the functions governing the triangle distribution
 *
 * 07mar2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include <math.h>
#include "tri.hpp"

// [[Rcpp::export]]
double d_triangle( double x, double a, double b, double c ){
  // Triangle distribution density on (a,b,c)
  if( x <= a || x >= c )
    return 0 ;
  double z = 2 / ( c - a ) ;
  // The peak of the triangle
  if( x <= b )
    return ( x - a ) / ( b - a ) * z ;
  if( x > b )
    return ( c - x ) / ( c - b ) * z ;
}

// [[Rcpp::export]]
double p_triangle( double x, double a, double b, double c ){
  // Triangle distribution CDF on (a,b,c)
  if( x <= a )
    return 0 ;
  if( x >= c )
    return 1 ;
  double z = 2 / ( c - a ) ;
  // The peak of the triangle
  if( x <= b )
    return 0.5 * z / ( b - a ) * ( x * ( x - 2 * a ) + pow( a, 2 ) ) ;
  if( x > b )
    return 1 - p_triangle( c - x, 0, c - b, c - a ) ;
}

// [[Rcpp::export]]
double d_triangle2( double x, double a, double b, double c ){
  if( x <= a )
    return 0 ;
  if( x >= c )
    return 0 ;
  double z = 2 / ( c - a ) ;
  // The peak of the triangle
  if( x < 0 )
    return 1 / ( b - a ) * z ;
  if( x >= 0 )
    return - 1 / ( c - b ) * z ;
}

// [[Rcpp::export]]
double q_triangle( double q, double a, double b, double c ){
// Triangle distribution inversion
  if( q <= ( b - a ) / ( c - a ) ){
    double disc = q * ( b - a ) * ( c - a )  ;
      // The discriminant
      // Because x^2 - 2ay + C = 0, where C = a^2 - disc
    return a + pow( disc, .5 ) ;
  }
  return( c - q_triangle( 1 - q, 0, c - b, c - a ) ) ;
      // If q is below the mode then call recursively
}
