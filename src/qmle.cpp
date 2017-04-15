/***********************************************************************************
 * qmle.cpp
 *
 * Computes quasi-maximum likelihood function for a regime-switching model
 *
 * 14apr2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#include "qmle.hpp"

// [[Rcpp::export]]
double rest_var_lhood( arma::mat Y, arma::vec par,
                       arma::uvec a_switch, arma::umat A_switch, arma::umat Sigma_switch,
                       arma::vec a_vals, arma::mat A_vals, arma::mat Sigma_vals ){
// Computes the (negative) likelihood of the restricted VAR of the form:
//    Y_t = a + A * Y_{t-1} + e_t,  e_t ~ N(0,Sigma)
// Where the integer matrices A_switch and Sigma_switch define the
// parameter restrictions (zeros for restricted to zero, one otherwise)

// TODO: ADD THE FIXED VALUES

  int m = Y.n_cols ;
  int n = Y.n_rows ;
      // Problem dimensions
  mat Yt = Y.cols(1,m-1) ;
  mat Yt1 = Y.cols(0,m-2) ;
      // Select the leading and lagged rows
  vec a = zeros(n) ;
  mat A = zeros(n,n) ;
  mat Sigma = zeros(n,n) ;
      // Initialize the matrives of parameters
  int counter = 0 ;
      // Counter
  for( int i = 0 ; i < n ; i ++ ){
    if( a_switch(i) > 0 ){
      a(i) = par(counter) ;
      counter++ ;
    }else{
      a(i) = a_vals(i) ;
    }
  }
  for( int j = 0 ; j < n ; j ++ ){
    for( int i = 0 ; i < n ; i ++ ){
      if( A_switch(i,j) > 0 ){
        A(i,j) = par(counter) ;
        counter++ ;
      }else{
        A(i,j) = A_vals(i,j) ;
      }
    }
  }
      // Copy parameters into A
  for( int i = 0 ; i < n ; i ++ ){
    for( int j = 0 ; j <= i ; j ++ ){
      if( Sigma_switch(i,j) > 0 ){
        Sigma(i,j) = par(counter) ;
        Sigma(j,i) = par(counter) ;
        counter++ ;
      }else{
        Sigma(i,j) = Sigma_vals(i,j) ;
        Sigma(j,i) = Sigma_vals(j,i) ;
      }
    }
  }
      // Copy parameters into Sigma
  mat e = Yt - A * Yt1 - a * ones<rowvec>(m-1) ;
      // The matrix of residuals
  mat Sigma_I = Sigma.i() ;
      // Inverse of Sigma
  mat e_t = e.t() ;
      // Inverse of Sigma
  double term = 0 ;
      // The final term in the likelihood
  for( int i = 0 ; i < m - 1 ; i ++ ){
    vec temp = e_t.row(i) * Sigma_I * e.col(i) ;
    term += temp(0) ;
  }
  // Rcout << "A = " << A <<std::endl ;
  // // Rcout << "Sigma:\n" << Sigma <<std::endl ;
  // // Rcout << "Sigma_I:\n" << Sigma_I <<std::endl ;
  // Rcout << "term = " << term <<std::endl ;
  // Rcout << "M_PI = " << M_PI <<std::endl ;
  // Rcout << "std::log( 2 * M_PI ) = " << std::log( 2 * M_PI ) <<std::endl ;
  // Rcout << "det(Sigma) = " << det(Sigma) <<std::endl ;
  // Rcout << "std::log(det(Sigma)) = " << std::log(det(Sigma)) <<std::endl ;

  double l = .5 * ( n * std::log( 2 * M_PI ) + std::log( det( Sigma ) ) )
                  + .5 * term / ( m - 1 ) ;
      // The negative log likelihood
  return l ;
}




