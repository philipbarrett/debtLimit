/***********************************************************************************
 * surp.hpp
 *
 * Interface to surp.cpp
 *
 * 20feb2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#ifndef SURP_HPP
#define SURP_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::vec v_surp( arma::vec d, arma::vec coeff, arma::vec G, bool tri ) ;
double surp( double d, arma::vec coeff, double G, bool tri ) ;
double surp_poly( double d, arma::vec coeff, double G ) ;
double surp_tri( double d, arma::vec coeff, double G ) ;
double d_surp_tri( double d, arma::vec coeff, double G ) ;

#endif
