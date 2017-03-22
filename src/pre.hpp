/***********************************************************************************
* pre.hpp
*
* Interface to pre.cpp
*
* 21marb2017
* Philip Barrett, DC
*
***********************************************************************************/

#ifndef PRE_HPP
#define PRE_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

double p_init_d_i( List params, arma::vec p, arma::vec d, arma::vec An,
                   arma::vec Bn, arma::vec Cn, arma::mat def, int i ) ;
arma::vec p_init_d( List params, arma::vec p, arma::vec d, arma::vec An,
                   arma::vec Bn, arma::vec Cn, arma::mat def ) ;

#endif
