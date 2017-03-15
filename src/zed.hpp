/***********************************************************************************
 * zed.hpp
 *
 * Interface to zed.cpp
 *
 * 20feb2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#ifndef ZED_HPP
#define ZED_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;
using namespace R ;

arma::vec zed( arma::vec p, arma::vec d, List params,
               arma::vec qd, arma::vec qe, arma::mat def ) ;
arma::mat zed_2( arma::vec p, arma::vec d, List params,
                 arma::vec qd, arma::vec qe, arma::mat def) ;
arma::mat zed_2_num( arma::vec p, arma::vec d, List params,
                     arma::vec qd, arma::vec qe, arma::mat def) ;
double zed_2_num_d_i( arma::vec p, arma::vec d, List params, int i,
                      arma::vec qd, arma::vec qe, arma::mat def ) ;

#endif
