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
               arma::vec An, arma::vec Cn, arma::mat def ) ;
arma::mat zed_2( arma::vec p, arma::vec d, List params,
                 arma::vec An, arma::vec Bn, arma::vec Cn, arma::mat def ) ;
arma::mat zed_2_num( arma::vec p, arma::vec d, List params,
                     arma::vec Cn, arma::vec An, arma::mat def) ;
double zed_2_num_d_i( arma::vec p, arma::vec d, List params, int i,
                      arma::vec An, arma::vec Cn, arma::mat def ) ;
arma::mat zed_2_ana( arma::vec p, arma::vec d, List params, arma::vec An,
                     arma::vec Bn, arma::vec Cn, arma::mat def) ;

#endif
