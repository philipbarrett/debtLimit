/***********************************************************************************
* q.hpp
*
* Interface to q.cpp
*
* 13mar2017
* Philip Barrett, DC
*
***********************************************************************************/

#ifndef Q_HPP
#define Q_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::vec q_fn( arma::vec R, arma::vec p, arma::mat trans, double lambda, double phi,
                int n, std::string cont_type="avg", arma::vec G = zeros(1),
                arma::vec q_e = zeros(1), arma::mat def=zeros(1,1) ) ;

#endif
