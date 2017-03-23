/***********************************************************************************
* ziter.hpp
*
* Interface to ziter.cpp
*
* 23mar2017
* Philip Barrett, DC
*
***********************************************************************************/

#ifndef ZITER_HPP
#define ZITER_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::mat ziter( arma::mat P, arma::vec d_bar, arma::mat QHat, arma::mat Q,
                 arma::vec d_grid, arma::vec e_grid,
                 bool tri, arma::mat D_prime_0, bool D_prime_0_flag,
                 List params, arma::vec An, arma::vec Cn,
                 int print_level, double tol, int maxit,
                 double q_tol, int q_maxit, double d_tol, int d_maxit ) ;

#endif
