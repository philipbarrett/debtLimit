/***********************************************************************************
 * qiter.hpp
 *
 * Interface to qiter.cpp
 *
 * 22mar2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#ifndef QITER_HPP
#define QITER_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::mat d_prime( int i_x, double d, arma::vec d_bar, double qhat, arma::mat Q,
                   arma::vec d_grid, arma::vec G, double lambda, arma::vec e_grid,
                   arma::vec coeff, bool tri,
                   arma::mat D_prime_0, bool D_prime_0_flag, bool verbose,
                   double tol, int maxit ) ;
arma::vec q_e( double d, arma::vec d_bar, double qhat, arma::mat Q,
               arma::vec d_grid, arma::vec G, double lambda, arma::vec e_grid,
               arma::vec coeff, bool tri, arma::mat D_prime_0, bool D_prime_0_flag,
               arma::mat trans, bool verbose, double tol, int maxit ) ;

#endif
