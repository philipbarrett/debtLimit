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
                   arma::mat D_prime_0, bool D_prime_0_flag, int print_level,
                   double tol, int maxit ) ;
arma::vec q_e( double d, arma::vec d_bar, arma::vec qhat, arma::mat Q,
               arma::vec d_grid, arma::vec G, double lambda, arma::vec e_grid,
               arma::vec coeff, bool tri, arma::mat D_prime_0, bool D_prime_0_flag,
               arma::mat trans, int print_level, double tol, int maxit ) ;
arma::vec q_hat_fn( double d, arma::vec p, arma::vec d_bar, arma::vec qhat, arma::mat Q,
                    arma::vec d_grid, arma::vec R, arma::vec G, double lambda, double phi,
                    arma::vec e_grid, arma::vec coeff, bool tri, arma::mat D_prime_0,
                    bool D_prime_0_flag, arma::mat trans, int print_level, double tol,
                    int maxit, double d_tol, int d_maxit ) ;
arma::mat q_hat_mat( arma::mat P, arma::vec d_bar, arma::mat QHat, arma::mat Q,
                     arma::vec d_grid, arma::vec R, arma::vec G, double lambda, double phi, arma::vec e_grid,
                     arma::vec coeff, bool tri, arma::mat D_prime_0, bool D_prime_0_flag,
                     arma::mat trans, int print_level, double tol, int maxit,
                     double d_tol, int d_maxit ) ;
arma::mat d_prime_mat( arma::vec d_bar, arma::mat QHat, arma::mat Q,
                       arma::vec d_grid, arma::vec G, double lambda, arma::vec e_grid,
                       arma::mat trans, arma::vec coeff, bool tri,
                       arma::mat D_prime_0, bool D_prime_0_flag, int print_level,
                       double tol, int maxit ) ;
arma::mat qe_mat( arma::vec d_bar, arma::mat QHat, arma::mat Q,
                  arma::vec d_grid, arma::vec G, double lambda, arma::vec e_grid,
                  arma::mat trans, arma::vec coeff, bool tri,
                  arma::mat D_prime_0, bool D_prime_0_flag, int print_level,
                  double tol, int maxit ) ;

#endif
