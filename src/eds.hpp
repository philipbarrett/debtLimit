/***********************************************************************************
 * eds.hpp
 *
 * Interface to eds.cpp
 *
 * 14jan2016
 * Philip Barrett, Chicago
 *
 ***********************************************************************************/

#ifndef EDS_HPP
#define EDS_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::mat keep_mat( arma::mat & X, arma::vec keep ) ;
    // Returns a submatrix of X where the rows are selected if the corresponding
    // entry of the keep vector is 1, and otherwise not.
arma::mat pcd_scaling( arma::mat & X ) ;
    // Reutrns the mean-zero, spherical unit-variance equivalent of a cloud of
    // points via PC decomposition
arma::vec eds_keep( arma::mat & X_norm, arma::vec & eps ) ;
    // Computes the vector of rows to keep in the EDS algorithm
arma::mat p_eps( arma::mat & X, arma::vec & eps ) ;
    // The EDS algorithm for a vector of distances epsilon
arma::mat p_eps_const( arma::mat & X, double eps ) ;
    // The EDS algorithm for fixed epsilon
arma::vec normal_kernel_density( arma::mat & X, double h ) ;
    // Coputes a normal kernal density given a bandwidth h
arma::vec almost_ergodic_indices( arma::mat & X, double delta, double h ) ;
    // Removes the points with cumulative probability mass less than delta
arma::mat almost_ergodic( arma::mat & X, double delta, double h ) ;
    // Removes the points with cumulative probability mass less than delta
arma::vec p_eps_cheap_idx( arma::mat & X, arma::vec & eps, double delta, double h ) ;
    // Indices of the cheap EDS algorithm with removal of low-probability points
arma::mat p_eps_cheap( arma::mat & X, arma::vec & eps, double delta, double h ) ;
    // The cheap EDS algorithm with removal of low-probability points
arma::vec p_eps_cheap_const_idx( arma::mat & X, double eps, double delta, double h ) ;
    // Indices wth constant epsilon
arma::mat p_eps_cheap_const( arma::mat & X, double eps, double delta, double h ) ;
    // With constant epsilon

#endif
