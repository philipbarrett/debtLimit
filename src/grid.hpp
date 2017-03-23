/***********************************************************************************
 * grid.hpp
 *
 * Interface to grid.cpp
 *
 * 22mar2017
 * Philip Barrett, DC
 *
 ***********************************************************************************/

#ifndef GRID_HPP
#define GRID_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::vec d_grid( arma::vec d, double x_sd, double x_sd_max, int n_pts ) ;
arma::vec e_grid_fn( double x_sd, int n_pts ) ;
arma::mat Q_init( arma::vec d_grid, arma::vec d, arma::vec R ) ;

#endif
