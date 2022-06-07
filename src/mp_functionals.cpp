#include <RcppArmadillo.h>
#include "matrix_functions.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' Phase-type mixed Poisson density
//'
//' Computes the density of a phase-type Mixed Poisson distribution with parameters
//' \code{alpha} and \code{S} at \code{x}.
//'
//' @param x Non-negative integer values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @return The density at \code{x}.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mp_density(Rcpp::NumericVector x, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector density(x.size());

  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;

  arma::mat I;
  I.eye(size(S));

  double max_val{max(x)};

  std::vector<arma::mat> vect = vector_of_powers(inv(I - S), max_val);

  arma::mat aux_mat(1,1);

  for (int k{0}; k < x.size(); ++k){
    aux_mat = alpha.t() * vect[x[k]] * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return density;
}



//' Phase-type mixed Poisson density
//'
//' Computes the density of a phase-type Mixed Poisson distribution with parameters
//' \code{alpha} and \code{S} at \code{x}.
//'
//' @param x Non-negative integer values.
//' @param ex Vector of covariates effect.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @return The density at \code{x}.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mp_density_cov(Rcpp::NumericVector x, Rcpp::NumericVector ex, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector density(x.size());

  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;

  arma::mat I;
  I.eye(size(S));

  arma::mat aux_mat(1,1);

  for (int k{0}; k < x.size(); ++k){
    aux_mat = alpha.t() * matrix_power(x[k], inv(ex[k] * I - S)) * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return density;
}
