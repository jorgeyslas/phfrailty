#include <RcppArmadillo.h>
#include "m_exp.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' Phase-type density
//'
//' Computes the density of a phase-type distribution with parameters \code{alpha} and \code{S} at \code{x}
//' @param x non-negative value
//' @param alpha vector of initial probabilities
//' @param S sub-intensity matrix
//' @return The density at \code{x}
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ph_density(Rcpp::NumericVector x, arma::vec alpha, arma::mat S) {

  Rcpp::NumericVector density(x.size());

  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;

  arma::mat aux_mat(1,1);

  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      density[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * x[k]) * exit_vect;
      density[k] = aux_mat(0,0);
    }
  }
  return density;
}


//' Phase-type cdf or tail
//'
//' Computes the cdf or tail of a phase-type distribution with parameters \code{alpha} and \code{S} at \code{x}
//' @param x non-negative values
//' @param alpha vector of initial probabilities
//' @param S sub-intensity matrix
//' @param lower_tail cdf or tail
//' @return The density at \code{x}
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ph_cdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, bool lower_tail = true) {

  Rcpp::NumericVector cdf(x.size());

  arma::mat e; e.ones(S.n_cols, 1);

  arma::mat aux_mat(1,1);

  for (int k = 0; k < x.size(); ++k){
    if (x[k] == 0) {
      aux_mat = alpha.t() * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
    else {
      aux_mat = alpha.t() * matrix_exponential(S * x[k]) * e;
      cdf[k] = 1.0 - aux_mat(0,0);
    }
  }
  if (lower_tail == true) {
    return cdf;
  }
  else {
    return (1 - cdf);
  }
}
