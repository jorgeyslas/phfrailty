#include <RcppArmadillo.h>
#include "matrix_functions.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Univariate case

//' Phase-type density
//'
//' Computes the density of a phase-type distribution with parameters
//' \code{alpha} and \code{S} at \code{x}.
//'
//' @param x Non-negative values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @return The density at \code{x}.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ph_density(Rcpp::NumericVector x, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector density(x.size());

  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;

  arma::mat aux_mat(1,1);

  for (int k{0}; k < x.size(); ++k){
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
//' Computes the cdf or tail of a phase-type distribution with parameters
//' \code{alpha} and \code{S} at \code{x}.
//'
//' @param x Non-negative values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @param lower_tail Cdf or tail.
//' @return The density at \code{x}.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ph_cdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, bool lower_tail = true) {
  Rcpp::NumericVector cdf(x.size());

  arma::mat e;
  e.ones(S.n_cols, 1);

  arma::mat aux_mat(1,1);

  for (int k{0}; k < x.size(); ++k) {
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


//' Laplace transform of a phase-type distribution
//'
//' Computes the Laplace transform at \code{r} of a phase-type distribution with
//'  parameters \code{alpha} and \code{S}.
//'
//' @param r Vector of real values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @return Laplace transform at \code{r}.
//' @export
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' S <- matrix(c(c(-1, 0, 0), c(1, -2, 0),c(0, 1, -5)), nrow = 3, ncol = 3)
//' ph_laplace(0.5, alpha, S)
// [[Rcpp::export]]
Rcpp::NumericVector ph_laplace(Rcpp::NumericVector r, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector laplace(r.size());

  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;

  arma::mat aux_mat(1,1);

  arma::mat identity_matrix;
  identity_matrix.eye(size(S));

  for (int i{0}; i < r.size(); ++i) {
    aux_mat = alpha.t() * inv(identity_matrix * r[i] +  S * (-1.0)) * exit_vect;
    laplace[i] = aux_mat(0,0);
  }
  return laplace;
}


//' Derivative of order n of the Laplace transform of a phase-type distribution
//' without the multiplying constant
//'
//' Computes the derivative of order n (without the multiplying constant) of the
//' Laplace transform at \code{r} of a phase-type distribution with parameters
//' \code{alpha} and \code{S}.
//'
//' @param r Vector of real values.
//' @param n An integer.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @return Laplace transform at \code{r}.
//' @export
//' @examples
//' alpha <- c(0.5, 0.3, 0.2)
//' S <- matrix(c(c(-1, 0, 0), c(1, -2, 0),c(0, 1, -5)), nrow = 3, ncol = 3)
//' ph_laplace_der_nocons(0.5, 2, alpha, S)
// [[Rcpp::export]]
Rcpp::NumericVector ph_laplace_der_nocons(Rcpp::NumericVector r, int n, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector laplace(r.size());

  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;

  arma::mat aux_mat(1,1);

  arma::mat identity_matrix;
  identity_matrix.eye(size(S));

  for (int i{0}; i < r.size(); ++i) {
    aux_mat = alpha.t() * matrix_power(n, inv(identity_matrix * r[i] +  S * (-1.0))) * exit_vect;
    laplace[i] = aux_mat(0,0);
  }
  return laplace;
}


// Multivariate case

//' Bivariate phase-type joint density
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return Joint density at \code{x}.
//' @export
//' @examples
//' alpha <- c(0.15, 0.85)
//' S11 <- matrix(c(c(-2, 9), c(0, -11)), nrow = 2, ncol = 2)
//' S12 <- matrix(c(c(2, 0), c(0, 2)), nrow = 2, ncol = 2)
//' S22 <- matrix(c(c(-1, 0), c(0.5, -5)), nrow = 2, ncol = 2)
//' x <- matrix(c(c(0.5, 1), c(2, 1.5)), ncol = 2)
//' bivph_density(x, alpha, S11, S12, S22)
// [[Rcpp::export]]
Rcpp::NumericVector bivph_density(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  long n{x.nrow()};

  Rcpp::NumericVector density(n);

  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = (S22 * (-1)) * e;

  arma::mat aux_mat(1,1);

  for (int k{0}; k < n; ++k) {
    aux_mat = alpha.t() * matrix_exponential(S11 * x(k,0)) * S12 * matrix_exponential(S22 * x(k,1)) * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return density;
}


//' Bivariate phase-type joint tail
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return Joint tail at \code{x}.
//' @export
//' @examples
//' alpha <- c(0.15, 0.85)
//' S11 <- matrix(c(c(-2, 9), c(0, -11)), nrow = 2, ncol = 2)
//' S12 <- matrix(c(c(2, 0), c(0, 2)), nrow = 2, ncol = 2)
//' S22 <- matrix(c(c(-1, 0), c(0.5, -5)), nrow = 2, ncol = 2)
//' x <- matrix(c(c(0.5, 1), c(2, 1.5)), ncol = 2)
//' bivph_tail(x, alpha, S11, S12, S22)
// [[Rcpp::export]]
Rcpp::NumericVector bivph_tail(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  long n{x.nrow()};

  Rcpp::NumericVector tail(n);

  arma::mat e;
  e.ones(S22.n_cols, 1);

  arma::mat aux_mat(1,1);

  for (int k{0}; k < n; ++k) {
    aux_mat = alpha.t() * inv(S11 * (-1)) * matrix_exponential(S11 * x(k,0)) * S12 * matrix_exponential(S22 * x(k,1)) * e;
    tail[k] = aux_mat(0,0);
  }
  return tail;
}


//' Bivariate phase-type joint Laplace
//'
//' @param r Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return Joint laplace at \code{r}.
//' @export
//' @examples
//' alpha <- c(0.15, 0.85)
//' S11 <- matrix(c(c(-2, 9), c(0, -11)), nrow = 2, ncol = 2)
//' S12 <- matrix(c(c(2, 0), c(0, 2)), nrow = 2, ncol = 2)
//' S22 <- matrix(c(c(-1, 0), c(0.5, -5)), nrow = 2, ncol = 2)
//' x <- matrix(c(c(0.5, 1), c(2, 1.5)), ncol = 2)
//' bivph_laplace(x, alpha, S11, S12, S22)
// [[Rcpp::export]]
Rcpp::NumericVector bivph_laplace(Rcpp::NumericMatrix r, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  long n{r.nrow()};

  Rcpp::NumericVector laplace(n);

  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = (S22 * (-1)) * e;

  arma::mat identity_matrix1;
  identity_matrix1.eye(size(S11));

  arma::mat identity_matrix2;
  identity_matrix2.eye(size(S22));

  arma::mat aux_mat(1,1);

  for (int k{0}; k < n; ++k) {
    aux_mat = alpha.t() * inv(identity_matrix1 * r(k,0) +  S11 * (-1.0)) * S12 * inv(identity_matrix2 * r(k,1) +  S22 * (-1.0)) * exit_vect;
    laplace[k] = aux_mat(0,0);
  }
  return laplace;
}


//' Derivative of order (n,m) of the joint Laplace of a bivariate phase-type
//'
//' @param r Matrix of values.
//' @param n Order of first component.
//' @param m Order of second component.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return Derivative of joint laplace at \code{r}, without multiplicative constants.
//' @export
//' @examples
//' alpha <- c(0.15, 0.85)
//' S11 <- matrix(c(c(-2, 9), c(0, -11)), nrow = 2, ncol = 2)
//' S12 <- matrix(c(c(2, 0), c(0, 2)), nrow = 2, ncol = 2)
//' S22 <- matrix(c(c(-1, 0), c(0.5, -5)), nrow = 2, ncol = 2)
//' x <- matrix(c(c(0.5, 1), c(2, 1.5)), ncol = 2)
//' bivph_laplace_der_nocons(x, 2, 3, alpha, S11, S12, S22)
// [[Rcpp::export]]
Rcpp::NumericVector bivph_laplace_der_nocons(Rcpp::NumericMatrix r, int n, int m, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  long r_size{r.nrow()};

  Rcpp::NumericVector laplace_der(r_size);

  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = (S22 * (-1)) * e;

  arma::mat identity_matrix1;
  identity_matrix1.eye(size(S11));

  arma::mat identity_matrix2;
  identity_matrix2.eye(size(S22));

  arma::mat aux_mat(1,1);

  for (int k{0}; k < r_size; ++k) {
    aux_mat = alpha.t() * matrix_power(n, inv(identity_matrix1 * r(k,0) +  S11 * (-1.0))) * S12 * matrix_power(m, inv(identity_matrix2 * r(k,1) +  S22 * (-1.0))) * exit_vect;
    laplace_der[k] = aux_mat(0,0);
  }
  return laplace_der;
}

