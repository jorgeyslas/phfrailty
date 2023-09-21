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



//' Phase-type mixed Poisson density with covariates
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
    aux_mat = alpha.t() * matrix_power(x[k] + 1, inv(ex[k] * I - S)) * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return density;
}



//' Phase-type mixed Poisson density with covariates
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
Rcpp::List mp_aux_density(Rcpp::NumericVector x, Rcpp::NumericVector ex, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector density(x.size());
  Rcpp::NumericVector dens_aux(x.size());

  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;

  arma::mat I;
  I.eye(size(S));

  arma::mat aux_mat(1,1);

  arma::mat temp_mat;

  for (int k{0}; k < x.size(); ++k){
    temp_mat = matrix_power(x[k] + 2, inv(ex[k] * I - S));
    aux_mat = alpha.t() * temp_mat * exit_vect;
    dens_aux[k] = aux_mat(0,0);
    aux_mat = alpha.t() * temp_mat * (ex[k] * I - S) * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return Rcpp::List::create(
    Rcpp::Named("density") = density,
    Rcpp::Named("dens_aux") = dens_aux
  );
}



// Multivariate case

//' Correlated Phase-type mixed Poisson density
//'
//' Computes the joint density of a correlated phase-type Mixed Poisson distribution
//' with parameters \code{alpha} and \code{S} at \code{x}.
//'
//' @param x Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return The density at \code{x}.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mp_cor_dens(Rcpp::NumericMatrix x, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  Rcpp::NumericVector density(x.nrow());

  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = (S22 * (-1)) * e;

  arma::mat I1;
  I1.eye(size(S11));

  arma::mat I2;
  I2.eye(size(S22));

  double max_val1{max(x.column(0))};
  std::vector<arma::mat> vect1 = vector_of_powers(inv(I1 - S11), max_val1);

  double max_val2{max(x.column(1))};
  std::vector<arma::mat> vect2 = vector_of_powers(inv(I2 - S22), max_val2);

  arma::mat aux_mat(1,1);

  for (int k{0}; k < x.nrow(); ++k){
    aux_mat = alpha.t() * vect1[x(k, 0)] * S12 * vect2[x(k, 1)] * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return density;
}


//' Correlated Phase-type mixed Poisson density with covariates
//'
//' Computes the joint density of a correlated phase-type Mixed Poisson distribution
//' with parameters \code{alpha} and \code{S} at \code{x}.
//'
//' @param x Matrix of values.
//' @param ex Matrix of covariates effect.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return The density at \code{x}.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mp_cor_dens_cov(Rcpp::NumericMatrix x, Rcpp::NumericMatrix ex, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  Rcpp::NumericVector density(x.nrow());

  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = (S22 * (-1)) * e;

  arma::mat I1;
  I1.eye(size(S11));

  arma::mat I2;
  I2.eye(size(S22));

  arma::mat aux_mat(1,1);

  for (int k{0}; k < x.nrow(); ++k){
    //aux_mat = alpha.t() *  matrix_power(x(k, 0) + 1, inv(ex(k, 0) * I1 - S11)) * S12 * matrix_power(x(k, 1) + 1, inv(ex(k, 1) * I2 - S22)) * exit_vect;
    aux_mat = alpha.t() *  matrix_power(x(k, 0) + 1, inv(I1 - S11 / ex(k, 0))) * S12 * matrix_power(x(k, 1) + 1, inv(I2 - S22 / ex(k, 1))) * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return density;
}



//' Correlated Phase-type mixed Poisson density with covariates
//'
//' Computes the joint density of a correlated phase-type Mixed Poisson distribution
//' with parameters \code{alpha} and \code{S} at \code{x}.
//'
//' @param x Matrix of values.
//' @param ex Matrix of covariates effect.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return The density at \code{x}.
//' @export
// [[Rcpp::export]]
Rcpp::List mp_cor_dens_aux(Rcpp::NumericMatrix x, Rcpp::NumericMatrix ex, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  Rcpp::NumericVector density(x.nrow());
  Rcpp::NumericVector dens_aux1(x.nrow());
  Rcpp::NumericVector dens_aux2(x.nrow());

  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = (S22 * (-1)) * e;

  arma::mat I1;
  I1.eye(size(S11));

  arma::mat I2;
  I2.eye(size(S22));

  arma::mat aux_mat(1,1);

  arma::mat temp_mat1;
  arma::mat temp_mat2;

  for (int k{0}; k < x.nrow(); ++k){
    //temp_mat1 = matrix_power(x(k, 0) + 2, inv(ex(k, 0) * I1 - S11));
    temp_mat1 = matrix_power(x(k, 0) + 2, inv(I1 - S11 / ex(k, 0)));
    //temp_mat2 = matrix_power(x(k, 1) + 2, inv(ex(k, 1) * I2 - S22));
    temp_mat2 = matrix_power(x(k, 1) + 2, inv(I2 - S22 / ex(k, 1)));
    //aux_mat = alpha.t() * temp_mat1 * S12 * temp_mat2 * (ex(k, 1) * I2 - S22)  * exit_vect;
    aux_mat = alpha.t() * temp_mat1 * S12 * temp_mat2 * (I2 - S22 / ex(k, 1))  * exit_vect;
    dens_aux1[k] = aux_mat(0,0);
    //aux_mat = alpha.t() * temp_mat1 * (ex(k, 0) * I1 - S11) * S12 * temp_mat2 * exit_vect;
    aux_mat = alpha.t() * temp_mat1 * (I1 - S11 / ex(k, 0)) * S12 * temp_mat2 * exit_vect;
    dens_aux2[k] = aux_mat(0,0);
    //aux_mat = alpha.t() * temp_mat1 * (ex(k, 0) * I1 - S11) * S12 * temp_mat2 * (ex(k, 1) * I2 - S22) * exit_vect;
    aux_mat = alpha.t() * temp_mat1 * (I1 - S11 / ex(k, 0)) * S12 * temp_mat2 * (I2 - S22 / ex(k, 1)) * exit_vect;
    density[k] = aux_mat(0,0);
  }
  return Rcpp::List::create(
    Rcpp::Named("density") = density,
    Rcpp::Named("dens_aux1") = dens_aux1,
    Rcpp::Named("dens_aux2") = dens_aux2
  );
}
