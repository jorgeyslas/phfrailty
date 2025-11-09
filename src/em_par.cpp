// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "matrix_functions.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// declarations you already have:
// arma::mat matrix_exponential(const arma::mat& A);
// arma::mat matrix_vanloan(const arma::mat& A, const arma::mat& B, const arma::mat& C);

//' EM for bivariate phase-type distributions using Pade for matrix exponential
//'
//' @param alpha Initial probabilities.
//' @param S11 Sub-intensity.
//' @param S12 A matrix.
//' @param S22 Sub-intensity.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' @return Fitted alpha, S11, S12 and S22 after one iteration.
//' @export
// [[Rcpp::export]]
void EMstep_bivph_omp(arma::vec &alpha,
                      arma::mat &S11,
                      arma::mat &S12,
                      arma::mat &S22,
                      const Rcpp::NumericMatrix &obs,
                      const Rcpp::NumericVector &weight) {

  const unsigned p1 = S11.n_rows;
  const unsigned p2 = S22.n_rows;
  const unsigned p  = p1 + p2;

  const std::size_t nObs = obs.nrow();
  const double eps = 1e-14;

  // Global sufficient statistics (to be reduced)
  arma::vec Bmean(p1, arma::fill::zeros);
  arma::vec Zmean(p,  arma::fill::zeros);
  arma::mat Nmean(p,  p + 1, arma::fill::zeros);
  double sum_weights = 0.0;

  // ----- Parallel section -----
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#pragma omp parallel
{
  // Thread-local accumulators
  arma::vec Bmean_loc(p1, arma::fill::zeros);
  arma::vec Zmean_loc(p,  arma::fill::zeros);
  arma::mat Nmean_loc(p,  p + 1, arma::fill::zeros);
  double sum_w_loc = 0.0;

  // Pre-allocated workspaces (thread-local)
  arma::mat g1(2 * p1, 2 * p1);
  arma::mat g2(2 * p2, 2 * p2);

  arma::mat aux_exp1(p1, p1, arma::fill::none);
  arma::mat cmatrix1(p1, p1, arma::fill::none);
  arma::mat aux_exp2(p2, p2, arma::fill::none);
  arma::mat cmatrix2(p2, p2, arma::fill::none);

  arma::vec ones_p2(p2, arma::fill::ones);
  arma::vec exit_vec = -S22 * ones_p2;   // read-only in E-step

  arma::rowvec alphaT = alpha.t();

#pragma omp for nowait
  for (std::size_t k = 0; k < nObs; ++k) {

    const double wgt = weight[k];
    sum_w_loc += wgt;

    const double t1 = obs(k, 0);
    const double t2 = obs(k, 1);

    // exp(S22 * t2)
    aux_exp2 = matrix_exponential(S22 * t2);

    // b1 = S12 * exp(S22 t2) * t2 * alpha^T
    arma::mat bmatrix1 = S12 * aux_exp2 * exit_vec * alphaT;

    // g1 = exp( [ S11  b1 ; 0  S11 ] * t1 )
    g1 = matrix_exponential(matrix_vanloan(S11, S11, bmatrix1) * t1);

    aux_exp1 = g1.submat(0,        0,        p1 - 1, p1 - 1);
    cmatrix1 = g1.submat(0,        p1,       p1 - 1, 2 * p1 - 1);

    // b2 = t2 * alpha^T * exp(S11 t1) * S12
    arma::mat bmatrix2 = exit_vec * alphaT * aux_exp1 * S12;

    // g2 = exp( [ S22  b2 ; 0  S22 ] * t2 )
    g2 = matrix_exponential(matrix_vanloan(S22, S22, bmatrix2) * t2);
    cmatrix2 = g2.submat(0, p2, p2 - 1, 2 * p2 - 1);

    // density
    const double dens = std::max(
      arma::as_scalar(alphaT * aux_exp1 * S12 * aux_exp2 * exit_vec), eps
    );
    const double w_over_d = wgt / dens;

    // helpers
    arma::vec    aux_matrix1 = aux_exp1 * S12 * aux_exp2 * exit_vec;     // p1×1
    arma::mat    aux_matrix2 = aux_exp2 * exit_vec * alphaT * aux_exp1;  // p2×p1
    arma::rowvec aux_matrix3 = alphaT * aux_exp1 * S12 * aux_exp2;       // 1×p2

    // Accumulate locals
    Bmean_loc += (alpha % aux_matrix1) * w_over_d;

    Zmean_loc.subvec(0, p1 - 1)   += cmatrix1.diag() * w_over_d;
    Zmean_loc.subvec(p1, p - 1)   += cmatrix2.diag() * w_over_d;

    arma::mat cm1T = cmatrix1.t();
    Nmean_loc.submat(0, 0, p1 - 1, p1 - 1) += (S11 % cm1T) * w_over_d;

    Nmean_loc.submat(0, p1, p1 - 1, p1 + p2 - 1) += (S12 % aux_matrix2.t()) * w_over_d;

    Nmean_loc.submat(p1, p, p - 1, p) += (aux_matrix3.t() % exit_vec) * w_over_d;

    arma::mat cm2T = cmatrix2.t();
    Nmean_loc.submat(p1, p1, p - 1, p1 + p2 - 1) += (S22 % cm2T) * w_over_d;
  }

  // Merge locals
#pragma omp critical
{
  Bmean  += Bmean_loc;
  Zmean  += Zmean_loc;
  Nmean  += Nmean_loc;
  sum_weights += sum_w_loc;
}
} // end parallel
#else
// Fallback: call your single-threaded EMstep_bivph here if desired
Rcpp::stop("OpenMP not available; compile with OpenMP to use EMstep_bivph_omp.");
#endif

// ----- M-step (same as before) -----
arma::vec exit_vec(p2,  arma::fill::zeros);

for (int i{0}; i < p1; ++i) {
  alpha[i] = Bmean[i] / sum_weights;
  if (alpha[i] < 0) {
    alpha[i] = 0;
  }
  S11(i,i) = 0;
  for (int j{0}; j < p1; ++j) {
    if (i != j) {
      S11(i,j) = Nmean(i,j) / Zmean(i,0);
      if (S11(i,j) < 0) {
        S11(i,j) = 0;
      }
      S11(i,i) -= S11(i,j);
    }
  }
  for (int j{0}; j < p2; ++j) {
    S12(i,j) = Nmean(i,j + p1) / Zmean(i,0);
    if (S12(i,j) < 0){
      S12(i,j) = 0;
    }
    S11(i,i) -= S12(i,j);
  }
}
for (int i{0}; i < p2; ++i) {
  exit_vec(i,0) = Nmean(i + p1,p) / Zmean(i + p1,0);
  if (exit_vec(i,0) < 0) {
    exit_vec(i,0) = 0;
  }
  S22(i,i) = -exit_vec(i,0);
  for (int j{0}; j < p2; ++j) {
    if (i != j) {
      S22(i,j) = Nmean(i + p1,j + p1) / Zmean(i + p1,0);
      if (S22(i,j) < 0){
        S22(i,j) = 0;
      }
      S22(i,i) -= S22(i,j);
    }
  }
}
}
