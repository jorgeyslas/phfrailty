#include <RcppArmadillo.h>
#include "matrix_functions.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Univariate case

//' EM for phase-type distributions using Pade for matrix exponential
//'
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' @return Fitted alpha and S after one iteration.
//' @export
// [[Rcpp::export]]
void EMstep(arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight) {
  unsigned p{S.n_rows};

  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;

  arma::mat Bmean = arma::zeros(p,1);
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);

  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);

  arma::mat aux_mat(1,1);

  arma::mat g(2 * p,2 * p);
  arma::mat s_prod_alpha(p,p);
  s_prod_alpha = exit_vect * alpha.t();

  g = matrix_vanloan(S, S, s_prod_alpha);

  double g_norm{inf_norm(g)};

  std::vector<arma::mat> the_vector;
  vector_of_matrices(the_vector, g, 6);

  arma::mat X(2 * p,2 * p);
  arma::mat d(2 * p,2 * p);

  const int q{6};
  int s{};
  double xmod{};
  double c{};

  double sum_weights{0.0};
  double density{0.0};

  //E-step
  for (int k{0}; k < obs.size(); ++k) {
    sum_weights += weight[k];

    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(g_norm  * obs[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = obs[k] / pow(2.0, s);
    c = 0.5;
    X = the_vector[1] * (c * xmod);

    g = the_vector[0] + X;
    d = the_vector[0] - X;

    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = the_vector[l] * (c * pow(xmod,l));
      g = g + X;
      if (pind) {
        d =  d + X;
      }
      else {
        d = d - X;
      }
      pind = !pind;
    }
    g = inv(d) * g;
    for (int l = 1; l <= s; ++l) {
      g = g * g;
    }

    // Separate matrix
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = g(i,j);
        cmatrix(i,j) = g(i,j + p);
      }
    }

    avector = alpha.t() * aux_exp;
    bvector = aux_exp * exit_vect;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);

    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * exit_vect(i,0) * weight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
      }
    }
  }

  // M step
  for (int i{0}; i < p; ++i) {
    alpha[i] = Bmean(i,0) / sum_weights;
    if (alpha[i] < 0) {
      alpha[i] = 0;
    }
    exit_vect(i,0) = Nmean(i,p) / Zmean(i,0);
    if (exit_vect(i,0) < 0) {
      exit_vect(i,0) = 0;
    }
    S(i,i) = -exit_vect(i,0);
    for (int j{0}; j < p; ++j) {
      if (i != j) {
        S(i,j) = Nmean(i,j) / Zmean(i,0);
        if (S(i,j) < 0) {
          S(i,j) = 0;
        }
        S(i,i) -= S(i,j);
      }
    }
  }
}


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
void EMstep_bivph(arma::vec & alpha, arma::mat & S11, arma::mat & S12, arma::mat & S22, const Rcpp::NumericMatrix & obs, const Rcpp::NumericVector & weight) {
  unsigned p1{S11.n_rows};
  unsigned p2{S22.n_rows};
  unsigned p{p1 + p2};

  double density{0};
  double aux;

  arma::mat Bmean = arma::zeros(p1,1);
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);

  arma::mat bmatrix1(p1,p1);
  arma::mat cmatrix1(p1,p1);
  arma::mat aux_exp1(p1,p1);

  arma::mat bmatrix2(p2,p2);
  arma::mat cmatrix2(p2,p2);
  arma::mat aux_exp2(p2,p2);

  arma::mat g1(2 * p1,2 * p1);
  arma::mat g2(2 * p2,2 * p2);

  arma::mat e;
  e.ones(p2, 1);
  arma::mat exit_vec = (S22 * (-1)) * e;

  arma::mat aux_matrix1(p1,1);
  arma::mat aux_matrix2(p2,p1);
  arma::mat aux_matrix3(1,p2);

  arma::mat aux_mat(1,1);

  double sum_weights{0.0};
  //E step
  for (int k{0}; k < obs.nrow(); ++k) {
    sum_weights += weight[k];

    aux_exp2 = matrix_exponential(S22 * obs(k,1));
    bmatrix1 = S12 * aux_exp2 * exit_vec * alpha.t();

    g1 = matrix_exponential(matrix_vanloan(S11, S11, bmatrix1) * obs(k,0));

    for (int i{0}; i < p1; ++i) {
      for (int j{0}; j < p1; ++j) {
        aux_exp1(i,j) = g1(i,j);
        cmatrix1(i,j) = g1(i,j + p1);
      }
    }

    bmatrix2 = exit_vec * alpha.t() * aux_exp1 * S12;

    g2 = matrix_exponential(matrix_vanloan(S22, S22, bmatrix2) * obs(k,1));

    for (int i{0}; i < p2; ++i) {
      for (int j{0}; j < p2; ++j) {
        cmatrix2(i,j) = g2(i,j + p2);
      }
    }

    aux_mat = alpha.t() * aux_exp1 * S12 * aux_exp2 * exit_vec;
    density = aux_mat(0,0);

    //E-step
    aux_matrix1 = aux_exp1 * S12 * aux_exp2 * exit_vec;
    aux_matrix2 =  aux_exp2 * exit_vec * alpha.t() * aux_exp1;

    for (int i{0}; i < p1; ++i) {
      aux = aux_matrix1(i,0);
      Bmean(i,0) += alpha[i] * aux * weight[k] / density;
      Zmean(i,0) += cmatrix1(i,i) * weight[k] / density;
      for (int j{0}; j < p1; ++j) {
        Nmean(i,j) += S11(i,j) * cmatrix1(j,i) * weight[k] / density;
      }
      for (int j{0}; j < p2; ++j) {
        aux = aux_matrix2(j,i);
        Nmean(i,j + p1) += S12(i,j) * aux * weight[k] / density;
      }
    }

    aux_matrix3 = alpha.t() * aux_exp1 * S12 * aux_exp2;

    for (int i{0}; i < p2; ++i) {
      Zmean(i + p1,0) += cmatrix2(i,i) * weight[k] / density;
      aux = aux_matrix3(0,i);
      Nmean(i + p1,p) += aux * exit_vec(i,0) * weight[k] / density;
      for (int j{0}; j < p2; ++j){
        Nmean(i + p1,j + p1) += S22(i,j) * cmatrix2(j,i) * weight[k] / density;
      }
    }
  }

  // M step
  for (int i{0}; i < p1; ++i) {
    alpha[i] = Bmean(i,0) / sum_weights;
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
void EMstep_bivph_gpt(arma::vec &alpha, arma::mat &S11, arma::mat &S12, arma::mat &S22,
                  const Rcpp::NumericMatrix &obs,
                  const Rcpp::NumericVector &weight) {

  const unsigned p1 = S11.n_rows;
  const unsigned p2 = S22.n_rows;
  const unsigned p  = p1 + p2;

  const std::size_t nObs = obs.nrow();
  const double eps = 1e-14;

  // Sufficient statistics
  arma::vec Bmean(p1, arma::fill::zeros);          // initial-state responsibilities
  arma::vec Zmean(p,  arma::fill::zeros);          // expected sojourn times
  arma::mat Nmean(p,  p + 1, arma::fill::zeros);   // expected jumps (last col: exits from block 2)

  // Pre-allocate workspaces
  arma::mat g1(2 * p1, 2 * p1);
  arma::mat g2(2 * p2, 2 * p2);

  arma::mat aux_exp1(p1, p1, arma::fill::none);    // exp(S11 * t1)
  arma::mat cmatrix1(p1, p1, arma::fill::none);    // ∫ exp(S11 * s) b exp(S11 * (t1-s)) ds (UR block)
  arma::mat aux_exp2(p2, p2, arma::fill::none);    // exp(S22 * t2)
  arma::mat cmatrix2(p2, p2, arma::fill::none);

  arma::vec ones_p2(p2, arma::fill::ones);
  arma::vec exit_vec = -S22 * ones_p2;            // t2 = -S22*1

  arma::rowvec alphaT = alpha.t();                // reuse

  double sum_weights = 0.0;

  // ----- E-step over observations -----
  // #pragma omp parallel  // optional: add OpenMP with local accumulators & critical section to merge
  for (std::size_t k = 0; k < nObs; ++k) {

    const double wgt = weight[k];
    sum_weights += wgt;

    const double t1 = obs(k, 0);
    const double t2 = obs(k, 1);

    // exp(S22 * t2)
    aux_exp2 = matrix_exponential(S22 * t2);

    // b1 = S12 * exp(S22 t2) * t2 * alpha^T
    arma::mat bmatrix1 = S12 * aux_exp2 * exit_vec * alphaT;

    // g1 = exp( [ S11  b1 ; 0  S11 ] * t1 )
    g1 = matrix_exponential(matrix_vanloan(S11, S11, bmatrix1) * t1);

    // Extract blocks without loops
    aux_exp1 = g1.submat(0,        0,        p1 - 1, p1 - 1);
    cmatrix1 = g1.submat(0,        p1,       p1 - 1, 2 * p1 - 1);

    // b2 = t2 * alpha^T * exp(S11 t1) * S12
    arma::mat bmatrix2 = exit_vec * alphaT * aux_exp1 * S12;

    // g2 = exp( [ S22  b2 ; 0  S22 ] * t2 )
    g2 = matrix_exponential(matrix_vanloan(S22, S22, bmatrix2) * t2);
    cmatrix2 = g2.submat(0, p2, p2 - 1, 2 * p2 - 1);

    // Density: alpha^T exp(S11 t1) S12 exp(S22 t2) t2
    const double density = std::max(
      arma::as_scalar(alphaT * aux_exp1 * S12 * aux_exp2 * exit_vec),
      eps
    );

    // Helpers reused in updates
    arma::vec  aux_matrix1 = aux_exp1 * S12 * aux_exp2 * exit_vec;          // p1×1
    arma::mat  aux_matrix2 = aux_exp2 * exit_vec * alphaT * aux_exp1;       // p2×p1
    arma::rowvec aux_matrix3 = alphaT * aux_exp1 * S12 * aux_exp2;          // 1×p2

    // ---- accumulate sufficient stats (vectorized where possible) ----

    // Block 1 (rows 0..p1-1)
    // Bmean: initial responsibilities (componentwise)
    Bmean += (alpha % aux_matrix1) * (wgt / density);

    // Zmean for block 1: take diagonal of cmatrix1
    Zmean.subvec(0, p1 - 1) += cmatrix1.diag() * (wgt / density);

    // Nmean for S11 off-diagonals: elementwise S11(i,j) * cmatrix1(j,i)
    // Use transpose once
    arma::mat cm1T = cmatrix1.t();  // p1×p1
    Nmean.submat(0, 0, p1 - 1, p1 - 1) += (S11 % cm1T) * (wgt / density);

    // Nmean for S12: elementwise S12(i,j) * aux_matrix2(j,i)
    Nmean.submat(0, p1, p1 - 1, p1 + p2 - 1) += (S12 % aux_matrix2.t()) * (wgt / density);

    // Block 2 (rows p1..p1+p2-1)
    Zmean.subvec(p1, p - 1) += cmatrix2.diag() * (wgt / density);

    // Exit counts (last column p): elementwise aux_matrix3(0,i) * exit_vec(i)
    Nmean.submat(p1, p, p - 1, p) += (aux_matrix3.t() % exit_vec) * (wgt / density);

    // Nmean for S22
    arma::mat cm2T = cmatrix2.t();
    Nmean.submat(p1, p1, p - 1, p1 + p2 - 1) += (S22 % cm2T) * (wgt / density);
  }

  // ----- M-step -----
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
