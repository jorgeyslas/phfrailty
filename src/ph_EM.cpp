#include <RcppArmadillo.h>
#include "matrix_functions.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Univariate case

//' EM using Pade for matrix exponential
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

  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;

  arma::mat Bmean = arma::zeros(p,1);
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);

  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);

  arma::mat aux_mat(1,1);

  arma::mat J(2 * p,2 * p);
  arma::mat t_prod_alpha(p,p);
  t_prod_alpha = t * alpha.t();

  J = matrix_VanLoan(S, S, t_prod_alpha);

  double JNorm{inf_norm(J)};

  std::vector<arma::mat> the_vector;
  vector_of_matrices(the_vector, J, 6);

  arma::mat X(2 * p,2 * p);
  arma::mat D(2 * p,2 * p);

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
    int ee{static_cast<int>(log2(JNorm  * obs[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = obs[k] / pow(2.0, s);
    c = 0.5;
    X = the_vector[1] * (c * xmod);

    J = the_vector[0] + X;
    D = the_vector[0] - X;

    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = the_vector[l] * (c * pow(xmod,l));
      J = J + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    J = inv(D) * J;
    for (int l = 1; l <= s; ++l) {
      J = J * J;
    }

    // Separate matrix
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }

    avector = alpha.t() * aux_exp;
    bvector = aux_exp * t;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);

    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * t(i,0) * weight[k] / density;
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
    t(i,0) = Nmean(i,p) / Zmean(i,0);
    if (t(i,0) < 0) {
      t(i,0) = 0;
    }
    S(i,i) = -t(i,0);
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



//' EM for bivariate PH using Pade for matrix exponential
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

  arma::mat J1(2 * p1,2 * p1);
  arma::mat J2(2 * p2,2 * p2);

  arma::mat e; e.ones(p2, 1);
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

    J1 = matrix_exponential(matrix_VanLoan(S11, S11, bmatrix1) * obs(k,0));

    for (int i{0}; i < p1; ++i) {
      for (int j{0}; j < p1; ++j) {
        aux_exp1(i,j) = J1(i,j);
        cmatrix1(i,j) = J1(i,j + p1);
      }
    }

    bmatrix2 = exit_vec * alpha.t() * aux_exp1 * S12;

    J2 = matrix_exponential(matrix_VanLoan(S22, S22, bmatrix2) * obs(k,1));

    for (int i{0}; i < p2; ++i) {
      for (int j{0}; j < p2; ++j) {
        cmatrix2(i,j) = J2(i,j + p2);
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

