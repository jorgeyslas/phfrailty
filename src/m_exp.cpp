#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' L_inf norm of a matrix
//'
//' Computes the inf norm of a matrix A, defined as
//' L-oo A =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
//' @param A a matrix
//' @return The norm
//'
// [[Rcpp::export]]
double inf_norm(arma::mat A) {
  double value{0.0};
  for (int i{0}; i < A.n_rows; ++i) {
    double row_sum{0.0};
    for (int j{0}; j < A.n_cols; ++j) {
      row_sum += abs(A(i,j));
    }
    value = std::max(value, row_sum);
  }
  return value;
}



//' Matrix exponential algorithm
//'
//' MATLAB's built-in algorithm - Pade approximation
//' @param A a matrix
//'
// [[Rcpp::export]]
arma::mat matrix_exponential(arma::mat A) {

  const int q{6};

  arma::mat matrixAuxiliar(A.n_rows,A.n_cols);
  arma::mat ExpM(A.n_rows,A.n_cols);

  double aNorm{inf_norm(A)};

  int ee{static_cast<int>(log2(aNorm)) + 1};

  int s{std::max(0, ee + 1)};

  double t{1.0 / pow(2.0, s)};

  arma::mat a2 = A * t;
  arma::mat x = a2;

  double c{0.5};

  ExpM.eye(size(A));

  ExpM = ExpM + (a2 * c);

  arma::mat d;
  d.eye(size(A));

  d = (d + a2 * (-c));

  int p{1};

  for (int k{2}; k <= q; ++k) {
    c = c * static_cast<double>(q - k + 1) / static_cast<double>(k * (2 * q - k + 1));

    x = (a2 * x);

    ExpM = (x * c) + ExpM;

    if (p) {
      d = (x * c) + d;
    }
    else {
      d = (x * (-c)) + d;
    }
    p = !p;
  }

  ExpM = inv(d) * ExpM;

  for (int k = 1; k <= s; ++k) {
    ExpM = ExpM * ExpM;
  }
  return(ExpM);
}
