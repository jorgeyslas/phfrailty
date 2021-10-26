#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' L_inf norm of a matrix
//'
//' Computes the inf norm of a matrix A, defined as
//' L-oo A =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
//' @param A A matrix.
//' @return The inf norm of A.
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
//' MATLAB's built-in algorithm - Pade approximation.
//' @param A A matrix.
//' @return exp(A).
//'
// [[Rcpp::export]]
arma::mat matrix_exponential(arma::mat A) {
  const int q{6};

  arma::mat matrixAuxiliar(A.n_rows,A.n_cols);
  arma::mat mexp(A.n_rows,A.n_cols);

  double a_norm{inf_norm(A)};

  int ee{static_cast<int>(log2(a_norm)) + 1};

  int s{std::max(0, ee + 1)};

  double t{1.0 / pow(2.0, s)};

  arma::mat a2 = A * t;
  arma::mat x = a2;

  double c{0.5};

  mexp.eye(size(A));

  mexp = mexp + (a2 * c);

  arma::mat d;
  d.eye(size(A));

  d = (d + a2 * (-c));

  int p{1};

  for (int k{2}; k <= q; ++k) {
    c = c * static_cast<double>(q - k + 1) / static_cast<double>(k * (2 * q - k + 1));

    x = (a2 * x);

    mexp = (x * c) + mexp;

    if (p) {
      d = (x * c) + d;
    }
    else {
      d = (x * (-c)) + d;
    }
    p = !p;
  }

  mexp = inv(d) * mexp;

  for (int k{1}; k <= s; ++k) {
    mexp = mexp * mexp;
  }
  return mexp;
}



//' Computes A^n
//'
//' @param A A matrix.
//' @param n An integer.
//' @return A^n.
//' @export
// [[Rcpp::export]]
arma::mat matrix_power(int n, arma::mat A) {
  if (n == 1) {
    return A;
  }
  else if (n == 2) {
    return A * A;
  }
  arma::mat previous_matrix = A * A;
  arma::mat new_matrix = A * previous_matrix;
  for (int i{4}; i <= n; ++i) {
    previous_matrix = new_matrix;
    new_matrix = A * previous_matrix;
  }
  return new_matrix;
}


//' Computes elements S^n until the value size
//'
//' @param the_vector A vector to save results.
//' @param S Sub-intensity matrix.
//' @param size_vec Size of vector.
//' @return Modified vector with the elements S^n.
// [[Rcpp::export]]
void vector_of_matrices(std::vector<arma::mat> & the_vector, const arma::mat & S, int size_vec) {
  arma::mat Id;
  Id.eye(size(S));

  the_vector.push_back(Id);

  for (int k{1}; k <= size_vec; ++k) {
    the_vector.push_back( S * the_vector[k - 1]);
  }
}


//' Creates the matrix  (A1, B1 ; 0, A2)
//'
//' @param A1 A matrix.
//' @param A2 A matrix.
//' @param B1 A matrix.
//' @return The matrix (A1, B1 ; 0, A2).
//'
// [[Rcpp::export]]
arma::mat matrix_VanLoan(const arma::mat & A1, const arma::mat & A2, const arma::mat & B1) {
  unsigned p1{A1.n_cols};
  unsigned p2{A2.n_cols};
  unsigned p{p1 + p2};

  arma::mat aux_matrix(p, p);

  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < p; ++j) {
      if (i < p1 && j < p1) {
        aux_matrix(i,j) = A1(i,j);
      }
      else if (i >= p1 && j < p1) {
        aux_matrix(i,j) = 0;
      }
      else if (i < p1 && j >= p1) {
        aux_matrix(i,j) = B1(i,j - p1);
      }
      else {
        aux_matrix(i,j) = A2(i - p1,j - p1);
      }
    }
  }
  return aux_matrix;
}
