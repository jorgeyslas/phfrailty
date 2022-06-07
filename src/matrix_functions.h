#include <RcppArmadillo.h>

double inf_norm(arma::mat A);
arma::mat matrix_exponential(arma::mat A);
arma::mat matrix_power(int n, arma::mat A);
void vector_of_matrices(std::vector<arma::mat> & the_vector, const arma::mat & S, int size_vec);
std::vector<arma::mat> vector_of_powers(const arma::mat & A, int vect_size);
arma::mat matrix_vanloan(const arma::mat & A1, const arma::mat & A2, const arma::mat & B1);
