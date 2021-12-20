#include <RcppArmadillo.h>

double inf_norm(arma::mat A);
arma::mat matrix_exponential(arma::mat A);
arma::mat matrix_power(int n, arma::mat A);
void vector_of_matrices(std::vector<arma::mat> & the_vector, const arma::mat & S, int size_vec);
arma::mat matrix_vanloan(const arma::mat & A1, const arma::mat & A2, const arma::mat & B1);
