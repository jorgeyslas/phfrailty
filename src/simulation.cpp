#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Embedded Markov chain of a sub-intensity matrix
//'
//' Returns the transition probabilities of the embedded Markov chain determined the sub-intensity matrix
//' @param S sub-intensity matrix
//' @return the embedded Markov chain
//'
// [[Rcpp::export]]
arma::mat embedded_mc(arma::mat S) {
  unsigned p{S.n_rows};
  arma::mat Q(p + 1, p + 1);

  arma::mat e; e.ones(S.n_cols, 1);
  arma::mat t = (S * (-1)) * e;

  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < p + 1; ++j) {
      if (j != i && j < p) {
        Q(i,j) = -1.0 * S(i,j) / S(i,i);
      }
      else if(j == p) {
        Q(i,j) = -1.0 * t(i,0) / S(i,i);
      }
    }
  }
  Q(p,p) = 1;

  return Q;
}

//' Cumulate matrix
//'
//' Creates a new matrix with entries the cumulated rows of \code{A}
//' @param A a matrix
//' @return the cumulated matrix
//'
// [[Rcpp::export]]
arma::mat cumulate_matrix(arma::mat A) {
  unsigned p1{A.n_rows};
  unsigned p2{A.n_cols};

  arma::mat cumulated(p1, p2);

  for (int i{0}; i < p1; ++i) {
    for (int j{0}; j < p2; ++j) {
      if (j == 0) {
        cumulated(i,j) = A(i,j);
      }
      else {
        cumulated(i,j) = cumulated(i,j - 1) + A(i,j);
      }
    }
  }
  return cumulated;
}

//' Cumulate vector
//'
//' Creates a new vector with entries the cumulated entries of \code{A}
//' @param A a vector
//' @return the cumulated vector
//'
// [[Rcpp::export]]
arma::vec cumulate_vector(arma::vec A) {
  unsigned p{A.size()};

  arma::vec cumulated(p);

  for (int i{0}; i < p; ++i) {
    if (i == 0) {
      cumulated[i] = A[i];
    }
    else {
      cumulated[i] = cumulated[i - 1] + A[i];
    }
  }
  return cumulated;
}

//' Initial state of Markov jump process
//'
//' Given the accumulated values of the initial probabilities \code{Pi} and a uniform value \code{u}, it returns the initial state of a Markov jump process
//' This corresponds to the states satisfying cum_pi_(k-1)<u<cum_pi_(k)
//' @param cum_pi a vector
//' @param u random value in (0,1)
//' @return initial state of the Markov jump process
//'
// [[Rcpp::export]]
long initial_state(arma::vec cum_pi, double u) {
  if (u <= cum_pi[0]) {
    return 0;
  }

  for( int i{1}; i < cum_pi.size(); ++i) {
    if (cum_pi[i - 1] < u && u <= cum_pi[i]) {
      return i;
    }
  }
  return 0;
}

//' New state in a Markov jump process
//'
//' Given a transition matrix \code{Q}, a uniform value \code{u}, and a previous state \code{k}, it returns the new state of a Markov jump process
//' @param previous_state previous state of the Markov jump process
//' @param cum_embedded_mc transition matrix
//' @param u random value in (0,1)
//' @return next state of the Markov jump process
//'
// [[Rcpp::export]]
long new_state(long previous_state, arma::mat cum_embedded_mc, double u) {
  if (u <= cum_embedded_mc(previous_state,0)) {
    return 0;
  }

  for (int i{1}; i < cum_embedded_mc.n_cols; ++i) {
    if (cum_embedded_mc(previous_state,i - 1) < u && u <= cum_embedded_mc(previous_state,i)) {
      return i;
    }
  }

  return 0;
}


//' Random phase-type
//'
//' Generates a sample of size \code{n} from a phase-type distribution with parameters \code{alpha} and \code{S}
//' @param n sample size
//' @param alpha vector of initial probabilities
//' @param S sub-intensity matrix
//' @return simulated sample
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericVector rphasetype(int n, arma::vec alpha, arma::mat S) {

  Rcpp::NumericVector sample(n);

  arma::mat cum_embedded_mc = cumulate_matrix(embedded_mc(S));
  arma::vec cum_pi = cumulate_vector(alpha);

  unsigned p{alpha.size()};
  long state{0};
  for (int i{0}; i < n; ++i) {
    double time{0.0};
    state = initial_state(cum_pi, Rcpp::runif(1)[0]);
    while (state != p) {
      time += log(1.0 - Rcpp::runif(1)[0]) / S(state,state);
      state = new_state(state, cum_embedded_mc, Rcpp::runif(1)[0]);
    }
    sample[i] = time;
  }
  return sample;
}
