# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Discretization of a univariate density
#'
#' Discretizates a univariate density function using Simpson's rule.
#'
#' @param density The density function.
#' @param parameters Parameters of the density function.
#' @param ini_point Initial value for the discretization.
#' @param truncation_point Max value for the discretization.
#' @param max_probability Maximum probability of a discrete point.
#' @param max_deltat Maximum size of interval.
#' @return List with values and weights.
discretizate_density <- function(density, parameters, ini_point, truncation_point, max_deltat, max_probability) {
    .Call(`_phfrailty_discretizate_density`, density, parameters, ini_point, truncation_point, max_deltat, max_probability)
}

#' Discretization of a bivariate density
#'
#' Discretizates a bivariate density function using Simpson's rule.
#'
#' @param density The density function.
#' @param parameters Parameters of the density function.
#' @param ini_point1 Max value for the discretization - first component.
#' @param truncation_point1 Max value for the discretization - first component.
#' @param max_deltat1 Maximum size of interval - first component.
#' @param ini_point2 Max value for the discretization - second component.
#' @param truncation_point2 Max value for the discretization - second component.
#' @param max_deltat2 Maximum size of interval - second component.
#' @param max_probability Maximum probability of a discrete point.
#' @return List with values and weights.
discretizate_bivdensity <- function(density, parameters, ini_point1, truncation_point1, max_deltat1, ini_point2, truncation_point2, max_deltat2, max_probability) {
    .Call(`_phfrailty_discretizate_bivdensity`, density, parameters, ini_point1, truncation_point1, max_deltat1, ini_point2, truncation_point2, max_deltat2, max_probability)
}

#' L_inf norm of a matrix
#'
#' Computes the inf norm of a matrix A, defined as
#'  L_inf A =  max (1 <= I <= M) sum (1 <= J <= N) abs (A(I,J)).
#'
#' @param A A matrix.
#' @return The inf norm of A.
#'
inf_norm <- function(A) {
    .Call(`_phfrailty_inf_norm`, A)
}

#' Matrix exponential
#'
#' MATLAB's built-in algorithm for matrix exponential - Pade approximation.
#'
#' @param A A matrix.
#' @return exp(A).
#'
matrix_exponential <- function(A) {
    .Call(`_phfrailty_matrix_exponential`, A)
}

#' Computes A^n
#'
#' @param A A matrix.
#' @param n An integer.
#' @return A^n.
#' @export
matrix_power <- function(n, A) {
    .Call(`_phfrailty_matrix_power`, n, A)
}

#' Computes elements S^n until the value size_vec
#'
#' @param the_vector A vector to save results.
#' @param S Sub-intensity matrix.
#' @param size_vec Size of vector.
#' @return Modified vector with the elements S^n.
vector_of_matrices <- function(the_vector, S, size_vec) {
    invisible(.Call(`_phfrailty_vector_of_matrices`, the_vector, S, size_vec))
}

#' Computes elements A^(n + 1) until the given size
#'
#' @param A A matrix.
#' @param vect_size Size of vector.
#'
vector_of_powers <- function(A, vect_size) {
    .Call(`_phfrailty_vector_of_powers`, A, vect_size)
}

#' Creates the matrix  (A1, B1 ; 0, A2)
#'
#' @param A1 A matrix.
#' @param A2 A matrix.
#' @param B1 A matrix.
#' @return The matrix (A1, B1 ; 0, A2).
#'
matrix_vanloan <- function(A1, A2, B1) {
    .Call(`_phfrailty_matrix_vanloan`, A1, A2, B1)
}

#' Phase-type mixed Poisson density
#'
#' Computes the density of a phase-type Mixed Poisson distribution with parameters
#' \code{alpha} and \code{S} at \code{x}.
#'
#' @param x Non-negative integer values.
#' @param alpha Vector of initial probabilities.
#' @param S Sub-intensity matrix.
#' @return The density at \code{x}.
#' @export
mp_density <- function(x, alpha, S) {
    .Call(`_phfrailty_mp_density`, x, alpha, S)
}

#' Phase-type mixed Poisson density with covariates
#'
#' Computes the density of a phase-type Mixed Poisson distribution with parameters
#' \code{alpha} and \code{S} at \code{x}.
#'
#' @param x Non-negative integer values.
#' @param ex Vector of covariates effect.
#' @param alpha Vector of initial probabilities.
#' @param S Sub-intensity matrix.
#' @return The density at \code{x}.
#' @export
mp_density_cov <- function(x, ex, alpha, S) {
    .Call(`_phfrailty_mp_density_cov`, x, ex, alpha, S)
}

#' Phase-type mixed Poisson density with covariates
#'
#' Computes the density of a phase-type Mixed Poisson distribution with parameters
#' \code{alpha} and \code{S} at \code{x}.
#'
#' @param x Non-negative integer values.
#' @param ex Vector of covariates effect.
#' @param alpha Vector of initial probabilities.
#' @param S Sub-intensity matrix.
#' @return The density at \code{x}.
#' @export
mp_aux_density <- function(x, ex, alpha, S) {
    .Call(`_phfrailty_mp_aux_density`, x, ex, alpha, S)
}

#' Correlated Phase-type mixed Poisson density
#'
#' Computes the joint density of a correlated phase-type Mixed Poisson distribution
#' with parameters \code{alpha} and \code{S} at \code{x}.
#'
#' @param x Matrix of values.
#' @param alpha Vector of initial probabilities.
#' @param S11 Sub-intensity matrix.
#' @param S12 Matrix.
#' @param S22 Sub-intensity matrix.
#' @return The density at \code{x}.
#' @export
mp_cor_dens <- function(x, alpha, S11, S12, S22) {
    .Call(`_phfrailty_mp_cor_dens`, x, alpha, S11, S12, S22)
}

#' Correlated Phase-type mixed Poisson density with covariates
#'
#' Computes the joint density of a correlated phase-type Mixed Poisson distribution
#' with parameters \code{alpha} and \code{S} at \code{x}.
#'
#' @param x Matrix of values.
#' @param ex Matrix of covariates effect.
#' @param alpha Vector of initial probabilities.
#' @param S11 Sub-intensity matrix.
#' @param S12 Matrix.
#' @param S22 Sub-intensity matrix.
#' @return The density at \code{x}.
#' @export
mp_cor_dens_cov <- function(x, ex, alpha, S11, S12, S22) {
    .Call(`_phfrailty_mp_cor_dens_cov`, x, ex, alpha, S11, S12, S22)
}

#' Correlated Phase-type mixed Poisson density with covariates
#'
#' Computes the joint density of a correlated phase-type Mixed Poisson distribution
#' with parameters \code{alpha} and \code{S} at \code{x}.
#'
#' @param x Matrix of values.
#' @param ex Matrix of covariates effect.
#' @param alpha Vector of initial probabilities.
#' @param S11 Sub-intensity matrix.
#' @param S12 Matrix.
#' @param S22 Sub-intensity matrix.
#' @return The density at \code{x}.
#' @export
mp_cor_dens_aux <- function(x, ex, alpha, S11, S12, S22) {
    .Call(`_phfrailty_mp_cor_dens_aux`, x, ex, alpha, S11, S12, S22)
}

#' EM for phase-type distributions using Pade for matrix exponential
#'
#' @param alpha Initial probabilities.
#' @param S Sub-intensity.
#' @param obs The observations.
#' @param weight The weights for the observations.
#' @return Fitted alpha and S after one iteration.
#' @export
EMstep <- function(alpha, S, obs, weight) {
    invisible(.Call(`_phfrailty_EMstep`, alpha, S, obs, weight))
}

#' EM for bivariate phase-type distributions using Pade for matrix exponential
#'
#' @param alpha Initial probabilities.
#' @param S11 Sub-intensity.
#' @param S12 A matrix.
#' @param S22 Sub-intensity.
#' @param obs The observations.
#' @param weight The weights for the observations.
#' @return Fitted alpha, S11, S12 and S22 after one iteration.
#' @export
EMstep_bivph <- function(alpha, S11, S12, S22, obs, weight) {
    invisible(.Call(`_phfrailty_EMstep_bivph`, alpha, S11, S12, S22, obs, weight))
}

#' Phase-type density
#'
#' Computes the density of a phase-type distribution with parameters
#' \code{alpha} and \code{S} at \code{x}.
#'
#' @param x Non-negative values.
#' @param alpha Vector of initial probabilities.
#' @param S Sub-intensity matrix.
#' @return The density at \code{x}.
#' @export
ph_density <- function(x, alpha, S) {
    .Call(`_phfrailty_ph_density`, x, alpha, S)
}

#' Phase-type cdf or tail
#'
#' Computes the cdf or tail of a phase-type distribution with parameters
#' \code{alpha} and \code{S} at \code{x}.
#'
#' @param x Non-negative values.
#' @param alpha Vector of initial probabilities.
#' @param S Sub-intensity matrix.
#' @param lower_tail Cdf or tail.
#' @return The density at \code{x}.
#' @export
ph_cdf <- function(x, alpha, S, lower_tail = TRUE) {
    .Call(`_phfrailty_ph_cdf`, x, alpha, S, lower_tail)
}

#' Laplace transform of a phase-type distribution
#'
#' Computes the Laplace transform at \code{r} of a phase-type distribution with
#'  parameters \code{alpha} and \code{S}.
#'
#' @param r Vector of real values.
#' @param alpha Vector of initial probabilities.
#' @param S Sub-intensity matrix.
#' @return Laplace transform at \code{r}.
#' @export
#' @examples
#' alpha <- c(0.5, 0.3, 0.2)
#' S <- matrix(c(c(-1, 0, 0), c(1, -2, 0),c(0, 1, -5)), nrow = 3, ncol = 3)
#' ph_laplace(0.5, alpha, S)
ph_laplace <- function(r, alpha, S) {
    .Call(`_phfrailty_ph_laplace`, r, alpha, S)
}

#' Derivative of order n of the Laplace transform of a phase-type distribution
#' without the multiplying constant
#'
#' Computes the derivative of order n (without the multiplying constant) of the
#' Laplace transform at \code{r} of a phase-type distribution with parameters
#' \code{alpha} and \code{S}.
#'
#' @param r Vector of real values.
#' @param n An integer.
#' @param alpha Vector of initial probabilities.
#' @param S Sub-intensity matrix.
#' @return Laplace transform at \code{r}.
#' @export
#' @examples
#' alpha <- c(0.5, 0.3, 0.2)
#' S <- matrix(c(c(-1, 0, 0), c(1, -2, 0),c(0, 1, -5)), nrow = 3, ncol = 3)
#' ph_laplace_der_nocons(0.5, 2, alpha, S)
ph_laplace_der_nocons <- function(r, n, alpha, S) {
    .Call(`_phfrailty_ph_laplace_der_nocons`, r, n, alpha, S)
}

#' Bivariate phase-type joint density
#'
#' @param x Matrix of values.
#' @param alpha Vector of initial probabilities.
#' @param S11 Sub-intensity matrix.
#' @param S12 Matrix.
#' @param S22 Sub-intensity matrix.
#' @return Joint density at \code{x}.
#' @export
#' @examples
#' alpha <- c(0.15, 0.85)
#' S11 <- matrix(c(c(-2, 9), c(0, -11)), nrow = 2, ncol = 2)
#' S12 <- matrix(c(c(2, 0), c(0, 2)), nrow = 2, ncol = 2)
#' S22 <- matrix(c(c(-1, 0), c(0.5, -5)), nrow = 2, ncol = 2)
#' x <- matrix(c(c(0.5, 1), c(2, 1.5)), ncol = 2)
#' bivph_density(x, alpha, S11, S12, S22)
bivph_density <- function(x, alpha, S11, S12, S22) {
    .Call(`_phfrailty_bivph_density`, x, alpha, S11, S12, S22)
}

#' Bivariate phase-type joint tail
#'
#' @param x Matrix of values.
#' @param alpha Vector of initial probabilities.
#' @param S11 Sub-intensity matrix.
#' @param S12 Matrix.
#' @param S22 Sub-intensity matrix.
#' @return Joint tail at \code{x}.
#' @export
#' @examples
#' alpha <- c(0.15, 0.85)
#' S11 <- matrix(c(c(-2, 9), c(0, -11)), nrow = 2, ncol = 2)
#' S12 <- matrix(c(c(2, 0), c(0, 2)), nrow = 2, ncol = 2)
#' S22 <- matrix(c(c(-1, 0), c(0.5, -5)), nrow = 2, ncol = 2)
#' x <- matrix(c(c(0.5, 1), c(2, 1.5)), ncol = 2)
#' bivph_tail(x, alpha, S11, S12, S22)
bivph_tail <- function(x, alpha, S11, S12, S22) {
    .Call(`_phfrailty_bivph_tail`, x, alpha, S11, S12, S22)
}

#' Bivariate phase-type joint Laplace
#'
#' @param r Matrix of values.
#' @param alpha Vector of initial probabilities.
#' @param S11 Sub-intensity matrix.
#' @param S12 Matrix.
#' @param S22 Sub-intensity matrix.
#' @return Joint laplace at \code{r}.
#' @export
#' @examples
#' alpha <- c(0.15, 0.85)
#' S11 <- matrix(c(c(-2, 9), c(0, -11)), nrow = 2, ncol = 2)
#' S12 <- matrix(c(c(2, 0), c(0, 2)), nrow = 2, ncol = 2)
#' S22 <- matrix(c(c(-1, 0), c(0.5, -5)), nrow = 2, ncol = 2)
#' x <- matrix(c(c(0.5, 1), c(2, 1.5)), ncol = 2)
#' bivph_laplace(x, alpha, S11, S12, S22)
bivph_laplace <- function(r, alpha, S11, S12, S22) {
    .Call(`_phfrailty_bivph_laplace`, r, alpha, S11, S12, S22)
}

#' Derivative of order (n,m) of the joint Laplace of a bivariate phase-type
#'
#' @param r Matrix of values.
#' @param n Order of first component.
#' @param m Order of second component.
#' @param alpha Vector of initial probabilities.
#' @param S11 Sub-intensity matrix.
#' @param S12 Matrix.
#' @param S22 Sub-intensity matrix.
#' @return Derivative of joint laplace at \code{r}, without multiplicative constants.
#' @export
#' @examples
#' alpha <- c(0.15, 0.85)
#' S11 <- matrix(c(c(-2, 9), c(0, -11)), nrow = 2, ncol = 2)
#' S12 <- matrix(c(c(2, 0), c(0, 2)), nrow = 2, ncol = 2)
#' S22 <- matrix(c(c(-1, 0), c(0.5, -5)), nrow = 2, ncol = 2)
#' x <- matrix(c(c(0.5, 1), c(2, 1.5)), ncol = 2)
#' bivph_laplace_der_nocons(x, 2, 3, alpha, S11, S12, S22)
bivph_laplace_der_nocons <- function(r, n, m, alpha, S11, S12, S22) {
    .Call(`_phfrailty_bivph_laplace_der_nocons`, r, n, m, alpha, S11, S12, S22)
}

#' Embedded Markov chain of a sub-intensity matrix
#'
#' Returns the transition probabilities of the embedded Markov chain determined
#'  the sub-intensity matrix.
#'
#' @param S Sub-intensity matrix.
#' @return The embedded Markov chain.
#'
embedded_mc <- function(S) {
    .Call(`_phfrailty_embedded_mc`, S)
}

#' Cumulate matrix
#'
#' Creates a new matrix with entries the cumulated rows of \code{A}.
#'
#' @param A A matrix.
#' @return The cumulated matrix.
#'
cumulate_matrix <- function(A) {
    .Call(`_phfrailty_cumulate_matrix`, A)
}

#' Cumulate vector
#'
#' Creates a new vector with entries the cumulated entries of \code{A}.
#'
#' @param A A vector.
#' @return The cumulated vector.
#'
cumulate_vector <- function(A) {
    .Call(`_phfrailty_cumulate_vector`, A)
}

#' Initial state of Markov jump process
#'
#' Given the accumulated values of the initial probabilities \code{Pi} and a
#' uniform value \code{u}, it returns the initial state of a Markov jump process.
#' This corresponds to the states satisfying cum_alpha_(k-1)<u<cum_alpha_(k).
#'
#' @param cum_alpha A vector.
#' @param u Random value in (0,1).
#' @return Initial state of the Markov jump process.
#'
initial_state <- function(cum_alpha, u) {
    .Call(`_phfrailty_initial_state`, cum_alpha, u)
}

#' New state in a Markov jump process
#'
#' Given a transition matrix \code{Q}, a uniform value \code{u}, and a previous
#' state \code{k}, it returns the new state of a Markov jump process.
#'
#' @param prev_state Previous state of the Markov jump process.
#' @param cum_embedded_mc Transition matrix.
#' @param u Random value in (0,1).
#' @return Next state of the Markov jump process.
#'
new_state <- function(prev_state, cum_embedded_mc, u) {
    .Call(`_phfrailty_new_state`, prev_state, cum_embedded_mc, u)
}

#' Simulate phase-type
#'
#' Generates a sample of size \code{n} from a phase-type distribution with
#' parameters \code{alpha} and \code{S}.
#'
#' @param n Sample size.
#' @param alpha Vector of initial probabilities.
#' @param S Sub-intensity matrix.
#' @return Simulated sample.
#' @export
#'
rphasetype <- function(n, alpha, S) {
    .Call(`_phfrailty_rphasetype`, n, alpha, S)
}

#' Simulate a MPH* random vector
#'
#' Generates a sample of size \code{n} from a MPH* distribution with parameters
#'  \code{alpha}, \code{S} and \code{R}.
#'
#' @param n Sample size.
#' @param alpha Initial probabilities.
#' @param S Sub-intensity matrix.
#' @param R Reward matrix.
#' @return The simulated sample.
#' @export
#' @examples
#' alpha <- c(0.5, 0.3, 0.2)
#' S <- matrix(c(c(-1, 0, 0), c(1, -2, 0), c(0, 1, -5)), nrow = 3, ncol = 3)
#' R <- matrix(c(c(1, 0, 0.8), c(0, 1, 0.2)), nrow = 3, ncol = 2)
#' n <- 10
#' rmph(n, alpha, S, R)
rmph <- function(n, alpha, S, R) {
    .Call(`_phfrailty_rmph`, n, alpha, S, R)
}

#' Clone a vector
#'
#' @param v A vector.
#' @return A clone of the vector.
#' @export
clone_vector <- function(v) {
    .Call(`_phfrailty_clone_vector`, v)
}

#' Clone a matrix
#'
#' @param m A matrix.
#' @return A clone of the matrix.
#' @export
clone_matrix <- function(m) {
    .Call(`_phfrailty_clone_matrix`, m)
}

#' Random structure of a phase-type
#'
#' Generates random parameters \code{alpha} and \code{S} of a phase-type
#' distribution of dimension \code{p} with chosen structure.
#'
#' @param p Dimension of the phase-type.
#' @param structure Type of structure: "general", "hyperexponential", "gerlang",
#'  "coxian" or "gcoxian".
#' @param scale_factor A factor that multiplies the sub-intensity matrix.
#' @return Random parameters \code{alpha} and \code{S} of a phase-type.
#'
random_structure <- function(p, structure = "general", scale_factor = 1) {
    .Call(`_phfrailty_random_structure`, p, structure, scale_factor)
}

#' Random structure of a bivariate phase-type
#'
#' Generates random parameters \code{alpha}, \code{S11}, \code{S12}, and \code{S22}
#' of a bivariate phase-type distribution of dimension \code{p  = p1 + p2}.
#'
#' @param p1 Dimension of the first block.
#' @param p2 Dimension of the second block.
#' @param scale_factor A factor that multiplies the sub-intensity matrix.
#' @return Random parameters  \code{alpha}, \code{S11}, \code{S12}, and \code{S22}
#'  of a bivariate phase-type.
#'
random_structure_bivph <- function(p1, p2, scale_factor = 1) {
    .Call(`_phfrailty_random_structure_bivph`, p1, p2, scale_factor)
}

#' Merges the matrices S11, S12 and S22 into a sub-intensity matrix
#'
#' @param S11 A sub-intensity matrix.
#' @param S12 A matrix.
#' @param S22 A sub-intensity matrix.
#' @return A sub-intensity matrix.
#'
merge_matrices <- function(S11, S12, S22) {
    .Call(`_phfrailty_merge_matrices`, S11, S12, S22)
}

