// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// inf_norm
double inf_norm(arma::mat A);
RcppExport SEXP _phfrailty_inf_norm(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(inf_norm(A));
    return rcpp_result_gen;
END_RCPP
}
// matrix_exponential
arma::mat matrix_exponential(arma::mat A);
RcppExport SEXP _phfrailty_matrix_exponential(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_exponential(A));
    return rcpp_result_gen;
END_RCPP
}
// ph_density
Rcpp::NumericVector ph_density(Rcpp::NumericVector x, arma::vec alpha, arma::mat S);
RcppExport SEXP _phfrailty_ph_density(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(ph_density(x, alpha, S));
    return rcpp_result_gen;
END_RCPP
}
// ph_cdf
Rcpp::NumericVector ph_cdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, bool lower_tail);
RcppExport SEXP _phfrailty_ph_cdf(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(ph_cdf(x, alpha, S, lower_tail));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_phfrailty_inf_norm", (DL_FUNC) &_phfrailty_inf_norm, 1},
    {"_phfrailty_matrix_exponential", (DL_FUNC) &_phfrailty_matrix_exponential, 1},
    {"_phfrailty_ph_density", (DL_FUNC) &_phfrailty_ph_density, 3},
    {"_phfrailty_ph_cdf", (DL_FUNC) &_phfrailty_ph_cdf, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_phfrailty(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
