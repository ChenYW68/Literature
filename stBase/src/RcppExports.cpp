// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Trace_Matrix
double Trace_Matrix(arma::colvec& h, arma::mat& S);
RcppExport SEXP _stBase_Trace_Matrix(SEXP hSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(Trace_Matrix(h, S));
    return rcpp_result_gen;
END_RCPP
}
// Trace_Muti
double Trace_Muti(arma::mat S1, arma::mat S2);
RcppExport SEXP _stBase_Trace_Muti(SEXP S1SEXP, SEXP S2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S2(S2SEXP);
    rcpp_result_gen = Rcpp::wrap(Trace_Muti(S1, S2));
    return rcpp_result_gen;
END_RCPP
}
// Solve
arma::mat Solve(arma::mat A);
RcppExport SEXP _stBase_Solve(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(Solve(A));
    return rcpp_result_gen;
END_RCPP
}
// theta2_fun_C
SEXP theta2_fun_C(SEXP m_r, SEXP n_r, SEXP theta2_r, SEXP ds_r, SEXP D_r, SEXP bandKernel_r, SEXP Q_r, SEXP S00_r, SEXP S01_r, SEXP theta1_r, SEXP sigTheta1_r, SEXP nThreads_r);
RcppExport SEXP _stBase_theta2_fun_C(SEXP m_rSEXP, SEXP n_rSEXP, SEXP theta2_rSEXP, SEXP ds_rSEXP, SEXP D_rSEXP, SEXP bandKernel_rSEXP, SEXP Q_rSEXP, SEXP S00_rSEXP, SEXP S01_rSEXP, SEXP theta1_rSEXP, SEXP sigTheta1_rSEXP, SEXP nThreads_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type m_r(m_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n_r(n_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type theta2_r(theta2_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ds_r(ds_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type D_r(D_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type bandKernel_r(bandKernel_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Q_r(Q_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type S00_r(S00_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type S01_r(S01_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type theta1_r(theta1_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sigTheta1_r(sigTheta1_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nThreads_r(nThreads_rSEXP);
    rcpp_result_gen = Rcpp::wrap(theta2_fun_C(m_r, n_r, theta2_r, ds_r, D_r, bandKernel_r, Q_r, S00_r, S01_r, theta1_r, sigTheta1_r, nThreads_r));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stBase_Trace_Matrix", (DL_FUNC) &_stBase_Trace_Matrix, 2},
    {"_stBase_Trace_Muti", (DL_FUNC) &_stBase_Trace_Muti, 2},
    {"_stBase_Solve", (DL_FUNC) &_stBase_Solve, 1},
    {"_stBase_theta2_fun_C", (DL_FUNC) &_stBase_theta2_fun_C, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_stBase(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
