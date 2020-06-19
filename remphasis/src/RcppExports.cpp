// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_mcem
List rcpp_mcem(const std::vector<double>& brts, const std::vector<double>& init_pars, int sample_size, int maxN, const std::string& plugin, int soc, int max_misssing, double max_lambda, const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, double xtol_rel, int num_threads, bool copy_trees);
RcppExport SEXP _remphasis_rcpp_mcem(SEXP brtsSEXP, SEXP init_parsSEXP, SEXP sample_sizeSEXP, SEXP maxNSEXP, SEXP pluginSEXP, SEXP socSEXP, SEXP max_misssingSEXP, SEXP max_lambdaSEXP, SEXP lower_boundSEXP, SEXP upper_boundSEXP, SEXP xtol_relSEXP, SEXP num_threadsSEXP, SEXP copy_treesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type brts(brtsSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type init_pars(init_parsSEXP);
    Rcpp::traits::input_parameter< int >::type sample_size(sample_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type maxN(maxNSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type plugin(pluginSEXP);
    Rcpp::traits::input_parameter< int >::type soc(socSEXP);
    Rcpp::traits::input_parameter< int >::type max_misssing(max_misssingSEXP);
    Rcpp::traits::input_parameter< double >::type max_lambda(max_lambdaSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type lower_bound(lower_boundSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type upper_bound(upper_boundSEXP);
    Rcpp::traits::input_parameter< double >::type xtol_rel(xtol_relSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type copy_trees(copy_treesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mcem(brts, init_pars, sample_size, maxN, plugin, soc, max_misssing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads, copy_trees));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_remphasis_rcpp_mcem", (DL_FUNC) &_remphasis_rcpp_mcem, 13},
    {NULL, NULL, 0}
};

void remphasis_init(DllInfo *dll);
RcppExport void R_init_remphasis(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    remphasis_init(dll);
}
