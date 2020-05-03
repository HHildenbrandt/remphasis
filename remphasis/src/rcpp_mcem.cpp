// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "emphasis.hpp"
#include "plugin.hpp"
using namespace Rcpp;


DataFrame unpack(const emphasis::tree_t& tree)
{
  NumericVector brts, n, t_ext;
  for (const emphasis::node_t& node : tree) {
    brts.push_back(node.brts);
    n.push_back(node.n);
    t_ext.push_back(node.t_ext);
  }
  return DataFrame::create(Named("brts") = brts, Named("n") = n, Named("t_ext") = t_ext);
}


// [[Rcpp::export(name = "em_cpp")]]
List rcpp_mcem(const NumericVector& brts_i,       
               const NumericVector& init_pars_i,      
               int sample_size,                       
               const std::string& plugin,             
               int soc = 2,                           
               int max_misssing = 10000,               
               double max_lambda = 500.0,             
               const NumericVector& lower_bound_i = NumericVector(),  
               const NumericVector& upper_bound_i = NumericVector(),  
               double xtol = 0.001,                     
               int num_threads = 0,
               bool copy_trees = false,
               bool verbose = true) 
{
  std::vector<double> brts = Rcpp::as<std::vector<double>>(brts_i);
  std::vector<double> init_pars = Rcpp::as<std::vector<double>>(init_pars_i);
  std::vector<double> lower_bound = Rcpp::as<std::vector<double>>(lower_bound_i);
  std::vector<double> upper_bound = Rcpp::as<std::vector<double>>(upper_bound_i);
  
  auto model = emphasis::create_plugin_model(plugin);
  init_pars.resize(model->nparams(), 0.0);
  if (verbose) {
    Rcout << "input has " << brts.size() << " speciations\n";
    Rcout << "number of parameters: " << model->nparams() << '\n';
    Rcout << "thread safe: " << (model->is_threadsafe() ? "true" : "false") << '\n';
    Rcout << "model description:\n" << model->description();
    Rcout << "\nrunning mcE_step with N = " << sample_size;
    Rcout << "\nInitial parameters: { ";
    for (auto p : init_pars) Rcout << p << " ";
    Rcout << "}\n";
  }
  auto mcem = emphasis::mcem(sample_size,
                             init_pars,
                             brts,
                             model.get(),
                             soc,
                             max_misssing,
                             max_lambda,
                             lower_bound,
                             upper_bound,
                             xtol,
                             emphasis::default_nlopt_algo,
                             num_threads);

  if(verbose) {
    Rcout << "valid trees: " << mcem.e.trees.size();
    if (mcem.e.trees.size() != static_cast<size_t>(sample_size)) {
      Rcout << "\n  reasons for rejection:";
      Rcout << "\n  * too big: " << mcem.e.rejected_overruns;
      Rcout << "\n  * insane lambda: " << mcem.e.rejected_lambda;
      Rcout << "\n  * zero weights: " << mcem.e.rejected_zero_weights;
    }
  }
  if (mcem.e.trees.empty()) throw std::runtime_error("no trees, no optimization");
  
  if(verbose) {
    Rcout << "\naverage branches for the valid trees: ";
    size_t b = 0; for (const auto& tree : mcem.e.trees) { b += tree.size(); }
    Rcout << static_cast<double>(b) / mcem.e.trees.size();
    Rcout << "\n\nnlopt result: " << mcem.m.opt << " (" << ((mcem.m.opt < 0) ? "fail" : "ok") << ")\n";
    Rcout << "augmentation runtime: " << mcem.e.elapsed << " ms\n";
    Rcout << "optimization runtime: " << mcem.m.elapsed << " ms\n";
    Rcout << "estimated parameters: { ";
    for (auto p : mcem.m.estimates) Rcout << p << " ";
    Rcout << "}\n";
  }
  
  List ret;
  if (copy_trees) {
    List trees;
    for (const emphasis::tree_t& tree : mcem.e.trees) {
      trees.push_back(unpack(tree));
    }
    ret["trees"] = trees;
  }
  ret["estimates"] = NumericVector(mcem.m.estimates.begin(), mcem.m.estimates.end());
  ret["nlopt"] = mcem.m.opt;
  ret["fhat"]  = mcem.e.fhat;
  ret["time"]  = mcem.e.elapsed + mcem.m.elapsed;
  return ret;
}
