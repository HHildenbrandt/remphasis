// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "emphasis.hpp"
#include "plugin.hpp"
#include "rinit.h"
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
List rcpp_mcem(const std::vector<double>& brts,       
               const std::vector<double>& init_pars,      
               int sample_size,                       
               const std::string& plugin,             
               int soc,
               int max_misssing,               
               double max_lambda,             
               const std::vector<double>& lower_bound,  
               const std::vector<double>& upper_bound,  
               double xtol,                     
               int num_threads,
               bool copy_trees) 
{
  try {
    auto model = emphasis::create_plugin_model(plugin);
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
                               num_threads);
    if (mcem.e.trees.empty()) {
      throw std::runtime_error("no trees, no optimization");
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
  catch (const std::exception &err) {
    forward_exception_to_r(err);
  }
  catch (...) {
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  return List{};
 }

