#include <memory>
#include <limits>
#include <algorithm>
#include "plugin.h"
#include "model_helpers.hpp"


using namespace emphasis;
using namespace detail;


namespace {

  constexpr int Npars = 4;
  constexpr const char* Description = R"descr(rpd5c model, dynamic link library.)descr";

  using reng_t = std::mt19937_64;   // we need doubles
  static thread_local reng_t reng_ = emphasis::detail::make_random_engine_low_entropy<reng_t>();

}


#ifdef __cplusplus
extern "C" {
#endif


  inline double speciation_rate(const double* pars, const emp_node_t& node)
  {
    const double lambda = pars[1] + pars[2] * node.n + pars[3] * node.pd / node.n;
    return std::max(0.0, lambda);
  }


  EMP_EXTERN(const char*) emp_description() { return Description; }
  EMP_EXTERN(bool) emp_is_threadsafe() { return true; }
  EMP_EXTERN(bool) emp_has_discrete_speciation_rate() { return true; }
  EMP_EXTERN(int) emp_nparams() { return Npars; }


  EMP_EXTERN(double) emp_nh_rate(void** state, double t, const double* pars, unsigned n, const emp_node_t* tree)
  {
    auto it = std::lower_bound(tree, tree + n, t, detail::node_less{});
    it = std::min(it, tree + n - 1);
    const double lambda = std::max(0.0, pars[1] + pars[2] * it->n + pars[3] * it->pd / it->n);
    return lambda * it->n * (1.0 - std::exp(-pars[0] * (tree[n - 1].brts - t)));
  }


  EMP_EXTERN(double) emp_extinction_time(void**, double t_speciation, const double* pars, unsigned n, const emp_node_t* tree)
  {
    return t_speciation + detail::trunc_exp(0, tree[n - 1].brts - t_speciation, pars[0], reng_);
  }


  EMP_EXTERN(double) emp_intensity(void** state, const double* pars, unsigned n, const emp_node_t* tree)
  {
    detail::mu_integral muint(pars[0], tree[n - 1].brts);
    double inte = 0;
    double prev_brts = 0;
    for (unsigned i = 0; i < n; ++i) {
      const auto& node = tree[i];
      const double lambda = speciation_rate(pars, node);
      inte += node.n * lambda * muint(prev_brts, node.brts);
      prev_brts = node.brts;
    }
    return inte;
  }


  EMP_EXTERN(double) emp_loglik(void** state, const double* pars, unsigned n, const emp_node_t* tree)
  {
    detail::log_sum log_sr{};
    int cex = 0;
    double inte = 0.0;
    double prev_brts = 0.0;
    for (unsigned i = 0; i < n; ++i) {
      const auto& node = tree[i];
      const double sr = speciation_rate(pars, node);
      if (detail::is_extinction(node)) {
        ++cex;
      }
      else {
        log_sr += sr;
      }
      // intensity.numerical2
      inte += (node.brts - prev_brts) * node.n * (sr + pars[0]);
      prev_brts = node.brts;
    }
    const double loglik = std::log(pars[0]) * cex + log_sr.result() - inte;
    return loglik;
  }


  EMP_EXTERN(double) emp_sampling_prob(void** state, const double* pars, unsigned n, const emp_node_t* tree)
  {
    double logg = -emp_intensity(nullptr, pars, n, tree);
    double tips = tree[0].n;
    for (unsigned i = 0; i < n; ++i) {
      const auto& node = tree[i];
      tips += detail::is_tip(node);
      if (detail::is_missing(node)) {
        const double lambda = speciation_rate(pars, node);
        const double lifespann = node.t_ext - node.brts;
        logg += std::log(node.n * pars[0] * lambda) - pars[0] * lifespann - std::log(node.n + tips);
      }
    }
    return logg;

    // old stuff
    //double logg = -emp_intensity(state, pars, n, tree);
    //tree_t mtree;
    //double nb = 0.0;
    //double no = tree[0].n;
    //double ne = 0.0;
    //struct nnn_t { 
    //  nnn_t() = default;
    //  nnn_t(double NB, double NO, double NE) : Nb(NB), No(NO), Ne(NE) {}
    //  double Nb = 0; double No = 0; double Ne = 0; 
    //};
    //
    //std::vector<nnn_t> nnn;
    //for (unsigned i = 0; i < n; ++i) {
    //  const auto& node = tree[i];
    //  const bool extinction = detail::is_extinction(node);
    //  const bool tip = detail::is_tip(node);
    //  const bool missing = !(extinction || tip);
    //  if (missing) {
    //    nnn.push_back(nnn_t(node.n, no, nb - ne));
    //    mtree.push_back(node);
    //  }
    //  if (extinction) ++ne;
    //  if (tip) ++no;
    //  if (missing) ++nb;
    //}
    //void* tmp = nullptr;
    //for (size_t i = 0; i < mtree.size(); ++i) {
    //  const auto lambda = speciation_rate(pars, mtree[i]);
    //  const auto lifespan = mtree[i].t_ext - mtree[i].brts;
    //  logg += std::log(nnn[i].Nb * pars[0] * lambda) - pars[0] * lifespan - std::log(2.0 * nnn[i].No + nnn[i].Ne);
    //}
    //return logg;
  }


  EMP_EXTERN(void) emp_lower_bound(double* pars)
  {
    pars[0] = 10e-9; pars[1] = 10e-9; pars[2] = pars[3] = -100.0;
  }


  EMP_EXTERN(void) emp_upper_bound(double* pars)
  {
    pars[0] = pars[1] = pars[2] = pars[3] = 100.0;
  }


#ifdef __cplusplus
}
#endif
