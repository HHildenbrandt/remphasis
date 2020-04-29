#include <memory>
#include <limits>
#include <algorithm>
#include "plugin.h"
#include "model_helpers.hpp"


using namespace emphasis;
using namespace detail;


namespace {

  constexpr int Npars = 4;
  constexpr const char* Description = R"descr(rpd5 model, dynamic link library.)descr";

  using reng_t = std::mt19937_64;   // we need doubles
  static thread_local reng_t reng_ = emphasis::detail::make_random_engine_low_entropy<reng_t>();

}


#ifdef __cplusplus
extern "C" {
#endif


  EMP_EXTERN(const char*) emp_description() { return Description; }
  EMP_EXTERN(bool) emp_is_threadsafe() { return true; }
  EMP_EXTERN(int) emp_nparams() { return Npars; }


  EMP_EXTERN(double) emp_speciation_rate(void** state, double t, const double* pars, unsigned n, const emp_node_t* tree)
  {
    auto it = std::lower_bound(tree, tree + n, t, detail::node_less{});
    it = std::max(it, tree + n - 1);
    const double lambda = pars[1] + pars[2] * it->n + pars[3] * it->pd / it->n;
    return std::max(0.0, lambda);
  }


  EMP_EXTERN(double) emp_speciation_rate_sum(void** state, double t, const double* pars, unsigned n, const emp_node_t* tree)
  {
    auto it = std::lower_bound(tree, tree + n, t, detail::node_less{});
    it = std::max(it, tree + n - 1);
    const double lambda = std::max(0.0, pars[1] + pars[2] * it->n + pars[3] * it->pd / it->n);
    return lambda * it->n * (1.0 - std::exp(-pars[0] * (tree[n - 1].brts - t)));
  }


  EMP_EXTERN(double) emp_extinction_time(void**, double t_speciation, const double* pars, unsigned n, const emp_node_t* tree)
  {
    return t_speciation + detail::trunc_exp(0, tree[n - 1].brts - t_speciation, pars[0], reng_);
  }


  inline double ind_rpd5(double x, double c1, double c2, double c3, double c4) 
  {
    return 0.5 * (c2 * x * x) + c1 * x - (c3 * std::exp(c4 * x) * (c2 * (c4 * x - 1) + c1 * c4)) / (c4 * c4);
  }


  EMP_EXTERN(double) emp_intensity(void** state, const double* pars, unsigned n, const emp_node_t* tree)
  {
    const double max_brts = tree[n - 1].brts;
    const double c2 = pars[3];
    const double c3 = std::exp(-pars[0] * max_brts);
    const double c4 = pars[0];
    double sum_inte = 0;
    double prev_brts = 0;
    double prev_pd = 0;
    for (unsigned i = 0; i < n; ++i) {
      const auto& node = tree[i];
      const double c1 = pars[1] + pars[2] * node.n + pars[3] * ((prev_pd - node.n * prev_brts) / node.n);
      double tmp0 = -c1 / c2;
      double tmp1 = node.brts;
      if ((pars[1] + pars[2] * node.n + pars[4] * prev_pd / node.n) < 0.0) std::swap(tmp0, tmp1);
      sum_inte += (ind_rpd5(tmp0, c1, c2, c3, c4) - ind_rpd5(tmp1, c1, c2, c3, c4)) * node.n;
      prev_brts = node.brts;
    }
    return sum_inte;
  }


  EMP_EXTERN(double) emp_loglik(void** state, const double* pars, unsigned n, const emp_node_t* tree)
  {
    double sum_inte = 0;
    auto sum_rho = detail::log_sum{};
    auto prev_sum_rho = detail::log_sum{};
    double prev_brts = 0;
    double prev_pd = 0;
    const bool z3 = pars[3] == 0;
    for (unsigned i = 0; i < n; ++i) {
      const auto& node = tree[i];
      const double wt = node.brts - prev_brts;
      const double pd2 = prev_pd + node.n * wt;
      const double lambda = pars[1] + pars[2] * node.n + pars[3] * pd2 / node.n;
      const auto to = detail::is_extinction(node) ? 0.0 : 1.0;
      sum_rho += lambda * to + pars[0] * (1.0 - to);
      double inte = node.n * (pars[0] * wt);
      if (!z3) {
        const double c1 = pars[1] + pars[2] * node.n + (pars[3] / node.n) * (prev_pd - prev_brts * node.n);
        const double r = -c1 / pars[0];
        double b0 = prev_brts;
        double b1 = node.brts;
        if ((prev_brts > r) && (node.brts < r)) {
          if (pars[0] > 0) b0 = r; else b1 = r;
        }
        inte += node.n * (c1 * (b1 - b0) + 0.5 * pars[3] * (b1 * b1 - b0 * b0));
      }
      sum_inte += inte;
      prev_sum_rho = sum_rho;
      prev_pd = node.pd;
    }
    double log_lik = prev_sum_rho.result() - sum_inte;
    return log_lik;
  }


  EMP_EXTERN(double) emp_sampling_prob(void** state, const double* pars, unsigned n, const emp_node_t* tree)
  {
    double logg = -emp_intensity(state, pars, n, tree);
    tree_t mtree;
    double nb = 0.0;
    double no = tree[0].n;
    double ne = 0.0;
    struct nnn_t { 
      nnn_t() = default;
      nnn_t(double NB, double NO, double NE) : Nb(NB), No(NO), Ne(NE) {}
      double Nb = 0; double No = 0; double Ne = 0; 
    };
    
    std::vector<nnn_t> nnn;
    for (unsigned i = 0; i < n; ++i) {
      const auto& node = tree[i];
      const bool extinction = detail::is_extinction(node);
      const bool tip = detail::is_tip(node);
      const bool missing = !(extinction || tip);
      if (missing) {
        nnn.push_back(nnn_t(node.n, no, nb - ne));
        mtree.push_back(node);
      }
      if (extinction) ++ne;
      if (tip) ++no;
      if (missing) ++nb;
    }
    void* tmp = nullptr;
    for (size_t i = 0; i < mtree.size(); ++i) {
      const auto lambda = emp_speciation_rate(&tmp, mtree[i].brts, pars, static_cast<unsigned>(mtree.size()), reinterpret_cast<const emp_node_t*>(mtree.data()));
      const auto lifespan = mtree[i].t_ext - mtree[i].brts;
      logg += std::log(nnn[i].Nb * pars[0] * lambda) - pars[0] * lifespan - std::log(2.0 * nnn[i].No + nnn[i].Ne);
    }
    return logg;
  }


  EMP_EXTERN(void) emp_lower_bound(double* pars)
  {
    pars[0] = 10e-9; pars[1] = 10e-9; pars[2] = -huge; pars[3] = -huge;
  }


  EMP_EXTERN(void) emp_upper_bound(double* pars)
  {
    pars[0] = pars[1] = pars[2] = pars[3] = huge;
  }


#ifdef __cplusplus
}
#endif
