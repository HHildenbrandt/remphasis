#include <memory>
#include <limits>
#include <algorithm>
#include <plugin.h>
#include <model_helpers.hpp>


using namespace emphasis;
using namespace detail;


namespace {

  constexpr int Npars = 3;
  constexpr const char* Description = R"descr(DDD model, dynamic link library.)descr";

  using reng_t = std::mt19937_64;   // we need doubles
  static thread_local reng_t reng_ = emphasis::detail::make_random_engine_low_entropy<reng_t>();
}


#ifdef __cplusplus
extern "C" {
#endif

  
  EMP_EXTERN(const char*) emp_description() { return Description; }
  EMP_EXTERN(bool) emp_is_threadsafe() { return true; }
  EMP_EXTERN(int) emp_nparams() { return Npars; }


  EMP_EXTERN(double) emp_extinction_time(void*, double t_speciation, const double* pars, unsigned n, const emp_node_t* tree)
  {
    return t_speciation + detail::trunc_exp(0, tree[n - 1].brts - t_speciation, pars[0], reng_);
  }


  EMP_EXTERN(double) emp_speciation_rate(void*, double t, const double* pars, unsigned n, const emp_node_t* tree)
  {
    auto it = std::lower_bound(tree, tree + n, t, detail::node_less{});
    it = std::max(it, tree + n - 1);
    const double N = it->n;
    const double lambda = std::max<double>(0.0, pars[1] + pars[2] * N);
    return lambda;
  }


  EMP_EXTERN(double) emp_speciation_rate_sum(void*, double t, const double* pars, unsigned n, const emp_node_t* tree)
  {
    auto it = std::lower_bound(tree, tree + n, t, detail::node_less{});
    it = std::max(it, tree + n - 1);
    const double N = it->n;
    const double lambda = std::max<double>(0.0, pars[1] + pars[2] * N);
    return lambda * N * (1.0 - std::exp(-pars[0] * (tree[n - 1].brts - t)));
  }


  EMP_EXTERN(double) emp_intensity(void*, const double* pars, unsigned n, const emp_node_t* tree)
  {
    double sum_sigma = 0;
    const double max_brts = tree[n - 1].brts;
    const double exp_max_term = std::exp(-pars[0] * max_brts) / pars[0];
    double exp_brts_m1_term = 1.0;
    double prev_brts = 0;
    for (unsigned i = 0; i < n; ++i) {
      const auto& node = tree[i];
      const double lambda = std::max(0.0, pars[1] + pars[2] * node.n);
      const double wt = node.brts - prev_brts;
      const double exp_brts_term = std::exp(pars[0] * node.brts);
      const double sigma = node.n * lambda * (wt - exp_max_term * (exp_brts_term - exp_brts_m1_term));
      sum_sigma += sigma;
      exp_brts_m1_term = exp_brts_term;
      prev_brts = node.brts;
    }
    return sum_sigma;
  }


  EMP_EXTERN(double) emp_sampling_prob(void*, const double* pars, unsigned n, const emp_node_t* tree)
  {
    double logg = -emp_intensity(nullptr, pars, n, tree);
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
        nnn.emplace_back(nnn_t(node.n, no, nb - ne));
        mtree.push_back(node);
      }
      if (extinction) ++ne;
      if (tip) ++no;
      if (missing) ++nb;
    }
    for (size_t i = 0; i < mtree.size(); ++i) {
      const auto lambda = emp_speciation_rate(nullptr, mtree[i].brts, pars, static_cast<unsigned>(mtree.size()), reinterpret_cast<const emp_node_t*>(mtree.data()));
      const auto lifespan = mtree[i].t_ext - mtree[i].brts;
      logg += std::log(nnn[i].Nb * pars[0] * lambda) - pars[0] * lifespan - std::log(2.0 * nnn[i].No + nnn[i].Ne);
    }
    return logg;
  }


  EMP_EXTERN(double) emp_loglik(void*, const double* pars, unsigned n, const emp_node_t* tree)
  {
    double sum_inte = 0;
    auto sum_rho = detail::log_sum{};
    auto prev_sum_rho = detail::log_sum{};
    double prev_brts = 0;
    for (unsigned i = 0; i < n; ++i) {
      const auto& node = tree[i];
      const double wt = node.brts - prev_brts;
      const double lambda = std::max<double>(0, pars[1] + pars[2] * node.n);
      sum_inte += node.n * (pars[0] + lambda) * wt;
      prev_sum_rho = sum_rho;
      const auto to = detail::is_extinction(node) ? 0.0 : 1.0;
      sum_rho += lambda * to + pars[0] * (1.0 - to);
      prev_brts = node.brts;
    }
    double log_lik = prev_sum_rho.result() - sum_inte;
    return log_lik;
  }


  EMP_EXTERN(void) emp_lower_bound(double* pars)
  {
    pars[0] = 10e-9; pars[1] = 10e-9; pars[2] = -1.0;
  }


  EMP_EXTERN(void) emp_upper_bound(double* pars)
  {
    pars[0] = pars[1] = pars[2] = huge;
  }

   
#ifdef __cplusplus
}
#endif
