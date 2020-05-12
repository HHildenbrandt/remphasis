#include <memory>
#include <limits>
#include <algorithm>
#include <plugin.h>
#include <model_helpers.hpp>


using namespace emphasis::detail;


namespace {

  using reng_t = std::mt19937_64;   // we need doubles
  static thread_local reng_t reng_ = make_random_engine<reng_t>();

}


inline double speciation_rate(const double* pars, const emp_node_t& node)
{
  const double lambda = std::max<double>(0.0, pars[1] + pars[2] * node.n);
  return lambda;
}


EMP_EXTERN(const char*) emp_description() { return "rpd1 model, dynamic link library."; }
EMP_EXTERN(bool) emp_is_threadsafe() { return true; }
EMP_EXTERN(int) emp_nparams() { return 3; }


EMP_EXTERN(double) emp_extinction_time(void**, double t_speciation, const double* pars, unsigned n, const emp_node_t* tree)
{
  return t_speciation + trunc_exp(0, tree[n - 1].brts - t_speciation, pars[0], reng_);
}


EMP_EXTERN(double) emp_nh_rate(void**, double t, const double* pars, unsigned n, const emp_node_t* tree)
{
  auto it = std::lower_bound(tree, tree + n, t, node_less{});
  it = std::min(it, tree + n - 1);
  const double N = it->n;
  const double lambda = speciation_rate(pars, *it);
  return lambda * N * (1.0 - std::exp(-pars[0] * (tree[n - 1].brts - t)));
}


EMP_EXTERN(double) emp_sampling_prob(void**, const double* pars, unsigned n, const emp_node_t* tree)
{
  mu_integral muint(pars[0], tree[n - 1].brts);
  double inte = 0;
  double logg = 0;
  double prev_brts = 0;
  double tips = tree[0].n;
  for (unsigned i = 0; i < n; ++i) {
    const auto& node = tree[i];
    const double lambda = speciation_rate(pars, node);
    inte += node.n * lambda * muint(prev_brts, node.brts);
    tips += is_tip(node);
    if (is_missing(node)) {
      const double lifespann = node.t_ext - node.brts;
      logg += std::log(node.n * pars[0] * lambda) - pars[0] * lifespann - std::log(node.n + tips);
    }
    prev_brts = node.brts;
  }
  return logg - inte;
}


EMP_EXTERN(double) emp_loglik(void**, const double* pars, unsigned n, const emp_node_t* tree)
{
  log_sum log_sr{};
  int cex = 0;
  double inte = 0.0;
  double prev_brts = 0.0;
  for (unsigned i = 0; i < n; ++i) {
    const auto& node = tree[i];
    const double sr = speciation_rate(pars, node);
    if (is_extinction(node)) {
      ++cex;
    }
    else {
      log_sr += sr;
    }
    inte += (node.brts - prev_brts) * node.n * (sr + pars[0]);
    prev_brts = node.brts;
  }
  const double loglik = std::log(pars[0]) * cex + log_sr.result() - inte;
  return loglik;
}


EMP_EXTERN(void) emp_lower_bound(double* pars)
{
  pars[0] = 10e-9; pars[1] = 10e-9; pars[2] = -1.0;
}


EMP_EXTERN(void) emp_upper_bound(double* pars)
{
  pars[0] = pars[1] = pars[2] = huge;
}
