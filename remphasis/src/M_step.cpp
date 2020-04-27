#ifdef _OPENMP
#   include <omp.h>
#else
#   define omp_set_num_threads(x) 
#endif
#include <thread>
#include <chrono>
#include <numeric>
#include <cmath>  
#include <algorithm>
#include <functional>
#ifdef EMP_BUILD_STANDALONE_CPP
# include <nlopt.h>
#else
# include <nloptrAPI.h>
#endif
#include "plugin.hpp"
#include "emphasis.hpp"


namespace emphasis {

  const int default_nlopt_algo = NLOPT_LN_SBPLX;


  M_step_t::M_step_t(M_step_t&& rhs) noexcept
  {
    this->operator=(std::move(rhs));
  }

  
  M_step_t& M_step_t::operator=(M_step_t&& rhs) noexcept
  {
    estimates = rhs.estimates;
    opt = rhs.opt;
    minf = rhs.minf;
    nlopt = rhs.nlopt; rhs.nlopt = nullptr;
    elapsed = rhs.elapsed;
    return *this;
  }


  M_step_t::~M_step_t()
  {
    if (nlopt) nlopt_destroy(nlopt);
  }


  namespace {

    struct nlopt_data
    {
      const std::vector<tree_t>& trees;
      const Model* model;
      const std::vector<double>& w;
    };


    double objective(unsigned int n, const double* x, double*, void* func_data)
    {
      auto psd = reinterpret_cast<const nlopt_data*>(func_data);
      double Q = 0.0;
      param_t pars(x, x + n);
#pragma omp parallel for if(psd->model->is_threadsafe()) schedule(dynamic) reduction(+:Q)
      for (int i = 0; i < static_cast<int>(psd->trees.size()); ++i) {
        const double loglik = psd->model->loglik(pars, psd->trees[i]);
        Q += loglik * psd->w[i];
      }
      return -Q;
    }

  }


  M_step_t M_step(const param_t& pars,
                  const std::vector<tree_t>& trees,          // augmented trees
                  const std::vector<double>& weights,
                  class Model* model,
                  const param_t& lower_bound, // overrides model.lower_bound
                  const param_t& upper_bound, // overrides model.upper.bound
                  double xtol,
                  int algo,
                  int num_threads)
  {
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    num_threads = std::max(1, std::min(num_threads, static_cast<int>(std::thread::hardware_concurrency())));
    omp_set_num_threads(num_threads);
    auto T0 = std::chrono::high_resolution_clock::now();
    nlopt_data sd{ trees, model, weights };
    nlopt_opt opt = nlopt_create(static_cast<nlopt_algorithm>(algo), static_cast<unsigned int>(model->nparams()));
    auto M = M_step_t{};
    if (opt != nullptr) {
      M.nlopt = opt;
      nlopt_set_xtol_rel(opt, xtol);
      auto lower = lower_bound.empty() ? model->lower_bound() : lower_bound;
      if (!lower.empty()) nlopt_set_lower_bounds(opt, lower.data());
      auto upper = upper_bound.empty() ? model->upper_bound() : upper_bound;
      if (!upper.empty()) nlopt_set_upper_bounds(opt, upper.data());
      nlopt_set_min_objective(opt, objective, &sd);
    }
    M.estimates = pars;
    M.opt = nlopt_optimize(opt, M.estimates.data(), &M.minf);
    auto T1 = std::chrono::high_resolution_clock::now();
    M.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
    return M;
  }

}

