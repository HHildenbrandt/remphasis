#ifdef _OPENMP
#   include <omp.h>
#else
#   define omp_set_num_threads(x) 
#endif
#include <thread>
#include <mutex>
#include <chrono>
#include <numeric>
#include <cmath>  
#include <algorithm>
#include <functional>
#include <nlopt.h>
#include "plugin.hpp"
#include "emphasis.hpp"
#include "nlopt_obj.hpp"


namespace emphasis {

  const int default_nlopt_algo = NLOPT_LN_SBPLX;


  namespace {

    struct nlopt_data
    {
      nlopt_data(const Model* M, const std::vector<tree_t>& Trees, const std::vector<double>& W)
        : model(M), state(Trees.size(), nullptr), trees(Trees), w(W)
      {
        for (size_t i = 0; i < trees.size(); ++i) {
          model->invalidate_state(&state[i]);
        }
      }

      ~nlopt_data()
      {
        auto empty = tree_t{};
        for (auto s : state) {
          model->free_state(&s);
        }
      }

      const Model* model;
      std::vector<void*> state;
      const std::vector<tree_t>& trees;
      const std::vector<double>& w;
    };


    double objective(unsigned int n, const double* x, double*, void* func_data)
    {
      auto psd = reinterpret_cast<nlopt_data*>(func_data);
      double Q = 0.0;
      param_t pars(x, x + n);
#pragma omp parallel for if(psd->model->is_threadsafe()) schedule(dynamic) reduction(+:Q)
      for (int i = 0; i < static_cast<int>(psd->trees.size()); ++i) {
        const double loglik = psd->model->loglik(&psd->state[i], pars, psd->trees[i]);
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
    num_threads = std::min(num_threads, static_cast<int>(std::thread::hardware_concurrency()));
    omp_set_num_threads(num_threads);
    auto T0 = std::chrono::high_resolution_clock::now();
    nlopt_data sd{ model, trees, weights };
    auto M = M_step_t{};
    auto nlopt = create_nlopt(pars.size());
    M.estimates = pars;
    nlopt->set_xtol_rel(xtol);
    auto lower = lower_bound.empty() ? model->lower_bound() : lower_bound;
    if (!lower.empty()) nlopt->set_lower_bounds(lower);
    auto upper = upper_bound.empty() ? model->upper_bound() : upper_bound;
    if (!upper.empty()) nlopt->set_upper_bounds(upper);
    nlopt->set_min_objective(objective, &sd);
    M.minf = nlopt->optimize(M.estimates);
    M.opt = nlopt->result();
    auto T1 = std::chrono::high_resolution_clock::now();
    M.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
    return M;
  }

}

