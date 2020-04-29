#ifdef _OPENMP
#   include <omp.h>
#else
#   define omp_get_max_threads() 1
#   define omp_get_num_threads() 1
#   define omp_set_num_threads(x) 
#endif
#include <mutex>
#include <atomic>
#include <algorithm>
#include <chrono>
#include <nlopt.h>
#include "emphasis.hpp"
#include "plugin.hpp"
#include "augment_tree.hpp"
#include "model_helpers.hpp"


namespace emphasis {

  mcem_t mcem(int N,      // number of augmented trees
              const param_t& pars,
              const brts_t& brts,
              class Model* model,
              int soc,
              int max_missing,
              double max_lambda,
              const param_t& lower_bound, // overrides model.lower_bound
              const param_t& upper_bound, // overrides model.upper.bound
              double xtol,
              int algo,
              int num_threads)
  {
    auto EM = mcem_t();
    EM.e = E_step(N, pars, brts, model, soc, max_missing, max_lambda, num_threads);
    // optimize
    if (!EM.e.trees.empty()) {
      EM.m = M_step(pars, EM.e.trees, EM.e.weights, model, lower_bound, upper_bound, xtol, algo, num_threads);
    }
    return EM;
  }

}
