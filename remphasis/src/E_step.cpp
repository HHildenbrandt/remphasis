#ifdef _OPENMP
#   include <omp.h>
#else
#   define omp_set_num_threads(x) 
#endif
#include <mutex>
#include <atomic>
#include <algorithm>
#include <numeric>
#include "emphasis.hpp"
#include "augment_tree.hpp"
#include "plugin.hpp"
#include "state_guard.hpp"
#include "model_helpers.hpp"


namespace emphasis {

  namespace detail {

    // this little addition reduces the load to memory allocator massively.
    tree_t thread_local pooled_tree;


    void inplace_cumsum_of_diff(brts_t& input)
    {
      double sum = 0.0;
      for (size_t i = 1; i < input.size(); ++i) {
        input[i - 1] = sum += (input[i - 1] - input[i]);
      }
      input.back() += sum;
    }


    tree_t create_tree(brts_t brts, double soc)
    {
      inplace_cumsum_of_diff(brts);
      tree_t tree;
      for (size_t i = 0; i < brts.size(); ++i) {
        tree.push_back({ brts[i], soc + i, t_ext_tip });
      }
      std::sort(tree.begin(), tree.end(), detail::node_less{});
      return(tree);
    }


    std::vector<double> log_w(param_t pars, 
                              const std::vector<tree_t>& trees, 
                              const Model* model)
    {
      std::vector<double> w(trees.size(), 0.0);
      std::mutex mutex;
      std::exception_ptr eptr = nullptr;
#pragma omp parallel for if(model->is_threadsafe()) schedule(static, 100)
      for (int i = 0; i < static_cast<int>(trees.size()); ++i) {
        try {
          state_guard state(model);
          state.invalidate_state();
          const auto logf = model->loglik(state, pars, trees[i]);
          const auto logg = model->sampling_prob(state, pars, trees[i]);
          w[i] = logf - logg;
        }
        catch (...) {
          std::lock_guard<std::mutex> _(mutex);
          eptr = std::current_exception();
        }
      }
      if (nullptr != eptr) {
        std::rethrow_exception(eptr);
      }
      return w;
    }
    
    void calculate_pd(tree_t& tree)
    {
      double sum = 0.0;
      double brts0 = 0.0;
      double ni = tree[0].n;
      for (auto& node : tree) {
        const double brts = node.brts;
        if (detail::is_missing(node)) {
          sum += (brts - brts0) * ni++;
          brts0 = brts;
        }
        node.pd = sum;
      }
    }

  }


  E_step_t E_step(int N,                      // number of augmented trees
                  const param_t& pars,
                  const brts_t& brts,
                  Model* model,
                  int soc,
                  int max_missing,
                  double max_lambda,
                  int num_threads)
  {
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    num_threads = std::min(num_threads, static_cast<int>(std::thread::hardware_concurrency()));
    omp_set_num_threads(num_threads);
    std::mutex mutex;
    std::atomic<bool> stop{ false };    // non-handled exception
    std::exception_ptr eptr = nullptr;
    tree_t init_tree = detail::create_tree(brts, static_cast<double>(soc));
    auto E = E_step_t{};
    auto T0 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for if(model->is_threadsafe()) schedule(dynamic, 1)
    for (int i = 0; i < N; ++i) {
      try {
        if (!stop) {
          // reuse tree from pool
          auto& pooled = detail::pooled_tree;
          emphasis::augment_tree(pars, init_tree, model, max_missing, max_lambda, pooled);
          detail::calculate_pd(pooled);
          std::lock_guard<std::mutex> _(mutex);
          E.trees.emplace_back(pooled.cbegin(), pooled.cend());
        }
      }
      catch (const augmentation_overrun&) {
        std::lock_guard<std::mutex> _(mutex);
        ++E.rejected_overruns;
      }
      catch (const augmentation_lambda&) {
        std::lock_guard<std::mutex> _(mutex);
        ++E.rejected_lambda;
      }
      catch (...) {
        std::lock_guard<std::mutex> _(mutex);
        eptr = std::current_exception();
        stop = true;    // catastrophic error; don't know how to handle this...
      }
    }
    if (nullptr != eptr) {
      std::rethrow_exception(eptr);
    }
    auto log_w = detail::log_w(pars, E.trees, model);
    const double max_log_w = *std::max_element(log_w.cbegin(), log_w.cend());
    double fhat = 0.0;
    // remove 'impossible' trees
    // calculate weights and fhat on the fly
    auto it = E.trees.begin();
    for (size_t i = 0; i < log_w.size(); ++i) {
      const auto w = std::exp(log_w[i] - max_log_w);
      if (w == 0.0) {
        it = E.trees.erase(it);
        ++E.rejected_zero_weights;
      }
      else {
        E.weights.push_back(w);
        ++it;
      }
	    fhat += std::exp(log_w[i]);
    }
    if (!E.trees.empty()) {
      E.fhat = std::log(fhat / static_cast<double>(log_w.size()));
    }
    auto T1 = std::chrono::high_resolution_clock::now();
    E.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
    return E;
  }

}
