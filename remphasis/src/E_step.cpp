#include <mutex>
#include <atomic>
#include <algorithm>
#include <numeric>
#include <tbb/tbb.h>
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
      tbb::parallel_for(tbb::blocked_range<size_t>(0, trees.size()), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i < r.end(); ++i) {
          state_guard state(model);
          state.invalidate_state();
          const auto logf = model->loglik(state, pars, trees[i]);
          const auto logg = model->sampling_prob(state, pars, trees[i]);
          w[i] = logf - logg;
        }
      });
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
                  int num_threads,
                  bool cont)
  {
    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
    std::mutex mutex;
    std::atomic<bool> stop{ false };    // non-handled exception
    tree_t init_tree = detail::create_tree(brts, static_cast<double>(soc));
    auto E = E_step_t{};
    auto T0 = std::chrono::high_resolution_clock::now();
    tbb::parallel_for(tbb::blocked_range<unsigned>(0, N, 1), [&](const tbb::blocked_range<unsigned>& r) {
      for (unsigned i = r.begin(); i < r.end(); ++i) {
        try {
          if (!stop) {
            // reuse tree from pool
            auto& pooled = detail::pooled_tree;
            emphasis::augment_tree(pars, init_tree, model, max_missing, max_lambda, pooled, cont);
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
      }
    });
    const auto log_w = detail::log_w(pars, E.trees, model);
    const double max_log_w = *std::max_element(log_w.cbegin(), log_w.cend());
    double sum_w = 0.0;
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
      sum_w += w;
    }
    E.fhat = std::log(sum_w / log_w.size()) + max_log_w;
    auto T1 = std::chrono::high_resolution_clock::now();
    E.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
    return E;
  }

}
