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


  E_step_t E_step(int N,               
                  int maxN,
                  const param_t& pars,
                  const brts_t& brts,
                  Model* model,
                  int soc,
                  int max_missing,
                  double max_lambda,
                  int num_threads)
  {
    if (!model->is_threadsafe()) num_threads = 1;
    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
    std::mutex mutex;
    std::atomic<bool> stop{ false };    // non-handled exception
    tree_t init_tree = detail::create_tree(brts, static_cast<double>(soc));
    auto E = E_step_t{};
    auto T0 = std::chrono::high_resolution_clock::now();
    tbb::parallel_for(tbb::blocked_range<unsigned>(0, maxN), [&](const tbb::blocked_range<unsigned>& r) {
      for (unsigned i = r.begin(); i < r.end(); ++i) {
        try {
          if (!stop) {
            // reuse tree from pool
            auto& pool_tree = detail::pooled_tree;
            emphasis::augment_tree(pars, init_tree, model, max_missing, max_lambda, pool_tree);
            double log_w = 0.0;
            {
              state_guard state(model);
              state.invalidate_state();
              const auto logf = model->loglik(state, pars, pool_tree);
              const auto logg = model->sampling_prob(state, pars, pool_tree);
              log_w = logf - logg;
            }
            if (std::isfinite(log_w) && (0.0 < std::exp(log_w))) {
              detail::calculate_pd(pool_tree);
              std::lock_guard<std::mutex> _(mutex);
              if (!stop) {
                E.trees.emplace_back(pool_tree.cbegin(), pool_tree.cend());
                E.weights.push_back(log_w);
                if (static_cast<int>(E.trees.size()) == N) {
                  stop = true;
                }
              }
            }
            else {
              std::lock_guard<std::mutex> _(mutex);
              ++E.rejected_zero_weights;
            }
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
    if (static_cast<int>(E.weights.size()) < N) {
      throw emphasis_error("maxN exceeded");
    }
    const double max_log_w = *std::max_element(E.weights.cbegin(), E.weights.cend());
    double sum_w = 0.0;
    for (size_t i = 0; i < E.weights.size(); ++i) {
      const double w = std::exp(E.weights[i] - max_log_w);
      sum_w += (E.weights[i] = w);
    }
    E.fhat = std::log(sum_w / E.weights.size()) + max_log_w;
    E.rejected = E.rejected_lambda + E.rejected_overruns + E.rejected_zero_weights;
    auto T1 = std::chrono::high_resolution_clock::now();
    E.elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
    return E;
  }

}
