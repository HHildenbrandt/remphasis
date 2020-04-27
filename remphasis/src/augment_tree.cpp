#include <cassert>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <atomic>
#include "plugin.hpp"
#include "augment_tree.hpp"
#include "model_helpers.hpp"


#ifdef EMP_AUG_DEBUG_SEED
static std::atomic<size_t> debug_seed{ EMP_AUG_DEBUG_SEED };
#endif


namespace emphasis {

  namespace {

    auto thread_local reng = detail::make_random_engine_low_entropy<std::default_random_engine>();


    double get_next_bt(const tree_t& tree, double cbt)
    {
      auto it = std::upper_bound(tree.cbegin(), tree.cend(), cbt, detail::node_less{});
      assert(it != tree.end());
      return it->brts;
    }


    // insert speciation node_t before t_spec,
    // inserts extinction node_t before t_ext
    // and tracks n
    void insert_species(double t_spec, double t_ext, tree_t& tree)
    {
      auto n_after = [](tree_t::iterator it) { 
        const auto to = detail::is_extinction(*it) ? -1.0 : 1.0;
        return it->n + to; 
      };
      tree.reserve(tree.size() + 2);   // keep iterators valid
      auto first = std::lower_bound(tree.begin(), tree.end(), t_spec, detail::node_less{});
      auto n = (first != tree.begin()) ? n_after(first - 1) : tree.front().n;
      first = tree.insert(first, node_t(t_spec, n, t_ext));
      // recalculate dirty range
      for (++first; first->brts < t_ext; ++first) {
        first->n = n_after(first - 1);
      }
      tree.insert(first, { t_ext, n_after(first - 1), t_ext_extinct });
    }


    void do_augment_tree(const param_t& pars, tree_t& tree, const Model& model, int max_missing, double max_lambda)
    {
      double cbt = 0;
      tree.reserve(5 * tree.size());    // just a guess, should cover most 'normal' cases
      int num_missing_branches = 0;
      const double b = tree.back().brts;
      double lambda2 = 0.0;
      bool dirty = true;
      while (cbt < b) {
        double next_bt = get_next_bt(tree, cbt);
        double lambda1 = (!dirty) ? lambda2 : model.speciation_rate_sum(cbt, pars, tree);
        lambda2 = model.speciation_rate_sum(next_bt, pars, tree);
        double lambda_max = std::max<double>(lambda1, lambda2);
        if (lambda_max > max_lambda) throw augmentation_lambda{};
        double u1 = std::uniform_real_distribution<>()(reng);
        double next_speciation_time = cbt - std::log(u1) / lambda_max;
        dirty = false;
        if (next_speciation_time < next_bt) {
          double u2 = std::uniform_real_distribution<>()(reng);
          double pt = model.speciation_rate_sum(next_speciation_time, pars, tree) / lambda_max;
          if (u2 < pt) {
            double extinction_time = model.extinction_time(next_speciation_time, pars, tree);
            insert_species(next_speciation_time, extinction_time, tree);
            num_missing_branches++;
            if (num_missing_branches > max_missing) {
              throw augmentation_overrun{};
            }
            dirty = true;   // tree changed
          }
          // else: tree unchanged; cbt <- next_bt; lambda1 <- lambda2
        }
        cbt = std::min(next_speciation_time, next_bt);
      }
    }

  } // namespace augment


  tree_t augment_tree(const param_t& pars, const tree_t& input_tree, Model* model, int max_missing, double max_lambda)
  {
    tree_t unpooled(input_tree);
    do_augment_tree(pars, unpooled, *model, max_missing, max_lambda);
    return unpooled;
  }


  void augment_tree(const param_t& pars, const tree_t& input_tree, Model* model, int max_missing, double max_lambda, tree_t& pooled)
  {
    pooled.resize(input_tree.size());
    std::copy(input_tree.cbegin(), input_tree.cend(), pooled.begin());
    do_augment_tree(pars, pooled, *model, max_missing, max_lambda);
  }


}
