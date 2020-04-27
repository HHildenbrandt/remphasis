#ifndef EMPHASIS_EMPHASIS_HPP_INCLUDED
#define EMPHASIS_EMPHASIS_HPP_INCLUDED

#include <stdexcept>
#include <limits>
#include <memory>
#include <vector>
#include "plugin.h"   // node_t 


namespace emphasis {

  extern const int default_nlopt_algo;  // resolves to NLOPT_SBPLX
  static constexpr int default_max_missing_branches = 10000;
  static constexpr double default_max_aug_lambda = 500.0;
  static constexpr double huge = std::numeric_limits<double>::max();


  class emphasis_error : public std::runtime_error
  {
  public:
    explicit emphasis_error(const std::string& what) : std::runtime_error(what) {}
    explicit emphasis_error(const char* what) : std::runtime_error(what) {}
  };


  using param_t = std::vector<double>;    // unspecific parameter set
  using brts_t = std::vector<double>;     // input tree

  // row in tree-matrix
  // for c++ this is wasteful but all double
  // makes it easier to represent a tree as
  // R-matrix
  struct node_t {
    node_t() = default;
    node_t(double B, double N, double T) : brts(B), n(N), t_ext(T) {}
    double brts = 0;
    double n = 0;         // n[i] = number of species in [time_i-1, time_i)
    double t_ext = 0;     // t_ext_tip for tips, t_ext_extinct extinction nodes
  };
  constexpr double t_ext_tip = 10e10;     // t_ext for present nodes
  constexpr double t_ext_extinct = 0.0;   // t_ext for extinction nodes


  // tree, sorted by node_t::time
  using tree_t = std::vector<node_t>;


  // results from e
  struct E_step_t
  {
    E_step_t(E_step_t&&) = default;
    E_step_t& operator=(E_step_t&&) = default;
    E_step_t(const E_step_t&) = delete;
    E_step_t& operator=(const E_step_t&) = delete;
    E_step_t() {};
    ~E_step_t() {};

    std::vector<tree_t> trees;          // augmented trees
    std::vector<double> weights;
    int rejected_overruns = 0;          // # trees rejected because overrun of missing branches
    int rejected_lambda = 0;            // # trees rejected because of lambda overrun
    int rejected_zero_weights = 0;      // # trees rejected because of zero-weight
    double elapsed = 0;                 // elapsed runtime [ms]
  };


  E_step_t E_step(int N,      // number of augmented trees
                  const param_t& pars,
                  const brts_t& brts,
                  class Model* model,
                  int soc = 2,
                  int max_missing = default_max_missing_branches,
                  double max_lambda = default_max_aug_lambda,
                  int num_threads = 0);


  // results from m
  struct M_step_t
  {
    M_step_t(M_step_t&&) noexcept;
    M_step_t& operator=(M_step_t&&) noexcept;
    M_step_t(const M_step_t&) = delete;
    M_step_t& operator=(const M_step_t&) = delete;
    M_step_t() {};
    ~M_step_t();

    param_t estimates;
    int opt = -1;                       // nlopt result
    double minf = 0.0;
    struct nlopt_opt_s* nlopt = nullptr;
    double elapsed = 0.0;               // elapsed runtime [ms]
  };


  M_step_t M_step(const param_t& pars,
                  const std::vector<tree_t>& trees,          // augmented trees
                  const std::vector<double>& weights,
                  class Model* model,
                  const param_t& lower_bound = {}, // overrides model.lower_bound
                  const param_t& upper_bound = {}, // overrides model.upper.bound
                  double xtol = 0,
                  int algo = default_nlopt_algo,
                  int num_threads = 0);


  // results from mcem
  struct mcem_t
  {
    mcem_t(mcem_t&&) = default;
    mcem_t& operator=(mcem_t&&) = default;
    mcem_t(const mcem_t&) = delete;
    mcem_t& operator=(const mcem_t&) = delete;
    mcem_t() {};
    ~mcem_t() {};

    E_step_t e;
    M_step_t m;
  };


  mcem_t mcem(int N,      // number of augmented trees
              const param_t& pars,
              const brts_t& brts,
              class Model* model,
              int soc = 2,
              int max_missing = default_max_missing_branches,
              double max_lambda = default_max_aug_lambda,
              const param_t& lower_bound = {}, // overrides model.lower_bound
              const param_t& upper_bound = {}, // overrides model.upper.bound
              double xtol = 0,
              int algo = default_nlopt_algo,
              int num_threads = 0);


  std::unique_ptr<class Model> create_plugin_model(const std::string& model_dll);

}

#endif
