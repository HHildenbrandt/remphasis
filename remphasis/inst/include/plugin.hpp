//
// C++-API for plug-in diversification model
//

#ifndef EMPHASIS_PLUGIN_HPP_INCLUDED
#define EMPHASIS_PLUGIN_HPP_INCLUDED

#include <memory>
#include <vector>
#include "plugin.h"   // emp_node_t


namespace emphasis {


  using param_t = std::vector<double>;                  // unspecific parameters
  using node_t = emp_node_t;                            // brts, n, t_ext, pd
  using tree_t = std::vector<node_t>;                   // tree, sorted by note_t::brts
  constexpr double t_ext_tip = emp_t_ext_tip;           // t_ext for present nodes
  constexpr double t_ext_extinct = emp_t_ext_extinct;   // t_ext for extinction nodes


  // abstract diversification model
  class Model
  {
  public:
    Model() = default;
    virtual ~Model() = default;

    // optional 
    virtual const char* description() const { return "not set"; }   // textual description of the model
    virtual bool is_threadsafe() const { return false; }            // is this implementation thread-save?

    // optional per-tree state handling
    virtual void free_state(void** state) const {}
    virtual void invalidate_state(void** state) const {}

    virtual int nparams() const = 0;              // number of parameters

    // diversification model, augmentation
    virtual double extinction_time(void** state, double t_speciations, const param_t& pars, const tree_t& tree) const = 0;
    virtual double speciation_rate(void** state, double t, const param_t& pars, const tree_t& tree) const = 0;
    virtual double speciation_rate_sum(void** state, double t, const param_t& pars, const tree_t& tree) const = 0;
    virtual double sampling_prob(void** state, const param_t& pars, const tree_t& tree) const = 0;
    virtual double intensity(void** state, const param_t& pars, const tree_t& tree) const = 0;

    // diversification model, optimization
    virtual double loglik(void** state, const param_t& pars, const tree_t& tree) const = 0;

    // optional hints for optimization step
    virtual param_t lower_bound() const { return param_t(); }
    virtual param_t upper_bound() const { return param_t(); }
  };

}

#endif
