//
// C++-API for plug-in diversification model
//

#ifndef EMPHASIS_PLUGIN_HPP_INCLUDED
#define EMPHASIS_PLUGIN_HPP_INCLUDED

#include "emphasis.hpp"


namespace emphasis {

  // abstract diversification model
  class Model
  {
  public:
    Model() = default;
    virtual ~Model() = default;

    virtual const char* description() const = 0;  // textual description of the model
    virtual bool is_threadsafe() const = 0;       // is this implementation thread-save?
    virtual int nparams() const = 0;              // number of parameters

    virtual double extinction_time(double t_speciations, const param_t& pars, const tree_t& tree) const = 0;
    virtual double speciation_rate(double t, const param_t& pars, const tree_t& tree) const = 0;
    virtual double speciation_rate_sum(double t, const param_t& pars, const tree_t& tree) const = 0;
    virtual double sampling_prob(const param_t& pars, const tree_t& tree) const = 0;
    virtual double intensity(const param_t& pars, const tree_t& tree) const = 0;
    virtual double loglik(const param_t& pars, const tree_t& tree) const = 0;

    // hints for optimization step
    virtual param_t lower_bound() const = 0;
    virtual param_t upper_bound() const = 0;
  };

}

#endif
