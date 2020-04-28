//
// C++-API for plug-in diversification model
//

#ifndef EMPHASIS_PLUGIN_HPP_INCLUDED
#define EMPHASIS_PLUGIN_HPP_INCLUDED

#include <memory>
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

    // state handling, default for stateless models
    virtual void free_state(void** state) const {}
    virtual void invalidate_state(void** state) const {}

    virtual double extinction_time(void** state, double t_speciations, const param_t& pars, const tree_t& tree) const = 0;
    virtual double speciation_rate(void** state, double t, const param_t& pars, const tree_t& tree) const = 0;
    virtual double speciation_rate_sum(void** state, double t, const param_t& pars, const tree_t& tree) const = 0;
    virtual double sampling_prob(void** state, const param_t& pars, const tree_t& tree) const = 0;
    virtual double intensity(void** state, const param_t& pars, const tree_t& tree) const = 0;
    virtual double loglik(void** state, const param_t& pars, const tree_t& tree) const = 0;

    // hints for optimization step
    virtual param_t lower_bound() const = 0;
    virtual param_t upper_bound() const = 0;
  };


  class state_guard
  {
  public:
    state_guard(const state_guard&) = delete;
    state_guard& operator=(const state_guard&) = delete;

    state_guard(const Model* model) 
      : model_(model), state_(nullptr) 
    {}

    ~state_guard() { model_->free_state(&state_); }
    operator void** () noexcept { return &state_; }
    void invalidate_state() { model_->invalidate_state(&state_); }

  private:
    const Model* model_;
    void* state_;
  };
}

#endif
