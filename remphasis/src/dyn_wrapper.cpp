#include <cmath>
#include "plugin.hpp"
#include "emphasis.hpp"
#include "model_helpers.hpp"
#include "dyn_wrapper.hpp"


#define emp_local_stringify(a) #a
#define emp_local_load_address(name, optional) \
name ## _ = (dynlib_.get_address<emp_ ## name ## _func>(emp_local_stringify(emp_ ## name), optional))


namespace emphasis {

  class dyn_model_t : public Model
  {
  public:
    dyn_model_t(const std::string& DLL)
    : dynlib_(DLL)
    {
      emp_local_load_address(description, true);
      emp_local_load_address(nparams, false);
      emp_local_load_address(is_threadsafe, true);
      emp_local_load_address(free_state, true);
      emp_local_load_address(invalidate_state, true);
      emp_local_load_address(extinction_time, false);
      emp_local_load_address(nh_rate, false);
      emp_local_load_address(sampling_prob, false);
      emp_local_load_address(intensity, false);
      emp_local_load_address(loglik, false);
      emp_local_load_address(lower_bound, true);
      emp_local_load_address(upper_bound, true);
    }

    ~dyn_model_t() override {}

    void free_state(void** state) const override {
      if (free_state_) free_state_(state);
    }

    void invalidate_state(void** state) const override {
      if (invalidate_state_) invalidate_state_(state);
    }

    const char* description() const override 
    { 
      return (description_) ? description_() : "not set"; 
    }
    
    bool is_threadsafe() const override 
    { 
      return (is_threadsafe_) ? is_threadsafe_() : false; 
    }
    
    int nparams() const override
    { 
      return nparams_(); 
    };

    double extinction_time(void** state, double t_speciation, const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(extinction_time_, state, t_speciation, pars, tree);
    }

    double nh_rate(void** state, double t, const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(nh_rate_, state, t, pars, tree);
    }

    double sampling_prob(void** state, const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(sampling_prob_, state, pars, tree);
    }

    double intensity(void** state, const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(intensity_, state, pars, tree);
    }

    double loglik(void** state, const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(loglik_, state, pars, tree);
    }

    param_t lower_bound() const override
    {
      if (lower_bound_) {
        param_t p(nparams());
        lower_bound_(p.data());
        return p;
      }
      return param_t();
    }

    param_t upper_bound() const override
    {
      if (upper_bound_) {
        param_t p(nparams());
        upper_bound_(p.data());
        return p;
      }
      return param_t();
    }

  private:
    static emp_free_state_func free_state_;
    static emp_invalidate_state_func invalidate_state_;
    static emp_extinction_time_func extinction_time_;
    static emp_nh_rate_func nh_rate_;
    static emp_sampling_prob_func sampling_prob_;
    static emp_intensity_func intensity_;
    static emp_loglik_func loglik_;
    static emp_description_func description_;
    static emp_is_threadsafe_func is_threadsafe_;
    static emp_nparams_func nparams_;
    static emp_lower_bound_func lower_bound_;
    static emp_upper_bound_func upper_bound_;
    dll::dynlib dynlib_;
  };

  emp_free_state_func dyn_model_t::free_state_ = nullptr;
  emp_invalidate_state_func dyn_model_t::invalidate_state_ = nullptr;
  emp_extinction_time_func dyn_model_t::extinction_time_ = nullptr;
  emp_nh_rate_func dyn_model_t::nh_rate_ = nullptr;
  emp_sampling_prob_func dyn_model_t::sampling_prob_ = nullptr;
  emp_intensity_func dyn_model_t::intensity_ = nullptr;
  emp_loglik_func dyn_model_t::loglik_ = nullptr;
  emp_description_func dyn_model_t::description_ = nullptr;
  emp_is_threadsafe_func dyn_model_t::is_threadsafe_ = nullptr;
  emp_nparams_func dyn_model_t::nparams_ = nullptr;
  emp_lower_bound_func dyn_model_t::lower_bound_ = nullptr;
  emp_upper_bound_func dyn_model_t::upper_bound_ = nullptr;


  std::unique_ptr<emphasis::Model> create_plugin_model(const std::string& DLL)
  {
    return std::unique_ptr<emphasis::dyn_model_t>(new emphasis::dyn_model_t(DLL));
  }

}

#undef emp_local_stringify
#undef emp_local_load_address
