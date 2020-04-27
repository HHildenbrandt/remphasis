#if defined(_WIN32)
#if !defined(WIN32_LEAN_AND_MEAN)
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
#else
#include <dlfcn.h>
#endif
#include <string>
#include "plugin.hpp"
#include "emphasis.hpp"
#include "model_helpers.hpp"


namespace dll {

#if defined(_WIN32)
  class dynlib
  {
  public:
    dynlib(const std::string& DLL)
    {
      hModule_ = LoadLibrary(DLL.c_str());
      if (NULL == hModule_) {
        throw emphasis::emphasis_error("Unable to load dynamic library");
      }
    }

    ~dynlib()
    {
      if (hModule_) {
        FreeLibrary(hModule_);
      }
    }

    template <typename FPTR>
    FPTR get_address(const char* fname, bool optional)
    {
      FPTR fptr = reinterpret_cast<FPTR>(GetProcAddress(hModule_, fname));
      if (!optional && (nullptr == fptr)) throw emphasis::emphasis_error("Can't load function address");
      return fptr;
    }

  private:
    HMODULE hModule_ = NULL;
  };

#else // _WIN32

  class dynlib
  {
  public:
    dynlib(const std::string& DLL)
    {
      hModule_ = dlopen(DLL.c_str(), RTLD_LAZY);
      if (nullptr == hModule_) {
        throw emphasis::emphasis_error("Unable to load dynamic library");
      }
    }

    ~dynlib()
    {
      if (hModule_) {
        dlclose(hModule_);
      }
    }

    template <typename FPTR>
    FPTR get_address(const char* fname, bool optional)
    {
      FPTR fptr = reinterpret_cast<FPTR>(dlsym(hModule_, fname));
      if (!optional && (nullptr == fptr)) throw emphasis::emphasis_error("Can't load function address");
      return fptr;
    }

  private:
    void* hModule_ = nullptr;
  };
#endif
}


#define emp_local_stringify(a) #a
#define emp_local_load_address(name, opt) \
name ## _ = (dynlib_.get_address<emp_ ## name ## _func>(emp_local_stringify(emp_ ## name), opt))


namespace emphasis {

  class dyn_model_t : public Model
  {
  public:
    dyn_model_t(const std::string& DLL)
    : dynlib_(DLL)
    {
      emp_local_load_address(extinction_time, false);
      emp_local_load_address(speciation_rate, false);
      emp_local_load_address(speciation_rate_sum, false);
      emp_local_load_address(sampling_prob, false);
      emp_local_load_address(intensity, false);
      emp_local_load_address(loglik, false);
      emp_local_load_address(lower_bound, false);
      emp_local_load_address(upper_bound, false);
      emp_local_load_address(nparams, false);
      emp_local_load_address(is_threadsafe, true);
      emp_local_load_address(description, false);
      emp_local_load_address(create, true);
      emp_local_load_address(destroy, true);
      if (create_) state_ = create_();
      if (nullptr != create_ && nullptr == state_) throw emphasis::emphasis_error("Can't create model");
    }

    ~dyn_model_t() override 
    {
      if (nullptr != state_) destroy_(state_);
    }

    const char* description() const override { return description_(); }
    bool is_threadsafe() const override { return is_threadsafe_(); }
    int nparams() const override { return nparams_(); };

    double extinction_time(double t_speciation, const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(extinction_time_, state_, t_speciation, pars, tree);
    }

    double speciation_rate(double t, const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(speciation_rate_, state_, t, pars, tree);
    }

    double speciation_rate_sum(double t, const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(speciation_rate_sum_, state_, t, pars, tree);
    }

    double sampling_prob(const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(sampling_prob_, state_, pars, tree);
    }

    double intensity(const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(intensity_, state_, pars, tree);
    }

    double loglik(const param_t& pars, const tree_t& tree) const override {
      return detail::wrap(loglik_, state_, pars, tree);
    }

    param_t lower_bound() const override
    {
      param_t p(nparams());
      lower_bound_(p.data());
      return p;
    }

    param_t upper_bound() const override
    {
      param_t p(nparams());
      upper_bound_(p.data());
      return p;
    }

  private:
    static void* state_;
    static emp_extinction_time_func extinction_time_;
    static emp_speciation_rate_func speciation_rate_;
    static emp_speciation_rate_sum_func speciation_rate_sum_;
    static emp_sampling_prob_func sampling_prob_;
    static emp_intensity_func intensity_;
    static emp_loglik_func loglik_;
    static emp_description_func description_;
    static emp_is_threadsafe_func is_threadsafe_;
    static emp_nparams_func nparams_;
    static emp_lower_bound_func lower_bound_;
    static emp_upper_bound_func upper_bound_;
    static emp_create_func create_;
    static emp_destroy_func destroy_;
    dll::dynlib dynlib_;
  };


  void* dyn_model_t::state_ = nullptr;
  emp_extinction_time_func dyn_model_t::extinction_time_ = nullptr;
  emp_speciation_rate_func dyn_model_t::speciation_rate_ = nullptr;
  emp_speciation_rate_sum_func dyn_model_t::speciation_rate_sum_ = nullptr;
  emp_sampling_prob_func dyn_model_t::sampling_prob_ = nullptr;
  emp_intensity_func dyn_model_t::intensity_ = nullptr;
  emp_loglik_func dyn_model_t::loglik_ = nullptr;
  emp_description_func dyn_model_t::description_ = nullptr;
  emp_is_threadsafe_func dyn_model_t::is_threadsafe_ = nullptr;
  emp_nparams_func dyn_model_t::nparams_ = nullptr;
  emp_lower_bound_func dyn_model_t::lower_bound_ = nullptr;
  emp_upper_bound_func dyn_model_t::upper_bound_ = nullptr;
  emp_create_func dyn_model_t::create_ = nullptr;
  emp_destroy_func dyn_model_t::destroy_ = nullptr;


  std::unique_ptr<emphasis::Model> create_plugin_model(const std::string& DLL)
  {
    return std::unique_ptr<emphasis::dyn_model_t>(new emphasis::dyn_model_t(DLL));
  }

}

#undef emp_local_stringify
#undef emp_local_load_address
