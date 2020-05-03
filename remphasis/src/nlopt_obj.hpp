#ifndef EMPHASIS_NLOPT_OBJ_HPP_INCLUDED
#define EMPHASIS_NLOPT_OBJ_HPP_INCLUDED

#include <memory>
#include <vector>


namespace emphasis {


  class NLopt_1
  {
  public:
    typedef double (*nlopt_func_fwd) (unsigned, const double*, double*, void*);

  public:
    virtual ~NLopt_1() {}
    virtual void set_xtol_rel(double) = 0;
    virtual void set_lower_bounds(double) = 0;
    virtual void set_upper_bounds(double) = 0;
    virtual void set_min_objective(nlopt_func_fwd, void*) = 0;
    virtual void set_max_objective(nlopt_func_fwd, void*) = 0;
    virtual double optimize(double&) = 0;
    virtual int result() = 0;     // static_cast<int>(nlopt_result)
  };


  class NLopt
  {
  public:
    typedef double (*nlopt_func_fwd) (unsigned, const double*, double*, void*);

  public:
    virtual ~NLopt() {}
    virtual void set_xtol_rel(double) = 0;
    virtual void set_lower_bounds(const std::vector<double>&) = 0;
    virtual void set_upper_bounds(const std::vector<double>&) = 0;
    virtual void set_min_objective(nlopt_func_fwd, void*) = 0;
    virtual void set_max_objective(nlopt_func_fwd, void*) = 0;
    virtual double optimize(std::vector<double>&) = 0;
    virtual int result() = 0;     // static_cast<int>(nlopt_result)
  };


  std::unique_ptr<NLopt_1> create_nlopt_1();
  std::unique_ptr<NLopt> create_nlopt(size_t nparams);

}

#endif
