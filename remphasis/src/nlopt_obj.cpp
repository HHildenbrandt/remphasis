#include <limits>
#ifdef EMP_BUILD_STANDALONE_CPP
# include <nlopt.h>
#else
# include <nloptrAPI.h>
#endif
#include "emphasis.hpp"
#include "nlopt_obj.hpp"


namespace emphasis {

  namespace {

    class linked_nlopt_1 : public NLopt_1
    {
    public:
      linked_nlopt_1()
        : lower_(-std::numeric_limits<double>::max()),
          upper_(+std::numeric_limits<double>::max())
      {
        nlopt_ = nlopt_create(NLOPT_LN_SBPLX, 1);
        if (nullptr == nlopt_) {
          throw emphasis_error("nlopt_create failed");
        }
      }

      ~linked_nlopt_1() override
      {
        if (nlopt_) nlopt_destroy(nlopt_);
      }

      void set_xtol_rel(double val) override
      {
        if (NLOPT_SUCCESS > (result_ = nlopt_set_xtol_rel(nlopt_, val))) {
          throw emphasis::emphasis_error("nlopt_set_xtol_failed");
        }
      }

      void set_lower_bounds(double val) override
      {
        lower_ = val;
        if (NLOPT_SUCCESS > (result_ = nlopt_set_lower_bounds(nlopt_, &lower_))) {
          throw emphasis_error("nlop_set_lower_bounds failed");
        }
      }

      void set_upper_bounds(double val) override
      {
        upper_ = val;
        if (NLOPT_SUCCESS > (result_ = nlopt_set_upper_bounds(nlopt_, &upper_))) {
          throw emphasis_error("nlopt_set_upper_bounds failed");
        }
      }

      void set_min_objective(nlopt_func dx, void* fdata) override
      {
        if (NLOPT_SUCCESS > (result_ = nlopt_set_min_objective(nlopt_, dx, fdata))) {
          throw emphasis_error("nlopt_set_min_objective failed");
        }
      }

      void set_max_objective(nlopt_func dx, void* fdata) override
      {
        if (NLOPT_SUCCESS > (result_ = nlopt_set_max_objective(nlopt_, dx, fdata))) {
          throw emphasis_error("nlopt_set_max_objective failed");
        }
      }

      double optimize(double& x) override
      {
        double fmin = 0.0;
        if (NLOPT_SUCCESS > (result_ = nlopt_optimize(nlopt_, &x, &fmin))) {
          throw emphasis_error("optimize failed");
        }
        return fmin;
      }

      int result() override { return static_cast<int>(result_); }

    private:
      nlopt_opt nlopt_ = nullptr;
      nlopt_result result_ = nlopt_result::NLOPT_FAILURE;
      double lower_, upper_;
    };


    class linked_nlopt : public NLopt
    {
    public:
      explicit linked_nlopt(size_t nparams)
        : lower_(nparams, -std::numeric_limits<double>::max()),
        upper_(nparams, +std::numeric_limits<double>::max())
      {
        nlopt_ = nlopt_create(NLOPT_LN_SBPLX, static_cast<unsigned>(nparams));
        if (nullptr == nlopt_) {
          throw emphasis_error("nlopt_create failed");
        }
      }

      ~linked_nlopt() override
      {
        if (nlopt_) nlopt_destroy(nlopt_);
      }

      void set_xtol_rel(double val) override
      {
        if (NLOPT_SUCCESS > (result_ = nlopt_set_xtol_rel(nlopt_, val))) {
          throw emphasis::emphasis_error("nlopt_set_xtol_failed");
        }
      }

      void set_lower_bounds(const std::vector<double>& val) override
      {
        lower_ = val;
        if (NLOPT_SUCCESS > (result_ = nlopt_set_lower_bounds(nlopt_, lower_.data()))) {
          throw emphasis_error("nlop_set_lower_bounds failed");
        }
      }

      void set_upper_bounds(const std::vector<double>& val) override
      {
        upper_ = val;
        if (NLOPT_SUCCESS > (result_ = nlopt_set_upper_bounds(nlopt_, upper_.data()))) {
          throw emphasis_error("nlopt_set_upper_bounds failed");
        }
      }

      void set_min_objective(nlopt_func dx, void* fdata) override
      {
        if (NLOPT_SUCCESS > (result_ = nlopt_set_min_objective(nlopt_, dx, fdata))) {
          throw emphasis_error("nlopt_set_min_objective failed");
        }
      }

      void set_max_objective(nlopt_func dx, void* fdata) override
      {
        if (NLOPT_SUCCESS > (result_ = nlopt_set_max_objective(nlopt_, dx, fdata))) {
          throw emphasis_error("nlopt_set_max_objective failed");
        }
      }

      double optimize(std::vector<double>& x) override
      {
        double fmin = 0.0;
        if (NLOPT_SUCCESS > (result_ = nlopt_optimize(nlopt_, x.data(), &fmin))) {
          throw emphasis_error("nlopt_optimize failed");
        }
        return fmin;
      }

      int result() override { return static_cast<int>(result_); }

    private:
      nlopt_opt nlopt_ = nullptr;
      nlopt_result result_ = nlopt_result::NLOPT_FAILURE;
      std::vector<double> lower_, upper_;
    };

  }



  std::unique_ptr<NLopt_1> create_nlopt_1()
  {
    return std::unique_ptr<NLopt_1>(new linked_nlopt_1());
  }


  std::unique_ptr<NLopt> create_nlopt(size_t nparams)
  {
    return std::unique_ptr<NLopt>(new linked_nlopt(nparams));
  }

}
