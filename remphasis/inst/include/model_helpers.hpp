#ifndef EMPHASIS_MODEL_HELPRES_HPP_INCLUDED
#define EMPHASIS_MODEL_HELPRES_HPP_INCLUDED

#include <limits>
#include <cmath>
#include <random>
#include <array>
#include <chrono>
#include <thread>
#include "plugin.hpp"


#ifndef EMPHASIS_LOGSUM_LOWER_TRESHOLD
#define EMPHASIS_LOGSUM_LOWER_TRESHOLD 10e-40
#endif

#ifndef EMPHASIS_LOGSUM_UPPER_TRESHOLD
#define EMPHASIS_LOGSUM_UPPER_TRESHOLD 10e+40
#endif


namespace emphasis {

  namespace detail {

    static constexpr double huge = std::numeric_limits<double>::max();


    // returns low-entropy 512 bit array for seed sequence
    // based on std::chrono::high_resolution_clock.
    // ripped from rndutils
    inline auto make_low_entropy_seed_array() noexcept->std::array<uint64_t, 8>
    {
      // the classic: time, advertised with nano-second resolution.
      const auto e1 = static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
      // different between invocations from different threads within one app: thread-id
      const auto tid = std::this_thread::get_id();
      const uint64_t e2{ std::hash<typename std::remove_const<decltype(tid)>::type>()(tid) };
      return std::array<uint64_t, 8>{ {
          e1, e2,
          0x000000003c10b019, 0x2bf820b4dd7c1a8a,
          0x9901cf90a40883da, 0x5a3686b2e1de6e51,
          0x000000cc0494d228, 0x000000cc04b66740
      }};
    }


    // random number generator from low-entropy seed sequence
    // ripped from rndutils
    template <typename URNG>
    inline auto make_random_engine() -> URNG
    {
      auto seed_array = make_low_entropy_seed_array();
      std::seed_seq sseq(seed_array.cbegin(), seed_array.cend());
      return URNG(sseq);
    }


    inline bool is_extinction(const node_t& node) { return node.t_ext == t_ext_extinct; }
    inline bool is_tip(const node_t& node) { return node.t_ext == t_ext_tip; }
    inline bool is_missing(const node_t& node) { return !(is_extinction(node) || is_tip(node)); }


    struct node_less
    {
      bool operator()(const node_t& a, const node_t& b) const noexcept { return a.brts < b.brts; };
      bool operator()(const node_t& a, double val) const noexcept { return a.brts < val; };
      bool operator()(double val, const node_t& a) const noexcept { return val < a.brts; };
    };


    class log_sum
    {
    public:
      double result() const 
      { 
        double r = std::log(prod_) + sum_;
        if (!std::isfinite(r)) {
          const double s = std::signbit(sum_) ? -1.0 : 1.0;
          r = s * std::numeric_limits<double>::infinity();
        }
        return r; 
      }

      void operator+=(double val)
      {
        if ((prod_ > EMPHASIS_LOGSUM_LOWER_TRESHOLD) && (prod_ < EMPHASIS_LOGSUM_UPPER_TRESHOLD)) {
          prod_ *= val;
        }
        else {
          sum_ += std::log(prod_) + std::log(val);
          prod_ = 1;
        }
      }

    private:
      double prod_ = 1;
      double sum_ = 0;
    };


    template <typename RENG>
    inline double trunc_exp(double lower, double upper, double rate, RENG& reng)
    {
      std::exponential_distribution<double> exp_dist(rate);
      double result = exp_dist(reng);
      while (result < lower || result > upper) {
        result = exp_dist(reng);
      }
      return result;
    }


    // Int_0_t1 (1-exp(-mu*(tm-t)))
    class mu_integral
    {
    public:
      mu_integral(double mu, double tm)
      : mu_(mu),
        s_(1.0 / (mu * std::exp(mu * tm)))
      {}

      double operator()(double t0, double t1)
      {
        const double expt0 = (pt1_ == t0) ? expt1_ : std::exp(mu_ * t0);
        expt1_ = std::exp(mu_ * t1);
        return (t1 - t0) - s_ * (expt1_ - expt0);
      }

    private:
      double pt1_ = -1.0;
      double expt1_ = 0;
      const double mu_ = 0;
      const double s_ = 0;
    };

  }

}

#endif
