#ifndef EMPHASIS_MODEL_HELPRES_HPP_INCLUDED
#define EMPHASIS_MODEL_HELPRES_HPP_INCLUDED

#include <cmath>
#include <random>
#include <array>
#include <chrono>
#include <thread>
#include "emphasis.hpp"


#ifndef EMPHASIS_LOGSUM_LOWER_TRESHOLD
#define EMPHASIS_LOGSUM_LOWER_TRESHOLD 10e-20
#endif

#ifndef EMPHASIS_LOGSUM_UPPER_TRESHOLD
#define EMPHASIS_LOGSUM_UPPER_TRESHOLD 10e20
#endif


namespace emphasis {

  namespace detail {


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
    // ripped from randutils
    template <typename URNG>
    inline auto make_random_engine_low_entropy() -> URNG
    {
      auto seed_array = detail::make_low_entropy_seed_array();
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
      double result() const { return std::log(prod_) + sum_; }

      void operator+=(double val)
      {
        if (true) { //(val > 0) {
          if ((prod_ < EMPHASIS_LOGSUM_LOWER_TRESHOLD) || (prod_ > EMPHASIS_LOGSUM_UPPER_TRESHOLD)) {
            sum_ += std::log(prod_) + std::log(val);
            prod_ = 1;
          }
          else {
            prod_ *= val;
          }
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


    template <typename Fun>
    inline double wrap(Fun&& fun, void** state, const param_t& pars, const tree_t& tree)
    {
      return fun(state, pars.data(), static_cast<unsigned>(tree.size()), reinterpret_cast<const emp_node_t*>(tree.data()));
    }


    template <typename Fun>
    inline double wrap(Fun&& fun, void** state, double t, const param_t& pars, const tree_t& tree)
    {
      return fun(state, t, pars.data(), static_cast<unsigned>(tree.size()), reinterpret_cast<const emp_node_t*>(tree.data()));
    }

  }

}

#endif
