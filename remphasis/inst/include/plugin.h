/***********************************************************/
/* C-API for plug-in diversification model                 */
/***********************************************************/


#ifndef EMPHASIS_PLUGIN_H_INCLUDED
#define EMPHASIS_PLUGIN_H_INCLUDED

#include <stddef.h>


#if (defined(_WIN32) || defined(__WIN32__))
#  if defined(__GNUC__)
#    define EMP_STDCALL __attribute__((stdcall))
#  elif defined(_MSC_VER) || defined(_ICC) || defined(_STDCALL_SUPPORTED)
#    define EMP_CALL __cdecl
#  else
#    define EMP_CALL
#  endif
#else
#  define EMP_CALL
#endif

#if (defined(_WIN32) || defined(__WIN32__)) && !defined(__LCC__)
#  if defined(EMP_DLL_IMPORT)
#    define EMP_EXTERN(T) extern __declspec(dllimport) T EMP_CALL
#  else
#    define EMP_EXTERN(T) extern __declspec(dllexport) T EMP_CALL
#  endif
#else
#  define EMP_EXTERN(T) extern T EMP_CALL
#endif


#ifdef __cplusplus

// exceptions must not cross DLL boundaries
# define emp_try try {
# define emp_catch(ret) catch (...) { return ret; }}
extern "C" {
#else
/*
   we didn't get into the pain to support longjump().
   returning NaN, Inf or nullptr as error-indicator is supported thought.
*/
# define emp_try
# define emp_catch(ret)
#endif /* __cplusplus */


  /* tree node */
  struct emp_node_t
  {
    double brts;
    double n;         /* n[i] = number of species in [time_i-1, time_i) */
    double t_ext;     /* +10e10 for present-day species. -1 for extinction nodes */
  };

  typedef const char* (*emp_description_func)();
  typedef bool (*emp_is_threadsafe_func)();
  typedef int (*emp_nparams_func)();
  typedef void (*emp_lower_bound_func)(double*);
  typedef void (*emp_upper_bound_func)(double*);

  typedef void (*emp_free_state_func)(void**);
  typedef void (*emp_invalidate_state_func)(void**);

  typedef double (*emp_extinction_time_func)(void**, double, const double*, unsigned, const emp_node_t*);
  typedef double (*emp_speciation_rate_func)(void**, double, const double*, unsigned, const emp_node_t*);
  typedef double (*emp_speciation_rate_sum_func)(void**, double, const double*, unsigned, const emp_node_t*);
  typedef double (*emp_sampling_prob_func)(void**, const double*, unsigned, const emp_node_t*);
  typedef double (*emp_intensity_func)(void**, const double*, unsigned, const emp_node_t*);
  typedef double (*emp_loglik_func)(void**, const double*, unsigned, const emp_node_t*);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
