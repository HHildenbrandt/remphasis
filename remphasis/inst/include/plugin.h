/***********************************************************/
/* C-API for plug-in diversification models                */
/***********************************************************/


#ifndef EMPHASIS_PLUGIN_H_INCLUDED
#define EMPHASIS_PLUGIN_H_INCLUDED

#include <stddef.h>


#if (defined(_WIN32) || defined(__WIN32__))
# if defined(__GNUC__)
#   define EMP_STDCALL __attribute__((stdcall))
# elif defined(_MSC_VER) || defined(_ICC) || defined(_STDCALL_SUPPORTED)
#   define EMP_CALL __cdecl
# else
#   define EMP_CALL
# endif
#else
# define EMP_CALL
#endif

#if (defined(_WIN32) || defined(__WIN32__)) && !defined(__LCC__)
# define EMP_EXTERN(T) extern "C" __declspec(dllexport) T EMP_CALL
#else
# define EMP_EXTERN(T) extern "C" T EMP_CALL
#endif


/* tree node */
struct emp_node_t
{
  double brts;
  double n;         /* n[i] = number of species in [time_i-1, time_i) */
  double t_ext;     /* emp_t_ext_tip for present-day species;  emp_t_ext_extinct for extinction nodes */
  double pd;
};


#define emp_t_ext_tip 10e10;     /* t_ext for present nodes */
#define emp_t_ext_extinct 0.0;   /* t_ext for extinction nodes */


typedef const char* (*emp_description_func)();
typedef bool (*emp_is_threadsafe_func)();
typedef int (*emp_nparams_func)();


/* optional per-tree state handling */
typedef void (*emp_free_state_func)(void**);
typedef void (*emp_invalidate_state_func)(void**);


/* diversification model */
typedef double (*emp_extinction_time_func)(void**, double, const double*, unsigned, const emp_node_t*);
typedef double (*emp_nh_rate_func)(void**, double, const double*, unsigned, const emp_node_t*);
typedef double (*emp_sampling_prob_func)(void**, const double*, unsigned, const emp_node_t*);
typedef double (*emp_loglik_func)(void**, const double*, unsigned, const emp_node_t*);


/* optional hints for optimizer */
typedef void (*emp_lower_bound_func)(double*);
typedef void (*emp_upper_bound_func)(double*);


#endif
