#include "rinit.h"



nlopt_opt(*remp_create)(nlopt_algorithm, unsigned) = NULL;
void(*remp_destroy)(nlopt_opt) = NULL;
nlopt_result(*remp_optimize)(nlopt_opt, double *, double *) = NULL;
nlopt_result(*remp_set_min_objective)(nlopt_opt, nlopt_func, void *) = NULL;
nlopt_result(*remp_set_max_objective)(nlopt_opt, nlopt_func, void *) = NULL;
nlopt_result(*remp_set_lower_bounds)(nlopt_opt, const double *) = NULL;
nlopt_result(*remp_set_lower_bounds1)(nlopt_opt, double) = NULL;
nlopt_result(*remp_set_upper_bounds)(nlopt_opt, const double *) = NULL;
nlopt_result(*remp_set_upper_bounds1)(nlopt_opt, double) = NULL;
nlopt_result(*remp_set_xtol_rel)(nlopt_opt, double) = NULL;
