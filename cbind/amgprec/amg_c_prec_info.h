#ifndef AMG_C_PREC_INFO_
#define AMG_C_PREC_INFO_

#include <stdint.h>
#include "psb_config.h"
#include "psb_types.h"

/* Feature test for the optional setters, queries and native-object destroy. */
#define AMG_C_PREC_NAME_LENGTH 80

/* Queries are local, do not perform MPI reductions, and do not build/rebuild
 * the hierarchy. Before hierarchy_build, nlevels is the number of configured
 * level slots, hierarchy_built is zero, and both complexity fields are zero.
 * After hierarchy_build, the complexity fields are the cached native values.
 * local_bytes is the native sizeof(global=.false.) count, excluding C handles.
 */
typedef struct AMG_C_PREC_INFO {
  psb_c_i_t nlevels;
  psb_c_i_t hierarchy_built;
  psb_c_i_t outer_sweeps;
  double operator_complexity;
  double average_coarsening_ratio;
  int64_t local_bytes;
} amg_c_prec_info;

/* The names describe the objects actually installed at this level, including
 * native fallback choices. Strings are NUL-terminated. A missing object is
 * reported as an empty string. residual_stop is 0/1 for Jacobi-family
 * smoothers and -1 when this property is not exposed by that smoother type.
 * This report is not a certificate that the complete AMG action is SPD.
 * cycle and coarse_matrix retain the native integer codes.
 */
typedef struct AMG_C_PREC_LEVEL_INFO {
  psb_c_i_t cycle;
  psb_c_i_t coarse_matrix;
  psb_c_i_t sweeps_pre;
  psb_c_i_t sweeps_post;
  psb_c_i_t residual_stop_pre;
  psb_c_i_t residual_stop_post;
  char smoother_pre[AMG_C_PREC_NAME_LENGTH];
  char solver_pre[AMG_C_PREC_NAME_LENGTH];
  char smoother_post[AMG_C_PREC_NAME_LENGTH];
  char solver_post[AMG_C_PREC_NAME_LENGTH];
} amg_c_prec_level_info;

#endif
