/* MPI regression for the C optional setters, resolved queries.
 * Uses a multilevel SPD/Hermitian tridiagonal example. No optional direct
 * solver package is needed: all levels use fixed L1-Jacobi sweeps.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "psb_base_cbind.h"
#include "amg_cbind.h"

static psb_c_ctxt context;
static psb_i_t rank_id;

static void fail(const char *what, int line)
{
  fprintf(stderr, "rank %lld: %s at line %d\n", (long long)rank_id, what, line);
  psb_c_print_errmsg();
  psb_c_abort(context);
  exit(EXIT_FAILURE);
}

#define REQUIRE(c) do { if (!(c)) fail(#c, __LINE__); } while (0)
#define CHECK(c) do { psb_i_t status_ = (c); if (status_ != 0) fail(#c, __LINE__); } while (0)

static void exercise_d(psb_c_dspmat *matrix, psb_c_descriptor *desc,
                       psb_i_t nlocal, psb_l_t first)
{
  amg_c_dprec *prec = amg_c_dprec_new();
  amg_c_prec_info info;
  amg_c_prec_level_info level;
  psb_c_dvector x = {0}, y = {0}, sum = {0};
  psb_c_dvector px = {0}, py = {0}, work = {0};
  psb_l_t *rows = (psb_l_t *)malloc((size_t)nlocal * sizeof(*rows));
  double *vx = (double *)malloc((size_t)nlocal * sizeof(*vx));
  double *vy = (double *)malloc((size_t)nlocal * sizeof(*vy));
  psb_i_t i, first_level, last_level, bad_level, idx;
  double xpy, ypx, xpx;
  double scale;

  REQUIRE(prec != NULL);
  REQUIRE(rows != NULL && vx != NULL && vy != NULL);
  for (i = 0; i < nlocal; ++i) {
    double g = (double)(first + i + 1);
    rows[i] = first + i;
    vx[i] = sin(0.1*g) + 0.5;
    vy[i] = cos(0.07*g) - 0.25;
  }
  CHECK(psb_c_dgeall(&x, desc));
  CHECK(psb_c_dgeall(&y, desc));
  CHECK(psb_c_dgeall(&sum, desc));
  CHECK(psb_c_dgeall(&px, desc));
  CHECK(psb_c_dgeall(&py, desc));
  CHECK(psb_c_dgeall(&work, desc));
  CHECK(psb_c_dgeins(nlocal, rows, vx, &x, desc));
  CHECK(psb_c_dgeins(nlocal, rows, vy, &y, desc));
  CHECK(psb_c_dgeasb(&x, desc));
  CHECK(psb_c_dgeasb(&y, desc));
  CHECK(psb_c_dgeasb(&sum, desc));
  CHECK(psb_c_dgeasb(&px, desc));
  CHECK(psb_c_dgeasb(&py, desc));
  CHECK(psb_c_dgeasb(&work, desc));
  free(rows);
  free(vx);
  free(vy);

  REQUIRE(amg_c_dprec_get_info(prec, &info) != 0);
  REQUIRE(info.nlevels == 0);
  REQUIRE(amg_c_dprec_get_info(NULL, &info) != 0);
  REQUIRE(amg_c_dprec_get_info(prec, NULL) != 0);

  {
    CHECK(amg_c_dprecinit(context, prec, "ML"));
    CHECK(amg_c_dprec_get_info(prec, &info));
    REQUIRE(info.nlevels > 0 && info.hierarchy_built == 0);
    REQUIRE(info.operator_complexity == 0);
    bad_level = 0;
    REQUIRE(amg_c_dprecseti_opt(prec, "SMOOTHER_SWEEPS", 1,
                                   &bad_level, NULL, NULL, NULL) != 0);
    REQUIRE(amg_c_dprecseti_opt(NULL, "SMOOTHER_SWEEPS", 1,
                                   NULL, NULL, NULL, NULL) != 0);
    REQUIRE(amg_c_dprecsetc_opt(prec, NULL, "VCYCLE",
                                   NULL, NULL, NULL, NULL) != 0);

    /* NULL options preserve the ordinary whole-hierarchy setter behavior. */
    CHECK(amg_c_dprecsetc_opt(prec, "ML_CYCLE", "VCYCLE", NULL, NULL, NULL, NULL));
    CHECK(amg_c_dprecsetc(prec, "PAR_AGGR_ALG", "DECOUPLED"));
    CHECK(amg_c_dprecsetc(prec, "AGGR_TYPE", "SOC1"));
    CHECK(amg_c_dprecsetc(prec, "AGGR_PROL", "SMOOTHED"));
    CHECK(amg_c_dprecsetc(prec, "SMOOTHER_TYPE", "L1-JACOBI"));
    CHECK(amg_c_dprecsetc(prec, "COARSE_SOLVE", "L1-JACOBI"));
    CHECK(amg_c_dprecseti_opt(prec, "MIN_COARSE_SIZE", 8, NULL, NULL, NULL, NULL));
    CHECK(amg_c_dprecseti(prec, "MIN_COARSE_SIZE_PER_PROCESS", 1));
    CHECK(amg_c_dprecseti(prec, "OUTER_SWEEPS", 1));
    CHECK(amg_c_dprecsetr_opt(prec, "AGGR_THRESH", 0.0, NULL, NULL, NULL, NULL));
    CHECK(amg_c_dhierarchy_build(matrix, desc, prec));
    CHECK(amg_c_dprec_get_info(prec, &info));
    REQUIRE(info.hierarchy_built == 1 && info.nlevels > 1);
    REQUIRE(isfinite(info.operator_complexity) && info.operator_complexity >= 1.0);
    REQUIRE(info.local_bytes > 0 && info.outer_sweeps == 1);
    first_level = 1;
    last_level = info.nlevels;
    bad_level = last_level + 1;
    REQUIRE(amg_c_dprec_get_level_info(prec, bad_level, &level) != 0);
    REQUIRE(amg_c_dprecseti_opt(prec, "SMOOTHER_SWEEPS", 2,
                                   &last_level, &first_level, NULL, NULL) != 0);

    /* Configure the actual built levels, including the coarsest one.
     * Explicit POST is separately allocated; its report must not alias PRE.
     */
    CHECK(amg_c_dprecsetc_opt(prec, "SMOOTHER_TYPE", "L1-JACOBI",
                                 &first_level, &last_level, "PRE", NULL));
    CHECK(amg_c_dprecsetc_opt(prec, "SMOOTHER_TYPE", "L1-JACOBI",
                                 &first_level, &last_level, "POST", NULL));
    CHECK(amg_c_dprecseti_opt(prec, "SMOOTHER_SWEEPS", 2,
                                 &first_level, &last_level, "PRE", NULL));
    CHECK(amg_c_dprecseti_opt(prec, "SMOOTHER_SWEEPS", 3,
                                 &first_level, &last_level, "POST", NULL));
    CHECK(amg_c_dprec_get_level_info(prec, 1, &level));
    REQUIRE(level.sweeps_pre == 2 && level.sweeps_post == 3);
    CHECK(amg_c_dprecseti_opt(prec, "SMOOTHER_SWEEPS", 2,
                                 &first_level, &last_level, "POST", NULL));
    CHECK(amg_c_dprecsetc_opt(prec, "SMOOTHER_STOP", "FALSE",
                                 &first_level, &last_level, "PRE", NULL));
    CHECK(amg_c_dprecsetc_opt(prec, "SMOOTHER_STOP", "FALSE",
                                 &first_level, &last_level, "POST", NULL));
    /* Exercise non-NULL idx forwarding on a named scalar parameter; a
     * parameter which does not consume an index retains native behavior.
     */
    idx = 1;
    CHECK(amg_c_dprecsetr_opt(prec, "SMOOTHER_STOPTOL", 1.e-12,
                                 &first_level, &last_level, "PRE", &idx));
    CHECK(amg_c_dsmoothers_build(matrix, desc, prec));

    for (i = 1; i <= last_level; ++i) {
      CHECK(amg_c_dprec_get_level_info(prec, i, &level));
      REQUIRE(level.sweeps_pre == 2 && level.sweeps_post == 2);
      REQUIRE(level.residual_stop_pre == 0 && level.residual_stop_post == 0);
      REQUIRE(strstr(level.solver_pre, "L1") != NULL);
      REQUIRE(strstr(level.solver_post, "L1") != NULL);
    }

    CHECK(amg_c_dprecapply(prec, &x, &px, desc));
    CHECK(amg_c_dprecapply(prec, &y, &py, desc));
    xpx = psb_c_dgedot(&x, &px, desc);
    xpy = psb_c_dgedot(&x, &py, desc);
    ypx = psb_c_dgedot(&y, &px, desc);
    REQUIRE(isfinite(xpx) && xpx > 0.0);

    scale = fmax(1.0, fmax(fabs(xpy), fabs(ypx)));
    REQUIRE(fabs(xpy - ypx) < 5.e-11 * scale);

    CHECK(psb_c_dgeaxpbyz(1.0, &x, 1.0, &y, &sum, desc));
    CHECK(amg_c_dprecapply(prec, &sum, &work, desc));
    CHECK(psb_c_dgeaxpby(-1.0, &px, 1.0, &work, desc));
    CHECK(psb_c_dgeaxpby(-1.0, &py, 1.0, &work, desc));
    scale = fmax(1.0, psb_c_dgenrm2(&px, desc) + psb_c_dgenrm2(&py, desc));
    REQUIRE(psb_c_dgenrm2(&work, desc) < 5.e-11 * scale);

    CHECK(psb_c_dvect_set_scal(&work, 19.0));
    CHECK(amg_c_dprecapply_opt(prec, &x, &work, desc, "N"));
    CHECK(psb_c_dgeaxpby(-1.0, &px, 1.0, &work, desc));
    REQUIRE(psb_c_dgenrm2(&work, desc) < 5.e-11 * scale);
    CHECK(psb_c_dvect_set_scal(&sum, 0.0));
    CHECK(psb_c_dvect_set_scal(&work, 19.0));
    CHECK(amg_c_dprecapply(prec, &sum, &work, desc));
    REQUIRE(psb_c_dgenrm2(&work, desc) < 5.e-13);

    CHECK(amg_c_dprecfree(prec));
    CHECK(amg_c_dprec_get_info(prec, &info));
    REQUIRE(info.nlevels == 0 && info.hierarchy_built == 0);
    REQUIRE(amg_c_dprec_get_level_info(prec, 1, &level) != 0);
  }

  CHECK(amg_c_dprec_delete(prec));
  CHECK(psb_c_dgefree(&x, desc));
  CHECK(psb_c_dgefree(&y, desc));
  CHECK(psb_c_dgefree(&sum, desc));
  CHECK(psb_c_dgefree(&px, desc));
  CHECK(psb_c_dgefree(&py, desc));
  CHECK(psb_c_dgefree(&work, desc));
  if (rank_id == 0) puts("d: optional setters, resolved queries, multilevel action passed");
}

static void exercise_z(psb_c_zspmat *matrix, psb_c_descriptor *desc,
                       psb_i_t nlocal, psb_l_t first)
{
  amg_c_zprec *prec = amg_c_zprec_new();
  amg_c_prec_info info;
  amg_c_prec_level_info level;
  psb_c_zvector x = {0}, y = {0}, sum = {0};
  psb_c_zvector px = {0}, py = {0}, work = {0};
  psb_l_t *rows = (psb_l_t *)malloc((size_t)nlocal * sizeof(*rows));
  psb_z_t *vx = (psb_z_t *)malloc((size_t)nlocal * sizeof(*vx));
  psb_z_t *vy = (psb_z_t *)malloc((size_t)nlocal * sizeof(*vy));
  psb_i_t i, first_level, last_level, bad_level, idx;
  psb_z_t xpy, ypx, xpx;
  double scale;

  REQUIRE(prec != NULL);
  REQUIRE(rows != NULL && vx != NULL && vy != NULL);
  for (i = 0; i < nlocal; ++i) {
    double g = (double)(first + i + 1);
    rows[i] = first + i;
    vx[i] = sin(0.1*g) + 0.5 + I*cos(0.13*g);
    vy[i] = cos(0.07*g) - 0.25 + I*sin(0.17*g);
  }
  CHECK(psb_c_zgeall(&x, desc));
  CHECK(psb_c_zgeall(&y, desc));
  CHECK(psb_c_zgeall(&sum, desc));
  CHECK(psb_c_zgeall(&px, desc));
  CHECK(psb_c_zgeall(&py, desc));
  CHECK(psb_c_zgeall(&work, desc));
  CHECK(psb_c_zgeins(nlocal, rows, vx, &x, desc));
  CHECK(psb_c_zgeins(nlocal, rows, vy, &y, desc));
  CHECK(psb_c_zgeasb(&x, desc));
  CHECK(psb_c_zgeasb(&y, desc));
  CHECK(psb_c_zgeasb(&sum, desc));
  CHECK(psb_c_zgeasb(&px, desc));
  CHECK(psb_c_zgeasb(&py, desc));
  CHECK(psb_c_zgeasb(&work, desc));
  free(rows);
  free(vx);
  free(vy);

  REQUIRE(amg_c_zprec_get_info(prec, &info) != 0);
  REQUIRE(info.nlevels == 0);
  REQUIRE(amg_c_zprec_get_info(NULL, &info) != 0);
  REQUIRE(amg_c_zprec_get_info(prec, NULL) != 0);

  {
    CHECK(amg_c_zprecinit(context, prec, "ML"));
    CHECK(amg_c_zprec_get_info(prec, &info));
    REQUIRE(info.nlevels > 0 && info.hierarchy_built == 0);
    REQUIRE(info.operator_complexity == 0);
    bad_level = 0;
    REQUIRE(amg_c_zprecseti_opt(prec, "SMOOTHER_SWEEPS", 1,
                                   &bad_level, NULL, NULL, NULL) != 0);
    REQUIRE(amg_c_zprecseti_opt(NULL, "SMOOTHER_SWEEPS", 1,
                                   NULL, NULL, NULL, NULL) != 0);
    REQUIRE(amg_c_zprecsetc_opt(prec, NULL, "VCYCLE",
                                   NULL, NULL, NULL, NULL) != 0);

    /* NULL options preserve the ordinary whole-hierarchy setter behavior. */
    CHECK(amg_c_zprecsetc_opt(prec, "ML_CYCLE", "VCYCLE", NULL, NULL, NULL, NULL));
    CHECK(amg_c_zprecsetc(prec, "PAR_AGGR_ALG", "DECOUPLED"));
    CHECK(amg_c_zprecsetc(prec, "AGGR_TYPE", "SOC1"));
    CHECK(amg_c_zprecsetc(prec, "AGGR_PROL", "SMOOTHED"));
    CHECK(amg_c_zprecsetc(prec, "SMOOTHER_TYPE", "L1-JACOBI"));
    CHECK(amg_c_zprecsetc(prec, "COARSE_SOLVE", "L1-JACOBI"));
    CHECK(amg_c_zprecseti_opt(prec, "MIN_COARSE_SIZE", 8, NULL, NULL, NULL, NULL));
    CHECK(amg_c_zprecseti(prec, "MIN_COARSE_SIZE_PER_PROCESS", 1));
    CHECK(amg_c_zprecseti(prec, "OUTER_SWEEPS", 1));
    CHECK(amg_c_zprecsetr_opt(prec, "AGGR_THRESH", 0.0, NULL, NULL, NULL, NULL));
    CHECK(amg_c_zhierarchy_build(matrix, desc, prec));
    CHECK(amg_c_zprec_get_info(prec, &info));
    REQUIRE(info.hierarchy_built == 1 && info.nlevels > 1);
    REQUIRE(isfinite(info.operator_complexity) && info.operator_complexity >= 1.0);
    REQUIRE(info.local_bytes > 0 && info.outer_sweeps == 1);
    first_level = 1;
    last_level = info.nlevels;
    bad_level = last_level + 1;
    REQUIRE(amg_c_zprec_get_level_info(prec, bad_level, &level) != 0);
    REQUIRE(amg_c_zprecseti_opt(prec, "SMOOTHER_SWEEPS", 2,
                                   &last_level, &first_level, NULL, NULL) != 0);

    /* Configure the actual built levels, including the coarsest one.
     * Explicit POST is separately allocated; its report must not alias PRE.
     */
    CHECK(amg_c_zprecsetc_opt(prec, "SMOOTHER_TYPE", "L1-JACOBI",
                                 &first_level, &last_level, "PRE", NULL));
    CHECK(amg_c_zprecsetc_opt(prec, "SMOOTHER_TYPE", "L1-JACOBI",
                                 &first_level, &last_level, "POST", NULL));
    CHECK(amg_c_zprecseti_opt(prec, "SMOOTHER_SWEEPS", 2,
                                 &first_level, &last_level, "PRE", NULL));
    CHECK(amg_c_zprecseti_opt(prec, "SMOOTHER_SWEEPS", 3,
                                 &first_level, &last_level, "POST", NULL));
    CHECK(amg_c_zprec_get_level_info(prec, 1, &level));
    REQUIRE(level.sweeps_pre == 2 && level.sweeps_post == 3);
    CHECK(amg_c_zprecseti_opt(prec, "SMOOTHER_SWEEPS", 2,
                                 &first_level, &last_level, "POST", NULL));
    CHECK(amg_c_zprecsetc_opt(prec, "SMOOTHER_STOP", "FALSE",
                                 &first_level, &last_level, "PRE", NULL));
    CHECK(amg_c_zprecsetc_opt(prec, "SMOOTHER_STOP", "FALSE",
                                 &first_level, &last_level, "POST", NULL));
    /* Exercise non-NULL idx forwarding on a named scalar parameter; a
     * parameter which does not consume an index retains native behavior.
     */
    idx = 1;
    CHECK(amg_c_zprecsetr_opt(prec, "SMOOTHER_STOPTOL", 1.e-12,
                                 &first_level, &last_level, "PRE", &idx));
    CHECK(amg_c_zsmoothers_build(matrix, desc, prec));

    for (i = 1; i <= last_level; ++i) {
      CHECK(amg_c_zprec_get_level_info(prec, i, &level));
      REQUIRE(level.sweeps_pre == 2 && level.sweeps_post == 2);
      REQUIRE(level.residual_stop_pre == 0 && level.residual_stop_post == 0);
      REQUIRE(strstr(level.solver_pre, "L1") != NULL);
      REQUIRE(strstr(level.solver_post, "L1") != NULL);
    }

    CHECK(amg_c_zprecapply(prec, &x, &px, desc));
    CHECK(amg_c_zprecapply(prec, &y, &py, desc));
    xpx = psb_c_zgedot(&x, &px, desc);
    xpy = psb_c_zgedot(&x, &py, desc);
    ypx = psb_c_zgedot(&y, &px, desc);
    REQUIRE(isfinite(creal(xpx)) && creal(xpx) > 0.0);
    REQUIRE(fabs(cimag(xpx)) < 5.e-11 * fmax(1.0, cabs(xpx)));
    scale = fmax(1.0, fmax(cabs(xpy), cabs(ypx)));
    REQUIRE(cabs(xpy - conj(ypx)) < 5.e-11 * scale);

    CHECK(psb_c_zgeaxpbyz(1.0, &x, 1.0, &y, &sum, desc));
    CHECK(amg_c_zprecapply(prec, &sum, &work, desc));
    CHECK(psb_c_zgeaxpby(-1.0, &px, 1.0, &work, desc));
    CHECK(psb_c_zgeaxpby(-1.0, &py, 1.0, &work, desc));
    scale = fmax(1.0, psb_c_zgenrm2(&px, desc) + psb_c_zgenrm2(&py, desc));
    REQUIRE(psb_c_zgenrm2(&work, desc) < 5.e-11 * scale);

    CHECK(psb_c_zvect_set_scal(&work, 19.0));
    CHECK(amg_c_zprecapply_opt(prec, &x, &work, desc, "N"));
    CHECK(psb_c_zgeaxpby(-1.0, &px, 1.0, &work, desc));
    REQUIRE(psb_c_zgenrm2(&work, desc) < 5.e-11 * scale);
    CHECK(psb_c_zvect_set_scal(&sum, 0.0));
    CHECK(psb_c_zvect_set_scal(&work, 19.0));
    CHECK(amg_c_zprecapply(prec, &sum, &work, desc));
    REQUIRE(psb_c_zgenrm2(&work, desc) < 5.e-13);

    CHECK(amg_c_zprecfree(prec));
    CHECK(amg_c_zprec_get_info(prec, &info));
    REQUIRE(info.nlevels == 0 && info.hierarchy_built == 0);
    REQUIRE(amg_c_zprec_get_level_info(prec, 1, &level) != 0);
  }

  CHECK(amg_c_zprec_delete(prec));
  CHECK(psb_c_zgefree(&x, desc));
  CHECK(psb_c_zgefree(&y, desc));
  CHECK(psb_c_zgefree(&sum, desc));
  CHECK(psb_c_zgefree(&px, desc));
  CHECK(psb_c_zgefree(&py, desc));
  CHECK(psb_c_zgefree(&work, desc));
  if (rank_id == 0) puts("z: optional setters, resolved queries, multilevel action passed");
}

int main(void)
{
  psb_c_descriptor desc = {0};
  psb_c_dspmat da = {0};
  psb_c_zspmat za = {0};
  amg_c_dprec *dp;
  amg_c_zprec *zp;
  psb_i_t nprocs, nlocal, j;
  const psb_l_t n = 256;
  psb_l_t first, last, row, rows[3], cols[3];
  double dv[3];
  psb_z_t zv[3];

  psb_c_init(&context);
  psb_c_info(context, &rank_id, &nprocs);
  REQUIRE(nprocs > 0 && nprocs <= 16);
  psb_c_set_index_base(0);
  psb_c_seterraction_ret();
  first = (n * rank_id) / nprocs;
  last = (n * (rank_id + 1)) / nprocs;
  nlocal = (psb_i_t)(last - first);
  CHECK(psb_c_cdall_nl(nlocal, context, &desc));
  CHECK(psb_c_dspall(&da, &desc));
  CHECK(psb_c_zspall(&za, &desc));
  for (row = first; row < last; ++row) {
    j = 0;
    if (row > 0) {
      rows[j] = row; cols[j] = row - 1; dv[j] = -1.0; zv[j] = -1.0; ++j;
    }
    rows[j] = row; cols[j] = row; dv[j] = 2.5; zv[j] = 2.5; ++j;
    if (row + 1 < n) {
      rows[j] = row; cols[j] = row + 1; dv[j] = -1.0; zv[j] = -1.0; ++j;
    }
    CHECK(psb_c_dspins(j, rows, cols, dv, &da, &desc));
    CHECK(psb_c_zspins(j, rows, cols, zv, &za, &desc));
  }
  CHECK(psb_c_cdasb(&desc));
  CHECK(psb_c_dspasb(&da, &desc));
  CHECK(psb_c_zspasb(&za, &desc));

  exercise_d(&da, &desc, nlocal, first);
  exercise_z(&za, &desc, nlocal, first);

  /* Heap constructors and deletion functions. */
  dp = amg_c_dprec_new();
  REQUIRE(dp != NULL);
  CHECK(amg_c_dprecinit(context, dp, "ML"));
  CHECK(amg_c_dprec_delete(dp));
  zp = amg_c_zprec_new();
  REQUIRE(zp != NULL);
  CHECK(amg_c_zprecinit(context, zp, "ML"));
  CHECK(amg_c_zprec_delete(zp));

  CHECK(psb_c_dspfree(&da, &desc));
  CHECK(psb_c_zspfree(&za, &desc));
  CHECK(psb_c_cdfree(&desc));
  if (rank_id == 0) puts("AMG C extensions: PASS");
  psb_c_exit(context);
  return EXIT_SUCCESS;
}
