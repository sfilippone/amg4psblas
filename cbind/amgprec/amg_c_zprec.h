#ifndef AMG_C_ZPREC_
#define AMG_C_ZPREC_

#include "amg_config.h"
#include "psb_base_cbind.h"
#include "psb_prec_cbind.h"
#include "psb_linsolve_cbind.h"
#include "amg_c_prec_info.h"

/* Object handle related routines */
/* Note:  amg_get_XXX_handle returns:  <= 0  unsuccessful */
/*                                     >0    valid handle */
#ifdef __cplusplus
extern "C"
{
#endif
  typedef struct AMG_C_ZPREC
  {
    void *dprec;
  } amg_c_zprec;

  amg_c_zprec *amg_c_zprec_new();
  psb_i_t amg_c_zprec_delete(amg_c_zprec *p);

  psb_i_t amg_c_zprecinit(psb_c_ctxt cctxt, amg_c_zprec *ph, const char *ptype);
  psb_i_t amg_c_zprecseti(amg_c_zprec *ph, const char *what, psb_i_t val);
  psb_i_t amg_c_zprecsetc(amg_c_zprec *ph, const char *what, const char *val);
  psb_i_t amg_c_zprecsetr(amg_c_zprec *ph, const char *what, double val);
  /* Optional arguments mirror the native set(ilev,ilmax,pos,idx) arguments.
   * NULL means absent, not zero. Levels and indexed native controls are
   * one-based; level 1 is the finest. Explicit ranges must be valid.
   * what, string val and non-NULL pos are NUL-terminated (at most 80 chars).
   * Set hierarchy options before hierarchy_build; set per-level smoother
   * options before smoothers_build. These calls do not rebuild anything.
   */
  psb_i_t amg_c_zprecseti_opt(amg_c_zprec *ph, const char *what, psb_i_t val,
                                const psb_i_t *ilev, const psb_i_t *ilmax,
                                const char *pos, const psb_i_t *idx);
  psb_i_t amg_c_zprecsetr_opt(amg_c_zprec *ph, const char *what, double val,
                                const psb_i_t *ilev, const psb_i_t *ilmax,
                                const char *pos, const psb_i_t *idx);
  psb_i_t amg_c_zprecsetc_opt(amg_c_zprec *ph, const char *what, const char *val,
                                const psb_i_t *ilev, const psb_i_t *ilmax,
                                const char *pos, const psb_i_t *idx);
  /* Local queries; ilev is one-based. See amg_c_prec_info.h. */
  psb_i_t amg_c_zprec_get_info(amg_c_zprec *ph, amg_c_prec_info *out);
  psb_i_t amg_c_zprec_get_level_info(amg_c_zprec *ph, psb_i_t ilev,
                                      amg_c_prec_level_info *out);
  psb_i_t amg_c_zprecbld(psb_c_zspmat *ah, psb_c_descriptor *cdh, amg_c_zprec *ph);
  psb_i_t amg_c_zhierarchy_build(psb_c_zspmat *ah, psb_c_descriptor *cdh, amg_c_zprec *ph);
  psb_i_t amg_c_zsmoothers_build(psb_c_zspmat *ah, psb_c_descriptor *cdh, amg_c_zprec *ph);
  psb_i_t amg_c_zsmoothers_build_opt(psb_c_zspmat *ah, psb_c_descriptor *cdh, amg_c_zprec *ph, const char *afmt, const char *chfmt);
  psb_i_t amg_c_zprecapply(amg_c_zprec *ph, psb_c_zvector *bh, psb_c_zvector *xh, psb_c_descriptor *cdh);
  psb_i_t amg_c_zprecapply_opt(amg_c_zprec *ph, psb_c_zvector *bh, psb_c_zvector *xh, psb_c_descriptor *cdh, const char *ctrans);
  psb_i_t amg_c_zprecfree(amg_c_zprec *ph);
  psb_i_t amg_c_zprecbld_opt(psb_c_zspmat *ah, psb_c_descriptor *cdh,
                             amg_c_zprec *ph, const char *afmt);

  psb_i_t amg_c_zdescr(amg_c_zprec *ph);
  psb_i_t amg_c_zallocate_wrk(amg_c_zprec *ph, const char *chfmt);

  psb_i_t amg_c_zkrylov(const char *method, psb_c_zspmat *ah, amg_c_zprec *ph,
                        psb_c_zvector *bh, psb_c_zvector *xh,
                        psb_c_descriptor *cdh, psb_c_SolverOptions *opt);

#ifdef __cplusplus
}
#endif

#endif
