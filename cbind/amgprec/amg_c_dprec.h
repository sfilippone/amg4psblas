#ifndef AMG_C_DPREC_
#define AMG_C_DPREC_

#include "amg_config.h"
#include "psb_base_cbind.h"
#include "psb_prec_cbind.h"
#include "psb_linsolve_cbind.h"
#include "amg_c_prec_info.h"

/* Object handle related routines */
/* Note:  amg_get_XXX_handle returns:  <= 0  unsuccessful */
/*                                     >0    valid handle */
#ifdef __cplusplus
extern "C" {
#endif
  typedef struct AMG_C_DPREC {
    void *dprec;
  } amg_c_dprec; 
  
  amg_c_dprec* amg_c_dprec_new();
  psb_i_t amg_c_dprec_delete(amg_c_dprec* p);
 
  psb_i_t amg_c_dprecinit(psb_c_ctxt cctxt, amg_c_dprec *ph, const char *ptype);
  psb_i_t amg_c_dprecseti(amg_c_dprec *ph, const char *what, psb_i_t val);
  psb_i_t amg_c_dprecsetc(amg_c_dprec *ph, const char *what, const char *val);
  psb_i_t amg_c_dprecsetr(amg_c_dprec *ph, const char *what, double val);
  /* Optional arguments mirror the native set(ilev,ilmax,pos,idx) arguments.
   * NULL means absent, not zero. Levels and indexed native controls are
   * one-based; level 1 is the finest. Explicit ranges must be valid.
   * what, string val and non-NULL pos are NUL-terminated (at most 80 chars).
   * Set hierarchy options before hierarchy_build; set per-level smoother
   * options before smoothers_build. These calls do not rebuild anything.
   */
  psb_i_t amg_c_dprecseti_opt(amg_c_dprec *ph, const char *what, psb_i_t val,
                                const psb_i_t *ilev, const psb_i_t *ilmax,
                                const char *pos, const psb_i_t *idx);
  psb_i_t amg_c_dprecsetr_opt(amg_c_dprec *ph, const char *what, double val,
                                const psb_i_t *ilev, const psb_i_t *ilmax,
                                const char *pos, const psb_i_t *idx);
  psb_i_t amg_c_dprecsetc_opt(amg_c_dprec *ph, const char *what, const char *val,
                                const psb_i_t *ilev, const psb_i_t *ilmax,
                                const char *pos, const psb_i_t *idx);
  /* Local queries; ilev is one-based. See amg_c_prec_info.h. */
  psb_i_t amg_c_dprec_get_info(amg_c_dprec *ph, amg_c_prec_info *out);
  psb_i_t amg_c_dprec_get_level_info(amg_c_dprec *ph, psb_i_t ilev,
                                      amg_c_prec_level_info *out);
  psb_i_t amg_c_dprecbld(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_dprec *ph);
  psb_i_t amg_c_dhierarchy_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_dprec *ph);
  psb_i_t amg_c_dsmoothers_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_dprec *ph);
  psb_i_t amg_c_dsmoothers_build_opt(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_dprec *ph, const char *afmt, const char *chfmt);
  psb_i_t amg_c_dprecapply(amg_c_dprec *ph, psb_c_dvector *bh, psb_c_dvector *xh, psb_c_descriptor *cdh);
  psb_i_t amg_c_dprecapply_opt(amg_c_dprec *ph, psb_c_dvector *bh, psb_c_dvector *xh, psb_c_descriptor *cdh, const char *ctrans);
  psb_i_t amg_c_dprecfree(amg_c_dprec *ph);
  psb_i_t amg_c_dprecbld_opt(psb_c_dspmat *ah, psb_c_descriptor *cdh, 
			  amg_c_dprec *ph, const char *afmt);
  psb_i_t amg_c_ddescr(amg_c_dprec *ph);
  psb_i_t amg_c_dallocate_wrk(amg_c_dprec *ph, const char *chfmt);

  psb_i_t amg_c_dkrylov(const char *method, psb_c_dspmat *ah, amg_c_dprec *ph, 
		  psb_c_dvector *bh, psb_c_dvector *xh,
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt);


#ifdef __cplusplus
}
#endif

#endif
