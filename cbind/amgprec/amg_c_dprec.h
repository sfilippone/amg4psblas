#ifndef AMG_C_DPREC_
#define AMG_C_DPREC_

#include "amg_const.h"
#include "psb_base_cbind.h"
#include "psb_prec_cbind.h"
#include "psb_krylov_cbind.h"

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
 
  psb_i_t amg_c_dprecinit(psb_i_t ictxt, amg_c_dprec *ph, const char *ptype);
  psb_i_t amg_c_dprecseti(amg_c_dprec *ph, const char *what, psb_i_t val);
  psb_i_t amg_c_dprecsetc(amg_c_dprec *ph, const char *what, const char *val);
  psb_i_t amg_c_dprecsetr(amg_c_dprec *ph, const char *what, double val);
  psb_i_t amg_c_dprecbld(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_dprec *ph);
  psb_i_t amg_c_dhierarchy_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_dprec *ph);
  psb_i_t amg_c_dsmoothers_build(psb_c_dspmat *ah, psb_c_descriptor *cdh, amg_c_dprec *ph);
  psb_i_t amg_c_dprecfree(amg_c_dprec *ph);
  psb_i_t amg_c_dprecbld_opt(psb_c_dspmat *ah, psb_c_descriptor *cdh, 
			  amg_c_dprec *ph, const char *afmt);
  psb_i_t amg_c_ddescr(amg_c_dprec *ph);

  psb_i_t amg_c_dkrylov(const char *method, psb_c_dspmat *ah, amg_c_dprec *ph, 
		  psb_c_dvector *bh, psb_c_dvector *xh,
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt);


#ifdef __cplusplus
}
#endif

#endif
